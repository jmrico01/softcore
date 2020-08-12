#version 450

#define TILE_SIZE 16
#define MAX_MATERIALS 32

#define FLOAT_BIG_POSITIVE 1e6
#define EPSILON 0.000001
#define BOUNCES 4
#define MIN_DIST 0.0001
#define MAX_DIST 20.0

#define QUAT_ONE vec4(0.0, 0.0, 0.0, 1.0)

layout(location = 0) in vec2 inUv;

layout(location = 0) out vec4 outColor;

struct Material
{
	vec3 albedo;
	float smoothness;
	vec3 emissionColor;
	float emission;
};

struct Triangle
{
	vec3 a, b, c;
	vec3 normal;
	uint materialIndex;
};

struct Bvh
{
	vec3 aabbMin;
	uint startTriangle;
	vec3 aabbMax;
	uint endTriangle;
	uint skip;
};

layout(binding = 0) uniform UBO
{
	vec3 cameraPos;
	uint seed;
	vec3 filmTopLeft;
	vec3 filmUnitOffsetX;
	vec3 filmUnitOffsetY;
	Material materials[MAX_MATERIALS];
} ubo;
layout(std140, binding = 1) readonly buffer Triangles
{
	Triangle triangles[];
};
layout(std140, binding = 2) readonly buffer Bvhs
{
	Bvh bvhs[];
};

const uint UINT32_MAX_VALUE = 0xffffffff;

// Reference https://en.wikipedia.org/wiki/Xorshift
uint XOrShift32(uint value)
{
	value ^= value << 13;
	value ^= value >> 17;
	value ^= value << 5;
	return value;
}

uint RandomUInt(inout uint state)
{
	state = XOrShift32(state);
	return state;
}

float RandomUnilateral(inout uint state)
{
	return float(RandomUInt(state)) / float(UINT32_MAX_VALUE);
}

float RandomBilateral(inout uint state)
{
	return 2.0 * RandomUnilateral(state) - 1.0;
}

vec3 Reciprocal(vec3 v)
{
	return vec3(1.0 / v.x, 1.0 / v.y, 1.0 / v.z);
}

vec4 QuatMultiply(vec4 q1, vec4 q2)
{
    vec4 result;
    result.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
    result.y = q1.w*q2.y + q1.y*q2.w + q1.z*q2.x - q1.x*q2.z;
    result.z = q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x;
    result.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
    return result;
}

vec4 QuatInverse(vec4 q)
{
	return vec4(-q.x, -q.y, -q.z, q.w);
}

vec3 QuatRotate(vec3 v, vec4 quat)
{
	// TODO Quat multiply with baked w=0
	const vec4 vQuat = vec4(v, 0.0);

	const vec4 qv = QuatMultiply(quat, vQuat);
	const vec4 qInv = QuatInverse(quat);
	const vec4 qvqInv = QuatMultiply(qv, qInv);

	return vec3(qvqInv.x, qvqInv.y, qvqInv.z);
}

vec4 QuatFromAngleUnitAxis(float angle, vec3 axis)
{
	const float cosHalfAngle = cos(angle / 2.0);
	const float sinHalfAngle = sin(angle / 2.0);

	return vec4(axis.x * sinHalfAngle, axis.y * sinHalfAngle, axis.z * sinHalfAngle, cosHalfAngle);
}

// Returns a quaternion for the rotation required to align unit X vector with unit vector reference
vec4 QuatRotationFromUnitX(vec3 reference)
{
	const float dotRefUnitX = dot(vec3(1.0, 0.0, 0.0), reference);
	if (dotRefUnitX > 0.99999f) {
		return QUAT_ONE;
	}
	else if (dotRefUnitX < -0.99999f) {
		// 180-degree rotation about y axis (perpendicular to x)
		// comptime angle = PI, cosHalfAngle = 0, sinHalfAngle = 1, axis = 0.0, 1.0, 0.0
		return vec4(0.0, 1.0, 0.0, 0.0);
	}

	const float angle = 1.0 + dotRefUnitX;
	const vec3 unitAxis = normalize(cross(vec3(1.0, 0.0, 0.0), reference));
	return QuatFromAngleUnitAxis(angle, unitAxis);
}

// NOTE remember this takes in the INVERSE ray direction!
bool RayAABBIntersection(vec3 rayOrigin, vec3 rayDirInv, vec3 aabbMin, vec3 aabbMax, out float tMin, out float tMax)
{
    tMin = -FLOAT_BIG_POSITIVE;
    tMax = FLOAT_BIG_POSITIVE;

    const float tX1 = (aabbMin.x - rayOrigin.x) * rayDirInv.x;
    const float tX2 = (aabbMax.x - rayOrigin.x) * rayDirInv.x;
    tMin = max(tMin, min(tX1, tX2));
    tMax = min(tMax, max(tX1, tX2));

    const float tY1 = (aabbMin.y - rayOrigin.y) * rayDirInv.y;
    const float tY2 = (aabbMax.y - rayOrigin.y) * rayDirInv.y;
    tMin = max(tMin, min(tY1, tY2));
    tMax = min(tMax, max(tY1, tY2));

    const float tZ1 = (aabbMin.z - rayOrigin.z) * rayDirInv.z;
    const float tZ2 = (aabbMax.z - rayOrigin.z) * rayDirInv.z;
    tMin = max(tMin, min(tZ1, tZ2));
    tMax = min(tMax, max(tZ1, tZ2));

    return tMax >= tMin;
}

bool RayTriangleIntersection(vec3 rayOrigin, vec3 rayDir, vec3 a, vec3 b, vec3 c, out float t)
{
	const vec3 ab = b - a;
	const vec3 ac = c - a;
	const vec3 h = cross(rayDir, ac);
	const float x = dot(ab, h);
	if (x > -EPSILON && x < EPSILON) {
		return false;
	}

	const float f = 1.0 / x;
	const vec3 s = rayOrigin - a;
	const float u = f * dot(s, h);
	if (u < 0.0 || u > 1.0) {
		return false;
	}

	const vec3 q = cross(s, ab);
	const float v = f * dot(rayDir, q);
	if (v < 0.0f || u + v > 1.0) {
		return false;
	}

	t = f * dot(ac, q);
	return true;
}

// TODO review inout hitMaterialIndex and other in/out params
void TraverseBvhs(vec3 rayOrigin, vec3 rayDir, vec3 inverseRayDir,
				  inout uint hitMaterialIndex, out vec3 hitNormal, inout float hitDist)
{
	uint i = 0;
	while (i < bvhs.length()) {
		const Bvh bvh = bvhs[i];

		float tMin, tMax;
		const bool aabbHit = RayAABBIntersection(rayOrigin, inverseRayDir, bvh.aabbMin, bvh.aabbMax, tMin, tMax);
		if (aabbHit && !(tMin < 0.0 && tMax < 0.0) && tMin < hitDist) {
			const uint tStart = bvh.startTriangle;
			const uint tEnd = bvh.endTriangle;
			for (uint i = tStart; i < tEnd; i++) {
				const vec3 normal = triangles[i].normal;
				const float dotNegRayNormal = dot(rayDir, normal);
				if (dotNegRayNormal <= 0.0) {
					float t;
					const bool triangleHit = RayTriangleIntersection(rayOrigin, rayDir,
																	 triangles[i].a, triangles[i].b, triangles[i].c, t);
					if (triangleHit && t > MIN_DIST && t < hitDist) {
						hitMaterialIndex = triangles[i].materialIndex;
						hitNormal = normal;
						hitDist = t;
					}
				}
			}

			i++;
		}
		else {
			i += bvh.skip;
		}
	}
}

#define SAMPLES 1

vec3 RaycastColor(vec3 rayOrigin, vec3 rayDir, uint randomSeed)
{
	uint randomState = randomSeed;
	float intensity = 1.0;
	vec3 color = vec3(0.0, 0.0, 0.0);

	float sampleWeight = 1.0 / SAMPLES;
	for (uint s = 0; s < SAMPLES; s++) {
		vec3 sampleRayOrigin = rayOrigin;
		vec3 sampleRayDir = rayDir;

		for (uint i = 0; i < BOUNCES; i++) {
			const vec3 inverseRayDir = Reciprocal(sampleRayDir);

			uint hitMaterialIndex = MAX_MATERIALS;
			vec3 hitNormal;
			float hitDist = MAX_DIST;

			TraverseBvhs(sampleRayOrigin, sampleRayDir, inverseRayDir, hitMaterialIndex, hitNormal, hitDist);

			if (hitMaterialIndex == MAX_MATERIALS) {
				break;
			}
			else {
				Material hitMaterial = ubo.materials[hitMaterialIndex];
				if (hitMaterial.emission > 0.0) {
					color += hitMaterial.emission * hitMaterial.emissionColor * intensity * sampleWeight;
					break;
				}
				else {
					intensity *= 0.5f;
					sampleRayOrigin += sampleRayDir * hitDist;

					const vec3 pureBounce = reflect(sampleRayDir, hitNormal);
					const vec3 randomHemisphereUnitX = normalize(vec3(
						RandomUnilateral(randomState),
						RandomBilateral(randomState),
						RandomBilateral(randomState)
					));
					const vec4 quatXToNormal = QuatRotationFromUnitX(hitNormal);
					const vec3 randomBounce = QuatRotate(randomHemisphereUnitX, quatXToNormal);

					sampleRayDir = normalize(mix(randomBounce, pureBounce, hitMaterial.smoothness));
				}
			}
		}
	}

	return color;
}

void main()
{
	const vec2 pixel = gl_FragCoord.xy;
	const vec3 filmOffsetX = ubo.filmUnitOffsetX * pixel.x;
	const vec3 filmOffsetY = ubo.filmUnitOffsetY * pixel.y;
	const vec3 filmPos = ubo.filmTopLeft + filmOffsetX + filmOffsetY;
	const vec3 rayDir = normalize(filmPos - ubo.cameraPos);

	uint seed = uint(pixel.x * 53829 + pixel.y * 17);
	seed += ubo.seed;
	//seed = uint(rayDir.x * 1203.0 + rayDir.y * 2383.0 + rayDir.z);
	const vec3 raycastColor = RaycastColor(ubo.cameraPos, rayDir, seed);

	outColor = vec4(raycastColor, 1.0);
}
