#version 450

// Vertex attributes
layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in vec3 inColor;

// Instance attributes
layout(location = 3) in mat4 model; // mat4 takes locations 3, 4, 5, 6
layout(location = 7) in vec3 color;
layout(location = 8) in vec3 collapseMid;
layout(location = 9) in float collapseT;

layout(location = 0) out vec3 outNormal;
layout(location = 1) out vec3 outColor;

layout(binding = 0) uniform UniformBufferObject {
    mat4 view;
    mat4 proj;
} ubo;

// TODO copied in lightmapMesh.frag
vec3 Lerp(vec3 a, vec3 b, float t)
{
	return a * (1.0 - t) + b * t;
}

float Random(vec2 uv)
{
	return fract(sin(dot(uv, vec2(12.9898, 78.233))) * 43758.5453123);
}

float CubicPulse(float xMax, float width, float x)
{
    x = abs(x - xMax);
    if (x > width) {
		return 0.0;
	}

    x /= width;
    return 1.0 - x * x * (3.0 - 2.0 * x);
}

void main()
{
	float randMagMax = 0.15;

	vec4 worldNormal = model * vec4(inNormal, 0.0);
    outNormal = normalize(worldNormal.xyz);
    outColor = inColor * color;

	float randOffset = 0.0;
	if (collapseT > 0.0) {
		vec2 seedUv = vec2(inPosition.x * inPosition.z * collapseT, inPosition.y * inPosition.z * collapseT);
		float rand = Random(seedUv);
		float randMag = max(CubicPulse(collapseT, 0.6, 0.5), 0.0) * randMagMax;
		randOffset = 2.0 * randMag * (rand - 0.5);
	}

	float t = clamp(collapseT + randOffset, 0.0, 1.0);
	t = t * t;
	vec3 collapsedPos = Lerp(inPosition, collapseMid, t);

	vec4 pos = ubo.proj * ubo.view * model * vec4(collapsedPos, 1.0);
    gl_Position = pos;
}
