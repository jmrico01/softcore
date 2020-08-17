#include "raytracer.h"

#include <intrin.h>
#include <time.h>

#include <stb_image_write.h>
#include <Tracy.hpp>
#include <TracyC.h>

#define SAMPLES_DIFFUSE 1
#define SAMPLES_SPECULAR 1
#define BOUNCES_DIFFUSE 1
#define BOUNCES_SPECULAR 1

struct ComputeBvh
{
	Vec3 aabbMin;
	uint32 startTriangle;
	Vec3 aabbMax;
	uint32 endTriangle;
    uint32 skip;
    uint32 pad[3];
};
static_assert(sizeof(ComputeBvh) % 16 == 0); // std140 layout

struct ComputeTriangle
{
    alignas(16) Vec3 a;
    alignas(16) Vec3 b;
    alignas(16) Vec3 c;
    alignas(16) Vec3 normal;
    uint32 materialIndex;
};
static_assert(sizeof(ComputeTriangle) % 16 == 0); // std140 layout

struct RandomSeries
{
    uint32 state;
};

// Reference http://www.burtleburtle.net/bob/hash/integer.html
uint32 WangHash(uint32 seed)
{
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
}

// Reference https://en.wikipedia.org/wiki/Xorshift
uint32 XOrShift32(RandomSeries* series)
{
    uint32 x = series->state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;

    series->state = x;
	return x;
}

uint32 RandomUInt32(RandomSeries* series)
{
    return XOrShift32(series);
}

int RandomInt32(RandomSeries* series, uint32 max)
{
    return (int)(RandomUInt32(series) % max);
}

float32 RandomUnilateral(RandomSeries* series)
{
    return (float32)RandomUInt32(series) / (float32)UINT32_MAX_VALUE;
}

float32 RandomBilateral(RandomSeries* series)
{
    return 2.0f * RandomUnilateral(series) - 1.0f;
}

Vec3 RandomInUnitSphereSurface(RandomSeries* series)
{
	const float32 angle = RandomUnilateral(series) * 2.0f * PI_F;
	const float32 z = RandomBilateral(series);
	const float32 r = Sqrt32(1.0f - z * z);
	return Vec3 { r * Cos32(angle), r * Sin32(angle), z };
}

float32 pow3(float32 x)
{
	return x * x * x;
}

float32 pow4(float32 x)
{
	return x * x * x * x;
}

Vec3 Reflect(Vec3 dir, Vec3 normal)
{
    return dir - 2.0f * Dot(dir, normal) * normal;
}

// Returns a quaternion for the rotation required to align unit X vector with unit vector reference
Quat QuatRotationFromUnitX(Vec3 reference)
{
    const float32 EPSILON = 0.000001f;

	const float dotUnitX = Dot(Vec3::unitX, reference);
	if (dotUnitX > (1.0 - EPSILON)) {
		return Quat::one;
	}
	else if (dotUnitX < (-1.0 + EPSILON)) {
		// 180-degree rotation about y axis (perpendicular to x)
		// comptime angle = PI, cosHalfAngle = 0, sinHalfAngle = 1, axis = 0.0, 1.0, 0.0
		return Quat { 0.0, 1.0, 0.0, 0.0 };
	}

	const Vec3 axis = Cross(Vec3::unitX, reference);
    const Quat q = { .x = axis.x, .y = axis.y, .z = axis.z, .w = 1.0f + dotUnitX };
	return Normalize(q);
}

Vec3 GetNewRayDirDiffuse(Vec3 hitNormal, RandomSeries* series)
{
	return Normalize(hitNormal + RandomInUnitSphereSurface(series));
}

Vec3 GetNewRayDirSpecular(Vec3 rayDir, Vec3 hitNormal, float32 smoothness, RandomSeries* series)
{
	// Gamma is a random variable between 0 and pi/2, that is picked from a distribution which is
	// lerp between uniform and gaussian w/ fancy params, based on parameter smoothness
	const float32 oneMinusS = 1.0f - smoothness;
	const float32 c = 1.0f / (pow3(oneMinusS)) - 1.0f;
	const float32 x = RandomUnilateral(series);
	const float32 gamma = Lerp(1.0f - x, Exp32(-c * x * x), smoothness) * PI_F / 2.0f;
	const float32 phi = RandomUnilateral(series) * 2.0f * PI_F;
	const Vec3 pureBounce = Reflect(rayDir, hitNormal);

    Quat quat = QuatFromAngleUnitAxis(gamma, Vec3::unitZ);
	quat = QuatFromAngleUnitAxis(phi, Vec3::unitX) * quat;
    // TODO hardcode unitX here?
	quat = QuatRotationFromUnitX(pureBounce) * quat;

	return quat * Vec3::unitX;
}

bool GetMaterial(const_string name, RaycastMaterial* material)
{
    if (StringEquals(name, ToString("Surfaces"))) {
        material->smoothness = 0.8f;
        material->albedo = Vec3 { 151, 162, 176 } / 255.0f;
        material->emission = 0.0f;
        material->emissionColor = Vec3::zero;
    }
    else if (StringEquals(name, ToString("LightHardGreen"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3::zero;
        material->emission = 1.5f;
        material->emissionColor = Vec3 { 0, 255, 129 } / 255.0f;
    }
    else if (StringEquals(name, ToString("LightSoftGreen"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3::zero;
        material->emission = 1.0f;
        material->emissionColor = Vec3 { 0, 255, 129 } / 255.0f;
    }
    else if (StringEquals(name, ToString("LightPink"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3::zero;
        material->emission = 1.0f;
        material->emissionColor = Vec3 { 237, 166, 255 } / 255.0f;
    }
    else if (StringEquals(name, ToString("LightRed"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3::zero;
        material->emission = 4.0f;
        material->emissionColor = Vec3 { 255, 0, 21 } / 255.0f;
    }
    else {
        return false;
    }

    return true;
}

Box CalculateAABBForTriangles(Array<RaycastTriangle> triangles)
{
    Box aabb = { .min = Vec3::one * 1e8, .max = -Vec3::one * 1e8 };

    for (uint32 i = 0; i < triangles.size; i++) {
        for (int k = 0; k < 3; k++) {
            const Vec3 v = triangles[i].pos[k];
            for (int e = 0; e < 3; e++) {
                aabb.min.e[e] = MinFloat32(v.e[e], aabb.min.e[e]);
                aabb.max.e[e] = MaxFloat32(v.e[e], aabb.max.e[e]);
            }
        }
    }

    return aabb;
}

uint32 CountTrianglesInsideBox(Array<RaycastTriangle> triangles, Box box)
{
    uint32 count = 0;
    for (uint32 i = 0; i < triangles.size; i++) {
        bool inside = false;
        for (int k = 0; k < 3; k++) {
            if (IsInsideInclusive(triangles[i].pos[k], box)) {
                // Only 1 vertex needs to be inside for the whole triangle to be inside
                inside = true;
                break;
            }
        }

        if (inside) {
            count++;
        }
    }

    return count;
}

Array<RaycastTriangle> GetTrianglesInsideBox(Array<RaycastTriangle> triangles, Box box, LinearAllocator* allocator)
{
    DynamicArray<RaycastTriangle, LinearAllocator> insideTriangles(allocator);
    for (uint32 i = 0; i < triangles.size; i++) {
        bool inside = false;
        for (int k = 0; k < 3; k++) {
            if (IsInsideInclusive(triangles[i].pos[k], box)) {
                // Only 1 vertex needs to be inside for the whole triangle to be inside
                inside = true;
                break;
            }
        }

        if (inside) {
            insideTriangles.Append(triangles[i]);
        }
    }

    return insideTriangles.ToArray();
}

bool DivideBvhTriangles(Array<RaycastTriangle> triangles, LinearAllocator* allocator,
                        Array<RaycastTriangle>* outTriangles1, Array<RaycastTriangle>* outTriangles2)
{
    const Box aabb = CalculateAABBForTriangles(triangles);
    const Vec3 aabbSize = aabb.max - aabb.min;
    const int longestDimension = (aabbSize.x >= aabbSize.y && aabbSize.x >= aabbSize.z) ? 0
        : ((aabbSize.y >= aabbSize.x && aabbSize.y >= aabbSize.z) ? 1 : 2);
    const float32 halfLongestDimension = aabbSize.e[longestDimension] / 2.0f;

    Box box1 = aabb;
    box1.max.e[longestDimension] -= halfLongestDimension;
    const uint32 countTriangles1 = CountTrianglesInsideBox(triangles, box1);
    if (countTriangles1 == triangles.size) {
        return false;
    }
    *outTriangles1 = GetTrianglesInsideBox(triangles, box1, allocator);

    Box box2 = aabb;
    box2.min.e[longestDimension] += halfLongestDimension;
    const uint32 countTriangles2 = CountTrianglesInsideBox(triangles, box2);
    if (countTriangles2 == triangles.size) {
        return false;
    }
    *outTriangles2 = GetTrianglesInsideBox(triangles, box2, allocator);

    return true;
}

bool FillRaycastBvhsRecursive(Array<RaycastTriangle> triangles, uint32 bvhMaxTriangles, uint32 depth, uint32* maxDepth,
                              LinearAllocator* allocator, DynamicArray<RaycastBvh, LinearAllocator>* outBvhs,
                              DynamicArray<RaycastTriangle, LinearAllocator>* outTriangles)
{
    if (depth > *maxDepth) {
        *maxDepth = depth;
    }

    const uint32 bvhIndex = outBvhs->size;
    outBvhs->Append();
    (*outBvhs)[bvhIndex].aabb = CalculateAABBForTriangles(triangles);

    if (triangles.size <= bvhMaxTriangles) {
        (*outBvhs)[bvhIndex].startTriangle = outTriangles->size;
        outTriangles->Append(triangles);
        (*outBvhs)[bvhIndex].endTriangle = outTriangles->size;
        (*outBvhs)[bvhIndex].skip = 1;
    }
    else {
        Array<RaycastTriangle> triangles1, triangles2;
        if (DivideBvhTriangles(triangles, allocator, &triangles1, &triangles2)) {
            const uint32 startBvh = outBvhs->size;
            if (!FillRaycastBvhsRecursive(triangles1, bvhMaxTriangles, depth + 1, maxDepth,
                                          allocator, outBvhs, outTriangles)) {
                return false;
            }
            if (!FillRaycastBvhsRecursive(triangles2, bvhMaxTriangles, depth + 1, maxDepth,
                                          allocator, outBvhs, outTriangles)) {
                return false;
            }

            (*outBvhs)[bvhIndex].startTriangle = 0;
            (*outBvhs)[bvhIndex].endTriangle = 0;
            (*outBvhs)[bvhIndex].skip = outBvhs->size - startBvh + 1;
        }
        else {
            (*outBvhs)[bvhIndex].startTriangle = outTriangles->size;
            outTriangles->Append(triangles);
            (*outBvhs)[bvhIndex].endTriangle = outTriangles->size;
            (*outBvhs)[bvhIndex].skip = 1;
        }
    }

    return true;
}

bool FillRaycastBvhs(Array<RaycastTriangle> triangles, uint32 bvhMaxTriangles, uint32* maxDepth,
                     LinearAllocator* allocator, DynamicArray<RaycastBvh, LinearAllocator>* outBvhs,
                     DynamicArray<RaycastTriangle, LinearAllocator>* outTriangles)
{
    *maxDepth = 0;
    return FillRaycastBvhsRecursive(triangles, bvhMaxTriangles, 0, maxDepth, allocator, outBvhs, outTriangles);
}

bool CreateRaycastGeometry(const LoadObjResult& obj, uint32 bvhMaxTriangles,
                           RaycastGeometry* geometry, LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    // Fill material info
    geometry->materialNames = allocator->NewArray<string>(obj.materials.size);
    if (geometry->materialNames.data == nullptr) {
        return false;
    }
    geometry->materials = allocator->NewArray<RaycastMaterial>(obj.materials.size);
    if (geometry->materials.data == nullptr) {
        return false;
    }

    for (uint32 i = 0; i < obj.materials.size; i++) {
        const_string materialName = obj.materials[i].name;
        geometry->materialNames[i] = allocator->NewArray<char>(materialName.size);
        MemCopy(geometry->materialNames[i].data, materialName.data, materialName.size);

        if (!GetMaterial(materialName, &geometry->materials[i])) {
            LOG_ERROR("Unrecognized material: %.*s\n", materialName.size, materialName.data);
            return false;
        }
    }

    // Fill mesh info
    geometry->meshes = allocator->NewArray<RaycastMesh>(obj.models.size);
    if (geometry->meshes.data == nullptr) {
        return false;
    }

    DynamicArray<RaycastBvh, LinearAllocator> bvhs(tempAllocator);
    DynamicArray<RaycastTriangle, LinearAllocator> triangles(tempAllocator);

    // Geometry stats
    uint32 trianglesTotal = 0;

    for (uint32 i = 0; i < obj.models.size; i++) {
        const_string name = obj.models[i].name;

        const uint32 numMeshTriangles = obj.models[i].triangles.size + obj.models[i].quads.size * 2;
        Array<RaycastTriangle> meshTriangles = tempAllocator->NewArray<RaycastTriangle>(numMeshTriangles);
        if (meshTriangles.data == nullptr) {
            return false;
        }
        trianglesTotal += numMeshTriangles;

        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                meshTriangles[j].pos[k] = t.v[k].pos;
            }
            meshTriangles[j].normal = normal;
            meshTriangles[j].materialIndex = t.materialIndex;
        }

        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const uint32 ind = obj.models[i].triangles.size + j * 2;
            const ObjQuad& q = obj.models[i].quads[j];
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                meshTriangles[ind].pos[k] = q.v[k].pos;
            }
            meshTriangles[ind].normal = normal;
            meshTriangles[ind].materialIndex = q.materialIndex;

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                meshTriangles[ind + 1].pos[k] = q.v[quadInd].pos;
            }
            meshTriangles[ind + 1].normal = normal;
            meshTriangles[ind + 1].materialIndex = q.materialIndex;
        }

        RaycastMesh& mesh = geometry->meshes[i];
        mesh.offset = Vec3::zero;
        mesh.inverseQuat = Quat::one;
        mesh.startBvh = bvhs.size;

        uint32 maxDepth;
        if (!FillRaycastBvhs(meshTriangles, bvhMaxTriangles, &maxDepth, tempAllocator, &bvhs, &triangles)) {
            LOG_ERROR("Failed to fill BVHs for mesh %.*s\n", name.size, name.data);
            return false;
        }

        mesh.endBvh = bvhs.size;
    }

    geometry->triangles = allocator->NewArray<RaycastTriangle>(triangles.size);
    geometry->triangles.CopyFrom(triangles.ToArray());

    geometry->bvhs = allocator->NewArray<RaycastBvh>(bvhs.size);
    geometry->bvhs.CopyFrom(bvhs.ToArray());

    // Report geometry stats
    LOG_INFO("BVH triangles / total triangles: %lu / %lu (%.02f %%)\n", triangles.size, trianglesTotal,
             (float32)triangles.size / (float32)trianglesTotal * 100.0f);

#if 0
    if (bvhMaxTriangles != UINT32_MAX_VALUE) {
        LOG_INFO("Total BVHs: %lu\n", totalBvhs);
        LOG_INFO("BVH triangles / total triangles: %lu / %lu (%.02f %%)\n", totalBvhTriangles, totalTriangles,
                 (float32)totalBvhTriangles / (float32)totalTriangles * 100.0f);

        LOG_INFO("Avg BVH triangles: %.02f\n", (float32)totalBvhTriangles / (float32)totalBvhs);
        DEBUG_ASSERT(maxBvhTrianglesInd != obj.models.size);
        const_string maxBvhTrianglesName = obj.models[maxBvhTrianglesInd].name;
        LOG_INFO("Max BVH triangles: %lu (mesh %.*s)\n",
                 maxBvhTriangles, maxBvhTrianglesName.size, maxBvhTrianglesName.data);

        LOG_INFO("Avg BVH nodes: %.02f\n", (float32)totalBvhs / (float32)obj.models.size);
        if (maxBvhNodesInd != obj.models.size) {
            const_string maxBvhNodesName = obj.models[maxBvhNodesInd].name;
            LOG_INFO("Max BVH nodes: %lu (mesh %.*s)\n",
                     maxBvhNodes, maxBvhNodesName.size, maxBvhNodesName.data);
        }
        if (maxBvhDepthInd != obj.models.size) {
            const_string maxBvhDepthName = obj.models[maxBvhDepthInd].name;
            LOG_INFO("Max BVH depth: %lu (mesh %.*s)\n",
                     maxBvhDepth, maxBvhDepthName.size, maxBvhDepthName.data);
        }
    }
#endif

    return true;
}

void TraverseMeshes(Vec3 rayOrigin, Vec3 rayDir, const RaycastGeometry& geometry, float32 minDist,
                    uint32* hitMaterialIndex, Vec3* hitNormal, float32* hitDist)
{
	for (uint32 m = 0; m < geometry.meshes.size; m++) {
		const RaycastMesh mesh = geometry.meshes[m];
		const Vec3 meshRayOrigin = mesh.inverseQuat * (rayOrigin + mesh.offset);
		const Vec3 meshRayDir = mesh.inverseQuat * rayDir;
		const Vec3 meshInverseRayDir = Reciprocal(meshRayDir);

		uint32 i = mesh.startBvh;
		while (i < mesh.endBvh) {
			const RaycastBvh bvh = geometry.bvhs[i];

			float32 tMin, tMax;
			const bool aabbHit = RayAABBIntersection(meshRayOrigin, meshInverseRayDir, bvh.aabb, &tMin, &tMax);
			if (aabbHit && !(tMin < 0.0 && tMax < 0.0) && tMin < *hitDist) {
				const uint32 tStart = bvh.startTriangle;
				const uint32 tEnd = bvh.endTriangle;
				for (uint32 j = tStart; j < tEnd; j++) {
                    const RaycastTriangle& triangle = geometry.triangles[j];
					const float32 dotNegRayNormal = Dot(meshRayDir, triangle.normal);
					if (dotNegRayNormal <= 0.0) {
						float32 t;
						const bool triangleHit = RayTriangleIntersection(meshRayOrigin, meshRayDir,
                                                                         triangle.pos[0],
                                                                         triangle.pos[1],
                                                                         triangle.pos[2],
                                                                         &t);
						if (triangleHit && t > minDist && t < *hitDist) {
							*hitMaterialIndex = triangle.materialIndex;
							*hitNormal = triangle.normal;
							*hitDist = t;
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
}

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, float32 minDist, float32 maxDist, const RaycastGeometry& geometry,
                  RandomSeries* series)
{
    UNREFERENCED_PARAMETER(series);
    ZoneScoped;

    const uint32 MAX_MATERIALS = geometry.materials.size;

	const float32 weightDiffuse = 0.5f;
	const float32 weightSpecular = 1.0f - weightDiffuse;
	const float32 smoothnessDiffuse = 0.2f;

    uint32 hitMaterialIndex = MAX_MATERIALS;
    Vec3 hitNormal;
	float hitDist = maxDist;
	TraverseMeshes(rayOrigin, rayDir, geometry, minDist, &hitMaterialIndex, &hitNormal, &hitDist);

    const RaycastMaterial* hitMaterial;
	if (hitMaterialIndex == MAX_MATERIALS) {
        return Vec3::zero;
	}
	else {
        hitMaterial = &geometry.materials[hitMaterialIndex];
		if (hitMaterial->emission > 0.0) {
			return hitMaterial->emission * hitMaterial->emissionColor;
		}
	}

	rayOrigin += rayDir * hitDist;
	Vec3 albedoMultDiffuse = hitMaterial->albedo;
	Vec3 albedoMultSpecular = hitMaterial->albedo;
	const float32 sampleWeightDiffuse = 1.0 / SAMPLES_DIFFUSE * weightDiffuse;
	const float32 sampleWeightSpecular = 1.0 / SAMPLES_SPECULAR * weightSpecular;
	const Vec3 sampleRayDirDiffuse = GetNewRayDirDiffuse(hitNormal, series);
	const Vec3 sampleRayDirSpecular = GetNewRayDirSpecular(rayDir, hitNormal, hitMaterial->smoothness, series);

	Vec3 color = Vec3::zero;

	// Diffuse
	for (uint32 s = 0; s < SAMPLES_DIFFUSE; s++) {
		Vec3 sampleRayOrigin = rayOrigin;
		Vec3 sampleRayDir = sampleRayDirDiffuse;

		for (uint32 i = 0; i < BOUNCES_DIFFUSE; i++) {
			hitMaterialIndex = MAX_MATERIALS;
			hitDist = maxDist;
			TraverseMeshes(sampleRayOrigin, sampleRayDir, geometry, minDist, &hitMaterialIndex, &hitNormal, &hitDist);

			if (hitMaterialIndex == MAX_MATERIALS) {
				break;
			}
			else {
				const RaycastMaterial& hitMaterialD = geometry.materials[hitMaterialIndex];
				if (hitMaterialD.emission > 0.0) {
					color += hitMaterialD.emission * sampleWeightDiffuse
                        * Multiply(hitMaterialD.emissionColor, albedoMultDiffuse);
                    break;
				}
				else if (i != BOUNCES_DIFFUSE - 1) {
					albedoMultDiffuse = Multiply(albedoMultDiffuse, hitMaterialD.albedo);
					sampleRayOrigin += sampleRayDir * hitDist;
					sampleRayDir = GetNewRayDirDiffuse(hitNormal, series);
				}
			}
		}
	}

	// Specular
	for (uint32 s = 0; s < SAMPLES_SPECULAR; s++) {
		Vec3 sampleRayOrigin = rayOrigin;
		Vec3 sampleRayDir = sampleRayDirSpecular;

		for (uint32 i = 0; i < BOUNCES_SPECULAR; i++) {
			hitMaterialIndex = MAX_MATERIALS;
			hitDist = maxDist;
			TraverseMeshes(sampleRayOrigin, sampleRayDir, geometry, minDist, &hitMaterialIndex, &hitNormal, &hitDist);

			if (hitMaterialIndex == MAX_MATERIALS) {
				break;
			}
			else {
				const RaycastMaterial& hitMaterialS = geometry.materials[hitMaterialIndex];
				if (hitMaterialS.emission > 0.0) {
					color += hitMaterialS.emission * sampleWeightSpecular
                        * Multiply(hitMaterialS.emissionColor, albedoMultSpecular);
					break;
				}
				else if (i != BOUNCES_SPECULAR - 1) {
					albedoMultSpecular = Multiply(albedoMultSpecular, hitMaterialS.albedo);
					sampleRayOrigin += sampleRayDir * hitDist;
					sampleRayDir = GetNewRayDirSpecular(rayDir, hitNormal, hitMaterialS.smoothness, series);
				}
			}
		}
	}

	return color;
}

#include "simd.cpp"

void TraverseMeshes_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, const RaycastGeometry& geometry, float32 minDist,
                      __m256i* hitMaterialIndex8, Vec3_8* hitNormal8, __m256* hitDist8)
{
    const __m256 minDist8 = _mm256_set1_ps(minDist);

	for (uint32 m = 0; m < geometry.meshes.size; m++) {
		const RaycastMesh mesh = geometry.meshes[m];
        const Vec3_8 meshOffset8 = Set1Vec3_8(mesh.offset);
        const Quat_8 meshInverseQuat8 = Set1Quat_8(mesh.inverseQuat);

		const Vec3_8 meshRayOrigin8 = Multiply_8(meshInverseQuat8, rayOrigin8 + meshOffset8);
		const Vec3_8 meshRayDir8 = Multiply_8(meshInverseQuat8, rayDir8);
		const Vec3_8 meshInverseRayDir8 = Reciprocal_8(meshRayDir8);

		uint32 i = mesh.startBvh;
		while (i < mesh.endBvh) {
			const RaycastBvh bvh = geometry.bvhs[i];

            __m256 tMin8, tMax8;
			const __m256 aabbHit8 = RayAABBIntersection_8(meshRayOrigin8, meshInverseRayDir8, bvh.aabb,
                                                          &tMin8, &tMax8);
            const __m256 hitPositive8 = _mm256_or_ps(_mm256_cmp_ps(tMin8, ZERO_8, _CMP_GE_OQ),
                                                     _mm256_cmp_ps(tMax8, ZERO_8, _CMP_GE_OQ));
            const __m256 hitCloser8 = _mm256_cmp_ps(tMin8, *hitDist8, _CMP_LT_OQ);

            if (AnyNonZero_8(aabbHit8) && AnyNonZero_8(hitPositive8) && AnyNonZero_8(hitCloser8)) {
                const uint32 tStart = bvh.startTriangle;
                const uint32 tEnd = bvh.endTriangle;
                for (uint32 j = tStart; j < tEnd; j++) {
                    const RaycastTriangle& triangle = geometry.triangles[j];
                    const Vec3_8 triangleNormal8 = Set1Vec3_8(triangle.normal);
                    const __m256 dotRayNormal8 = Dot_8(meshRayDir8, triangleNormal8);
                    const __m256 rayNormalSameDir8 = _mm256_cmp_ps(dotRayNormal8, ZERO_8, _CMP_LE_OQ);
                    if (AnyNonZero_8(rayNormalSameDir8)) {
                        __m256 t8;
                        const __m256 triangleHit8 = RayTriangleIntersection_8(meshRayOrigin8, meshRayDir8,
                                                                              triangle.pos[0],
                                                                              triangle.pos[1],
                                                                              triangle.pos[2],
                                                                              &t8);
                        const __m256 triangleCloser8 = _mm256_and_ps(_mm256_cmp_ps(t8, minDist8, _CMP_GT_OQ),
                                                                     _mm256_cmp_ps(t8, *hitDist8, _CMP_LT_OQ));
                        const __m256 triangleHitCloser8 = _mm256_and_ps(triangleHit8, triangleCloser8);
                        const __m256i triangleHitCloser8i = _mm256_castps_si256(triangleHitCloser8);

                        const __m256i materialIndex8 = _mm256_set1_epi32(triangle.materialIndex);
                        const Vec3_8 normal8 = Set1Vec3_8(triangle.normal);
                        *hitMaterialIndex8 = _mm256_blendv_epi8(*hitMaterialIndex8, materialIndex8, triangleHitCloser8i);
                        hitNormal8->x = _mm256_blendv_ps(hitNormal8->x, normal8.x, triangleHitCloser8);
                        hitNormal8->y = _mm256_blendv_ps(hitNormal8->y, normal8.y, triangleHitCloser8);
                        hitNormal8->z = _mm256_blendv_ps(hitNormal8->z, normal8.z, triangleHitCloser8);
                        *hitDist8 = _mm256_blendv_ps(*hitDist8, t8, triangleHitCloser8);
                    }
                }

                i++;
            }
            else {
                i += bvh.skip;
            }
        }
    }
}

Vec3_8 RaycastColor_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, float32 minDist, float32 maxDist,
                      const RaycastGeometry& geometry, RandomSeries* series)
{
    UNREFERENCED_PARAMETER(series);
    ZoneScoped;

    const __m256i MAX_MATERIALS_8 = _mm256_set1_epi32(geometry.materials.size);

    const __m256 weightDiffuse8 = _mm256_set1_ps(0.5f);
    const __m256 weightSpecular8 = ONE_8 - weightDiffuse8;
    const __m256 smoothnessDiffuse = _mm256_set1_ps(0.2f);

    __m256i hitMaterialIndex8 = MAX_MATERIALS_8;
    Vec3_8 hitNormal8;
    __m256 hitDist8 = _mm256_set1_ps(maxDist);
    TraverseMeshes_8(rayOrigin8, rayDir8, geometry, minDist, &hitMaterialIndex8, &hitNormal8, &hitDist8);

    const __m256i notHit8 = _mm256_cmpeq_epi32(hitMaterialIndex8, MAX_MATERIALS_8);
    Vec3_8 color = Set1Vec3_8(Vec3::zero);
    color.x = _mm256_blendv_ps(ONE_8, ZERO_8, _mm256_castsi256_ps(notHit8));

    return color;

#if 0
    const RaycastMaterial* hitMaterial;
    if (hitMaterialIndex == MAX_MATERIALS) {
        return Vec3::zero;
    }
    else {
        hitMaterial = &geometry.materials[hitMaterialIndex];
        if (hitMaterial->emission > 0.0) {
            return hitMaterial->emission * hitMaterial->emissionColor;
        }
    }

    rayOrigin += rayDir * hitDist;
    Vec3 albedoMultDiffuse = hitMaterial->albedo;
    Vec3 albedoMultSpecular = hitMaterial->albedo;
    const float32 sampleWeightDiffuse = 1.0 / SAMPLES_DIFFUSE * weightDiffuse;
    const float32 sampleWeightSpecular = 1.0 / SAMPLES_SPECULAR * weightSpecular;
    const Vec3 sampleRayDirDiffuse = GetNewRayDirDiffuse(hitNormal, series);
    const Vec3 sampleRayDirSpecular = GetNewRayDirSpecular(rayDir, hitNormal, hitMaterial->smoothness, series);

    Vec3 color = Vec3::zero;

    // Diffuse
    for (uint32 s = 0; s < SAMPLES_DIFFUSE; s++) {
        Vec3 sampleRayOrigin = rayOrigin;
        Vec3 sampleRayDir = sampleRayDirDiffuse;

        for (uint32 i = 0; i < BOUNCES_DIFFUSE; i++) {
            hitMaterialIndex = MAX_MATERIALS;
            hitDist = maxDist;
            TraverseMeshes(sampleRayOrigin, sampleRayDir, geometry, minDist, &hitMaterialIndex, &hitNormal, &hitDist);

            if (hitMaterialIndex == MAX_MATERIALS) {
                break;
            }
            else {
                const RaycastMaterial& hitMaterialD = geometry.materials[hitMaterialIndex];
                if (hitMaterialD.emission > 0.0) {
                    color += hitMaterialD.emission * sampleWeightDiffuse
                        * Multiply(hitMaterialD.emissionColor, albedoMultDiffuse);
                    break;
                }
                else if (i != BOUNCES_DIFFUSE - 1) {
                    albedoMultDiffuse = Multiply(albedoMultDiffuse, hitMaterialD.albedo);
                    sampleRayOrigin += sampleRayDir * hitDist;
                    sampleRayDir = GetNewRayDirDiffuse(hitNormal, series);
                }
            }
        }
    }

    // Specular
    for (uint32 s = 0; s < SAMPLES_SPECULAR; s++) {
        Vec3 sampleRayOrigin = rayOrigin;
        Vec3 sampleRayDir = sampleRayDirSpecular;

        for (uint32 i = 0; i < BOUNCES_SPECULAR; i++) {
            hitMaterialIndex = MAX_MATERIALS;
            hitDist = maxDist;
            TraverseMeshes(sampleRayOrigin, sampleRayDir, geometry, minDist, &hitMaterialIndex, &hitNormal, &hitDist);

            if (hitMaterialIndex == MAX_MATERIALS) {
                break;
            }
            else {
                const RaycastMaterial& hitMaterialS = geometry.materials[hitMaterialIndex];
                if (hitMaterialS.emission > 0.0) {
                    color += hitMaterialS.emission * sampleWeightSpecular
                        * Multiply(hitMaterialS.emissionColor, albedoMultSpecular);
                    break;
                }
                else if (i != BOUNCES_SPECULAR - 1) {
                    albedoMultSpecular = Multiply(albedoMultSpecular, hitMaterialS.albedo);
                    sampleRayOrigin += sampleRayDir * hitDist;
                    sampleRayDir = GetNewRayDirSpecular(rayDir, hitNormal, hitMaterialS.smoothness, series);
                }
            }
        }
    }

    return color;
#endif
}

struct RaycastThreadWorkCommon
{
    const RaycastGeometry* geometry;
    uint32 seed;
    Vec3 filmTopLeft;
    Vec3 filmUnitOffsetX;
    Vec3 filmUnitOffsetY;
    Vec3 cameraPos;
    uint32 width, height;
    float32 minDist;
    float32 maxDist;
};

struct RaycastThreadWork
{
    const RaycastThreadWorkCommon* common;
    Vec3* colorHdr;
    uint32 minX, maxX, minY, maxY;
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(threadIndex);
    UNREFERENCED_PARAMETER(queue);
    ZoneScoped;

    RaycastThreadWork* work = (RaycastThreadWork*)data;
    const RaycastThreadWorkCommon& common = *work->common;
    const uint32 minX = work->minX;
    const uint32 maxX = work->maxX;
    const uint32 minY = work->minY;
    const uint32 maxY = work->maxY;

    const __m256 minX8 = _mm256_set1_ps((float32)minX);
    const Vec3_8 filmUnitOffsetX8 = Set1Vec3_8(common.filmUnitOffsetX);
    const Vec3_8 filmUnitOffsetY8 = Set1Vec3_8(common.filmUnitOffsetY);
    const Vec3_8 filmTopLeft8 = Set1Vec3_8(common.filmTopLeft);
    const Vec3_8 cameraPos8 = Set1Vec3_8(common.cameraPos);

    for (uint32 y = minY; y < maxY; y++) {
#if 1
        const __m256 y8 = _mm256_set1_ps((float32)y);
        const Vec3_8 filmOffsetY8 = filmUnitOffsetY8 * y8;

        uint32 n = (maxX - minX) / 8;
        for (uint32 i = 0; i < n; i++) {
            const __m256 baseX8 = _mm256_set1_ps((float32)i * 8.0f);
            const __m256 x8 = _mm256_set_ps(7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f, 0.0f) + minX8 + baseX8;

            const Vec3_8 filmOffsetX8 = filmUnitOffsetX8 * x8;
            const Vec3_8 filmPos8 = filmTopLeft8 + filmOffsetX8 + filmOffsetY8;
            const Vec3_8 rayDir8 = Normalize_8(filmPos8 - cameraPos8);
            RandomSeries series = { .state = 0 }; // TODO init
            series.state += common.seed;

            const Vec3_8 raycastColor8 = RaycastColor_8(cameraPos8, rayDir8, common.minDist, common.maxDist,
                                                        *common.geometry, &series);
            const uint32 basePixelIndex = y * common.width + minX + i * 8;
            StoreVec3_8(raycastColor8, work->colorHdr + basePixelIndex);
        }

#else
        const Vec3 filmOffsetY = common.filmUnitOffsetY * (float32)y;

        for (uint32 x = minX; x < maxX; x++) {
            const uint32 pixelIndex = y * common.width + x;

            const Vec3 filmOffsetX = common.filmUnitOffsetX * (float32)x;
            const Vec3 filmPos = common.filmTopLeft + filmOffsetX + filmOffsetY;
            const Vec3 rayDir = Normalize(filmPos - common.cameraPos);
            RandomSeries series = { .state = WangHash(pixelIndex) };
            series.state += common.seed;

            const Vec3 raycastColor = RaycastColor(common.cameraPos, rayDir, common.minDist, common.maxDist,
                                                   *common.geometry, &series);

            work->colorHdr[pixelIndex] = raycastColor;
        }
#endif
    }
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, CanvasState* canvas, uint32* pixels,
                    LinearAllocator* allocator, AppWorkQueue* queue)
{
    ZoneScoped;
    const uint32 numPixels = width * height;

    {
        ZoneScopedN("Clear");
        // TODO make GPU clear pixels? eh, maybe not
        MemSet(canvas->colorHdr.data, 0, numPixels * sizeof(Vec3));
    }

    {
        ZoneScopedN("ProduceWork");

        const Quat inverseCameraRot = Inverse(cameraRot);
        const Vec3 cameraUp = inverseCameraRot * Vec3::unitZ;
        const Vec3 cameraForward = inverseCameraRot * Vec3::unitX;
        const Vec3 cameraLeft = inverseCameraRot * Vec3::unitY;

        const float32 filmDist = 1.0f;
        const float32 filmHeight = tanf(fov / 2.0f) * 2.0f;
        const float32 filmWidth = filmHeight * (float32)width / (float32)height;

        const Vec3 filmTopLeft = cameraPos + cameraForward * filmDist
            + (filmWidth / 2.0f) * cameraLeft + (filmHeight / 2.0f) * cameraUp;
        const Vec3 filmUnitOffsetX = -cameraLeft * filmWidth / (float32)(width - 1);
        const Vec3 filmUnitOffsetY = -cameraUp * filmHeight / (float32)(height - 1);

        RaycastThreadWorkCommon workCommon = {
            .geometry = &geometry,
            .seed = (uint32)rand(),
            .filmTopLeft = filmTopLeft,
            .filmUnitOffsetX = filmUnitOffsetX,
            .filmUnitOffsetY = filmUnitOffsetY,
            .cameraPos = cameraPos,
            .width = width,
            .height = height,
            .minDist = 0.0f,
            .maxDist = 20.0f,
        };

        const uint32 WORK_TILE_SIZE = 16;
        static_assert(WORK_TILE_SIZE % 8 == 0); // AVX width
        const uint32 wholeTilesX = width / WORK_TILE_SIZE;
        const uint32 wholeTilesY = height / WORK_TILE_SIZE;
        const uint32 numWorkEntries = wholeTilesX * wholeTilesY;
        Array<RaycastThreadWork> workEntries = allocator->NewArray<RaycastThreadWork>(numWorkEntries);

        for (uint32 tileY = 0; tileY < wholeTilesY; tileY++) {
            for (uint32 tileX = 0; tileX < wholeTilesX; tileX++) {
                const uint32 workIndex = tileY * wholeTilesX + tileX;
                RaycastThreadWork* work = &workEntries[workIndex];
                work->common = &workCommon;
                work->colorHdr = canvas->colorHdr.data,
                work->minX = tileX * WORK_TILE_SIZE;
                work->maxX = (tileX + 1) * WORK_TILE_SIZE;
                work->minY = tileY * WORK_TILE_SIZE;
                work->maxY = (tileY + 1) * WORK_TILE_SIZE;

                if (!TryAddWork(queue, &RaycastThreadProc, work)) {
                    CompleteAllWork(queue, 0);
                    tileX--;
                    continue;
                }
            }
        }
    }

    {
        ZoneScopedN("FinishWork");
        CompleteAllWork(queue, 0);
    }

    {
        ZoneScopedN("TranslateHdr");

        float32 maxColor = 1.0f;
        for (uint32 i = 0; i < numPixels; i++) {
            const Vec3 colorHdr = canvas->colorHdr[i];
            maxColor = MaxFloat32(maxColor, colorHdr.r);
            maxColor = MaxFloat32(maxColor, colorHdr.g);
            maxColor = MaxFloat32(maxColor, colorHdr.b);
        }

        for (uint32 i = 0; i < numPixels; i++) {
            //const Vec3 colorNormalized = canvas->colorHdr[i] / maxColor;
            const Vec3 colorNormalized = Vec3 {
                ClampFloat32(canvas->colorHdr[i].r, 0.0f, 1.0f),
                ClampFloat32(canvas->colorHdr[i].g, 0.0f, 1.0f),
                ClampFloat32(canvas->colorHdr[i].b, 0.0f, 1.0f)
            };
            // TODO gamma correction?
            const uint8 r = (uint8)(colorNormalized.r * 255.0f);
            const uint8 g = (uint8)(colorNormalized.g * 255.0f);
            const uint8 b = (uint8)(colorNormalized.b * 255.0f);
            pixels[i] = ((uint32)0xff << 24) + ((uint32)b << 16) + ((uint32)g << 8) + r;
        }
    }
}

bool LoadRaytracePipeline(const VulkanWindow& window, VkCommandPool commandPool, uint32 width, uint32 height,
                          const RaycastGeometry& geometry, LinearAllocator* allocator, VulkanRaytracePipeline* pipeline)
{
    QueueFamilyInfo queueFamilyInfo = GetQueueFamilyInfo(window.surface, window.physicalDevice, allocator);
    if (!queueFamilyInfo.hasComputeFamily) {
        LOG_ERROR("Device has no compute family\n");
        return false;
    }

    // compute queue and command pool
    {
        VkDeviceQueueCreateInfo queueCreateInfo = {};
        queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
        queueCreateInfo.pNext = nullptr;
        queueCreateInfo.queueFamilyIndex = queueFamilyInfo.computeFamilyIndex;
        queueCreateInfo.queueCount = 1;
        vkGetDeviceQueue(window.device, queueFamilyInfo.computeFamilyIndex, 0, &pipeline->queue);

        VkCommandPoolCreateInfo commandPoolCreateInfo = {};
        commandPoolCreateInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        commandPoolCreateInfo.queueFamilyIndex = queueFamilyInfo.computeFamilyIndex;
        commandPoolCreateInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

        if (vkCreateCommandPool(window.device, &commandPoolCreateInfo, nullptr,
                                &pipeline->commandPool) != VK_SUCCESS) {
            LOG_ERROR("vkCreateCommandPool failed\n");
            return false;
        }
    }

    // meshes, triangles, BVHs
    {
        SCOPED_ALLOCATOR_RESET(*allocator);

        pipeline->meshes.Clear();
        for (uint32 i = 0; i < geometry.meshes.size; i++) {
            const RaycastMesh& srcMesh = geometry.meshes[i];
            ComputeMesh* newMesh = pipeline->meshes.Append();
            newMesh->offset = srcMesh.offset;
            newMesh->inverseQuat = Vec4 {
                srcMesh.inverseQuat.x, srcMesh.inverseQuat.y, srcMesh.inverseQuat.z, srcMesh.inverseQuat.w
            };
            newMesh->startBvh = srcMesh.startBvh;
            newMesh->endBvh = srcMesh.endBvh;
        }

        Array<ComputeTriangle> triangles = allocator->NewArray<ComputeTriangle>(geometry.triangles.size);
        for (uint32 i = 0; i < geometry.triangles.size; i++) {
            const RaycastTriangle& srcTriangle = geometry.triangles[i];
            triangles[i].a = srcTriangle.pos[0];
            triangles[i].b = srcTriangle.pos[1];
            triangles[i].c = srcTriangle.pos[2];
            triangles[i].normal = srcTriangle.normal;
            triangles[i].materialIndex = srcTriangle.materialIndex;
        }

        Array<ComputeBvh> bvhs = allocator->NewArray<ComputeBvh>(geometry.bvhs.size);
        for (uint32 i = 0; i < geometry.bvhs.size; i++) {
            const RaycastBvh& srcBvh = geometry.bvhs[i];
            bvhs[i].aabbMin = srcBvh.aabb.min;
            bvhs[i].aabbMax = srcBvh.aabb.max;
            bvhs[i].startTriangle = srcBvh.startTriangle;
            bvhs[i].endTriangle = srcBvh.endTriangle;
            bvhs[i].skip = srcBvh.skip;
        }

        const VkDeviceSize triangleBufferSize = triangles.size * sizeof(ComputeTriangle);
        const VkDeviceSize bvhBufferSize = bvhs.size * sizeof(ComputeBvh);
        pipeline->numTriangles = triangles.size;
        pipeline->numBvhs = bvhs.size;

        const VkDeviceSize stagingBufferSize = MaxUInt64(triangleBufferSize, bvhBufferSize);
        VulkanBuffer stagingBuffer;
        if (!CreateVulkanBuffer(stagingBufferSize,
                                VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &stagingBuffer)) {
            LOG_ERROR("CreateVulkanBuffer failed\n");
            return false;
        }
        defer(DestroyVulkanBuffer(window.device, &stagingBuffer));
        void* data;

        // Create triangle storage buffer
        if (!CreateVulkanBuffer(triangleBufferSize,
                                VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                                window.device, window.physicalDevice, &pipeline->triangles)) {
            LOG_ERROR("CreateVulkanBuffer failed\n");
            return false;
        }

        // Copy triangle data to staging buffer
        vkMapMemory(window.device, stagingBuffer.memory, 0, triangleBufferSize, 0, &data);
        MemCopy(data, triangles.data, triangleBufferSize);
        vkUnmapMemory(window.device, stagingBuffer.memory);

        // Copy triangle data to final storage buffer
        {
            SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, commandPool, window.graphicsQueue);
            CopyBuffer(commandBuffer, stagingBuffer.buffer, pipeline->triangles.buffer, triangleBufferSize);
        }

        // Create bvh storage buffer
        if (!CreateVulkanBuffer(bvhBufferSize,
                                VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                                window.device, window.physicalDevice, &pipeline->bvhs)) {
            LOG_ERROR("CreateVulkanBuffer failed\n");
            return false;
        }

        // Copy bvh data to staging buffer
        vkMapMemory(window.device, stagingBuffer.memory, 0, bvhBufferSize, 0, &data);
        MemCopy(data, bvhs.data, bvhBufferSize);
        vkUnmapMemory(window.device, stagingBuffer.memory);

        // Copy bvh data to final storage buffer
        {
            SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, commandPool, window.graphicsQueue);
            CopyBuffer(commandBuffer, stagingBuffer.buffer, pipeline->bvhs.buffer, bvhBufferSize);
        }
    }

    // uniform buffer
    {
        const VkDeviceSize uniformBufferSize = sizeof(ComputeUbo);
        if (!CreateVulkanBuffer(uniformBufferSize,
                                VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &pipeline->uniform)) {
            LOG_ERROR("CreateBuffer failed for vertex buffer\n");
            return false;
        }
    }

    // image
    {
        const VkFormat format = VK_FORMAT_R16G16B16A16_SFLOAT;
        VkFormatProperties formatProperties;
        vkGetPhysicalDeviceFormatProperties(window.physicalDevice, format, &formatProperties);
        if ((formatProperties.optimalTilingFeatures & VK_FORMAT_FEATURE_STORAGE_IMAGE_BIT) == 0) {
            LOG_ERROR("Storage image unsupported for format %d\n", format);
            return false;
        }

        if (!CreateImage(window.device, window.physicalDevice, width, height, format,
                         VK_IMAGE_TILING_OPTIMAL,
                         VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT,
                         VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                         &pipeline->image.image, &pipeline->image.memory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        if (!CreateImageView(window.device, pipeline->image.image, format, VK_IMAGE_ASPECT_COLOR_BIT,
                             &pipeline->image.view)) {
            LOG_ERROR("CreateImageView failed\n");
            return false;
        }

        SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, commandPool, window.graphicsQueue);
        TransitionImageLayout(commandBuffer, pipeline->image.image,
                              VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    }

    // descriptor set layout
    {
        VkDescriptorSetLayoutBinding layoutBindings[4] = {};

        layoutBindings[0].binding = 0;
        layoutBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        layoutBindings[0].descriptorCount = 1;
        layoutBindings[0].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        layoutBindings[0].pImmutableSamplers = nullptr;

        layoutBindings[1].binding = 1;
        layoutBindings[1].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        layoutBindings[1].descriptorCount = 1;
        layoutBindings[1].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        layoutBindings[1].pImmutableSamplers = nullptr;

        layoutBindings[2].binding = 2;
        layoutBindings[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        layoutBindings[2].descriptorCount = 1;
        layoutBindings[2].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        layoutBindings[2].pImmutableSamplers = nullptr;

        layoutBindings[3].binding = 3;
        layoutBindings[3].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        layoutBindings[3].descriptorCount = 1;
        layoutBindings[3].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        layoutBindings[3].pImmutableSamplers = nullptr;

        VkDescriptorSetLayoutCreateInfo layoutCreateInfo = {};
        layoutCreateInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layoutCreateInfo.bindingCount = C_ARRAY_LENGTH(layoutBindings);
        layoutCreateInfo.pBindings = layoutBindings;

        if (vkCreateDescriptorSetLayout(window.device, &layoutCreateInfo, nullptr,
                                        &pipeline->descriptorSetLayout) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorSetLayout failed\n");
            return false;
        }
    }

    // descriptor pool
    {
        VkDescriptorPoolSize poolSizes[3] = {};

        poolSizes[0].type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        poolSizes[0].descriptorCount = 1;

        poolSizes[1].type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        poolSizes[1].descriptorCount = 1;

        poolSizes[2].type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        poolSizes[2].descriptorCount = 2;

        VkDescriptorPoolCreateInfo poolInfo = {};
        poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        poolInfo.poolSizeCount = C_ARRAY_LENGTH(poolSizes);
        poolInfo.pPoolSizes = poolSizes;
        poolInfo.maxSets = 1;

        if (vkCreateDescriptorPool(window.device, &poolInfo, nullptr, &pipeline->descriptorPool) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorPool failed\n");
            return false;
        }
    }

    // descriptor set
    {
        VkDescriptorSetAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        allocInfo.descriptorPool = pipeline->descriptorPool;
        allocInfo.descriptorSetCount = 1;
        allocInfo.pSetLayouts = &pipeline->descriptorSetLayout;

        if (vkAllocateDescriptorSets(window.device, &allocInfo, &pipeline->descriptorSet) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateDescriptorSets failed\n");
            return false;
        }

        VkWriteDescriptorSet descriptorWrites[4] = {};

        const VkDescriptorImageInfo imageInfo = {
            .sampler = VK_NULL_HANDLE,
            .imageView = pipeline->image.view,
            .imageLayout = VK_IMAGE_LAYOUT_GENERAL
        };
        descriptorWrites[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        descriptorWrites[0].dstSet = pipeline->descriptorSet;
        descriptorWrites[0].dstBinding = 0;
        descriptorWrites[0].dstArrayElement = 0;
        descriptorWrites[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        descriptorWrites[0].descriptorCount = 1;
        descriptorWrites[0].pImageInfo = &imageInfo;

        const VkDescriptorBufferInfo uniformBufferInfo = {
            .buffer = pipeline->uniform.buffer,
            .offset = 0,
            .range = sizeof(ComputeUbo),
        };
        descriptorWrites[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        descriptorWrites[1].dstSet = pipeline->descriptorSet;
        descriptorWrites[1].dstBinding = 1;
        descriptorWrites[1].dstArrayElement = 0;
        descriptorWrites[1].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        descriptorWrites[1].descriptorCount = 1;
        descriptorWrites[1].pBufferInfo = &uniformBufferInfo;

        const VkDescriptorBufferInfo triangleBufferInfo = {
            .buffer = pipeline->triangles.buffer,
            .offset = 0,
            .range = pipeline->numTriangles * sizeof(ComputeTriangle),
        };
        descriptorWrites[2].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        descriptorWrites[2].dstSet = pipeline->descriptorSet;
        descriptorWrites[2].dstBinding = 2;
        descriptorWrites[2].dstArrayElement = 0;
        descriptorWrites[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        descriptorWrites[2].descriptorCount = 1;
        descriptorWrites[2].pBufferInfo = &triangleBufferInfo;

        const VkDescriptorBufferInfo bvhBufferInfo = {
            .buffer = pipeline->bvhs.buffer,
            .offset = 0,
            .range = pipeline->numBvhs * sizeof(ComputeBvh),
        };
        descriptorWrites[3].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        descriptorWrites[3].dstSet = pipeline->descriptorSet;
        descriptorWrites[3].dstBinding = 3;
        descriptorWrites[3].dstArrayElement = 0;
        descriptorWrites[3].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        descriptorWrites[3].descriptorCount = 1;
        descriptorWrites[3].pBufferInfo = &bvhBufferInfo;

        vkUpdateDescriptorSets(window.device, C_ARRAY_LENGTH(descriptorWrites), descriptorWrites, 0, nullptr);
    }

    // pipeline
    {
        const Array<uint8> shaderCode = LoadEntireFile(ToString("data/shaders/raytracer.comp.spv"), allocator);
        if (shaderCode.data == nullptr) {
            LOG_ERROR("Failed to load compute shader code\n");
            return false;
        }

        VkShaderModule shaderModule;
        if (!CreateShaderModule(shaderCode, window.device, &shaderModule)) {
            LOG_ERROR("Failed to create compute shader module\n");
            return false;
        }
        defer(vkDestroyShaderModule(window.device, shaderModule, nullptr));

        VkPipelineShaderStageCreateInfo shaderStageCreateInfo = {};
        shaderStageCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        shaderStageCreateInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        shaderStageCreateInfo.module = shaderModule;
        shaderStageCreateInfo.pName = "main";

        VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo = {};
        pipelineLayoutCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipelineLayoutCreateInfo.setLayoutCount = 1;
        pipelineLayoutCreateInfo.pSetLayouts = &pipeline->descriptorSetLayout;
        pipelineLayoutCreateInfo.pushConstantRangeCount = 0;
        pipelineLayoutCreateInfo.pPushConstantRanges = nullptr;

        if (vkCreatePipelineLayout(window.device, &pipelineLayoutCreateInfo, nullptr,
                                   &pipeline->pipelineLayout) != VK_SUCCESS) {
            LOG_ERROR("vkCreatePipelineLayout failed\n");
            return false;
        }

        VkComputePipelineCreateInfo pipelineCreateInfo = {};
        pipelineCreateInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
        pipelineCreateInfo.pNext = nullptr;
        pipelineCreateInfo.flags = 0;
        pipelineCreateInfo.stage = shaderStageCreateInfo;
        pipelineCreateInfo.layout = pipeline->pipelineLayout;

        if (vkCreateComputePipelines(window.device, VK_NULL_HANDLE, 1, &pipelineCreateInfo, nullptr,
                                     &pipeline->pipeline) != VK_SUCCESS) {
            LOG_ERROR("vkCreateComputePipelines failed\n");
            return false;
        }
    }

    // command buffer
    {
        VkCommandBufferAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        allocInfo.pNext = nullptr;
        allocInfo.commandPool = pipeline->commandPool;
        allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        allocInfo.commandBufferCount = 1;

        if (vkAllocateCommandBuffers(window.device, &allocInfo, &pipeline->commandBuffer) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateCommandBuffers failed\n");
            return false;
        }

        VkCommandBufferBeginInfo beginInfo = {};
        beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        beginInfo.flags = 0;
        beginInfo.pInheritanceInfo = nullptr;

        if (vkBeginCommandBuffer(pipeline->commandBuffer, &beginInfo) != VK_SUCCESS) {
            LOG_ERROR("vkBeginCommandBuffer failed\n");
            return false;
        }

        vkCmdBindPipeline(pipeline->commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline->pipeline);
        vkCmdBindDescriptorSets(pipeline->commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                                pipeline->pipelineLayout, 0, 1, &pipeline->descriptorSet, 0, 0);

        const uint32 batchSize = VulkanRaytracePipeline::BATCH_SIZE;
        vkCmdDispatch(pipeline->commandBuffer, width / batchSize, height / batchSize, 1);

        if (vkEndCommandBuffer(pipeline->commandBuffer) != VK_SUCCESS) {
            LOG_ERROR("vkEndCommandBuffer failed\n");
            return false;
        }
    }

    // fence
    {
        VkFenceCreateInfo fenceCreateInfo = {};
        fenceCreateInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fenceCreateInfo.flags = 0;

        if (vkCreateFence(window.device, &fenceCreateInfo, nullptr, &pipeline->fence) != VK_SUCCESS) {
            LOG_ERROR("vkCreateFence failed\n");
            return false;
        }
    }

    return true;
}

void UnloadRaytracePipeline(VkDevice device, VulkanRaytracePipeline* pipeline)
{
    vkDestroyFence(device, pipeline->fence, nullptr);
    vkDestroyPipeline(device, pipeline->pipeline, nullptr);

    vkDestroyPipelineLayout(device, pipeline->pipelineLayout, nullptr);
    vkDestroyDescriptorPool(device, pipeline->descriptorPool, nullptr);
    vkDestroyDescriptorSetLayout(device, pipeline->descriptorSetLayout, nullptr);

    DestroyVulkanImage(device, &pipeline->image);
    DestroyVulkanBuffer(device, &pipeline->uniform);

    DestroyVulkanBuffer(device, &pipeline->bvhs);
    DestroyVulkanBuffer(device, &pipeline->triangles);

    vkDestroyCommandPool(device, pipeline->commandPool, nullptr);
}
