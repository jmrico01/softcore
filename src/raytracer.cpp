#include "raytracer.h"

#include <intrin.h>
#include <time.h>

#include <stb_image_write.h>
#include <Tracy.hpp>
#include <TracyC.h>

struct RandomSeries
{
    uint32 state;
};

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

bool FillRaycastMeshBvh(RaycastMeshBvh* bvh, Box minBox, uint32 bvhMaxTriangles, bool firstPass,
                        uint32 depth, uint32* maxDepth, Array<RaycastTriangle> triangles,
                        LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    DEBUG_ASSERT(bvh != nullptr);

    if (firstPass) {
        *maxDepth = 0;
    }
    else if (depth > *maxDepth) {
        *maxDepth = depth;
    }

    DynamicArray<RaycastTriangle, LinearAllocator> insideTriangles(tempAllocator);

    bvh->aabb.min = Vec3::one * 1e8;
    bvh->aabb.max = -Vec3::one * 1e8;
    for (uint32 i = 0; i < triangles.size; i++) {
        bool inside = false;
        for (int k = 0; k < 3; k++) {
            if (IsInsideInclusive(triangles[i].pos[k], minBox)) {
                // Only 1 vertex needs to be inside for the whole triangle to be inside
                inside = true;
                break;
            }
        }

        if (inside) {
            for (int k = 0; k < 3; k++) {
                const Vec3 v = triangles[i].pos[k];
                for (int e = 0; e < 3; e++) {
                    bvh->aabb.min.e[e] = MinFloat32(bvh->aabb.min.e[e], v.e[e]);
                    bvh->aabb.max.e[e] = MaxFloat32(bvh->aabb.max.e[e], v.e[e]);
                }
            }

            insideTriangles.Append(triangles[i]);
        }
    }

    // TODO examine the second case. # of triangles should decrease reasonably each partition,
    // otherwise we're duplicating checks
    if (insideTriangles.size <= bvhMaxTriangles || (!firstPass && insideTriangles.size == triangles.size)) {
        bvh->child1 = nullptr;
        bvh->child2 = nullptr;
        bvh->triangles = allocator->NewArray<RaycastTriangle>(insideTriangles.size);
        if (bvh->triangles.data == nullptr) {
            return false;
        }
        bvh->triangles.CopyFrom(insideTriangles.ToArray());
    }
    else {
        const Vec3 boxSize = bvh->aabb.max - bvh->aabb.min;
        const int longestDimension = (boxSize.x >= boxSize.y && boxSize.x >= boxSize.z) ? 0
            : ((boxSize.y >= boxSize.x && boxSize.y >= boxSize.z) ? 1 : 2);
        const float32 halfLongestDimension = boxSize.e[longestDimension] / 2.0f;

        Box minBox1 = bvh->aabb;
        minBox1.max.e[longestDimension] -= halfLongestDimension;
        Box minBox2 = bvh->aabb;
        minBox2.min.e[longestDimension] += halfLongestDimension;

        bvh->child1 = allocator->New<RaycastMeshBvh>();
        if (bvh->child1 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBvh(bvh->child1, minBox1, bvhMaxTriangles, false, depth + 1, maxDepth,
                                insideTriangles.ToArray(), allocator, tempAllocator)) {
            return false;
        }

        bvh->child2 = allocator->New<RaycastMeshBvh>();
        if (bvh->child2 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBvh(bvh->child2, minBox2, bvhMaxTriangles, false, depth + 1, maxDepth,
                                insideTriangles.ToArray(), allocator, tempAllocator)) {
            return false;
        }
    }

    return true;
}

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, uint32 bvhMaxTriangles,
                                      LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    RaycastGeometry geometry = {};

    geometry.materialNames = allocator->NewArray<string>(obj.materials.size);
    if (geometry.materialNames.data == nullptr) {
        return geometry;
    }
    geometry.materials = allocator->NewArray<RaycastMaterial>(obj.materials.size);
    if (geometry.materials.data == nullptr) {
        return geometry;
    }

    for (uint32 i = 0; i < obj.materials.size; i++) {
        const_string materialName = obj.materials[i].name;
        geometry.materialNames[i] = allocator->NewArray<char>(materialName.size);
        MemCopy(geometry.materialNames[i].data, materialName.data, materialName.size);

        if (!GetMaterial(materialName, &geometry.materials[i])) {
            LOG_ERROR("Unrecognized material: %.*s\n", materialName.size, materialName.data);
            return geometry;
        }
    }

    geometry.meshes = allocator->NewArray<RaycastMesh>(obj.models.size);
    if (geometry.meshes.data == nullptr) {
        return geometry;
    }

    uint32 totalTriangles = 0;
    uint32 totalBvhs = 0;
    uint32 maxBvhTrianglesInd = obj.models.size;
    uint32 maxBvhTriangles = 0;
    uint32 totalBvhTriangles = 0;
    uint32 maxBvhDepthInd = obj.models.size;
    uint32 maxBvhDepth = 0;
    uint32 maxBvhNodesInd = obj.models.size;
    uint32 maxBvhNodes = 0;

    for (uint32 i = 0; i < obj.models.size; i++) {
        const_string name = obj.models[i].name;
        RaycastMesh& mesh = geometry.meshes[i];

        const uint32 numTriangles = obj.models[i].triangles.size + obj.models[i].quads.size * 2;
        mesh.numTriangles = numTriangles;
        totalTriangles += numTriangles;

        Array<RaycastTriangle> triangles = tempAllocator->NewArray<RaycastTriangle>(numTriangles);
        if (triangles.data == nullptr) {
            LOG_ERROR("Failed to allocate triangles for raycast mesh %.*s\n", name.size, name.data);
            geometry.meshes.data = nullptr;
            return geometry;
        }

        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                triangles[j].pos[k] = t.v[k].pos;
            }
            triangles[j].normal = normal;
            triangles[j].materialIndex = t.materialIndex;
        }
        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const uint32 ind = obj.models[i].triangles.size + j * 2;
            const ObjQuad& q = obj.models[i].quads[j];
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                triangles[ind].pos[k] = q.v[k].pos;
            }
            triangles[ind].normal = normal;
            triangles[ind].materialIndex = q.materialIndex;

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                triangles[ind + 1].pos[k] = q.v[quadInd].pos;
            }
            triangles[ind + 1].normal = normal;
            triangles[ind + 1].materialIndex = q.materialIndex;
        }

        const Box veryBigBox = { .min = -Vec3::one * 1e8, .max = Vec3::one * 1e8 };
        uint32 maxDepth;
        if (!FillRaycastMeshBvh(&mesh.bvh, veryBigBox, bvhMaxTriangles, true, 0, &maxDepth,
                                triangles, allocator, tempAllocator)) {
            LOG_ERROR("Failed to fill raycast mesh box for mesh %.*s\n", name.size, name.data);
            geometry.meshes.data = nullptr;
            return geometry;
        }

        if (maxDepth > maxBvhDepth) {
            maxBvhDepthInd = i;
            maxBvhDepth = maxDepth;
        }

        uint32 numBvhNodes = 0;
        DynamicArray<RaycastMeshBvh*, LinearAllocator> bvhStack(tempAllocator);
        bvhStack.Append(&mesh.bvh);
        while (bvhStack.size > 0) {
            RaycastMeshBvh* bvh = bvhStack[bvhStack.size - 1];
            bvhStack.RemoveLast();

            if (bvh->child1 == nullptr) {
                numBvhNodes++;
                totalBvhTriangles += bvh->triangles.size;
                if (bvh->triangles.size > maxBvhTriangles) {
                    maxBvhTrianglesInd = i;
                    maxBvhTriangles = bvh->triangles.size;
                }
            }
            else {
                bvhStack.Append(bvh->child1);
                bvhStack.Append(bvh->child2);
            }
        }

        totalBvhs += numBvhNodes;
        if (numBvhNodes > maxBvhNodes) {
            maxBvhNodesInd = i;
            maxBvhNodes = numBvhNodes;
        }
    }

    // Report BVH stats
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

    return geometry;
}

bool HitTriangles(Vec3 rayOrigin, Vec3 rayDir, float32 minDist, Array<RaycastTriangle> triangles,
                  uint32* hitMaterialIndex, Vec3* hitNormal, float32* hitDist)
{
    bool hit = false;

    for (uint32 i = 0; i < triangles.size; i++) {
        const RaycastTriangle& triangle = triangles[i];
        const float32 dotNegRayNormal = Dot(rayDir, triangle.normal);
        if (dotNegRayNormal > 0.0f) {
            continue;
        }

        float32 t;
        const bool tIntersect = RayTriangleIntersection(rayOrigin, rayDir,
                                                        triangle.pos[0], triangle.pos[1], triangle.pos[2], &t);
        UNREFERENCED_PARAMETER(minDist);
        if (tIntersect && t > minDist && t < *hitDist) {
            *hitMaterialIndex = triangle.materialIndex;
            *hitNormal = triangle.normal;
            *hitDist = t;
            hit = true;
        }
    }

    return hit;
}

const uint32 BVH_STACK_SIZE = 4096;

bool TraverseMeshBvh(Vec3 rayOrigin, Vec3 rayDir, Vec3 inverseRayDir, float32 minDist,
                     const RaycastMeshBvh* bvhStack[BVH_STACK_SIZE],
                     uint32* hitMaterialIndex, Vec3* hitNormal, float32* hitDist)
{
    uint32 n = 1;

    bool hit = false;
    while (n > 0) {
        const RaycastMeshBvh* bvh = bvhStack[--n];

        float32 tMin, tMax;
        const bool intersect = RayAxisAlignedBoxIntersection(rayOrigin, inverseRayDir, bvh->aabb, &tMin, &tMax);
        if (!intersect || (tMin < 0.0f && tMax < 0.0f) || tMin > *hitDist) {
            continue;
        }

        if (bvh->child1 == nullptr) {
            if (HitTriangles(rayOrigin, rayDir, minDist, bvh->triangles, hitMaterialIndex, hitNormal, hitDist)) {
                hit = true;
            }
        }
        else {
            bvhStack[n++] = bvh->child1;
            bvhStack[n++] = bvh->child2;
        }
    }

    return hit;
}

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, uint32 bounces, float32 minDist, float32 maxDist,
                  const RaycastMeshBvh* bvhStack[BVH_STACK_SIZE], const RaycastGeometry& geometry, RandomSeries* series)
{
    ZoneScoped;

    float32 intensity = 1.0f;
    Vec3 color = Vec3::zero;

    for (uint32 b = 0; b < bounces; b++) {
        // TODO I think this is causing artifacts, because we could be dividing by zero
        const Vec3 inverseRayDir = Reciprocal(rayDir);

        uint32 hitMaterialIndex = geometry.materials.size;
        Vec3 hitNormal = Vec3::zero;
        float32 hitDist = maxDist;
        for (uint32 i = 0; i < geometry.meshes.size; i++) {
            const RaycastMesh& mesh = geometry.meshes[i];
            bvhStack[0] = &mesh.bvh;
            TraverseMeshBvh(rayOrigin, rayDir, inverseRayDir, minDist, bvhStack,
                            &hitMaterialIndex, &hitNormal, &hitDist);
        }

        if (hitMaterialIndex == geometry.materials.size) {
            break;
        }

        const RaycastMaterial& hitMaterial = geometry.materials[hitMaterialIndex];
        if (hitMaterial.emission > 0.0f) {
            color = hitMaterial.emission * hitMaterial.emissionColor * intensity;
            break;
        }
        else {
            intensity *= 0.5f;
            rayOrigin += rayDir * hitDist;

            const Quat xToNormal = QuatRotBetweenVectors(Vec3::unitX, hitNormal);
            const Vec3 pureBounce = rayDir - 2.0f * Dot(rayDir, hitNormal) * hitNormal;
            const Vec3 randomBounce = xToNormal * NormalizeOrZero(Vec3 {
                                                                      RandomUnilateral(series),
                                                                      RandomBilateral(series),
                                                                      RandomBilateral(series),
                                                                  });
            rayDir = NormalizeOrZero(Lerp(randomBounce, pureBounce, hitMaterial.smoothness));
        }
    }

    return color;
}

struct RaycastThreadWorkCommon
{
    const RaycastGeometry* geometry;
    Vec3 filmTopLeft;
    Vec3 filmUnitOffsetX;
    Vec3 filmUnitOffsetY;
    Vec3 cameraPos;
    uint32 width, height;
    uint32 bounces;
    float32 fill;
    float32 minDist;
    float32 maxDist;
};

struct RaycastThreadWork
{
    const RaycastThreadWorkCommon* common;
    Vec3* colorHdr;
    RandomSeries randomSeries;
    uint32 minX, maxX, minY, maxY;
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(threadIndex);
    UNREFERENCED_PARAMETER(queue);
    ZoneScoped;

    // TODO have some guarantee that this stack will be big enough
    const RaycastMeshBvh* bvhStack[BVH_STACK_SIZE];

    RaycastThreadWork* work = (RaycastThreadWork*)data;
    const RaycastThreadWorkCommon& common = *work->common;
    RandomSeries randomSeries = work->randomSeries;
    const uint32 minX = work->minX;
    const uint32 maxX = work->maxX;
    const uint32 minY = work->minY;
    const uint32 maxY = work->maxY;

    for (uint32 y = minY; y < maxY; y++) {
        for (uint32 x = minX; x < maxX; x++) {
            if (RandomUnilateral(&randomSeries) < common.fill) {
                const Vec3 filmOffsetX = common.filmUnitOffsetX * (float32)x;
                const Vec3 filmOffsetY = common.filmUnitOffsetY * (float32)y;
                const Vec3 filmPos = common.filmTopLeft + filmOffsetX + filmOffsetY;
                const Vec3 rayDir = Normalize(filmPos - common.cameraPos);

                const Vec3 raycastColor = RaycastColor(common.cameraPos, rayDir, common.bounces,
                                                       common.minDist, common.maxDist,
                                                       bvhStack, *common.geometry, &randomSeries);

                const uint32 pixelIndex = y * common.width + x;
                work->colorHdr[pixelIndex] = raycastColor;
            }
        }
    }
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    const uint8* materialIndices, uint32 width, uint32 height, CanvasState* canvas, uint32* pixels,
                    LinearAllocator* allocator, AppWorkQueue* queue)
{
    UNREFERENCED_PARAMETER(materialIndices);
    ZoneScoped;

    const uint32 numPixels = width * height;

    {
        ZoneScopedN("PixelDecay");

#if 0
        for (uint32 i = 0; i < numPixels; i++) {
            uint8* decay = &canvas->decay.data[i];
            if (*decay != 0xff) {
                (*decay) += 1;
                if (*decay >= canvas->decayFrames) {
                    canvas->colorHdr[i] = Vec3::zero;
                    *decay = 0xff;
                }
            }
        }
#endif

        MemSet(canvas->colorHdr.data, 0, numPixels * sizeof(Vec3));
    }

    TracyCZoneN(zoneProduceWork, "ProduceWork", true);

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
        .filmTopLeft = filmTopLeft,
        .filmUnitOffsetX = filmUnitOffsetX,
        .filmUnitOffsetY = filmUnitOffsetY,
        .cameraPos = cameraPos,
        .width = width,
        .height = height,
        .bounces = canvas->bounces,
        .fill = canvas->screenFill,
        .minDist = 0.0f,
        .maxDist = 20.0f,
    };

    RandomSeries series;
    const uint32 seed = (uint32)(cameraPos.x * 1000.0f + cameraPos.y * 1000.0f + cameraPos.z * 1000.0f);
    series.state = seed;
    series.state = (uint32)rand();

    const uint32 WORK_TILE_SIZE = 16;
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
            work->randomSeries.state = RandomUInt32(&series);
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

    TracyCZoneEnd(zoneProduceWork);

    {
        ZoneScopedN("FinishWork");
        CompleteAllWork(queue, 0);
    }

    TracyCZoneN(zoneHdrTranslate, "TranslateHdr", true);

    float32 maxColor = 1.0f;
    for (uint32 i = 0; i < numPixels; i++) {
        const Vec3 colorHdr = canvas->colorHdr[i];
        maxColor = MaxFloat32(maxColor, colorHdr.r);
        maxColor = MaxFloat32(maxColor, colorHdr.g);
        maxColor = MaxFloat32(maxColor, colorHdr.b);
    }

    for (uint32 i = 0; i < numPixels; i++) {
        const Vec3 colorNormalized = canvas->colorHdr[i] / maxColor;
        // TODO gamma correction?
        const uint8 r = (uint8)(colorNormalized.r * 255.0f);
        const uint8 g = (uint8)(colorNormalized.g * 255.0f);
        const uint8 b = (uint8)(colorNormalized.b * 255.0f);
        pixels[i] = ((uint32)0xff << 24) + ((uint32)b << 16) + ((uint32)g << 8) + r;
    }

    TracyCZoneEnd(zoneHdrTranslate);
}

uint32 CountSubBvhs(const RaycastMeshBvh* root, LinearAllocator* allocator)
{
    SCOPED_ALLOCATOR_RESET(*allocator);

    uint32 count = 0;
    DynamicArray<const RaycastMeshBvh*, LinearAllocator> bvhStack(allocator);
    bvhStack.Append(root);

    while (bvhStack.size > 0) {
        const RaycastMeshBvh* bvh = bvhStack[bvhStack.size - 1];
        bvhStack.RemoveLast();
        count++;

        if (bvh->child1 != nullptr) {
            bvhStack.Append(bvh->child1);
            bvhStack.Append(bvh->child2);
        }
    }

    return count - 1;
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

    // triangles, BVHs, meshes
    {
        SCOPED_ALLOCATOR_RESET(*allocator);

        DynamicArray<ComputeTriangle, LinearAllocator> triangles(allocator);
        DynamicArray<ComputeBvh, LinearAllocator> bvhs(allocator);
        pipeline->meshes.Clear();
        DynamicArray<const RaycastMeshBvh*, LinearAllocator> bvhStack(allocator);

        for (uint32 i = 0; i < geometry.meshes.size; i++) {
            const RaycastMesh& mesh = geometry.meshes[i];

            ComputeMesh* newMesh = pipeline->meshes.Append();
            newMesh->offset = Vec3::zero;
            newMesh->quat = Vec4 {
                Quat::one.x, Quat::one.y, Quat::one.z, Quat::one.w
            };
            newMesh->startBvh = bvhs.size;

            bvhStack.Clear();
            bvhStack.Append(&mesh.bvh);
            while (bvhStack.size > 0) {
                const RaycastMeshBvh* bvh = bvhStack[bvhStack.size - 1];
                bvhStack.RemoveLast();

                ComputeBvh* newBvh = bvhs.Append();
                newBvh->aabbMin = bvh->aabb.min;
                newBvh->aabbMax = bvh->aabb.max;
                newBvh->skip = CountSubBvhs(bvh, allocator) + 1;

                if (bvh->child1 == nullptr) {
                    newBvh->startTriangle = triangles.size;

                    for (uint32 j = 0; j < bvh->triangles.size; j++) {
                        const RaycastTriangle& srcTriangle = bvh->triangles[j];
                        ComputeTriangle* dstTriangle = triangles.Append();
                        dstTriangle->a = srcTriangle.pos[0];
                        dstTriangle->b = srcTriangle.pos[1];
                        dstTriangle->c = srcTriangle.pos[2];
                        dstTriangle->normal = srcTriangle.normal;
                        dstTriangle->materialIndex = srcTriangle.materialIndex;
                    }

                    newBvh->endTriangle = triangles.size;
                }
                else {
                    newBvh->startTriangle = 0;
                    newBvh->endTriangle = 0;

                    // pre-order traversal
                    bvhStack.Append(bvh->child2);
                    bvhStack.Append(bvh->child1);
                }
            }

            newMesh->endBvh = bvhs.size;
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
        fenceCreateInfo.flags = 0;//VK_FENCE_CREATE_SIGNALED_BIT;

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
