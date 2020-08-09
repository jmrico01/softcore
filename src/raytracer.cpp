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

bool GetMaterial(const_string name, RaycastMaterial* material)
{
    if (StringEquals(name, ToString("Surfaces"))) {
        material->smoothness = 0.8f;
        material->albedo = Vec3 { 1.0f, 1.0f, 1.0f };
        material->emission = 0.0f;
        material->emissionColor = Vec3::zero;
    }
    else if (StringEquals(name, ToString("LightHardGreen"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3 { 0.0f, 1.0f, 0.42f };
        material->emission = 1.5f;
        material->emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };
    }
    else if (StringEquals(name, ToString("LightSoftGreen"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3 { 0.0f, 1.0f, 0.42f };
        material->emission = 1.0f;
        material->emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };
    }
    else if (StringEquals(name, ToString("LightPink"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3 { 0.9294f, 0.651f, 1.0f };
        material->emission = 1.0f;
        material->emissionColor = Vec3 { 0.9294f, 0.651f, 1.0f };
    }
    else if (StringEquals(name, ToString("LightRed"))) {
        material->smoothness = 0.0f;
        material->albedo = Vec3 { 1.0f, 0.0f, 0.082f };
        material->emission = 1.5f;
        material->emissionColor = Vec3 { 1.0f, 0.0f, 0.082f };
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
            if (IsInside(triangles[i].pos[k], minBox)) {
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

    geometry.materials = allocator->NewArray<RaycastMaterial>(obj.materials.size);
    if (geometry.materials.data == nullptr) {
        return geometry;
    }

    for (uint32 i = 0; i < obj.materials.size; i++) {
        const_string materialName = obj.materials[i].name;
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

bool TraverseMeshBox(Vec3 rayOrigin, Vec3 rayDir, Vec3 inverseRayDir, float32 minDist,
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
            TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, minDist, bvhStack,
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
