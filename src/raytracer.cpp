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

struct DebugTimer
{
    static bool initialized;
    static uint64 win32Freq;

    uint64 cycles;
    uint64 win32Time;
};

bool DebugTimer::initialized = false;
uint64 DebugTimer::win32Freq;

DebugTimer StartDebugTimer()
{
    if (!DebugTimer::initialized) {
        DebugTimer::initialized = true;
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        DebugTimer::win32Freq = freq.QuadPart;
    }

    DebugTimer timer;
    LARGE_INTEGER win32Time;
    QueryPerformanceCounter(&win32Time);
    timer.win32Time = win32Time.QuadPart;
    timer.cycles = __rdtsc();
    return timer;
}

void StopDebugTimer(DebugTimer* timer)
{
    LARGE_INTEGER win32End;
    QueryPerformanceCounter(&win32End);
    timer->cycles = __rdtsc() - timer->cycles;
    timer->win32Time = win32End.QuadPart - timer->win32Time;
}

void StopAndPrintDebugTimer(DebugTimer* timer)
{
    StopDebugTimer(timer);
    const float32 win32Time = (float32)timer->win32Time / DebugTimer::win32Freq * 1000.0f;
    LOG_INFO("Timer: %.03fms | %llu MC\n", win32Time, timer->cycles / 1000000);
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

bool FillRaycastMeshBvh(RaycastMeshBvh* bvh, Box minBox, uint32 boxMaxTriangles, bool firstPass, 
                        Array<RaycastTriangle> triangles, LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    DEBUG_ASSERT(bvh != nullptr);

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
    if (insideTriangles.size <= boxMaxTriangles || (!firstPass && insideTriangles.size == triangles.size)) {
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
        if (!FillRaycastMeshBvh(bvh->child1, minBox1, boxMaxTriangles, false, insideTriangles.ToArray(),
                                allocator, tempAllocator)) {
            return false;
        }

        bvh->child2 = allocator->New<RaycastMeshBvh>();
        if (bvh->child2 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBvh(bvh->child2, minBox2, boxMaxTriangles, false, insideTriangles.ToArray(),
                                allocator, tempAllocator)) {
            return false;
        }
    }

    return true;
}

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, uint32 boxMaxTriangles,
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

    const Vec3 vertexColor = Vec3::zero;

    for (uint32 i = 0; i < obj.models.size; i++) {
        RaycastMesh& mesh = geometry.meshes[i];

        const uint32 numTriangles = obj.models[i].triangles.size + obj.models[i].quads.size * 2;
        mesh.numTriangles = numTriangles;

        Array<RaycastTriangle> triangles = tempAllocator->NewArray<RaycastTriangle>(numTriangles);
        if (triangles.data == nullptr) {
            LOG_ERROR("Failed to allocate triangles for raycast mesh %lu\n", i);
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
        if (!FillRaycastMeshBvh(&mesh.bvh, veryBigBox, boxMaxTriangles, true, triangles, allocator, tempAllocator)) {
            LOG_ERROR("Failed to fill raycast mesh box for mesh %lu\n", i);
            geometry.meshes.data = nullptr;
            return geometry;
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
                                                        triangle.pos[0], triangle.pos[1], triangle.pos[2],
                                                        &t);
        if (tIntersect && t > minDist && t < *hitDist) {
            *hitMaterialIndex = triangle.materialIndex;
            *hitNormal = triangle.normal;
            *hitDist = t;
            hit = true;
        }
    }

    return hit;
}

bool TraverseMeshBox(Vec3 rayOrigin, Vec3 rayDir, Vec3 inverseRayDir, float32 minDist,
                     const RaycastMeshBvh* bvhStack[], uint32* hitMaterialIndex, Vec3* hitNormal, float32* hitDist)
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
                  const RaycastGeometry& geometry, RandomSeries* series)
{
    ZoneScoped;

    // TODO have some guarantee that this stack will be big enough
    const RaycastMeshBvh* bvhStack[4096];

    Vec3 color = Vec3::zero;
    for (uint32 b = 0; b < bounces; b++) {
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
            color = hitMaterial.emission * hitMaterial.emissionColor;
            break;
        }
        else {
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
    uint32 bounces;
    float32 minDist;
    float32 maxDist;
};

struct RaycastThreadWorkState
{
    Array<RandomSeries> threadRandomSeries;
};

struct RaycastThreadWork
{
    static const uint32 PIXELS_PER_WORK_UNIT = 64;

    const RaycastThreadWorkCommon* common;
    RaycastThreadWorkState* state;
    Vec2Int pixels[PIXELS_PER_WORK_UNIT];
    Vec3 outputColors[PIXELS_PER_WORK_UNIT];
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(queue);

    RaycastThreadWork* work = (RaycastThreadWork*)data;
    for (uint32 i = 0; i < work->PIXELS_PER_WORK_UNIT; i++) {
        const Vec3 filmOffsetX = work->common->filmUnitOffsetX * (float32)work->pixels[i].x;
        const Vec3 filmOffsetY = work->common->filmUnitOffsetY * (float32)work->pixels[i].y;
        const Vec3 filmPos = work->common->filmTopLeft + filmOffsetX + filmOffsetY;
        const Vec3 rayDir = Normalize(filmPos - work->common->cameraPos);
        work->outputColors[i] = RaycastColor(work->common->cameraPos, rayDir, work->common->bounces,
                                             work->common->minDist, work->common->maxDist,
                                             *(work->common->geometry), &work->state->threadRandomSeries[threadIndex]);
    }
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    const uint8* materialIndices, uint32 width, uint32 height, CanvasState* canvas,
                    LinearAllocator* allocator, AppWorkQueue* queue)
{
    UNREFERENCED_PARAMETER(materialIndices);
    ZoneScoped;

    TracyCZoneN(zonePixelDecay, "PixelDecay", true);

#if 0
    for (uint32 i = 0; i < width * height; i++) {
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
    MemSet(canvas->colorHdr.data, 0, width * height * sizeof(Vec3));

    TracyCZoneEnd(zonePixelDecay);

    TracyCZoneN(zoneProduceWork, "ProduceWork", true);

    const Quat inverseCameraRot = Inverse(cameraRot);
    const Vec3 cameraUp = inverseCameraRot * Vec3::unitZ;
    const Vec3 cameraForward = inverseCameraRot * Vec3::unitX;
    const Vec3 cameraLeft = inverseCameraRot * Vec3::unitY;

    const float32 filmDist = 1.0f;
    const float32 filmHeight = tanf(fov / 2.0f) * 2.0f;
    const float32 filmWidth = filmHeight * (float32)width / (float32)height;

    const Vec3 filmTopLeft = cameraPos + cameraForward * filmDist
        + (filmWidth / 2.0f + canvas->test2) * cameraUp + (filmHeight / 2.0f + canvas->test3) * cameraLeft;
    const Vec3 filmUnitOffsetX = -cameraLeft * filmWidth / (float32)(width - 1);
    const Vec3 filmUnitOffsetY = -cameraUp * filmHeight / (float32)(height - 1);

    const RaycastThreadWorkCommon workCommon = {
        .geometry = &geometry,
        .filmTopLeft = filmTopLeft,
        .filmUnitOffsetX = filmUnitOffsetX,
        .filmUnitOffsetY = filmUnitOffsetY,
        .cameraPos = cameraPos,
        .bounces = canvas->bounces,
        .minDist = 0.0f,
        .maxDist = 20.0f,
    };

    RandomSeries series;
    const uint32 seed = (uint32)(cameraPos.x * 1000.0f + cameraPos.y * 1000.0f + cameraPos.z * 1000.0f);
    series.state = seed;
    // series.state = (uint32)rand();

    const uint32 numThreads = GetCpuCount();
    RaycastThreadWorkState workState = {};
    workState.threadRandomSeries = allocator->NewArray<RandomSeries>(numThreads);
    for (uint32 i = 0; i < workState.threadRandomSeries.size; i++) {
        // TODO these don't guarantee same results, because threads will pick work entries in undefined order
        // Need a more specific, dedicated thread pool system for this task
        workState.threadRandomSeries[i].state = RandomUInt32(&series);
    }

    uint32 NUM_RAYS_PER_FRAME = (uint32)(canvas->screenFill * (float32)width * (float32)height);
    NUM_RAYS_PER_FRAME = RoundUpToPowerOfTwo(NUM_RAYS_PER_FRAME, RaycastThreadWork::PIXELS_PER_WORK_UNIT);

    const uint32 numWorkEntries = NUM_RAYS_PER_FRAME / RaycastThreadWork::PIXELS_PER_WORK_UNIT;
    Array<RaycastThreadWork> workEntries = allocator->NewArray<RaycastThreadWork>(numWorkEntries);
    for (uint32 i = 0; i < workEntries.size; i++) {
        workEntries[i].common = &workCommon;
        workEntries[i].state = &workState;
        for (uint32 j = 0; j < RaycastThreadWork::PIXELS_PER_WORK_UNIT; j++) {
            workEntries[i].pixels[j].x = RandomInt32(&series, width);
            workEntries[i].pixels[j].y = RandomInt32(&series, height);
        }

        if (!TryAddWork(queue, &RaycastThreadProc, &workEntries[i])) {
            CompleteAllWork(queue, 0);
            i--;
            continue;
        }
    }

    TracyCZoneEnd(zoneProduceWork);

    {
        ZoneScopedN("FinishWork");
        CompleteAllWork(queue, 0);
    }

    TracyCZoneN(zoneBlend, "Blend", true);

    const float32 neighborIntensity = 0.5f;
    const float32 prevWeight = 0.5f;
    for (uint32 i = 0; i < workEntries.size; i++) {
        for (uint32 j = 0; j < RaycastThreadWork::PIXELS_PER_WORK_UNIT; j++) {
            const Vec2Int pixel = workEntries[i].pixels[j];
            const Vec3 color = workEntries[i].outputColors[j];

            const uint32 index = pixel.y * width + pixel.x;
            const Vec3 prevColor = canvas->colorHdr[index];
            canvas->colorHdr[index] = Lerp(color, prevColor + color, prevWeight);
            canvas->decay[index] = 0;

            if (pixel.x > 0) {
                const uint32 indexLeft = index - 1;
                canvas->colorHdr[indexLeft] += color * neighborIntensity;
                canvas->decay[indexLeft] = 0;
            }
            if ((uint32)pixel.x < width - 1) {
                const uint32 indexRight = index + 1;
                canvas->colorHdr[indexRight] += color * neighborIntensity;
                canvas->decay[indexRight] = 0;
            }
            if (pixel.y > 0) {
                const uint32 indexTop = index - width;
                canvas->colorHdr[indexTop] += color * neighborIntensity;
                canvas->decay[indexTop] = 0;
            }
            if ((uint32)pixel.y < height - 1) {
                const uint32 indexBottom = index + width;
                canvas->colorHdr[indexBottom] += color * neighborIntensity;
                canvas->decay[indexBottom] = 0;
            }
        }
    }

    TracyCZoneEnd(zoneBlend);
}
