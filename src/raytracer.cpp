#include "raytracer.h"

#include <intrin.h>
#include <time.h>

#include <stb_image_write.h>
#ifdef TRACY_ENABLE
#include <Tracy.hpp>
#endif

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
        const_string name = obj.materials[i].name;
        RaycastMaterial& material = geometry.materials[i];

        if (StringEquals(name, ToString("Surfaces")) || StringEquals(name, ToString("None"))) {
            material.smoothness = 0.8f;
            material.albedo = Vec3 { 1.0f, 1.0f, 1.0f };
            material.emission = 0.0f;
            material.emissionColor = Vec3::zero;
        }
        else if (StringEquals(name, ToString("LightHardGreen"))) {
            material.smoothness = 0.0f;
            material.albedo = Vec3 { 0.0f, 1.0f, 0.42f };
            material.emission = 1.5f;
            material.emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };
        }
        else if (StringEquals(name, ToString("LightSoftGreen"))) {
            material.smoothness = 0.0f;
            material.albedo = Vec3 { 0.0f, 1.0f, 0.42f };
            material.emission = 1.0f;
            material.emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };
        }
        else if (StringEquals(name, ToString("LightPink"))) {
            material.smoothness = 0.0f;
            material.albedo = Vec3 { 0.9294f, 0.651f, 1.0f };
            material.emission = 1.0f;
            material.emissionColor = Vec3 { 0.9294f, 0.651f, 1.0f };
        }
        else if (StringEquals(name, ToString("LightRed"))) {
            material.smoothness = 0.0f;
            material.albedo = Vec3 { 1.0f, 0.0f, 0.082f };
            material.emission = 1.5f;
            material.emissionColor = Vec3 { 1.0f, 0.0f, 0.082f };
        }
        else {
            LOG_ERROR("Unrecognized material: %.*s\n", name.size, name.data);
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

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, float32 minDist, float32 maxDist, const RaycastGeometry& geometry)
{
    ZoneScoped;
    const uint32 SAMPLES = 8;
    const uint32 BOUNCES = 4;

    // TODO have some guarantee that this stack will be big enough
    const RaycastMeshBvh* bvhStack[4096];

    const float32 sampleWeight = 1.0f / (float32)SAMPLES;

    Vec3 color = Vec3::zero;
    for (uint32 s = 0; s < SAMPLES; s++) {
        for (uint32 b = 0; b < BOUNCES; b++) {
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
                color += sampleWeight * hitMaterial.emission * hitMaterial.emissionColor;
                break;
            }
            else {
                rayOrigin += rayDir * hitDist;

                const Quat xToNormal = QuatRotBetweenVectors(Vec3::unitX, hitNormal);
                const Vec3 pureBounce = rayDir - 2.0f * Dot(rayDir, hitNormal) * hitNormal;
                const Vec3 randomBounce = xToNormal * NormalizeOrZero(Vec3 {
                                                                          RandFloat32( 0.0f, 1.0f),
                                                                          RandFloat32(-1.0f, 1.0f),
                                                                          RandFloat32(-1.0f, 1.0f)
                                                                      });
                rayDir = NormalizeOrZero(Lerp(randomBounce, pureBounce, hitMaterial.smoothness));
            }
        }
    }

    return color;
}

struct RaycastThreadWorkCommon
{
    Vec3 filmTopLeft;
    Vec3 filmUnitOffsetX;
    Vec3 filmUnitOffsetY;
    Vec3 cameraPos;
    float32 minDist;
    float32 maxDist;
};

struct RaycastThreadWork
{
    static const uint32 PIXELS_PER_WORK_UNIT = 64;

    const RaycastThreadWorkCommon* common;
    const RaycastGeometry* geometry;
    Vec2Int pixels[PIXELS_PER_WORK_UNIT];
    Vec3 outputColors[PIXELS_PER_WORK_UNIT];
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(threadIndex);
    UNREFERENCED_PARAMETER(queue);

    RaycastThreadWork* work = (RaycastThreadWork*)data;
    for (uint32 i = 0; i < work->PIXELS_PER_WORK_UNIT; i++) {
        const Vec3 filmOffsetX = work->common->filmUnitOffsetX * (float32)work->pixels[i].x;
        const Vec3 filmOffsetY = work->common->filmUnitOffsetY * (float32)work->pixels[i].y;
        const Vec3 filmPos = work->common->filmTopLeft + filmOffsetX + filmOffsetY;
        const Vec3 rayDir = Normalize(filmPos - work->common->cameraPos);
        work->outputColors[i] = RaycastColor(work->common->cameraPos, rayDir, 0.0f, 20.0f, *work->geometry);
    }
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, CanvasState* canvas, LinearAllocator* allocator, AppWorkQueue* queue)
{
    ZoneScoped;
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

    const Quat inverseCameraRot = Inverse(cameraRot);
    const Vec3 cameraUp = inverseCameraRot * Vec3::unitZ;
    const Vec3 cameraForward = inverseCameraRot * Vec3::unitX;
    const Vec3 cameraLeft = inverseCameraRot * Vec3::unitY;

    const float32 filmDist = 1.0f;
    const float32 filmWidth = 2.0f;
    const float32 filmHeight = filmWidth * (float32)height / (float32)width;

    const Vec3 filmTopLeft = cameraPos + cameraForward * filmDist
        + filmWidth / 2.0f * cameraUp + filmHeight / 2.0f * cameraLeft;
    const Vec3 filmUnitOffsetX = -cameraLeft * filmWidth / (float32)width;
    const Vec3 filmUnitOffsetY = -cameraUp * filmHeight / (float32)height;

    const RaycastThreadWorkCommon workCommon = {
        .filmTopLeft = filmTopLeft,
        .filmUnitOffsetX = filmUnitOffsetX,
        .filmUnitOffsetY = filmUnitOffsetY,
        .cameraPos = cameraPos,
        .minDist = 0.0f,
        .maxDist = 20.0f
    };

    // The idea is that, over 1 second, we'll get a whole frame-pixels' worth of rays
    uint32 NUM_RAYS_PER_FRAME = (uint32)(canvas->screenFill * (float32)width * (float32)height);
    NUM_RAYS_PER_FRAME = RoundUpToPowerOfTwo(NUM_RAYS_PER_FRAME, RaycastThreadWork::PIXELS_PER_WORK_UNIT);

    const uint32 numWorkEntries = NUM_RAYS_PER_FRAME / RaycastThreadWork::PIXELS_PER_WORK_UNIT;
    Array<RaycastThreadWork> workEntries = allocator->NewArray<RaycastThreadWork>(numWorkEntries);
    for (uint32 i = 0; i < workEntries.size; i++) {
        workEntries[i].common = &workCommon;
        workEntries[i].geometry = &geometry;
        for (uint32 j = 0; j < RaycastThreadWork::PIXELS_PER_WORK_UNIT; j++) {
            workEntries[i].pixels[j].x = RandInt(width);
            workEntries[i].pixels[j].y = RandInt(height);
        }

        if (!TryAddWork(queue, &RaycastThreadProc, &workEntries[i])) {
            CompleteAllWork(queue, 0);
            i--;
            continue;
        }
    }

    CompleteAllWork(queue, 0);

    const float32 neighborIntensity = 0.5f;
    for (uint32 i = 0; i < workEntries.size; i++) {
        for (uint32 j = 0; j < RaycastThreadWork::PIXELS_PER_WORK_UNIT; j++) {
            const uint32 pixelInd = workEntries[i].pixels[j].y * width + workEntries[i].pixels[j].x;
            canvas->colorHdr[pixelInd] = workEntries[i].outputColors[j];
            canvas->decay[pixelInd] = 0;
        }
    }
}
