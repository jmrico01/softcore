#include "raytracer.h"

#include <intrin.h>
#include <stb_image_write.h>
#include <time.h>

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

bool FillRaycastMeshBoxRecursive(RaycastMeshBox* box, Box minBox, uint32 boxMaxTriangles, Array<RaycastTriangle> triangles, 
                                 LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    DEBUG_ASSERT(box != nullptr);

    DynamicArray<RaycastTriangle, LinearAllocator> insideTriangles(tempAllocator);

    box->aabb.min = Vec3::one * 1e8;
    box->aabb.max = -Vec3::one * 1e8;
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
                    box->aabb.min.e[e] = MinFloat32(box->aabb.min.e[e], v.e[e]);
                    box->aabb.max.e[e] = MaxFloat32(box->aabb.max.e[e], v.e[e]);
                }
            }

            insideTriangles.Append(triangles[i]);
        }
    }

    // TODO examine the second case. # of triangles should decrease reasonably each partition,
    // otherwise we're duplicating checks
    if (insideTriangles.size <= boxMaxTriangles || insideTriangles.size == triangles.size) {
        box->child1 = nullptr;
        box->child2 = nullptr;
        box->triangles = allocator->NewArray<RaycastTriangle>(insideTriangles.size);
        if (box->triangles.data == nullptr) {
            return false;
        }
        box->triangles.CopyFrom(insideTriangles.ToArray());
    }
    else {
        const Vec3 boxSize = box->aabb.max - box->aabb.min;
        const int longestDimension = (boxSize.x >= boxSize.y && boxSize.x >= boxSize.z) ? 0
            : ((boxSize.y >= boxSize.x && boxSize.y >= boxSize.z) ? 1 : 2);
        const float32 halfLongestDimension = boxSize.e[longestDimension] / 2.0f;

        Box minBox1 = box->aabb;
        minBox1.max.e[longestDimension] -= halfLongestDimension;
        Box minBox2 = box->aabb;
        minBox2.min.e[longestDimension] += halfLongestDimension;

        box->child1 = allocator->New<RaycastMeshBox>();
        if (box->child1 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBoxRecursive(box->child1, minBox1, boxMaxTriangles, insideTriangles.ToArray(),
                                         allocator, tempAllocator)) {
            return false;
        }

        box->child2 = allocator->New<RaycastMeshBox>();
        if (box->child2 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBoxRecursive(box->child2, minBox2, boxMaxTriangles, insideTriangles.ToArray(),
                                         allocator, tempAllocator)) {
            return false;
        }
    }

    return true;
}

bool FillRaycastMeshBox(RaycastMeshBox* box, uint32 boxMaxTriangles, Array<RaycastTriangle> triangles,
                        LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    DEBUG_ASSERT(box != nullptr);

    // Calculate AABB
    box->aabb.min = Vec3::one * 1e8;
    box->aabb.max = -Vec3::one * 1e8;
    for (uint32 i = 0; i < triangles.size; i++) {
        for (int k = 0; k < 3; k++) {
            const Vec3 v = triangles[i].pos[k];
            for (int e = 0; e < 3; e++) {
                box->aabb.min.e[e] = MinFloat32(box->aabb.min.e[e], v.e[e]);
                box->aabb.max.e[e] = MaxFloat32(box->aabb.max.e[e], v.e[e]);
            }
        }
    }

    if (triangles.size <= boxMaxTriangles) {
        box->child1 = nullptr;
        box->child2 = nullptr;
        box->triangles = allocator->NewArray<RaycastTriangle>(triangles.size);
        if (box->triangles.data == nullptr) {
            return false;
        }
        box->triangles.CopyFrom(triangles);
    }
    else {
        const Vec3 boxSize = box->aabb.max - box->aabb.min;
        const int longestDimension = (boxSize.x >= boxSize.y && boxSize.x >= boxSize.z) ? 0
            : ((boxSize.y >= boxSize.x && boxSize.y >= boxSize.z) ? 1 : 2);
        const float32 halfLongestDimension = boxSize.e[longestDimension] / 2.0f;

        Box minBox1 = box->aabb;
        minBox1.max.e[longestDimension] -= halfLongestDimension;
        Box minBox2 = box->aabb;
        minBox2.min.e[longestDimension] += halfLongestDimension;

        box->child1 = allocator->New<RaycastMeshBox>();
        if (box->child1 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBoxRecursive(box->child1, minBox1, boxMaxTriangles, triangles, allocator, tempAllocator)) {
            return false;
        }

        box->child2 = allocator->New<RaycastMeshBox>();
        if (box->child2 == nullptr) {
            return false;
        }
        if (!FillRaycastMeshBoxRecursive(box->child2, minBox2, boxMaxTriangles, triangles, allocator, tempAllocator)) {
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

        if (!FillRaycastMeshBox(&mesh.box, boxMaxTriangles, triangles, allocator, tempAllocator)) {
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

bool TraverseMeshBox(Vec3 rayOrigin, Vec3 rayDir, Vec3 inverseRayDir, float32 minDist, float32 maxDist, 
                     const RaycastMeshBox& box, uint32* hitMaterialIndex, Vec3* hitNormal, float32* hitDist)
{
    // TODO get tMin and tMax out of AABB intersect, can discard rays by distance
    float32 tMin, tMax;
    bool intersect = RayAxisAlignedBoxIntersection(rayOrigin, inverseRayDir, box.aabb, &tMin, &tMax);
    if (!intersect || (tMin < 0.0f && tMax < 0.0f) || tMin > maxDist) {
        return false;
    }

    if (box.child1 == nullptr) {
        return HitTriangles(rayOrigin, rayDir, minDist, box.triangles, hitMaterialIndex, hitNormal, hitDist);
    }

    bool hit = false;

    // TODO remove this recursion? need thread-safe allocators though
    if (TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, minDist, maxDist, *box.child1,
                        hitMaterialIndex, hitNormal, hitDist)) {
        hit = true;
    }

    if (TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, minDist, maxDist, *box.child2,
                        hitMaterialIndex, hitNormal, hitDist)) {
        hit = true;
    }

    return hit;
}

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, float32 minDist, float32 maxDist, const RaycastGeometry& geometry)
{
    const uint32 SAMPLES = 8;
    const uint32 BOUNCES = 4;

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

                TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, minDist, maxDist, mesh.box,
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

struct RaycastThreadWork
{
    Vec3 cameraPos;
    Vec3 rayDir;
    const RaycastGeometry* geometry;
    Vec3 outputColor;
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(queue);

    RaycastThreadWork* work = (RaycastThreadWork*)data;

    work->outputColor = RaycastColor(work->cameraPos, work->rayDir, 0.0f, 20.0f, *work->geometry);
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, CanvasState* canvas, LinearAllocator* allocator, AppWorkQueue* queue)
{
    // TODO proper timing & profiling
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

    // The idea is that, over 1 second, we'll get a whole frame-pixels' worth of rays
    const uint32 NUM_RAYS_PER_FRAME = (uint32)(canvas->screenFill * (float32)width * (float32)height);
    Array<Vec2Int> randomPixels = allocator->NewArray<Vec2Int>(NUM_RAYS_PER_FRAME);
    for (uint32 i = 0; i < randomPixels.size; i++) {
        randomPixels[i].x = RandInt(width);
        randomPixels[i].y = RandInt(height);
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

    Array<RaycastThreadWork> workEntries = allocator->NewArray<RaycastThreadWork>(randomPixels.size);

    for (uint32 i = 0; i < randomPixels.size; i++) {
        const Vec3 filmOffsetX = -cameraLeft * (float32)randomPixels[i].x / (float32)width * filmWidth;
        const Vec3 filmOffsetY = -cameraUp * (float32)randomPixels[i].y / (float32)height * filmHeight;
        const Vec3 filmPos = filmTopLeft + filmOffsetX + filmOffsetY;
        const Vec3 rayDir = Normalize(filmPos - cameraPos);

        workEntries[i].cameraPos = cameraPos;
        workEntries[i].rayDir = rayDir;
        workEntries[i].geometry = &geometry;

        if (!TryAddWork(queue, &RaycastThreadProc, &workEntries[i])) {
            CompleteAllWork(queue);
            i--;
            continue;
        }
    }

    CompleteAllWork(queue);

    for (uint32 i = 0; i < randomPixels.size; i++) {
        const uint32 pixelInd = randomPixels[i].y * width + randomPixels[i].x;
        canvas->colorHdr[pixelInd] = workEntries[i].outputColor;
        canvas->decay[pixelInd] = 0;
    }
}
