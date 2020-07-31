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

void FillRaycastMeshBoxRecursive(RaycastMeshBox* box, Box minBox, uint32 boxMaxTriangles, Array<RaycastTriangle> triangles, 
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
        FillRaycastMeshBoxRecursive(box->child1, minBox1, boxMaxTriangles, insideTriangles.ToArray(),
                                    allocator, tempAllocator);

        box->child2 = allocator->New<RaycastMeshBox>();
        FillRaycastMeshBoxRecursive(box->child2, minBox2, boxMaxTriangles, insideTriangles.ToArray(),
                                    allocator, tempAllocator);
    }
}

void FillRaycastMeshBox(RaycastMeshBox* box, uint32 boxMaxTriangles, Array<RaycastTriangle> triangles,
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
        FillRaycastMeshBoxRecursive(box->child1, minBox1, boxMaxTriangles, triangles, allocator, tempAllocator);

        box->child2 = allocator->New<RaycastMeshBox>();
        FillRaycastMeshBoxRecursive(box->child2, minBox2, boxMaxTriangles, triangles, allocator, tempAllocator);
    }
}

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, uint32 boxMaxTriangles,
                                      LinearAllocator* allocator, LinearAllocator* tempAllocator)
{
    RaycastGeometry geometry = {};

    geometry.materials = allocator->NewArray<RaycastMaterial>(4);
    if (geometry.materials.data == nullptr) {
        return geometry;
    }

    // Light cubes
    geometry.materials[0].smoothness = 0.0f;
    geometry.materials[0].albedo = Vec3 { 1.0f, 0.0f, 0.082f };
    geometry.materials[0].emission = 4.0f;
    geometry.materials[0].emissionColor = Vec3 { 1.0f, 0.0f, 0.0f };

    // Light cone
    geometry.materials[1].smoothness = 0.0f;
    geometry.materials[1].albedo = Vec3 { 0.0f, 1.0f, 0.42f };
    geometry.materials[1].emission = 4.0f;
    geometry.materials[1].emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };

    // Light disk (overhead)
    geometry.materials[2].smoothness = 0.0f;
    geometry.materials[2].albedo = Vec3 { 0.9294f, 0.651f, 1.0f };
    geometry.materials[2].emission = 1.0f;
    geometry.materials[2].emissionColor = Vec3 { 0.9294f, 0.651f, 1.0f };

    // Surfaces (everything else)
    geometry.materials[3].smoothness = 0.5f;
    geometry.materials[3].albedo = Vec3 { 1.0f, 1.0f, 1.0f };
    geometry.materials[3].emission = 0.0f;
    geometry.materials[3].emissionColor = Vec3::zero;

    geometry.meshes = allocator->NewArray<RaycastMesh>(obj.models.size);
    if (geometry.meshes.data == nullptr) {
        return geometry;
    }

    const Vec3 vertexColor = Vec3::zero;

    for (uint32 i = 0; i < obj.models.size; i++) {
        RaycastMesh& mesh = geometry.meshes[i];

        // For light cones min scene
        if (i == 2 || i == 3 || i == 4) { // cubes
            mesh.materialIndex = 0;
        }
        else if (i == obj.models.size - 3) { // light cone
            mesh.materialIndex = 1;
        }
        else if (i == obj.models.size - 1) { // light overhead
            mesh.materialIndex = 2;
        }
        else { // everything else
            mesh.materialIndex = 3;
        }

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
        }
        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const uint32 ind = obj.models[i].triangles.size + j * 2;
            const ObjQuad& q = obj.models[i].quads[j];
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                triangles[ind].pos[k] = q.v[k].pos;
            }
            triangles[ind].normal = normal;

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                triangles[ind + 1].pos[k] = q.v[quadInd].pos;
            }
            triangles[ind + 1].normal = normal;
        }

        FillRaycastMeshBox(&mesh.box, boxMaxTriangles, triangles, allocator, tempAllocator);
    }

    return geometry;
}

bool HitTriangles(Vec3 rayOrigin, Vec3 rayDir, Array<RaycastTriangle> triangles, Vec3* hitNormal, float32* hitDist)
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
        if (tIntersect && t > 0.0f && t < *hitDist) {
            *hitNormal = triangle.normal;
            *hitDist = t;
            hit = true;
        }
    }

    return hit;
}

bool TraverseMeshBox(Vec3 rayOrigin, Vec3 rayDir, Vec3 inverseRayDir, const RaycastMeshBox& box,
                     Vec3* hitNormal, float32* hitDist)
{
    float32 t;
    bool intersect = RayAxisAlignedBoxIntersection(rayOrigin, inverseRayDir, box.aabb, &t);
    if (!intersect) {
        return false;
    }

    if (box.child1 == nullptr) {
        return HitTriangles(rayOrigin, rayDir, box.triangles, hitNormal, hitDist);
    }

    bool hit = false;

    // TODO remove this recursion? need thread-safe allocators though
    if (TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, *box.child1, hitNormal, hitDist)) {
        hit = true;
    }
    if (TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, *box.child2, hitNormal, hitDist)) {
        hit = true;
    }

    return hit;
}

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, const RaycastGeometry& geometry)
{
    const uint32 SAMPLES = 8;
    const uint32 BOUNCES = 4;

    const float32 sampleWeight = 1.0f / (float32)SAMPLES;

    Vec3 color = Vec3::zero;
    for (uint32 s = 0; s < SAMPLES; s++) {
        for (uint32 b = 0; b < BOUNCES; b++) {
            Vec3 inverseRayDir = Reciprocal(rayDir);

            uint32 hitMaterialInd = geometry.materials.size;
            Vec3 hitNormal = Vec3::zero;
            float32 hitDist = 1e8;
            for (uint32 i = 0; i < geometry.meshes.size; i++) {
                const RaycastMesh& mesh = geometry.meshes[i];

#if 1
                if (TraverseMeshBox(rayOrigin, rayDir, inverseRayDir, mesh.box, &hitNormal, &hitDist)) {
                    hitMaterialInd = mesh.materialIndex;
                }
#else
                float32 tAabb;
                bool intersect = RayAxisAlignedBoxIntersection(rayOrigin, inverseRayDir, mesh.box.aabb, &tAabb);
                if (!intersect) {
                    continue;
                }

                if (HitTriangles(rayOrigin, rayDir, mesh.triangles, &hitNormal, &hitDist)) {
                    hitMaterialInd = mesh.materialIndex;
                }
#endif
            }

            if (hitMaterialInd == geometry.materials.size) {
                break;
            }

            const RaycastMaterial& hitMaterial = geometry.materials[hitMaterialInd];
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

    work->outputColor = RaycastColor(work->cameraPos, work->rayDir, *work->geometry);
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
