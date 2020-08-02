#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

struct RaycastMaterial
{
    float32 smoothness;
    Vec3 albedo;
    float32 emission;
    Vec3 emissionColor;
};

struct RaycastTriangle
{
    Vec3 pos[3];
    Vec3 normal;
    uint32 materialIndex;
};

struct RaycastMeshBvh
{
    Box aabb;
    RaycastMeshBvh* child1;
    RaycastMeshBvh* child2;
    Array<RaycastTriangle> triangles;
};

struct RaycastMesh
{
    RaycastMeshBvh bvh;
    uint32 numTriangles;
};

struct RaycastGeometry
{
    Array<RaycastMaterial> materials;
    Array<RaycastMesh> meshes;
};

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, uint32 boxMaxTriangles,
                                      LinearAllocator* allocator, LinearAllocator* tempAllocator);

struct CanvasState
{
    static const uint32 MAX_WIDTH = 3840;
    static const uint32 MAX_HEIGHT = 2160;
    static const uint32 MAX_PIXELS = MAX_WIDTH * MAX_HEIGHT;

    float32 screenFill;
    uint8 decayFrames;
    uint32 bounces;
    uint32 samples;

    StaticArray<Vec3, MAX_PIXELS> colorHdr;
    StaticArray<uint8, MAX_PIXELS> decay;
};

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, CanvasState* canvas, LinearAllocator* allocator, AppWorkQueue* queue);
