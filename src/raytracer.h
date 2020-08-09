#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

const uint32 BOX_MAX_TRIANGLES = 32;

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

bool GetMaterial(const_string name, RaycastMaterial* material);

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, uint32 boxMaxTriangles,
                                      LinearAllocator* allocator, LinearAllocator* tempAllocator);

struct CanvasState
{
    static const uint32 MAX_WIDTH = 3840;
    static const uint32 MAX_HEIGHT = 2160;
    static const uint32 MAX_PIXELS = MAX_WIDTH * MAX_HEIGHT;

    uint32 prevSeed;

    float32 screenFill;
    uint8 decayFrames;
    uint32 bounces;

    float32 test1;
    float32 test2;
    float32 test3;

    StaticArray<Vec3, MAX_PIXELS> colorHdr;
    StaticArray<uint8, MAX_PIXELS> decay;
};

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    const uint8* materialIndices, uint32 width, uint32 height, CanvasState* canvas, uint32* pixels,
                    LinearAllocator* allocator, AppWorkQueue* queue);
