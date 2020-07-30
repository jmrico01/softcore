#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

struct RaycastTriangle
{
    Vec3 pos[3];
    Vec3 normal;
};

struct RaycastMesh
{
    Box aabb;
    float32 smoothness;
    Vec3 albedo;
    float32 emission;
    Vec3 emissionColor;
    Array<RaycastTriangle> triangles;
};

struct RaycastGeometry
{
    Array<RaycastMesh> meshes;
};

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, LinearAllocator* allocator);

bool GenerateLightmaps(const LoadObjResult& obj, uint32 bounces, AppWorkQueue* queue, LinearAllocator* allocator,
                       const_string lightmapDirPath);

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, uint32* pixels, LinearAllocator* allocator, AppWorkQueue* queue);
