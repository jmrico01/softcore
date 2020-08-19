#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

const uint32 BOX_MAX_TRIANGLES = 32;

struct RaycastMaterial
{
    Vec3 albedo;
    float32 smoothness;
    Vec3 emissionColor;
    float32 emission;
};

struct RaycastTriangle
{
    Vec3 pos[3];
    Vec3 normal;
    uint32 materialIndex;
};

struct RaycastBvh
{
    Box aabb;
    uint32 startTriangle;
    uint32 endTriangle;
    uint32 skip;
};

struct RaycastMesh
{
    Quat inverseQuat;
    Vec3 offset;
    uint32 startBvh;
    uint32 endBvh;
};

struct RaycastGeometry
{
    Array<string> materialNames;
    Array<RaycastMaterial> materials;
    Array<RaycastTriangle> triangles;
    Array<RaycastBvh> bvhs;
    Array<RaycastMesh> meshes;
};

struct ComputeMaterial
{
    Vec3 albedo;
    float32 smoothness;
    Vec3 emissionColor;
    float32 emission;
};

struct ComputeMesh
{
    Vec4 inverseQuat;
	Vec3 offset;
    uint32 startBvh;
    uint32 endBvh;
    uint32 pad[3];
};
static_assert(sizeof(ComputeMesh) % 16 == 0); // std140 layout

struct ComputeUbo
{
    const static uint32 MAX_MATERIALS = 16;
    const static uint32 MAX_MESHES = 128;

	alignas(16) Vec3 cameraPos;
    uint32 seed;
	alignas(16) Vec3 filmTopLeft;
    uint32 numMeshes;
	alignas(16) Vec3 filmUnitOffsetX;
    float32 frac;
	alignas(16) Vec3 filmUnitOffsetY;
	alignas(16) ComputeMaterial materials[MAX_MATERIALS];
    alignas(16) ComputeMesh meshes[MAX_MESHES];
};

bool CreateRaycastGeometry(const LoadObjResult& obj, uint32 boxMaxTriangles,
                           RaycastGeometry* geometry, LinearAllocator* allocator, LinearAllocator* tempAllocator);

struct CanvasState
{
    static const uint32 MAX_WIDTH = 3840;
    static const uint32 MAX_HEIGHT = 2160;
    static const uint32 MAX_PIXELS = MAX_WIDTH * MAX_HEIGHT;

    float32 test1;
    float32 test2;
    float32 test3;

    StaticArray<Vec3, MAX_PIXELS> colorHdr;
    StaticArray<uint8, MAX_PIXELS> decay;
};

struct VulkanRaytracePipeline
{
    const static uint32 BATCH_SIZE = 16;

    VkQueue queue;
    VkCommandPool commandPool;

    VkImage imageCpu;
    VkDeviceMemory imageCpuMemory;
    VulkanImage imageGpu;

    uint32 numTriangles;
    VulkanBuffer triangles;
    uint32 numBvhs;
    VulkanBuffer bvhs;
    FixedArray<ComputeMesh, ComputeUbo::MAX_MESHES> meshes;
    VulkanBuffer uniform;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    VkDescriptorSet descriptorSet;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;

    VkCommandBuffer commandBuffer;
    VkFence fence;
};

bool LoadRaytracePipeline(const VulkanWindow& window, VkCommandPool commandPool, uint32 width, uint32 height,
                          const RaycastGeometry& geometry, LinearAllocator* allocator, VulkanRaytracePipeline* pipeline);
void UnloadRaytracePipeline(VkDevice device, VulkanRaytracePipeline* pipeline);

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, float32 frac, CanvasState* canvas, uint32* pixels,
                    LinearAllocator* allocator, AppWorkQueue* queue);
