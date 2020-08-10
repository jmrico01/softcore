#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

const uint32 BOX_MAX_TRIANGLES = 32;

struct ComputeMaterial
{
    Vec3 albedo;
    float32 smoothness;
    Vec3 emissionColor;
    float32 emission;
};

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

struct ComputeUbo
{
    const static uint32 MAX_MATERIALS = 32;

	alignas(16) Vec3 cameraPos;
    uint32 numTopLevelBvhs;
	alignas(16) Vec3 filmTopLeft;
	alignas(16) Vec3 filmUnitOffsetX;
	alignas(16) Vec3 filmUnitOffsetY;
	alignas(16) ComputeMaterial materials[MAX_MATERIALS];
};

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

struct VulkanRaytracePipeline
{
    const static uint32 BATCH_SIZE = 16;

    VkQueue computeQueue;
    VkCommandPool computeCommandPool;

    VulkanImage computeImage;

    uint32 numTriangles;
    VulkanBuffer computeTriangles;
    uint32 numBvhs;
    VulkanBuffer computeBvhs;
    VulkanBuffer computeUniform;

    VkDescriptorSetLayout computeDescriptorSetLayout;
    VkDescriptorPool computeDescriptorPool;
    VkDescriptorSet computeDescriptorSet;

    VkPipelineLayout computePipelineLayout;
    VkPipeline computePipeline;

    VkCommandBuffer computeCommandBuffer;
    VkFence computeFence;
};

bool LoadRaytracePipeline(const VulkanWindow& window, VkCommandPool commandPool, uint32 width, uint32 height,
                          const RaycastGeometry& geometry, LinearAllocator* allocator, VulkanRaytracePipeline* pipeline);
void UnloadRaytracePipeline(VkDevice device, VulkanRaytracePipeline* pipeline);

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, float32 fov, const RaycastGeometry& geometry,
                    const uint8* materialIndices, uint32 width, uint32 height, CanvasState* canvas, uint32* pixels,
                    LinearAllocator* allocator, AppWorkQueue* queue);
