#pragma once

#include <km_common/km_array.h>
#include <km_common/km_math.h>
#include <km_common/vulkan/km_vulkan_core.h>
#include <km_common/vulkan/km_vulkan_util.h>

enum class MeshId
{
    TILE,
    MOB,

    COUNT
};

struct VulkanLightmapMeshPipeline
{
    static const uint32 MAX_MESHES = 64;
    FixedArray<uint32, MAX_MESHES> meshTriangleEndInds;
    VulkanBuffer vertexBuffer;

    static const uint32 MAX_LIGHTMAPS = 64;
    FixedArray<VulkanImage, MAX_LIGHTMAPS> lightmaps;
    VkSampler lightmapSampler;

    VulkanBuffer uniformBuffer;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    FixedArray<VkDescriptorSet, MAX_LIGHTMAPS> descriptorSets;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;
};

struct VulkanMeshPipeline
{
    static const uint32 MAX_MESHES = (uint32)MeshId::COUNT;
    static const uint32 MAX_INSTANCES = 128 * 128 * 128;

    StaticArray<uint32, MAX_MESHES> numVertices;
    VulkanBuffer vertexBuffer;
    VulkanBuffer instanceBuffer;
    VulkanBuffer uniformBuffer;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    VkDescriptorSet descriptorSet;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;
};

struct VulkanMeshInstanceData
{
    Mat4 model;
    Vec3 color;
    Vec3 collapseMid;
    float32 collapseT;
};

struct VulkanMeshRenderState
{
    using MeshInstanceData = FixedArray<VulkanMeshInstanceData, VulkanMeshPipeline::MAX_INSTANCES>;
    StaticArray<MeshInstanceData, VulkanMeshPipeline::MAX_MESHES> meshInstanceData;
};

void PushMesh(MeshId meshId, Mat4 model, Vec3 color, Vec3 collapseMid, float32 collapseT,
              VulkanMeshRenderState* renderState);

void ResetMeshRenderState(VulkanMeshRenderState* renderState);

void UploadAndSubmitMeshDrawCommands(VkDevice device, VkCommandBuffer commandBuffer,
                                     const VulkanMeshPipeline& meshPipeline, const VulkanMeshRenderState& renderState,
                                     Mat4 view, Mat4 proj, LinearAllocator* allocator);

void UploadAndSubmitLightmapMeshDrawCommands(VkDevice device, VkCommandBuffer commandBuffer,
                                             const VulkanLightmapMeshPipeline& lightmapMeshPipeline,
                                             Mat4 model, Mat4 view, Mat4 proj);

bool LoadMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, LinearAllocator* allocator,
                               VulkanMeshPipeline* meshPipeline);
void UnloadMeshPipelineSwapchain(VkDevice device, VulkanMeshPipeline* meshPipeline);

bool LoadMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                            VulkanMeshPipeline* meshPipeline);
void UnloadMeshPipelineWindow(VkDevice device, VulkanMeshPipeline* meshPipeline);

bool LoadLightmapMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, LinearAllocator* allocator, VulkanLightmapMeshPipeline* lightmapMeshPipeline);
void UnloadLightmapMeshPipelineSwapchain(VkDevice device, VulkanLightmapMeshPipeline* lightmapMeshPipeline);

bool LoadLightmapMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                                    VulkanLightmapMeshPipeline* lightmapMeshPipeline);
void UnloadLightmapMeshPipelineWindow(VkDevice device, VulkanLightmapMeshPipeline* lightmapMeshPipeline);