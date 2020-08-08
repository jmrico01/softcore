#pragma once

#include <km_common/km_array.h>
#include <km_common/km_math.h>
#include <km_common/vulkan/km_vulkan_core.h>
#include <km_common/vulkan/km_vulkan_util.h>

struct VulkanMeshPipeline
{
    uint32 numVertices;
    VulkanBuffer vertexBuffer;
    VulkanBuffer uniformBuffer;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    VkDescriptorSet descriptorSet;

    VkCommandBuffer commandBuffer;

    VkRenderPass renderPass;

    VulkanImage colorImage;
    VulkanImage materialIndexImage;
    VkImage materialIndexImageDst;
    VkDeviceMemory materialIndexImageDstMemory;

    VulkanImage depthImage;

    VkFramebuffer framebuffer;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;

    VkFence fence;
};

struct MeshUniformBufferObject
{
    alignas(16) Mat4 view;
    alignas(16) Mat4 proj;
};

struct VulkanCompositePipeline
{
    const static uint32 TRIANGLE_GEOMETRY_ATLAS_SIZE = 2048;
    const static uint32 TRIANGLE_MATERIALS_ATLAS_SIZE = TRIANGLE_GEOMETRY_ATLAS_SIZE / 4;

    VulkanBuffer vertexBuffer;
    VulkanBuffer uniformBuffer;

    VulkanImage raytracedImage;
    VkSampler sampler;

    VulkanImage triangleGeometry;
    VulkanImage triangleMaterialInds;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    VkDescriptorSet descriptorSet;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;
};

struct CompositeMaterial
{
    alignas(16) Vec3 albedo;
    alignas(16) Vec3 emissionColor;
    float32 smoothness;
    float32 emission;
};

struct CompositeUniformBufferObject
{
    const static uint32 MAX_MATERIALS = 8;

    alignas(16) CompositeMaterial materials[MAX_MATERIALS];
    alignas(16) Vec3 filmTopLeft;
    alignas(16) Vec3 filmUnitOffsetX;
    alignas(16) Vec3 filmUnitOffsetY;
    alignas(16) Vec3 cameraPos;
    uint32 numTriangles;
};

bool LoadMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, VkCommandPool commandPool,
                               LinearAllocator* allocator, VulkanMeshPipeline* meshPipeline);
void UnloadMeshPipelineSwapchain(VkDevice device, VulkanMeshPipeline* meshPipeline);

bool LoadMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, const LoadObjResult& obj,
                            LinearAllocator* allocator, VulkanMeshPipeline* meshPipeline);
void UnloadMeshPipelineWindow(VkDevice device, VulkanMeshPipeline* meshPipeline);

bool LoadCompositePipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain,
                                    VkImageView imageView, LinearAllocator* allocator,
                                    VulkanCompositePipeline* compositePipeline);
void UnloadCompositePipelineSwapchain(VkDevice device, VulkanCompositePipeline* compositePipeline);

bool LoadCompositePipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, const LoadObjResult& obj,
                                 LinearAllocator* allocator, VulkanCompositePipeline* compositePipeline);
void UnloadCompositePipelineWindow(VkDevice device, VulkanCompositePipeline* compositePipeline);
