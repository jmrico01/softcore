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
    VulkanBuffer vertexBuffer;
    VkSampler sampler;

    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    VkDescriptorSet descriptorSet;

    VkPipelineLayout pipelineLayout;
    VkPipeline pipeline;
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

bool LoadCompositePipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                                 VulkanCompositePipeline* compositePipeline);
void UnloadCompositePipelineWindow(VkDevice device, VulkanCompositePipeline* compositePipeline);
