#include "mesh.h"

const VkFormat MATERIAL_INDEX_FORMAT = VK_FORMAT_R8_UINT;

#pragma pack(push, 1)
struct VulkanMeshVertex
{
    Vec3 pos;
    Vec3 normal;
    uint8 materialIndex;
};
#pragma pack(pop)

using VulkanMeshTriangle = StaticArray<VulkanMeshVertex, 3>;

struct VulkanMeshGeometry
{
    Array<uint32> meshEndInds;
    Array<VulkanMeshTriangle> triangles;
};

internal bool ObjToVulkanMeshGeometry(const LoadObjResult& obj, LinearAllocator* allocator, VulkanMeshGeometry* geometry)
{
    geometry->meshEndInds = allocator->NewArray<uint32>(obj.models.size);
    if (geometry->meshEndInds.data == nullptr) {
        return false;
    }

    uint32 totalTriangles = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        totalTriangles += obj.models[i].triangles.size + obj.models[i].quads.size * 2;
    }
    geometry->triangles = allocator->NewArray<VulkanMeshTriangle>(totalTriangles);
    if (geometry->triangles.data == nullptr) {
        return false;
    }

    const Vec3 color = Vec3::one;

    uint32 endInd = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const uint32 tInd = endInd + j;
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry->triangles[tInd][k].pos = t.v[k].pos;
                geometry->triangles[tInd][k].normal = normal;
                geometry->triangles[tInd][k].materialIndex = t.materialIndex;
            }
        }
        endInd += obj.models[i].triangles.size;

        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const ObjQuad& q = obj.models[i].quads[j];
            const uint32 tInd = endInd + j * 2;
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry->triangles[tInd][k].pos = q.v[k].pos;
                geometry->triangles[tInd][k].normal = normal;
                geometry->triangles[tInd][k].materialIndex = q.materialIndex;
            }

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                geometry->triangles[tInd + 1][k].pos = q.v[quadInd].pos;
                geometry->triangles[tInd + 1][k].normal = normal;
                geometry->triangles[tInd + 1][k].materialIndex = q.materialIndex;
            }
        }
        endInd += obj.models[i].quads.size * 2;

        geometry->meshEndInds[i] = endInd;
    }

    return true;
}

bool LoadMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, VkCommandPool commandPool,
                               LinearAllocator* allocator, VulkanMeshPipeline* meshPipeline)
{
    VkFormatProperties formatProperties;
    vkGetPhysicalDeviceFormatProperties(window.physicalDevice, MATERIAL_INDEX_FORMAT, &formatProperties);
    if ((formatProperties.linearTilingFeatures & VK_FORMAT_FEATURE_TRANSFER_SRC_BIT) == 0) {
        LOG_ERROR("Material index format not eligible as transfer src for linear tiling\n");
        return false;
    }
    if ((formatProperties.optimalTilingFeatures & VK_FORMAT_FEATURE_TRANSFER_DST_BIT) == 0) {
        LOG_ERROR("Material index format not eligible as transfer dst for optimal tiling\n");
        return false;
    }
    if ((formatProperties.optimalTilingFeatures & VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BIT) == 0) {
        LOG_ERROR("Material index format not eligible as color attachment for optimal tiling\n");
        return false;
    }
    if ((formatProperties.bufferFeatures & VK_FORMAT_FEATURE_VERTEX_BUFFER_BIT) == 0) {
        LOG_ERROR("Material index format not eligible as vertex buffer input\n");
        return false;
    }

    // Create render pass
    {
        VkAttachmentDescription attachments[2] = {};

        // TODO possibly save a GPU-local color buffer for later compositing with ray traced image

        // Attachment 0 - material index
        attachments[0].format = MATERIAL_INDEX_FORMAT;
        attachments[0].samples = VK_SAMPLE_COUNT_1_BIT;
        attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        attachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
        attachments[0].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        attachments[0].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        attachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        attachments[0].finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;

        // Attachment 1 - depth
        attachments[1].format = VK_FORMAT_D32_SFLOAT;
        attachments[1].samples = VK_SAMPLE_COUNT_1_BIT;
        attachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
        attachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        attachments[1].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
        attachments[1].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
        attachments[1].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        attachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

        VkAttachmentReference colorAttachmentRefs[1] = {};
        colorAttachmentRefs[0].attachment = 0;
        colorAttachmentRefs[0].layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

        VkAttachmentReference depthAttachmentRef = {};
        depthAttachmentRef.attachment = 1;
        depthAttachmentRef.layout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

        VkSubpassDescription subpass = {};
        subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
        subpass.colorAttachmentCount = C_ARRAY_LENGTH(colorAttachmentRefs);
        subpass.pColorAttachments = colorAttachmentRefs;
        subpass.pDepthStencilAttachment = &depthAttachmentRef;

        VkSubpassDependency dependency = {};
        dependency.srcSubpass = VK_SUBPASS_EXTERNAL;
        dependency.dstSubpass = 0;
        dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.srcAccessMask = 0;
        dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

        VkRenderPassCreateInfo renderPassCreateInfo = {};
        renderPassCreateInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
        renderPassCreateInfo.attachmentCount = C_ARRAY_LENGTH(attachments);
        renderPassCreateInfo.pAttachments = attachments;
        renderPassCreateInfo.subpassCount = 1;
        renderPassCreateInfo.pSubpasses = &subpass;
        renderPassCreateInfo.dependencyCount = 1;
        renderPassCreateInfo.pDependencies = &dependency;

        if (vkCreateRenderPass(window.device, &renderPassCreateInfo, nullptr, &meshPipeline->renderPass) != VK_SUCCESS) {
            LOG_ERROR("vkCreateRenderPass failed\n");
            return false;
        }
    }

    // Create color attachment and dst image
    {
        if (!CreateImage(window.device, window.physicalDevice,
                         swapchain.extent.width, swapchain.extent.height,
                         MATERIAL_INDEX_FORMAT, VK_IMAGE_TILING_OPTIMAL,
                         VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
                         VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                         &meshPipeline->materialIndexImage.image, &meshPipeline->materialIndexImage.memory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        if (!CreateImageView(window.device, meshPipeline->materialIndexImage.image, MATERIAL_INDEX_FORMAT,
                             VK_IMAGE_ASPECT_COLOR_BIT, &meshPipeline->materialIndexImage.view)) {
            LOG_ERROR("CreateImageView failed\n");
            return false;
        }

        if (!CreateImage(window.device, window.physicalDevice,
                         swapchain.extent.width, swapchain.extent.height,
                         MATERIAL_INDEX_FORMAT, VK_IMAGE_TILING_LINEAR,
                         VK_IMAGE_USAGE_TRANSFER_DST_BIT,
                         VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                         &meshPipeline->materialIndexImageDst, &meshPipeline->materialIndexImageDstMemory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, commandPool, window.graphicsQueue);
        TransitionImageLayout(commandBuffer, meshPipeline->materialIndexImageDst,
                              VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL);
    }

    // Create depth attachment
    {
        if (!CreateImage(window.device, window.physicalDevice,
                         swapchain.extent.width, swapchain.extent.height,
                         VK_FORMAT_D32_SFLOAT, VK_IMAGE_TILING_OPTIMAL,
                         VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
                         VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                         &meshPipeline->depthImage.image, &meshPipeline->depthImage.memory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        if (!CreateImageView(window.device, meshPipeline->depthImage.image, VK_FORMAT_D32_SFLOAT,
                             VK_IMAGE_ASPECT_DEPTH_BIT, &meshPipeline->depthImage.view)) {
            LOG_ERROR("CreateImageView failed\n");
            return false;
        }
    }

    // Create framebuffers
    {
        const VkImageView attachments[] = {
            meshPipeline->materialIndexImage.view,
            meshPipeline->depthImage.view,
        };

        VkFramebufferCreateInfo framebufferCreateInfo = {};
        framebufferCreateInfo.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
        framebufferCreateInfo.renderPass = meshPipeline->renderPass;
        framebufferCreateInfo.attachmentCount = C_ARRAY_LENGTH(attachments);
        framebufferCreateInfo.pAttachments = attachments;
        framebufferCreateInfo.width = swapchain.extent.width;
        framebufferCreateInfo.height = swapchain.extent.height;
        framebufferCreateInfo.layers = 1;

        if (vkCreateFramebuffer(window.device, &framebufferCreateInfo, nullptr,
                                &meshPipeline->framebuffer) != VK_SUCCESS) {
            LOG_ERROR("vkCreateFramebuffer failed\n");
            return false;
        }
    }

    // Create pipeline
    {
        const Array<uint8> vertShaderCode = LoadEntireFile(ToString("data/shaders/mesh.vert.spv"), allocator);
        if (vertShaderCode.data == nullptr) {
            LOG_ERROR("Failed to load vertex shader code\n");
            return false;
        }
        const Array<uint8> fragShaderCode = LoadEntireFile(ToString("data/shaders/mesh.frag.spv"), allocator);
        if (fragShaderCode.data == nullptr) {
            LOG_ERROR("Failed to load fragment shader code\n");
            return false;
        }

        VkShaderModule vertShaderModule;
        if (!CreateShaderModule(vertShaderCode, window.device, &vertShaderModule)) {
            LOG_ERROR("Failed to create vertex shader module\n");
            return false;
        }
        defer(vkDestroyShaderModule(window.device, vertShaderModule, nullptr));

        VkShaderModule fragShaderModule;
        if (!CreateShaderModule(fragShaderCode, window.device, &fragShaderModule)) {
            LOG_ERROR("Failed to create fragment shader module\n");
            return false;
        }
        defer(vkDestroyShaderModule(window.device, fragShaderModule, nullptr));

        VkPipelineShaderStageCreateInfo vertShaderStageCreateInfo = {};
        vertShaderStageCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        vertShaderStageCreateInfo.stage = VK_SHADER_STAGE_VERTEX_BIT;
        vertShaderStageCreateInfo.module = vertShaderModule;
        vertShaderStageCreateInfo.pName = "main";
        // vertShaderStageCreateInfo.pSpecializationInfo is useful for setting shader constants

        VkPipelineShaderStageCreateInfo fragShaderStageCreateInfo = {};
        fragShaderStageCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        fragShaderStageCreateInfo.stage = VK_SHADER_STAGE_FRAGMENT_BIT;
        fragShaderStageCreateInfo.module = fragShaderModule;
        fragShaderStageCreateInfo.pName = "main";

        VkPipelineShaderStageCreateInfo shaderStages[] = { vertShaderStageCreateInfo, fragShaderStageCreateInfo };

        VkVertexInputBindingDescription bindingDescription = {};
        bindingDescription.binding = 0;
        bindingDescription.stride = sizeof(VulkanMeshVertex);
        bindingDescription.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

        VkVertexInputAttributeDescription attributeDescriptions[3] = {};

        attributeDescriptions[0].binding = 0;
        attributeDescriptions[0].location = 0;
        attributeDescriptions[0].format = VK_FORMAT_R32G32B32_SFLOAT;
        attributeDescriptions[0].offset = offsetof(VulkanMeshVertex, pos);

        attributeDescriptions[1].binding = 0;
        attributeDescriptions[1].location = 1;
        attributeDescriptions[1].format = VK_FORMAT_R32G32B32_SFLOAT;
        attributeDescriptions[1].offset = offsetof(VulkanMeshVertex, normal);

        attributeDescriptions[2].binding = 0;
        attributeDescriptions[2].location = 2;
        attributeDescriptions[2].format = MATERIAL_INDEX_FORMAT;
        attributeDescriptions[2].offset = offsetof(VulkanMeshVertex, materialIndex);

        VkPipelineVertexInputStateCreateInfo vertexInputCreateInfo = {};
        vertexInputCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
        vertexInputCreateInfo.vertexBindingDescriptionCount = 1;
        vertexInputCreateInfo.pVertexBindingDescriptions = &bindingDescription;
        vertexInputCreateInfo.vertexAttributeDescriptionCount = C_ARRAY_LENGTH(attributeDescriptions);
        vertexInputCreateInfo.pVertexAttributeDescriptions = attributeDescriptions;

        VkPipelineInputAssemblyStateCreateInfo inputAssemblyCreateInfo = {};
        inputAssemblyCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO;
        inputAssemblyCreateInfo.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
        inputAssemblyCreateInfo.primitiveRestartEnable = VK_FALSE;

        VkViewport viewport = {};
        viewport.x = 0.0f;
        viewport.y = 0.0f;
        viewport.width = (float32)swapchain.extent.width;
        viewport.height = (float32)swapchain.extent.height;
        viewport.minDepth = 0.0f;
        viewport.maxDepth = 1.0f;

        VkRect2D scissor = {};
        scissor.offset = { 0, 0 };
        scissor.extent = swapchain.extent;

        VkPipelineViewportStateCreateInfo viewportStateCreateInfo = {};
        viewportStateCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO;
        viewportStateCreateInfo.viewportCount = 1;
        viewportStateCreateInfo.pViewports = &viewport;
        viewportStateCreateInfo.scissorCount = 1;
        viewportStateCreateInfo.pScissors = &scissor;

        VkPipelineRasterizationStateCreateInfo rasterizerCreateInfo = {};
        rasterizerCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO;
        rasterizerCreateInfo.depthClampEnable = VK_FALSE;
        rasterizerCreateInfo.rasterizerDiscardEnable = VK_FALSE;
        rasterizerCreateInfo.polygonMode = VK_POLYGON_MODE_FILL;
        rasterizerCreateInfo.lineWidth = 1.0f;
        rasterizerCreateInfo.cullMode = VK_CULL_MODE_BACK_BIT;
        rasterizerCreateInfo.frontFace = VK_FRONT_FACE_CLOCKWISE;
        rasterizerCreateInfo.depthBiasEnable = VK_FALSE;
        rasterizerCreateInfo.depthBiasConstantFactor = 0.0f;
        rasterizerCreateInfo.depthBiasClamp = 0.0f;
        rasterizerCreateInfo.depthBiasSlopeFactor = 0.0f;

        VkPipelineMultisampleStateCreateInfo multisampleCreateInfo = {};
        multisampleCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO;
        multisampleCreateInfo.sampleShadingEnable = VK_FALSE;
        multisampleCreateInfo.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;
        multisampleCreateInfo.minSampleShading = 1.0f;
        multisampleCreateInfo.pSampleMask = nullptr;
        multisampleCreateInfo.alphaToCoverageEnable = VK_FALSE;
        multisampleCreateInfo.alphaToOneEnable = VK_FALSE;

        VkPipelineColorBlendAttachmentState colorBlendAttachments[1] = {};
        colorBlendAttachments[0].blendEnable = VK_FALSE;

        VkPipelineColorBlendStateCreateInfo colorBlendingCreateInfo = {};
        colorBlendingCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
        colorBlendingCreateInfo.logicOpEnable = VK_FALSE;
        colorBlendingCreateInfo.logicOp = VK_LOGIC_OP_COPY;
        colorBlendingCreateInfo.attachmentCount = C_ARRAY_LENGTH(colorBlendAttachments);
        colorBlendingCreateInfo.pAttachments = colorBlendAttachments;
        colorBlendingCreateInfo.blendConstants[0] = 0.0f;
        colorBlendingCreateInfo.blendConstants[1] = 0.0f;
        colorBlendingCreateInfo.blendConstants[2] = 0.0f;
        colorBlendingCreateInfo.blendConstants[3] = 0.0f;

        VkPipelineDepthStencilStateCreateInfo depthStencilCreateInfo = {};
        depthStencilCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO;
        depthStencilCreateInfo.depthTestEnable = VK_TRUE;
        depthStencilCreateInfo.depthWriteEnable = VK_TRUE;
        depthStencilCreateInfo.depthCompareOp = VK_COMPARE_OP_LESS;
        depthStencilCreateInfo.depthBoundsTestEnable = VK_FALSE;
        depthStencilCreateInfo.minDepthBounds = 0.0f; // disabled
        depthStencilCreateInfo.maxDepthBounds = 1.0f; // disabled
        depthStencilCreateInfo.stencilTestEnable = VK_FALSE;

        VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo = {};
        pipelineLayoutCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
        pipelineLayoutCreateInfo.setLayoutCount = 1;
        pipelineLayoutCreateInfo.pSetLayouts = &meshPipeline->descriptorSetLayout;
        pipelineLayoutCreateInfo.pushConstantRangeCount = 0;
        pipelineLayoutCreateInfo.pPushConstantRanges = nullptr;

        if (vkCreatePipelineLayout(window.device, &pipelineLayoutCreateInfo, nullptr,
                                   &meshPipeline->pipelineLayout) != VK_SUCCESS) {
            LOG_ERROR("vkCreatePipelineLayout failed\n");
            return false;
        }

        VkGraphicsPipelineCreateInfo pipelineCreateInfo = {};
        pipelineCreateInfo.sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO;
        pipelineCreateInfo.stageCount = C_ARRAY_LENGTH(shaderStages);
        pipelineCreateInfo.pStages = shaderStages;
        pipelineCreateInfo.pVertexInputState = &vertexInputCreateInfo;
        pipelineCreateInfo.pInputAssemblyState = &inputAssemblyCreateInfo;
        pipelineCreateInfo.pViewportState = &viewportStateCreateInfo;
        pipelineCreateInfo.pRasterizationState = &rasterizerCreateInfo;
        pipelineCreateInfo.pMultisampleState = &multisampleCreateInfo;
        pipelineCreateInfo.pDepthStencilState = &depthStencilCreateInfo;
        pipelineCreateInfo.pColorBlendState = &colorBlendingCreateInfo;
        pipelineCreateInfo.pDynamicState = nullptr;
        pipelineCreateInfo.layout = meshPipeline->pipelineLayout;
        pipelineCreateInfo.renderPass = meshPipeline->renderPass;
        pipelineCreateInfo.subpass = 0;
        pipelineCreateInfo.basePipelineHandle = VK_NULL_HANDLE;
        pipelineCreateInfo.basePipelineIndex = -1;

        if (vkCreateGraphicsPipelines(window.device, VK_NULL_HANDLE, 1, &pipelineCreateInfo, nullptr,
                                      &meshPipeline->pipeline) != VK_SUCCESS) {
            LOG_ERROR("vkCreateGraphicsPipeline failed\n");
            return false;
        }
    }

    // Record command buffer
    {
        VkCommandBufferBeginInfo beginInfo = {};
        beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        beginInfo.flags = 0;
        beginInfo.pInheritanceInfo = nullptr;

        if (vkBeginCommandBuffer(meshPipeline->commandBuffer, &beginInfo) != VK_SUCCESS) {
            LOG_ERROR("vkBeginCommandBuffer failed\n");
            return false;
        }

        const VkClearValue clearValues[] = {
            { 0, 0, 0, 0 },
            { 1.0f, 0 }
        };

        VkRenderPassBeginInfo renderPassInfo = {};
        renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
        renderPassInfo.renderPass = meshPipeline->renderPass;
        renderPassInfo.framebuffer = meshPipeline->framebuffer;
        renderPassInfo.renderArea.offset = { 0, 0 };
        renderPassInfo.renderArea.extent = swapchain.extent;
        renderPassInfo.clearValueCount = C_ARRAY_LENGTH(clearValues);
        renderPassInfo.pClearValues = clearValues;

        vkCmdBeginRenderPass(meshPipeline->commandBuffer, &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

        vkCmdBindPipeline(meshPipeline->commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, meshPipeline->pipeline);

        const VkDeviceSize offsets[] = { 0 };
        vkCmdBindVertexBuffers(meshPipeline->commandBuffer, 0, 1, &meshPipeline->vertexBuffer.buffer, offsets);

        vkCmdBindDescriptorSets(meshPipeline->commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS,
                                meshPipeline->pipelineLayout, 0, 1,
                                &meshPipeline->descriptorSet, 0, nullptr);

        vkCmdDraw(meshPipeline->commandBuffer, meshPipeline->numVertices, 1, 0, 0);

        vkCmdEndRenderPass(meshPipeline->commandBuffer);

        VkImageSubresourceLayers imageSubresourceLayers;
        imageSubresourceLayers.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageSubresourceLayers.mipLevel = 0;
        imageSubresourceLayers.baseArrayLayer = 0;
        imageSubresourceLayers.layerCount = 1;

        VkImageCopy imageCopy;
        imageCopy.srcSubresource = imageSubresourceLayers;
        imageCopy.srcOffset.x = 0;
        imageCopy.srcOffset.y = 0;
        imageCopy.srcOffset.z = 0;
        imageCopy.dstSubresource = imageSubresourceLayers;
        imageCopy.dstOffset.x = 0;
        imageCopy.dstOffset.y = 0;
        imageCopy.dstOffset.z = 0;
        imageCopy.extent.width = swapchain.extent.width;
        imageCopy.extent.height = swapchain.extent.height;
        imageCopy.extent.depth = 1;

        vkCmdCopyImage(meshPipeline->commandBuffer,
                       meshPipeline->materialIndexImage.image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                       meshPipeline->materialIndexImageDst, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                       1, &imageCopy);

        if (vkEndCommandBuffer(meshPipeline->commandBuffer) != VK_SUCCESS) {
            LOG_ERROR("vkEndCommandBuffer failed\n");
        }
    }

    return true;
}

void UnloadMeshPipelineSwapchain(VkDevice device, VulkanMeshPipeline* meshPipeline)
{
    vkDestroyPipeline(device, meshPipeline->pipeline, nullptr);
    vkDestroyPipelineLayout(device, meshPipeline->pipelineLayout, nullptr);

    vkDestroyFramebuffer(device, meshPipeline->framebuffer, nullptr);

    DestroyVulkanImage(device, &meshPipeline->depthImage);

    vkDestroyImage(device, meshPipeline->materialIndexImageDst, nullptr);
    vkFreeMemory(device, meshPipeline->materialIndexImageDstMemory, nullptr);
    DestroyVulkanImage(device, &meshPipeline->materialIndexImage);

    vkDestroyRenderPass(device, meshPipeline->renderPass, nullptr);
}

bool LoadMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                            VulkanMeshPipeline* meshPipeline)
{
    // Create vertex buffers
    {
        SCOPED_ALLOCATOR_RESET(*allocator);

        const_string objPath = ToString("data/models/interior-1.obj");

        LoadObjResult obj;
        if (!LoadObj(objPath, Vec3::zero, 1.0f, &obj, allocator)) {
            LOG_ERROR("Failed to load obj file %.*s\n", objPath.size, objPath.data);
            return false;
        }

        VulkanMeshGeometry geometry;
        if (!ObjToVulkanMeshGeometry(obj, allocator, &geometry)) {
            LOG_ERROR("Failed to load Vulkan geometry from obj\n");
            return false;
        }

        const Array<VulkanMeshVertex> vertices = {
            .size = geometry.triangles.size * 3,
            .data = &geometry.triangles[0][0]
        };
        meshPipeline->numVertices = vertices.size;

        const VkDeviceSize vertexBufferSize = vertices.size * sizeof(VulkanMeshVertex);

        VulkanBuffer stagingBuffer;
        if (!CreateVulkanBuffer(vertexBufferSize,
                                VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &stagingBuffer)) {
            LOG_ERROR("CreateBuffer failed for staging buffer\n");
            return false;
        }

        // Copy vertex data from CPU into memory-mapped staging buffer
        void* data;
        vkMapMemory(window.device, stagingBuffer.memory, 0, vertexBufferSize, 0, &data);
        MemCopy(data, vertices.data, vertexBufferSize);
        vkUnmapMemory(window.device, stagingBuffer.memory);

        if (!CreateVulkanBuffer(vertexBufferSize,
                                VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                                window.device, window.physicalDevice, &meshPipeline->vertexBuffer)) {
            LOG_ERROR("CreateBuffer failed for vertex buffer\n");
            return false;
        }

        // Copy vertex data from staging buffer into GPU vertex buffer
        {
            SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, commandPool, window.graphicsQueue);
            CopyBuffer(commandBuffer, stagingBuffer.buffer, meshPipeline->vertexBuffer.buffer, vertexBufferSize);
        }

        DestroyVulkanBuffer(window.device, &stagingBuffer);
    }

    // Create uniform buffer
    {
        VkDeviceSize uniformBufferSize = sizeof(MeshUniformBufferObject);
        if (!CreateVulkanBuffer(uniformBufferSize,
                                VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &meshPipeline->uniformBuffer)) {
            LOG_ERROR("CreateBuffer failed for uniform buffer\n");
            return false;
        }
    }

    // Create descriptor set layout
    {
        VkDescriptorSetLayoutBinding uboLayoutBinding = {};
        uboLayoutBinding.binding = 0;
        uboLayoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        uboLayoutBinding.descriptorCount = 1;
        uboLayoutBinding.stageFlags = VK_SHADER_STAGE_VERTEX_BIT;
        uboLayoutBinding.pImmutableSamplers = nullptr;

        VkDescriptorSetLayoutCreateInfo layoutCreateInfo = {};
        layoutCreateInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layoutCreateInfo.bindingCount = 1;
        layoutCreateInfo.pBindings = &uboLayoutBinding;

        if (vkCreateDescriptorSetLayout(window.device, &layoutCreateInfo, nullptr,
                                        &meshPipeline->descriptorSetLayout) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorSetLayout failed\n");
            return false;
        }
    }

    // Create descriptor pool
    {
        VkDescriptorPoolSize poolSize = {};
        poolSize.type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        poolSize.descriptorCount = 1;

        VkDescriptorPoolCreateInfo poolInfo = {};
        poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        poolInfo.poolSizeCount = 1;
        poolInfo.pPoolSizes = &poolSize;
        poolInfo.maxSets = 1;

        if (vkCreateDescriptorPool(window.device, &poolInfo, nullptr, &meshPipeline->descriptorPool) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorPool failed\n");
            return false;
        }
    }

    // Create descriptor set
    {
        VkDescriptorSetAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        allocInfo.descriptorPool = meshPipeline->descriptorPool;
        allocInfo.descriptorSetCount = 1;
        allocInfo.pSetLayouts = &meshPipeline->descriptorSetLayout;

        if (vkAllocateDescriptorSets(window.device, &allocInfo, &meshPipeline->descriptorSet) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateDescriptorSets failed\n");
            return false;
        }

        VkWriteDescriptorSet descriptorWrite = {};

        VkDescriptorBufferInfo bufferInfo = {};
        bufferInfo.buffer = meshPipeline->uniformBuffer.buffer;
        bufferInfo.offset = 0;
        bufferInfo.range = sizeof(MeshUniformBufferObject);

        descriptorWrite.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        descriptorWrite.dstSet = meshPipeline->descriptorSet;
        descriptorWrite.dstBinding = 0;
        descriptorWrite.dstArrayElement = 0;
        descriptorWrite.descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        descriptorWrite.descriptorCount = 1;
        descriptorWrite.pBufferInfo = &bufferInfo;

        vkUpdateDescriptorSets(window.device, 1, &descriptorWrite, 0, nullptr);
    }

    // Allocate command buffer
    {
        VkCommandBufferAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        allocInfo.commandPool = commandPool;
        allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        allocInfo.commandBufferCount = 1;

        if (vkAllocateCommandBuffers(window.device, &allocInfo, &meshPipeline->commandBuffer) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateCommandBuffers failed\n");
            return false;
        }
    }

    // Create fence
    {
        VkFenceCreateInfo fenceCreateInfo = {};
        fenceCreateInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fenceCreateInfo.flags = 0;

        if (vkCreateFence(window.device, &fenceCreateInfo, nullptr, &meshPipeline->fence) != VK_SUCCESS) {
            LOG_ERROR("vkCreateFence failed\n");
            return false;
        }
    }

    return true;
}

void UnloadMeshPipelineWindow(VkDevice device, VulkanMeshPipeline* meshPipeline)
{
    vkDestroyFence(device, meshPipeline->fence, nullptr);

    vkDestroyDescriptorPool(device, meshPipeline->descriptorPool, nullptr);
    vkDestroyDescriptorSetLayout(device, meshPipeline->descriptorSetLayout, nullptr);

    DestroyVulkanBuffer(device, &meshPipeline->uniformBuffer);
    DestroyVulkanBuffer(device, &meshPipeline->vertexBuffer);
}
