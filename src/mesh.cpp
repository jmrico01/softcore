#include "mesh.h"

struct VulkanMeshVertex
{
    Vec3 pos;
    Vec3 normal;
    Vec3 color;
};

using VulkanMeshTriangle = StaticArray<VulkanMeshVertex, 3>;

struct VulkanMeshGeometry
{
    bool valid;
    Array<uint32> meshEndInds;
    Array<VulkanMeshTriangle> triangles;
};

struct MeshUniformBufferObject
{
    alignas(16) Mat4 view;
    alignas(16) Mat4 proj;
};

struct VulkanLightmapMeshVertex
{
    Vec3 pos;
    Vec3 normal;
    Vec3 color;
    Vec2 uv;
    float32 lightmapWeight;
};

using VulkanLightmapMeshTriangle = StaticArray<VulkanLightmapMeshVertex, 3>;

struct VulkanLightmapMeshGeometry
{
    bool valid;
    Array<uint32> meshEndInds;
    Array<VulkanLightmapMeshTriangle> triangles;
};

struct LightmapMeshUniformBufferObject
{
    alignas(16) Mat4 model;
    alignas(16) Mat4 view;
    alignas(16) Mat4 proj;
};

internal VulkanMeshGeometry ObjToVulkanMeshGeometry(const LoadObjResult& obj, LinearAllocator* allocator)
{
    VulkanMeshGeometry geometry;
    geometry.valid = false;
    geometry.meshEndInds = allocator->NewArray<uint32>(obj.models.size);
    if (geometry.meshEndInds.data == nullptr) {
        return geometry;
    }

    uint32 totalTriangles = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        totalTriangles += obj.models[i].triangles.size + obj.models[i].quads.size * 2;
    }
    geometry.triangles = allocator->NewArray<VulkanMeshTriangle>(totalTriangles);
    if (geometry.triangles.data == nullptr) {
        return geometry;
    }

    const Vec3 color = Vec3::one;

    uint32 endInd = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const uint32 tInd = endInd + j;
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry.triangles[tInd][k].pos = t.v[k].pos;
                geometry.triangles[tInd][k].normal = normal;
                geometry.triangles[tInd][k].color = color;
            }
        }
        endInd += obj.models[i].triangles.size;

        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const ObjQuad& q = obj.models[i].quads[j];
            const uint32 tInd = endInd + j * 2;
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry.triangles[tInd][k].pos = q.v[k].pos;
                geometry.triangles[tInd][k].normal = normal;
                geometry.triangles[tInd][k].color = color;
            }

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                geometry.triangles[tInd + 1][k].pos = q.v[quadInd].pos;
                geometry.triangles[tInd + 1][k].normal = normal;
                geometry.triangles[tInd + 1][k].color = color;
            }
        }
        endInd += obj.models[i].quads.size * 2;

        geometry.meshEndInds[i] = endInd;
    }

    geometry.valid = true;
    return geometry;
}

internal VulkanLightmapMeshGeometry ObjToVulkanLightmapMeshGeometry(const LoadObjResult& obj, LinearAllocator* allocator)
{
    VulkanLightmapMeshGeometry geometry;
    geometry.valid = false;
    geometry.meshEndInds = allocator->NewArray<uint32>(obj.models.size);
    if (geometry.meshEndInds.data == nullptr) {
        return geometry;
    }

    uint32 totalTriangles = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        totalTriangles += obj.models[i].triangles.size + obj.models[i].quads.size * 2;
    }
    geometry.triangles = allocator->NewArray<VulkanLightmapMeshTriangle>(totalTriangles);
    if (geometry.triangles.data == nullptr) {
        return geometry;
    }

    const Vec3 vertexColor = Vec3::zero;

    uint32 endInd = 0;
    for (uint32 i = 0; i < obj.models.size; i++) {
        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const uint32 tInd = endInd + j;
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry.triangles[tInd][k].pos = t.v[k].pos;
                geometry.triangles[tInd][k].normal = normal;
                geometry.triangles[tInd][k].color = vertexColor;
                geometry.triangles[tInd][k].uv = t.v[k].uv;
            }
        }
        endInd += obj.models[i].triangles.size;

        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const ObjQuad& q = obj.models[i].quads[j];
            const uint32 tInd = endInd + j * 2;
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                geometry.triangles[tInd][k].pos = q.v[k].pos;
                geometry.triangles[tInd][k].normal = normal;
                geometry.triangles[tInd][k].color = vertexColor;
                geometry.triangles[tInd][k].uv = q.v[k].uv;
            }

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                geometry.triangles[tInd + 1][k].pos = q.v[quadInd].pos;
                geometry.triangles[tInd + 1][k].normal = normal;
                geometry.triangles[tInd + 1][k].color = vertexColor;
                geometry.triangles[tInd + 1][k].uv = q.v[quadInd].uv;
            }
        }
        endInd += obj.models[i].quads.size * 2;

        geometry.meshEndInds[i] = endInd;
    }

    geometry.valid = true;
    return geometry;
}

void PushMesh(MeshId meshId, Mat4 model, Vec3 color, Vec3 collapseMid, float32 collapseT,
              VulkanMeshRenderState* renderState)
{
    const uint32 meshIndex = (uint32)meshId;
    DEBUG_ASSERT(meshIndex < renderState->meshInstanceData.SIZE);

    VulkanMeshInstanceData* instanceData = renderState->meshInstanceData[meshIndex].Append();
    instanceData->model = model;
    instanceData->color = color;
    instanceData->collapseMid = collapseMid;
    instanceData->collapseT = collapseT;
}

void ResetMeshRenderState(VulkanMeshRenderState* renderState)
{
    for (uint32 i = 0; i < renderState->meshInstanceData.SIZE; i++) {
        renderState->meshInstanceData[i].Clear();
    }
}

void UploadAndSubmitMeshDrawCommands(VkDevice device, VkCommandBuffer commandBuffer,
                                     const VulkanMeshPipeline& meshPipeline, const VulkanMeshRenderState& renderState,
                                     Mat4 view, Mat4 proj, LinearAllocator* allocator)
{
    Array<uint32> meshNumInstances = allocator->NewArray<uint32>(meshPipeline.MAX_MESHES);
    uint32 totalNumInstances = 0;
    for (uint32 i = 0; i < meshPipeline.MAX_MESHES; i++) {
        meshNumInstances[i] = renderState.meshInstanceData[i].size;
        totalNumInstances += meshNumInstances[i];
    }

    if (totalNumInstances > meshPipeline.MAX_INSTANCES) {
        LOG_ERROR("Too many mesh instances: %lu, max %lu\n", totalNumInstances, meshPipeline.MAX_INSTANCES);
        // TODO what to do here?
        DEBUG_PANIC("too many meshes!\n");
    }

    Array<VulkanMeshInstanceData> instanceData = allocator->NewArray<VulkanMeshInstanceData>(totalNumInstances);
    uint32 instanceOffset = 0;
    for (uint32 i = 0; i < meshPipeline.MAX_MESHES; i++) {
        const uint32 instances = meshNumInstances[i];
        MemCopy(instanceData.data + instanceOffset,
                renderState.meshInstanceData[i].data, instances * sizeof(VulkanMeshInstanceData));
        instanceOffset += instances;
    }

    if (instanceData.size > 0) {
        void* data;
        const uint32 bufferSize = instanceData.size * sizeof(VulkanMeshInstanceData);
        vkMapMemory(device, meshPipeline.instanceBuffer.memory, 0, bufferSize, 0, &data);
        MemCopy(data, instanceData.data, bufferSize);
        vkUnmapMemory(device, meshPipeline.instanceBuffer.memory);
    }

    const MeshUniformBufferObject ubo = { .view = view, .proj = proj };

    void* data;
    vkMapMemory(device, meshPipeline.uniformBuffer.memory, 0, sizeof(MeshUniformBufferObject), 0, &data);
    MemCopy(data, &ubo, sizeof(ubo));
    vkUnmapMemory(device, meshPipeline.uniformBuffer.memory);

    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, meshPipeline.pipeline);

    const VkBuffer vertexBuffers[] = {
        meshPipeline.vertexBuffer.buffer,
        meshPipeline.instanceBuffer.buffer
    };
    const VkDeviceSize offsets[] = { 0, 0 };
    vkCmdBindVertexBuffers(commandBuffer, 0, C_ARRAY_LENGTH(vertexBuffers), vertexBuffers, offsets);

    vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, meshPipeline.pipelineLayout, 0, 1,
                            &meshPipeline.descriptorSet, 0, nullptr);

    uint32 startVertex = 0;
    uint32 startInstance = 0;
    for (uint32 i = 0; i < meshPipeline.MAX_MESHES; i++) {
        if (meshNumInstances[i] > 0) {
            vkCmdDraw(commandBuffer, meshPipeline.numVertices[i], meshNumInstances[i], startVertex, startInstance);
            startInstance += meshNumInstances[i];
        }
        startVertex += meshPipeline.numVertices[i];
    }
}

void UploadAndSubmitLightmapMeshDrawCommands(VkDevice device, VkCommandBuffer commandBuffer,
                                             const VulkanLightmapMeshPipeline& lightmapMeshPipeline,
                                             Mat4 model, Mat4 view, Mat4 proj)
{
    const LightmapMeshUniformBufferObject ubo = {
        .model = model,
        .view = view,
        .proj = proj
    };

    void* data;
    vkMapMemory(device, lightmapMeshPipeline.uniformBuffer.memory, 0, sizeof(LightmapMeshUniformBufferObject), 0, &data);
    MemCopy(data, &ubo, sizeof(ubo));
    vkUnmapMemory(device, lightmapMeshPipeline.uniformBuffer.memory);

    vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, lightmapMeshPipeline.pipeline);

    const VkBuffer vertexBuffers[] = { lightmapMeshPipeline.vertexBuffer.buffer };
    const VkDeviceSize offsets[] = { 0 };
    vkCmdBindVertexBuffers(commandBuffer, 0, C_ARRAY_LENGTH(vertexBuffers), vertexBuffers, offsets);

    uint32 startTriangleInd = 0;
    for (uint32 i = 0; i < lightmapMeshPipeline.meshTriangleEndInds.size; i++) {
        vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, lightmapMeshPipeline.pipelineLayout, 0, 1,
                                &lightmapMeshPipeline.descriptorSets[i], 0, nullptr);

        const uint32 numTriangles = lightmapMeshPipeline.meshTriangleEndInds[i] - startTriangleInd;
        vkCmdDraw(commandBuffer, numTriangles * 3, 1, startTriangleInd * 3, 0);

        startTriangleInd = lightmapMeshPipeline.meshTriangleEndInds[i];
    }
}

bool LoadMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, LinearAllocator* allocator,
                               VulkanMeshPipeline* meshPipeline)
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

    VkVertexInputBindingDescription bindingDescriptions[2] = {};
    VkVertexInputAttributeDescription attributeDescriptions[10] = {};

    // Per-vertex attribute bindings
    bindingDescriptions[0].binding = 0;
    bindingDescriptions[0].stride = sizeof(VulkanMeshVertex);
    bindingDescriptions[0].inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

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
    attributeDescriptions[2].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[2].offset = offsetof(VulkanMeshVertex, color);

    // Per-instance attribute bindings
    bindingDescriptions[1].binding = 1;
    bindingDescriptions[1].stride = sizeof(VulkanMeshInstanceData);
    bindingDescriptions[1].inputRate = VK_VERTEX_INPUT_RATE_INSTANCE;

    attributeDescriptions[3].binding = 1;
    attributeDescriptions[3].location = 3;
    attributeDescriptions[3].format = VK_FORMAT_R32G32B32A32_SFLOAT;
    attributeDescriptions[3].offset = 0 * sizeof(float32);
    attributeDescriptions[4].binding = 1;
    attributeDescriptions[4].location = 4;
    attributeDescriptions[4].format = VK_FORMAT_R32G32B32A32_SFLOAT;
    attributeDescriptions[4].offset = 4 * sizeof(float32);
    attributeDescriptions[5].binding = 1;
    attributeDescriptions[5].location = 5;
    attributeDescriptions[5].format = VK_FORMAT_R32G32B32A32_SFLOAT;
    attributeDescriptions[5].offset = 8 * sizeof(float32);
    attributeDescriptions[6].binding = 1;
    attributeDescriptions[6].location = 6;
    attributeDescriptions[6].format = VK_FORMAT_R32G32B32A32_SFLOAT;
    attributeDescriptions[6].offset = 12 * sizeof(float32);

    attributeDescriptions[7].binding = 1;
    attributeDescriptions[7].location = 7;
    attributeDescriptions[7].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[7].offset = offsetof(VulkanMeshInstanceData, color);

    attributeDescriptions[8].binding = 1;
    attributeDescriptions[8].location = 8;
    attributeDescriptions[8].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[8].offset = offsetof(VulkanMeshInstanceData, collapseMid);

    attributeDescriptions[9].binding = 1;
    attributeDescriptions[9].location = 9;
    attributeDescriptions[9].format = VK_FORMAT_R32_SFLOAT;
    attributeDescriptions[9].offset = offsetof(VulkanMeshInstanceData, collapseT);

    VkPipelineVertexInputStateCreateInfo vertexInputCreateInfo = {};
    vertexInputCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO;
    vertexInputCreateInfo.vertexBindingDescriptionCount = C_ARRAY_LENGTH(bindingDescriptions);
    vertexInputCreateInfo.pVertexBindingDescriptions = bindingDescriptions;
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

    VkPipelineColorBlendAttachmentState colorBlendAttachment = {};
    colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    colorBlendAttachment.blendEnable = VK_FALSE;
    colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;
    colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ZERO;
    colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
    colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
    colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
    colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

    VkPipelineColorBlendStateCreateInfo colorBlendingCreateInfo = {};
    colorBlendingCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
    colorBlendingCreateInfo.logicOpEnable = VK_FALSE;
    colorBlendingCreateInfo.logicOp = VK_LOGIC_OP_COPY;
    colorBlendingCreateInfo.attachmentCount = 1;
    colorBlendingCreateInfo.pAttachments = &colorBlendAttachment;
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
    pipelineCreateInfo.renderPass = swapchain.renderPass;
    pipelineCreateInfo.subpass = 0;
    pipelineCreateInfo.basePipelineHandle = VK_NULL_HANDLE;
    pipelineCreateInfo.basePipelineIndex = -1;

    if (vkCreateGraphicsPipelines(window.device, VK_NULL_HANDLE, 1, &pipelineCreateInfo, nullptr,
                                  &meshPipeline->pipeline) != VK_SUCCESS) {
        LOG_ERROR("vkCreateGraphicsPipeline failed\n");
        return false;
    }

    return true;
}

void UnloadMeshPipelineSwapchain(VkDevice device, VulkanMeshPipeline* meshPipeline)
{
    vkDestroyPipeline(device, meshPipeline->pipeline, nullptr);
    vkDestroyPipelineLayout(device, meshPipeline->pipelineLayout, nullptr);
}

bool LoadMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                            VulkanMeshPipeline* meshPipeline)
{
    // Create vertex buffers
    {
        ALLOCATOR_SCOPE_RESET(*allocator);
        DynamicArray<VulkanMeshVertex, LinearAllocator> vertices(allocator);

        const Vec3 color = Vec3::one;
        const StaticArray<VulkanMeshVertex, 6> tileVertices = {
            .data = {
                { { 0.0f, 0.0f, 0.0f }, Vec3::unitZ, color },
                { { 0.0f, 1.0f, 0.0f }, Vec3::unitZ, color },
                { { 1.0f, 1.0f, 0.0f }, Vec3::unitZ, color },
                { { 1.0f, 1.0f, 0.0f }, Vec3::unitZ, color },
                { { 1.0f, 0.0f, 0.0f }, Vec3::unitZ, color },
                { { 0.0f, 0.0f, 0.0f }, Vec3::unitZ, color },
            }
        };

        static_assert((uint32)MeshId::TILE == 0);
        meshPipeline->numVertices[0] = tileVertices.SIZE;
        vertices.Append(tileVertices.ToArray());

        struct ObjMesh {
            MeshId id;
            const_string path;
        };

        const ObjMesh objMeshes[] = {
            { MeshId::MOB, ToString("data/models/enemy1.obj") },
        };

        for (uint32 i = 0; i < C_ARRAY_LENGTH(objMeshes); i++) {
            LoadObjResult obj;
            if (!LoadObj(objMeshes[i].path, &obj, allocator)) {
                LOG_ERROR("Failed to load rockie1 obj file\n");
                return false;
            }

            VulkanMeshGeometry geometry = ObjToVulkanMeshGeometry(obj, allocator);
            if (!geometry.valid) {
                LOG_ERROR("Failed to load Vulkan geometry from obj\n");
                return false;
            }

            const Array<VulkanMeshVertex> objVertices = {
                .size = geometry.triangles.size * 3,
                .data = &geometry.triangles[0][0]
            };

            const uint32 meshIndex = (uint32)objMeshes[i].id;
            meshPipeline->numVertices[meshIndex] = objVertices.size;
            vertices.Append(objVertices);
        }

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
        CopyBuffer(window.device, commandPool, window.graphicsQueue, stagingBuffer.buffer,
                   meshPipeline->vertexBuffer.buffer, vertexBufferSize);

        DestroyVulkanBuffer(window.device, &stagingBuffer);
    }

    // Create instance buffer
    {
        const VkDeviceSize bufferSize = meshPipeline->MAX_INSTANCES * sizeof(VulkanMeshInstanceData);

        if (!CreateVulkanBuffer(bufferSize,
                                VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &meshPipeline->instanceBuffer)) {
            LOG_ERROR("CreateBuffer failed for instance buffer\n");
            return false;
        }
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

    return true;
}

void UnloadMeshPipelineWindow(VkDevice device, VulkanMeshPipeline* meshPipeline)
{
    vkDestroyDescriptorPool(device, meshPipeline->descriptorPool, nullptr);
    vkDestroyDescriptorSetLayout(device, meshPipeline->descriptorSetLayout, nullptr);

    DestroyVulkanBuffer(device, &meshPipeline->uniformBuffer);
    DestroyVulkanBuffer(device, &meshPipeline->instanceBuffer);
    DestroyVulkanBuffer(device, &meshPipeline->vertexBuffer);
}

bool LoadLightmapMeshPipelineSwapchain(const VulkanWindow& window, const VulkanSwapchain& swapchain, LinearAllocator* allocator, VulkanLightmapMeshPipeline* lightmapMeshPipeline)
{
    const Array<uint8> vertShaderCode = LoadEntireFile(ToString("data/shaders/lightmapMesh.vert.spv"), allocator);
    if (vertShaderCode.data == nullptr) {
        LOG_ERROR("Failed to load vertex shader code\n");
        return false;
    }
    const Array<uint8> fragShaderCode = LoadEntireFile(ToString("data/shaders/lightmapMesh.frag.spv"), allocator);
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
    bindingDescription.stride = sizeof(VulkanLightmapMeshVertex);
    bindingDescription.inputRate = VK_VERTEX_INPUT_RATE_VERTEX;

    VkVertexInputAttributeDescription attributeDescriptions[5] = {};
    attributeDescriptions[0].binding = 0;
    attributeDescriptions[0].location = 0;
    attributeDescriptions[0].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[0].offset = offsetof(VulkanLightmapMeshVertex, pos);

    attributeDescriptions[1].binding = 0;
    attributeDescriptions[1].location = 1;
    attributeDescriptions[1].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[1].offset = offsetof(VulkanLightmapMeshVertex, normal);

    attributeDescriptions[2].binding = 0;
    attributeDescriptions[2].location = 2;
    attributeDescriptions[2].format = VK_FORMAT_R32G32B32_SFLOAT;
    attributeDescriptions[2].offset = offsetof(VulkanLightmapMeshVertex, color);

    attributeDescriptions[3].binding = 0;
    attributeDescriptions[3].location = 3;
    attributeDescriptions[3].format = VK_FORMAT_R32G32_SFLOAT;
    attributeDescriptions[3].offset = offsetof(VulkanLightmapMeshVertex, uv);

    attributeDescriptions[4].binding = 0;
    attributeDescriptions[4].location = 4;
    attributeDescriptions[4].format = VK_FORMAT_R32_SFLOAT;
    attributeDescriptions[4].offset = offsetof(VulkanLightmapMeshVertex, lightmapWeight);

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

    VkPipelineColorBlendAttachmentState colorBlendAttachment = {};
    colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    colorBlendAttachment.blendEnable = VK_FALSE;
    colorBlendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;
    colorBlendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ZERO;
    colorBlendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
    colorBlendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
    colorBlendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
    colorBlendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;

    VkPipelineColorBlendStateCreateInfo colorBlendingCreateInfo = {};
    colorBlendingCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO;
    colorBlendingCreateInfo.logicOpEnable = VK_FALSE;
    colorBlendingCreateInfo.logicOp = VK_LOGIC_OP_COPY;
    colorBlendingCreateInfo.attachmentCount = 1;
    colorBlendingCreateInfo.pAttachments = &colorBlendAttachment;
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
    pipelineLayoutCreateInfo.pSetLayouts = &lightmapMeshPipeline->descriptorSetLayout;
    pipelineLayoutCreateInfo.pushConstantRangeCount = 0;
    pipelineLayoutCreateInfo.pPushConstantRanges = nullptr;

    if (vkCreatePipelineLayout(window.device, &pipelineLayoutCreateInfo, nullptr,
                               &lightmapMeshPipeline->pipelineLayout) != VK_SUCCESS) {
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
    pipelineCreateInfo.layout = lightmapMeshPipeline->pipelineLayout;
    pipelineCreateInfo.renderPass = swapchain.renderPass;
    pipelineCreateInfo.subpass = 0;
    pipelineCreateInfo.basePipelineHandle = VK_NULL_HANDLE;
    pipelineCreateInfo.basePipelineIndex = -1;

    if (vkCreateGraphicsPipelines(window.device, VK_NULL_HANDLE, 1, &pipelineCreateInfo, nullptr,
                                  &lightmapMeshPipeline->pipeline) != VK_SUCCESS) {
        LOG_ERROR("vkCreateGraphicsPipeline failed\n");
        return false;
    }

    return true;
}

void UnloadLightmapMeshPipelineSwapchain(VkDevice device, VulkanLightmapMeshPipeline* lightmapMeshPipeline)
{
    vkDestroyPipeline(device, lightmapMeshPipeline->pipeline, nullptr);
    vkDestroyPipelineLayout(device, lightmapMeshPipeline->pipelineLayout, nullptr);
}

bool LoadLightmapMeshPipelineWindow(const VulkanWindow& window, VkCommandPool commandPool, LinearAllocator* allocator,
                                    VulkanLightmapMeshPipeline* lightmapMeshPipeline)
{
    // Load vulkan vertex geometry
    VulkanLightmapMeshGeometry geometry;
    {
        LoadObjResult obj;
        if (!LoadObj(ToString("data/models/light-cones.obj"), &obj, allocator)) {
            LOG_ERROR("Failed to load reference scene .obj\n");
            return false;
        }

        geometry = ObjToVulkanLightmapMeshGeometry(obj, allocator);
        if (!geometry.valid) {
            LOG_ERROR("Failed to load Vulkan geometry from obj\n");
            return false;
        }

        // Set per-vertex lightmap weights based on triangle areas
        for (uint32 i = 0; i < geometry.triangles.size; i++) {
            VulkanLightmapMeshTriangle& t = geometry.triangles[i];
            const float32 area = TriangleArea(t[0].pos, t[1].pos, t[2].pos);
            const float32 weight = ClampFloat32(SmoothStep(0.0f, 0.005f, area), 0.0f, 1.0f);
            for (int j = 0; j < 3; j++) {
                t[j].lightmapWeight = weight;
            }
        }

        // Load vertex colors from lightmap data
        uint32 startInd = 0;
        for (uint32 i = 0; i < geometry.meshEndInds.size; i++) {
            const_string filePath = AllocPrintf(allocator, "data/lightmaps/%llu.v", i);
            Array<uint8> vertexColors = LoadEntireFile(filePath, allocator);
            if (vertexColors.data == nullptr) {
                LOG_ERROR("Failed to load vertex colors for mesh %lu\n", i);
                return false;
            }

            const bool sizeEvenVec3 = vertexColors.size % sizeof(Vec3) == 0;
            const bool sizeEvenTriangles = (vertexColors.size / sizeof(Vec3)) % 3 == 0;
            if (!sizeEvenVec3 || !sizeEvenTriangles) {
                LOG_ERROR("Incorrect format for vertex colors at %.*s, mesh %lu\n", filePath.size, filePath.data, i);
                return false;
            }

            const uint32 expectedColors = (geometry.meshEndInds[i] - startInd) * 3;
            const uint32 numColors = vertexColors.size / sizeof(Vec3);
            if (expectedColors != numColors) {
                LOG_ERROR("Mismatched number of vertex colors, expected %lu, got %lu\n", expectedColors, numColors);
                return false;
            }

            const Array<Vec3> colors = {
                .size = vertexColors.size / sizeof(Vec3),
                .data = (Vec3*)vertexColors.data
            };
            for (uint32 j = startInd; j < geometry.meshEndInds[i]; j++) {
                const uint32 colorInd = (j - startInd) * 3;
                geometry.triangles[j][0].color = colors[colorInd];
                geometry.triangles[j][1].color = colors[colorInd + 1];
                geometry.triangles[j][2].color = colors[colorInd + 2];
            }

            startInd = geometry.meshEndInds[i];
        }

        // Save mesh triangle end inds to VulkanApp structure for draw commands to use
        for (uint32 i = 0; i < geometry.meshEndInds.size; i++) {
            lightmapMeshPipeline->meshTriangleEndInds.Append(geometry.meshEndInds[i]);
        }
    }

    // Create vertex buffer
    {
        const VkDeviceSize vertexBufferSize = geometry.triangles.size * 3 * sizeof(VulkanLightmapMeshVertex);

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
        MemCopy(data, geometry.triangles.data, vertexBufferSize);
        vkUnmapMemory(window.device, stagingBuffer.memory);

        if (!CreateVulkanBuffer(vertexBufferSize,
                                VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                                window.device, window.physicalDevice, &lightmapMeshPipeline->vertexBuffer)) {
            LOG_ERROR("CreateBuffer failed for vertex buffer\n");
            return false;
        }

        // Copy vertex data from staging buffer into GPU vertex buffer
        CopyBuffer(window.device, commandPool, window.graphicsQueue, stagingBuffer.buffer,
                   lightmapMeshPipeline->vertexBuffer.buffer, vertexBufferSize);

        DestroyVulkanBuffer(window.device, &stagingBuffer);
    }

    // Create lightmaps
    {
        for (uint32 i = 0; i < geometry.meshEndInds.size; i++) {
            const char* filePath = ToCString(AllocPrintf(allocator, "data/lightmaps/%llu.png", i), allocator);
            int width, height, channels;
            unsigned char* imageData = stbi_load(filePath, &width, &height, &channels, 0);
            if (imageData == NULL) {
                LOG_ERROR("Failed to load lightmap: %s\n", filePath);
                return false;
            }
            defer(stbi_image_free(imageData));

            VulkanImage* lightmapImage = lightmapMeshPipeline->lightmaps.Append();
            if (!LoadVulkanImage(window.device, window.physicalDevice, window.graphicsQueue, commandPool,
                                 width, height, channels, (const uint8*)imageData, lightmapImage)) {
                LOG_ERROR("Failed to Vulkan image for lightmap %s\n", filePath);
                return false;
            }
        }
    }

    // Create lightmap sampler
    {
        VkSamplerCreateInfo createInfo = {};
        createInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
        createInfo.magFilter = LIGHTMAP_TEXTURE_FILTER;
        createInfo.minFilter = LIGHTMAP_TEXTURE_FILTER;
        createInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        createInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        createInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_REPEAT;
        createInfo.anisotropyEnable = VK_FALSE;
        createInfo.maxAnisotropy = 1.0f;
        createInfo.borderColor = VK_BORDER_COLOR_INT_OPAQUE_BLACK;
        createInfo.unnormalizedCoordinates = VK_FALSE;
        createInfo.compareEnable = VK_FALSE;
        createInfo.compareOp = VK_COMPARE_OP_ALWAYS;
        createInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
        createInfo.mipLodBias = 0.0f;
        createInfo.minLod = 0.0f;
        createInfo.maxLod = 0.0f;

        if (vkCreateSampler(window.device, &createInfo, nullptr, &lightmapMeshPipeline->lightmapSampler) != VK_SUCCESS) {
            LOG_ERROR("vkCreateSampler failed\n");
            return false;
        }
    }

    // Create uniform buffer
    {
        VkDeviceSize uniformBufferSize = sizeof(LightmapMeshUniformBufferObject);
        if (!CreateVulkanBuffer(uniformBufferSize,
                                VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                window.device, window.physicalDevice, &lightmapMeshPipeline->uniformBuffer)) {
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

        VkDescriptorSetLayoutBinding samplerLayoutBinding = {};
        samplerLayoutBinding.binding = 1;
        samplerLayoutBinding.descriptorCount = 1;
        samplerLayoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        samplerLayoutBinding.pImmutableSamplers = nullptr;
        samplerLayoutBinding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;

        const VkDescriptorSetLayoutBinding bindings[] = { uboLayoutBinding, samplerLayoutBinding };

        VkDescriptorSetLayoutCreateInfo layoutCreateInfo = {};
        layoutCreateInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
        layoutCreateInfo.bindingCount = C_ARRAY_LENGTH(bindings);
        layoutCreateInfo.pBindings = bindings;

        if (vkCreateDescriptorSetLayout(window.device, &layoutCreateInfo, nullptr,
                                        &lightmapMeshPipeline->descriptorSetLayout) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorSetLayout failed\n");
            return false;
        }
    }

    // Create descriptor pool
    {
        VkDescriptorPoolSize poolSizes[2] = {};
        poolSizes[0].type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
        poolSizes[0].descriptorCount = geometry.meshEndInds.size;
        poolSizes[1].type = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
        poolSizes[1].descriptorCount = geometry.meshEndInds.size;

        VkDescriptorPoolCreateInfo poolInfo = {};
        poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        poolInfo.poolSizeCount = C_ARRAY_LENGTH(poolSizes);
        poolInfo.pPoolSizes = poolSizes;
        poolInfo.maxSets = geometry.meshEndInds.size;

        if (vkCreateDescriptorPool(window.device, &poolInfo, nullptr, &lightmapMeshPipeline->descriptorPool) != VK_SUCCESS) {
            LOG_ERROR("vkCreateDescriptorPool failed\n");
            return false;
        }
    }

    // Create descriptor set
    {
        FixedArray<VkDescriptorSetLayout, lightmapMeshPipeline->MAX_LIGHTMAPS> layouts;
        layouts.Clear();
        for (uint32 i = 0; i < geometry.meshEndInds.size; i++) {
            layouts.Append(lightmapMeshPipeline->descriptorSetLayout);
        }

        VkDescriptorSetAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        allocInfo.descriptorPool = lightmapMeshPipeline->descriptorPool;
        allocInfo.descriptorSetCount = geometry.meshEndInds.size;
        allocInfo.pSetLayouts = layouts.data;

        if (vkAllocateDescriptorSets(window.device, &allocInfo, lightmapMeshPipeline->descriptorSets.data) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateDescriptorSets failed\n");
            return false;
        }
        lightmapMeshPipeline->descriptorSets.size = geometry.meshEndInds.size;

        for (uint32 i = 0; i < geometry.meshEndInds.size; i++) {
            VkWriteDescriptorSet descriptorWrites[2] = {};

            VkDescriptorBufferInfo bufferInfo = {};
            bufferInfo.buffer = lightmapMeshPipeline->uniformBuffer.buffer;
            bufferInfo.offset = 0;
            bufferInfo.range = sizeof(LightmapMeshUniformBufferObject);

            descriptorWrites[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrites[0].dstSet = lightmapMeshPipeline->descriptorSets[i];
            descriptorWrites[0].dstBinding = 0;
            descriptorWrites[0].dstArrayElement = 0;
            descriptorWrites[0].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
            descriptorWrites[0].descriptorCount = 1;
            descriptorWrites[0].pBufferInfo = &bufferInfo;

            VkDescriptorImageInfo imageInfo = {};
            imageInfo.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            imageInfo.imageView = lightmapMeshPipeline->lightmaps[i].view;
            imageInfo.sampler = lightmapMeshPipeline->lightmapSampler;

            descriptorWrites[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrites[1].dstSet = lightmapMeshPipeline->descriptorSets[i];
            descriptorWrites[1].dstBinding = 1;
            descriptorWrites[1].dstArrayElement = 0;
            descriptorWrites[1].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
            descriptorWrites[1].descriptorCount = 1;
            descriptorWrites[1].pImageInfo = &imageInfo;

            vkUpdateDescriptorSets(window.device, C_ARRAY_LENGTH(descriptorWrites), descriptorWrites, 0, nullptr);
        }
    }

    return true;
}

void UnloadLightmapMeshPipelineWindow(VkDevice device, VulkanLightmapMeshPipeline* lightmapMeshPipeline)
{
    lightmapMeshPipeline->descriptorSets.Clear();
    vkDestroyDescriptorPool(device, lightmapMeshPipeline->descriptorPool, nullptr);
    vkDestroyDescriptorSetLayout(device, lightmapMeshPipeline->descriptorSetLayout, nullptr);

    DestroyVulkanBuffer(device, &lightmapMeshPipeline->uniformBuffer);

    vkDestroySampler(device, lightmapMeshPipeline->lightmapSampler, nullptr);

    for (uint32 i = 0; i < lightmapMeshPipeline->lightmaps.size; i++) {
        vkDestroyImageView(device, lightmapMeshPipeline->lightmaps[i].view, nullptr);
        vkDestroyImage(device, lightmapMeshPipeline->lightmaps[i].image, nullptr);
        vkFreeMemory(device, lightmapMeshPipeline->lightmaps[i].memory, nullptr);
    }
    lightmapMeshPipeline->lightmaps.Clear();

    DestroyVulkanBuffer(device, &lightmapMeshPipeline->vertexBuffer);

    lightmapMeshPipeline->meshTriangleEndInds.Clear();
}
