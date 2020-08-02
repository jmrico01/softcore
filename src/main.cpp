#include "main.h"

#include <stdio.h>

#include <stb_image.h>
#include <Tracy.hpp>

#include <km_common/km_array.h>
#include <km_common/km_defines.h>
#include <km_common/km_load_obj.h>
#include <km_common/km_os.h>
#include <km_common/km_string.h>
#include <km_common/app/km_app.h>

#define ENABLE_THREADS 1

#define BEGIN_INTERIOR 1

// Required for platform main
const char* WINDOW_NAME = "softcore";
const int WINDOW_START_WIDTH  = 1024;
const int WINDOW_START_HEIGHT = 768;
const bool WINDOW_LOCK_CURSOR = true;
const uint64 PERMANENT_MEMORY_SIZE = MEGABYTES(256);
const uint64 TRANSIENT_MEMORY_SIZE = GIGABYTES(1);

internal AppState* GetAppState(AppMemory* memory)
{
    static_assert(sizeof(AppState) < PERMANENT_MEMORY_SIZE);
    DEBUG_ASSERT(sizeof(AppState) < memory->permanent.size);

    AppState* appState = (AppState*)memory->permanent.data;
    appState->arena = {
        .size = memory->permanent.size - sizeof(AppState),
        .data = memory->permanent.data + sizeof(AppState),
    };
    return appState;
}

internal TransientState* GetTransientState(AppMemory* memory)
{
    static_assert(sizeof(TransientState) < TRANSIENT_MEMORY_SIZE);
    DEBUG_ASSERT(sizeof(TransientState) < memory->transient.size);

    TransientState* transientState = (TransientState*)memory->transient.data;
    transientState->scratch = {
        .size = memory->transient.size - sizeof(TransientState),
        .data = memory->transient.data + sizeof(TransientState),
    };
    return transientState;
}

internal bool LoadScene(const_string scene, AppState* appState, LinearAllocator* allocator)
{
    const_string scenePath = AllocPrintf(allocator, "data/models/%.*s.obj", scene.size, scene.data);

    LoadObjResult obj;
    if (!LoadObj(scenePath, Vec3::zero, 1.0f, &obj, allocator)) {
        LOG_ERROR("Failed to load obj file \"%.*s\" for scene %.*s\n", scenePath.size, scenePath.data,
                  scene.size, scene.data);
        return false;
    }

    appState->arenaAllocator.Clear();

    const uint32 boxMaxTriangles = 32;
    RaycastGeometry geometry = CreateRaycastGeometry(obj, boxMaxTriangles, &appState->arenaAllocator, allocator);
    if (geometry.meshes.data == nullptr) {
        LOG_ERROR("Failed to construct raycast geometry from obj for scene %.*s\n", scene.size, scene.data);
        return false;
    }

    uint32 totalTriangles = 0;
    for (uint32 i = 0; i < geometry.meshes.size; i++) {
        totalTriangles += geometry.meshes[i].numTriangles;
    }
    LOG_INFO("Loaded raycast geometry for scene \"%.*s\" - %lu meshes, %lu total triangles\n",
             scene.size, scene.data, geometry.meshes.size, totalTriangles);
    appState->raycastGeometry = geometry;

    return true;
}

APP_UPDATE_AND_RENDER_FUNCTION(AppUpdateAndRender)
{
    UNREFERENCED_PARAMETER(audio);
    FrameMark;

    AppState* appState = GetAppState(memory);
    TransientState* transientState = GetTransientState(memory);

    const Vec2Int screenSize = {
        (int)vulkanState.swapchain.extent.width,
        (int)vulkanState.swapchain.extent.height
    };
    DEBUG_ASSERT(screenSize.x < CanvasState::MAX_WIDTH);
    DEBUG_ASSERT(screenSize.y < CanvasState::MAX_HEIGHT);

    const float32 CAMERA_HEIGHT = 1.7f;

    // Initialize memory if necessary
    if (!memory->initialized) {
        ZoneScopedN("InitMemory");

        appState->arena = {
            .size = memory->permanent.size - sizeof(AppState),
            .data = memory->permanent.data + sizeof(AppState),
        };
        appState->arenaAllocator = LinearAllocator(appState->arena);

#if BEGIN_INTERIOR
        // Scene - interiors
        appState->cameraPos = Vec3 { -1.15f, -2.21f, 0.20f };
        appState->cameraAngles = Vec2 { -1.25f, 0.03f };
#else
        // Scene - light cones
        appState->cameraPos = Vec3 { 3.74f, -5.21f, 2.49f };
        appState->cameraAngles = Vec2 { -2.26f, -0.77f };
#endif

        appState->canvas.screenFill = 1.0f / 60.0f;
        appState->canvas.decayFrames = 255;
        appState->canvas.samples = 8;
        appState->canvas.bounces = 4;

        MemSet(appState->canvas.colorHdr.data, 0, CanvasState::MAX_PIXELS * sizeof(Vec3));
        MemSet(appState->canvas.decay.data, 0, CanvasState::MAX_PIXELS * sizeof(uint8));

        {
            LinearAllocator allocator(transientState->scratch);

#if BEGIN_INTERIOR
            LoadScene(ToString("interior-1"), appState, &allocator);
#else
            LoadScene(ToString("light-cones"), appState, &allocator);
#endif
        }

        // Debug views 
        appState->debugView = false;
        LockCursor(false);

        appState->inputScreenFillState.value = appState->canvas.screenFill;
        appState->inputDecayFramesState.Initialize(appState->canvas.decayFrames);
        appState->inputSamplesState.Initialize(appState->canvas.samples);
        appState->inputBouncesState.Initialize(appState->canvas.bounces);

        memory->initialized = true;
    }

    // Reset frame state
    {
        ResetSpriteRenderState(&transientState->frameState.spriteRenderState);
        ResetTextRenderState(&transientState->frameState.textRenderState);
    }

    appState->elapsedTime += deltaTime;

    if (KeyPressed(input, KM_KEY_G)) {
        appState->debugView = !appState->debugView;
        if (!appState->debugView) {
            // LockCursor(true);
        }
    }

    if (KeyPressed(input, KM_KEY_ESCAPE)) {
        bool cursorLocked = IsCursorLocked();
        LockCursor(!cursorLocked);
    }

    const float32 cameraSensitivity = 2.0f;
    const Vec2 mouseDeltaFrac = {
        (float32)input.mouseDelta.x / (float32)screenSize.x,
        (float32)input.mouseDelta.y / (float32)screenSize.y
    };
    if (IsCursorLocked() || MouseDown(input, KM_MOUSE_RIGHT)) {
        appState->cameraAngles.x += mouseDeltaFrac.x * cameraSensitivity;
        appState->cameraAngles.y -= mouseDeltaFrac.y * cameraSensitivity;
    }

    appState->cameraAngles.x = ModFloat32(appState->cameraAngles.x, PI_F * 2.0f);
    appState->cameraAngles.y = ClampFloat32(appState->cameraAngles.y, -PI_F / 2.0f, PI_F / 2.0f);

    const Quat cameraRotYaw = QuatFromAngleUnitAxis(appState->cameraAngles.x, Vec3::unitZ);
    const Quat cameraRotPitch = QuatFromAngleUnitAxis(appState->cameraAngles.y, Vec3::unitY);

    const Quat cameraRotYawInv = Inverse(cameraRotYaw);
    const Vec3 cameraForwardXY = cameraRotYawInv * Vec3::unitX;
    const Vec3 cameraRightXY = cameraRotYawInv * -Vec3::unitY;
    const Vec3 cameraUp = Vec3::unitZ;

    Vec3 velocity = Vec3::zero;
    if (KeyDown(input, KM_KEY_W)) {
        velocity += cameraForwardXY;
    }
    if (KeyDown(input, KM_KEY_S)) {
        velocity -= cameraForwardXY;
    }
    if (KeyDown(input, KM_KEY_A)) {
        velocity -= cameraRightXY;
    }
    if (KeyDown(input, KM_KEY_D)) {
        velocity += cameraRightXY;
    }
    if (KeyDown(input, KM_KEY_SPACE)) {
        velocity += cameraUp;
    }
    if (KeyDown(input, KM_KEY_CTRL)) {
        velocity -= cameraUp;
    }

    if (velocity != Vec3::zero) {
        float32 speed = 0.7f;
        if (KeyDown(input, KM_KEY_SHIFT)) {
            speed *= 2.0f;
        }

        appState->cameraPos += velocity * speed * deltaTime;
    }

    const Quat cameraRot = cameraRotPitch * cameraRotYaw;
    const Mat4 cameraRotMat4 = UnitQuatToMat4(cameraRot);

    const Quat inverseCameraRot = Inverse(cameraRot);
    const Vec3 cameraForward = inverseCameraRot * Vec3::unitX;

    // Transforms world-view camera (+X forward, +Z up) to Vulkan camera (+Z forward, -Y up)
    const Quat baseCameraRot = QuatFromAngleUnitAxis(-PI_F / 2.0f, Vec3::unitY)
        * QuatFromAngleUnitAxis(PI_F / 2.0f, Vec3::unitX);
    const Mat4 baseCameraRotMat4 = UnitQuatToMat4(baseCameraRot);

    const Mat4 view = baseCameraRotMat4 * cameraRotMat4 * Translate(-appState->cameraPos);

    const float32 aspect = (float32)screenSize.x / (float32)screenSize.y;
    const float32 nearZ = 0.1f;
    const float32 farZ = 100.0f;
    const Mat4 proj = Perspective(PI_F / 4.0f, aspect, nearZ, farZ);

    // Debug views
    if (appState->debugView) {
        LinearAllocator allocator(transientState->scratch);

        const VulkanFontFace& fontNormal = appState->fontFaces[(uint32)FontId::OCR_A_REGULAR_18];
        // const VulkanFontFace& fontTitle = appState->fontFaces[(uint32)FontId::OCR_A_REGULAR_24];

        const Vec4 backgroundColor = Vec4 { 0.0f, 0.0f, 0.0f, 0.5f };
        const Vec4 inputTextColor = { 1.0f, 1.0f, 1.0f, 1.0f };
        const Vec2Int panelBorderSize = Vec2Int { 6, 8 };
        const int panelPosMargin = 50;

        // General debug info
        const Vec2Int panelDebugInfoPos = { panelPosMargin, panelPosMargin };

        Panel panelDebugInfo(&allocator);
        panelDebugInfo.Begin(input, &fontNormal, PanelFlag::GROW_DOWNWARDS, panelDebugInfoPos, 0.0f);

        const int fps = (int)(1.0f / deltaTime);
        const float32 frameMs = deltaTime * 1000.0f;
        const_string frameTiming = AllocPrintf(&allocator, "%d fps | %.03f ms", fps, frameMs);
        panelDebugInfo.Text(frameTiming);
        const_string resolution = AllocPrintf(&allocator, "%d x %d", screenSize.x, screenSize.y);
        panelDebugInfo.Text(resolution);

        panelDebugInfo.Text(string::empty);

        const_string cameraPosString = AllocPrintf(&allocator, "%.02f, %.02f, %.02f POS",
                                                   appState->cameraPos.x, appState->cameraPos.y, appState->cameraPos.z);
        panelDebugInfo.Text(cameraPosString);
        const_string cameraAnglesString = AllocPrintf(&allocator, "%.02f, %.02f ROT",
                                                      appState->cameraAngles.x, appState->cameraAngles.y);
        panelDebugInfo.Text(cameraAnglesString);

        panelDebugInfo.Text(string::empty);

        bool cursorLocked = IsCursorLocked();
        if (panelDebugInfo.Checkbox(&cursorLocked, ToString("cam lock (ESC)"))) {
            LockCursor(cursorLocked);
        }

        panelDebugInfo.Text(string::empty);

        panelDebugInfo.Text(ToString("Screen fill"));
        if (panelDebugInfo.SliderFloat(&appState->inputScreenFillState, 0.0f, 0.1f)) {
            appState->canvas.screenFill = appState->inputScreenFillState.value;
        }

        panelDebugInfo.Text(ToString("Decay frames (0 - 255)"));
        if (panelDebugInfo.InputInt(&appState->inputDecayFramesState)) {
            const int value = appState->inputDecayFramesState.value;
            if (value > 0 && value <= 255) {
                appState->canvas.decayFrames = (uint8)appState->inputDecayFramesState.value;
            }
        }

        panelDebugInfo.Text(ToString("Samples (at least 1)"));
        if (panelDebugInfo.InputInt(&appState->inputSamplesState)) {
            const int value = appState->inputSamplesState.value;
            if (value > 0) {
                appState->canvas.samples = appState->inputSamplesState.value;
            }
        }

        panelDebugInfo.Text(ToString("Bounces (at least 1)"));
        if (panelDebugInfo.InputInt(&appState->inputBouncesState)) {
            const int value = appState->inputBouncesState.value;
            if (value > 0) {
                appState->canvas.bounces = appState->inputBouncesState.value;
            }
        }

        panelDebugInfo.Text(string::empty);

        const Array<string> modelFiles = ListDir(ToString("data/models"), &allocator);
        if (modelFiles.data != nullptr) {
            DynamicArray<string, LinearAllocator> scenes(&allocator);
            for (uint32 i = 0; i < modelFiles.size; i++) {
                const uint32 find = SubstringSearch(modelFiles[i], ToString(".obj"));
                if (find != modelFiles[i].size) {
                    scenes.Append(modelFiles[i].SliceTo(find));
                }
            }

            panelDebugInfo.Text(ToString("Load scene"));
            if (panelDebugInfo.Dropdown(&appState->inputSceneDropdownState, scenes.ToArray())) {
                const_string scene = scenes[appState->inputSceneDropdownState.selected];
                if (LoadScene(scene, appState, &allocator)) {
                    LOG_INFO("Loaded scene %.*s\n", scene.size, scene.data);
                }
                else {
                    LOG_ERROR("Failed to load scene %.*s\n", scene.size, scene.data);
                }
            }
        }

        panelDebugInfo.Draw(panelBorderSize, Vec4::one, backgroundColor, screenSize,
                            &transientState->frameState.spriteRenderState, &transientState->frameState.textRenderState);
    }

    // Ray traced rendering
    {
        ZoneScopedN("RayTraceRender");

        LinearAllocator allocator(transientState->scratch);

        const uint32 numPixels = screenSize.x * screenSize.y;
        RaytraceRender(appState->cameraPos, cameraRot, appState->raycastGeometry,
                       screenSize.x, screenSize.y, &appState->canvas, &allocator, queue);

        const uint32 numBytes = numPixels * 4;
        void* imageMemory;
        vkMapMemory(vulkanState.window.device, appState->vulkanAppState.imageMemory, 0, numBytes, 0, &imageMemory);

        float32 maxColor = 1.0f;
        for (uint32 i = 0; i < numPixels; i++) {
            const Vec3 colorHdr = appState->canvas.colorHdr[i];
            maxColor = MaxFloat32(maxColor, colorHdr.r);
            maxColor = MaxFloat32(maxColor, colorHdr.g);
            maxColor = MaxFloat32(maxColor, colorHdr.b);
        }

        uint32* pixels = (uint32*)imageMemory;
        for (uint32 i = 0; i < numPixels; i++) {
            const Vec3 colorNormalized = appState->canvas.colorHdr[i] / maxColor;
            // TODO gamma correction?
            const uint8 r = (uint8)(colorNormalized.r * 255.0f);
            const uint8 g = (uint8)(colorNormalized.g * 255.0f);
            const uint8 b = (uint8)(colorNormalized.b * 255.0f);
            pixels[i] = ((uint32)r << 16) + ((uint32)g << 8) + b;
        }
        // MemCopy(imageMemory, appState->pixels.data, numBytes);
        vkUnmapMemory(vulkanState.window.device, appState->vulkanAppState.imageMemory);
    }

    // ================================================================================================
    // Vulkan rendering ===============================================================================
    // ================================================================================================

    VkCommandBuffer buffer = appState->vulkanAppState.commandBuffer;
    VkFence fence = appState->vulkanAppState.fence;

    {
        ZoneScopedN("WaitForFence");

        // TODO revisit this. should the platform coordinate something like this in some other way?
        // swapchain image acquisition timings seem to be kind of sloppy tbh, so this might be the best way.
        if (vkWaitForFences(vulkanState.window.device, 1, &fence, VK_TRUE, UINT64_MAX) != VK_SUCCESS) {
            LOG_ERROR("vkWaitForFences didn't return success for fence %lu\n", swapchainImageIndex);
        }
        if (vkResetFences(vulkanState.window.device, 1, &fence) != VK_SUCCESS) {
            LOG_ERROR("vkResetFences didn't return success for fence %lu\n", swapchainImageIndex);
        }
    }

    if (vkResetCommandBuffer(buffer, 0) != VK_SUCCESS) {
        LOG_ERROR("vkResetCommandBuffer failed\n");
    }

    VkCommandBufferBeginInfo beginInfo = {};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = 0;
    beginInfo.pInheritanceInfo = nullptr;

    if (vkBeginCommandBuffer(buffer, &beginInfo) != VK_SUCCESS) {
        LOG_ERROR("vkBeginCommandBuffer failed\n");
    }

    TransitionImageLayout(buffer, appState->vulkanAppState.image,
                          VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL);
    TransitionImageLayout(buffer, vulkanState.swapchain.images[swapchainImageIndex],
                          VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL);

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
    imageCopy.extent.width = screenSize.x;
    imageCopy.extent.height = screenSize.y;
    imageCopy.extent.depth = 1;

    vkCmdCopyImage(buffer,
                   appState->vulkanAppState.image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                   vulkanState.swapchain.images[swapchainImageIndex], VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                   1, &imageCopy);

    const VkClearValue clearValues[] = {
        { 0.0f, 0.0f, 0.0f, 1.0f },
        { 1.0f, 0 }
    };

    VkRenderPassBeginInfo renderPassInfo = {};
    renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    renderPassInfo.renderPass = vulkanState.swapchain.renderPass;
    renderPassInfo.framebuffer = vulkanState.swapchain.framebuffers[swapchainImageIndex];
    renderPassInfo.renderArea.offset = { 0, 0 };
    renderPassInfo.renderArea.extent = vulkanState.swapchain.extent;
    renderPassInfo.clearValueCount = C_ARRAY_LENGTH(clearValues);
    renderPassInfo.pClearValues = clearValues;

    vkCmdBeginRenderPass(buffer, &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

    // Sprites
    {
        LinearAllocator allocator(transientState->scratch);
        UploadAndSubmitSpriteDrawCommands(vulkanState.window.device, buffer, appState->vulkanAppState.spritePipeline,
                                          transientState->frameState.spriteRenderState, &allocator);
    }

    // Text
    {
        LinearAllocator allocator(transientState->scratch);
        UploadAndSubmitTextDrawCommands(vulkanState.window.device, buffer, appState->vulkanAppState.textPipeline,
                                        transientState->frameState.textRenderState, &allocator);
    }

    vkCmdEndRenderPass(buffer);

    if (vkEndCommandBuffer(buffer) != VK_SUCCESS) {
        LOG_ERROR("vkEndCommandBuffer failed\n");
    }

    const VkSemaphore waitSemaphores[] = { vulkanState.window.imageAvailableSemaphore };
    const VkPipelineStageFlags waitStages[] = { VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT };

    const VkSemaphore signalSemaphores[] = { vulkanState.window.renderFinishedSemaphore };

    VkSubmitInfo submitInfo = {};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.waitSemaphoreCount = C_ARRAY_LENGTH(waitSemaphores);
    submitInfo.pWaitSemaphores = waitSemaphores;
    submitInfo.pWaitDstStageMask = waitStages;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &buffer;
    submitInfo.signalSemaphoreCount = C_ARRAY_LENGTH(signalSemaphores);
    submitInfo.pSignalSemaphores = signalSemaphores;

    if (vkQueueSubmit(vulkanState.window.graphicsQueue, 1, &submitInfo, fence) != VK_SUCCESS) {
        LOG_ERROR("Failed to submit draw command buffer\n");
    }

    return true;
}

APP_LOAD_VULKAN_SWAPCHAIN_STATE_FUNCTION(AppLoadVulkanSwapchainState)
{
    LOG_INFO("Loading Vulkan swapchain-dependent app state\n");

    const VulkanWindow& window = vulkanState.window;
    const VulkanSwapchain& swapchain = vulkanState.swapchain;

    VulkanAppState* app = &(GetAppState(memory)->vulkanAppState);
    TransientState* transientState = GetTransientState(memory);
    LinearAllocator allocator(transientState->scratch);

    // Create image
    if (!CreateImage(window.device, window.physicalDevice, swapchain.extent.width, swapchain.extent.height,
                     VK_FORMAT_B8G8R8A8_SRGB, VK_IMAGE_TILING_LINEAR, VK_IMAGE_USAGE_TRANSFER_SRC_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                     &app->image, &app->imageMemory)) {
        LOG_ERROR("CreateImage failed\n");
        return false;
    }

    if (!LoadSpritePipelineSwapchain(window, swapchain, &allocator, &app->spritePipeline)) {
        LOG_ERROR("Failed to load swapchain-dependent Vulkan sprite pipeline\n");
        return false;
    }

    if (!LoadTextPipelineSwapchain(window, swapchain, &allocator, &app->textPipeline)) {
        LOG_ERROR("Failed to load swapchain-dependent Vulkan text pipeline\n");
        return false;
    }

    return true;
}

APP_UNLOAD_VULKAN_SWAPCHAIN_STATE_FUNCTION(AppUnloadVulkanSwapchainState)
{
    LOG_INFO("Unloading Vulkan swapchain-dependent app state\n");

    const VkDevice& device = vulkanState.window.device;
    VulkanAppState* app = &(GetAppState(memory)->vulkanAppState);

    UnloadTextPipelineSwapchain(device, &app->textPipeline);
    UnloadSpritePipelineSwapchain(device, &app->spritePipeline);

    vkDestroyImage(device, app->image, nullptr);
    vkFreeMemory(device, app->imageMemory, nullptr);
}

APP_LOAD_VULKAN_WINDOW_STATE_FUNCTION(AppLoadVulkanWindowState)
{
    LOG_INFO("Loading Vulkan window-dependent app state\n");

    const VulkanWindow& window = vulkanState.window;

    AppState* appState = GetAppState(memory);
    VulkanAppState* app = &appState->vulkanAppState;
    TransientState* transientState = GetTransientState(memory);
    LinearAllocator allocator(transientState->scratch);

    // Create command pool
    {
        QueueFamilyInfo queueFamilyInfo = GetQueueFamilyInfo(window.surface, window.physicalDevice, &allocator);

        VkCommandPoolCreateInfo poolCreateInfo = {};
        poolCreateInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
        poolCreateInfo.queueFamilyIndex = queueFamilyInfo.graphicsFamilyIndex;
        poolCreateInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

        if (vkCreateCommandPool(window.device, &poolCreateInfo, nullptr, &app->commandPool) != VK_SUCCESS) {
            LOG_ERROR("vkCreateCommandPool failed\n");
            return false;
        }
    }

    // Create command buffer
    {
        VkCommandBufferAllocateInfo allocInfo = {};
        allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
        allocInfo.commandPool = app->commandPool;
        allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        allocInfo.commandBufferCount = 1;

        if (vkAllocateCommandBuffers(window.device, &allocInfo, &app->commandBuffer) != VK_SUCCESS) {
            LOG_ERROR("vkAllocateCommandBuffers failed\n");
            return false;
        }
    }

    // Create fence
    {
        VkFenceCreateInfo fenceCreateInfo = {};
        fenceCreateInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
        fenceCreateInfo.flags = VK_FENCE_CREATE_SIGNALED_BIT;

        if (vkCreateFence(window.device, &fenceCreateInfo, nullptr, &app->fence) != VK_SUCCESS) {
            LOG_ERROR("vkCreateFence failed\n");
            return false;
        }
    }

    const bool spritePipeline = LoadSpritePipelineWindow(window, app->commandPool, &allocator, &app->spritePipeline);
    if (!spritePipeline) {
        LOG_ERROR("Failed to load window-dependent Vulkan sprite pipeline\n");
        return false;
    }

    const bool textPipeline = LoadTextPipelineWindow(window, app->commandPool, &allocator, &app->textPipeline);
    if (!textPipeline) {
        LOG_ERROR("Failed to load window-dependent Vulkan text pipeline\n");
        return false;
    }

    // Sprites
    {
        const char* spriteFilePaths[] = {
            "data/sprites/pixel.png",
            "data/sprites/jon.png",
            "data/sprites/rock.png"
        };
        static_assert(C_ARRAY_LENGTH(spriteFilePaths) == (uint32)SpriteId::COUNT);

        for (uint32 i = 0; i < C_ARRAY_LENGTH(spriteFilePaths); i++) {
            int width, height, channels;
            unsigned char* imageData = stbi_load(spriteFilePaths[i], &width, &height, &channels, 0);
            if (imageData == NULL) {
                DEBUG_PANIC("Failed to load sprite: %s\n", spriteFilePaths[i]);
            }
            defer(stbi_image_free(imageData));

            uint8* vulkanImageData = (uint8*)imageData;
            if (channels == 3) {
                LOG_ERROR("Image %s with 3 channels, converting to RGBA\n", spriteFilePaths[i]);

                vulkanImageData = allocator.New<uint8>(width * height * 4);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        const int imageDataInd = (y * width + x) * 3;
                        vulkanImageData[imageDataInd + 0] = imageData[imageDataInd + 0];
                        vulkanImageData[imageDataInd + 1] = imageData[imageDataInd + 1];
                        vulkanImageData[imageDataInd + 2] = imageData[imageDataInd + 2];
                        vulkanImageData[imageDataInd + 3] = 255;
                    }
                }

                channels = 4;
            }

            VulkanImage sprite;
            if (!LoadVulkanImage(window.device, window.physicalDevice, app->commandPool, window.graphicsQueue,
                                 width, height, channels, vulkanImageData, &sprite)) {
                DEBUG_PANIC("Failed to Vulkan image for sprite %s\n", spriteFilePaths[i]);
            }

            uint32 spriteIndex;
            if (!RegisterSprite(window.device, &app->spritePipeline, sprite, &spriteIndex)) {
                DEBUG_PANIC("Failed to register sprite %s\n", spriteFilePaths[i]);
            }
            DEBUG_ASSERT(spriteIndex == i);
        }
    }

    // Fonts
    {
        struct FontData {
            const_string filePath;
            uint32 height;
        };
        const FontData fontData[] = {
            { ToString("data/fonts/ocr-a/regular.ttf"), 18 },
            { ToString("data/fonts/ocr-a/regular.ttf"), 24 },
        };
        static_assert(C_ARRAY_LENGTH(fontData) == (uint32)FontId::COUNT);

        FT_Library ftLibrary;
        FT_Error error = FT_Init_FreeType(&ftLibrary);
        if (error) {
            DEBUG_PANIC("FreeType init error: %d\n", error);
        }

        for (uint32 i = 0; i < C_ARRAY_LENGTH(fontData); i++) {
            LoadFontFaceResult fontFaceResult;
            if (!LoadFontFace(ftLibrary, fontData[i].filePath, fontData[i].height, &allocator, &fontFaceResult)) {
                DEBUG_PANIC("Failed to load font face at %.*s\n", fontData[i].filePath.size, fontData[i].filePath.data);
            }

            if (!RegisterFont(window.device, window.physicalDevice, window.graphicsQueue, app->commandPool,
                              &app->textPipeline, fontFaceResult, &appState->fontFaces[i])) {
                DEBUG_PANIC("Failed to register font %lu\n", i);
            }
        }
    }

    return true;
}

APP_UNLOAD_VULKAN_WINDOW_STATE_FUNCTION(AppUnloadVulkanWindowState)
{
    LOG_INFO("Unloading Vulkan window-dependent app state\n");

    const VkDevice& device = vulkanState.window.device;
    VulkanAppState* app = &(GetAppState(memory)->vulkanAppState);

    UnloadTextPipelineWindow(device, &app->textPipeline);
    UnloadSpritePipelineWindow(device, &app->spritePipeline);

    vkDestroyFence(device, app->fence, nullptr);
    vkDestroyCommandPool(device, app->commandPool, nullptr);
}

#include "imgui.cpp"
#include "raytracer.cpp"

#include <km_common/km_array.cpp>
#include <km_common/km_container.cpp>
#include <km_common/km_load_font.cpp>
#include <km_common/km_load_obj.cpp>
#include <km_common/km_memory.cpp>
#include <km_common/km_os.cpp>
#include <km_common/km_string.cpp>

#include <km_common/app/km_app.cpp>
#include <km_common/app/km_input.cpp>
#include <km_common/app/km_log.cpp>

#include <km_common/vulkan/km_vulkan_core.cpp>
#include <km_common/vulkan/km_vulkan_sprite.cpp>
#include <km_common/vulkan/km_vulkan_text.cpp>
#include <km_common/vulkan/km_vulkan_util.cpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#undef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#undef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_SPRINTF_IMPLEMENTATION
#include <stb_sprintf.h>
#undef STB_SPRINTF_IMPLEMENTATION
