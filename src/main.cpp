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

#define DISCRETE_GPU 1

struct ComputeTriangle
{
    alignas(16) Vec3 a;
    alignas(16) Vec3 b;
    alignas(16) Vec3 c;
    alignas(16) Vec3 normal;
    uint32 materialIndex;
};

struct ComputeMaterial
{
    Vec3 albedo;
    float32 smoothness;
    Vec3 emissionColor;
    float32 emission;
};

struct ComputeUbo
{
    const static uint32 MAX_MATERIALS = 32;

	alignas(16) Vec3 cameraPos;
	alignas(16) Vec3 filmTopLeft;
	alignas(16) Vec3 filmUnitOffsetX;
	alignas(16) Vec3 filmUnitOffsetY;
	alignas(16) ComputeMaterial materials[MAX_MATERIALS];
};

struct StartSceneInfo
{
    const_string scene;
    Vec3 pos;
    Vec2 angles;
};

const StartSceneInfo START_SCENE_INFOS[] = {
    { // front view
        .scene = ToString("interior-1"),
        .pos = Vec3 { 0.64f, -2.35f, 0.44f },
        .angles = Vec2 { -2.04f, -0.11f },
    },
    { // buggy center-box view
        .scene = ToString("interior-1"),
        .pos = Vec3 { 0.65f, 0.87f, 0.66f },
        .angles = Vec2 { 2.25f, -0.54f },
    },
    {
        .scene = ToString("interior-2"),
        .pos = Vec3 { -0.73f, -4.79f, 0.45f },
        .angles = Vec2 { 4.88f, 0.20f },
    },
    {
        .scene = ToString("light-cones"),
        .pos = Vec3 { 9.06f, -0.54f, 3.25f },
        .angles = Vec2 { -3.26f, -0.29f },
    },
    {
        .scene = ToString("simple"),
        .pos = Vec3 { 1.95f, -0.73f, 2.00f },
        .angles = Vec2 { 3.27f, -0.29f },
    },
};

const StartSceneInfo START_SCENE_INFO = START_SCENE_INFOS[0];

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

internal bool LoadSceneObj(const_string scene, LinearAllocator* allocator, LoadObjResult* obj)
{
    const_string scenePath = AllocPrintf(allocator, "data/models/%.*s.obj", scene.size, scene.data);
    return LoadObj(scenePath, Vec3::zero, 1.0f, obj, allocator);
}

internal bool LoadScene(const_string scene, AppState* appState, LinearAllocator* allocator,
                        const VulkanWindow& window, const VulkanSwapchain& swapchain, VkCommandPool commandPool)
{
    LoadObjResult obj;
    if (!LoadSceneObj(scene, allocator, &obj)) {
        LOG_ERROR("Failed to load obj file for scene %.*s\n", scene.size, scene.data);
        return false;
    }

    appState->arenaAllocator.Clear();

    RaycastGeometry geometry = CreateRaycastGeometry(obj, BOX_MAX_TRIANGLES, &appState->arenaAllocator, allocator);
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

    VulkanMeshPipeline* meshPipeline = &appState->vulkanAppState.meshPipeline;
    VulkanCompositePipeline* compositePipeline = &appState->vulkanAppState.compositePipeline;

    UnloadCompositePipelineSwapchain(window.device, compositePipeline);
    UnloadMeshPipelineSwapchain(window.device, meshPipeline);
    UnloadCompositePipelineWindow(window.device, compositePipeline);
    UnloadMeshPipelineWindow(window.device, meshPipeline);

    if (!LoadMeshPipelineWindow(window, commandPool, obj, allocator, meshPipeline)) {
        LOG_ERROR("Failed to reload window-dependent Vulkan mesh pipeline\n");
        return false;
    }
    if (!LoadCompositePipelineWindow(window, commandPool, obj, allocator, compositePipeline)) {
        LOG_ERROR("Failed to reload window-dependent Vulkan mesh pipeline\n");
        return false;
    }
    if (!LoadMeshPipelineSwapchain(window, swapchain, commandPool, allocator, meshPipeline)) {
        LOG_ERROR("Failed to reload window-dependent Vulkan mesh pipeline\n");
        return false;
    }

    // const VkImageView rasterizedImageView = meshPipeline->colorImage.view;
    const VkImageView rasterizedImageView = appState->vulkanAppState.computeImage.view;
    if (!LoadCompositePipelineSwapchain(window, swapchain, rasterizedImageView, allocator, compositePipeline)) {
        LOG_ERROR("Failed to reload window-dependent Vulkan mesh pipeline\n");
        return false;
    }

    return true;
}

APP_UPDATE_AND_RENDER_FUNCTION(AppUpdateAndRender)
{
    UNREFERENCED_PARAMETER(audio);
    ZoneScoped;

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

        appState->cameraPos = START_SCENE_INFO.pos;
        appState->cameraAngles = START_SCENE_INFO.angles;

        appState->canvas.screenFill = 0.5f;
        appState->canvas.decayFrames = 2;
        appState->canvas.bounces = 4;

        appState->canvas.test1 = 0.0f;
        appState->canvas.test2 = 0.0f;
        appState->canvas.test3 = 0.0f;

        MemSet(appState->canvas.colorHdr.data, 0, CanvasState::MAX_PIXELS * sizeof(Vec3));
        MemSet(appState->canvas.decay.data, 0, CanvasState::MAX_PIXELS * sizeof(uint8));

        {
            LinearAllocator allocator(transientState->scratch);

            LoadScene(START_SCENE_INFO.scene, appState, &allocator, vulkanState.window, vulkanState.swapchain,
                      appState->vulkanAppState.commandPool);
        }

        // Debug views 
        appState->debugView = false;
        LockCursor(false);

        appState->inputScreenFillState.value = appState->canvas.screenFill;
        appState->inputDecayFramesState.Initialize(appState->canvas.decayFrames);
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
    if (KeyDown(input, KM_KEY_E)) {
        velocity += cameraUp;
    }
    if (KeyDown(input, KM_KEY_Q)) {
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

    // Transforms world-view camera (+X forward, +Z up) to Vulkan camera (+Z forward, -Y up)
    const Quat baseCameraRot = QuatFromAngleUnitAxis(-PI_F / 2.0f, Vec3::unitY)
        * QuatFromAngleUnitAxis(PI_F / 2.0f, Vec3::unitX);
    const Mat4 baseCameraRotMat4 = UnitQuatToMat4(baseCameraRot);

    const Mat4 view = baseCameraRotMat4 * cameraRotMat4 * Translate(-appState->cameraPos);

    const float32 fov = PI_F / 4.0f;
    const float32 aspect = (float32)screenSize.x / (float32)screenSize.y;
    const float32 nearZ = 0.1f;
    const float32 farZ = 100.0f;
    const Mat4 proj = Perspective(fov, aspect, nearZ, farZ);

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
        if (panelDebugInfo.SliderFloat(&appState->inputScreenFillState, 0.0f, 1.0f)) {
            appState->canvas.screenFill = appState->inputScreenFillState.value;
        }

        panelDebugInfo.Text(ToString("Decay frames (0 - 255)"));
        if (panelDebugInfo.InputInt(&appState->inputDecayFramesState)) {
            const int value = appState->inputDecayFramesState.value;
            if (value > 0 && value <= 255) {
                appState->canvas.decayFrames = (uint8)appState->inputDecayFramesState.value;
            }
        }

        panelDebugInfo.Text(ToString("Bounces (at least 1)"));
        if (panelDebugInfo.InputInt(&appState->inputBouncesState)) {
            const int value = appState->inputBouncesState.value;
            if (value > 0) {
                appState->canvas.bounces = appState->inputBouncesState.value;
                LOG_INFO("Set bounces: %d\n", appState->canvas.bounces);
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
                if (LoadScene(scene, appState, &allocator, vulkanState.window, vulkanState.swapchain,
                              appState->vulkanAppState.commandPool)) {
                    LOG_INFO("Loaded scene %.*s\n", scene.size, scene.data);
                }
                else {
                    LOG_ERROR("Failed to load scene %.*s\n", scene.size, scene.data);
                }
            }
        }

        panelDebugInfo.Text(string::empty);
        static PanelSliderState inputTest1 = {
            .value = appState->canvas.test1
        };
        static PanelSliderState inputTest2 = {
            .value = appState->canvas.test2
        };
        static PanelSliderState inputTest3 = {
            .value = appState->canvas.test3
        };
        if (panelDebugInfo.SliderFloat(&inputTest1, 0.0f, 1.0f)) {
            appState->canvas.test1 = inputTest1.value;
        }
        if (panelDebugInfo.SliderFloat(&inputTest2, -0.2f, 0.0f)) {
            appState->canvas.test2 = inputTest2.value;
        }
        if (panelDebugInfo.SliderFloat(&inputTest3, 0.0f, 0.2f)) {
            appState->canvas.test3 = inputTest3.value;
        }

        panelDebugInfo.Draw(panelBorderSize, Vec4::one, backgroundColor, screenSize,
                            &transientState->frameState.spriteRenderState, &transientState->frameState.textRenderState);
    }

    // ================================================================================================
    // Rendering ======================================================================================
    // ================================================================================================

    {
        ZoneScopedN("PreRasterize");

        const VulkanMeshPipeline& meshPipeline = appState->vulkanAppState.meshPipeline;

        const MeshUniformBufferObject ubo = { .view = view, .proj = proj };
        void* data;
        vkMapMemory(vulkanState.window.device, meshPipeline.uniformBuffer.memory, 0,
                    sizeof(MeshUniformBufferObject), 0, &data);
        MemCopy(data, &ubo, sizeof(ubo));
        vkUnmapMemory(vulkanState.window.device, meshPipeline.uniformBuffer.memory);

        VkSubmitInfo submitInfo = {};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submitInfo.waitSemaphoreCount = 0;
        submitInfo.pWaitSemaphores = nullptr;
        submitInfo.pWaitDstStageMask = nullptr;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &meshPipeline.commandBuffer;
        submitInfo.signalSemaphoreCount = 0;
        submitInfo.pSignalSemaphores = nullptr;

        if (vkQueueSubmit(vulkanState.window.graphicsQueue, 1, &submitInfo, meshPipeline.fence) != VK_SUCCESS) {
            LOG_ERROR("Failed to submit mesh command buffer\n");
        }

        // TODO right now, CPU takes long enough to raytrace that there's no chance of this work
        // overlapping with the compositing pass. not guaranteed though, make sure to have some semaphores at least.
        // If we need the material indices, also need the fence below
        {
            ZoneScopedN("Wait");

            if (vkWaitForFences(vulkanState.window.device, 1, &meshPipeline.fence, VK_TRUE, UINT64_MAX) != VK_SUCCESS) {
                LOG_ERROR("vkWaitForFences didn't return success for fence %lu\n", swapchainImageIndex);
            }
            if (vkResetFences(vulkanState.window.device, 1, &meshPipeline.fence) != VK_SUCCESS) {
                LOG_ERROR("vkResetFences didn't return success for fence %lu\n", swapchainImageIndex);
            }
        }
    }

#if 0
    // Ray traced rendering
    {
        ZoneScopedN("RaytracedRendering");

        LinearAllocator allocator(transientState->scratch);

        uint8* materialIndices = allocator.New<uint8>(screenSize.x * screenSize.y);
#if !DISCRETE_GPU
        {
            ZoneScopedN("ReadMaterialIndices");

            const uint32 materialIndicesSize = screenSize.x * screenSize.y * sizeof(uint8);
            const VkDeviceMemory& materialIndexMemory = appState->vulkanAppState.meshPipeline.materialIndexImageDstMemory;

            void* data;
            vkMapMemory(vulkanState.window.device, materialIndexMemory, 0, materialIndicesSize, 0, &data);
            MemCopy(materialIndices, data, materialIndicesSize);
            vkUnmapMemory(vulkanState.window.device, materialIndexMemory);
        }
#endif

        const uint32 numBytes = screenSize.x * screenSize.y * 4;
        void* imageMemory;
        vkMapMemory(vulkanState.window.device, appState->vulkanAppState.imageMemory, 0, numBytes, 0, &imageMemory);
        uint32* pixels = (uint32*)imageMemory;

        RaytraceRender(appState->cameraPos, cameraRot, fov, appState->raycastGeometry, materialIndices,
                       screenSize.x, screenSize.y, &appState->canvas, pixels, &allocator, queue);

        vkUnmapMemory(vulkanState.window.device, appState->vulkanAppState.imageMemory);
    }
#else
    UNREFERENCED_PARAMETER(queue);
#endif

    // Compute pipeline raytracing
    {
        const VulkanAppState& vulkanAppState = appState->vulkanAppState;

        LinearAllocator allocator(transientState->scratch);

        const RaycastGeometry& geometry = appState->raycastGeometry;

        ComputeUbo* computeUbo = allocator.New<ComputeUbo>();
        {
            const uint32 width = screenSize.x;
            const uint32 height = screenSize.y;
            const Vec3 cameraPos = appState->cameraPos;

            // TODO copy-pasted from raytracer.cpp for now
            const Quat inverseCameraRot = Inverse(cameraRot);
            const Vec3 cameraUp2 = inverseCameraRot * Vec3::unitZ;
            const Vec3 cameraForward2 = inverseCameraRot * Vec3::unitX;
            const Vec3 cameraLeft2 = inverseCameraRot * Vec3::unitY;

            const float32 filmDist = 1.0f;
            const float32 filmHeight = tanf(fov / 2.0f) * 2.0f;
            const float32 filmWidth = filmHeight * (float32)width / (float32)height;

            computeUbo->cameraPos = cameraPos;
            computeUbo->filmTopLeft = cameraPos + cameraForward2 * filmDist
                + (filmWidth / 2.0f + -0.136f) * cameraUp2 + (filmHeight / 2.0f + 0.136f) * cameraLeft2;
            computeUbo->filmUnitOffsetX = -cameraLeft2 * filmWidth / (float32)(width - 1);
            computeUbo->filmUnitOffsetY = -cameraUp2 * filmHeight / (float32)(height - 1);
        }

        uint32 numMaterials = 0;
        for (uint32 i = 0; i < geometry.materials.size; i++) {
            const RaycastMaterial& srcMaterial = geometry.materials[i];
            ComputeMaterial* dstMaterial = &computeUbo->materials[numMaterials++];

            dstMaterial->albedo = srcMaterial.albedo;
            dstMaterial->smoothness = srcMaterial.smoothness;
            dstMaterial->emissionColor = srcMaterial.emissionColor;
            dstMaterial->emission = srcMaterial.emission;

            if (numMaterials >= ComputeUbo::MAX_MATERIALS) {
                LOG_ERROR("Too many materials, hit limit of %lu\n", ComputeUbo::MAX_MATERIALS);
                break;
            }
        }

        void* data;
        vkMapMemory(vulkanState.window.device, vulkanAppState.computeUniform.memory, 0, sizeof(ComputeUbo), 0, &data);
        MemCopy(data, computeUbo, sizeof(ComputeUbo));
        vkUnmapMemory(vulkanState.window.device, vulkanAppState.computeUniform.memory);

		VkSubmitInfo submitInfo = {};
        submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
        submitInfo.waitSemaphoreCount = 0;
        submitInfo.pWaitSemaphores = nullptr;
        submitInfo.pWaitDstStageMask = nullptr;
		submitInfo.commandBufferCount = 1;
		submitInfo.pCommandBuffers = &vulkanAppState.computeCommandBuffer;
        submitInfo.signalSemaphoreCount = 0;
        submitInfo.pSignalSemaphores = nullptr;

		if (vkQueueSubmit(vulkanAppState.computeQueue, 1, &submitInfo, vulkanAppState.computeFence) != VK_SUCCESS) {
            LOG_ERROR("vkQueueSubmit failed\n");
        }

        if (vkWaitForFences(vulkanState.window.device, 1, &vulkanAppState.computeFence,
                            VK_TRUE, UINT64_MAX) != VK_SUCCESS) {
            LOG_ERROR("vkWaitForFences failed\n");
        }
        if (vkResetFences(vulkanState.window.device, 1, &vulkanAppState.computeFence) != VK_SUCCESS) {
            LOG_ERROR("vkResetFences failed\n");
        }
    }

    {
        LinearAllocator allocator(transientState->scratch);

        const VulkanCompositePipeline& compositePipeline = appState->vulkanAppState.compositePipeline;
        const RaycastGeometry& geometry = appState->raycastGeometry;

        CompositeUniformBufferObject* compositeUbo = allocator.New<CompositeUniformBufferObject>();
        {
            const uint32 width = screenSize.x;
            const uint32 height = screenSize.y;
            const Vec3 cameraPos = appState->cameraPos;

            // TODO copy-pasted from raytracer.cpp for now
            const Quat inverseCameraRot = Inverse(cameraRot);
            const Vec3 cameraUp2 = inverseCameraRot * Vec3::unitZ;
            const Vec3 cameraForward2 = inverseCameraRot * Vec3::unitX;
            const Vec3 cameraLeft2 = inverseCameraRot * Vec3::unitY;

            const float32 filmDist = 1.0f;
            const float32 filmHeight = tanf(fov / 2.0f) * 2.0f;
            const float32 filmWidth = filmHeight * (float32)width / (float32)height;

            compositeUbo->filmTopLeft = cameraPos + cameraForward2 * filmDist
                + (filmWidth / 2.0f + -0.136f) * cameraUp2 + (filmHeight / 2.0f + 0.136f) * cameraLeft2;
            compositeUbo->filmUnitOffsetX = -cameraLeft2 * filmWidth / (float32)(width - 1);
            compositeUbo->filmUnitOffsetY = -cameraUp2 * filmHeight / (float32)(height - 1);
            compositeUbo->cameraPos = cameraPos;
        }

        uint32 numMaterials = 0;
        for (uint32 i = 0; i < geometry.materials.size; i++) {
            const RaycastMaterial& srcMaterial = geometry.materials[i];
            CompositeMaterial* dstMaterial = &compositeUbo->materials[numMaterials++];

            dstMaterial->albedo = srcMaterial.albedo;
            dstMaterial->smoothness = srcMaterial.smoothness;
            dstMaterial->emissionColor = srcMaterial.emissionColor;
            dstMaterial->emission = srcMaterial.emission;

            if (numMaterials >= CompositeUniformBufferObject::MAX_MATERIALS) {
                break;
            }
        }

        // TODO this is only for counting triangles now. should do this when creating the pipeline
        // number of triangles could easily be made a shader constant
        DynamicArray<const RaycastMeshBvh*, LinearAllocator> bvhStack(&allocator);
        compositeUbo->numTriangles = 0;
        for (uint32 i = 0; i < geometry.meshes.size; i++) {
            bvhStack.Clear();
            bvhStack.Append(&geometry.meshes[i].bvh);
            while (bvhStack.size > 0) {
                const RaycastMeshBvh* bvh = bvhStack[bvhStack.size - 1];
                bvhStack.RemoveLast();

                if (bvh->child1 == nullptr) {
                    compositeUbo->numTriangles += bvh->triangles.size;
                }
                else {
                    bvhStack.Append(bvh->child1);
                    bvhStack.Append(bvh->child2);
                }
            }
        }

        void* data;
        vkMapMemory(vulkanState.window.device, compositePipeline.uniformBuffer.memory, 0,
                    sizeof(CompositeUniformBufferObject), 0, &data);
        MemCopy(data, compositeUbo, sizeof(CompositeUniformBufferObject));
        vkUnmapMemory(vulkanState.window.device, compositePipeline.uniformBuffer.memory);
    }

    // Vulkan rendering
    VkCommandBuffer buffer = appState->vulkanAppState.commandBuffer;
    VkFence fence = appState->vulkanAppState.fence;

    {
        ZoneScopedN("WaitForFinalRenderFence");

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

    TransitionImageLayout(buffer, appState->vulkanAppState.computeImage.image,
                          VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL);

    TransitionImageLayout(buffer, appState->vulkanAppState.compositePipeline.raytracedImage.image,
                          VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL);
    vkCmdCopyImage(buffer,
                   appState->vulkanAppState.image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                   appState->vulkanAppState.compositePipeline.raytracedImage.image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                   1, &imageCopy);
    TransitionImageLayout(buffer, appState->vulkanAppState.compositePipeline.raytracedImage.image,
                          VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL);

    const VkClearValue clearValues[] = {
        { 0.0f, 0.0f, 0.0f, 0.0f },
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

    // Composite
    {
        const VulkanCompositePipeline& compositePipeline = appState->vulkanAppState.compositePipeline;

        vkCmdBindPipeline(buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, compositePipeline.pipeline);

        const VkBuffer vertexBuffers[] = { compositePipeline.vertexBuffer.buffer };
        const VkDeviceSize offsets[] = { 0 };
        vkCmdBindVertexBuffers(buffer, 0, C_ARRAY_LENGTH(vertexBuffers), vertexBuffers, offsets);

        vkCmdBindDescriptorSets(buffer, VK_PIPELINE_BIND_POINT_GRAPHICS, compositePipeline.pipelineLayout, 0,
                                1, &compositePipeline.descriptorSet, 0, nullptr);

        vkCmdDraw(buffer, 6, 1, 0, 0);
    }

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

    // TODO uhhhh
    TransitionImageLayout(buffer, appState->vulkanAppState.computeImage.image,
                          VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);

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

    // Create images
    {
        if (!CreateImage(window.device, window.physicalDevice, swapchain.extent.width, swapchain.extent.height,
                         VK_FORMAT_R8G8B8A8_UINT, VK_IMAGE_TILING_LINEAR, VK_IMAGE_USAGE_TRANSFER_SRC_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                         &app->image, &app->imageMemory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, app->commandPool, window.graphicsQueue);
        TransitionImageLayout(commandBuffer, app->image,
                              VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL);

#if 0
        if (!CreateImage(window.device, window.physicalDevice, swapchain.extent.width, swapchain.extent.height,
                         VK_FORMAT_R8G8B8A8_UINT, VK_IMAGE_TILING_OPTIMAL,
                         VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT,
                         VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                         &app->imageSampled.image, &app->imageSampled.memory)) {
            LOG_ERROR("CreateImage failed\n");
            return false;
        }

        if (!CreateImageView(window.device, app->imageSampled.image, VK_FORMAT_R8G8B8A8_UINT, VK_IMAGE_ASPECT_COLOR_BIT,
                             &app->imageSampled.view)) {
            LOG_ERROR("CreateImageView failed\n");
            return false;
        }
#endif
    }

    if (!LoadMeshPipelineSwapchain(window, swapchain, app->commandPool, &allocator, &app->meshPipeline)) {
        LOG_ERROR("Failed to load swapchain-dependent Vulkan mesh pipeline\n");
        return false;
    }

    // const VkImageView rasterizedImageView = app->meshPipeline.colorImage.view;
    const VkImageView rasterizedImageView = app->computeImage.view;
    if (!LoadCompositePipelineSwapchain(window, swapchain, rasterizedImageView, &allocator, &app->compositePipeline)) {
        LOG_ERROR("Failed to load swapchain-dependent Vulkan composite pipeline\n");
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
    UnloadCompositePipelineSwapchain(device, &app->compositePipeline);
    UnloadMeshPipelineSwapchain(device, &app->meshPipeline);

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

    LoadObjResult obj;
    if (!LoadSceneObj(START_SCENE_INFO.scene, &allocator, &obj)) {
        LOG_ERROR("Failed to load obj file for scene %.*s\n", START_SCENE_INFO.scene.size, START_SCENE_INFO.scene.data);
        return false;
    }

    RaycastGeometry geometry = CreateRaycastGeometry(obj, UINT32_MAX_VALUE, &allocator, &allocator);
    if (geometry.meshes.data == nullptr) {
        LOG_ERROR("Failed to load geometry from obj\n");
        return false;
    }

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

    // Create compute resources
    {
        const uint32 imageWidth = WINDOW_START_WIDTH;
        const uint32 imageHeight = WINDOW_START_HEIGHT;

        // compute queue and command pool
        {
            QueueFamilyInfo qfi = GetQueueFamilyInfo(window.surface, window.physicalDevice, &allocator);
            if (!qfi.hasComputeFamily) {
                LOG_ERROR("Device has no compute family\n");
                return false;
            }

            VkDeviceQueueCreateInfo queueCreateInfo = {};
            queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
            queueCreateInfo.pNext = nullptr;
            queueCreateInfo.queueFamilyIndex = qfi.computeFamilyIndex;
            queueCreateInfo.queueCount = 1;
            vkGetDeviceQueue(window.device, qfi.computeFamilyIndex, 0, &app->computeQueue);

            VkCommandPoolCreateInfo commandPoolCreateInfo = {};
            commandPoolCreateInfo.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
            commandPoolCreateInfo.queueFamilyIndex = qfi.computeFamilyIndex;
            commandPoolCreateInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

            if (vkCreateCommandPool(window.device, &commandPoolCreateInfo, nullptr,
                                    &app->computeCommandPool) != VK_SUCCESS) {
                LOG_ERROR("vkCreateCommandPool failed\n");
                return false;
            }
        }

        // triangles
        {
            DynamicArray<ComputeTriangle, LinearAllocator> triangles(&allocator);
            for (uint32 i = 0; i < geometry.meshes.size; i++) {
                const RaycastMesh& mesh = geometry.meshes[i];
                if (mesh.numTriangles != mesh.bvh.triangles.size) {
                    LOG_ERROR("BVH structure, expected flat\n");
                    return false;
                }

                for (uint32 j = 0; j < mesh.bvh.triangles.size; j++) {
                    const RaycastTriangle& srcTriangle = mesh.bvh.triangles[j];
                    ComputeTriangle* dstTriangle = triangles.Append();
                    dstTriangle->a = srcTriangle.pos[0];
                    dstTriangle->b = srcTriangle.pos[1];
                    dstTriangle->c = srcTriangle.pos[2];
                    dstTriangle->normal = srcTriangle.normal;
                    dstTriangle->materialIndex = srcTriangle.materialIndex;
                }
            }

            const VkDeviceSize bufferSize = triangles.size * sizeof(ComputeTriangle);
            app->numTriangles = triangles.size;

            VulkanBuffer stagingBuffer;
            if (!CreateVulkanBuffer(bufferSize,
                                    VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                    window.device, window.physicalDevice, &stagingBuffer)) {
                LOG_ERROR("CreateBuffer failed for staging buffer\n");
                return false;
            }

            // Copy vertex data from CPU into memory-mapped staging buffer
            void* data;
            vkMapMemory(window.device, stagingBuffer.memory, 0, bufferSize, 0, &data);
            MemCopy(data, triangles.data, bufferSize);
            vkUnmapMemory(window.device, stagingBuffer.memory);

            if (!CreateVulkanBuffer(bufferSize,
                                    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                    VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                                    window.device, window.physicalDevice, &app->computeTriangles)) {
                LOG_ERROR("CreateBuffer failed for vertex buffer\n");
                return false;
            }

            // Copy vertex data from staging buffer into GPU vertex buffer
            {
                SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, app->commandPool, window.graphicsQueue);
                CopyBuffer(commandBuffer, stagingBuffer.buffer, app->computeTriangles.buffer, bufferSize);
            }

            DestroyVulkanBuffer(window.device, &stagingBuffer);
        }

        // uniform buffer
        {
            const VkDeviceSize uniformBufferSize = sizeof(ComputeUbo);
            if (!CreateVulkanBuffer(uniformBufferSize,
                                    VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                    window.device, window.physicalDevice, &app->computeUniform)) {
                LOG_ERROR("CreateBuffer failed for vertex buffer\n");
                return false;
            }
        }

        // image
        {
            const VkFormat format = VK_FORMAT_R8G8B8A8_UNORM;
            VkFormatProperties formatProperties;
            vkGetPhysicalDeviceFormatProperties(window.physicalDevice, format, &formatProperties);
            if ((formatProperties.optimalTilingFeatures & VK_FORMAT_FEATURE_STORAGE_IMAGE_BIT) == 0) {
                LOG_ERROR("Storage image unsupported for format %d\n", format);
                return false;
            }

            if (!CreateImage(window.device, window.physicalDevice, imageWidth, imageHeight, format,
                             VK_IMAGE_TILING_OPTIMAL,
                             VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT,
                             VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                             &app->computeImage.image, &app->computeImage.memory)) {
                LOG_ERROR("CreateImage failed\n");
                return false;
            }

            if (!CreateImageView(window.device, app->computeImage.image, format, VK_IMAGE_ASPECT_COLOR_BIT,
                                 &app->computeImage.view)) {
                LOG_ERROR("CreateImageView failed\n");
                return false;
            }

            SCOPED_VK_COMMAND_BUFFER(commandBuffer, window.device, app->commandPool, window.graphicsQueue);
            TransitionImageLayout(commandBuffer, app->computeImage.image,
                                  VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
        }

        // sampler
        {
            VkSamplerCreateInfo createInfo = {};
            createInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
            createInfo.magFilter = VK_FILTER_NEAREST;
            createInfo.minFilter = VK_FILTER_NEAREST;
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

            if (vkCreateSampler(window.device, &createInfo, nullptr, &app->computeSampler) != VK_SUCCESS) {
                LOG_ERROR("vkCreateSampler failed\n");
                return false;
            }
        }

        // descriptor set layout
        {
            VkDescriptorSetLayoutBinding layoutBindings[3] = {};
            layoutBindings[0].binding = 0;
            layoutBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            layoutBindings[0].descriptorCount = 1;
            layoutBindings[0].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
            layoutBindings[0].pImmutableSamplers = nullptr;

            layoutBindings[1].binding = 1;
            layoutBindings[1].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
            layoutBindings[1].descriptorCount = 1;
            layoutBindings[1].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
            layoutBindings[1].pImmutableSamplers = nullptr;

            layoutBindings[2].binding = 2;
            layoutBindings[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            layoutBindings[2].descriptorCount = 1;
            layoutBindings[2].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
            layoutBindings[2].pImmutableSamplers = nullptr;

            VkDescriptorSetLayoutCreateInfo layoutCreateInfo = {};
            layoutCreateInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
            layoutCreateInfo.bindingCount = C_ARRAY_LENGTH(layoutBindings);
            layoutCreateInfo.pBindings = layoutBindings;

            if (vkCreateDescriptorSetLayout(window.device, &layoutCreateInfo, nullptr,
                                            &app->computeDescriptorSetLayout) != VK_SUCCESS) {
                LOG_ERROR("vkCreateDescriptorSetLayout failed\n");
                return false;
            }
        }

        // descriptor pool
        {
            VkDescriptorPoolSize poolSizes[3] = {};

            poolSizes[0].type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            poolSizes[0].descriptorCount = 1;

            poolSizes[1].type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
            poolSizes[1].descriptorCount = 1;

            poolSizes[2].type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            poolSizes[2].descriptorCount = 1;

            VkDescriptorPoolCreateInfo poolInfo = {};
            poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
            poolInfo.poolSizeCount = C_ARRAY_LENGTH(poolSizes);
            poolInfo.pPoolSizes = poolSizes;
            poolInfo.maxSets = 1;

            if (vkCreateDescriptorPool(window.device, &poolInfo, nullptr, &app->computeDescriptorPool) != VK_SUCCESS) {
                LOG_ERROR("vkCreateDescriptorPool failed\n");
                return false;
            }
        }

        // descriptor set
        {
            VkDescriptorSetAllocateInfo allocInfo = {};
            allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
            allocInfo.descriptorPool = app->computeDescriptorPool;
            allocInfo.descriptorSetCount = 1;
            allocInfo.pSetLayouts = &app->computeDescriptorSetLayout;

            if (vkAllocateDescriptorSets(window.device, &allocInfo, &app->computeDescriptorSet) != VK_SUCCESS) {
                LOG_ERROR("vkAllocateDescriptorSets failed\n");
                return false;
            }

            VkWriteDescriptorSet descriptorWrites[3] = {};

            const VkDescriptorImageInfo imageInfo = {
                .sampler = VK_NULL_HANDLE,
                .imageView = app->computeImage.view,
                .imageLayout = VK_IMAGE_LAYOUT_GENERAL
            };
            descriptorWrites[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrites[0].dstSet = app->computeDescriptorSet;
            descriptorWrites[0].dstBinding = 0;
            descriptorWrites[0].dstArrayElement = 0;
            descriptorWrites[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            descriptorWrites[0].descriptorCount = 1;
            descriptorWrites[0].pImageInfo = &imageInfo;

            const VkDescriptorBufferInfo uniformBufferInfo = {
                .buffer = app->computeUniform.buffer,
                .offset = 0,
                .range = sizeof(ComputeUbo),
            };
            descriptorWrites[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrites[1].dstSet = app->computeDescriptorSet;
            descriptorWrites[1].dstBinding = 1;
            descriptorWrites[1].dstArrayElement = 0;
            descriptorWrites[1].descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER;
            descriptorWrites[1].descriptorCount = 1;
            descriptorWrites[1].pBufferInfo = &uniformBufferInfo;

            const VkDescriptorBufferInfo triangleBufferInfo = {
                .buffer = app->computeTriangles.buffer,
                .offset = 0,
                .range = app->numTriangles * sizeof(ComputeTriangle),
            };
            descriptorWrites[2].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            descriptorWrites[2].dstSet = app->computeDescriptorSet;
            descriptorWrites[2].dstBinding = 2;
            descriptorWrites[2].dstArrayElement = 0;
            descriptorWrites[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            descriptorWrites[2].descriptorCount = 1;
            descriptorWrites[2].pBufferInfo = &triangleBufferInfo;

            vkUpdateDescriptorSets(window.device, C_ARRAY_LENGTH(descriptorWrites), descriptorWrites, 0, nullptr);
        }

        // pipeline
        {
            const Array<uint8> shaderCode = LoadEntireFile(ToString("data/shaders/raytrace.comp.spv"), &allocator);
            if (shaderCode.data == nullptr) {
                LOG_ERROR("Failed to load compute shader code\n");
                return false;
            }

            VkShaderModule shaderModule;
            if (!CreateShaderModule(shaderCode, window.device, &shaderModule)) {
                LOG_ERROR("Failed to create compute shader module\n");
                return false;
            }
            defer(vkDestroyShaderModule(window.device, shaderModule, nullptr));

            VkPipelineShaderStageCreateInfo shaderStageCreateInfo = {};
            shaderStageCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
            shaderStageCreateInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
            shaderStageCreateInfo.module = shaderModule;
            shaderStageCreateInfo.pName = "main";

            VkPipelineLayoutCreateInfo pipelineLayoutCreateInfo = {};
            pipelineLayoutCreateInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
            pipelineLayoutCreateInfo.setLayoutCount = 1;
            pipelineLayoutCreateInfo.pSetLayouts = &app->computeDescriptorSetLayout;
            pipelineLayoutCreateInfo.pushConstantRangeCount = 0;
            pipelineLayoutCreateInfo.pPushConstantRanges = nullptr;

            if (vkCreatePipelineLayout(window.device, &pipelineLayoutCreateInfo, nullptr,
                                       &app->computePipelineLayout) != VK_SUCCESS) {
                LOG_ERROR("vkCreatePipelineLayout failed\n");
                return false;
            }

            VkComputePipelineCreateInfo pipelineCreateInfo = {};
            pipelineCreateInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
            pipelineCreateInfo.pNext = nullptr;
            pipelineCreateInfo.flags = 0;
            pipelineCreateInfo.stage = shaderStageCreateInfo;
            pipelineCreateInfo.layout = app->computePipelineLayout;

            if (vkCreateComputePipelines(window.device, VK_NULL_HANDLE, 1, &pipelineCreateInfo, nullptr,
                                         &app->computePipeline) != VK_SUCCESS) {
                LOG_ERROR("vkCreateComputePipelines failed\n");
                return false;
            }
        }

        // command buffer
        {
            VkCommandBufferAllocateInfo allocInfo = {};
            allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
            allocInfo.pNext = nullptr;
            allocInfo.commandPool = app->computeCommandPool;
            allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
            allocInfo.commandBufferCount = 1;

            if (vkAllocateCommandBuffers(window.device, &allocInfo, &app->computeCommandBuffer) != VK_SUCCESS) {
                LOG_ERROR("vkAllocateCommandBuffers failed\n");
                return false;
            }

            VkCommandBufferBeginInfo beginInfo = {};
            beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
            beginInfo.flags = 0;
            beginInfo.pInheritanceInfo = nullptr;

            if (vkBeginCommandBuffer(app->computeCommandBuffer, &beginInfo) != VK_SUCCESS) {
                LOG_ERROR("vkBeginCommandBuffer failed\n");
                return false;
            }

            vkCmdBindPipeline(app->computeCommandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, app->computePipeline);
            vkCmdBindDescriptorSets(app->computeCommandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                                    app->computePipelineLayout, 0, 1, &app->computeDescriptorSet, 0, 0);

            const uint32 batchSize = VulkanAppState::BATCH_SIZE;
            vkCmdDispatch(app->computeCommandBuffer, imageWidth / batchSize, imageHeight / batchSize, 1);

            if (vkEndCommandBuffer(app->computeCommandBuffer) != VK_SUCCESS) {
                LOG_ERROR("vkEndCommandBuffer failed\n");
                return false;
            }
        }

        // fence
        {
            VkFenceCreateInfo fenceCreateInfo = {};
            fenceCreateInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
            fenceCreateInfo.flags = 0;//VK_FENCE_CREATE_SIGNALED_BIT;

            if (vkCreateFence(window.device, &fenceCreateInfo, nullptr, &app->computeFence) != VK_SUCCESS) {
                LOG_ERROR("vkCreateFence failed\n");
                return false;
            }
        }
    }

    const bool meshPipeline = LoadMeshPipelineWindow(window, app->commandPool, obj, &allocator, &app->meshPipeline);
    if (!meshPipeline) {
        LOG_ERROR("Failed to load window-dependent Vulkan mesh pipeline\n");
        return false;
    }

    const bool compositePipeline = LoadCompositePipelineWindow(window, app->commandPool, obj, &allocator,
                                                               &app->compositePipeline);
    if (!compositePipeline) {
        LOG_ERROR("Failed to load window-dependent Vulkan composite pipeline\n");
        return false;
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
    UnloadCompositePipelineWindow(device, &app->compositePipeline);
    UnloadMeshPipelineWindow(device, &app->meshPipeline);

    vkDestroyFence(device, app->fence, nullptr);
    vkDestroyCommandPool(device, app->commandPool, nullptr);
}

#include "imgui.cpp"
#include "raytracer.cpp"
#include "render.cpp"

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
