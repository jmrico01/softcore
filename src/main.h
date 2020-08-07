#pragma once

// #define TRACY_ENABLE
#include <TracyClient.cpp>

#ifdef TRACY_ENABLE
#pragma message("PROFILING ENABLED")
#endif

#include <km_common/app/km_log.h>

#include <km_common/km_array.h>
#include <km_common/km_debug.h>
#include <km_common/vulkan/km_vulkan_core.h>
#include <km_common/vulkan/km_vulkan_sprite.h>
#include <km_common/vulkan/km_vulkan_text.h>

#include "imgui.h"
#include "raytracer.h"
#include "render.h"

enum class SpriteId
{
    PIXEL,

    COUNT
};

enum class FontId
{
    OCR_A_REGULAR_18,
    OCR_A_REGULAR_24,

    COUNT
};

struct VulkanAppState
{
    VkCommandPool commandPool;
    VkCommandBuffer commandBuffer;
    VkFence fence;

    VkImage image;
    VkDeviceMemory imageMemory;

    VulkanMeshPipeline meshPipeline;
    VulkanCompositePipeline compositePipeline;

    VulkanSpritePipeline<(uint32)SpriteId::COUNT> spritePipeline;
    VulkanTextPipeline<(uint32)FontId::COUNT> textPipeline;
};

struct AppState
{
    LargeArray<uint8> arena;
    LinearAllocator arenaAllocator;

    VulkanAppState vulkanAppState;
    VulkanFontFace fontFaces[FontId::COUNT];

    float32 elapsedTime;

    Vec3 cameraPos;
    Vec2 cameraAngles;

    CanvasState canvas;
    RaycastGeometry raycastGeometry;

    // Debug
    bool debugView;
    PanelSliderState inputScreenFillState;
    PanelInputIntState inputDecayFramesState;
    PanelInputIntState inputBouncesState;
    PanelDropdownState inputSceneDropdownState;
};

struct FrameState
{
    // VulkanMeshRenderState meshRenderState;
    VulkanSpriteRenderState<(uint32)SpriteId::COUNT> spriteRenderState;
    VulkanTextRenderState<(uint32)FontId::COUNT> textRenderState;
};

struct TransientState
{
    FrameState frameState;
    LargeArray<uint8> scratch;
};
