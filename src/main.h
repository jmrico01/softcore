#pragma once

#include <km_common/app/km_log.h>

#include <km_common/km_array.h>
#include <km_common/km_debug.h>
#include <km_common/vulkan/km_vulkan_core.h>
#include <km_common/vulkan/km_vulkan_sprite.h>
#include <km_common/vulkan/km_vulkan_text.h>

#include "imgui.h"
#include "level.h"
#include "mesh.h"

enum class SpriteId
{
    PIXEL,
    JON,
    ROCK,

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

    VulkanSpritePipeline<(uint32)SpriteId::COUNT> spritePipeline;
    VulkanTextPipeline<(uint32)FontId::COUNT> textPipeline;

    VulkanMeshPipeline meshPipeline;
    VulkanLightmapMeshPipeline lightmapMeshPipeline;
};

struct AppState
{
    static const uint32 MAX_MOBS = 1024;

    VulkanAppState vulkanAppState;
    VulkanFontFace fontFaces[FontId::COUNT];

    float32 elapsedTime;

    Vec3 cameraPos;
    Vec2 cameraAngles;

    LevelData levelData;

    // Debug
    bool debugView;

    bool noclip;
    Vec3 noclipPos;

    bool blockEditor;
    PanelSliderState sliderBlockSize;
    PanelDropdownState loadLevelDropdownState;
};

struct FrameState
{
    VulkanSpriteRenderState<(uint32)SpriteId::COUNT> spriteRenderState;
    VulkanTextRenderState<(uint32)FontId::COUNT> textRenderState;

    VulkanMeshRenderState meshRenderState;
};

struct TransientState
{
    FrameState frameState;
    LargeArray<uint8> scratch;
};
