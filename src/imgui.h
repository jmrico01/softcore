#pragma once

#include <km_common/km_array.h>
#include <km_common/km_string.h>
#include <km_common/app/km_input.h>

enum class PanelRenderCommandType
{
	RECT,
	TEXT
};

using PanelRenderCommandFlags = uint8;
namespace PanelRenderCommandFlag
{
constexpr PanelRenderCommandFlags OVERRIDE_COLOR = 1 << 0;
constexpr PanelRenderCommandFlags OVERRIDE_ALPHA = 1 << 1;
}

struct PanelRenderCommandRect
{
	Vec2Int size;
};

struct PanelRenderCommandText
{
    string text; // NOTE we don't own this data
	const VulkanFontFace* fontFace;
};

struct PanelRenderCommand
{
	PanelRenderCommandType type;
    PanelRenderCommandFlags flags;
    Vec2Int position;
    float32 depth;
    Vec2 anchor;
    Vec4 color;
	union
	{
		PanelRenderCommandRect commandRect;
		PanelRenderCommandText commandText;
	};
};

using PanelFlags = uint32;
namespace PanelFlag
{
constexpr PanelFlags GROW_DOWNWARDS = 1 << 0;
constexpr PanelFlags MINIMIZED      = 1 << 1;
}

struct PanelDropdownState
{
    uint32 selected;
    bool opened;
};

struct PanelSliderState
{
    float32 value;
    bool drag;
};

struct PanelInputTextState
{
    static const uint32 MAX_LENGTH = 256;
    bool focused;
    FixedArray<char, MAX_LENGTH> text;
};

struct PanelInputIntState
{
    bool valid;
    int value;
    PanelInputTextState textState;

    void Initialize(int value);
};

struct Panel
{
    PanelFlags flags;
    Vec2Int position;
	Vec2Int positionCurrent;
    float32 anchorX;
	Vec2Int size;
	DynamicArray<PanelRenderCommand, LinearAllocator> renderCommands;
    LinearAllocator* allocator;
	const AppInput* input;
    const VulkanFontFace* fontFaceDefault;

    Panel(LinearAllocator* allocator);

    // NOTE: position is needed on Panel::Begin to properly handle mouse input for all components
    void Begin(const AppInput& input, const VulkanFontFace* fontFace, PanelFlags flags, Vec2Int position,
               float32 anchorX);

    bool TitleBar(const_string text, bool* minimized = nullptr, Vec4 color = Vec4::zero,
                  const VulkanFontFace* fontFace = nullptr);

	void Text(const_string text, Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);

	bool Button(const_string text, Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);
	bool Checkbox(bool* value, const_string text, Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);

	bool Dropdown(PanelDropdownState* state, const Array<string>& values, Vec4 color = Vec4::zero,
                  const VulkanFontFace* fontFace = nullptr);

	bool SliderFloat(PanelSliderState* state, float32 min, float32 max, const_string text = const_string::empty,
                     Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);

    bool InputText(PanelInputTextState* state, Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);
    bool InputInt(PanelInputIntState* state, Vec4 color = Vec4::zero, const VulkanFontFace* fontFace = nullptr);

    template <uint32 S1, uint32 S2>
        void Draw(Vec2Int borderSize, Vec4 defaultColor, Vec4 backgroundColor, Vec2Int screenSize,
                  VulkanSpriteRenderState<S1>* spriteRenderState, VulkanTextRenderState<S2>* textRenderState);

    // Unimplemented ...
    int SliderInt();
	float32 InputFloat();
	Vec2 InputVec2();
	Vec3 InputVec3();
	Vec4 InputVec4();
	Vec4 InputColor();
};
