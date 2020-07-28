#include "imgui.h"

const float32 MARGIN_FRACTION = 1.2f;
// TODO fix, panel items depth-cull against regular world geometry, so these values have to be really small for now
const float32 BACKGROUND_DEPTH = 0.1f;
const float32 PANEL_ITEM_DEPTH = 0.05f;
const float32 PANEL_TEXT_DEPTH = 0.0f;

// Add text render command to panel. Returns size of text in panel
Vec2Int AddText(Panel* panel, Vec2Int offset, const_string text, Vec4 color, const VulkanFontFace* fontFace,
                PanelRenderCommand** command = nullptr)
{
    DEBUG_ASSERT(panel != nullptr);

    const VulkanFontFace* fontToUse = fontFace == nullptr ? panel->fontFaceDefault : fontFace;
    const Vec2Int size = {
        (int)((float32)GetTextWidth(*fontToUse, text) * MARGIN_FRACTION),
        (int)((float32)fontToUse->height * MARGIN_FRACTION)
    };

    PanelRenderCommand* newCommand = panel->renderCommands.Append();
    newCommand->type = PanelRenderCommandType::TEXT;
    newCommand->flags = color == Vec4::zero ? 0 : PanelRenderCommandFlag::OVERRIDE_COLOR;
    offset.x = Lerp(offset.x, -offset.x, panel->anchorX);
    newCommand->position = offset + panel->positionCurrent;
    if (panel->flags & PanelFlag::GROW_DOWNWARDS) {
        newCommand->position.y += size.y;
    }
    newCommand->color = color;
    newCommand->anchor = { panel->anchorX, 0.0f };
    newCommand->commandText.text = ToNonConstString(text);
    newCommand->commandText.fontFace = fontToUse;

    if (command != nullptr) {
        *command = newCommand;
    }
    return size;
}

// Add rect render command to panel. Returns whether the mouse is inside ("hovering") the rect
bool AddRect(Panel* panel, Vec2Int offset, Vec2Int size, Vec4 color, bool interactiveAlpha,
             PanelRenderCommand** command = nullptr)
{
    DEBUG_ASSERT(panel != nullptr);

    const float32 PRESSED_ALPHA = 1.0f;
    const float32 HOVERED_ALPHA = 0.65f;
    const float32 IDLE_ALPHA = 0.3f;

    const Vec2Int rectCenterRelative = Vec2Int {
        Lerp(size.x / 2, -size.x / 2, panel->anchorX),
        panel->flags & PanelFlag::GROW_DOWNWARDS ? size.y / 2 : -size.y / 2
    };

    const Vec2Int rectCenter = offset + panel->positionCurrent + rectCenterRelative;
    const RectInt rectRect = {
        .min = rectCenter - size / 2,
        .max = rectCenter + size / 2
    };
    const bool hovered = IsInside(panel->input->mousePos, rectRect);
    const bool pressed = hovered && MouseDown(*panel->input, KM_MOUSE_LEFT);

    PanelRenderCommand* newCommand = panel->renderCommands.Append();
    newCommand->type = PanelRenderCommandType::RECT;
    if (interactiveAlpha) {
        newCommand->flags = PanelRenderCommandFlag::OVERRIDE_ALPHA;
    }
    if (color != Vec4::zero) {
        newCommand->flags |= PanelRenderCommandFlag::OVERRIDE_COLOR;
    }
    newCommand->position = rectCenter;
    newCommand->anchor = Vec2 { 0.5f, 0.5f };
    newCommand->color = color;
    if (interactiveAlpha) {
        newCommand->color.a = pressed ? PRESSED_ALPHA : hovered ? HOVERED_ALPHA : IDLE_ALPHA;
    }
    newCommand->commandRect.size = size;

    if (command != nullptr) {
        *command = newCommand;
    }
    return hovered;
}

void PanelInputIntState::Initialize(int value)
{
    string text = { .size = textState.MAX_LENGTH, .data = textState.text.data };
    valid = SizedPrintf(&text, "%d", value);
    if (valid) {
        textState.text.size = text.size;
        this->value = value;
    }
}

Panel::Panel(LinearAllocator* allocator)
: allocator(allocator), renderCommands(allocator)
{
}

// TODO assert that this is called before anything else?
void Panel::Begin(const AppInput& input, const VulkanFontFace* fontFace, PanelFlags flags, Vec2Int position,
                  float32 anchorX)
{
    DEBUG_ASSERT(fontFace != nullptr);

    this->flags = flags;
    this->position = position;
    this->positionCurrent = position;
    this->anchorX = anchorX;
    this->size = Vec2Int::zero;
    this->input = &input;
    this->fontFaceDefault = fontFace;
}

bool Panel::TitleBar(const_string text, bool* minimized, Vec4 color, const VulkanFontFace* fontFace)
{
    DEBUG_ASSERT(renderCommands.size == 0);

    bool changed = false;
    if (minimized == nullptr) {
        Text(text, color, fontFace);
    }
    else {
        bool notMinimized = !(*minimized);
        changed = Checkbox(&notMinimized, text, color, fontFace);
        *minimized = !notMinimized;
        if (*minimized) {
            flags |= PanelFlag::MINIMIZED;
        }
    }

    Text(const_string::empty, color, fontFaceDefault);

    return changed;
}

void Panel::Text(const_string text, Vec4 color, const VulkanFontFace* fontFace)
{
    if (flags & PanelFlag::MINIMIZED) {
        return;
    }

    const Vec2Int textSize = AddText(this, Vec2Int::zero, text, color, fontFace);

    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += textSize.y;
    }
    else {
        positionCurrent.y -= textSize.y;
    }
    size.x = MaxInt(size.x, textSize.x);
    size.y += textSize.y;
}

bool Panel::Button(const_string text, Vec4 color, const VulkanFontFace* fontFace)
{
    if (flags & PanelFlag::MINIMIZED) {
        return false;
    }

    PanelRenderCommand* textCommand;
    const Vec2Int textSize = AddText(this, Vec2Int::zero, text, color, fontFace, &textCommand);

    const Vec2Int rectOffset = Vec2Int {
        0, (int)((float32)textCommand->commandText.fontFace->height * (MARGIN_FRACTION - 1.0f))
    };
    const bool rectHovered = AddRect(this, rectOffset, textSize, color, true);

    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += textSize.y;
    }
    else {
        positionCurrent.y -= textSize.y;
    }
    size.x = MaxInt(size.x, textSize.x);
    size.y += textSize.y;

    return rectHovered && MousePressed(*input, KM_MOUSE_LEFT);
}

bool Panel::Checkbox(bool* value, const_string text, Vec4 color, const VulkanFontFace* fontFace)
{
    DEBUG_ASSERT(value != nullptr);
    if (flags & PanelFlag::MINIMIZED) {
        return false;
    }

    const float32 CHECKED_ALPHA   = 1.0f;
    const float32 HOVERED_ALPHA   = 0.65f;
    const float32 UNCHECKED_ALPHA = 0.3f;
    bool valueChanged = false;

    const VulkanFontFace* fontToUse = fontFace == nullptr ? fontFaceDefault : fontFace;
    const int boxSize = (int)(fontToUse->height / MARGIN_FRACTION);
    const Vec2Int boxSize2 = { boxSize, boxSize };

    PanelRenderCommand* rectCommand;
    const bool rectHovered = AddRect(this, Vec2Int::zero, boxSize2, color, false, &rectCommand);
    if (rectHovered && MouseReleased(*input, KM_MOUSE_LEFT)) {
        *value = !(*value);
        valueChanged = true;
    }
    rectCommand->flags = PanelRenderCommandFlag::OVERRIDE_ALPHA;
    rectCommand->color.a = *value ? CHECKED_ALPHA : (rectHovered ? HOVERED_ALPHA : UNCHECKED_ALPHA);

    PanelRenderCommand* textCommand;
    const Vec2Int textOffset = { (int)((float32)boxSize * MARGIN_FRACTION), 0 };
    const Vec2Int textSize = AddText(this, textOffset, text, color, fontFace, &textCommand);
    if (flags & PanelFlag::GROW_DOWNWARDS) {
        textCommand->position.y -= textSize.y - boxSize;
    }

    const int sizeX = boxSize + textSize.x;
    const int sizeY = MaxInt(boxSize, textSize.y);
    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += sizeY;
    }
    else {
        positionCurrent.y -= sizeY;
    }
    size.x = MaxInt(size.x, sizeX);
    size.y += sizeY;

    return valueChanged;
}

bool Panel::Dropdown(PanelDropdownState* state, const Array<string>& values, Vec4 color, const VulkanFontFace* fontFace)
{
    DEBUG_ASSERT(state != nullptr);
    DEBUG_ASSERT(state->selected < values.size);
    if (flags & PanelFlag::MINIMIZED) {
        return false;
    }

    PanelRenderCommand* textCommand;
    const Vec2Int textSize = AddText(this, Vec2Int::zero, values[state->selected], color, fontFace, &textCommand);
    if (state->opened) {
        AddText(this, Vec2Int { textSize.x, 0 }, ToString(" |"), color, fontFace);
    }
    else {
        AddText(this, Vec2Int { textSize.x, 0 }, ToString(" -"), color, fontFace);
    }

    const Vec2Int rectOffset = Vec2Int {
        0, (int)((float32)textCommand->commandText.fontFace->height * (MARGIN_FRACTION - 1.0f))
    };
    const bool rectHovered = AddRect(this, rectOffset, textSize, color, true);

    if (MousePressed(*input, KM_MOUSE_LEFT)) {
        if (rectHovered) {
            state->opened = !state->opened;
        }
    }

    bool changed = false;
    if (state->opened) {
        Vec2Int offset = Vec2Int { 0, textSize.y };
        for (uint32 i = 0; i < values.size; i++) {
            if (i == state->selected) continue;

            const Vec2Int valueTextSize = AddText(this, offset, values[i], color, fontFace, &textCommand);
            const bool valueRectHovered = AddRect(this, offset + rectOffset, valueTextSize, color, true);
            if (valueRectHovered && MouseReleased(*input, KM_MOUSE_LEFT)) {
                state->selected = i;
                state->opened = false;
                changed = true;
                break;
            }
            offset.y += valueTextSize.y;
        }
    }

    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += textSize.y;
    }
    else {
        positionCurrent.y -= textSize.y;
    }
    size.x = MaxInt(size.x, textSize.x);
    size.y += textSize.y;

    return changed;
}

bool Panel::SliderFloat(PanelSliderState* state, float32 min, float32 max, const_string text, Vec4 color,
                        const VulkanFontFace* fontFace)
{
    UNREFERENCED_PARAMETER(text);

    DEBUG_ASSERT(state != nullptr);
    if (flags & PanelFlag::MINIMIZED) {
        return false;
    }

    const Vec2Int sliderBarSize = Vec2Int { 200, 5 };
    const Vec2Int sliderSize = Vec2Int { 10, 30 };
    const int totalHeight = 40;

    const float32 SLIDER_DRAGGED_ALPHA = 1.0f;
    const float32 SLIDER_HOVER_ALPHA = 0.8f;
    const float32 SLIDER_ALPHA = 0.6f;
    const float32 SLIDER_BAR_ALPHA = 0.4f;

    float32 sliderT = (state->value - min) / (max - min);

    bool valueChanged = false;

    const Vec2Int sliderOffset = { (int)((float32)sliderBarSize.x * sliderT), 0 };
    PanelRenderCommand* sliderRectCommand;
    const bool sliderHovered = AddRect(this, sliderOffset, sliderSize, color, false, &sliderRectCommand);
    sliderRectCommand->flags = PanelRenderCommandFlag::OVERRIDE_ALPHA;
    sliderRectCommand->color.a = state->drag ? SLIDER_DRAGGED_ALPHA : sliderHovered ? SLIDER_HOVER_ALPHA : SLIDER_ALPHA;

    if (!state->drag && sliderHovered && MouseDown(*input, KM_MOUSE_LEFT)) {
        state->drag = true;
    }
    if (state->drag && !MouseDown(*input, KM_MOUSE_LEFT)) {
        state->drag = false;
    }
    if (state->drag) {
        int newSliderX = ClampInt(input->mousePos.x - positionCurrent.x, 0, sliderBarSize.x);
        float32 newSliderT = (float32)newSliderX / sliderBarSize.x;
        state->value = newSliderT * (max - min) + min;
        valueChanged = true;
    }

    int sliderBarOffsetY = totalHeight / 2 - sliderBarSize.y;
    const Vec2Int sliderBarOffset = {
        0, flags & PanelFlag::GROW_DOWNWARDS ? sliderBarOffsetY : -sliderBarOffsetY
    };
    PanelRenderCommand* sliderBarRectCommand;
    AddRect(this, sliderBarOffset, sliderBarSize, color, false, &sliderBarRectCommand);
    sliderBarRectCommand->flags = PanelRenderCommandFlag::OVERRIDE_ALPHA;
    sliderBarRectCommand->color.a = SLIDER_BAR_ALPHA;

    const_string valueText = AllocPrintf(allocator, "%.03f", state->value);
    const Vec2Int textSize = AddText(this, Vec2Int { sliderBarSize.x + 10, 0 }, valueText, color, fontFace);

    const int sizeX = sliderBarSize.x + textSize.x;
    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += totalHeight;
    }
    else {
        positionCurrent.y -= totalHeight;
    }
    size.x = MaxInt(size.x, sizeX);
    size.y += totalHeight;

    return valueChanged;
}

bool Panel::InputText(PanelInputTextState* state, Vec4 color, const VulkanFontFace* fontFace)
{
    DEBUG_ASSERT(state != nullptr);
    if (flags & PanelFlag::MINIMIZED) {
        state->focused = false;
        return false;
    }

    bool inputChanged = false;
    if (state->focused) {
        for (uint64 i = 0; i < input->keyboardStringLen; i++) {
            char c = input->keyboardString[i];
            if (c == 8) { // backspace
                if (state->text.size > 0) {
                    state->text.RemoveLast();
                }
                inputChanged = true;
                continue;
            }
            if (c == '\n') { // enter
                state->focused = false;
                break;
            }

            if (state->text.size >= state->MAX_LENGTH) {
                continue;
            }
            inputChanged = true;
            state->text.Append(c);
        }
    }

    const float32 PRESSED_ALPHA = 0.4f;
    const float32 HOVERED_ALPHA = 0.3f;
    const float32 IDLE_ALPHA = 0.1f;
    const int BOX_MIN_WIDTH = 20;

    PanelRenderCommand* textCommand;
    const Vec2Int textSize = AddText(this, Vec2Int::zero, state->text.ToArray(), color, fontFace, &textCommand);

    const Vec2Int rectOffset = Vec2Int {
        0, (int)((float32)textCommand->commandText.fontFace->height * (MARGIN_FRACTION - 1.0f))
    };
    const Vec2Int rectSize = { MaxInt(textSize.x, BOX_MIN_WIDTH), textSize.y };
    PanelRenderCommand* rectCommand;
    const bool rectHovered = AddRect(this, rectOffset, rectSize, color, false, &rectCommand);
    rectCommand->flags = PanelRenderCommandFlag::OVERRIDE_ALPHA;
    // TODO add these alpha values as args to AddRect
    const bool pressed = rectHovered && MouseDown(*input, KM_MOUSE_LEFT);
    rectCommand->color.a = pressed ? PRESSED_ALPHA : rectHovered ? HOVERED_ALPHA : IDLE_ALPHA;

    if (MousePressed(*input, KM_MOUSE_LEFT)) {
        state->focused = rectHovered;
    }

    if (flags & PanelFlag::GROW_DOWNWARDS) {
        positionCurrent.y += rectSize.y;
    }
    else {
        positionCurrent.y -= rectSize.y;
    }
    size.x = MaxInt(size.x, rectSize.x);
    size.y += rectSize.y;

    return inputChanged;
}

bool Panel::InputInt(PanelInputIntState* state, Vec4 color, const VulkanFontFace* fontFace)
{
    bool inputChanged = false;
    if (InputText(&state->textState, color, fontFace)) {
        state->valid = StringToIntBase10(state->textState.text.ToArray(), &state->value);
        if (state->valid) {
            inputChanged = true;
        }
    }

    if (!state->valid) {
        // TODO a little dangerous / weird?
        renderCommands[renderCommands.size - 1].flags |= PanelRenderCommandFlag::OVERRIDE_COLOR;
        renderCommands[renderCommands.size - 1].color.r = 1.0f;
        renderCommands[renderCommands.size - 1].color.g = 0.0f;
        renderCommands[renderCommands.size - 1].color.b = 0.0f;
    }

    return inputChanged;
}

template <uint32 S1, uint32 S2>
void Panel::Draw(Vec2Int borderSize, Vec4 defaultColor, Vec4 backgroundColor, Vec2Int screenSize,
                 VulkanSpriteRenderState<S1>* spriteRenderState, VulkanTextRenderState<S2>* textRenderState)
{
    // Draw background
    const Vec2Int backgroundOffset = {
        Lerp(-borderSize.x, borderSize.x, anchorX),
        flags & PanelFlag::GROW_DOWNWARDS ? -borderSize.y : borderSize.y
    };
    const Vec2Int backgroundSize = size + borderSize * 2;
    const Vec2 backgroundAnchor = { anchorX, flags & PanelFlag::GROW_DOWNWARDS ? 0.0f : 1.0f };
    PushSprite((uint32)SpriteId::PIXEL, position + backgroundOffset, backgroundSize, BACKGROUND_DEPTH,
               backgroundAnchor, backgroundColor, screenSize, spriteRenderState);

    for (uint32 i = 0; i < renderCommands.size; i++) {
        const PanelRenderCommand& command = renderCommands[i];
        Vec4 color = command.flags & PanelRenderCommandFlag::OVERRIDE_COLOR ? command.color : defaultColor;
        if (command.flags & PanelRenderCommandFlag::OVERRIDE_ALPHA) {
            color.a = command.color.a;
        }

        switch (command.type) {
            case PanelRenderCommandType::RECT: {
                PushSprite((uint32)SpriteId::PIXEL, command.position, command.commandRect.size, PANEL_ITEM_DEPTH,
                           command.anchor, color, screenSize, spriteRenderState);
            } break;
            case PanelRenderCommandType::TEXT: {
                PushText(*command.commandText.fontFace, command.commandText.text, command.position, PANEL_TEXT_DEPTH,
                         command.anchor.x, color, screenSize, textRenderState);
            } break;
        }
    }
}
