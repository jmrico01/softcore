#pragma once

#include <km_common/km_load_obj.h>
#include <km_common/km_memory.h>
#include <km_common/app/km_app.h>

// MODEL INDICES:
// 0, 1, 2 - some rocks
// 3 - box, left
// 4 - head
// 5 - box, right
// 6 - some other rock
// 7 - walls

#define RESTRICT_LIGHTING 1
const int MODELS_TO_LIGHT[] = {
    3
};

#define RESTRICT_OCCLUSION 0
#define MODEL_TO_OCCLUDE 3

#define RESTRICT_WALL 1
const uint64 PLANE_LEFT  = 0;
const uint64 PLANE_BACK  = 1;
const uint64 PLANE_RIGHT = 2;
const uint64 PLANE_FLOOR = 3;
#define PLANE_TO_LIGHT PLANE_FLOOR

const uint32 LIGHTMAP_NUM_BOUNCES = 1;
const VkFilter LIGHTMAP_TEXTURE_FILTER = VK_FILTER_LINEAR;

const float32 RESOLUTION_PER_WORLD_UNIT = 64.0f;
const uint32 NUM_HEMISPHERE_SAMPLES = 64;
const uint32 SAMPLES_PER_GROUP = 8; // AVX width, won't compile if != 8
static_assert(NUM_HEMISPHERE_SAMPLES % SAMPLES_PER_GROUP == 0);
const uint32 NUM_HEMISPHERE_SAMPLE_GROUPS = NUM_HEMISPHERE_SAMPLES / SAMPLES_PER_GROUP;

bool GenerateLightmaps(const LoadObjResult& obj, uint32 bounces, AppWorkQueue* queue, LinearAllocator* allocator,
                       const_string lightmapDirPath);
