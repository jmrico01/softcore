#pragma once

#include <km_common/km_array.h>

#include "mesh.h"

constexpr Vec3Int BLOCKS_SIZE = { 128, 128, 32 };
constexpr Vec3Int BLOCK_ORIGIN = { BLOCKS_SIZE.x / 2, BLOCKS_SIZE.y / 2, 4 };

enum class GridTemplateId
{
    FLATGRASS,
    TUNNELS
};

enum class BlockId
{
    NONE = 0,
    SIDEWALK,
    STREET,
    BUILDING,

    COUNT
};

struct Block
{
    BlockId id;
};

using BlockGrid = StaticArray<StaticArray<StaticArray<Block, BLOCKS_SIZE.x>, BLOCKS_SIZE.y>, BLOCKS_SIZE.z>;

struct BlockRenderInfo
{
    MeshId meshId;
    Mat4 model;
    Vec3 color;
};

struct Mob
{
    Vec3 pos;
    float32 yaw;
    Box hitbox;
    float32 collapseT;
    bool collapsed;
};

struct BlockUpdate
{
    Vec3Int index;
    Block block;
};

struct LevelData
{
    static const uint32 NUM_BLOCKS = BLOCKS_SIZE.x * BLOCKS_SIZE.y * BLOCKS_SIZE.z;
    static const uint32 MAX_MOBS = 1024;

    GridTemplateId gridTemplateId;
    float32 blockSize;
    BlockGrid grid;
    FixedArray<BlockRenderInfo, NUM_BLOCKS> gridRenderInfo;

    FixedArray<Mob, MAX_MOBS> mobs;
    uint32 collapsingMobIndex;
};

bool IsWalkable(const LevelData& levelData, Vec3Int blockIndex);

void SubmitBlockUpdates(Array<BlockUpdate> updates, LevelData* levelData);

bool LoadLevel(const_string levelName, LevelData* levelData, LinearAllocator* allocator);
bool SaveLevel(const_string levelName, const LevelData& levelData, LinearAllocator* allocator);
Array<string> GetSavedLevels(LinearAllocator* allocator);

void UpdateGridRenderInfo(const BlockGrid& blockGrid, float32 blockSize, Vec3Int blockOrigin,
                          FixedArray<BlockRenderInfo, LevelData::NUM_BLOCKS>* renderInfo);
