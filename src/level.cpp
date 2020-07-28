#include "level.h"

string GetLevelFilePath(const_string levelName, LinearAllocator* allocator)
{
    string path = AllocPrintf(allocator, "data/levels/%.*s.blockgrid", levelName.size, levelName.data);
    DEBUG_ASSERT(path.data != nullptr);
    return path;
}

bool IsWalkable(const LevelData& levelData, Vec3Int blockIndex)
{
    if (blockIndex.x < 0 || blockIndex.y < 0 || blockIndex.z < 0) {
        return false;
    }

    const Block& block = levelData.grid[blockIndex.z][blockIndex.y][blockIndex.x];
    return block.id == BlockId::NONE;
}

void SubmitBlockUpdates(Array<BlockUpdate> updates, LevelData* levelData)
{
    for (uint32 i = 0; i < updates.size; i++) {
        const Vec3Int index = updates[i].index;
        levelData->grid[index.z][index.y][index.x] = updates[i].block;
    }

    UpdateGridRenderInfo(levelData->grid, levelData->blockSize, BLOCK_ORIGIN, &levelData->gridRenderInfo);
}

bool LoadLevel(const_string levelName, LevelData* levelData, LinearAllocator* allocator)
{
    string path = GetLevelFilePath(levelName, allocator);
    Array<uint8> data = LoadEntireFile(path, allocator);
    if (data.data == nullptr) {
        return false;
    }
    if (data.size != sizeof(Vec3Int) + sizeof(BlockGrid)) {
        return false;
    }

    const Vec3Int* fileBlockGridSize = (Vec3Int*)data.data;
    if ((uint32)fileBlockGridSize->x != levelData->grid[0][0].SIZE
        || (uint32)fileBlockGridSize->y != levelData->grid[0].SIZE
        || (uint32)fileBlockGridSize->z != levelData->grid.SIZE) {
        return false;
    }

    MemCopy(&levelData->grid, data.data + sizeof(Vec3Int), sizeof(BlockGrid));

    UpdateGridRenderInfo(levelData->grid, levelData->blockSize, BLOCK_ORIGIN, &levelData->gridRenderInfo);

    return true;
}

bool SaveLevel(const_string levelName, const LevelData& levelData, LinearAllocator* allocator)
{
    const Vec3Int blockGridSize = { (int)levelData.grid[0][0].SIZE, (int)levelData.grid[0].SIZE, (int)levelData.grid.SIZE };

    Array<uint8> data = allocator->NewArray<uint8>(sizeof(Vec3Int) + sizeof(BlockGrid));
    MemCopy(data.data, &blockGridSize, sizeof(Vec3Int));
    MemCopy(data.data + sizeof(Vec3Int), &levelData.grid, sizeof(BlockGrid));

    string path = GetLevelFilePath(levelName, allocator);
    return WriteFile(path, data, false);
}

Array<string> GetSavedLevels(LinearAllocator* allocator)
{
    const Array<string> levelFiles = ListDir(ToString("data/levels"), allocator);
    if (levelFiles.data == nullptr) {
        return { .size = 0, .data = nullptr };
    }

    DynamicArray<string, LinearAllocator> validLevelFiles(allocator);
    for (uint32 i = 0; i < levelFiles.size; i++) {
        const uint32 find = SubstringSearch(levelFiles[i], ToString(".blockgrid"));
        if (find != levelFiles[i].size) {
            validLevelFiles.Append(levelFiles[i].SliceTo(find));
        }
    }

    return validLevelFiles.ToArray();
}

void UpdateGridRenderInfo(const BlockGrid& blockGrid, float32 blockSize, Vec3Int blockOrigin,
                          FixedArray<BlockRenderInfo, LevelData::NUM_BLOCKS>* renderInfo)
{
    const Mat4 scale = Scale(blockSize);
    DEFINE_BLOCK_ROTATION_MATRICES;

    renderInfo->Clear();

    const uint32 sizeZ = blockGrid.SIZE;
    for (uint32 z = 0; z < sizeZ; z++) {
        const uint32 sizeY = blockGrid[z].SIZE;
        for (uint32 y = 0; y < sizeY; y++) {
            const uint32 sizeX = blockGrid[z][y].SIZE;
            for (uint32 x = 0; x < sizeX; x++) {
                const Block block = blockGrid[z][y][x];

                bool drawTop    = z != sizeZ - 1 && blockGrid[z + 1][y][x].id == BlockId::NONE;
                bool drawLeft   = y != sizeY - 1 && blockGrid[z][y + 1][x].id == BlockId::NONE;
                bool drawFront  = x != sizeX - 1 && blockGrid[z][y][x + 1].id == BlockId::NONE;
                bool drawBottom = z != 0 && blockGrid[z - 1][y][x].id == BlockId::NONE;
                bool drawRight  = y != 0 && blockGrid[z][y - 1][x].id == BlockId::NONE;
                bool drawBack   = x != 0 && blockGrid[z][y][x - 1].id == BlockId::NONE;

                Vec3 color = Vec3::zero;
                switch (block.id) {
                    case BlockId::NONE: {
                        drawTop = false;
                        drawLeft = false;
                        drawFront = false;
                        drawBottom = false;
                        drawRight = false;
                        drawBack = false;
                    } break;
                    case BlockId::SIDEWALK: {
                        color = Vec3::one * 0.7f;
                    } break;
                    case BlockId::STREET: {
                        color = Vec3::one * 0.5f;
                    } break;
                    case BlockId::BUILDING: {
                        color = Vec3::one * 0.8f;
                    } break;
                }

                const Vec3 pos = BlockIndexToWorldPos(Vec3Int { (int)x, (int)y, (int)z }, blockSize, blockOrigin);
                const Mat4 posTransform = Translate(pos);
                if (drawTop) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * top;
                    info->color = color;
                }
                if (drawBottom) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * bottom;
                    info->color = color;
                }
                if (drawLeft) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * left;
                    info->color = color;
                }
                if (drawRight) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * right;
                    info->color = color;
                }
                if (drawFront) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * front;
                    info->color = color;
                }
                if (drawBack) {
                    BlockRenderInfo* info = renderInfo->Append();
                    info->meshId = MeshId::TILE;
                    info->model = posTransform * scale * back;
                    info->color = color;
                }
            }
        }
    }
}

void GenerateCityBlocks(uint32 streetSize, uint32 sidewalkSize, uint32 buildingSize, uint32 buildingHeight,
                        LinearAllocator* allocator, BlockGrid* blockGrid)
{
    MemSet(blockGrid, 0, sizeof(BlockGrid));

    const uint32 cityBlockSize = streetSize + sidewalkSize * 2 + buildingSize;

    Array<uint8> streetMask = allocator->NewArray<uint8>(cityBlockSize);
    for (uint32 i = 0; i < cityBlockSize; i++) {
        streetMask[i] = i < streetSize;
    }
    Array<uint8> sidewalkMask = allocator->NewArray<uint8>(cityBlockSize);
    for (uint32 i = 0; i < cityBlockSize; i++) {
        sidewalkMask[i] = (i >= streetSize && i < (streetSize + sidewalkSize))
            || (i >= (streetSize + sidewalkSize + buildingSize) && i < (streetSize + 2 * sidewalkSize + buildingSize));
    }
    Array<uint8> buildingMask = allocator->NewArray<uint8>(cityBlockSize);
    for (uint32 i = 0; i < cityBlockSize; i++) {
        buildingMask[i] = i >= (streetSize + sidewalkSize) && i < (streetSize + sidewalkSize + buildingSize);
    }

    static_assert(BLOCK_ORIGIN.z > 0);
    for (uint32 z = 0; z < blockGrid->SIZE; z++) {
        for (uint32 y = 0; y < (*blockGrid)[z].SIZE; y++) {
            for (uint32 x = 0; x < (*blockGrid)[z][y].SIZE; x++) {
                const uint32 cityBlockX = x % cityBlockSize;
                const uint32 cityBlockY = y % cityBlockSize;

                const uint8 building = buildingMask[cityBlockX] * buildingMask[cityBlockY];
                if ((int)z == BLOCK_ORIGIN.z - 1) {
                    const uint8 street = streetMask[cityBlockX] + streetMask[cityBlockY];
                    const uint8 sidewalk = sidewalkMask[cityBlockX] + sidewalkMask[cityBlockY];
                    if (street) {
                        (*blockGrid)[z][y][x].id = BlockId::STREET;
                    }
                    else if (sidewalk) {
                        (*blockGrid)[z][y][x].id = BlockId::SIDEWALK;
                    }
                    else if (building) {
                        (*blockGrid)[z][y][x].id = BlockId::BUILDING;
                    }
                }
                else if (building && z - BLOCK_ORIGIN.z < buildingHeight) {
                    (*blockGrid)[z][y][x].id = BlockId::BUILDING;
                }
            }
        }
    }
}
