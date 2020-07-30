#include "raytracer.h"

#include <intrin.h>
#include <stb_image_write.h>
#include <time.h>

struct DebugTimer
{
    static bool initialized;
    static uint64 win32Freq;

    uint64 cycles;
    uint64 win32Time;
};

bool DebugTimer::initialized = false;
uint64 DebugTimer::win32Freq;

DebugTimer StartDebugTimer()
{
    if (!DebugTimer::initialized) {
        DebugTimer::initialized = true;
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        DebugTimer::win32Freq = freq.QuadPart;
    }

    DebugTimer timer;
    LARGE_INTEGER win32Time;
    QueryPerformanceCounter(&win32Time);
    timer.win32Time = win32Time.QuadPart;
    timer.cycles = __rdtsc();
    return timer;
}

void StopDebugTimer(DebugTimer* timer)
{
    LARGE_INTEGER win32End;
    QueryPerformanceCounter(&win32End);
    timer->cycles = __rdtsc() - timer->cycles;
    timer->win32Time = win32End.QuadPart - timer->win32Time;
}

void StopAndPrintDebugTimer(DebugTimer* timer)
{
    StopDebugTimer(timer);
    const float32 win32Time = (float32)timer->win32Time / DebugTimer::win32Freq * 1000.0f;
    LOG_INFO("Timer: %.03fms | %llu MC\n", win32Time, timer->cycles / 1000000);
}

RaycastGeometry CreateRaycastGeometry(const LoadObjResult& obj, LinearAllocator* allocator)
{
    RaycastGeometry geometry;
    geometry.meshes = allocator->NewArray<RaycastMesh>(obj.models.size);
    if (geometry.meshes.data == nullptr) {
        return geometry;
    }

    const Vec3 vertexColor = Vec3::zero;

    for (uint32 i = 0; i < obj.models.size; i++) {
        RaycastMesh& mesh = geometry.meshes[i];

        // For light cones min scene
        if (i == 2 || i == 3 || i == 4) { // cubes
            mesh.smoothness = 0.0f;
            mesh.albedo = Vec3 { 1.0f, 0.0f, 0.0f };
            mesh.emission = 2.0f;
            mesh.emissionColor = Vec3 { 1.0f, 0.0f, 0.0f };
        }
        else if (i == obj.models.size - 3) { // light cone
            mesh.smoothness = 0.0f;
            mesh.albedo = Vec3 { 0.0f, 1.0f, 0.42f };
            mesh.emission = 4.0f;
            mesh.emissionColor = Vec3 { 0.0f, 1.0f, 0.42f };
        }
        else if (i == obj.models.size - 1) { // light overhead
            mesh.smoothness = 0.0f;
            mesh.albedo = Vec3 { 0.9294f, 0.651f, 1.0f };
            mesh.emission = 1.0f;
            mesh.emissionColor = Vec3 { 0.9294f, 0.651f, 1.0f };
        }
        else { // everything else
            mesh.smoothness = 0.8f;
            mesh.albedo = Vec3 { 1.0f, 1.0f, 1.0f };
            mesh.emission = 0.0f;
            mesh.emissionColor = Vec3::zero;
        }

        const uint32 numTriangles = obj.models[i].triangles.size + obj.models[i].quads.size * 2;
        mesh.triangles = allocator->NewArray<RaycastTriangle>(numTriangles);
        if (mesh.triangles.data == nullptr) {
            LOG_ERROR("Failed to allocate triangles for raycast mesh %lu\n", i);
            geometry.meshes.data = nullptr;
            return geometry;
        }

        // Fill in triangle geometry data
        for (uint32 j = 0; j < obj.models[i].triangles.size; j++) {
            const ObjTriangle& t = obj.models[i].triangles[j];
            const Vec3 normal = CalculateTriangleUnitNormal(t.v[0].pos, t.v[1].pos, t.v[2].pos);

            for (int k = 0; k < 3; k++) {
                mesh.triangles[j].pos[k] = t.v[k].pos;
            }
            mesh.triangles[j].normal = normal;
        }
        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const uint32 ind = obj.models[i].triangles.size + j * 2;
            const ObjQuad& q = obj.models[i].quads[j];
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                mesh.triangles[ind].pos[k] = q.v[k].pos;
            }
            mesh.triangles[ind].normal = normal;

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                mesh.triangles[ind + 1].pos[k] = q.v[quadInd].pos;
            }
            mesh.triangles[ind + 1].normal = normal;
        }

        float32 surfaceArea = 0.0f;
        for (uint32 j = 0; j < mesh.triangles.size; j++) {
            const RaycastTriangle& t = mesh.triangles[j];
            surfaceArea += TriangleArea(t.pos[0], t.pos[1], t.pos[2]);
        }

        // Calculate AABB
        mesh.aabb.min = Vec3::one * 1e8;
        mesh.aabb.max = -Vec3::one * 1e8;
        for (uint32 j = 0; j < mesh.triangles.size; j++) {
            for (int k = 0; k < 3; k++) {
                const Vec3 v = mesh.triangles[j].pos[k];
                for (int e = 0; e < 3; e++) {
                    mesh.aabb.min.e[e] = MinFloat32(mesh.aabb.min.e[e], v.e[e]);
                    mesh.aabb.max.e[e] = MaxFloat32(mesh.aabb.max.e[e], v.e[e]);
                }
            }
        }
    }

    return geometry;
}

Vec3 RaycastColor(Vec3 rayOrigin, Vec3 rayDir, const RaycastGeometry& geometry)
{
    const uint32 SAMPLES = 8;
    const uint32 BOUNCES = 4;

    const float32 sampleWeight = 1.0f / (float32)SAMPLES;

    Vec3 color = Vec3::zero;
    for (uint32 s = 0; s < SAMPLES; s++) {
        for (uint32 b = 0; b < BOUNCES; b++) {
            Vec3 inverseRayDir = Reciprocal(rayDir);

            uint32 closestMeshInd = geometry.meshes.size;
            uint32 closestTriangleInd = 0;
            float32 closestTriangleDist = 1e8;
            for (uint32 i = 0; i < geometry.meshes.size; i++) {
                const RaycastMesh& mesh = geometry.meshes[i];
                // TODO also compare min distance with closestTriangleDist ?
                float32 tAabb;
                bool intersect = RayAxisAlignedBoxIntersection(rayOrigin, inverseRayDir, mesh.aabb, &tAabb);
                if (!intersect || tAabb >= closestTriangleDist) {
                    continue;
                }

                for (uint32 j = 0; j < mesh.triangles.size; j++) {
                    const RaycastTriangle& triangle = mesh.triangles[j];
                    float32 t;
                    const bool tIntersect = RayTriangleIntersection(rayOrigin, rayDir,
                                                                    triangle.pos[0], triangle.pos[1], triangle.pos[2],
                                                                    &t);
                    if (tIntersect && t > 0.0f && t < closestTriangleDist) {
                        closestMeshInd = i;
                        closestTriangleInd = j;
                        closestTriangleDist = t;
                    }
                }
            }

            if (closestMeshInd == geometry.meshes.size) {
                break;
            }

            const RaycastMesh& hitMesh = geometry.meshes[closestMeshInd];
            if (hitMesh.emission > 0.0f) {
                color += sampleWeight * hitMesh.emissionColor;
                break;
            }
            else {
                rayOrigin += rayDir * closestTriangleDist;

                const Vec3 normal = hitMesh.triangles[closestTriangleInd].normal;
                const Vec3 pureBounce = rayDir - 2.0f * Dot(rayDir, normal) * normal;
                const Vec3 randomBounce = NormalizeOrZero(Vec3 {
                                                              RandFloat32(-1.0f, 1.0f),
                                                              RandFloat32(-1.0f, 1.0f),
                                                              RandFloat32(-1.0f, 1.0f)
                                                          });
                rayDir = NormalizeOrZero(Lerp(randomBounce, pureBounce, hitMesh.smoothness));
            }
        }
    }

    return color;
}

struct RaycastThreadWork
{
    Vec3 cameraPos;
    Vec3 rayDir;
    const RaycastGeometry* geometry;
    uint32* pixel;
};

APP_WORK_QUEUE_CALLBACK_FUNCTION(RaycastThreadProc)
{
    UNREFERENCED_PARAMETER(queue);

    const RaycastThreadWork* work = (RaycastThreadWork*)data;

    const Vec3 color = RaycastColor(work->cameraPos, work->rayDir, *work->geometry);
    const uint8 r = (uint8)(color.r * 255.0f);
    const uint8 g = (uint8)(color.g * 255.0f);
    const uint8 b = (uint8)(color.b * 255.0f);
    *(work->pixel) = ((uint32)r << 16) + ((uint32)g << 8) + b;
}

void RaytraceRender(Vec3 cameraPos, Quat cameraRot, const RaycastGeometry& geometry,
                    uint32 width, uint32 height, uint32* pixels, LinearAllocator* allocator, AppWorkQueue* queue)
{
    UNREFERENCED_PARAMETER(allocator);

    const uint8 RAY_DECAY_FRAMES = 255;

    // NOTE this assumes little-endian
    uint8* pixelChannels = (uint8*)pixels;
    for (uint32 i = 0; i < width * height; i++) {
        uint8* alpha = pixelChannels + 3;
        if (*alpha != 0xff) {
            (*alpha) += 1;
            if (*alpha >= RAY_DECAY_FRAMES) {
                *(pixelChannels + 0) = 0;
                *(pixelChannels + 1) = 0;
                *(pixelChannels + 2) = 0;
                *alpha = 0xff;
            }
        }
        pixelChannels += 4;
    }

    // The idea is that, over 1 second, we'll get a whole frame-pixels' worth of rays
    const uint32 NUM_RAYS_PER_FRAME = (uint32)(1.0f / 60.0f * width * height);
    Array<Vec2Int> randomPixels = allocator->NewArray<Vec2Int>(NUM_RAYS_PER_FRAME);
    for (uint32 i = 0; i < randomPixels.size; i++) {
        randomPixels[i].x = RandInt(width);
        randomPixels[i].y = RandInt(height);
    }

    const Quat inverseCameraRot = Inverse(cameraRot);
    const Vec3 cameraUp = inverseCameraRot * Vec3::unitZ;
    const Vec3 cameraForward = inverseCameraRot * Vec3::unitX;
    const Vec3 cameraLeft = inverseCameraRot * Vec3::unitY;

    const float32 filmDist = 1.0f;
    const float32 filmWidth = 2.0f;
    const float32 filmHeight = filmWidth * (float32)height / (float32)width;

    const Vec3 filmTopLeft = cameraPos + cameraForward * filmDist
        + filmWidth / 2.0f * cameraUp + filmHeight / 2.0f * cameraLeft;

    Array<RaycastThreadWork> workEntries = allocator->NewArray<RaycastThreadWork>(randomPixels.size);

    for (uint32 i = 0; i < randomPixels.size; i++) {
        const Vec3 filmOffsetX = -cameraLeft * (float32)randomPixels[i].x / (float32)width * filmWidth;
        const Vec3 filmOffsetY = -cameraUp * (float32)randomPixels[i].y / (float32)height * filmHeight;
        const Vec3 filmPos = filmTopLeft + filmOffsetX + filmOffsetY;
        const Vec3 rayDir = Normalize(filmPos - cameraPos);

        workEntries[i].cameraPos = cameraPos;
        workEntries[i].rayDir = rayDir;
        workEntries[i].geometry = &geometry;
        workEntries[i].pixel = pixels + randomPixels[i].y * width + randomPixels[i].x;

        if (!TryAddWork(queue, &RaycastThreadProc, &workEntries[i])) {
            CompleteAllWork(queue);
            i--;
            continue;
        }
    }

    CompleteAllWork(queue);
}

#if 0
internal void GenerateHemisphereSamples(Array<Vec3> samples)
{
    for (uint32 i = 0; i < samples.size; i++) {
        Vec3 dir;
        do {
            dir.x = RandFloat32();
            dir.y = RandFloat32(-1.0f, 1.0f);
            dir.z = RandFloat32(-1.0f, 1.0f);
        } while (MagSq(dir) > 1.0f);

        samples[i] = Normalize(dir);
    }
}

struct SampleGroup
{
    StaticArray<Vec3, SAMPLES_PER_GROUP> group;
};

// Provides a metric that is lower the closer together 2 unit direction vectors are (can be negative)
internal float32 HemisphereDirCloseness(Vec3 dir1, Vec3 dir2)
{
    return -Dot(dir1, dir2);
}

internal float32 HemisphereGroupCloseness(const Array<Vec3>& samples, int groupIndices[SAMPLES_PER_GROUP])
{
    float32 closeness = 0.0f;
    for (uint32 i = 0; i < 8; i++) {
        for (uint32 j = i + 1; j < 8; j++) {
            closeness += HemisphereDirCloseness(samples[groupIndices[i]], samples[groupIndices[j]]);
        }
    }
    return closeness;
}

internal bool GenerateHemisphereSampleGroups(Array<SampleGroup> sampleGroups, LinearAllocator* allocator)
{
    ALLOCATOR_SCOPE_RESET(*allocator);

    Array<Vec3> samples = {
        .size = sampleGroups.size * SAMPLES_PER_GROUP,
        .data = &sampleGroups[0].group[0]
    };
    GenerateHemisphereSamples(samples);

    DEBUG_ASSERT(samples.size % 8 == 0);
    const uint32 numGroups = samples.size / 8;

    Array<int> indices = allocator->NewArray<int>(samples.size);
    if (indices.data == nullptr) {
        return false;
    }
    for (uint32 i = 0; i < samples.size; i++) {
        indices[i] = i;
    }

    const unsigned int seed = (unsigned int)time(NULL);
    srand(seed);

    const uint32 ITERATIONS = 10000;
    Array<int> minIndices = allocator->NewArray<int>(samples.size);
    if (minIndices.data == nullptr) {
        return false;
    }
    float32 minCloseness = 1e8;
    for (uint32 n = 0; n < ITERATIONS; n++) {
        indices.Shuffle();
        float32 totalCloseness = 0.0f;
        for (uint32 i = 0; i < numGroups; i++) {
            totalCloseness += HemisphereGroupCloseness(samples, &indices[i * 8]);
        }
        if (totalCloseness < minCloseness) {
            minIndices.CopyFrom(indices);
            minCloseness = totalCloseness;
        }
    }

    Array<Vec3> samplesCopy = allocator->NewArray<Vec3>(samples.size);
    if (samplesCopy.data == nullptr) {
        return false;
    }
    samplesCopy.CopyFrom(samples);
    for (uint32 i = 0; i < samples.size; i++) {
        samples[i] = samplesCopy[minIndices[i]];
    }

    // Debug purposes only, log closeness
    for (uint32 i = 0; i < samples.size; i++) {
        indices[i] = i;
    }
    float32 totalCloseness = 0.0f;
    for (uint32 i = 0; i < numGroups; i++) {
        totalCloseness += HemisphereGroupCloseness(samples, &indices[i * 8]);
    }
    LOG_INFO("hemisphere closeness: %f\n", totalCloseness);

    return true;
}

internal Vec3 RaycastColor(Array<SampleGroup> sampleGroups, Vec3 pos, Vec3 normal, const RaycastGeometry& geometry)
{
    // For reference: left wall is at  Y =  1.498721
    //                right wall is at Y = -1.544835
    const LightRect LIGHT_RECTS[] = {
        {
            .origin = { 4.0f, 1.498721f - 0.005f, 2.24f },
            .width = { -2.0f, 0.0f, 0.0f },
            .height = { 0.0f, 0.0f, -2.2f },
            .color = Vec3 { 1.0f, 0.0f, 0.0f },
            .intensity = 2.0f
        },
        {
            .origin = { 2.0f, -1.544835f + 0.005f, 2.24f },
            .width = { 2.0f, 0.0f, 0.0f },
            .height = { 0.0f, 0.0f, -2.2f },
            .color = Vec3 { 0.0f, 0.0f, 1.0f },
            .intensity = 2.0f
        },
    };

    const float32 MATERIAL_REFLECTANCE = 0.3f;

    const uint32 numSamples = sampleGroups.size * SAMPLES_PER_GROUP;

    // NOTE this will do unknown-ish things with "up" direction
    const Quat xToNormalRot = QuatRotBetweenVectors(Vec3::unitX, normal);

    const float32 largeFloat = 1e8;
    const __m256 largeFloat8 = _mm256_set1_ps(largeFloat);
    const __m256 zero8 = _mm256_setzero_ps();
    const Vec3_8 pos8 = Set1Vec3_8(pos);
    const Quat_8 xToNormalRot8 = Set1Quat_8(xToNormalRot);
    const float32 offset = 0.001f;
    const __m256 offset8 = _mm256_set1_ps(offset);
    const float32 sampleContribution = 1.0f / (float32)numSamples;

    static_assert(SAMPLES_PER_GROUP == 8);

    Vec3 outputColor = Vec3::zero;
    for (uint32 m = 0; m < sampleGroups.size; m++) {
        const Vec3_8 sample8 = SetVec3_8(sampleGroups[m].group);
        const Vec3_8 sampleNormal8 = Multiply_8(xToNormalRot8, sample8);
        const Vec3_8 sampleNormalInv8 = Inverse_8(sampleNormal8);
        const Vec3_8 originOffset8 = Add_8(pos8, Multiply_8(sampleNormal8, offset8));

        __m256i closestMeshInd8 = _mm256_set1_epi32(geometry.meshes.size);
        __m256i closestTriangleInd8 = _mm256_undefined_si256();
        __m256 closestTriangleDist8 = largeFloat8;
        for (uint32 i = 0; i < geometry.meshes.size; i++) {
#if RESTRICT_LIGHTING && RESTRICT_OCCLUSION
            if (i != MODEL_TO_OCCLUDE) continue;
#endif
            const RaycastMesh& mesh = geometry.meshes[i];
            const __m256i meshInd8 = _mm256_set1_epi32(i);
            // TODO also return min distance, and compare with closestTriangleDist8 ?
            const __m256 intersect8 = RayAxisAlignedBoxIntersection_8(originOffset8, sampleNormalInv8,
                                                                      mesh.aabb.min, mesh.aabb.max);
            const int allZero = _mm256_testc_ps(zero8, intersect8);
            if (allZero) {
                continue;
            }

            for (uint32 j = 0; j < mesh.triangles.size; j++) {
                const RaycastTriangle& triangle = mesh.triangles[j];
                const __m256i triangleInd8 = _mm256_set1_epi32(j);
                __m256 t8;
                const __m256 tIntersect8 = RayTriangleIntersection_8(originOffset8, sampleNormal8,
                                                                     triangle.pos[0], triangle.pos[1], triangle.pos[2],
                                                                     &t8);

                const __m256 closerMask8 = _mm256_and_ps(_mm256_cmp_ps(t8, closestTriangleDist8, _CMP_LT_OQ), tIntersect8);
                closestTriangleDist8 = _mm256_blendv_ps(closestTriangleDist8, t8, closerMask8);

                const __m256i closerMask8i = _mm256_castps_si256(closerMask8);
                closestMeshInd8 = _mm256_blendv_epi8(closestMeshInd8, meshInd8, closerMask8i);
                closestTriangleInd8 = _mm256_blendv_epi8(closestTriangleInd8, triangleInd8, closerMask8i);
            }
        }

        __m256i closestLightInd8 = _mm256_set1_epi32(C_ARRAY_LENGTH(LIGHT_RECTS));
        __m256 closestLightDist8 = largeFloat8;
        for (int l = 0; l < C_ARRAY_LENGTH(LIGHT_RECTS); l++) {
            const __m256i lightInd = _mm256_set1_epi32(l);
            const Vec3_8 lightRectOrigin8 = Set1Vec3_8(LIGHT_RECTS[l].origin);
            const Vec3_8 lightRectNormal8 = Set1Vec3_8(Normalize(Cross(LIGHT_RECTS[l].width, LIGHT_RECTS[l].height)));

            const Vec3_8 lightWidth8 = Set1Vec3_8(LIGHT_RECTS[l].width);
            const Vec3_8 lightHeight8 = Set1Vec3_8(LIGHT_RECTS[l].height);
            const __m256 lightRectWidth8 = Mag_8(lightWidth8);
            const Vec3_8 lightRectUnitWidth8 = Divide_8(lightWidth8, lightRectWidth8);
            const __m256 lightRectHeight8 = Mag_8(lightHeight8);
            const Vec3_8 lightRectUnitHeight8 = Divide_8(lightHeight8, lightRectHeight8);

            __m256 t8;
            const __m256 pIntersect8 = RayPlaneIntersection_8(pos8, sampleNormal8,
                                                              lightRectOrigin8, lightRectNormal8, &t8);

            // Pixels are lit only when 0.0f <= t < closestTriangleDist
            __m256 lit8 = _mm256_and_ps(pIntersect8, _mm256_cmp_ps(zero8, t8, _CMP_LE_OQ));
            lit8 = _mm256_and_ps(lit8, _mm256_cmp_ps(t8, closestTriangleDist8, _CMP_LT_OQ));

            const Vec3_8 intersect8 = Add_8(pos8, Multiply_8(sampleNormal8, t8));
            const Vec3_8 rectOriginToIntersect8 = Subtract_8(intersect8, lightRectOrigin8);

            const __m256 projWidth8 = Dot_8(rectOriginToIntersect8, lightRectUnitWidth8);
            lit8 = _mm256_and_ps(lit8, _mm256_cmp_ps(zero8, projWidth8, _CMP_LE_OQ));
            lit8 = _mm256_and_ps(lit8, _mm256_cmp_ps(projWidth8, lightRectWidth8, _CMP_LE_OQ));

            const __m256 projHeight8 = Dot_8(rectOriginToIntersect8, lightRectUnitHeight8);
            lit8 = _mm256_and_ps(lit8, _mm256_cmp_ps(zero8, projHeight8, _CMP_LE_OQ));
            lit8 = _mm256_and_ps(lit8, _mm256_cmp_ps(projHeight8, lightRectHeight8, _CMP_LE_OQ));

            closestLightDist8 = _mm256_blendv_ps(closestLightDist8, t8, lit8);
            const __m256i lit8i = _mm256_castps_si256(lit8);
            closestLightInd8 = _mm256_blendv_epi8(closestLightInd8, lightInd, lit8i);
        }

#if 0
        // TODO wow... there's definitely a better way to do this... right?
        const int32 lightInds[8] = {
            _mm256_extract_epi32(closestLightInd8, 0),
            _mm256_extract_epi32(closestLightInd8, 1),
            _mm256_extract_epi32(closestLightInd8, 2),
            _mm256_extract_epi32(closestLightInd8, 3),
            _mm256_extract_epi32(closestLightInd8, 4),
            _mm256_extract_epi32(closestLightInd8, 5),
            _mm256_extract_epi32(closestLightInd8, 6),
            _mm256_extract_epi32(closestLightInd8, 7),
        };
        const int32 meshInds[8] = {
            _mm256_extract_epi32(closestMeshInd8, 0),
            _mm256_extract_epi32(closestMeshInd8, 1),
            _mm256_extract_epi32(closestMeshInd8, 2),
            _mm256_extract_epi32(closestMeshInd8, 3),
            _mm256_extract_epi32(closestMeshInd8, 4),
            _mm256_extract_epi32(closestMeshInd8, 5),
            _mm256_extract_epi32(closestMeshInd8, 6),
            _mm256_extract_epi32(closestMeshInd8, 7),
        };
        const int32 triangleInds[8] = {
            _mm256_extract_epi32(closestTriangleInd8, 0),
            _mm256_extract_epi32(closestTriangleInd8, 1),
            _mm256_extract_epi32(closestTriangleInd8, 2),
            _mm256_extract_epi32(closestTriangleInd8, 3),
            _mm256_extract_epi32(closestTriangleInd8, 4),
            _mm256_extract_epi32(closestTriangleInd8, 5),
            _mm256_extract_epi32(closestTriangleInd8, 6),
            _mm256_extract_epi32(closestTriangleInd8, 7),
        };
        for (int i = 0; i < 8; i++) {
            if (lightInds[i] != C_ARRAY_LENGTH(LIGHT_RECTS)) {
                const float32 lightIntensity = LIGHT_RECTS[lightInds[i]].intensity;
                outputColor += lightIntensity * sampleContribution * LIGHT_RECTS[lightInds[i]].color;
            }
            else if ((uint32)meshInds[i] != geometry.meshes.size) {
                const RaycastMesh& mesh = geometry.meshes[meshInds[i]];
                const Lightmap& lightmap = mesh.lightmap;
                const int squareSize = lightmap.squareSize;
                const RaycastTriangle& triangle = mesh.triangles[triangleInds[i]];
                const Vec3 sampleNormal = xToNormalRot * sampleGroups[m].group[i];
                const Vec3 originOffset = pos + sampleNormal * offset;

                Vec3 b;
                const bool result = BarycentricCoordinates(originOffset, sampleNormal,
                                                           triangle.pos[0], triangle.pos[1], triangle.pos[2], &b);
                const Vec2 uv = triangle.uvs[0] * b.x + triangle.uvs[1] * b.y + triangle.uvs[2] * b.z;
                const Vec2Int pixel = { (int)(uv.x * squareSize), (int)(uv.y * squareSize) };
                if (0 <= pixel.x && pixel.x < squareSize && 0 <= pixel.y && pixel.y < squareSize) {
                    const uint32 pixelValue = lightmap.pixels[pixel.y * squareSize + pixel.x];
                    uint32 pixelR = pixelValue & 0xff;
                    uint32 pixelG = (pixelValue >> 8) & 0xff;
                    uint32 pixelB = (pixelValue >> 16) & 0xff;
                    const Vec3 pixelColor = {
                        (float32)pixelR / 255.0f,
                        (float32)pixelG / 255.0f,
                        (float32)pixelB / 255.0f
                    };
                    // TODO adjust color based on material properties, e.g. material should absorb some light
                    float32 weight = MATERIAL_REFLECTANCE;
                    outputColor += weight * sampleContribution * pixelColor;
                }
            }
        }
#endif
    }

    outputColor.r = ClampFloat32(outputColor.r, 0.0f, 1.0f);
    outputColor.g = ClampFloat32(outputColor.g, 0.0f, 1.0f);
    outputColor.b = ClampFloat32(outputColor.b, 0.0f, 1.0f);
    return outputColor;
}
#endif
