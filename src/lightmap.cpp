#include "lightmap.h"

#include <intrin.h>
#include <stb_image_write.h>
#include <time.h>

// SIMD helpers ------------------------------------------------------------------------

struct Vec3_8
{
    __m256 x, y, z;
};

struct Quat_8
{
    __m256 x, y, z, w;
};

Vec3_8 Set1Vec3_8(Vec3 v)
{
    return Vec3_8 { _mm256_set1_ps(v.x), _mm256_set1_ps(v.y), _mm256_set1_ps(v.z) };
}

Vec3_8 SetVec3_8(const StaticArray<Vec3, 8>& vs)
{
    return Vec3_8 {
        .x = _mm256_set_ps(vs[7].x, vs[6].x, vs[5].x, vs[4].x, vs[3].x, vs[2].x, vs[1].x, vs[0].x),
        .y = _mm256_set_ps(vs[7].y, vs[6].y, vs[5].y, vs[4].y, vs[3].y, vs[2].y, vs[1].y, vs[0].y),
        .z = _mm256_set_ps(vs[7].z, vs[6].z, vs[5].z, vs[4].z, vs[3].z, vs[2].z, vs[1].z, vs[0].z),
    };
}

Vec3_8 Add_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = _mm256_add_ps(v1.x, v2.x),
        .y = _mm256_add_ps(v1.y, v2.y),
        .z = _mm256_add_ps(v1.z, v2.z),
    };
}

Vec3_8 Subtract_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = _mm256_sub_ps(v1.x, v2.x),
        .y = _mm256_sub_ps(v1.y, v2.y),
        .z = _mm256_sub_ps(v1.z, v2.z),
    };
}

Vec3_8 Multiply_8(Vec3_8 v, __m256 s)
{
    return Vec3_8 {
        .x = _mm256_mul_ps(v.x, s),
        .y = _mm256_mul_ps(v.y, s),
        .z = _mm256_mul_ps(v.z, s),
    };
}

Vec3_8 Multiply_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = _mm256_mul_ps(v1.x, v2.x),
        .y = _mm256_mul_ps(v1.y, v2.y),
        .z = _mm256_mul_ps(v1.z, v2.z),
    };
}

Vec3_8 Divide_8(Vec3_8 v, __m256 s)
{
    const __m256 sInv = _mm256_rcp_ps(s);
    return Vec3_8 {
        .x = _mm256_mul_ps(v.x, sInv),
        .y = _mm256_mul_ps(v.y, sInv),
        .z = _mm256_mul_ps(v.z, sInv),
    };
}

__m256 MagSq_8(Vec3_8 v)
{
    return _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v.x, v.x), _mm256_mul_ps(v.y, v.y)), _mm256_mul_ps(v.z, v.z));
}
__m256 Mag_8(Vec3_8 v)
{
    return _mm256_sqrt_ps(MagSq_8(v));
}

Vec3_8 Normalize_8(Vec3_8 v)
{
    const __m256 mag = Mag_8(v);
    return Divide_8(v, mag);
}

__m256 Dot_8(Vec3_8 v1, Vec3_8 v2)
{
    return _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(v1.x, v2.x), _mm256_mul_ps(v1.y, v2.y)),
                         _mm256_mul_ps(v1.z, v2.z));
}

Vec3_8 Cross_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = _mm256_sub_ps(_mm256_mul_ps(v1.y, v2.z), _mm256_mul_ps(v1.z, v2.y)),
        .y = _mm256_sub_ps(_mm256_mul_ps(v1.z, v2.x), _mm256_mul_ps(v1.x, v2.z)),
        .z = _mm256_sub_ps(_mm256_mul_ps(v1.x, v2.y), _mm256_mul_ps(v1.y, v2.x)),
    };
}

Vec3_8 Inverse_8(Vec3_8 v)
{
    return Vec3_8 {
        .x = _mm256_rcp_ps(v.x),
        .y = _mm256_rcp_ps(v.y),
        .z = _mm256_rcp_ps(v.z),
    };
}

Quat_8 Set1Quat_8(Quat q)
{
    return Quat_8 {
        .x = _mm256_set1_ps(q.x),
        .y = _mm256_set1_ps(q.y),
        .z = _mm256_set1_ps(q.z),
        .w = _mm256_set1_ps(q.w)
    };
}

Quat_8 Multiply_8(Quat_8 q1, Quat_8 q2)
{
    return Quat_8 {
        .x = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(q1.w, q2.x), _mm256_mul_ps(q1.x, q2.w)),
                           _mm256_sub_ps(_mm256_mul_ps(q1.y, q2.z), _mm256_mul_ps(q1.z, q2.y))),
        .y = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(q1.w, q2.y), _mm256_mul_ps(q1.y, q2.w)),
                           _mm256_sub_ps(_mm256_mul_ps(q1.z, q2.x), _mm256_mul_ps(q1.x, q2.z))),
        .z = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(q1.w, q2.z), _mm256_mul_ps(q1.z, q2.w)),
                           _mm256_sub_ps(_mm256_mul_ps(q1.x, q2.y), _mm256_mul_ps(q1.y, q2.x))),
        .w = _mm256_sub_ps(_mm256_sub_ps(_mm256_mul_ps(q1.w, q2.w), _mm256_mul_ps(q1.x, q2.x)),
                           _mm256_add_ps(_mm256_mul_ps(q1.y, q2.y), _mm256_mul_ps(q1.z, q2.z))),
    };
}

Quat_8 Inverse_8(Quat_8 q)
{
    const __m256 zero8 = _mm256_setzero_ps();
    return Quat_8 {
        .x = _mm256_sub_ps(zero8, q.x),
        .y = _mm256_sub_ps(zero8, q.y),
        .z = _mm256_sub_ps(zero8, q.z),
        .w = q.w
    };
}

Vec3_8 Multiply_8(Quat_8 q, Vec3_8 v)
{
    const __m256 zero8 = _mm256_setzero_ps();
    Quat_8 vQuat = { v.x, v.y, v.z, zero8 };
    Quat_8 qv = Multiply_8(q, vQuat);

    Quat_8 qInv = Inverse_8(q);
    Quat_8 qvqInv = Multiply_8(qv, qInv);

    return Vec3_8 { qvqInv.x, qvqInv.y, qvqInv.z };
}

__m256 RayPlaneIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, Vec3_8 planeOrigin8, Vec3_8 planeNormal8, __m256* t8)
{
    const __m256 zero8 = _mm256_setzero_ps();

    const __m256 dotDirNormal8 = Dot_8(rayDir8, planeNormal8);
    // Set mask when dot is non-zero (otherwise, ray direction is perpendicular to plane normal, so no intersection)
    const __m256 result8 = _mm256_cmp_ps(dotDirNormal8, zero8, _CMP_NEQ_OQ);

    const __m256 invDotDirNormal8 = _mm256_rcp_ps(dotDirNormal8);
    *t8 = _mm256_mul_ps(Dot_8(Subtract_8(planeOrigin8, rayOrigin8), planeNormal8), invDotDirNormal8);
    return result8;
}

__m256 RayAxisAlignedBoxIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDirInv8, Vec3 boxMin, Vec3 boxMax)
{
    const Vec3_8 boxMin8 = Set1Vec3_8(boxMin);
    const Vec3_8 boxMax8 = Set1Vec3_8(boxMax);

    __m256 tMin = _mm256_set1_ps(-INFINITY);
    __m256 tMax = _mm256_set1_ps(INFINITY);

    const __m256 tX1 = _mm256_mul_ps(_mm256_sub_ps(boxMin8.x, rayOrigin8.x), rayDirInv8.x);
    const __m256 tX2 = _mm256_mul_ps(_mm256_sub_ps(boxMax8.x, rayOrigin8.x), rayDirInv8.x);
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tX1, tX2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tX1, tX2));

    const __m256 tY1 = _mm256_mul_ps(_mm256_sub_ps(boxMin8.y, rayOrigin8.y), rayDirInv8.y);
    const __m256 tY2 = _mm256_mul_ps(_mm256_sub_ps(boxMax8.y, rayOrigin8.y), rayDirInv8.y);
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tY1, tY2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tY1, tY2));

    const __m256 tZ1 = _mm256_mul_ps(_mm256_sub_ps(boxMin8.z, rayOrigin8.z), rayDirInv8.z);
    const __m256 tZ2 = _mm256_mul_ps(_mm256_sub_ps(boxMax8.z, rayOrigin8.z), rayDirInv8.z);
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tZ1, tZ2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tZ1, tZ2));

    // NOTE: doing an ordered (O) and non-signaling (Q) compare for greater than or equals here
    // This means that if there's a NaN value, the comparison will return false, but no exception will be triggered
    __m256 result8 = _mm256_cmp_ps(tMax, tMin, _CMP_GE_OQ);
    return result8;
}

__m256 RayTriangleIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, Vec3 a, Vec3 b, Vec3 c, __m256* t8)
{
    const __m256 zero8 = _mm256_setzero_ps();
    const __m256 one8 = _mm256_set1_ps(1.0f);
    const float32 epsilon = 0.000001f;
    const __m256 epsilon8 = _mm256_set1_ps(epsilon);
    const __m256 negEpsilon8 = _mm256_set1_ps(-epsilon);

    const Vec3_8 a8 = Set1Vec3_8(a);
    const Vec3_8 b8 = Set1Vec3_8(b);
    const Vec3_8 c8 = Set1Vec3_8(c);

    const Vec3_8 ab8 = Subtract_8(b8, a8);
    const Vec3_8 ac8 = Subtract_8(c8, a8);
    const Vec3_8 h8 = Cross_8(rayDir8, ac8);
    const __m256 x8 = Dot_8(ab8, h8);
    // Result mask is set when x < -EPSILON || x > EPSILON
    __m256 result8 = _mm256_or_ps(_mm256_cmp_ps(x8, negEpsilon8, _CMP_LT_OQ), _mm256_cmp_ps(x8, epsilon8, _CMP_GT_OQ));

    const __m256 f8 = _mm256_rcp_ps(x8);
    const Vec3_8 s8 = Subtract_8(rayOrigin8, a8);
    const __m256 u8 = _mm256_mul_ps(f8, Dot_8(s8, h8));
    // Result mask is set when 0.0f <= u <= 1.0f
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(zero8, u8, _CMP_LE_OQ));
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(u8, one8, _CMP_LE_OQ));

    const Vec3_8 q8 = Cross_8(s8, ab8);
    const __m256 v8 = _mm256_mul_ps(f8, Dot_8(rayDir8, q8));
    // Result mask is set when 0.0f <= v && u + v <= 1.0f
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(zero8, v8, _CMP_LE_OQ));
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(_mm256_add_ps(u8, v8), one8, _CMP_LE_OQ));

    *t8 = _mm256_mul_ps(f8, Dot_8(ac8, q8));
    // Result mask is set when t8 >= 0.0f (otherwise, intersection point is behind the ray origin)
    // NOTE if t is 0, intersection is a line (I think)
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(*t8, zero8, _CMP_GE_OQ));
    return result8;
}

// -------------------------------------------------------------------------------------

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

struct LightRect
{
    Vec3 origin;
    Vec3 width;
    Vec3 height;
    Vec3 color;
    float32 intensity;
};

struct Lightmap
{
    uint32 squareSize;
    uint32* pixels;
};

struct RaycastTriangle
{
    Vec3 pos[3];
    // TODO move these somewhere else. I think only pos is used in the hot raycasting loops
    Vec3 color[3];
    Vec3 normal;
    Vec2 uvs[3];
};

struct RaycastMesh
{
    Vec3 min, max;
    Lightmap lightmap;
    Array<RaycastTriangle> triangles;
};

struct RaycastGeometry
{
    Array<RaycastMesh> meshes;
};

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
                mesh.triangles[j].uvs[k] = t.v[k].uv;
                mesh.triangles[j].color[k] = vertexColor;
            }
            mesh.triangles[j].normal = normal;
        }
        for (uint32 j = 0; j < obj.models[i].quads.size; j++) {
            const uint32 ind = obj.models[i].triangles.size + j * 2;
            const ObjQuad& q = obj.models[i].quads[j];
            const Vec3 normal = CalculateTriangleUnitNormal(q.v[0].pos, q.v[1].pos, q.v[2].pos);

            for (int k = 0; k < 3; k++) {
                mesh.triangles[ind].pos[k] = q.v[k].pos;
                mesh.triangles[ind].uvs[k] = q.v[k].uv;
                mesh.triangles[ind].color[k] = vertexColor;
            }
            mesh.triangles[ind].normal = normal;

            for (int k = 0; k < 3; k++) {
                const uint32 quadInd = (k + 2) % 4;
                mesh.triangles[ind + 1].pos[k] = q.v[quadInd].pos;
                mesh.triangles[ind + 1].uvs[k] = q.v[quadInd].uv;
                mesh.triangles[ind + 1].color[k] = vertexColor;
            }
            mesh.triangles[ind + 1].normal = normal;
        }

        float32 surfaceArea = 0.0f;
        for (uint32 j = 0; j < mesh.triangles.size; j++) {
            const RaycastTriangle& t = mesh.triangles[j];
            surfaceArea += TriangleArea(t.pos[0], t.pos[1], t.pos[2]);
        }

        // Calculate AABB
        mesh.min = Vec3::one * 1e8;
        mesh.max = -Vec3::one * 1e8;
        for (uint32 j = 0; j < mesh.triangles.size; j++) {
            for (int k = 0; k < 3; k++) {
                const Vec3 v = mesh.triangles[j].pos[k];
                for (int e = 0; e < 3; e++) {
                    mesh.min.e[e] = MinFloat32(mesh.min.e[e], v.e[e]);
                    mesh.max.e[e] = MaxFloat32(mesh.max.e[e], v.e[e]);
                }
            }
        }

        // Allocate lightmap
        const uint32 size = (uint32)(sqrt(surfaceArea) * RESOLUTION_PER_WORLD_UNIT);
        const uint32 squareSize = RoundUpToPowerOfTwo(MinInt(size, 1024));
        mesh.lightmap.squareSize = squareSize;
        mesh.lightmap.pixels = allocator->New<uint32>(squareSize * squareSize);
        if (mesh.lightmap.pixels == nullptr) {
            LOG_ERROR("Failed to allocate %dx%d pixels for lightmap %lu\n", squareSize, squareSize, i);
            geometry.meshes.data = nullptr;
            return geometry;
        }

        MemSet(mesh.lightmap.pixels, 0, squareSize * squareSize * sizeof(uint32));
    }

    return geometry;
}

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
                                                                      mesh.min, mesh.max);
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
    }

    outputColor.r = ClampFloat32(outputColor.r, 0.0f, 1.0f);
    outputColor.g = ClampFloat32(outputColor.g, 0.0f, 1.0f);
    outputColor.b = ClampFloat32(outputColor.b, 0.0f, 1.0f);
    return outputColor;
}

struct WorkLightmapRasterizeRowCommon
{
    Array<SampleGroup> hemisphereSampleGroups;
    const RaycastGeometry* geometry;
    uint32 meshInd;
    Lightmap* lightmap;
};

struct WorkLightmapRasterizeRow
{
    const WorkLightmapRasterizeRowCommon* common;
    uint32 triangleInd;
    int minPixelX, maxPixelX, pixelY;
};

void LightmapRasterizeRow(const RaycastGeometry& geometry, uint32 meshInd, uint32 triangleInd,
                          int minPixelX, int maxPixelX, int pixelY, Array<SampleGroup> hemisphereSampleGroups,
                          Lightmap* lightmap)
{
    const RaycastTriangle& triangle = geometry.meshes[meshInd].triangles[triangleInd];
    const uint32 squareSize = lightmap->squareSize;
    const float32 uvY = (float32)pixelY / squareSize;

    for (int x = minPixelX; x < maxPixelX; x++) {
        const float32 uvX = (float32)x / squareSize;
        const Vec3 bC = BarycentricCoordinates(Vec2 { uvX, uvY }, triangle.uvs[0], triangle.uvs[1], triangle.uvs[2]);
        const Vec3 p = triangle.pos[0] * bC.x + triangle.pos[1] * bC.y + triangle.pos[2] * bC.z;

        const Vec3 raycastColor = RaycastColor(hemisphereSampleGroups, p, triangle.normal, geometry);
        uint8 r = (uint8)(raycastColor.r * 255.0f);
        uint8 g = (uint8)(raycastColor.g * 255.0f);
        uint8 b = (uint8)(raycastColor.b * 255.0f);
        uint8 a = 0xff;
        lightmap->pixels[pixelY * squareSize + x] = (a << 24) + (b << 16) + (g << 8) + r;
    }
}

uint32 bounce_ = 0;

void ThreadLightmapRasterizeRow(AppWorkQueue* queue, void* data)
{
    WorkLightmapRasterizeRow* workData = (WorkLightmapRasterizeRow*)data;

    const uint32 remaining = queue->entriesTotal - queue->entriesComplete;
    if (remaining % 100 == 0) {
        const RaycastGeometry& geometry = *workData->common->geometry;
        const int32 meshInd = workData->common->meshInd;
        LOG_INFO("%d rows in queue | bounce %lu, mesh %lu, triangle %lu/%lu, row %d (%d pixels)\n",
                 remaining, bounce_, meshInd, workData->triangleInd, geometry.meshes[meshInd].triangles.size,
                 workData->pixelY, workData->maxPixelX - workData->minPixelX);
    }

    LightmapRasterizeRow(*workData->common->geometry, workData->common->meshInd, workData->triangleInd,
                         workData->minPixelX, workData->maxPixelX, workData->pixelY, 
                         workData->common->hemisphereSampleGroups, workData->common->lightmap);
}

#if 0
internal bool CalculateLightmapForMesh(const RaycastGeometry& geometry, uint32 meshInd, AppWorkQueue* queue,
                                       LinearAllocator* allocator, Lightmap* lightmap)
{
    const int LIGHTMAP_PIXEL_MARGIN = 1;

    Array<SampleGroup> hemisphereSampleGroups = allocator->NewArray<SampleGroup>(NUM_HEMISPHERE_SAMPLE_GROUPS);
    if (hemisphereSampleGroups.data == nullptr) {
        LOG_ERROR("Failed to allocate hemisphere sample groups\n");
        return false;
    }
    if (!GenerateHemisphereSampleGroups(hemisphereSampleGroups, allocator)) {
        LOG_ERROR("Failed to generate hemisphere sample groups\n");
        return false;
    }

    const WorkLightmapRasterizeRowCommon workCommon = {
        .hemisphereSampleGroups = hemisphereSampleGroups,
        .geometry = &geometry,
        .meshInd = meshInd,
        .lightmap = lightmap
    };
    auto allocatorState = allocator->SaveState();

    const RaycastMesh& mesh = geometry.meshes[meshInd];
    const uint32 squareSize = mesh.lightmap.squareSize;

    for (uint32 i = 0; i < mesh.triangles.size; i++) {
#if RESTRICT_LIGHTING && RESTRICT_WALL
        if (i / 2 != PLANE_TO_LIGHT) continue;
#endif
        const RaycastTriangle& triangle = mesh.triangles[i];
        const Vec2 minUv = {
            MinFloat32(triangle.uvs[0].x, MinFloat32(triangle.uvs[1].x, triangle.uvs[2].x)),
            MinFloat32(triangle.uvs[0].y, MinFloat32(triangle.uvs[1].y, triangle.uvs[2].y))
        };
        const Vec2 maxUv = {
            MaxFloat32(triangle.uvs[0].x, MaxFloat32(triangle.uvs[1].x, triangle.uvs[2].x)),
            MaxFloat32(triangle.uvs[0].y, MaxFloat32(triangle.uvs[1].y, triangle.uvs[2].y))
        };
        const int minPixelY = MaxInt((int)(minUv.y * squareSize) - LIGHTMAP_PIXEL_MARGIN, 0);
        const int maxPixelY = MinInt((int)(maxUv.y * squareSize) + LIGHTMAP_PIXEL_MARGIN, squareSize);
        for (int y = minPixelY; y < maxPixelY; y++) {
            const float32 uvY = (float32)y / squareSize;
            const float32 t0 = (uvY - triangle.uvs[0].y) / (triangle.uvs[1].y - triangle.uvs[0].y);
            const float32 x0 = triangle.uvs[0].x + (triangle.uvs[1].x - triangle.uvs[0].x) * t0;
            const float32 t1 = (uvY - triangle.uvs[1].y) / (triangle.uvs[2].y - triangle.uvs[1].y);
            const float32 x1 = triangle.uvs[1].x + (triangle.uvs[2].x - triangle.uvs[1].x) * t1;
            const float32 t2 = (uvY - triangle.uvs[2].y) / (triangle.uvs[0].y - triangle.uvs[2].y);
            const float32 x2 = triangle.uvs[2].x + (triangle.uvs[0].x - triangle.uvs[2].x) * t2;
            float32 minX = maxUv.x;
            if (x0 >= minUv.x && x0 < minX) {
                minX = x0;
            }
            if (x1 >= minUv.x && x1 < minX) {
                minX = x1;
            }
            if (x2 >= minUv.x && x2 < minX) {
                minX = x2;
            }
            float32 maxX = minUv.x;
            if (x0 <= maxUv.x && x0 > maxX) {
                maxX = x0;
            }
            if (x1 <= maxUv.x && x1 > maxX) {
                maxX = x1;
            }
            if (x2 <= maxUv.x && x2 > maxX) {
                maxX = x2;
            }
            int minPixelX = MaxInt((int)(minX * squareSize) - LIGHTMAP_PIXEL_MARGIN, 0);
            int maxPixelX = MinInt((int)(maxX * squareSize) + LIGHTMAP_PIXEL_MARGIN, squareSize);
            int numPixels = maxPixelX - minPixelX;
            if (numPixels < 0) {
                // This happens because we give Y pixels an extra margin
                // TODO yeah, maybe this is OK
                minPixelX = (minPixelX + maxPixelX) / 2;
                maxPixelX = minPixelX + 1;
                numPixels = 1;
            }

            WorkLightmapRasterizeRow* work = allocator->New<WorkLightmapRasterizeRow>();
            if (work == nullptr) {
                CompleteAllWork(queue);
                allocator->LoadState(allocatorState);
                work = allocator->New<WorkLightmapRasterizeRow>();
                if (work == nullptr) {
                    LOG_ERROR("Failed to allocate work entry after allocator flush\n");
                    return false;
                }
            }
            work->common = &workCommon;
            work->triangleInd = i;
            work->minPixelX = minPixelX;
            work->maxPixelX = maxPixelX;
            work->pixelY = y;
            if (!TryAddWork(queue, ThreadLightmapRasterizeRow, work)) {
                CompleteAllWork(queue);
                DEBUG_ASSERT(TryAddWork(queue, ThreadLightmapRasterizeRow, work));
            }
        }
    }

    CompleteAllWork(queue);

    return true;
}
#endif

bool LightMeshVertices(const RaycastGeometry& geometry, uint32 meshInd, LinearAllocator* allocator,
                       Array<Vec3> vertexColors)
{
    Array<SampleGroup> hemisphereSampleGroups = allocator->NewArray<SampleGroup>(NUM_HEMISPHERE_SAMPLE_GROUPS);
    if (hemisphereSampleGroups.data == nullptr) {
        LOG_ERROR("Failed to allocate hemisphere sample groups\n");
        return false;
    }
    if (!GenerateHemisphereSampleGroups(hemisphereSampleGroups, allocator)) {
        LOG_ERROR("Failed to generate hemisphere sample groups\n");
        return false;
    }


    for (uint32 i = 0; i < geometry.meshes[meshInd].triangles.size; i++) {
        const RaycastTriangle& t = geometry.meshes[meshInd].triangles[i];
        for (int j = 0; j < 3; j++) {
            const Vec3 dir = t.normal;
            const Vec3 raycastColor = RaycastColor(hemisphereSampleGroups, t.pos[j], dir, geometry);
            vertexColors[i * 3 + j] = raycastColor;
        }
    }

    return true;
}

bool GenerateLightmaps(const LoadObjResult& obj, uint32 bounces, AppWorkQueue* queue, LinearAllocator* allocator,
                       const_string lightmapDirPath)
{
    RaycastGeometry geometry = CreateRaycastGeometry(obj, allocator);
    if (geometry.meshes.data == nullptr) {
        LOG_ERROR("Failed to construct raycast geometry from obj\n");
        return false;
    }

    uint32 totalTriangles = 0;
    for (uint32 i = 0; i < geometry.meshes.size; i++) {
        totalTriangles += geometry.meshes[i].triangles.size;
    }

    LOG_INFO("Generating lightmaps for %lu meshes, %lu total triangles, %lu bounces\n",
             geometry.meshes.size, totalTriangles, bounces);

    DebugTimer lightmapTimer = StartDebugTimer();

    for (uint32 b = 0; b < bounces; b++) {
        LOG_INFO("Bounce %lu\n", b);
        bounce_ = b; // NOTE for logging purposes only

        ALLOCATOR_SCOPE_RESET(*allocator);

        Array<Lightmap> lightmaps = allocator->NewArray<Lightmap>(geometry.meshes.size);
        if (lightmaps.data == nullptr) {
            LOG_ERROR("Failed to allocate lightmaps array in bounce %lu\n", b);
            return false;
        }

        Array<Array<Vec3>> meshVertexColors = allocator->NewArray<Array<Vec3>>(geometry.meshes.size);
        if (meshVertexColors.data == nullptr) {
            LOG_ERROR("Failed to allocate vertex color arrays in bounce %lu\n", b);
            return false;
        }

        for (uint32 i = 0; i < geometry.meshes.size; i++) {
            const uint32 squareSize = geometry.meshes[i].lightmap.squareSize;
            lightmaps[i].squareSize = squareSize,
            lightmaps[i].pixels = allocator->New<uint32>(squareSize * squareSize);
            if (lightmaps[i].pixels == nullptr) {
                LOG_ERROR("Failed to allocate pixels for lightmap, bounce %lu, mesh %lu\n", b, i);
                return false;
            }
            MemSet(lightmaps[i].pixels, 0, squareSize * squareSize * sizeof(uint32));

            meshVertexColors[i] = allocator->NewArray<Vec3>(geometry.meshes[i].triangles.size * 3);
            if (meshVertexColors[i].data == nullptr) {
                LOG_ERROR("Failed to allocate vertex colors, bounce %lu, mesh %lu\n", b, i);
                return false;
            }
        }

#if RESTRICT_LIGHTING
        for (uint32 m = 0, i = MODELS_TO_LIGHT[m]; m < C_ARRAY_LENGTH(MODELS_TO_LIGHT); i = MODELS_TO_LIGHT[++m])
#else
        for (uint32 i = 0; i < geometry.meshes.size; i++)
#endif
        {
            LOG_INFO("Lighting mesh %lu\n", i);

            ALLOCATOR_SCOPE_RESET(*allocator);

            UNREFERENCED_PARAMETER(queue);
#if 0
            // Calculate lightmap for mesh and save to file
            if (!CalculateLightmapForMesh(geometry, i, queue, allocator, &lightmaps[i])) {
                LOG_ERROR("Failed to compute lightmap for mesh %lu, bounce %lu\n", i, b);
                return false;
            }
            const char* lightmapFilePath = ToCString(AllocPrintf(allocator, "%.*s/%d.png",
                                                                 lightmapDirPath.size, lightmapDirPath.data, i),
                                                     allocator);
            if (!stbi_write_png(lightmapFilePath, lightmaps[i].squareSize, lightmaps[i].squareSize, 4,
                                lightmaps[i].pixels, 0)) {
                LOG_ERROR("Failed to save lightmap file to %s for mesh %lu, bounce %lu\n",
                          lightmapFilePath, i, b);
                return false;
            }
#endif

            // Calculate vertex light for mesh and save to file
            if (!LightMeshVertices(geometry, i, allocator, meshVertexColors[i])) {
                LOG_ERROR("Failed to light vertices for mesh %lu, bounce %lu\n", i, b);
                return false;
            }
            const Array<uint8> vertexColorData = {
                .size = meshVertexColors[i].size * sizeof(Vec3),
                .data = (uint8*)meshVertexColors[i].data
            };
            string verticesFilePath = AllocPrintf(allocator, "%.*s/%d.v", lightmapDirPath.size, lightmapDirPath.data, i);
            if (!WriteFile(verticesFilePath, vertexColorData, false)) {
                LOG_ERROR("Failed to write light vertices to %.*s for mesh %lu, bounce %lu\n",
                          verticesFilePath.size, verticesFilePath.data, i, b);
                return false;
            }
        }

        if (b != bounces - 1) {
            // Copy calculated lightmaps to RaycastGeometry for the next bounce
            for (uint32 i = 0; i < geometry.meshes.size; i++) {
                const uint32 squareSize = lightmaps[i].squareSize;
                MemCopy(geometry.meshes[i].lightmap.pixels, lightmaps[i].pixels, squareSize * squareSize * sizeof(uint32));
            }
        }
    }

    StopAndPrintDebugTimer(&lightmapTimer);

    return true;
}
