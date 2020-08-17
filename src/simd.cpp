// SIMD helpers ------------------------------------------------------------------------

const __m256 ZERO_8 = _mm256_setzero_ps();
const __m256 ONE_8 = _mm256_set1_ps(1.0f);

bool AllZero_8(__m256 x)
{
    return _mm256_testc_ps(ZERO_8, x) == 1;
}

bool AnyNonZero_8(__m256 x)
{
    return _mm256_testc_ps(ZERO_8, x) != 1;
}

inline __m256 operator+(__m256 a, __m256 b)
{
    return _mm256_add_ps(a, b);
}
inline __m256 operator-(__m256 a, __m256 b)
{
    return _mm256_sub_ps(a, b);
}
inline __m256 operator*(__m256 a, __m256 b)
{
    return _mm256_mul_ps(a, b);
}

inline __m256 operator-(__m256 x)
{
    return ZERO_8 - x;
}

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

void StoreVec3_8(Vec3_8 v, Vec3* mem)
{
    // TODO probably slow? but I really don't wanna do a shuffle-y implementation
    alignas(32) float32 x[8];
    alignas(32) float32 y[8];
    alignas(32) float32 z[8];

    _mm256_store_ps(x, v.x);
    _mm256_store_ps(y, v.y);
    _mm256_store_ps(z, v.z);

    for (uint32 k = 0; k < 8; k++) {
        mem[k].x = x[k];
        mem[k].y = y[k];
        mem[k].z = z[k];
    }
}

inline Vec3_8 operator+(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = v1.x + v2.x,
        .y = v1.y + v2.y,
        .z = v1.z + v2.z,
    };
}
inline Vec3_8 operator-(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = v1.x - v2.x,
        .y = v1.y - v2.y,
        .z = v1.z - v2.z,
    };
}

inline Vec3_8 operator*(Vec3_8 v, __m256 s)
{
    return Vec3_8 {
        .x = v.x * s,
        .y = v.y * s,
        .z = v.z * s,
    };
}
inline Vec3_8 operator/(Vec3_8 v, __m256 s)
{
    return Vec3_8 {
        .x = v.x * s,
        .y = v.y * s,
        .z = v.z * s,
    };
}

inline Vec3_8 operator-(Vec3_8 v)
{
    return { -v.x, -v.y, -v.z };
}

Vec3_8 Multiply_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = v1.x * v2.x,
        .y = v1.y * v2.y,
        .z = v1.z * v2.z,
    };
}

__m256 MagSq_8(Vec3_8 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
__m256 Mag_8(Vec3_8 v)
{
    return _mm256_sqrt_ps(MagSq_8(v));
}

Vec3_8 Normalize_8(Vec3_8 v)
{
    const __m256 mag = Mag_8(v);
    return v / mag;
}

__m256 Dot_8(Vec3_8 v1, Vec3_8 v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec3_8 Cross_8(Vec3_8 v1, Vec3_8 v2)
{
    return Vec3_8 {
        .x = v1.y * v2.z - v1.z * v2.y,
        .y = v1.z * v2.x - v1.x * v2.z,
        .z = v1.x * v2.y - v1.y * v2.x,
    };
}

Vec3_8 Reciprocal_8(Vec3_8 v)
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
        .w = _mm256_set1_ps(q.w),
    };
}

Quat_8 Multiply_8(Quat_8 q1, Quat_8 q2)
{
    return Quat_8 {
        .x = q1.w * q2.x - q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
        .y = q1.w * q2.y - q1.y * q2.w + q1.z * q2.x - q1.x * q2.z,
        .z = q1.w * q2.z - q1.z * q2.w + q1.x * q2.y - q1.y * q2.x,
        .w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
    };
}

Quat_8 Inverse_8(Quat_8 q)
{
    return Quat_8 {
        .x = ZERO_8 - q.x,
        .y = ZERO_8 - q.y,
        .z = ZERO_8 - q.z,
        .w = q.w
    };
}

Vec3_8 Multiply_8(Quat_8 q, Vec3_8 v)
{
    const Quat_8 vQuat = { v.x, v.y, v.z, ZERO_8 };
    const Quat_8 qv = Multiply_8(q, vQuat);

    const Quat_8 qInv = Inverse_8(q);
    const Quat_8 qvqInv = Multiply_8(qv, qInv);

    return Vec3_8 { qvqInv.x, qvqInv.y, qvqInv.z };
}

__m256 RayPlaneIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, Vec3_8 planeOrigin8, Vec3_8 planeNormal8, __m256* t8)
{
    const __m256 dotDirNormal8 = Dot_8(rayDir8, planeNormal8);
    // Set mask when dot is non-zero (otherwise, ray direction is perpendicular to plane normal, so no intersection)
    const __m256 result8 = _mm256_cmp_ps(dotDirNormal8, ZERO_8, _CMP_NEQ_OQ);

    const __m256 invDotDirNormal8 = _mm256_rcp_ps(dotDirNormal8);
    *t8 = _mm256_mul_ps(Dot_8(planeOrigin8 - rayOrigin8, planeNormal8), invDotDirNormal8);
    return result8;
}

__m256 RayAABBIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDirInv8, Box aabb, __m256* tMinOut, __m256* tMaxOut)
{
    const Vec3_8 aabbMin8 = Set1Vec3_8(aabb.min);
    const Vec3_8 aabbMax8 = Set1Vec3_8(aabb.max);

    __m256 tMin = _mm256_set1_ps(-INFINITY);
    __m256 tMax = _mm256_set1_ps(INFINITY);

    const __m256 tX1 = (aabbMin8.x - rayOrigin8.x) * rayDirInv8.x;
    const __m256 tX2 = (aabbMax8.x - rayOrigin8.x) * rayDirInv8.x;
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tX1, tX2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tX1, tX2));

    const __m256 tY1 = (aabbMin8.y - rayOrigin8.y) * rayDirInv8.y;
    const __m256 tY2 = (aabbMax8.y - rayOrigin8.y) * rayDirInv8.y;
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tY1, tY2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tY1, tY2));

    const __m256 tZ1 = (aabbMin8.z - rayOrigin8.z) * rayDirInv8.z;
    const __m256 tZ2 = (aabbMax8.z - rayOrigin8.z) * rayDirInv8.z;
    tMin = _mm256_max_ps(tMin, _mm256_min_ps(tZ1, tZ2));
    tMax = _mm256_min_ps(tMax, _mm256_max_ps(tZ1, tZ2));

    *tMinOut = tMin;
    *tMaxOut = tMax;
    // NOTE: doing an ordered (O) and non-signaling (Q) compare for greater than or equals here
    // This means that if there's a NaN value, the comparison will return false, but no exception will be triggered
    return _mm256_cmp_ps(tMax, tMin, _CMP_GE_OQ);
}

__m256 RayTriangleIntersection_8(Vec3_8 rayOrigin8, Vec3_8 rayDir8, Vec3 a, Vec3 b, Vec3 c, __m256* t8)
{
    const float32 epsilon = 0.000001f;
    const __m256 epsilon8 = _mm256_set1_ps(epsilon);
    const __m256 negEpsilon8 = _mm256_set1_ps(-epsilon);

    const Vec3_8 a8 = Set1Vec3_8(a);
    const Vec3_8 b8 = Set1Vec3_8(b);
    const Vec3_8 c8 = Set1Vec3_8(c);

    const Vec3_8 ab8 = b8 - a8;
    const Vec3_8 ac8 = c8 - a8;
    const Vec3_8 h8 = Cross_8(rayDir8, ac8);
    const __m256 x8 = Dot_8(ab8, h8);
    // Result mask is set when x < -EPSILON || x > EPSILON
    __m256 result8 = _mm256_or_ps(_mm256_cmp_ps(x8, negEpsilon8, _CMP_LT_OQ), _mm256_cmp_ps(x8, epsilon8, _CMP_GT_OQ));

    const __m256 f8 = _mm256_rcp_ps(x8);
    const Vec3_8 s8 = rayOrigin8 - a8;
    const __m256 u8 = f8 * Dot_8(s8, h8);
    // Result mask is set when 0.0f <= u <= 1.0f
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(ZERO_8, u8, _CMP_LE_OQ));
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(u8, ONE_8, _CMP_LE_OQ));

    const Vec3_8 q8 = Cross_8(s8, ab8);
    const __m256 v8 = f8 * Dot_8(rayDir8, q8);
    // Result mask is set when 0.0f <= v && u + v <= 1.0f
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(ZERO_8, v8, _CMP_LE_OQ));
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(u8 + v8, ONE_8, _CMP_LE_OQ));

    *t8 = f8 * Dot_8(ac8, q8);
    // Result mask is set when t8 >= 0.0f (otherwise, intersection point is behind the ray origin)
    // NOTE if t is 0, intersection is a line (I think)
    result8 = _mm256_and_ps(result8, _mm256_cmp_ps(*t8, ZERO_8, _CMP_GE_OQ));
    return result8;
}

// -------------------------------------------------------------------------------------
