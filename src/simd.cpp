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
