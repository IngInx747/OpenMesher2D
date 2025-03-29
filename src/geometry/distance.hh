#ifndef EUCLIDEAN_DISTANCE_HH
#define EUCLIDEAN_DISTANCE_HH

#include "math_def.hh"
#include "vector_n.hh"

inline double hausdorff_distance(const VecN<Vec2, 2> &s, const Vec2 &p)
{
    const auto &a = s[0];
    const auto &b = s[1];
    const auto d0 = p - a;
    const auto d1 = p - b;
    const auto e = (b - a).normalized();
    return dot(d0, e) < 0 ? norm(d0)
        :  dot(d1,-e) < 0 ? norm(d1)
        :  fabs(cross(d0, e));
}

inline double hausdorff_distance(const VecN<Vec3, 2> &s, const Vec3 &p)
{
    const auto &a = s[0];
    const auto &b = s[1];
    const auto d0 = p - a;
    const auto d1 = p - b;
    const auto e = (b - a).normalized();
    return dot(d0, e) < 0 ? norm(d0)
        :  dot(d1,-e) < 0 ? norm(d1)
        :  norm(cross(d0, e));
}

inline double hausdorff_distance(const VecN<Vec2, 3> &t, const Vec2 &p)
{
    const auto &a = t[0];
    const auto &b = t[1];
    const auto &c = t[2];
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const auto na = (Vec2 { bc[1], -bc[0] }).normalized();
    const auto nb = (Vec2 { ca[1], -ca[0] }).normalized();
    const auto nc = (Vec2 { ab[1], -ab[0] }).normalized();
    if (dot(ab, ap) <= 0 && dot(ca, ap) >= 0) return norm(ap);
    if (dot(bc, bp) <= 0 && dot(ab, bp) >= 0) return norm(bp);
    if (dot(ca, cp) <= 0 && dot(bc, cp) >= 0) return norm(cp);
    if (dot(nc, ap) >= 0 && dot(ab, ap) >= 0 && dot(ab, bp) <= 0) return dot(nc, ap);
    if (dot(na, bp) >= 0 && dot(bc, bp) >= 0 && dot(bc, cp) <= 0) return dot(na, bp);
    if (dot(nb, cp) >= 0 && dot(ca, cp) >= 0 && dot(ca, ap) <= 0) return dot(nb, cp);
    return 0; // inside triangle
}

inline double hausdorff_distance(const VecN<Vec3, 3> &t, const Vec3 &p)
{
    const auto &a = t[0];
    const auto &b = t[1];
    const auto &c = t[2];
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const auto n =  cross(ab,bc).normalized();
    const auto na = cross(bc, n).normalized();
    const auto nb = cross(ca, n).normalized();
    const auto nc = cross(ab, n).normalized();
    if (dot(ab, ap) <= 0 && dot(ca, ap) >= 0) return norm(ap);
    if (dot(bc, bp) <= 0 && dot(ab, bp) >= 0) return norm(bp);
    if (dot(ca, cp) <= 0 && dot(bc, cp) >= 0) return norm(cp);
    if (dot(nc, ap) >= 0 && dot(ab, ap) >= 0 && dot(ab, bp) <= 0) return norm(cross(ab, ap)) / norm(ab);
    if (dot(na, bp) >= 0 && dot(bc, bp) >= 0 && dot(bc, cp) <= 0) return norm(cross(bc, bp)) / norm(bc);
    if (dot(nb, cp) >= 0 && dot(ca, cp) >= 0 && dot(ca, ap) <= 0) return norm(cross(ca, cp)) / norm(ca);
    return fabs(dot(ap, n)); // inside triangle
}

#endif