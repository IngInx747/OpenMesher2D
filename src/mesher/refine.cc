#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "delaunay.hh"
#include "triangle.hh"
#include "segment.hh"
#include "mesh.hh"
#include "topology.hh"
#include "search.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// Attributes
////////////////////////////////////////////////////////////////

static inline bool is_exterior(const TriMesh &mesh, const Fh &fh)
{
    return is_hidden(mesh, fh);
}

static inline bool is_exterior(const TriMesh &mesh, const Eh &eh)
{
    return is_hidden(mesh, eh);
}

static inline bool is_exterior(const TriMesh &mesh, const Hh &hh)
{
    return mesh.is_boundary(hh) || is_hidden(mesh, mesh.face_handle(hh));
}

static inline bool is_segment(const TriMesh &mesh, const Eh &eh)
{
    return is_sharp(mesh, eh);
}

static inline bool is_segment(const TriMesh &mesh, const Hh &hh)
{
    return is_sharp(mesh, mesh.edge_handle(hh));
}

static inline bool is_segment(const TriMesh &mesh, const Vh &vh)
{
    for (Eh eh : mesh.ve_range(vh)) if (is_segment(mesh, eh)) return true;
    return false; // is the vertex on a segment
}

static inline bool is_endian(const TriMesh &mesh, const Vh &vh)
{
    return is_sharp(mesh, vh); // is the vertex from an original segment
}

struct OnSegment
{
    inline bool operator()(const TriMesh &mesh, const Hh &hh) const
    { return is_segment(mesh, hh); }
};

static inline Vh segment_head(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = hh; OnSegment pred {};
    while (!is_endian(mesh, mesh.from_vertex_handle(hi)))
    { if ((hi = prev(mesh, pred, hi)) == hh) break; }
    return mesh.from_vertex_handle(hi);
}

static inline Vh segment_tail(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = hh; OnSegment pred {};
    while (!is_endian(mesh, mesh.to_vertex_handle(hi)))
    { if ((hi = next(mesh, pred, hi)) == hh) break; }
    return mesh.to_vertex_handle(hi);
}

static void mark_endians(TriMesh &mesh)
{
    for (Hh hh : mesh.halfedges()) if (is_sharp(mesh, mesh.edge_handle(hh)))
    { set_sharp(mesh, mesh.to_vertex_handle(hh), true); }
}

////////////////////////////////////////////////////////////////
/// Delaunay
////////////////////////////////////////////////////////////////

//   1  
//  / \ 
// 2---0
//  \ / 
//   3  
static inline bool is_delaunay(const TriMesh &mesh, const Eh &eh)
{
    Hh hh0 = mesh.halfedge_handle(eh, 0);
    Hh hh1 = mesh.halfedge_handle(eh, 1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(mesh.next_halfedge_handle(hh0)));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u3 = get_xy(mesh, mesh.to_vertex_handle(mesh.next_halfedge_handle(hh1)));
    return fuzzy_delaunay(u0, u1, u2, u3);
}

struct EuclideanDelaunay
{
    inline bool operator()(const TriMesh &mesh, const Eh &eh) const
    { return is_sharp(mesh, eh) || is_delaunay(mesh, eh); } // If true, do not flip
};

////////////////////////////////////////////////////////////////
/// Containers
////////////////////////////////////////////////////////////////

template <class T>
struct unique_vector
{
    inline void push_back(const T &v)
    { if (std::find(vs.begin(), vs.end(), v) == vs.end()) { vs.push_back(v); } }

    inline void pop_back()
    { vs.pop_back(); }

    inline void clear()
    { vs.clear(); }

    inline const std::vector<T> &vector() const
    { return vs; }

protected:

    std::vector<T> vs;
};

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

static inline bool is_inside(const TriMesh &mesh, const Fh &fh, const Vec2 &u)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    const auto loc = exact_locate(u0, u1, u2, u);
    return loc != TRI_LOC::OUT;
}

static inline TRI_LOC locate(const TriMesh &mesh, const Fh &fh, const Vec2 &u, Hh &hh)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto loc = exact_locate(u0, u1, u2, u);
    if (loc == TRI_LOC::E0) { hh = hh0; }
    if (loc == TRI_LOC::E1) { hh = hh1; }
    if (loc == TRI_LOC::E2) { hh = hh2; }
    if (loc == TRI_LOC::V0) { hh = hh1; }
    if (loc == TRI_LOC::V1) { hh = hh2; }
    if (loc == TRI_LOC::V2) { hh = hh0; }
    return loc;
}

static inline void split(TriMesh &mesh, Hh hh, Vh vh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);

    mesh.split_edge_copy(mesh.edge_handle(hh), vh);

    // OpenMesh copies property to all 4 new edges even 2 of
    // which are not sharp. Unsharp the 2 edges accordingly.
    for (Hh hh : mesh.voh_range(vh))
    if (mesh.to_vertex_handle(hh) != vh0)
    if (mesh.to_vertex_handle(hh) != vh1)
    { set_sharp(mesh, mesh.edge_handle(hh), false); }

    // One of 2 new diagonal edges can be in the exterior
    // region and hence should be hidden.
    for (Eh eh : mesh.ve_range(vh))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 0)))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 1)))
    { set_hidden(mesh, eh, true); }
}

static inline void split(TriMesh &mesh, Fh fh, Vh vh)
{
    mesh.split_copy(fh, vh);
}

static inline bool is_collapsable(TriMesh &mesh, const Hh &hhc)
{
    Vh vhc = mesh.from_vertex_handle(hhc);
    Vh vh0 = mesh.to_vertex_handle  (hhc);
    const auto u0 = get_xy(mesh, vh0);

    for (Hh hh : mesh.voh_range(vhc)) // check n-2 triangles
    {
        Vh vh1 = mesh.to_vertex_handle(hh);
        Vh vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
        if (vh1 == vh0 || vh2 == vh0) continue;
        const auto u1 = get_xy(mesh, vh1);
        const auto u2 = get_xy(mesh, vh2);
        const int r = orientation(u0, u1, u2);
        if (r <= 0) return false;
    }

    return true;
}

static inline Vh collapse(TriMesh &mesh, Vh vh)
{
    assert(!mesh.is_boundary(vh));

    Hh hhc {}; // to collapse

    for (Hh hh : mesh.voh_range(vh)) if (is_collapsable(mesh, hh))
    { hhc = hh; break; }

    if (!hhc.is_valid()) return Vh {};

    Vh vi = mesh.to_vertex_handle(hhc);
    mesh.collapse(hhc);
    return vi;
}

static Fh search_primitive(const TriMesh &mesh, const Vec2 &u1, const Fh &fh0, Hh &hh)
{
    if (is_inside(mesh, fh0, u1)) return fh0;

    PrimitivePlow pp(mesh);

    init(pp, fh0, u1); // setup the plow

    const int max_n_iter = (int)mesh.n_edges(); int n_iter {};

    for (pp.next(); n_iter < max_n_iter; pp.next(), ++n_iter)
    {
        if (pp.status() == PLOW_STATUS::EDGE) // record segments
        {
            if (is_segment(mesh, pp.halfedge_handle()))
            { hh = pp.halfedge_handle(); break; }
        }

        Fh fh {}; // check if the target is reached

        if (pp.status() == PLOW_STATUS::EDGE)
        {
            fh = mesh.opposite_face_handle(pp.halfedge_handle());
        }
        else if (pp.status() == PLOW_STATUS::VERT)
        {
            fh = mesh.face_handle(mesh.halfedge_handle(pp.vertex_handle()));
        }
        if (fh.is_valid() && is_inside(mesh, fh, u1))
        {
            return fh;
        }

        if (pp.status() == PLOW_STATUS::MISS) // searching lost in vain, for some reasons
        {
            break;
        }
    }

    return Fh {};
}

////////////////////////////////////////////////////////////////
/// Primitive union
////////////////////////////////////////////////////////////////

struct Primitive // either a triangle or a segment
{
    Vh vh0, vh1, vh2;
};

template<>
struct std::hash<Primitive>
{
    size_t operator()(const Primitive &primitive) const noexcept
    {
        size_t h0 = hash<Vh>{}(primitive.vh0);
        size_t h1 = hash<Vh>{}(primitive.vh1);
        size_t h2 = hash<Vh>{}(primitive.vh2);
        return (((h2 << 1) ^ h1) << 1) ^ h0;
    }
};

template<>
struct std::equal_to<Primitive>
{
    bool operator()(const Primitive &lhs, const Primitive &rhs) const noexcept
    {
        return
            lhs.vh0 == rhs.vh0 &&
            lhs.vh1 == rhs.vh1 &&
            lhs.vh2 == rhs.vh2;
        }
};

static inline bool is_triangle(const Primitive &primitive)
{
    return primitive.vh2.is_valid();
}

static inline bool is_segment(const Primitive &primitive)
{
    return !primitive.vh2.is_valid();
}

static inline Fh get_triangle(const TriMesh &mesh, const Primitive &primitive)
{
    if (is_deleted(mesh, primitive.vh0) ||
        is_deleted(mesh, primitive.vh1) ||
        is_deleted(mesh, primitive.vh2))
        return Fh {};

    Hh hh = mesh.find_halfedge(primitive.vh0, primitive.vh1);
    if (!hh.is_valid()) return Fh {};

    Vh vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
    if (vh2 != primitive.vh2) return Fh {};

    return mesh.face_handle(hh);
}

static inline Hh get_segment(const TriMesh &mesh, const Primitive &primitive)
{
    if (is_deleted(mesh, primitive.vh0) ||
        is_deleted(mesh, primitive.vh1))
        return Hh {};

    return mesh.find_halfedge(primitive.vh0, primitive.vh1);
}

static inline Primitive make_primitive(const TriMesh &mesh, const Fh &fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    Vh vh0 = mesh.to_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle(hi);
    Vh vh2 = mesh.to_vertex_handle(hj);
    return { vh0, vh1, vh2 };
}

static inline Primitive make_primitive(const TriMesh &mesh, const Hh &hh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);
    return { vh0, vh1, Vh {} };
}

////////////////////////////////////////////////////////////////
/// Constants
////////////////////////////////////////////////////////////////

static inline double apex_squared_cosine(const double min_angle)
{
    const double cs0 = cos(min_angle); // cos(t)
    const double cs1 = cs0*cs0*2. - 1; // cos(pi - 2t)
    return cs1*cs1;
}

static inline double offcenter_height(const double cs)
{
    // formular adapted from Triangle lib [J.R.Shewchuk]
    return .475 * sqrt((1. + cs)/(1. - cs));
}

////////////////////////////////////////////////////////////////
/// Segment encroachment
////////////////////////////////////////////////////////////////

static inline bool is_encroached(const TriMesh &mesh, const Hh &hh)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    const auto d0 = u0 - u2;
    const auto d1 = u1 - u2;
    return (dot(d0, d1) < 0.); // a.b < 0
}

static inline bool is_encroached(const TriMesh &mesh, const Hh &hh, const double cs2)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    const auto d0 = u0 - u2;
    const auto d1 = u1 - u2;
    const double d00 = dot(d0, d0); // a^2
    const double d11 = dot(d1, d1); // b^2
    const double d01 = dot(d0, d1); // a.b
    return (d01 < 0.) && (d01*d01 >= d00*d11*cs2); // a.b < 0 and |a.b| > (a*b)|cos(t0)|
}

struct Encroachment
{
    Encroachment(const double min_angle): cs2(apex_squared_cosine(min_angle)) {}

    inline bool operator()(const TriMesh &mesh, const Hh &hh) const // If true, split it
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, cs2); }

    const double cs2; // (cos(pi - 2t))^2
};

////////////////////////////////////////////////////////////////
/// Triangle quality
////////////////////////////////////////////////////////////////

static inline Hh the_shortest_edge(const TriMesh &mesh, const Fh &fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hi));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hj));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh));
    const double d0 = dot(u2 - u1, u2 - u1); // len2(hh)
    const double d1 = dot(u0 - u2, u0 - u2); // len2(hi)
    const double d2 = dot(u1 - u0, u1 - u0); // len2(hj)
    return (d1<d0) ? (d1<d2 ? hi : hj) : (d0<d2 ? hh : hj);
}

static inline double the_shortest_len2(const TriMesh &mesh, const Fh &fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hi));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hj));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh));
    const double d0 = dot(u2 - u1, u2 - u1); // len2(hh)
    const double d1 = dot(u0 - u2, u0 - u2); // len2(hi)
    const double d2 = dot(u1 - u0, u1 - u0); // len2(hj)
    return (d1<d0) ? (d1<d2 ? d1 : d2) : (d0<d2 ? d0 : d2);
}

static inline bool is_subtending_input_angle(const TriMesh &mesh, const Hh &hh)
{
    OnSegment pred {};

    // try getting two segments between which the base is
    Hh hh0 = next(mesh, pred, hh);
    Hh hh1 = prev(mesh, pred, hh);

    // one or both ends of the base not lying on any segment
    if (!is_segment(mesh, hh0) || !is_segment(mesh, hh1)) return false;

    Vh vh01 = segment_head(mesh, hh0);
    Vh vh00 = segment_tail(mesh, hh0);
    Vh vh11 = segment_head(mesh, hh1);
    Vh vh10 = segment_tail(mesh, hh1);
    Vh vhc {}; // the common end

    if (vh01 == vh10) vhc = vh01;
    if (vh00 == vh11) vhc = vh00;

    // two segment not sharing a common end
    if (!vhc.is_valid()) return false;

    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto uc = get_xy(mesh, vhc);
    const auto d0 = u0 - uc;
    const auto d1 = u1 - uc;
    const double rd = dot(d0, d0) / dot(d1, d1);

    // check if two ends lie on a circle centering the common end
    return fabs(rd - 1.0) < 1e-3;
}

static inline bool is_bad_triangle(
    const TriMesh &mesh, const Fh &fh,
    const double max_cos2,
    const double max_len2,
    const double max_twoa)
{
    const Hh hh[3] {
        mesh.halfedge_handle(fh),
        mesh.next_halfedge_handle(mesh.halfedge_handle(fh)),
        mesh.prev_halfedge_handle(mesh.halfedge_handle(fh)),
    };

    //      1      //
    //  0  / \  2  //
    //    /   \    //
    //   2-----0   //
    //      1      //

    const Vec2 u[3] {
        get_xy(mesh, mesh.to_vertex_handle(hh[1])),
        get_xy(mesh, mesh.to_vertex_handle(hh[2])),
        get_xy(mesh, mesh.to_vertex_handle(hh[0])),
    };

    const Vec2 d[3] {
        u[2] - u[1],
        u[0] - u[2],
        u[1] - u[0],
    };

    const double twoa = fabs(cross(-d[1], d[2])); // area*2

    // check area upper bound
    if (twoa > max_twoa) return true;

    const double dd[3] {
        dot(d[0], d[0]), // len^2 of edge 0
        dot(d[1], d[1]), // len^2 of edge 1
        dot(d[2], d[2]), // len^2 of edge 2
    };

    // check length upper bound
    if (dd[0] > max_len2 ||
        dd[1] > max_len2 ||
        dd[2] > max_len2 ) return true;

    const double dp[3] {
        dot(-d[1], d[2]), // dot product at corner 0
        dot(-d[2], d[0]), // dot product at corner 1
        dot(-d[0], d[1]), // dot product at corner 2
    };

    // use the shortest edge as base
    double ddb = dd[0]; int i = 0;
    if (ddb > dd[1]) { ddb = dd[1]; i = 1; }
    if (ddb > dd[2]) { ddb = dd[2]; i = 2; }
    int j = (i + 1) % 3; // leg 1
    int k = (i + 2) % 3; // leg 2

    const double cs2 = (dp[i]*dp[i])/(dd[j]*dd[k]); // cos^2 of the smallest corner angle

    // check angle lower bound
    if (cs2 > max_cos2) {

    // skip if the smallest edge subtends an input angle
    if (!is_segment(mesh, hh[i]) && is_subtending_input_angle(mesh, hh[i])) return false;

    return true; }

    // the triangle passes all tests, hence no refinement is needed
    return false;
}

struct BadTriangle
{
    BadTriangle(const double min_angle, const double max_length, const double max_area):
    max_cos2(cos(min_angle)*cos(min_angle)),
    max_len2(max_length*max_length),
    max_twoa(max_area*2.0) {}

    inline bool operator()(const TriMesh &mesh, const Fh &fh) const // If true, split it
    { return is_bad_triangle(mesh, fh, max_cos2, max_len2, max_twoa); }

    const double max_cos2; // lower bound of corner angle
    const double max_len2; // upper bound of edge length
    const double max_twoa; // upper bound of triangle area
};

////////////////////////////////////////////////////////////////
/// Refinement
////////////////////////////////////////////////////////////////

static inline Vec2 circumcenter(const TriMesh &mesh, const Fh &fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hi));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hj));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh));
    return circumcenter(u0, u1, u2);
}

static inline Vec2 offcenter(const TriMesh &mesh, const Hh &hh, const double h)
{
    // off-center [A. Ungor]
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto du = u1 - u0;
    const Vec2 dv { -du[1], du[0] };
    return u0 + du*.5 + dv*h;
}

static inline Vec2 splitting_position(const TriMesh &mesh, const Fh &fh, const double h)
{
    Hh hh = the_shortest_edge(mesh, fh);
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto uc = circumcenter(mesh, fh);
    const auto uh = offcenter(mesh, hh, h);
    const auto um = (u0 + u1) * .5;
    const double dc = dot(uc - um, uc - um);
    const double dh = dot(uh - um, uh - um);
    return (dc>dh) ? uh : uc;
}

static inline Vec2 splitting_position(const TriMesh &mesh, const Hh &hh)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    double t = 0.5;

    /// The concentric circle algorithm is adapted from Triangle
    ///   lib to handle small angle in the input. [J.R.Shewchuk]

    Hh hi = mesh.opposite_halfedge_handle(hh);

    bool is_acute_v0 = 
        !is_exterior(mesh, hh) && is_segment(mesh, mesh.prev_halfedge_handle(hh)) ||
        !is_exterior(mesh, hi) && is_segment(mesh, mesh.next_halfedge_handle(hi)) ;

    bool is_acute_v1 = 
        !is_exterior(mesh, hh) && is_segment(mesh, mesh.next_halfedge_handle(hh)) ||
        !is_exterior(mesh, hi) && is_segment(mesh, mesh.prev_halfedge_handle(hi)) ;

    if (is_acute_v0 || is_acute_v1)
    {
        const double l = norm(u1 - u0);
        double ep = 1; // nearest power of 2

        // find the ratio that splits the segment most evenly
        while (ep*3.0 < l) { ep *= 2.0; }
        while (ep*1.5 > l) { ep *= 0.5; }

        t = ep / l;
        t = is_acute_v1 ? 1 - t : t;
    }

    return u0*(1-t) + u1*t;
}

static int refine_segments(TriMesh &mesh, const Encroachment &encroached)
{
    int n_new_vertices {};

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    // Before enqueueing encroached segments, free vertices in
    // the diametral circle of the segment should be deleted.
    // [TODO]

    std::deque<Primitive> segments {}; // primitives to split

    for (Eh eh : mesh.edges()) { if (is_segment(mesh, eh))
    for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh))
    { segments.push_back(make_primitive(mesh, hh)); } } }

    while (!segments.empty())
    {
        auto primitive = segments.front(); segments.pop_front();

        if (is_segment(primitive))
        {
            // restore the segment
            Hh ho = get_segment(mesh, primitive);

            // the segment was gone during refining
            if (!ho.is_valid()) continue;

            // find a position for splitting
            const auto u = splitting_position(mesh, ho);

            // split the segment
            Vh vo = mesh.new_vertex({ 0,0,0 });
            set_xy(mesh, vo, u);
            split(mesh, ho, vo);

            // affected edges
            unique_vector<Eh> eas {};

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record affected edges while testing Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
            { if (!is_exterior(mesh, eh)) eas.push_back(eh); }

            // check encroachment of two new segments
            for (Eh eh : mesh.ve_range(vo)) { if (is_segment(mesh, eh)) eas.push_back(eh); }

            // check encroachment of affected segments
            for (Eh eh : eas.vector()) { if (!is_deleted(mesh, eh)) if (is_segment(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh))
            { segments.push_back(make_primitive(mesh, hh)); } } }

            ++n_new_vertices;
        }
    }

    return 0;
}

static int refine_interior(TriMesh &mesh, const BadTriangle &bad_triangle, const Encroachment &encroached)
{
    using Priority = float;

    struct Entry { Primitive primitive; Priority priority; };

    const auto less = [](const Entry &a, const Entry &b)
    { return a.priority < b.priority; };

    // Key is the length of the shortest edge. The shorter,
    //   the higher the order of the triangle in the queue.
    // Use negative value to place the shorter ones on top.
    const auto calc_priority = [](const TriMesh &mesh, const Fh &fh)
    { return static_cast<Priority>(-the_shortest_len2(mesh, fh)); };

    // constant of off-center [A. Ungor]
    const double off_h = offcenter_height(sqrt(bad_triangle.max_cos2));

    int n_new_vertices {};

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    // segments to split
    std::deque<Primitive> segments {};

    // triangles to refine
    std::priority_queue<Entry, std::deque<Entry>, decltype(less)> triangles(less);

    // check quality of all triangles at initialization
    for (Fh fh : mesh.faces()) if (bad_triangle(mesh, fh))
    { triangles.push({ make_primitive(mesh, fh), calc_priority(mesh, fh) }); }

    while (true)
    {
        if (!segments.empty())
        {
            auto primitive = segments.front(); segments.pop_front();

            // restore the segment
            Hh ho = get_segment(mesh, primitive);

            // the segment was gone during refining
            if (!ho.is_valid()) continue;

            // affected edges
            unique_vector<Hh> es {};

            // remove free vertices within the diametral circle of the segment before splitting
            for (const auto &segment : { primitive, Primitive { primitive.vh1, primitive.vh0, Vh {} } }) while (true)
            {
                Hh hh = get_segment(mesh, segment);
                if (!hh.is_valid()) break;

                Vh vh = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
                if (is_segment(mesh, vh)) break;

                // stop when the segment is no longer encroached within the diametral circle
                if (!is_encroached(mesh, hh)) break;

                // remove the apex vertex
                Vh vc = collapse(mesh, vh);
                if (!vc.is_valid()) break;

                // edges that are potentially non-Delaunay
                std::vector<Eh> ehs {}; for (Eh eh : mesh.ve_range(vc)) { ehs.push_back(eh); }
                delaunifier.reset(); delaunifier.enqueue(ehs.data(), (int)ehs.size());

                // record edges while maintaining Delaunayhood
                for (Eh eh = delaunifier.flip(); eh.is_valid(); eh = delaunifier.flip()) {
                for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh)) es.push_back(hh); }
            }

            // the segment was gone during last step (won't happen)
            if (!(ho = get_segment(mesh, primitive)).is_valid()) continue;

            // find a position for splitting
            const auto u = splitting_position(mesh, ho);

            // split the segment
            Vh vo = mesh.new_vertex({ 0,0,0 });
            set_xy(mesh, vo, u);
            split(mesh, ho, vo);

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record affected edges while maintaining Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next()) {
            for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh)) es.push_back(hh); }

            // record two new segments
            for (Eh eh : mesh.ve_range(vo)) { if (is_segment(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh)) es.push_back(hh); }

            // check encroachment of affected segments
            for (Hh hh : es.vector()) { if (!is_deleted(mesh, hh) && is_segment(mesh, hh) && encroached(mesh, hh))
            { segments.push_back(make_primitive(mesh, hh)); } }

            // check quality of affected triangles
            for (Hh hh : es.vector()) { if (!is_deleted(mesh, hh) && !is_exterior(mesh, hh)) { Fh fh = mesh.face_handle(hh);
            if (bad_triangle(mesh, fh)) { triangles.push({ make_primitive(mesh, fh), calc_priority(mesh, fh) }); } } }

            ++n_new_vertices;
        }

        else if (!triangles.empty())
        {            
            auto primitive = triangles.top().primitive; triangles.pop();

            // restore the triangle
            Fh fr = get_triangle(mesh, primitive);

            // the triangle was gone during refining
            if (!fr.is_valid()) continue;

            // find the circumcenter of the triangle
            const auto u = splitting_position(mesh, fr, off_h);

            // search for the primitive at which the circumcenter locates
            Hh sh; Fh fo = search_primitive(mesh, u, fr, sh);

            if (sh.is_valid()) // some segments are in the way hence encroached
            {
                segments.push_back(make_primitive(mesh, sh));
                continue; // abort splitting
            }

            // skip any point out of domain (won't happen as some segments must
            // be encroached and splitting must have been aborted.)
            if (!fo.is_valid() || is_exterior(mesh, fo)) continue;

            // at which part of the triangle the point locates
            Hh ho {}; const auto loc = locate(mesh, fo, u, ho);

            // skip any vertex with invalid location (won't happen)
            if (loc == TRI_LOC::OUT) continue;

            // skip any duplicated vertices (hardly happens, or Delaunayhood is
            // not maintained as the circumcircle contains a point in the mesh.)
            if ((int)loc & (int)TRI_LOC::VS) continue;

            // allocate a new vertex
            Vh vo = mesh.new_vertex({ 0,0,0 });
            set_xy(mesh, vo, u);

            // insert the vertex into the triangle or onto the edge
            if (loc == TRI_LOC::IN) split(mesh, fo, vo);
            else                    split(mesh, ho, vo);

            // affected edges
            unique_vector<Hh> es {};

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record affected edges while maintaining Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next()) {
            for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh)) es.push_back(hh); }

            // check encroachment over affected segments
            std::vector<Hh> hhs {}; for (Hh hh : es.vector()) {
            if (!is_deleted(mesh, hh) && is_segment(mesh, hh) && encroached(mesh, hh)) hhs.push_back(hh); }

            // record encroached segments
            if (!hhs.empty()) { for (Hh hh : hhs) segments.push_back(make_primitive(mesh, hh)); }

            // undo vertex insertion if some segments are encroached
            while (!hhs.empty())
            {
                Vh vc = collapse(mesh, vo);
                if (!vc.is_valid()) break;

                // In the extreme concave cases will the undo-insertion fail.
                // One can ease the concavity by unflipping around its 1-ring.
                // [TODO]

                // edges that are potentially non-Delaunay
                std::vector<Eh> ehs {}; for (Eh eh : mesh.ve_range(vc)) { ehs.push_back(eh); }
                delaunifier.reset(); delaunifier.enqueue(ehs.data(), (int)ehs.size());

                // record affected edges while maintaining Delaunayhood
                for (Eh eh = delaunifier.flip(); eh.is_valid(); eh = delaunifier.flip()) {
                for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh)) es.push_back(hh); }

                break;
            }

            // check quality of affected triangles
            for (Hh hh : es.vector()) { if (!is_deleted(mesh, hh) && !is_exterior(mesh, hh)) { Fh fh = mesh.face_handle(hh);
            if (bad_triangle(mesh, fh)) { triangles.push({ make_primitive(mesh, fh), calc_priority(mesh, fh) }); } } }

            ++n_new_vertices;
        }

        else break;
    }

    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int refine(TriMesh &mesh, const double min_angle, const double max_length, const double max_area)
{
    int err {};

    Encroachment encroached(min_angle);

    BadTriangle bad_triangle(min_angle, max_length, max_area);

    mark_endians(mesh);

    refine_segments(mesh, encroached);

    refine_interior(mesh, bad_triangle, encroached);

    mesh.garbage_collection();

    return err;
}