#include "delaunay.hh"
#include "triangle.hh"
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

static inline bool is_on_segment(const TriMesh &mesh, const Vh &vh)
{
    for (Eh eh : mesh.ve_range(vh)) if (is_segment(mesh, eh)) return true;
    return false; // check if the vertex lies on a segment
}

static inline bool is_segment(const TriMesh &mesh, const Vh &vh)
{
    return is_sharp(mesh, vh) || is_on_segment(mesh, vh);
}

static inline void set_segment(TriMesh &mesh, const Vh &vh, const bool value)
{
    set_sharp(mesh, vh, value);
}

static void mark_segment_vertices(TriMesh &mesh)
{
    for (Hh hh : mesh.halfedges()) if (is_segment(mesh, hh))
    { set_segment(mesh, mesh.to_vertex_handle(hh), true); }
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
/// Locations
////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

static inline void split(TriMesh &mesh, Hh hh, Vh vh)
{
    split_copy(mesh, mesh.edge_handle(hh), vh);

    // mark as on-segment
    set_segment(mesh, vh, true);

    // One of 2 new diagonal edges can be in the exterior
    // region and hence should be hidden.
    for (Eh eh : mesh.ve_range(vh))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 0)))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 1)))
    { set_hidden(mesh, eh, true); }
}

static inline void split(TriMesh &mesh, Fh fh, Vh vh)
{
    bool hidden = mesh.status(fh).hidden();

    mesh.split_copy(fh, vh);

    // Suppose no vertex will fall into exterior
    // region during smoothing.
    assert(!hidden);

    // The point by chance can be in the exterior
    // region and incident edges should be hidden.
    if (hidden) for (Eh eh : mesh.ve_range(vh))
    { mesh.status(eh).set_hidden(true); }

    if (hidden)
    { mesh.status(vh).set_hidden(true); }
}

static inline Fh search_triangle(const TriMesh &mesh, const Vec2 &u, Fh fh = Fh {})
{
    if (!fh.is_valid()) for (Fh fi : mesh.faces()) { fh = fi; break; }

    //return search_triangle_zigzag(mesh, u, fh);
    return search_triangle_linear(mesh, u, fh);
}

static inline int make_delaunay(TriMesh &mesh, const Vh &vh)
{
    auto delaunifier = make_flipper(mesh, EuclideanDelaunay {});

    const int max_n_flip = (int)mesh.n_edges() * 50;

    for (Hh hh : mesh.voh_range(vh)) if (!mesh.is_boundary(hh))
    {
        delaunifier.enqueue(mesh.edge_handle(hh));
        delaunifier.enqueue(mesh.edge_handle(mesh.next_halfedge_handle(hh)));
    }

    delaunifier.flip_all(max_n_flip);
    delaunifier.clear();

    return 0;
}

static inline bool is_movable(const TriMesh &mesh, const Vh &vh, const Vec2 &u)
{
    const auto u0 = get_xy(mesh, vh);

    for (Hh hh : mesh.voh_range(vh)) if (!mesh.is_boundary(hh))
    {
        const auto u1 = get_xy(mesh, hh);
        const auto u2 = get_xy(mesh, mesh.next_halfedge_handle(hh));
        const auto o0 = orientation(u0, u1, u2);
        const auto o1 = orientation(u , u1, u2);
        if (o0 != o1) return false;
    }

    return true;
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

    // remove the vertex from the topology
    Vh vi = mesh.to_vertex_handle(hhc);
    mesh.collapse(hhc);

    // reuse the vertex in the memory
    mesh.status(vh).set_deleted(false);

    return vi;
}

////////////////////////////////////////////////////////////////
/// Frameworks
////////////////////////////////////////////////////////////////

struct Relocatable
{
    virtual Vec2 operator()(const TriMesh&, Vh) const = 0;
};

static inline int force_move(TriMesh &mesh, const Vh &vh, const Vec2 &u)
{
    // the movement preserves the topology
    if (is_movable(mesh, vh, u))
    {
        set_xy(mesh, vh, u);
        make_delaunay(mesh, vh);
        return 0;
    }

    // the movement screws up the topology, hence
    // we need to re-triangulate around the vertex

    // first, remove the vertex from the topoology
    Vh vi = collapse(mesh, vh);

    // collapsing fails in extreme concave cases
    if (!vi.is_valid()) return 1;

    // improve local triangulation
    make_delaunay(mesh, vi);

    // start with a face nearby
    Fh fh {}; for (Fh fi : mesh.vf_range(vi))
    if (!is_exterior(mesh, fi)) { fh = fi; break;}

    // find the triangle where the point locates
    fh = search_triangle(mesh, u, fh);

    // skip any point out of domain
    if (!fh.is_valid() || is_exterior(mesh, fh)) return 1;

    // at which part of the triangle the point locates
    Hh hh {}; const auto loc = locate(mesh, fh, u, hh);

    // skip any invalid location
    if (loc == TRI_LOC::OUT) return 1;

    // skip any overlapping case
    if ((int)loc & (int)TRI_LOC::VS) return 1;

    // update the position of the vertex
    set_xy(mesh, vh, u);

    // insert back to the topoology
    if (loc == TRI_LOC::IN) split(mesh, fh, vh);
    else                    split(mesh, hh, vh);

    // improve local triangulation
    make_delaunay(mesh, vh);

    return 0;
}

static int forward_smoothing(TriMesh &mesh, const Relocatable &pos_eval, const int max_num_iter)
{
    for (int iter = 0; iter < max_num_iter; ++iter)
    { // begin of iteration

    for (Vh vh : mesh.vertices())
    if (!is_segment(mesh, vh))
    if (!mesh.is_boundary(vh))
    {
        // calculate new positions
        const auto u = pos_eval(mesh, vh);

        // update positions
        force_move(mesh, vh, u);
    }

    } // end of iteration

    mesh.garbage_collection();

    return 0;
}

static int deferred_smoothing(TriMesh &mesh, const Relocatable &pos_eval, const int max_num_iter)
{
    mesh.request_vertex_texcoords2D();

    for (int iter = 0; iter < max_num_iter; ++iter)
    { // begin of iteration

    // calculate new positions
    for (Vh vh : mesh.vertices())
    if (!is_segment(mesh, vh))
    if (!mesh.is_boundary(vh))
    { set_uv(mesh, vh, pos_eval(mesh, vh)); }

    // update positions
    for (Vh vh : mesh.vertices())
    if (!is_segment(mesh, vh))
    if (!mesh.is_boundary(vh))
    { force_move(mesh, vh, get_uv(mesh, vh)); }

    } // end of iteration

    mesh.release_vertex_texcoords2D();
    mesh.garbage_collection();

    return 0;
}

////////////////////////////////////////////////////////////////
/// Algorithms
////////////////////////////////////////////////////////////////

static inline Vec2 delta_laplacian(const TriMesh &mesh, const Vh &vh)
{
    const auto u0 = get_xy(mesh, vh);
    Vec2 du { 0 }; double w {};

    for (Vh vi : mesh.vv_range(vh))
    {
        const auto u = get_xy(mesh, vi);
        double l = norm(u - u0);
        du += (u-u0)*l; w += l; // attractive force
    }

    return du / w;
}

int laplacian_smoothing(TriMesh &mesh, const double step, const int max_num_iter)
{
    struct Relocator : public Relocatable
    {
        Relocator(const double s): lambda(s) {}

        virtual Vec2 operator()(const TriMesh &mesh, Vh vh) const final
        {
            const auto u0 = get_xy(mesh, vh);
            const auto du = delta_laplacian(mesh, vh);
            return u0 + du * lambda;
        }

        const double lambda;
    };

    Relocator relocator(step);

    mark_segment_vertices(mesh);

    return forward_smoothing(mesh, relocator, max_num_iter);
}

static inline Vec2 circumcenter(const TriMesh &mesh, const Fh &fh)
{
    auto hh = mesh.halfedge_handle(fh);
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    return circumcenter(u0, u1, u2);
}

static inline Vec2 centroid(const Vec2 us[], const int n)
{
    Vec2 uc { 0 }; double area {};

    for (int i = 0; i < n; ++i)
    {
        const auto &u0 = us[i];
        const auto &u1 = us[(i+1)%n];
        const double a = cross(u0, u1);
        uc += (u0 + u1) * a;
        area += a;
    }

    return uc / (area*3);
}

static inline Vec2 delta_local_CVT(const TriMesh &mesh, const Vh &vh)
{
    Vec2 us[256]; int n {};

    for (auto fadj : mesh.vf_range(vh)) if (n < 256)
    {
        const auto u = circumcenter(mesh, fadj);
        us[n++] = u;
    }

    if (n > 128) return { 0,0 };

    const auto u0 = get_xy(mesh, vh);
    const auto u1 = centroid(us, n);

    return u1 - u0;
}

int local_CVT_smoothing(TriMesh &mesh, const double step, const int max_num_iter)
{
    struct Relocator : public Relocatable
    {
        Relocator(const double s): lambda(s) {}

        virtual Vec2 operator()(const TriMesh &mesh, Vh vh) const final
        {
            const auto u0 = get_xy(mesh, vh);
            const auto du = delta_local_CVT(mesh, vh);
            return u0 + du * lambda;
        }

        const double lambda;
    };

    Relocator relocator(step);

    mark_segment_vertices(mesh);

    return deferred_smoothing(mesh, relocator, max_num_iter);
}

static inline Vec2 delta_local_ODT(const TriMesh &mesh, const Vh &vh)
{
    const auto u0 = get_xy(mesh, vh);
    Vec2 du { 0 }; double area {};

    for (Hh hh : mesh.voh_range(vh))
    {
        const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh));
        const auto u2 = get_xy(mesh, mesh.to_vertex_handle(mesh.next_halfedge_handle(hh)));
        const auto d1 = u1 - u0;
        const auto d2 = u2 - u0;
        const auto dt = u2 - u1;
        const auto n = Vec2 { dt[1], -dt[0] } * 0.5;
        const double s = (dot(d1, d1) + dot(d2, d2));
        const double a = fabs(cross(d1, d2));
        du += n*s;
        area += a;
    }

    return du / area;
}

int local_ODT_smoothing(TriMesh &mesh, const double step, const int max_num_iter)
{
    struct Relocator : public Relocatable
    {
        Relocator(const double s): lambda(s) {}

        virtual Vec2 operator()(const TriMesh &mesh, Vh vh) const final
        {
            const auto u0 = get_xy(mesh, vh);
            const auto du = delta_local_ODT(mesh, vh);
            return u0 + du * lambda;
        }

        const double lambda;
    };

    Relocator relocator(step);

    return deferred_smoothing(mesh, relocator, max_num_iter);
}
