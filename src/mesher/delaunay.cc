#include "delaunay.hh"
#include "mesh.hh"
#include "topology.hh"

using namespace OpenMesh;

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

int make_delaunay(TriMesh &mesh)
{
    auto delaunifier = make_flipper(mesh, EuclideanDelaunay {});

    const int max_n_flip = (int)mesh.n_edges() * 50;

    delaunifier.enqueue_all(); int n_flip {};

#if 1 // non-stop flipping until Delaunayhood is met
    n_flip = delaunifier.flip_all(max_n_flip);

#elif 1 // test Delaunayhood while checking each flipped edge
    for (Eh eh = delaunifier.flip(); eh.is_valid() && n_flip < max_n_flip; eh = delaunifier.flip(), ++n_flip)
    {} // do something each flipping

#else // test Delaunayhood while checking each encountered edge
    for (Eh eh = delaunifier.next(); eh.is_valid() && n_flip < max_n_flip; eh = delaunifier.next(), ++n_flip)
    {} // do something each encountering

#endif

    delaunifier.clear();

    return 0;
}
