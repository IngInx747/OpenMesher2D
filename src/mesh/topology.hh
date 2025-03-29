#ifndef MESH_TOPOLOGY_HH
#define MESH_TOPOLOGY_HH

#include <queue>
#include <unordered_map>
#include <unordered_set>

#include "mesh.hh"

////////////////////////////////////////////////////////////////
/// Delaunay
////////////////////////////////////////////////////////////////

namespace OpenMesh
{

inline bool is_rotate_ok(const PolyConnectivity &mesh, const HalfedgeHandle &hh)
{
    // boundary edges cannot be rotated
    if (mesh.is_boundary(mesh.edge_handle(hh))) return false;

    const auto oh = mesh.opposite_halfedge_handle(hh);

    // check if the rotated edge is already present in the mesh
    const auto vh0 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
    const auto vh1 = mesh.to_vertex_handle(mesh.next_halfedge_handle(oh));

    if (vh0 == vh1) return false; // this is generally a bad sign !!!

    for (PolyConnectivity::ConstVertexVertexIter vvi(mesh, vh0); vvi.is_valid(); ++vvi)
        if (*vvi == vh1) return false;

    return true;
}

inline void rotate(PolyConnectivity &mesh, const HalfedgeHandle &hh)
{
    // CAUTION : Rotating a halfedge may result in
    // a non-manifold mesh, hence check for yourself
    // whether this operation is allowed or not!
    assert(!mesh.is_boundary(mesh.edge_handle(hh)));
    assert(is_rotate_ok(mesh, hh));

    auto ha0 = hh;
    auto hb0 = mesh.opposite_halfedge_handle(hh);

    auto ha1 = mesh.next_halfedge_handle(ha0);
    auto ha2 = mesh.next_halfedge_handle(ha1);
    auto ha3 = mesh.prev_halfedge_handle(ha0);

    auto hb1 = mesh.next_halfedge_handle(hb0);
    auto hb2 = mesh.next_halfedge_handle(hb1);
    auto hb3 = mesh.prev_halfedge_handle(hb0);

    auto va0 = mesh.to_vertex_handle(ha0);
    auto vb0 = mesh.to_vertex_handle(hb0);

    auto va1 = mesh.to_vertex_handle(ha1);
    auto vb1 = mesh.to_vertex_handle(hb1);

    auto fa = mesh.face_handle(ha0);
    auto fb = mesh.face_handle(hb0);

    mesh.set_vertex_handle(ha0, va1);
    mesh.set_vertex_handle(hb0, vb1);

    mesh.set_next_halfedge_handle(ha0, ha2);
    mesh.set_next_halfedge_handle(hb1, ha0);
    mesh.set_next_halfedge_handle(ha3, hb1);

    mesh.set_next_halfedge_handle(hb0, hb2);
    mesh.set_next_halfedge_handle(ha1, hb0);
    mesh.set_next_halfedge_handle(hb3, ha1);

    mesh.set_face_handle(ha1, fb);
    mesh.set_face_handle(hb1, fa);

    mesh.set_halfedge_handle(fa, ha0);
    mesh.set_halfedge_handle(fb, hb0);

    if (mesh.halfedge_handle(va0) == hb0)
        mesh.set_halfedge_handle(va0, ha1);
    if (mesh.halfedge_handle(vb0) == ha0)
        mesh.set_halfedge_handle(vb0, hb1);
}

inline void flip(PolyConnectivity &mesh, const EdgeHandle &eh)
{
    rotate(mesh, mesh.halfedge_handle(eh, 0));
}

inline void flip(TriConnectivity &mesh, const EdgeHandle &eh)
{
    mesh.flip(eh);
}

} // namespace OpenMesh

template <class MeshT, class DelaunayT>
class Flipper
{
public:

    Flipper(MeshT &mesh, const DelaunayT &is_delaunay)
    : mesh(mesh), is_delaunay(is_delaunay) {}

    /// Clear all enqueued edges
    void clear();

    /// Enqueue all edges in the mesh
    void enqueue_all();

    /// Enqueue an edge
    void enqueue(const Eh&);

    /// Enqueue edges
    void enqueue(const Eh[], const int);

    /// Enqueue edges with a rule
    template <class PredicateT> void enqueue(const PredicateT&);

    /// Test Delaunayhood (and flip) one edge at a time.
    Eh next();

    /// Keep testing Delaunayhood until an edge is flipped.
    /// The flipped edge or a null handle will be returned.
    Eh flip();

    // Keep testing Delaunayhood and flipping non-Delaunay
    // edges until no more edges to flip or it exceeds the 
    // maxinum number of flippings.
    int flip_all(const int max_n_flip);

    size_t size() const { return frontier.size(); }

protected:

    void enqueue_adjacent(const Eh&);

    inline bool is_enqueued(const Eh &eh) const
    { return mesh.status(eh).tagged2(); }

    inline void set_enqueued(const Eh &eh, const bool val)
    { mesh.status(eh).set_tagged2(val); }

protected:

    MeshT &mesh;

    const DelaunayT &is_delaunay;

    std::deque<Eh> frontier;
};

template <class MeshT, class DelaunayT>
inline void Flipper<MeshT, DelaunayT>::clear()
{
    for (Eh eh : frontier) { set_enqueued(eh, false); }
    frontier.clear();
}

template <class MeshT, class DelaunayT>
inline void Flipper<MeshT, DelaunayT>::enqueue_all()
{
    for (Eh eh : mesh.edges()) if (!is_enqueued(eh))
    {
        frontier.push_back(eh);
        set_enqueued(eh, true);
    }
}

template <class MeshT, class DelaunayT>
inline void Flipper<MeshT, DelaunayT>::enqueue(const Eh &eh)
{
    if (!is_enqueued(eh))
    {
        frontier.push_back(eh);
        set_enqueued(eh, true);
    }
}

template <class MeshT, class DelaunayT>
inline void Flipper<MeshT, DelaunayT>::enqueue(const Eh ehs[], const int ne)
{
    for (int i = 0; i < ne; ++i) if (!is_enqueued(ehs[i]))
    {
        frontier.push_back(ehs[i]);
        set_enqueued(ehs[i], true);
    }
}

template <class MeshT, class DelaunayT>
template <class PredicateT>
inline void Flipper<MeshT, DelaunayT>::enqueue(const PredicateT &predicate)
{
    for (Eh eh : mesh.edges())
    if (!is_enqueued(eh))
    if (predicate(mesh, eh))
    {
        frontier.push_back(eh);
        set_enqueued(eh, true);
    }
}

template <class MeshT, class DelaunayT>
inline void Flipper<MeshT, DelaunayT>::enqueue_adjacent(const Eh &eh_base)
{
    const Eh ehs[4] {
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh_base, 0))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh_base, 0))),
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh_base, 1))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh_base, 1)))
    };

    for (const Eh &eh : ehs) if (!is_enqueued(eh))
    {
        frontier.push_back(eh);
        set_enqueued(eh, true);
    }
}

template <class MeshT, class DelaunayT>
inline Eh Flipper<MeshT, DelaunayT>::next()
{
    if (frontier.empty()) return Eh {};

    Eh eh = frontier.front();
    set_enqueued(eh, false);
    frontier.pop_front();

    if (mesh.is_boundary(eh)) return eh;
    if (is_delaunay(mesh,eh)) return eh;

    OpenMesh::flip(mesh, eh);
    enqueue_adjacent(eh);

    return eh;
}

template <class MeshT, class DelaunayT>
inline Eh Flipper<MeshT, DelaunayT>::flip()
{
    while (!frontier.empty())
    {
        Eh eh = frontier.front();
        set_enqueued(eh, false);
        frontier.pop_front();

        if (mesh.is_boundary(eh)) continue;
        if (is_delaunay(mesh,eh)) continue;

        OpenMesh::flip(mesh, eh);
        enqueue_adjacent(eh);

        return eh;
    }

    return Eh {};
}

template <class MeshT, class DelaunayT>
inline int Flipper<MeshT, DelaunayT>::flip_all(const int max_n_flip)
{
    int n_flip {};

    for ( ; !frontier.empty() && n_flip < max_n_flip; )
    {
        Eh eh = frontier.front();
        set_enqueued(eh, false);
        frontier.pop_front();

        if (mesh.is_boundary(eh)) continue;
        if (is_delaunay(mesh,eh)) continue;

        OpenMesh::flip(mesh, eh);
        enqueue_adjacent(eh);
        ++n_flip;
    }

    return n_flip;
}

template <class MeshT, class DelaunayT>
inline Flipper<MeshT, DelaunayT> make_flipper(MeshT &mesh, const DelaunayT &is_delaunay)
{ return Flipper<MeshT, DelaunayT>(mesh, is_delaunay); }

////////////////////////////////////////////////////////////////
/// Modification
////////////////////////////////////////////////////////////////

namespace OpenMesh
{

inline void split_copy(TriConnectivity &mesh, EdgeHandle _eh, VertexHandle _vh)
{
    const VertexHandle v0 = mesh.to_vertex_handle(mesh.halfedge_handle(_eh, 0));
    const VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(_eh, 1));

    const size_t nf = mesh.n_faces();

    // Split the halfedge (handle will be preserved)
    mesh.split(_eh, _vh);

    // Copy the properties of the original edge to the corresponding edges
    for (auto voh_it = mesh.voh_iter(_vh); voh_it.is_valid(); ++voh_it)
    {
        EdgeHandle eh = mesh.edge_handle(*voh_it);
        VertexHandle vh = mesh.to_vertex_handle(*voh_it);

        // only copy to the edge that is split from the original one
        if (vh == v0 || vh == v1)
            mesh.copy_all_properties(_eh, eh, true);
    }

    for (VertexHandle vh : {v0, v1})
    {
        // get the halfedge pointing from new vertex to old vertex
        const HalfedgeHandle hh = mesh.find_halfedge(_vh, vh);

        if (!mesh.is_boundary(hh)) // for boundaries there are no faces whose properties need to be copied
        {
            FaceHandle fh0 = mesh.face_handle(hh);
            FaceHandle fh1 = mesh.opposite_face_handle(mesh.prev_halfedge_handle(hh));

            if (static_cast<size_t>(fh0.idx()) >= nf) // is fh0 the new face?
                std::swap(fh0, fh1);

            // copy properties from old face to new face
            mesh.copy_all_properties(fh0, fh1, true);
        }
    }
}

} // namespace OpenMesh

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

template <class MeshT>
inline int get_boundaries(const MeshT&, std::vector<Hh>&);

#endif