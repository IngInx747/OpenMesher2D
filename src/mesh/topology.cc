#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "mesh.hh"

using namespace OpenMesh;

template <class MeshT>
inline int get_boundaries(const MeshT &mesh, std::vector<Hh> &boundaries)
{
    std::unordered_set<Hh> visit {};

    for (auto hdge : mesh.halfedges())
    if (hdge.is_boundary())
    if (!visit.count(hdge))
    {
        boundaries.push_back(hdge);
        for (auto hbnd : hdge.loop())
            visit.insert(hbnd);
    }

    return (int)boundaries.size();
}

template
int get_boundaries(const TriMesh&, std::vector<Hh>&);

#if 0

template <class MeshT>
inline bool is_topological_disk(const MeshT &mesh)
{
    std::queue<Fh> frontier;
    std::unordered_set<Fh> visited {};

    Fh fh0 = mesh.face_handle(0);
    frontier.push(fh0);
    visited.insert(fh0);

    while (!frontier.empty())
    {
        auto fh = frontier.front(); frontier.pop();

        for (auto hdge : mesh.fh_range(fh))
        if (!hdge.edge().is_boundary())
        if (!is_marked(mesh, hdge.edge()))
        {
            auto fadj = hdge.opp().face();

            if (!visited.count(fadj))
            {
                frontier.push(fadj);
                visited.insert(fadj);
            }
        }
    }

    return visited.size() == mesh.n_faces();
}

template
bool is_topological_disk(const TriMesh &mesh);

template <class MeshT>
inline void prune_branch_edges(MeshT &mesh)
{
    std::queue<Vh> frontier; // valence-1 vertices
    auto valences = makeTemporaryProperty<Vh, int>(mesh); // incident sharp edges number

    for (auto edge : mesh.edges()) if (is_marked(mesh, edge))
    { valences[edge.v0()] += 1; valences[edge.v1()] += 1; }

    // 1. Compute the valence of each vertex, and record all valence-1 vertices.
    // This step just fetches current valence-1 vertices. During pruning, more
    // valence-1 vertices will be generated.
    for (auto vert : mesh.vertices())
    if (valences[vert] == 1) frontier.push(vert);

    // 2. Remove the segments which attached to valence-1 vertices.
    // After getting first batch of valence-1 vertices, prune and find
    // new valence-1 vertices recursively.
    while (!frontier.empty())
    {
        auto vh = frontier.front(); frontier.pop();

#if 1 // ignore the branches incident to a marked vertex
        if (is_marked(mesh, vh)) continue;
#endif

        for (auto hdge : mesh.voh_range(vh))
        {
            auto edge = hdge.edge();
            auto vadj = hdge.to();

            if (is_marked(mesh, edge))
            {
                valences[vadj] -= 1; // update valence of adjacent vertex
                set_marked(mesh, edge, false); // remove branch edge
                if (valences[vadj] == 1) frontier.push(vadj); // update frontier
            }
        }
    }
}

template
void prune_branch_edges(TriMesh &mesh);

template <class MeshT>
inline void remove_exact_loops(MeshT &mesh)
{
    auto visited = makeTemporaryProperty<Fh, bool>(mesh);

    for (auto face : mesh.faces()) if (!visited[face])
    {
        // Tsunami algorithm
        // Flood the entire surface and keep breaking exact loops.
        // Any untouched face outside of the current patch is a
        // potential spot for flooding in next iteration.
        std::queue<std::pair<Eh, Fh>> wave_front;
        wave_front.push({ Eh {}, face });

        while (!wave_front.empty())
        {
            auto ehc = wave_front.front().first;
            auto fhc = wave_front.front().second;
            wave_front.pop();

            // the face is under water
            if (visited[fhc]) continue;

            // break a random fence in the loop
            if (ehc.is_valid()) set_marked(mesh, ehc, false);

            // flood the local patch and record wave fronts
            std::queue<Fh> frontier;
            frontier.push(fhc); visited[fhc] = true;

            while (!frontier.empty())
            {
                auto fh = frontier.front(); frontier.pop();

                for (auto hdge : mesh.fh_range(fh))
                if (!hdge.edge().is_boundary())
                {
                    auto edge = hdge.edge();
                    auto fadj = hdge.opp().face();

                    if (is_marked(mesh, edge))
                    {
                        wave_front.push({ edge, fadj });
                    }
                    else if (!visited[fadj])
                    {
                        frontier.push(fadj); visited[fadj] = true;
                    }
                }
            }
        }
    }

    prune_branch_edges(mesh);
}

template
void remove_exact_loops(TriMesh &mesh);

template <class MeshT>
inline void remove_local_spikes(MeshT &mesh)
{}

///   *      *          *------*
///    \    /   ----->   
///     \  /    (switch  
///      \/   sharp edges) 
template <>
static void remove_local_spikes(TriMesh &mesh)
{
    for (auto face : mesh.faces())
    {
        int ns {}; // number of sharp edges
        for (auto edge : face.edges())
            if (is_marked(mesh, edge)) ++ns;

        if (ns == 2)
        {
            auto hdge = face.halfedge();
            while (is_marked(mesh, hdge.edge()))
                hdge = hdge.next();
            auto vert = hdge.next().to();

            if (vert.is_boundary()) continue;
            if (is_marked(mesh, vert)) continue;

            int valence {};
            for (auto edge : vert.edges())
                if (is_marked(mesh, edge)) ++valence;
            if (valence != 2) continue;

            for (auto edge : face.edges())
                set_marked(mesh, edge, !is_marked(mesh, edge));
        }
    }
}

template <class MeshT>
inline void generate_cut_graph(MeshT &mesh)
{
    auto sharped = makeTemporaryProperty<Eh, bool>(mesh);
    auto visited = makeTemporaryProperty<Fh, bool>(mesh);

    for (auto edge : mesh.edges())
        sharped[edge] = true;

    for (auto face : mesh.faces()) if (!visited[face])
    {
        std::queue<Fh> frontier;

        // traverse the connecting component
        frontier.push(face);
        visited[face] = true;

        // Flood over the surface till its waves confront.
        // The flood does not go across a pre-marked edge,
        // thus preserves prescribed layout on the mesh.
        while (!frontier.empty())
        {
            auto fh = frontier.front(); frontier.pop();

            for (auto hdge : mesh.fh_range(fh))
            if (!hdge.edge().is_boundary())    // is cliff on the other side?
            if (!is_marked(mesh, hdge.edge())) // is blocked by user constraint?
            {
                auto edge = hdge.edge();
                auto fadj = hdge.opp().face();

                if (!visited[fadj])
                {
                    frontier.push(fadj);
                    sharped[edge] = false;
                    visited[fadj] = true;
                }
            }
        }
    }

    for (auto edge : mesh.edges()) if (sharped[edge])
        set_marked(mesh, edge, true);

    prune_branch_edges (mesh);
    remove_local_spikes(mesh);
}

template
void generate_cut_graph(TriMesh &mesh);

/// @brief Virtually slice the mesh by re-indexing the halfedges
/// @param mesh with all the edges to be sliced marked
/// @return Number of vertices of the sliced mesh
template <class MeshT>
inline int slice_virtual(MeshT &mesh)
{
    auto h_i = getOrMakeProperty<Hh, int>(mesh, var_h_index());

    const auto is_seam = [](const MeshT &mesh, const Hh &hh) -> bool {
        return is_marked(mesh, mesh.edge_handle(hh)) && !mesh.is_boundary(hh) ||
            mesh.is_boundary(mesh.opposite_halfedge_handle(hh)); };

    int nv {}; // number of vertices of the sliced mesh

    // index with ordinary vertex order
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++nv)
        for (auto hdge : (*v_it).incoming_halfedges()) h_i[hdge] = nv;

    // duplicate vertex along the seams
    for (auto vert : mesh.vertices())
    {
        // next counter-clockwise sharp halfedge (exclusion)
        auto hdge1 = vert.halfedge().prev().opp();
        do { if (is_seam(mesh, hdge1)) break; }
        while ((hdge1 = hdge1.prev().opp()) != vert.halfedge());

        // next clockwise sharp halfedge (inclusion)
        auto hdge2 = vert.halfedge();
        do { if (is_seam(mesh, hdge2)) break; }
        while ((hdge2 = hdge2.opp().next()) != vert.halfedge());

        //if (hdge1 == hdge2) continue; // valence <= 1
        for (auto hdge = hdge1.prev(); hdge != hdge2.prev(); hdge = hdge.opp().prev()) //if (!hdge.is_boundary())
        { h_i[hdge] = nv; if (is_seam(mesh, hdge)) ++nv; } // traverse valence-1 fans (ccw incoming halfedges)
    }

    // most-clockwise incoming halfedges on the boundary are always misindexed
    for (auto hdge : mesh.halfedges()) if (hdge.is_boundary())
        h_i[hdge] = h_i[hdge.opp().prev()];

    return nv;
}

template
int slice_virtual(TriMesh &mesh);

#endif
