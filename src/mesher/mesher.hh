#ifndef MESHER2D_HH
#define MESHER2D_HH

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mesh.hh"

int make_delaunay(TriMesh&);

int triangulate(
    const std::vector<Vec2>&,
    const std::vector<Int2>&,
    TriMesh&,
    std::unordered_map<Vh, Vh>         &duplicated_vertices,
    std::vector<std::pair<Int2, Int2>> &intersecting_segments,
    std::vector<std::pair<Int2, int>>  &overlapping_vertices);

int triangulate(
    const std::vector<Vec2>&,
    const std::vector<Int2>&,
    TriMesh&,
    std::unordered_map<Vh, Vh> &duplicated_vertices);
    
int triangulate(
    const std::vector<Vec2>&,
    const std::vector<Int2>&,
    TriMesh&);

Hh restore_segment(TriMesh&, const Vh &vertex0, const Vh &vertex1);

int hide_exterior_region(TriMesh&, const std::vector<Vec2> &seeds);

int refine(TriMesh&, const double min_angle, const double max_length, const double max_area);

int laplacian_smoothing(TriMesh&, const double step, const int max_num_iter);

int local_CVT_smoothing(TriMesh&, const double step, const int max_num_iter);

#endif