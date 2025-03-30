#include "math_def.hh"
#include "mesh.hh"
#include "mesher.hh"
#include "mesh.io.hh"
#include "main.hh"

using namespace OpenMesh;

enum SMOOTHING { NONE, LAPLACIAN, CVT, ODT };

int main(const int argc, const char **argv)
{
    int err {};
    std::string filename, prefix, path;
    if (argc < 2) { printf("No file provided.\n"); return 1; }

    double min_angle = 26.0;
    double max_length = 1e5;
    double max_area = .5e10;
    int smooth_type = NONE;
    int smooth_iter = 20;
    double smooth_step = 1.;

    const char *arg {};
    if ((arg = get_value(argv, argv + argc, "--min-angle")))  { min_angle = atof(arg); }
    if ((arg = get_value(argv, argv + argc, "--max-length"))) { max_length = atof(arg); }
    if ((arg = get_value(argv, argv + argc, "--max-area")))   { max_area = atof(arg); }
    if ((arg = get_value(argv, argv + argc, "--smooth-type"))) { smooth_type = atoi(arg); }
    if ((arg = get_value(argv, argv + argc, "--smooth-iter"))) { smooth_iter = atoi(arg); }
    if ((arg = get_value(argv, argv + argc, "--smooth-step"))) { smooth_step = atof(arg); }

    filename.append(argv[1]);
    prefix = filename.substr(0, filename.find_last_of("."));
    path = filename.substr(0, filename.find_last_of("/\\"));

    std::vector<Vec2> vs {};
    std::vector<Int2> es {};
    std::vector<Vec2> ss {};
    if (read_poly(vs, es, ss, filename.c_str()) != 0)
    { printf("Cannot load poly file.\n"); return 1; }

    TriMesh mesh;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_name())() = prefix;
    getOrMakeProperty<Mh, std::string>(mesh, var_m_path())() = path;

    err = triangulate(vs, es, mesh);
    mesh.delete_isolated_vertices(); // remove dups
    save_mesh(mesh, (prefix + ".CDT.mesh").c_str());
    if (err) { printf("Segments intersecting.\n"); return err; }

    hide_exterior_region(mesh, ss);
    save_mesh(mesh, (prefix + ".CDT.mesh").c_str());

    err = refine(mesh, radian(min_angle), max_length, max_area);
    save_mesh(mesh, (prefix + ".refined.mesh").c_str());
    if (err) { printf("Refinement failed.\n"); return err; }

    if (smooth_type == LAPLACIAN)
    { err = laplacian_smoothing(mesh, smooth_step, smooth_iter); }
    else if (smooth_type == CVT)
    { err = local_CVT_smoothing(mesh, smooth_step, smooth_iter); }
    else if (smooth_type == ODT)
    { err = local_ODT_smoothing(mesh, smooth_step, smooth_iter); }
    if (smooth_type)
    { save_mesh(mesh, (prefix + ".smooth.mesh").c_str()); }
    if (err) { printf("Smoothing failed.\n"); return err; }

    return err;
}