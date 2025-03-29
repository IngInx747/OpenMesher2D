#include <string>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.hh"

using namespace OpenMesh;

enum IO_ERR_TYPE
{
    NO_ERROR,
    INTERNAL,
    CANNOT_OPEN,
    UNSUPPORTED_FORMAT,
};

static const char *__io_err_msg[] = {
    "",
    "internal error",
    "cannot open file",
    "unsupported format"
};

////////////////////////////////////////////////////////////////
/// Utils
////////////////////////////////////////////////////////////////

inline std::string get_file_extension(const char *filename)
{
    std::string s { filename };
    auto found = s.find_last_of(".");
    return s.substr(found + 1);
}

////////////////////////////////////////////////////////////////
/// Raw data
////////////////////////////////////////////////////////////////

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char *filename,
    const int offset = 1,
    const std::streamsize prec = 17);

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char *filename,
    const int offset,
    const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::defaultfloat << std::setprecision(prec);

    for (int i = 0; i < nv; ++i)
        out << "v "
            << vs[i*3 + 0] << " "
            << vs[i*3 + 1] << " "
            << vs[i*3 + 2] << " "
            << "\n";

    for (int i = 0; i < nf; ++i)
        out << "f "
            << fs[i*3 + 0] + offset << " "
            << fs[i*3 + 1] + offset << " "
            << fs[i*3 + 2] + offset << " "
            << "\n";

    for (int i = 0; i < ne; ++i)
        out << "l "
            << es[i*2 + 0] + offset << " "
            << es[i*2 + 1] + offset << " "
            << "\n";

    return 0;
}

int read_node(std::vector<Vec2> &vs, const char *filename)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    int ver {}, dim {}, nv {}, ne {};

    in >> nv;
    in >> dim;
    in.ignore(kInf, '\n');

    vs.reserve(nv);

    for (int i = 0; i < nv; ++i)
    {
        int vid {}; in >> vid;
        Vec2 p {}; in >> p[0] >> p[1];
        vs.push_back(p);
        // optional: process vertex tag
        in.ignore(kInf, '\n');
    }

    return IO_ERR_TYPE::NO_ERROR;
}

int save_node(const double* vs, const int nv, const char* filename, const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::defaultfloat << std::setprecision(prec);

    out << nv << " 2\n";

    for (int i = 0; i < nv; ++i)
        out << i + 1 << " "
        << vs[i*2 + 0] << " "
        << vs[i*2 + 1] << " "
        << "\n";

    return 0;
}

int read_poly(
    std::vector<Vec2> &vs,
    std::vector<Int2> &es,
    std::vector<Vec2> &ss,
    const char *filename,
    const int offset)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    int ver {}, dim {}, nv {}, ne {}, ns {};

    in >> nv;
    in >> dim;
    in.ignore(kInf, '\n');

    vs.reserve(nv);

    for (int i = 0; i < nv; ++i)
    {
        int vid {}; in >> vid;
        Vec2 p {}; in >> p[0] >> p[1];
        vs.push_back(p);
        // optional: process vertex tag
        in.ignore(kInf, '\n');
    }

    in >> ne;
    in.ignore(kInf, '\n');

    es.reserve(ne);

    for (int i = 0; i < ne; ++i)
    {
        int eid {}; in >> eid;
        Int2 vv; in >> vv[0] >> vv[1];
        es.emplace_back(vv[0] - offset, vv[1] - offset);
        // optional: process edge tag
        in.ignore(kInf, '\n');
    }

    in >> ns;
    in.ignore(kInf, '\n');

    ss.reserve(ns);

    for (int i = 0; i < ns; ++i)
    {
        int vid {}; in >> vid;
        Vec2 p {}; in >> p[0] >> p[1];
        ss.push_back(p);
        // optional: process vertex tag
        in.ignore(kInf, '\n');
    }

    return IO_ERR_TYPE::NO_ERROR;
}

int save_poly(
    const double *vs, const int nv,
    const int    *es, const int ne,
    const double *ss, const int ns,
    const char* filename,
    const int offset,
    const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::defaultfloat << std::setprecision(prec);

    out << nv << " 2\n";

    for (int i = 0; i < nv; ++i)
        out << i + 1 << " "
            << vs[i*2 + 0] << " "
            << vs[i*2 + 1] << " "
            << -1 << "\n";

    out << ne << " 1\n";

    for (int i = 0; i < ne; ++i)
        out << i + 1 << " "
            << es[i*2 + 0] + offset << " "
            << es[i*2 + 1] + offset << " "
            << -1 << "\n";

    out << ns << "\n";

    for (int i = 0; i < ns; ++i)
        out << i + 1 << " "
            << ss[i*2 + 0] << " "
            << ss[i*2 + 1] << " "
            << "\n";

    return 0;
}

int read_mesh(
    std::vector<Vec2> &vs,
    std::vector<Int3> &fs,
    std::vector<Int2> &es,
    const char *filename,
    const int offset)
{
    // [TODO]

    return IO_ERR_TYPE::UNSUPPORTED_FORMAT;
}

int save_mesh(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char* filename,
    const int offset,
    const std::streamsize prec)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    out << std::defaultfloat << std::setprecision(prec);

    out << "MeshVersionFormatted 1\n\nDimension\n2\n\n";

    out << "Vertices\n" << nv << "\n";

    for (int i = 0; i < nv; ++i)
        out << vs[i*2 + 0] << " "
            << vs[i*2 + 1] << " "
            << " 0\n";

    out << "Triangles\n" << nf << "\n";

    for (int i = 0; i < nf; ++i)
        out << fs[i*3 + 0] + offset << " "
            << fs[i*3 + 1] + offset << " "
            << fs[i*3 + 2] + offset << " "
            << " 1\n";

    out << "Edges\n" << ne << "\n";

    for (int i = 0; i < ne; ++i)
    {
        out << es[i*2 + 0] + offset << " "
            << es[i*2 + 1] + offset << " "
            << " -1\n";
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

////////////////////////////////////////////////////////////////
/// Mesh IO
////////////////////////////////////////////////////////////////

template <class MeshT>
inline int reindex(const MeshT &mesh, std::unordered_map<Vh, int> &v2i)
{
    int nv {};

    for (Vh vh : mesh.vertices())
    { v2i[vh] = nv++; }

    return nv;
}

template <class MeshT>
inline int reindex(const MeshT &mesh, std::unordered_map<Fh, int> &f2i)
{
    int nf {};

    for (Fh fh : mesh.faces())
    { f2i[fh] = nf++; }

    return nf;
}

template <class MeshT>
inline int reindex(const MeshT &mesh, std::unordered_map<Eh, int> &e2i)
{
    int ne {};

    for (Eh eh : mesh.edges())
    { e2i[eh] = ne++; }

    return ne;
}

template <class MeshT>
inline int read_mesh_builtin(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    IO::Options opt;

    opt += IO::Options::VertexTexCoord;
    mesh.request_vertex_texcoords2D();

    opt += IO::Options::VertexNormal;
    mesh.request_vertex_normals();

    opt += IO::Options::VertexColor;
    mesh.request_vertex_colors();

    opt += IO::Options::FaceNormal;
    mesh.request_face_normals();

    opt += IO::Options::FaceColor;
    mesh.request_face_colors();

    if (!IO::read_mesh(mesh, filename, opt)) err = IO_ERR_TYPE::INTERNAL;

    if (!opt.vertex_has_texcoord()) mesh.release_vertex_texcoords2D();

    if (!opt.vertex_has_normal()) mesh.release_vertex_normals();

    if (!opt.vertex_has_color()) mesh.release_vertex_colors();

    if (!opt.face_has_normal()) mesh.release_face_normals();

    if (!opt.face_has_color()) mesh.release_face_colors();

    return err;
}

template <class MeshT>
inline int save_mesh_builtin(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    IO::Options opt;

    if (mesh.has_vertex_texcoords2D()) opt += IO::Options::VertexTexCoord;

    if (mesh.has_vertex_normals()) opt += IO::Options::VertexNormal;

    if (mesh.has_vertex_colors()) opt += IO::Options::VertexColor;

    if (mesh.has_face_normals()) opt += IO::Options::FaceNormal;

    if (mesh.has_face_colors()) opt += IO::Options::FaceColor;

    if (mesh.has_halfedge_texcoords2D()) opt += IO::Options::FaceTexCoord;

    if (!IO::write_mesh(mesh, filename, opt, 17)) err = IO_ERR_TYPE::INTERNAL;

    return err;
}

template <class MeshT>
inline int read_mesh_dot_mesh(MeshT &mesh, const char *filename)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    std::vector<Int3> fs;
    std::vector<Int2> es;
    int ver {}, dim {}, nv {}, nf {}, ne {};

    while (in)
    {
        std::string buf {};
        in >> buf;

        if (buf.compare("MeshVersionFormatted") == 0) // version
        {
            in >> ver;
        }
        else if (buf.compare("Dimension") == 0) // dimension
        {
            in >> dim;
        }
        else if (buf.compare("Vertices") == 0) // vertices
        {
            in >> nv;
            in.ignore(kInf, '\n');
            for (int i = 0; i < nv; ++i)
            {
                TriMesh::Point p { 0,0,0 };
                for (int j = 0; j < dim; ++j)
                    in >> p[j];
                mesh.new_vertex(p);
                // optional: process vertex tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Triangles") == 0) // faces
        {
            in >> nf;
            in.ignore(kInf, '\n');
            fs.reserve(nf);
            for (int i = 0; i < nf; ++i)
            {
                int ids[3];
                in >> ids[0] >> ids[1] >> ids[2];
                fs.emplace_back(ids[0]-1, ids[1]-1, ids[2]-1);
                // optional: process facet tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Edges") == 0) // edge
        {
            in >> ne;
            in.ignore(kInf, '\n');
            es.reserve(ne);
            for (int i = 0; i < ne; ++i)
            {
                int ids[2];
                in >> ids[0] >> ids[1];
                es.emplace_back(ids[0]-1, ids[1]-1);
                // optional: process edge tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("End") == 0)
        {}
    }

    for (const auto &face : fs)
    {
        int i0 = face[0];
        int i1 = face[1];
        int i2 = face[2];
        mesh.add_face({
            mesh.vertex_handle(i0),
            mesh.vertex_handle(i1),
            mesh.vertex_handle(i2)
        });
    }

    for (const auto &edge : es)
    {
        int i0 = edge[0];
        int i1 = edge[1];
        auto hdge = mesh.find_halfedge(mesh.vertex_handle(i0), mesh.vertex_handle(i1));
        if (hdge.is_valid()) set_sharp(mesh, hdge.edge(), true);
    }

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_mesh_dot_mesh(MeshT &mesh, const char *filename)
{
    return IO_ERR_TYPE::UNSUPPORTED_FORMAT;
}

template <>
int save_mesh_dot_mesh(const TriMesh &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    std::unordered_map<Vh, int> v2i {};
    std::unordered_map<Fh, int> f2i {};

    const int nv = reindex(mesh, v2i);
    const int nf = reindex(mesh, f2i);

    out << std::defaultfloat << std::setprecision(17);

    out << "MeshVersionFormatted 1\n\nDimension\n3\n\n";

    out << "Vertices\n" << nv << "\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "0\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vert] << "\n";
        }
    }

    out << "Triangles\n" << nf << "\n";

    if (!hasProperty<Fh, int>(mesh, var_f_label()))
    {
        for (auto face : mesh.faces())
        {
            for (auto vert : face.vertices())
                out << v2i[vert] + 1 << " ";
            out << "1\n"; // arbitrary positive tag
        }
    }
    else
    {
        auto f_tag = getProperty<Fh, int>(mesh, var_f_label());

        for (auto face : mesh.faces())
        {
            for (auto vert : face.vertices())
                out << v2i[vert] + 1 << " ";
            out << f_tag[face] << "\n";
        }
    }

    int ne_sharp {}; for (Eh eh : mesh.edges()) { if (is_sharp(mesh, eh)) ++ne_sharp; }

    out << "Edges\n" << ne_sharp << "\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto edge : mesh.edges()) if (is_sharp(mesh, edge))
        {
            auto vert0 = edge.v0(), vert1 = edge.v1();
            out << v2i[vert0] + 1 << " "
                << v2i[vert1] + 1 << " -1\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto edge : mesh.edges()) if (is_sharp(mesh, edge))
        {
            auto vert0 = edge.v0(), vert1 = edge.v1();
            out << v2i[vert0] + 1 << " "
                << v2i[vert1] + 1 << " "
                << e_tag[edge] << "\n";
        }
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_mesh_dot_node(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    std::unordered_map<Vh, int> v2i {};
    const int nv = reindex(mesh, v2i);

    out << std::defaultfloat << std::setprecision(17);

    out << nv << " 3\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (Vh vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << " "
                << p[0] << " "
                << p[1] << " "
                << p[2] << " "
                << "0\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (Vh vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << " "
                << p[0] << " "
                << p[1] << " "
                << p[2] << " "
                << v_tag[vh] << "\n";
        }
    }

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_mesh_dot_poly(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    std::unordered_map<Vh, int> v2i {};
    const int nv = reindex(mesh, v2i);

    out << std::defaultfloat << std::setprecision(17);

    out << nv << " 3 0 1\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "-1\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vh] << "\n";
        }
    }

    int ne {}; for (Eh eh : mesh.edges()) { ++ne; }

    //int ne_sharp {}; for (Eh eh : mesh.edges()) { if (is_sharp(mesh, eh)) ++ne_sharp; }

    out << ne << " 1\n";

    ne = 0;

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto eh : mesh.edges())// if (is_sharp(mesh, eh))
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << (++ne) << " "
                << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << -1 << "\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto eh : mesh.edges())// if (is_sharp(mesh, eh))
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << (++ne) + 1 << " "
                << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << e_tag[eh] << "\n";
        }
    }

    out << "0\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
int read_mesh(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        err = read_mesh_builtin(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = read_mesh_dot_mesh(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "read_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template <class MeshT>
int save_mesh(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        err = save_mesh_builtin(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = save_mesh_dot_mesh(mesh, filename);
    }
    else if (ext.compare("node") == 0)
    {
        err = save_mesh_dot_node(mesh, filename);
    }
    else if (ext.compare("poly") == 0)
    {
        err = save_mesh_dot_poly(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "save_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template
int read_mesh(TriMesh&, const char*);

template
int read_mesh(PolyMesh&, const char*);

template
int save_mesh(const TriMesh&, const char*);

template
int save_mesh(const PolyMesh&, const char*);

////////////////////////////////////////////////////////////////
/// Poly IO
////////////////////////////////////////////////////////////////

template <class MeshT>
inline int read_poly_dot_poly(MeshT &mesh, const char *filename)
{
    return IO_ERR_TYPE::UNSUPPORTED_FORMAT;
}

template <class MeshT>
inline int save_poly_obj(const MeshT &mesh, const char *filename)
{
    std::unordered_map<Vh, int> v2i {};
    std::unordered_map<Eh, int> e2i {};

    const int nv = reindex(mesh, v2i);
    const int ne = reindex(mesh, e2i);

    std::vector<int> es;

    for (auto eh : mesh.edges())
    {
        Hh hh = mesh.halfedge_handle(eh, 0);
        es.push_back(v2i[mesh.from_vertex_handle(hh)] + 1);
        es.push_back(v2i[mesh.to_vertex_handle  (hh)] + 1);
    }

    return save_obj(
        mesh.points()->data(), nv,
        nullptr, 0,
        es.data(), ne,
        filename);
}

template <class MeshT>
inline int save_poly_dot_poly(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    std::unordered_map<Vh, int> v2i {};
    std::unordered_map<Eh, int> e2i {};

    const int nv = reindex(mesh, v2i);
    const int ne = reindex(mesh, e2i);

    out << std::defaultfloat << std::setprecision(17);

    out << nv << " 3 0 1\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "-1\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vh : mesh.vertices())
        {
            const auto &p = mesh.point(vh);
            out << v2i[vh] + 1 << "  "
                << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vh] << "\n";
        }
    }

    out << ne << " 1\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << e2i[eh] + 1 << " "
                << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << -1 << "\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << e2i[eh] + 1 << " "
                << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << e_tag[eh] << "\n";
        }
    }

    out << "0\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
inline int save_poly_dot_mesh(const MeshT &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    std::unordered_map<Vh, int> v2i {};
    std::unordered_map<Eh, int> e2i {};

    const int nv = reindex(mesh, v2i);
    const int ne = reindex(mesh, e2i);

    out << std::defaultfloat << std::setprecision(17);

    out << "MeshVersionFormatted 1\n\nDimension\n3\n\n";

    out << "Vertices\n" << nv << "\n";

    if (!hasProperty<Vh, int>(mesh, var_v_label()))
    {
        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << "0\n"; // arbitrary tag
        }
    }
    else
    {
        auto v_tag = getProperty<Vh, int>(mesh, var_v_label());

        for (auto vert : mesh.vertices())
        {
            const auto &p = mesh.point(vert);
            out << p[0] << "  "
                << p[1] << "  "
                << p[2] << "  "
                << v_tag[vert] << "\n";
        }
    }

    out << "Triangles\n" << 0 << "\n";

    out << "Edges\n" << ne << "\n";

    if (!hasProperty<Eh, int>(mesh, var_e_label()))
    {
        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << "-1\n";
        }
    }
    else
    {
        auto e_tag = getProperty<Eh, int>(mesh, var_e_label());

        for (auto eh : mesh.edges())
        {
            Hh hh = mesh.halfedge_handle(eh, 0);
            out << v2i[mesh.from_vertex_handle(hh)] + 1 << " "
                << v2i[mesh.to_vertex_handle  (hh)] + 1 << " "
                << e_tag[eh] << "\n";
        }
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

template <class MeshT>
int read_poly(MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("poly") == 0)
    {
        err = read_poly_dot_poly(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "read_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template <class MeshT>
int save_poly(const MeshT &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };

    const auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0)
    {
        err = save_poly_obj(mesh, filename);
    }
    else if (ext.compare("poly") == 0)
    {
        err = save_poly_dot_poly(mesh, filename);
    }
    else if (ext.compare("mesh") == 0)
    {
        err = save_poly_dot_mesh(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err) fprintf(stderr, "save_mesh error %d: %s\n", err, __io_err_msg[err]);

    return err;
}

template
int save_poly(const PolyMesh&, const char*);

template
int save_poly(const TriMesh&, const char*);

////////////////////////////////////////////////////////////////
/// Primitives IO
////////////////////////////////////////////////////////////////

template <class MeshT>
inline int save_marked_face_centroids(const MeshT &mesh, const char *filename, const double offset)
{
    std::vector<Vec3> ps {};
    int np {};

    for (const auto face : mesh.faces()) if (is_marked(mesh, face))
    {
        const auto p = mesh.calc_face_centroid(face);
        const auto n = mesh.calc_normal(face) * offset;
        ps.push_back(p + n);
        ++np;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        nullptr,                  0,
        filename);
}

template
int save_marked_face_centroids(const TriMesh&, const char*, const double);

template
int save_marked_face_centroids(const PolyMesh&, const char*, const double);

template <class MeshT>
int save_marked_vertices(const MeshT &mesh, const char *filename, const double offset)
{
    std::vector<Vec3> ps {};
    int np {};

    for (const auto vert : mesh.vertices()) if (is_marked(mesh, vert))
    {
        const auto p = mesh.point(vert);
        const auto n = mesh.calc_normal(vert) * offset;
        ps.push_back(p + n);
        ++np;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        nullptr,                  0,
        filename);
}

template
int save_marked_vertices(const TriMesh&, const char*, const double);

template
int save_marked_vertices(const PolyMesh&, const char*, const double);
