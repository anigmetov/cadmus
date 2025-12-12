#include <bit>

#include "logger.h"
#include "profile.h"
#include "cube.h"

namespace cadmus {

    GridRef Cube::global_domain_ {nullptr, Point(), true};

    Int get_vertex_part(Int cube_id)
    {
        // extract last CADMUS_DIM bits of id
        return cube_id >> CADMUS_DIM;
    }

    Point Cube::get_vertex() const
    {
        // extract last CADMUS_DIM bits of id
        Int v_idx = ::cadmus::get_vertex_part(id_);
        return global_domain_.vertex(v_idx);
    }

    dim_type get_dim(Int cube_id)
    {
        dim_type res = 0;
        for(int d = 0; d < CADMUS_DIM; ++d) {
            if (cube_id & (1 << d))
                res++;
        }
        return res;
    }

    dim_type Cube::dim() const
    {
        return get_dim(id_);
    }

    std::vector<Point> get_cube_vertices(Int cube_id)
    {
        std::vector<Point> vertices;

        Int vertex_id = get_vertex_part(cube_id);

        Point vertex = Cube::global_domain_.vertex(vertex_id);

        // displacements = standard basis vectors whose bits are set in cube_id
        std::vector<diy::Point<int, CADMUS_DIM>> disps;

        for(dim_type d = 0; d < CADMUS_DIM; ++d) {
            if (cube_id & (1 << d)) {
                Point pt = Point::zero();
                pt[d] = 1;
                disps.push_back(pt);
            }
        }

        int n_vertices = (1 << disps.size());

        // vertex is sum of basis vectors belonging to an arbitrary subset of disps
        // we need to iterate over all subsets, vertex_idx = subset index
        for(int vertex_idx = 0; vertex_idx < n_vertices; ++vertex_idx) {
            Point v = vertex;
            for(int disp_idx = 0; disp_idx < disps.size(); ++disp_idx) {
                if (vertex_idx & (1 << disp_idx)) {
                    v += disps[disp_idx];
                }
            }
            vertices.push_back(v);
        }

        return vertices;
    }

    std::vector<Int> Cube::boundary() const
    {
        auto logger = spd::get("console");

        std::vector<Int> result;

        std::vector<int> disps;

        for(int d = 0; d < CADMUS_DIM; ++d) {
            if (id() & (1 << d)) {
                disps.push_back(d);
            }
        }

//        if (logger) logger->debug("compute boundary: cube_id = {:b}, displacements = {}", id(), container_to_string(disps));
        // add faces from current vertex

        for(int k = 0; k < disps.size(); ++k) {
            assert(id() & (1 << disps[k]));
            // unset the disps[k] bit
            Int face_id = id() & ~(1 << disps[k]);
//            if (logger) logger->debug("compute boundary: cube_id = {:b}, disps[{}] = {}, face_id = {:b}, dim = {}, face dim = {}", id(), k, disps[k], face_id, get_dim(id()), get_dim(face_id));
            assert(::cadmus::get_vertex_part(id()) == ::cadmus::get_vertex_part(face_id));
            assert(get_dim(id()) == get_dim(face_id) + 1);
            result.push_back(face_id);
        }

        // pick one displacement vector, move vertex along it
        // add face defined by the remaining vectors
        for(int k = 0; k < disps.size(); ++k) {
            Point p = get_vertex();
            p[disps[k]] += 1;
            Int face_id = global_domain_.index(p) << CADMUS_DIM;
            // set the bits corresponding to other displacements
            for(int j = 0; j < disps.size(); ++j) {
                if (j == k)
                    continue;
                face_id |= (1 << disps[j]);
            }
//            if (logger) logger->debug("compute boundary: other vertex, cube_id = {:b}, disps[{}] = {}, face_id = {:b}, p = {}, get_vertex = {}, dim ={}, face dim = {}", id(), k, disps[k], face_id, p, get_vertex(), get_dim(id()), get_dim(face_id));
//            assert(::cadmus::get_vertex_part(id()) = ::cadmus::get_vertex_part(face_id));
            assert(get_dim(id()) == get_dim(face_id) + 1);
            result.push_back(face_id);
        }

        return result;
    }

    UidValue Cube::uid_value() const
    {
        return {id(), value()};
    }

    Real Cube::value() const
    {
        return value_;
    }

    Int Cube::sorted_id() const
    {
        return sorted_id_;
    }

    void Cube::set_sorted_id(Int new_sorted_id)
    {
        sorted_id_ = new_sorted_id;
    }

    Int Cube::id() const
    {
        return id_;
    }

//std::string Cube::repr() const
//{
//    std::stringstream ss;
//    ss << "Cube(id_=" << id_ << ", sorted_id_ = " << sorted_id_ << ", vertices_=(";
//
//    for(size_t i = 0 ; i < vertices_.size() - 1 ; ++i)
//        ss << vertices_[i] << ", ";
//
//    ss << vertices_[vertices_.size() - 1] << "), ";
//
//    ss << ", u_globals = (";
//    for(size_t i = 0 ; i < u_globals.size() - 1 ; ++i)
//        ss << "(" << u_globals[i] << "), ";
//    ss << "(" << u_globals[u_globals.size() - 1] << ")), ";
//
//    ss << "value_=" << value_ << ")";
//
//    return ss.str();
//}
//
    std::ostream& operator<<(std::ostream& out, const Cube& s)
    {
        out << "Cube(id_=" << s.id_ << ", sorted_id_ = " << s.sorted_id_ << ", vertices=(";

        auto vs = get_cube_vertices(s.id());

        for(size_t i = 0; i < vs.size() - 1; ++i)
            out << vs[i] << ", ";

        out << vs[vs.size() - 1] << "), ";
        out << "value_=" << s.value_ << ")";

        return out;
    }

    bool Cube::operator==(const Cube& other) const
    {
        return id() == other.id();
    }

    bool Cube::operator!=(const Cube& other) const
    {
        return !(*this == other);
    }

} // namespace cadmus

