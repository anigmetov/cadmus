#include "profile.h"
#include "simplex.h"

namespace cadmus {

Simplex::Simplex(const IdxVector& _vertices, Real _value)
        :
        vertices_(_vertices), value_(_value)
{
    if (vertices_.empty())
        throw std::runtime_error("Empty simplex not allowed");

    std::sort(vertices_.begin(), vertices_.end());
}

Simplex::Simplex(Int _id, const IdxVector& _vertices, Real _value, std::vector<diy::Point<int, CADMUS_DIM>> _u_globals)
        :
        vertices_(_vertices), value_(_value), id_(_id), u_globals(_u_globals)
{
    if (vertices_.empty())
        throw std::runtime_error("Empty simplex not allowed");

    if (vertices_.size() > 1)
        std::sort(vertices_.begin(), vertices_.end());
}

dim_type Simplex::dim() const
{
    return static_cast<dim_type>(vertices_.size()) - 1;
}

Real Simplex::value() const
{
    return value_;
}

UidValue Simplex::uid_value() const
{
    return {id_, value()};
}

std::vector<IdxVector> Simplex::boundary() const
{
    std::vector<IdxVector> bdry;
    bdry.reserve(vertices_.size());

    for(size_t i = 0 ; i < vertices_.size() ; ++i) {
        IdxVector tau;
        tau.reserve(vertices_.size() - 1);

        for(size_t j = 0 ; j < vertices_.size() ; ++j)
            if (j != i)
                tau.push_back(vertices_[j]);

        // vertices_ is sorted -> tau is sorted automatically

        bdry.push_back(tau);
    }

    return bdry;
}

Int Simplex::sorted_id() const
{
    return sorted_id_;
}

void Simplex::set_sorted_id(Int new_sorted_id)
{
    sorted_id_ = new_sorted_id;
}

Int Simplex::id() const
{
    return id_;
}

IdxVector Simplex::vertices() const
{
    return vertices_;
}

std::string Simplex::repr() const
{
    std::stringstream ss;
    ss << "Simplex(id_=" << id_ << ", sorted_id_ = " << sorted_id_ << ", vertices_=(";

    for(size_t i = 0 ; i < vertices_.size() - 1 ; ++i)
        ss << vertices_[i] << ", ";

    ss << vertices_[vertices_.size() - 1] << "), ";

    ss << ", u_globals = (";
    for(size_t i = 0 ; i < u_globals.size() - 1 ; ++i)
        ss << "(" << u_globals[i] << "), ";
    ss << "(" << u_globals[u_globals.size() - 1] << ")), ";

    ss << "value_=" << value_ << ")";

    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const Simplex& s)
{
    out << "Simplex(id_=" << s.id_ << ", sorted_id_ = " << s.sorted_id_ << ", vertices_=(";

    for(size_t i = 0 ; i < s.vertices_.size() - 1 ; ++i)
        out << s.vertices_[i] << ", ";

    out << s.vertices_[s.vertices_.size() - 1] << "), ";

    out << ", u_globals = (";
    for(size_t i = 0 ; i < s.u_globals.size() - 1 ; ++i)
        out << "(" << s.u_globals[i] << "), ";
    out << "(" << s.u_globals[s.u_globals.size() - 1] << ")), ";

    out << "value_=" << s.value_ << ")";

    return out;
}

bool Simplex::operator==(const Simplex& other) const
{
    return sorted_id_ == other.sorted_id_ and vertices_ == other.vertices_ and value_ == other.value_;
}

bool Simplex::operator!=(const Simplex& other) const
{
    return !(*this == other);
}

} // namespace cadmus

