#pragma once

#include <limits>
#include <exception>
#include <cassert>
#include <algorithm>
#include <ostream>

#include "types.h"
#include "logger.h"

namespace cadmus {

class Simplex {
public:
    using Vertices = IdxVector;
    Simplex() = default;

    Simplex(const IdxVector& _vertices, Real _value);

    Simplex(Int _id, const IdxVector& _vertices, Real _value, std::vector<diy::Point<int, CADMUS_DIM>> );

    dim_type dim() const;
    Real value() const;
    Int sorted_id() const;
    void set_sorted_id(Int new_sorted_id);
    Int id() const;
    IdxVector vertices() const;
    std::vector<IdxVector> boundary() const;
    UidValue uid_value() const;

    std::vector<diy::Point<int, CADMUS_DIM>> u_globals;

    bool operator==(const Simplex& other) const;
    bool operator!=(const Simplex& other) const;

    std::string repr() const;

    friend std::ostream& operator<<(std::ostream&, const Simplex&);
private:
    static constexpr Int k_invalid_id = Int(-1);

    Int id_ {k_invalid_id};
    Int sorted_id_ {k_invalid_id};
    Vertices vertices_;
    Real value_ {std::numeric_limits<Real>::max()};

};

using Simplices = std::vector<Simplex>;

} // namespace cadmus
