#pragma once

#include "types.h"

namespace cadmus {

    class Cube {
    public:
        Cube() = default;

        Cube(Int _id, Real _value) : id_(_id), value_(_value) {};

        static GridRef global_domain_;

        dim_type dim() const;
        Real value() const;
        Int sorted_id() const;
        void set_sorted_id(Int new_sorted_id);
        Int id() const;
        std::vector<Int> boundary() const;
        UidValue uid_value() const;

        bool operator==(const Cube& other) const;
        bool operator!=(const Cube& other) const;

        std::string repr() const;

        friend std::ostream& operator<<(std::ostream&, const Cube&);

        Point get_vertex() const;
    private:
        static constexpr Int k_invalid_id = Int(-1);

        Int id_ {k_invalid_id};
        Int sorted_id_ {k_invalid_id};
        Real value_ {-std::numeric_limits<Real>::max()};
    };

    using Cubes = std::vector<Cube>;

    std::vector<Point> get_cube_vertices(Int id);
}