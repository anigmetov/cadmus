#pragma once

#include <vector>
#include <unordered_map>

#include "logger.h"
#include "profile.h"
#include "types.h"
#include "cube.h"
#include "matrix_with_values.h"

namespace cadmus {
struct CubicalFiltration {

    static constexpr Int k_invalid_id = std::numeric_limits<Int>::max();

    Cubes cells_;

    std::unordered_map<Int, Int> id_to_sorted_id_;

    CubicalFiltration() = default;
    CubicalFiltration(Cubes&& simplices);
    CubicalFiltration(const CubicalFiltration&) = delete;
    CubicalFiltration(CubicalFiltration&&) = default;
    CubicalFiltration& operator=(const CubicalFiltration&) = delete;
    CubicalFiltration& operator=(CubicalFiltration&&) = default;


    void init();

    size_t size() const { return cells_.size(); }

    void sanity_check();

    // will clear lookup tables
    MatrixWithValues boundary_matrix();
    MatrixWithValues coboundary_matrix();

    friend std::ostream& operator<<(std::ostream& out, const CubicalFiltration& fil);

};
}