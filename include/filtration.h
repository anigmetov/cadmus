#pragma once

#include <vector>
#include <unordered_map>

#include "logger.h"
#include "profile.h"
#include "types.h"
#include "simplex.h"
#include "matrix_with_values.h"

namespace cadmus {
struct SimplicialFiltration {

    static constexpr Int k_invalid_id = std::numeric_limits<Int>::max();

    using SimplexVertices = Simplex::Vertices;

    Simplices cells_;

    // TODO: replace with unordered_map, need a good hash for vector of vertices
    std::map<SimplexVertices, Int> vertices_to_sorted_id_;
    std::unordered_map<Int, Int> id_to_sorted_id_;
    std::vector<Int> sorted_id_to_id_;

    SimplicialFiltration() = default;
    SimplicialFiltration(Simplices&& simplices);
    SimplicialFiltration(const SimplicialFiltration&) = delete;
    SimplicialFiltration(SimplicialFiltration&&) = default;
    SimplicialFiltration& operator=(const SimplicialFiltration&) = delete;
    SimplicialFiltration& operator=(SimplicialFiltration&&) = default;


    void init();

    size_t size() const { return cells_.size(); }

    void sanity_check();

    MatrixWithValues boundary_matrix() const;
    MatrixWithValues coboundary_matrix() const;

    friend std::ostream& operator<<(std::ostream& out, const SimplicialFiltration& fil);

};
}