#include <tuple>

#include "types.h"
#include "profile.h"
#include "matrix.h"
#include "cubical_filtration.h"

namespace cadmus {

    CubicalFiltration::CubicalFiltration(Cubes&& cubes)
            :
            cells_(std::move(cubes))
    {
        init();
    }

    void CubicalFiltration::init()
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

//        print_memory_usage(-1, fmt::format("CubicalFiltration::init, enter, size = {}", cells_.size()));

        // sort cubes by UidValue
        // Sorting is done by dimension
        std::sort(cells_.begin(), cells_.end(),
                [this](const Cube& a, const Cube& b) { return a.uid_value() < b.uid_value(); });

        // assign sorted_id to all cubes
        for(Int sorted_id = 0; sorted_id < cells_.size(); ++sorted_id)
            cells_[sorted_id].set_sorted_id(sorted_id);

        // fill in lookup tables
        id_to_sorted_id_.clear();

#ifndef NDEBUG
        std::unordered_set<Int> sigma_ids;
#endif

        for(const auto& sigma: cells_) {
#ifndef NDEBUG
            assert(sigma_ids.count(sigma.id()) == 0);
            sigma_ids.insert(sigma.id());
#endif
            id_to_sorted_id_[sigma.id()] = sigma.sorted_id();
        }

//        print_memory_usage(-1, fmt::format("CubicalFiltration::init, after filling in lookup tables, size {}", cells_.size()));
        assert(id_to_sorted_id_.size() == cells_.size());
    }

    MatrixWithValues CubicalFiltration::boundary_matrix()
    {
        CALI_CXX_MARK_FUNCTION;
        auto logger = spd::get("console");
        Matrix d;

        // boundary matrix is created columnwise

        d.stores_columns = true;
        d.stores_rows = false;
        d.n_columns = d.n_rows = size();

        d.columns.reserve(size());

        for(auto& sigma: cells_) {
            Column col;

            if (sigma.dim() > 0) {
                col.reserve(1 << sigma.dim());

                for(auto& tau_id: sigma.boundary()) {
                    if (logger) logger->trace("sigma = {}, looking for simplex tau = {}", sigma, tau_id);
                    auto tau_sorted_id = id_to_sorted_id_.at(tau_id);
                    col.push_back(tau_sorted_id);
                }

                std::sort(col.begin(), col.end());
            }

            d.columns.push_back(std::move(col));
        }

        id_to_sorted_id_ = decltype(id_to_sorted_id_)();

//        print_memory_usage(-1, fmt::format("CubicalFiltration::boundary_matrix, after filling d"));

        UidValues row_values, col_values;

        row_values.reserve(size());

        for(const Cube& sigma: cells_) {
            row_values.push_back(sigma.uid_value());
        }

//        print_memory_usage(-1, fmt::format("CubicalFiltration::boundary_matrix, after filling in row_values"));

        // for now, we compute boundary matrix in all dimensions, it is a square matrix
        col_values = row_values;

        return {std::move(d), std::move(col_values), std::move(row_values)};;
    }

    MatrixWithValues CubicalFiltration::coboundary_matrix()
    {
        CALI_CXX_MARK_FUNCTION;

        MatrixWithValues bm = boundary_matrix();

//        print_memory_usage(-1, fmt::format("CubicalFiltration::coboundary_matrix, after getting bm"));

        // column values of anti-transpose are the reverse of row values of D
        UidValues col_values = std::move(*bm.row_values);
        std::reverse(col_values.begin(), col_values.end());

        // row values of anti-transpose are the reverse of column values of D
        // could just copy col_values, but this will break if we want one dimension only
        UidValues row_values = std::move(*bm.column_values);
        std::reverse(row_values.begin(), row_values.end());

//        print_memory_usage(-1, fmt::format("CubicalFiltration::coboundary_matrix, after reverting row_values"));

        bool fill_columns = true, fill_rows = false;

        Matrix d_perp = bm.matrix.antitranspose(fill_columns, fill_rows);

//        print_memory_usage(-1, fmt::format("CubicalFiltration::coboundary_matrix, after getting d_perp"));

        return {std::move(d_perp), std::move(col_values), std::move(row_values)};
    }

    std::ostream& operator<<(std::ostream& out, const CubicalFiltration& fil)
    {
        out << "CubicalFiltration(" << "cubes = [\n";
        for(auto&& sigma: fil.cells_) {
            out << sigma << "\n";
        }
        out << "]\n";

        out << "id_to_sorted_id_ = {\n";
        for(const auto[k, v]: fil.id_to_sorted_id_) {
            out << k << " -> " << v << "\n";
        }
        out << "}\n";

        out << "}\n)";

        return out;
    }
}
