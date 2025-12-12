#include <tuple>

#include "types.h"
#include "profile.h"
#include "matrix.h"
#include "filtration.h"

namespace cadmus {

    SimplicialFiltration::SimplicialFiltration(Simplices&& simplices)
            :
            cells_(std::move(simplices))
    {
        init();
    }

    void SimplicialFiltration::init()
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

        // sort simplices by UidValue
        // Sorting is done by dimension
        std::sort(cells_.begin(), cells_.end(),
                [this](const Simplex& a, const Simplex& b) { return a.uid_value() < b.uid_value(); });

        // assign sorted_id to all simplices
        for(Int sorted_id = 0; sorted_id < cells_.size(); ++sorted_id)
            cells_[sorted_id].set_sorted_id(sorted_id);

        // fill in lookup tables
        id_to_sorted_id_.clear();
        vertices_to_sorted_id_.clear();
        sorted_id_to_id_ = std::vector<Int>(size(), k_invalid_id);

        std::unordered_set<Int> sigma_ids;

        for(const auto& sigma: cells_) {
            assert(sigma_ids.count(sigma.id()) == 0);

            sigma_ids.insert(sigma.id());

            id_to_sorted_id_[sigma.id()] = sigma.sorted_id();

            sorted_id_to_id_[sigma.sorted_id()] = sigma.id();

            vertices_to_sorted_id_[sigma.vertices()] = sigma.sorted_id();
        }

        assert(id_to_sorted_id_.size() == cells_.size() and vertices_to_sorted_id_.size() == cells_.size());
    }

    MatrixWithValues SimplicialFiltration::boundary_matrix() const
    {
        CALI_CXX_MARK_FUNCTION;
        auto logger = spd::get("console");
        Matrix d;

        // boundary matrix is created columnwise

        d.stores_columns = true;
        d.stores_rows = false;
        d.n_columns = d.n_rows = size();

        for(auto& sigma: cells_) {
            Column col;

            if (sigma.dim() > 0) {

                for(auto& tau_vertices: sigma.boundary()) {
                    if (logger) logger->trace("sigma = {}, looking for simplex tau = {}", sigma, tau_vertices);
                    auto tau_sorted_id = vertices_to_sorted_id_.at(tau_vertices);
                    col.push_back(tau_sorted_id);
                }

                std::sort(col.begin(), col.end());
            }

            d.columns.push_back(std::move(col));
        }

        UidValues row_values, col_values;

        row_values.reserve(size());

        for(const Simplex& sigma: cells_) {
            row_values.push_back(sigma.uid_value());
        }

        // for now, we compute boundary matrix in all dimensions, it is a square matrix
        col_values = row_values;

        return {std::move(d), std::move(col_values), std::move(row_values)};;
    }

    MatrixWithValues SimplicialFiltration::coboundary_matrix() const
    {
        CALI_CXX_MARK_FUNCTION;

        MatrixWithValues bm = boundary_matrix();

        // column values of anti-transpose are the reverse of row values of D
        UidValues col_values = std::move(*bm.row_values);
        std::reverse(col_values.begin(), col_values.end());

        // row values of anti-transpose are the reverse of column values of D
        // could just copy col_values, but this will break if we want one dimension only
        UidValues row_values = std::move(*bm.column_values);
        std::reverse(row_values.begin(), row_values.end());

        bool fill_columns = true, fill_rows = false;

        Matrix d_perp = bm.matrix.antitranspose(fill_columns, fill_rows);

        return {std::move(d_perp), std::move(col_values), std::move(row_values)};
    }

    std::ostream& operator<<(std::ostream& out, const SimplicialFiltration& fil)
    {
        out << "SimplicialFiltration(" << "simplices = [\n";
        for(auto&& sigma: fil.cells_) {
            out << sigma << "\n";
        }
        out << "]\n";

        out << "id_to_sorted_id_ = {\n";
        for(const auto[k, v]: fil.id_to_sorted_id_) {
            out << k << " -> " << v << "\n";
        }
        out << "}\n";

        out << "sorted_id_to_id_ = [\n";
        for(size_t i = 0; i < fil.sorted_id_to_id_.size(); ++i) {
            out << i << " : " << fil.sorted_id_to_id_[i] << ",\n";
        }
        out << "]\n";

        out << "vertices_to_sorted_id_ = {\n";
        for(auto&&[k, v]: fil.vertices_to_sorted_id_) {
            out << "[";
            for(int vertex_idx = 0; vertex_idx < k.size(); ++vertex_idx) {
                out << k[vertex_idx];
                if (vertex_idx + 1 < k.size())
                    out << ", ";
            }
            out << "] -> " << v << "\n";
        }
        out << "}\n)";

        return out;
    }
}
