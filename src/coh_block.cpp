#include <chrono>

#include <diy/algorithms.hpp>

#include "coh_block.h"

namespace std {
    std::ostream& operator<<(std::ostream& out, const std::map<int, int>& m)
    {
        out << "map{";
        for(auto&[k, v]: m) {
            out << k << " -> " << v << ", ";
        }
        out << "}";
        return out;
    }
}

namespace cadmus {

    void CohBlock::count_simplices_in_cube()
    {
        auto logger = spd::get("console");

        total_n_cells_per_cube_ = 0;

        for(dim_type d = 0; d <= CADMUS_DIM; ++d) {
            total_n_cells_per_cube_ += get_fr_displacements(d).size();
            if (d == 0)
                n_cells_in_lower_dimension_[d] = 0;
            else
                n_cells_in_lower_dimension_[d] = n_cells_in_lower_dimension_.at(d - 1) + get_fr_displacements(d - 1).size();
        }

        if (logger) logger->trace("total_n_cells_per_cube = {}, n_cells_in_lower_dimension = {}", total_n_cells_per_cube_, n_cells_in_lower_dimension_);
    }

    void CohBlock::init(FilType fil_type, bool negate, bool clearing_opt, bool swap_reduction)
    {
        CALI_CXX_MARK_FUNCTION;
        negate_ = negate;
        stats_.clearing_opt = clearing_opt_ = clearing_opt;
        stats_.swap_reduction = swap_reduction_ = swap_reduction;

        // get all values in
        for(size_t i = 0; i < fab_.size(); ++i) {
            all_values_.push_back(fab_.data()[i]);
        }

//        print_memory_usage(stats_.gid, "init, after copying to all_values_");

        // all indexing at this point is local
        init_coboundary_matrix_local(fil_type);

//        print_memory_usage(stats_.gid, "init, after init_coboundary_matrix");

        get_min_max_cells();

//        print_memory_usage(stats_.gid, "init, after get_min_max_cells");

        // in Release, we don't need filtration after that point,
        // clear it to save memory
#ifdef NDEBUG
        if (fil_type == FilType::Freudenthal)
            freudenthal_fil_ = SimplicialFiltration();
        else if (fil_type == FilType::Cubical)
            cubical_fil_ = CubicalFiltration();
#endif

//        print_memory_usage(stats_.gid, "init, after clearing filtration");
    }


    void CohBlock::get_min_max_cells()
    {
        if (freudenthal_fil_.size() > 0) {
            for(const Simplex& sigma: freudenthal_fil_.cells_) {
                dim_type d = sigma.dim();
                if (min_values_.find(d) == min_values_.end()) {
                    min_values_[d] = sigma.uid_value();
                    max_values_[d] = sigma.uid_value();
                } else {
                    min_values_[d] = std::min(sigma.uid_value(), min_values_[d]);
                    max_values_[d] = std::max(sigma.uid_value(), max_values_[d]);
                }
            }
        } else {
            if (cubical_fil_.size() == 0)
                throw std::runtime_error("get_min_max_cells: both filtrations empty");
             for(const Cube& sigma: cubical_fil_.cells_) {
                dim_type d = sigma.dim();
                if (min_values_.find(d) == min_values_.end()) {
                    min_values_[d] = sigma.uid_value();
                    max_values_[d] = sigma.uid_value();
                } else {
                    min_values_[d] = std::min(sigma.uid_value(), min_values_[d]);
                    max_values_[d] = std::max(sigma.uid_value(), max_values_[d]);
                }
            }
        }

        auto logger = spd::get("console");

        for(const auto&[dim, val]: min_values_) {
            bool is_interior = union_uids_.count(val.uid) == 0;
            logger->info("gid = {}, dim = {}, min simplex value = {}, is_interior = {}", stats_.gid, dim, val, is_interior);
        }

        for(const auto&[dim, val]: max_values_) {
            bool is_interior = union_uids_.count(val.uid) == 0;
            logger->info("gid = {}, dim = {}, max simplex value = {}, is_interior = {}", stats_.gid, dim, val, is_interior);
        }

        global_min_values_ = min_values_;
        global_max_values_ = max_values_;
    }

    void CohBlock::init_coboundary_matrix_local(FilType fil_type)
    {
        CALI_CXX_MARK_FUNCTION;
        Timer timer;

        UidValues uid_values_interior, uid_values_union;

        auto logger = spd::get("console");

        MatrixWithValues d_total;

        if (fil_type == FilType::Cubical) {
            logger->warn("Using cubical");
            cubical_fil_ = get_cubical_filtration(CADMUS_DIM);
            d_total = cubical_fil_.coboundary_matrix();
        } else if (fil_type == FilType::Freudenthal) {
            logger->warn("Using simplicial Freudenthal");
            freudenthal_fil_ = get_freudenthal_filtration(CADMUS_DIM);
            d_total = freudenthal_fil_.coboundary_matrix();
        } else {
            throw std::runtime_error("unknown filtration type");
        }

//        print_memory_usage(stats_.gid, "init_coboundary_matrix_local, after getting d_total");

        if (logger) logger->debug("CohBlock::init, freudenthal_fil.coboundary_matrix OK, {}", d_total.debug_info());
        assert(d_total.sanity_check(true));

        auto is_interior_pred = [&d_total, this](size_t i) {
          return union_uids_.find((*d_total.column_values)[i].uid) == union_uids_.end();
        };

        for(int i = 0; i < d_total.n_columns(); ++i) {
            if (is_interior_pred(i)) {
                uid_values_interior.emplace_back((*d_total.column_values)[i]);
            } else {
                uid_values_union.emplace_back((*d_total.column_values)[i]);
            }
        }

        if (logger) logger->debug("CohBlock::init, uid_values OK");

        // split into 3 parts
        // TODO: refactor as method

        Matrix i_to_i, u_to_i, u_to_u;

        i_to_i.stores_rows = u_to_i.stores_rows = u_to_u.stores_rows = false;
        i_to_i.stores_columns = u_to_i.stores_columns = u_to_u.stores_columns = true;

        std::unordered_map<size_t, size_t> i_lookup, u_lookup;

        {
            size_t i_lookup_index = 0, u_lookup_index = 0;
            for(size_t i = 0; i < d_total.n_rows(); ++i) {
                if (is_interior_pred(i)) {
                    i_lookup[i] = i_lookup_index++;
                } else {
                    u_lookup[i] = u_lookup_index++;
                }
            }
        }

//        print_memory_usage(stats_.gid, "init_coboundary_matrix_local, after getting i_lookup, u_lookup");

        for(size_t col_idx = 0; col_idx < d_total.n_columns(); ++col_idx) {
            Column& current_col = d_total.matrix.columns.at(col_idx);
            if (is_interior_pred(col_idx)) {
                // coboundary of interior simplex is interior
                for(auto& x: current_col) {
                    x = i_lookup.at(x);
                }

                i_to_i.columns.emplace_back(std::move(current_col));
            } else {
                // split coboundary of union simplex in 2 parts
                // insert empty columns
                u_to_i.columns.emplace_back();
                u_to_u.columns.emplace_back();
                for(auto x: current_col) {
                    if (is_interior_pred(x)) {
                        u_to_i.columns.back().push_back(i_lookup.at(x));
                    } else {
                        u_to_u.columns.back().push_back(u_lookup.at(x));
                    }
                }
            }
        }

        i_to_i.n_columns = i_to_i.columns.size();
        u_to_i.n_columns = u_to_i.columns.size();
        i_to_i.n_rows = u_to_i.n_rows = i_lookup.size();

        u_to_u.n_columns = u_to_u.columns.size();
        u_to_u.n_rows = u_lookup.size();

        SPUidValues sp_uid_values_interior = std::make_shared<UidValues>(std::move(uid_values_interior));
        SPUidValues sp_uid_values_union    = std::make_shared<UidValues>(std::move(uid_values_union));

        d_interior_ = MatrixWithValues(std::move(i_to_i), sp_uid_values_interior, sp_uid_values_interior);
        d_union_to_interior_ = MatrixWithValues(std::move(u_to_i), sp_uid_values_union, sp_uid_values_interior);
        d_union_to_union_ = MatrixWithValues(std::move(u_to_u), sp_uid_values_union, sp_uid_values_union);

//        print_memory_usage(stats_.gid, "init_coboundary_matrix_local, after getting d_interior_, d_union_to_interior_, d_union_to_union_");

        assert(d_interior_.sanity_check(true));
        assert(d_union_to_interior_.sanity_check(true));
        assert(d_union_to_union_.sanity_check(true));

        if (logger) logger->debug("d_interior_: {}", d_interior_.debug_info());
        if (logger) logger->debug("d_union_to_interior_: {}", d_union_to_interior_.debug_info());
        if (logger) logger->debug("d_union_to_union_: {}", d_union_to_union_.debug_info());

        // clear union_uids_, not used after this point
        union_uids_ = decltype(union_uids_)();

        stats_.local_init = timer.elapsed_reset();
    }

    void CohBlock::reduce_interior()
    {
        CALI_CXX_MARK_FUNCTION;

        Timer timer;

        auto logger = spd::get("console");

        std::vector<Int> pivots(d_interior_.n_rows(), -1);

        if (logger) logger->debug("started reducing interior matrix, n_rows = {} = {}, n_columns = {} = {}", d_interior_.n_rows(), d_interior_.matrix.rows.size(), d_interior_.n_columns(), d_interior_.matrix.columns.size());

        assert(d_interior_.sanity_check(true));

        for(Int col_idx = 0; col_idx < d_interior_.matrix.n_columns; ++col_idx) {
            Column& reduced_col = d_interior_.matrix.columns[col_idx];
            // TODO: apparent pairs
            // Clearing optimization works: if there is a low in a certain row,
            // we know for sure that the row simplex is positive

            if (clearing_opt_) {
                if (not d_interior_.matrix.columns[col_idx].is_zero()) {
                    // simplex col_idx is pivot -> col_idx is positive -> its column is 0
                    if (pivots[col_idx] >= 0) {
                        logger->debug("reduce_interior: cleared {}", (*d_interior_.column_values)[col_idx]);
                        d_interior_.matrix.columns[col_idx].clear();
                        stats_.interior_cleared++;
                        continue;
                    }
                }
            }

            if (logger) logger->trace("reducing column {}", col_idx);

            load_to_cache(reduced_col);

            while(not addition_cache_.empty()) {

                if (logger) logger->trace("col_idx = {}, low = {}, low value = {}", col_idx, d_interior_.matrix.columns[col_idx].low(), (*d_interior_.row_values)[d_interior_.matrix.columns[col_idx].low()]);

                Int& pivot = pivots[*addition_cache_.rbegin()];

                if (logger) logger->trace("col_idx = {}, low = {}, pivot = {}", col_idx, d_interior_.matrix.columns[col_idx].low(), pivot);

                if (pivot == -1) {
                    pivot = col_idx;
                    if (logger) logger->trace("set pivot to col_idx = {}, pivot = {}", col_idx, pivot);
                    load_from_cache(reduced_col);
                    break;
                } else {
                    add_to_cached(d_interior_.matrix.columns[pivot]);
                    if (logger) logger->trace("AFTER addition column= {}", d_interior_.matrix.columns[col_idx]);
                    if (logger) logger->trace("col_idx = {}, added pivot, now low = {}", col_idx, d_interior_.matrix.columns[col_idx].low());
                }
            } // reduction loop

            if (addition_cache_.empty()) {
                reduced_col.clear();
            }

            if (logger) logger->trace("done reducing column {}", col_idx);
        } // loop over columns

        assert(d_union_to_interior_.sanity_check(true));
        assert(d_union_to_union_.sanity_check(true));
        assert(d_interior_.sanity_check(true));

        if (logger) logger->trace("reduction interior done");

#ifndef NDEBUG
        // debug: print local diagram
        {
            std::unordered_set<size_t> paired;

            for(size_t col_idx = 0; col_idx < d_interior_.n_columns(); ++col_idx) {
                if (d_interior_.matrix.columns[col_idx].is_zero())
                    continue;
                size_t birth_idx = col_idx;
                size_t death_idx = d_interior_.matrix.columns[col_idx].low();

                paired.insert(death_idx);

                auto birth_v = (*d_interior_.column_values)[birth_idx];
                auto death_v = (*d_interior_.column_values)[death_idx];

                if (freudenthal_fil_.size() > 0) {
                    Simplex birth_simplex = *std::find_if(freudenthal_fil_.cells_.begin(), freudenthal_fil_.cells_.end(), [birth_v](const Simplex& x) { return x.id() == birth_v.uid; });
                    Simplex death_simplex = *std::find_if(freudenthal_fil_.cells_.begin(), freudenthal_fil_.cells_.end(), [death_v](const Simplex& x) { return x.id() == death_v.uid; });
                    if (birth_v.value != death_v.value)
                        if (logger) logger->debug("LOCAL PAIRING: ({}, {}), simplices {} <-> {}", birth_v.value, death_v.value, birth_simplex, death_simplex);
                } else {
                    Cube birth_simplex = *std::find_if(cubical_fil_.cells_.begin(), cubical_fil_.cells_.end(), [birth_v](const auto& x) { return x.id() == birth_v.uid; });
                    Cube death_simplex = *std::find_if(cubical_fil_.cells_.begin(), cubical_fil_.cells_.end(), [death_v](const auto& x) { return x.id() == death_v.uid; });
                    if (birth_v.value != death_v.value)
                        if (logger) logger->debug("LOCAL PAIRING: ({}, {}), simplices {} <-> {}", birth_v.value, death_v.value, birth_simplex, death_simplex);
                }

            }

            for(size_t col_idx = 0; col_idx < d_interior_.n_columns(); ++col_idx) {
                if (d_interior_.matrix.columns[col_idx].is_zero() and paired.count(col_idx) == 0) {
                    size_t birth_idx = col_idx;

                    auto birth_v = (*d_interior_.column_values)[birth_idx];
                    if (freudenthal_fil_.size() > 0) {
                        Simplex birth_simplex = *std::find_if(freudenthal_fil_.cells_.begin(), freudenthal_fil_.cells_.end(), [birth_v](const Simplex& x) { return x.id() == birth_v.uid; });
                        if (logger) logger->debug("LOCAL PAIRING: ({}, INF), birth simplex {}", birth_v.value, birth_simplex);
                    } else {
                        Cube birth_cell = *std::find_if(cubical_fil_.cells_.begin(), cubical_fil_.cells_.end(), [birth_v](const auto& x) { return x.id() == birth_v.uid; });
                        if (logger) logger->debug("LOCAL PAIRING: ({}, INF), birth simplex {}", birth_v.value, birth_cell);
                    }
                }
            }
        }
#endif

        stats_.interior_reduction = timer.elapsed_reset();
    }

    void CohBlock::sparsify()
    {
        CALI_CXX_MARK_FUNCTION;

        Timer timer;

        auto logger = spd::get("console");

        if (logger) logger->debug("started sparsification");

        d_interior_.convert_to_rows(true);
        d_union_to_interior_.convert_to_rows(false);

        assert(d_interior_.stores_rows() and d_union_to_interior_.stores_rows());

        std::vector<std::pair<size_t, size_t>> non_zero_cols_with_lows;

        for(size_t col_idx = 0; col_idx < d_interior_.n_columns(); ++col_idx) {
            if (not d_interior_.matrix.columns[col_idx].is_zero() and d_interior_.matrix.columns[col_idx].size() != 1) {
                non_zero_cols_with_lows.emplace_back(col_idx, d_interior_.matrix.columns[col_idx].low());
            }
        }

        std::sort(non_zero_cols_with_lows.begin(), non_zero_cols_with_lows.end(), [](const auto a, const auto b) { return a.second > b.second; });

        // we start with column with the biggest low, so pivot_row always contains only one non-zero
        // thus all the remaining columns remain unchanged and we don't need to update them
        for(auto[column_idx, pivot_row]: non_zero_cols_with_lows) {
            Column& col = d_interior_.matrix.columns[column_idx];

            for(auto reduced_row: col) {
                if (reduced_row < pivot_row) {
                    // here we don't use cache, because there is only one addtion -- it's not reduction
                    // column in d_interior.matrix is updated at once, we know the final answer: only the lowest element survives
                    d_union_to_interior_.matrix.add_to_row(reduced_row, pivot_row, 1);
                }
            }

            // move the lowest entry to front
            col.front() = col.back();
            col.resize(1);

        }

        if (logger) logger->debug("sparsified, filling in sparsifed_d_");

        sparsified_d_.reserve(d_interior_.n_columns());

        for(size_t col_idx = 0; col_idx < d_interior_.matrix.columns.size(); ++col_idx) {
            UidValue col_value = (*d_interior_.column_values)[col_idx];

            if (d_interior_.matrix.columns[col_idx].is_zero()) {
                //sparsified_d_.emplace_back(col_value);
            } else {
                UidValue row_value = (*d_interior_.row_values)[d_interior_.matrix.columns[col_idx].low()];
                sparsified_d_.emplace_back(col_value, row_value);
            }
        }


        // convert u_to_i back to column format
        d_union_to_interior_.matrix.clear_columns();
        d_union_to_interior_.matrix.convert_to_columns(false);

        assert(d_union_to_interior_.sanity_check(true));
        assert(d_union_to_union_.sanity_check(true));

        if (logger) logger->info("sparsify: sparsified_d_.size = {} (pid = {})", sparsified_d_.size(), getpid());

        // sparsify d_union_to_interior matrix
        {
            size_t di_idx = 0;
            size_t total_before {0}, total_after {0};

            std::unordered_set<size_t> to_remove;

            for(size_t dui_idx = 0; dui_idx < d_union_to_interior_.n_columns(); ++dui_idx) {

                Column& col = d_union_to_interior_.matrix.columns[dui_idx];

                while(di_idx < d_interior_.n_columns() and (*d_interior_.column_values)[di_idx] > (*d_union_to_interior_.column_values)[dui_idx]) {
                    if (not d_interior_.matrix.columns[di_idx].is_zero()) {
                        to_remove.insert(d_interior_.matrix.columns[di_idx].low());
                    }
                    di_idx++;
                }

                if (to_remove.empty())
                    continue;

                Column new_col;
                new_col.reserve(col.size());

                std::copy_if(col.begin(), col.end(), std::back_inserter(new_col), [&to_remove](size_t x) { return to_remove.count(x) == 0; });

                total_before += col.size();
                total_after += new_col.size();

                col = std::move(new_col);
            }
            if (logger) logger->info("sparsify: processed d_union_to_interior, for affected columns: total_before = {}, total_after = {} (pid = {})", total_before, total_after, getpid());
        }

        d_interior_.matrix.clear();

        stats_.sparsify = timer.elapsed_reset();

        int ui_columns = 0;
        for(size_t col_idx = 0; col_idx < d_union_to_interior_.n_columns(); ++col_idx)
            if (not d_union_to_interior_.matrix.columns[col_idx].is_zero())
                ui_columns++;

        int uu_columns;
        for(size_t col_idx = 0; col_idx < d_union_to_union_.n_columns(); ++col_idx)
            if (not d_union_to_union_.matrix.columns[col_idx].is_zero())
                uu_columns++;

        if (logger)
            logger->info("sparsify exit: ultrasparse columns: {}, UI columns: {}, non-zero {}, UU columns: {}, non-zero {}", sparsified_d_.size(), d_union_to_interior_.n_columns(), ui_columns, d_union_to_union_.n_columns(), uu_columns);
    }

    bool CohBlock::check_sparsified_d(bool check_sorting)
    {
        auto logger = spd::get("console");

        if (logger) logger->debug("enter check_sparsified_d");

        for(auto se: sparsified_d_) {
            if (not((0 <= se.row_uid_value.dim()) and (se.row_uid_value.dim() <= CADMUS_DIM))) {
                if (logger) logger->critical("check_sparsified_d: bad entry {} large row dim", se);
                return false;
            }

            if (not((0 <= se.col_uid_value.dim()) and (se.col_uid_value.dim() <= CADMUS_DIM))) {
                if (logger) logger->critical("check_sparsified_d: bad entry {} large col dim", se);
                return false;
            }

            if (se.col_uid_value.dim() + 1 != se.row_uid_value.dim()) {
                if (logger) logger->critical("check_sparsified_d: bad entry {} dim mismatch", se);
                return false;
            }

            if (se.col_uid_value.value > se.row_uid_value.value) {
                if (logger) logger->critical("check_sparsified_d: bad entry {} value order", se);
                return false;
            }

        }
        if (check_sorting) {
            if (not std::is_sorted(sparsified_d_.begin(), sparsified_d_.end(), [](const SparseEntry& x, const SparseEntry& y) { return x.col_uid_value > y.col_uid_value; })) {
                if (logger) logger->critical("check_sparsified_d: not sorted");
                return false;
            }
        }

        return true;
    }

    void CohBlock::add_to_cached(const Column& pivot_col)
    {
        for(size_t i: pivot_col) {
            auto insertion_res = addition_cache_.insert(i);
            // if the element was already there, remove it
            if (not insertion_res.second) {
                addition_cache_.erase(insertion_res.first);
            }
        }
    }

    void CohBlock::cached_drop_low()
    {
        addition_cache_.erase(std::prev(addition_cache_.end()));
    }

    void CohBlock::load_to_cache(const Column& col)
    {
        addition_cache_.clear();
        addition_cache_.insert(col.cbegin(), col.cend());
    }

    void CohBlock::load_from_cache(Column& col)
    {
        col.data.assign(addition_cache_.begin(), addition_cache_.end());
    }

    PersistencePairs CohBlock::read_off_persistence()
    {
        PersistencePairs result;
        for(dim_type dim = 0; dim < CADMUS_DIM; ++dim) {
            global_r_[dim].read_off_persistence(negate_, global_min_vertex_uid_, result);
        }
        return result;
    }

    void CohBlock::save_diagrams(std::string fname_out, int gid)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

        bool save_as_text = true;

        auto pers_pairs = read_off_persistence();

        if (pers_pairs.empty()) {
            logger->trace("no persistence points on this rank, exiting save_diagrams");
            return;
        }

        if (save_as_text) {
            fname_out = fname_out + "_" + std::to_string(gid);
            std::ofstream f(fname_out);

            if (not f.good())
                throw std::runtime_error(std::string("cannot open file ") + fname_out);

            f << "dim; birth; death\n";

            for(auto p: pers_pairs) {
                f << p.dim << "; " << p.birth << "; " << p.death << "\n";
            }

            for(auto p: pers_pairs) {
                logger->trace("found pair, column: {}, low: {}", p.column_id, p.low_id);
            }

            f.close();
        } else {
            throw std::runtime_error("not implemented");
        }
    }

    std::unordered_map<int, ValueMatrix> CohBlock::split_matrix_by_gids(MatrixWithValues& m)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");
        std::unordered_map<int, MatrixWithValues> result;

        for(size_t col_idx = 0; col_idx < m.n_columns(); ++col_idx) {

            if (m.matrix.columns[col_idx].is_zero()) {
                if (logger) logger->trace("split_matrix_by_gids, col_idx = {} column is ZERO, skip", col_idx);
                continue;
            }

            int to_gid = gid_by_col_value((*m.column_values)[col_idx]);

            if (logger) logger->trace("split_matrix_by_gids, my_gid = {}, col_val = {}, col_idx = {}, to_gid = {}", stats_.gid,  (*m.column_values)[col_idx], col_idx, to_gid);

            if (!result.count(to_gid)) {
                Matrix empty = Matrix::create_col_matrix({}, m.n_rows());
                SPUidValues empty_col_vals = std::make_shared<UidValues>(UidValues());
                result[to_gid] = MatrixWithValues(std::move(empty), empty_col_vals, m.row_values);
            }

            result[to_gid].append_column(m.matrix.columns[col_idx], (*m.column_values)[col_idx]);
        }

        std::unordered_map<int, ValueMatrix> result_;

        for(auto&[k, v]: result) {
            result_[k] = v;
            if (logger) logger->trace("split_matrix_by_gids, gid = {}, n_columns = {}, final result size = {}", k, v.n_columns(), result_[k].size());
        }

        return result_;
    }

    std::unordered_map<int, UltraSparseMatrixWithValues> CohBlock::split_ultrasparse_matrix_by_gids()
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

        std::unordered_map<int, UltraSparseMatrixWithValues> result;

        for(const SparseEntry& se: sparsified_d_) {
            int to_gid = gid_by_row_value(se.col_uid_value);
            if (logger) logger->trace("split_ultrasparse_matrix_by_gids, my gid = {}, to_gid = {}, col = {}, row = {}", stats_.gid, to_gid, se.col_uid_value, se.row_uid_value);
            result[to_gid].push_back(se);
        }

        sparsified_d_.clear();
        sparsified_d_.shrink_to_fit();

        return result;
    }

    void CohBlock::rearrange_matrix_send(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

        auto gid_to_ui_matrix = split_matrix_by_gids(d_union_to_interior_);
        auto gid_to_uu_matrix = split_matrix_by_gids(d_union_to_union_);

        // ultrasparse part
        auto gid_to_us_matrix = split_ultrasparse_matrix_by_gids();

        std::map<int, RearrangeMessage> gid_to_message;

        for(auto&[gid, ui_matrix]: gid_to_ui_matrix) {
            gid_to_message[gid].has_ui_matrix = true;
            gid_to_message[gid].ui_matrix = std::move(ui_matrix);
        }

        for(auto&[gid, uu_matrix]: gid_to_uu_matrix) {
            gid_to_message[gid].has_uu_matrix = true;
            gid_to_message[gid].uu_matrix = std::move(uu_matrix);
        }

        for(auto&[gid, us_matrix]: gid_to_us_matrix) {
            gid_to_message[gid].has_us_matrix = true;
            gid_to_message[gid].us_matrix = std::move(us_matrix);
        }

        if (logger) logger->info("rearrange_send, gid_to_message.size = {}", gid_to_message.size());

        for(auto&[gid, msg]: gid_to_message) {
            if (gid == cp.gid()) {
                // these columns remain at this gid, do not send them
                // at the receive phase we merge the received matrices into u_to_i/u_local, local_sparse_
                if (logger) logger->trace("rearrange_send, gid = {}, stays here ui.size = {} uu.size ={}, us.size = {}", gid, msg.ui_matrix.size(), msg.uu_matrix.size(), msg.us_matrix.size());
                u_to_i_local_ = std::move(msg.ui_matrix);
                u_to_u_local_ = std::move(msg.uu_matrix);
                local_sparse_ = std::move(msg.us_matrix);
            } else {
                enqueue(cp, assigner, gid, msg);
            }
        }
        if (logger) logger->trace("rearrange_send exit");
    }

    void CohBlock::rearrange_matrix_receive(const diy::Master::ProxyWithLink& cp)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");
        std::vector<int> incoming_gids;

        if (logger) logger->trace("rearrange_matrix_receive start");

        cp.incoming(incoming_gids);

        if (logger) logger->trace("rearrange_matrix_receive incoming.size = {}, ui.size = {} uu.size ={}, us.size = {}", incoming_gids.size(), u_to_i_local_.size(), u_to_u_local_.size(), local_sparse_.size());

        for(int from_gid: incoming_gids) {
            if (cp.incoming(from_gid).size()) {
                RearrangeMessage msg;
                cp.dequeue(from_gid, msg);

                if (logger)
                    logger->info("rearrange_matrix_receive my_gid = {}, from_gid = {}, has_ui = {}, ui.size = {}, has_uu = {}, uu.size = {}", cp.gid(), from_gid, msg.has_ui_matrix, msg.ui_matrix.size(), msg.has_uu_matrix,
                            msg.uu_matrix.size());

                if (msg.has_ui_matrix)
                    u_to_i_local_ = ValueMatrix::merge(u_to_i_local_, msg.ui_matrix);

                if (msg.has_uu_matrix)
                    u_to_u_local_ = ValueMatrix::merge(u_to_u_local_, msg.uu_matrix);

                if (msg.has_us_matrix)
                    local_sparse_ = merge(local_sparse_, msg.us_matrix);

                if (logger) logger->debug("rearrange_matrix_receive, after receiving from_gid = {}, ui.size = {} uu.size ={}, us.size = {}", from_gid, u_to_i_local_.size(), u_to_u_local_.size(), local_sparse_.size());
            }
        }

        if (logger) logger->trace("rearrange_matrix_receive loop done, u_to_u_local = {}, u_to_i_local = {}", u_to_u_local_.size(), u_to_i_local_.size());


        // after we merged all u_to_u and u_to_i, combine them to get final matrix
        global_r_ = ValueMatrix::merge_and_split(u_to_u_local_, u_to_i_local_);
        if (logger) logger->info("rearrange_matrix_receive merge done, size {}", global_r_.size());

        for(dim_type dim = 0; dim < CADMUS_DIM; ++dim) {
            global_r_[dim].gid_ = cp.gid();
        }

        append_sparse(global_r_, local_sparse_);
        if (logger) logger->info("rearrange_matrix_receive append_sparse done, size = {}", global_r_.size());

        // after this point we only need global_r_, local matrices can be cleared
        u_to_u_local_ = ValueMatrix();
        u_to_i_local_ = ValueMatrix();
        local_sparse_ = UltraSparseMatrixWithValues();

        // for clearing optimization
        std::unordered_set<UidValue> positive_uids;

        for(dim_type dim = 0; dim < CADMUS_DIM; ++dim) {
            if (dim > 0)
                global_r_[dim].reduce_local(global_min_vertex_uid_, stats_, &global_r_[dim - 1].local_pivots_, clearing_opt_, dim);
            else
                global_r_[dim].reduce_local(global_min_vertex_uid_, stats_, nullptr, false, dim);
        }

        if (logger) logger->info("rearrange_matrix_receive local reduction done");

    }


    void CohBlock::send_columns(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, dim_type dim, int round)
    {
        CALI_CXX_MARK_FUNCTION;

        timer_total_.reset();

        auto logger = spd::get("console");

        std::map<int, UidValues> gid_to_locations;

        if (round == 0) {
            // iterate over all columns
            for(const auto& [col_val, col] : global_r_[dim]) {
                int to_gid = gid_by_row_value(col.low());
                if (logger) logger->debug("send_columns, round = {}, my gid = {}, to_gid = {}, low = {}", round, cp.gid(), to_gid, col.low());
                if (to_gid != cp.gid()) {
                    // column belongs to another gid
                    gid_to_locations[to_gid].push_back(col_val);
                }
            }
        } else {
            // iterate only over new_lows (filled in in receive_columns)
            for(UidValue new_low : new_lows_) {
                int to_gid = gid_by_row_value(new_low);
                if (logger) logger->debug("send_columns, round = {}, my gid = {}, to_gid = {}, low = {}", round, cp.gid(), to_gid, new_low);
                 if (to_gid != cp.gid()) {
                    // column belongs to another gid
                    gid_to_locations[to_gid].push_back(global_r_[dim].column_by_low(new_low));
                }
            }
        }

        int n_sent = 0;

        // enqueue columns to gids that are going to reduce it and remove them from current rank
        for(const auto&[reducer_gid, locs]: gid_to_locations) {
            if (logger) logger->debug("send_columns, my gid = {}, reducer_gid = {}, lows = {}", cp.gid(), reducer_gid, locs.size());
            auto submatrix = global_r_[dim].extract_submatrix(locs);
            enqueue(cp, assigner, reducer_gid, submatrix);
            n_sent++;
        }

        if (logger) logger->trace("send_columns, my gid = {}, n_sent = {}", cp.gid(), n_sent);

        cp.collectives()->clear();
        cp.all_reduce(n_sent, std::plus<int>());

        double elapsed = timer_total_.elapsed_reset();
        stats_.send_columns_per_dim[dim] += elapsed;
        stats_.final_reduction_per_dim[dim] += elapsed;
    }


    void CohBlock::receive_columns(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, dim_type dim)
    {
        CALI_CXX_MARK_FUNCTION;

        timer_total_.elapsed_reset();

        auto logger = spd::get("console");

        std::vector<int> incoming_gids;

        cp.incoming(incoming_gids);

        UidValue min_dense_col_val;

        bool received = false;

        if (logger) logger->info("receive_columns_ enter, global_r_ in dimension {} size {}", dim, global_r_[dim].size());

        for(int gid: incoming_gids) {
            if (cp.incoming(gid).size()) {
                ValueMatrix rc;
                cp.dequeue(gid, rc);

                if (received) {
                    min_dense_col_val = std::min(rc.first_value(), min_dense_col_val, std::greater());
                } else {
                    min_dense_col_val = rc.first_value();
                    received = true;
                }

                if (logger) logger->debug("receive_columns_, received from gid {}, global_r_ before append {}, rc.size = {}, min_dense_col_val = {}", gid, global_r_.size(), rc.size(), min_dense_col_val);
                global_r_[dim].append_mv(std::move(rc));
                if (logger) logger->debug("receive_columns_, received from gid {}, global_r_ after append {}, rc.size = {}, min_dense_col_val = {}", gid, global_r_.size(), rc.size(), min_dense_col_val);
            }
        }

        if (not received) {
            double elapsed = timer_total_.elapsed_reset();
            new_lows_ = UidValues();
            stats_.final_reduction_per_dim[dim] += elapsed;
            stats_.receive_columns_per_dim[dim] += elapsed;
            return;
        }

        Timer tmr;

        logger->debug("gid = {}, receive_columns, dim = {}, matrix address = {}", stats_.gid, dim, fmt::ptr(&global_r_[dim]));

        new_lows_ = global_r_[dim].reduce_starting_from(min_dense_col_val, global_min_vertex_uid_);

        stats_.reduce_starting_from_per_dim[dim] += tmr.elapsed_reset();

        double elapsed = timer_total_.elapsed_reset();
        stats_.final_reduction_per_dim[dim] += elapsed;
        stats_.receive_columns_per_dim[dim] += elapsed;
    }


    void print_cob(const MatrixWithValues& d, const SimplicialFiltration& fil)
    {
        auto logger = spd::get("console");

        for(size_t col_idx = 0; col_idx < d.n_columns(); ++col_idx) {
            Int uid = (*d.column_values)[col_idx].uid;
            Simplex sigma = *std::find_if(fil.cells_.begin(), fil.cells_.end(), [uid](const Simplex& x) { return x.id() == uid; });
            if (logger) logger->info("------------------");
            if (logger) logger->info("Coboundary of simplex {}:", sigma);
            for(auto row_idx: d.matrix.columns[col_idx]) {
                Int tau_uid = (*d.row_values)[row_idx].uid;
                auto tau = *std::find_if(fil.cells_.begin(), fil.cells_.end(), [tau_uid](const Simplex& x) { return x.id() == tau_uid; });
                if (logger) logger->info("row_idx = {}, tau = {}", row_idx, tau);
            }
            if (logger) logger->info("------------------");
        }
    }

} // namespace cadmus


