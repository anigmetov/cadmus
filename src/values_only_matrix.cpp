#include <algorithm>

#include "profile.h"
#include "values_only_matrix.h"

namespace cadmus {

    std::ostream& operator<<(std::ostream& out, const ValueColumn& col)
    {
        out << "ValueColumn(" << container_to_string(col.data_) << ")";
        return out;
    }

    ValueColumn ValueColumn::merge(ValueColumn& a, ValueColumn& b)
    {
        UidValues merged;

        std::set_union(a.data_.begin(), a.data_.end(), b.data_.begin(), b.data_.end(),
                std::back_inserter(merged), std::greater<>());

        return ValueColumn(std::move(merged));
    }

    bool ValueColumn::sanity_check(ThrowInSanCheck throw_in_san_check)
    {
        if (not std::is_sorted(data_.begin(), data_.end(), std::greater<>())) {
            if (throw_in_san_check == ThrowInSanCheck::yes)
                throw std::runtime_error("unsorted value column");
            else
                return false;
        }
        return true;
    }

    std::set<UidValue, std::greater<>> ValueColumn::addition_cache_;

    ValueMatrix::ValueMatrix(MatrixWithValues& mv)
    {
        if (!mv.stores_columns())
            throw std::runtime_error("not implemented");

        for(size_t col_idx = 0; col_idx < mv.n_columns(); ++col_idx) {

            if (mv.matrix.columns[col_idx].is_zero())
                continue;

            auto col_value = (*mv.column_values)[col_idx];
            ValueColumn& col = columns_[col_value];
            for(size_t entry: mv.matrix.columns[col_idx]) {
                col.push_back((*mv.row_values)[entry]);
            }

            mv.matrix.columns[col_idx].clear();
            mv.matrix.columns[col_idx].shrink_to_fit();

        }
    }

    ValueMatrix::DataType::iterator ValueMatrix::erase(ValueMatrix::DataType::iterator iter)
    {
        auto logger = spd::get("console");

        if (not iter->second.is_zero()) {
//            UidValue low_val = iter->second.low();
//            UidValue col_val = iter->first;
            size_t lp_erased = local_pivots_.erase(iter->second.low());
            logger->trace("ValueMatrix::erase by iter, low = {}, erased from local_pivots: {}", iter->second.low(), lp_erased);
        }

        auto result = columns_.erase(iter);

        logger->trace("ValueMatrix::erase by iter, column value = {}, result == end: {}", iter->first, result == columns_.end());

        return result;
    }

    ValueMatrix::DataType::size_type ValueMatrix::erase(ValueMatrix::DataType::key_type col_val, bool not_found_ok)
    {
        auto logger = spd::get("console");

        if (not_found_ok and columns_.count(col_val) == 0)
            return 0;

        if (not columns_.at(col_val).is_zero()) {
            local_pivots_.erase(columns_.at(col_val).low());
        }

        return columns_.erase(col_val);
    }

    ValueMatrix::DataType::size_type ValueMatrix::erase(ValueMatrix::DataType::key_type col_val, UidValue low_val)
    {
        auto logger = spd::get("console");

        size_t lp_erased = local_pivots_.erase(low_val);
        size_t col_erased = columns_.erase(col_val);

        logger->trace("ValueMatrix::erase: col_val = {}, low_val = {}, erased from local_pivots: {}, erased from columns: {}", col_val, low_val, lp_erased, col_erased);

        return col_erased;
    }

    ValueMatrix ValueMatrix::merge(ValueMatrix& a, ValueMatrix& b)
    {
        ValueMatrix result;

        for(auto&[col_val, a_column]: a.columns_) {
            auto b_iter = b.columns_.find(col_val);
            if (b_iter == b.columns_.end())
                result.columns_[col_val] = std::move(a_column);
            else
                result.columns_[col_val] = ValueColumn::merge(a_column, b_iter->second);
        }

        for(auto&[col_val, b_column]: b.columns_) {
            // if this column value was in a, we have already processed it
            if (a.columns_.count(col_val))
                continue;
            assert(result.columns_.count(col_val) == 0);
            result.columns_[col_val] = std::move(b_column);
        }

        return result;
    }

    ValueMatrices ValueMatrix::merge_and_split(ValueMatrix& a, ValueMatrix& b)
    {
        ValueMatrices result;

        for(auto&[col_val, a_column]: a.columns_) {
            auto b_iter = b.columns_.find(col_val);
            if (b_iter == b.columns_.end())
                result[col_val.dim()].columns_[col_val] = std::move(a_column);
            else
                result[col_val.dim()].columns_[col_val] = ValueColumn::merge(a_column, b_iter->second);
        }

        for(auto&[col_val, b_column]: b.columns_) {
            // if this column value was in a, we have already processed it
            if (a.columns_.count(col_val))
                continue;
            assert(result[col_val.dim()].columns_.count(col_val) == 0);
            result[col_val.dim()].columns_[col_val] = std::move(b_column);
        }

        return result;
    }

    // no check if matrix is reduced
    void ValueMatrix::read_off_persistence(bool negate, Int global_min_vertex_uid, PersistencePairs& result) const
    {
        auto logger = spd::get("console");

        for(auto[col_val, col]: columns_) {
            UidValue low_val = col.low();

            if (col_val.value == low_val.value) {
                continue;
            }

            if (col_val.uid == global_min_vertex_uid) {
                // we use cohomology only: dimension comes from column, like birth value
//                if (negate)
//                    result.emplace_back(col_val.dim, -col_val.value, -std::numeric_limits<Real>::infinity(), col_val, low_val);
//                else
//                    result.emplace_back(col_val.dim, col_val.value, std::numeric_limits<Real>::infinity(), col_val, low_val);
                continue;
            } else {
                // we use cohomology only: dimension comes from column, like birth value
                if (negate)
                    result.emplace_back(col_val.dim(), -col_val.value, -low_val.value, col_val, low_val);
                else
                    result.emplace_back(col_val.dim(), col_val.value, low_val.value, col_val, low_val);
            }
        }
    }

    void append_sparse(ValueMatrices& matrices, const UltraSparseMatrixWithValues& sparse_)
    {
        for(auto[col_val, row_val]: sparse_) {
            assert(matrices[col_val.dim()].columns_.count(col_val) == 0);
            matrices[col_val.dim()].columns_[col_val] = ValueColumn({row_val});
        }
    }

    void ValueMatrix::append_mv(ValueMatrix&& m)
    {
        for(auto&&[col_val, col]: m) {
#ifndef NDEBUG
            if (columns_.count(col_val)) {
                std::cerr << "ERROR IN append_mv, col_val = " << col_val << std::endl;
            }
#endif
                assert(columns_.count(col_val) == 0);
//            assert(col_val > columns_.rbegin()->first);
            columns_[col_val] = std::move(col);
        }
    }

    void ValueColumn::sparsify_reduced()
    {
        if (data_.size() > 1) {
            data_[0] = data_.back();
            data_.resize(1);
            data_.shrink_to_fit();
        }
    }

    void ValueMatrix::sparsify_reduced()
    {
        for(auto& [col_val, col] : columns_) {
            col.sparsify_reduced();
        }
    }

    UidValues ValueMatrix::reduce_starting_from(UidValue start_val, Int skip_uid)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");

        Timer timer;

        UidValues new_lows;

        if (logger) logger->trace("reduce_starting_from, gid = {}, size {}, address = {}", gid_, size(), fmt::ptr(this));

        assert(columns_.count(start_val) == 1);

#ifdef CADMUS_PER_COL_STATS
        std::map<double, UidValue, std::greater<>> top_cols;
        std::map<int, UidValue, std::greater<>> top_cols_additions;
#endif

        std::unordered_set<UidValue> to_delete;
        for(auto col_iter = columns_.find(start_val); col_iter != columns_.end();) {

            timer.reset();

            auto next_col_iter = std::next(col_iter);

            UidValue col_val = col_iter->first;

            if (col_val.uid == skip_uid) {
                col_iter = next_col_iter;
                continue;
            }

            ValueColumn* reduced = &col_iter->second;

            reduced->load_to_cache();

            bool seen = false;
            int n_adds = 0;

            while(not ValueColumn::is_cache_zero()) {
                UidValue low_val = ValueColumn::cache_low();
                auto pivot_iter = local_pivots_.find(low_val);
                if (pivot_iter != local_pivots_.end()) {
                    UidValue pivot_loc = pivot_iter->second;
                    assert(pivot_loc.dim() == col_val.dim());
                    if (pivot_loc > col_val) {
                        logger->trace("reduce_starting_from, pivot to the left, low = {}, pivot_loc = {}, col_val = {}", low_val, pivot_loc, col_val);
                        // pivot to the left from reduced
//                        if (columns_[pivot_loc].is_zero()) {
//                            throw std::runtime_error("bad pivot: is zero");
//                        }
                        n_adds++;
                        ValueColumn::add_to_cache(columns_[pivot_loc]);
                        assert(ValueColumn::is_cache_zero() or ValueColumn::cache_low() > low_val);
                    } else if (pivot_loc < col_val) {
                        // pivot to the right, mark reduced column as pivot and jump right
#ifdef CADMUS_PER_COL_STATS
                        double elapsed = timer.elapsed_reset();
                        if (elapsed > 0.1) {
                            top_cols[elapsed] = col_val;
                            top_cols_additions[n_adds] = col_val;
                        }
#endif
                        logger->trace("reduce_starting_from, pivot to the right, setting col_val = {} as pivot for low = {}, will jump right to reduce pivot_loc = {}", col_val, low_val, pivot_loc);
                        reduced->load_from_cache();
                        set_local_pivot(low_val, col_val);
                        col_iter = columns_.find(pivot_loc);
                        col_val = col_iter->first;
                        // write current column back
                        // make reduced point to the right pivot and load it to cache
                        reduced = &col_iter->second;
                        reduced->load_to_cache();
                    } else {
                        // pivot_loc == col_val: this can happen if we already reduced this column by right jump
                        // it's fine, the column has been reduced, just go to next one
                        // set seen flag, so that we don't include it in new_lows again
                        logger->trace("reduce_starting_from,  col_val = {}  = pivot_loc = {}, setting seen to true and breaking", col_val, pivot_loc);
                        seen = true;
                        break;
                    }
                } else {
#ifdef CADMUS_PER_COL_STATS
                    double elapsed = timer.elapsed_reset();
                    if (elapsed > 1) {
                        top_cols[elapsed] = col_val;
                        top_cols_additions[n_adds] = col_val;
                    }
#endif
                    break;
                }
            }

            reduced->load_from_cache();

            logger->trace("reduce_starting_from,  col_val = {}  after reduction loop", col_val);
            if (reduced->is_zero()) {
                // do not delete now: we may invalidate next_col_iter
                to_delete.insert(col_val);
                assert(local_pivots_.count(col_val) == 0);
//                local_pivots_.erase(col_val);
                logger->trace("reduce_starting_from, marked col_val = {} to delete, ", col_val);
            } else if (not seen) {
                assert(std::find(new_lows.begin(), new_lows.end(), reduced->low()) == new_lows.end());
                logger->trace("reduce_starting_from, adding to new_lows {}, setting col_val {} as pivot", reduced->low(), col_val);
                new_lows.push_back(reduced->low());
                set_local_pivot(reduced->low(), col_val);
            }

            // finally, go to the next column saved in next_col_iter
            // do not advance col_iter, since we may have jumped left
            col_iter = next_col_iter;
        }

        // delete all zeroed columns
        for(UidValue cv: to_delete) {
            logger->trace("reduce_starting_from, erasing column {}", cv);
            erase(cv);
        }

        if (logger) logger->info("reduce_starting_from, gid = {}, new_lows.size = {}", gid_, new_lows.size());

#ifdef CADMUS_PER_COL_STATS
        {
            double total_top_elapsed = 0;
            int n_printed = 0;
            for(auto&[elapsed, col]: top_cols) {
                total_top_elapsed += elapsed;
                if (n_printed < 10)
                    logger->warn("REDUCE_STARTING_FROM TOP column {} took {} sec to reduce", col, elapsed);
                n_printed++;
            }
            if (total_top_elapsed > 10)
                logger->warn("REDUCE_STARTING_FROM TOP total_top_elapsed = {}", total_top_elapsed);
        }

        {
            long int total_top_additions = 0;
            int n_printed = 0;
            for(auto&[adds, col]: top_cols_additions) {
                total_top_additions += adds;
                if (n_printed < 10)
                    logger->warn("REDUCE_STARTING_FROM TOP ADDS column {} took {} additions to reduce", col, adds);
                n_printed++;
            }
            logger->warn("REDUCE_STARTING_FROM TOP total_top_additions = {}", total_top_additions);
        }
#endif

        assert(is_reduced(skip_uid));
        assert(is_local_pivots_correct(skip_uid));

        return new_lows;
    }

    void ValueColumn::load_to_cache() const
    {
        addition_cache_.clear();
        addition_cache_.insert(data_.begin(), data_.end());
    }

    void ValueColumn::load_from_cache()
    {
        data_.assign(addition_cache_.begin(), addition_cache_.end());
    }

    void ValueColumn::add_to_cache(const ValueColumn& pivot_col)
    {
        assert(not is_cache_zero() and not pivot_col.is_zero() and cache_low() == pivot_col.low());

        for(UidValue v : pivot_col) {
            auto insertion_res = addition_cache_.insert(v);
            // if the element was already there, remove it
            if (not insertion_res.second) {
                addition_cache_.erase(insertion_res.first);
            }
        }
    }

    UidValue ValueColumn::cache_low()
    {
        return *addition_cache_.rbegin();
    }

    bool ValueColumn::is_cache_zero()
    {
        return addition_cache_.empty();
    }

    void ValueMatrix::reduce_local(Int skip_uid, Stats& stats, PivotsType* prev_dim_pivots, bool clearing, dim_type dim)
    {
        CALI_CXX_MARK_FUNCTION;

        Timer timer;

        auto logger = spd::get("console");

        if (columns_.empty()) {
            if (logger) logger->trace("reduce_local, empty matrix, return");
            return;
        }

        if (logger) logger->trace("reduce_local, size {}", size());

#ifdef CADMUS_PER_COL_STATS
        Timer column_timer;
        Timer addition_timer;
        double elapsed_in_add = 0;

        std::vector<std::tuple<UidValue, double, int, size_t, size_t>> elapsed_per_col;
        elapsed_per_col.reserve(columns_.size());
#endif

        for(auto col_iter = columns_.begin(); col_iter != columns_.end();) {
#ifdef CADMUS_PER_COL_STATS
            column_timer.reset();
#endif

            int additions = 0;

            UidValue col_val = col_iter->first;

            logger->trace("ValueMatrix::reduce_local gid = {}, started on column {} of size {}, skip uid = {}", stats.gid, col_val, col_iter->second.size(), skip_uid);

            if (col_val.uid == skip_uid) {
                col_iter++;
                continue;
            }

            ValueColumn& reduced = col_iter->second;
            const size_t orig_size = reduced.size();

            if (clearing) {
                if (prev_dim_pivots and prev_dim_pivots->count(col_iter->first)) {
                    logger->trace("ValueMatrix::reduce_local gid = {}, cleared low {} column {}", stats.gid, reduced.low(), col_iter->first);
                    stats.local_cleared++;
                    // do not call ValueMatrix::erase here: it will attempt to remove the entry from local_pivots_
                    // we don't need that: we have not processed this column yet,
                    // so it is not in local_pivots_
                    assert(local_pivots_.count(reduced.low()) == 0 or local_pivots_.at(reduced.low()) != col_val);
                    col_iter = columns_.erase(col_iter);
                    continue;
                }
            }

            reduced.load_to_cache();

            while(not ValueColumn::is_cache_zero()) {
                UidValue low_val = ValueColumn::cache_low();
                logger->trace("ValueMatrix::reduce_local gid = {}, reduction loop for col_val {}, low = {}", stats.gid, col_val, low_val);
                auto pivot_iter = local_pivots_.find(low_val);
                if (pivot_iter != local_pivots_.end()) {

                    UidValue pivot_loc = pivot_iter->second;

                    assert(pivot_loc > col_val);
#ifdef CADMUS_PER_COL_STATS
                    addition_timer.reset();
#endif

                    logger->trace("ValueMatrix::reduce_local gid = {}, low_val = {}, pivot_loc = {}, col_val = {}", stats.gid, low_val, pivot_loc, col_val);
                    ValueColumn::add_to_cache(columns_[pivot_loc]);

#ifdef CADMUS_PER_COL_STATS
                    elapsed_in_add += addition_timer.elapsed();
#endif
                    additions++;

                } else {
#ifdef CADMUS_PER_COL_STATS
                    double elapsed_sec = column_timer.elapsed_reset();
                    elapsed_per_col.emplace_back(col_val, elapsed_sec, additions, ValueColumn::addition_cache_.size(), orig_size);
#endif
                    logger->trace("ValueMatrix::reduce_local gid = {}, low_val = {}, done reducing this column", stats.gid, low_val);
                    break;
                }
            }

            reduced.load_from_cache();

            if (reduced.is_zero()) {
                logger->debug("ValueMatrix::reduce_local gid = {}, column {} reduced to zero", stats.gid, col_val);
                col_iter = erase(col_iter);
            } else {
                logger->trace("ValueMatrix::reduce_local gid ={}, setting column = {} as pivot for {}", stats.gid, col_val, reduced.low());
                set_local_pivot(reduced.low(), col_val);
                col_iter++;
            }
        }

#ifdef CADMUS_PER_COL_STATS
        {
            std::sort(elapsed_per_col.begin(), elapsed_per_col.end(), [](auto x, auto y) { return std::get<1>(x) > std::get<1>(y); });
            double total_col_elapsed = std::accumulate(elapsed_per_col.begin(), elapsed_per_col.end(), 0.0, [](auto x, auto y) { return x + std::get<1>(y); });
            logger->warn("REDUCE_LOCAL TOP total_col_elapsed = {}, elapsed in additions = {}", total_col_elapsed, elapsed_in_add);
            if (elapsed_in_add > 500) {
                int n_printed = 0;
                for(auto[col_val, elapsed, additions, final_size, orig_size]: elapsed_per_col) {
                    if (n_printed < 300) {
                        logger->warn("REDUCE_LOCAL TOP column {} took {} sec to reduce, {} additions, final size {}, original size {}", col_val, elapsed, additions, final_size, orig_size);
                        n_printed++;
                    } else {
                        break;
                    }
                }
            }
        }
#endif

        stats.local_reduction_per_dim[dim] += timer.elapsed();

        assert(is_reduced(skip_uid));
        assert(is_local_pivots_correct(skip_uid));
    }

    bool ValueMatrix::is_reduced(Int skip_uid) const
    {
        std::unordered_set<UidValue> lows;
        for(const auto& col_iter : columns_) {
            if (col_iter.first.uid == skip_uid)
                continue;
            if (not col_iter.second.is_zero()) {
                auto ins_res = lows.insert(col_iter.second.low());
                if (not ins_res.second)
                    return false;
            }
        }
        return true;
    }

    bool ValueMatrix::is_local_pivots_correct(Int skip_uid) const
    {
        auto logger = spd::get("console");
        size_t n_non_zero_cols = 0;
        for(const auto& col_iter : columns_) {
            if (col_iter.first.uid == skip_uid)
                continue;
            if (not col_iter.second.is_zero()) {
                n_non_zero_cols++;
                if (local_pivots_.count(col_iter.second.low()) == 0) {
                    logger->critical("local_pivots_ is incorrect, gid = {}, low = {}, col_val = {}, NO entry in local_pivots", gid_, col_iter.second.low(), col_iter.first);
                    return false;
                }
                if (local_pivots_.at(col_iter.second.low()) != col_iter.first) {
                    logger->critical("local_pivots_ is incorrect, gid = {}, low = {}, col_val = {}, entry in local_pivots: {}", gid_, col_iter.second.low(), col_iter.first, local_pivots_.at(col_iter.second.low()));
                    return false;
                }
            }
        }
        if (n_non_zero_cols != local_pivots_.size()) {
            logger->critical("local_pivots_ is incorrect, gid = {}, n_non_zero_cols = {}, local_pivots_.size = {}", gid_, n_non_zero_cols, local_pivots_.size());
        }
        return n_non_zero_cols == local_pivots_.size();
    }

    UidValue ValueMatrix::column_by_low(const UidValue& low)
    {
        auto iter = local_pivots_.find(low);
#ifndef NDEBUG
        if (iter  == local_pivots_.end()) {
            auto logger = spd::get("console");
            logger->critical("not found in pivots: low = {}", low);
        }
#endif

        assert(iter != local_pivots_.end());
        return iter->second;
    }

    ValueMatrix ValueMatrix::extract_submatrix(const UidValues& col_vals)
    {
        CALI_CXX_MARK_FUNCTION;

        ValueMatrix submatrix;

        for(UidValue col_val: col_vals) {
            UidValue low_val = columns_.at(col_val).low();
            submatrix.columns_[col_val] = std::move(columns_.at(col_val));
            // after move, column in columns_ is zero, but we still must delete from local_pivots_
            erase(col_val, low_val);
        }

        return submatrix;
    }

    void ValueMatrix::set_local_pivot(const UidValue& low_val, const UidValue& col_val)
    {
        auto logger = spd::get("console");

        assert(low_val.dim() == col_val.dim() + 1);
        assert(columns_.at(col_val).low() == low_val);
        // if we jump to reduce the column right from us, it is already set as pivot, so count can be 1
        // but in this case the new pivot must be to the left
        assert(local_pivots_.count(low_val) == 0 or local_pivots_.at(low_val) < col_val);

        logger->trace("set_local_pivot: low_val = {}, col_val = {}", low_val, col_val);

        local_pivots_[low_val] = col_val;
    }

    size_t ValueMatrix::total_column_size_in_bytes() const
    {
        size_t result = 0;

        for(const auto&[col_val, col]: columns_)
            result += col.size();

        result *= sizeof(UidValue);

        return result;
    }

    size_t ValueMatrix::total_column_size_in_megabytes() const
    {
        return total_column_size_in_bytes() / 1024 / 1024;
    }

    size_t ValueMatrices::size() const
    {
        return std::accumulate(matrices_.begin(), matrices_.end(), 0, [](auto s, const auto& m) { return s + m.size(); });
    }

    size_t ValueMatrices::total_column_size_in_megabytes() const
    {
        return std::accumulate(matrices_.begin(), matrices_.end(), 0, [](auto s, const auto& m) { return s + m.total_column_size_in_megabytes(); });
    }

    void ValueMatrices::append_mv(ValueMatrices&& m)
    {
        for(size_t dim = 0; dim < CADMUS_DIM; ++dim) {
            matrices_[dim].append_mv(std::move(m.matrices_[dim]));
        }
    }

} // namespace cadmus

namespace diy {
    void Serialization<cadmus::ValueColumn>::save(BinaryBuffer& bb, const cadmus::ValueColumn& col)
    {
        diy::save(bb, col.data_);
    }

    void Serialization<cadmus::ValueColumn>::load(BinaryBuffer& bb, cadmus::ValueColumn& col)
    {
        diy::load(bb, col.data_);
    }

    void Serialization<cadmus::ValueMatrix>::save(BinaryBuffer& bb, const cadmus::ValueMatrix& matrix)
    {
        auto logger = spd::get("console");
//        diy::save(bb, matrix.columns_);
        diy::save(bb, matrix.columns_.size());
        for(const auto&[col_val, col]: matrix.columns_) {
            diy::save(bb, col_val);
            diy::save(bb, col);
        }
    }

    void Serialization<cadmus::ValueMatrix>::load(BinaryBuffer& bb, cadmus::ValueMatrix& matrix)
    {
//        diy::load(bb, matrix.columns_);
        matrix.columns_.clear();
        size_t sz;
        diy::load(bb, sz);
        for(size_t i = 0; i < sz; ++i) {
            cadmus::UidValue col_val;
            cadmus::ValueColumn col;
            diy::load(bb, col_val);
            diy::load(bb, col);
            matrix.columns_[col_val] = std::move(col);
        }
    }

    void Serialization<cadmus::ValueMatrices>::save(BinaryBuffer& bb, const cadmus::ValueMatrices& m)
    {
        for(cadmus::dim_type dim = 0; dim < CADMUS_DIM; ++dim)
            diy::save(bb, m.matrices_[dim]);
    }

    void Serialization<cadmus::ValueMatrices>::load(BinaryBuffer& bb, cadmus::ValueMatrices& m)
    {
        for(cadmus::dim_type dim = 0; dim < CADMUS_DIM; ++dim)
            diy::load(bb, m.matrices_[dim]);
    }

}