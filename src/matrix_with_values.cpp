#include <iterator>
#include <algorithm>

#include "profile.h"
#include "matrix_with_values.h"

namespace cadmus {

    std::ostream& operator<<(std::ostream& out, const UltraSparseMatrixWithValues& m)
    {
        out << "UltraSparseMatrixWithValues[";
        for(size_t i = 0; i < m.size(); ++i) {
            out << m[i];
            if (i + 1 < m.size())
                out << ", ";
        }
        out << "]";
        return out;
    }

// assume same row range
    UltraSparseMatrixWithValues merge(const UltraSparseMatrixWithValues& a, const UltraSparseMatrixWithValues& b)
    {
        UltraSparseMatrixWithValues result;
        result.reserve(a.size() + b.size());
        // we can use either merge or set_union here: interior blocks do not intersect,
        // so each cell appears only once
        std::merge(a.cbegin(), a.cend(), b.cbegin(), b.cend(), std::back_inserter(result), std::greater<SparseEntry>());

#ifndef NDEBUG
        UltraSparseMatrixWithValues result_1;
        result_1.reserve(a.size() + b.size());
        std::set_union(a.cbegin(), a.cend(), b.cbegin(), b.cend(), std::back_inserter(result_1), std::greater<SparseEntry>());
        assert(result == result_1);
#endif

        return result;
    }

    bool MatrixWithValues::sanity_check(bool cohomology)
    {
        auto logger = spd::get("console");
        matrix.sanity_check();

        if (n_rows() != row_values->size()) {
            if (logger) logger->critical("sanity check: bad row_values, {}", debug_info());
            throw std::runtime_error(debug_info());
            return false;
        }

        if (n_columns() != column_values->size()) {
            if (logger) logger->critical("sanity check: bad column_values, {}", debug_info());
            throw std::runtime_error(debug_info());
            return false;
        }

        if (cohomology) {
            if (not std::is_sorted(row_values->rbegin(), row_values->rend())) {
                if (logger) logger->critical("sanity check: unsorted row_values, {}", debug_info());

                for(int i = 0; i < row_values->size() - 1; ++i) {
                    if ((*row_values)[i] <= (*row_values)[i + 1]) { ;
                        if (logger) logger->critical("sanity check: {}, {}, {}", i, (*row_values)[i], (*row_values)[i + 1]);
                    }
                }

                throw std::runtime_error(debug_info());
                return false;
            }
        } else {
            if (not std::is_sorted(row_values->begin(), row_values->end())) {
                if (logger) logger->critical("sanity check: unsorted row_values, {}", debug_info());
                throw std::runtime_error(debug_info());
                return false;
            }
        }

        if (cohomology) {
            if (not std::is_sorted(column_values->rbegin(), column_values->rend())) {
                if (logger) logger->critical("sanity check: unsorted column_values, {}", debug_info());
                throw std::runtime_error(debug_info());
                return false;
            }
        } else {
            if (not std::is_sorted(column_values->begin(), column_values->end())) {
                if (logger) logger->critical("sanity check: unsorted column_values, {}", debug_info());
                throw std::runtime_error(debug_info());
                return false;
            }
        }

        if (stores_columns()) {
            for(size_t col_idx = 0; col_idx < n_columns(); ++col_idx) {
                dim_type dim_col = column_values->at(col_idx).dim();
                for(auto& row_idx: matrix.columns.at(col_idx)) {
                    dim_type dim_row = row_values->at(row_idx).dim();

                    if (cohomology and (dim_row != dim_col + 1) or (not cohomology and (dim_row != dim_col - 1))) {
                        if (logger) logger->critical("sanity check: dimension mismatch, cohomology = {}, {}, col_idx = {}, row_idx = {}, {}, {}", cohomology, debug_info(), col_idx, row_idx, (*column_values)[col_idx], (*row_values)[row_idx]);
                        throw std::runtime_error(debug_info());
                        return false;
                    }
                }
            }
        }

        if (stores_rows()) {
            for(size_t row_idx = 0; row_idx < n_rows(); ++row_idx) {
                dim_type dim_row = row_values->at(row_idx).dim();
                for(auto& col_idx: matrix.columns.at(row_idx)) {
                    dim_type dim_col = column_values->at(col_idx).dim();
                    if (cohomology and (dim_row != dim_col + 1) or (not cohomology and (dim_row != dim_col - 1))) {
                        logger->critical("sanity check: dimension mismatch, cohomology = {}, {}, col_idx = {}, row_idx = {}, {}, {}", cohomology, debug_info(), col_idx, row_idx, (*column_values)[col_idx], (*row_values)[row_idx]);
                        throw std::runtime_error(debug_info());
                        return false;
                    }
                }
            }
        }
        return true;
    }

    MatrixWithValues MatrixWithValues::merge(MatrixWithValues& a, MatrixWithValues& b)
    {
        if (a.stores_rows()) {
            // we only need to merge matrices that are either in row or column order
            if (a.stores_columns() or not b.stores_rows() or b.stores_columns())
                throw std::runtime_error("not implemented");

            return merge_row_matrix(a, b);
        }

        if (a.stores_columns()) {
            // we only need to merge matrices that are either in row or column order
            if (a.stores_rows() or not b.stores_columns() or b.stores_rows())
                throw std::runtime_error("not implemented");

            return merge_column_matrix(a, b);
        }

        throw std::runtime_error("bad parameters");
    }

    void MatrixWithValues::append_column(Column& col, Value val)
    {
        assert(matrix.stores_columns and not matrix.stores_rows);
        matrix.append_column(col);
        column_values->push_back(val);
    }

    void apply_lookup(Column& col, const std::vector<size_t>& lookup_table)
    {
        for(auto& x: col) {
            assert(x < lookup_table.size());
            x = lookup_table[x];
        }
    }

    MatrixWithValues MatrixWithValues::merge_column_matrix(MatrixWithValues& a, MatrixWithValues& b)
    {
        CALI_CXX_MARK_FUNCTION;

        auto logger = spd::get("console");
        logger->debug("in merge_column_matrix, a = {}, b = {}", a.debug_info(), b.debug_info());

        Values union_row_values;
        union_row_values.reserve(a.n_rows() + b.n_rows());

        std::set_union(a.row_values->cbegin(), a.row_values->cend(), b.row_values->cbegin(), b.row_values->cend(), std::back_inserter(union_row_values), std::greater<UidValue>());

        assert(std::is_sorted(union_row_values.begin(), union_row_values.end(), std::greater<>()));

        // index i in columns of a becomes a_lookup_table[i] in the merged matrix, same for b
        std::vector<size_t> a_lookup_table, b_lookup_table;

        size_t a_idx = 0, b_idx = 0;

        for(size_t union_idx = 0; union_idx < union_row_values.size(); ++union_idx) {
            assert(a_idx == a.row_values->size() or a.row_values->at(a_idx) <= union_row_values.at(union_idx));

            if (a_idx < a.row_values->size() and (*a.row_values)[a_idx] == union_row_values[union_idx]) {
                assert(a_lookup_table.size() == a_idx);

                a_lookup_table.push_back(union_idx);

                a_idx++;
            }

            assert(b_idx == b.row_values->size() or b.row_values->at(b_idx) <= union_row_values.at(union_idx));

            if (b_idx < b.row_values->size() and (*b.row_values)[b_idx] == union_row_values[union_idx]) {

                assert(b_lookup_table.size() == b_idx);

                b_lookup_table.push_back(union_idx);

                b_idx++;
            }
        }

#ifdef NDEBUG
        // row_values of a and b are used only in asserts below,
        // free their memory
        a.row_values->clear();
        a.row_values->shrink_to_fit();
        b.row_values->clear();
        b.row_values->shrink_to_fit();
#endif


        MatrixWithValues result;
        result.matrix.n_rows = union_row_values.size();
        result.row_values = std::make_shared<UidValues>(std::move(union_row_values));
        result.matrix.stores_rows = false;
        result.matrix.stores_columns = true;

        size_t a_col_idx = 0, b_col_idx = 0;

        while(a_col_idx < a.n_columns() && b_col_idx < b.n_columns()) {
            if ((*a.column_values)[a_col_idx] > (*b.column_values)[b_col_idx]) {
                // value in A is greater -> it comes first in cohomology, take A column and adjust indices
                result.append_column(a.matrix.columns[a_col_idx], (*a.column_values)[a_col_idx]);

                apply_lookup(result.matrix.columns.back(), a_lookup_table);

                // debug: verify that dimensions and values match
                assert(std::all_of(result.matrix.columns.back().begin(), result.matrix.columns.back().end(),
                        [&result, &a, &a_lookup_table, a_col_idx](auto x) {
                          if (result.row_values->at(x).dim() != a.column_values->at(a_col_idx).dim() + 1)
                              return false;
                          return result.row_values->at(x) == a.row_values->at(std::find(a_lookup_table.begin(), a_lookup_table.end(), x) - a_lookup_table.begin());
                        }));

                // point to the next column in A
                a_col_idx++;
            } else if ((*a.column_values)[a_col_idx] < (*b.column_values)[b_col_idx]) {
                // value in B is greater -> it comes first in cohomology, take B column and adjust indices
                result.append_column(b.matrix.columns[b_col_idx], (*b.column_values)[b_col_idx]);
                apply_lookup(result.matrix.columns.back(), b_lookup_table);

                // debug: verify that dimensions and values match
                assert(std::all_of(result.matrix.columns.back().begin(), result.matrix.columns.back().end(),
                        [&result, &b, &b_lookup_table, b_col_idx](auto x) {
                          if (result.row_values->at(x).dim() != b.column_values->at(b_col_idx).dim() + 1)
                              return false;
                          return result.row_values->at(x) == b.row_values->at(std::find(b_lookup_table.begin(), b_lookup_table.end(), x) - b_lookup_table.begin());
                        }));

                b_col_idx++;
            } else {
                // A and B columns contain (partial) coboundary of the same simplex
                assert(a.column_values->at(a_col_idx) == b.column_values->at(b_col_idx));
                // first, adjust indices in both columns
                Column& a_col = a.matrix.columns[a_col_idx];
                Column& b_col = b.matrix.columns[b_col_idx];

                apply_lookup(a_col, a_lookup_table);
                apply_lookup(b_col, b_lookup_table);

                // now merge them, duplicates are OK - just take one of them
                Column result_column = Column::merge_union(a_col, b_col);

                result.append_column(result_column, (*a.column_values)[a_col_idx]);

                // verify that dimensions match
                assert(std::all_of(result.matrix.columns.back().begin(), result.matrix.columns.back().end(),
                        [&result, &a, &b, a_col_idx, &a_lookup_table, &b_lookup_table](size_t x) {
                          if (result.row_values->at(x).dim() != a.column_values->at(a_col_idx).dim() + 1)
                              return false;
                          auto merged_value = result.row_values->at(x);

                          size_t aaa_idx = std::find(a_lookup_table.begin(), a_lookup_table.end(), x) - a_lookup_table.begin();
                          size_t bbb_idx = std::find(b_lookup_table.begin(), b_lookup_table.end(), x) - b_lookup_table.begin();


                          if (aaa_idx < a_lookup_table.size()) {
                              auto a_value = a.row_values->at(aaa_idx);
                              if (merged_value != a_value)
                                  return false;
                          }

                          if (bbb_idx < b_lookup_table.size()) {
                              auto b_value = b.row_values->at(bbb_idx);
                              if (merged_value != b_value)
                                  return false;
                          }
                          return (aaa_idx < a_lookup_table.size() or bbb_idx < b_lookup_table.size());

                        }));

                a_col_idx++;
                b_col_idx++;
            }
        }

        // process remaining columns in a
        while(a_col_idx < a.n_columns()) {
            result.append_column(a.matrix.columns[a_col_idx], (*a.column_values)[a_col_idx]);

            apply_lookup(result.matrix.columns.back(), a_lookup_table);

            // debug: verify that dimensions and values match
            assert(std::all_of(result.matrix.columns.back().begin(), result.matrix.columns.back().end(),
                    [&result, &a, &a_lookup_table, a_col_idx](auto x) {
                      if (result.row_values->at(x).dim() != a.column_values->at(a_col_idx).dim() + 1)
                          return false;
                      return result.row_values->at(x) == a.row_values->at(std::find(a_lookup_table.begin(), a_lookup_table.end(), x) - a_lookup_table.begin());
                    }));

            a_col_idx++;
        }

        // process remaining columns in b
        while(b_col_idx < b.n_columns()) {
            result.append_column(b.matrix.columns[b_col_idx], (*b.column_values)[b_col_idx]);
            apply_lookup(result.matrix.columns.back(), b_lookup_table);

            // debug: verify that dimensions and values match
            assert(std::all_of(result.matrix.columns.back().begin(), result.matrix.columns.back().end(),
                    [&result, &b, &b_lookup_table, b_col_idx](auto x) {
                      if (result.row_values->at(x).dim() != b.column_values->at(b_col_idx).dim() + 1)
                          return false;
                      return result.row_values->at(x) == b.row_values->at(std::find(b_lookup_table.begin(), b_lookup_table.end(), x) - b_lookup_table.begin());
                    }));

            b_col_idx++;
        }

        return result;
    }

    MatrixWithValues MatrixWithValues::merge_row_matrix(MatrixWithValues&, MatrixWithValues&)
    {
        throw std::runtime_error("not implemented");
    }

    std::string MatrixWithValues::debug_info() const
    {
        std::stringstream ss;
        ss << "MatrixWithValues{";
        ss << "column_values->size = " << column_values->size();
        ss << ", Mb = " << column_values->size() * sizeof(UidValue) / 1024 / 1024;
        if (not column_values->empty())
            ss << ", column_values = [" << column_values->front() << ", ..., " << column_values->back() << "]";
        ss << ", row_values->size = " << row_values->size();
        ss << ", Mb = " << row_values->size() * sizeof(UidValue) / 1024 / 1024;
        if (not row_values->empty())
            ss << ", row_values = [" << row_values->front() << ", ..., " << row_values->back() << "]";
        ss << ", " << matrix.debug_info();
        ss << "}";
        return ss.str();
    }

    size_t MatrixWithValues::size_in_bytes() const
    {
        size_t result;
        result += column_values->size() * sizeof(UidValue);
        result += row_values->size() * sizeof(UidValue);
        return result + matrix.size_in_bytes_total();
    }

    float MatrixWithValues::size_in_gigabytes() const
    {
        return size_in_bytes() / 1024 / 1024 / 1024;
    }

    std::ostream& operator<<(std::ostream& out, const MatrixWithValues& mv)
    {
        out << "MatrixWithValues(matrix = ";
        out << mv.matrix;
        out << ", column_values = " << mv.column_values;
        out << ", row_values = " << mv.row_values;
        out << ")";

        return out;
    }
} // namespace cadmus


namespace diy {

    void Serialization<cadmus::MatrixWithValues>::save(BinaryBuffer& bb, const cadmus::MatrixWithValues& m)
    {
        CALI_CXX_MARK_FUNCTION;
        diy::save(bb, m.matrix);
        diy::save(bb, *m.column_values);
        diy::save(bb, *m.row_values);
    }

    void Serialization<cadmus::MatrixWithValues>::load(BinaryBuffer& bb, cadmus::MatrixWithValues& m)
    {
        CALI_CXX_MARK_FUNCTION;
        diy::load(bb, m.matrix);
        cadmus::UidValues col_values;
        diy::load(bb, col_values);
        m.column_values = std::make_shared<cadmus::UidValues>(std::move(col_values));
        cadmus::UidValues row_values;
        diy::load(bb, row_values);
        m.row_values = std::make_shared<cadmus::UidValues>(std::move(row_values));
    }

} // namespace diy
