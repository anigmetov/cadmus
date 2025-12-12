#pragma once

#include <cassert>

#include <diy/serialization.hpp>

#include "matrix.h"

namespace cadmus {


using UltraSparseMatrixWithValues = std::vector<SparseEntry>;

UltraSparseMatrixWithValues merge(const UltraSparseMatrixWithValues& a, const UltraSparseMatrixWithValues& b);

struct MatrixWithValues {
    using Value = UidValue;
    using Values  = UidValues;
    using SPValues = SPUidValues;

    Matrix matrix;

    SPValues column_values {nullptr};
    SPValues row_values {nullptr};

    MatrixWithValues() = default;
    MatrixWithValues(const MatrixWithValues&) = default;
    MatrixWithValues(MatrixWithValues&&) = default;

    MatrixWithValues& operator=(const MatrixWithValues&) = default;
    MatrixWithValues& operator=(MatrixWithValues&&) = default;

    MatrixWithValues(Matrix&& matrix, Values&& column_values, Values&& row_values) : matrix(std::move(matrix)),
                                                                                     column_values(std::make_shared<UidValues>(std::move(column_values))),
                                                                                     row_values(std::make_shared<UidValues>(std::move(row_values))) {}

    MatrixWithValues(Matrix&& matrix, SPValues column_values, SPValues row_values) : matrix(std::move(matrix)),
                                                                                     column_values(column_values),
                                                                                     row_values(row_values) {}


    bool stores_rows() const { return matrix.stores_rows; }
    bool stores_columns() const { return matrix.stores_columns; }

    void convert_to_rows(bool keep_columns) { matrix.convert_to_rows(keep_columns); }
    void convert_to_columns(bool keep_rows) { matrix.convert_to_columns(keep_rows); }

    [[nodiscard]] size_t n_columns() const { return matrix.n_columns; }
    [[nodiscard]] size_t n_rows() const { return matrix.n_rows; }

    void append_column(Column& col, Value val);

    static MatrixWithValues merge(MatrixWithValues&, MatrixWithValues&);
    static MatrixWithValues merge_row_matrix(MatrixWithValues&, MatrixWithValues&);
    static MatrixWithValues merge_column_matrix(MatrixWithValues&, MatrixWithValues&);

    std::string debug_info() const;

    size_t size_in_bytes() const;
    float size_in_gigabytes() const;
    inline size_t max_column_length() const { return matrix.max_column_length(); }


    bool sanity_check(bool is_cohomology);

    friend std::ostream& operator<<(std::ostream& out, const MatrixWithValues& m);
};

}

namespace diy {

template<>
struct Serialization<cadmus::MatrixWithValues>
{
    static void save(BinaryBuffer& bb, const cadmus::MatrixWithValues& m);
    static void load(BinaryBuffer& bb, cadmus::MatrixWithValues& m);
};
}