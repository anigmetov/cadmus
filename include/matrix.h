#pragma once

#include <vector>
#include <initializer_list>
#include <ostream>
#include <unordered_map>

#include "logger.h"
#include "types.h"

namespace cadmus {


struct Column {
    static constexpr size_t k_invalid_index = std::numeric_limits<size_t>::max();
    // types
    using Index = std::size_t;
    using Data = std::vector<Index>;

    // for now
    using value_type = Index;

    // constructor
    Column() = default;
    Column(const Column&) = default;
    Column(Column&&) = default;

    template<class Int>
    Column(std::initializer_list<Int> init)
    {
        for (Int value : init) {
            data.push_back(value);
        }
    }

    Column& operator=(const Column&) = default;
    Column& operator=(Column&&) = default;

    friend bool operator==(const Column& a, const Column& b) { return a.data == b.data; }
    friend bool operator!=(const Column& a, const Column& b) { return not(a== b); }

    explicit Column(Index i) : data({i}) {}

    // data
    Data data;

    // methods
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }

    auto begin() const { return data.cbegin(); }
    auto end() const { return data.cend(); }

    auto cbegin() const { return data.cbegin(); }
    auto cend() const { return data.cend(); }

    auto rbegin() { return data.rbegin(); }
    auto rend() { return data.rend(); }

    auto rbegin() const { return data.rbegin(); }
    auto rend() const { return data.rend(); }

    void push_back(Index v) { data.push_back(v); }
    void reserve(size_t n) { data.reserve(n); }
    size_t size() const { return data.size(); }

    void resize(size_t n) { data.resize(n); }

    auto& front()       { return data.front(); }
    auto  front() const { return data.front(); }
    auto& back()        { return data.back(); }
    auto  back() const  { return data.back(); }

    bool is_zero() const { return data.size() == 0; }
    void clear() { data.clear(); }
    void shrink_to_fit() { data.shrink_to_fit(); }
    void clear_and_shrink() { clear(); shrink_to_fit(); }

    size_t low() const { return is_zero() ? k_invalid_index : data.back(); }

    void drop_low() { assert(not is_zero()); data.pop_back(); }

    void ultrasparsify()
    {
        assert(not is_zero());
        data[0] = data.back();
        data.resize(1);
        data.shrink_to_fit();
    }

    void sanity_check();

    size_t size_in_bytes() const { return data.size() * sizeof(Data::value_type); }

    inline friend void swap(Column& a, Column& b)
    {
        std::swap(a.data, b.data);
    }

    static Column merge_union(const Column& a, const Column& b);

    friend std::ostream& operator<<(std::ostream& out, const Column& a);
};


void add(const Column& reduced, const Column& pivot, Column& result, int alpha);

struct Matrix {
    // both column and row are stored as vectors,
    // just an alias for readability
    using Row = Column;
    using Columns = std::vector<Column>;
    using Rows = std::vector<Row>;

    // data
    Columns columns;
    Rows rows;
    bool stores_rows{false};
    bool stores_columns{false};

    size_t n_rows {0};
    size_t n_columns{0};

    // methods
    Matrix() = default;
    Matrix(const Matrix&) = default;
    Matrix(Matrix&&) noexcept = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&&) = default;

    Matrix(const Rows& rows, const Columns& columns, int n_rows, int n_columns);

    static Matrix create_col_matrix(const Columns& columns, size_t n_rows);

    // will move from col
    void append_column(Column& col);


    friend bool operator==(const Matrix& a, const Matrix& b);
    friend bool operator!=(const Matrix& a, const Matrix& b);


    void add_to_column(size_t reduced_idx, size_t pivot_idx, Int alpha);
    void add_to_row(size_t reduced_idx, size_t pivot_idx, Int alpha);

    void convert_to_rows(bool keep_columns);
    void convert_to_columns(bool keep_rows);
    void clear_columns();
    void clear_rows();

    void swap_columns(size_t i, size_t j);

    void clear();

    Matrix antitranspose(bool fill_columns, bool fill_rows);


    size_t size_in_bytes_columns() const;
    size_t size_in_bytes_rows() const;
    size_t size_in_bytes_total() const;
    size_t max_column_length() const;

    template<class P>
    std::pair<Matrix, Matrix> split_columns(P pred) const
    {
        if (not stores_columns or stores_rows) {
            throw std::runtime_error("not implemented");
        }

        Matrix first, second;

        first.stores_columns = second.stores_columns = true;
        first.stores_rows = second.stores_rows = false;

        std::unordered_map<size_t, size_t> first_lookup, second_lookup;

        size_t new_first = 0, new_second = 0;

        for(size_t i = 0; i < n_rows; ++i) {
            if (pred(i))
                first_lookup[i] = new_first++;
            else
                second_lookup[i] = new_second++;
        }

        first.n_rows = first_lookup.size();
        second.n_rows = second_lookup.size();

        for(size_t col_idx = 0; col_idx < columns.size(); ++col_idx) {
            if (pred(col_idx)) {
                first.columns.push_back(columns[col_idx]);
                for(auto& x : first.columns.back())
                    x = first_lookup.at(x);
            } else {
                second.columns.push_back(columns[col_idx]);
                for(auto& x : second.columns.back())
                    x = second_lookup.at(x);
            }
        }

        first.n_columns = first.columns.size();
        second.n_columns = second.columns.size();

        return {first, second};
    }


    template<class P>
    std::tuple<Matrix, Matrix, std::unordered_map<size_t, size_t>, std::unordered_map<size_t, size_t>> split_rows(P pred) const
    {
        if (not stores_columns or stores_rows) {
            throw std::runtime_error("not implemented");
        }

        Matrix first, second;

        first.stores_columns = second.stores_columns = true;
        first.stores_rows = second.stores_rows = false;

        // initialize with empty columns
        first.columns = second.columns = Columns(n_columns);

        std::unordered_map<size_t, size_t> first_lookup, second_lookup;

        size_t new_first = 0, new_second = 0;

        for(size_t i = 0; i < n_rows; ++i) {
            if (pred(i))
                first_lookup[i] = new_first++;
            else
                second_lookup[i] = new_second++;
        }

        first.n_rows = first_lookup.size();
        second.n_rows = second_lookup.size();

        for(size_t col_idx = 0; col_idx < columns.size(); ++col_idx) {
            for(auto entry : columns[col_idx]) {
                if (pred(entry)) {
                    first.columns[col_idx].push_back(first_lookup.at(entry));
                } else {
                    second.columns[col_idx].push_back(second_lookup.at(entry));
                }
            }
        }

        first.n_columns = first.columns.size();
        second.n_columns = second.columns.size();

        return {first, second, first_lookup, second_lookup};
    }

    void sanity_check();

    std::string debug_info() const;

    friend std::ostream& operator<<(std::ostream& out, const Matrix& a);

private:
    // to avoid many memory allocations, always put addition results in this column
    // and then swap
    Column addition_result_cache_;
};

} // namespace cadmus


namespace diy {

template<>
struct Serialization<cadmus::Column>
{
    static void save(BinaryBuffer& bb, const cadmus::Column& c);
    static void load(BinaryBuffer& bb, cadmus::Column& c);
};

template<>
struct Serialization<cadmus::Matrix>
{
    static void save(BinaryBuffer& bb, const cadmus::Matrix& m);
    static void load(BinaryBuffer& bb, cadmus::Matrix& m);
};
} // namespace diy