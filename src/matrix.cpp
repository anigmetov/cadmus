#include "profile.h"
#include "matrix.h"

namespace cadmus {

void Column::sanity_check()
{
    for(size_t i = 0 ; i < size() ; ++i) {
        if (i > 0 and data[i] <= data[i - 1])
            throw std::runtime_error("unsorted column");
    }
}

Column Column::merge_union(const Column& a, const Column& b)
{
    Column result;
    result.data.reserve(a.size() + b.size());
    std::set_union(a.cbegin(), a.cend(), b.cbegin(), b.cend(), std::back_inserter(result.data));
    return result;
}

std::ostream& operator<<(std::ostream& out, const Column& a)
{
    out << "Column(size=" << a.size() << " ";
    out << a.data;
    out << ")";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
    out << "Matrix(n_rows=" << m.n_rows << ", n_columns = " << m.n_columns;
    out << ", rows.size = " << m.rows.size() << ", columns.size = " << m.columns.size();
    out << ", stores_rows = " << m.stores_rows << ", stores_columns = " << m.stores_columns;
    out << ")";
    return out;
}

void add(const Column& reduced, const Column& pivot, Column& result, int alpha)
{
    auto reduced_iter = reduced.cbegin();
    auto pivot_iter = pivot.cbegin();

    while(reduced_iter != reduced.cend() and pivot_iter != pivot.cend()) {
        if (*reduced_iter < *pivot_iter) {
            result.push_back(*reduced_iter);
            reduced_iter++;
        } else if (*reduced_iter > *pivot_iter) {
            result.push_back(*pivot_iter);
            pivot_iter++;
        } else {
            reduced_iter++;
            pivot_iter++;
        }
    }

    while(reduced_iter != reduced.cend()) {
        result.push_back(*reduced_iter);
        reduced_iter++;
    }

    while(pivot_iter != pivot.cend()) {
        result.push_back(*pivot_iter);
        pivot_iter++;
    }
}

Matrix::Matrix(const Rows& rows, const Columns& columns, int n_rows, int n_columns)
        :columns(columns), rows(rows), n_columns(n_columns), n_rows(n_rows)
{
    if (n_columns == 0 or n_rows == 0)
        throw std::runtime_error("currently empty matrix is not supported");

    if (rows.empty() and columns.empty())
        throw std::runtime_error("currently empty matrix is not supported");

    stores_rows = not rows.empty();
    stores_columns = not columns.empty();
}

Matrix Matrix::create_col_matrix(const Columns& columns, size_t n_rows)
{
    Matrix result;

    result.stores_rows = false;
    result.stores_columns = true;
    result.columns = columns;
    result.n_rows = n_rows;

    return result;
}

void Matrix::append_column(Column& col)
{
    n_columns++;
    columns.emplace_back(std::move(col));
}

void Matrix::add_to_column(size_t reduced_idx, size_t pivot_idx, Int alpha)
{
    Column& reduced = columns[reduced_idx];
    Column& pivot = columns[pivot_idx];

    addition_result_cache_.clear();
    addition_result_cache_.reserve(std::max(reduced.size(), pivot.size()));

    add(reduced, pivot, addition_result_cache_, alpha);

    reduced = std::move(addition_result_cache_);
}

void Matrix::swap_columns(size_t i, size_t j)
{
    std::swap(columns[i], columns[j]);
}


void Matrix::add_to_row(size_t reduced_idx, size_t pivot_idx, Int alpha)
{
    auto logger = spd::get("console");
    if (logger) logger->trace("add_to_row: reduced_idx = {}, pivot_idx = {}, n_rows = {}, rows.size = {}", reduced_idx, pivot_idx, n_rows, rows.size());

    Row& reduced = rows[reduced_idx];
    Row& pivot = rows[pivot_idx];

    if (logger) logger->debug("add_to_row: reduced_idx = {}, pivot_idx = {}, reduced.size = {}, pivot.size = {}", reduced_idx, pivot_idx, reduced.size(), pivot.size());

    addition_result_cache_.reserve(std::max(reduced.size(), pivot.size()));

    add(reduced, pivot, addition_result_cache_, alpha);

    if (logger) logger->debug("add_to_row: reduced_idx = {}, pivot_idx = {}, reduced.size = {}, pivot.size = {}, result.size = {}", reduced_idx, pivot_idx, reduced.size(), pivot.size(), addition_result_cache_.size());

    reduced = std::move(addition_result_cache_);
}

Matrix Matrix::antitranspose(bool fill_columns, bool fill_rows)
{
    CALI_CXX_MARK_FUNCTION;

    auto logger = spd::get("console");

    if (not(stores_columns and not stores_rows))
        throw std::runtime_error("not implmented, should not happen for D");

    Matrix result;
    result.n_rows = n_columns;
    result.n_columns = n_rows;

    result.stores_columns = fill_columns;
    result.stores_rows = fill_rows;

    if (result.stores_columns)
        result.columns = std::vector<Column>(result.n_columns);

    if (result.stores_rows)
        result.rows = std::vector<Column>(result.n_rows);

    if (stores_columns) {
        for(size_t col_idx = 0 ; col_idx < n_columns ; ++col_idx) {
            auto& col = columns.at(n_columns - 1 - col_idx).data;

            for(auto rev_iter = col.rbegin() ; rev_iter != col.rend() ; ++rev_iter) {
                size_t row_idx = n_rows - 1 - *rev_iter;

                if (logger) logger->trace("col_idx = {}", col_idx);

                if (result.stores_columns)
                    result.columns.at(row_idx).push_back(col_idx);

                if (result.stores_rows)
                    result.rows.at(col_idx).push_back(row_idx);
            }
            col.clear();
            col.shrink_to_fit();
        }
    } else if (stores_rows) {
         for(size_t row_idx = 0 ; row_idx < n_columns ; ++row_idx) {

            const auto& row = rows.at(n_rows - 1 - row_idx);

            for(auto c = row.rbegin() ; c != row.rend() ; ++c) {
                size_t col_idx = n_columns - 1 - *c;

                if (logger) logger->trace("row_idx = {}", row_idx);

                if (result.stores_columns)
                    result.columns.at(col_idx).push_back(row_idx);

                if (result.stores_rows)
                    result.rows.at(row_idx).push_back(col_idx);
            }
        }
    }

    if (logger) logger->trace("antitranspose OK");

    return result;
}

void Matrix::sanity_check()
{
    if (not stores_rows and not rows.empty())
        throw std::runtime_error("rows must be empty");

    if (stores_rows and n_rows != rows.size())
        throw std::runtime_error("bad n_rows");

    for(auto& row: rows)
        row.sanity_check();

    if (not stores_columns and not columns.empty())
        throw std::runtime_error("columns must be empty");

    if (stores_columns and n_columns != columns.size())
        throw std::runtime_error("bad n_columns");

    for(auto& column: columns)
        column.sanity_check();

    if (stores_rows and stores_columns) {
        std::set<std::pair<size_t, size_t>> row_entries, column_entries;

        for(size_t i = 0 ; i < n_rows ; ++i)
            for(size_t j: rows[i])
                row_entries.emplace(i, j);

        for(size_t j = 0 ; j < n_columns ; ++j)
            for(size_t i: columns[j])
                column_entries.emplace(i, j);

        if (row_entries != column_entries)
            throw std::runtime_error("inconsistency between rows and columns");
    }
}

void Matrix::clear()
{
    clear_columns();
    clear_rows();
}

void Matrix::clear_columns()
{
    stores_columns = false;
    columns.clear();
}

void Matrix::clear_rows()
{
    stores_rows = false;
    rows.clear();
}

void Matrix::convert_to_rows(bool keep_columns)
{
    if (not stores_columns)
        throw std::runtime_error("cannot convert to rows: columns empty");

    if (stores_rows or not rows.empty())
        throw std::runtime_error("refuset to convert to rows: rows already filled");

    stores_rows = true;

    rows = std::vector<Row>(n_rows);

    for(size_t col_idx = 0 ; col_idx < columns.size() ; ++col_idx) {
        for(auto row_idx: columns[col_idx]) {
            rows.at(row_idx).push_back(col_idx);
        }
    }

    if (not keep_columns) {
        clear_columns();
    }
}

void Matrix::convert_to_columns(bool keep_rows)
{
    if (not stores_rows)
        throw std::runtime_error("cannot convert to columns: rows empty");

    if (stores_columns or not columns.empty())
        throw std::runtime_error("refuset to convert to rows: rows already filled");

    stores_columns = true;

    columns = std::vector<Row>(n_columns);

    for(size_t row_idx = 0 ; row_idx < rows.size() ; ++row_idx) {
        for(auto col_idx: rows[row_idx]) {
            columns.at(col_idx).push_back(row_idx);
        }
    }

    if (not keep_rows) {
        clear_rows();
    }
}

std::string Matrix::debug_info() const
{
    std::stringstream ss;
    ss << "Matrix(stores_columns = " << stores_columns;
    ss << ", stores_rows = " << stores_rows;
    ss << ", n_columns = " << n_columns;
    ss << ", n_rows = " << n_rows;
    ss << ", columns.size = " << columns.size();
    ss << ", rows.size = " << rows.size();
    ss << ")";
    return ss.str();
}

size_t Matrix::size_in_bytes_columns() const
{
    size_t result = 0;
    if (stores_columns) {
        result += columns.size() * sizeof(Column);
        result += std::accumulate(columns.cbegin(), columns.cend(), 0LL, [](size_t x, const Column& col) { return x + col.size_in_bytes(); });
    }
    return result;
}

size_t Matrix::max_column_length() const
{
    if (columns.empty())
        return 0;
    else
        return std::max_element(columns.cbegin(), columns.cend(), [](const Column& a, const Column& b) { return a.size() < b.size(); })->size();
}

size_t Matrix::size_in_bytes_rows() const
{
    size_t result = 0;
    if (stores_rows) {
        result += rows.size() * sizeof(Column);
        result += std::accumulate(rows.cbegin(), rows.cend(), 0LL, [](size_t x, const Column& col) { return x + col.size_in_bytes(); });
    }
    return result;
}

size_t Matrix::size_in_bytes_total() const
{
    return size_in_bytes_columns() + size_in_bytes_rows();
}

bool operator==(const Matrix& a, const Matrix& b)
{
    if (a.stores_rows != b.stores_rows)
        return false;

    if (a.stores_columns != b.stores_columns)
        return false;

    if (a.n_rows != b.n_rows)
        return false;

    if (a.n_columns != b.n_columns)
        return false;

    if (a.rows != b.rows)
        return false;

    if (a.columns != b.columns)
        return false;

    return true;
}

bool operator!=(const Matrix& a, const Matrix& b)
{
    return !(a == b);
}

} // namespace cadmus

namespace diy {

void Serialization<cadmus::Column>::save(BinaryBuffer& bb, const cadmus::Column& c)
{
    diy::save(bb, c.data);
}

void Serialization<cadmus::Column>::load(BinaryBuffer& bb, cadmus::Column& c)
{
    diy::load(bb, c.data);
}

void Serialization<cadmus::Matrix>::save(BinaryBuffer& bb, const cadmus::Matrix& m)
{
    CALI_CXX_MARK_FUNCTION;
    diy::save(bb, m.n_rows);
    diy::save(bb, m.n_columns);
    diy::save(bb, m.stores_rows);
    diy::save(bb, m.stores_columns);
    diy::save(bb, m.rows);
    diy::save(bb, m.columns);
}

void Serialization<cadmus::Matrix>::load(BinaryBuffer& bb, cadmus::Matrix& m)
{
    CALI_CXX_MARK_FUNCTION;
    diy::load(bb, m.n_rows);
    diy::load(bb, m.n_columns);
    diy::load(bb, m.stores_rows);
    diy::load(bb, m.stores_columns);
    diy::load(bb, m.rows);
    diy::load(bb, m.columns);
}

} // namespace diy
