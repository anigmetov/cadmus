#pragma once

#include <vector>
#include <map>
#include <cassert>

#include "types.h"
#include "matrix_with_values.h"

namespace cadmus {

    class ValueColumn {
    public:
        // data must be sorted
        UidValues data_;
    public:
        ValueColumn() = default;
        ValueColumn(const ValueColumn&) = default;
        ValueColumn(ValueColumn&&) = default;
        ValueColumn& operator=(const ValueColumn&) = delete;
        ValueColumn& operator=(ValueColumn&&) = default;
        ValueColumn(const UidValues& data) : data_(data) { assert(sanity_check(ThrowInSanCheck::no)); }
        ValueColumn(UidValues&& data) : data_(std::move(data)) { assert(sanity_check(ThrowInSanCheck::no)); }

        void push_back(const UidValue& x) { data_.push_back(x); }

        auto begin() { return data_.begin(); }
        auto end() { return data_.end(); }
        auto begin() const { return data_.cbegin(); }
        auto end() const { return data_.cend(); }

        bool is_zero() const { return data_.empty(); }
        // will throw if data_ is zero
        UidValue low() const { return data_.back(); }

        UidValue&       operator[](size_t idx)       { return data_[idx]; }
        const UidValue& operator[](size_t idx) const { return data_[idx]; }

        void drop_low() { data_.pop_back(); }

        size_t size() const { return data_.size(); }

        // merge data from other into this
        static ValueColumn merge(ValueColumn& a, ValueColumn& b);

        bool sanity_check(ThrowInSanCheck throw_in_san_check);

        friend std::ostream& operator<<(std::ostream&, const ValueColumn&);

        friend diy::Serialization<ValueColumn>;

        static std::set<UidValue, std::greater<>> addition_cache_;

        void load_to_cache() const;
        void load_from_cache();
        static UidValue cache_low();
        static bool is_cache_zero();

        static void add_to_cache(const ValueColumn& pivot_col);
        void sparsify_reduced();
    };

    class ValueMatrix;

    class ValueMatrices;

    struct ValueMatrix {
    public:
        using DataType = std::map<UidValue, ValueColumn, std::greater<UidValue>>;
        using PivotsType = std::unordered_map<UidValue, UidValue>;

        ValueMatrix() = default;
        ValueMatrix(const ValueMatrix&) = default;
        ValueMatrix(ValueMatrix&&) = default;
        ValueMatrix& operator=(const ValueMatrix&) = default;
        ValueMatrix& operator=(ValueMatrix&&) = default;
        ValueMatrix(MatrixWithValues& mv);
        bool sanity_check(ThrowInSanCheck throw_in_san_check);
        static ValueMatrices merge_and_split(ValueMatrix& a, ValueMatrix& b);
        static ValueMatrix merge(ValueMatrix& a, ValueMatrix& b);

        auto begin()       { return columns_.begin(); }
        auto begin() const { return columns_.cbegin(); }
        auto end()         { return columns_.end(); }
        auto end() const   { return columns_.cend(); }

        auto find(const UidValue& x) { return columns_.find(x); }
        auto find(const UidValue& x) const { return columns_.find(x); }
        auto count(const UidValue& x) const { return columns_.count(x); }

        auto size() const { return columns_.size(); }

        DataType::iterator  erase(DataType::iterator iter);
        DataType::iterator  erase(DataType::const_iterator iter);
        DataType::size_type erase(DataType::key_type col_val, bool not_found_ok = false);
        DataType::size_type erase(DataType::key_type col_val, UidValue low_val);

        ValueColumn& operator[](const UidValue& x)             { return columns_[x]; }

        ValueColumn& at(const UidValue& x)             { return columns_.at(x); }
        const ValueColumn& at(const UidValue& x) const { return columns_.at(x); }

        void insert(DataType::value_type& x)       { columns_.insert(std::move(x)); }
        void insert(const DataType::value_type& x) { columns_.insert(x); }

        void reduce_local(Int skip_uid, Stats& stats, PivotsType* prev_dim_pivots, bool clearing, dim_type dim);

        UidValue column_by_low(const UidValue& low);

        ValueMatrix extract_submatrix(const UidValues& col_vals);

        void append_mv(ValueMatrix&& m);

        UidValues reduce_starting_from(UidValue start_val, Int skip_uid);

        // appends to result
        void read_off_persistence(bool negate, Int global_min_vertex_uid, PersistencePairs& result) const;

        size_t total_column_size_in_bytes() const;
        size_t total_column_size_in_megabytes() const;

        UidValue first_value() const { return columns_.begin()->first; }

        bool is_reduced(Int skip_uid) const;
        bool is_local_pivots_correct(Int skip_uid) const;

        int gid_ {-1};

        friend diy::Serialization<ValueMatrix>;
        // members
        // columns are ordered in reverse order
        DataType columns_;

        PivotsType local_pivots_;
        // for debug
        std::unordered_set<UidValue> seen_cols_;
        // methods

        void set_local_pivot(const UidValue& low_val, const UidValue& col_val);
        unsigned long erase_column(const UidValue& col_val);
        void sparsify_reduced();
    };

    void append_sparse(ValueMatrices& matrices, const UltraSparseMatrixWithValues& sparse_);

    struct ValueMatrices {
        std::array<ValueMatrix, CADMUS_DIM> matrices_;

        ValueMatrix&       operator[](size_t idx)       { return matrices_[idx]; }
        const ValueMatrix& operator[](size_t idx) const { return matrices_[idx]; }

        size_t size() const;
        size_t total_column_size_in_megabytes() const;

        void append_mv(ValueMatrices&& m);
    };


} // namespace cadmus


namespace diy {

template<>
struct Serialization<cadmus::ValueColumn>
{
    static void save(BinaryBuffer& bb, const cadmus::ValueColumn&);
    static void load(BinaryBuffer& bb, cadmus::ValueColumn&);
};

template<>
struct Serialization<cadmus::ValueMatrix>
{
    static void save(BinaryBuffer& bb, const cadmus::ValueMatrix&);
    static void load(BinaryBuffer& bb, cadmus::ValueMatrix&);
};

template<>
struct Serialization<cadmus::ValueMatrices>
{
    static void save(BinaryBuffer& bb, const cadmus::ValueMatrices&);
    static void load(BinaryBuffer& bb, cadmus::ValueMatrices&);
};

} // namespace diy
