#pragma once

#include <utility>
#include <ostream>
#include <vector>
#include <limits>
#include <memory>

#include <diy/master.hpp>
#include <diy/types.hpp>
#include <diy/proxy.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi.hpp>
#include <diy/grid.hpp>
#include <diy/io/numpy.hpp>
#include <diy/point.hpp>

#include "hash_combine.h"

//#ifndef CADMUS_DIM
//#define CADMUS_DIM 3
//#endif

namespace cadmus {

    using Int = long long int;
    using dim_type = int;

    using Real = float;
    using IdxVector = std::vector<Int>;

    inline constexpr size_t k_invalid_id = std::numeric_limits<size_t>::max();
    inline constexpr Real k_plus_inf = std::numeric_limits<Real>::infinity();
    inline constexpr Real k_minus_inf = -std::numeric_limits<Real>::infinity();

    enum class ThrowInSanCheck {
        yes,
        no
    };

    struct UidValue {
        Int uid {-1};
        dim_type dim() const;

        Real value {std::numeric_limits<Real>::max()};

        UidValue() = default;
        UidValue(const UidValue&) = default;
        UidValue(UidValue&&) = default;
        UidValue& operator=(const UidValue&) = default;
        UidValue& operator=(UidValue&&) = default;

        UidValue(Int uid, Real value)
                :uid(uid), value(value) { }

        static UidValue smallest()
        {
            return UidValue(std::numeric_limits<dim_type>::max(), -std::numeric_limits<Real>::max());
        }

        inline bool operator==(const UidValue& other) const
        {
            return uid == other.uid and value == other.value;
        }

        inline bool operator!=(const UidValue& other) const
        {
            return !(*this == other);
        }

        inline bool operator<(const UidValue& other) const
        {
            if (dim() > other.dim())
                return true;

            if (dim() == other.dim() and value < other.value)
                return true;

            if (dim() == other.dim() and value == other.value and uid < other.uid)
                return true;

            return false;
        }

        inline bool operator<=(const UidValue& other) const
        {
            return (*this == other) or (*this < other);
        }

        inline bool operator>(const UidValue& other) const
        {
            return !(*this <= other);
        }

        inline bool operator>=(const UidValue& other) const
        {
            return !(*this < other);
        }

        friend std::ostream& operator<<(std::ostream& out, const UidValue& x);
    };

    struct SparseEntry {
        UidValue col_uid_value;
        UidValue row_uid_value;

        SparseEntry() = default;
        SparseEntry(const SparseEntry&) = default;
        SparseEntry(SparseEntry&&) = default;
        SparseEntry& operator=(const SparseEntry&) = default;
        SparseEntry& operator=(SparseEntry&&) = default;

        SparseEntry(const UidValue& col_uid_val, const UidValue& row_uid_val)
                :col_uid_value(col_uid_val), row_uid_value(row_uid_val) { }

        // for zero columns, row_uid_value is not initialized
        inline bool operator==(const SparseEntry& other) const
        {
            return col_uid_value == other.col_uid_value and row_uid_value == other.row_uid_value;
        }

        inline bool operator!=(const SparseEntry& other) const
        {
            return !(*this == other);
        }

        inline bool operator<(const SparseEntry& other) const
        {
            return (col_uid_value < other.col_uid_value) or (col_uid_value == other.col_uid_value and row_uid_value < other.row_uid_value);
        }

        inline bool operator<=(const SparseEntry& other) const
        {
            return (*this == other) or (*this < other);
        }

        inline bool operator>(const SparseEntry& other) const
        {
            return !(*this <= other);
        }

        inline bool operator>=(const SparseEntry& other) const
        {
            return !(*this < other);
        }

        friend std::ostream& operator<<(std::ostream& out, const SparseEntry& x);
    };

    using UidsSet = std::unordered_set<Int>;
    using UidValues = std::vector<UidValue>;
    using SPUidValues = std::shared_ptr<UidValues>;

    std::ostream& operator<<(std::ostream& out, const UidValues& x);

    struct PeristencePair {
        dim_type dim;
        Real birth;
        Real death;
        UidValue column_id;
        UidValue low_id;

        PeristencePair(dim_type dim, Real birth, Real death, UidValue cid, UidValue rid)
                :dim(dim), birth(birth), death(death), column_id(cid), low_id(rid) { }

    };

    using PersistencePairs = std::vector<PeristencePair>;

// diy types
    using Bounds = diy::DiscreteBounds;
    using Decomposer = diy::RegularDecomposer<Bounds>;
// NB: diy::Point is fixed-dim, subclass of std::array, but inside Bounds Point refers to DynamicPoint
    using DynamicPoint = Bounds::Point;
    using Point = diy::Point<DynamicPoint::Coordinate, CADMUS_DIM>;
    using Link = diy::RegularLink<Bounds>;
    using WrapVector = Decomposer::BoolVector;
    using Shape = diy::Point<Int, CADMUS_DIM>;
    using Grid = diy::Grid<Real, CADMUS_DIM>;
    using GridRef = diy::GridRef<Real, CADMUS_DIM>;
    using BlockID = diy::BlockID;

    using Vertex = Grid::Vertex;
    using VertexVec = std::vector<Vertex>;
    using VertexVecVec = std::vector<VertexVec>;

};

namespace diy {

    template<>
    struct Serialization<cadmus::UidValue> {
        static void save(BinaryBuffer& bb, const cadmus::UidValue&);
        static void load(BinaryBuffer& bb, cadmus::UidValue&);
    };

    template<>
    struct Serialization<::cadmus::SparseEntry> {
        static void save(BinaryBuffer& bb, const cadmus::SparseEntry&);
        static void load(BinaryBuffer& bb, cadmus::SparseEntry&);
    };

} // namespace diy

namespace std {

    template<>
    struct hash<cadmus::UidValue> {
        std::size_t operator()(const cadmus::UidValue& v) const
        {
            std::size_t seed = 0;
            cadmus::hash_combine(seed, std::hash<decltype(v.uid)>()(v.uid));
            cadmus::hash_combine(seed, std::hash<decltype(v.value)>()(v.value));
            return seed;
        }
    };

    template<>
    struct hash<cadmus::SparseEntry> {
        std::size_t operator()(const cadmus::SparseEntry& x) const
        {
            std::size_t seed = 0;
            cadmus::hash_combine(seed, std::hash<decltype(x.col_uid_value)>()(x.col_uid_value));
            cadmus::hash_combine(seed, std::hash<decltype(x.row_uid_value)>()(x.row_uid_value));
            //cadmus::hash_combine(seed, std::hash<decltype(x.is_empty)>()(x.is_empty));
            return seed;
        }
    };

}

