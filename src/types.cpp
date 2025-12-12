#include "types.h"

namespace cadmus {

dim_type UidValue::dim() const
{
    dim_type res = 0;
    for(int d = 0; d < CADMUS_DIM; ++d) {
        if (uid & (1 << d))
            res++;
    }
    return res;
}

std::ostream& operator<<(std::ostream& out, const UidValue& x)
{
    out << "UidValue(value=" << x.value << ", dim=" << x.dim() << ", uid=" << x.uid << ")";
    return out;
}

std::ostream& operator<<(std::ostream& out, const SparseEntry& x)
{
    out << "SparseEntry(col = " << x.col_uid_value;
    out << ", row = " << x.row_uid_value << ")";
    return out;
}

std::ostream& operator<<(std::ostream& out, const UidValues& x)
{
    out << "[";
    for(size_t i = 0 ; i < x.size() ; ++i) {
        out << x[i];
        if (i + 1 < x.size())
            out << ", ";
    }
    out << "]";
    return out;
}

} // namespace cadmus

namespace diy {

void Serialization<cadmus::UidValue>::save(BinaryBuffer& bb, const cadmus::UidValue& x)
{
    diy::save(bb, x.value);
    diy::save(bb, x.uid);
}

void Serialization<cadmus::UidValue>::load(BinaryBuffer& bb, cadmus::UidValue& x)
{
    diy::load(bb, x.value);
    diy::load(bb, x.uid);
}

void Serialization<cadmus::SparseEntry>::save(BinaryBuffer& bb, const cadmus::SparseEntry& x)
{
    diy::save(bb, x.col_uid_value);
    diy::save(bb, x.row_uid_value);
}

void Serialization<cadmus::SparseEntry>::load(BinaryBuffer& bb, cadmus::SparseEntry& x)
{
    diy::load(bb, x.col_uid_value);
    diy::load(bb, x.row_uid_value);
}

}

