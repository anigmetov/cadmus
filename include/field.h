#pragma once

#include <vector>


namespace cadmus {

template<int P>
class Field {
};

template<>
class Field<2> {
    using Entry = size_t;

    static Entry one(size_t pos) { return pos; }
    static Entry minus_one(size_t pos) { return pos; }

};


}