#include "profile.h"
#include "messages.h"

namespace cadmus {
    void PivotMessage::push_back(UidValue row_value)
    {
        lows.emplace_back(row_value);
    }
}

