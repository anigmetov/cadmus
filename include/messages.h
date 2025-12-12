#pragma once

#include <vector>
#include <utility>

#include <diy/serialization.hpp>

#include "types.h"
#include "profile.h"
#include "values_only_matrix.h"

namespace cadmus {

// helper struct to redistribute matrix across ranks
    struct RearrangeMessage {
        bool has_uu_matrix {false};
        bool has_ui_matrix {false};
        bool has_us_matrix {false};
        ValueMatrix ui_matrix {};
        ValueMatrix uu_matrix {};
        UltraSparseMatrixWithValues us_matrix {};
        RearrangeMessage() = default;
        RearrangeMessage(const RearrangeMessage&) = default;
        RearrangeMessage(RearrangeMessage&&) = default;
        RearrangeMessage& operator=(const RearrangeMessage&) = default;
        RearrangeMessage& operator=(RearrangeMessage&&) = default;
    };

    struct MainLoopMessageBit {
        int reducer_gid;
        UidValue low;
        MainLoopMessageBit(int reducer_gid, UidValue low) : reducer_gid(reducer_gid), low(low) {}
        MainLoopMessageBit() = default;
        MainLoopMessageBit(const MainLoopMessageBit&) = default;
        MainLoopMessageBit(MainLoopMessageBit&&) = default;
        MainLoopMessageBit& operator=(const MainLoopMessageBit&) = default;
        MainLoopMessageBit& operator=(MainLoopMessageBit&&) = default;
    };

    struct MainLoopMessage {
        std::vector<MainLoopMessageBit> message;

        size_t size() const { return message.size(); }

        void emplace_back(int reducer_gid, UidValue low)
        {
            message.emplace_back(reducer_gid, low);
        };

    };

// helper struct for distributed reduction
    struct LowInfo {
        bool is_ultrasparse;
        UidValue col_value;
        UidValue low;
    };

    struct PivotMessage {
        std::vector<UidValue> lows;

        void push_back(UidValue);

        PivotMessage() = default;
        PivotMessage(const PivotMessage&) = default;
        PivotMessage(PivotMessage&&) = default;
        PivotMessage& operator=(const PivotMessage&) = default;
        PivotMessage& operator=(PivotMessage&&) = default;
    };

    using SendTo = std::pair<int, UidValue>;
    using SendToMessage = std::vector<SendTo>;

} //namespace cadmus

namespace diy {

    template<>
    struct Serialization<cadmus::RearrangeMessage> {
        static void save(BinaryBuffer& bb, const cadmus::RearrangeMessage& msg)
        {
            auto logger = spd::get("console");

            if (logger) logger->trace("RearrangeMessage save, start");
            diy::save(bb, msg.has_ui_matrix);
            diy::save(bb, msg.has_uu_matrix);
            diy::save(bb, msg.has_us_matrix);

            if (logger) logger->trace("RearrangeMessage save, flags ok");

            if (msg.has_ui_matrix)
                diy::save(bb, msg.ui_matrix);

            if (logger) logger->trace("RearrangeMessage save, ui ok");

            if (msg.has_uu_matrix)
                diy::save(bb, msg.uu_matrix);

            if (logger) logger->trace("RearrangeMessage save, uu ok");

            if (msg.has_us_matrix)
                diy::save(bb, msg.us_matrix);

            if (logger) logger->trace("RearrangeMessage save, us ok");
        }

        static void load(BinaryBuffer& bb, cadmus::RearrangeMessage& msg)
        {
            diy::load(bb, msg.has_ui_matrix);
            diy::load(bb, msg.has_uu_matrix);
            diy::load(bb, msg.has_us_matrix);

            if (msg.has_ui_matrix)
                diy::load(bb, msg.ui_matrix);

            if (msg.has_uu_matrix)
                diy::load(bb, msg.uu_matrix);

            if (msg.has_us_matrix)
                diy::load(bb, msg.us_matrix);
        }

    };

    template<>
    struct Serialization<cadmus::MainLoopMessageBit> {
        static void save(BinaryBuffer& bb, const cadmus::MainLoopMessageBit& msg)
        {
            diy::save(bb, msg.reducer_gid);
            diy::save(bb, msg.low);
        }

        static void load(BinaryBuffer& bb, cadmus::MainLoopMessageBit& msg)
        {
            diy::load(bb, msg.reducer_gid);
            diy::load(bb, msg.low);
        }
    };

    template<>
    struct Serialization<cadmus::MainLoopMessage> {
        static void save(BinaryBuffer& bb, const cadmus::MainLoopMessage& msg)
        {
            diy::save(bb, msg.message);
        }

        static void load(BinaryBuffer& bb, cadmus::MainLoopMessage& msg)
        {
            diy::load(bb, msg.message);
        }
    };

    template<>
    struct Serialization<cadmus::PivotMessage> {
        static void save(BinaryBuffer& bb, const cadmus::PivotMessage& msg)
        {
            CALI_CXX_MARK_FUNCTION;
            diy::save(bb, msg.lows);
        }

        static void load(BinaryBuffer& bb, cadmus::PivotMessage& msg)
        {
            CALI_CXX_MARK_FUNCTION;
            diy::load(bb, msg.lows);
        }
    };
} // namespace diy


