#pragma once

#include <vector>
#include <unordered_map>
#include <map>
#include <array>
#include <unordered_set>
#include <set>
#include <ostream>
#include <iostream>
#include <sstream>
#include <memory>
#include <string>

template<class T>
std::string container_to_string(const std::vector<T>& cont)
{
    std::stringstream ss;
    ss << "[";
    for(const auto& x : cont) {
        ss << x << ", ";
    }
    ss << "]";
    return ss.str();
}

template<class T>
std::string container_to_string(const std::unordered_set<T>& cont)
{
    std::stringstream ss;
    ss << "{";
    for(const auto& x : cont) {
        ss << x << ", ";
    }
    ss << "}";
    return ss.str();
}

template<class T, class C>
std::string container_to_string(const std::set<T, C>& cont)
{
    std::stringstream ss;
    ss << "{";
    for(const auto& x : cont) {
        ss << x << ", ";
    }
    ss << "}";
    return ss.str();
}

#ifdef CADMUS_USE_SPDLOG

#include "spdlog/spdlog.h"
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace spd=spdlog;

#else

namespace spd {

    namespace level {
        enum level_enum {
            trace,
            debug,
            info,
            warn,
            err,
            critical,
            off,
        };
    }

    struct logger {

        template<class... Args>
        static void set_pattern(Args... args) { };

        template<class... Args>
        static void critical(Args... args) { };

        template<class... Args>
        static void error(Args... args) { };

        template<class... Args>
        static void warn(Args... args) { };

        template<class... Args>
        static void info(Args... args) { };

        template<class... Args>
        static void debug(Args... args) { };

        template<class... Args>
        static void trace(Args... args) { };

        template<class T>
        static void set_level(T) {}

        template<class T>
        static void flush_on(T) {}

    };

    inline logger* stderr_color_mt(std::string level) { return nullptr; }
    inline logger* get(std::string level) { return nullptr; }


} // namespace spd
#endif

inline spd::level::level_enum log_level_from_string(std::string s)
{
    if (s == "trace")
        return spd::level::trace;
    else if (s == "debug")
        return spd::level::debug;
    else if (s == "info")
        return spd::level::info;
    else if (s == "warn")
        return spd::level::warn;
    else if (s == "err")
        return spd::level::err;
    else if (s == "critical")
        return spd::level::critical;
    else if (s == "off")
        return spd::level::off;
    else {
        std::cerr << "WARNING: unknown logging level " << s << ", setting to info" << std::endl;
        return spd::level::info;
    }
}


namespace std {
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
    {
        os << "[";

        for(size_t i = 0; i < v.size(); ++i) {
            os << v[i];
            if (i + 1 < v.size())
                os << ", ";
        }

        os << "]";
        return os;
    }

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& v)
    {
        os << "USet{";

        for(auto v_iter = v.cbegin(); v_iter != v.cend(); ++v_iter) {
            os << *v_iter;
            auto next_iter = v_iter;
            std::advance(next_iter, 1);
            if (next_iter != v.cend())
                os << ", ";
        }

        os << "}";
        return os;
    }

    template<typename Int, size_t DD>
    std::ostream& operator<<(std::ostream& os, const std::array<Int, DD>& p)
    {
        os << "(";

        for(size_t i = 0; i < DD; ++i) {
            os << p[i];
            if (i + 1 < DD)
                os << ", ";
        }

        os << ")";
        return os;
    }
}
