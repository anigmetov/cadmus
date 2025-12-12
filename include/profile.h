#pragma once

#ifdef CADMUS_USE_CALIPER
#include <caliper/cali.h>
#else
#define CALI_CXX_MARK_FUNCTION
#define CALI_MARK_BEGIN(x) x;
#define CALI_MARK_END(x) x;
#endif


#include <sys/resource.h>
#include <sstream>
#include <ostream>
#include <iostream>
#include <chrono>

#include <spdlog/fmt/fmt.h>

#include "types.h"
#include "logger.h"

namespace cadmus {
    // all timings are in seconds
    struct Stats {
        int gid;
        int n_blocks {0};
        long int n_simplices {0};
        long int n_columns_after_rearrange {0};
        long int final_n_columns {0};
        size_t final_columns_in_mb {0};
        int interior_cleared {0};
        int local_cleared {0};
        int global_cleared {0};
        bool clearing_opt;
        bool swap_reduction;
        double local_init {0};
        double interior_reduction {0};

        double sparsify {0};
        double exchange {0};
        double rearrange_send {0};
        double rearrange_recv {0};
        double final_sort {0};
        double total_wo_reduction {0};

        std::array<double, CADMUS_DIM> local_reduction_per_dim {};
        std::array<double, CADMUS_DIM> final_reduction_per_dim {};
        std::array<double, CADMUS_DIM> exchange_pivots_send_per_dim {};
        std::array<double, CADMUS_DIM> exchange_pivots_recv_per_dim {};
        std::array<double, CADMUS_DIM> request_sending_columns_per_dim {};
        std::array<double, CADMUS_DIM> send_columns_per_dim {};
        std::array<double, CADMUS_DIM> receive_columns_per_dim {};
        std::array<double, CADMUS_DIM> reduce_starting_from_per_dim {};
        std::array<double, CADMUS_DIM> update_global_pivots_receive_per_dim {};

        double local_reduction()              const { return std::accumulate(local_reduction_per_dim.begin(), local_reduction_per_dim.end(), 0.0); };
        double final_reduction()              const { return std::accumulate(final_reduction_per_dim.begin(), final_reduction_per_dim.end(), 0.0); };
        double exchange_pivots_send()         const { return std::accumulate(exchange_pivots_send_per_dim.begin(), exchange_pivots_send_per_dim.end(), 0.0); };
        double exchange_pivots_recv()         const { return std::accumulate(exchange_pivots_recv_per_dim.begin(), exchange_pivots_recv_per_dim.end(), 0.0); };
        double request_sending_columns()      const { return std::accumulate(request_sending_columns_per_dim.begin(), request_sending_columns_per_dim.end(), 0.0); };
        double send_columns()                 const { return std::accumulate(send_columns_per_dim.begin(), send_columns_per_dim.end(), 0.0); };
        double receive_columns()              const { return std::accumulate(receive_columns_per_dim.begin(), receive_columns_per_dim.end(), 0.0); };
        double reduce_starting_from()         const { return std::accumulate(reduce_starting_from_per_dim.begin(), reduce_starting_from_per_dim.end(), 0.0); };
        double update_global_pivots_receive() const { return std::accumulate(update_global_pivots_receive_per_dim.begin(), update_global_pivots_receive_per_dim.end(), 0.0); };

        double total() const { return total_wo_reduction; }

        friend std::ostream& operator<<(std::ostream&, const Stats&);
        std::string as_csv() const;

        std::string per_dimension() const;
    };

    struct Timer {
        std::chrono::high_resolution_clock timer_;
        std::chrono::time_point<decltype(timer_)> start_;

        Timer()
                :start_(timer_.now()) { }

        inline void reset() { start_ = timer_.now(); }

        inline double elapsed() const
        {
            std::chrono::duration<double> res = timer_.now() - start_;
            return res.count();
        }

        inline double elapsed_reset()
        {
            std::chrono::duration<double> res = timer_.now() - start_;
            reset();
            return res.count();
        }
    };

    // collect all the output in one string first,
    // works better under MPI
    inline std::ostream& operator<<(std::ostream& out, const Stats& s)
    {
        auto rep = fmt::format(
                "gid = {}, n_blocks = {}, n_columns_after_rearrange = {}, final_n_columns = {}, size_in_mb = {}, interior_cleared = {}, local_cleared = {}, global_cleared = {}, local_init = {:.2f}, interior_reduction = {:.2f}, sparsify = {:.2f}, diy::sort = {:.2f}, local_reduction = {:.2f},  final_reduction = {:.2f}, reduce_starting_from = {:.2f}, rearrange_send = {:.2f}, rearrange_recv = {:.2f}, exchange_pivots_send = {:.2f}, exchange_pivots_recv = {:.2f}",
                s.gid, s.n_blocks, s.n_columns_after_rearrange, s.final_n_columns, s.final_columns_in_mb, s.interior_cleared, s.local_cleared, s.global_cleared, s.local_init, s.interior_reduction, s.sparsify, s.final_sort,
                s.local_reduction(), s.final_reduction(), s.reduce_starting_from(), s.rearrange_send, s.rearrange_recv, s.exchange_pivots_send(), s.exchange_pivots_recv());
        out << rep;
        return out;
    }

    inline std::string Stats::per_dimension() const
    {
        std::stringstream ss;
        for(dim_type dim = 0; dim < CADMUS_DIM; ++dim) {
            ss << fmt::format(
                    "DIM = {}, gid={}, local_reduction={:.2f}, pivots send={:.2f}, receive={:.2f}, request sending cols={:.2f}, send columns={:.2f}, receive columns={:.2f}, reduce_starting_from={:.2f}, update global pivots={:.2f}\n",
                    dim, gid, local_reduction_per_dim[dim], exchange_pivots_send_per_dim[dim], exchange_pivots_recv_per_dim[dim], request_sending_columns_per_dim[dim],
                    send_columns_per_dim[dim], receive_columns_per_dim[dim], reduce_starting_from_per_dim[dim], update_global_pivots_receive_per_dim[dim]);
        }
        ss << "\n";

        return ss.str();
    }

    // collect all the output in one string first,
    // works better under MPI
    inline std::string Stats::as_csv() const
    {
        std::stringstream ss;
        ss << "FINAL_STATS_GID; n_blocks; n_simplices; clearing_opt; swap_reduction; local_init; interior_reduction; sparsify; exchange; final_sort; final_reduction; reduce_starting_from; total\n";
        ss << gid << "; " << n_blocks << "; " << n_simplices << "; ";
        if (clearing_opt)
            ss << "true;  ";
        else
            ss << "false; ";

        if (swap_reduction)
            ss << "true;  ";
        else
            ss << "false; ";

        ss << local_init << "; " << interior_reduction << "; " << sparsify << "; " << exchange << "; " << final_sort << "; " << final_reduction() << "; " << reduce_starting_from() << ";" << total() << "\n";
        return ss.str();
    }

    struct ColReductionStat {
        UidValue col_value;
        double elapsed;
        size_t original_size;
        size_t final_size;
        int n_additions;

        ColReductionStat(const UidValue& _col_value, double _elapsed, size_t _orig_size, size_t _fin_size, int _n_adds)
                :col_value(_col_value), elapsed(_elapsed), original_size(_orig_size), final_size(_fin_size), n_additions(_n_adds) { }

        ColReductionStat(const ColReductionStat&) = default;
        ColReductionStat& operator=(const ColReductionStat&) = default;

        bool operator<(const ColReductionStat& other) const
        {
            return std::tie(elapsed, n_additions, final_size, original_size, col_value) < std::tie(other.elapsed, other.n_additions, other.final_size, other.original_size, other.col_value);
        }

        bool operator<=(const ColReductionStat& other) const
        {
            return (*this < other) or (*this == other);
        }

        bool operator>(const ColReductionStat& other) const
        {
            return not(*this <= other);
        }

        bool operator>=(const ColReductionStat& other) const
        {
            return not(*this < other);
        }

        bool operator==(const ColReductionStat& other) const
        {
            return std::tie(col_value, elapsed, original_size, final_size, n_additions) == std::tie(other.col_value, other.elapsed, other.original_size, other.final_size, other.n_additions);
        }

        bool operator!=(const ColReductionStat& other) const
        {
            return not(*this == other);
        }
    };

    inline std::ostream& operator<<(std::ostream& out, const ColReductionStat& s)
    {
        auto rep = fmt::format("Column {} took {} sec to reduce, original size = {}, final_size = {}, n_additions = {}", s.col_value, s.elapsed, s.original_size, s.final_size, s.n_additions);
        out << rep;
        return out;
    }

    inline void print_memory_usage(int rank, std::string prefix = "")
    {
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long long mem_in_kb = usage.ru_maxrss;
        long long mem_in_mb = mem_in_kb / 1024;
        double mem_in_gb = static_cast<double>(mem_in_mb) / 1024.0;
        std::stringstream ss;
        if (not prefix.empty())
            ss << prefix << ", rank " << rank << ", max mem = " << mem_in_mb << " Mb = " << mem_in_gb << " Gb";
        else
            ss << "Rank " << rank << ", max mem = " << mem_in_mb << " Mb = " << mem_in_gb << " Gb";
        auto logger = spd::get("console");
        if (logger)
            logger->warn(ss.str());
        else
            std::cerr << ss.str() << std::endl;
    }

};
