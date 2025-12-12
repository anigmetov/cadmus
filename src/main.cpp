#include <iostream>
#include <vector>
#include <unistd.h>

#include <diy/master.hpp>
#include <diy/proxy.hpp>
#include <diy/decomposition.hpp>
#include <diy/mpi.hpp>
#include <diy/algorithms.hpp>
#include <diy/io/numpy.hpp>
#include <diy/mpi/operations.hpp>

#include <opts/opts.h>

#include "logger.h"
#include "profile.h"
#include "types.h"
#include "coh_block.h"


//constexpr int D = 3;
using Real = cadmus::Real;
using Int = cadmus::Int;

using namespace cadmus;


Decomposer read_from_npy_file(std::string fname_in,
        diy::mpi::communicator& world,
        int nblocks,
        diy::Master& master,
        diy::ContiguousAssigner& assigner,
        diy::MemoryBuffer& header,
        Bounds& domain,
        WrapVector wrap,
        bool negate)
{
    auto logger = spd::get("console");

    diy::mpi::io::file in(world, fname_in, diy::mpi::io::file::rdonly);
    diy::io::NumPy reader(in);
    reader.read_header();

    if (reader.word_size() != sizeof(Real)) {
        throw std::runtime_error("Type mismatch");
    }

    if (reader.shape().size() != CADMUS_DIM)
        throw std::runtime_error("compiled with wrong dimension, cannot read file");

    Int cube_part = 1 << CADMUS_DIM;
    if (std::accumulate(reader.shape().begin(), reader.shape().end(), 1, [](auto x, auto y) { return x * y; }) > std::numeric_limits<Int>::max() / cube_part)
        throw std::runtime_error("data too big: cannot fit vertex id in the first bits");

    domain.min = {0, 0, 0, 0};
    domain.max = {0, 0, 0, 0};
    DynamicPoint one = DynamicPoint::one(CADMUS_DIM);

    for(unsigned i = 0; i < CADMUS_DIM; ++i) {
        domain.max[i] = reader.shape()[i] - 1;
    }

    if (logger) logger->debug("read domain, domain = {}", domain);

    Decomposer decomposer(CADMUS_DIM,
            domain,
            nblocks,
            Decomposer::BoolVector {true, true, true},   // share_face
            wrap,
            Decomposer::CoordinateVector {0, 0, 0});

    decomposer.decompose(world.rank(), assigner, [&master, &wrap, &reader, one, &logger, negate, nblocks](int gid,
            const Bounds& core,
            const Bounds& bounds,
            const Bounds& domain,
            const Link& link) {
      auto* b = new cadmus::CohBlock;

      b->stats_.n_blocks = nblocks;

      auto shape_4d = bounds.max - bounds.min + one;
      Shape shape(&shape_4d[0]);
      bool c_order = true;
      b->fab_storage_ = decltype(b->fab_storage_)(shape, c_order);
      b->fab_ = decltype(b->fab_)(b->fab_storage_.data(), shape, c_order);
      b->core = core;

      Point glob_domain_shape;
      for(unsigned i = 0; i < CADMUS_DIM; ++i)
          glob_domain_shape[i] = domain.max[i] + 1;
      b->global_domain_ = GridRef(nullptr, glob_domain_shape);

      if (logger) logger->debug("core = {}, bounds = {}, domain = {}, wrap = {}", core, bounds, glob_domain_shape, wrap);
      reader.read(core, b->fab_.data());

      if (negate) {
          b->fab_ /= -1;
      }

      if (logger) logger->debug("read core, core = {}", core);
      if (logger) logger->trace("link: neighbors.size = {}, ", link.neighbors().size());
      for(auto i = 0; i < link.neighbors().size(); ++i)
          if (logger) logger->trace("link: i = {}, gid = {}, proc = {}", i, link.neighbors()[i].gid, link.neighbors()[i].proc);

      master.add(gid, b, link);
    });

    return decomposer;
}

int main(int argc, char** argv)
{

    Timer timer;
    Timer timer_local;

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    using opts::Option;
    using opts::PosOption;

    opts::Options ops;

    int nthreads = 1;
    int nblocks = -1;
    int num_samples = 10;
    int in_memory = -1; // max # blocks to store in memory
    std::string log_level = "info";

    bool negate {false}, help, clearing {false}, swap_reduction {false}, freudenthal {false};
    ops
            >> Option('n', "negate", negate, "negate the function")
            >> Option('s', "swap", swap_reduction, "use swap_reduction optimization")
            >> Option('f', "freudenthal", freudenthal, "use simplicial (Freudenthal) complexes")
            >> Option('c', "clear", clearing, "use clearing_opt optimization")
            >> Option('a', "num-samples", num_samples, "number of samples for sample sort")
            >> Option('b', "blocks", nblocks, "number of DIY blocks")
            >> Option('l', "logger", log_level, "level for the logger output (trace, debug, info, ...)")
            >> Option('h', "help", help, "show help message");


    std::string fname_in, fname_out;
    if (!ops.parse(argc, argv) || help ||
            !(ops >> PosOption(fname_in) >> PosOption(fname_out))) {
        std::cout << "Usage: " << argv[0] << " [options] INFILE OUTFILE\n\n";
        std::cout << "Distributed cohomology\n\n";
        std::cout << ops << std::endl;
        return 1;
    }

    if (nblocks < 0)
        nblocks = world.size();

    // Create a logger with a console sink
    auto logger = spd::stderr_color_mt("console");
    if (logger) logger->set_pattern("[%H:%M:%S.%e] [PID %P] [%l] %v");
    // Set the logger's logger level

    if (logger) logger->set_level(log_level_from_string(log_level));
    if (logger) logger->flush_on(log_level_from_string(log_level));

    if (logger) logger->warn("world.size = {}, nblocks = {}", world.size(), nblocks);

    CohBlock::FilType fil_type = freudenthal ? CohBlock::FilType::Freudenthal : CohBlock::FilType::Cubical;

    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);

    int k = 2; // use binary tree for reduction

    if (logger) logger->info("Started, rank = {}, pid = {}, node = {}", world.rank(), getpid(), std::string(processor_name));
    if (logger and world.rank() == 0) logger->info("negate = {}, swap = {}, clearing_opt = {}", negate, swap_reduction, clearing);

    diy::Master master(world, nthreads, in_memory, &cadmus::CohBlock::create, &cadmus::CohBlock::destroy);

    bool is_rank_0 = world.rank() == 0;

    if (logger)
        logger->debug("Created master, negate = {}.", negate);
    else {
        if (is_rank_0) std::cerr << "created master" << std::endl;
    }

    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain(CADMUS_DIM);
    cadmus::WrapVector wrap;

    auto decomposer = read_from_npy_file(fname_in, world, nblocks, master, assigner, header, domain, wrap, negate);

    // set static variable in Cube class to global domain
    {
        Point glob_domain_shape;
        for(unsigned i = 0; i < CADMUS_DIM; ++i)
            glob_domain_shape[i] = domain.max[i] + 1;
        cadmus::Cube::global_domain_ = GridRef(nullptr, glob_domain_shape);
    }

//    print_memory_usage(world.rank(), "Right after reading file ");

    if (logger)
        logger->info("Data read");
    else {
        if (is_rank_0) std::cerr << "Data read" << std::endl;
    }


    // set up matrix
    master.foreach([negate, clearing, swap_reduction, fil_type](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->stats_.gid = cp.gid();
      b->init(fil_type, negate, clearing, swap_reduction);
    });

    print_memory_usage(world.rank(), "After init");

    world.barrier();

    CALI_MARK_BEGIN("reduction_kernel");

    timer.reset();

    master.foreach([](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->timer_total_.reset();
      b->reduce_interior();
      b->sparsify();
    });

//    print_memory_usage(world.rank(), "After local reduction and sparsification");

    // determine global minimum

    bool contiguous = true;

    diy::RegularMergePartners partners(decomposer, k, contiguous);

    diy::reduce(master, assigner, partners, [](CohBlock* b, const diy::ReduceProxy& rp, const diy::RegularMergePartners& ps) {
      for(int i = 0; i < rp.in_link().size(); ++i) {
          BlockID from_bid = rp.in_link().target(i);
          if (from_bid.gid == rp.gid())
              continue;
          for(dim_type d = 0; d < CADMUS_DIM; ++d) {
              UidValue in_min_value, in_max_value;
              rp.dequeue(from_bid, in_min_value);
              rp.dequeue(from_bid, in_max_value);
              b->global_min_values_[d] = std::min(b->global_min_values_[d], in_min_value);
              b->global_max_values_[d] = std::max(b->global_max_values_[d], in_max_value);
          }
      }

      assert(rp.out_link().size() <= 1);

      if (rp.out_link().size() == 1) {
          BlockID to_bid = rp.out_link().target(0);

          for(dim_type d = 0; d < CADMUS_DIM; ++d) {
              rp.enqueue(to_bid, b->global_min_values_[d]);
              rp.enqueue(to_bid, b->global_max_values_[d]);
          }
      }

    });

    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      if (cp.gid() == 0) {
          auto logger = spd::get("console");
          for(dim_type d = 0; d < CADMUS_DIM; ++d) {
              logger->info("In dimension {}, global min = {}, max = {}", d, b->global_min_values_[d], b->global_max_values_[d]);
          }
      }
    });

    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      if (cp.gid() == 0) {
          for(dim_type d = 0; d < CADMUS_DIM; ++d) {
              for(int to_gid = 1; to_gid < nblocks; ++to_gid) {
                  enqueue(cp, assigner, to_gid, b->global_min_values_[d]);
                  enqueue(cp, assigner, to_gid, b->global_max_values_[d]);
              }
          }
      }
    });

    master.exchange(true);

    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      if (cp.gid() != 0) {
          int from_gid = 0;
          for(dim_type d = 0; d < CADMUS_DIM; ++d) {
              cp.dequeue(from_gid, b->global_min_values_[d]);
              cp.dequeue(from_gid, b->global_max_values_[d]);
          }
      }
      b->global_min_vertex_uid_ = b->global_min_values_[0].uid;
    });

//    print_memory_usage(world.rank(), "After global minimum");

    timer_local.reset();

    diy::sort<CohBlock, Real, std::greater<Real>>(master, assigner,
            &CohBlock::all_values_,
            &CohBlock::boundaries_, num_samples, std::greater<Real>(), k);

//    print_memory_usage(world.rank(), "After diy::sort");

    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->stats_.final_sort = timer_local.elapsed();
    });

    timer_local.reset();

    // redistribute matrix across ranks
    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->rearrange_matrix_send(cp, assigner);
      b->stats_.rearrange_send = timer_local.elapsed();
    });

    timer_local.reset();

//    print_memory_usage(world.rank(), "After rearrange_send");

    master.exchange(true);

//    print_memory_usage(world.rank(), "After exchange");

    master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->rearrange_matrix_receive(cp);
      b->stats_.rearrange_recv = timer_local.elapsed();
      b->stats_.n_columns_after_rearrange = b->global_r_.size();
    });

    timer_local.reset();

//    print_memory_usage(world.rank(), "After rearrange_recv");

    // after that matrix was converted to hybrid format (ultrasparse + CSR), all entries are indexed by UidValue

    Timer timer_dim;

    for(dim_type dim = 0; dim < CADMUS_DIM; ++dim) {
        timer_dim.reset();

        timer_local.reset();

        // main reduction loop

        bool remote = true;

        int round = 0;

        while(true) {

            // in the first round, all columns will be redistributed by low
            // after that, send columns will only send new lows
            master.collectives().clear();
            master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) { b->send_columns(cp, assigner, dim, round); });
            master.exchange(remote);

            // blocks write the number of columns they need to send to local variables and use collective all_reduce
            int n_not_done = master.proxy(master.loaded_block()).read<int>();
            if (logger)
                logger->info("dim = {}, round = {}, n_not_done = {}", dim, round, n_not_done);
            else if (world.rank() == 0) {
                std::cerr << fmt::format("dim = {}, round = {}, n_not_done = {}", dim, round, n_not_done) << std::endl;
            }

            if (n_not_done == 0)
                break;

            master.collectives().clear();
            master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) { b->receive_columns(cp, assigner, dim); });
            master.exchange(remote);

            if (logger)
                logger->info("send_columns done");
            else {
                if (is_rank_0) std::cerr << "send_columns done, dim = " << dim << ", round = " << round << std::endl;
            }
            round++;
        }

//        print_memory_usage(world.rank(), fmt::format("After global reduction loop in dimension {}", dim));

//        master.foreach([&](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) { b->global_r_[dim].sparsify_reduced(); });

//        print_memory_usage(world.rank(), fmt::format("After sparsify_reduced in dimension {}", dim));

        CALI_MARK_BEGIN("global_clearing");

        if (clearing and dim < CADMUS_DIM - 1) {
            // globally distribute positive uids
            // and perform clearing optimization
            diy::all_to_all(master, assigner,
                    [&](cadmus::CohBlock* b, const diy::ReduceProxy& rp) {
                      if (rp.in_link().size() == 0) {
                          std::map<int, UidValues> gid_to_positive_vals;

                          for(auto& col : b->global_r_[dim]) {
                              UidValue positive_uid_value = col.second.low();
                              int dest_gid = b->gid_by_row_value(positive_uid_value);
                              gid_to_positive_vals[dest_gid].push_back(positive_uid_value);
                          }

                          for(auto&[dest_gid, positive_values] : gid_to_positive_vals) {
                              diy::BlockID dest = rp.out_link().target(dest_gid);
                              assert(dest.gid == dest_gid);
                              rp.enqueue(dest, positive_values);
                          }
                      } else {
                          for(int i = 0; i < rp.in_link().size(); ++i) {
                              int gid = rp.in_link().target(i).gid;
                              assert(gid == i);
                              if (rp.incoming(gid).size()) {
                                  UidValues positive_values;
                                  rp.dequeue(gid, positive_values);

                                  for(UidValue pos_val : positive_values) {
                                      size_t n_erased = b->global_r_[dim + 1].erase(pos_val, true);
                                      b->stats_.global_cleared += n_erased;
                                      if (n_erased) {
                                          logger->debug("global clearing: gid = {}, cleared {}, global cleared = {}", b->stats_.gid, pos_val, b->stats_.global_cleared);
                                      }
                                  }
                              }
                          }
                      }
                    }, // end of all-to-all lambda
                    k);
        }

        print_memory_usage(world.rank(), fmt::format("After global clearing in dimension {}", dim));

        CALI_MARK_END("global_clearing");
    }

    CALI_MARK_END("reduction_kernel");

    if (logger)
        logger->warn("Done reducing, TIME={}", timer.elapsed());
    else if (is_rank_0)
        std::cerr << "Done reducing, TIME=" << timer.elapsed() << std::endl;

    print_memory_usage(world.rank(), "After reduction");

    master.foreach([logger](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      b->stats_.final_n_columns = b->global_r_.size();
      b->stats_.final_columns_in_mb = b->global_r_.total_column_size_in_megabytes();
      if (logger) {
          logger->warn("{}", b->stats_);
      } else {
          std::cerr << b->stats_ << std::endl;
      }
    });

    master.foreach([logger](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) {
      if (logger) {
          logger->warn("{}", b->stats_.per_dimension());
      } else {
          std::cerr << b->stats_.per_dimension() << std::endl;
      }
    });

    if (fname_out != "/dev/null")
        master.foreach([fname_out](cadmus::CohBlock* b, const diy::Master::ProxyWithLink& cp) { b->save_diagrams(fname_out, cp.gid()); });
    return 0;
}
