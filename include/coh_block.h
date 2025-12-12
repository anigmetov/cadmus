#pragma once

#include <utility>
#include <tuple>
#include <vector>

#include <diy/serialization.hpp>

#include "logger.h"
#include "types.h"
#include "messages.h"
#include "profile.h"
#include "simplex.h"
#include "cube.h"
#include "filtration.h"
#include "cubical_filtration.h"
#include "matrix_with_values.h"
#include "values_only_matrix.h"


namespace cadmus {

    VertexVecVec get_fr_displacements(size_t d);

    struct CohBlock {
        enum class FilType {
            Cubical,
            Freudenthal
        };

//    using Field = Field<2>;

        using GidToPivotMsg = std::map<int, PivotMessage>;
        using Gids = std::vector<int>;

        // data
        bool negate_ {false};
        bool wrap_ {false};

        GridRef fab_;
        Grid fab_storage_;

        // used for converting Vertex to global Index and back
        // doesn't store any data (data_ is nullptr)
        GridRef global_domain_;

        Bounds core {CADMUS_DIM};

        // for simplex i (in filtration order) return the corresponding index (C order)
        std::vector<size_t> simplex_to_global;

        // components of coboundary matrix:
        // columns of interior cells (coboundary automatically in the interior)
        MatrixWithValues d_interior_;
        // columns of K_i \cap union cells, rows of interior cells
        MatrixWithValues d_union_to_interior_;
        // columns of K_i \cap union cells, rows of K_i \cap union cells
        MatrixWithValues d_union_to_union_;

        // used in final reduction only
        ValueMatrices global_r_;
        UidValues new_lows_;

        // only one of the two is used
        std::vector<Simplex> simplices;
        std::vector<Cube> cubes;

        std::unordered_set<Int> union_uids_;

        UltraSparseMatrixWithValues sparsified_d_;

        bool is_root_ {false};

        int n_blocks_ {0};

        bool clearing_opt_ {true};
        bool swap_reduction_ {true};

        std::set<size_t> addition_cache_;

        Stats stats_;

        // for rearrange only
        ValueMatrix u_to_u_local_;
        ValueMatrix u_to_i_local_;
        UltraSparseMatrixWithValues local_sparse_;

        Timer timer_total_;

        // helpers to get unique id
        int total_n_cells_per_cube_;   // total number of simplices of Freudenthal triangulation in a single cube
        std::map<int, int> n_cells_in_lower_dimension_; // for each key d, stores total number of simplices of dimension < d in a single cube

        CohBlock()
                :
                fab_(nullptr, Shape()),
                fab_storage_(Shape()),
                global_domain_(nullptr, Shape())
        {
            count_simplices_in_cube();
        }

        CohBlock(Real* data, const Shape& shape, const Shape& global_domain_shape, int n_blocks)
                :
                fab_(data, shape),
                fab_storage_(shape),
                global_domain_(nullptr, global_domain_shape),
                is_root_(false),
                n_blocks_(n_blocks)
        {
            count_simplices_in_cube();
        }

        ~CohBlock()
        {
            if (fab_.data() != fab_storage_.data())
                delete[] fab_.data();
        }

        static void* create()
        {
            return new CohBlock();
        }

        static void destroy(void* b)
        {
            delete static_cast<CohBlock*>(b);
        }

        void count_simplices_in_cube();

        void init(FilType fil_type, bool negate, bool clearing_opt, bool swap_reduction);

        bool is_union(const Vertex& v_local) const;

        void init_coboundary_matrix_local(FilType fil_type);

        void reduce_interior();

        void sparsify();

        // lower-star utils

        SimplicialFiltration get_freudenthal_filtration(dim_type top_d = CADMUS_DIM);
        CubicalFiltration get_cubical_filtration(dim_type top_d = CADMUS_DIM);

        void add_freudenthal_simplices(dim_type d);

        void add_freudenthal_simplices_from_vertex(const Vertex& v, dim_type d, const VertexVecVec& disps);

        Int global_id(const IdxVector& vertices, Point v, size_t deltas_idx) const;

        // cubical LS utils

        std::tuple<bool, bool, Real> get_cube_validity_interior_and_value(Int cube_id);

        void save_diagrams(std::string fname_out, int gid);
        Vertex from_local_to_global(const Vertex& v);
        Vertex from_global_to_local(const Vertex& v);

        // distributed algorithm routines
        std::vector<Real> all_values_;
        std::vector<Real> boundaries_;

        int gid_by_row_value(const UidValue& val) const;
        int gid_by_col_value(const UidValue& val) const;

        std::unordered_map<int, ValueMatrix> split_matrix_by_gids(MatrixWithValues& m);
        std::unordered_map<int, UltraSparseMatrixWithValues> split_ultrasparse_matrix_by_gids();

        void rearrange_matrix_send(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner);
        void rearrange_matrix_receive(const diy::Master::ProxyWithLink& cp);

        void send_columns(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, dim_type dim, int round);
        void receive_columns(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, dim_type dim);

        SimplicialFiltration freudenthal_fil_;
        CubicalFiltration    cubical_fil_;

        void add_to_cached(const Column& pivot_col);
        void load_to_cache(const Column& col);
        void load_from_cache(Column& col);

        void cached_drop_low();

        bool check_sparsified_d(bool check_sorting);

        std::map<dim_type, UidValue> min_values_;
        std::map<dim_type, UidValue> max_values_;
        std::map<dim_type, UidValue> global_min_values_;
        std::map<dim_type, UidValue> global_max_values_;
        Int                          global_min_vertex_uid_ { -1 };


        PersistencePairs read_off_persistence();

        void get_min_max_cells();

    };

    void print_cob(const MatrixWithValues& d, const SimplicialFiltration& fil);

    template<class T>
    void enqueue(const diy::Master::ProxyWithLink& cp, const diy::Assigner& assigner, int gid, const T& msg)
    {
        int proc = assigner.rank(gid);
        diy::BlockID bid {gid, proc};
        cp.enqueue(bid, msg);
    }

} // namespace cadmus