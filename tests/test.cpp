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

#include <catch2/catch_test_macros.hpp>

#include "matrix.h"
#include "logger.h"
#include "coh_block.h"
#include "profile.h"
#include "types.h"


TEST_CASE("Merge-union", "[column]")
{
    cadmus::Column col_1;
    cadmus::Column col_2;
    cadmus::Column merged_col;
    cadmus::Column correct_merged_col;

    SECTION("empty cols") {
        merged_col = cadmus::Column::merge_union(col_1, col_2);
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("one empty col-1") {
        col_2.push_back(1);
        col_2.push_back(3);
        merged_col = cadmus::Column::merge_union(col_1, col_2);
        correct_merged_col = col_2;
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("one empty col-2") {
        col_1.push_back(1);
        col_1.push_back(3);
        merged_col = cadmus::Column::merge_union(col_1, col_2);
        correct_merged_col = col_1;
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("simple merge-no duplicates") {
        col_1.push_back(1);
        col_1.push_back(3);

        col_2.push_back(2);
        col_2.push_back(7);

        correct_merged_col.push_back(1);
        correct_merged_col.push_back(2);
        correct_merged_col.push_back(3);
        correct_merged_col.push_back(7);

        merged_col = cadmus::Column::merge_union(col_1, col_2);
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("simple merge-duplicates") {
        col_1.push_back(1);
        col_1.push_back(3);

        col_2.push_back(1);
        col_2.push_back(2);
        col_2.push_back(3);
        col_2.push_back(7);

        correct_merged_col.push_back(1);
        correct_merged_col.push_back(2);
        correct_merged_col.push_back(3);
        correct_merged_col.push_back(7);

        merged_col = cadmus::Column::merge_union(col_1, col_2);
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("simple merge-duplicates") {
        col_2.push_back(1);
        col_2.push_back(3);

        col_1.push_back(1);
        col_1.push_back(2);
        col_1.push_back(3);
        col_1.push_back(7);

        correct_merged_col.push_back(1);
        correct_merged_col.push_back(2);
        correct_merged_col.push_back(3);
        correct_merged_col.push_back(7);

        merged_col = cadmus::Column::merge_union(col_1, col_2);
        REQUIRE(merged_col == correct_merged_col);
    }

    SECTION("simple merge-duplicates") {
        col_1.push_back(1);
        col_1.push_back(2);
        col_1.push_back(3);
        col_1.push_back(7);

        correct_merged_col = col_2 = col_1;
        merged_col = cadmus::Column::merge_union(col_1, col_2);
        REQUIRE(merged_col == correct_merged_col);
    }

}

TEST_CASE("Matrix", "[matrix]")
{
    cadmus::Column col_1;
    cadmus::Column col_2;
    cadmus::Matrix m;
    cadmus::Matrix m_perp;
    cadmus::Matrix correct_m_perp;
    bool fill_columns, fill_rows;

//    SECTION("empty matrix") {
//        m_perp = m.antitranspose(fill_columns, fill_rows);
//        REQUIRE(m_perp == correct_m_perp);
//    }

    SECTION("simple matrix") {
        fill_rows = fill_columns = true;
        m.n_columns = 2;
        m.n_rows = 2;

        /*
         * 1 1   0 1
         * 1 0   1 1
         */

        col_1 = cadmus::Column({0, 1});
        col_2 = cadmus::Column({0});

        m = cadmus::Matrix({}, {col_1, col_2}, 2, 2);
        m_perp = m.antitranspose(fill_columns, fill_rows);
        correct_m_perp = cadmus::Matrix({{1}, {0, 1}}, {{1}, {0, 1}}, 2, 2);

        REQUIRE(m_perp == correct_m_perp);
    }

    SECTION("simple matrix-non-symmetric") {
        fill_rows = fill_columns = true;

        /*
         * x 0 0    0 1 1 0
         * 1 0 1    0 0 0 0
         * 1 0 1    0 1 1 x
         * 0 0 0
         */

        m = cadmus::Matrix({}, {{0, 1, 2}, {}, {1, 2}}, 4, 3);
        m_perp = m.antitranspose(fill_columns, fill_rows);
        correct_m_perp = cadmus::Matrix(/*rows*/ {{1, 2}, {}, {1, 2, 3}},
                /*cols*/ {{}, {0, 2}, {0, 2}, {2}}, 3, 4);

        REQUIRE(m_perp == correct_m_perp);
    }
}

TEST_CASE("Matrix cols rows", "[matrix]")
{
    SECTION("simple matrix-non-symmetric 1 col") {
        /*
         * 1
         * 0
         * 1
         * 0
         */

        int n_rows = 4, n_cols = 1;

        auto m_rows = cadmus::Matrix({{0}, {}, {0}, {}}, {}, n_rows, n_cols);
        auto m_cols = cadmus::Matrix({}, {{0, 2}}, n_rows, n_cols);
        auto m_full = cadmus::Matrix({{0}, {}, {0}, {}}, {{0, 2}}, n_rows, n_cols);

        auto orig_m_rows = m_rows;
        auto orig_m_cols = m_cols;

        m_rows.convert_to_columns(false);
        m_cols.convert_to_rows(false);

        REQUIRE(m_rows == orig_m_cols);
        REQUIRE(m_cols == orig_m_rows);

        m_rows.convert_to_rows(false);
        m_cols.convert_to_columns(false);

        REQUIRE(m_rows == orig_m_rows);
        REQUIRE(m_cols == orig_m_cols);

        m_rows.convert_to_columns(true);
        m_cols.convert_to_rows(true);

        REQUIRE(m_rows == m_full);
        REQUIRE(m_cols == m_full);

    }

    SECTION("simple matrix-non-symmetric 3 cols") {
        /*
         * 1 0 0
         * 0 0 0
         * 1 0 1
         * 0 0 1
         */

        int n_rows = 4, n_cols = 3;

        auto m_rows = cadmus::Matrix({{0}, {}, {0, 2}, {2}}, {}, n_rows, n_cols);
        auto m_cols = cadmus::Matrix({}, {{0, 2}, {}, {2, 3}}, n_rows, n_cols);
        auto m_full = cadmus::Matrix({{0}, {}, {0, 2}, {2}}, {{0, 2}, {}, {2, 3}}, n_rows, n_cols);

        auto orig_m_rows = m_rows;
        auto orig_m_cols = m_cols;

        m_rows.convert_to_columns(false);
        m_cols.convert_to_rows(false);

        REQUIRE(m_rows == orig_m_cols);
        REQUIRE(m_cols == orig_m_rows);

        m_rows.convert_to_rows(false);
        m_cols.convert_to_columns(false);

        REQUIRE(m_rows == orig_m_rows);
        REQUIRE(m_cols == orig_m_cols);

        m_rows.convert_to_columns(true);
        m_cols.convert_to_rows(true);

        REQUIRE(m_rows == m_full);
        REQUIRE(m_cols == m_full);

    }
}

TEST_CASE("CohBlock", "[sparsify]")
{
    auto logger = spd::get("console");
    if (!logger)
        logger = spd::stderr_color_mt("console");

    // Set the logger's log level
    logger->set_level(spd::level::debug);
    logger->flush_on(spd::level::trace);

    cadmus::CohBlock b = cadmus::CohBlock();

    SECTION("simple sparsification-1") {
        /* i_to_i       u_to_i
         * 1 1 1      0 1 0
         * 0 1 0      1 1 1
         * 0 0 1      0 1 0
         *
         * after sparsification
         *
         * 1 0 0     1 1 1
         * 0 1 0     1 1 1
         * 0 0 1     0 1 0
         *
         * after adding columns of i_to_i to u_to_i to kill more ones
         * column 0 to columns 1, 2
         * column 1 to column 2
         *
         * 1 0 0     0 0 0
         * 0 1 0     1 1 0
         * 0 0 1     0 1 0
          *
         */

        // UidValue ctor: {uid, dim, value}
        cadmus::UidValues int_col_values = {{10, 0, 50}, {11, 0, 40}, {12, 0, 30}};
        cadmus::UidValues u_col_values = {{13, 0, 45}, {14, 0, 44}, {15, 0, 23}};
        cadmus::UidValues int_row_values = {{8, 1, 20}, {7, 1, 10}, {6, 1, 0}};

        auto m_i_to_i = cadmus::Matrix({}, {{0}, {0, 1}, {0, 2}}, int_row_values.size(), int_col_values.size());
        auto m_u_to_i = cadmus::Matrix({}, {{1}, {0, 1, 2}, {1}}, int_row_values.size(), u_col_values.size());
        auto s_u_to_i = cadmus::Matrix({}, {{1}, {1, 2}, {}}, int_row_values.size(), u_col_values.size());

        cadmus::UltraSparseMatrixWithValues correct_sparsified_d = {{int_col_values[0], int_row_values[0]}, {int_col_values[1], int_row_values[1]}, {int_col_values[2], int_row_values[2]}};

        cadmus::SPUidValues pint_col_values = std::make_shared<cadmus::UidValues>(int_col_values);
        cadmus::SPUidValues pu_col_values = std::make_shared<cadmus::UidValues>(u_col_values);
        cadmus::SPUidValues pint_row_values = std::make_shared<cadmus::UidValues>(int_row_values);

        b.d_interior_ = cadmus::MatrixWithValues(std::move(m_i_to_i), pint_col_values, pint_row_values);
        b.d_union_to_interior_ = cadmus::MatrixWithValues(std::move(m_u_to_i), pu_col_values, pint_row_values);;

        b.sparsify();

        for(auto col : b.d_union_to_interior_.matrix.columns) {
            logger->debug("col {}", col);
        }

        REQUIRE(b.d_union_to_interior_.matrix == s_u_to_i);
        REQUIRE(b.sparsified_d_ == correct_sparsified_d);
    }


    SECTION("simple sparsification-1") {
        /* i_to_i       u_to_i
         * 1 1 1 0      0 1 0 1
         * 1 1 0 0      1 1 1 0
         * 1 0 0 0      0 1 0 1
         *
         * sequence of additions: add row 2 to rows 1 and 0, then add row 1 to row 0
         *
         * after sparsification
         *
         * 0 0 1 0     1 0 1 1
         * 0 1 0 0     1 0 1 1
         * 1 0 0 0     0 1 0 1
         *
         * after adding columns of i_to_id to u_to_i to kill more ones
         *
         * 0 0 1 0     1 0 1 0
         * 0 1 0 0     1 0 0 0
         * 1 0 0 0     0 0 0 0
         *
         *
         */

        // UidValue ctor: {uid, dim, value}
        cadmus::UidValues int_col_values = {{10, 0, 50}, {11, 0, 40}, {12, 0, 30}, {103, 0, 25}};
        cadmus::UidValues u_col_values = {{13, 0, 45}, {14, 0, 44}, {15, 0, 33}, {101, 0, 27}};
        cadmus::UidValues int_row_values = {{8, 1, 20}, {7, 1, 10}, {6, 1, 0}};

        auto m_i_to_i = cadmus::Matrix({}, {{0, 1, 2}, {0, 1}, {0}, {}}, int_row_values.size(), int_col_values.size());
        auto m_u_to_i = cadmus::Matrix({}, {{1}, {0, 1, 2}, {1}, {0, 2}}, int_row_values.size(), u_col_values.size());
        auto s_u_to_i = cadmus::Matrix({}, {{0, 1}, {}, {0}, {}}, int_row_values.size(), u_col_values.size());

        cadmus::UltraSparseMatrixWithValues correct_sparsified_d = {{int_col_values[0], int_row_values[2]}, {int_col_values[1], int_row_values[1]}, {int_col_values[2], int_row_values[0]}};

        cadmus::SPUidValues pint_col_values = std::make_shared<cadmus::UidValues>(int_col_values);
        cadmus::SPUidValues pu_col_values = std::make_shared<cadmus::UidValues>(u_col_values);
        cadmus::SPUidValues pint_row_values = std::make_shared<cadmus::UidValues>(int_row_values);


        b.d_interior_ = cadmus::MatrixWithValues(std::move(m_i_to_i), pint_col_values, pint_row_values);
        b.d_union_to_interior_ = cadmus::MatrixWithValues(std::move(m_u_to_i), pu_col_values, pint_row_values);;

        b.sparsify();

        for(auto col : b.d_union_to_interior_.matrix.columns) {
            logger->debug("col {}", col);
        }

        REQUIRE(b.sparsified_d_ == correct_sparsified_d);
        REQUIRE(b.d_union_to_interior_.matrix == s_u_to_i);
    }

}


//TEST_CASE("Cube", "[2D]")
//{
//    auto logger = spd::get("console");
//    if (!logger)
//        logger = spd::stderr_color_mt("console");
//
//    // Set the logger's log level
//    logger->set_level(spd::level::debug);
//    logger->flush_on(spd::level::trace);
//
//    cadmus::Point global_shape { 3, 4 };
//    cadmus::DynamicPoint core_min {0, 0, 0, 0};
//    cadmus::DynamicPoint core_max {3, 4, 0, 0};
//    cadmus::Bounds core {core_min, core_max};
//
//    std::vector<cadmus::Real> global_data { 1,  2,  3,  4,
//                                            5,  6,  7,  8,
//                                            9, 10, 11, 12 };
//
//    cadmus::GridRef global_domain {nullptr, global_shape, true };
//
//    cadmus::Cube::global_domain_ = global_domain;
//
//    cadmus::CohBlock b = cadmus::CohBlock();
//    b.core = core;
//    b.fab_ = cadmus::GridRef { global_data.data(), global_shape };
//
//    SECTION("cubical filtration-1") {
//
//    }
//
//}