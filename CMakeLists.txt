cmake_minimum_required(VERSION 3.10)
project(cadmus)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE)

option(cadmus_caliper        "Enable profiling with Caliper"     OFF)
option(cadmus_spdlog         "Enable logging with spdlog"        OFF)
option(cadmus_per_col_stats  "Detailed timing for each column"   OFF)
option(cadmus_march          "Enable march=native in Release"    ON)
option(cadmus_backward       "Enable Backward"    OFF)

find_path(DIY_INCLUDE_DIR diy/master.hpp)

find_package(MPI REQUIRED)

find_package(Boost REQUIRED)

if(cadmus_backward)
    add_subdirectory(extern/backward-cpp)
endif()

set(libraries MPI::MPI_C)

if(cadmus_caliper)
    find_package(Caliper REQUIRED)
    add_compile_definitions (CADMUS_USE_CALIPER)
    set(libraries ${libraries} caliper)
endif()

if(cadmus_spdlog)
    add_compile_definitions (CADMUS_USE_SPDLOG)
    if(cadmus_per_col_stats)
        add_compile_definitions (CADMUS_PER_COL_STATS)
    endif()
else()
    remove_definitions (-DDIY_USE_SPDLOG)
endif()

if (cadmus_march)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        include(CheckCXXCompilerFlag)

        CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
        endif()
    endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/extern/opts/include
                    ${DIY_INCLUDE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/extern/)

add_subdirectory(src)

find_package(Catch2 3)

if (Catch2_FOUND)
    add_subdirectory(tests)
endif()
