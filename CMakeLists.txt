cmake_minimum_required(VERSION 3.22)
project(kin6d LANGUAGES Fortran C)
include(CTest)
include(FetchContent)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
set(FORTPLOT_BUILD_FORTARRAY_ADAPTER OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG 3abf42c9a3c1c7026e1eda5e4aea8d607f2b8f1e)
FetchContent_MakeAvailable(fortplot)
# The first slice needs only FortNum's independently maintained interval modules.
# Compile them unchanged, from the pinned upstream tree; no copied numerics.
FetchContent_Declare(fortnum
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    GIT_TAG 1a52e13c5e28b09d083947a6ec7e8a1dcf9d502e)
FetchContent_GetProperties(fortnum)
if(NOT fortnum_POPULATED)
    FetchContent_Populate(fortnum)
endif()
add_library(kin6d_intervals
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_rounding.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_interval.f90)
target_compile_options(kin6d_intervals PRIVATE -fno-fast-math -ffp-contract=off)
add_library(kin6d src/generated/kin6d_leaves.f90 src/laboratory.f90)
target_link_libraries(kin6d PUBLIC kin6d_intervals)
target_compile_options(kin6d PRIVATE -Wall -Wextra -fno-fast-math -ffp-contract=off)
add_executable(kin6d_demo app/kin6d_demo.f90)
target_link_libraries(kin6d_demo PRIVATE kin6d fortplot)
add_executable(test_laboratory test/test_laboratory.f90)
target_link_libraries(test_laboratory PRIVATE kin6d)
add_test(NAME characteristic_and_steady_oracles COMMAND test_laboratory)
set(KIN6D_FORTSYM_BUILD "" CACHE PATH "Optional FortSym build tree for regeneration")
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_laboratory gen/gen_laboratory.f90)
    target_include_directories(gen_laboratory PRIVATE ${KIN6D_FORTSYM_BUILD}/include)
    target_link_libraries(gen_laboratory PRIVATE
        ${KIN6D_FORTSYM_BUILD}/lib/libfortsym.so)
endif()
# Explicit performance audit; not part of the ordinary correctness gate.
add_executable(benchmark_streaming app/benchmark_streaming.f90)
target_link_libraries(benchmark_streaming PRIVATE kin6d)
