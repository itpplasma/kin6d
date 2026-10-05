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
    GIT_TAG e589583a790e1b2dec662c10cca57e03739a0f3d)
FetchContent_MakeAvailable(fortplot)
# The first slice needs only FortNum's independently maintained interval modules.
# Compile them unchanged, from the pinned upstream tree; no copied numerics.
FetchContent_Declare(fortnum
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    GIT_TAG aa15a16e50f5c2e619a9a65a5100b2a994342949)
FetchContent_GetProperties(fortnum)
if(NOT fortnum_POPULATED)
    FetchContent_Populate(fortnum)
endif()
add_library(kin6d_intervals
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_rounding.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_interval.f90)
target_compile_options(kin6d_intervals PRIVATE -fno-fast-math -ffp-contract=off)
# Generated derivative families retain the primal ABI even for constant partials.
# Unused primal arguments are intentional; keep other diagnostics enabled.
set_source_files_properties(src/generated/kin6d_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_library(kin6d src/generated/kin6d_leaves.f90 src/laboratory.f90
    src/generated/kin6d_poisson_leaves.f90 src/poisson_laboratory.f90
    src/generated/kin6d_fp_leaves.f90 src/fp_laboratory.f90
    src/generated/kin6d_transport_leaves.f90 src/transport_laboratory.f90
    src/transport_certificate.f90 src/transport_sampling.f90)
set_source_files_properties(
    src/generated/kin6d_poisson_leaves.f90 src/generated/kin6d_fp_leaves.f90
    src/generated/kin6d_transport_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d PUBLIC kin6d_intervals)
target_compile_options(kin6d PRIVATE -Wall -Wextra -fno-fast-math -ffp-contract=off)
add_executable(kin6d_demo app/kin6d_demo.f90)
target_link_libraries(kin6d_demo PRIVATE kin6d fortplot)
add_executable(steady_laboratory_demo app/steady_laboratory_demo.f90)
target_link_libraries(steady_laboratory_demo PRIVATE kin6d fortplot)
add_executable(transport_demo app/transport_demo.f90)
target_link_libraries(transport_demo PRIVATE kin6d fortplot)
add_executable(formal_krook_certificate app/formal_krook_certificate.f90)
target_link_libraries(formal_krook_certificate PRIVATE kin6d)
add_executable(test_laboratory test/test_laboratory.f90)
target_link_libraries(test_laboratory PRIVATE kin6d)
add_test(NAME characteristic_and_steady_oracles COMMAND test_laboratory)
add_executable(test_poisson test/test_poisson.f90)
target_link_libraries(test_poisson PRIVATE kin6d)
add_test(NAME poisson_and_coupled_steady_oracles COMMAND test_poisson)
add_executable(test_fp test/test_fp.f90)
target_link_libraries(test_fp PRIVATE kin6d)
add_test(NAME bounded_fp_steady_oracles COMMAND test_fp)
add_executable(test_transport test/test_transport.f90)
target_link_libraries(test_transport PRIVATE kin6d)
add_test(NAME general_transport_and_tangent_oracles COMMAND test_transport)
add_executable(test_transport_sampling test/test_transport_sampling.f90)
target_link_libraries(test_transport_sampling PRIVATE kin6d)
add_test(NAME cached_transport_band_oracles COMMAND test_transport_sampling)
find_program(KIN6D_LAKE lake)
option(KIN6D_CHECK_LEAN "Check continuous Lean theorems and computed certificate" OFF)
if(KIN6D_CHECK_LEAN)
    if(NOT KIN6D_LAKE)
        message(FATAL_ERROR "Lean/Lake is required; see proof/README.md for setup")
    endif()
    add_test(NAME lean_continuous_certificates
        COMMAND bash ${CMAKE_SOURCE_DIR}/proof/check.sh
            $<TARGET_FILE:formal_krook_certificate>)
    set_tests_properties(lean_continuous_certificates PROPERTIES
        TIMEOUT 300 PROCESSORS 2)
endif()
set(KIN6D_FORTSYM_BUILD "" CACHE PATH "Optional FortSym build tree for regeneration")
set(KIN6D_FORTSYM_REVISION "93c8b2ada93a2f58e2c36349560e9c52f2c9f79a")
if(KIN6D_FORTSYM_BUILD)
    file(STRINGS "${KIN6D_FORTSYM_BUILD}/CMakeCache.txt" fortsym_source
        REGEX "^CMAKE_HOME_DIRECTORY:INTERNAL=")
    string(REPLACE "CMAKE_HOME_DIRECTORY:INTERNAL=" "" fortsym_source
        "${fortsym_source}")
    execute_process(COMMAND git -C "${fortsym_source}" rev-parse HEAD
        OUTPUT_VARIABLE fortsym_revision OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE fortsym_revision_status)
    if(NOT fortsym_revision_status EQUAL 0 OR
            NOT fortsym_revision STREQUAL KIN6D_FORTSYM_REVISION)
        message(FATAL_ERROR "Regeneration requires FortSym ${KIN6D_FORTSYM_REVISION}")
    endif()
    foreach(generator laboratory poisson fp transport)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_include_directories(gen_${generator} PRIVATE
            ${KIN6D_FORTSYM_BUILD}/include)
        target_link_libraries(gen_${generator} PRIVATE
            ${KIN6D_FORTSYM_BUILD}/lib/libfortsym.so)
    endforeach()
endif()
# Explicit performance audit; not part of the ordinary correctness gate.
add_executable(benchmark_streaming app/benchmark_streaming.f90)
target_link_libraries(benchmark_streaming PRIVATE kin6d)
add_executable(benchmark_fp app/benchmark_fp.f90)
target_link_libraries(benchmark_fp PRIVATE kin6d)
add_executable(benchmark_transport app/benchmark_transport.f90)
target_link_libraries(benchmark_transport PRIVATE kin6d)
