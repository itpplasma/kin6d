# Opt-in CPU precision experiment. Production residuals/estimators are untouched.
option(KIN6D_BUILD_CPU_PRECISION_BENCHMARK "Build experimental CPU precision benchmark" OFF)
if(NOT KIN6D_BUILD_CPU_PRECISION_BENCHMARK)
    return()
endif()
if(NOT DEFINED KIN6D_CPU_PRECISION_ROOT)
    set(KIN6D_CPU_PRECISION_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
endif()
add_executable(benchmark_cpu_precision
    ${KIN6D_CPU_PRECISION_ROOT}/app/benchmark_cpu_precision.f90
    ${KIN6D_CPU_PRECISION_ROOT}/src/generated/kin6d_gs_cpu_precision.f90
    ${KIN6D_CPU_PRECISION_ROOT}/src/generated/kin6d_map_cpu_precision.f90)
# The numerical experiment uses published runtime APIs only.
if(NOT TARGET fortgen_runtime)
    include(FetchContent)
    set(FORTGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(kin6d_cpu_fortgen
        GIT_REPOSITORY https://github.com/lazy-fortran/fortgen.git
        GIT_TAG d7946a1adad695c5f68823a9f0650a9174eb32fb)
    FetchContent_MakeAvailable(kin6d_cpu_fortgen)
endif()
target_link_libraries(benchmark_cpu_precision PRIVATE
    OpenMP::OpenMP_Fortran fortgen_runtime)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Keep the large Hessian leaf inside its SIMD loop. Default inlining
    # heuristics otherwise leave this loop scalar (see benchmark recipe).
    target_compile_options(benchmark_cpu_precision PRIVATE
        -finline-limit=100000)
endif()
# Native behavioral checks are explicit because this benchmark is opt-in.
foreach(kernel gs map)
    add_test(NAME cpu_precision_${kernel}_comparison
        COMMAND benchmark_cpu_precision ${kernel} 4096 2)
    set_tests_properties(cpu_precision_${kernel}_comparison PROPERTIES
        PROCESSORS 2
        ENVIRONMENT "OMP_NUM_THREADS=2;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1")
endforeach()
if(NOT TARGET kin6d_mhd)
    set_tests_properties(cpu_precision_gs_comparison cpu_precision_map_comparison
        PROPERTIES LABELS "mhd;performance")
    return()
endif()
add_executable(benchmark_equilibrium_cpu
    app/benchmark_equilibrium_cpu.f90 test/support/ideal_map_cases.f90)
target_link_libraries(benchmark_equilibrium_cpu PRIVATE kin6d_mhd)
# The oracle is also compiled by an existing test executable. Keep this
# executable's support modules private to avoid duplicate module producers.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cpu-profile-modules)
set_target_properties(benchmark_equilibrium_cpu PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cpu-profile-modules)
target_include_directories(benchmark_equilibrium_cpu PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/cpu-profile-modules)
# Reuse the production precision/refinement APIs. Compiling older source
# slices here would duplicate module providers and reintroduce trampolines.
add_library(kin6d_cpu_precision_framework INTERFACE)
target_link_libraries(kin6d_cpu_precision_framework INTERFACE
    kin6d_ode fortnum fortsparse fortgen_runtime OpenMP::OpenMP_Fortran)
add_executable(benchmark_cpu_solve app/benchmark_cpu_solve.f90)
target_link_libraries(benchmark_cpu_solve PRIVATE
    kin6d_equilibrium kin6d_cpu_precision_framework)
add_test(NAME cpu_precision_gs_ir COMMAND benchmark_cpu_solve 24 2)
set_tests_properties(cpu_precision_gs_ir PROPERTIES PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
