cmake_minimum_required(VERSION 3.24)
if(NOT CMAKE_CUDA_ARCHITECTURES)
    message(FATAL_ERROR "Set CMAKE_CUDA_ARCHITECTURES explicitly for this prototype")
endif()
project(kin6d_gpu_prototype LANGUAGES Fortran CXX CUDA)
# Standalone experiment; it defines no production backend/configuration policy.
find_package(OpenMP REQUIRED COMPONENTS Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
add_executable(benchmark_equilibrium_gpu benchmark.f90 gpu.cu
    generated/ideal_map_gpu_cpu.f90
    ../../src/generated/kin6d_ideal_map_leaves.f90)
target_link_libraries(benchmark_equilibrium_gpu PRIVATE OpenMP::OpenMP_Fortran)
set_target_properties(benchmark_equilibrium_gpu PROPERTIES CUDA_STANDARD 17)
target_compile_options(benchmark_equilibrium_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-O3;-fno-fast-math;-ffp-contract=off>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3;--ptxas-options=-v>)
enable_testing()
add_test(NAME gpu_equilibrium_leaf_oracles COMMAND benchmark_equilibrium_gpu 257 3)
set_tests_properties(gpu_equilibrium_leaf_oracles PROPERTIES LABELS "mhd;fast"
    PROCESSORS 4 ENVIRONMENT
    "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
