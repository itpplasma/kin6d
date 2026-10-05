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
    GIT_TAG 87ff7a43809e868e0153f400a21cdb31b9a92d84)
FetchContent_MakeAvailable(fortplot)
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ftrampoline-impl=heap" KIN6D_HAVE_HEAP_TRAMPOLINES)
if(KIN6D_HAVE_HEAP_TRAMPOLINES)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ftrampoline-impl=heap>)
endif()
# Compile FortNum's interval and ODE modules unchanged from the pinned tree.
# The application owns physical laws; generic numerics remain upstream.
FetchContent_Declare(fortnum
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    GIT_TAG 015bb8b1490d99d55493b3936e0dd6adf77bc293)
FetchContent_GetProperties(fortnum)
if(NOT fortnum_POPULATED)
    FetchContent_Populate(fortnum)
endif()
add_library(kin6d_intervals
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_rounding.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_interval.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_verified_midpoint_kernel.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_verified_quadrature.f90)
target_compile_options(kin6d_intervals PRIVATE -fno-fast-math -ffp-contract=off)
add_library(kin6d_ode
    ${fortnum_SOURCE_DIR}/src/fortnum_kinds.f90
    ${fortnum_SOURCE_DIR}/src/roots/fortnum_roots.f90
    ${fortnum_SOURCE_DIR}/src/quadrature/fortnum_quadrature.f90
    ${fortnum_SOURCE_DIR}/src/fortnum_status.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_rk54_ck_tableau.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_cash_karp.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_geometric.f90)
target_compile_options(kin6d_ode PRIVATE -fno-fast-math -ffp-contract=off)
# Generated derivative families retain the primal ABI even for constant partials.
# Unused primal arguments are intentional; keep other diagnostics enabled.
set_source_files_properties(src/generated/kin6d_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_library(kin6d src/generated/kin6d_leaves.f90 src/laboratory.f90
    src/generated/kin6d_species_leaves.f90 src/species.f90
    src/generated/kin6d_coulomb_leaves.f90 src/ewald.f90 src/nbody.f90
    src/generated/kin6d_poisson_leaves.f90 src/poisson_laboratory.f90
    src/generated/kin6d_fp_leaves.f90 src/fp_laboratory.f90
    src/generated/kin6d_transport_leaves.f90 src/transport_laboratory.f90
    src/transport_certificate.f90 src/transport_sampling.f90
    src/generated/kin6d_streaming_leaves.f90 src/streaming_comparison.f90
    src/streaming_particle_density.f90
    src/generated/kin6d_streaming_moments_leaves.f90 src/streaming_moments.f90
    src/generated/kin6d_debye_sphere_leaves.f90 src/debye_sphere_reference.f90
    src/generated/kin6d_slab_oscillation_leaves.f90 src/slab_oscillation.f90
    src/generated/kin6d_shielding_dynamics_leaves.f90 src/shielding_dynamics.f90
    src/generated/kin6d_interaction_leaves.f90 src/interaction_laboratory.f90
    src/generated/kin6d_cloud_leaves.f90 src/microscopic_cloud.f90)
set_source_files_properties(
    src/generated/kin6d_species_leaves.f90 src/generated/kin6d_coulomb_leaves.f90
    src/generated/kin6d_poisson_leaves.f90 src/generated/kin6d_fp_leaves.f90
    src/generated/kin6d_transport_leaves.f90 src/generated/kin6d_streaming_leaves.f90
    src/generated/kin6d_streaming_moments_leaves.f90
    src/generated/kin6d_debye_sphere_leaves.f90
    src/generated/kin6d_slab_oscillation_leaves.f90
    src/generated/kin6d_shielding_dynamics_leaves.f90
    src/generated/kin6d_interaction_leaves.f90
    src/generated/kin6d_cloud_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
foreach(quantum_family fermi_gas quantum_fock quantum_evolution quantum_correlations
        quantum_reference_gallery classical_pair)
    target_sources(kin6d PRIVATE
        src/generated/kin6d_${quantum_family}_leaves.f90 src/${quantum_family}.f90)
    set_source_files_properties(src/generated/kin6d_${quantum_family}_leaves.f90
        PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
endforeach()
target_link_libraries(kin6d PUBLIC kin6d_intervals kin6d_ode)
target_compile_options(kin6d PRIVATE -Wall -Wextra -fno-fast-math -ffp-contract=off)
add_executable(kin6d_demo app/kin6d_demo.f90)
target_link_libraries(kin6d_demo PRIVATE kin6d fortplot)
add_executable(steady_laboratory_demo app/steady_laboratory_demo.f90)
target_link_libraries(steady_laboratory_demo PRIVATE kin6d fortplot)
add_executable(transport_demo app/transport_demo.f90)
target_link_libraries(transport_demo PRIVATE kin6d fortplot)
add_executable(streaming_comparison_demo app/streaming_comparison_demo.f90)
target_link_libraries(streaming_comparison_demo PRIVATE kin6d fortplot)
foreach(example streaming_moments debye_sphere slab_oscillation shielding_dynamics interaction)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE kin6d fortplot)
endforeach()
add_executable(microscopic_cloud_demo app/microscopic_cloud_demo.f90)
target_link_libraries(microscopic_cloud_demo PRIVATE kin6d fortplot)
foreach(example quantum_reference classical_pair)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE kin6d fortplot)
endforeach()
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
add_executable(test_transport_linear_bands test/test_transport_linear_bands.f90)
target_link_libraries(test_transport_linear_bands PRIVATE kin6d)
add_test(NAME continuous_linear_band_oracles COMMAND test_transport_linear_bands)
add_executable(test_streaming_comparison test/test_streaming_comparison.f90)
target_link_libraries(test_streaming_comparison PRIVATE kin6d)
add_test(NAME streaming_representation_oracles COMMAND test_streaming_comparison)
add_executable(test_streaming_particle_density test/test_streaming_particle_density.f90)
target_link_libraries(test_streaming_particle_density PRIVATE kin6d)
add_test(NAME streaming_particle_density_oracles COMMAND test_streaming_particle_density)
foreach(example streaming_moments debye_sphere slab_oscillation shielding_dynamics interaction)
    add_executable(test_${example} test/test_${example}.f90)
    target_link_libraries(test_${example} PRIVATE kin6d)
    add_test(NAME ${example}_oracles COMMAND test_${example})
endforeach()
add_executable(test_nbody test/test_nbody.f90)
target_link_libraries(test_nbody PRIVATE kin6d)
add_test(NAME nbody_coulomb_species_oracles COMMAND test_nbody)
add_executable(test_cloud test/test_cloud.f90)
target_link_libraries(test_cloud PRIVATE kin6d)
add_test(NAME microscopic_cloud_oracles COMMAND test_cloud)
foreach(quantum_family fermi_gas quantum_fock quantum_evolution quantum_correlations
        quantum_reference_gallery classical_pair)
    add_executable(test_${quantum_family} test/test_${quantum_family}.f90)
    target_link_libraries(test_${quantum_family} PRIVATE kin6d)
    add_test(NAME ${quantum_family}_oracles COMMAND test_${quantum_family})
endforeach()
add_executable(test_quantum_cutoff test/test_quantum_cutoff.f90)
target_link_libraries(test_quantum_cutoff PRIVATE kin6d)
add_test(NAME quantum_cutoff_oracles COMMAND test_quantum_cutoff)
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
set(KIN6D_FORTSYM_REVISION "a6142712af998b727c245fb6cfbb3745b296ec24")
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
    configure_file(gen/kin6d_fortsym_revision.inc.in
        ${CMAKE_BINARY_DIR}/gen_include/kin6d_fortsym_revision.inc @ONLY)
    add_library(kin6d_gen_support STATIC gen/kin6d_gen_support.f90
        gen/kin6d_symbolic_physics.f90)
    target_include_directories(kin6d_gen_support PUBLIC
        ${KIN6D_FORTSYM_BUILD}/include ${CMAKE_BINARY_DIR}/gen_include)
    target_link_libraries(kin6d_gen_support PUBLIC
        ${KIN6D_FORTSYM_BUILD}/lib/libfortsym${CMAKE_SHARED_LIBRARY_SUFFIX})
    foreach(generator laboratory poisson fp transport streaming streaming_moments
            debye_sphere slab_oscillation shielding_dynamics interaction cloud
            fermi_gas quantum_fock quantum_evolution quantum_correlations
            quantum_reference_gallery classical_pair coulomb species)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE kin6d_gen_support)
    endforeach()
endif()
# Explicit performance audit; not part of the ordinary correctness gate.
add_executable(benchmark_streaming app/benchmark_streaming.f90)
target_link_libraries(benchmark_streaming PRIVATE kin6d)
add_executable(benchmark_streaming_comparison app/benchmark_streaming_comparison.f90)
target_link_libraries(benchmark_streaming_comparison PRIVATE kin6d)
add_executable(benchmark_fp app/benchmark_fp.f90)
target_link_libraries(benchmark_fp PRIVATE kin6d)
add_executable(benchmark_transport app/benchmark_transport.f90)
target_link_libraries(benchmark_transport PRIVATE kin6d)
add_executable(benchmark_slab_oscillation app/benchmark_slab_oscillation.f90)
target_link_libraries(benchmark_slab_oscillation PRIVATE kin6d)

add_executable(benchmark_fermi_gas app/benchmark_fermi_gas.f90)
target_link_libraries(benchmark_fermi_gas PRIVATE kin6d)
add_executable(benchmark_quantum_cutoff app/benchmark_quantum_cutoff.f90)
target_link_libraries(benchmark_quantum_cutoff PRIVATE kin6d)

# Q7: static screening hierarchy and 1X1V Wigner-Poisson (appended block).
add_library(kin6d_fft
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft8_kernel.f90
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft.f90)
target_link_libraries(kin6d_fft PUBLIC kin6d_ode)
target_compile_options(kin6d_fft PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(kin6d PRIVATE
    src/generated/kin6d_screening_leaves.f90 src/screening_hierarchy.f90
    src/generated/kin6d_wigner_leaves.f90 src/wigner_dispersion.f90
    src/wigner_poisson.f90)
set_source_files_properties(src/generated/kin6d_screening_leaves.f90
    src/generated/kin6d_wigner_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d PUBLIC kin6d_fft)
foreach(q7 screening wigner)
    add_executable(${q7}_demo app/${q7}_demo.f90)
    target_link_libraries(${q7}_demo PRIVATE kin6d fortplot)
    add_executable(test_${q7} test/test_${q7}.f90)
    target_link_libraries(test_${q7} PRIVATE kin6d)
    add_test(NAME ${q7}_q7_oracles COMMAND test_${q7})
    if(KIN6D_FORTSYM_BUILD)
        add_executable(gen_${q7} gen/gen_${q7}.f90)
        target_link_libraries(gen_${q7} PRIVATE kin6d_gen_support)
    endif()
endforeach()
