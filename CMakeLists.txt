cmake_minimum_required(VERSION 3.22)
project(kin6d LANGUAGES Fortran C)
include(CTest)
include(FetchContent)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
# Dependency pins: one source used by FetchContent and by result provenance.
set(KIN6D_FORTPLOT_REVISION "67834eb4568720226ca97a94e02267f48542613d")
set(KIN6D_FORTNUM_REVISION "80b73401c973d1451ee5e802775a7feda41981e1")
set(KIN6D_FORTIO_REVISION 38add03dc691a7cd4715cd3e611ad42422e46091)
set(FORTPLOT_BUILD_FORTARRAY_ADAPTER OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG ${KIN6D_FORTPLOT_REVISION})
FetchContent_MakeAvailable(fortplot)
FetchContent_Declare(fortio
    GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
    GIT_TAG ${KIN6D_FORTIO_REVISION})
FetchContent_MakeAvailable(fortio)
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ftrampoline-impl=heap" KIN6D_HAVE_HEAP_TRAMPOLINES)
if(KIN6D_HAVE_HEAP_TRAMPOLINES)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ftrampoline-impl=heap>)
endif()
# Compile FortNum's interval and ODE modules unchanged from the pinned tree.
# The application owns physical laws; generic numerics remain upstream.
FetchContent_Declare(fortnum
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    GIT_TAG ${KIN6D_FORTNUM_REVISION})
FetchContent_GetProperties(fortnum)
if(NOT fortnum_POPULATED)
    FetchContent_Populate(fortnum)
endif()
add_library(kin6d_intervals
    ${fortnum_SOURCE_DIR}/src/fortnum_status.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_rounding.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_interval.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_verified_midpoint_kernel.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_verified_quadrature.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_ode_residual_certificate.f90)
target_compile_options(kin6d_intervals PRIVATE -fno-fast-math -ffp-contract=off)
add_library(kin6d_ode
    ${fortnum_SOURCE_DIR}/src/fortnum_kinds.f90
    ${fortnum_SOURCE_DIR}/src/roots/fortnum_roots.f90
    ${fortnum_SOURCE_DIR}/src/quadrature/fortnum_quadrature.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_rk54_ck_tableau.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_cash_karp.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_geometric.f90
    ${fortnum_SOURCE_DIR}/src/rng/fortnum_rng.f90
    ${fortnum_SOURCE_DIR}/src/interp/fortnum_polynomial.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_legendre_recurrence_kernel.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_legendre.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_krylov.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_symmetric_eigen.f90)
target_compile_options(kin6d_ode PRIVATE -fno-fast-math -ffp-contract=off)
target_link_libraries(kin6d_ode PUBLIC kin6d_intervals)
# Generated derivative families retain the primal ABI even for constant partials.
# Unused primal arguments are intentional; keep other diagnostics enabled.
set_source_files_properties(src/generated/kin6d_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_library(kin6d src/generated/kin6d_leaves.f90 src/laboratory.f90 src/sort_keys.f90
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
        quantum_reference_gallery classical_pair jellium_cell)
    target_sources(kin6d PRIVATE
        src/generated/kin6d_${quantum_family}_leaves.f90 src/${quantum_family}.f90)
    set_source_files_properties(src/generated/kin6d_${quantum_family}_leaves.f90
        PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
endforeach()
# Q1/Q1b partial-wave encounter: banded symmetric eigensolver from LAPACK,
# the external numerical dependency FortNum already declares.
find_package(LAPACK REQUIRED)
target_sources(kin6d PRIVATE src/generated/kin6d_coulomb_encounter_leaves.f90
    src/radial_dvr.f90 src/coulomb_encounter.f90 src/coulomb_kepler.f90)
set_source_files_properties(src/generated/kin6d_coulomb_encounter_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d PUBLIC kin6d_intervals kin6d_ode LAPACK::LAPACK)
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
foreach(example quantum_reference classical_pair coulomb_encounter jellium_cell)
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
        quantum_reference_gallery classical_pair jellium_cell)
    add_executable(test_${quantum_family} test/test_${quantum_family}.f90)
    target_link_libraries(test_${quantum_family} PRIVATE kin6d)
    add_test(NAME ${quantum_family}_oracles COMMAND test_${quantum_family})
endforeach()
add_executable(test_coulomb_encounter test/test_coulomb_encounter.f90)
target_link_libraries(test_coulomb_encounter PRIVATE kin6d)
add_test(NAME coulomb_encounter_oracles COMMAND test_coulomb_encounter)
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
            quantum_reference_gallery classical_pair coulomb species coulomb_encounter
            jellium_cell)
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

# Campaign Q3/Q8: finite-N ideal fermion ensembles and relativistic Fermi EOS.
if(KIN6D_FORTSYM_BUILD)
    foreach(generator relativistic_fermi fermion_ensemble)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE kin6d_gen_support)
    endforeach()
endif()
target_sources(kin6d PRIVATE src/generated/kin6d_relativistic_fermi_leaves.f90)
set_source_files_properties(src/generated/kin6d_relativistic_fermi_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_relativistic_fermi_gas test/test_relativistic_fermi_gas.f90)
target_link_libraries(test_relativistic_fermi_gas PRIVATE kin6d)
add_test(NAME relativistic_fermi_gas_oracles COMMAND test_relativistic_fermi_gas)
add_executable(test_fermi_gas_cold test/test_fermi_gas_cold.f90)
target_link_libraries(test_fermi_gas_cold PRIVATE kin6d)
add_test(NAME fermi_gas_cold_sommerfeld COMMAND test_fermi_gas_cold)
target_sources(kin6d PRIVATE src/generated/kin6d_fermion_ensemble_leaves.f90
    src/fermion_ensemble.f90)
set_source_files_properties(src/generated/kin6d_fermion_ensemble_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_fermion_ensemble test/test_fermion_ensemble.f90)
target_link_libraries(test_fermion_ensemble PRIVATE kin6d)
add_test(NAME fermion_ensemble_oracles COMMAND test_fermion_ensemble)
add_executable(fermion_ensemble_demo app/fermion_ensemble_demo.f90)
target_link_libraries(fermion_ensemble_demo PRIVATE kin6d fortplot)
add_executable(relativistic_fermi_demo app/relativistic_fermi_demo.f90)
target_link_libraries(relativistic_fermi_demo PRIVATE kin6d fortplot)

# Campaign F5: namelist configuration in, NetCDF-4/HDF5 result with embedded
# literal/resolved configuration and provenance out (src/result_io.f90).
# Physics stays in kin6d; only kin6d_io and applications see FortIO.
set(kin6d_build_info_inc ${CMAKE_BINARY_DIR}/result_include/kin6d_build_info.inc)
string(TOUPPER "${CMAKE_BUILD_TYPE}" kin6d_build_type_upper)
add_custom_target(kin6d_build_info
    COMMAND ${CMAKE_COMMAND} -DKIN6D_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DKIN6D_INFO_OUTPUT=${kin6d_build_info_inc}
        -DKIN6D_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        "-DKIN6D_BUILD_FLAGS=${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${kin6d_build_type_upper}}"
        "-DKIN6D_BUILD_COMPILER=${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}"
        -DKIN6D_FORTSYM_REVISION=${KIN6D_FORTSYM_REVISION}
        -DKIN6D_FORTNUM_REVISION=${KIN6D_FORTNUM_REVISION}
        -DKIN6D_FORTPLOT_REVISION=${KIN6D_FORTPLOT_REVISION}
        -DKIN6D_FORTIO_REVISION=${KIN6D_FORTIO_REVISION}
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_build_info.cmake
    BYPRODUCTS ${kin6d_build_info_inc} VERBATIM)
add_library(kin6d_io src/result_io.f90 src/periodic_pair_config.f90)
add_dependencies(kin6d_io kin6d_build_info)
target_include_directories(kin6d_io PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(kin6d_io PUBLIC kin6d fortio::fortio)
target_compile_options(kin6d_io PRIVATE -Wall -Wextra)
# Every application reads --config FILE and writes output/<name>.nc.
foreach(application kin6d_demo steady_laboratory_demo transport_demo
        streaming_comparison_demo streaming_moments_demo debye_sphere_demo
        slab_oscillation_demo shielding_dynamics_demo interaction_demo
        microscopic_cloud_demo quantum_reference_demo classical_pair_demo
        fermion_ensemble_demo relativistic_fermi_demo formal_krook_certificate
        benchmark_streaming benchmark_streaming_comparison benchmark_fp
        benchmark_transport benchmark_slab_oscillation benchmark_fermi_gas
        benchmark_quantum_cutoff jellium_cell_demo coulomb_encounter_demo)
    target_link_libraries(${application} PRIVATE kin6d_io)
endforeach()
add_executable(test_result_io test/test_result_io.f90)
target_link_libraries(test_result_io PRIVATE kin6d_io)
add_test(NAME result_io_contract
    COMMAND ${CMAKE_COMMAND} -DPROBE=$<TARGET_FILE:test_result_io>
        -DCLASSICAL=$<TARGET_FILE:classical_pair_demo>
        -DQUANTUM=$<TARGET_FILE:quantum_reference_demo>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/result_io_check
        -DFORTIO_PIN=${KIN6D_FORTIO_REVISION}
        -P ${CMAKE_SOURCE_DIR}/test/result_io_check.cmake)
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
    target_link_libraries(${q7}_demo PRIVATE kin6d kin6d_io fortplot)
    add_executable(test_${q7} test/test_${q7}.f90)
    target_link_libraries(test_${q7} PRIVATE kin6d)
    add_test(NAME ${q7}_q7_oracles COMMAND test_${q7})
    if(KIN6D_FORTSYM_BUILD)
        add_executable(gen_${q7} gen/gen_${q7}.f90)
        target_link_libraries(gen_${q7} PRIVATE kin6d_gen_support)
    endif()
endforeach()
# Q5: Hartree/HF/Mermin/TDHF mean-field models of the Q4 jellium cell (appended).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_mean_field_cell gen/gen_mean_field_cell.f90)
    target_link_libraries(gen_mean_field_cell PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d PRIVATE src/generated/kin6d_mean_field_cell_leaves.f90
    src/mean_field_cell.f90)
set_source_files_properties(src/generated/kin6d_mean_field_cell_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_mean_field_cell test/test_mean_field_cell.f90)
target_link_libraries(test_mean_field_cell PRIVATE kin6d)
add_test(NAME mean_field_cell_oracles COMMAND test_mean_field_cell)
add_executable(mean_field_cell_demo app/mean_field_cell_demo.f90)
target_link_libraries(mean_field_cell_demo PRIVATE kin6d kin6d_io fortplot)
