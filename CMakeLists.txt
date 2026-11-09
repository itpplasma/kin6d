cmake_minimum_required(VERSION 3.22)
project(milca LANGUAGES Fortran C)

# gfortran < 13 mis-optimizes an allocatable left-hand side of `x = matmul(...)`
# at -O1 and above: the inlined matmul does not reallocate x and writes past
# the heap block (minimal reproducer: gfortran 12.2 -O2, valgrind). Calling the
# library matmul avoids it; newer compilers are not affected.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 13)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-finline-matmul-limit=0>)
endif()
include(CTest)
include(FetchContent)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
# OpenMP is on by default; MILCA_ENABLE_OPENMP=OFF gives a serial build. The
# thread count and determinism mode are written to the application provenance.
option(MILCA_ENABLE_OPENMP "Build with OpenMP threading" ON)
if(MILCA_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
# Dependency pins: one source used by FetchContent and by result provenance.
set(MILCA_FORTPLOT_REVISION "35a9cd26ab20cf50ba6f7347d7daf13566e0dc91")
set(MILCA_FORTNUM_REVISION "742d49a73517a9be61ebc9ca18fc0a4c63f2f08c")
set(MILCA_FORTIO_REVISION 38add03dc691a7cd4715cd3e611ad42422e46091)
set(FORTPLOT_BUILD_FORTARRAY_ADAPTER OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fortplot
    GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
    GIT_TAG ${MILCA_FORTPLOT_REVISION})
FetchContent_MakeAvailable(fortplot)
FetchContent_Declare(fortio
    GIT_REPOSITORY https://github.com/lazy-fortran/fortio.git
    GIT_TAG ${MILCA_FORTIO_REVISION})
FetchContent_MakeAvailable(fortio)
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ftrampoline-impl=heap" MILCA_HAVE_HEAP_TRAMPOLINES)
if(MILCA_HAVE_HEAP_TRAMPOLINES)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ftrampoline-impl=heap>)
endif()
# Use the complete pinned numerical target, including precision modules needed
# by FortSparse. The compatibility layer names preserve MILCA's layer graph;
# FortNum owns module compilation, strict arithmetic and OpenMP settings.
set(FORTNUM_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FORTNUM_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_Declare(fortnum
    GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
    GIT_TAG ${MILCA_FORTNUM_REVISION})
FetchContent_MakeAvailable(fortnum)
add_library(milca_intervals INTERFACE)
target_link_libraries(milca_intervals INTERFACE fortnum)
# The pinned FortNum target does not yet register its geometric ODE module.
# Retain this existing consumer leaf unchanged until its upstream CMake list
# includes it; every other numerical module comes from the complete target.
add_library(milca_ode
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_geometric.f90)
target_link_libraries(milca_ode PUBLIC fortnum)
target_compile_options(milca_ode PRIVATE -fno-fast-math -ffp-contract=off)
# MAKEGRID coils reader with circuit groups (itpplasma/libneo#469) and the ten-column
# Boozer .bc reader (itpplasma/libneo#509) from libneo main.
# Only its two dependency-free modules are compiled; libneo itself is not configured.
set(MILCA_LIBNEO_REVISION "d7c0cc08d3ae85467d81ee7d5c84f3ec42e4f745")
# SOURCE_SUBDIR names a directory without CMakeLists.txt, so MakeAvailable only
# populates; unlike a direct FetchContent_Populate it honours dependency
# providers and FETCHCONTENT_SOURCE_DIR_LIBNEO in disconnected builds.
FetchContent_Declare(libneo
    GIT_REPOSITORY https://github.com/itpplasma/libneo.git
    GIT_TAG ${MILCA_LIBNEO_REVISION}
    SOURCE_SUBDIR milca-populate-only)
FetchContent_MakeAvailable(libneo)
add_library(milca_libneo_coils
    ${libneo_SOURCE_DIR}/src/field/biotsavart.f90
    ${libneo_SOURCE_DIR}/src/field/makegrid_coils.f90
    ${libneo_SOURCE_DIR}/src/field/bc_file.f90)
# Generated derivative families retain the primal ABI even for constant partials.
# Unused primal arguments are intentional; keep other diagnostics enabled.
set_source_files_properties(src/generated/milca_leaves.f90
    src/generated/milca_vrmhd_leaves.f90
    src/generated/milca_vrmhd_section_leaves.f90
    src/generated/milca_pp_fusion_leaves.f90
    src/generated/milca_pp_fusion_eft_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
# Layered libraries (ROADMAP "Build structure and fast development loop"):
# L0 milca_core; L1 milca_fem2d, milca_pairs; L2 model classes; `milca` is the
# umbrella. Each source names its layer; the milca_layering test checks the use
# graph against this partition.
include(cmake/milca_layers.cmake)
target_sources(milca_mhd PRIVATE src/generated/milca_vrmhd_leaves.f90
    src/mhd/vrmhd_parameters.f90
    src/generated/milca_vrmhd_leaf_ids.f90
    src/generated/milca_vrmhd_dispatch.f90
    src/generated/milca_vrmhd_section_leaves.f90
    src/mhd/vrmhd_point.f90 src/mhd/vrmhd_section.f90
    src/mhd/vrmhd_boundary.f90 src/mhd/vrmhd_stationary.f90
    src/mhd/vrmhd_drive.f90 src/mhd/vrmhd_evolve.f90 src/mhd/vrmhd_fourier.f90
    src/mhd/vrmhd_fourier_nonlinear.f90)
add_library(milca_vrmhd_wall_test STATIC test/support/vrmhd_wall_oracle.f90)
target_link_libraries(milca_vrmhd_wall_test PUBLIC milca_mhd)
add_library(milca_vrmhd_mms_test STATIC test/support/vrmhd_mms_leaves.f90)
target_link_libraries(milca_vrmhd_mms_test PUBLIC milca_mhd)
add_executable(test_vrmhd_mms test/test_vrmhd_mms.f90)
target_link_libraries(test_vrmhd_mms PRIVATE milca_vrmhd_mms_test)
add_executable(test_vrmhd_stationary test/test_vrmhd_stationary.f90)
target_link_libraries(test_vrmhd_stationary PRIVATE milca_vrmhd_mms_test milca_vrmhd_wall_test
    fortnum)
add_executable(test_vrmhd_gs test/test_vrmhd_gs.f90
    src/generated/milca_vrmhd_gs_leaves.f90)
target_link_libraries(test_vrmhd_gs PRIVATE milca_vrmhd_wall_test)
add_test(NAME vrmhd_gs_consistency COMMAND test_vrmhd_gs)
set_tests_properties(vrmhd_gs_consistency PROPERTIES LABELS rm2 PROCESSORS 8
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
foreach(gate boundary derivatives harmonic gauge solve)
    add_test(NAME vrmhd_stationary_${gate} COMMAND test_vrmhd_stationary ${gate})
    set_tests_properties(vrmhd_stationary_${gate} PROPERTIES LABELS rm2
        PROCESSORS 8 ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
endforeach()
set_tests_properties(vrmhd_stationary_gauge PROPERTIES
    TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=8")
add_test(NAME vrmhd_section_manufactured_consistency COMMAND test_vrmhd_mms)
set_tests_properties(vrmhd_section_manufactured_consistency PROPERTIES LABELS rm2
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_library(milca_vrmhd_slab_test STATIC test/support/vrmhd_periodic_slab.f90)
target_link_libraries(milca_vrmhd_slab_test PUBLIC milca_mhd)
add_executable(test_vrmhd_evolve test/test_vrmhd_evolve.f90)
target_link_libraries(test_vrmhd_evolve PRIVATE milca_vrmhd_slab_test
    milca_spectrum_support)
foreach(gate energy decay waves steady)
    add_test(NAME vrmhd_evolve_${gate} COMMAND test_vrmhd_evolve ${gate})
    set_tests_properties(vrmhd_evolve_${gate} PROPERTIES LABELS rm2-evolve
        PROCESSORS 8 TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_vrmhd_island test/test_vrmhd_island.f90
    test/support/vrmhd_island_leaves.f90)
target_link_libraries(test_vrmhd_island PRIVATE milca_vrmhd_slab_test)
foreach(family equal mixed)
    foreach(gate boundary equilibrium width ledger graded)
        add_test(NAME vrmhd_island_${gate}_${family}
            COMMAND test_vrmhd_island ${gate} ${family})
        set_tests_properties(vrmhd_island_${gate}_${family} PROPERTIES LABELS rm-b5
            PROCESSORS 4 TIMEOUT 300
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()
add_test(NAME vrmhd_island_induction_mixed
    COMMAND test_vrmhd_island induction mixed 8 8)
set_tests_properties(vrmhd_island_induction_mixed PROPERTIES LABELS rm-b5
    PROCESSORS 4 TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_island_preconditioner COMMAND test_vrmhd_island preconditioner mixed)
set_tests_properties(vrmhd_island_preconditioner PROPERTIES LABELS rm-b5
    PROCESSORS 4 TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_island_decimal_steps
    COMMAND test_vrmhd_island probe mixed 4 4 0.01 0.1 8 0.01)
set_tests_properties(vrmhd_island_decimal_steps PROPERTIES LABELS rm-b5
    PROCESSORS 4 TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_island_graded_input COMMAND test_vrmhd_island
    width mixed 4 6 0.01 0.125 4 0.01 supernodal
    ${CMAKE_SOURCE_DIR}/test/fixtures/vrmhd_island/radial_mesh.nml)
set_tests_properties(vrmhd_island_graded_input PROPERTIES LABELS rm-b5
    PROCESSORS 4 TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_vrmhd_section test/test_vrmhd_section.f90)
target_link_libraries(test_vrmhd_section PRIVATE milca_mhd)
add_executable(test_vrmhd_energy test/test_vrmhd_energy.f90)
target_link_libraries(test_vrmhd_energy PRIVATE milca_vrmhd_wall_test)
add_test(NAME vrmhd_section_energy COMMAND test_vrmhd_energy)
set_tests_properties(vrmhd_section_energy PROPERTIES LABELS rm2
    TIMEOUT 60 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_section_mixed_derivatives COMMAND test_vrmhd_section mixed)
set_tests_properties(vrmhd_section_mixed_derivatives PROPERTIES LABELS rm2-mixed
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
foreach(degree 1 2)
    add_test(NAME vrmhd_section_mixed_energy_p${degree}
        COMMAND test_vrmhd_energy mixed${degree})
    add_test(NAME vrmhd_section_mixed_mms_p${degree}
        COMMAND test_vrmhd_mms mixed${degree})
    set_tests_properties(vrmhd_section_mixed_energy_p${degree}
        vrmhd_section_mixed_mms_p${degree} PROPERTIES LABELS rm2-mixed
        TIMEOUT 900 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(gate geometry derivatives volume closure)
    add_test(NAME vrmhd_section_${gate} COMMAND test_vrmhd_section ${gate})
    set_tests_properties(vrmhd_section_${gate} PROPERTIES LABELS rm2
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
add_library(milca_vrmhd_mms_fourier_test STATIC
    test/support/vrmhd_mms_fourier_leaves.f90
    test/support/vrmhd_mms_nonlinear_leaves.f90)
target_link_libraries(milca_vrmhd_mms_fourier_test PUBLIC milca_mhd)
add_executable(test_vrmhd_fourier_nonlinear test/test_vrmhd_fourier_nonlinear.f90)
target_link_libraries(test_vrmhd_fourier_nonlinear PRIVATE milca_mhd)
foreach(gate coupling derivatives energy)
    add_test(NAME vrmhd_fourier_nonlinear_${gate}
        COMMAND test_vrmhd_fourier_nonlinear ${gate})
    set_tests_properties(vrmhd_fourier_nonlinear_${gate} PROPERTIES LABELS rm3b
        PROCESSORS 4 TIMEOUT 900
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_vrmhd_harmonics test/test_vrmhd_harmonics.f90)
target_link_libraries(test_vrmhd_harmonics PRIVATE milca_vrmhd_mms_fourier_test
    milca_spectrum_support)
foreach(gate exactness reduction decoupling manufactured spectra accumulation)
    add_test(NAME vrmhd_harmonics_${gate} COMMAND test_vrmhd_harmonics ${gate})
    set_tests_properties(vrmhd_harmonics_${gate} PROPERTIES LABELS rm3
        PROCESSORS 8 ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(gate exactness reduction decoupling manufactured)
    add_test(NAME vrmhd_harmonics_mixed_${gate}
        COMMAND test_vrmhd_harmonics mixed_${gate})
    set_tests_properties(vrmhd_harmonics_mixed_${gate} PROPERTIES LABELS rm3-mixed
        PROCESSORS 4 TIMEOUT 900
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_test(NAME vrmhd_harmonics_spectra_fine COMMAND test_vrmhd_harmonics spectra_fine)
set_tests_properties(vrmhd_harmonics_spectra_fine PROPERTIES PROCESSORS 8
    TIMEOUT 3600 ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_harmonics_helical COMMAND test_vrmhd_harmonics helical)
set_tests_properties(vrmhd_harmonics_helical PROPERTIES PROCESSORS 4 TIMEOUT 900
    LABELS rm3b ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
foreach(family equal mixed)
    add_test(NAME vrmhd_fourier_nonlinear_manufactured_${family}
        COMMAND test_vrmhd_harmonics nonlinear_${family})
    set_tests_properties(vrmhd_fourier_nonlinear_manufactured_${family} PROPERTIES
        PROCESSORS 4 TIMEOUT 900 LABELS rm3b
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(family equal mixed)
    add_test(NAME vrmhd_helical_nonlinear_manufactured_${family}
        COMMAND test_vrmhd_harmonics helical_nonlinear_${family})
    set_tests_properties(vrmhd_helical_nonlinear_manufactured_${family} PROPERTIES
        PROCESSORS 4 TIMEOUT 900 LABELS rm3b
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_vrmhd_leaves test/test_vrmhd_leaves.f90
    test/support/vrmhd_leaf_dispatch.f90)
target_link_libraries(test_vrmhd_leaves PRIVATE milca_mhd)
foreach(gate derivatives physics energy covariance)
    add_test(NAME vrmhd_${gate} COMMAND test_vrmhd_leaves ${gate})
    set_tests_properties(vrmhd_${gate} PROPERTIES LABELS rm1)
endforeach()
target_sources(milca_core PRIVATE src/generated/milca_leaves.f90 src/core/laboratory.f90
    src/core/sort_keys.f90 src/core/real_identity.f90
    src/generated/milca_species_leaves.f90 src/core/species.f90)
target_sources(milca_pairs PRIVATE src/generated/milca_coulomb_leaves.f90 src/pairs/ewald.f90
    src/generated/milca_ocp_leaves.f90 src/pairs/near_encounter.f90 src/pairs/pme.f90
    src/pairs/nbody.f90 src/generated/milca_interaction_leaves.f90
    src/pairs/interaction_laboratory.f90)
target_sources(milca_classical PRIVATE src/classical/ocp.f90 src/classical/pair_correlation.f90
    src/generated/milca_cloud_leaves.f90 src/classical/microscopic_cloud.f90)
target_sources(milca_continuum PRIVATE src/continuum/vlasov_response.f90
    src/generated/milca_poisson_leaves.f90 src/continuum/poisson_laboratory.f90
    src/generated/milca_fp_leaves.f90 src/continuum/fp_laboratory.f90
    src/generated/milca_fp_relaxation_leaves.f90 src/continuum/fp_relaxation.f90
    src/generated/milca_dk_pitch_leaves.f90 src/continuum/dk_pitch.f90
    src/generated/milca_dk_mono_leaves.f90 src/continuum/dk_mono.f90
    src/generated/milca_dk_hp_leaves.f90 src/continuum/dk_hp.f90
    src/generated/milca_dk_ntv_leaves.f90 src/continuum/dk_ntv.f90
    src/generated/milca_dk_collop_leaves.f90 src/continuum/dk_collop.f90
    src/generated/milca_thermal_screening_leaves.f90 src/continuum/screening_moments.f90
    src/generated/milca_transport_leaves.f90 src/continuum/transport_laboratory.f90
    src/continuum/transport_certificate.f90 src/continuum/transport_sampling.f90
    src/generated/milca_streaming_leaves.f90 src/continuum/streaming_comparison.f90
    src/continuum/streaming_particle_density.f90
    src/generated/milca_streaming_moments_leaves.f90 src/continuum/streaming_moments.f90
    src/generated/milca_debye_sphere_leaves.f90 src/continuum/debye_sphere_reference.f90
    src/generated/milca_slab_oscillation_leaves.f90 src/continuum/slab_oscillation.f90
    src/generated/milca_shielding_dynamics_leaves.f90 src/continuum/shielding_dynamics.f90
    src/generated/milca_homma_leaves.f90 src/continuum/rosenbluth_reference.f90)
set_source_files_properties(
    src/generated/milca_species_leaves.f90 src/generated/milca_coulomb_leaves.f90
    src/generated/milca_ocp_leaves.f90
    src/generated/milca_poisson_leaves.f90 src/generated/milca_fp_leaves.f90
    src/generated/milca_transport_leaves.f90 src/generated/milca_streaming_leaves.f90
    src/generated/milca_streaming_moments_leaves.f90
    src/generated/milca_debye_sphere_leaves.f90
    src/generated/milca_slab_oscillation_leaves.f90
    src/generated/milca_shielding_dynamics_leaves.f90
    src/generated/milca_interaction_leaves.f90
    src/generated/milca_cloud_leaves.f90 src/generated/milca_homma_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
foreach(quantum_family fermi_gas quantum_fock quantum_evolution quantum_correlations
        quantum_reference_gallery classical_pair jellium_cell)
    set(family_layer milca_quantum)
    if(quantum_family STREQUAL "classical_pair")
        set(family_layer milca_classical)
    endif()
    string(REPLACE "milca_" "" family_directory "${family_layer}")
    target_sources(${family_layer} PRIVATE
        src/generated/milca_${quantum_family}_leaves.f90
        src/${family_directory}/${quantum_family}.f90)
    set_source_files_properties(src/generated/milca_${quantum_family}_leaves.f90
        PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
endforeach()
# Q1/Q1b partial-wave encounter: banded symmetric eigensolver from LAPACK,
# the external numerical dependency FortNum already declares.
find_package(LAPACK REQUIRED)
target_sources(milca_pairs PRIVATE src/generated/milca_coulomb_encounter_leaves.f90
    src/pairs/radial_dvr.f90)
target_sources(milca_quantum PRIVATE src/generated/milca_pp_fusion_leaves.f90
    src/quantum/pp_fusion.f90
    src/generated/milca_pp_fusion_eft_leaves.f90 src/quantum/pp_fusion_eft.f90
    src/quantum/coulomb_encounter.f90 src/quantum/coulomb_kepler.f90)
set_source_files_properties(src/generated/milca_coulomb_encounter_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(milca_core PUBLIC LAPACK::LAPACK)
option(MILCA_BUILD_GLL_SPIKE "Build the axisymmetric matrix-free GLL map experiment" OFF)
# FortFEM main (block_csc assembly, Hiptmair-Xu, circulant BEM with per-harmonic
# Steklov form #144, exact FEEC tables, split stretched BEM fans #145,
# second-kind Nedelec edges #148, exact moment gradients #149 and the
# boundary trace projection #152 and the toroidal period #154).
set(MILCA_FORTFEM_REVISION "cb7dfd809221a85cd33fcfffc87c5069a66d885f")
FetchContent_Declare(fortfem
    GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
    GIT_TAG ${MILCA_FORTFEM_REVISION})
FetchContent_GetProperties(fortfem)
if(NOT fortfem_POPULATED)
    FetchContent_Populate(fortfem)
endif()
# Narrow, unchanged upstream geometry/meshing slice; no private FEM geometry.
file(GLOB MILCA_TRIANGLE_SOURCES CONFIGURE_DEPENDS
    ${fortfem_SOURCE_DIR}/src/triangle_compat/*.f90)
add_library(milca_spatial ${MILCA_TRIANGLE_SOURCES}
    ${fortfem_SOURCE_DIR}/src/utils/kinds.f90
    ${fortfem_SOURCE_DIR}/src/assembly/element_to_csr.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/delaunay_types.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/robust_predicates.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/geometric_predicates.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/point_location.f90
    ${fortfem_SOURCE_DIR}/src/basis/basis_p1_2d.f90
    ${fortfem_SOURCE_DIR}/src/basis/basis_p2_2d.f90
    ${fortfem_SOURCE_DIR}/src/geometry/geometry_jet_2d.f90
    ${fortfem_SOURCE_DIR}/src/elements/tetra_piola_maps.f90
    ${fortfem_SOURCE_DIR}/src/quadrature/gauss_quadrature_2d.f90
    ${fortfem_SOURCE_DIR}/src/quadrature/triangle_duffy_quadrature.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_reference_basis_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_lagrange_product_jets.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_reference_scalar_intervals.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_triangle_level_geometry.f90
    ${fortfem_SOURCE_DIR}/src/quadrature/triangle_level_quadrature.f90
    ${fortfem_SOURCE_DIR}/src/bem/axisymmetric_gs_kernel.f90
    ${fortfem_SOURCE_DIR}/src/bem/axisymmetric_gs_boundary_operators_2d.f90
    ${fortfem_SOURCE_DIR}/src/estimators/bem_residual_energy.f90)
target_link_libraries(milca_spatial PUBLIC milca_ode milca_intervals)
# The upstream BEM assembly (panel-pair colouring) is threaded with OpenMP.
if(MILCA_ENABLE_OPENMP)
    target_link_libraries(milca_spatial PUBLIC OpenMP::OpenMP_Fortran)
endif()
set_source_files_properties(
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_reference_scalar_intervals.f90
    PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
target_sources(milca_fem2d PRIVATE src/fem2d/gs_load_batch.f90
    src/fem2d/gs_timers.f90 src/fem2d/gs_options.f90 src/fem2d/gs_csr.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_symmetry.f90
    src/generated/milca_gs_symmetry_leaves.f90
    src/generated/milca_grad_shafranov_leaves.f90
    src/generated/milca_grad_shafranov_batch.f90 src/equilibrium/grad_shafranov.f90)
if(MILCA_ENABLE_OPENMP)
    target_link_libraries(milca_core PUBLIC OpenMP::OpenMP_Fortran)
endif()
add_library(milca_gs_reference_test STATIC
    test/generated/milca_gs_reference_error_leaves.f90
    test/support/gs_reference_error.f90)
target_link_libraries(milca_gs_reference_test PUBLIC milca)
target_compile_options(milca_gs_reference_test PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(milca_equilibrium PRIVATE src/generated/milca_gs_coarea_leaves.f90 src/equilibrium/gs_coarea.f90)
target_sources(milca_equilibrium PRIVATE src/generated/milca_boozer_file_leaves.f90
    src/equilibrium/boozer_file.f90 src/generated/milca_boozer_orbit_leaves.f90
    src/equilibrium/boozer_orbit.f90)
target_sources(milca_equilibrium PRIVATE src/generated/milca_gs_inverse_leaves.f90
    src/equilibrium/gs_inverse.f90 src/equilibrium/gs_prescribed_q.f90)
set_source_files_properties(src/generated/milca_gs_inverse_leaves.f90
    src/generated/milca_gs_coarea_leaves.f90
    src/generated/milca_boozer_file_leaves.f90
    src/generated/milca_dk_pitch_leaves.f90
    src/generated/milca_dk_mono_leaves.f90
    src/generated/milca_dk_hp_leaves.f90
    src/generated/milca_dk_ntv_leaves.f90
    src/generated/milca_dk_collop_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
set_source_files_properties(src/generated/milca_grad_shafranov_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(milca_fem2d PUBLIC milca_spatial)
add_executable(test_real_identity test/test_real_identity.f90)
target_link_libraries(test_real_identity PRIVATE milca)
add_test(NAME real_value_identity_edge_cases COMMAND test_real_identity)
add_executable(test_gs_mirror_mesh test/test_gs_mirror_mesh.f90)
target_link_libraries(test_gs_mirror_mesh PRIVATE milca)
add_test(NAME gs_mirror_mesh COMMAND test_gs_mirror_mesh)
set_tests_properties(gs_mirror_mesh PROPERTIES TIMEOUT 60 PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(test_grad_shafranov test/test_grad_shafranov.f90)
target_link_libraries(test_grad_shafranov PRIVATE milca)
add_test(NAME grad_shafranov_analytic_and_derivative_oracles
    COMMAND test_grad_shafranov)
add_executable(test_gs_profile_quadrature test/test_gs_profile_quadrature.f90)
target_link_libraries(test_gs_profile_quadrature PRIVATE milca)
add_test(NAME grad_shafranov_profile_cut_quadrature_oracles
    COMMAND test_gs_profile_quadrature)
add_executable(test_gs_profile_domain test/test_gs_profile_domain.f90)
target_link_libraries(test_gs_profile_domain PRIVATE milca)
foreach(domain valid zero negative nan)
    add_test(NAME grad_shafranov_profile_mu0_${domain}
        COMMAND test_gs_profile_domain ${domain})
    add_test(NAME grad_shafranov_p2_profile_mu0_${domain}
        COMMAND test_gs_profile_domain ${domain} 2)
endforeach()
# Production sparse factorizations use project-owned FortSparse CPU factors.
# The opt-in fixed-source GS prototype requests its native FP64 SPD CUDA path.
option(MILCA_ENABLE_CUDA_GS_SOLVE "Build experimental host-staged CUDA GS solve" OFF)
if(MILCA_ENABLE_CUDA_GS_SOLVE)
    if(NOT CMAKE_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "CUDA GS prototype requires explicit CMAKE_CUDA_ARCHITECTURES")
    endif()
    enable_language(CUDA)
    set(FORTSPARSE_ENABLE_NATIVE_CUDA ON CACHE BOOL "" FORCE)
endif()
set(MILCA_FORTSPARSE_REVISION "f5ab3e0884123794521d4202f515ad80002fab9a")
set(FORTSPARSE_ENABLE_SUPERLU OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_ENABLE_UMFPACK_IPC OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_ENABLE_MUMPS OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FORTSPARSE_INSTALL_HELPER OFF CACHE BOOL "" FORCE)
if(NOT TARGET fortsparse)
    FetchContent_Declare(fortsparse
        GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
        GIT_TAG ${MILCA_FORTSPARSE_REVISION}
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortsparse)
endif()
get_target_property(_milca_sparse_definitions fortsparse COMPILE_DEFINITIONS)
if(NOT "FORTSPARSE_HAVE_AMD" IN_LIST _milca_sparse_definitions)
    message(FATAL_ERROR "MILCA native sparse solves require AMD graph ordering. "
        "Install AMD and configure FortSparse with FORTSPARSE_ENABLE_AMD=ON.")
endif()
foreach(_external FORTSPARSE_HAVE_SUPERLU FORTSPARSE_HAVE_MUMPS)
    if(_external IN_LIST _milca_sparse_definitions)
        message(FATAL_ERROR "MILCA native-only build received an external "
            "sparse backend: ${_external}")
    endif()
endforeach()
if(MILCA_ENABLE_OPENMP)
    target_link_libraries(fortsparse PRIVATE OpenMP::OpenMP_Fortran)
endif()
add_library(milca_sparse_direct
    ${fortfem_SOURCE_DIR}/src/topology/tree_cotree_gauge.f90
    ${fortfem_SOURCE_DIR}/src/solvers/fortfem_sparse_direct.f90)
target_link_libraries(milca_sparse_direct PUBLIC milca_spatial fortsparse)
target_link_libraries(milca_fem2d PUBLIC milca_sparse_direct)
# The local-theta NTV kernel factors its (theta, l) system with FortSparse.
target_link_libraries(milca_continuum PUBLIC fortsparse)
target_sources(milca_fem2d PRIVATE src/fem2d/gs_direct.f90)
message(STATUS "MILCA sparse direct backend: native FortSparse CPU")
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_coils.f90 src/equilibrium/gs_coupling.f90
    src/equilibrium/gs_free_boundary.f90 src/equilibrium/gs_coil_inverse.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/equilibrium_harness.f90)
add_executable(test_equilibrium_harness test/test_equilibrium_harness.f90)
target_link_libraries(test_equilibrium_harness PRIVATE milca)
add_test(NAME equilibrium_d8_registry_contract COMMAND test_equilibrium_harness)
# N1: nonlinear 3D frozen-in energy map (fixed boundary), DESIGN 3D lane.
target_sources(milca_mhd PRIVATE src/generated/milca_ideal_map_leaves.f90
    src/mhd/ideal_map.f90 src/mhd/map_3d_equilibrium.f90
    src/mhd/equilibrium_3d_driver.f90
    src/generated/milca_map_3d_readback_leaves.f90 src/mhd/ideal_seed.f90
    src/generated/milca_ideal_field_angle_leaves.f90
    src/generated/milca_ideal_seed_3d_leaves.f90 src/mhd/ideal_seed_3d.f90
    src/mhd/nested_seed_input.f90)
set_source_files_properties(src/generated/milca_ideal_map_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
# N4: 3D free boundary of the map lane (Gamma_c shell, per-n exterior).
target_sources(milca_mhd PRIVATE src/generated/milca_vacuum_shell_leaves.f90
    src/generated/milca_fold_barrier_leaves.f90 src/mhd/ideal_free_boundary.f90)
set_source_files_properties(src/generated/milca_vacuum_shell_leaves.f90
    src/generated/milca_fold_barrier_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_sources(milca_equilibrium PRIVATE
    src/generated/milca_nested_oracle_leaves.f90 src/equilibrium/nested_oracles.f90)
set_source_files_properties(src/generated/milca_nested_oracle_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
if(MILCA_ENABLE_OPENMP)
    add_executable(benchmark_nested_oracles test/benchmark_nested_oracles.f90)
    target_link_libraries(benchmark_nested_oracles PRIVATE milca OpenMP::OpenMP_Fortran)
    add_test(NAME equilibrium_nested_oracle_batch COMMAND benchmark_nested_oracles 1024)
    set_tests_properties(equilibrium_nested_oracle_batch PROPERTIES PROCESSORS 8)
endif()
add_executable(benchmark_nested_reference EXCLUDE_FROM_ALL benchmarks/nested/compare.f90)
target_link_libraries(benchmark_nested_reference PRIVATE milca_equilibrium)
add_executable(test_nested_boundary_input test/test_nested_boundary_input.f90)
target_link_libraries(test_nested_boundary_input PRIVATE milca)
foreach(family 1 2)
    if(family EQUAL 1)
        set(nested_config issan-integer)
    else()
        set(nested_config issan-sheared)
    endif()
    add_test(NAME equilibrium_nested_boundary_input_${family} COMMAND ${CMAKE_COMMAND}
        -DAPP=$<TARGET_FILE:nested_oracle_export>
        -DCHECK=$<TARGET_FILE:test_nested_boundary_input> -DFAMILY=${family}
        -DCONFIG=${CMAKE_SOURCE_DIR}/examples/nested-oracles/${nested_config}.nml
        -DWORK=${CMAKE_BINARY_DIR}/nested-export-${family}
        -P ${CMAKE_SOURCE_DIR}/test/check_nested_oracle_export.cmake)
    set_tests_properties(equilibrium_nested_boundary_input_${family} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
endforeach()
add_executable(test_nested_oracles test/test_nested_oracles.f90)
target_link_libraries(test_nested_oracles PRIVATE milca)
add_test(NAME equilibrium_nested_exact_oracles COMMAND test_nested_oracles)
set_tests_properties(equilibrium_nested_exact_oracles PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=4" PROCESSORS 4)
target_sources(milca_mhd PRIVATE src/generated/milca_landreman_leaves.f90
    src/mhd/landreman_benchmark.f90)
set_source_files_properties(src/generated/milca_landreman_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_landreman_free_boundary test/test_landreman_free_boundary.f90)
target_link_libraries(test_landreman_free_boundary PRIVATE milca)
add_test(NAME ideal_fb_landreman_axisymmetric COMMAND test_landreman_free_boundary)
set_tests_properties(ideal_fb_landreman_axisymmetric PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
# E3.3: the same cases through the full-surface R2 BEM exterior (24 rotated
# sectors, block-circulant assembly, per-harmonic Steklov form), compared
# with the per-n exterior. bem3d16 is the chain=16 level only; the two-level
# gate with its refinement checks runs nightly (ci/nightly-jobs.json).
add_test(NAME ideal_fb_landreman_bem3d16_axisymmetric
    COMMAND test_landreman_free_boundary bem3d16)
add_test(NAME ideal_fb_landreman_bem3d_axisymmetric
    COMMAND test_landreman_free_boundary bem3d)
set_tests_properties(ideal_fb_landreman_bem3d16_axisymmetric PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_tests_properties(ideal_fb_landreman_bem3d_axisymmetric PROPERTIES
    TIMEOUT 1800 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_ideal_free_boundary test/test_ideal_free_boundary.f90
    test/generated/ideal_fb_kink_leaves.f90)
target_link_libraries(test_ideal_free_boundary PRIVATE milca)
add_test(NAME ideal_fb_schur_spectrum COMMAND test_ideal_free_boundary spectrum)
set_tests_properties(ideal_fb_schur_spectrum PROPERTIES ENVIRONMENT
    "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4 TIMEOUT 120)
add_test(NAME ideal_fb_external_kink COMMAND test_ideal_free_boundary kink)
set_tests_properties(ideal_fb_external_kink PROPERTIES ENVIRONMENT
    "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4 TIMEOUT 120)
add_test(NAME ideal_fb_derivatives COMMAND test_ideal_free_boundary derivatives)
add_test(NAME ideal_fb_vertical_newton COMMAND test_ideal_free_boundary vertical)
foreach(k 2 3)
    add_test(NAME ideal_fb_vacuum_exact_n0_p${k} COMMAND test_ideal_free_boundary vacuum ${k})
endforeach()
add_test(NAME ideal_fb_vacuum_gamma_c_independence_n3
    COMMAND test_ideal_free_boundary vacuum3d 3)
add_test(NAME ideal_fb_slide_gauge_curl_free COMMAND test_ideal_free_boundary slide)
add_test(NAME ideal_fb_hypotheses_bordering COMMAND test_ideal_free_boundary hypotheses)
add_test(NAME ideal_fb_negative_curvature COMMAND test_ideal_free_boundary curvature)
set_tests_properties(ideal_fb_derivatives ideal_fb_vertical_newton
    ideal_fb_slide_gauge_curl_free ideal_fb_hypotheses_bordering
    ideal_fb_negative_curvature
    ideal_fb_vacuum_exact_n0_p2
    ideal_fb_vacuum_exact_n0_p3 ideal_fb_vacuum_gamma_c_independence_n3 PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4")
# One compiled ideal_map_cases module shared by the ideal map test programs.
add_library(milca_ideal_map_cases_test STATIC test/support/ideal_map_cases.f90)
target_link_libraries(milca_ideal_map_cases_test PUBLIC milca)
add_executable(test_ideal_map_axisymmetric test/test_ideal_map_axisymmetric.f90)
target_link_libraries(test_ideal_map_axisymmetric PRIVATE milca_ideal_map_cases_test)
# Parallel-only Levenberg-Marquardt term (docs/ideal-map-3d.md).
add_executable(test_ideal_map_parallel test/test_ideal_map_parallel.f90)
target_link_libraries(test_ideal_map_parallel PRIVATE milca_ideal_map_cases_test)
add_test(NAME ideal_map_parallel_mass COMMAND test_ideal_map_parallel unit 32 2)
add_test(NAME ideal_map_parallel_damping_p2
    COMMAND test_ideal_map_parallel compare 32 2)
set_tests_properties(ideal_map_parallel_mass ideal_map_parallel_damping_p2
    PROPERTIES TIMEOUT 600 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_ideal_map_hessian test/test_ideal_map_hessian.f90
    test/support/ideal_map_hessian_oracle.f90)
target_link_libraries(test_ideal_map_hessian PRIVATE milca)
add_test(NAME ideal_map_hessian_sum_factorisation_oracle COMMAND test_ideal_map_hessian)
set_tests_properties(ideal_map_hessian_sum_factorisation_oracle PROPERTIES TIMEOUT 120
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
foreach(degree RANGE 4 6)
    add_test(NAME ideal_map_hessian_p${degree}_oracle
        COMMAND test_ideal_map_hessian ${degree})
    set_tests_properties(ideal_map_hessian_p${degree}_oracle PROPERTIES TIMEOUT 120
        PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_equilibrium_3d test/test_equilibrium_3d.f90)
target_link_libraries(test_equilibrium_3d PRIVATE milca)
# Arguments: TAU MPOL NTOR NODES DEGREE and B/p/label gates at about twice the
# observed sampled errors recorded in docs/ideal-seed-3d.md; the P3 case also
# solves the coarse P2 level and requires every error to decrease.
add_test(NAME equilibrium_d8_landreman_3d_integer
    COMMAND test_equilibrium_3d 0 2 2 8 2 0.075 0.19 0.35)
add_test(NAME equilibrium_d8_landreman_3d_p3
    COMMAND test_equilibrium_3d 0 4 4 16 3 2.3e-3 2.2e-2 3.3e-2 2 2 8 2)
set_tests_properties(equilibrium_d8_landreman_3d_p3 PROPERTIES PROCESSORS 4 TIMEOUT 180)
add_test(NAME equilibrium_d8_issan_3d_integer
    COMMAND test_equilibrium_3d 0.4 2 2 8 2 0.055 0.14 0.26)
set_tests_properties(equilibrium_d8_landreman_3d_integer equilibrium_d8_issan_3d_integer
    PROPERTIES PROCESSORS 4 TIMEOUT 60)
add_executable(test_map_3d_equilibrium test/test_map_3d_equilibrium.f90)
target_link_libraries(test_map_3d_equilibrium PRIVATE milca)
add_test(NAME equilibrium_map_3d_readback COMMAND test_map_3d_equilibrium)
set_tests_properties(equilibrium_map_3d_readback PROPERTIES PROCESSORS 4 TIMEOUT 60)
add_executable(test_ideal_seed_3d test/test_ideal_seed_3d.f90)
target_link_libraries(test_ideal_seed_3d PRIVATE milca)
add_test(NAME ideal_seed_3d_boundary_profiles COMMAND test_ideal_seed_3d)
set_tests_properties(ideal_seed_3d_boundary_profiles PROPERTIES PROCESSORS 4 TIMEOUT 60)
add_executable(test_nested_seed_input test/test_nested_seed_input.f90)
target_link_libraries(test_nested_seed_input PRIVATE milca)
add_executable(test_ideal_map_field_angle test/test_ideal_map_field_angle.f90)
target_link_libraries(test_ideal_map_field_angle PRIVATE milca)
add_test(NAME ideal_map_field_angle COMMAND test_ideal_map_field_angle)
set_tests_properties(ideal_map_field_angle PROPERTIES PROCESSORS 4 TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME ideal_seed_3d_nested_input COMMAND test_nested_seed_input)
set_tests_properties(ideal_seed_3d_nested_input PROPERTIES PROCESSORS 4 TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME ideal_seed_boundary_label
    COMMAND test_ideal_map_axisymmetric seedboundary)
add_test(NAME ideal_map_derivatives_and_rotation_null
    COMMAND test_ideal_map_axisymmetric derivatives)
add_test(NAME ideal_map_n3_leakage COMMAND test_ideal_map_axisymmetric leakage)
foreach(k 2 3)
    add_test(NAME ideal_map_n0_solovev_a3_p${k}
        COMMAND test_ideal_map_axisymmetric solovev ${k})
    add_test(NAME ideal_map_n0_shaped_solovev_p${k}
        COMMAND test_ideal_map_axisymmetric shaped ${k})
endforeach()
# Shaped Solovev ladders at degrees 4 - 6 against the exact flux (no 2D GS
# oracle above P3), and the Lagrange jets and curved meshes for degrees 2 .. 8.
add_test(NAME ideal_map_triangle_lagrange_p2_p8
    COMMAND test_ideal_map_axisymmetric basis)
foreach(k 4 5 6)
    add_test(NAME ideal_map_n0_shaped_solovev_p${k}
        COMMAND test_ideal_map_axisymmetric shapedho ${k})
    set_tests_properties(ideal_map_n0_shaped_solovev_p${k} PROPERTIES
        TIMEOUT 1800 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
endforeach()
# Three-mesh scientific regression ladders reserve four threads.
set_tests_properties(ideal_map_n3_leakage ideal_map_n0_solovev_a3_p2
    ideal_map_n0_solovev_a3_p3 ideal_map_n0_shaped_solovev_p2 PROPERTIES
    TIMEOUT 900 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
# Nested seed with offset axis must return, at order k + 1.
add_test(NAME ideal_map_n0_seed_return_a3_p2 COMMAND test_ideal_map_axisymmetric seed 2)
add_test(NAME ideal_map_n0_seed_return_shaped_p2
    COMMAND test_ideal_map_axisymmetric seedshaped 2)
# P3 ladder: transversal gauge weight 1e-2, initial damping 1e-6, stop at
# 1e-11 x the seed gradient (docs/ideal-map-3d.md). Shaped P3 is open.
add_test(NAME ideal_map_n0_seed_return_a3_p3
    COMMAND test_ideal_map_axisymmetric seed 3 v 7 -1 1 1e-6 1e-2 1e-11)
set_tests_properties(ideal_map_n0_seed_return_a3_p3
    PROPERTIES TIMEOUT 2400 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_seed_return_shaped_p2
    PROPERTIES TIMEOUT 2400 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
# Nested solve: omega eliminated exactly (no gauge form), Newton on W_red(nu)
# with the Schur-complement step; LM damping 1e-4 on nu only (1e-6 reached a
# second stationary point with higher W on the P2 ladder).
add_test(NAME ideal_map_reduced_omega_derivatives
    COMMAND test_ideal_map_axisymmetric reduced)
# Field-line slides g b_pol (and, at P3, the r^2-weighted translation) are
# near-nulls of W_red that vanish with h; the nested axis basis removes them.
# omega = f(chi) is a near-null of W at fixed nu (probes 7, 8), removed by
# the omega gauge moments.
foreach(k 2 3)
    add_test(NAME ideal_map_nested_slide_null_p${k}
        COMMAND test_ideal_map_axisymmetric gauge ${k})
endforeach()
add_test(NAME ideal_map_n0_seed_nested_a3_p2
    COMMAND test_ideal_map_axisymmetric seed 2 v 7 -1 1 1e-4 -1 1e-11 n)
add_test(NAME ideal_map_n0_seed_nested_a3_p3
    COMMAND test_ideal_map_axisymmetric seed 3 v 7 -1 1 1e-4 -1 1e-11 n)
add_test(NAME ideal_map_n0_seed_nested_shaped_p3
    COMMAND test_ideal_map_axisymmetric seedshaped 3 v 7 -1 1 1e-4 -1 1e-11 n)
set_tests_properties(ideal_map_n0_seed_nested_a3_p2 ideal_map_n0_seed_nested_a3_p3
    PROPERTIES TIMEOUT 2400 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_seed_nested_shaped_p3
    PROPERTIES TIMEOUT 4800 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
add_test(NAME ideal_map_factor_cache COMMAND test_ideal_map_axisymmetric factor 2)
add_test(NAME ideal_map_masked_factor_step COMMAND test_ideal_map_axisymmetric factorstep 2)
add_test(NAME ideal_map_krylov_step COMMAND test_ideal_map_axisymmetric krylov 4)
set_tests_properties(ideal_map_krylov_step
    PROPERTIES PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
add_test(NAME ideal_map_lagged_hessian COMMAND test_ideal_map_axisymmetric lagged 4)
set_tests_properties(ideal_map_lagged_hessian
    PROPERTIES PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
add_test(NAME ideal_map_continuation COMMAND test_ideal_map_axisymmetric continuation 4)
set_tests_properties(ideal_map_continuation
    PROPERTIES PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
add_test(NAME ideal_map_n0_normalised_residual_a3_p2
    COMMAND test_ideal_map_axisymmetric residual 2)
set_tests_properties(ideal_map_n0_normalised_residual_a3_p2
    PROPERTIES TIMEOUT 900 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_seed_return_a3_p2
    PROPERTIES TIMEOUT 1800 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_shaped_solovev_p3 PROPERTIES
    TIMEOUT 2400 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
# gs_coupling and gs_free_boundary threads come from the project-wide MILCA_ENABLE_OPENMP option.
add_executable(test_gs_coils test/test_gs_coils.f90)
target_link_libraries(test_gs_coils PRIVATE milca)
add_test(NAME gs_coils_independent_oracles COMMAND test_gs_coils)
add_executable(test_gs_freegs test/test_gs_freegs.f90)
target_link_libraries(test_gs_freegs PRIVATE milca_equilibrium)
add_test(NAME gs_free_boundary_freegs_public_parity
    COMMAND test_gs_freegs ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/freegs)
set_tests_properties(gs_free_boundary_freegs_public_parity PROPERTIES
    TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_library(milca_gs_coupling_test STATIC test/support/gs_coupling_fixture.f90)
target_link_libraries(milca_gs_coupling_test PUBLIC milca)
foreach(coupling_test vacuum fixed_current structure)
    add_executable(test_gs_coupling_${coupling_test} test/test_gs_coupling_${coupling_test}.f90)
    target_link_libraries(test_gs_coupling_${coupling_test} PRIVATE milca_gs_coupling_test)
endforeach()
add_test(NAME gs_coupling_t1_vacuum_coil_flux COMMAND test_gs_coupling_vacuum)
set_tests_properties(gs_coupling_t1_vacuum_coil_flux PROPERTIES TIMEOUT 900)
add_test(NAME gs_coupling_t2_fixed_current_cubature COMMAND test_gs_coupling_fixed_current)
set_tests_properties(gs_coupling_t2_fixed_current_cubature PROPERTIES TIMEOUT 900)
add_test(NAME gs_coupling_structure_symmetry_energy COMMAND test_gs_coupling_structure)
set_tests_properties(gs_coupling_structure_symmetry_energy PROPERTIES TIMEOUT 900)
add_executable(test_gs_steklov_cache test/test_gs_steklov_cache.f90)
target_link_libraries(test_gs_steklov_cache PRIVATE milca)
add_test(NAME gs_steklov_cache_relabelled_curve COMMAND test_gs_steklov_cache)
add_executable(test_gs_fb_interface test/test_gs_fb_interface.f90)
target_link_libraries(test_gs_fb_interface PRIVATE milca)
add_test(NAME gs_free_boundary_interface_flow_quadrature COMMAND test_gs_fb_interface)
add_library(milca_gs_fb_test STATIC test/support/gs_fb_fixture.f90)
target_link_libraries(milca_gs_fb_test PUBLIC milca_gs_coupling_test)
add_executable(test_gs_fb_jacobian test/test_gs_fb_jacobian.f90)
target_link_libraries(test_gs_fb_jacobian PRIVATE milca_gs_fb_test)
add_test(NAME gs_free_boundary_t7_jacobian_finite_differences COMMAND test_gs_fb_jacobian)
add_test(NAME gs_free_boundary_t7_native_lu COMMAND test_gs_fb_jacobian)
set_tests_properties(gs_free_boundary_t7_native_lu PROPERTIES
    TIMEOUT 600 PROCESSORS 2
    ENVIRONMENT "MILCA_FB_SOLVER=lu;OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
set_tests_properties(gs_free_boundary_t7_jacobian_finite_differences PROPERTIES TIMEOUT 900)
add_executable(test_gs_fb_law_tangent test/test_gs_fb_law_tangent.f90)
target_link_libraries(test_gs_fb_law_tangent PRIVATE milca_gs_fb_test)
add_test(NAME gs_free_boundary_law_tangent_finite_differences COMMAND test_gs_fb_law_tangent)
add_executable(test_gs_fb_manufactured test/test_gs_fb_manufactured.f90)
target_link_libraries(test_gs_fb_manufactured PRIVATE milca_gs_fb_test)
add_library(milca_gs_ci_test STATIC test/support/gs_ci_fixture.f90)
target_link_libraries(milca_gs_ci_test PUBLIC milca_gs_fb_test)
add_executable(test_gs_coil_inverse test/test_gs_coil_inverse.f90)
target_link_libraries(test_gs_coil_inverse PRIVATE milca_gs_ci_test)
add_test(NAME gs_coil_inverse_t6_twin COMMAND test_gs_coil_inverse)
set_tests_properties(gs_coil_inverse_t6_twin PROPERTIES TIMEOUT 1800)
add_executable(test_gs_fb_refit test/test_gs_fb_refit.f90)
target_link_libraries(test_gs_fb_refit PRIVATE milca_gs_fb_test)
add_test(NAME gs_free_boundary_p8_refit_loop COMMAND test_gs_fb_refit)
set_tests_properties(gs_free_boundary_p8_refit_loop PROPERTIES TIMEOUT 1800)
add_executable(test_gs_fb_greens test/test_gs_fb_greens.f90)
target_link_libraries(test_gs_fb_greens PRIVATE milca_gs_fb_test)
if(OpenMP_Fortran_FOUND)
    target_compile_options(test_gs_fb_greens PRIVATE ${OpenMP_Fortran_FLAGS})
endif()
add_test(NAME gs_free_boundary_t4_green_function_cross_check COMMAND test_gs_fb_greens)
set_tests_properties(gs_free_boundary_t4_green_function_cross_check PROPERTIES TIMEOUT 1800)
add_test(NAME gs_free_boundary_t3_manufactured_newton COMMAND test_gs_fb_manufactured)
foreach(start reference rp rm zp zm el)
    add_test(NAME equilibrium_d8_psi_offset_${start}
        COMMAND test_gs_fb_manufactured offset ${start})
    set_tests_properties(equilibrium_d8_psi_offset_${start} PROPERTIES
        TIMEOUT 600 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
set_tests_properties(gs_free_boundary_t3_manufactured_newton PROPERTIES TIMEOUT 1800)
add_executable(test_gs_coupled_majorant test/test_gs_coupled_majorant.f90)
target_link_libraries(test_gs_coupled_majorant PRIVATE milca_gs_fb_test)
add_test(NAME gs_free_boundary_p10_coupled_majorant COMMAND test_gs_coupled_majorant)
set_tests_properties(gs_free_boundary_p10_coupled_majorant PROPERTIES TIMEOUT 1800)
add_test(NAME gs_free_boundary_equilibrated_quick
    COMMAND test_gs_coupled_majorant quick)
set_tests_properties(gs_free_boundary_equilibrated_quick PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_bem_residual test/test_gs_bem_residual.f90)
target_link_libraries(test_gs_bem_residual PRIVATE milca_equilibrium)
add_test(NAME gs_bem_single_layer_reference_residual COMMAND test_gs_bem_residual)
set_tests_properties(gs_bem_single_layer_reference_residual PROPERTIES TIMEOUT 30
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_fb_large_aspect test/test_gs_fb_large_aspect.f90)
target_link_libraries(test_gs_fb_large_aspect PRIVATE milca)
add_test(NAME gs_free_boundary_t5_large_aspect_shafranov COMMAND test_gs_fb_large_aspect)
set_tests_properties(gs_free_boundary_t5_large_aspect_shafranov PROPERTIES TIMEOUT 1800)
target_sources(milca_fem2d PRIVATE src/fem2d/curved_triangle.f90
    src/generated/milca_segment_leaves.f90 src/fem2d/triangle_basis.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_majorant.f90 src/equilibrium/gs_profile_law.f90 src/equilibrium/gs_nonlinear.f90
    src/equilibrium/gs_equilibrated.f90 src/generated/milca_gs_error_leaves.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_quadrature_budget.f90)
set_source_files_properties(
    src/equilibrium/gs_quadrature_budget.f90
    src/equilibrium/gs_error_control.f90
    src/generated/milca_gs_error_leaves.f90
    src/generated/milca_grad_shafranov_leaves.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_rt0_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_rt0_majorant.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_bernstein_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_bernstein_majorant.f90
    PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
add_executable(test_gs_checked_majorant test/test_gs_checked_majorant.f90)
target_link_libraries(test_gs_checked_majorant PRIVATE milca_equilibrium)
add_test(NAME grad_shafranov_checked_bdm_majorant COMMAND test_gs_checked_majorant)
set_tests_properties(grad_shafranov_checked_bdm_majorant PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")
add_executable(test_gs_checked_source_remainder test/test_gs_checked_source_remainder.f90)
target_link_libraries(test_gs_checked_source_remainder PRIVATE milca_equilibrium)
target_compile_options(test_gs_checked_source_remainder PRIVATE
    -fno-fast-math -ffp-contract=off)
add_test(NAME grad_shafranov_checked_source_remainder
    COMMAND test_gs_checked_source_remainder)
set_tests_properties(grad_shafranov_checked_source_remainder PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_test(NAME grad_shafranov_checked_source_study
    COMMAND test_gs_checked_source_remainder study)
set_tests_properties(grad_shafranov_checked_source_study PROPERTIES
    TIMEOUT 300 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_property(SOURCE src/generated/milca_gs_error_leaves.f90 APPEND PROPERTY
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_gs_error_leaves test/test_gs_error_leaves.f90)
target_link_libraries(test_gs_error_leaves PRIVATE milca_equilibrium)
add_test(NAME grad_shafranov_error_leaf_oracles COMMAND test_gs_error_leaves)
set_tests_properties(grad_shafranov_error_leaf_oracles PROPERTIES TIMEOUT 20)
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_error_control.f90)
add_executable(test_gs_error_control test/test_gs_error_control.f90)
target_link_libraries(test_gs_error_control PRIVATE milca)
add_test(NAME gs_fixed_stability_and_linear_goal COMMAND test_gs_error_control)
add_executable(test_gs_verified_stability test/test_gs_verified_stability.f90)
target_link_libraries(test_gs_verified_stability PRIVATE milca_equilibrium)
add_test(NAME gs_verified_nonlinear_coercivity COMMAND test_gs_verified_stability)
set_tests_properties(gs_verified_nonlinear_coercivity PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
set_tests_properties(gs_fixed_stability_and_linear_goal PROPERTIES
    TIMEOUT 30 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_equilibrated_majorant test/test_gs_equilibrated_majorant.f90)
target_link_libraries(test_gs_equilibrated_majorant PRIVATE milca)
add_test(NAME grad_shafranov_equilibrated_flux_majorant
    COMMAND test_gs_equilibrated_majorant)
set_tests_properties(grad_shafranov_equilibrated_flux_majorant PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_gs_p2 test/test_gs_p2.f90)
target_link_libraries(test_gs_p2 PRIVATE milca)
add_test(NAME grad_shafranov_p2_physical_oracles COMMAND test_gs_p2)
add_executable(test_gs_p2_gauge test/test_gs_p2_gauge.f90)
target_link_libraries(test_gs_p2_gauge PRIVATE milca)
add_test(NAME grad_shafranov_p2_gauge_oracles COMMAND test_gs_p2_gauge)
add_executable(test_gs_p2_refinement test/test_gs_p2_refinement.f90)
target_link_libraries(test_gs_p2_refinement PRIVATE milca)
add_test(NAME grad_shafranov_p2_polynomial_refinement COMMAND test_gs_p2_refinement)
set_tests_properties(grad_shafranov_p2_polynomial_refinement PROPERTIES
    FIXTURES_SETUP gs_p2_solovev_export)
add_executable(test_boozer_file test/test_boozer_file.f90)
target_link_libraries(test_boozer_file PRIVATE milca)
# NEO-2 (e932b64c) on libneo's public tok02 file at s = 0.50781: iota,
# diota/ds, <|grad s|> [1/cm] and <B^2> = avbhat2*bmod0^2 [T^2].
add_test(NAME equilibrium_boozer_file_geometry COMMAND test_boozer_file
    ${libneo_SOURCE_DIR}/python/tests/circ.bc 0.50781 0.76312999999999
    -0.5561487528261847 0.031219209564037516 3.9494899)
milca_test_labels(equilibrium_boozer_file_geometry fast CLASSES equilibrium
    TARGETS test_boozer_file)
add_executable(test_gs_coarea test/test_gs_coarea.f90)
target_link_libraries(test_gs_coarea PRIVATE milca)
add_test(NAME grad_shafranov_p2_native_coarea_oracles COMMAND test_gs_coarea)
add_executable(test_gs_inverse_cache_reuse test/test_gs_inverse_cache_reuse.f90)
target_link_libraries(test_gs_inverse_cache_reuse PRIVATE milca)
add_test(NAME grad_shafranov_inverse_cache_reuse_oracles COMMAND test_gs_inverse_cache_reuse)
add_executable(test_gs_inverse_rhs_batch test/test_gs_inverse_rhs_batch.f90)
target_link_libraries(test_gs_inverse_rhs_batch PRIVATE milca)
add_test(NAME grad_shafranov_inverse_rhs_batch_oracles COMMAND test_gs_inverse_rhs_batch)
add_executable(test_gs_inverse_boundary_response test/test_gs_inverse_boundary_response.f90)
target_link_libraries(test_gs_inverse_boundary_response PRIVATE milca)
add_test(NAME grad_shafranov_inverse_boundary_response_oracles
    COMMAND test_gs_inverse_boundary_response)
add_executable(test_gs_inverse_cache_equation test/test_gs_inverse_cache_equation.f90)
target_link_libraries(test_gs_inverse_cache_equation PRIVATE milca)
add_test(NAME grad_shafranov_inverse_cache_equation_oracle COMMAND test_gs_inverse_cache_equation)
add_executable(test_gs_inverse_source_cut test/test_gs_inverse_source_cut.f90)
target_link_libraries(test_gs_inverse_source_cut PRIVATE milca)
add_test(NAME grad_shafranov_inverse_source_cut_oracles COMMAND test_gs_inverse_source_cut)
add_executable(test_gs_inverse_mean test/test_gs_inverse_mean.f90)
target_link_libraries(test_gs_inverse_mean PRIVATE milca)
add_test(NAME grad_shafranov_profile_mean_gauge_oracles COMMAND test_gs_inverse_mean)
add_executable(test_gs_inverse test/test_gs_inverse.f90)
target_link_libraries(test_gs_inverse PRIVATE milca)
add_executable(test_gs_prescribed_q test/test_gs_prescribed_q.f90)
target_link_libraries(test_gs_prescribed_q PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_p3_convergence COMMAND test_gs_prescribed_q)
set_tests_properties(grad_shafranov_prescribed_q_p3_convergence PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_variable_p test/test_gs_prescribed_q_variable_p.f90)
target_link_libraries(test_gs_prescribed_q_variable_p PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_variable_p COMMAND test_gs_prescribed_q_variable_p)
set_tests_properties(grad_shafranov_prescribed_q_variable_p PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_intervals test/test_gs_prescribed_q_intervals.f90)
target_link_libraries(test_gs_prescribed_q_intervals PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_jacobian_richardson
         COMMAND test_gs_prescribed_q_intervals jacobian)
add_test(NAME grad_shafranov_prescribed_q_exact_vs_secant
         COMMAND test_gs_prescribed_q_intervals secant)
add_test(NAME grad_shafranov_prescribed_q_knots64_uniform
         COMMAND test_gs_prescribed_q_intervals mms 64 0 8)
add_test(NAME grad_shafranov_prescribed_q_knots64_graded
         COMMAND test_gs_prescribed_q_intervals mms 64 1 8)
set_tests_properties(grad_shafranov_prescribed_q_jacobian_richardson
    grad_shafranov_prescribed_q_exact_vs_secant PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_tests_properties(grad_shafranov_prescribed_q_knots64_uniform
    grad_shafranov_prescribed_q_knots64_graded PROPERTIES
    TIMEOUT 3600 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_pedestal test/test_gs_prescribed_q_pedestal.f90)
target_link_libraries(test_gs_prescribed_q_pedestal PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_pedestal COMMAND test_gs_prescribed_q_pedestal)
set_tests_properties(grad_shafranov_prescribed_q_pedestal PROPERTIES
    TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_edge_mode test/test_gs_prescribed_q_edge_mode.f90)
target_link_libraries(test_gs_prescribed_q_edge_mode PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_edge_mode
         COMMAND test_gs_prescribed_q_edge_mode)
set_tests_properties(grad_shafranov_prescribed_q_edge_mode PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_circle test/test_gs_prescribed_q_circle.f90)
target_link_libraries(test_gs_prescribed_q_circle PRIVATE milca)
add_test(NAME grad_shafranov_prescribed_q_e1 COMMAND test_gs_prescribed_q_circle 1)
add_test(NAME grad_shafranov_prescribed_q_e2 COMMAND test_gs_prescribed_q_circle 2)
set_tests_properties(grad_shafranov_prescribed_q_e1 grad_shafranov_prescribed_q_e2
    PROPERTIES TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME grad_shafranov_phase1_inverse_app COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_phase1> -DDECK=$<TARGET_FILE:test_gs_prescribed_q_circle>
    -DWORK=${CMAKE_BINARY_DIR}/phase1_inverse_app
    -P ${CMAKE_SOURCE_DIR}/test/gs_phase1_inverse.cmake)
add_test(NAME grad_shafranov_phase1_inverse_pressure_law_app COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_phase1> -DDECK=$<TARGET_FILE:test_gs_prescribed_q_circle>
    -DWORK=${CMAKE_BINARY_DIR}/phase1_inverse_pressure_law_app
    -P ${CMAKE_SOURCE_DIR}/test/gs_phase1_inverse_pressure_law.cmake)
add_test(NAME grad_shafranov_constant_q_inverse_unit COMMAND test_gs_inverse)
add_test(NAME grad_shafranov_constant_q_inverse_refinement COMMAND test_gs_inverse fine)
set_tests_properties(grad_shafranov_constant_q_inverse_unit grad_shafranov_constant_q_inverse_refinement PROPERTIES TIMEOUT 120)
set_tests_properties(grad_shafranov_p2_physical_oracles PROPERTIES
    FIXTURES_SETUP gs_native_mesh_export)
add_executable(grad_shafranov_demo app/grad_shafranov_demo.f90)
target_link_libraries(grad_shafranov_demo PRIVATE milca milca_io fortplot)
add_test(NAME grad_shafranov_inverse_app_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:grad_shafranov_demo> -DWORK=${CMAKE_BINARY_DIR}/gs_inverse_app_check
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_inverse_app.cmake)
set_tests_properties(grad_shafranov_inverse_app_contract PROPERTIES
    TIMEOUT 300 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_p1_weak_current test/test_gs_p1_weak_current.f90)
target_link_libraries(test_gs_p1_weak_current PRIVATE milca)
add_test(NAME grad_shafranov_p1_weak_current COMMAND test_gs_p1_weak_current)
foreach(mutation omit-jump omit-cylindrical wrong-current-sign)
    add_test(NAME grad_shafranov_p1_rejects_${mutation}
        COMMAND test_gs_p1_weak_current ${mutation})
    set_tests_properties(grad_shafranov_p1_rejects_${mutation} PROPERTIES WILL_FAIL TRUE)
endforeach()
add_executable(gold_hoyle_demo app/gold_hoyle_demo.f90)
add_executable(gs_reference_error app/gs_reference_error.f90)
target_link_libraries(gs_reference_error PRIVATE milca_gs_reference_test)
configure_file(test/fixtures/gs_reference_rejection.in
    ${CMAKE_BINARY_DIR}/gs_reference_rejection.in COPYONLY)
add_test(NAME grad_shafranov_p1_reference_accepts_mesh_v1
    COMMAND gs_reference_error p1_native_oracle_mesh.dat gs_reference_rejection.in)
add_test(NAME grad_shafranov_p1_reference_rejects_mesh_v2
    COMMAND ${CMAKE_COMMAND} -DGS_REFERENCE_ERROR=$<TARGET_FILE:gs_reference_error>
        -DMESH=${CMAKE_BINARY_DIR}/p2_native_oracle_mesh.dat
        -DPARAMETERS=${CMAKE_BINARY_DIR}/gs_reference_rejection.in
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_gs_p2_rejection.cmake)
set_tests_properties(grad_shafranov_p1_reference_accepts_mesh_v1
    grad_shafranov_p1_reference_rejects_mesh_v2 PROPERTIES
    FIXTURES_REQUIRED gs_native_mesh_export)
add_executable(test_gs_reference_error test/test_gs_reference_error.f90)
target_link_libraries(test_gs_reference_error PRIVATE milca_gs_reference_test)
add_test(NAME grad_shafranov_reference_cell_enclosures COMMAND test_gs_reference_error)
add_executable(test_gs_p2_reference_error test/test_gs_p2_reference_error.f90)
target_link_libraries(test_gs_p2_reference_error PRIVATE milca_gs_reference_test)
add_test(NAME grad_shafranov_p2_reference_cell_enclosures COMMAND test_gs_p2_reference_error)
configure_file(test/fixtures/gs_p2_reference.in
    ${CMAKE_BINARY_DIR}/gs_p2_reference.in COPYONLY)
add_test(NAME grad_shafranov_p2_reference_rejects_curved_metadata
    COMMAND ${CMAKE_COMMAND} -DGS_REFERENCE_ERROR=$<TARGET_FILE:gs_reference_error>
        -DMESH=${CMAKE_BINARY_DIR}/p2_curved_metadata_mesh.dat
        -DPARAMETERS=${CMAKE_BINARY_DIR}/gs_p2_reference.in
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_gs_p2_geometry_rejection.cmake)
set_tests_properties(grad_shafranov_p2_reference_rejects_curved_metadata PROPERTIES
    FIXTURES_REQUIRED gs_p2_solovev_export)

target_link_libraries(gold_hoyle_demo PRIVATE milca milca_io)
foreach(layer IN LISTS MILCA_LAYER_LIBRARIES)
    target_compile_options(${layer} PRIVATE -Wall -Wextra -fno-fast-math
        -ffp-contract=off -frecursive)
endforeach()
add_executable(milca_demo app/milca_demo.f90)
find_package(OpenMP COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
    # Inverse Jacobian columns and other independent loops use OpenMP.
    target_link_libraries(milca_core PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_link_libraries(milca_demo PRIVATE milca fortplot)
add_executable(steady_laboratory_demo app/steady_laboratory_demo.f90)
target_link_libraries(steady_laboratory_demo PRIVATE milca fortplot)
add_executable(transport_demo app/transport_demo.f90)
target_link_libraries(transport_demo PRIVATE milca fortplot)
add_executable(streaming_comparison_demo app/streaming_comparison_demo.f90)
target_link_libraries(streaming_comparison_demo PRIVATE milca fortplot)
foreach(example streaming_moments debye_sphere slab_oscillation shielding_dynamics interaction)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE milca fortplot)
endforeach()
add_executable(microscopic_cloud_demo app/microscopic_cloud_demo.f90)
target_link_libraries(microscopic_cloud_demo PRIVATE milca fortplot)
foreach(example quantum_reference classical_pair coulomb_encounter jellium_cell)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE milca fortplot)
endforeach()
# Presentation-only palette and number formatting for the gallery demos.
add_library(milca_demo_style app/demo_style.f90)
target_compile_options(milca_demo_style PRIVATE -Wall -Wextra)
foreach(application milca_demo steady_laboratory_demo transport_demo
        streaming_comparison_demo streaming_moments_demo debye_sphere_demo
        slab_oscillation_demo shielding_dynamics_demo interaction_demo
        microscopic_cloud_demo classical_pair_demo)
    target_link_libraries(${application} PRIVATE milca_demo_style)
endforeach()
add_executable(test_demo_style test/test_demo_style.f90)
target_link_libraries(test_demo_style PRIVATE milca_demo_style)
add_test(NAME demo_style_numbers COMMAND test_demo_style)
add_executable(formal_krook_certificate app/formal_krook_certificate.f90)
target_link_libraries(formal_krook_certificate PRIVATE milca)
add_executable(test_laboratory test/test_laboratory.f90)
target_link_libraries(test_laboratory PRIVATE milca)
add_test(NAME characteristic_and_steady_oracles COMMAND test_laboratory)
add_executable(test_poisson test/test_poisson.f90)
target_link_libraries(test_poisson PRIVATE milca)
add_test(NAME poisson_and_coupled_steady_oracles COMMAND test_poisson)
add_executable(test_fp test/test_fp.f90)
target_link_libraries(test_fp PRIVATE milca)
add_test(NAME bounded_fp_steady_oracles COMMAND test_fp)
add_executable(test_rosenbluth_reference test/test_rosenbluth_reference.f90)
target_link_libraries(test_rosenbluth_reference PRIVATE milca)
foreach(check maxwellian l1_brute invariants slow_limit)
    add_test(NAME rosenbluth_${check}_oracle COMMAND test_rosenbluth_reference ${check})
endforeach()
# Homogeneous FP relaxation (B04): the same executable, one check per name.
add_executable(test_ntv_axisymmetric test/test_ntv_axisymmetric.f90)
target_link_libraries(test_ntv_axisymmetric PRIVATE milca)
foreach(check onsager momentum local collisionless)
    add_test(NAME ntv_axisymmetric_${check} COMMAND test_ntv_axisymmetric ${check})
    milca_test_labels(ntv_axisymmetric_${check} fast CLASSES continuum
        TARGETS test_ntv_axisymmetric)
endforeach()
add_executable(boozer_orbit_scan app/boozer_orbit_scan.f90)
target_link_libraries(boozer_orbit_scan PRIVATE milca)
add_executable(test_boozer_orbit test/test_boozer_orbit.f90)
target_link_libraries(test_boozer_orbit PRIVATE milca)
foreach(check deeply_trapped thin_bounce invariants pitch_map shear potential pitch_bracket)
    add_test(NAME boozer_orbit_${check} COMMAND test_boozer_orbit ${check})
    milca_test_labels(boozer_orbit_${check} fast CLASSES equilibrium
        TARGETS test_boozer_orbit)
endforeach()
add_executable(test_dk_ntv test/test_dk_ntv.f90)
target_link_libraries(test_dk_ntv PRIVATE milca)
add_test(NAME dk_ntv_jvp COMMAND test_dk_ntv jvp)
milca_test_labels(dk_ntv_jvp fast CLASSES continuum TARGETS test_dk_ntv)
add_test(NAME dk_ntv_assembly COMMAND test_dk_ntv assembly)
milca_test_labels(dk_ntv_assembly fast CLASSES continuum TARGETS test_dk_ntv)
add_test(NAME dk_ntv_tangents COMMAND test_dk_ntv tangents)
milca_test_labels(dk_ntv_tangents fast CLASSES continuum TARGETS test_dk_ntv)
add_test(NAME dk_ntv_rotation_limit COMMAND test_dk_ntv rotation_limit)
milca_test_labels(dk_ntv_rotation_limit fast CLASSES continuum TARGETS test_dk_ntv)
add_test(NAME dk_ntv_one_over_nu_limit COMMAND test_dk_ntv one_over_nu_limit)
milca_test_labels(dk_ntv_one_over_nu_limit fast CLASSES continuum TARGETS test_dk_ntv)
foreach(check ql_dpar ql_f10 ql_oracle ql_collisionless f10_provider reuse)
    add_test(NAME dk_ntv_${check} COMMAND test_dk_ntv ${check})
    milca_test_labels(dk_ntv_${check} fast CLASSES continuum TARGETS test_dk_ntv)
endforeach()
add_test(NAME dk_ntv_line COMMAND test_dk_ntv line)
milca_test_labels(dk_ntv_line fast CLASSES continuum TARGETS test_dk_ntv)
add_executable(test_dk_collop test/test_dk_collop.f90)
target_link_libraries(test_dk_collop PRIVATE milca)
foreach(check maxwellian_null rates jvp zero_jvp fp_moments conservation self_adjoint conservation_pair
        exchange_pair)
    add_test(NAME dk_collop_${check} COMMAND test_dk_collop ${check})
    milca_test_labels(dk_collop_${check} fast CLASSES continuum TARGETS test_dk_collop)
endforeach()
add_test(NAME dk_collop_neo2_matrix COMMAND test_dk_collop neo2_matrix
    ${CMAKE_SOURCE_DIR}/test/data/neo2_collop_aug30835.txt)
milca_test_labels(dk_collop_neo2_matrix fast CLASSES continuum TARGETS test_dk_collop)
add_executable(test_dk_mono test/test_dk_mono.f90)
target_link_libraries(test_dk_mono PRIVATE milca)
foreach(check collisionless banana pfirsch_schlueter plateau onsager jvp estimate
        factor_api coupling quasisymmetric_3d)
    add_test(NAME dk_mono_${check} COMMAND test_dk_mono ${check})
    milca_test_labels(dk_mono_${check} fast CLASSES continuum TARGETS test_dk_mono)
endforeach()
# Against NEO-2-PAR in Lorentz mode on libneo's public tok02 file (circ.bc).
add_test(NAME dk_mono_neo2_tok02 COMMAND test_dk_mono neo2_tok02
    ${libneo_SOURCE_DIR}/python/tests/circ.bc)
milca_test_labels(dk_mono_neo2_tok02 fast CLASSES continuum equilibrium
    TARGETS test_dk_mono)
add_executable(test_dk_hp test/test_dk_hp.f90)
target_link_libraries(test_dk_hp PRIVATE milca)
foreach(check reference onsager collisionless banana pfirsch_schlueter estimate
        ntv_reference ntv_superbanana_plateau ntv_one_over_nu ntv_rotation)
    add_test(NAME dk_hp_${check} COMMAND test_dk_hp ${check})
    milca_test_labels(dk_hp_${check} fast CLASSES continuum TARGETS test_dk_hp)
endforeach()
add_executable(test_dk_pitch test/test_dk_pitch.f90)
target_link_libraries(test_dk_pitch PRIVATE milca)
add_test(NAME dk_pitch_lorentz_operator COMMAND test_dk_pitch)
milca_test_labels(dk_pitch_lorentz_operator fast CLASSES continuum TARGETS test_dk_pitch)
add_executable(test_fp_relaxation test/test_fp_relaxation.f90)
target_link_libraries(test_fp_relaxation PRIVATE milca)
foreach(check maxwellian mass slow_rate convergence_v convergence_dt)
    add_test(NAME fp_relaxation_${check} COMMAND test_fp_relaxation ${check})
    milca_test_labels(fp_relaxation_${check} fast CLASSES continuum
        TARGETS test_fp_relaxation)
endforeach()
# First-order thermal-force moments (B05-B07): the same executable, one check per name.
add_executable(test_screening_moments test/test_screening_moments.f90)
target_link_libraries(test_screening_moments PRIVATE milca)
foreach(check maxwellian zero_density rate0 strong_field flr_ratio convergence_v
        oracle_rosenbluth oracle_poisson temperature_derivative transient steady_transient)
    add_test(NAME screening_moments_${check} COMMAND test_screening_moments ${check})
    milca_test_labels(screening_moments_${check} fast CLASSES continuum
        TARGETS test_screening_moments)
endforeach()
add_executable(test_transport test/test_transport.f90)
target_link_libraries(test_transport PRIVATE milca)
add_test(NAME general_transport_and_tangent_oracles COMMAND test_transport)
add_executable(test_transport_sampling test/test_transport_sampling.f90)
target_link_libraries(test_transport_sampling PRIVATE milca)
add_test(NAME cached_transport_band_oracles COMMAND test_transport_sampling)
add_executable(test_transport_linear_bands test/test_transport_linear_bands.f90)
target_link_libraries(test_transport_linear_bands PRIVATE milca)
add_test(NAME continuous_linear_band_oracles COMMAND test_transport_linear_bands)
add_executable(test_streaming_comparison test/test_streaming_comparison.f90)
target_link_libraries(test_streaming_comparison PRIVATE milca)
add_test(NAME streaming_representation_oracles COMMAND test_streaming_comparison)
add_executable(test_streaming_particle_density test/test_streaming_particle_density.f90)
target_link_libraries(test_streaming_particle_density PRIVATE milca)
add_test(NAME streaming_particle_density_oracles COMMAND test_streaming_particle_density)
foreach(example streaming_moments debye_sphere slab_oscillation shielding_dynamics interaction)
    add_executable(test_${example} test/test_${example}.f90)
    target_link_libraries(test_${example} PRIVATE milca)
    add_test(NAME ${example}_oracles COMMAND test_${example})
endforeach()
add_executable(test_nbody test/test_nbody.f90)
target_link_libraries(test_nbody PRIVATE milca)
add_test(NAME nbody_coulomb_species_oracles COMMAND test_nbody)
add_executable(test_pme test/test_pme.f90)
target_link_libraries(test_pme PRIVATE milca)
add_test(NAME pme_direct_ewald_oracles COMMAND test_pme)
add_executable(test_ocp test/test_ocp.f90)
target_link_libraries(test_ocp PRIVATE milca)
add_test(NAME ocp_vlasov_debye_oracles COMMAND test_ocp)
add_executable(test_cloud test/test_cloud.f90)
target_link_libraries(test_cloud PRIVATE milca)
add_test(NAME microscopic_cloud_oracles COMMAND test_cloud)
foreach(quantum_family fermi_gas quantum_fock quantum_evolution quantum_correlations
        quantum_reference_gallery classical_pair jellium_cell)
    add_executable(test_${quantum_family} test/test_${quantum_family}.f90)
    target_link_libraries(test_${quantum_family} PRIVATE milca)
    add_test(NAME ${quantum_family}_oracles COMMAND test_${quantum_family})
endforeach()
add_executable(pp_fusion_export app/pp_fusion_export.f90)
add_dependencies(pp_fusion_export milca_build_info)
target_include_directories(pp_fusion_export PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(pp_fusion_export PRIVATE milca milca_io)

add_test(NAME pp_fusion_export_contract
    COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:pp_fusion_export>
        -DWORK=${CMAKE_BINARY_DIR}/pp-fusion-export-test
        -P ${CMAKE_SOURCE_DIR}/test/check_pp_fusion_export.cmake)
set_tests_properties(pp_fusion_export_contract PROPERTIES PROCESSORS 4)

add_executable(nested_oracle_export app/nested_oracle_export.f90)
target_link_libraries(nested_oracle_export PRIVATE milca milca_io)
add_executable(nested_equilibrium app/nested_equilibrium.f90)
target_link_libraries(nested_equilibrium PRIVATE milca milca_io)
add_executable(test_nested_equilibrium_app test/test_nested_equilibrium_app.f90)
target_link_libraries(test_nested_equilibrium_app PRIVATE milca_io)
# PARVMEC (ns 101, mpol 8) on the same inputs: raxis, wb and, for iota, the
# exact volume 2 pi^2 R a^2; GVEC gives raxis 3.03683 m for the iota case.
foreach(case IN ITEMS "circular-tokamak|3.0366;0.0750949;59.2176264"
        "circular-tokamak-current|3.04266;0.0755050")
    string(REPLACE "|" ";" parts "${case}")
    list(GET parts 0 name)
    list(SUBLIST parts 1 -1 expected)
    string(REPLACE ";" " " expected "${expected}")
    add_test(NAME equilibrium_nested_${name} COMMAND ${CMAKE_COMMAND}
        -DAPP=$<TARGET_FILE:nested_equilibrium>
        -DCHECK=$<TARGET_FILE:test_nested_equilibrium_app>
        -DCONFIG=${CMAKE_SOURCE_DIR}/examples/nested-equilibrium/${name}.nml
        -DWORK=${CMAKE_BINARY_DIR}/nested-equilibrium-${name}
        "-DEXPECTED=${expected}"
        -P ${CMAKE_SOURCE_DIR}/test/check_nested_equilibrium.cmake)
    set_tests_properties(equilibrium_nested_${name} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
endforeach()
add_executable(pp_fusion_eft_export app/pp_fusion_eft_export.f90)
add_dependencies(pp_fusion_eft_export milca_build_info)
target_include_directories(pp_fusion_eft_export PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(pp_fusion_eft_export PRIVATE milca milca_io)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pp_fusion_eft_export_contract)
add_test(NAME pp_fusion_eft_export_contract
    COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:pp_fusion_eft_export>
        -DWORK=${CMAKE_BINARY_DIR}/pp_fusion_eft_export_contract
        -P ${CMAKE_SOURCE_DIR}/test/check_pp_fusion_eft_export.cmake)
set_tests_properties(pp_fusion_eft_export_contract PROPERTIES PROCESSORS 4)

add_executable(test_pp_fusion_eft test/test_pp_fusion_eft.f90)
target_link_libraries(test_pp_fusion_eft PRIVATE milca)
add_test(NAME pp_fusion_eft_oracles COMMAND test_pp_fusion_eft
    ${CMAKE_SOURCE_DIR}/test/reference/pp_fusion_eft.dat)
set_tests_properties(pp_fusion_eft_oracles PROPERTIES PROCESSORS 4)

add_executable(test_pp_fusion test/test_pp_fusion.f90)
target_link_libraries(test_pp_fusion PRIVATE milca)
add_test(NAME pp_fusion_oracles COMMAND test_pp_fusion
    ${CMAKE_SOURCE_DIR}/test/reference/pp_fusion_ode.dat)
set_tests_properties(pp_fusion_oracles PROPERTIES PROCESSORS 4)

add_executable(test_coulomb_encounter test/test_coulomb_encounter.f90)
target_link_libraries(test_coulomb_encounter PRIVATE milca)
add_test(NAME coulomb_encounter_oracles COMMAND test_coulomb_encounter)
add_executable(test_quantum_cutoff test/test_quantum_cutoff.f90)
target_link_libraries(test_quantum_cutoff PRIVATE milca)
add_test(NAME quantum_cutoff_oracles COMMAND test_quantum_cutoff)
find_program(MILCA_LAKE lake)
option(MILCA_CHECK_LEAN "Check continuous Lean theorems and computed certificate" OFF)
if(MILCA_CHECK_LEAN)
    if(NOT MILCA_LAKE)
        message(FATAL_ERROR "Lean/Lake is required; see proof/README.md for setup")
    endif()
    add_test(NAME lean_continuous_certificates
        COMMAND bash ${CMAKE_SOURCE_DIR}/proof/check.sh
            $<TARGET_FILE:formal_krook_certificate>)
    set_tests_properties(lean_continuous_certificates PROPERTIES
        TIMEOUT 300 PROCESSORS 2)
endif()
set(MILCA_FORTSYM_BUILD "" CACHE PATH "Optional FortSym build tree for regeneration")
set(MILCA_FORTSYM_REVISION "b9ea95c5a0475b016b1e149875a9bec9d8562533")
if(MILCA_FORTSYM_BUILD)
    file(STRINGS "${MILCA_FORTSYM_BUILD}/CMakeCache.txt" fortsym_source
        REGEX "^CMAKE_HOME_DIRECTORY:INTERNAL=")
    string(REPLACE "CMAKE_HOME_DIRECTORY:INTERNAL=" "" fortsym_source
        "${fortsym_source}")
    execute_process(COMMAND git -C "${fortsym_source}" rev-parse HEAD
        OUTPUT_VARIABLE fortsym_revision OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE fortsym_revision_status)
    if(NOT fortsym_revision_status EQUAL 0 OR
            NOT fortsym_revision STREQUAL MILCA_FORTSYM_REVISION)
        message(FATAL_ERROR "Regeneration requires FortSym ${MILCA_FORTSYM_REVISION}")
    endif()
    configure_file(gen/milca_fortsym_revision.inc.in
        ${CMAKE_BINARY_DIR}/gen_include/milca_fortsym_revision.inc @ONLY)
    file(STRINGS "${CMAKE_SOURCE_DIR}/gen/ideal-mhd-revision.txt"
        MILCA_IDEAL_MHD_REVISION)
    FetchContent_Declare(milca_shared_mhd
        URL "https://github.com/itpplasma/GLISS/archive/${MILCA_IDEAL_MHD_REVISION}.tar.gz"
        SOURCE_SUBDIR gen/codegen_sources_only)
    FetchContent_MakeAvailable(milca_shared_mhd)
    set(shared_mhd_gen "${milca_shared_mhd_SOURCE_DIR}/gen/ideal_mhd")
    add_library(milca_gen_support STATIC gen/milca_gen_support.f90
        gen/milca_symbolic_physics.f90 gen/milca_symbolic_ideal.f90
        ${shared_mhd_gen}/ideal_mhd_symbolic.f90)
    target_include_directories(milca_gen_support PUBLIC
        ${MILCA_FORTSYM_BUILD}/include ${CMAKE_BINARY_DIR}/gen_include)
    target_link_libraries(milca_gen_support PUBLIC
        ${MILCA_FORTSYM_BUILD}/lib/libfortsym${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_executable(test_gs_p1_current_identities test/test_gs_p1_current_identities.f90)
    target_link_libraries(test_gs_p1_current_identities PRIVATE milca_gen_support)
    add_test(NAME grad_shafranov_p1_exact_current_antiderivatives
        COMMAND test_gs_p1_current_identities)
    add_executable(gen_vrmhd_island gen/gen_vrmhd_island.f90)
    target_link_libraries(gen_vrmhd_island PRIVATE milca_gen_support)
    add_test(NAME vrmhd_island_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_island>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/test/support/vrmhd_island_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-island-regeneration
        -DGENERATED_DIR=test/support -DLEAF=vrmhd_island_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_island_regeneration PROPERTIES LABELS rm-b5)
    add_executable(gen_vrmhd_mms gen/gen_vrmhd_mms.f90)
    target_link_libraries(gen_vrmhd_mms PRIVATE milca_gen_support)
    add_test(NAME vrmhd_mms_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_mms>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/test/support/vrmhd_mms_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-mms-regeneration
        -DGENERATED_DIR=test/support -DLEAF=vrmhd_mms_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_mms_regeneration PROPERTIES LABELS rm2)
    add_executable(gen_vrmhd_mms_fourier gen/gen_vrmhd_mms_fourier.f90)
    target_link_libraries(gen_vrmhd_mms_fourier PRIVATE milca_gen_support)
    add_test(NAME vrmhd_mms_fourier_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_mms_fourier>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/test/support/vrmhd_mms_fourier_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-mms-fourier-regeneration
        -DGENERATED_DIR=test/support -DLEAF=vrmhd_mms_fourier_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_mms_fourier_regeneration PROPERTIES LABELS rm3)
    add_test(NAME vrmhd_mms_nonlinear_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_mms_fourier>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/test/support/vrmhd_mms_nonlinear_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-mms-nonlinear-regeneration
        -DGENERATED_DIR=test/support -DLEAF=vrmhd_mms_nonlinear_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_mms_nonlinear_regeneration PROPERTIES LABELS rm3b)
    add_executable(gen_vrmhd_gs gen/gen_vrmhd_gs.f90)
    target_link_libraries(gen_vrmhd_gs PRIVATE milca_gen_support)
    add_test(NAME vrmhd_gs_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_gs>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_vrmhd_gs_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-gs-regeneration
        -DLEAF=milca_vrmhd_gs_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_gs_regeneration PROPERTIES LABELS rm2)
    add_executable(gen_vrmhd_section gen/gen_vrmhd_section.f90)
    target_link_libraries(gen_vrmhd_section PRIVATE milca_gen_support)
    add_test(NAME vrmhd_section_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_section>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_vrmhd_section_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-section-regeneration
        -DLEAF=milca_vrmhd_section_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    set_tests_properties(vrmhd_section_regeneration PROPERTIES LABELS rm2)
    add_executable(gen_gs_symmetry gen/gen_gs_symmetry.f90)
    target_link_libraries(gen_gs_symmetry PRIVATE milca_gen_support)
    add_test(NAME grad_shafranov_rm0_symbolic
        COMMAND ${CMAKE_COMMAND}
            -DGENERATOR=$<TARGET_FILE:gen_gs_symmetry>
            -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_gs_symmetry_leaves.f90
            -DWORK=${CMAKE_BINARY_DIR}/rm0_regeneration
            -DLEAF=milca_gs_symmetry_leaves.f90
            -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    add_executable(gen_ideal_fb_kink gen/gen_ideal_fb_kink.f90)
    target_link_libraries(gen_ideal_fb_kink PRIVATE milca_gen_support)
    add_executable(test_gs_symbolic test/test_gs_symbolic.f90)
    target_link_libraries(test_gs_symbolic PRIVATE milca_gen_support)
    add_test(NAME grad_shafranov_symbolic_identities COMMAND test_gs_symbolic)
    add_executable(test_gs_periodic_cylinder_symbolic
        test/test_gs_periodic_cylinder_symbolic.f90)
    target_link_libraries(test_gs_periodic_cylinder_symbolic PRIVATE milca_gen_support)
    add_test(NAME grad_shafranov_periodic_cylinder_exact_identities
        COMMAND test_gs_periodic_cylinder_symbolic)
    add_executable(test_gs_variational_identities test/test_gs_variational_identities.f90)
    target_link_libraries(test_gs_variational_identities PRIVATE milca_gen_support)
    add_executable(test_gs_forward_majorant_identities
        test/test_gs_forward_majorant_identities.f90)
    target_link_libraries(test_gs_forward_majorant_identities PRIVATE milca_gen_support)
    add_executable(test_gs_coupled_majorant_identities
        test/test_gs_coupled_majorant_identities.f90)
    target_link_libraries(test_gs_coupled_majorant_identities PRIVATE milca_gen_support)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/gs_coupled_majorant_evidence)
    add_test(NAME grad_shafranov_coupled_majorant_exact_identities
        COMMAND test_gs_coupled_majorant_identities)
    set_tests_properties(grad_shafranov_coupled_majorant_exact_identities PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gs_coupled_majorant_evidence)
    set(gs_majorant_check_dir ${CMAKE_BINARY_DIR}/gs_forward_majorant_evidence)
    file(MAKE_DIRECTORY ${gs_majorant_check_dir}/native)
    add_test(NAME grad_shafranov_forward_majorant_exact_identities
        COMMAND test_gs_forward_majorant_identities)
    set_tests_properties(grad_shafranov_forward_majorant_exact_identities PROPERTIES
        WORKING_DIRECTORY ${gs_majorant_check_dir}/native)
    foreach(premise jump boundary quadrature)
        file(MAKE_DIRECTORY ${gs_majorant_check_dir}/${premise})
        add_test(NAME grad_shafranov_forward_majorant_rejects_${premise}_omission
            COMMAND test_gs_forward_majorant_identities ${premise})
        set_tests_properties(
            grad_shafranov_forward_majorant_rejects_${premise}_omission PROPERTIES
            WILL_FAIL TRUE WORKING_DIRECTORY ${gs_majorant_check_dir}/${premise})
    endforeach()
    set(gs_leaf_source ${CMAKE_CURRENT_SOURCE_DIR}/src/generated/milca_grad_shafranov_leaves.f90)
    set(gs_check_dir ${CMAKE_BINARY_DIR}/gs_variational_evidence)
    file(MAKE_DIRECTORY ${gs_check_dir}/native ${gs_check_dir}/wrong_field
        ${gs_check_dir}/wrong_weight)
    add_test(NAME grad_shafranov_variational_exact_identities
        COMMAND test_gs_variational_identities ${gs_leaf_source})
    set_tests_properties(grad_shafranov_variational_exact_identities PROPERTIES
        WORKING_DIRECTORY ${gs_check_dir}/native)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${gs_leaf_source})
    file(READ ${gs_leaf_source} gs_leaf_text)
    string(REPLACE "b_r = -psi_z/r" "b_r = psi_z/r"
        gs_wrong_field "${gs_leaf_text}")
    file(WRITE ${gs_check_dir}/wrong_field.f90 "${gs_wrong_field}")
    add_test(NAME grad_shafranov_rejects_wrong_field_orientation
        COMMAND test_gs_variational_identities ${gs_check_dir}/wrong_field.f90)
    set_tests_properties(grad_shafranov_rejects_wrong_field_orientation PROPERTIES
        WILL_FAIL TRUE WORKING_DIRECTORY ${gs_check_dir}/wrong_field)
    string(REPLACE "source = ffprime/r + mu0*pprime*r"
        "source = ffprime/r + mu0*pprime/r" gs_wrong_weight "${gs_leaf_text}")
    file(WRITE ${gs_check_dir}/wrong_weight.f90 "${gs_wrong_weight}")
    add_test(NAME grad_shafranov_rejects_wrong_weak_source_weight
        COMMAND test_gs_variational_identities ${gs_check_dir}/wrong_weight.f90)
    set_tests_properties(grad_shafranov_rejects_wrong_weak_source_weight PROPERTIES
        WILL_FAIL TRUE WORKING_DIRECTORY ${gs_check_dir}/wrong_weight)
    foreach(generator vrmhd_leaves gs_error grad_shafranov ideal_map ideal_gll ideal_xi
            ideal_response_error equilibrium_core curve_chart ideal_seed_3d map_3d_readback
            ideal_field_angle nonnested islands metriplectic vacuum_shell fold_barrier
            landreman nested_oracles
            response_strip_jets
            response_fluid
            gs_coarea gs_inverse boozer_file boozer_orbit dk_pitch dk_mono dk_hp dk_ntv dk_collop ntv_boozer
            gs_reference_error laboratory poisson fp fp_relaxation transport
            streaming streaming_moments
            debye_sphere slab_oscillation shielding_dynamics interaction cloud homma
            fermi_gas quantum_fock quantum_evolution quantum_correlations
            quantum_reference_gallery classical_pair coulomb species coulomb_encounter
            pp_fusion pp_fusion_eft
            jellium_cell ocp island_seed thermal_screening rm_oracles)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE milca_gen_support)
    endforeach()
    add_test(NAME dk_pitch_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_dk_pitch>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_dk_pitch_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/dk-pitch-regeneration
        -DLEAF=milca_dk_pitch_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(dk_pitch_leaf_regeneration fast CLASSES continuum
        TARGETS gen_dk_pitch)
    add_test(NAME dk_mono_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_dk_mono>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_dk_mono_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/dk-mono-regeneration
        -DLEAF=milca_dk_mono_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(dk_mono_leaf_regeneration fast CLASSES continuum
        TARGETS gen_dk_mono)
    add_test(NAME dk_hp_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_dk_hp>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_dk_hp_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/dk-hp-regeneration
        -DLEAF=milca_dk_hp_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(dk_hp_leaf_regeneration fast CLASSES continuum
        TARGETS gen_dk_hp)
    add_test(NAME boozer_orbit_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_boozer_orbit>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_boozer_orbit_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/boozer-orbit-regeneration
        -DLEAF=milca_boozer_orbit_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(boozer_orbit_leaf_regeneration fast CLASSES equilibrium
        TARGETS gen_boozer_orbit)
    add_test(NAME ntv_boozer_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ntv_boozer>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ntv_boozer_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ntv-boozer-regeneration
        -DLEAF=milca_ntv_boozer_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(ntv_boozer_leaf_regeneration fast CLASSES apps
        TARGETS gen_ntv_boozer)
    add_test(NAME dk_ntv_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_dk_ntv>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_dk_ntv_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/dk-ntv-regeneration
        -DLEAF=milca_dk_ntv_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(dk_ntv_leaf_regeneration fast CLASSES continuum
        TARGETS gen_dk_ntv)
    add_test(NAME dk_collop_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_dk_collop>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_dk_collop_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/dk-collop-regeneration
        -DLEAF=milca_dk_collop_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(dk_collop_leaf_regeneration fast CLASSES continuum
        TARGETS gen_dk_collop)
    add_test(NAME boozer_file_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_boozer_file>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_boozer_file_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/boozer-file-regeneration
        -DLEAF=milca_boozer_file_leaves.f90
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(boozer_file_leaf_regeneration fast CLASSES equilibrium
        TARGETS gen_boozer_file)
    add_test(NAME rm_oracle_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_rm_oracles>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_rm_oracle_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/rm-oracle-regeneration
        -P ${CMAKE_SOURCE_DIR}/test/check_rm_oracle_regeneration.cmake)
    set_tests_properties(rm_oracle_regeneration PROPERTIES LABELS rm-oracle)
    add_test(NAME vrmhd_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_vrmhd_leaves>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_vrmhd_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/vrmhd-leaf-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_vrmhd_regeneration.cmake)
    set_tests_properties(vrmhd_leaf_regeneration PROPERTIES LABELS rm1)
    add_test(NAME island_leaf_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_islands>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_island_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/island-leaf-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_island_regeneration.cmake)
    milca_test_labels(island_leaf_regeneration fast CLASSES equilibrium
        TARGETS gen_islands)
    add_executable(gen_ideal_mhd_batches ${shared_mhd_gen}/gen_ideal_mhd_batches.f90)
    target_link_libraries(gen_ideal_mhd_batches PRIVATE milca_gen_support)
    add_test(NAME ideal_mhd_batch_exact_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_mhd_batches>
        -DREVISION=${MILCA_FORTSYM_REVISION}
        -DSOURCE=${CMAKE_SOURCE_DIR}/src/generated
        -DWORK=${CMAKE_BINARY_DIR}/ideal-mhd-batch-regeneration
        -P ${shared_mhd_gen}/check_batches.cmake)
    milca_test_labels(ideal_mhd_batch_exact_regeneration fast CLASSES mhd
        TARGETS gen_ideal_mhd_batches)
    add_executable(gen_ideal_mhd_coordinates
        ${shared_mhd_gen}/gen_ideal_mhd_coordinates.f90)
    target_link_libraries(gen_ideal_mhd_coordinates PRIVATE milca_gen_support)
    add_test(NAME ideal_mhd_coordinates_exact_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_mhd_coordinates>
        -DSOURCE=${CMAKE_SOURCE_DIR}/src/generated/ideal_mhd_energy_generated.f90
        -DREVISION=${MILCA_FORTSYM_REVISION}
        -DWORK=${CMAKE_BINARY_DIR}/ideal-mhd-coordinate-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_mhd_coordinates.cmake)
    add_test(NAME ideal_xi_exact_identity_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_xi>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_xi_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-xi-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_xi_regeneration.cmake)
    add_test(NAME nonnested_exact_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_nonnested>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_nonnested_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/nonnested-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_nonnested_regeneration.cmake)
    milca_test_labels(nonnested_exact_regeneration fast CLASSES mhd TARGETS gen_nonnested)
    add_test(NAME metriplectic_identities_and_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_metriplectic>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_metriplectic_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/metriplectic-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_metriplectic_regeneration.cmake)
    milca_test_labels(metriplectic_identities_and_regeneration fast CLASSES mhd
        TARGETS gen_metriplectic)
    add_test(NAME island_seed_exact_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_island_seed>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_island_seed_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/island-seed-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_island_seed_regeneration.cmake)
    milca_test_labels(island_seed_exact_regeneration fast CLASSES mhd
        TARGETS gen_island_seed)
    add_test(NAME equilibrium_nested_oracle_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_nested_oracles>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_nested_oracle_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/nested-oracle-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_nested_oracle_regeneration.cmake)
    add_test(NAME ideal_field_angle_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_field_angle>
        -DLEAF=milca_ideal_field_angle_leaves.f90
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_field_angle_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-field-angle-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    add_test(NAME homma_leaves_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_homma>
        -DLEAF=milca_homma_leaves.f90
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_homma_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/homma-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    add_test(NAME fp_relaxation_leaves_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_fp_relaxation>
        -DLEAF=milca_fp_relaxation_leaves.f90
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_fp_relaxation_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/fp-relaxation-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    add_test(NAME thermal_screening_leaves_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_thermal_screening>
        -DLEAF=milca_thermal_screening_leaves.f90
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_thermal_screening_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/thermal-screening-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    milca_test_labels(thermal_screening_leaves_regeneration fast CLASSES continuum
        TARGETS gen_thermal_screening)
    add_test(NAME ideal_map_leaves_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_map>
        -DLEAF=milca_ideal_map_leaves.f90
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_map_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-map-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_leaf_regeneration.cmake)
    add_test(NAME ideal_seed_3d_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ideal_seed_3d>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_ideal_seed_3d_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/ideal-seed-3d-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_ideal_seed_3d_regeneration.cmake)
    add_test(NAME equilibrium_map_3d_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_map_3d_readback>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_map_3d_readback_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/map-3d-readback-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_map_3d_readback_regeneration.cmake)
    add_executable(check_response_strip_current gen/check_response_strip_current.f90)
    target_link_libraries(check_response_strip_current PRIVATE milca_gen_support)
    add_test(NAME response_strip_current_exact_identities
        COMMAND check_response_strip_current)
    add_executable(check_homma_derivations gen/check_homma_derivations.f90)
    target_link_libraries(check_homma_derivations PRIVATE milca_gen_support)
    add_test(NAME check_homma_derivations COMMAND check_homma_derivations)
    add_test(NAME equilibrium_curve_chart_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_curve_chart>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/milca_curve_chart_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/curve-chart-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_check_curve_chart_regeneration.cmake)
endif()
# Explicit performance audit; not part of the ordinary correctness gate.
add_executable(benchmark_streaming app/benchmark_streaming.f90)
target_link_libraries(benchmark_streaming PRIVATE milca)
add_executable(benchmark_streaming_comparison app/benchmark_streaming_comparison.f90)
target_link_libraries(benchmark_streaming_comparison PRIVATE milca)
add_executable(homma_reference app/homma_reference.f90)
add_executable(benchmark_fp app/benchmark_fp.f90)
target_link_libraries(benchmark_fp PRIVATE milca)
add_executable(fp_relaxation app/fp_relaxation.f90)
target_link_libraries(fp_relaxation PRIVATE milca milca_io)
add_executable(screening_moments app/screening_moments.f90)
target_link_libraries(screening_moments PRIVATE milca milca_io)
add_executable(benchmark_transport app/benchmark_transport.f90)
target_link_libraries(benchmark_transport PRIVATE milca)
add_executable(benchmark_slab_oscillation app/benchmark_slab_oscillation.f90)
target_link_libraries(benchmark_slab_oscillation PRIVATE milca)

add_executable(benchmark_fermi_gas app/benchmark_fermi_gas.f90)
target_link_libraries(benchmark_fermi_gas PRIVATE milca)
add_executable(benchmark_quantum_cutoff app/benchmark_quantum_cutoff.f90)
target_link_libraries(benchmark_quantum_cutoff PRIVATE milca)
add_executable(benchmark_ocp app/benchmark_ocp.f90)
target_link_libraries(benchmark_ocp PRIVATE milca)
# The periodic OCP demos run independent ensemble members on OpenMP threads;
# the library is reentrant (no module state) and compiled with automatic
# locals (-frecursive) so that concurrent calls do not share scratch.
foreach(example ocp_oscillation ocp_shielding)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE milca fortplot)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(${example}_demo PRIVATE OpenMP::OpenMP_Fortran)
    endif()
endforeach()

# Campaign Q3/Q8: finite-N ideal fermion ensembles and relativistic Fermi EOS.
if(MILCA_FORTSYM_BUILD)
    foreach(generator relativistic_fermi fermion_ensemble)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE milca_gen_support)
    endforeach()
endif()
target_sources(milca_quantum PRIVATE src/generated/milca_relativistic_fermi_leaves.f90)
set_source_files_properties(src/generated/milca_relativistic_fermi_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_relativistic_fermi_gas test/test_relativistic_fermi_gas.f90)
target_link_libraries(test_relativistic_fermi_gas PRIVATE milca)
add_test(NAME relativistic_fermi_gas_oracles COMMAND test_relativistic_fermi_gas)
add_executable(test_fermi_gas_cold test/test_fermi_gas_cold.f90)
target_link_libraries(test_fermi_gas_cold PRIVATE milca)
add_test(NAME fermi_gas_cold_sommerfeld COMMAND test_fermi_gas_cold)
target_sources(milca_quantum PRIVATE src/generated/milca_fermion_ensemble_leaves.f90
    src/quantum/fermion_ensemble.f90)
set_source_files_properties(src/generated/milca_fermion_ensemble_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_fermion_ensemble test/test_fermion_ensemble.f90)
target_link_libraries(test_fermion_ensemble PRIVATE milca)
add_test(NAME fermion_ensemble_oracles COMMAND test_fermion_ensemble)
add_executable(fermion_ensemble_demo app/fermion_ensemble_demo.f90)
target_link_libraries(fermion_ensemble_demo PRIVATE milca fortplot)
add_executable(relativistic_fermi_demo app/relativistic_fermi_demo.f90)
target_link_libraries(relativistic_fermi_demo PRIVATE milca fortplot)

# Campaign F5: namelist configuration in, NetCDF-4/HDF5 result with embedded
# literal/resolved configuration and provenance out (src/apps/result_io.f90).
# Physics stays in milca; only milca_io and applications see FortIO.
set(milca_build_info_inc ${CMAKE_BINARY_DIR}/result_include/milca_build_info.inc)
string(TOUPPER "${CMAKE_BUILD_TYPE}" milca_build_type_upper)
add_custom_target(milca_build_info
    COMMAND ${CMAKE_COMMAND} -DMILCA_SOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DMILCA_INFO_OUTPUT=${milca_build_info_inc}
        -DMILCA_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        "-DMILCA_BUILD_FLAGS=${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${milca_build_type_upper}}"
        "-DMILCA_BUILD_COMPILER=${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}"
        -DMILCA_FORTSYM_REVISION=${MILCA_FORTSYM_REVISION}
        -DMILCA_FORTNUM_REVISION=${MILCA_FORTNUM_REVISION}
        -DMILCA_FORTPLOT_REVISION=${MILCA_FORTPLOT_REVISION}
        -DMILCA_FORTIO_REVISION=${MILCA_FORTIO_REVISION}
        -DMILCA_FORTFEM_REVISION=${MILCA_FORTFEM_REVISION}
        -DMILCA_FORTSPARSE_REVISION=${MILCA_FORTSPARSE_REVISION}
        -DMILCA_LIBNEO_REVISION=${MILCA_LIBNEO_REVISION}
        -P ${CMAKE_SOURCE_DIR}/cmake/milca_build_info.cmake
    BYPRODUCTS ${milca_build_info_inc} VERBATIM)
add_library(milca_io src/apps/result_io.f90 src/apps/periodic_pair_config.f90
    src/apps/animation_export.f90)
add_dependencies(milca_io milca_build_info)
# NTV field lines from Boozer files: composes the equilibrium reader with the
# continuum kernel, so it lives in the application layer.
add_library(milca_ntv_geometry src/apps/ntv_boozer_line.f90 src/apps/ntv_orbit_geometry.f90
    src/generated/milca_ntv_boozer_leaves.f90 src/apps/ntv_device.f90)
target_link_libraries(milca_ntv_geometry PUBLIC milca)
if(MILCA_ENABLE_OPENMP)
    target_link_libraries(milca_ntv_geometry PUBLIC OpenMP::OpenMP_Fortran)
endif()
add_executable(ntv_profile app/ntv_profile.f90)
target_link_libraries(ntv_profile PRIVATE milca_ntv_geometry)
add_executable(ntv_axisymmetric_profile app/ntv_axisymmetric_profile.f90)
target_link_libraries(ntv_axisymmetric_profile PRIVATE milca_ntv_geometry)
add_executable(ntv_indicator_profile app/ntv_indicator_profile.f90)
target_link_libraries(ntv_indicator_profile PRIVATE milca_ntv_geometry)
add_executable(test_ntv_orbit_geometry test/test_ntv_orbit_geometry.f90)
target_link_libraries(test_ntv_orbit_geometry PRIVATE milca_ntv_geometry)
add_test(NAME ntv_orbit_geometry_thin COMMAND test_ntv_orbit_geometry)
milca_test_labels(ntv_orbit_geometry_thin fast CLASSES apps TARGETS test_ntv_orbit_geometry)
add_executable(test_ntv_device test/test_ntv_device.f90)
target_link_libraries(test_ntv_device PRIVATE milca_ntv_geometry)
add_test(NAME ntv_device COMMAND test_ntv_device)
milca_test_labels(ntv_device fast CLASSES apps TARGETS test_ntv_device)
add_executable(test_ntv_boozer_line test/test_ntv_boozer_line.f90)
target_link_libraries(test_ntv_boozer_line PRIVATE milca_ntv_geometry)
add_test(NAME ntv_boozer_line COMMAND test_ntv_boozer_line)
milca_test_labels(ntv_boozer_line fast CLASSES apps TARGETS test_ntv_boozer_line)
target_include_directories(milca_io PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(milca_io PUBLIC milca fortio::fortio)
target_compile_options(milca_io PRIVATE -Wall -Wextra)
# Every application reads --config FILE and writes output/<name>.nc.
foreach(application milca_demo steady_laboratory_demo transport_demo
        streaming_comparison_demo streaming_moments_demo debye_sphere_demo
        slab_oscillation_demo shielding_dynamics_demo interaction_demo
        microscopic_cloud_demo quantum_reference_demo classical_pair_demo
        fermion_ensemble_demo relativistic_fermi_demo formal_krook_certificate
        benchmark_streaming benchmark_streaming_comparison benchmark_fp
        benchmark_transport benchmark_slab_oscillation benchmark_fermi_gas
        benchmark_quantum_cutoff jellium_cell_demo coulomb_encounter_demo fp_relaxation
        benchmark_ocp ocp_oscillation_demo ocp_shielding_demo homma_reference
        screening_moments)
    target_link_libraries(${application} PRIVATE milca_io)
endforeach()
add_executable(test_result_io test/test_result_io.f90)
target_link_libraries(test_result_io PRIVATE milca_io)
add_test(NAME result_io_contract
    COMMAND ${CMAKE_COMMAND} -DPROBE=$<TARGET_FILE:test_result_io>
        -DCLASSICAL=$<TARGET_FILE:classical_pair_demo>
        -DQUANTUM=$<TARGET_FILE:quantum_reference_demo>
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/result_io_check
        -DFORTIO_PIN=${MILCA_FORTIO_REVISION}
        -DFORTSPARSE_PIN=${MILCA_FORTSPARSE_REVISION}
        -DLIBNEO_PIN=${MILCA_LIBNEO_REVISION}
        -P ${CMAKE_SOURCE_DIR}/test/result_io_check.cmake)
add_executable(test_animation_export test/test_animation_export.f90)
target_link_libraries(test_animation_export PRIVATE milca_io)
add_test(NAME animation_export_contract COMMAND test_animation_export
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# Q7: static screening hierarchy and 1X1V Wigner-Poisson (appended block).
add_library(milca_fft
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft8_kernel.f90
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft.f90)
target_link_libraries(milca_fft PUBLIC milca_ode)
target_compile_options(milca_fft PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(milca_quantum PRIVATE
    src/generated/milca_screening_leaves.f90 src/quantum/screening_hierarchy.f90)
# Wigner-Poisson is a phase-space continuum solver sharing the streaming leaves.
target_sources(milca_continuum PRIVATE
    src/generated/milca_wigner_leaves.f90 src/continuum/wigner_dispersion.f90
    src/continuum/wigner_poisson.f90)
set_source_files_properties(src/generated/milca_screening_leaves.f90
    src/generated/milca_wigner_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(milca_core PUBLIC milca_fft)
foreach(q7 screening wigner)
    add_executable(${q7}_demo app/${q7}_demo.f90)
    target_link_libraries(${q7}_demo PRIVATE milca milca_io fortplot)
    add_executable(test_${q7} test/test_${q7}.f90)
    target_link_libraries(test_${q7} PRIVATE milca)
    add_test(NAME ${q7}_q7_oracles COMMAND test_${q7})
    if(MILCA_FORTSYM_BUILD)
        add_executable(gen_${q7} gen/gen_${q7}.f90)
        target_link_libraries(gen_${q7} PRIVATE milca_gen_support)
    endif()
endforeach()
# Q5: Hartree/HF/Mermin/TDHF mean-field models of the Q4 jellium cell (appended).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_mean_field_cell gen/gen_mean_field_cell.f90)
    target_link_libraries(gen_mean_field_cell PRIVATE milca_gen_support)
endif()
target_sources(milca_quantum PRIVATE src/generated/milca_mean_field_cell_leaves.f90
    src/quantum/mean_field_cell.f90)
set_source_files_properties(src/generated/milca_mean_field_cell_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_mean_field_cell test/test_mean_field_cell.f90)
target_link_libraries(test_mean_field_cell PRIVATE milca)
add_test(NAME mean_field_cell_oracles COMMAND test_mean_field_cell)
add_executable(mean_field_cell_demo app/mean_field_cell_demo.f90)
target_link_libraries(mean_field_cell_demo PRIVATE milca milca_io fortplot)
# QG6: quantum binary-collision kernel and hybrid classical/quantum encounters
# (appended block).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_quantum_collision gen/gen_quantum_collision.f90)
    target_link_libraries(gen_quantum_collision PRIVATE milca_gen_support)
endif()
target_sources(milca_pairs PRIVATE src/generated/milca_quantum_collision_leaves.f90
    src/pairs/quantum_collision_kernel.f90)
set_source_files_properties(src/generated/milca_quantum_collision_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_sources(milca_quantum PRIVATE src/quantum/hybrid_encounter.f90)
add_executable(test_quantum_collision test/test_quantum_collision.f90)
target_link_libraries(test_quantum_collision PRIVATE milca)
add_test(NAME quantum_collision_qg6_oracles COMMAND test_quantum_collision)
add_executable(test_hybrid_encounter test/test_hybrid_encounter.f90)
target_link_libraries(test_hybrid_encounter PRIVATE milca)
add_test(NAME hybrid_encounter_qg6_oracles COMMAND test_hybrid_encounter)
add_executable(hybrid_encounter_demo app/hybrid_encounter_demo.f90)
target_link_libraries(hybrid_encounter_demo PRIVATE milca milca_io fortplot)

# Campaign Q2: Gaussian/coherent-sum packet models measured against the Q1
# partial-wave reference (src/quantum/packet_models.f90, docs/packet-models.md).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_packet_models gen/gen_packet_models.f90)
    target_link_libraries(gen_packet_models PRIVATE milca_gen_support)
endif()
target_sources(milca_quantum PRIVATE src/generated/milca_packet_models_leaves.f90
    src/quantum/gaussian_pair.f90 src/quantum/packet_models.f90)
set_source_files_properties(src/generated/milca_packet_models_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_packet_models test/test_packet_models.f90)
target_link_libraries(test_packet_models PRIVATE milca)
add_test(NAME packet_models_oracles COMMAND test_packet_models)
add_executable(packet_models_demo app/packet_models_demo.f90)
target_link_libraries(packet_models_demo PRIVATE milca fortplot milca_io)

# Campaign H1b: quantum vs classical finite-N pair/three-body correlations
# (docs/quantum-classical-correlations.md).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_quantum_classical gen/gen_quantum_classical.f90)
    target_link_libraries(gen_quantum_classical PRIVATE milca_gen_support)
endif()
target_sources(milca_pairs PRIVATE src/generated/milca_quantum_classical_leaves.f90
    src/pairs/pair_law.f90)
target_sources(milca_classical PRIVATE src/classical/classical_cell.f90)
target_sources(milca_quantum PRIVATE src/quantum/quantum_pair_cell.f90
    src/quantum/quantum_classical_map.f90)
set_source_files_properties(src/generated/milca_quantum_classical_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_quantum_classical test/test_quantum_classical.f90)
target_link_libraries(test_quantum_classical PRIVATE milca)
add_test(NAME quantum_classical_h1b_oracles COMMAND test_quantum_classical)
add_executable(correlation_map_quantum_demo app/correlation_map_quantum_demo.f90)
target_link_libraries(correlation_map_quantum_demo PRIVATE milca milca_io
    milca_demo_style fortplot)
# Campaign H1a: classical OCP correlation-hierarchy map (Metropolis MC with
# incremental Ewald, g3/cumulant estimators, YBG closures, HNC; appended).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_correlation_hierarchy gen/gen_correlation_hierarchy.f90)
    target_link_libraries(gen_correlation_hierarchy PRIVATE milca_gen_support)
endif()
target_sources(milca_core PRIVATE src/core/series_statistics.f90)
target_sources(milca_classical PRIVATE src/generated/milca_correlation_hierarchy_leaves.f90
    src/classical/ocp_monte_carlo.f90 src/classical/correlation_hierarchy.f90 src/classical/ocp_hnc.f90)
set_source_files_properties(src/generated/milca_correlation_hierarchy_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_correlation_hierarchy test/test_correlation_hierarchy.f90)
target_link_libraries(test_correlation_hierarchy PRIVATE milca)
add_test(NAME correlation_hierarchy_h1a_oracles COMMAND test_correlation_hierarchy)
add_executable(correlation_map_classical_demo app/correlation_map_types.f90
    app/correlation_map_md.f90 app/correlation_map_plots.f90
    app/correlation_map_classical_demo.f90)
target_link_libraries(correlation_map_classical_demo PRIVATE milca milca_io fortplot
    milca_demo_style)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(correlation_map_classical_demo PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Error-control surveys are explicit performance runs, never timing tests.
add_library(error_benchmark_support app/error_benchmark_support.f90)
foreach(survey continuum quantum microscopic modes)
    add_executable(benchmark_error_${survey} app/benchmark_error_${survey}.f90)
    target_link_libraries(benchmark_error_${survey} PRIVATE milca error_benchmark_support)
endforeach()
add_executable(test_error_modes test/test_error_modes.f90)
target_link_libraries(test_error_modes PRIVATE milca)
add_test(NAME error_modes_oracles COMMAND test_error_modes)
foreach(fixture verified_matvec sampled_norm)
    add_executable(test_error_${fixture} ${fortnum_SOURCE_DIR}/test/test_${fixture}.f90)
    target_link_libraries(test_error_${fixture} PRIVATE milca_ode)
    add_test(NAME error_${fixture}_oracles COMMAND test_error_${fixture})
endforeach()

add_executable(test_curved_triangle test/test_curved_triangle.f90)
target_link_libraries(test_curved_triangle PRIVATE milca)
add_test(NAME grad_shafranov_curved_geometry COMMAND test_curved_triangle)
add_executable(test_gs_curved_convergence test/test_gs_curved_convergence.f90)
target_link_libraries(test_gs_curved_convergence PRIVATE milca)
add_test(NAME grad_shafranov_curved_convergence COMMAND test_gs_curved_convergence)
add_test(NAME grad_shafranov_p3_curved_convergence COMMAND test_gs_curved_convergence 3)
add_executable(test_gs_p3 test/test_gs_p3.f90)
add_executable(test_gs_nonlinear test/test_gs_nonlinear.f90)
add_executable(test_gs_profile_newton test/test_gs_profile_newton.f90)
target_link_libraries(test_gs_profile_newton PRIVATE milca)
add_test(NAME grad_shafranov_cubic_profile_oracles COMMAND test_gs_profile_newton)
add_test(NAME grad_shafranov_profile_quadrature_work COMMAND test_gs_profile_newton work)
set_tests_properties(grad_shafranov_cubic_profile_oracles
    grad_shafranov_profile_quadrature_work PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_link_libraries(test_gs_nonlinear PRIVATE milca)
add_test(NAME grad_shafranov_p2_nonlinear_convergence COMMAND test_gs_nonlinear 2)
add_test(NAME grad_shafranov_p3_nonlinear_convergence COMMAND test_gs_nonlinear 3)
add_test(NAME grad_shafranov_p3_spline_knot_convergence COMMAND test_gs_nonlinear 3 knots)
set_tests_properties(grad_shafranov_p2_nonlinear_convergence
    grad_shafranov_p3_nonlinear_convergence
    grad_shafranov_p3_spline_knot_convergence PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_link_libraries(test_gs_p3 PRIVATE milca)
add_test(NAME grad_shafranov_p3_physical_oracles COMMAND test_gs_p3)
add_executable(test_gs_threads test/test_gs_threads.f90)
target_link_libraries(test_gs_threads PRIVATE milca)
if(OpenMP_Fortran_FOUND)
    target_compile_options(test_gs_threads PRIVATE ${OpenMP_Fortran_FLAGS})
endif()
add_test(NAME grad_shafranov_p3_thread_invariance COMMAND test_gs_threads)
set_tests_properties(grad_shafranov_p3_thread_invariance PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_native_adapter test/test_gs_native_adapter.f90)
target_link_libraries(test_gs_native_adapter PRIVATE milca_fem2d)
add_test(NAME gs_native_adapter_contract COMMAND test_gs_native_adapter)
set_tests_properties(gs_native_adapter_contract PROPERTIES
    TIMEOUT 5 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_direct test/test_gs_direct.f90)
target_link_libraries(test_gs_direct PRIVATE milca)
add_test(NAME grad_shafranov_csr_direct_newton_oracles COMMAND test_gs_direct)
set_tests_properties(grad_shafranov_csr_direct_newton_oracles PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_device_selection test/test_gs_device_selection.f90)
target_link_libraries(test_gs_device_selection PRIVATE milca_equilibrium)
add_test(NAME grad_shafranov_device_selection_contract
    COMMAND test_gs_device_selection)
set_tests_properties(grad_shafranov_device_selection_contract PROPERTIES
    TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
# Affine-profile Picard fixed point through the public GS route: the CPU
# route always runs; the CUDA guard below adds the CPU/GPU pair.
add_executable(test_gs_cuda_picard test/test_gs_cuda_picard.f90)
target_link_libraries(test_gs_cuda_picard PRIVATE milca_equilibrium)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_gs_cuda_picard PRIVATE -fno-fast-math -ffp-contract=off)
endif()
add_test(NAME grad_shafranov_affine_picard_cpu COMMAND test_gs_cuda_picard cpu)
set_tests_properties(grad_shafranov_affine_picard_cpu PROPERTIES TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
if(MILCA_ENABLE_CUDA_GS_SOLVE)
    add_executable(test_gs_cuda_fixed_source test/test_gs_cuda_fixed_source.f90)
    target_link_libraries(test_gs_cuda_fixed_source PRIVATE milca_nn_solovev_test)
    add_executable(test_gs_cuda_public_solve test/test_gs_cuda_public_solve.f90)
    target_link_libraries(test_gs_cuda_public_solve PRIVATE milca_nn_solovev_test)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(test_gs_cuda_fixed_source PRIVATE
            -fno-fast-math -ffp-contract=off)
        target_compile_options(test_gs_cuda_public_solve PRIVATE
            -fno-fast-math -ffp-contract=off)
    endif()
    add_test(NAME grad_shafranov_cuda_fixed_source_a3
        COMMAND test_gs_cuda_fixed_source a3 16 0)
    add_test(NAME grad_shafranov_cuda_affine_strict_audit
        COMMAND test_gs_cuda_fixed_source affine 16 0)
    add_test(NAME grad_shafranov_cuda_public_fixed_source
        COMMAND test_gs_cuda_public_solve)
    add_test(NAME grad_shafranov_cuda_affine_picard
        COMMAND test_gs_cuda_picard pair)
    set_tests_properties(grad_shafranov_cuda_fixed_source_a3
        grad_shafranov_cuda_affine_strict_audit
        grad_shafranov_cuda_public_fixed_source
        grad_shafranov_cuda_affine_picard PROPERTIES
        TIMEOUT 120 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endif()
add_executable(test_gs_solver_work test/test_gs_solver_work.f90)
target_link_libraries(test_gs_solver_work PRIVATE milca)
add_test(NAME grad_shafranov_shaped_solver_work COMMAND test_gs_solver_work)
add_executable(test_gs_axis_refinement test/test_gs_axis_refinement.f90)
target_link_libraries(test_gs_axis_refinement PRIVATE milca)
add_test(NAME grad_shafranov_axis_refinement COMMAND test_gs_axis_refinement)
add_executable(test_gs_profile_law_readback test/test_gs_synthetic_law_readback.f90)
target_link_libraries(test_gs_profile_law_readback PRIVATE milca)
add_test(NAME grad_shafranov_profile_law_readback
    COMMAND test_gs_profile_law_readback)
set_tests_properties(grad_shafranov_curved_convergence
    grad_shafranov_p3_curved_convergence grad_shafranov_axis_refinement
    grad_shafranov_shaped_solver_work PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_sources(milca_equilibrium PRIVATE src/equilibrium/gs_readback.f90
    src/equilibrium/gs_axis_recovery.f90 src/equilibrium/gs_field_readback.f90
    src/equilibrium/gs_ampere_loop.f90)
add_executable(test_gs_curved_readback test/test_gs_curved_readback.f90)
target_link_libraries(test_gs_curved_readback PRIVATE milca)
add_test(NAME grad_shafranov_curved_readback COMMAND test_gs_curved_readback)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/p3_readback_test)
add_test(NAME grad_shafranov_p3_curved_readback COMMAND test_gs_curved_readback 3)
set_tests_properties(grad_shafranov_p3_curved_readback PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/p3_readback_test
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
set_tests_properties(grad_shafranov_curved_readback PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_ampere_loop test/test_gs_ampere_loop.f90)
target_link_libraries(test_gs_ampere_loop PRIVATE milca)
foreach(degree 2 3)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/ampere_p${degree}_test)
    add_test(NAME grad_shafranov_p${degree}_ampere_loop
        COMMAND test_gs_ampere_loop ${degree})
    set_tests_properties(grad_shafranov_p${degree}_ampere_loop PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/ampere_p${degree}_test
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_gs_weak_ampere_loop test/test_gs_weak_ampere_loop.f90)
target_link_libraries(test_gs_weak_ampere_loop PRIVATE milca)
foreach(degree 2 3)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/weak_ampere_p${degree}_test)
    add_test(NAME grad_shafranov_p${degree}_weak_ampere_loop
        COMMAND test_gs_weak_ampere_loop ${degree})
    set_tests_properties(grad_shafranov_p${degree}_weak_ampere_loop PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/weak_ampere_p${degree}_test
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(gs_phase1 app/gs_phase1.f90)
target_link_libraries(gs_phase1 PRIVATE milca_io fortplot)
add_executable(gs_load_bench app/gs_load_bench.f90)
target_link_libraries(gs_load_bench PRIVATE milca)
add_executable(gs_free app/gs_free.f90)
target_link_libraries(gs_free PRIVATE milca_io fortplot)
add_executable(gs_free_3d app/gs_free_3d.f90)
target_link_libraries(gs_free_3d PRIVATE milca)
# NL-3a on the gs_free fixture (N = 0): reproduce the 2D free boundary and
# return from +-15 mm R/Z and +2 % elongation offset starts.
# One registered test per start case (see milca_test_labels.cmake).
foreach(case r0 rp rm zp zm el gzm)
    add_test(NAME gs_free_3d_nl3a_fixture_return_${case} COMMAND ${CMAKE_COMMAND}
        -DAPP=$<TARGET_FILE:gs_free_3d> -DGS_FREE=$<TARGET_FILE:gs_free>
        -DFIXTURES=${CMAKE_SOURCE_DIR}/test/fixtures/gs_free
        -DWORK=${CMAKE_BINARY_DIR}/gs_free_3d_return_${case} -DCASE=${case}
        -DDISTANCE_LIMIT=1.0e-4 -DDEFECT_LIMIT=1.0e-4
        -P ${CMAKE_SOURCE_DIR}/test/gs_free_3d_return.cmake)
    set_tests_properties(gs_free_3d_nl3a_fixture_return_${case} PROPERTIES
        TIMEOUT 3600 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
# Discrete hypotheses (FB-A2, FB-D2) at converged fixture states: A_h + S_h > 0,
# [grad C; M] of full row rank and one negative direction of S_red (vertical).
add_test(NAME gs_free_3d_fixture_hypotheses COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_free_3d> -DGS_FREE=$<TARGET_FILE:gs_free>
    -DFIXTURES=${CMAKE_SOURCE_DIR}/test/fixtures/gs_free -DMODE=hypotheses
    -DWORK=${CMAKE_BINARY_DIR}/gs_free_3d_hypotheses -DDISTANCE_LIMIT=1.0e-4
    -DDEFECT_LIMIT=1.0e-4 -P ${CMAKE_SOURCE_DIR}/test/gs_free_3d_return.cmake)
set_tests_properties(gs_free_3d_fixture_hypotheses PROPERTIES TIMEOUT 3600
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_test(NAME gs_free_app_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_free> -DFIXTURES=${CMAKE_SOURCE_DIR}/test/fixtures/gs_free
    -DWORK=${CMAKE_BINARY_DIR}/gs_free_app
    -P ${CMAKE_SOURCE_DIR}/test/gs_free_app.cmake)
set_tests_properties(gs_free_app_contract PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
# The comparison producer is optional; the native solver needs no Python.

add_executable(test_gs_readback_order test/test_gs_readback_order.f90)
target_link_libraries(test_gs_readback_order PRIVATE milca)
add_test(NAME grad_shafranov_p2_readback_order COMMAND test_gs_readback_order 2)
add_test(NAME grad_shafranov_p3_readback_order COMMAND test_gs_readback_order 3)

add_executable(test_gs_readback_free test/test_gs_readback_free.f90)
target_link_libraries(test_gs_readback_free PRIVATE milca)
add_test(NAME grad_shafranov_p2_readback_free COMMAND test_gs_readback_free 2)
add_test(NAME grad_shafranov_p3_readback_free COMMAND test_gs_readback_free 3)
set_tests_properties(grad_shafranov_p2_readback_free
    grad_shafranov_p3_readback_free PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
    TIMEOUT 300)

add_executable(test_gs_readback_work test/test_gs_readback_work.f90)
target_link_libraries(test_gs_readback_work PRIVATE milca)
add_test(NAME grad_shafranov_readback_work COMMAND test_gs_readback_work)
set_tests_properties(grad_shafranov_p2_readback_order
    grad_shafranov_p3_readback_order grad_shafranov_readback_work PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
    TIMEOUT 300)
add_test(NAME grad_shafranov_stream_readback
    COMMAND ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:gs_phase1>
    -DWORK=${CMAKE_BINARY_DIR}/gs_stream_test
    -P ${CMAKE_SOURCE_DIR}/test/gs_phase1_stream.cmake)
set_tests_properties(grad_shafranov_stream_readback PROPERTIES TIMEOUT 30)

add_test(NAME grad_shafranov_production_estimator
    COMMAND ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:gs_phase1>
    -DWORK=${CMAKE_BINARY_DIR}/gs_estimator_test
    -P ${CMAKE_SOURCE_DIR}/test/gs_phase1_estimator.cmake)
set_tests_properties(grad_shafranov_production_estimator PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")

# Common equilibrium core (docs/equilibrium-core.md): state, chart, adapters.
target_sources(milca_mhd PRIVATE src/generated/ideal_mhd_energy_generated.f90
    src/mhd/ideal_xi.F90
    src/generated/milca_ideal_xi_leaves.f90)
# Generated derivatives retain the primal ABI for constant partials.
set_source_files_properties(src/generated/milca_ideal_xi_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_ideal_mhd_coordinates test/test_ideal_mhd_coordinates.f90)
target_link_libraries(test_ideal_mhd_coordinates PRIVATE milca_mhd)
add_test(NAME ideal_mhd_coordinate_analytic_limits COMMAND test_ideal_mhd_coordinates)
add_executable(test_ideal_xi test/test_ideal_xi.f90)
target_link_libraries(test_ideal_xi PRIVATE milca_mhd)
add_test(NAME ideal_xi_hessian_and_analytic_limits COMMAND test_ideal_xi)
set_tests_properties(ideal_xi_hessian_and_analytic_limits PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2" PROCESSORS 2)
target_sources(milca_mhd PRIVATE src/mhd/map_equilibrium.f90)
target_sources(milca_equilibrium PRIVATE src/generated/milca_equilibrium_core_leaves.f90
    src/generated/milca_curve_chart_leaves.f90
    src/equilibrium/equilibrium_chart.f90 src/equilibrium/equilibrium_state.f90 src/equilibrium/gs_equilibrium.f90
    src/equilibrium/equilibrium_residual.f90 src/equilibrium/equilibrium_readback.f90
    src/equilibrium/equilibrium_certificate.f90 src/equilibrium/gs_certificate.f90
    src/equilibrium/equilibrium_boundary.f90 src/equilibrium/equilibrium_geometry.f90)
set_source_files_properties(src/generated/milca_equilibrium_core_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_equilibrium_boundary test/test_equilibrium_boundary.f90)
target_link_libraries(test_equilibrium_boundary PRIVATE milca)
add_test(NAME equilibrium_core_d6_exterior_per_n
    COMMAND test_equilibrium_boundary exterior)
set_tests_properties(equilibrium_core_d6_exterior_per_n PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=4" PROCESSORS 4 TIMEOUT 600)
add_test(NAME equilibrium_core_d6_coil_interface
    COMMAND test_equilibrium_boundary coils)
add_executable(test_equilibrium_geometry test/test_equilibrium_geometry.f90)
target_link_libraries(test_equilibrium_geometry PRIVATE milca)
foreach(mode solovev circle drift state)
    add_test(NAME equilibrium_core_d7_geometry_${mode}
        COMMAND test_equilibrium_geometry ${mode})
    set_tests_properties(equilibrium_core_d7_geometry_${mode} PROPERTIES
        TIMEOUT 300 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_equilibrium_core test/test_equilibrium_core.f90)
target_link_libraries(test_equilibrium_core PRIVATE milca)
add_test(NAME equilibrium_core_chart_jet COMMAND test_equilibrium_core chart)
add_executable(test_curve_chart test/test_curve_chart.f90)
set_source_files_properties(test/test_curve_chart.f90 PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(test_curve_chart PRIVATE milca_equilibrium)
foreach(mode circle jets tube periodic invalid)
    add_test(NAME equilibrium_curve_chart_${mode} COMMAND test_curve_chart ${mode})
endforeach()
add_executable(test_curve_chart_solver test/test_curve_chart_solver.f90)
target_link_libraries(test_curve_chart_solver PRIVATE milca_mhd)
foreach(mode landreman figure8)
    add_test(NAME ideal_map_curve_chart_${mode} COMMAND test_curve_chart_solver ${mode})
    set_tests_properties(ideal_map_curve_chart_${mode} PROPERTIES TIMEOUT 120
        PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(k 2 3)
    add_test(NAME equilibrium_core_d7_readback_p${k}
        COMMAND test_equilibrium_core readback ${k})
    add_test(NAME equilibrium_core_d3_solovev_p${k}
        COMMAND test_equilibrium_core residual ${k})
    set_tests_properties(equilibrium_core_d3_solovev_p${k} PROPERTIES TIMEOUT 300
        PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_test(NAME equilibrium_core_d3_upper_bound_p2
    COMMAND test_equilibrium_core residual_bound)
set_tests_properties(equilibrium_core_d3_upper_bound_p2 PROPERTIES TIMEOUT 60
    PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
# Synthetic Solov'ev inputs exercise the cubic-profile Newton, P3, direct
# factor, current constraint, readback and majorant path on a coarse mesh.
add_test(NAME gs_k6_mini_regression COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_phase1> -DCHECK=$<TARGET_FILE:test_equilibrium_core>
    -DFIXTURES=${CMAKE_SOURCE_DIR}/test/fixtures/gs_synthetic
    -DWORK=${CMAKE_BINARY_DIR}/gs_k6_mini
    -P ${CMAKE_SOURCE_DIR}/test/gs_k6_mini.cmake)
set_tests_properties(gs_k6_mini_regression PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

# K-C1/L3: 3D polyline coils, per-n Fourier projection, vacuum b_mn (L3).
target_sources(milca_equilibrium PRIVATE src/equilibrium/coils_3d.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/makegrid_circuits.f90)
target_link_libraries(milca_equilibrium PUBLIC milca_libneo_coils)
target_sources(milca_mhd PRIVATE src/mhd/vacuum_surfaces.f90)
add_executable(test_coils_3d test/test_coils_3d.f90)
target_link_libraries(test_coils_3d PRIVATE milca)
add_test(NAME coils_3d_biot_savart_fourier_oracles COMMAND test_coils_3d)
set_tests_properties(coils_3d_biot_savart_fourier_oracles PROPERTIES TIMEOUT 120)
add_executable(test_vacuum_surfaces test/test_vacuum_surfaces.f90)
target_link_libraries(test_vacuum_surfaces PRIVATE milca)
add_test(NAME vacuum_surfaces_sfl_harmonics_oracles COMMAND test_vacuum_surfaces)
set_tests_properties(vacuum_surfaces_sfl_harmonics_oracles PROPERTIES TIMEOUT 300)
add_executable(coil_modes app/coil_modes.f90)
target_link_libraries(coil_modes PRIVATE milca)
add_executable(gs_vacuum_3d app/gs_vacuum_3d.f90)
target_link_libraries(gs_vacuum_3d PRIVATE milca)

# L4 (first stage): per-n FEEC complex, toroidal-mode BEM and their magnetic
# coupling from FortFEM; flux-aligned meshes and the per-n Ampere operator
# with the exact exterior.
add_library(milca_feec_n
    ${fortfem_SOURCE_DIR}/src/operators/direct_fourier_transform.f90
    ${fortfem_SOURCE_DIR}/src/boundary/boundary_types.f90
    ${fortfem_SOURCE_DIR}/src/mesh/triangle_io.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/bowyer_watson.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/constrained_delaunay.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/mesh_refinement.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/triangulation_fortran.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/triangulator.f90
    ${fortfem_SOURCE_DIR}/src/mesh/mesh_2d.f90
    ${fortfem_SOURCE_DIR}/src/geometry/meridian_curved_mesh.f90
    ${fortfem_SOURCE_DIR}/src/elements/polynomial_candidate_jets.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_lagrange_arbitrary_order.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_triangle_feec_tables.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_triangle_second_gradient.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_triangle_bernstein_jets.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_bernstein_jets.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_nedelec_arbitrary_order.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_rt_arbitrary_order.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_nedelec_second_kind.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_bdm_arbitrary_order.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_affine_map_products.f90
    ${fortfem_SOURCE_DIR}/src/geometry/triangle_affine_map.f90
    ${fortfem_SOURCE_DIR}/src/interpolation/triangle_vector_interpolation.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_triangle_piola_products.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_piola_maps.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_mapped_oscillation_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_equilibrated_flux_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_rt0_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/equilibrated_flux_topology.f90
    ${fortfem_SOURCE_DIR}/src/estimators/equilibrated_flux_2d.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_meridian_hcurl_equilibration_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/equilibrated_hcurl_meridian_n.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_rt0_majorant.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_bernstein_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_bernstein_majorant.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_feec_operators.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_feec_n.f90
    ${fortfem_SOURCE_DIR}/src/solvers/hcurl_auxiliary_space.f90
    ${fortfem_SOURCE_DIR}/src/assembly/meridian_hcurl_auxiliary.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_scalar.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_meridian_fourier_equilibration_products.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_equilibration.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_curl.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_apply.f90
    ${fortfem_SOURCE_DIR}/src/bem/toroidal_mode_kernel.f90
    ${fortfem_SOURCE_DIR}/src/bem/toroidal_mode_boundary_operators_2d.f90
    ${fortfem_SOURCE_DIR}/src/bem/meridian_magnetic_coupling.f90)
target_link_libraries(milca_feec_n PUBLIC
    milca_spatial milca_fft milca_sparse_direct fortsparse LAPACK::LAPACK)
if(MILCA_ENABLE_OPENMP)
    target_link_libraries(milca_feec_n PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_sources(milca_mhd PRIVATE src/mhd/mhd_mesh.f90 src/mhd/mhd_ampere_n.f90
    src/generated/milca_ideal_response_error_leaves.f90
    src/mhd/mhd_response_n.f90 src/mhd/mhd_ideal_n.f90 src/mhd/mhd_closure_n.f90
    src/mhd/mhd_layer_n.f90)
target_link_libraries(milca_fem2d PUBLIC milca_feec_n)
target_sources(milca_mhd PRIVATE src/mhd/mhd_perturbation_n.f90)
add_executable(test_mhd_perturbation test/test_mhd_perturbation.f90)
target_link_libraries(test_mhd_perturbation PRIVATE milca)
add_test(NAME mhd_perturbation_l7_xi_lagrangian_charts COMMAND test_mhd_perturbation)
set_tests_properties(mhd_perturbation_l7_xi_lagrangian_charts PROPERTIES TIMEOUT 300)
add_executable(test_mhd_ampere_error test/test_mhd_ampere_error.f90)
target_link_libraries(test_mhd_ampere_error PRIVATE milca)
add_test(NAME mhd_ampere_n_error_status COMMAND test_mhd_ampere_error)
set_tests_properties(mhd_ampere_n_error_status PROPERTIES TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ampere_sheet_error test/test_mhd_ampere_sheet_error.f90)
target_link_libraries(test_mhd_ampere_sheet_error PRIVATE milca)
add_test(NAME mhd_ampere_n_sheet_coupled_equilibrated
    COMMAND test_mhd_ampere_sheet_error)
set_tests_properties(mhd_ampere_n_sheet_coupled_equilibrated PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_meridian_edge_remainder
    ${fortfem_SOURCE_DIR}/test/test_meridian_edge_remainder.f90)
target_link_libraries(test_meridian_edge_remainder PRIVATE milca)
add_test(NAME mhd_ampere_n_sheet_interpolation_remainder
    COMMAND test_meridian_edge_remainder)
set_tests_properties(mhd_ampere_n_sheet_interpolation_remainder PROPERTIES TIMEOUT 30
    PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ampere_n test/test_mhd_ampere_n.f90)
target_link_libraries(test_mhd_ampere_n PRIVATE milca)
add_test(NAME mhd_ampere_n_exterior_manufactured COMMAND test_mhd_ampere_n)
set_tests_properties(mhd_ampere_n_exterior_manufactured PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_shielding test/test_mhd_shielding.f90)
target_link_libraries(test_mhd_shielding PRIVATE milca)
add_test(NAME mhd_shielding_cylinder_limit_reciprocity COMMAND test_mhd_shielding)
set_tests_properties(mhd_shielding_cylinder_limit_reciprocity PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(gs_response_3d app/gs_response_3d.f90)
target_link_libraries(gs_response_3d PRIVATE milca)
add_test(NAME gs_response_3d_ampere_indicator_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_response_3d> -DWORK=${CMAKE_BINARY_DIR}/gs_response_3d_check
    -DEQDSK=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_response_3d_app.cmake)
set_tests_properties(gs_response_3d_ampere_indicator_contract PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ideal_n test/test_mhd_ideal_n.f90)
target_link_libraries(test_mhd_ideal_n PRIVATE milca)
add_test(NAME mhd_ideal_n_mde_force_free_identity COMMAND test_mhd_ideal_n)
add_test(NAME mhd_ideal_n_final_state_error COMMAND test_mhd_ideal_n error)
set_tests_properties(mhd_ideal_n_final_state_error PROPERTIES TIMEOUT 60
    ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_test(NAME gs_response_3d_ideal_error_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_response_3d>
    -DWORK=${CMAKE_BINARY_DIR}/gs_response_3d_ideal_check
    -DEQDSK=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_response_3d_ideal_error.cmake)
# ~45 s serial (complex Ampere LU) on a 2026 workstation, ~60 s on sCluster.
set_tests_properties(gs_response_3d_ideal_error_contract PROPERTIES TIMEOUT 240
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
set_tests_properties(mhd_ideal_n_mde_force_free_identity PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_shield_torus test/test_mhd_shield_torus.f90)
target_link_libraries(test_mhd_shield_torus PRIVATE milca)
add_test(NAME mhd_shield_torus_nrad_convergence COMMAND test_mhd_shield_torus
    ${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk)
set_tests_properties(mhd_shield_torus_nrad_convergence PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_flux_mesh_axis test/test_flux_mesh_axis.f90)
target_link_libraries(test_flux_mesh_axis PRIVATE milca)
add_test(NAME flux_mesh_axis_interpolation_order COMMAND test_flux_mesh_axis)
set_tests_properties(flux_mesh_axis_interpolation_order PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(test_mhd_closure test/test_mhd_closure.f90)
target_link_libraries(test_mhd_closure PRIVATE milca)
add_test(NAME response_layer_all_rational_setup COMMAND test_mhd_closure setup)
set_tests_properties(response_layer_all_rational_setup PROPERTIES TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_test(NAME gs_response_3d_frobenius_all_rational COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_response_3d>
    -DWORK=${CMAKE_BINARY_DIR}/gs_response_3d_frobenius_check
    -DEQDSK=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk
    -P ${CMAKE_SOURCE_DIR}/test/check_response_all_rational.cmake)
set_tests_properties(gs_response_3d_frobenius_all_rational PROPERTIES TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_test(NAME mhd_closure_tearing_delta_prime COMMAND test_mhd_closure)
set_tests_properties(mhd_closure_tearing_delta_prime PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
# R-B2 increment 1: strip coupled to the FEEC torus (response-unified.md, R-B2).
# Gated at nrad 64 (smallest passing nrad of the (b) ladder), strip half-width
# xd 0.01 plus an (a)-only run at xd 0.02. Physical strip current is checked
# at both widths; interface-sheet sum and its width order remain diagnostics
# with their existing gates (GMRES tolerance 1e-7).
add_test(NAME response_strip_torus_coupling COMMAND test_mhd_closure strip)
set_tests_properties(response_strip_torus_coupling PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
# R-B1: resolved resistive layer strip, periodic cylinder (response-unified.md (B)).
target_sources(milca_mhd PRIVATE src/mhd/response_strip.f90)
target_sources(milca_mhd PRIVATE src/generated/milca_response_fluid_leaves.f90
    src/mhd/response_fluid.f90 src/mhd/response_toroidal_strip.f90)
add_executable(test_response_toroidal_strip test/test_response_toroidal_strip.f90)
target_link_libraries(test_response_toroidal_strip PRIVATE milca)
add_executable(test_response_interfaces test/test_response_synthetic_interfaces.f90)
target_link_libraries(test_response_interfaces PRIVATE milca)
add_test(NAME response_toroidal_strip_two_ring_traces
    COMMAND test_response_interfaces traces)
add_test(NAME response_cut_total_pressure_balance
    COMMAND test_response_interfaces cut)
add_test(NAME response_toroidal_strip_edge_merged_domains
    COMMAND test_response_interfaces domains)
set_tests_properties(response_toroidal_strip_two_ring_traces
    response_cut_total_pressure_balance response_toroidal_strip_edge_merged_domains
    PROPERTIES TIMEOUT 120 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=PASSIVE")
add_test(NAME response_toroidal_strip_pressure_action
    COMMAND test_response_toroidal_strip action)
add_test(NAME response_toroidal_strip_pressure_solve
    COMMAND test_response_toroidal_strip solve)
add_test(NAME response_toroidal_strip_radial_convergence
    COMMAND test_response_toroidal_strip convergence)
set_tests_properties(response_toroidal_strip_pressure_action
    response_toroidal_strip_pressure_solve response_toroidal_strip_radial_convergence
    PROPERTIES TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=PASSIVE")
add_executable(test_response_fluid test/test_response_fluid.f90)
target_link_libraries(test_response_fluid PRIVATE milca)
add_test(NAME response_fluid_full_pressure_oracles COMMAND test_response_fluid)
set_tests_properties(response_fluid_full_pressure_oracles PROPERTIES
    TIMEOUT 120 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OMP_WAIT_POLICY=PASSIVE")
add_executable(test_response_strip test/test_response_strip.f90)
target_link_libraries(test_response_strip PRIVATE milca)
add_test(NAME response_strip_cylinder_layer COMMAND test_response_strip)
add_library(milca_rm_oracles STATIC src/generated/milca_rm_oracle_leaves.f90
    test/support/rm_mhd_oracles.f90)
target_link_libraries(milca_rm_oracles PUBLIC fortnum)
add_library(milca_spectrum_support STATIC test/support/vrmhd_structural_spectrum.f90
    test/support/vrmhd_spectrum_support.f90)
target_link_libraries(milca_spectrum_support PUBLIC milca_mhd milca_rm_oracles)
add_executable(test_vrmhd_spectra test/test_vrmhd_spectra.f90)
target_link_libraries(test_vrmhd_spectra PRIVATE milca_spectrum_support)
add_executable(test_vrmhd_spectral_folding test/test_vrmhd_spectral_folding.f90)
target_link_libraries(test_vrmhd_spectral_folding PRIVATE milca_spectrum_support)
add_test(NAME vrmhd_spectral_folding COMMAND test_vrmhd_spectral_folding)
set_tests_properties(vrmhd_spectral_folding PROPERTIES LABELS rm2)
foreach(mode spurious ambiguous)
    add_test(NAME vrmhd_resolved_power_rejects_${mode}
        COMMAND test_vrmhd_spectral_folding ${mode})
    set_tests_properties(vrmhd_resolved_power_rejects_${mode}
        PROPERTIES LABELS rm2 WILL_FAIL TRUE)
endforeach()
add_executable(test_vrmhd_fourier test/test_vrmhd_fourier.f90)
target_link_libraries(test_vrmhd_fourier PRIVATE milca_spectrum_support)
add_test(NAME vrmhd_fourier_content COMMAND test_vrmhd_fourier)
set_tests_properties(vrmhd_fourier_content PROPERTIES LABELS rm2)
add_test(NAME vrmhd_section_spectra COMMAND test_vrmhd_spectra)
set_tests_properties(vrmhd_section_spectra PROPERTIES LABELS rm2 PROCESSORS 8
    TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_section_spectra_fine COMMAND test_vrmhd_spectra fine)
set_tests_properties(vrmhd_section_spectra_fine PROPERTIES PROCESSORS 8
    TIMEOUT 7200 ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
add_executable(vrmhd_pairing_study test/vrmhd_pairing_study.f90)
target_link_libraries(vrmhd_pairing_study PRIVATE milca_spectrum_support)
add_test(NAME vrmhd_pairing_measurement COMMAND vrmhd_pairing_study smoke)
set_tests_properties(vrmhd_pairing_measurement PROPERTIES PROCESSORS 8 TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
add_test(NAME vrmhd_pairing_benchmark COMMAND vrmhd_pairing_study
    ${CMAKE_CURRENT_BINARY_DIR}/vrmhd-pairings.csv)
set_tests_properties(vrmhd_pairing_benchmark PROPERTIES PROCESSORS 8 TIMEOUT 3600
    LABELS benchmark ENVIRONMENT "OMP_NUM_THREADS=8;OPENBLAS_NUM_THREADS=1")
add_executable(test_rm_oracles test/test_rm_oracles.f90)
target_link_libraries(test_rm_oracles PRIVATE milca_rm_oracles milca_mhd)
foreach(case waves rates outer newcomb convergence saturation table)
    add_test(NAME rm_oracle_${case} COMMAND test_rm_oracles ${case})
    set_tests_properties(rm_oracle_${case} PROPERTIES TIMEOUT 60
        LABELS rm-oracle ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
set_tests_properties(response_strip_cylinder_layer PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(response_strip_ladder app/response_strip_ladder.f90)
target_link_libraries(response_strip_ladder PRIVATE milca)
add_test(NAME mhd_resistive_layer_rotation COMMAND test_mhd_closure resistive)
set_tests_properties(mhd_resistive_layer_rotation PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_newcomb test/test_mhd_newcomb.f90)
target_link_libraries(test_mhd_newcomb PRIVATE milca)
add_test(NAME mhd_newcomb_screw_pinch_limit COMMAND test_mhd_newcomb)
set_tests_properties(mhd_newcomb_screw_pinch_limit PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")

# NN1: closed free-B states and the vacuum / axisymmetric-pressure restrictions.
target_sources(milca_mhd PRIVATE src/generated/milca_nonnested_leaves.f90
    src/mhd/equilibrium_nonnested.f90 src/mhd/nonnested_residual.f90)
set_source_files_properties(src/generated/milca_nonnested_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")

add_library(milca_nn_solovev_test STATIC test/support/nn_solovev_fixture.f90)
target_link_libraries(milca_nn_solovev_test PUBLIC milca)
add_executable(test_nonnested_solovev test/test_nonnested_solovev.f90)
target_link_libraries(test_nonnested_solovev PRIVATE milca_nn_solovev_test)
foreach(k 2 3)
    add_test(NAME nonnested_n0_solovev_p${k} COMMAND test_nonnested_solovev ${k})
    set_tests_properties(nonnested_n0_solovev_p${k} PROPERTIES TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
# Equilibrated energy bound of the N = 0 restriction against the Solov'ev flux.
target_sources(milca_mhd PRIVATE src/mhd/nonnested_error.f90)
add_executable(test_nonnested_error test/test_nonnested_error.f90)
target_link_libraries(test_nonnested_error PRIVATE milca_nn_solovev_test)
add_test(NAME nonnested_n0_equilibrated_solovev COMMAND test_nonnested_error)
set_tests_properties(nonnested_n0_equilibrated_solovev PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")
add_executable(test_nonnested_fourier_error
    ${fortfem_SOURCE_DIR}/test/test_meridian_fourier_equilibration_slow.f90)
target_link_libraries(test_nonnested_fourier_error PRIVATE milca_feec_n)
add_test(NAME nonnested_transport_fourier_equilibrated
    COMMAND test_nonnested_fourier_error)
set_tests_properties(nonnested_transport_fourier_equilibrated PROPERTIES
    TIMEOUT 60 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

add_executable(test_nonnested_transport_error test/test_nonnested_transport_error.f90)
target_link_libraries(test_nonnested_transport_error PRIVATE milca)
add_test(NAME nonnested_transport_equilibrated_n0
    COMMAND test_nonnested_transport_error)
set_tests_properties(nonnested_transport_equilibrated_n0 PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")

add_executable(test_nonnested_anisotropy test/test_nonnested_anisotropy.f90)
target_link_libraries(test_nonnested_anisotropy PRIVATE milca)
add_test(NAME nonnested_transport_anisotropy_pollution
    COMMAND test_nonnested_anisotropy pollution)
add_test(NAME nonnested_transport_anisotropy_correction
    COMMAND test_nonnested_anisotropy correction)
set_tests_properties(nonnested_transport_anisotropy_pollution
    nonnested_transport_anisotropy_correction PROPERTIES TIMEOUT 60
    PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

target_sources(milca_mhd PRIVATE src/mhd/island_diagnostics.f90)
target_sources(milca_equilibrium PRIVATE src/generated/milca_island_leaves.f90
    src/equilibrium/coil_quadrature_field.f90)
add_executable(test_island_diagnostics test/test_island_diagnostics.f90)
add_executable(test_island_driven test/support/island_outer_reference.f90
    test/test_island_driven.f90)
target_link_libraries(test_island_driven PRIVATE milca_mhd)
add_test(NAME island_driven_outer_reference COMMAND test_island_driven outer)
add_test(NAME island_driven_stable_outer_reference
    COMMAND test_island_driven stable_outer)
add_test(NAME island_zero_beta_ohmic COMMAND test_island_driven base 3 4 1)
set_tests_properties(island_zero_beta_ohmic PROPERTIES TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
# The NN2 driven and stable gates (test_island_driven driven|stable 3 8 1)
# are open studies, not registered: docs/island-driven-cylinder.md.
target_link_libraries(test_island_diagnostics PRIVATE milca_mhd)
add_test(NAME island_diagnostics_pendulum COMMAND test_island_diagnostics)
set_tests_properties(island_diagnostics_pendulum PROPERTIES TIMEOUT 60 PROCESSORS 2
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")

add_library(milca_nn_vacuum_test STATIC test/support/nn_vacuum_fixture.f90
    test/support/island_a_fixture.f90)
add_executable(test_island_vacuum_oracle test/test_island_vacuum_oracle.f90)
target_link_libraries(test_island_vacuum_oracle PRIVATE milca_nn_vacuum_test)
add_test(NAME island_coil_quadrature COMMAND test_island_vacuum_oracle quadrature)
add_executable(test_island_direct_coil_split test/test_island_direct_coil_split.f90)
target_link_libraries(test_island_direct_coil_split PRIVATE milca_nn_vacuum_test)
add_test(NAME island_direct_coil_split COMMAND test_island_direct_coil_split)
set_tests_properties(island_direct_coil_split PROPERTIES TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME island_vacuum_seed_p3 COMMAND test_island_vacuum_oracle seed 3 16)
set_tests_properties(island_vacuum_seed_p3 PROPERTIES TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME island_coarse_feec_map COMMAND test_island_vacuum_oracle coarse 3 8)
add_test(NAME island_coarse_feec_map_p2 COMMAND test_island_vacuum_oracle coarse 2 16)
foreach(degree 2 3)
    add_test(NAME island_vacuum_width_phase_p${degree}
        COMMAND test_island_vacuum_oracle ladder ${degree} 4 3)
    set_tests_properties(island_vacuum_width_phase_p${degree} PROPERTIES TIMEOUT 3600
        PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
set_tests_properties(island_coarse_feec_map island_coarse_feec_map_p2 PROPERTIES TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
target_link_libraries(milca_nn_vacuum_test PUBLIC milca)

add_executable(nn1_spec_case test/nn1_spec_case.f90)
target_link_libraries(nn1_spec_case PRIVATE milca)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
    add_test(NAME nn1_spec_case_smoke
        COMMAND ${CMAKE_COMMAND}
            -DEXE=$<TARGET_FILE:nn1_spec_case>
            -DDATA=${CMAKE_SOURCE_DIR}/test/data/nn1_spec
            -DCHECK=${CMAKE_SOURCE_DIR}/test/check_nn1_spec_case.py
            -DPY=${Python3_EXECUTABLE} -DOUT=${CMAKE_CURRENT_BINARY_DIR}/nn1_spec_smoke
            -P ${CMAKE_SOURCE_DIR}/test/run_nn1_spec_case.cmake)
    set_tests_properties(nn1_spec_case_smoke PROPERTIES TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=4")
endif()

add_executable(test_nonnested_vacuum test/test_nonnested_vacuum.f90)
target_link_libraries(test_nonnested_vacuum PRIVATE milca_nn_vacuum_test)
# NN1 per-mode solve ladder: Cholesky versus Hiptmair-Xu CG (benchmarks/nn1_hx).
add_executable(nn1_hx_ladder benchmarks/nn1_hx/nn1_hx_ladder.f90)
target_link_libraries(nn1_hx_ladder PRIVATE milca_nn_vacuum_test)
add_executable(test_nonnested_hx test/test_nonnested_hx.f90)
target_link_libraries(test_nonnested_hx PRIVATE milca_nn_vacuum_test)
add_test(NAME nonnested_hx_vacuum_p2 COMMAND test_nonnested_hx 2 4 8)
add_test(NAME nonnested_hx_vacuum_p3 COMMAND test_nonnested_hx 3 4 8)
set_tests_properties(nonnested_hx_vacuum_p2 nonnested_hx_vacuum_p3 PROPERTIES
    TIMEOUT 120 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME nonnested_vacuum_helical_coil_island COMMAND test_nonnested_vacuum)
add_executable(test_nonnested_vacuum_mode_residual
    test/test_nonnested_vacuum_mode_residual.f90)
target_link_libraries(test_nonnested_vacuum_mode_residual PRIVATE milca_nn_vacuum_test)
add_test(NAME nonnested_vacuum_mode_residual_p3
    COMMAND test_nonnested_vacuum_mode_residual)
add_test(NAME nonnested_vacuum_mode_residual_p3_nr32
    COMMAND test_nonnested_vacuum_mode_residual 32)
set_tests_properties(nonnested_vacuum_mode_residual_p3 PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4")
set_tests_properties(nonnested_vacuum_mode_residual_p3_nr32 PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_tests_properties(nonnested_vacuum_helical_coil_island PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")

add_executable(test_nonnested_leaves test/test_nonnested_leaves.f90)
target_link_libraries(test_nonnested_leaves PRIVATE milca)
add_test(NAME nonnested_symbolic_derivatives_and_d3_identity
    COMMAND test_nonnested_leaves)

# NN2 (first increment): anisotropic pressure transport for a given B.
target_sources(milca_mhd PRIVATE src/mhd/nonnested_transport.f90)
add_executable(test_nonnested_transport test/test_nonnested_transport.f90)
target_link_libraries(test_nonnested_transport PRIVATE milca)
add_test(NAME nonnested_transport_manufactured_3d
    COMMAND test_nonnested_transport manufactured)
set_tests_properties(nonnested_transport_manufactured_3d PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
foreach(k 2 3 4)
    add_test(NAME nonnested_transport_sovinec_chi_ratio_p${k}
        COMMAND test_nonnested_transport sovinec ${k})
    set_tests_properties(nonnested_transport_sovinec_chi_ratio_p${k} PROPERTIES
        TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()

# NN2 oracle: Hudson & Breslau (2008) anisotropic heat transport in a chaotic
# field (docs/equilibrium-nonnested-nn2.md).
add_executable(test_nonnested_hudson_breslau test/test_nonnested_hudson_breslau.f90)
target_link_libraries(test_nonnested_hudson_breslau PRIVATE milca)
add_test(NAME nonnested_hudson_breslau_integrable_exact_p3
    COMMAND test_nonnested_hudson_breslau integrable 3)
set_tests_properties(nonnested_hudson_breslau_integrable_exact_p3 PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")
add_test(NAME nonnested_hudson_breslau_chaos_responds_p3
    COMMAND test_nonnested_hudson_breslau chaos 3 7 48 3)
set_tests_properties(nonnested_hudson_breslau_chaos_responds_p3 PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# NN2 (second increment): axisymmetric coupled transport / Ampere / Ohm Newton.
target_sources(milca_mhd PRIVATE src/mhd/nonnested_profile.f90
    src/mhd/island_closure.f90
    src/mhd/nonnested_coupled.f90
    src/mhd/island_case.f90 src/mhd/island_report.f90
    src/mhd/island_provider.f90 src/mhd/frozen_current.f90)
add_executable(test_frozen_current test/test_frozen_current.f90)
target_link_libraries(test_frozen_current PRIVATE milca_mhd)
add_test(NAME frozen_current_ring_p3 COMMAND test_frozen_current ring 3 8)
add_test(NAME frozen_current_helical_p2 COMMAND test_frozen_current helical 2)
add_test(NAME frozen_current_helical_p3 COMMAND test_frozen_current helical 3)
add_test(NAME frozen_current_exterior_p3
    COMMAND test_frozen_current exterior 3 4 4 1.5e-3)
add_test(NAME frozen_current_exterior_refined_p3
    COMMAND test_frozen_current exterior 3 8 6 2e-4)
add_test(NAME frozen_current_source_p3 COMMAND test_frozen_current source 3 4 1e-6)
set_tests_properties(frozen_current_ring_p3 frozen_current_helical_p2
    frozen_current_helical_p3 frozen_current_exterior_p3
    frozen_current_exterior_refined_p3 frozen_current_source_p3 PROPERTIES
    TIMEOUT 1800 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_island_closure test/test_island_closure.f90)
target_link_libraries(test_island_closure PRIVATE milca_nn_solovev_test)
foreach(case declaration drive report)
    add_test(NAME island_closure_${case} COMMAND test_island_closure ${case})
    set_tests_properties(island_closure_${case} PROPERTIES TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
endforeach()
add_executable(test_nonnested_actions test/test_nonnested_actions.f90)
target_link_libraries(test_nonnested_actions PRIVATE milca_nn_solovev_test)
foreach(k 2 3)
    foreach(case fixed floating voltage)
        add_test(NAME nonnested_actions_${case}_p${k}
            COMMAND test_nonnested_actions ${case} ${k})
        set_tests_properties(nonnested_actions_${case}_p${k} PROPERTIES
            TIMEOUT 180 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()
add_test(NAME nonnested_actions_invalid COMMAND test_nonnested_actions invalid 2)
set_tests_properties(nonnested_actions_invalid PROPERTIES TIMEOUT 180 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_nonnested_profile test/test_nonnested_profile.f90)
target_link_libraries(test_nonnested_profile PRIVATE milca_nn_solovev_test)
foreach(k 2 3)
    foreach(case manufactured nested reproduce reproduce3d sampling)
        add_test(NAME nonnested_fitted_${case}_p${k}
            COMMAND test_nonnested_profile ${case} ${k})
        set_tests_properties(nonnested_fitted_${case}_p${k} PROPERTIES
            TIMEOUT 600 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()
add_executable(test_island_fitted test/test_island_fitted.f90)
target_link_libraries(test_island_fitted PRIVATE milca_nn_solovev_test)
foreach(case readiness heat current)
    add_test(NAME island_fitted_${case} COMMAND test_island_fitted ${case})
    set_tests_properties(island_fitted_${case} PROPERTIES TIMEOUT 180 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_library(milca_g6_comparison_test STATIC
    test/support/g6_sampled_force.f90 test/support/g6_external_comparison.f90)
target_compile_options(milca_g6_comparison_test PRIVATE
    -fno-fast-math -ffp-contract=off)
add_executable(test_g6_sampled_force test/test_g6_sampled_force.f90)
target_link_libraries(test_g6_sampled_force PRIVATE milca_g6_comparison_test)
add_test(NAME nonnested_g6_external_sampled_force COMMAND test_g6_sampled_force)
set_tests_properties(nonnested_g6_external_sampled_force PROPERTIES
    TIMEOUT 30 PROCESSORS 1 ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(test_g6_external_comparison test/test_g6_external_comparison.f90)
target_link_libraries(test_g6_external_comparison PRIVATE milca_g6_comparison_test)
add_test(NAME nonnested_g6_external_comparison_contract
    COMMAND test_g6_external_comparison)
set_tests_properties(nonnested_g6_external_comparison_contract PROPERTIES
    TIMEOUT 30 PROCESSORS 1 ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(test_nonnested_coupled test/test_nonnested_coupled.f90
    test/support/g6_solovev_oracle.f90)
target_link_libraries(test_nonnested_coupled PRIVATE milca_nn_solovev_test
    milca_g6_comparison_test)
foreach(k 2 3)
    add_test(NAME nonnested_coupled_solovev_p${k}
        COMMAND test_nonnested_coupled solovev ${k} 1e6)
    add_test(NAME nonnested_coupled_offset_return_p${k}
        COMMAND test_nonnested_coupled offset ${k} 1e6)
    add_test(NAME nonnested_coupled_fourier_reduction_p${k}
        COMMAND test_nonnested_coupled reduction ${k} 1e6)
    set_tests_properties(nonnested_coupled_solovev_p${k}
        nonnested_coupled_offset_return_p${k}
        nonnested_coupled_fourier_reduction_p${k} PROPERTIES
        TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
add_executable(test_nonnested_island test/test_nonnested_island.f90)
target_link_libraries(test_nonnested_island PRIVATE milca_nn_vacuum_test)
add_test(NAME island_coarse_nn2_map COMMAND test_nonnested_island map 3 8 1 1e6)
set_tests_properties(island_coarse_nn2_map PROPERTIES TIMEOUT 600 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME nonnested_island_flattening_p3
    COMMAND test_nonnested_island ratio 3 8 1)
add_test(NAME nonnested_island_d3_ladder_p2
    COMMAND test_nonnested_island ladder 2 4 1 1e4 3)
set_tests_properties(nonnested_island_flattening_p3 nonnested_island_d3_ladder_p2
    PROPERTIES TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4")
# Closure A+D provider: implicit JVP/VJP and nested limit (island case).
add_executable(test_island_provider test/test_island_provider.f90)
target_link_libraries(test_island_provider PRIVATE milca_nn_vacuum_test)
# Study driver (flattening scan, Poincare data); run by scripts, not CTest.
add_executable(study_island_closure_a test/study_island_closure_a.f90)
target_link_libraries(study_island_closure_a PRIVATE milca_nn_vacuum_test)
foreach(case tangent nested vacuum)
    add_test(NAME island_provider_${case}_p2 COMMAND test_island_provider ${case} 2 4)
    set_tests_properties(island_provider_${case}_p2 PROPERTIES TIMEOUT 900
        PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
# Landreman family 1 (exact 3D, uniform iota) through the coupled solve.
add_executable(test_nonnested_landreman test/test_nonnested_landreman.f90)
target_link_libraries(test_nonnested_landreman PRIVATE milca)
add_test(NAME nonnested_coupled_landreman_family1_p2
    COMMAND test_nonnested_landreman 2 1 3 2)
set_tests_properties(nonnested_coupled_landreman_family1_p2
    PROPERTIES TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4")

# G0: retain the published P3 high-ratio closed-line oracle and its gates.
add_test(NAME nonnested_coupled_closed_line_p3
    COMMAND test_nonnested_landreman 3 1e6 3 4)
set_tests_properties(nonnested_coupled_closed_line_p3 PROPERTIES
    TIMEOUT 3600 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1" PROCESSORS 4)
# NN1 on a phi-dependent (mapped) domain: docs/nn1-mapped-3d.md.
target_sources(milca_mhd PRIVATE src/mhd/nonnested_mapped.f90
    src/mhd/nonnested_newton_mr.f90
    src/mhd/nonnested_beltrami.f90 src/mhd/nonnested_beltrami_mr.f90
    src/mhd/nonnested_relax.f90 src/mhd/nonnested_force.f90
    src/generated/milca_metriplectic_leaves.f90)
set_source_files_properties(src/generated/milca_metriplectic_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_nonnested_beltrami_spec test/test_nonnested_beltrami_spec.f90)
target_link_libraries(test_nonnested_beltrami_spec PRIVATE milca)
set(MILCA_BELTRAMI_SPEC_SAMPLES "" CACHE PATH "Optional native SPEC sample directory")
if(MILCA_BELTRAMI_SPEC_SAMPLES)
    foreach(case axis helical)
        if(case STREQUAL "axis")
            set(helical_amplitude 0)
            set(modes 0)
        else()
            set(helical_amplitude 0.1)
            set(modes 2)
        endif()
        foreach(rings 2 4 8)
            add_test(NAME nonnested_beltrami_spec_${case}_${rings}
                COMMAND test_nonnested_beltrami_spec
                ${MILCA_BELTRAMI_SPEC_SAMPLES}/${case}.samples
                ${helical_amplitude} ${rings} 2 ${modes}
                ${MILCA_BELTRAMI_SPEC_SAMPLES}/candidate_${case}_${rings}_${modes}.dat)
            set_tests_properties(nonnested_beltrami_spec_${case}_${rings}
                PROPERTIES TIMEOUT 900 PROCESSORS 4
                ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
        endforeach()
    endforeach()
    foreach(modes 4 6)
        add_test(NAME nonnested_beltrami_spec_modes_${modes}
            COMMAND test_nonnested_beltrami_spec
            ${MILCA_BELTRAMI_SPEC_SAMPLES}/helical.samples 0.1 4 2 ${modes}
            ${MILCA_BELTRAMI_SPEC_SAMPLES}/candidate_helical_4_${modes}.dat)
        set_tests_properties(nonnested_beltrami_spec_modes_${modes} PROPERTIES
            TIMEOUT 900 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endif()
add_executable(test_nonnested_beltrami_pairing
    ${fortfem_SOURCE_DIR}/test/check.f90
    ${fortfem_SOURCE_DIR}/test/meridian_test_meshes.f90
    ${fortfem_SOURCE_DIR}/test/test_meridian_fourier_curl.f90)
target_link_libraries(test_nonnested_beltrami_pairing PRIVATE milca_feec_n)
add_test(NAME nonnested_beltrami_pairing COMMAND test_nonnested_beltrami_pairing)
set_tests_properties(nonnested_beltrami_pairing PROPERTIES
    TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_nonnested_beltrami test/test_nonnested_beltrami.f90)
target_link_libraries(test_nonnested_beltrami PRIVATE milca)
add_test(NAME nonnested_beltrami_resonance
    COMMAND test_nonnested_beltrami resonance 2)
set_tests_properties(nonnested_beltrami_resonance PROPERTIES
    TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME nonnested_beltrami_constraints
    COMMAND test_nonnested_beltrami constraints 2)
set_tests_properties(nonnested_beltrami_constraints PROPERTIES
    TIMEOUT 60 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
foreach(k 2 3)
    foreach(case manufactured lundquist vacuum)
        add_test(NAME nonnested_beltrami_${case}_p${k}
            COMMAND test_nonnested_beltrami ${case} ${k})
        set_tests_properties(nonnested_beltrami_${case}_p${k} PROPERTIES
            TIMEOUT 600 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()

# Matrix-free MINRES / Newton-MR Taylor closure (opt-in), against the
# factored path: docs/nonnested-newton-mr.md.
foreach(case constraints vacuum)
    add_test(NAME nonnested_beltrami_${case}_mr
        COMMAND test_nonnested_beltrami ${case} 2 mr)
    set_tests_properties(nonnested_beltrami_${case}_mr PROPERTIES
        TIMEOUT 60 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
foreach(k 2 3)
    foreach(case manufactured lundquist)
        add_test(NAME nonnested_beltrami_${case}_p${k}_mr
            COMMAND test_nonnested_beltrami ${case} ${k} mr)
        set_tests_properties(nonnested_beltrami_${case}_p${k}_mr PROPERTIES
            TIMEOUT 600 PROCESSORS 4
            ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    endforeach()
endforeach()
add_executable(test_nonnested_beltrami_mr test/test_nonnested_beltrami_mr.f90)
target_link_libraries(test_nonnested_beltrami_mr PRIVATE milca)
add_test(NAME nonnested_beltrami_mr_agreement
    COMMAND test_nonnested_beltrami_mr agree 2 2 2)
add_test(NAME nonnested_beltrami_mr_block_ladder
    COMMAND test_nonnested_beltrami_mr ladder block 2 2 4.5 0.1 0 2 4 8)
set_tests_properties(nonnested_beltrami_mr_agreement nonnested_beltrami_mr_block_ladder
    PROPERTIES TIMEOUT 120 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

add_executable(test_nonnested_relax test/test_nonnested_relax.f90)
target_link_libraries(test_nonnested_relax PRIVATE milca_nn_vacuum_test)
foreach(case taylor ideal)
    add_test(NAME nonnested_relax_${case} COMMAND test_nonnested_relax ${case})
    set_tests_properties(nonnested_relax_${case} PROPERTIES
        TIMEOUT 300 LABELS "mhd;fast" PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
# Harmonic deflation, potential-velocity and Leray force, residual floor:
# docs/nonnested-force.md.
add_executable(test_nonnested_force test/test_nonnested_force.f90)
target_link_libraries(test_nonnested_force PRIVATE milca_nn_solovev_test
    milca_nn_vacuum_test)
foreach(case hodge solovev beltrami shell sheet)
    add_test(NAME nonnested_force_${case} COMMAND test_nonnested_force ${case})
    set_tests_properties(nonnested_force_${case} PROPERTIES
        TIMEOUT 600 PROCESSORS 2
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
endforeach()
add_test(NAME nonnested_relax_vacuum_island COMMAND test_nonnested_relax island)
set_tests_properties(nonnested_relax_vacuum_island PROPERTIES
    TIMEOUT 1800 LABELS "mhd;slow" PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

# Energy-criterion island seeding (docs/island-seed.md).
target_sources(milca_mhd PRIVATE src/generated/milca_island_seed_leaves.f90
    src/mhd/island_seed.f90)
set_source_files_properties(src/generated/milca_island_seed_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_island_seed test/test_island_seed.f90)
target_link_libraries(test_island_seed PRIVATE milca)
foreach(case seed none ladder leaves)
    add_test(NAME island_seed_${case} COMMAND test_island_seed ${case})
endforeach()
set_tests_properties(island_seed_seed island_seed_ladder PROPERTIES
    TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_tests_properties(island_seed_none island_seed_leaves PROPERTIES
    TIMEOUT 300 PROCESSORS 1
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

add_executable(test_nonnested_mapped test/test_nonnested_mapped.f90)
target_link_libraries(test_nonnested_mapped PRIVATE milca)
add_test(NAME nonnested_mapped_identity_regression
    COMMAND test_nonnested_mapped identity)
add_test(NAME nonnested_mapped_ellipse_h_ladder
    COMMAND test_nonnested_mapped hladder 10 4 8 16)
add_test(NAME nonnested_mapped_ellipse_mode_ladder
    COMMAND test_nonnested_mapped modes 2 1 3)
add_test(NAME nonnested_mapped_zero_forcing_residual
    COMMAND test_nonnested_mapped zero_forcing)
set_tests_properties(nonnested_mapped_identity_regression
    nonnested_mapped_ellipse_h_ladder nonnested_mapped_ellipse_mode_ladder
    nonnested_mapped_zero_forcing_residual
    PROPERTIES TIMEOUT 900 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
foreach(mapped_case point_oracles cartesian_identity cartesian_quadrature cartesian_vacuum)
    add_test(NAME nonnested_mapped_${mapped_case}
        COMMAND test_nonnested_mapped ${mapped_case})
    set_tests_properties(nonnested_mapped_${mapped_case} PROPERTIES
        PROCESSORS 4 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()

# Matrix-free Newton-MR for mapped NN1: docs/nonnested-newton-mr.md.
add_executable(test_nonnested_newton_mr test/test_nonnested_newton_mr.f90)
target_link_libraries(test_nonnested_newton_mr PRIVATE milca)
add_test(NAME nonnested_newton_mr_agreement
    COMMAND test_nonnested_newton_mr agree 4 4)
add_test(NAME nonnested_newton_mr_lumped_ladder
    COMMAND test_nonnested_newton_mr ladder lumped 4 2 4 8)
set_tests_properties(nonnested_newton_mr_agreement nonnested_newton_mr_lumped_ladder
    PROPERTIES TIMEOUT 900 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

# Precise-QA vacuum benchmark of mapped NN1: docs/nn1-qa-benchmark.md.
add_executable(test_nonnested_qa test/test_nonnested_qa.f90)
target_link_libraries(test_nonnested_qa PRIVATE milca)
add_test(NAME nonnested_qa_oracle_check
    COMMAND test_nonnested_qa ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/nn1_qa check)
set_tests_properties(nonnested_qa_oracle_check
    PROPERTIES TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")

# Electromagnetic wave at a plasma edge: 1D Yee/Mur fields coupled to
# relativistic Boris macroparticles (docs/em-slab-reflection.md; appended).
if(MILCA_FORTSYM_BUILD)
    add_executable(gen_em_slab gen/gen_em_slab.f90)
    target_link_libraries(gen_em_slab PRIVATE milca_gen_support)
endif()
target_sources(milca_continuum PRIVATE src/generated/milca_em_slab_leaves.f90
    src/continuum/em_slab.f90)
set_source_files_properties(src/generated/milca_em_slab_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_em_slab test/test_em_slab.f90)
target_link_libraries(test_em_slab PRIVATE milca)
foreach(em_case oracles vacuum reflection transmission refinement)
    add_test(NAME em_slab_${em_case} COMMAND test_em_slab ${em_case})
    set_tests_properties(em_slab_${em_case} PROPERTIES
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
add_executable(em_slab_demo app/em_slab_demo.f90)
target_link_libraries(em_slab_demo PRIVATE milca milca_io)

# R-B2 jets: independent of outer coupling and cut-boundary orchestration.
target_sources(milca_mhd PRIVATE src/generated/milca_response_strip_jets_leaves.f90
    src/mhd/response_strip_jets.f90)
set_source_files_properties(src/generated/milca_response_strip_jets_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_response_strip_jets test/test_response_strip_jets.f90)
target_link_libraries(test_response_strip_jets PRIVATE milca)
add_test(NAME response_strip_equilibrium_jets COMMAND test_response_strip_jets)
set_tests_properties(response_strip_equilibrium_jets PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive"
    PROCESSORS 2 TIMEOUT 60)
set(MILCA_STRIP_EQDSK_TEST "" CACHE FILEPATH
    "Optional imported equilibrium for independent strip-jet difference checks")
if(MILCA_STRIP_EQDSK_TEST)
    add_test(NAME response_strip_tc24_jets
        COMMAND test_response_strip_jets ${MILCA_STRIP_EQDSK_TEST})
    set_tests_properties(response_strip_tc24_jets PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1" TIMEOUT 60)
endif()
target_sources(milca_mhd PRIVATE src/mhd/response_strip_preconditioner.f90)
add_executable(test_response_strip_preconditioner
    test/test_response_strip_preconditioner.f90)
target_link_libraries(test_response_strip_preconditioner PRIVATE milca)
if(MILCA_STRIP_EQDSK_TEST)
    add_executable(test_response_strip_tc24_preconditioner
        test/test_response_strip_tc24_preconditioner.f90)
    target_link_libraries(test_response_strip_tc24_preconditioner PRIVATE milca)
    add_test(NAME response_strip_tc24_preconditioner
        COMMAND test_response_strip_tc24_preconditioner ${MILCA_STRIP_EQDSK_TEST})
    set_tests_properties(response_strip_tc24_preconditioner PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive"
        PROCESSORS 2 TIMEOUT 180)
endif()
foreach(check oracles negative limits scaling)
    add_test(NAME response_strip_preconditioner_${check}
        COMMAND test_response_strip_preconditioner ${check})
    set_tests_properties(response_strip_preconditioner_${check} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive"
    PROCESSORS 2 TIMEOUT 180)
endforeach()

include(cmake/response_calibration.cmake)

# Narrow umbrella links to the used layers, label tests, add the per-class
# milca_check_<class> targets and the layering test (cmake/milca_layers.cmake).
add_executable(test_gs_symmetry test/test_gs_symmetry.f90)
target_link_libraries(test_gs_symmetry PRIVATE milca_equilibrium)
foreach(gate rotation translation screw manufactured limit derivatives guards)
    add_test(NAME grad_shafranov_rm0_${gate} COMMAND test_gs_symmetry ${gate})
endforeach()

include(cmake/milca_cpu_precision.cmake)

add_executable(test_map_equilibrium test/test_map_equilibrium.f90)
target_link_libraries(test_map_equilibrium PRIVATE milca_mhd)
foreach(k 2 3)
    add_test(NAME equilibrium_core_map_state_p${k} COMMAND test_map_equilibrium ${k})
    set_tests_properties(equilibrium_core_map_state_p${k} PROPERTIES
        PROCESSORS 8 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()

option(MILCA_BUILD_CUDA_PRECISION_BENCHMARK
    "Build experimental CUDA equilibrium precision benchmark" OFF)
if(MILCA_BUILD_CUDA_PRECISION_BENCHMARK)
    add_subdirectory(benchmarks/equilibrium-gpu)
endif()
include(cmake/milca_s0.cmake)
include(cmake/milca_mhd_batches.cmake)
include(cmake/milca_e21.cmake)
include(cmake/milca_separatrix.cmake)
include(cmake/milca_current_sheet.cmake)
# Existing PR47 R2 consumer, recovered on current providers; no selector adoption.
target_sources(milca_feec_n PRIVATE
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_parametric_bem_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_bem_high_order_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_exterior_harmonic_products.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_harmonic_field_products.f90
    ${fortfem_SOURCE_DIR}/src/bem/parametric_laplace_bem_3d.f90
    ${fortfem_SOURCE_DIR}/src/bem/parametric_bem_high_order.f90
    ${fortfem_SOURCE_DIR}/src/bem/bem_circulant.f90
    ${fortfem_SOURCE_DIR}/src/bem/exterior_harmonic_fields.f90)
target_sources(milca_equilibrium PRIVATE src/equilibrium/equilibrium_boundary_3d.f90)
target_sources(milca_mhd PRIVATE src/mhd/ideal_free_boundary_bem3d.f90)
add_executable(test_ideal_fb_surface test/test_ideal_fb_surface.f90)
target_link_libraries(test_ideal_fb_surface PRIVATE milca)
foreach(k 2 3)
    foreach(geometry cylinder circle helix)
        add_test(NAME ideal_fb_chart_surface_${geometry}_p${k}
            COMMAND test_ideal_fb_surface ${geometry} ${k})
    endforeach()
endforeach()
add_test(NAME ideal_fb_chart_surface_invalid COMMAND test_ideal_fb_surface invalid)
add_test(NAME ideal_fb_chart_surface_circulant COMMAND test_ideal_fb_surface circulant)
set_tests_properties(ideal_fb_chart_surface_circulant PROPERTIES
    TIMEOUT 120 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_test(NAME ideal_fb_exterior_frozen_derivatives
    COMMAND test_ideal_free_boundary frozen-exterior)
set_tests_properties(ideal_fb_exterior_frozen_derivatives PROPERTIES
    TIMEOUT 300 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
# Circuit currents as parameters of the nested-map free-boundary residual.
target_sources(milca_mhd PRIVATE src/mhd/ideal_fb_circuits.f90)
add_executable(test_ideal_fb_circuits test/test_ideal_fb_circuits.f90)
target_link_libraries(test_ideal_fb_circuits PRIVATE milca)
add_test(NAME ideal_fb_circuits_derivatives COMMAND test_ideal_fb_circuits derivatives)
add_test(NAME ideal_fb_circuits_exterior3d COMMAND test_ideal_fb_circuits exterior3d)
add_executable(test_ideal_fb_traces3d test/test_ideal_fb_traces3d.f90)
target_link_libraries(test_ideal_fb_traces3d PRIVATE milca)
add_test(NAME ideal_fb_traces3d_axisymmetric COMMAND test_ideal_fb_traces3d axisymmetric)
add_test(NAME ideal_fb_traces3d_helical COMMAND test_ideal_fb_traces3d helical)
add_test(NAME ideal_fb_traces3d_makegrid COMMAND test_ideal_fb_traces3d makegrid)
add_test(NAME ideal_fb_traces3d_curve_hessian COMMAND test_ideal_fb_traces3d curve-hessian)
add_test(NAME ideal_fb_exterior3d_periodic COMMAND test_ideal_fb_traces3d exterior)
set_tests_properties(ideal_fb_exterior3d_periodic PROPERTIES
    TIMEOUT 1200 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
# Synthetic nfp = 3 free-boundary case against VMEC++ (benchmarks/fb3d-synthetic).
add_executable(fb3d_synthetic app/fb3d_synthetic.f90)
target_link_libraries(fb3d_synthetic PRIVATE milca)
add_test(NAME ideal_fb_landreman_axisymmetric_circuits
    COMMAND test_landreman_free_boundary circuits)
add_test(NAME ideal_fb_landreman_helical_circuits
    COMMAND test_landreman_free_boundary helical)
set_tests_properties(ideal_fb_circuits_derivatives ideal_fb_circuits_exterior3d
    PROPERTIES TIMEOUT 300 PROCESSORS 4
    ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
set_tests_properties(ideal_fb_landreman_axisymmetric_circuits
    ideal_fb_landreman_helical_circuits PROPERTIES
    TIMEOUT 600 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")

if(MILCA_BUILD_GLL_SPIKE)
    if(NOT MILCA_ENABLE_OPENMP)
        message(FATAL_ERROR "The GLL spike requires MILCA_ENABLE_OPENMP=ON")
    endif()
    add_subdirectory(benchmarks/gll)
    add_test(NAME ideal_gll_coupled_block_oracle COMMAND test_gll_block_preconditioner)
    set_tests_properties(ideal_gll_coupled_block_oracle PROPERTIES
        TIMEOUT 60 PROCESSORS 1
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
    foreach(p 4 6 8)
        add_test(NAME ideal_gll_derivatives_p${p}
            COMMAND benchmark_ideal_gll ${p} 1 3)
        set_tests_properties(ideal_gll_derivatives_p${p} PROPERTIES
            TIMEOUT 120 PROCESSORS 2
            ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive")
    endforeach()
endif()
include(cmake/milca_libneo.cmake)
set(MILCA_TIAGO_SOURCE_DIR "" CACHE PATH "TIAGO source for the optional island inference benchmark")
if(MILCA_TIAGO_SOURCE_DIR)
    include(benchmarks/tiago-provider/CMakeLists.txt)
endif()
milca_layers_finalize()

# Test resource declaration. Every registered test gets an explicit OpenMP
# thread count, a matching PROCESSORS reservation, and OMP_WAIT_POLICY=passive.
# Without this, concurrently running tests each start one spinning OpenMP thread
# per core, oversubscribe the machine under `ctest -j`, and starve each other
# until a wall-clock limit trips. A test keeps the thread count it declared via
# ENVIRONMENT (or PROCESSORS); the rest get MILCA_TEST_DEFAULT_THREADS.
set(MILCA_TEST_DEFAULT_THREADS 2 CACHE STRING
    "OpenMP threads (and ctest PROCESSORS) for tests that do not declare their own")
get_directory_property(_milca_all_tests TESTS)
foreach(_t IN LISTS _milca_all_tests)
    get_test_property(${_t} ENVIRONMENT _env)
    get_test_property(${_t} PROCESSORS _procs)
    if(NOT _env)
        set(_env "")
    endif()
    set(_threads "")
    foreach(_item IN LISTS _env)
        if(_item MATCHES "^OMP_NUM_THREADS=([0-9]+)$")
            set(_threads "${CMAKE_MATCH_1}")
        endif()
    endforeach()
    if(NOT _threads)
        if(_procs)
            set(_threads "${_procs}")
        else()
            set(_threads "${MILCA_TEST_DEFAULT_THREADS}")
        endif()
        list(APPEND _env "OMP_NUM_THREADS=${_threads}")
    endif()
    if(NOT "${_env}" MATCHES "(^|;)OPENBLAS_NUM_THREADS=")
        list(APPEND _env "OPENBLAS_NUM_THREADS=1")
    endif()
    if(NOT "${_env}" MATCHES "(^|;)OMP_WAIT_POLICY=")
        list(APPEND _env "OMP_WAIT_POLICY=passive")
    endif()
    set_tests_properties(${_t} PROPERTIES ENVIRONMENT "${_env}" PROCESSORS "${_threads}")
endforeach()

include(cmake/milca_estimator_benchmark.cmake)

include(cmake/milca_stack_targets.cmake)
milca_stack_targets()

# Must follow the last add_executable: wraps test commands of all targets.
include(cmake/milca_test_output_cap.cmake)
