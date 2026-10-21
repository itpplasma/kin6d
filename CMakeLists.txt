cmake_minimum_required(VERSION 3.22)
project(kin6d LANGUAGES Fortran C)
include(CTest)
include(FetchContent)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
# OpenMP is on by default; KIN6D_ENABLE_OPENMP=OFF gives a serial build. The
# thread count and determinism mode are written to the application provenance.
option(KIN6D_ENABLE_OPENMP "Build with OpenMP threading" ON)
if(KIN6D_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
# Dependency pins: one source used by FetchContent and by result provenance.
set(KIN6D_FORTPLOT_REVISION "759957d2b196fc7382ea1c370958f7adaa2fc753")
set(KIN6D_FORTNUM_REVISION "0a27a7f84249fc91e42118a1a8ae7f8573de1bd4")
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
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_ode_residual_certificate.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_verified_matvec.f90)
target_compile_options(kin6d_intervals PRIVATE -fno-fast-math -ffp-contract=off)
add_library(kin6d_ode
    ${fortnum_SOURCE_DIR}/src/fortnum_kinds.f90
    ${fortnum_SOURCE_DIR}/src/roots/fortnum_roots.f90
    ${fortnum_SOURCE_DIR}/src/roots/fortnum_multiroot.f90
    ${fortnum_SOURCE_DIR}/src/quadrature/fortnum_quadrature.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_rk54_ck_tableau.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_cash_karp.f90
    ${fortnum_SOURCE_DIR}/src/ode/fortnum_ode_geometric.f90
    ${fortnum_SOURCE_DIR}/src/rng/fortnum_rng.f90
    ${fortnum_SOURCE_DIR}/src/verified/fortnum_sampled_norm.f90
    ${fortnum_SOURCE_DIR}/src/interp/fortnum_polynomial.f90
    ${fortnum_SOURCE_DIR}/src/interp/fortnum_bspline.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_legendre_recurrence_kernel.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_legendre.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_coulomb_series_kernel.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_coulomb.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_gamma.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_elliptic.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_hypergeom_2f1_term_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_toroidal_order_kernel.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_toroidal.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_faddeeva.f90
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_riccati_bessel.f90
    ${fortnum_SOURCE_DIR}/src/fortnum_ad_backend.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_det2_jvp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_det2_vjp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_det3_jvp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_det3_vjp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_inv2_jvp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_inv2_vjp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_inv3_jvp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortnum_inv3_vjp_kernel.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortad/fortnum_det2_jvp_fortad.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortad/fortnum_det2_vjp_fortad.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortad/fortnum_det3_jvp_fortad.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortad/fortnum_det3_vjp_fortad.f90
    ${fortnum_SOURCE_DIR}/src/generated/fortad/fortnum_inv2_jvp_fortad.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_linalg.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_parallel_vector.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_krylov.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_symmetric_eigen.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_bvls.f90
    ${fortnum_SOURCE_DIR}/src/linalg/fortnum_gsvd.f90)
target_compile_options(kin6d_ode PRIVATE -fno-fast-math -ffp-contract=off)
target_link_libraries(kin6d_ode PUBLIC kin6d_intervals)
if(KIN6D_ENABLE_OPENMP)
    # These numerical leaves are also called from threaded KIN6D operators.
    # Compile them for concurrent calls, including runtime-checked Debug.
    target_link_libraries(kin6d_ode PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(kin6d_intervals PUBLIC OpenMP::OpenMP_Fortran)
endif()
# Generated derivative families retain the primal ABI even for constant partials.
# Unused primal arguments are intentional; keep other diagnostics enabled.
set_source_files_properties(src/generated/kin6d_leaves.f90
    src/generated/kin6d_pp_fusion_leaves.f90
    src/generated/kin6d_pp_fusion_eft_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
# Layered libraries (ROADMAP "Build structure and fast development loop"):
# L0 kin6d_core; L1 kin6d_fem2d, kin6d_pairs; L2 model classes; `kin6d` is the
# umbrella. Each source names its layer; the kin6d_layering test checks the use
# graph against this partition.
include(cmake/kin6d_layers.cmake)
target_sources(kin6d_core PRIVATE src/generated/kin6d_leaves.f90 src/core/laboratory.f90
    src/core/sort_keys.f90 src/core/real_identity.f90
    src/generated/kin6d_species_leaves.f90 src/core/species.f90)
target_sources(kin6d_pairs PRIVATE src/generated/kin6d_coulomb_leaves.f90 src/pairs/ewald.f90
    src/generated/kin6d_ocp_leaves.f90 src/pairs/near_encounter.f90 src/pairs/pme.f90
    src/pairs/nbody.f90 src/generated/kin6d_interaction_leaves.f90
    src/pairs/interaction_laboratory.f90)
target_sources(kin6d_classical PRIVATE src/classical/ocp.f90 src/classical/pair_correlation.f90
    src/generated/kin6d_cloud_leaves.f90 src/classical/microscopic_cloud.f90)
target_sources(kin6d_continuum PRIVATE src/continuum/vlasov_response.f90
    src/generated/kin6d_poisson_leaves.f90 src/continuum/poisson_laboratory.f90
    src/generated/kin6d_fp_leaves.f90 src/continuum/fp_laboratory.f90
    src/generated/kin6d_transport_leaves.f90 src/continuum/transport_laboratory.f90
    src/continuum/transport_certificate.f90 src/continuum/transport_sampling.f90
    src/generated/kin6d_streaming_leaves.f90 src/continuum/streaming_comparison.f90
    src/continuum/streaming_particle_density.f90
    src/generated/kin6d_streaming_moments_leaves.f90 src/continuum/streaming_moments.f90
    src/generated/kin6d_debye_sphere_leaves.f90 src/continuum/debye_sphere_reference.f90
    src/generated/kin6d_slab_oscillation_leaves.f90 src/continuum/slab_oscillation.f90
    src/generated/kin6d_shielding_dynamics_leaves.f90 src/continuum/shielding_dynamics.f90)
set_source_files_properties(
    src/generated/kin6d_species_leaves.f90 src/generated/kin6d_coulomb_leaves.f90
    src/generated/kin6d_ocp_leaves.f90
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
    set(family_layer kin6d_quantum)
    if(quantum_family STREQUAL "classical_pair")
        set(family_layer kin6d_classical)
    endif()
    string(REPLACE "kin6d_" "" family_directory "${family_layer}")
    target_sources(${family_layer} PRIVATE
        src/generated/kin6d_${quantum_family}_leaves.f90
        src/${family_directory}/${quantum_family}.f90)
    set_source_files_properties(src/generated/kin6d_${quantum_family}_leaves.f90
        PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
endforeach()
# Q1/Q1b partial-wave encounter: banded symmetric eigensolver from LAPACK,
# the external numerical dependency FortNum already declares.
find_package(LAPACK REQUIRED)
target_sources(kin6d_pairs PRIVATE src/generated/kin6d_coulomb_encounter_leaves.f90
    src/pairs/radial_dvr.f90)
target_sources(kin6d_quantum PRIVATE src/generated/kin6d_pp_fusion_leaves.f90
    src/quantum/pp_fusion.f90
    src/generated/kin6d_pp_fusion_eft_leaves.f90 src/quantum/pp_fusion_eft.f90
    src/quantum/coulomb_encounter.f90 src/quantum/coulomb_kepler.f90)
set_source_files_properties(src/generated/kin6d_coulomb_encounter_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d_core PUBLIC LAPACK::LAPACK)
set(KIN6D_FORTFEM_REVISION "72a2c44bac62d61f08541359c79590b1dfcf99ab")
FetchContent_Declare(fortfem
    GIT_REPOSITORY https://github.com/lazy-fortran/fortfem.git
    GIT_TAG ${KIN6D_FORTFEM_REVISION})
FetchContent_GetProperties(fortfem)
if(NOT fortfem_POPULATED)
    FetchContent_Populate(fortfem)
endif()
# Narrow, unchanged upstream geometry/meshing slice; no private FEM geometry.
file(GLOB KIN6D_TRIANGLE_SOURCES CONFIGURE_DEPENDS
    ${fortfem_SOURCE_DIR}/src/triangle_compat/*.f90)
add_library(kin6d_spatial ${KIN6D_TRIANGLE_SOURCES}
    ${fortfem_SOURCE_DIR}/src/utils/kinds.f90
    ${fortfem_SOURCE_DIR}/src/assembly/element_to_csr.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/delaunay_types.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/robust_predicates.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/geometric_predicates.f90
    ${fortfem_SOURCE_DIR}/src/triangulation/point_location.f90
    ${fortfem_SOURCE_DIR}/src/basis/basis_p1_2d.f90
    ${fortfem_SOURCE_DIR}/src/basis/basis_p2_2d.f90
    ${fortfem_SOURCE_DIR}/src/geometry/geometry_jet_2d.f90
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
target_link_libraries(kin6d_spatial PUBLIC kin6d_ode kin6d_intervals)
# The upstream BEM assembly (panel-pair colouring) is threaded with OpenMP.
if(KIN6D_ENABLE_OPENMP)
    target_link_libraries(kin6d_spatial PUBLIC OpenMP::OpenMP_Fortran)
endif()
set_source_files_properties(
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_reference_scalar_intervals.f90
    PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
target_sources(kin6d_fem2d PRIVATE src/fem2d/gs_load_batch.f90
    src/fem2d/gs_timers.f90 src/fem2d/gs_options.f90 src/fem2d/gs_csr.f90)
target_sources(kin6d_equilibrium PRIVATE src/generated/kin6d_grad_shafranov_leaves.f90
    src/generated/kin6d_grad_shafranov_batch.f90 src/equilibrium/grad_shafranov.f90)
if(KIN6D_ENABLE_OPENMP)
    target_link_libraries(kin6d_core PUBLIC OpenMP::OpenMP_Fortran)
endif()
add_library(kin6d_gs_reference_test STATIC
    test/generated/kin6d_gs_reference_error_leaves.f90
    test/support/gs_reference_error.f90)
target_link_libraries(kin6d_gs_reference_test PUBLIC kin6d)
target_compile_options(kin6d_gs_reference_test PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(kin6d_equilibrium PRIVATE src/generated/kin6d_gs_coarea_leaves.f90 src/equilibrium/gs_coarea.f90)
target_sources(kin6d_equilibrium PRIVATE src/generated/kin6d_gs_inverse_leaves.f90
    src/equilibrium/gs_inverse.f90 src/equilibrium/gs_prescribed_q.f90)
set_source_files_properties(src/generated/kin6d_gs_inverse_leaves.f90
    src/generated/kin6d_gs_coarea_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
set_source_files_properties(src/generated/kin6d_grad_shafranov_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d_fem2d PUBLIC kin6d_spatial)
add_executable(test_real_identity test/test_real_identity.f90)
target_link_libraries(test_real_identity PRIVATE kin6d)
add_test(NAME real_value_identity_edge_cases COMMAND test_real_identity)
add_executable(test_grad_shafranov test/test_grad_shafranov.f90)
target_link_libraries(test_grad_shafranov PRIVATE kin6d)
add_test(NAME grad_shafranov_analytic_and_derivative_oracles
    COMMAND test_grad_shafranov)
add_executable(test_gs_profile_quadrature test/test_gs_profile_quadrature.f90)
target_link_libraries(test_gs_profile_quadrature PRIVATE kin6d)
add_test(NAME grad_shafranov_profile_cut_quadrature_oracles
    COMMAND test_gs_profile_quadrature)
add_executable(test_gs_profile_domain test/test_gs_profile_domain.f90)
target_link_libraries(test_gs_profile_domain PRIVATE kin6d)
foreach(domain valid zero negative nan)
    add_test(NAME grad_shafranov_profile_mu0_${domain}
        COMMAND test_gs_profile_domain ${domain})
    add_test(NAME grad_shafranov_p2_profile_mu0_${domain}
        COMMAND test_gs_profile_domain ${domain} 2)
endforeach()
# Sparse direct factorisations go through FortSparse. Its AUTO backend picks
# per matrix between serial SuperLU (small systems) and threaded MUMPS (OpenMP,
# AMD or METIS ordering) for every general direct solve (GS stiffness Cholesky,
# ideal-map and coupled Newton LDL^T/LU); like OpenMP, MUMPS is required when
# enabled, so configure fails if it is not found.
# KIN6D_ENABLE_MUMPS=OFF selects serial SuperLU for the Newton systems and the
# iterative GS path. FortFEM's sparse_direct wrapper (SuperLU) stays the
# fallback for asymmetric systems.
# MUMPS is a sequential OpenMP scivision/mumps install (MUMPSConfig.cmake on
# CMAKE_PREFIX_PATH or MUMPS_DIR) or a distribution package (libdmumps,
# sequential or MPI). Detected by file, not find_package: fortsparse runs the
# config itself and the config may define its imported targets only once.
option(KIN6D_ENABLE_MUMPS "Direct sparse solves through fortsparse MUMPS" ON)
if(KIN6D_ENABLE_MUMPS)
    find_file(KIN6D_MUMPS_CONFIG NAMES MUMPSConfig.cmake
        HINTS ${MUMPS_DIR} ENV MUMPS_DIR
        PATH_SUFFIXES cmake lib/cmake lib/cmake/MUMPS)
    if(KIN6D_MUMPS_CONFIG)
        file(STRINGS "${KIN6D_MUMPS_CONFIG}" _kin6d_mumps_omp
            REGEX "set\\(MUMPS_openmp ON\\)")
        set(_kin6d_mumps_found "${KIN6D_MUMPS_CONFIG}")
        # The config resolves METIS/GKlib with plain find_library; without its
        # prefix first a system libmetis.so is mixed with the install's static
        # GKlib, which corrupts the heap inside METIS_NodeND.
        get_filename_component(_kin6d_mumps_prefix "${KIN6D_MUMPS_CONFIG}" DIRECTORY)
        foreach(_kin6d_level RANGE 2)
            get_filename_component(_kin6d_name "${_kin6d_mumps_prefix}" NAME)
            if(_kin6d_name MATCHES "^(MUMPS|cmake|lib|lib64)$")
                get_filename_component(_kin6d_mumps_prefix
                    "${_kin6d_mumps_prefix}" DIRECTORY)
            endif()
        endforeach()
        list(PREPEND CMAKE_PREFIX_PATH "${_kin6d_mumps_prefix}")
    else()
        find_library(KIN6D_MUMPS_LIBRARY NAMES dmumps
            HINTS ${MUMPS_DIR} ${MUMPS_ROOT} ENV MUMPS_DIR PATH_SUFFIXES lib lib64)
        if(NOT KIN6D_MUMPS_LIBRARY)
            message(FATAL_ERROR "MUMPS not found. Point MUMPS_DIR or "
                "CMAKE_PREFIX_PATH at an OpenMP MUMPS installation, or configure "
                "with -DKIN6D_ENABLE_MUMPS=OFF to use SuperLU.")
        endif()
        # An OpenMP build references the GNU OpenMP runtime entry points.
        file(STRINGS "${KIN6D_MUMPS_LIBRARY}" _kin6d_mumps_omp
            REGEX "^GOMP_parallel" LIMIT_COUNT 1)
        set(_kin6d_mumps_found "${KIN6D_MUMPS_LIBRARY}")
    endif()
    if(_kin6d_mumps_omp)
        message(STATUS "KIN6D sparse direct backend: MUMPS with OpenMP "
            "(${_kin6d_mumps_found})")
    else()
        message(WARNING "MUMPS at ${_kin6d_mumps_found} is not built with "
            "OpenMP: factorisations run on one thread. Prefer an OpenMP MUMPS "
            "via MUMPS_DIR.")
    endif()
    set(FORTSPARSE_ENABLE_MUMPS ON CACHE BOOL "" FORCE)
else()
    message(STATUS "KIN6D sparse direct backend: SuperLU (KIN6D_ENABLE_MUMPS=OFF)")
endif()
if(NOT TARGET fortsparse)
    set(FORTSPARSE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(FORTSPARSE_INSTALL_HELPER OFF CACHE BOOL "" FORCE)
    # FortSparse's fetched solvers share the selected BLAS; a prebuilt system
    # SuperLU may carry an incompatible reference-BLAS dependency.
    FetchContent_Declare(fortsparse
        GIT_REPOSITORY https://github.com/lazy-fortran/fortsparse.git
        GIT_TAG 09fc2b406aecf4f5ddda30e53aef8cc8bdb8b0c5
        GIT_SHALLOW FALSE)
    FetchContent_MakeAvailable(fortsparse)
endif()
if(KIN6D_ENABLE_OPENMP)
    # Block-RHS callers enter sparse solves concurrently even without MUMPS.
    target_link_libraries(fortsparse PRIVATE OpenMP::OpenMP_Fortran)
endif()
add_library(kin6d_sparse_direct
    ${fortfem_SOURCE_DIR}/src/topology/tree_cotree_gauge.f90
    ${fortfem_SOURCE_DIR}/src/solvers/fortfem_sparse_direct.f90)
target_link_libraries(kin6d_sparse_direct PUBLIC kin6d_spatial fortsparse)
target_link_libraries(kin6d_fem2d PUBLIC kin6d_sparse_direct)
if(KIN6D_ENABLE_MUMPS)
    target_sources(kin6d_fem2d PRIVATE src/fem2d/gs_direct.f90)
    # Follow the selected FortSparse backend, including sequential mpiseq
    # installations that do not provide a Fortran MPI module.
    set_source_files_properties(src/fem2d/gs_direct.f90 PROPERTIES Fortran_PREPROCESS ON)
    get_target_property(_kin6d_sparse_definitions fortsparse COMPILE_DEFINITIONS)
    if("FORTSPARSE_MUMPS_MPI" IN_LIST _kin6d_sparse_definitions)
        set_property(SOURCE src/fem2d/gs_direct.f90 APPEND PROPERTY
            COMPILE_DEFINITIONS KIN6D_MUMPS_MPI)
        set(_kin6d_mumps_mpi ON)
    endif()
else()
    target_sources(kin6d_fem2d PRIVATE src/fem2d/gs_direct_stub.f90)
endif()
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/gs_coils.f90 src/equilibrium/gs_coupling.f90
    src/equilibrium/gs_free_boundary.f90 src/equilibrium/gs_coil_inverse.f90)
# N1: nonlinear 3D frozen-in energy map (fixed boundary), DESIGN 3D lane.
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_ideal_map_leaves.f90
    src/mhd/ideal_map.f90 src/mhd/ideal_seed.f90)
set_source_files_properties(src/generated/kin6d_ideal_map_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
# N4: 3D free boundary of the map lane (Gamma_c shell, per-n exterior).
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_vacuum_shell_leaves.f90
    src/mhd/ideal_free_boundary.f90)
set_source_files_properties(src/generated/kin6d_vacuum_shell_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_landreman_leaves.f90
    src/mhd/landreman_benchmark.f90)
set_source_files_properties(src/generated/kin6d_landreman_leaves.f90 PROPERTIES
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_landreman_free_boundary test/test_landreman_free_boundary.f90)
target_link_libraries(test_landreman_free_boundary PRIVATE kin6d)
add_test(NAME ideal_fb_landreman_axisymmetric COMMAND test_landreman_free_boundary)
set_tests_properties(ideal_fb_landreman_axisymmetric PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_ideal_free_boundary test/test_ideal_free_boundary.f90)
target_link_libraries(test_ideal_free_boundary PRIVATE kin6d)
add_test(NAME ideal_fb_derivatives COMMAND test_ideal_free_boundary derivatives)
add_test(NAME ideal_fb_vertical_newton COMMAND test_ideal_free_boundary vertical)
foreach(k 2 3)
    add_test(NAME ideal_fb_vacuum_exact_n0_p${k} COMMAND test_ideal_free_boundary vacuum ${k})
endforeach()
add_test(NAME ideal_fb_vacuum_gamma_c_independence_n3
    COMMAND test_ideal_free_boundary vacuum3d 3)
add_test(NAME ideal_fb_slide_gauge_curl_free COMMAND test_ideal_free_boundary slide)
add_test(NAME ideal_fb_hypotheses_bordering COMMAND test_ideal_free_boundary hypotheses)
set_tests_properties(ideal_fb_derivatives ideal_fb_vertical_newton
    ideal_fb_slide_gauge_curl_free ideal_fb_hypotheses_bordering
    ideal_fb_vacuum_exact_n0_p2
    ideal_fb_vacuum_exact_n0_p3 ideal_fb_vacuum_gamma_c_independence_n3 PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_ideal_map_axisymmetric test/test_ideal_map_axisymmetric.f90
    test/support/ideal_map_cases.f90)
target_link_libraries(test_ideal_map_axisymmetric PRIVATE kin6d)
add_test(NAME ideal_map_derivatives_and_rotation_null
    COMMAND test_ideal_map_axisymmetric derivatives)
add_test(NAME ideal_map_n3_leakage COMMAND test_ideal_map_axisymmetric leakage)
foreach(k 2 3)
    add_test(NAME ideal_map_n0_solovev_a3_p${k}
        COMMAND test_ideal_map_axisymmetric solovev ${k})
    add_test(NAME ideal_map_n0_shaped_solovev_p${k}
        COMMAND test_ideal_map_axisymmetric shaped ${k})
endforeach()
# Three-mesh ladders: 6-15 s each with 4 threads; shaped P3 about 30 s with
# MUMPS (one factor thread, see factor_threads), about 1 min with SuperLU.
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
add_test(NAME ideal_map_n0_normalised_residual_a3_p2
    COMMAND test_ideal_map_axisymmetric residual 2)
set_tests_properties(ideal_map_n0_normalised_residual_a3_p2
    PROPERTIES TIMEOUT 900 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_seed_return_a3_p2
    PROPERTIES TIMEOUT 1800 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
set_tests_properties(ideal_map_n0_shaped_solovev_p3 PROPERTIES
    TIMEOUT 2400 PROCESSORS 4 ENVIRONMENT OMP_NUM_THREADS=4)
# gs_coupling and gs_free_boundary threads come from the project-wide KIN6D_ENABLE_OPENMP option.
add_executable(test_gs_coils test/test_gs_coils.f90)
target_link_libraries(test_gs_coils PRIVATE kin6d)
add_test(NAME gs_coils_independent_oracles COMMAND test_gs_coils)
add_library(kin6d_gs_coupling_test STATIC test/support/gs_coupling_fixture.f90)
target_link_libraries(kin6d_gs_coupling_test PUBLIC kin6d)
foreach(coupling_test vacuum fixed_current structure)
    add_executable(test_gs_coupling_${coupling_test} test/test_gs_coupling_${coupling_test}.f90)
    target_link_libraries(test_gs_coupling_${coupling_test} PRIVATE kin6d_gs_coupling_test)
endforeach()
add_test(NAME gs_coupling_t1_vacuum_coil_flux COMMAND test_gs_coupling_vacuum)
set_tests_properties(gs_coupling_t1_vacuum_coil_flux PROPERTIES TIMEOUT 900)
add_test(NAME gs_coupling_t2_fixed_current_cubature COMMAND test_gs_coupling_fixed_current)
set_tests_properties(gs_coupling_t2_fixed_current_cubature PROPERTIES TIMEOUT 900)
add_test(NAME gs_coupling_structure_symmetry_energy COMMAND test_gs_coupling_structure)
set_tests_properties(gs_coupling_structure_symmetry_energy PROPERTIES TIMEOUT 900)
add_executable(test_gs_steklov_cache test/test_gs_steklov_cache.f90)
target_link_libraries(test_gs_steklov_cache PRIVATE kin6d)
add_test(NAME gs_steklov_cache_relabelled_curve COMMAND test_gs_steklov_cache)
add_executable(test_gs_fb_interface test/test_gs_fb_interface.f90)
target_link_libraries(test_gs_fb_interface PRIVATE kin6d)
add_test(NAME gs_free_boundary_interface_flow_quadrature COMMAND test_gs_fb_interface)
add_library(kin6d_gs_fb_test STATIC test/support/gs_fb_fixture.f90)
target_link_libraries(kin6d_gs_fb_test PUBLIC kin6d_gs_coupling_test)
add_executable(test_gs_fb_jacobian test/test_gs_fb_jacobian.f90)
target_link_libraries(test_gs_fb_jacobian PRIVATE kin6d_gs_fb_test)
add_test(NAME gs_free_boundary_t7_jacobian_finite_differences COMMAND test_gs_fb_jacobian)
set_tests_properties(gs_free_boundary_t7_jacobian_finite_differences PROPERTIES TIMEOUT 900)
add_executable(test_gs_fb_manufactured test/test_gs_fb_manufactured.f90)
target_link_libraries(test_gs_fb_manufactured PRIVATE kin6d_gs_fb_test)
add_library(kin6d_gs_ci_test STATIC test/support/gs_ci_fixture.f90)
target_link_libraries(kin6d_gs_ci_test PUBLIC kin6d_gs_fb_test)
add_executable(test_gs_coil_inverse test/test_gs_coil_inverse.f90)
target_link_libraries(test_gs_coil_inverse PRIVATE kin6d_gs_ci_test)
add_test(NAME gs_coil_inverse_t6_twin COMMAND test_gs_coil_inverse)
set_tests_properties(gs_coil_inverse_t6_twin PROPERTIES TIMEOUT 1800)
add_executable(test_gs_fb_refit test/test_gs_fb_refit.f90)
target_link_libraries(test_gs_fb_refit PRIVATE kin6d_gs_fb_test)
add_test(NAME gs_free_boundary_p8_refit_loop COMMAND test_gs_fb_refit)
set_tests_properties(gs_free_boundary_p8_refit_loop PROPERTIES TIMEOUT 1800)
add_executable(test_gs_fb_greens test/test_gs_fb_greens.f90)
target_link_libraries(test_gs_fb_greens PRIVATE kin6d_gs_fb_test)
if(OpenMP_Fortran_FOUND)
    target_compile_options(test_gs_fb_greens PRIVATE ${OpenMP_Fortran_FLAGS})
endif()
add_test(NAME gs_free_boundary_t4_green_function_cross_check COMMAND test_gs_fb_greens)
set_tests_properties(gs_free_boundary_t4_green_function_cross_check PROPERTIES TIMEOUT 1800)
add_test(NAME gs_free_boundary_t3_manufactured_newton COMMAND test_gs_fb_manufactured)
set_tests_properties(gs_free_boundary_t3_manufactured_newton PROPERTIES TIMEOUT 1800)
add_executable(test_gs_coupled_majorant test/test_gs_coupled_majorant.f90)
target_link_libraries(test_gs_coupled_majorant PRIVATE kin6d_gs_fb_test)
add_test(NAME gs_free_boundary_p10_coupled_majorant COMMAND test_gs_coupled_majorant)
set_tests_properties(gs_free_boundary_p10_coupled_majorant PROPERTIES TIMEOUT 1800)
add_test(NAME gs_free_boundary_equilibrated_quick
    COMMAND test_gs_coupled_majorant quick)
set_tests_properties(gs_free_boundary_equilibrated_quick PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_bem_residual test/test_gs_bem_residual.f90)
target_link_libraries(test_gs_bem_residual PRIVATE kin6d_equilibrium)
add_test(NAME gs_bem_single_layer_reference_residual COMMAND test_gs_bem_residual)
set_tests_properties(gs_bem_single_layer_reference_residual PROPERTIES TIMEOUT 30
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_fb_large_aspect test/test_gs_fb_large_aspect.f90)
target_link_libraries(test_gs_fb_large_aspect PRIVATE kin6d)
add_test(NAME gs_free_boundary_t5_large_aspect_shafranov COMMAND test_gs_fb_large_aspect)
set_tests_properties(gs_free_boundary_t5_large_aspect_shafranov PROPERTIES TIMEOUT 1800)
target_sources(kin6d_fem2d PRIVATE src/fem2d/curved_triangle.f90
    src/generated/kin6d_segment_leaves.f90 src/fem2d/triangle_basis.f90)
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/gs_majorant.f90 src/equilibrium/gs_profile_law.f90 src/equilibrium/gs_nonlinear.f90
    src/equilibrium/gs_equilibrated.f90 src/generated/kin6d_gs_error_leaves.f90)
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/gs_quadrature_budget.f90)
set_source_files_properties(
    src/equilibrium/gs_quadrature_budget.f90
    src/equilibrium/gs_error_control.f90
    src/generated/kin6d_gs_error_leaves.f90
    src/generated/kin6d_grad_shafranov_leaves.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_rt0_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_rt0_majorant.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_verified_bernstein_products.f90
    ${fortfem_SOURCE_DIR}/src/estimators/verified_bernstein_majorant.f90
    PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
add_executable(test_gs_checked_majorant test/test_gs_checked_majorant.f90)
target_link_libraries(test_gs_checked_majorant PRIVATE kin6d_equilibrium)
add_test(NAME grad_shafranov_checked_bdm_majorant COMMAND test_gs_checked_majorant)
set_tests_properties(grad_shafranov_checked_bdm_majorant PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")
add_executable(test_gs_checked_source_remainder test/test_gs_checked_source_remainder.f90)
target_link_libraries(test_gs_checked_source_remainder PRIVATE kin6d_equilibrium)
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
set_property(SOURCE src/generated/kin6d_gs_error_leaves.f90 APPEND PROPERTY
    COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_gs_error_leaves test/test_gs_error_leaves.f90)
target_link_libraries(test_gs_error_leaves PRIVATE kin6d_equilibrium)
add_test(NAME grad_shafranov_error_leaf_oracles COMMAND test_gs_error_leaves)
set_tests_properties(grad_shafranov_error_leaf_oracles PROPERTIES TIMEOUT 20)
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/gs_error_control.f90)
add_executable(test_gs_error_control test/test_gs_error_control.f90)
target_link_libraries(test_gs_error_control PRIVATE kin6d)
add_test(NAME gs_fixed_stability_and_linear_goal COMMAND test_gs_error_control)
add_executable(test_gs_verified_stability test/test_gs_verified_stability.f90)
target_link_libraries(test_gs_verified_stability PRIVATE kin6d_equilibrium)
add_test(NAME gs_verified_nonlinear_coercivity COMMAND test_gs_verified_stability)
set_tests_properties(gs_verified_nonlinear_coercivity PROPERTIES
    TIMEOUT 30 PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
set_tests_properties(gs_fixed_stability_and_linear_goal PROPERTIES
    TIMEOUT 30 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_equilibrated_majorant test/test_gs_equilibrated_majorant.f90)
target_link_libraries(test_gs_equilibrated_majorant PRIVATE kin6d)
add_test(NAME grad_shafranov_equilibrated_flux_majorant
    COMMAND test_gs_equilibrated_majorant)
set_tests_properties(grad_shafranov_equilibrated_flux_majorant PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_gs_p2 test/test_gs_p2.f90)
target_link_libraries(test_gs_p2 PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_physical_oracles COMMAND test_gs_p2)
add_executable(test_gs_p2_gauge test/test_gs_p2_gauge.f90)
target_link_libraries(test_gs_p2_gauge PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_gauge_oracles COMMAND test_gs_p2_gauge)
add_executable(test_gs_p2_refinement test/test_gs_p2_refinement.f90)
target_link_libraries(test_gs_p2_refinement PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_polynomial_refinement COMMAND test_gs_p2_refinement)
set_tests_properties(grad_shafranov_p2_polynomial_refinement PROPERTIES
    FIXTURES_SETUP gs_p2_solovev_export)
add_executable(test_gs_coarea test/test_gs_coarea.f90)
target_link_libraries(test_gs_coarea PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_native_coarea_oracles COMMAND test_gs_coarea)
add_executable(test_gs_inverse_cache_reuse test/test_gs_inverse_cache_reuse.f90)
target_link_libraries(test_gs_inverse_cache_reuse PRIVATE kin6d)
add_test(NAME grad_shafranov_inverse_cache_reuse_oracles COMMAND test_gs_inverse_cache_reuse)
add_executable(test_gs_inverse_rhs_batch test/test_gs_inverse_rhs_batch.f90)
target_link_libraries(test_gs_inverse_rhs_batch PRIVATE kin6d)
add_test(NAME grad_shafranov_inverse_rhs_batch_oracles COMMAND test_gs_inverse_rhs_batch)
add_executable(test_gs_inverse_boundary_response test/test_gs_inverse_boundary_response.f90)
target_link_libraries(test_gs_inverse_boundary_response PRIVATE kin6d)
add_test(NAME grad_shafranov_inverse_boundary_response_oracles
    COMMAND test_gs_inverse_boundary_response)
add_executable(test_gs_inverse_cache_equation test/test_gs_inverse_cache_equation.f90)
target_link_libraries(test_gs_inverse_cache_equation PRIVATE kin6d)
add_test(NAME grad_shafranov_inverse_cache_equation_oracle COMMAND test_gs_inverse_cache_equation)
add_executable(test_gs_inverse_source_cut test/test_gs_inverse_source_cut.f90)
target_link_libraries(test_gs_inverse_source_cut PRIVATE kin6d)
add_test(NAME grad_shafranov_inverse_source_cut_oracles COMMAND test_gs_inverse_source_cut)
add_executable(test_gs_inverse_mean test/test_gs_inverse_mean.f90)
target_link_libraries(test_gs_inverse_mean PRIVATE kin6d)
add_test(NAME grad_shafranov_profile_mean_gauge_oracles COMMAND test_gs_inverse_mean)
add_executable(test_gs_inverse test/test_gs_inverse.f90)
target_link_libraries(test_gs_inverse PRIVATE kin6d)
add_executable(test_gs_prescribed_q test/test_gs_prescribed_q.f90)
target_link_libraries(test_gs_prescribed_q PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_p3_convergence COMMAND test_gs_prescribed_q)
set_tests_properties(grad_shafranov_prescribed_q_p3_convergence PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_variable_p test/test_gs_prescribed_q_variable_p.f90)
target_link_libraries(test_gs_prescribed_q_variable_p PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_variable_p COMMAND test_gs_prescribed_q_variable_p)
set_tests_properties(grad_shafranov_prescribed_q_variable_p PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_intervals test/test_gs_prescribed_q_intervals.f90)
target_link_libraries(test_gs_prescribed_q_intervals PRIVATE kin6d)
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
target_link_libraries(test_gs_prescribed_q_pedestal PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_pedestal COMMAND test_gs_prescribed_q_pedestal)
set_tests_properties(grad_shafranov_prescribed_q_pedestal PROPERTIES
    TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_edge_mode test/test_gs_prescribed_q_edge_mode.f90)
target_link_libraries(test_gs_prescribed_q_edge_mode PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_edge_mode
         COMMAND test_gs_prescribed_q_edge_mode)
set_tests_properties(grad_shafranov_prescribed_q_edge_mode PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_prescribed_q_circle test/test_gs_prescribed_q_circle.f90)
target_link_libraries(test_gs_prescribed_q_circle PRIVATE kin6d)
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
target_link_libraries(grad_shafranov_demo PRIVATE kin6d kin6d_io fortplot)
add_test(NAME grad_shafranov_inverse_app_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:grad_shafranov_demo> -DWORK=${CMAKE_BINARY_DIR}/gs_inverse_app_check
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_inverse_app.cmake)
set_tests_properties(grad_shafranov_inverse_app_contract PROPERTIES
    TIMEOUT 300 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_p1_weak_current test/test_gs_p1_weak_current.f90)
target_link_libraries(test_gs_p1_weak_current PRIVATE kin6d)
add_test(NAME grad_shafranov_p1_weak_current COMMAND test_gs_p1_weak_current)
foreach(mutation omit-jump omit-cylindrical wrong-current-sign)
    add_test(NAME grad_shafranov_p1_rejects_${mutation}
        COMMAND test_gs_p1_weak_current ${mutation})
    set_tests_properties(grad_shafranov_p1_rejects_${mutation} PROPERTIES WILL_FAIL TRUE)
endforeach()
add_executable(gold_hoyle_demo app/gold_hoyle_demo.f90)
add_executable(gs_reference_error app/gs_reference_error.f90)
target_link_libraries(gs_reference_error PRIVATE kin6d_gs_reference_test)
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
target_link_libraries(test_gs_reference_error PRIVATE kin6d_gs_reference_test)
add_test(NAME grad_shafranov_reference_cell_enclosures COMMAND test_gs_reference_error)
add_executable(test_gs_p2_reference_error test/test_gs_p2_reference_error.f90)
target_link_libraries(test_gs_p2_reference_error PRIVATE kin6d_gs_reference_test)
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

target_link_libraries(gold_hoyle_demo PRIVATE kin6d kin6d_io)
foreach(layer IN LISTS KIN6D_LAYER_LIBRARIES)
    target_compile_options(${layer} PRIVATE -Wall -Wextra -fno-fast-math
        -ffp-contract=off -frecursive)
endforeach()
add_executable(kin6d_demo app/kin6d_demo.f90)
find_package(OpenMP COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
    # Inverse Jacobian columns and other independent loops use OpenMP.
    target_link_libraries(kin6d_core PUBLIC OpenMP::OpenMP_Fortran)
endif()
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
# Presentation-only palette and number formatting for the gallery demos.
add_library(kin6d_demo_style app/demo_style.f90)
target_compile_options(kin6d_demo_style PRIVATE -Wall -Wextra)
foreach(application kin6d_demo steady_laboratory_demo transport_demo
        streaming_comparison_demo streaming_moments_demo debye_sphere_demo
        slab_oscillation_demo shielding_dynamics_demo interaction_demo
        microscopic_cloud_demo classical_pair_demo)
    target_link_libraries(${application} PRIVATE kin6d_demo_style)
endforeach()
add_executable(test_demo_style test/test_demo_style.f90)
target_link_libraries(test_demo_style PRIVATE kin6d_demo_style)
add_test(NAME demo_style_numbers COMMAND test_demo_style)
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
add_executable(test_pme test/test_pme.f90)
target_link_libraries(test_pme PRIVATE kin6d)
add_test(NAME pme_direct_ewald_oracles COMMAND test_pme)
add_executable(test_ocp test/test_ocp.f90)
target_link_libraries(test_ocp PRIVATE kin6d)
add_test(NAME ocp_vlasov_debye_oracles COMMAND test_ocp)
add_executable(test_cloud test/test_cloud.f90)
target_link_libraries(test_cloud PRIVATE kin6d)
add_test(NAME microscopic_cloud_oracles COMMAND test_cloud)
foreach(quantum_family fermi_gas quantum_fock quantum_evolution quantum_correlations
        quantum_reference_gallery classical_pair jellium_cell)
    add_executable(test_${quantum_family} test/test_${quantum_family}.f90)
    target_link_libraries(test_${quantum_family} PRIVATE kin6d)
    add_test(NAME ${quantum_family}_oracles COMMAND test_${quantum_family})
endforeach()
add_executable(pp_fusion_export app/pp_fusion_export.f90)
add_dependencies(pp_fusion_export kin6d_build_info)
target_include_directories(pp_fusion_export PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(pp_fusion_export PRIVATE kin6d kin6d_io)

add_test(NAME pp_fusion_export_contract
    COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:pp_fusion_export>
        -DWORK=${CMAKE_BINARY_DIR}/pp-fusion-export-test
        -P ${CMAKE_SOURCE_DIR}/test/check_pp_fusion_export.cmake)
set_tests_properties(pp_fusion_export_contract PROPERTIES PROCESSORS 4)

add_executable(pp_fusion_eft_export app/pp_fusion_eft_export.f90)
add_dependencies(pp_fusion_eft_export kin6d_build_info)
target_include_directories(pp_fusion_eft_export PRIVATE ${CMAKE_BINARY_DIR}/result_include)
target_link_libraries(pp_fusion_eft_export PRIVATE kin6d kin6d_io)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pp_fusion_eft_export_contract)
add_test(NAME pp_fusion_eft_export_contract
    COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:pp_fusion_eft_export>
        -DWORK=${CMAKE_BINARY_DIR}/pp_fusion_eft_export_contract
        -P ${CMAKE_SOURCE_DIR}/test/check_pp_fusion_eft_export.cmake)
set_tests_properties(pp_fusion_eft_export_contract PROPERTIES PROCESSORS 4)

add_executable(test_pp_fusion_eft test/test_pp_fusion_eft.f90)
target_link_libraries(test_pp_fusion_eft PRIVATE kin6d)
add_test(NAME pp_fusion_eft_oracles COMMAND test_pp_fusion_eft
    ${CMAKE_SOURCE_DIR}/test/reference/pp_fusion_eft.dat)
set_tests_properties(pp_fusion_eft_oracles PROPERTIES PROCESSORS 4)

add_executable(test_pp_fusion test/test_pp_fusion.f90)
target_link_libraries(test_pp_fusion PRIVATE kin6d)
add_test(NAME pp_fusion_oracles COMMAND test_pp_fusion
    ${CMAKE_SOURCE_DIR}/test/reference/pp_fusion_ode.dat)
set_tests_properties(pp_fusion_oracles PROPERTIES PROCESSORS 4)

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
set(KIN6D_FORTSYM_REVISION "435f8975ab5cec1624f8f76a26d9ea57e3b06195")
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
    add_executable(test_gs_p1_current_identities test/test_gs_p1_current_identities.f90)
    target_link_libraries(test_gs_p1_current_identities PRIVATE kin6d_gen_support)
    add_test(NAME grad_shafranov_p1_exact_current_antiderivatives
        COMMAND test_gs_p1_current_identities)
    add_executable(test_gs_symbolic test/test_gs_symbolic.f90)
    target_link_libraries(test_gs_symbolic PRIVATE kin6d_gen_support)
    add_test(NAME grad_shafranov_symbolic_identities COMMAND test_gs_symbolic)
    add_executable(test_gs_periodic_cylinder_symbolic
        test/test_gs_periodic_cylinder_symbolic.f90)
    target_link_libraries(test_gs_periodic_cylinder_symbolic PRIVATE kin6d_gen_support)
    add_test(NAME grad_shafranov_periodic_cylinder_exact_identities
        COMMAND test_gs_periodic_cylinder_symbolic)
    add_executable(test_gs_variational_identities test/test_gs_variational_identities.f90)
    target_link_libraries(test_gs_variational_identities PRIVATE kin6d_gen_support)
    add_executable(test_gs_forward_majorant_identities
        test/test_gs_forward_majorant_identities.f90)
    target_link_libraries(test_gs_forward_majorant_identities PRIVATE kin6d_gen_support)
    add_executable(test_gs_coupled_majorant_identities
        test/test_gs_coupled_majorant_identities.f90)
    target_link_libraries(test_gs_coupled_majorant_identities PRIVATE kin6d_gen_support)
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
    set(gs_leaf_source ${CMAKE_CURRENT_SOURCE_DIR}/src/generated/kin6d_grad_shafranov_leaves.f90)
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
    foreach(generator gs_error grad_shafranov ideal_map ideal_response_error
            equilibrium_core curve_chart nonnested vacuum_shell landreman
            response_strip_jets
            response_fluid
            gs_coarea gs_inverse
            gs_reference_error laboratory poisson fp transport
            streaming streaming_moments
            debye_sphere slab_oscillation shielding_dynamics interaction cloud
            fermi_gas quantum_fock quantum_evolution quantum_correlations
            quantum_reference_gallery classical_pair coulomb species coulomb_encounter
            pp_fusion pp_fusion_eft
            jellium_cell ocp)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE kin6d_gen_support)
    endforeach()
    add_test(NAME equilibrium_curve_chart_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_curve_chart>
        -DEXPECTED=${CMAKE_SOURCE_DIR}/src/generated/kin6d_curve_chart_leaves.f90
        -DWORK=${CMAKE_BINARY_DIR}/curve-chart-regeneration
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_check_curve_chart_regeneration.cmake)
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
add_executable(benchmark_ocp app/benchmark_ocp.f90)
target_link_libraries(benchmark_ocp PRIVATE kin6d)
# The periodic OCP demos run independent ensemble members on OpenMP threads;
# the library is reentrant (no module state) and compiled with automatic
# locals (-frecursive) so that concurrent calls do not share scratch.
foreach(example ocp_oscillation ocp_shielding)
    add_executable(${example}_demo app/${example}_demo.f90)
    target_link_libraries(${example}_demo PRIVATE kin6d fortplot)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(${example}_demo PRIVATE OpenMP::OpenMP_Fortran)
    endif()
endforeach()

# Campaign Q3/Q8: finite-N ideal fermion ensembles and relativistic Fermi EOS.
if(KIN6D_FORTSYM_BUILD)
    foreach(generator relativistic_fermi fermion_ensemble)
        add_executable(gen_${generator} gen/gen_${generator}.f90)
        target_link_libraries(gen_${generator} PRIVATE kin6d_gen_support)
    endforeach()
endif()
target_sources(kin6d_quantum PRIVATE src/generated/kin6d_relativistic_fermi_leaves.f90)
set_source_files_properties(src/generated/kin6d_relativistic_fermi_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_relativistic_fermi_gas test/test_relativistic_fermi_gas.f90)
target_link_libraries(test_relativistic_fermi_gas PRIVATE kin6d)
add_test(NAME relativistic_fermi_gas_oracles COMMAND test_relativistic_fermi_gas)
add_executable(test_fermi_gas_cold test/test_fermi_gas_cold.f90)
target_link_libraries(test_fermi_gas_cold PRIVATE kin6d)
add_test(NAME fermi_gas_cold_sommerfeld COMMAND test_fermi_gas_cold)
target_sources(kin6d_quantum PRIVATE src/generated/kin6d_fermion_ensemble_leaves.f90
    src/quantum/fermion_ensemble.f90)
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
# literal/resolved configuration and provenance out (src/apps/result_io.f90).
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
        -DKIN6D_FORTFEM_REVISION=${KIN6D_FORTFEM_REVISION}
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_build_info.cmake
    BYPRODUCTS ${kin6d_build_info_inc} VERBATIM)
add_library(kin6d_io src/apps/result_io.f90 src/apps/periodic_pair_config.f90
    src/apps/animation_export.f90)
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
        benchmark_quantum_cutoff jellium_cell_demo coulomb_encounter_demo
        benchmark_ocp ocp_oscillation_demo ocp_shielding_demo)
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
add_executable(test_animation_export test/test_animation_export.f90)
target_link_libraries(test_animation_export PRIVATE kin6d_io)
add_test(NAME animation_export_contract COMMAND test_animation_export
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# Q7: static screening hierarchy and 1X1V Wigner-Poisson (appended block).
add_library(kin6d_fft
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft8_kernel.f90
    ${fortnum_SOURCE_DIR}/src/fft/fortnum_fft.f90)
target_link_libraries(kin6d_fft PUBLIC kin6d_ode)
target_compile_options(kin6d_fft PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(kin6d_quantum PRIVATE
    src/generated/kin6d_screening_leaves.f90 src/quantum/screening_hierarchy.f90)
# Wigner-Poisson is a phase-space continuum solver sharing the streaming leaves.
target_sources(kin6d_continuum PRIVATE
    src/generated/kin6d_wigner_leaves.f90 src/continuum/wigner_dispersion.f90
    src/continuum/wigner_poisson.f90)
set_source_files_properties(src/generated/kin6d_screening_leaves.f90
    src/generated/kin6d_wigner_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d_core PUBLIC kin6d_fft)
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
target_sources(kin6d_quantum PRIVATE src/generated/kin6d_mean_field_cell_leaves.f90
    src/quantum/mean_field_cell.f90)
set_source_files_properties(src/generated/kin6d_mean_field_cell_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_mean_field_cell test/test_mean_field_cell.f90)
target_link_libraries(test_mean_field_cell PRIVATE kin6d)
add_test(NAME mean_field_cell_oracles COMMAND test_mean_field_cell)
add_executable(mean_field_cell_demo app/mean_field_cell_demo.f90)
target_link_libraries(mean_field_cell_demo PRIVATE kin6d kin6d_io fortplot)
# QG6: quantum binary-collision kernel and hybrid classical/quantum encounters
# (appended block).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_quantum_collision gen/gen_quantum_collision.f90)
    target_link_libraries(gen_quantum_collision PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d_pairs PRIVATE src/generated/kin6d_quantum_collision_leaves.f90
    src/pairs/quantum_collision_kernel.f90)
set_source_files_properties(src/generated/kin6d_quantum_collision_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_sources(kin6d_quantum PRIVATE src/quantum/hybrid_encounter.f90)
add_executable(test_quantum_collision test/test_quantum_collision.f90)
target_link_libraries(test_quantum_collision PRIVATE kin6d)
add_test(NAME quantum_collision_qg6_oracles COMMAND test_quantum_collision)
add_executable(test_hybrid_encounter test/test_hybrid_encounter.f90)
target_link_libraries(test_hybrid_encounter PRIVATE kin6d)
add_test(NAME hybrid_encounter_qg6_oracles COMMAND test_hybrid_encounter)
add_executable(hybrid_encounter_demo app/hybrid_encounter_demo.f90)
target_link_libraries(hybrid_encounter_demo PRIVATE kin6d kin6d_io fortplot)

# Campaign Q2: Gaussian/coherent-sum packet models measured against the Q1
# partial-wave reference (src/quantum/packet_models.f90, docs/packet-models.md).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_packet_models gen/gen_packet_models.f90)
    target_link_libraries(gen_packet_models PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d_quantum PRIVATE src/generated/kin6d_packet_models_leaves.f90
    src/quantum/gaussian_pair.f90 src/quantum/packet_models.f90)
set_source_files_properties(src/generated/kin6d_packet_models_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_packet_models test/test_packet_models.f90)
target_link_libraries(test_packet_models PRIVATE kin6d)
add_test(NAME packet_models_oracles COMMAND test_packet_models)
add_executable(packet_models_demo app/packet_models_demo.f90)
target_link_libraries(packet_models_demo PRIVATE kin6d fortplot kin6d_io)

# Campaign H1b: quantum vs classical finite-N pair/three-body correlations
# (docs/quantum-classical-correlations.md).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_quantum_classical gen/gen_quantum_classical.f90)
    target_link_libraries(gen_quantum_classical PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d_pairs PRIVATE src/generated/kin6d_quantum_classical_leaves.f90
    src/pairs/pair_law.f90)
target_sources(kin6d_classical PRIVATE src/classical/classical_cell.f90)
target_sources(kin6d_quantum PRIVATE src/quantum/quantum_pair_cell.f90
    src/quantum/quantum_classical_map.f90)
set_source_files_properties(src/generated/kin6d_quantum_classical_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_quantum_classical test/test_quantum_classical.f90)
target_link_libraries(test_quantum_classical PRIVATE kin6d)
add_test(NAME quantum_classical_h1b_oracles COMMAND test_quantum_classical)
add_executable(correlation_map_quantum_demo app/correlation_map_quantum_demo.f90)
target_link_libraries(correlation_map_quantum_demo PRIVATE kin6d kin6d_io
    kin6d_demo_style fortplot)
# Campaign H1a: classical OCP correlation-hierarchy map (Metropolis MC with
# incremental Ewald, g3/cumulant estimators, YBG closures, HNC; appended).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_correlation_hierarchy gen/gen_correlation_hierarchy.f90)
    target_link_libraries(gen_correlation_hierarchy PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d_core PRIVATE src/core/series_statistics.f90)
target_sources(kin6d_classical PRIVATE src/generated/kin6d_correlation_hierarchy_leaves.f90
    src/classical/ocp_monte_carlo.f90 src/classical/correlation_hierarchy.f90 src/classical/ocp_hnc.f90)
set_source_files_properties(src/generated/kin6d_correlation_hierarchy_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_correlation_hierarchy test/test_correlation_hierarchy.f90)
target_link_libraries(test_correlation_hierarchy PRIVATE kin6d)
add_test(NAME correlation_hierarchy_h1a_oracles COMMAND test_correlation_hierarchy)
add_executable(correlation_map_classical_demo app/correlation_map_types.f90
    app/correlation_map_md.f90 app/correlation_map_plots.f90
    app/correlation_map_classical_demo.f90)
target_link_libraries(correlation_map_classical_demo PRIVATE kin6d kin6d_io fortplot
    kin6d_demo_style)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(correlation_map_classical_demo PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Error-control surveys are explicit performance runs, never timing tests.
add_library(error_benchmark_support app/error_benchmark_support.f90)
foreach(survey continuum quantum microscopic modes)
    add_executable(benchmark_error_${survey} app/benchmark_error_${survey}.f90)
    target_link_libraries(benchmark_error_${survey} PRIVATE kin6d error_benchmark_support)
endforeach()
add_executable(test_error_modes test/test_error_modes.f90)
target_link_libraries(test_error_modes PRIVATE kin6d)
add_test(NAME error_modes_oracles COMMAND test_error_modes)
foreach(fixture verified_matvec sampled_norm)
    add_executable(test_error_${fixture} ${fortnum_SOURCE_DIR}/test/test_${fixture}.f90)
    target_link_libraries(test_error_${fixture} PRIVATE kin6d_ode)
    add_test(NAME error_${fixture}_oracles COMMAND test_error_${fixture})
endforeach()

add_executable(test_curved_triangle test/test_curved_triangle.f90)
target_link_libraries(test_curved_triangle PRIVATE kin6d)
add_test(NAME grad_shafranov_curved_geometry COMMAND test_curved_triangle)
add_executable(test_gs_curved_convergence test/test_gs_curved_convergence.f90)
target_link_libraries(test_gs_curved_convergence PRIVATE kin6d)
add_test(NAME grad_shafranov_curved_convergence COMMAND test_gs_curved_convergence)
add_test(NAME grad_shafranov_p3_curved_convergence COMMAND test_gs_curved_convergence 3)
add_executable(test_gs_p3 test/test_gs_p3.f90)
add_executable(test_gs_nonlinear test/test_gs_nonlinear.f90)
add_executable(test_gs_profile_newton test/test_gs_profile_newton.f90)
target_link_libraries(test_gs_profile_newton PRIVATE kin6d)
add_test(NAME grad_shafranov_cubic_profile_oracles COMMAND test_gs_profile_newton)
add_test(NAME grad_shafranov_profile_quadrature_work COMMAND test_gs_profile_newton work)
set_tests_properties(grad_shafranov_cubic_profile_oracles
    grad_shafranov_profile_quadrature_work PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_link_libraries(test_gs_nonlinear PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_nonlinear_convergence COMMAND test_gs_nonlinear 2)
add_test(NAME grad_shafranov_p3_nonlinear_convergence COMMAND test_gs_nonlinear 3)
add_test(NAME grad_shafranov_p3_spline_knot_convergence COMMAND test_gs_nonlinear 3 knots)
set_tests_properties(grad_shafranov_p2_nonlinear_convergence
    grad_shafranov_p3_nonlinear_convergence
    grad_shafranov_p3_spline_knot_convergence PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_link_libraries(test_gs_p3 PRIVATE kin6d)
add_test(NAME grad_shafranov_p3_physical_oracles COMMAND test_gs_p3)
add_executable(test_gs_threads test/test_gs_threads.f90)
target_link_libraries(test_gs_threads PRIVATE kin6d)
if(OpenMP_Fortran_FOUND)
    target_compile_options(test_gs_threads PRIVATE ${OpenMP_Fortran_FLAGS})
endif()
add_test(NAME grad_shafranov_p3_thread_invariance COMMAND test_gs_threads)
set_tests_properties(grad_shafranov_p3_thread_invariance PROPERTIES
    TIMEOUT 600 ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_direct test/test_gs_direct.f90)
target_link_libraries(test_gs_direct PRIVATE kin6d)
add_test(NAME grad_shafranov_csr_direct_newton_oracles COMMAND test_gs_direct)
if(_kin6d_mumps_mpi)
    add_executable(test_gs_direct_mpi test/test_gs_direct_mpi.f90)
    target_link_libraries(test_gs_direct_mpi PRIVATE kin6d)
    foreach(_mode automatic external finalized)
        add_test(NAME gs_direct_mpi_${_mode} COMMAND test_gs_direct_mpi ${_mode})
        set_tests_properties(gs_direct_mpi_${_mode} PROPERTIES
            TIMEOUT 15 ENVIRONMENT "OMP_NUM_THREADS=1")
    endforeach()
endif()
set_tests_properties(grad_shafranov_csr_direct_newton_oracles PROPERTIES
    SKIP_RETURN_CODE 77
    TIMEOUT 600 ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_solver_work test/test_gs_solver_work.f90)
target_link_libraries(test_gs_solver_work PRIVATE kin6d)
add_test(NAME grad_shafranov_shaped_solver_work COMMAND test_gs_solver_work)
add_executable(test_gs_axis_refinement test/test_gs_axis_refinement.f90)
target_link_libraries(test_gs_axis_refinement PRIVATE kin6d)
add_test(NAME grad_shafranov_axis_refinement COMMAND test_gs_axis_refinement)
add_executable(test_gs_profile_law_readback test/test_gs_profile_law_readback.f90)
target_link_libraries(test_gs_profile_law_readback PRIVATE kin6d)
add_test(NAME grad_shafranov_profile_law_readback
    COMMAND test_gs_profile_law_readback)
set_tests_properties(grad_shafranov_curved_convergence
    grad_shafranov_p3_curved_convergence grad_shafranov_axis_refinement
    grad_shafranov_shaped_solver_work PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/gs_readback.f90
    src/equilibrium/gs_axis_recovery.f90 src/equilibrium/gs_field_readback.f90
    src/equilibrium/gs_ampere_loop.f90)
add_executable(test_gs_curved_readback test/test_gs_curved_readback.f90)
target_link_libraries(test_gs_curved_readback PRIVATE kin6d)
add_test(NAME grad_shafranov_curved_readback COMMAND test_gs_curved_readback)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/p3_readback_test)
add_test(NAME grad_shafranov_p3_curved_readback COMMAND test_gs_curved_readback 3)
set_tests_properties(grad_shafranov_p3_curved_readback PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/p3_readback_test
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
set_tests_properties(grad_shafranov_curved_readback PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_executable(test_gs_ampere_loop test/test_gs_ampere_loop.f90)
target_link_libraries(test_gs_ampere_loop PRIVATE kin6d)
foreach(degree 2 3)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/ampere_p${degree}_test)
    add_test(NAME grad_shafranov_p${degree}_ampere_loop
        COMMAND test_gs_ampere_loop ${degree})
    set_tests_properties(grad_shafranov_p${degree}_ampere_loop PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/ampere_p${degree}_test
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_gs_weak_ampere_loop test/test_gs_weak_ampere_loop.f90)
target_link_libraries(test_gs_weak_ampere_loop PRIVATE kin6d)
foreach(degree 2 3)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/weak_ampere_p${degree}_test)
    add_test(NAME grad_shafranov_p${degree}_weak_ampere_loop
        COMMAND test_gs_weak_ampere_loop ${degree})
    set_tests_properties(grad_shafranov_p${degree}_weak_ampere_loop PROPERTIES
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/weak_ampere_p${degree}_test
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(gs_phase1 app/gs_phase1.f90)
target_link_libraries(gs_phase1 PRIVATE kin6d_io fortplot)
add_executable(gs_load_bench app/gs_load_bench.f90)
target_link_libraries(gs_load_bench PRIVATE kin6d)
add_executable(gs_free app/gs_free.f90)
target_link_libraries(gs_free PRIVATE kin6d_io fortplot)
add_executable(gs_free_3d app/gs_free_3d.f90)
target_link_libraries(gs_free_3d PRIVATE kin6d)
# NL-3a on the gs_free fixture (N = 0): reproduce the 2D free boundary and
# return from +-15 mm R/Z and +2 % elongation offset starts.
# One registered test per start case (see kin6d_test_labels.cmake).
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
target_link_libraries(test_gs_readback_order PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_readback_order COMMAND test_gs_readback_order 2)
add_test(NAME grad_shafranov_p3_readback_order COMMAND test_gs_readback_order 3)

add_executable(test_gs_readback_free test/test_gs_readback_free.f90)
target_link_libraries(test_gs_readback_free PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_readback_free COMMAND test_gs_readback_free 2)
add_test(NAME grad_shafranov_p3_readback_free COMMAND test_gs_readback_free 3)
set_tests_properties(grad_shafranov_p2_readback_free
    grad_shafranov_p3_readback_free PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
    TIMEOUT 300)

add_executable(test_gs_readback_work test/test_gs_readback_work.f90)
target_link_libraries(test_gs_readback_work PRIVATE kin6d)
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
target_sources(kin6d_mhd PRIVATE src/mhd/map_equilibrium.f90)
target_sources(kin6d_equilibrium PRIVATE src/generated/kin6d_equilibrium_core_leaves.f90
    src/generated/kin6d_curve_chart_leaves.f90
    src/equilibrium/equilibrium_chart.f90 src/equilibrium/equilibrium_state.f90 src/equilibrium/gs_equilibrium.f90
    src/equilibrium/equilibrium_residual.f90 src/equilibrium/equilibrium_readback.f90
    src/equilibrium/equilibrium_certificate.f90 src/equilibrium/gs_certificate.f90
    src/equilibrium/equilibrium_boundary.f90 src/equilibrium/equilibrium_geometry.f90)
set_source_files_properties(src/generated/kin6d_equilibrium_core_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_equilibrium_boundary test/test_equilibrium_boundary.f90)
target_link_libraries(test_equilibrium_boundary PRIVATE kin6d)
add_test(NAME equilibrium_core_d6_exterior_per_n
    COMMAND test_equilibrium_boundary exterior)
set_tests_properties(equilibrium_core_d6_exterior_per_n PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=4" PROCESSORS 4 TIMEOUT 600)
add_test(NAME equilibrium_core_d6_coil_interface
    COMMAND test_equilibrium_boundary coils)
add_executable(test_equilibrium_geometry test/test_equilibrium_geometry.f90)
target_link_libraries(test_equilibrium_geometry PRIVATE kin6d)
foreach(mode solovev circle drift state)
    add_test(NAME equilibrium_core_d7_geometry_${mode}
        COMMAND test_equilibrium_geometry ${mode})
    set_tests_properties(equilibrium_core_d7_geometry_${mode} PROPERTIES
        TIMEOUT 300 PROCESSORS 4
        ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()
add_executable(test_equilibrium_core test/test_equilibrium_core.f90)
target_link_libraries(test_equilibrium_core PRIVATE kin6d)
add_test(NAME equilibrium_core_chart_jet COMMAND test_equilibrium_core chart)
add_executable(test_curve_chart test/test_curve_chart.f90)
set_source_files_properties(test/test_curve_chart.f90 PROPERTIES Fortran_PREPROCESS ON)
target_link_libraries(test_curve_chart PRIVATE kin6d_equilibrium)
foreach(mode circle jets tube periodic invalid)
    add_test(NAME equilibrium_curve_chart_${mode} COMMAND test_curve_chart ${mode})
endforeach()
add_executable(test_curve_chart_solver test/test_curve_chart_solver.f90)
target_link_libraries(test_curve_chart_solver PRIVATE kin6d_mhd)
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
# K6-mini: the accepted TC24 K6 law and path (cubic-profile Newton, P3,
# direct factor, current constraint, readback, majorant) on a coarse mesh.
# Every core migration step must reproduce the fixture to 1e-12 relative.
add_test(NAME gs_k6_mini_regression COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_phase1> -DCHECK=$<TARGET_FILE:test_equilibrium_core>
    -DFIXTURES=${CMAKE_SOURCE_DIR}/test/fixtures/k6_mini
    -DWORK=${CMAKE_BINARY_DIR}/gs_k6_mini
    -P ${CMAKE_SOURCE_DIR}/test/gs_k6_mini.cmake)
set_tests_properties(gs_k6_mini_regression PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

# K-C1/L3: 3D polyline coils, per-n Fourier projection, vacuum b_mn (L3).
target_sources(kin6d_equilibrium PRIVATE src/equilibrium/coils_3d.f90)
target_sources(kin6d_mhd PRIVATE src/mhd/vacuum_surfaces.f90)
add_executable(test_coils_3d test/test_coils_3d.f90)
target_link_libraries(test_coils_3d PRIVATE kin6d)
add_test(NAME coils_3d_biot_savart_fourier_oracles COMMAND test_coils_3d)
set_tests_properties(coils_3d_biot_savart_fourier_oracles PROPERTIES TIMEOUT 120)
add_executable(test_vacuum_surfaces test/test_vacuum_surfaces.f90)
target_link_libraries(test_vacuum_surfaces PRIVATE kin6d)
add_test(NAME vacuum_surfaces_sfl_harmonics_oracles COMMAND test_vacuum_surfaces)
set_tests_properties(vacuum_surfaces_sfl_harmonics_oracles PROPERTIES TIMEOUT 300)
add_executable(coil_modes app/coil_modes.f90)
target_link_libraries(coil_modes PRIVATE kin6d)
add_executable(gs_vacuum_3d app/gs_vacuum_3d.f90)
target_link_libraries(gs_vacuum_3d PRIVATE kin6d)

# L4 (first stage): per-n FEEC complex, toroidal-mode BEM and their magnetic
# coupling from FortFEM; flux-aligned meshes and the per-n Ampere operator
# with the exact exterior.
add_library(kin6d_feec_n
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
    ${fortfem_SOURCE_DIR}/src/elements/triangle_nedelec_arbitrary_order.f90
    ${fortfem_SOURCE_DIR}/src/elements/triangle_rt_arbitrary_order.f90
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
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_scalar.f90
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_meridian_fourier_equilibration_products.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_equilibration.f90
    ${fortfem_SOURCE_DIR}/src/assembly/assembly_meridian_fourier_curl.f90
    ${fortfem_SOURCE_DIR}/src/bem/toroidal_mode_kernel.f90
    ${fortfem_SOURCE_DIR}/src/bem/toroidal_mode_boundary_operators_2d.f90
    ${fortfem_SOURCE_DIR}/src/bem/meridian_magnetic_coupling.f90)
target_link_libraries(kin6d_feec_n PUBLIC kin6d_spatial kin6d_fft fortsparse)
if(KIN6D_ENABLE_OPENMP)
    target_link_libraries(kin6d_feec_n PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_sources(kin6d_mhd PRIVATE src/mhd/mhd_mesh.f90 src/mhd/mhd_ampere_n.f90
    src/generated/kin6d_ideal_response_error_leaves.f90
    src/mhd/mhd_response_n.f90 src/mhd/mhd_ideal_n.f90 src/mhd/mhd_closure_n.f90
    src/mhd/mhd_layer_n.f90)
target_link_libraries(kin6d_fem2d PUBLIC kin6d_feec_n)
target_sources(kin6d_mhd PRIVATE src/mhd/mhd_perturbation_n.f90)
add_executable(test_mhd_perturbation test/test_mhd_perturbation.f90)
target_link_libraries(test_mhd_perturbation PRIVATE kin6d)
add_test(NAME mhd_perturbation_l7_xi_lagrangian_charts COMMAND test_mhd_perturbation)
set_tests_properties(mhd_perturbation_l7_xi_lagrangian_charts PROPERTIES TIMEOUT 300)
add_executable(test_mhd_ampere_error test/test_mhd_ampere_error.f90)
target_link_libraries(test_mhd_ampere_error PRIVATE kin6d)
add_test(NAME mhd_ampere_n_error_status COMMAND test_mhd_ampere_error)
set_tests_properties(mhd_ampere_n_error_status PROPERTIES TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ampere_sheet_error test/test_mhd_ampere_sheet_error.f90)
target_link_libraries(test_mhd_ampere_sheet_error PRIVATE kin6d)
add_test(NAME mhd_ampere_n_sheet_coupled_equilibrated
    COMMAND test_mhd_ampere_sheet_error)
set_tests_properties(mhd_ampere_n_sheet_coupled_equilibrated PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_meridian_edge_remainder
    ${fortfem_SOURCE_DIR}/test/test_meridian_edge_remainder.f90)
target_link_libraries(test_meridian_edge_remainder PRIVATE kin6d)
add_test(NAME mhd_ampere_n_sheet_interpolation_remainder
    COMMAND test_meridian_edge_remainder)
set_tests_properties(mhd_ampere_n_sheet_interpolation_remainder PROPERTIES TIMEOUT 30
    PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ampere_n test/test_mhd_ampere_n.f90)
target_link_libraries(test_mhd_ampere_n PRIVATE kin6d)
add_test(NAME mhd_ampere_n_exterior_manufactured COMMAND test_mhd_ampere_n)
set_tests_properties(mhd_ampere_n_exterior_manufactured PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_shielding test/test_mhd_shielding.f90)
target_link_libraries(test_mhd_shielding PRIVATE kin6d)
add_test(NAME mhd_shielding_cylinder_limit_reciprocity COMMAND test_mhd_shielding)
set_tests_properties(mhd_shielding_cylinder_limit_reciprocity PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(gs_response_3d app/gs_response_3d.f90)
target_link_libraries(gs_response_3d PRIVATE kin6d)
add_test(NAME gs_response_3d_ampere_indicator_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_response_3d> -DWORK=${CMAKE_BINARY_DIR}/gs_response_3d_check
    -DEQDSK=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_response_3d_app.cmake)
set_tests_properties(gs_response_3d_ampere_indicator_contract PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_executable(test_mhd_ideal_n test/test_mhd_ideal_n.f90)
target_link_libraries(test_mhd_ideal_n PRIVATE kin6d)
add_test(NAME mhd_ideal_n_mde_force_free_identity COMMAND test_mhd_ideal_n)
add_test(NAME mhd_ideal_n_final_state_error COMMAND test_mhd_ideal_n error)
set_tests_properties(mhd_ideal_n_final_state_error PROPERTIES TIMEOUT 60
    ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
add_test(NAME gs_response_3d_ideal_error_contract COMMAND ${CMAKE_COMMAND}
    -DAPP=$<TARGET_FILE:gs_response_3d>
    -DWORK=${CMAKE_BINARY_DIR}/gs_response_3d_ideal_check
    -DEQDSK=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk
    -P ${CMAKE_SOURCE_DIR}/test/check_gs_response_3d_ideal_error.cmake)
set_tests_properties(gs_response_3d_ideal_error_contract PROPERTIES TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
set_tests_properties(mhd_ideal_n_mde_force_free_identity PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_shield_torus test/test_mhd_shield_torus.f90)
target_link_libraries(test_mhd_shield_torus PRIVATE kin6d)
add_test(NAME mhd_shield_torus_nrad_convergence COMMAND test_mhd_shield_torus
    ${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/shield_torus_a10.eqdsk)
set_tests_properties(mhd_shield_torus_nrad_convergence PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_flux_mesh_axis test/test_flux_mesh_axis.f90)
target_link_libraries(test_flux_mesh_axis PRIVATE kin6d)
add_test(NAME flux_mesh_axis_interpolation_order COMMAND test_flux_mesh_axis)
set_tests_properties(flux_mesh_axis_interpolation_order PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(test_mhd_closure test/test_mhd_closure.f90)
target_link_libraries(test_mhd_closure PRIVATE kin6d)
add_test(NAME mhd_closure_tearing_delta_prime COMMAND test_mhd_closure)
set_tests_properties(mhd_closure_tearing_delta_prime PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
# R-B2 increment 1: strip coupled to the FEEC torus (response-unified.md, R-B2).
# Gated at nrad 64 (smallest passing nrad of the (b) ladder), strip half-width
# xd 0.01 plus an (a)-only run at xd 0.02 for the xd^2 order of gate (a)
# (~30 s with MUMPS, ~100 s with SuperLU; GMRES tol 1e-7 on both).
add_test(NAME response_strip_torus_coupling COMMAND test_mhd_closure strip)
set_tests_properties(response_strip_torus_coupling PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
# R-B1: resolved resistive layer strip, periodic cylinder (response-unified.md (B)).
target_sources(kin6d_mhd PRIVATE src/mhd/response_strip.f90)
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_response_fluid_leaves.f90
    src/mhd/response_fluid.f90 src/mhd/response_toroidal_strip.f90)
add_executable(test_response_toroidal_strip test/test_response_toroidal_strip.f90)
target_link_libraries(test_response_toroidal_strip PRIVATE kin6d)
add_executable(test_response_interfaces test/test_response_interfaces.f90)
target_link_libraries(test_response_interfaces PRIVATE kin6d)
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
target_link_libraries(test_response_fluid PRIVATE kin6d)
add_test(NAME response_fluid_full_pressure_oracles COMMAND test_response_fluid)
set_tests_properties(response_fluid_full_pressure_oracles PROPERTIES
    TIMEOUT 120 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OMP_WAIT_POLICY=PASSIVE")
add_executable(test_response_strip test/test_response_strip.f90)
target_link_libraries(test_response_strip PRIVATE kin6d)
add_test(NAME response_strip_cylinder_layer COMMAND test_response_strip)
set_tests_properties(response_strip_cylinder_layer PROPERTIES TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")
add_executable(response_strip_ladder app/response_strip_ladder.f90)
target_link_libraries(response_strip_ladder PRIVATE kin6d)
add_test(NAME mhd_resistive_layer_rotation COMMAND test_mhd_closure resistive)
set_tests_properties(mhd_resistive_layer_rotation PROPERTIES TIMEOUT 1800
    ENVIRONMENT "OMP_NUM_THREADS=4")
add_executable(test_mhd_newcomb test/test_mhd_newcomb.f90)
target_link_libraries(test_mhd_newcomb PRIVATE kin6d)
add_test(NAME mhd_newcomb_screw_pinch_limit COMMAND test_mhd_newcomb)
set_tests_properties(mhd_newcomb_screw_pinch_limit PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")

# NN1: closed free-B states and the vacuum / axisymmetric-pressure restrictions.
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_nonnested_leaves.f90
    src/mhd/equilibrium_nonnested.f90 src/mhd/nonnested_residual.f90)
set_source_files_properties(src/generated/kin6d_nonnested_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")

add_library(kin6d_nn_solovev_test STATIC test/support/nn_solovev_fixture.f90)
target_link_libraries(kin6d_nn_solovev_test PUBLIC kin6d)
add_executable(test_nonnested_solovev test/test_nonnested_solovev.f90)
target_link_libraries(test_nonnested_solovev PRIVATE kin6d_nn_solovev_test)
foreach(k 2 3)
    add_test(NAME nonnested_n0_solovev_p${k} COMMAND test_nonnested_solovev ${k})
    set_tests_properties(nonnested_n0_solovev_p${k} PROPERTIES TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
# Equilibrated energy bound of the N = 0 restriction against the Solov'ev flux.
target_sources(kin6d_mhd PRIVATE src/mhd/nonnested_error.f90)
add_executable(test_nonnested_error test/test_nonnested_error.f90)
target_link_libraries(test_nonnested_error PRIVATE kin6d_nn_solovev_test)
add_test(NAME nonnested_n0_equilibrated_solovev COMMAND test_nonnested_error)
set_tests_properties(nonnested_n0_equilibrated_solovev PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")
add_executable(test_nonnested_fourier_error
    ${fortfem_SOURCE_DIR}/test/test_meridian_fourier_equilibration.f90)
target_link_libraries(test_nonnested_fourier_error PRIVATE kin6d_feec_n)
add_test(NAME nonnested_transport_fourier_equilibrated
    COMMAND test_nonnested_fourier_error)
set_tests_properties(nonnested_transport_fourier_equilibrated PROPERTIES
    TIMEOUT 60 PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

add_executable(test_nonnested_transport_error test/test_nonnested_transport_error.f90)
target_link_libraries(test_nonnested_transport_error PRIVATE kin6d)
add_test(NAME nonnested_transport_equilibrated_n0
    COMMAND test_nonnested_transport_error)
set_tests_properties(nonnested_transport_equilibrated_n0 PROPERTIES TIMEOUT 60
    PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=2")

add_executable(test_nonnested_anisotropy test/test_nonnested_anisotropy.f90)
target_link_libraries(test_nonnested_anisotropy PRIVATE kin6d)
add_test(NAME nonnested_transport_anisotropy_pollution
    COMMAND test_nonnested_anisotropy pollution)
add_test(NAME nonnested_transport_anisotropy_correction
    COMMAND test_nonnested_anisotropy correction)
set_tests_properties(nonnested_transport_anisotropy_pollution
    nonnested_transport_anisotropy_correction PROPERTIES TIMEOUT 60
    PROCESSORS 4 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")

add_library(kin6d_nn_vacuum_test STATIC test/support/nn_vacuum_fixture.f90)
target_link_libraries(kin6d_nn_vacuum_test PUBLIC kin6d)

add_executable(nn1_spec_case test/nn1_spec_case.f90)
target_link_libraries(nn1_spec_case PRIVATE kin6d)
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
target_link_libraries(test_nonnested_vacuum PRIVATE kin6d_nn_vacuum_test)
add_test(NAME nonnested_vacuum_helical_coil_island COMMAND test_nonnested_vacuum)
set_tests_properties(nonnested_vacuum_helical_coil_island PROPERTIES TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=4")

add_executable(test_nonnested_leaves test/test_nonnested_leaves.f90)
target_link_libraries(test_nonnested_leaves PRIVATE kin6d)
add_test(NAME nonnested_symbolic_derivatives_and_d3_identity
    COMMAND test_nonnested_leaves)

# NN2 (first increment): anisotropic pressure transport for a given B.
target_sources(kin6d_mhd PRIVATE src/mhd/nonnested_transport.f90)
add_executable(test_nonnested_transport test/test_nonnested_transport.f90)
target_link_libraries(test_nonnested_transport PRIVATE kin6d)
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
target_link_libraries(test_nonnested_hudson_breslau PRIVATE kin6d)
add_test(NAME nonnested_hudson_breslau_integrable_exact_p3
    COMMAND test_nonnested_hudson_breslau integrable 3)
set_tests_properties(nonnested_hudson_breslau_integrable_exact_p3 PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")
add_test(NAME nonnested_hudson_breslau_chaos_responds_p3
    COMMAND test_nonnested_hudson_breslau chaos 3 7 48 3)
set_tests_properties(nonnested_hudson_breslau_chaos_responds_p3 PROPERTIES
    TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# NN2 (second increment): axisymmetric coupled transport / Ampere / Ohm Newton.
target_sources(kin6d_mhd PRIVATE src/mhd/nonnested_coupled.f90)
add_executable(test_nonnested_coupled test/test_nonnested_coupled.f90)
target_link_libraries(test_nonnested_coupled PRIVATE kin6d_nn_solovev_test)
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
target_link_libraries(test_nonnested_island PRIVATE kin6d_nn_vacuum_test)
add_test(NAME nonnested_island_flattening_p3
    COMMAND test_nonnested_island ratio 3 8 1)
add_test(NAME nonnested_island_d3_ladder_p2
    COMMAND test_nonnested_island ladder 2 4 1 1e4 3)
set_tests_properties(nonnested_island_flattening_p3 nonnested_island_d3_ladder_p2
    PROPERTIES TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4")
# Landreman family 1 (exact 3D, uniform iota) through the coupled solve.
add_executable(test_nonnested_landreman test/test_nonnested_landreman.f90)
target_link_libraries(test_nonnested_landreman PRIVATE kin6d)
add_test(NAME nonnested_coupled_landreman_family1_p2
    COMMAND test_nonnested_landreman 2 1 3 2)
set_tests_properties(nonnested_coupled_landreman_family1_p2
    PROPERTIES TIMEOUT 1800 ENVIRONMENT "OMP_NUM_THREADS=4")

# NN1 on a phi-dependent (mapped) domain: docs/nn1-mapped-3d.md.
target_sources(kin6d_mhd PRIVATE src/mhd/nonnested_mapped.f90)
add_executable(test_nonnested_mapped test/test_nonnested_mapped.f90)
target_link_libraries(test_nonnested_mapped PRIVATE kin6d)
add_test(NAME nonnested_mapped_identity_regression
    COMMAND test_nonnested_mapped identity)
add_test(NAME nonnested_mapped_ellipse_h_ladder
    COMMAND test_nonnested_mapped hladder 10 4 8 16)
add_test(NAME nonnested_mapped_ellipse_mode_ladder
    COMMAND test_nonnested_mapped modes 2 1 3)
set_tests_properties(nonnested_mapped_identity_regression
    nonnested_mapped_ellipse_h_ladder nonnested_mapped_ellipse_mode_ladder
    PROPERTIES TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4")

# Precise-QA vacuum benchmark of mapped NN1: docs/nn1-qa-benchmark.md.
add_executable(test_nonnested_qa test/test_nonnested_qa.f90)
target_link_libraries(test_nonnested_qa PRIVATE kin6d)
add_test(NAME nonnested_qa_oracle_check
    COMMAND test_nonnested_qa ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/nn1_qa check)
set_tests_properties(nonnested_qa_oracle_check
    PROPERTIES TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")

# Electromagnetic wave at a plasma edge: 1D Yee/Mur fields coupled to
# relativistic Boris macroparticles (docs/em-slab-reflection.md; appended).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_em_slab gen/gen_em_slab.f90)
    target_link_libraries(gen_em_slab PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d_continuum PRIVATE src/generated/kin6d_em_slab_leaves.f90
    src/continuum/em_slab.f90)
set_source_files_properties(src/generated/kin6d_em_slab_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_em_slab test/test_em_slab.f90)
target_link_libraries(test_em_slab PRIVATE kin6d)
foreach(em_case oracles vacuum reflection transmission refinement)
    add_test(NAME em_slab_${em_case} COMMAND test_em_slab ${em_case})
    set_tests_properties(em_slab_${em_case} PROPERTIES
        TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=4")
endforeach()
add_executable(em_slab_demo app/em_slab_demo.f90)
target_link_libraries(em_slab_demo PRIVATE kin6d kin6d_io)

# R-B2 jets: independent of outer coupling and cut-boundary orchestration.
target_sources(kin6d_mhd PRIVATE src/generated/kin6d_response_strip_jets_leaves.f90
    src/mhd/response_strip_jets.f90)
set_source_files_properties(src/generated/kin6d_response_strip_jets_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
add_executable(test_response_strip_jets test/test_response_strip_jets.f90)
target_link_libraries(test_response_strip_jets PRIVATE kin6d)
add_test(NAME response_strip_equilibrium_jets COMMAND test_response_strip_jets)
set_tests_properties(response_strip_equilibrium_jets PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1;OMP_WAIT_POLICY=passive"
    PROCESSORS 2 TIMEOUT 60)
set(KIN6D_STRIP_EQDSK_TEST "" CACHE FILEPATH
    "Optional imported equilibrium for independent strip-jet difference checks")
if(KIN6D_STRIP_EQDSK_TEST)
    add_test(NAME response_strip_tc24_jets
        COMMAND test_response_strip_jets ${KIN6D_STRIP_EQDSK_TEST})
    set_tests_properties(response_strip_tc24_jets PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1" TIMEOUT 60)
endif()
target_sources(kin6d_mhd PRIVATE src/mhd/response_strip_preconditioner.f90)
add_executable(test_response_strip_preconditioner
    test/test_response_strip_preconditioner.f90)
target_link_libraries(test_response_strip_preconditioner PRIVATE kin6d)
if(KIN6D_STRIP_EQDSK_TEST)
    add_executable(test_response_strip_tc24_preconditioner
        test/test_response_strip_tc24_preconditioner.f90)
    target_link_libraries(test_response_strip_tc24_preconditioner PRIVATE kin6d)
    add_test(NAME response_strip_tc24_preconditioner
        COMMAND test_response_strip_tc24_preconditioner ${KIN6D_STRIP_EQDSK_TEST})
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
# kin6d_check_<class> targets and the layering test (cmake/kin6d_layers.cmake).

add_executable(test_map_equilibrium test/test_map_equilibrium.f90)
target_link_libraries(test_map_equilibrium PRIVATE kin6d_mhd)
foreach(k 2 3)
    add_test(NAME equilibrium_core_map_state_p${k} COMMAND test_map_equilibrium ${k})
    set_tests_properties(equilibrium_core_map_state_p${k} PROPERTIES
        PROCESSORS 8 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
endforeach()

kin6d_layers_finalize()

# Test resource declaration. Every registered test gets an explicit OpenMP
# thread count, a matching PROCESSORS reservation, and OMP_WAIT_POLICY=passive.
# Without this, concurrently running tests each start one spinning OpenMP thread
# per core, oversubscribe the machine under `ctest -j`, and starve each other
# until a wall-clock limit trips. A test keeps the thread count it declared via
# ENVIRONMENT (or PROCESSORS); the rest get KIN6D_TEST_DEFAULT_THREADS.
set(KIN6D_TEST_DEFAULT_THREADS 2 CACHE STRING
    "OpenMP threads (and ctest PROCESSORS) for tests that do not declare their own")
get_directory_property(_kin6d_all_tests TESTS)
foreach(_t IN LISTS _kin6d_all_tests)
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
            set(_threads "${KIN6D_TEST_DEFAULT_THREADS}")
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

include(cmake/kin6d_estimator_benchmark.cmake)
