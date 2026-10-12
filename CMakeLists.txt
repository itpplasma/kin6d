cmake_minimum_required(VERSION 3.22)
project(kin6d LANGUAGES Fortran C)
include(CTest)
include(FetchContent)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
# Dependency pins: one source used by FetchContent and by result provenance.
set(KIN6D_FORTPLOT_REVISION "759957d2b196fc7382ea1c370958f7adaa2fc753")
set(KIN6D_FORTNUM_REVISION "d053752a4160c9188c0a12c1e55ed4408788d1d4")
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
    ${fortnum_SOURCE_DIR}/src/special/fortnum_special_gamma.f90
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
    src/generated/kin6d_coulomb_leaves.f90 src/ewald.f90
    src/generated/kin6d_ocp_leaves.f90 src/near_encounter.f90 src/pme.f90
    src/nbody.f90 src/ocp.f90 src/vlasov_response.f90 src/pair_correlation.f90
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
set(KIN6D_FORTFEM_REVISION "de4245c5b15a87424b53b754c3189f6a25f43fc7")
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
    ${fortfem_SOURCE_DIR}/src/quadrature/triangle_level_quadrature.f90)
target_link_libraries(kin6d_spatial PUBLIC kin6d_ode kin6d_intervals)
set_source_files_properties(
    ${fortfem_SOURCE_DIR}/src/generated/fortfem_reference_scalar_intervals.f90
    PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
target_sources(kin6d PRIVATE src/generated/kin6d_grad_shafranov_leaves.f90
    src/grad_shafranov.f90)
add_library(kin6d_gs_reference_test STATIC
    test/generated/kin6d_gs_reference_error_leaves.f90
    test/support/gs_reference_error.f90)
target_link_libraries(kin6d_gs_reference_test PUBLIC kin6d)
target_compile_options(kin6d_gs_reference_test PRIVATE -fno-fast-math -ffp-contract=off)
target_sources(kin6d PRIVATE src/generated/kin6d_gs_coarea_leaves.f90 src/gs_coarea.f90)
target_sources(kin6d PRIVATE src/generated/kin6d_gs_inverse_leaves.f90
    src/gs_inverse.f90 src/gs_prescribed_q.f90)
set_source_files_properties(src/generated/kin6d_gs_inverse_leaves.f90
    src/generated/kin6d_gs_coarea_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
set_source_files_properties(src/generated/kin6d_grad_shafranov_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_link_libraries(kin6d PUBLIC kin6d_spatial)
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
target_sources(kin6d PRIVATE src/curved_triangle.f90 src/triangle_basis.f90
    src/gs_majorant.f90 src/gs_profile_law.f90 src/gs_nonlinear.f90)
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
set_tests_properties(grad_shafranov_prescribed_q_p3_convergence PROPERTIES TIMEOUT 300)
add_executable(test_gs_prescribed_q_variable_p test/test_gs_prescribed_q_variable_p.f90)
target_link_libraries(test_gs_prescribed_q_variable_p PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_variable_p COMMAND test_gs_prescribed_q_variable_p)
set_tests_properties(grad_shafranov_prescribed_q_variable_p PROPERTIES TIMEOUT 900)
add_executable(test_gs_prescribed_q_circle test/test_gs_prescribed_q_circle.f90)
target_link_libraries(test_gs_prescribed_q_circle PRIVATE kin6d)
add_test(NAME grad_shafranov_prescribed_q_e1 COMMAND test_gs_prescribed_q_circle 1)
add_test(NAME grad_shafranov_prescribed_q_e2 COMMAND test_gs_prescribed_q_circle 2)
set_tests_properties(grad_shafranov_prescribed_q_e1 grad_shafranov_prescribed_q_e2
    PROPERTIES TIMEOUT 900)
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
target_compile_options(kin6d PRIVATE -Wall -Wextra -fno-fast-math -ffp-contract=off
    -frecursive)
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
set(KIN6D_FORTSYM_REVISION "f749e65a4c00f014d88d0a4786a0689662b6a28f")
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
    foreach(generator grad_shafranov gs_coarea gs_inverse gs_reference_error
            laboratory poisson fp transport
            streaming streaming_moments
            debye_sphere slab_oscillation shielding_dynamics interaction cloud
            fermi_gas quantum_fock quantum_evolution quantum_correlations
            quantum_reference_gallery classical_pair coulomb species coulomb_encounter
            jellium_cell ocp)
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
add_executable(benchmark_ocp app/benchmark_ocp.f90)
target_link_libraries(benchmark_ocp PRIVATE kin6d)
# The periodic OCP demos run independent ensemble members on OpenMP threads;
# the library is reentrant (no module state) and compiled with automatic
# locals (-frecursive) so that concurrent calls do not share scratch.
find_package(OpenMP COMPONENTS Fortran)
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
        -DKIN6D_FORTFEM_REVISION=${KIN6D_FORTFEM_REVISION}
        -P ${CMAKE_SOURCE_DIR}/cmake/kin6d_build_info.cmake
    BYPRODUCTS ${kin6d_build_info_inc} VERBATIM)
add_library(kin6d_io src/result_io.f90 src/periodic_pair_config.f90
    src/animation_export.f90)
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
# QG6: quantum binary-collision kernel and hybrid classical/quantum encounters
# (appended block).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_quantum_collision gen/gen_quantum_collision.f90)
    target_link_libraries(gen_quantum_collision PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d PRIVATE src/generated/kin6d_quantum_collision_leaves.f90
    src/quantum_collision_kernel.f90)
set_source_files_properties(src/generated/kin6d_quantum_collision_leaves.f90
    PROPERTIES COMPILE_OPTIONS "-Wno-unused-dummy-argument")
target_sources(kin6d PRIVATE src/hybrid_encounter.f90)
add_executable(test_quantum_collision test/test_quantum_collision.f90)
target_link_libraries(test_quantum_collision PRIVATE kin6d)
add_test(NAME quantum_collision_qg6_oracles COMMAND test_quantum_collision)
add_executable(test_hybrid_encounter test/test_hybrid_encounter.f90)
target_link_libraries(test_hybrid_encounter PRIVATE kin6d)
add_test(NAME hybrid_encounter_qg6_oracles COMMAND test_hybrid_encounter)
add_executable(hybrid_encounter_demo app/hybrid_encounter_demo.f90)
target_link_libraries(hybrid_encounter_demo PRIVATE kin6d kin6d_io fortplot)

# Campaign Q2: Gaussian/coherent-sum packet models measured against the Q1
# partial-wave reference (src/packet_models.f90, docs/packet-models.md).
if(KIN6D_FORTSYM_BUILD)
    add_executable(gen_packet_models gen/gen_packet_models.f90)
    target_link_libraries(gen_packet_models PRIVATE kin6d_gen_support)
endif()
target_sources(kin6d PRIVATE src/generated/kin6d_packet_models_leaves.f90
    src/gaussian_pair.f90 src/packet_models.f90)
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
target_sources(kin6d PRIVATE src/generated/kin6d_quantum_classical_leaves.f90
    src/pair_law.f90 src/classical_cell.f90 src/quantum_pair_cell.f90
    src/quantum_classical_map.f90)
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
target_sources(kin6d PRIVATE src/generated/kin6d_correlation_hierarchy_leaves.f90
    src/ocp_monte_carlo.f90 src/series_statistics.f90 src/correlation_hierarchy.f90
    src/ocp_hnc.f90)
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
add_executable(test_gs_solver_work test/test_gs_solver_work.f90)
target_link_libraries(test_gs_solver_work PRIVATE kin6d)
add_test(NAME grad_shafranov_shaped_solver_work COMMAND test_gs_solver_work)
add_executable(test_gs_axis_refinement test/test_gs_axis_refinement.f90)
target_link_libraries(test_gs_axis_refinement PRIVATE kin6d)
add_test(NAME grad_shafranov_axis_refinement COMMAND test_gs_axis_refinement)
set_tests_properties(grad_shafranov_curved_convergence
    grad_shafranov_p3_curved_convergence grad_shafranov_axis_refinement
    grad_shafranov_shaped_solver_work PROPERTIES
    TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
target_sources(kin6d PRIVATE src/gs_readback.f90 src/gs_axis_recovery.f90
    src/gs_field_readback.f90)
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
add_executable(gs_phase1 app/gs_phase1.f90)
target_link_libraries(gs_phase1 PRIVATE kin6d_io fortplot)
# The comparison producer is optional; the native solver needs no Python.

add_executable(test_gs_readback_order test/test_gs_readback_order.f90)
target_link_libraries(test_gs_readback_order PRIVATE kin6d)
add_test(NAME grad_shafranov_p2_readback_order COMMAND test_gs_readback_order 2)
add_test(NAME grad_shafranov_p3_readback_order COMMAND test_gs_readback_order 3)

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
