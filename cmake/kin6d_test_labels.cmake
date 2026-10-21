# Class and tier labels of every registered test (ROADMAP "Build structure").
# kin6d_test_labels(<test> <tier> CLASSES <class>... TARGETS <executable>...)
# Class: the highest layer whose modules the test executables use (mhd implies
# equilibrium); apps for application-only tests. Tier: slow above 60 s in a
# full run, fast otherwise (regression and performance are also accepted).
# kin6d_check_<class> builds the TARGETS of its tests. Tests that are not
# configured (FortSym, Lean) are skipped.
kin6d_test_labels(real_value_identity_edge_cases fast CLASSES core
    TARGETS test_real_identity)
kin6d_test_labels(grad_shafranov_analytic_and_derivative_oracles fast CLASSES equilibrium
    TARGETS test_grad_shafranov)
kin6d_test_labels(grad_shafranov_profile_cut_quadrature_oracles fast CLASSES equilibrium
    TARGETS test_gs_profile_quadrature)
kin6d_test_labels(grad_shafranov_profile_mu0_valid fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_p2_profile_mu0_valid fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_profile_mu0_zero fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_p2_profile_mu0_zero fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_profile_mu0_negative fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_p2_profile_mu0_negative fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_profile_mu0_nan fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(grad_shafranov_p2_profile_mu0_nan fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
kin6d_test_labels(ideal_fb_derivatives fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_vertical_newton fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_vacuum_exact_n0_p2 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_vacuum_exact_n0_p3 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_vacuum_gamma_c_independence_n3 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_slide_gauge_curl_free fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_fb_hypotheses_bordering fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
kin6d_test_labels(ideal_map_derivatives_and_rotation_null fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_curve_chart_landreman fast CLASSES mhd
    TARGETS test_curve_chart_solver)
kin6d_test_labels(ideal_map_curve_chart_figure8 fast CLASSES mhd
    TARGETS test_curve_chart_solver)
kin6d_test_labels(ideal_map_n3_leakage fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_solovev_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_shaped_solovev_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_solovev_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_shaped_solovev_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_return_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_return_shaped_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_return_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_reduced_omega_derivatives fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_nested_slide_null_p2 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_nested_slide_null_p3 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_nested_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_nested_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_seed_nested_shaped_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(ideal_map_n0_normalised_residual_a3_p2 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
kin6d_test_labels(gs_coils_independent_oracles fast CLASSES equilibrium
    TARGETS test_gs_coils)
kin6d_test_labels(gs_coupling_t1_vacuum_coil_flux fast CLASSES equilibrium
    TARGETS test_gs_coupling_vacuum)
kin6d_test_labels(gs_coupling_t2_fixed_current_cubature fast CLASSES equilibrium
    TARGETS test_gs_coupling_fixed_current)
kin6d_test_labels(gs_coupling_structure_symmetry_energy fast CLASSES equilibrium
    TARGETS test_gs_coupling_structure)
kin6d_test_labels(gs_steklov_cache_relabelled_curve fast CLASSES equilibrium
    TARGETS test_gs_steklov_cache)
kin6d_test_labels(gs_free_boundary_interface_flow_quadrature fast CLASSES equilibrium
    TARGETS test_gs_fb_interface)
kin6d_test_labels(gs_free_boundary_t7_jacobian_finite_differences fast CLASSES equilibrium
    TARGETS test_gs_fb_jacobian)
kin6d_test_labels(gs_coil_inverse_t6_twin fast CLASSES equilibrium
    TARGETS test_gs_coil_inverse)
kin6d_test_labels(gs_free_boundary_p8_refit_loop slow CLASSES equilibrium
    TARGETS test_gs_fb_refit)
kin6d_test_labels(gs_free_boundary_t4_green_function_cross_check fast CLASSES equilibrium
    TARGETS test_gs_fb_greens)
kin6d_test_labels(gs_free_boundary_t3_manufactured_newton fast CLASSES equilibrium
    TARGETS test_gs_fb_manufactured)
kin6d_test_labels(gs_free_boundary_p10_coupled_majorant slow CLASSES equilibrium
    TARGETS test_gs_coupled_majorant)
kin6d_test_labels(gs_free_boundary_t5_large_aspect_shafranov fast CLASSES equilibrium
    TARGETS test_gs_fb_large_aspect)
kin6d_test_labels(grad_shafranov_equilibrated_flux_majorant fast CLASSES equilibrium
    TARGETS test_gs_equilibrated_majorant)
kin6d_test_labels(grad_shafranov_checked_bdm_majorant fast CLASSES equilibrium
    TARGETS test_gs_checked_majorant)
kin6d_test_labels(grad_shafranov_checked_source_remainder fast CLASSES equilibrium
    TARGETS test_gs_checked_source_remainder)
kin6d_test_labels(grad_shafranov_checked_source_study slow CLASSES equilibrium
    TARGETS test_gs_checked_source_remainder)
kin6d_test_labels(grad_shafranov_p2_physical_oracles fast CLASSES equilibrium
    TARGETS test_gs_p2)
kin6d_test_labels(grad_shafranov_p2_gauge_oracles fast CLASSES equilibrium
    TARGETS test_gs_p2_gauge)
kin6d_test_labels(grad_shafranov_p2_polynomial_refinement fast CLASSES equilibrium
    TARGETS test_gs_p2_refinement)
kin6d_test_labels(grad_shafranov_p2_native_coarea_oracles fast CLASSES equilibrium
    TARGETS test_gs_coarea)
kin6d_test_labels(grad_shafranov_inverse_cache_reuse_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_cache_reuse)
kin6d_test_labels(grad_shafranov_inverse_rhs_batch_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_rhs_batch)
kin6d_test_labels(grad_shafranov_inverse_boundary_response_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_boundary_response)
kin6d_test_labels(grad_shafranov_inverse_cache_equation_oracle fast CLASSES equilibrium
    TARGETS test_gs_inverse_cache_equation)
kin6d_test_labels(grad_shafranov_inverse_source_cut_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_source_cut)
kin6d_test_labels(grad_shafranov_profile_mean_gauge_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_mean)
kin6d_test_labels(grad_shafranov_prescribed_q_p3_convergence fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q)
kin6d_test_labels(grad_shafranov_prescribed_q_variable_p fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_variable_p)
kin6d_test_labels(grad_shafranov_prescribed_q_jacobian_richardson fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
kin6d_test_labels(grad_shafranov_prescribed_q_exact_vs_secant fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
kin6d_test_labels(grad_shafranov_prescribed_q_knots64_uniform slow CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
kin6d_test_labels(grad_shafranov_prescribed_q_knots64_graded slow CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
kin6d_test_labels(grad_shafranov_prescribed_q_pedestal fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_pedestal)
kin6d_test_labels(grad_shafranov_prescribed_q_edge_mode fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_edge_mode)
kin6d_test_labels(grad_shafranov_prescribed_q_e1 fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_circle)
kin6d_test_labels(grad_shafranov_prescribed_q_e2 fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_circle)
kin6d_test_labels(grad_shafranov_phase1_inverse_app fast CLASSES equilibrium
    TARGETS gs_phase1 test_gs_prescribed_q_circle)
kin6d_test_labels(grad_shafranov_phase1_inverse_pressure_law_app fast CLASSES equilibrium
    TARGETS gs_phase1 test_gs_prescribed_q_circle)
kin6d_test_labels(grad_shafranov_constant_q_inverse_unit fast CLASSES equilibrium
    TARGETS test_gs_inverse)
kin6d_test_labels(grad_shafranov_constant_q_inverse_refinement fast CLASSES equilibrium
    TARGETS test_gs_inverse)
kin6d_test_labels(grad_shafranov_inverse_app_contract fast CLASSES equilibrium
    TARGETS grad_shafranov_demo)
kin6d_test_labels(grad_shafranov_p1_weak_current fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
kin6d_test_labels(grad_shafranov_p1_rejects_omit-jump fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
kin6d_test_labels(grad_shafranov_p1_rejects_omit-cylindrical fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
kin6d_test_labels(grad_shafranov_p1_rejects_wrong-current-sign fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
kin6d_test_labels(grad_shafranov_p1_reference_accepts_mesh_v1 fast CLASSES fem2d
    TARGETS gs_reference_error)
kin6d_test_labels(grad_shafranov_p1_reference_rejects_mesh_v2 fast CLASSES fem2d
    TARGETS gs_reference_error)
kin6d_test_labels(grad_shafranov_reference_cell_enclosures fast CLASSES fem2d
    TARGETS test_gs_reference_error)
kin6d_test_labels(grad_shafranov_p2_reference_cell_enclosures fast CLASSES fem2d
    TARGETS test_gs_p2_reference_error)
kin6d_test_labels(grad_shafranov_p2_reference_rejects_curved_metadata fast CLASSES fem2d
    TARGETS gs_reference_error)
kin6d_test_labels(demo_style_numbers fast CLASSES apps TARGETS test_demo_style)
kin6d_test_labels(characteristic_and_steady_oracles fast CLASSES core
    TARGETS test_laboratory)
kin6d_test_labels(poisson_and_coupled_steady_oracles fast CLASSES continuum
    TARGETS test_poisson)
kin6d_test_labels(bounded_fp_steady_oracles fast CLASSES continuum TARGETS test_fp)
kin6d_test_labels(general_transport_and_tangent_oracles fast CLASSES continuum
    TARGETS test_transport)
kin6d_test_labels(cached_transport_band_oracles fast CLASSES continuum
    TARGETS test_transport_sampling)
kin6d_test_labels(continuous_linear_band_oracles fast CLASSES continuum
    TARGETS test_transport_linear_bands)
kin6d_test_labels(streaming_representation_oracles fast CLASSES continuum
    TARGETS test_streaming_comparison)
kin6d_test_labels(streaming_particle_density_oracles fast CLASSES continuum
    TARGETS test_streaming_particle_density)
kin6d_test_labels(streaming_moments_oracles fast CLASSES continuum
    TARGETS test_streaming_moments)
kin6d_test_labels(debye_sphere_oracles fast CLASSES continuum TARGETS test_debye_sphere)
kin6d_test_labels(slab_oscillation_oracles fast CLASSES continuum
    TARGETS test_slab_oscillation)
kin6d_test_labels(shielding_dynamics_oracles fast CLASSES continuum
    TARGETS test_shielding_dynamics)
kin6d_test_labels(interaction_oracles fast CLASSES pairs TARGETS test_interaction)
kin6d_test_labels(nbody_coulomb_species_oracles fast CLASSES pairs TARGETS test_nbody)
kin6d_test_labels(pme_direct_ewald_oracles fast CLASSES pairs TARGETS test_pme)
kin6d_test_labels(ocp_vlasov_debye_oracles fast CLASSES classical continuum
    TARGETS test_ocp)
kin6d_test_labels(microscopic_cloud_oracles fast CLASSES classical TARGETS test_cloud)
kin6d_test_labels(fermi_gas_oracles fast CLASSES quantum TARGETS test_fermi_gas)
kin6d_test_labels(quantum_fock_oracles fast CLASSES quantum TARGETS test_quantum_fock)
kin6d_test_labels(quantum_evolution_oracles fast CLASSES quantum
    TARGETS test_quantum_evolution)
kin6d_test_labels(quantum_correlations_oracles fast CLASSES quantum
    TARGETS test_quantum_correlations)
kin6d_test_labels(quantum_reference_gallery_oracles fast CLASSES quantum
    TARGETS test_quantum_reference_gallery)
kin6d_test_labels(classical_pair_oracles fast CLASSES classical
    TARGETS test_classical_pair)
kin6d_test_labels(jellium_cell_oracles fast CLASSES quantum TARGETS test_jellium_cell)
kin6d_test_labels(coulomb_encounter_oracles fast CLASSES quantum
    TARGETS test_coulomb_encounter)
kin6d_test_labels(quantum_cutoff_oracles fast CLASSES quantum
    TARGETS test_quantum_cutoff)
kin6d_test_labels(relativistic_fermi_gas_oracles fast CLASSES quantum
    TARGETS test_relativistic_fermi_gas)
kin6d_test_labels(fermi_gas_cold_sommerfeld fast CLASSES quantum
    TARGETS test_fermi_gas_cold)
kin6d_test_labels(fermion_ensemble_oracles fast CLASSES quantum
    TARGETS test_fermion_ensemble)
kin6d_test_labels(result_io_contract fast CLASSES classical quantum
    TARGETS test_result_io classical_pair_demo quantum_reference_demo)
kin6d_test_labels(animation_export_contract fast CLASSES apps
    TARGETS test_animation_export)
kin6d_test_labels(screening_q7_oracles fast CLASSES quantum TARGETS test_screening)
kin6d_test_labels(wigner_q7_oracles fast CLASSES continuum TARGETS test_wigner)
kin6d_test_labels(mean_field_cell_oracles fast CLASSES quantum
    TARGETS test_mean_field_cell)
kin6d_test_labels(quantum_collision_qg6_oracles fast CLASSES pairs
    TARGETS test_quantum_collision)
kin6d_test_labels(hybrid_encounter_qg6_oracles fast CLASSES quantum
    TARGETS test_hybrid_encounter)
kin6d_test_labels(packet_models_oracles fast CLASSES quantum TARGETS test_packet_models)
kin6d_test_labels(quantum_classical_h1b_oracles fast CLASSES classical quantum
    TARGETS test_quantum_classical)
kin6d_test_labels(correlation_hierarchy_h1a_oracles fast CLASSES classical
    TARGETS test_correlation_hierarchy)
kin6d_test_labels(error_modes_oracles fast CLASSES continuum quantum
    TARGETS test_error_modes)
kin6d_test_labels(error_verified_matvec_oracles fast CLASSES core
    TARGETS test_error_verified_matvec)
kin6d_test_labels(error_sampled_norm_oracles fast CLASSES core
    TARGETS test_error_sampled_norm)
kin6d_test_labels(grad_shafranov_curved_geometry fast CLASSES fem2d
    TARGETS test_curved_triangle)
kin6d_test_labels(grad_shafranov_curved_convergence fast CLASSES equilibrium
    TARGETS test_gs_curved_convergence)
kin6d_test_labels(grad_shafranov_p3_curved_convergence fast CLASSES equilibrium
    TARGETS test_gs_curved_convergence)
kin6d_test_labels(grad_shafranov_cubic_profile_oracles fast CLASSES equilibrium
    TARGETS test_gs_profile_newton)
kin6d_test_labels(grad_shafranov_profile_quadrature_work fast CLASSES equilibrium
    TARGETS test_gs_profile_newton)
kin6d_test_labels(grad_shafranov_p2_nonlinear_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
kin6d_test_labels(grad_shafranov_p3_nonlinear_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
kin6d_test_labels(grad_shafranov_p3_spline_knot_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
kin6d_test_labels(grad_shafranov_p3_physical_oracles fast CLASSES equilibrium
    TARGETS test_gs_p3)
kin6d_test_labels(grad_shafranov_p3_thread_invariance fast CLASSES equilibrium
    TARGETS test_gs_threads)
kin6d_test_labels(grad_shafranov_csr_direct_newton_oracles fast CLASSES equilibrium
    TARGETS test_gs_direct)
foreach(_mode automatic external finalized)
    kin6d_test_labels(gs_direct_mpi_${_mode} fast CLASSES fem2d
        TARGETS test_gs_direct_mpi)
endforeach()
kin6d_test_labels(grad_shafranov_shaped_solver_work fast CLASSES equilibrium
    TARGETS test_gs_solver_work)
kin6d_test_labels(grad_shafranov_axis_refinement fast CLASSES equilibrium
    TARGETS test_gs_axis_refinement)
kin6d_test_labels(grad_shafranov_profile_law_readback fast CLASSES equilibrium
    TARGETS test_gs_profile_law_readback)
kin6d_test_labels(grad_shafranov_curved_readback fast CLASSES equilibrium
    TARGETS test_gs_curved_readback)
kin6d_test_labels(grad_shafranov_p3_curved_readback fast CLASSES equilibrium
    TARGETS test_gs_curved_readback)
foreach(case r0 rp rm zp zm el gzm)
    kin6d_test_labels(gs_free_3d_nl3a_fixture_return_${case} slow CLASSES mhd
        TARGETS gs_free_3d gs_free)
endforeach()
kin6d_test_labels(gs_free_3d_fixture_hypotheses slow CLASSES mhd
    TARGETS gs_free_3d gs_free)
kin6d_test_labels(gs_free_app_contract fast CLASSES equilibrium TARGETS gs_free)
kin6d_test_labels(grad_shafranov_p2_readback_order fast CLASSES equilibrium
    TARGETS test_gs_readback_order)
kin6d_test_labels(grad_shafranov_p3_readback_order fast CLASSES equilibrium
    TARGETS test_gs_readback_order)
kin6d_test_labels(grad_shafranov_p2_readback_free fast CLASSES equilibrium
    TARGETS test_gs_readback_free)
kin6d_test_labels(grad_shafranov_p3_readback_free fast CLASSES equilibrium
    TARGETS test_gs_readback_free)
kin6d_test_labels(grad_shafranov_readback_work fast CLASSES equilibrium
    TARGETS test_gs_readback_work)
kin6d_test_labels(grad_shafranov_stream_readback fast CLASSES equilibrium
    TARGETS gs_phase1)
kin6d_test_labels(equilibrium_core_d6_exterior_per_n fast CLASSES equilibrium
    TARGETS test_equilibrium_boundary)
kin6d_test_labels(equilibrium_core_d6_coil_interface fast CLASSES equilibrium
    TARGETS test_equilibrium_boundary)
kin6d_test_labels(equilibrium_core_d7_geometry_solovev fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
kin6d_test_labels(equilibrium_core_d7_geometry_circle fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
kin6d_test_labels(equilibrium_core_d7_geometry_drift fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
kin6d_test_labels(equilibrium_core_d7_geometry_state fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
kin6d_test_labels(equilibrium_core_chart_jet fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(equilibrium_core_d7_readback_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(equilibrium_core_d3_solovev_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(equilibrium_core_d3_upper_bound_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(gs_response_3d_ampere_indicator_contract fast CLASSES mhd
    TARGETS gs_response_3d)
kin6d_test_labels(equilibrium_core_d7_readback_p3 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(equilibrium_core_d3_solovev_p3 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
kin6d_test_labels(gs_k6_mini_regression slow CLASSES equilibrium
    TARGETS gs_phase1 test_equilibrium_core)
kin6d_test_labels(coils_3d_biot_savart_fourier_oracles fast CLASSES equilibrium
    TARGETS test_coils_3d)
kin6d_test_labels(vacuum_surfaces_sfl_harmonics_oracles fast CLASSES mhd
    TARGETS test_vacuum_surfaces)
kin6d_test_labels(mhd_perturbation_l7_xi_lagrangian_charts fast CLASSES mhd
    TARGETS test_mhd_perturbation)
kin6d_test_labels(mhd_ampere_n_error_status fast CLASSES mhd
    TARGETS test_mhd_ampere_error)
kin6d_test_labels(mhd_ampere_n_sheet_coupled_equilibrated fast CLASSES mhd
    TARGETS test_mhd_ampere_sheet_error)
kin6d_test_labels(mhd_ampere_n_sheet_interpolation_remainder fast CLASSES mhd
    TARGETS test_meridian_edge_remainder)
kin6d_test_labels(mhd_ampere_n_exterior_manufactured fast CLASSES mhd
    TARGETS test_mhd_ampere_n)
kin6d_test_labels(mhd_shielding_cylinder_limit_reciprocity fast CLASSES mhd
    TARGETS test_mhd_shielding)
kin6d_test_labels(mhd_ideal_n_mde_force_free_identity fast CLASSES mhd
    TARGETS test_mhd_ideal_n)
kin6d_test_labels(mhd_ideal_n_final_state_error fast CLASSES mhd
    TARGETS test_mhd_ideal_n)
kin6d_test_labels(gs_response_3d_ideal_error_contract fast CLASSES apps mhd
    TARGETS gs_response_3d)
kin6d_test_labels(mhd_shield_torus_nrad_convergence fast CLASSES mhd
    TARGETS test_mhd_shield_torus)
kin6d_test_labels(flux_mesh_axis_interpolation_order fast CLASSES mhd
    TARGETS test_flux_mesh_axis)
kin6d_test_labels(mhd_closure_tearing_delta_prime slow CLASSES mhd
    TARGETS test_mhd_closure)
kin6d_test_labels(response_strip_torus_coupling slow CLASSES mhd
    TARGETS test_mhd_closure)
kin6d_test_labels(response_toroidal_strip_two_ring_traces fast CLASSES mhd
    TARGETS test_response_interfaces)
kin6d_test_labels(response_cut_total_pressure_balance fast CLASSES mhd
    TARGETS test_response_interfaces)
kin6d_test_labels(response_toroidal_strip_edge_merged_domains fast CLASSES mhd
    TARGETS test_response_interfaces)
kin6d_test_labels(response_toroidal_strip_pressure_action fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
kin6d_test_labels(response_toroidal_strip_pressure_solve fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
kin6d_test_labels(response_toroidal_strip_radial_convergence fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
kin6d_test_labels(response_fluid_full_pressure_oracles fast CLASSES mhd
    TARGETS test_response_fluid)
kin6d_test_labels(response_strip_cylinder_layer fast CLASSES mhd
    TARGETS test_response_strip)
kin6d_test_labels(mhd_resistive_layer_rotation slow CLASSES mhd
    TARGETS test_mhd_closure)
kin6d_test_labels(mhd_newcomb_screw_pinch_limit fast CLASSES mhd
    TARGETS test_mhd_newcomb)
kin6d_test_labels(nonnested_n0_solovev_p2 fast CLASSES mhd
    TARGETS test_nonnested_solovev)
kin6d_test_labels(nonnested_n0_solovev_p3 fast CLASSES mhd
    TARGETS test_nonnested_solovev)
kin6d_test_labels(nonnested_n0_equilibrated_solovev fast CLASSES mhd
    TARGETS test_nonnested_error)
kin6d_test_labels(nonnested_transport_equilibrated_n0 fast CLASSES mhd
    TARGETS test_nonnested_transport_error)
kin6d_test_labels(nn1_spec_case_smoke fast CLASSES mhd TARGETS nn1_spec_case)
kin6d_test_labels(nonnested_vacuum_helical_coil_island slow CLASSES mhd
    TARGETS test_nonnested_vacuum)
kin6d_test_labels(nonnested_symbolic_derivatives_and_d3_identity fast CLASSES mhd
    TARGETS test_nonnested_leaves)
kin6d_test_labels(nonnested_transport_manufactured_3d slow CLASSES mhd
    TARGETS test_nonnested_transport)
kin6d_test_labels(nonnested_transport_sovinec_chi_ratio_p2 fast CLASSES mhd
    TARGETS test_nonnested_transport)
kin6d_test_labels(nonnested_transport_sovinec_chi_ratio_p3 fast CLASSES mhd
    TARGETS test_nonnested_transport)
kin6d_test_labels(nonnested_transport_sovinec_chi_ratio_p4 fast CLASSES mhd
    TARGETS test_nonnested_transport)
kin6d_test_labels(nonnested_hudson_breslau_integrable_exact_p3 fast CLASSES mhd
    TARGETS test_nonnested_hudson_breslau)
kin6d_test_labels(nonnested_hudson_breslau_chaos_responds_p3 fast CLASSES mhd
    TARGETS test_nonnested_hudson_breslau)
kin6d_test_labels(nonnested_coupled_solovev_p2 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_coupled_offset_return_p2 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_coupled_fourier_reduction_p2 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_coupled_solovev_p3 slow CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_coupled_offset_return_p3 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_coupled_fourier_reduction_p3 slow CLASSES mhd
    TARGETS test_nonnested_coupled)
kin6d_test_labels(nonnested_island_flattening_p3 slow CLASSES mhd
    TARGETS test_nonnested_island)
kin6d_test_labels(nonnested_island_d3_ladder_p2 slow CLASSES mhd
    TARGETS test_nonnested_island)
kin6d_test_labels(nonnested_coupled_landreman_family1_p2 slow CLASSES mhd
    TARGETS test_nonnested_landreman)
kin6d_test_labels(nonnested_mapped_identity_regression fast CLASSES mhd
    TARGETS test_nonnested_mapped)
kin6d_test_labels(nonnested_mapped_ellipse_h_ladder slow CLASSES mhd
    TARGETS test_nonnested_mapped)
kin6d_test_labels(nonnested_mapped_ellipse_mode_ladder fast CLASSES mhd
    TARGETS test_nonnested_mapped)
kin6d_test_labels(nonnested_qa_oracle_check fast CLASSES mhd TARGETS test_nonnested_qa)
kin6d_test_labels(lean_continuous_certificates fast CLASSES continuum
    TARGETS formal_krook_certificate)
kin6d_test_labels(grad_shafranov_p1_exact_current_antiderivatives fast CLASSES equilibrium
    TARGETS test_gs_p1_current_identities)
kin6d_test_labels(grad_shafranov_symbolic_identities fast CLASSES equilibrium
    TARGETS test_gs_symbolic)
kin6d_test_labels(grad_shafranov_periodic_cylinder_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_periodic_cylinder_symbolic)
kin6d_test_labels(grad_shafranov_coupled_majorant_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_coupled_majorant_identities)
kin6d_test_labels(grad_shafranov_forward_majorant_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
kin6d_test_labels(grad_shafranov_forward_majorant_rejects_jump_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
kin6d_test_labels(grad_shafranov_forward_majorant_rejects_boundary_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
kin6d_test_labels(grad_shafranov_forward_majorant_rejects_quadrature_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
kin6d_test_labels(grad_shafranov_variational_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)
kin6d_test_labels(grad_shafranov_rejects_wrong_field_orientation fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)
kin6d_test_labels(grad_shafranov_rejects_wrong_weak_source_weight fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)

kin6d_test_labels(grad_shafranov_error_leaf_oracles fast CLASSES equilibrium
    TARGETS test_gs_error_leaves)
kin6d_test_labels(gs_free_boundary_equilibrated_quick fast CLASSES equilibrium
    TARGETS test_gs_coupled_majorant)
kin6d_test_labels(gs_bem_single_layer_reference_residual fast CLASSES equilibrium
    TARGETS test_gs_bem_residual)
kin6d_test_labels(gs_fixed_stability_and_linear_goal fast CLASSES equilibrium
    TARGETS test_gs_error_control)
kin6d_test_labels(gs_verified_nonlinear_coercivity fast CLASSES equilibrium
    TARGETS test_gs_verified_stability)
kin6d_test_labels(grad_shafranov_production_estimator fast CLASSES apps
    TARGETS gs_phase1)

# Saved-field Ampere current diagnostics.
foreach(degree 2 3)
    kin6d_test_labels(grad_shafranov_p${degree}_ampere_loop fast CLASSES equilibrium
        TARGETS test_gs_ampere_loop)
    kin6d_test_labels(grad_shafranov_p${degree}_weak_ampere_loop fast CLASSES equilibrium
        TARGETS test_gs_weak_ampere_loop)
endforeach()

# Electromagnetic plane-wave teaching case.
foreach(mode oracles vacuum reflection transmission refinement)
    kin6d_test_labels(em_slab_${mode} fast CLASSES continuum TARGETS test_em_slab)
endforeach()

kin6d_test_labels(response_strip_equilibrium_jets fast CLASSES mhd
    TARGETS test_response_strip_jets)
kin6d_test_labels(response_strip_tc24_jets fast CLASSES mhd
    TARGETS test_response_strip_jets)
kin6d_test_labels(response_strip_tc24_preconditioner slow CLASSES mhd
    TARGETS test_response_strip_tc24_preconditioner)
foreach(check oracles negative limits)
    kin6d_test_labels(response_strip_preconditioner_${check} fast CLASSES mhd
        TARGETS test_response_strip_preconditioner)
endforeach()
kin6d_test_labels(response_strip_preconditioner_scaling performance CLASSES mhd
    TARGETS test_response_strip_preconditioner)

kin6d_test_labels(response_spot_calibration_algebra fast CLASSES mhd
    TARGETS test_response_calibration)
kin6d_test_labels(response_unit_sheet_spot_cylinder fast CLASSES mhd
    TARGETS test_response_sheet_cylinder)
kin6d_test_labels(response_qualification_frozen_rules fast CLASSES mhd
    TARGETS test_response_qualification)
kin6d_test_labels(response_qualification_external_adapter fast CLASSES apps
    TARGETS response_qualify)
kin6d_test_labels(response_spot_resolved_single_ring fast CLASSES mhd
    TARGETS test_response_strip_readback)
kin6d_test_labels(response_calibration_deck_validation fast CLASSES apps
    TARGETS gs_response_calibration)
kin6d_test_labels(response_calibration_deck_rejects_zero_mode fast CLASSES apps
    TARGETS gs_response_calibration)

kin6d_test_labels(nonnested_transport_anisotropy_pollution fast CLASSES mhd
    TARGETS test_nonnested_anisotropy)
kin6d_test_labels(nonnested_transport_anisotropy_correction fast CLASSES mhd
    TARGETS test_nonnested_anisotropy)
kin6d_test_labels(nonnested_transport_fourier_equilibrated fast CLASSES mhd
    TARGETS test_nonnested_fourier_error)

kin6d_test_labels(pp_fusion_oracles fast CLASSES quantum TARGETS test_pp_fusion)
kin6d_test_labels(pp_fusion_export_contract fast CLASSES quantum TARGETS pp_fusion_export)

kin6d_test_labels(pp_fusion_eft_oracles fast CLASSES quantum TARGETS test_pp_fusion_eft)
kin6d_test_labels(pp_fusion_eft_export_contract fast CLASSES quantum TARGETS pp_fusion_eft_export)

foreach(mode circle jets tube periodic invalid)
    kin6d_test_labels(equilibrium_curve_chart_${mode} fast CLASSES equilibrium
        TARGETS test_curve_chart)
endforeach()
if(TARGET gen_curve_chart)
    kin6d_test_labels(equilibrium_curve_chart_regeneration fast CLASSES equilibrium
        TARGETS gen_curve_chart)
endif()
kin6d_test_labels(ideal_fb_landreman_axisymmetric fast CLASSES mhd
    TARGETS test_landreman_free_boundary)

foreach(k 2 3)
    kin6d_test_labels(equilibrium_core_map_state_p${k} fast CLASSES mhd
        TARGETS test_map_equilibrium)
endforeach()
