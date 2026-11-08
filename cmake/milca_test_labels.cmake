foreach(k 2 3)
    foreach(geometry cylinder circle helix)
        milca_test_labels(ideal_fb_chart_surface_${geometry}_p${k} fast
            CLASSES mhd equilibrium TARGETS test_ideal_fb_surface)
    endforeach()
endforeach()
milca_test_labels(ideal_fb_chart_surface_invalid fast CLASSES mhd equilibrium
    TARGETS test_ideal_fb_surface)
milca_test_labels(ideal_fb_chart_surface_circulant fast CLASSES mhd equilibrium
    TARGETS test_ideal_fb_surface)
milca_test_labels(ideal_fb_exterior_frozen_derivatives fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_circuits_derivatives fast CLASSES mhd
    TARGETS test_ideal_fb_circuits)
milca_test_labels(ideal_fb_circuits_exterior3d fast CLASSES mhd
    TARGETS test_ideal_fb_circuits)
milca_test_labels(ideal_fb_traces3d_axisymmetric fast CLASSES mhd
    TARGETS test_ideal_fb_traces3d)
milca_test_labels(ideal_fb_traces3d_helical fast CLASSES mhd
    TARGETS test_ideal_fb_traces3d)
milca_test_labels(ideal_fb_traces3d_makegrid fast CLASSES mhd
    TARGETS test_ideal_fb_traces3d)
milca_test_labels(ideal_fb_traces3d_curve_hessian fast CLASSES mhd
    TARGETS test_ideal_fb_traces3d)
milca_test_labels(ideal_fb_exterior3d_periodic slow CLASSES mhd
    TARGETS test_ideal_fb_traces3d)
milca_test_labels(ideal_fb_landreman_axisymmetric_circuits fast CLASSES mhd
    TARGETS test_landreman_free_boundary)
milca_test_labels(ideal_fb_landreman_helical_circuits fast CLASSES mhd
    TARGETS test_landreman_free_boundary)
foreach(k 2 3)
    foreach(case fixed floating voltage)
        milca_test_labels(nonnested_actions_${case}_p${k} fast CLASSES mhd
            TARGETS test_nonnested_actions)
    endforeach()
endforeach()
milca_test_labels(nonnested_actions_invalid fast CLASSES mhd
    TARGETS test_nonnested_actions)
# Class and tier labels of every registered test (ROADMAP "Build structure").
foreach(case waves rates outer newcomb convergence saturation table)
    milca_test_labels(rm_oracle_${case} fast CLASSES mhd TARGETS test_rm_oracles)
endforeach()
milca_test_labels(rm_oracle_regeneration fast CLASSES mhd TARGETS gen_rm_oracles)
# milca_test_labels(<test> <tier> CLASSES <class>... TARGETS <executable>...)
# Class: the highest layer whose modules the test executables use (mhd implies
# equilibrium); apps for application-only tests. Tier: slow above 60 s in a
# full run, fast otherwise (regression and performance are also accepted).
# milca_check_<class> builds the TARGETS of its tests. Tests that are not
# configured (FortSym, Lean) are skipped.
milca_test_labels(grad_shafranov_device_selection_contract fast CLASSES equilibrium
    TARGETS test_gs_device_selection)
milca_test_labels(grad_shafranov_cuda_public_fixed_source fast CLASSES equilibrium
    TARGETS test_gs_cuda_public_solve)
milca_test_labels(grad_shafranov_cuda_affine_picard fast CLASSES equilibrium
    TARGETS test_gs_cuda_picard)
milca_test_labels(grad_shafranov_affine_picard_cpu fast CLASSES equilibrium
    TARGETS test_gs_cuda_picard)
milca_test_labels(grad_shafranov_cuda_fixed_source_a3 fast CLASSES equilibrium
    TARGETS test_gs_cuda_fixed_source)
milca_test_labels(grad_shafranov_cuda_affine_strict_audit fast CLASSES equilibrium
    TARGETS test_gs_cuda_fixed_source)
milca_test_labels(gs_free_boundary_freegs_public_parity fast CLASSES equilibrium
    TARGETS test_gs_freegs)
milca_test_labels(s0_alfven_wave_convergence fast CLASSES mhd TARGETS test_s0)
milca_test_labels(s0_shear_continuum_targeting fast CLASSES mhd TARGETS test_s0)
milca_test_labels(s0_physical_metric_contracts fast CLASSES mhd TARGETS test_s0)
milca_test_labels(s0_exact_identities_and_regeneration fast CLASSES mhd TARGETS gen_s0)
milca_test_labels(s0_gs_space_regeneration fast CLASSES mhd TARGETS gen_s0_gs_space)
milca_test_labels(s0_gs_space_jets_and_derivatives fast CLASSES mhd
    TARGETS test_s0_gs_space)
milca_test_labels(s0_gs_kink_radial_refinement fast CLASSES mhd TARGETS test_s0_gs_kink)
milca_test_labels(s0_gs_kink_serial_refinement fast CLASSES mhd TARGETS test_s0_gs_kink)
milca_test_labels(s0_gs_kink_aspect_ratio fast CLASSES mhd TARGETS test_s0_gs_kink)
milca_test_labels(s0_gs_kink_equilibrium_refinement fast CLASSES mhd
    TARGETS test_s0_gs_kink)
milca_test_labels(s0_cylinder_oracle_limits fast CLASSES mhd TARGETS test_s0_gs_kink)
milca_test_labels(real_value_identity_edge_cases fast CLASSES core
    TARGETS test_real_identity)
milca_test_labels(gs_mirror_mesh fast CLASSES equilibrium
    TARGETS test_gs_mirror_mesh)
milca_test_labels(grad_shafranov_analytic_and_derivative_oracles fast CLASSES equilibrium
    TARGETS test_grad_shafranov)
milca_test_labels(grad_shafranov_profile_cut_quadrature_oracles fast CLASSES equilibrium
    TARGETS test_gs_profile_quadrature)
milca_test_labels(grad_shafranov_profile_mu0_valid fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_p2_profile_mu0_valid fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_profile_mu0_zero fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_p2_profile_mu0_zero fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_profile_mu0_negative fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_p2_profile_mu0_negative fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_profile_mu0_nan fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(grad_shafranov_p2_profile_mu0_nan fast CLASSES equilibrium
    TARGETS test_gs_profile_domain)
milca_test_labels(ideal_fb_schur_spectrum fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_external_kink fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_derivatives fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_vertical_newton fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_vacuum_exact_n0_p2 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_vacuum_exact_n0_p3 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_vacuum_gamma_c_independence_n3 fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_slide_gauge_curl_free fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_hypotheses_bordering fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_fb_negative_curvature fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(ideal_map_derivatives_and_rotation_null fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_factor_cache fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_masked_factor_step fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_krylov_step fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_lagged_hessian fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_continuation fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_hessian_sum_factorisation_oracle fast CLASSES mhd
    TARGETS test_ideal_map_hessian)
foreach(degree RANGE 4 6)
    milca_test_labels(ideal_map_hessian_p${degree}_oracle fast CLASSES mhd
        TARGETS test_ideal_map_hessian)
endforeach()
milca_test_labels(ideal_map_curve_chart_landreman fast CLASSES mhd
    TARGETS test_curve_chart_solver)
milca_test_labels(ideal_map_curve_chart_figure8 fast CLASSES mhd
    TARGETS test_curve_chart_solver)
milca_test_labels(ideal_map_n3_leakage fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_solovev_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_shaped_solovev_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_solovev_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_triangle_lagrange_p2_p8 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_shaped_solovev_p4 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_shaped_solovev_p5 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_shaped_solovev_p6 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_shaped_solovev_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_return_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_return_shaped_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_return_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_reduced_omega_derivatives fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_nested_slide_null_p2 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_nested_slide_null_p3 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_nested_a3_p2 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_nested_a3_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_n0_seed_nested_shaped_p3 slow CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(ideal_map_parallel_mass fast CLASSES mhd
    TARGETS test_ideal_map_parallel)
milca_test_labels(ideal_map_parallel_damping_p2 fast CLASSES mhd
    TARGETS test_ideal_map_parallel)
milca_test_labels(ideal_map_n0_normalised_residual_a3_p2 fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)
milca_test_labels(gs_coils_independent_oracles fast CLASSES equilibrium
    TARGETS test_gs_coils)
milca_test_labels(gs_coupling_t1_vacuum_coil_flux fast CLASSES equilibrium
    TARGETS test_gs_coupling_vacuum)
milca_test_labels(gs_coupling_t2_fixed_current_cubature fast CLASSES equilibrium
    TARGETS test_gs_coupling_fixed_current)
milca_test_labels(gs_coupling_structure_symmetry_energy fast CLASSES equilibrium
    TARGETS test_gs_coupling_structure)
milca_test_labels(gs_steklov_cache_relabelled_curve fast CLASSES equilibrium
    TARGETS test_gs_steklov_cache)
milca_test_labels(gs_free_boundary_interface_flow_quadrature fast CLASSES equilibrium
    TARGETS test_gs_fb_interface)
milca_test_labels(gs_free_boundary_t7_jacobian_finite_differences fast CLASSES equilibrium
    TARGETS test_gs_fb_jacobian)
milca_test_labels(gs_free_boundary_law_tangent_finite_differences fast CLASSES equilibrium
    TARGETS test_gs_fb_law_tangent)
milca_test_labels(gs_coil_inverse_t6_twin fast CLASSES equilibrium
    TARGETS test_gs_coil_inverse)
milca_test_labels(gs_free_boundary_p8_refit_loop slow CLASSES equilibrium
    TARGETS test_gs_fb_refit)
milca_test_labels(gs_free_boundary_t4_green_function_cross_check fast CLASSES equilibrium
    TARGETS test_gs_fb_greens)
milca_test_labels(gs_free_boundary_t3_manufactured_newton fast CLASSES equilibrium
    TARGETS test_gs_fb_manufactured)
milca_test_labels(gs_free_boundary_p10_coupled_majorant slow CLASSES equilibrium
    TARGETS test_gs_coupled_majorant)
milca_test_labels(gs_free_boundary_t5_large_aspect_shafranov fast CLASSES equilibrium
    TARGETS test_gs_fb_large_aspect)
milca_test_labels(grad_shafranov_equilibrated_flux_majorant fast CLASSES equilibrium
    TARGETS test_gs_equilibrated_majorant)
milca_test_labels(grad_shafranov_checked_bdm_majorant fast CLASSES equilibrium
    TARGETS test_gs_checked_majorant)
milca_test_labels(grad_shafranov_checked_source_remainder fast CLASSES equilibrium
    TARGETS test_gs_checked_source_remainder)
milca_test_labels(grad_shafranov_checked_source_study slow CLASSES equilibrium
    TARGETS test_gs_checked_source_remainder)
milca_test_labels(grad_shafranov_p2_physical_oracles fast CLASSES equilibrium
    TARGETS test_gs_p2)
milca_test_labels(grad_shafranov_p2_gauge_oracles fast CLASSES equilibrium
    TARGETS test_gs_p2_gauge)
milca_test_labels(grad_shafranov_p2_polynomial_refinement fast CLASSES equilibrium
    TARGETS test_gs_p2_refinement)
milca_test_labels(grad_shafranov_p2_native_coarea_oracles fast CLASSES equilibrium
    TARGETS test_gs_coarea)
milca_test_labels(grad_shafranov_inverse_cache_reuse_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_cache_reuse)
milca_test_labels(grad_shafranov_inverse_rhs_batch_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_rhs_batch)
milca_test_labels(grad_shafranov_inverse_boundary_response_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_boundary_response)
milca_test_labels(grad_shafranov_inverse_cache_equation_oracle fast CLASSES equilibrium
    TARGETS test_gs_inverse_cache_equation)
milca_test_labels(grad_shafranov_inverse_source_cut_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_source_cut)
milca_test_labels(grad_shafranov_profile_mean_gauge_oracles fast CLASSES equilibrium
    TARGETS test_gs_inverse_mean)
milca_test_labels(grad_shafranov_prescribed_q_p3_convergence fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q)
milca_test_labels(grad_shafranov_prescribed_q_variable_p fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_variable_p)
milca_test_labels(grad_shafranov_prescribed_q_jacobian_richardson fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
milca_test_labels(grad_shafranov_prescribed_q_exact_vs_secant fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
milca_test_labels(grad_shafranov_prescribed_q_knots64_uniform slow CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
milca_test_labels(grad_shafranov_prescribed_q_knots64_graded slow CLASSES equilibrium
    TARGETS test_gs_prescribed_q_intervals)
milca_test_labels(grad_shafranov_prescribed_q_pedestal fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_pedestal)
milca_test_labels(grad_shafranov_prescribed_q_edge_mode fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_edge_mode)
milca_test_labels(grad_shafranov_prescribed_q_e1 fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_circle)
milca_test_labels(grad_shafranov_prescribed_q_e2 fast CLASSES equilibrium
    TARGETS test_gs_prescribed_q_circle)
milca_test_labels(grad_shafranov_phase1_inverse_app fast CLASSES equilibrium
    TARGETS gs_phase1 test_gs_prescribed_q_circle)
milca_test_labels(grad_shafranov_phase1_inverse_pressure_law_app fast CLASSES equilibrium
    TARGETS gs_phase1 test_gs_prescribed_q_circle)
milca_test_labels(grad_shafranov_constant_q_inverse_unit fast CLASSES equilibrium
    TARGETS test_gs_inverse)
milca_test_labels(grad_shafranov_constant_q_inverse_refinement fast CLASSES equilibrium
    TARGETS test_gs_inverse)
milca_test_labels(grad_shafranov_inverse_app_contract fast CLASSES equilibrium
    TARGETS grad_shafranov_demo)
milca_test_labels(grad_shafranov_p1_weak_current fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
milca_test_labels(grad_shafranov_p1_rejects_omit-jump fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
milca_test_labels(grad_shafranov_p1_rejects_omit-cylindrical fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
milca_test_labels(grad_shafranov_p1_rejects_wrong-current-sign fast CLASSES equilibrium
    TARGETS test_gs_p1_weak_current)
milca_test_labels(grad_shafranov_p1_reference_accepts_mesh_v1 fast CLASSES fem2d
    TARGETS gs_reference_error)
milca_test_labels(grad_shafranov_p1_reference_rejects_mesh_v2 fast CLASSES fem2d
    TARGETS gs_reference_error)
milca_test_labels(grad_shafranov_reference_cell_enclosures fast CLASSES fem2d
    TARGETS test_gs_reference_error)
milca_test_labels(grad_shafranov_p2_reference_cell_enclosures fast CLASSES fem2d
    TARGETS test_gs_p2_reference_error)
milca_test_labels(grad_shafranov_p2_reference_rejects_curved_metadata fast CLASSES fem2d
    TARGETS gs_reference_error)
milca_test_labels(demo_style_numbers fast CLASSES apps TARGETS test_demo_style)
milca_test_labels(characteristic_and_steady_oracles fast CLASSES core
    TARGETS test_laboratory)
milca_test_labels(poisson_and_coupled_steady_oracles fast CLASSES continuum
    TARGETS test_poisson)
milca_test_labels(bounded_fp_steady_oracles fast CLASSES continuum TARGETS test_fp)
milca_test_labels(general_transport_and_tangent_oracles fast CLASSES continuum
    TARGETS test_transport)
milca_test_labels(cached_transport_band_oracles fast CLASSES continuum
    TARGETS test_transport_sampling)
milca_test_labels(continuous_linear_band_oracles fast CLASSES continuum
    TARGETS test_transport_linear_bands)
milca_test_labels(streaming_representation_oracles fast CLASSES continuum
    TARGETS test_streaming_comparison)
milca_test_labels(streaming_particle_density_oracles fast CLASSES continuum
    TARGETS test_streaming_particle_density)
milca_test_labels(streaming_moments_oracles fast CLASSES continuum
    TARGETS test_streaming_moments)
milca_test_labels(debye_sphere_oracles fast CLASSES continuum TARGETS test_debye_sphere)
foreach(check maxwellian l1_brute invariants slow_limit)
    milca_test_labels(rosenbluth_${check}_oracle fast CLASSES continuum
        TARGETS test_rosenbluth_reference)
endforeach()
milca_test_labels(slab_oscillation_oracles fast CLASSES continuum
    TARGETS test_slab_oscillation)
milca_test_labels(shielding_dynamics_oracles fast CLASSES continuum
    TARGETS test_shielding_dynamics)
milca_test_labels(interaction_oracles fast CLASSES pairs TARGETS test_interaction)
milca_test_labels(nbody_coulomb_species_oracles fast CLASSES pairs TARGETS test_nbody)
milca_test_labels(pme_direct_ewald_oracles fast CLASSES pairs TARGETS test_pme)
milca_test_labels(ocp_vlasov_debye_oracles fast CLASSES classical continuum
    TARGETS test_ocp)
milca_test_labels(microscopic_cloud_oracles fast CLASSES classical TARGETS test_cloud)
milca_test_labels(fermi_gas_oracles fast CLASSES quantum TARGETS test_fermi_gas)
milca_test_labels(quantum_fock_oracles fast CLASSES quantum TARGETS test_quantum_fock)
milca_test_labels(quantum_evolution_oracles fast CLASSES quantum
    TARGETS test_quantum_evolution)
milca_test_labels(quantum_correlations_oracles fast CLASSES quantum
    TARGETS test_quantum_correlations)
milca_test_labels(quantum_reference_gallery_oracles fast CLASSES quantum
    TARGETS test_quantum_reference_gallery)
milca_test_labels(classical_pair_oracles fast CLASSES classical
    TARGETS test_classical_pair)
milca_test_labels(jellium_cell_oracles fast CLASSES quantum TARGETS test_jellium_cell)
milca_test_labels(coulomb_encounter_oracles fast CLASSES quantum
    TARGETS test_coulomb_encounter)
milca_test_labels(quantum_cutoff_oracles fast CLASSES quantum
    TARGETS test_quantum_cutoff)
milca_test_labels(relativistic_fermi_gas_oracles fast CLASSES quantum
    TARGETS test_relativistic_fermi_gas)
milca_test_labels(fermi_gas_cold_sommerfeld fast CLASSES quantum
    TARGETS test_fermi_gas_cold)
milca_test_labels(fermion_ensemble_oracles fast CLASSES quantum
    TARGETS test_fermion_ensemble)
milca_test_labels(result_io_contract fast CLASSES classical quantum
    TARGETS test_result_io classical_pair_demo quantum_reference_demo)
milca_test_labels(animation_export_contract fast CLASSES apps
    TARGETS test_animation_export)
milca_test_labels(screening_q7_oracles fast CLASSES quantum TARGETS test_screening)
milca_test_labels(wigner_q7_oracles fast CLASSES continuum TARGETS test_wigner)
milca_test_labels(mean_field_cell_oracles fast CLASSES quantum
    TARGETS test_mean_field_cell)
milca_test_labels(quantum_collision_qg6_oracles fast CLASSES pairs
    TARGETS test_quantum_collision)
milca_test_labels(hybrid_encounter_qg6_oracles fast CLASSES quantum
    TARGETS test_hybrid_encounter)
milca_test_labels(packet_models_oracles fast CLASSES quantum TARGETS test_packet_models)
milca_test_labels(quantum_classical_h1b_oracles fast CLASSES classical quantum
    TARGETS test_quantum_classical)
milca_test_labels(correlation_hierarchy_h1a_oracles fast CLASSES classical
    TARGETS test_correlation_hierarchy)
milca_test_labels(error_modes_oracles fast CLASSES continuum quantum
    TARGETS test_error_modes)
milca_test_labels(error_verified_matvec_oracles fast CLASSES core
    TARGETS test_error_verified_matvec)
milca_test_labels(error_sampled_norm_oracles fast CLASSES core
    TARGETS test_error_sampled_norm)
milca_test_labels(grad_shafranov_curved_geometry fast CLASSES fem2d
    TARGETS test_curved_triangle)
milca_test_labels(grad_shafranov_curved_convergence fast CLASSES equilibrium
    TARGETS test_gs_curved_convergence)
milca_test_labels(grad_shafranov_p3_curved_convergence fast CLASSES equilibrium
    TARGETS test_gs_curved_convergence)
milca_test_labels(grad_shafranov_cubic_profile_oracles fast CLASSES equilibrium
    TARGETS test_gs_profile_newton)
milca_test_labels(grad_shafranov_profile_quadrature_work fast CLASSES equilibrium
    TARGETS test_gs_profile_newton)
milca_test_labels(grad_shafranov_p2_nonlinear_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
milca_test_labels(grad_shafranov_p3_nonlinear_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
milca_test_labels(grad_shafranov_p3_spline_knot_convergence fast CLASSES equilibrium
    TARGETS test_gs_nonlinear)
milca_test_labels(grad_shafranov_p3_physical_oracles fast CLASSES equilibrium
    TARGETS test_gs_p3)
milca_test_labels(grad_shafranov_p3_thread_invariance fast CLASSES equilibrium
    TARGETS test_gs_threads)
milca_test_labels(grad_shafranov_csr_direct_newton_oracles fast CLASSES equilibrium
    TARGETS test_gs_direct)
foreach(_mode automatic external finalized)
    milca_test_labels(gs_direct_mpi_${_mode} fast CLASSES fem2d
        TARGETS test_gs_direct_mpi)
endforeach()
milca_test_labels(grad_shafranov_shaped_solver_work fast CLASSES equilibrium
    TARGETS test_gs_solver_work)
milca_test_labels(grad_shafranov_axis_refinement fast CLASSES equilibrium
    TARGETS test_gs_axis_refinement)
milca_test_labels(grad_shafranov_profile_law_readback fast CLASSES equilibrium
    TARGETS test_gs_profile_law_readback)
milca_test_labels(grad_shafranov_curved_readback fast CLASSES equilibrium
    TARGETS test_gs_curved_readback)
milca_test_labels(grad_shafranov_p3_curved_readback fast CLASSES equilibrium
    TARGETS test_gs_curved_readback)
foreach(case r0 rp rm zp zm el gzm)
    milca_test_labels(gs_free_3d_nl3a_fixture_return_${case} slow CLASSES mhd
        TARGETS gs_free_3d gs_free)
endforeach()
milca_test_labels(gs_free_3d_fixture_hypotheses slow CLASSES mhd
    TARGETS gs_free_3d gs_free)
milca_test_labels(gs_free_app_contract fast CLASSES equilibrium TARGETS gs_free)
milca_test_labels(grad_shafranov_p2_readback_order fast CLASSES equilibrium
    TARGETS test_gs_readback_order)
milca_test_labels(grad_shafranov_p3_readback_order fast CLASSES equilibrium
    TARGETS test_gs_readback_order)
milca_test_labels(grad_shafranov_p2_readback_free fast CLASSES equilibrium
    TARGETS test_gs_readback_free)
milca_test_labels(grad_shafranov_p3_readback_free fast CLASSES equilibrium
    TARGETS test_gs_readback_free)
milca_test_labels(grad_shafranov_readback_work fast CLASSES equilibrium
    TARGETS test_gs_readback_work)
milca_test_labels(grad_shafranov_stream_readback fast CLASSES equilibrium
    TARGETS gs_phase1)
milca_test_labels(equilibrium_core_d6_exterior_per_n fast CLASSES equilibrium
    TARGETS test_equilibrium_boundary)
milca_test_labels(equilibrium_core_d6_coil_interface fast CLASSES equilibrium
    TARGETS test_equilibrium_boundary)
milca_test_labels(equilibrium_core_d7_geometry_solovev fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
milca_test_labels(equilibrium_core_d7_geometry_circle fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
milca_test_labels(equilibrium_core_d7_geometry_drift fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
milca_test_labels(equilibrium_core_d7_geometry_state fast CLASSES equilibrium
    TARGETS test_equilibrium_geometry)
milca_test_labels(equilibrium_core_chart_jet fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(equilibrium_core_d7_readback_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(equilibrium_core_d3_solovev_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(equilibrium_core_d3_upper_bound_p2 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(gs_response_3d_ampere_indicator_contract fast CLASSES mhd
    TARGETS gs_response_3d)
milca_test_labels(equilibrium_core_d7_readback_p3 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(equilibrium_core_d3_solovev_p3 fast CLASSES equilibrium
    TARGETS test_equilibrium_core)
milca_test_labels(gs_k6_mini_regression slow CLASSES equilibrium
    TARGETS gs_phase1 test_equilibrium_core)
milca_test_labels(coils_3d_biot_savart_fourier_oracles fast CLASSES equilibrium
    TARGETS test_coils_3d)
milca_test_labels(vacuum_surfaces_sfl_harmonics_oracles fast CLASSES mhd
    TARGETS test_vacuum_surfaces)
milca_test_labels(mhd_perturbation_l7_xi_lagrangian_charts fast CLASSES mhd
    TARGETS test_mhd_perturbation)
milca_test_labels(mhd_ampere_n_error_status fast CLASSES mhd
    TARGETS test_mhd_ampere_error)
milca_test_labels(mhd_ampere_n_sheet_coupled_equilibrated fast CLASSES mhd
    TARGETS test_mhd_ampere_sheet_error)
milca_test_labels(mhd_ampere_n_sheet_interpolation_remainder fast CLASSES mhd
    TARGETS test_meridian_edge_remainder)
milca_test_labels(mhd_ampere_n_exterior_manufactured fast CLASSES mhd
    TARGETS test_mhd_ampere_n)
milca_test_labels(mhd_shielding_cylinder_limit_reciprocity fast CLASSES mhd
    TARGETS test_mhd_shielding)
milca_test_labels(mhd_ideal_n_mde_force_free_identity fast CLASSES mhd
    TARGETS test_mhd_ideal_n)
milca_test_labels(mhd_ideal_n_final_state_error fast CLASSES mhd
    TARGETS test_mhd_ideal_n)
milca_test_labels(gs_response_3d_ideal_error_contract fast CLASSES apps mhd
    TARGETS gs_response_3d)
milca_test_labels(mhd_shield_torus_nrad_convergence fast CLASSES mhd
    TARGETS test_mhd_shield_torus)
milca_test_labels(flux_mesh_axis_interpolation_order fast CLASSES mhd
    TARGETS test_flux_mesh_axis)
milca_test_labels(mhd_closure_tearing_delta_prime slow CLASSES mhd
    TARGETS test_mhd_closure)
milca_test_labels(response_layer_all_rational_setup fast CLASSES mhd
    TARGETS test_mhd_closure)
milca_test_labels(gs_response_3d_frobenius_all_rational slow CLASSES apps mhd
    TARGETS gs_response_3d)
milca_test_labels(response_strip_torus_coupling slow CLASSES mhd
    TARGETS test_mhd_closure)
milca_test_labels(response_toroidal_strip_two_ring_traces fast CLASSES mhd
    TARGETS test_response_interfaces)
milca_test_labels(response_cut_total_pressure_balance fast CLASSES mhd
    TARGETS test_response_interfaces)
milca_test_labels(response_toroidal_strip_edge_merged_domains fast CLASSES mhd
    TARGETS test_response_interfaces)
milca_test_labels(response_toroidal_strip_pressure_action fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
milca_test_labels(response_toroidal_strip_pressure_solve fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
milca_test_labels(response_toroidal_strip_radial_convergence fast CLASSES mhd
    TARGETS test_response_toroidal_strip)
milca_test_labels(response_fluid_full_pressure_oracles fast CLASSES mhd
    TARGETS test_response_fluid)
milca_test_labels(response_strip_cylinder_layer fast CLASSES mhd
    TARGETS test_response_strip)
if(TARGET check_response_strip_current)
    milca_test_labels(response_strip_current_exact_identities fast CLASSES mhd
        TARGETS check_response_strip_current)
endif()
if(TARGET check_homma_derivations)
    milca_test_labels(check_homma_derivations fast CLASSES classical
        TARGETS check_homma_derivations)
endif()
milca_test_labels(mhd_resistive_layer_rotation slow CLASSES mhd
    TARGETS test_mhd_closure)
milca_test_labels(mhd_newcomb_screw_pinch_limit fast CLASSES mhd
    TARGETS test_mhd_newcomb)
milca_test_labels(nonnested_n0_solovev_p2 fast CLASSES mhd
    TARGETS test_nonnested_solovev)
milca_test_labels(nonnested_n0_solovev_p3 fast CLASSES mhd
    TARGETS test_nonnested_solovev)
milca_test_labels(nonnested_n0_equilibrated_solovev fast CLASSES mhd
    TARGETS test_nonnested_error)
milca_test_labels(nonnested_transport_equilibrated_n0 fast CLASSES mhd
    TARGETS test_nonnested_transport_error)
milca_test_labels(nn1_spec_case_smoke fast CLASSES mhd TARGETS nn1_spec_case)
milca_test_labels(nonnested_vacuum_helical_coil_island slow CLASSES mhd
    TARGETS test_nonnested_vacuum)
milca_test_labels(nonnested_vacuum_mode_residual_p3 fast CLASSES mhd
    TARGETS test_nonnested_vacuum_mode_residual)
milca_test_labels(nonnested_vacuum_mode_residual_p3_nr32 slow CLASSES mhd
    TARGETS test_nonnested_vacuum_mode_residual)
milca_test_labels(island_diagnostics_pendulum fast CLASSES mhd
    TARGETS test_island_diagnostics)
milca_test_labels(island_driven_outer_reference fast CLASSES mhd
    TARGETS test_island_driven)
milca_test_labels(island_driven_stable_outer_reference fast CLASSES mhd
    TARGETS test_island_driven)
milca_test_labels(island_zero_beta_ohmic fast CLASSES mhd
    TARGETS test_island_driven)
milca_test_labels(island_coil_quadrature fast CLASSES mhd
    TARGETS test_island_vacuum_oracle)
milca_test_labels(island_direct_coil_split fast CLASSES mhd
    TARGETS test_island_direct_coil_split)
milca_test_labels(island_vacuum_seed_p3 fast CLASSES mhd
    TARGETS test_island_vacuum_oracle)
milca_test_labels(island_coarse_feec_map slow CLASSES mhd
    TARGETS test_island_vacuum_oracle)
milca_test_labels(island_coarse_feec_map_p2 slow CLASSES mhd
    TARGETS test_island_vacuum_oracle)
milca_test_labels(island_coarse_nn2_map slow CLASSES mhd
    TARGETS test_nonnested_island)
foreach(degree 2 3)
    milca_test_labels(island_vacuum_width_phase_p${degree} slow CLASSES mhd
        TARGETS test_island_vacuum_oracle)
endforeach()
foreach(k 2 3)
    milca_test_labels(nonnested_hx_vacuum_p${k} fast CLASSES mhd
        TARGETS test_nonnested_hx)
endforeach()
milca_test_labels(nonnested_symbolic_derivatives_and_d3_identity fast CLASSES mhd
    TARGETS test_nonnested_leaves)
milca_test_labels(nonnested_transport_manufactured_3d slow CLASSES mhd
    TARGETS test_nonnested_transport)
milca_test_labels(nonnested_transport_sovinec_chi_ratio_p2 fast CLASSES mhd
    TARGETS test_nonnested_transport)
milca_test_labels(nonnested_transport_sovinec_chi_ratio_p3 fast CLASSES mhd
    TARGETS test_nonnested_transport)
milca_test_labels(nonnested_transport_sovinec_chi_ratio_p4 fast CLASSES mhd
    TARGETS test_nonnested_transport)
milca_test_labels(nonnested_hudson_breslau_integrable_exact_p3 fast CLASSES mhd
    TARGETS test_nonnested_hudson_breslau)
milca_test_labels(nonnested_hudson_breslau_chaos_responds_p3 fast CLASSES mhd
    TARGETS test_nonnested_hudson_breslau)
foreach(_fc ring_p3 helical_p2 helical_p3 exterior_p3 source_p3)
    milca_test_labels(frozen_current_${_fc} fast CLASSES mhd
        TARGETS test_frozen_current)
endforeach()
milca_test_labels(frozen_current_exterior_refined_p3 slow CLASSES mhd
    TARGETS test_frozen_current)
milca_test_labels(nonnested_coupled_solovev_p2 slow CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_coupled_offset_return_p2 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_coupled_fourier_reduction_p2 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_coupled_solovev_p3 slow CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_coupled_offset_return_p3 fast CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_coupled_fourier_reduction_p3 slow CLASSES mhd
    TARGETS test_nonnested_coupled)
milca_test_labels(nonnested_island_flattening_p3 slow CLASSES mhd
    TARGETS test_nonnested_island)
milca_test_labels(nonnested_island_d3_ladder_p2 slow CLASSES mhd
    TARGETS test_nonnested_island)
milca_test_labels(nonnested_coupled_landreman_family1_p2 slow CLASSES mhd
    TARGETS test_nonnested_landreman)
milca_test_labels(nonnested_mapped_identity_regression fast CLASSES mhd
    TARGETS test_nonnested_mapped)
milca_test_labels(nonnested_mapped_ellipse_h_ladder slow CLASSES mhd
    TARGETS test_nonnested_mapped)
milca_test_labels(nonnested_mapped_ellipse_mode_ladder fast CLASSES mhd
    TARGETS test_nonnested_mapped)
milca_test_labels(nonnested_mapped_zero_forcing_residual fast CLASSES mhd
    TARGETS test_nonnested_mapped)
foreach(mapped_case point_oracles cartesian_identity cartesian_quadrature cartesian_vacuum)
    milca_test_labels(nonnested_mapped_${mapped_case} fast CLASSES mhd
        TARGETS test_nonnested_mapped)
endforeach()
milca_test_labels(nonnested_newton_mr_agreement fast CLASSES mhd
    TARGETS test_nonnested_newton_mr)
milca_test_labels(nonnested_newton_mr_lumped_ladder fast CLASSES mhd
    TARGETS test_nonnested_newton_mr)
milca_test_labels(nonnested_qa_oracle_check fast CLASSES mhd TARGETS test_nonnested_qa)
foreach(case hodge solovev beltrami shell sheet)
    milca_test_labels(nonnested_force_${case} fast CLASSES mhd
        TARGETS test_nonnested_force)
endforeach()
# Non-nested Beltrami equilibria (docs/nn1-mapped-3d.md). The SPEC sample tests
# are configured only when MILCA_BELTRAMI_SPEC_SAMPLES is set, so that executable
# is recorded for milca_check_mhd directly instead of through a test.
set_property(GLOBAL APPEND PROPERTY MILCA_CHECK_TARGETS_mhd
    test_nonnested_beltrami_spec)
foreach(case axis helical)
    foreach(rings 2 4 8)
        milca_test_labels(nonnested_beltrami_spec_${case}_${rings} slow CLASSES mhd
            TARGETS test_nonnested_beltrami_spec)
    endforeach()
endforeach()
foreach(modes 4 6)
    milca_test_labels(nonnested_beltrami_spec_modes_${modes} slow CLASSES mhd
        TARGETS test_nonnested_beltrami_spec)
endforeach()
milca_test_labels(nonnested_beltrami_pairing fast CLASSES mhd
    TARGETS test_nonnested_beltrami_pairing)
milca_test_labels(nonnested_beltrami_resonance fast CLASSES mhd
    TARGETS test_nonnested_beltrami)
milca_test_labels(nonnested_beltrami_constraints fast CLASSES mhd
    TARGETS test_nonnested_beltrami)
milca_test_labels(nonnested_beltrami_constraints_mr fast CLASSES mhd
    TARGETS test_nonnested_beltrami)
milca_test_labels(nonnested_beltrami_vacuum_mr fast CLASSES mhd
    TARGETS test_nonnested_beltrami)
foreach(k 2 3)
    foreach(case manufactured lundquist)
        milca_test_labels(nonnested_beltrami_${case}_p${k}_mr fast CLASSES mhd
            TARGETS test_nonnested_beltrami)
    endforeach()
endforeach()
milca_test_labels(nonnested_beltrami_mr_agreement fast CLASSES mhd
    TARGETS test_nonnested_beltrami_mr)
milca_test_labels(nonnested_beltrami_mr_block_ladder fast CLASSES mhd
    TARGETS test_nonnested_beltrami_mr)
foreach(k 2 3)
    foreach(case manufactured lundquist vacuum)
        milca_test_labels(nonnested_beltrami_${case}_p${k} slow CLASSES mhd
            TARGETS test_nonnested_beltrami)
    endforeach()
endforeach()
milca_test_labels(lean_continuous_certificates fast CLASSES continuum
    TARGETS formal_krook_certificate)
milca_test_labels(grad_shafranov_p1_exact_current_antiderivatives fast CLASSES equilibrium
    TARGETS test_gs_p1_current_identities)
milca_test_labels(grad_shafranov_symbolic_identities fast CLASSES equilibrium
    TARGETS test_gs_symbolic)
milca_test_labels(grad_shafranov_periodic_cylinder_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_periodic_cylinder_symbolic)
milca_test_labels(grad_shafranov_coupled_majorant_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_coupled_majorant_identities)
milca_test_labels(grad_shafranov_forward_majorant_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
milca_test_labels(grad_shafranov_forward_majorant_rejects_jump_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
milca_test_labels(grad_shafranov_forward_majorant_rejects_boundary_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
milca_test_labels(grad_shafranov_forward_majorant_rejects_quadrature_omission fast CLASSES equilibrium
    TARGETS test_gs_forward_majorant_identities)
milca_test_labels(grad_shafranov_variational_exact_identities fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)
milca_test_labels(grad_shafranov_rejects_wrong_field_orientation fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)
milca_test_labels(grad_shafranov_rejects_wrong_weak_source_weight fast CLASSES equilibrium
    TARGETS test_gs_variational_identities)

milca_test_labels(grad_shafranov_error_leaf_oracles fast CLASSES equilibrium
    TARGETS test_gs_error_leaves)
milca_test_labels(gs_free_boundary_equilibrated_quick fast CLASSES equilibrium
    TARGETS test_gs_coupled_majorant)
milca_test_labels(gs_bem_single_layer_reference_residual fast CLASSES equilibrium
    TARGETS test_gs_bem_residual)
milca_test_labels(gs_fixed_stability_and_linear_goal fast CLASSES equilibrium
    TARGETS test_gs_error_control)
milca_test_labels(gs_verified_nonlinear_coercivity fast CLASSES equilibrium
    TARGETS test_gs_verified_stability)
milca_test_labels(grad_shafranov_production_estimator fast CLASSES apps
    TARGETS gs_phase1)

# Saved-field Ampere current diagnostics.
foreach(degree 2 3)
    milca_test_labels(grad_shafranov_p${degree}_ampere_loop fast CLASSES equilibrium
        TARGETS test_gs_ampere_loop)
    milca_test_labels(grad_shafranov_p${degree}_weak_ampere_loop fast CLASSES equilibrium
        TARGETS test_gs_weak_ampere_loop)
endforeach()

# Electromagnetic plane-wave teaching case.
foreach(mode oracles vacuum reflection transmission refinement)
    milca_test_labels(em_slab_${mode} fast CLASSES continuum TARGETS test_em_slab)
endforeach()

milca_test_labels(response_strip_equilibrium_jets fast CLASSES mhd
    TARGETS test_response_strip_jets)
milca_test_labels(response_strip_tc24_jets fast CLASSES mhd
    TARGETS test_response_strip_jets)
milca_test_labels(response_strip_tc24_preconditioner slow CLASSES mhd
    TARGETS test_response_strip_tc24_preconditioner)
foreach(check oracles negative limits)
    milca_test_labels(response_strip_preconditioner_${check} fast CLASSES mhd
        TARGETS test_response_strip_preconditioner)
endforeach()
milca_test_labels(response_strip_preconditioner_scaling performance CLASSES mhd
    TARGETS test_response_strip_preconditioner)

milca_test_labels(response_spot_calibration_algebra fast CLASSES mhd
    TARGETS test_response_calibration)
milca_test_labels(response_unit_sheet_spot_cylinder fast CLASSES mhd
    TARGETS test_response_sheet_cylinder)
milca_test_labels(response_qualification_frozen_rules fast CLASSES mhd
    TARGETS test_response_qualification)
milca_test_labels(response_qualification_external_adapter fast CLASSES apps
    TARGETS response_qualify)
milca_test_labels(response_spot_resolved_single_ring fast CLASSES mhd
    TARGETS test_response_strip_readback)
milca_test_labels(response_calibration_deck_validation fast CLASSES apps
    TARGETS gs_response_calibration)
milca_test_labels(response_calibration_deck_rejects_zero_mode fast CLASSES apps
    TARGETS gs_response_calibration)

milca_test_labels(nonnested_transport_anisotropy_pollution fast CLASSES mhd
    TARGETS test_nonnested_anisotropy)
milca_test_labels(nonnested_transport_anisotropy_correction fast CLASSES mhd
    TARGETS test_nonnested_anisotropy)
milca_test_labels(nonnested_transport_fourier_equilibrated fast CLASSES mhd
    TARGETS test_nonnested_fourier_error)

milca_test_labels(gs_native_adapter_contract fast CLASSES fem2d
    TARGETS test_gs_native_adapter)

milca_test_labels(gs_free_boundary_t7_native_lu fast CLASSES equilibrium
    TARGETS test_gs_fb_jacobian)
milca_test_labels(pp_fusion_oracles fast CLASSES quantum TARGETS test_pp_fusion)
milca_test_labels(pp_fusion_export_contract fast CLASSES quantum TARGETS pp_fusion_export)

milca_test_labels(pp_fusion_eft_oracles fast CLASSES quantum TARGETS test_pp_fusion_eft)
milca_test_labels(pp_fusion_eft_export_contract fast CLASSES quantum TARGETS pp_fusion_eft_export)

foreach(mode circle jets tube periodic invalid)
    milca_test_labels(equilibrium_curve_chart_${mode} fast CLASSES equilibrium
        TARGETS test_curve_chart)
endforeach()
if(TARGET gen_curve_chart)
    milca_test_labels(equilibrium_curve_chart_regeneration fast CLASSES equilibrium
        TARGETS gen_curve_chart)
endif()
milca_test_labels(ideal_fb_landreman_axisymmetric fast CLASSES mhd
    TARGETS test_landreman_free_boundary)
milca_test_labels(ideal_fb_landreman_bem3d16_axisymmetric slow CLASSES mhd equilibrium
    TARGETS test_landreman_free_boundary)
milca_test_labels(ideal_fb_landreman_bem3d_axisymmetric slow CLASSES mhd equilibrium
    TARGETS test_landreman_free_boundary)

milca_test_labels(cpu_precision_gs_comparison performance CLASSES mhd
    TARGETS benchmark_cpu_precision)
milca_test_labels(cpu_precision_map_comparison performance CLASSES mhd
    TARGETS benchmark_cpu_precision)
milca_test_labels(cpu_precision_gs_ir performance CLASSES mhd
    TARGETS benchmark_cpu_solve)
# D8/M6: same native solver and manufactured oracle as the T3 ladder.
foreach(start reference rp rm zp zm el)
    milca_test_labels(equilibrium_d8_psi_offset_${start} fast CLASSES equilibrium
        TARGETS test_gs_fb_manufactured)
endforeach()
milca_test_labels(equilibrium_d8_registry_contract fast CLASSES equilibrium
    TARGETS test_equilibrium_harness)
milca_test_labels(equilibrium_nested_exact_oracles fast CLASSES equilibrium
    TARGETS test_nested_oracles)

milca_test_labels(equilibrium_nested_oracle_batch fast CLASSES equilibrium
    TARGETS benchmark_nested_oracles)
milca_test_labels(equilibrium_nested_oracle_regeneration fast CLASSES equilibrium
    TARGETS gen_nested_oracles)
foreach(family 1 2)
    milca_test_labels(equilibrium_nested_boundary_input_${family} fast CLASSES apps
        TARGETS nested_oracle_export test_nested_boundary_input)
endforeach()
foreach(name circular-tokamak circular-tokamak-current)
    milca_test_labels(equilibrium_nested_${name} fast CLASSES apps
        TARGETS nested_equilibrium test_nested_equilibrium_app)
endforeach()
foreach(k 2 3)
    milca_test_labels(equilibrium_core_map_state_p${k} fast CLASSES mhd
        TARGETS test_map_equilibrium)
endforeach()

milca_test_labels(ideal_xi_hessian_and_analytic_limits fast CLASSES mhd
    TARGETS test_ideal_xi)
if(TARGET gen_ideal_xi)
    milca_test_labels(ideal_xi_exact_identity_and_regeneration fast CLASSES mhd
        TARGETS gen_ideal_xi)
endif()
milca_test_labels(ideal_seed_boundary_label fast CLASSES mhd
    TARGETS test_ideal_map_axisymmetric)

milca_test_labels(ideal_seed_3d_boundary_profiles fast CLASSES mhd
    TARGETS test_ideal_seed_3d)
milca_test_labels(ideal_map_field_angle fast CLASSES mhd
    TARGETS test_ideal_map_field_angle)
milca_test_labels(ideal_seed_3d_nested_input fast CLASSES mhd
    TARGETS test_nested_seed_input)
if(TARGET gen_homma)
    milca_test_labels(homma_leaves_regeneration fast CLASSES continuum
        TARGETS gen_homma)
endif()
if(TARGET gen_fp_relaxation)
    milca_test_labels(fp_relaxation_leaves_regeneration fast CLASSES continuum
        TARGETS gen_fp_relaxation)
endif()
if(TARGET gen_ideal_field_angle)
    milca_test_labels(ideal_field_angle_regeneration fast CLASSES mhd
        TARGETS gen_ideal_field_angle)
    milca_test_labels(ideal_map_leaves_regeneration fast CLASSES mhd
        TARGETS gen_ideal_map)
endif()
if(TARGET gen_ideal_seed_3d)
    milca_test_labels(ideal_seed_3d_regeneration fast CLASSES mhd
        TARGETS gen_ideal_seed_3d)
endif()

milca_test_labels(ideal_mhd_coordinate_analytic_limits fast CLASSES mhd
    TARGETS test_ideal_mhd_coordinates)
if(TARGET gen_ideal_mhd_coordinates)
    milca_test_labels(ideal_mhd_coordinates_exact_regeneration fast CLASSES mhd
        TARGETS gen_ideal_mhd_coordinates)
endif()

milca_test_labels(island_closure_declaration fast CLASSES mhd TARGETS test_island_closure)
milca_test_labels(island_closure_report fast CLASSES mhd TARGETS test_island_closure)
milca_test_labels(island_closure_drive fast CLASSES mhd TARGETS test_island_closure)
milca_test_labels(nonnested_coupled_closed_line_p3 slow CLASSES mhd
    TARGETS test_nonnested_landreman)

milca_test_labels(nonnested_g6_external_comparison_contract fast CLASSES mhd
    TARGETS test_g6_external_comparison)

milca_test_labels(nonnested_g6_external_sampled_force fast CLASSES mhd
    TARGETS test_g6_sampled_force)

foreach(k 2 3)
    foreach(case manufactured nested reproduce reproduce3d sampling)
        milca_test_labels(nonnested_fitted_${case}_p${k} fast CLASSES mhd
            TARGETS test_nonnested_profile)
    endforeach()
endforeach()
foreach(case readiness heat current)
    milca_test_labels(island_fitted_${case} fast CLASSES mhd
        TARGETS test_island_fitted)
endforeach()
foreach(p 4 6 8)
    milca_test_labels(ideal_gll_derivatives_p${p} fast CLASSES mhd
        TARGETS benchmark_ideal_gll)
endforeach()
milca_test_labels(ideal_gll_coupled_block_oracle fast CLASSES mhd
    TARGETS test_gll_block_preconditioner)
milca_test_labels(island_provider_tangent_p2 slow CLASSES mhd
    TARGETS test_island_provider)
milca_test_labels(island_provider_nested_p2 fast CLASSES mhd
    TARGETS test_island_provider)
milca_test_labels(island_provider_vacuum_p2 fast CLASSES mhd
    TARGETS test_island_provider)
foreach(case vacuum pressure)
    milca_test_labels(tiago_island_${case} fast CLASSES mhd
        TARGETS test_tiago_island)
endforeach()
milca_test_labels(tiago_island_direct slow CLASSES mhd
    TARGETS test_tiago_island)
foreach(case seed none ladder leaves)
    milca_test_labels(island_seed_${case} fast CLASSES mhd TARGETS test_island_seed)
endforeach()
milca_test_labels(s0_toroidal_continuum_gap fast CLASSES mhd TARGETS test_s0_tae)
milca_test_labels(s0_external_kink_rigid_growth fast CLASSES mhd
    TARGETS test_ideal_free_boundary)
milca_test_labels(equilibrium_map_3d_readback fast CLASSES mhd
    TARGETS test_map_3d_equilibrium)
if(TARGET gen_map_3d_readback)
    milca_test_labels(equilibrium_map_3d_regeneration fast CLASSES mhd
        TARGETS gen_map_3d_readback)
endif()

milca_test_labels(equilibrium_d8_landreman_3d_integer fast CLASSES mhd
    TARGETS test_equilibrium_3d)
milca_test_labels(equilibrium_d8_issan_3d_integer fast CLASSES mhd
    TARGETS test_equilibrium_3d)

milca_test_labels(equilibrium_d8_landreman_3d_p3 fast CLASSES mhd
    TARGETS test_equilibrium_3d)

foreach(gate rotation translation screw manufactured limit derivatives guards)
    milca_test_labels(grad_shafranov_rm0_${gate} fast CLASSES equilibrium
        TARGETS test_gs_symmetry)
    set_property(TEST grad_shafranov_rm0_${gate} APPEND PROPERTY LABELS rm0)
endforeach()
milca_test_labels(grad_shafranov_rm0_symbolic fast CLASSES equilibrium
    TARGETS gen_gs_symmetry)
if(TEST grad_shafranov_rm0_symbolic)
    set_property(TEST grad_shafranov_rm0_symbolic APPEND PROPERTY LABELS rm0)
endif()
foreach(gate derivatives physics energy covariance)
    milca_test_labels(vrmhd_${gate} fast CLASSES mhd TARGETS test_vrmhd_leaves)
endforeach()
milca_test_labels(vrmhd_leaf_regeneration fast CLASSES mhd TARGETS gen_vrmhd_leaves)

foreach(gate geometry derivatives volume closure)
    milca_test_labels(vrmhd_section_${gate} fast CLASSES mhd TARGETS test_vrmhd_section)
endforeach()

milca_test_labels(vrmhd_section_manufactured_consistency fast CLASSES mhd TARGETS test_vrmhd_mms)
milca_test_labels(vrmhd_section_spectra slow CLASSES mhd TARGETS test_vrmhd_spectra)
milca_test_labels(vrmhd_section_spectra_fine slow CLASSES mhd TARGETS test_vrmhd_spectra)
milca_test_labels(vrmhd_spectral_folding fast CLASSES mhd TARGETS test_vrmhd_spectral_folding)
foreach(mode spurious ambiguous)
    milca_test_labels(vrmhd_resolved_power_rejects_${mode} fast CLASSES mhd
        TARGETS test_vrmhd_spectral_folding)
endforeach()
milca_test_labels(vrmhd_fourier_content fast CLASSES mhd TARGETS test_vrmhd_fourier)
milca_test_labels(vrmhd_section_energy fast CLASSES mhd TARGETS test_vrmhd_energy)
milca_test_labels(vrmhd_section_mixed_derivatives fast CLASSES mhd TARGETS test_vrmhd_section)
foreach(degree 1 2)
    milca_test_labels(vrmhd_section_mixed_energy_p${degree} fast CLASSES mhd
        TARGETS test_vrmhd_energy)
    milca_test_labels(vrmhd_section_mixed_mms_p${degree} fast CLASSES mhd
        TARGETS test_vrmhd_mms)
endforeach()
foreach(gate boundary derivatives harmonic solve)
    milca_test_labels(vrmhd_stationary_${gate} fast CLASSES mhd
        TARGETS test_vrmhd_stationary)
endforeach()
milca_test_labels(vrmhd_stationary_gauge slow CLASSES mhd
    TARGETS test_vrmhd_stationary)
milca_test_labels(vrmhd_gs_consistency fast CLASSES mhd TARGETS test_vrmhd_gs)
milca_test_labels(vrmhd_gs_regeneration fast CLASSES mhd TARGETS gen_vrmhd_gs)

milca_test_labels(vrmhd_section_regeneration fast CLASSES mhd TARGETS gen_vrmhd_section)

milca_test_labels(vrmhd_mms_regeneration fast CLASSES mhd TARGETS gen_vrmhd_mms)

milca_test_labels(vrmhd_pairing_measurement fast CLASSES mhd TARGETS vrmhd_pairing_study)
milca_test_labels(vrmhd_pairing_benchmark performance CLASSES mhd TARGETS vrmhd_pairing_study)
foreach(gate energy decay waves steady)
    milca_test_labels(vrmhd_evolve_${gate} fast CLASSES mhd
        TARGETS test_vrmhd_evolve)
endforeach()

foreach(gate exactness reduction decoupling manufactured accumulation)
    milca_test_labels(vrmhd_harmonics_${gate} fast CLASSES mhd TARGETS test_vrmhd_harmonics)
endforeach()
milca_test_labels(vrmhd_harmonics_spectra slow CLASSES mhd TARGETS test_vrmhd_harmonics)
foreach(gate exactness reduction decoupling)
    milca_test_labels(vrmhd_harmonics_mixed_${gate} fast CLASSES mhd
        TARGETS test_vrmhd_harmonics)
endforeach()
milca_test_labels(vrmhd_harmonics_mixed_manufactured slow CLASSES mhd
    TARGETS test_vrmhd_harmonics)
milca_test_labels(vrmhd_harmonics_spectra_fine slow CLASSES mhd
    TARGETS test_vrmhd_harmonics)
milca_test_labels(vrmhd_mms_fourier_regeneration fast CLASSES mhd
    TARGETS gen_vrmhd_mms_fourier)

foreach(gate coupling energy)
    milca_test_labels(vrmhd_fourier_nonlinear_${gate} fast CLASSES mhd
        TARGETS test_vrmhd_fourier_nonlinear)
endforeach()

milca_test_labels(vrmhd_fourier_nonlinear_manufactured_equal fast
    CLASSES mhd TARGETS test_vrmhd_harmonics)
milca_test_labels(vrmhd_fourier_nonlinear_manufactured_mixed slow
    CLASSES mhd TARGETS test_vrmhd_harmonics)
milca_test_labels(vrmhd_mms_nonlinear_regeneration fast CLASSES mhd
    TARGETS gen_vrmhd_mms_fourier)

milca_test_labels(vrmhd_harmonics_helical fast CLASSES mhd
    TARGETS test_vrmhd_harmonics)

milca_test_labels(vrmhd_helical_nonlinear_manufactured_equal fast
    CLASSES mhd TARGETS test_vrmhd_harmonics)
milca_test_labels(vrmhd_helical_nonlinear_manufactured_mixed slow
    CLASSES mhd TARGETS test_vrmhd_harmonics)

milca_test_labels(vrmhd_fourier_nonlinear_derivatives slow CLASSES mhd
    TARGETS test_vrmhd_fourier_nonlinear)

foreach(family equal mixed)
    foreach(gate boundary equilibrium width ledger graded)
        milca_test_labels(vrmhd_island_${gate}_${family} fast CLASSES mhd
            TARGETS test_vrmhd_island)
    endforeach()
endforeach()
milca_test_labels(vrmhd_island_induction_mixed fast CLASSES mhd
    TARGETS test_vrmhd_island)
milca_test_labels(vrmhd_island_regeneration fast CLASSES mhd
    TARGETS gen_vrmhd_island)

milca_test_labels(vrmhd_island_preconditioner fast CLASSES mhd
    TARGETS test_vrmhd_island)
milca_test_labels(vrmhd_island_decimal_steps fast CLASSES mhd
    TARGETS test_vrmhd_island)
milca_test_labels(vrmhd_island_graded_input fast CLASSES mhd
    TARGETS test_vrmhd_island)
