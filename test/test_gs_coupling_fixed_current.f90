program test_gs_coupling_fixed_current
    !! T2 (DESIGN 5.1): a prescribed toroidal current density j(R, Z),
    !! independent of psi and supported inside the fitted Gamma_p chain.
    !! The coupled solution of a(psi, v) + <S psi, v> = int mu0 j v (+ coils)
    !! must equal psi_p + sum_k I_k g_k, with psi_p = mu0 int G j from an
    !! independent tensor cubature (polar about the evaluation point inside
    !! the current, elliptic coordinates outside). Compared at points inside
    !! Gamma_c, including points inside the current, and on Gamma_c (all trace
    !! nodes, edge midpoints included). Also checked: Ampere's law, the total
    !! energy against mu0 int j psi_p dA, and the energy balance of the
    !! solution.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_mu0
    use kin6d_grad_shafranov, only: gs_state_t
    use kin6d_gs_coils, only: gs_circuit_set_t, gs_circuit_flux_vectors
    use kin6d_gs_coupling, only: gs_coupling_t, gs_coupling_build, gs_coupled_apply, &
                     gs_coupling_coil_load, gs_coupling_source_load, gs_coupled_solve, &
                                 gs_coupling_energies, GS_COUPLING_OK
    use gs_coupling_fixture, only: pi, gp_a, gp_b, plasma_current, j0, build_state, &
                      build_coils, density, plasma_flux_ref, plasma_energy_ref, fe_value
    implicit none

    real(dp), parameter :: currents(2) = [1.2e5_dp, -0.8e5_dp]
    real(dp), parameter :: samples(2, 7) = reshape([3.0_dp, 0.0_dp, 3.4_dp, 0.5_dp, &
                     3.5_dp, 0.2_dp, 2.5_dp, -0.9_dp, 3.9_dp, -0.6_dp, 2.2_dp, 1.0_dp, &
                                                    3.8_dp, 1.3_dp], [2, 7])
    real(dp) :: energy_ref
    integer :: failures

    failures = 0
    call check_oracle()
    energy_ref = plasma_energy_ref(2)
    call sweep(2, [24, 48, 96, 192], 1, .false., 3.0_dp)
    call sweep(3, [24, 48, 96], 4, .false., 3.5_dp)
    call sweep(3, [48], 4, .true., 0.0_dp)
    call panel_refinement()
    if (failures /= 0) error stop 'gs_coupling fixed-current oracle failures'
    print '(a)', 'gs_coupling fixed current (T2) oracles passed'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            failures = failures + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine require

    subroutine check_oracle()
        !! The cubature oracle itself: two resolutions agree, and the current
        !! it integrates is the stated one (Ampere, via the flux far away is
        !! covered by the energy test).
        real(dp) :: coarse, fine, total
        integer :: k

        do k = 1, size(samples, 2)
            coarse = plasma_flux_ref(samples(1, k), samples(2, k), 1)
            fine = plasma_flux_ref(samples(1, k), samples(2, k), 2)
            call require(abs(coarse - fine) < 1.0e-10_dp * abs(fine), &
                         'cubature oracle converged at sample')
        end do
        coarse = plasma_energy_ref(1)
        call require(abs(coarse - plasma_energy_ref(2)) < 1.0e-8_dp * abs(coarse), &
                     'energy oracle converged')
        ! int j dA = j0 pi a b / 4 = Ip.
        total = j0 * pi * gp_a * gp_b / 4.0_dp
        call require(abs(total - plasma_current) < 1.0e-9_dp * plasma_current, &
                     'current density normalisation')
    end subroutine check_oracle

    subroutine solve_case(n_outer, order, panels, with_coils, errors, iterations)
        !! errors = [rel max trace-node error, rel max sample error, rel Ampere
        !! error, rel energy error, rel energy balance, rel residual].
        integer, intent(in) :: n_outer, order, panels
        logical, intent(in) :: with_coils
        real(dp), intent(out) :: errors(6)
        integer, intent(out) :: iterations
        type(gs_state_t) :: state
        type(gs_circuit_set_t) :: set
        type(gs_coupling_t) :: coupling
        real(dp), allocatable :: load(:, :), source(:), rhs(:), psi(:), action(:)
        real(dp), allocatable :: ref_trace(:), coil_flux(:, :)
        real(dp) :: scale, value, ref, interior, exterior, ampere
        integer :: status, i, k
        logical :: found

        call build_state(state, n_outer, order, status)
        call require(status == 0, 'two-loop mesh')
        call gs_coupling_build(state, coupling, status, panels_per_edge=panels)
        call require(status == GS_COUPLING_OK, 'coupling build')
        call gs_coupling_source_load(state, density, source, status)
        call require(status == GS_COUPLING_OK, 'source load')
        rhs = source
        if (with_coils) then
            call build_coils(set)
            call gs_coupling_coil_load(state, coupling, set, load, status)
            call require(status == GS_COUPLING_OK, 'coil load')
            rhs = rhs + matmul(load, currents)
        end if
        call gs_coupled_solve(coupling, rhs, psi, status, 1.0e-11_dp, iterations)
        call require(status == GS_COUPLING_OK, 'coupled CG converged')
        allocate (action(size(psi)), ref_trace(coupling % n_trace))
        call gs_coupled_apply(coupling, psi, action)
        errors(6) = norm2(rhs - action) / norm2(rhs)

        do i = 1, coupling % n_trace
            k = coupling % trace_node(i)
           ref_trace(i) = plasma_flux_ref(state % points(1, k), state % points(2, k), 1)
        end do
        if (with_coils) then
            allocate (coil_flux(coupling % n_trace, 2))
           call gs_circuit_flux_vectors(set, state % points(1, coupling % trace_node), &
                                    state % points(2, coupling % trace_node), coil_flux)
            ref_trace = ref_trace + matmul(coil_flux, currents)
        end if
        scale = maxval(abs(ref_trace))
        errors(1) = maxval(abs(psi(coupling % trace_node) - ref_trace)) / scale
        errors(2) = 0.0_dp
        do i = 1, size(samples, 2)
            call fe_value(state, psi, samples(1, i), samples(2, i), value, found)
            call require(found, 'sample point located')
            ref = plasma_flux_ref(samples(1, i), samples(2, i), 1)
            if (with_coils) then
                call gs_circuit_flux_vectors(set, samples(1:1, i), samples(2:2, i), &
                                             coil_flux(1:1, :))
                ref = ref + dot_product(coil_flux(1, :), currents)
            end if
            errors(2) = max(errors(2), abs(value - ref) / scale)
        end do

        ! Ampere: the exterior flux through Gamma_c is mu0 Ip (the coil flux
        ! has none: its source-free interior field carries no net flux).
        call gs_coupling_energies(coupling, psi, interior, exterior)
        if (.not. with_coils) then
           ampere = sum(matmul(coupling % s_trace, psi(coupling % trace_node))) / gs_mu0
            errors(3) = abs(ampere - plasma_current) / plasma_current
            ! Whole-plane energy of the plasma flux, mu0 int j psi dA.
            errors(4) = abs(interior + exterior - energy_ref) / energy_ref
        else
            errors(3) = 0.0_dp
            errors(4) = 0.0_dp
        end if
        ! Energy balance of the discrete solution: psi^T A psi = rhs^T psi.
        errors(5) = abs(interior + exterior - dot_product(rhs, psi)) &
                    / abs(dot_product(rhs, psi))
    end subroutine solve_case

    subroutine sweep(order, levels, panels, with_coils, min_rate)
        integer, intent(in) :: order, levels(:), panels
        logical, intent(in) :: with_coils
        real(dp), intent(in) :: min_rate
        real(dp) :: errors(6, size(levels)), rate(max(1, size(levels) - 1))
        integer :: i, n, iterations

        n = size(levels)
        do i = 1, n
         call solve_case(levels(i), order, panels, with_coils, errors(:, i), iterations)
           print '(a,i1,a,l1,a,i4,a,6es10.2,a,i5)', 'P', order, ' coils ', with_coils, &
                ' edges ', levels(i), '  trace pt ampere energy balance res ', &
                errors(:, i), '  cg ', iterations
            call require(errors(6, i) < 1.0e-9_dp, 'coupled residual')
          call require(errors(5, i) < 1.0e-8_dp, 'energy balance psi^T A psi = f^T psi')
            if (.not. with_coils) call require(errors(3, i) < 1.0e-7_dp, &
                                               'Ampere law: exterior flux = mu0 Ip')
        end do
        if (n < 2) then
            call require(errors(1, 1) < 1.0e-5_dp .and. errors(2, 1) < 1.0e-5_dp, &
                         'plasma plus coils against cubature plus direct coil flux')
            return
        end if
        rate = log(errors(1, :n - 1) / errors(1, 2:)) / log(2.0_dp)
        print '(a,i1,a,10f7.2)', 'P', order, ' observed trace rates ', rate
        rate = log(errors(2, :n - 1) / errors(2, 2:)) / log(2.0_dp)
        print '(a,i1,a,10f7.2)', 'P', order, ' observed sample rates ', rate
        call require(all(errors(1, 2:) < errors(1, :n - 1)), 'trace error decreases')
        call require(all(errors(2, 2:) < errors(2, :n - 1)), 'sample error decreases')
        ! Mean trace rate over the sweep and sample rate on the finest pair.
  call require(log(errors(1, 1)/errors(1, n))/log(2.0_dp)/(n - 1) > min_rate - 1.0_dp, &
                     'mean trace rate')
        ! The sample error is the maximum over a few points of a signed
        ! mesh-realisation dependent error: with the CG tolerance varied from 1e-9
        ! to 1e-12 the sample values are unchanged, but a flag that perturbs the
        ! triangulation (the axis test shows 2616 vs 2622 cells) moved the finest-mesh
        ! error by a factor 1.8 (single-pair rate 2.9 to 3.7). The rate is
        ! therefore the slope over the two finest pairs.
        call require(log(errors(2, n - 2)/errors(2, n))/log(2.0_dp)/2 > min_rate, &
                     'sample rate over the two finest pairs')
      call require(errors(4, n) < 1.0e-6_dp, 'whole-plane energy against mu0 int j psi')
    end subroutine sweep

    subroutine panel_refinement()
        !! Fixed P3 mesh: the BEM trace is the P2 interpolant of the P3 trace on
        !! `panels` panels per edge, and its error falls like panels**3 until
        !! the finite element error is reached.
        integer, parameter :: panels(3) = [1, 2, 4]
        real(dp) :: errors(6, 3)
        integer :: i, iterations

        do i = 1, 3
            call solve_case(96, 3, panels(i), .false., errors(:, i), iterations)
            print '(a,i2,a,6es10.2)', 'P3 edges 96 panels/edge ', panels(i), &
                '  trace pt ampere energy balance res ', errors(:, i)
        end do
        call require(errors(1, 2) < 0.1_dp*errors(1, 1), 'trace error falls 10x, 1 to 2 panels')
        call require(errors(1, 3) < 0.25_dp*errors(1, 2), 'trace error falls 4x, 2 to 4 panels')
    end subroutine panel_refinement

end program test_gs_coupling_fixed_current
