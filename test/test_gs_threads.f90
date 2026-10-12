program test_gs_threads
    ! Thread invariance of the element-parallel P3 load, the cached majorant
    ! recovery and the contour readback. The load gathers per-element results in
    ! element order, so it must be bitwise independent of the thread count; the
    ! other results are sums that agree to a few ulp.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_set_num_threads
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, solve_gs
    use kin6d_gs_profile_law, only: gs_profile_law_t
    use kin6d_gs_nonlinear, only: gs_profile_jacobian_t, gs_cubic_load
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    use kin6d_gs_majorant, only: gs_functional_majorant
    implicit none
    integer, parameter :: counts(3) = [1, 3, 8]
    type(gs_state_t) :: state, uncached
    type(gs_profile_law_t) :: law
    type(gs_profile_jacobian_t) :: reference, trial
    type(gs_readback_t) :: base, other
    real(dp), allocatable :: rhs_ref(:), rhs(:)
    real(dp) :: contour(2, 32), integral_ref, integral, flux, balance
    real(dp) :: majorant_ref, majorant, majorant_uncached
    integer :: i, status, n

    do i = 1, 32
        call boundary(real(i - 1, dp)/32, contour(:, i))
    end do
    call initialize_gs(state, contour, 0.06_dp, status, 3, curve=boundary)
    if (status /= 0) error stop 'curved P3 mesh'
    call solve_gs(state, 1.0_dp, -0.3_dp, -0.2_dp, status)
    if (status /= 0) error stop 'Solovev reference'
    n = size(state%psi)
    allocate (rhs_ref(n), rhs(n))

    ! Normalised cubic law: many cut cells and the adaptive quadrature.
    allocate (law%knots(2), law%pp(4, 1), law%ff(4, 1))
    law%knots = [0.0_dp, 1.0_dp]
    law%pp(:, 1) = [-0.3_dp, -0.08_dp, -0.02_dp, -0.01_dp]
    law%ff = 0
    law%f_edge = 5; law%p_edge = 0; law%normalized = .true.
    law%axis = minval(state%psi)
    call omp_set_num_threads(1)
    call gs_cubic_load(state, 1.0_dp, law, rhs_ref, reference, status, 1.0e-10_dp, &
                       source_integral=integral_ref)
    if (status /= 0) error stop 'serial cubic load'
    do i = 2, size(counts)
        call omp_set_num_threads(counts(i))
        call gs_cubic_load(state, 1.0_dp, law, rhs, trial, status, 1.0e-10_dp, &
                           source_integral=integral)
        if (status /= 0) error stop 'parallel cubic load'
        if (any(rhs /= rhs_ref)) error stop 'load depends on thread count'
        if (any(trial%mass /= reference%mass)) error stop 'mass depends on threads'
        if (any(trial%axis_load /= reference%axis_load)) &
            error stop 'axis load depends on thread count'
        if (integral /= integral_ref) error stop 'source integral depends on threads'
        if (any(trial%mval /= reference%mval)) error stop 'CSR mass values'
    end do

    ! Cached majorant recovery: the first call builds the mesh-only system, later
    ! calls reuse it. All agree with the uncached evaluation to the recovery
    ! tolerance, at any thread count.
    law%normalized = .false.; law%knots = [-10.0_dp, 0.0_dp]
    law%pp = 0; law%ff(:, 1) = [0.3_dp, 0.05_dp, 0.0_dp, 0.0_dp]
    uncached = state
    nullify (uncached%recovery)
    call omp_set_num_threads(1)
    call gs_functional_majorant(uncached, 1.0_dp, 0.0_dp, 0.0_dp, majorant_uncached, &
                                flux, balance, status, law=law)
    if (status /= 4) error stop 'uncached majorant'
    call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, majorant_ref, &
                                flux, balance, status, law=law)
    if (status /= 4) error stop 'majorant'
    if (.not. state%recovery%ready) error stop 'recovery system not cached'
    if (abs(majorant_ref - majorant_uncached) > 1.0e-9_dp*majorant_ref) &
        error stop 'cached majorant differs from uncached'
    do i = 2, size(counts)
        call omp_set_num_threads(counts(i))
        call gs_functional_majorant(state, 1.0_dp, 0.0_dp, 0.0_dp, majorant, &
                                    flux, balance, status, law=law)
        if (status /= 4) error stop 'threaded majorant'
        if (abs(majorant - majorant_ref) > 1.0e-9_dp*majorant_ref) &
            error stop 'majorant depends on thread count'
    end do

    ! Contour readback over levels, with the reduced volume and flux integrals.
    law%pp = 0; law%ff = 0
    law%pp(1, :) = -0.3_dp; law%ff(1, :) = -0.2_dp
    call omp_set_num_threads(1)
    call read_gs_equilibrium(state, 5.0_dp, -0.2_dp, 65, 160, base, status, law)
    if (status /= 0) error stop 'serial readback'
    do i = 2, size(counts)
        call omp_set_num_threads(counts(i))
        call read_gs_equilibrium(state, 5.0_dp, -0.2_dp, 65, 160, other, status, law)
        if (status /= 0) error stop 'parallel readback'
        if (maxval(abs(other%q - base%q)) > 1.0e-13_dp*maxval(abs(base%q))) &
            error stop 'q depends on thread count'
        if (maxval(abs(other%s_tor - base%s_tor)) > 1.0e-13_dp) &
            error stop 's_tor depends on thread count'
        if (abs(other%volume - base%volume) > 1.0e-13_dp*base%volume) &
            error stop 'volume depends on thread count'
        if (abs(other%phi_edge - base%phi_edge) > 1.0e-13_dp*abs(base%phi_edge)) &
            error stop 'toroidal flux depends on thread count'
    end do
    print *, 'PASS: P3 load, cached majorant and readback are thread invariant'
contains
    subroutine boundary(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [3 + cos(2*acos(-1.0_dp)*t), 1.4_dp*sin(2*acos(-1.0_dp)*t)]
    end subroutine
end program
