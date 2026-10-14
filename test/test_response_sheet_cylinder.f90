program test_response_sheet_cylinder
    !! Independent periodic-cylinder scalar-potential oracle. A unit outward
    !! n x grad Lambda sheet obeys [potential]=1 and continuous B_normal.
    !! I_m from its convergent series, K_m from its integral representation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
                                       meridian_point_basis_t, meridian_point_basis_init
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
                                  ampere_n_curl, mu0
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
                                    ring_flux_harmonic, ring_sheet_load, ring_area
    use kin6d_response_calibration, only: spot_weights, spot_jump, &
                                          equivalent_sheet_field
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), rs = 0.5_dp
    integer, parameter :: m = 3, n = 2, nf = 256
    real(dp) :: r0, nodes(5, 2), w(5, 2), th(nf), radius, k, z, ip, kp
    real(dp) :: self_ref, jump_ref, area, errs(2, 2), potential, at
    real(dp), allocatable :: labels(:), extra(:)
    complex(dp), allocatable :: f(:), a(:), field(:)
    complex(dp) :: flux(5, 2), self, jump, b
    type(flux_mesh_t) :: fm
    type(meridian_feec_t) :: fs
    type(meridian_point_basis_t) :: pb
    type(ampere_n_t) :: amp
    type(ring_angles_t) :: ra
    integer :: lev, side, i, j, js, st
    th = [(2*pi*(i - 1)/nf, i=1, nf)]
    call ring_angles_set(real(m, dp)/n, th, real(m, dp)/n*th, ra)
    do side = 1, 2
        do i = 1, 5
            nodes(i, side) = rs**2 + merge(-0.01_dp, 0.01_dp, side == 1) &
                             + (i - 3)*0.002_dp
        end do
    end do
    extra = [rs, sqrt(nodes(:, 1)), sqrt(nodes(:, 2))]
    call ring_labels_merge(32, extra, 0.0_dp, labels)
    call spot_weights(nodes, rs**2, n, 1.0_dp, 0.02_dp, w, st)
    if (st /= 0) error stop 'spot weights'
    do lev = 1, 2
        r0 = 5.0_dp*lev
        call flux_mesh_build(labels, 8, 512, 2, circle, fm, st)
        if (st /= 0) error stop 'mesh'
        call meridian_feec_build(fm%mm, 2, fs, st)
        if (st /= 0) error stop 'FEEC'
        call ampere_n_build(fs, n, amp, st)
        if (st /= 0) error stop 'Ampere'
        call meridian_point_basis_init(fs, pb, st)
        if (st /= 0) error stop 'point basis'
        allocate (f(amp%nn), a(amp%nn))
        f = (0.0_dp, 0.0_dp)
        js = minloc(abs(labels - rs), 1)
        call ring_sheet_load(fs, pb, fm, js, ra, n, m, (1.0_dp, 0.0_dp), f)
        call ampere_n_solve(amp, f, a, st)
        if (st /= 0) error stop 'unit sheet solve'
        call ampere_n_curl(amp, a, field)
        self = ring_flux_harmonic(fs, pb, fm, js, ra, n, field, m)
        area = ring_area(fm, js)
        do side = 1, 2
            do i = 1, 5
                j = minloc(abs(labels**2 - nodes(i, side)), 1)
                flux(i, side) = ring_flux_harmonic(fs, pb, fm, j, ra, n, field, m)
            end do
        end do
        jump = spot_jump(w, flux)
        k = real(n, dp)/r0; z = k*rs
        ip = 0.5_dp*(bessi(m - 1, z) + bessi(m + 1, z))
        kp = -0.5_dp*(bessk(m - 1, z) + bessk(m + 1, z))
        self_ref = 2*mu0*r0*z*z*ip*kp
        jump_ref = 0.0_dp
        do side = 1, 2
            at = rs**2 + merge(-0.01_dp, 0.01_dp, side == 1)
            radius = sqrt(at)
            if (side == 1) then
                potential = z*kp*bessi(m, k*radius)
            else
                potential = z*ip*bessk(m, k*radius)
            end if
            jump_ref = jump_ref + merge(-1.0_dp, 1.0_dp, side == 1)* &
                       mu0*r0*(k*k + m*m/at)*potential
        end do
        errs(1, lev) = abs(self/self_ref - 1)
        errs(2, lev) = abs(jump/jump_ref - 1)
        call equivalent_sheet_field(cmplx(0.7_dp, -0.4_dp, dp)*jump, jump, &
                                    self, area, b, st)
        print '(a,f6.1,a,2es12.3,a,2es14.6)', 'R0=', r0, &
            ' self/jump relative errors=', errs(:, lev), ' singflx=', b
        if (st /= 0) error stop 'normalization'
        if (abs(b + cmplx(0.7_dp, -0.4_dp, dp)*self_ref/(r0*rs)) &
            > 1.0e-3_dp*abs(self_ref/(r0*rs))) error stop 'sign/peak/area'
        deallocate (f, a, field)
    end do
    if (maxval(errs(:, 2)) > 1.0e-3_dp) error stop 'cylinder absolute oracle'
    if (errs(1, 2) > 0.4_dp*errs(1, 1)) error stop 'self flux torus limit convergence'
    ! Differentiating sampled FEEC fluxes has a radial-error floor; do not
    ! assert an aspect-ratio order after its curvature correction is hidden.
    if (maxval(errs(2, :)) > 1.0e-5_dp) error stop 'spot jump absolute oracle'
    print '(a)', 'unit-sheet spot cylinder calibration: PASS'
contains
    subroutine circle(u, theta, point)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: point(2)
        point = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine circle
    real(dp) function bessi(order, x) result(value)
        integer, intent(in) :: order
        real(dp), intent(in) :: x
        real(dp) :: term
        integer :: j
        term = (0.5_dp*x)**order
        do j = 2, order
            term = term/j
        end do
        value = term
        do j = 1, 100
            term = term*(0.5_dp*x)**2/(j*(j + order))
            value = value + term
            if (abs(term) < epsilon(value)*abs(value)) exit
        end do
    end function bessi
    real(dp) function bessk(order, x) result(value)
        integer, intent(in) :: order
        real(dp), intent(in) :: x
        real(dp), parameter :: h = 0.001_dp
        real(dp) :: t
        integer :: j
        value = 0.5_dp*exp(-x)
        do j = 1, ceiling(acosh(100.0_dp/x)/h)
            t = j*h
            value = value + exp(-x*cosh(t))*cosh(order*t)
        end do
        value = value*h
    end function bessk
end program test_response_sheet_cylinder
