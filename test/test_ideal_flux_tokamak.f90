program test_ideal_flux_tokamak
    ! P1e.3 toroidal checks on an exact Solov'ev equilibrium (mu0 = rho = 1):
    ! psi = s0 [(R^2 - R0^2)^2/8 + R^2 Z^2/(2 kappa^2)], p' = -s0 (1 + 1/kappa^2),
    ! F = F0. The PEST and Hamada charts must agree with the grid field, satisfy
    ! force balance, and give the same spectrum.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsparse, only: csc_z_t, fortsparse_status_t, FORTSPARSE_OK
    use milca_vacuum_surfaces, only: psi_grid_t, psi_grid_init, angle_pest, &
        angle_hamada
    use milca_ideal_flux, only: flux_space_t, flux_space_init, flux_assemble, &
        flux_point_t
    use milca_ideal_flux_leaves, only: flux_jets_array
    use milca_s0_spectrum, only: s0_pencil_t, s0_report_t, s0_prepare, s0_solve
    use fortnum_eigen, only: eigen_options_t, EIGEN_OK, EIGEN_NEAREST
    use milca_ideal_flux_tokamak, only: tokamak_chart_t, tokamak_chart_init, &
        tokamak_position
    implicit none
    real(dp), parameter :: r0 = 3, s0 = 1, kappa = 1.5_dp, psib = 3, &
        psin_wall = 0.8_dp
    real(dp) :: f0
    type(psi_grid_t) :: grid
    type(tokamak_chart_t) :: pest, hamada, geometric
    character(32) :: mode
    integer :: ok
    mode = 'all'
    if (command_argument_count() > 0) call get_command_argument(1, mode)
    f0 = 16
    call solovev_grid(grid)
    call tokamak_chart_init(pest, grid, psin_wall, angle_pest, 24, 24, 96, 512, &
        1.0_dp, ok)
    if (ok /= 0) error stop 'PEST chart'
    call tokamak_chart_init(hamada, grid, psin_wall, angle_hamada, 24, 24, 96, &
        512, 1.0_dp, ok)
    if (ok /= 0) error stop 'Hamada chart'
    call tokamak_chart_init(geometric, grid, psin_wall, 0, 24, 24, 96, 512, &
        1.0_dp, ok)
    if (ok /= 0) error stop 'geometric chart'
    if (mode == 'all' .or. mode == 'field') then
        call jets_check(pest, 'PEST  ')
        call jets_check(hamada, 'Hamada')
        call jets_check(geometric, 'geom. ')
        call field_and_balance(pest, 'PEST  ')
        call field_and_balance(hamada, 'Hamada')
        call field_and_balance(geometric, 'geom. ')
    end if
    if (mode == 'all' .or. mode == 'invariance') call invariance()
    print '(a)', 'PASS ideal flux tokamak ('//trim(mode)//')'
contains
    subroutine solovev_grid(g)
        type(psi_grid_t), intent(out) :: g
        integer, parameter :: nr = 641, nz = 801, nf = 101
        real(dp) :: psirz(nr, nz), fpol(nf), pres(nf), r, z, rleft, rdim, zdim
        integer :: i, j, status
        rleft = 1.2_dp; rdim = 3.6_dp; zdim = 5.0_dp
        do j = 1, nz
            z = -zdim/2 + zdim*(j - 1)/(nz - 1)
            do i = 1, nr
                r = rleft + rdim*(i - 1)/(nr - 1)
                psirz(i, j) = s0*((r*r - r0*r0)**2/8 + r*r*z*z/(2*kappa**2))
            end do
        end do
        fpol = f0
        do i = 1, nf
            pres(i) = s0*(1 + 1/kappa**2)*psib*(1 - real(i - 1, dp)/(nf - 1))
        end do
        call psi_grid_init(g, psirz, rleft, rdim, 0.0_dp, zdim, r0, 0.0_dp, psib, &
            1.0_dp, fpol, status, pres)
        if (status /= 0) error stop 'Solovev grid'
    end subroutine solovev_grid

    subroutine field_and_balance(chart, label)
        type(tokamak_chart_t), intent(in) :: chart
        character(*), intent(in) :: label
        type(flux_point_t) :: pt
        real(dp) :: jac, b(3), g(3, 3), q(3), t(12, 12), bref(3), cur(3), force(3)
        real(dp) :: errb, errf, rho, theta
        real(dp), parameter :: radii(5) = [0.03_dp, 0.1_dp, 0.3_dp, 0.6_dp, 0.95_dp]
        integer :: i, j
        errb = 0; errf = 0
        do j = 1, 7
            do i = 1, 5
                rho = radii(i); theta = 0.9_dp*j
                call chart%evaluate(rho, theta, 0.0_dp, pt)
                call flux_jets_array(pt%xu, pt%xuu, pt%psi1, pt%psi2, pt%iota, &
                    pt%iota1, pt%pp, pt%lam, jac, b, g, q, t)
                bref = solovev_field(tokamak_position(chart, rho, theta, 0.0_dp))
                errb = max(errb, norm2(b - bref)/norm2(bref))
                cur = [g(3, 2) - g(2, 3), g(1, 3) - g(3, 1), g(2, 1) - g(1, 2)]
                force = [cur(2)*b(3) - cur(3)*b(2), cur(3)*b(1) - cur(1)*b(3), &
                    cur(1)*b(2) - cur(2)*b(1)] - q
                errf = max(errf, norm2(force)/norm2(q))
            end do
        end do
        print '(a,a,2es10.2)', label, ' field error, force residual', errb, errf
        if (errb > 1e-6_dp) error stop 'chart field differs from the exact field'
        if (errf > 3e-3_dp) error stop 'force balance'
    end subroutine field_and_balance

    subroutine jets_check(chart, label)
        ! Chart jets against central differences of the chart position.
        type(tokamak_chart_t), intent(in) :: chart
        character(*), intent(in) :: label
        type(flux_point_t) :: pt, pp_, pm
        real(dp) :: u(3), h, e1, e2, xp(3), xm(3)
        integer :: k, i
        h = 1e-5_dp; e1 = 0; e2 = 0
        do i = 1, 6
            u = [0.2_dp + 0.13_dp*i, 0.8_dp*i, 0.3_dp]
            call chart%evaluate(u(1), u(2), u(3), pt)
            do k = 1, 3
                xp = cart(tokamak_position(chart, u(1) + h*d(k, 1), &
                    u(2) + h*d(k, 2), u(3) + h*d(k, 3)))
                xm = cart(tokamak_position(chart, u(1) - h*d(k, 1), &
                    u(2) - h*d(k, 2), u(3) - h*d(k, 3)))
                e1 = max(e1, maxval(abs((xp - xm)/(2*h) - pt%xu(:, k))))
                call chart%evaluate(u(1) + h*d(k, 1), u(2) + h*d(k, 2), &
                    u(3) + h*d(k, 3), pp_)
                call chart%evaluate(u(1) - h*d(k, 1), u(2) - h*d(k, 2), &
                    u(3) - h*d(k, 3), pm)
                e2 = max(e2, maxval(abs((pp_%xu - pm%xu)/(2*h) - pt%xuu(:, :, k))))
            end do
        end do
        print '(a,a,2es10.2)', label, ' jet errors x_u, x_uu', e1, e2
        if (e1 > 1e-8_dp .or. e2 > 1e-7_dp) error stop 'chart jets'
    end subroutine jets_check

    pure real(dp) function d(k, j)
        integer, intent(in) :: k, j
        d = merge(1.0_dp, 0.0_dp, k == j)
    end function d

    pure function cart(x) result(c)
        real(dp), intent(in) :: x(3)
        real(dp) :: c(3)
        c = [x(1)*cos(x(3)), x(1)*sin(x(3)), x(2)]
    end function cart

    function solovev_field(x) result(b)
        ! Exact B = (1/R)(psi_Z e_R - psi_R e_Z) + (F0/R) e_phi at (R, Z, phi).
        real(dp), intent(in) :: x(3)
        real(dp) :: b(3), r, z, pr, pz
        r = x(1); z = x(2)
        pr = s0*(r*(r*r - r0*r0)/2 + r*z*z/kappa**2)
        pz = s0*r*r*z/kappa**2
        b = pz/r*[cos(x(3)), sin(x(3)), 0.0_dp] &
            + f0/r*[-sin(x(3)), cos(x(3)), 0.0_dp] + [0.0_dp, 0.0_dp, -pr/r]
    end function solovev_field

    subroutine invariance()
        ! The m = 1 internal kink (q crosses 1 at rho ~ 0.6) in two different
        ! straight-field-line charts. Converged with m = -10..28 and 32 cells:
        ! omega^2 = -1.30328e-3 (PEST) and -1.30331e-3 (Hamada).
        real(dp) :: lp, lh, cp, ch
        call lowest_mode(pest, 16, lp, cp, -8, 18)
        call lowest_mode(hamada, 16, lh, ch, -8, 18)
        print '(a,2es16.8,a,2f4.0)', 'internal kink PEST, Hamada', lp, lh, &
            ' unstable', cp, ch
        if (cp < 1 .or. ch < 1) error stop 'internal kink missing'
        if (abs(lp/lh - 1) > 3e-3_dp) error stop 'PEST and Hamada disagree'
    end subroutine invariance

    subroutine lowest_mode(chart, ne, lowest, count, mlo, mhi)
        ! Sparse LDL^H inertia below -1e-6 and the eigenvalue nearest -1e-3.
        type(tokamak_chart_t), intent(in) :: chart
        integer, intent(in) :: ne, mlo, mhi
        real(dp), intent(out) :: lowest, count
        integer :: nthq
        type(flux_space_t) :: space
        type(csc_z_t) :: stiffness, mass
        type(fortsparse_status_t) :: status
        type(s0_pencil_t) :: pencil
        type(s0_report_t) :: report
        type(eigen_options_t) :: options
        complex(dp), allocatable :: start(:), vectors(:, :), omega2(:)
        real(dp), allocatable :: residuals(:), edges(:)
        integer :: ok, j
        nthq = 64
        edges = [(real(j, dp)/ne, j = 0, ne)]
        call flux_space_init(space, edges, 3, [(j, j = mlo, mhi)], &
            [(1, j = mlo, mhi)], 1, nthq, 1, ok)
        call flux_assemble(space, chart, 5.0_dp/3, 1.0_dp, stiffness, mass, status)
        if (status%code /= FORTSPARSE_OK) error stop 'assembly'
        call s0_prepare(pencil, stiffness, mass, cmplx(-1e-6_dp, 0, dp), report)
        count = report%below
        call s0_prepare(pencil, stiffness, mass, cmplx(-2e-3_dp, 0, dp), report)
        allocate(start(space%ndof), vectors(space%ndof, 1), omega2(1), residuals(1))
        do j = 1, space%ndof
            start(j) = cmplx(sin(1.3_dp*j), cos(0.7_dp*j), dp)
        end do
        options%max_basis = 24; options%max_restarts = 60
        options%which = EIGEN_NEAREST
        call s0_solve(pencil, start, 1e-10_dp, options, omega2, vectors, &
            residuals, report)
        lowest = real(omega2(1))
        if (report%status /= EIGEN_OK) lowest = huge(1.0_dp)
    end subroutine lowest_mode

end program test_ideal_flux_tokamak
