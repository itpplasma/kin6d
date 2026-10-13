module nn_mapped_fixture
    ! Closed-form vacuum oracle: B = F grad phi + B_v e_Z + sum_m grad u_m,
    ! u_m = c_m R0 (R/R0)^m cos(m phi) = Re of a holomorphic function of
    ! x + i y, so curl B = 0 and div B = 0 everywhere. Potential
    ! A = (0, B_v R/2, -F ln R + sum_m c_m R0 (R/R0)^m sin(m phi)).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3, minor = 1, ftor = 3
    real(dp), parameter :: bv = 0.1_dp, coef(3) = [0.05_dp, 0.03_dp, 0.02_dp]
    ! Rotating ellipse: semi-axes e1, e2 (equal area), turning with angle phi
    ! (the cross-section has period pi, two field periods).
    real(dp), parameter :: e1 = 1.25_dp, e2 = 0.8_dp
contains
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + minor*u*cos(theta), minor*u*sin(theta)]
    end subroutine

    subroutine identity_map(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        y = x + 0*phi
    end subroutine

    subroutine ellipse_map(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        real(dp) :: q(2, 2), d(2)
        q = reshape([cos(phi), sin(phi), -sin(phi), cos(phi)], [2, 2])
        d = matmul(transpose(q), x - [r0, 0.0_dp])
        y = [r0, 0.0_dp] + matmul(q, [e1*d(1), e2*d(2)])
    end subroutine

    subroutine potential(x, a)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3)
        integer :: m
        a = [0.0_dp, bv*x(1)/2, -ftor*log(x(1))]
        do m = 1, size(coef)
            a(3) = a(3) + coef(m)*r0*(x(1)/r0)**m*sin(m*x(3))
        end do
    end subroutine

    subroutine exact_field(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        integer :: m
        b = [0.0_dp, ftor/x(1), bv]
        do m = 1, size(coef)
            b(1) = b(1) + coef(m)*m*(x(1)/r0)**(m - 1)*cos(m*x(3))
            b(2) = b(2) - coef(m)*m*(x(1)/r0)**(m - 1)*sin(m*x(3))
        end do
    end subroutine
end module nn_mapped_fixture

program test_nonnested_mapped
    ! Usage: identity | hladder N nr... | modes nr N...
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_mapped_fixture
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build, nn_mapped_data, &
        nn_mapped_solve, nn_mapped_error, nn_mapped_report_t
    implicit none
    character(32) :: mode, arg
    integer, allocatable :: list(:)
    integer :: i, fixed_arg
    real(dp), allocatable :: errs(:)
    call get_command_argument(1, mode)
    if (command_argument_count() > 1) then
        call get_command_argument(2, arg)
        read (arg, *) fixed_arg
        allocate (list(command_argument_count() - 2))
        do i = 1, size(list)
            call get_command_argument(i + 2, arg)
            read (arg, *) list(i)
        end do
    end if
    select case (trim(mode))
    case ('identity')
        call identity_regression()
    case ('hladder')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call ellipse_case(list(i), fixed_arg, errs(i))
            if (i > 1) then
                write (*, '(a,i0,a,f6.3)') 'order nr=', list(i), ': ', &
                    log(errs(i - 1)/errs(i))/log(real(list(i), dp)/list(i - 1))
            end if
        end do
        if (size(errs) > 1) then
            if (log(errs(size(errs) - 1)/errs(size(errs))) &
                /log(real(list(size(list)), dp)/list(size(list) - 1)) < 1.7_dp) &
                error stop 'mapped h order below 1.7 for P2'
        end if
    case ('modes')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call ellipse_case(fixed_arg, list(i), errs(i))
        end do
        if (size(errs) > 1) then
            if (errs(size(errs)) > 0.5_dp*errs(1)) error stop 'mode ladder stalls'
        end if
    case default
        error stop 'usage: identity | hladder N nr... | modes nr N...'
    end select
contains
    subroutine make_field(nr, modes, field)
        integer, intent(in) :: nr, modes
        type(free_b_field_t), intent(out) :: field
        type(flux_mesh_t) :: fm
        real(dp) :: rings(nr)
        integer :: j, st
        do j = 1, nr
            rings(j) = real(j, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 4096, 2, circle, fm, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(fm, 2, modes, [r0, 0.0_dp], minor, ftor, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
    end subroutine

    subroutine identity_regression()
        ! Identity phi-dependence: the coupled FGMRES solve must reproduce
        ! the per-n NN1 solves for the same traces to rounding.
        integer, parameter :: nr = 2, modes = 3, nphi = 4*modes + 4
        type(free_b_field_t), target :: mapped, ref
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        type(nn_report_t) :: r1
        real(dp), allocatable :: data(:)
        complex(dp), allocatable :: bc(:), seed(:)
        real(dp) :: diff, err, dev, eh, ee, err_ref
        integer :: st, n, n1, nn, n0
        call make_field(nr, modes, mapped)
        call make_field(nr, modes, ref)
        call nn_map_build(mapped%fs, 0, 1, identity_map, map, st)
        if (st /= 0) error stop 'identity map'
        call nn_mapped_data(mapped, map, nphi, potential, data, st)
        if (st /= 0) error stop 'identity data'
        call nn_mapped_solve(mapped, map, nphi, data, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'identity mapped solve'
        n1 = mapped%fs%sp(1)%ndof; nn = mapped%fs%nn; n0 = mapped%fs%n0
        bc = cmplx(data(nn + 1:n1), 0.0_dp, dp)
        seed = 0*bc
        call nn_solve(ref, 0, bc, seed, r1, st)
        if (st /= 0) error stop 'NN1 n=0'
        do n = 1, modes
            bc = cmplx(data((2*n - 1)*n1 + 1:(2*n - 1)*n1 + nn), &
                -data(2*n*n1 + 1:2*n*n1 + nn), dp)/2
            seed = 0*bc
            call nn_solve(ref, n, bc, seed, r1, st)
            if (st /= 0) error stop 'NN1 mode'
        end do
        diff = maxval(abs(mapped%b - ref%b))/maxval(abs(ref%b))
        call nn_mapped_error(mapped, map, 2*nphi + 1, exact_field, err, dev, eh, ee, st)
        if (st /= 0) error stop 'identity error'
        call nn_mapped_error(ref, map, 2*nphi + 1, exact_field, err_ref, dev, eh, ee, st)
        if (st /= 0) error stop 'identity error ref'
        write (*, '(a,es10.3,a,i0,a,es10.3,a,es10.3,a,es10.3)') &
            'identity: max|b-b_NN1|/max|b|=', diff, ' iterations=', &
            report%iterations, ' divergence=', report%divergence, &
            ' field error=', err/dev, ' NN1=', err_ref/dev
        if (diff > 1.0e-10_dp) error stop 'identity map does not reproduce NN1'
        if (report%iterations > 3) error stop 'identity preconditioner not exact'
        if (report%divergence > 1.0e-12_dp) error stop 'identity divergence'
        if (abs(err - err_ref) > 1.0e-9_dp*dev) error stop 'identity error mismatch'
    end subroutine

    subroutine ellipse_case(nr, modes, rel)
        integer, intent(in) :: nr, modes
        real(dp), intent(out) :: rel
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: data(:)
        real(dp) :: err, dev, eh, ee
        integer :: st, nphi, c0, c1, rate
        nphi = 4*modes + 4
        call system_clock(c0, rate)
        call make_field(nr, modes, field)
        call nn_map_build(field%fs, 2, 8, ellipse_map, map, st)
        if (st /= 0) error stop 'ellipse map'
        call nn_mapped_data(field, map, nphi, potential, data, st)
        if (st /= 0) error stop 'ellipse data'
        call nn_mapped_solve(field, map, nphi, data, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'ellipse solve'
        call nn_mapped_error(field, map, 2*nphi + 1, exact_field, err, dev, eh, ee, st)
        if (st /= 0) error stop 'ellipse error'
        call system_clock(c1)
        rel = err/dev
        write (*, '(a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,i0,a,es9.2,a,es9.2,a,f8.2)') &
            'ellipse nr=', nr, ' N=', modes, ' ndof=', report%ndof, &
            ' field error=', rel, ' energy rel=', abs(report%energy - ee)/ee, &
            ' fgmres=', report%iterations, ' residual=', report%residual, &
            ' div=', report%divergence, ' s=', real(c1 - c0, dp)/rate
        if (report%divergence > 1.0e-12_dp) error stop 'mapped divergence identity'
        if (report%residual > 1.0e-9_dp) error stop 'mapped residual'
    end subroutine
end program test_nonnested_mapped
