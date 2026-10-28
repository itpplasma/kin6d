module newton_mr_fixture
    ! The mapped vacuum oracle of test_nonnested_mapped (rotating ellipse,
    ! closed-form vacuum B) and a peak-memory probe.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 3, minor = 1, ftor = 3
    real(dp), parameter :: bv = 0.1_dp, coef(3) = [0.05_dp, 0.03_dp, 0.02_dp]
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

    integer function peak_kib() result(kib)
        ! VmHWM of this process (Linux), -1 if unavailable.
        character(256) :: line
        integer :: u, ios
        kib = -1
        open (newunit=u, file='/proc/self/status', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') read (line(7:), *) kib
        end do
        close (u)
    end function
end module newton_mr_fixture


program test_nonnested_newton_mr
    ! Usage: agree nr N | ladder pc N nr... |
    !        path factored|mr|assembled nr N
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use newton_mr_fixture
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build, nn_mapped_data, &
        nn_mapped_solve, nn_mapped_error, nn_mapped_report_t, map_metric_t
    use kin6d_nonnested_newton_mr, only: nn_mr_options_t, nn_mr_report_t, &
        nn_mapped_solve_mr, NN_MR_PC_NONE, NN_MR_PC_JACOBI, NN_MR_PC_LUMPED
    use kin6d_newton_mr, only: MR_FORCING_FIXED
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl_curl
    use fortfem_meridian_fourier_apply, only: meridian_fourier_op_t, &
        meridian_fourier_op_init, meridian_fourier_curl_curl_apply
    use fortfem_meridian_fourier_curl, only: meridian_fourier_curl
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, fortsparse_status_t
    implicit none
    character(32) :: mode, arg
    integer, allocatable :: list(:)
    integer :: i
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('agree')
        call read_ints(2, list)
        call agreement(list(1), list(2))
    case ('ladder')
        call get_command_argument(2, arg)
        call read_ints(3, list)
        call ladder(trim(arg), list(1), list(2:))
    case ('path')
        call get_command_argument(2, arg)
        call read_ints(3, list)
        call one_path(trim(arg), list(1), list(2))
    case default
        error stop 'usage: agree nr N | ladder pc N nr... | path p nr N'
    end select
contains
    subroutine read_ints(first, values)
        integer, intent(in) :: first
        integer, allocatable, intent(out) :: values(:)
        character(32) :: word
        allocate (values(command_argument_count() - first + 1))
        do i = 1, size(values)
            call get_command_argument(first + i - 1, word)
            read (word, *) values(i)
        end do
    end subroutine

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

    subroutine setup(nr, modes, ellipse, field, map, data, nphi)
        integer, intent(in) :: nr, modes
        logical, intent(in) :: ellipse
        type(free_b_field_t), intent(out), target :: field
        type(nn_map_t), intent(out), target :: map
        real(dp), allocatable, intent(out) :: data(:)
        integer, intent(out) :: nphi
        integer :: st
        nphi = 4*modes + 4
        call make_field(nr, modes, field)
        if (ellipse) then
            call nn_map_build(field%fs, 2, 8, ellipse_map, map, st)
        else
            call nn_map_build(field%fs, 0, 1, identity_map, map, st)
        end if
        if (st /= 0) error stop 'map'
        call nn_mapped_data(field, map, nphi, potential, data, st)
        if (st /= 0) error stop 'data'
    end subroutine

    subroutine agreement(nr, modes)
        ! (1) The matrix-free Hessian action equals the assembled S.
        ! (2) Newton-MR (lumped preconditioner, matrix-free) and the default
        ! assembled FGMRES path give the same B, energy and oracle error.
        ! (3) Identity map: Newton-MR reproduces the per-mode NN1 solves.
        integer, intent(in) :: nr, modes
        type(free_b_field_t), target :: f_ref, f_mr, f_nn1
        type(nn_map_t), target :: map
        type(map_metric_t) :: metric
        type(meridian_fourier_op_t) :: op
        type(csc_t) :: s, c, ct
        type(fortsparse_status_t) :: fst
        type(nn_mapped_report_t) :: rep_ref, rep_mr
        type(nn_mr_options_t) :: opt
        type(nn_mr_report_t) :: mr
        type(nn_report_t) :: r1
        real(dp), allocatable :: data(:), u(:), y1(:), y2(:)
        complex(dp), allocatable :: bc(:), seed(:)
        real(dp) :: err_ref, err_mr, dev, eh, ee, diff, de
        integer :: st, nphi, k, n, n1, nn
        call setup(nr, modes, .true., f_ref, map, data, nphi)
        metric%fs => f_ref%fs; metric%map => map
        call meridian_fourier_curl_curl(f_ref%fs, modes, nphi, metric, s, st)
        if (st /= 0) error stop 'assembled S'
        call meridian_fourier_curl(f_ref%fs, modes, c, st)
        call csc_transpose(c, ct, fst)
        call meridian_fourier_op_init(f_ref%fs, 2, modes, nphi, metric, op, st)
        if (st /= 0) error stop 'op init'
        u = [(sin(0.37_dp*k) + 0.1_dp*cos(1.3_dp*k), k=1, s%ncol)]
        allocate (y1(size(u)), y2(size(u)))
        y1 = csc_matvec(s, u)
        call meridian_fourier_curl_curl_apply(op, c, ct, u, y2)
        diff = norm2(y1 - y2)/norm2(y1)
        write (*, '(a,es11.3)') 'apply: |S_mf u - S u|/|S u| =', diff
        if (diff > 1.0e-12_dp) error stop 'matrix-free apply differs from assembled S'
        call nn_mapped_solve(f_ref, map, nphi, data, 1.0e-13_dp, rep_ref, st)
        if (st /= 0) error stop 'reference solve'
        f_mr = f_ref
        opt%newton%tolerance = 1.0e-12_dp
        call nn_mapped_solve_mr(f_mr, map, nphi, data, opt, rep_mr, mr, st)
        if (st /= 0) error stop 'Newton-MR solve'
        call nn_mapped_error(f_ref, map, 2*nphi + 1, exact_field, err_ref, dev, eh, ee, st)
        call nn_mapped_error(f_mr, map, 2*nphi + 1, exact_field, err_mr, dev, eh, ee, st)
        diff = maxval(abs(f_mr%b - f_ref%b))/maxval(abs(f_ref%b))
        de = abs(rep_mr%energy - rep_ref%energy)/abs(rep_ref%energy)
        write (*, '(a,i0,a,i0,a,es10.3,a,es10.3,a,2es10.3,a,2es10.3)') 'agree nr=', nr, &
            ' N=', modes, ' max|b_mr-b_ref|/max|b|=', diff, ' energy rel diff=', de, &
            ' field error ref/mr=', err_ref/dev, err_mr/dev, ' residual ref/mr=', &
            rep_ref%residual, rep_mr%residual
        write (*, '(a,i0,a,i0,a,i0,a,i0)') '  fgmres=', rep_ref%iterations, &
            ' newton=', mr%newton%iterations, ' minres=', mr%newton%minres_iterations, &
            ' npc=', mr%newton%npc_steps
        if (diff > 1.0e-9_dp) error stop 'Newton-MR B differs from the factored path'
        if (de > 1.0e-11_dp) error stop 'Newton-MR energy differs'
        if (abs(err_mr - err_ref) > 1.0e-9_dp*dev) error stop 'oracle error differs'
        if (rep_mr%residual > 1.0e-11_dp) error stop 'Newton-MR residual'
        if (rep_mr%divergence > 1.0e-12_dp) error stop 'Newton-MR divergence'
        ! Identity map: the per-n NN1 solves of the same traces.
        call setup(nr, modes, .false., f_mr, map, data, nphi)
        call make_field(nr, modes, f_nn1)
        call nn_mapped_solve_mr(f_mr, map, nphi, data, opt, rep_mr, mr, st)
        if (st /= 0) error stop 'identity Newton-MR'
        n1 = f_mr%fs%sp(1)%ndof; nn = f_mr%fs%nn
        bc = cmplx(data(nn + 1:n1), 0.0_dp, dp); seed = 0*bc
        call nn_solve(f_nn1, 0, bc, seed, r1, st)
        if (st /= 0) error stop 'NN1 n=0'
        do n = 1, modes
            bc = cmplx(data((2*n - 1)*n1 + 1:(2*n - 1)*n1 + nn), &
                -data(2*n*n1 + 1:2*n*n1 + nn), dp)/2
            seed = 0*bc
            call nn_solve(f_nn1, n, bc, seed, r1, st)
            if (st /= 0) error stop 'NN1 mode'
        end do
        diff = maxval(abs(f_mr%b - f_nn1%b))/maxval(abs(f_nn1%b))
        write (*, '(a,es10.3,a,i0,a,i0)') 'identity: max|b_mr-b_NN1|/max|b|=', diff, &
            ' newton=', mr%newton%iterations, ' minres=', mr%newton%minres_iterations
        if (diff > 1.0e-9_dp) error stop 'identity: Newton-MR differs from NN1'
        if (mr%newton%minres_iterations > 3*mr%newton%iterations) &
            error stop 'identity: lumped preconditioner not exact'
    end subroutine

    subroutine ladder(pc, modes, rings)
        ! MINRES iterations versus mesh size for one preconditioner (a solve
        ! that does not converge within max_minres is reported, not fatal).
        character(*), intent(in) :: pc
        integer, intent(in) :: modes, rings(:)
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: rep
        type(nn_mr_options_t) :: opt
        type(nn_mr_report_t) :: mr
        real(dp), allocatable :: data(:)
        integer :: k, st, nphi, its(size(rings))
        select case (pc)
        case ('none'); opt%preconditioner = NN_MR_PC_NONE
        case ('jacobi'); opt%preconditioner = NN_MR_PC_JACOBI
        case ('lumped'); opt%preconditioner = NN_MR_PC_LUMPED
        case default; error stop 'ladder: none | jacobi | lumped'
        end select
        ! One MINRES solve to the relative residual 1e-8 (fixed forcing,
        ! one Newton step): its iteration count measures the preconditioner.
        opt%newton%tolerance = 1.0e-8_dp
        opt%newton%forcing = MR_FORCING_FIXED
        opt%newton%forcing_max = 1.0e-8_dp
        opt%newton%max_iterations = 1
        opt%newton%max_minres = 5000
        do k = 1, size(rings)
            call setup(rings(k), modes, .true., field, map, data, nphi)
            call nn_mapped_solve_mr(field, map, nphi, data, opt, rep, mr, st)
            write (*, '(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,es9.2,a,f8.2,a,f8.2)') &
                'ladder pc=', pc, ' N=', modes, ' nr=', rings(k), ' ndof=', rep%ndof, &
                ' status=', mr%newton%status, ' minres=', &
                mr%newton%minres_iterations, ' gradient=', &
                mr%newton%gradient/mr%newton%initial_gradient, ' pc_s=', &
                mr%preconditioner_s, ' minres_s=', mr%newton%minres_s
            if (pc == 'lumped' .and. st /= 0) error stop 'ladder solve'
            its(k) = mr%newton%minres_iterations
        end do
        if (pc == 'lumped' .and. size(rings) > 1) then
            ! Mesh independence of the lumped preconditioner.
            if (its(size(its)) > 2*its(1) + 4) error stop 'lumped iterations grow with h'
        end if
    end subroutine

    subroutine one_path(path, nr, modes)
        ! Wall time and peak memory of one solve path in a fresh process.
        character(*), intent(in) :: path
        integer, intent(in) :: nr, modes
        type(free_b_field_t), target :: field
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: rep
        type(nn_mr_options_t) :: opt
        type(nn_mr_report_t) :: mr
        real(dp), allocatable :: data(:)
        real(dp) :: err, dev, eh, ee
        integer :: st, nphi, kib0
        integer(8) :: c0, c1, rate
        call setup(nr, modes, .true., field, map, data, nphi)
        kib0 = peak_kib()
        call system_clock(c0, rate)
        select case (path)
        case ('factored')
            call nn_mapped_solve(field, map, nphi, data, 1.0e-12_dp, rep, st)
        case ('mr', 'assembled')
            opt%matrix_free = path == 'mr'
            opt%newton%tolerance = 1.0e-12_dp
            call nn_mapped_solve_mr(field, map, nphi, data, opt, rep, mr, st)
        case default
            error stop 'path: factored | mr | assembled'
        end select
        call system_clock(c1)
        if (st /= 0) error stop 'path solve'
        call nn_mapped_error(field, map, 2*nphi + 1, exact_field, err, dev, eh, ee, st)
        write (*, '(a,a,a,i0,a,i0,a,i0,a,f8.3,a,i0,a,i0,a,i0,a,es10.3,a,es9.2)') &
            'path=', path, ' nr=', nr, ' N=', modes, ' ndof=', rep%ndof, ' solve_s=', &
            real(c1 - c0, dp)/rate, ' peak_MiB=', peak_kib()/1024, ' before_MiB=', &
            kib0/1024, ' iterations=', rep%iterations, ' field error=', err/dev, &
            ' residual=', rep%residual
        if (path /= 'factored') write (*, '(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,i0)') &
            '  operator_s=', mr%operator_s, ' pc_s=', mr%preconditioner_s, &
            ' minres_s=', mr%newton%minres_s, ' evaluate_s=', mr%newton%evaluate_s, &
            ' operator_MiB=', int(mr%operator_bytes/1048576)
    end subroutine
end program test_nonnested_newton_mr
