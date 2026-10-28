module beltrami_mr_fixture
    ! The synthetic helical Taylor case of test_nonnested_beltrami_spec
    ! (R0 = 6, unit circular section, rotating elliptic deformation of
    ! amplitude epsilon in 2 phi, toroidal flux 1, closed wall), without
    ! the external SPEC samples, and a peak-memory probe.
    use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), r0 = 6.0_dp
    real(dp) :: epsilon = 0.1_dp
contains
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine

    subroutine deform(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        real(dp) :: c, s
        c = epsilon*cos(2*phi); s = epsilon*sin(2*phi)
        y = [r0 + (1 + c)*(x(1) - r0) + s*x(2), s*(x(1) - r0) + (1 - c)*x(2)]
    end subroutine

    subroutine helical_case(k, nr, modes, field, map)
        integer, intent(in) :: k, nr, modes
        type(free_b_field_t), intent(out) :: field
        type(nn_map_t), intent(out) :: map
        type(flux_mesh_t) :: mesh
        real(dp), allocatable :: rings(:)
        integer :: i, st
        allocate (rings(nr)); rings = [(real(i, dp)/nr, i=1, nr)]
        call flux_mesh_build(rings, 8, 256, k, circle, mesh, st)
        if (st /= 0) error stop 'mesh'
        call nn_build(mesh, k, modes, [r0, 0.0_dp], 1.0_dp, 0.0_dp, 1.0_dp, field, st)
        if (st /= 0) error stop 'field'
        call nn_map_build(field%fs, 2, 9, deform, map, st)
        if (st /= 0) error stop 'map'
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

    real(dp) function wall()
        integer(i8) :: count, rate
        call system_clock(count, rate)
        wall = real(count, dp)/real(rate, dp)
    end function
end module beltrami_mr_fixture

program test_nonnested_beltrami_mr
    ! Matrix-free MINRES / Newton-MR Taylor closure against the factored
    ! path on the mapped helical case. Usage:
    !   agree k nr N               fields, energy, helicity, mu tangent and the
    !                              helicity constraint at mu = 1.1 and 4.5
    !   ladder pc k N mu eps shift nr...
    !                              MINRES iterations (pc: none|jacobi|block) at
    !                              deformation eps and preconditioner shift
    !                              (0: automatic)
    !   path factored|mr k nr N mu wall time and peak RSS of one solve
    !   inertia k nr N mu          negative eigenvalues of the A block at mu
    !                              (factored inertia at mu minus that at 0)
    use beltrami_mr_fixture
    use kin6d_nonnested_beltrami, only: beltrami_system_t, beltrami_report_t, &
        beltrami_prepare, beltrami_solve, beltrami_constraint, beltrami_mr_options_t, &
        BELTRAMI_PC_NONE, BELTRAMI_PC_JACOBI, BELTRAMI_PC_BLOCK
    implicit none
    character(32) :: mode, arg
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('agree')
        call agree()
    case ('ladder')
        call ladder()
    case ('path')
        call path()
    case ('inertia')
        call inertia()
    case default
        error stop 'usage: agree | ladder | path'
    end select
contains
    integer function iarg(i)
        integer, intent(in) :: i
        call get_command_argument(i, arg); read (arg, *) iarg
    end function

    real(dp) function rarg(i)
        integer, intent(in) :: i
        call get_command_argument(i, arg); read (arg, *) rarg
    end function

    subroutine prepare(field, map, modes, mr, system)
        type(free_b_field_t), intent(in) :: field
        type(nn_map_t), intent(in) :: map
        integer, intent(in) :: modes
        type(beltrami_mr_options_t), optional, intent(in) :: mr
        type(beltrami_system_t), intent(out) :: system
        real(dp), allocatable :: data(:)
        integer :: st
        allocate (data((2*modes + 1)*field%fs%sp(1)%ndof)); data = 0
        call beltrami_prepare(field, 1.0_dp, data, 4*modes + 9, system, st, &
            map=map, nfp=2, mr=mr)
        if (st /= 0) error stop 'prepare'
    end subroutine

    subroutine agree()
        ! Factored vs MINRES to the MINRES tolerance 1e-11. mu = 4.5 lies
        ! beyond the first Taylor eigenvalue of this section: the A block is
        ! itself indefinite there (see the printed inertia).
        type(free_b_field_t) :: field, fm
        type(nn_map_t) :: map
        type(beltrami_system_t) :: sf, sm
        type(beltrami_report_t) :: rf, rm
        type(beltrami_mr_options_t) :: mr
        real(dp), allocatable :: af(:), am(:), df(:), dm(:)
        real(dp) :: mus(2), db, de, dh, dd, target
        integer :: k, nr, modes, i
        k = iarg(2); nr = iarg(3); modes = iarg(4)
        call helical_case(k, nr, modes, field, map)
        mr%tolerance = 1.0e-11_dp
        call prepare(field, map, modes, system=sf)
        call prepare(field, map, modes, mr, sm)
        mus = [1.1_dp, 4.5_dp]
        do i = 1, 2
            fm = field
            call beltrami_solve(sf, mus(i), field, rf, potential=af, dmu=df)
            call beltrami_solve(sm, mus(i), fm, rm, potential=am, dmu=dm)
            if (rf%status /= 0 .or. rm%status /= 0) error stop 'agree solve'
            db = maxval(abs(field%b - fm%b))/maxval(abs(field%b))
            de = abs(rf%energy - rm%energy)/abs(rf%energy)
            dh = abs(rf%helicity - rm%helicity)/abs(rf%helicity)
            dd = norm2(df - dm)/norm2(df)
            print '(a,f4.1,a,3i7)', 'mu ', mus(i), ' factored inertia (-,0,+) ', &
                rf%negative, rf%zero, rf%positive
            print '(a,4es11.3)', '  rel. diff b / energy / helicity / dA/dmu ', &
                db, de, dh, dd
            print '(a,2es11.3,a,i5)', '  residual factored / MINRES ', rf%residual, &
                rm%residual, '  MINRES its ', rm%minres_iterations
            print '(a,2es11.3)', '  divergence factored / MINRES ', rf%divergence, &
                rm%divergence
            if (db > 1.0e-8_dp .or. de > 1.0e-9_dp .or. dh > 1.0e-8_dp .or. &
                dd > 1.0e-7_dp) error stop 'factored/MINRES agreement'
            if (rm%divergence > 2.0e-11_dp) error stop 'MINRES div curl identity'
            if (i == 1) target = rf%helicity
        end do
        ! Helicity constraint: factored bisection vs Newton-MR from mid-bracket
        ! (mu0 = 1.2) on the same target, recovering mu = 1.1.
        call beltrami_constraint(sf, 'helicity', target, [0.5_dp, 1.9_dp], &
            1.0e-10_dp, field, rf)
        fm = field
        call beltrami_constraint(sm, 'helicity', target, [0.5_dp, 1.9_dp], &
            1.0e-10_dp, fm, rm)
        print '(a,2es22.14)', 'helicity mu factored / Newton-MR ', rf%mu, rm%mu
        print '(a,i3,a,i5,a,es11.3)', '  Newton-MR steps ', rm%newton_iterations, &
            ' MINRES its ', rm%minres_iterations, ' |H-H0|/|H0| ', &
            abs(rm%helicity - target)/abs(target)
        if (rf%status /= 0 .or. rm%status /= 0) error stop 'constraint status'
        if (abs(rf%mu - 1.1_dp) > 1.0e-8_dp .or. abs(rm%mu - 1.1_dp) > 1.0e-8_dp) &
            error stop 'helicity multiplier'
        if (rm%newton_iterations < 2) error stop 'Newton-MR took no step'
    end subroutine

    subroutine ladder()
        ! MINRES iterations of one fixed-mu solve to relative residual 1e-8.
        type(free_b_field_t) :: field
        type(nn_map_t) :: map
        type(beltrami_system_t) :: system
        type(beltrami_report_t) :: report
        type(beltrami_mr_options_t) :: mr
        integer :: k, modes, i, nr, lo, hi, its(16), n
        real(dp) :: mu
        character(16) :: pc
        call get_command_argument(2, pc)
        k = iarg(3); modes = iarg(4); mu = rarg(5); epsilon = rarg(6)
        mr%shift = rarg(7)
        select case (trim(pc))
        case ('none')
            mr%preconditioner = BELTRAMI_PC_NONE
        case ('jacobi')
            mr%preconditioner = BELTRAMI_PC_JACOBI
        case default
            mr%preconditioner = BELTRAMI_PC_BLOCK
        end select
        mr%tolerance = 1.0e-8_dp; mr%max_minres = 20000
        n = command_argument_count() - 7
        do i = 1, n
            nr = iarg(7 + i)
            call helical_case(k, nr, modes, field, map)
            call prepare(field, map, modes, mr, system)
            call beltrami_solve(system, mu, field, report)
            its(i) = report%minres_iterations
            print '(a,a,a,i0,a,i0,a,f4.1,a,i3,a,i8,a,i6,a,es10.2,a,i0)', 'ladder ', &
                trim(pc), ' P', k, ' N=', modes, ' mu=', mu, ' nr=', nr, &
                ' V1 dofs=', (2*modes + 1)*field%fs%sp(1)%ndof, ' its=', its(i), &
                ' res=', report%residual, ' status=', report%status
        end do
        if (trim(pc) == 'block') then
            lo = minval(its(:n)); hi = maxval(its(:n))
            if (hi > 2*lo + 10) error stop 'block preconditioner not mesh independent'
        end if
    end subroutine

    subroutine inertia()
        type(free_b_field_t) :: field
        type(nn_map_t) :: map
        type(beltrami_system_t) :: system
        type(beltrami_report_t) :: r0, r1
        integer :: k, nr, modes
        real(dp) :: mu
        k = iarg(2); nr = iarg(3); modes = iarg(4); mu = rarg(5)
        call helical_case(k, nr, modes, field, map)
        call prepare(field, map, modes, system=system)
        call beltrami_solve(system, 0.0_dp, field, r0)
        call beltrami_solve(system, mu, field, r1)
        print '(a,i0,a,i0,a,f4.1,a,i0)', 'inertia nr=', nr, ' N=', modes, ' mu=', &
            mu, ' negative A-block eigenvalues=', r1%negative - r0%negative
    end subroutine

    subroutine path()
        ! One prepare + solve; wall time and peak RSS (run in a fresh process).
        type(free_b_field_t) :: field
        type(nn_map_t) :: map
        type(beltrami_system_t) :: system
        type(beltrami_report_t) :: report
        type(beltrami_mr_options_t) :: mr
        character(16) :: which
        integer :: k, nr, modes, base
        real(dp) :: mu, t0, t1
        call get_command_argument(2, which)
        k = iarg(3); nr = iarg(4); modes = iarg(5); mu = rarg(6)
        call helical_case(k, nr, modes, field, map)
        base = peak_kib()
        t0 = wall()
        if (trim(which) == 'mr') then
            call prepare(field, map, modes, mr, system)
        else
            call prepare(field, map, modes, system=system)
        end if
        t1 = wall()
        call beltrami_solve(system, mu, field, report)
        if (report%status /= 0) error stop 'path solve'
        print '(a,a,a,i0,a,i0,a,i8,a,f8.2,a,f8.2,a,i7,a,i7,a,i5,a,es10.2,a,es14.6)', &
            'path ', trim(which), ' nr=', nr, ' N=', modes, ' V1 dofs=', &
            (2*modes + 1)*field%fs%sp(1)%ndof, ' prepare_s=', t1 - t0, &
            ' solve_s=', wall() - t1, ' peak_MiB=', peak_kib()/1024, &
            ' before_MiB=', base/1024, ' its=', report%minres_iterations, &
            ' res=', report%residual, ' energy=', report%energy
    end subroutine
end program test_nonnested_beltrami_mr
