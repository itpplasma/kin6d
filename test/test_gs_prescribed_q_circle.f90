program test_gs_prescribed_q_circle
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, sample_gs, &
                                    write_gs_mesh, read_gs_mesh
    use kin6d_gs_prescribed_q, only: solve_gs_prescribed_q
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_primitives, &
                                    write_gs_profile_law, read_gs_profile_law
    use kin6d_gs_readback, only: gs_readback_t, read_gs_equilibrium
    implicit none
    real(dp), parameter :: r0 = 6.2_dp, a = r0/10, b0 = 5.3_dp
    real(dp), parameter :: mu0 = 1.2566370614359173e-6_dp
    type(gs_state_t) :: state, negative, loaded
    type(gs_profile_law_t) :: law, reverse_law, loaded_law
    type(gs_readback_t) :: readback
    real(dp), allocatable :: boundary(:, :)
    real(dp) :: s(4), q(4), phi, pprime, qerror, h, p, f
    real(dp) :: values(4, 160, 3), point(2), psi, br, bz, bt, error(4, 2), rates(4)
    integer :: mesh, n, i, status, family, k, evaluations, unit
    integer(int64) :: start, finish, clock_rate
    character(8) :: arg
    character(32) :: path
    call get_command_argument(1, arg)
    read (arg, *) family
    s = [0.0_dp, 1.0_dp/3, 2.0_dp/3, 1.0_dp]; q = 1.5_dp
    phi = acos(-1.0_dp)*a*a*b0
    pprime = 0
    if (family == 2) pprime = -2*b0/(1.5_dp*mu0*r0*r0)
    if (command_argument_count() == 2) then
        open (newunit=unit, file='boundary.dat', status='replace')
        write (unit, *) 512
        do i = 1, 512
            call circle(real(i - 1, dp)/512, point)
            write (unit, '(2es25.16)') point
        end do
        close (unit)
        open (newunit=unit, file='q.dat', status='replace')
        write (unit, *) size(s)
        do i = 1, size(s)
            write (unit, '(2es25.16)') s(i), q(i)
        end do
        close (unit)
        open (newunit=unit, file='input.nml', status='replace')
        write (unit, '(a)') "&gs_phase1 boundary_file='boundary.dat', q_file='q.dat',"
        write (unit, '(a)') ' degree=3, boundary_nodes=25, inverse_intervals=2,'
        write (unit, '(a)') ' radial_points=51, nonlinear_tolerance=1e-10,'
        write (unit, '(a,es25.16,a,es25.16,a)') ' phi_target=', phi, &
            ', pprime=', pprime, ' /'
        close (unit)
        stop
    end if
    do mesh = 1, 3
        n = nint(24*1.5_dp**(mesh - 1)) + 1
        allocate (boundary(2, n))
        do i = 1, n
            call circle(real(i - 1, dp)/n, boundary(:, i))
        end do
        h = 2*acos(-1.0_dp)*a/n
        call initialize_gs(state, boundary, h*h/2, status, 3, curve=circle)
        if (status /= 0) error stop 'circle mesh'
        call system_clock(start, clock_rate)
        call solve_gs_prescribed_q(state, mu0, pprime, s, q, phi, law, status, &
                             intervals=1, radial_points=2*n + 1, tolerance=1.0e-10_dp, &
                                   q_error=qerror, evaluations=evaluations)
        print *, 'E, n, status, q, evaluations:', family, n, status, qerror, evaluations
        call system_clock(finish)
        print *, 'inverse solve seconds:', real(finish - start, dp)/clock_rate
        if (real(finish - start, dp)/clock_rate > 300) error stop 'solve time limit'
        if (status /= 0) error stop 'constant q solve'
        call read_gs_equilibrium(state, law%f_edge, 0.0_dp, 4*n + 1, 160, &
                                 readback, status, law)
        if (status /= 0) error stop 'independent circle readback'
        if (abs(readback%phi_edge/phi - 1) > 1.0e-9_dp) error stop 'physical Phi'
        do i = 1, 160
            call circle(mod(i*0.6180339887498948_dp, 1.0_dp), point)
            point = [r0, 0.0_dp] + sqrt(real(i, dp)/180)*(point - [r0, 0.0_dp])
            call sample_gs(state, point(1), point(2), 1.0_dp, psi, br, bz, bt, status)
            if (status /= 0) error stop 'sample circle'
            call gs_law_primitives(law, psi, p, f, status)
            if (status /= 0) error stop 'circle F'
            values(:, i, mesh) = [psi, br, bz, f/point(1)]
        end do
        if (mesh == 1) then
            negative = state
            call solve_gs_prescribed_q(negative, mu0, pprime, s, -q, -phi, &
                              reverse_law, status, intervals=1, radial_points=2*n + 1, &
                                       tolerance=1.0e-10_dp)
            if (status /= 0) error stop 'negative q and Phi'
            if (maxval(abs(negative%psi - state%psi)) > 1.0e-9_dp) &
                error stop 'F reversal changed poloidal state'
            write (path, '(a,i0,a)') 'inverse_circle_', family, '.dat'
            call write_gs_mesh(negative, trim(path), status)
            if (status /= 0) error stop 'inverse mesh export'
            call read_gs_mesh(loaded, trim(path), status)
            if (status /= 0) error stop 'inverse mesh import'
            call write_gs_profile_law(trim(path)//'.law', reverse_law)
            call read_gs_profile_law(trim(path)//'.law', loaded_law, status)
            if (status /= 0) error stop 'inverse source import'
            call sample_gs(loaded, point(1), point(2), 1.0_dp, psi, br, bz, bt, status)
            if (status /= 0) error stop 'signed mesh readback'
            call gs_law_primitives(loaded_law, psi, p, f, status)
            if (status /= 0) error stop 'signed F readback'
            if (abs(f/point(1) + values(4, 160, mesh)) > 1.0e-8_dp) &
                error stop 'signed toroidal field round trip'
            call read_gs_equilibrium(negative, reverse_law%f_edge, 0.0_dp, &
                                     4*n + 1, 160, readback, status, reverse_law)
            if (status /= 0 .or. any(readback%q >= 0)) error stop 'signed q round trip'
            call solve_gs_prescribed_q(negative, mu0, pprime, s, -q, phi, &
                                       reverse_law, status)
            if (status == 0) error stop 'inconsistent q/Phi sign accepted'
            if (family == 1) then
                call solve_gs_prescribed_q(negative, mu0, pprime, s, q, phi, &
                              reverse_law, status, intervals=1, radial_points=2*n + 1, &
                                           q_tolerance=1.0e-8_dp)
                if (status /= 9) error stop 'underresolved q profile accepted'
                call read_gs_equilibrium(negative, reverse_law%f_edge, 0.0_dp, &
                                         4*n + 1, 160, readback, status, reverse_law)
                if (status /= 0) error stop 'rejected inverse lost its candidate'
                if (abs(readback%phi_edge/phi - 1) > 1.0e-9_dp) &
                    error stop 'rejected inverse lost physical normalization'
            end if
        end if
        deallocate (boundary)
    end do
    do k = 1, 4
        error(k, 1) = norm2(values(k, :, 2) - values(k, :, 1))
        error(k, 2) = norm2(values(k, :, 3) - values(k, :, 2))
    end do
    rates = log(error(:, 1)/error(:, 2))/log(1.5_dp)
    print *, 'self convergence psi BR BZ Bphi:', rates
    if (rates(1) < 3.3_dp .or. minval(rates(2:3)) < 2.5_dp) &
        error stop 'constant-q P3 field self convergence'
    if (rates(4) < 2.5_dp) error stop 'constant-q F self convergence'
contains
    subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [r0 + a*cos(2*acos(-1.0_dp)*t), a*sin(2*acos(-1.0_dp)*t)]
    end subroutine
end program
