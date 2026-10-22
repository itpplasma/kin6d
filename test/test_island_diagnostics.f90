module island_diagnostic_fixture
    ! Independent Hamiltonian oracle in I=r^2/2, theta:
    ! H=I/2+s(I-Is)^2/2-eps cos(2 theta-phi-alpha).
    ! In the rotating frame this is an exact pendulum. This is a diagnostic
    ! oracle, not an NN2 solution or a tokamak equilibrium.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, EQ_BPHI
    implicit none
    real(dp), parameter :: centre(2) = [10.0_dp, 0.0_dp], is = 0.25_dp, shear = 2
    type, extends(equilibrium_field_t) :: pendulum_t
        real(dp) :: eps = 1.0e-3_dp, alpha = 0.6_dp
        logical :: zero_toroidal = .false.
    contains
        procedure :: sample
    end type
contains
    subroutine sample(field, points, values, statuses)
        class(pendulum_t), intent(in) :: field
        real(dp), intent(in) :: points(:, :)
        real(dp), intent(out) :: values(:, :)
        integer, intent(out) :: statuses(:)
        real(dp) :: x, z, r, theta, phi, vr, vt
        integer :: j
        values = 0; statuses = 0
        do j = 1, size(points, 2)
            x = points(1, j) - centre(1); z = points(2, j)
            r = sqrt(x*x + z*z)
            if (r < 0.2_dp .or. r > 1.3_dp) then
                statuses(j) = 1
                cycle
            end if
            theta = atan2(z, x); phi = points(3, j)
            vr = -2*field%eps*sin(2*theta - phi - field%alpha)/r
            vt = 0.5_dp + shear*(r*r/2 - is)
            values(EQ_BR, j) = (vr*cos(theta) - r*vt*sin(theta))/points(1, j)
            values(EQ_BZ, j) = (vr*sin(theta) + r*vt*cos(theta))/points(1, j)
            if (.not. field%zero_toroidal) values(EQ_BPHI, j) = 1
        end do
    end subroutine
end module

program test_island_diagnostics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use island_diagnostic_fixture, only: pendulum_t, centre, is, shear
    use kin6d_island_diagnostics, only: island_trace_options_t, island_point_t, &
                       island_width_t, island_trace, island_periodic_point, island_width
    implicit none
    type(pendulum_t) :: field
    type(island_trace_options_t) :: opt
    type(island_point_t) :: op, xp
    type(island_width_t) :: width
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: r, theta, start(2), finish(2), expected, exact_edges(2), delta
    real(dp) :: widths(3), eps(3), exponent, qseed(2)
    integer :: st, j
    r = sqrt(2*is); theta = field%alpha/2
    start = centre + r*[cos(theta), sin(theta)]
    opt%root_tol = 1.0e-8_dp
    call island_trace(field, start, 0.0_dp, 1.0_dp, opt, finish, st)
    if (st /= 0) error stop 'trace'
    if (norm2(finish - (2*centre - start)) > 1.0e-7_dp) error stop 'one turn'
    call island_periodic_point(field, start + [0.002_dp, -0.003_dp], &
                               0.0_dp, 2, centre, 2, 1, opt, op)
    if (op%status /= 0) error stop 'two-turn O root'
    expected = (1 - cos(4*pi*2*sqrt(shear*field%eps)))/2
    print '(a,4es13.5)', 'O position/closure/residue error/det defect ', &
        norm2(op%x - start), op%closure, abs(op%residue - expected), &
        abs(op%determinant - 1)
    if (norm2(op%x - start) > 1.0e-6_dp) error stop 'O position'
    if (abs(op%phase - field%alpha) > 1.0e-6_dp) error stop 'phase'
    if (abs(op%residue - expected) > 1.0e-4_dp) error stop 'Greene residue'
    if (abs(op%determinant - 1) > 1.0e-4_dp) error stop 'flux preservation'
    theta = (pi + field%alpha)/2
    start = centre + r*[cos(theta), sin(theta)]
    call island_periodic_point(field, start, 0.0_dp, 2, centre, 2, 1, opt, xp)
    if (xp%status /= 0 .or. xp%residue >= 0) error stop 'X root'
    expected = (1 - cosh(4*pi*2*sqrt(shear*field%eps)))/2
    if (abs(xp%residue - expected) > 1.0e-4_dp) error stop 'hyperbolic residue'
    eps = [1.0e-3_dp, 2.5e-4_dp, 6.25e-5_dp]
    do j = 1, 3
        field%eps = eps(j)
        delta = 2*sqrt(field%eps/shear)
        exact_edges = sqrt(2*[is - delta, is + delta])
        theta = field%alpha/2
        call island_width(field, [r - 0.18_dp, r + 0.18_dp], r, centre, &
                          [cos(theta), sin(theta)], 0.0_dp, &
                          64.0_dp*sqrt(eps(1)/eps(j)), 2, 1, opt, &
                          action_label, 'I = r^2/2 [m^2]', width, scans=16)
        if (width%status /= 0) error stop 'width diagnostic'
        widths(j) = width%metres
        print '(a,5es13.5)', 'eps/width/error/trace/duration ', eps(j), widths(j), &
            abs(widths(j) - (exact_edges(2) - exact_edges(1))), &
            width%trace_error, width%duration_error
        if (maxval(abs(width%edge - exact_edges)) > 2*width%scan_error) &
            error stop 'pendulum separatrix width'
        if (abs(width%label_width - 2*delta) > 2*width%scan_error) &
            error stop 'declared radial label'
    end do
    exponent = log(widths(1)/widths(3))/log(eps(1)/eps(3))
    print '(a,f10.6)', 'pendulum width exponent ', exponent
    if (abs(exponent - 0.5_dp) > 0.05_dp) error stop 'pendulum scaling'
    ! A symmetry seed is not a full-map root unless the field is reversible
    ! at this section. A failed fallback must retain unresolved status.
    field%alpha = 0.6_dp; qseed = centre + [r, 0.0_dp]
    opt%max_steps = 1
    call island_periodic_point(field, qseed, 0.0_dp, 2, centre, 2, 1, opt, op, &
                               symmetry_bracket=[10.6_dp, 10.8_dp])
    if (op%status == 0) error stop 'failed trace accepted as root'
    opt%max_steps = 200000
    field%zero_toroidal = .true.
    call island_trace(field, start, 0.0_dp, 1.0_dp, opt, finish, st)
    if (st == 0) error stop 'zero toroidal field accepted'
    field%zero_toroidal = .false.
    call island_trace(field, centre, 0.0_dp, 1.0_dp, opt, finish, st)
    if (st == 0) error stop 'domain failure ignored'
    print '(a)', 'PASS: reusable N-turn island diagnostics against pendulum oracle'
contains
    function action_label(x) result(label)
        real(dp), intent(in) :: x(2)
        real(dp) :: label
        label = sum((x - centre)**2)/2
    end function
end program
