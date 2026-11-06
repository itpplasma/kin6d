!> Boozer `.bc` geometry: closed forms on a synthetic file, refusal outside
!! the tabulated surfaces, JVP/VJP pairing of the metric leaf, and agreement
!! with NEO-2 on a public tokamak file when its path and NEO-2's values are
!! given (args: file s iota diota_ds av_nabla_stor[1/cm] avb2[T^2]).
program test_boozer_file
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_boozer_file
    use kin6d_boozer_file_leaves, only: boozer_metric, boozer_metric_jvp, &
        boozer_metric_vjp
    implicit none
    real(dp), parameter :: r0 = 1.7_dp, a = 0.5_dp, b0 = 2.3_dp, eps = 0.1_dp
    integer :: failures = 0

    call synthetic_file()
    call metric_pairing()
    if (command_argument_count() >= 6) call neo2_comparison()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'
contains

    subroutine synthetic_file()
        character(*), parameter :: name = 'synthetic.bc'
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(boozer_point_t) :: p
        real(dp) :: s, theta, grad_s, b2, delta
        integer :: unit, k, status, i
        real(dp), parameter :: ss(9) = [(0.1_dp*k, k=1, 9)]

        ! R = r0 + a s cos(theta), Z = a s sin(theta): s = r/a, |grad s| = 1/a.
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 9 1 -1.0 0.5 1.7'
        do k = 1, 9
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') ss(k), 0.9_dp - 0.3_dp*ss(k), 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, b0
            write (unit, '(2i5,4es24.16)') 1, 0, a*ss(k), a*ss(k), 0.0_dp, &
                -eps*b0
        end do
        close (unit)
        call load_boozer_file(name, bf)

        do i = 1, 2
            s = merge(0.5_dp, 0.43_dp, i == 1)
            call boozer_surface(bf, s, surf, status)
            call check('surface status', status == BOOZER_OK)
            call close_to('iota', surf%iota, 0.9_dp - 0.3_dp*s, 1e-14_dp)
            call close_to('diota_ds', surf%diota_ds, -0.3_dp, 1e-13_dp)
            theta = 0.7_dp
            p = boozer_point(surf, theta, 0.3_dp)
            call close_to('R', p%r, r0 + a*s*cos(theta), 1e-14_dp)
            call close_to('|grad s|', p%grad_s, 1/a, 1e-13_dp)
            call close_to('jacobian', p%jacobian, -a*a*s*p%r, 1e-13_dp)
            call boozer_average(surf, 16, 1, grad_s, b2)
            delta = a*s/r0
            call close_to('<|grad s|>', grad_s, 1/a, 1e-13_dp)
            call close_to('<B^2>', b2, b0**2*(1 - eps*delta + eps**2/2), 1e-13_dp)
        end do

        call boozer_surface(bf, 0.05_dp, surf, status)
        call check('refuse s below the file', status == BOOZER_OUTSIDE)
        call boozer_surface(bf, 0.95_dp, surf, status)
        call check('refuse s above the file', status == BOOZER_OUTSIDE)
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine synthetic_file

    subroutine metric_pairing()
        real(dp) :: x(10), v(10), g(10), w(2), jv(2), lhs, rhs
        integer :: k
        call random_seed(put=[(17 + k, k=1, 64)])
        call random_number(x); call random_number(v); call random_number(w)
        x(1) = x(1) + 1
        call boozer_metric_jvp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), v(1), v(2), v(3), v(4), &
            v(5), v(6), v(7), v(8), v(9), v(10), jv(1), jv(2))
        call boozer_metric_vjp(x(1), x(2), x(3), x(4), x(5), x(6), x(7), &
            x(8), x(9), x(10), w(1), w(2), g(1), g(2), &
            g(3), g(4), g(5), g(6), g(7), g(8), g(9), g(10))
        lhs = dot_product(w, jv)
        rhs = dot_product(g, v)
        call close_to('metric JVP/VJP pairing', lhs, rhs, 1e-13_dp)
    end subroutine metric_pairing

    subroutine neo2_comparison()
        character(512) :: file
        real(dp) :: ref(5), grad_s, b2
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        integer :: k, status
        character(64) :: arg

        call get_command_argument(1, file)
        do k = 1, 5
            call get_command_argument(k + 1, arg)
            read (arg, *) ref(k)
        end do
        call load_boozer_file(trim(file), bf)
        call boozer_surface(bf, ref(1), surf, status)
        call check('NEO-2 surface status', status == BOOZER_OK)
        call boozer_average(surf, 128, 1, grad_s, b2)
        print '(a,4es10.2)', 'relative deviation from NEO-2 (iota, diota/ds, <|grad s|>, <B^2>):', &
            surf%iota/ref(2) - 1, surf%diota_ds/ref(3) - 1, grad_s/(100*ref(4)) - 1, b2/ref(5) - 1
        call close_to('NEO-2 iota', surf%iota, ref(2), 1e-6_dp)
        call close_to('NEO-2 diota/ds', surf%diota_ds, ref(3), 1e-4_dp)
        ! The radial derivatives of the file's coefficients are uncertain by
        ! about 1e-3 (shifting the four-point stencil by one surface moves
        ! <|grad s|> by 4e-4 to 9e-4); NEO-2 uses a global spline instead.
        call close_to('NEO-2 <|grad s|>', grad_s/100, ref(4), 1e-3_dp)
        call close_to('NEO-2 <B^2>', b2, ref(5), 1e-4_dp)
    end subroutine neo2_comparison

    subroutine close_to(label, got, want, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, want, rtol
        if (abs(got - want) > rtol*abs(want)) then
            print '(a,2es24.16)', 'FAIL '//label//': ', got, want
            failures = failures + 1
        end if
    end subroutine close_to

    subroutine check(label, ok)
        character(*), intent(in) :: label
        logical, intent(in) :: ok
        if (.not. ok) then
            print '(a)', 'FAIL '//label
            failures = failures + 1
        end if
    end subroutine check
end program test_boozer_file
