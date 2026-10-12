program test_equilibrium_core
    ! Common equilibrium core (docs/equilibrium-core.md).
    !   chart              jet of the cylindrical chart against finite differences
    !   compare WORK FIX   K6-mini regression: gs_phase1 outputs in WORK against
    !                      the fixture FIX to 1e-12 relative
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    implicit none
    character(512) :: mode, work, fixtures
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('chart')
        call test_chart()
    case ('compare')
        call get_command_argument(2, work)
        call get_command_argument(3, fixtures)
        call compare_k6_mini(trim(work), trim(fixtures))
    case default
        error stop 'mode: chart | compare WORK FIXTURES'
    end select
contains
    subroutine test_chart()
        ! Second and third q-derivatives are the derivatives of the lower ones.
        type(cylindrical_chart_t) :: chart
        real(dp) :: q(3), x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        real(dp) :: xp(3), d1p(3, 3), d2p(3, 3, 3), d3p(3, 3, 3, 3)
        real(dp) :: xm(3), d1m(3, 3), d2m(3, 3, 3), d3m(3, 3, 3, 3), dq(3)
        real(dp) :: error, h
        integer :: k, sample
        error = 0
        h = 1.0e-5_dp
        do sample = 1, 4
            q = [6.2_dp + 0.3_dp*sample, -0.5_dp + 0.4_dp*sample, 0.7_dp*sample]
            call chart%jet(q, x, d1, d2, d3)
            if (maxval(abs(x - [q(1)*cos(q(3)), q(1)*sin(q(3)), q(2)])) > 1.0e-14_dp) &
                error stop 'cylindrical chart position'
            do k = 1, 3
                dq = 0; dq(k) = h
                call chart%jet(q + dq, xp, d1p, d2p, d3p)
                call chart%jet(q - dq, xm, d1m, d2m, d3m)
                error = max(error, maxval(abs((xp - xm)/(2*h) - d1(:, k))), &
                    maxval(abs((d1p - d1m)/(2*h) - d2(:, :, k))), &
                    maxval(abs((d2p - d2m)/(2*h) - d3(:, :, :, k))))
            end do
        end do
        print '(a,es10.3)', 'cylindrical chart jet finite-difference defect: ', error
        if (error > 1.0e-8_dp) error stop 'cylindrical chart jet inconsistent'
        print '(a)', 'PASS: chart jet'
    end subroutine

    subroutine compare_k6_mini(work, fixtures)
        character(*), intent(in) :: work, fixtures
        real(dp), parameter :: tolerance = 1.0e-12_dp
        real(dp) :: worst
        worst = 0
        call compare_summary(work//'/summary.dat', &
            fixtures//'/reference_summary.dat', worst)
        call compare_table(work//'/profiles.dat', &
            fixtures//'/reference_profiles.dat', 7, 1, worst)
        call compare_table(work//'/field.dat', &
            fixtures//'/reference_field.dat', 4, 0, worst)
        print '(a,es10.3)', 'K6-mini worst relative deviation: ', worst
        if (worst > tolerance) error stop 'K6-mini regression exceeds 1e-12'
        print '(a)', 'PASS: K6-mini regression'
    end subroutine

    subroutine compare_summary(path, reference, worst)
        ! Every key of the reference must be present; values relative.
        character(*), intent(in) :: path, reference
        real(dp), intent(inout) :: worst
        character(64) :: key, other
        character(256) :: line
        real(dp) :: expected, value, deviation
        integer :: unit, ref_unit, stat
        logical :: found
        open (newunit=ref_unit, file=reference, status='old', action='read')
        do
            read (ref_unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (len_trim(line) == 0) cycle
            read (line, *) key, expected
            found = .false.
            open (newunit=unit, file=path, status='old', action='read')
            do
                read (unit, '(a)', iostat=stat) line
                if (stat /= 0) exit
                read (line, *) other
                if (trim(other) /= trim(key)) cycle
                read (line, *) other, value
                found = .true.
                exit
            end do
            close (unit)
            if (.not. found) then
                print '(2a)', 'missing summary key: ', trim(key)
                error stop 'K6-mini summary incomplete'
            end if
            deviation = abs(value - expected)/max(abs(expected), tiny(1.0_dp))
            print '(a24,2es25.16e3,es10.3)', key, expected, value, deviation
            worst = max(worst, deviation)
        end do
        close (ref_unit)
    end subroutine

    subroutine compare_table(path, reference, columns, header, worst)
        ! Columnwise normwise deviation max|a - b|/max|b|: the fast (default)
        ! reductions leave 1e-14 run-to-run noise in the flux coefficients,
        ! which a pointwise ratio amplifies where a component crosses zero.
        character(*), intent(in) :: path, reference
        integer, intent(in) :: columns, header
        real(dp), intent(inout) :: worst
        real(dp), allocatable :: a(:, :), b(:, :)
        real(dp) :: deviation
        integer :: j
        call read_table(path, columns, header, a)
        call read_table(reference, columns, header, b)
        if (size(a, 2) /= size(b, 2)) error stop 'K6-mini table length differs'
        do j = 1, columns
            deviation = maxval(abs(a(j, :) - b(j, :)))/max(maxval(abs(b(j, :))), &
                tiny(1.0_dp))
            print '(a,a,i0,es10.3)', trim(reference(index(reference, '/', &
                back=.true.) + 1:)), ' column ', j, deviation
            worst = max(worst, deviation)
        end do
    end subroutine

    subroutine read_table(path, columns, header, table)
        character(*), intent(in) :: path
        integer, intent(in) :: columns, header
        real(dp), allocatable, intent(out) :: table(:, :)
        real(dp) :: row(columns)
        integer :: unit, stat, n, i
        open (newunit=unit, file=path, status='old', action='read')
        do i = 1, header
            read (unit, *)
        end do
        n = 0
        do
            read (unit, *, iostat=stat) row
            if (stat /= 0) exit
            n = n + 1
        end do
        rewind (unit)
        do i = 1, header
            read (unit, *)
        end do
        allocate (table(columns, n))
        do i = 1, n
            read (unit, *) table(:, i)
        end do
        close (unit)
        if (n == 0) error stop 'empty K6-mini table'
    end subroutine
end program test_equilibrium_core
