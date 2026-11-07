!> NTV field line of a Boozer `.bc` surface (kin6d_ntv_boozer_line): closed
!! forms on a synthetic file of circular surfaces with a toroidal current and
!! magnetic shear, and refusal of a file with n /= 0 rows. With arguments
!! (file s n_theta) it prints the line of a given file instead.
program test_ntv_boozer_line
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_boozer_file, only: boozer_file_t, boozer_surface_t, load_boozer_file, &
        boozer_surface, BOOZER_OK
    use kin6d_ntv_boozer_line, only: ntv_boozer_line_t, ntv_boozer_line, NTV_LINE_OK, &
        NTV_LINE_REFUSED
    implicit none
    real(dp), parameter :: r0 = 1.7_dp, a = 0.5_dp, b0 = 2.3_dp, eps = 0.1_dp
    real(dp), parameter :: jpol = -1.2e7_dp, itor_s = -3.0e5_dp, flux = -1.3_dp
    real(dp), parameter :: twopi = 8*atan(1.0_dp), mu0 = 2*twopi*1e-7_dp
    integer :: failures = 0

    if (command_argument_count() >= 3) then
        call print_line()
        stop
    end if
    call circular()
    call refuse_nonaxisymmetric()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'
contains

    subroutine write_file(name, n_row)
        character(*), intent(in) :: name
        integer, intent(in) :: n_row
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a,es24.16,a)') ' 1 0 9 1 ', flux, ' 0.5 1.7'
        do k = 1, 9
            s = 0.1_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.9_dp - 0.3_dp*s, jpol, itor_s*s, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, r0, 0.0_dp, 0.0_dp, b0
            write (unit, '(2i5,4es24.16)') 1, n_row, a*sqrt(s), a*sqrt(s), 0.0_dp, &
                -eps*b0*sqrt(s)
        end do
        close (unit)
    end subroutine write_file

    subroutine circular()
        character(*), parameter :: name = 'ntv_line_circular.bc'
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(ntv_boozer_line_t) :: line
        real(dp), parameter :: s = 0.43_dp
        real(dp) :: theta(5), stream(5), mirror(5), drift(5), drift2(5), source(5), &
            measure(5), bmod(5), iota, g, i, psi, want
        integer :: status, unit, k
        call write_file(name, 0)
        call load_boozer_file(name, bf, status)
        call boozer_surface(bf, s, surf, status)
        call check('surface', status == BOOZER_OK)
        call ntv_boozer_line(surf, line, status)
        call check('line', status == NTV_LINE_OK)
        theta = [(0.3_dp + 1.1_dp*k, k=1, 5)]
        call line%sample(theta, stream, mirror, drift, drift2, source, measure)
        iota = 0.9_dp - 0.3_dp*s
        g = mu0*jpol/twopi
        i = mu0*itor_s*s/twopi
        psi = flux/twopi
        bmod = b0*(1 - eps*sqrt(s)*cos(theta))
        do k = 1, 5
            call close_to('b.grad theta', stream(k), iota*bmod(k)/(g + iota*i), 1e-12_dp)
            call close_to('mirror', mirror(k), iota/(g + iota*i)*b0*eps*sqrt(s) &
                *sin(theta(k))/2, 1e-12_dp)
            ! dB/ds = -eps b0 cos(theta)/(2 sqrt(s))
            want = eps*b0*cos(theta(k))/(2*sqrt(s))/(iota*psi*bmod(k))
            call close_to('drift', drift(k), want, 1e-12_dp)
            want = 2*(mu0*itor_s/twopi - g*(-0.3_dp)/iota**2)/(psi*(g + iota*i))
            call close_to('drift2', drift2(k), want, 1e-12_dp)
            call close_to('source', source(k), 2*b0/(3*bmod(k)), 1e-13_dp)
            call close_to('measure', measure(k), 2*(b0/bmod(k))**3, 1e-13_dp)
        end do
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine circular

    !> theta, b, b.grad theta, mirror, drift, drift2 and the flux-Jacobian
    !> average <B^2> = 1/<1/B^2> of a surface of a file.
    subroutine print_line()
        character(512) :: file, arg
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(ntv_boozer_line_t) :: line
        real(dp) :: s
        real(dp), allocatable :: theta(:), stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:), b(:)
        integer :: n, status, k
        call get_command_argument(1, file)
        call get_command_argument(2, arg); read (arg, *) s
        call get_command_argument(3, arg); read (arg, *) n
        call load_boozer_file(trim(file), bf, status)
        call boozer_surface(bf, s, surf, status)
        if (status /= BOOZER_OK) error stop 'surface outside the file'
        call ntv_boozer_line(surf, line, status)
        if (status /= NTV_LINE_OK) error stop 'line refused'
        allocate (theta(n), stream(n), mirror(n), drift(n), drift2(n), source(n), &
            measure(n), b(n))
        theta = [(twopi*(k - 1)/n, k=1, n)]
        call line%sample(theta, stream, mirror, drift, drift2, source, measure)
        b = 2/(3*source)
        print '(a,3es24.16)', 'iota bref <B^2>_flux ', surf%iota, line%bref, &
            line%bref**2/(sum(1/b**2)/n)
        do k = 1, n
            print '(7es24.16)', theta(k), b(k), stream(k), mirror(k), drift(k), &
                drift2(k), measure(k)
        end do
    end subroutine print_line

    subroutine refuse_nonaxisymmetric()
        character(*), parameter :: name = 'ntv_line_n1.bc'
        type(boozer_file_t) :: bf
        type(boozer_surface_t) :: surf
        type(ntv_boozer_line_t) :: line
        integer :: status, unit
        call write_file(name, 1)
        call load_boozer_file(name, bf, status)
        call boozer_surface(bf, 0.43_dp, surf, status)
        call ntv_boozer_line(surf, line, status)
        call check('refuse n /= 0', status == NTV_LINE_REFUSED)
        open (newunit=unit, file=name)
        close (unit, status='delete')
    end subroutine refuse_nonaxisymmetric

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
end program test_ntv_boozer_line
