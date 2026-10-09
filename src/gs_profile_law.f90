module kin6d_gs_profile_law
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use kin6d_grad_shafranov_leaves, only: gs_cubic, gs_cubic_jvp, &
                          gs_cubic_primitive, gs_normalized_flux, gs_normalized_flux_jvp
    implicit none
    private
    public :: gs_profile_law_t, read_gs_profile_law, write_gs_profile_law
    public :: gs_law_sources, gs_law_primitives, gs_law_coordinate
    public :: gs_law_kink
    type :: gs_profile_law_t
        ! Coefficients in ascending powers of (coordinate - knots(i)).
        ! Sources are physical dp/dpsi and d(F^2/2)/dpsi, also for normalized s.
        real(dp), allocatable :: knots(:), pp(:, :), ff(:, :)
        real(dp) :: p_edge = 0, f_edge = 1, axis = -1
        logical :: normalized = .false.
    contains
        procedure :: valid => law_valid
    end type
contains
    logical function gs_law_kink(law, psi) result(kink)
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi
        real(dp) :: x, dx, da, value, left_slope, right_slope
        integer :: i, field
        call gs_law_coordinate(law, psi, x, dx, da)
        kink = .false.
        do i = 2, size(law%knots) - 1
            if (abs(x - law%knots(i)) > 0) cycle
            do field = 1, 2
                if (field == 1) then
                    call polynomial(law%pp(:, i - 1), x - law%knots(i - 1), &
                                    value, left_slope)
                    right_slope = law%pp(2, i)
                else
                    call polynomial(law%ff(:, i - 1), x - law%knots(i - 1), &
                                    value, left_slope)
                    right_slope = law%ff(2, i)
                end if
                if (abs(left_slope - right_slope) > 64*epsilon(1.0_dp) &
                    *max(1.0_dp, abs(left_slope), abs(right_slope))) kink = .true.
            end do
        end do
    end function

    logical function law_valid(law) result(valid)
        class(gs_profile_law_t), intent(in) :: law
        real(dp) :: left, right, h
        integer :: i, k, n
        valid = .false.
        if (.not. allocated(law%knots)) return
        if (.not. allocated(law%pp) .or. .not. allocated(law%ff)) return
        n = size(law%knots)
        if (n < 2) return
        if (any(shape(law%pp) /= [4, n - 1])) return
        if (any(shape(law%ff) /= [4, n - 1])) return
        if (.not. all(ieee_is_finite(law%knots))) return
        if (.not. all(ieee_is_finite(law%pp))) return
        if (.not. all(ieee_is_finite(law%ff))) return
        if (.not. all(ieee_is_finite([law%p_edge, law%f_edge, law%axis]))) return
        if (any(law%knots(2:) <= law%knots(:n - 1))) return
        if (abs(law%f_edge) <= tiny(1.0_dp) .or. law%p_edge < 0) return
        if (law%normalized) then
            if (abs(law%axis) <= tiny(1.0_dp)) return
            if (abs(law%knots(1)) > 0 .or. abs(law%knots(n) - 1) > 0) return
        end if
        ! Moving source-cut terms cancel only for continuous laws.
        do i = 1, n - 2
            h = law%knots(i + 1) - law%knots(i)
            do k = 1, 2
                if (k == 1) then
                    call polynomial(law%pp(:, i), h, left)
                    right = law%pp(1, i + 1)
                else
                    call polynomial(law%ff(:, i), h, left)
                    right = law%ff(1, i + 1)
                end if
                if (abs(left - right) > 128*epsilon(1.0_dp)*max(1.0_dp, &
                                                          abs(left), abs(right))) return
            end do
        end do
        valid = .true.
    end function

    subroutine read_gs_profile_law(path, law, status)
        character(*), intent(in) :: path
        type(gs_profile_law_t), intent(out) :: law
        integer, intent(out) :: status
        integer :: unit, n, mode, i, ios
        status = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read (unit, *, iostat=ios) n, mode
        if (ios == 0) then
            if (n < 2 .or. n > 100000 .or. mode < 0 .or. mode > 1) ios = 1
        end if
        if (ios == 0) then
            law%normalized = mode == 1
            allocate (law%knots(n), law%pp(4, n - 1), law%ff(4, n - 1))
            read (unit, *, iostat=ios) law%p_edge, law%f_edge, law%axis
        end if
        if (ios == 0) read (unit, *, iostat=ios) law%knots
        if (ios == 0) then
            do i = 1, n - 1
                read (unit, *, iostat=ios) law%pp(:, i), law%ff(:, i)
                if (ios /= 0) exit
            end do
        end if
        close (unit)
        if (ios /= 0) return
        if (law%valid()) status = 0
    end subroutine

    subroutine write_gs_profile_law(path, law)
        character(*), intent(in) :: path
        type(gs_profile_law_t), intent(in) :: law
        integer :: unit, i
        open (newunit=unit, file=path, status='replace')
        write (unit, *) size(law%knots), merge(1, 0, law%normalized)
        write (unit, '(3es25.16e3)') law%p_edge, law%f_edge, law%axis
        write (unit, '(*(es25.16e3,1x))') law%knots
        do i = 1, size(law%knots) - 1
            write (unit, '(8es25.16e3)') law%pp(:, i), law%ff(:, i)
        end do
        close (unit)
    end subroutine

    subroutine gs_law_coordinate(law, psi, x, dx, da)
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: x, dx, da
        x = psi; dx = 1; da = 0
        if (.not. law%normalized) return
        call gs_normalized_flux(psi, law%axis, x)
        call gs_normalized_flux_jvp(psi, law%axis, 1.0_dp, 0.0_dp, dx)
        call gs_normalized_flux_jvp(psi, law%axis, 0.0_dp, 1.0_dp, da)
    end subroutine

    subroutine gs_law_sources(law, psi, pp, ff, dpp, dff, app, aff)
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: pp, ff
        real(dp), intent(out), optional :: dpp, dff, app, aff
        real(dp) :: x, dx, da, p, f
        integer :: i
        call gs_law_coordinate(law, psi, x, dx, da)
        i = interval(law%knots, x)
        x = x - law%knots(i)
        if (present(dpp) .or. present(dff) .or. present(app) .or. present(aff)) then
            call polynomial(law%pp(:, i), x, pp, p)
            call polynomial(law%ff(:, i), x, ff, f)
        else
            call polynomial(law%pp(:, i), x, pp)
            call polynomial(law%ff(:, i), x, ff)
        end if
        if (present(dpp)) dpp = p*dx
        if (present(dff)) dff = f*dx
        if (present(app)) app = p*da
        if (present(aff)) aff = f*da
    end subroutine

    subroutine gs_law_primitives(law, psi, pressure, f, status)
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: pressure, f
        integer, intent(out) :: status
        real(dp) :: x, dx, da, edge, scale, squared
        call gs_law_coordinate(law, psi, x, dx, da)
        edge = 0; scale = 1
        if (law%normalized) then
            edge = 1; scale = -law%axis
        end if
        pressure = law%p_edge + scale*(primitive(law, law%pp, x) &
                                       - primitive(law, law%pp, edge))
        squared = law%f_edge**2 + 2*scale*(primitive(law, law%ff, x) &
                                           - primitive(law, law%ff, edge))
        status = 1
        f = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. all(ieee_is_finite([pressure, squared]))) return
        if (pressure < -1.0e-10_dp*max(1.0_dp, law%p_edge)) return
        if (squared <= 0) return
        f = sign(sqrt(squared), law%f_edge)
        status = 0
    end subroutine

    integer function interval(knots, x) result(i)
        real(dp), intent(in) :: knots(:), x
        integer :: lo, hi, mid
        lo = 1; hi = size(knots)
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (x < knots(mid)) then
                hi = mid
            else
                lo = mid
            end if
        end do
        i = lo
    end function

    real(dp) function primitive(law, coefficients, x) result(value)
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(in) :: coefficients(:, :), x
        real(dp) :: part, width
        integer :: i, last
        last = interval(law%knots, x)
        value = 0
        do i = 1, last
            width = x - law%knots(i)
            if (i < last) width = law%knots(i + 1) - law%knots(i)
            call gs_cubic_primitive(width, coefficients(1, i), coefficients(2, i), &
                                    coefficients(3, i), coefficients(4, i), part)
            value = value + part
        end do
    end function

    subroutine polynomial(c, x, value, slope)
        real(dp), intent(in) :: c(4), x
        real(dp), intent(out) :: value
        real(dp), intent(out), optional :: slope
        call gs_cubic(x, c(1), c(2), c(3), c(4), value)
        if (present(slope)) call gs_cubic_jvp(x, c(1), c(2), c(3), c(4), &
                                          1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, slope)
    end subroutine
end module
