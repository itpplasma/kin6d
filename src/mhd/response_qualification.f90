module kin6d_response_qualification
    !! Frozen TC24 RULE.md gates. These are empirical qualification tests,
    !! never deterministic certificates or fitted ideal-limit estimates.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    integer, parameter, public :: qualification_axes = 8
    character(12), parameter, public :: qualification_axis_names(8) = &
        [character(12) :: 'radial', 'core', 'grading', 'harmonic', &
        'poloidal', 'extent', 'resistivity', 'frequency']
    type, public :: response_record_t
        complex(dp) :: b = (0.0_dp, 0.0_dp), b_vac = (0.0_dp, 0.0_dp)
        complex(dp) :: total_flux = (0.0_dp, 0.0_dp), current = (0.0_dp, 0.0_dp)
        real(dp) :: residual = huge(1.0_dp), balance = huge(1.0_dp)
        real(dp) :: extent(2) = 0.0_dp
        logical :: available = .false., coupled = .false., calibrated = .false.
        logical :: full_pressure = .false., all_rationals = .false.
        logical :: interface_test = .false., whole_edge = .false.
    end type response_record_t
    type, public :: qualification_t
        logical :: gates(8) = .false., passed = .false.
        real(dp) :: changes(2, 8) = huge(1.0_dp)
        real(dp) :: gpec_ratio = -1.0_dp, gpec_phase = 0.0_dp
        real(dp) :: mars_ratio = -1.0_dp
        integer :: verdict = 0 !! 0 no verdict, 1 GPEC, 2 MARS
    end type qualification_t
    public :: qualify_response, record_admissible, ladder_gate
contains
    logical function record_admissible(row) result(ok)
        type(response_record_t), intent(in) :: row
        ok = .false.
        if (.not. row%available) return
        if (.not. row%coupled .or. .not. row%calibrated) return
        if (.not. row%full_pressure .or. .not. row%all_rationals) return
        if (.not. all(ieee_is_finite([real(row%b), aimag(row%b), &
            real(row%b_vac), aimag(row%b_vac), real(row%total_flux), &
            aimag(row%total_flux), real(row%current), aimag(row%current), &
            row%residual, row%balance, row%extent]))) return
        if (abs(row%b_vac) <= tiny(1.0_dp)) return
        if (row%residual < 0.0_dp .or. row%balance < 0.0_dp) return
        if (row%residual > 1.0e-8_dp .or. row%balance > 1.0e-6_dp) return
        if (any(row%extent <= 0.0_dp)) return
        ok = .true.
    end function record_admissible

    subroutine ladder_gate(rows, axis, ok, changes)
        type(response_record_t), intent(in) :: rows(3)
        integer, intent(in) :: axis
        logical, intent(out) :: ok
        real(dp), intent(out) :: changes(2)
        real(dp) :: floor, spread, denom
        integer :: i, j
        ok = .false.; changes = huge(1.0_dp)
        if (axis < 1 .or. axis > 8) return
        if (.not. all([(record_admissible(rows(i)), i=1, 3)])) return
        ! The declared drive is unchanged across the ladder. The finer-row
        ! vacuum floor is retained in each successive relative difference.
        do i = 1, 2
            floor = 1.0e-3_dp*abs(rows(i + 1)%b_vac)
            changes(i) = abs(rows(i + 1)%b - rows(i)%b)/ &
                max(abs(rows(i + 1)%b), floor)
        end do
        if (axis == 6) then
            spread = 0.0_dp
            do i = 1, 2
                do j = i + 1, 3
                    denom = max(abs(rows(j)%b), 1.0e-3_dp*abs(rows(j)%b_vac))
                    spread = max(spread, abs(rows(j)%b - rows(i)%b)/denom)
                end do
            end do
            changes(2) = spread
            if (spread > 0.05_dp) return
            ! Capped identical extents cannot count as an extent experiment.
            ! Each repeated pair needs an independently passed interface test
            ! or a whole edge-strip solve on both members of that pair.
            do i = 1, 2
                do j = i + 1, 3
                    if (all(abs(rows(i)%extent - rows(j)%extent) <= &
                        64*epsilon(denom)*max(rows(i)%extent, rows(j)%extent))) then
                        if (.not. (rows(i)%interface_test .and. &
                            rows(j)%interface_test) .and. &
                            .not. (rows(i)%whole_edge .and. rows(j)%whole_edge)) return
                    end if
                end do
            end do
            ok = .true.
        else if (axis >= 7) then
            ok = changes(2) <= 0.01_dp .and. changes(2) <= changes(1)
        else
            ok = changes(2) <= 0.05_dp .and. &
                (changes(2) <= changes(1) .or. all(changes <= 0.005_dp))
        end if
    end subroutine ladder_gate

    subroutine qualify_response(rows, gpec, mars_amplitude, mars_admitted, result)
        type(response_record_t), intent(in) :: rows(3, 8)
        complex(dp), intent(in) :: gpec
        real(dp), intent(in) :: mars_amplitude
        logical, intent(in) :: mars_admitted
        type(qualification_t), intent(out) :: result
        logical :: gpec_ok, mars_ok
        complex(dp) :: b
        integer :: axis
        result = qualification_t()
        do axis = 1, 8
            call ladder_gate(rows(:, axis), axis, result%gates(axis), &
                result%changes(:, axis))
        end do
        b = rows(3, 1)%b
        if (ieee_is_finite(real(gpec)) .and. ieee_is_finite(aimag(gpec))) then
            if (abs(gpec) > tiny(1.0_dp) .and. rows(3, 1)%available) then
                if (ieee_is_finite(real(b)) .and. ieee_is_finite(aimag(b))) then
                    result%gpec_ratio = abs(b)/abs(gpec)
                    result%gpec_phase = atan2(aimag(b*conjg(gpec)), real(b*conjg(gpec)))
                end if
            end if
        end if
        result%passed = all(result%gates) .and. result%gpec_ratio >= 0.0_dp
        if (.not. result%passed) return
        ! Until the jump-based MARS normalization is independently admitted,
        ! code comparisons remain diagnostic and no exclusive verdict exists.
        if (.not. mars_admitted) return
        if (.not. ieee_is_finite(mars_amplitude)) return
        if (mars_amplitude <= tiny(1.0_dp)) return
        result%mars_ratio = abs(b)/mars_amplitude
        gpec_ok = abs(result%gpec_ratio - 1) <= 0.10_dp
        mars_ok = abs(result%mars_ratio - 1) <= 0.10_dp
        if (gpec_ok .and. .not. mars_ok) result%verdict = 1
        if (mars_ok .and. .not. gpec_ok) result%verdict = 2
    end subroutine qualify_response
end module kin6d_response_qualification
