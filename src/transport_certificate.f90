module kin6d_transport_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, mag
    use kin6d_transport_laboratory, only: transport_state_t, transport_valid, &
                                          sample_transport_coefficients_interval
    use kin6d_transport_leaves, only: transport_residual_cell_interval, &
               transport_residual_box_interval, transport_certificate_radius_interval, &
                 transport_tangent_radius_interval, transport_tangent_source_interval, &
                             transport_goal_band_interval, transport_scalar_add_interval, &
                             transport_abs_value
    implicit none
    private
    public :: transport_certificate_t, certify_transport, &
              transport_certificate_accepts, certify_transport_tangent, &
              transport_observable_band, certify_transport_sample, &
              certify_transport_tangent_sample

    type :: transport_certificate_t
        real(dp) :: residual_bound = huge(0.0_dp)
        real(dp) :: state_radius = huge(0.0_dp)
        real(dp) :: tangent_residual_bound = huge(0.0_dp)
        real(dp) :: tangent_radius = huge(0.0_dp)
        real(dp) :: nu = 0.0_dp
        integer :: nodes = 0, modes = 0
        logical :: valid = .false., tangent_valid = .false., strict = .false.
    end type
contains
    subroutine certify_transport(state, certificate, ok, strict)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(out) :: certificate
        logical, intent(out) :: ok
        logical, intent(in), optional :: strict
        type(interval_t) :: radius

        certificate = transport_certificate_t()
        ok = .false.
        if (present(strict)) certificate%strict = strict
        if (.not. transport_valid(state)) return
        call residual_bound(state, state%cosine, state%sine, state%source_cos, &
                            state%source_sin, 0.0_dp, certificate%residual_bound, ok, &
                            certificate%strict)
        if (.not. ok) return
        ok = .false.
        call transport_certificate_radius_interval(interval(state%nu), &
                         interval(certificate%residual_bound), interval(0.0_dp), radius)
        if (.not. finite_interval(radius)) return
        if (radius%hi < 0.0_dp) return
        certificate%state_radius = radius%hi
        certificate%nu = state%nu
        certificate%nodes = size(state%v)
        certificate%modes = size(state%cosine, 1) - 1
        certificate%valid = .true.
        ok = .true.
    end subroutine

    logical function transport_certificate_accepts(state, radius, strict) result(ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: radius
        logical, intent(in), optional :: strict
        type(transport_certificate_t) :: certificate
        logical :: prepared

        ok = .false.
        if (.not. ieee_is_finite(radius)) return
        if (radius < 0.0_dp) return
        call certify_transport(state, certificate, prepared, strict)
        if (.not. prepared) return
        ok = radius >= certificate%state_radius
    end function

    subroutine certify_transport_tangent(state, dnu, dsource_cos, dsource_sin, &
                                         dcos, dsin, certificate, ok, strict)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: dnu, dsource_cos(0:, :), dsource_sin(0:, :)
        real(dp), intent(in) :: dcos(0:, :), dsin(0:, :)
        type(transport_certificate_t), intent(out) :: certificate
        logical, intent(out) :: ok
        logical, intent(in), optional :: strict
        type(interval_t) :: radius
        real(dp) :: dnu_magnitude

        call certify_transport(state, certificate, ok, strict)
        if (.not. ok) return
        ok = .false.
        if (.not. ieee_is_finite(dnu)) return
        if (.not. valid_coefficients(state, dsource_cos, dsource_sin)) return
        if (.not. valid_coefficients(state, dcos, dsin)) return
        call residual_bound(state, dcos, dsin, dsource_cos, dsource_sin, dnu, &
                            certificate%tangent_residual_bound, ok, certificate%strict)
        if (.not. ok) return
        ok = .false.
        call transport_abs_value(dnu, dnu_magnitude)
        call transport_tangent_radius_interval(interval(state%nu), &
                     interval(certificate%tangent_residual_bound), interval(dnu_magnitude), &
                           interval(certificate%state_radius), interval(0.0_dp), radius)
        if (.not. finite_interval(radius)) return
        if (radius%hi < 0.0_dp) return
        certificate%tangent_radius = radius%hi
        certificate%tangent_valid = .true.
        ok = .true.
    end subroutine

    logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value

        ok = .false.
        if (.not. all(ieee_is_finite([value%lo, value%hi]))) return
        if (value%lo > value%hi) return
        ok = .true.
    end function

    logical function valid_coefficients(state, cosine, sine) result(ok)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: cosine(0:, :), sine(0:, :)

        ok = .false.
        if (any(shape(cosine) /= shape(state%cosine))) return
        if (any(shape(sine) /= shape(state%sine))) return
        if (.not. all(ieee_is_finite(cosine))) return
        if (.not. all(ieee_is_finite(sine))) return
        if (any(abs(sine(0, :)) > 0.0_dp)) return
        ok = .true.
    end function

    subroutine residual_bound(state, cosine, sine, source_cos, source_sin, dnu, &
                              bound, ok, strict)
        type(transport_state_t), intent(in) :: state
        real(dp), intent(in) :: cosine(0:, :), sine(0:, :)
        real(dp), intent(in) :: source_cos(0:, :), source_sin(0:, :), dnu
        real(dp), intent(out) :: bound
        logical, intent(out) :: ok
        logical, intent(in) :: strict
        type(interval_t) :: gc0, gc1, gs0, gs1, bc0, bc1, bc2, bs0, bs1, bs2
        type(interval_t) :: cell, next, rc, rs, vr, cr, sr, gr, hr
        real(dp) :: cosine_bound, sine_bound
        integer :: j, mode

        bound = huge(0.0_dp)
        ok = .false.
        bound = 0.0_dp
        do j = 1, size(state%v) - 1
            cell = interval(0.0_dp)
            do mode = 0, size(cosine, 1) - 1
                gc0 = interval(source_cos(mode, j))
                gc1 = interval(source_cos(mode, j + 1))
                gs0 = interval(source_sin(mode, j))
                gs1 = interval(source_sin(mode, j + 1))
                if (abs(dnu) > 0.0_dp) then
                    call transport_tangent_source_interval(gc0, &
                                   interval(state%cosine(mode, j)), interval(dnu), next)
                    gc0 = next
                    call transport_tangent_source_interval(gc1, &
                               interval(state%cosine(mode, j + 1)), interval(dnu), next)
                    gc1 = next
                    call transport_tangent_source_interval(gs0, &
                                     interval(state%sine(mode, j)), interval(dnu), next)
                    gs0 = next
                    call transport_tangent_source_interval(gs1, &
                                 interval(state%sine(mode, j + 1)), interval(dnu), next)
                    gs1 = next
                end if
                if (.not. finite_interval(gc0)) return
                if (.not. finite_interval(gc1)) return
                if (.not. finite_interval(gs0)) return
                if (.not. finite_interval(gs1)) return
                if (strict) then
                    vr = interval(state%v(j), state%v(j + 1))
                    cr = interval(min(cosine(mode, j), cosine(mode, j + 1)), &
                                  max(cosine(mode, j), cosine(mode, j + 1)))
                    sr = interval(min(sine(mode, j), sine(mode, j + 1)), &
                                  max(sine(mode, j), sine(mode, j + 1)))
                    gr = interval(min(gc0%lo, gc1%lo), max(gc0%hi, gc1%hi))
                    hr = interval(min(gs0%lo, gs1%lo), max(gs0%hi, gs1%hi))
                    call transport_residual_box_interval(vr, interval(state%nu), &
                                       interval(real(mode, dp)), cr, sr, gr, hr, rc, rs)
                    if (.not. finite_interval(rc)) return
                    if (.not. finite_interval(rs)) return
                    cosine_bound = mag(rc)
                    sine_bound = mag(rs)
                else
                    call transport_residual_cell_interval(interval(state%v(j)), &
                                         interval(state%v(j + 1)), interval(state%nu), &
                                  interval(real(mode, dp)), interval(cosine(mode, j)), &
                               interval(cosine(mode, j + 1)), interval(sine(mode, j)), &
                                      interval(sine(mode, j + 1)), gc0, gc1, gs0, gs1, &
                                                          bc0, bc1, bc2, bs0, bs1, bs2)
                    if (.not. finite_interval(bc0)) return
                    if (.not. finite_interval(bc1)) return
                    if (.not. finite_interval(bc2)) return
                    if (.not. finite_interval(bs0)) return
                    if (.not. finite_interval(bs1)) return
                    if (.not. finite_interval(bs2)) return
                    cosine_bound = max(mag(bc0), mag(bc1), mag(bc2))
                    sine_bound = max(mag(bs0), mag(bs1), mag(bs2))
                end if
                call transport_scalar_add_interval(cell, interval(cosine_bound), next)
                cell = next
                if (mode > 0) then
                    call transport_scalar_add_interval(cell, interval(sine_bound), next)
                    cell = next
                end if
                if (.not. finite_interval(cell)) return
            end do
            bound = max(bound, cell%hi)
        end do
        ok = .true.
    end subroutine

    subroutine transport_observable_band(candidate, weight_norm, radius, &
                                         lower, upper, ok)
        type(interval_t), intent(in) :: candidate, weight_norm
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: lo, hi

        lower = -huge(0.0_dp); upper = huge(0.0_dp)
        ok = .false.
        if (.not. finite_interval(candidate)) return
        if (.not. finite_interval(weight_norm)) return
        if (weight_norm%lo < 0.0_dp) return
        if (.not. ieee_is_finite(radius)) return
        if (radius < 0.0_dp) return
     call transport_goal_band_interval(candidate, weight_norm, interval(radius), lo, hi)
        if (.not. finite_interval(lo)) return
        if (.not. finite_interval(hi)) return
        lower = lo%lo; upper = hi%hi
        ok = .true.
    end subroutine

    logical function prepared_layout(state, certificate) result(ok)
        ! The caller keeps the prepared source/candidate packet immutable.
        ! These O(1) checks do not revalidate its entire scientific content.
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate

        ok = .false.
        if (.not. certificate%valid) return
        if (.not. allocated(state%v)) return
        if (.not. allocated(state%cosine)) return
        if (.not. allocated(state%sine)) return
        if (size(state%v) /= certificate%nodes) return
        if (size(state%cosine, 1) /= certificate%modes + 1) return
        if (size(state%cosine, 2) /= certificate%nodes) return
        if (any(shape(state%sine) /= shape(state%cosine))) return
        if (.not. ieee_is_finite(state%nu)) return
        if (state%nu < certificate%nu .or. state%nu > certificate%nu) return
        ok = .true.
    end function

  subroutine certify_transport_sample(state, certificate, x, velocity, lower, upper, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate
        type(interval_t), intent(in) :: x, velocity
        real(dp), intent(out) :: lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: candidate

        lower = -huge(0.0_dp); upper = huge(0.0_dp)
        ok = .false.
        if (.not. prepared_layout(state, certificate)) return
        call sample_transport_coefficients_interval(state, state%cosine, state%sine, &
                                                    x, velocity, candidate, ok)
        if (.not. ok) return
        call transport_observable_band(candidate, interval(1.0_dp), &
                                       certificate%state_radius, lower, upper, ok)
    end subroutine

    subroutine certify_transport_tangent_sample(state, certificate, dcos, dsin, &
                                                x, velocity, lower, upper, ok)
        type(transport_state_t), intent(in) :: state
        type(transport_certificate_t), intent(in) :: certificate
        real(dp), intent(in) :: dcos(0:, :), dsin(0:, :)
        type(interval_t), intent(in) :: x, velocity
        real(dp), intent(out) :: lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: candidate

        lower = -huge(0.0_dp); upper = huge(0.0_dp)
        ok = .false.
        if (.not. prepared_layout(state, certificate)) return
        if (.not. certificate%tangent_valid) return
        call sample_transport_coefficients_interval(state, dcos, dsin, x, velocity, &
                                                    candidate, ok)
        if (.not. ok) return
        call transport_observable_band(candidate, interval(1.0_dp), &
                                       certificate%tangent_radius, lower, upper, ok)
    end subroutine
end module
