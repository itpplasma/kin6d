!> C0: two-body Coulomb encounter in the relative coordinate r = x1 - x2.
!> The reduced problem is one particle of the reduced mass mu and charge q1
!> moving in the field of a fixed charge q2 (kin6d_nbody, isolated, one
!> source); relative acceleration (g/mu) r/|r|^3 with g from kin6d_species.
!> The candidate is the shared KDK step; the Picard flow tube below encloses
!> the exact flow independently of it.
module kin6d_interaction_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, ipoint, interval, iadd, ipowi, idiv
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_species, only: units_t, species_t, pair_scales_t, pair_scales
    use kin6d_nbody, only: nbody_system_t, nbody_isolated, nbody_add_source, &
                           nbody_step
    use kin6d_coulomb_leaves, only: coulomb_pair_interval
    use kin6d_interaction_leaves
    implicit none
    private
    public :: interaction_state_t, init_interaction, step_interaction, &
              sample_interaction, neutral_contact_trajectory, interaction_coupling
    type :: interaction_state_t
        private
        real(dp) :: y(6) = 0
        !> Relative acceleration coupling g/mu as a rigorous enclosure.
        type(interval_t) :: coupling
        type(nbody_system_t) :: relative
        type(interval_t) :: rb(3), ub(3)
        logical :: initialized = .false., certified = .false.
    end type
contains
    !> Pair (first, second) at relative position r and velocity u. This named
    !> model admits only repulsive pairs (g > 0).
    subroutine init_interaction(state, r, u, units, first, second, ok)
        type(interaction_state_t), intent(out) :: state
        real(dp), intent(in) :: r(3), u(3)
        type(units_t), intent(in) :: units
        type(species_t), intent(in) :: first, second
        logical, intent(out) :: ok
        type(pair_scales_t) :: scales
        ok = .false.
        if (.not. all(ieee_is_finite(r))) return
        if (.not. all(ieee_is_finite(u))) return
        if (maxval(abs(r)) <= 0) return
        if (.not. ieee_is_finite(sum(r*r))) return
        scales = pair_scales(units, first, second)
        if (.not. (scales%coupling > 0) .or. .not. (scales%reduced_mass > 0)) return
        call nbody_isolated(state%relative, units, &
                            [species_t('relative', first%charge, &
                                       scales%reduced_mass)], ok)
        if (.not. ok) return
        call nbody_add_source(state%relative, [0.0_dp, 0.0_dp, 0.0_dp], &
                              second%charge, ok)
        if (.not. ok) return
        state%coupling = idiv(ipoint(scales%coupling), ipoint(scales%reduced_mass))
        state%y = [r, u]
        state%rb = ipoint(r); state%ub = ipoint(u)
        state%initialized = .true.; state%certified = .true.; ok = .true.
    end subroutine

    !> Relative acceleration coupling g/mu (midpoint of its enclosure).
    pure real(dp) function interaction_coupling(state) result(c)
        type(interaction_state_t), intent(in) :: state
        c = (state%coupling%lo + state%coupling%hi)/2
    end function

    subroutine step_interaction(state, h, ok, certify)
        type(interaction_state_t), intent(inout) :: state
        real(dp), intent(in) :: h
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        real(dp) :: y(6)
        type(fortnum_status_t) :: status
        logical :: checked
        ok = .false.
        if (.not. state%initialized) return
        if (.not. ieee_is_finite(h)) return
        if (h <= 0) return
        y = state%y
        call nbody_step(state%relative, y, h, status)
        if (status%code /= FORTNUM_OK) return
        if (.not. all(ieee_is_finite(y))) return
        checked = .true.
        if (present(certify)) checked = certify
        if (.not. checked) state%certified = .false.
        if (state%certified) call validated_step(state, h)
        state%y = y; ok = .true.
    end subroutine

    subroutine validated_step(state, h)
        type(interaction_state_t), intent(inout) :: state
        real(dp), intent(in) :: h
        type(interval_t) :: yr(3), yu(3), nr(3), nu(3), a(3), a0(3), j(3), pot
        type(interval_t) :: ih, span, norm2
        real(dp) :: margin
        integer :: k, it
        logical :: included
        ih = ipoint(h); span = interval(0.0_dp, h)
        yr = state%rb; yu = state%ub
        call acceleration(yr, state%coupling, a, pot)
        do k = 1, 3
            call interaction_affine_interval(state%rb(k), state%ub(k), span, nr(k))
            call interaction_affine_interval(state%ub(k), a(k), span, nu(k))
            yr(k) = hull(yr(k), nr(k)); yu(k) = hull(yu(k), nu(k))
        end do
        included = .false.
        do it = 1, 24
            do k = 1, 3
                margin = max(1.e-13_dp, 0.05_dp*(yr(k)%hi - yr(k)%lo))
                yr(k) = iadd(yr(k), interval(-margin, margin))
                margin = max(1.e-13_dp, 0.05_dp*(yu(k)%hi - yu(k)%lo))
                yu(k) = iadd(yu(k), interval(-margin, margin))
            end do
            norm2 = iadd(iadd(ipowi(yr(1), 2), ipowi(yr(2), 2)), ipowi(yr(3), 2))
            if (norm2%lo <= 0) exit
            if (.not. valid_box(yr) .or. .not. valid_box(yu)) exit
            call acceleration(yr, state%coupling, a, pot)
            do k = 1, 3
                call interaction_affine_interval(state%rb(k), yu(k), span, nr(k))
                call interaction_affine_interval(state%ub(k), a(k), span, nu(k))
            end do
            if (.not. valid_box(nr) .or. .not. valid_box(nu)) exit
            included = all(nr%lo > yr%lo) .and. all(nr%hi < yr%hi) .and. &
                       all(nu%lo > yu%lo) .and. all(nu%hi < yu%hi)
            if (included) exit
            do k = 1, 3
                yr(k) = hull(yr(k), nr(k)); yu(k) = hull(yu(k), nu(k))
            end do
        end do
        if (.not. included) then
            state%certified = .false.; return
        end if
        ! Picard self inclusion on a nonsingular box gives a complete flow tube.
        ! Taylor's integral remainder uses acceleration and force-rate over it.
        call acceleration(state%rb, state%coupling, a0, pot)
        call interaction_force_rate_interval(yr(1), yr(2), yr(3), &
                                             state%coupling, yu(1), yu(2), yu(3), &
                                             j(1), j(2), j(3))
        do k = 1, 3
            call interaction_taylor_position_interval(state%rb(k), state%ub(k), &
                                                      a(k), ih, nr(k))
            call interaction_taylor_velocity_interval(state%ub(k), a0(k), &
                                                      j(k), ih, nu(k))
        end do
        if (.not. valid_box(nr) .or. .not. valid_box(nu)) then
            state%certified = .false.; return
        end if
        state%rb = nr; state%ub = nu
    end subroutine

    !> Relative acceleration (g/mu) r/|r|^3 enclosed over the box r.
    subroutine acceleration(r, coupling, a, potential)
        type(interval_t), intent(in) :: r(3), coupling
        type(interval_t), intent(out) :: a(3), potential
        call coulomb_pair_interval(r(1), r(2), r(3), coupling, &
                                   a(1), a(2), a(3), potential)
    end subroutine

    subroutine sample_interaction(state, r, u, rb, ub, ok, bounded)
        type(interaction_state_t), intent(in) :: state
        real(dp), intent(out) :: r(3), u(3)
        type(interval_t), intent(out) :: rb(3), ub(3)
        logical, intent(out) :: ok, bounded
        ok = state%initialized; bounded = state%certified
        r = state%y(:3); u = state%y(4:); rb = state%rb; ub = state%ub
        ! Bounds are meaningful only when bounded is true. Failed continuation
        ! remains a candidate trajectory, never a certified output.
    end subroutine

    subroutine neutral_contact_trajectory(r0, u0, radius, t, r, u, ok)
        real(dp), intent(in) :: r0(3), u0(3), radius, t
        real(dp), intent(out) :: r(3), u(3)
        logical, intent(out) :: ok
        real(dp) :: hit, discriminant, gap2, speed2, radial
        ok = .false.; r = 0; u = 0
        if (.not. all(ieee_is_finite(r0))) return
        if (.not. all(ieee_is_finite(u0))) return
        if (.not. ieee_is_finite(radius) .or. .not. ieee_is_finite(t)) return
        if (radius <= 0 .or. t < 0) return
        call interaction_contact_geometry(r0(1), r0(2), r0(3), u0(1), u0(2), &
                                          u0(3), &
                                          radius, gap2, speed2, radial, &
                                          discriminant)
        if (.not. ieee_is_finite(discriminant)) return
        if (gap2 <= 0 .or. speed2 <= 0) return
        if (discriminant <= 0 .or. radial >= 0) then
            call interaction_affine(r0, u0, t, r); u = u0
        else
            call interaction_contact_time(r0(1), r0(2), r0(3), u0(1), u0(2), u0(3), &
                                          radius, hit)
            if (t <= hit) then
                call interaction_affine(r0, u0, t, r); u = u0
            else
                call interaction_contact_after(r0(1), r0(2), r0(3), u0(1), &
                                               u0(2), u0(3), &
                                               radius, t, r(1), r(2), r(3), u(1), &
                                               u(2), u(3))
            end if
        end if
        ok = all(ieee_is_finite(r)) .and. all(ieee_is_finite(u))
    end subroutine

    pure function hull(a, b) result(c)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: c
        c = interval(min(a%lo, b%lo), max(a%hi, b%hi))
    end function
    pure logical function valid_box(a)
        type(interval_t), intent(in) :: a(:)
        valid_box = all(ieee_is_finite(a%lo)) .and. &
                    all(ieee_is_finite(a%hi)) .and. all(a%lo <= a%hi)
    end function
end module kin6d_interaction_laboratory
