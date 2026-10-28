module island_outer_reference
    ! Cylindrical, zero-beta, single-helicity, stationary, fully reconnected
    ! constant-psi outer problem. See the V1 section of the NN2 notes.
    ! FortNum controls ODE error. The Frobenius cut delta is refined separately.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_island_leaves, only: cylindrical_outer_rhs
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
                           ode_integrate
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    public :: outer_response
contains
    subroutine outer_response(qaxis, qedge, delta, response, tearing, drive, status, &
                              radii, radial)
        real(dp), intent(in) :: qaxis, qedge, delta
        real(dp), intent(out) :: response, tearing, drive
        integer, intent(out) :: status
        ! Signed f(r)/f_b; response remains the separate resonant f_s/f_b.
        ! Positive radii below rmin and unresolved cut points except r_s fail.
        real(dp), intent(in), optional :: radii(:)
        real(dp), intent(out), optional :: radial(:)
        real(dp) :: rs, qc, coeff, ymin(2), ymax(2), small(2), psis, b1, b2, r
        real(dp), parameter :: rmin = 1.0e-4_dp
        integer :: i
        qc = qedge - qaxis
        response = 0; tearing = 0; drive = 0; status = 1
        if (present(radial)) radial = 0
        if (present(radii) .neqv. present(radial)) return
        if (present(radii)) then
            if (size(radial) /= size(radii)) return
            if (.not. all(ieee_is_finite(radii))) return
            if (any(radii < 0) .or. any(radii > 1)) return
        end if
        if (qaxis >= 2 .or. qedge <= 2 .or. delta <= 0) return
        rs = sqrt((2 - qaxis)/qc)
        ! Residue of the simple pole in J'/(r(1/q-1/2)).
        coeff = 2*qaxis/rs
        call integrate(rmin, rs - delta, [1.0_dp, 2/rmin], ymin)
        if (status /= 0) return
        ! psi=psi_s(1+A x log|x|)+L x+O(x^2 log|x|).
        ! Equal finite parts L_+=L_- mean zero sheet current. The omitted
        ! terms are tested by the cut ladder, not claimed rigorously bounded.
        psis = (ymin(1) + delta*ymin(2))/(1 + coeff*delta)
        ymax = [2*psis - ymin(1), ymin(2)]/psis
        call integrate(rs + delta, 1.0_dp, ymax, small)
        if (status /= 0) return
        b1 = small(1)
        ! The small Frobenius branch has psi_s=0 and unit finite derivative.
        call integrate(rs + delta, 1.0_dp, [delta, 1.0_dp], small)
        if (status /= 0) return
        b2 = small(1)
        if (min(abs(b1), abs(b2)) <= tiny(1.0_dp)) then
            status = 3
            return
        end if
        tearing = -b1/b2; drive = 1/b2; response = 1/b1
        if (present(radii)) then
            do i = 1, size(radii)
                r = radii(i)
                if (abs(r - rs) <= 32*epsilon(rs)) then
                    radial(i) = response
                else if (r == 0) then
                    radial(i) = 0
                else if (r == rmin) then
                    radial(i) = 1/(psis*b1)
                else if (r >= rmin .and. r <= rs - delta) then
                    call integrate(rmin, r, [1.0_dp, 2/rmin], small)
                    if (status /= 0) return
                    radial(i) = small(1)/(psis*b1)
                else if (r == rs + delta) then
                    radial(i) = ymax(1)/b1
                else if (r >= rs + delta) then
                    call integrate(rs + delta, r, ymax, small)
                    if (status /= 0) return
                    radial(i) = small(1)/b1
                else
                    status = 4
                    return
                end if
            end do
        end if
    contains
        subroutine integrate(a, b, initial, final)
            real(dp), intent(in) :: a, b, initial(2)
            real(dp), intent(out) :: final(2)
            type(ode_problem_t) :: problem
            type(ode_workspace_t) :: work
            type(ode_solution_t) :: sol
            type(fortnum_status_t) :: stat
            problem%rhs => rhs
            problem%t0 = a; problem%t1 = b; problem%y0 = initial
            problem%rtol = 1.0e-11_dp; problem%atol = 1.0e-13_dp
            problem%hmax = 0.005_dp
            call ode_integrate(problem, work, sol, stat)
            status = stat%code
            final = initial
            if (status == 0) final = sol%y(:, size(sol%t))
        end subroutine
        subroutine rhs(r, y, dy, ctx)
            real(dp), intent(in) :: r, y(:)
            real(dp), intent(out) :: dy(:)
            class(*), intent(in), optional :: ctx
            call cylindrical_outer_rhs(r, y(1), y(2), qaxis, qc, 2.0_dp, 1.0_dp, &
                                       dy(1), dy(2))
        end subroutine
    end subroutine
end module
