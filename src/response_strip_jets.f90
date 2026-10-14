module kin6d_response_strip_jets
    !! Imported TC24 candidate jets, in physical (R,phi,Z) at phi=0.
    !! u=sqrt(psin), theta=geometric angle about the supplied magnetic axis.
    !! Smooth spatial jets differentiate the represented spline and implicit
    !! level equation, never the root iterations. Input-data derivatives and
    !! continuous error enclosures are unavailable. No pressure is subtracted.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_eval, profiles_of_psin, &
                                     flux_ray_point
    use kin6d_mhd_ideal_n, only: equilibrium_point, ideal_options_t
    use kin6d_mhd_ampere_n, only: mu0
    use kin6d_response_fluid, only: response_fluid_equilibrium_t
    use kin6d_response_strip_jets_leaves, only: strip_field_jets, strip_map_jets
    implicit none
    private
    type, public :: strip_equilibrium_jets_t
        real(dp) :: x(2) = 0, jac(2, 2) = 0, inverse(2, 2) = 0
        !! hmap(:,1:3)=X_uu,X_u_theta,X_theta_theta.
        real(dp) :: hmap(2, 3) = 0, db(3, 3) = 0, curl_current(3) = 0
        real(dp) :: label_defect = huge(1.0_dp), ampere_defect = huge(1.0_dp)
        type(response_fluid_equilibrium_t) :: fluid
    end type
    public :: strip_equilibrium_point, strip_equilibrium_grid
contains
    subroutine strip_equilibrium_point(g, u, theta, rho, eta, jets, info)
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: u, theta, rho, eta
        type(strip_equilibrium_jets_t), intent(out) :: jets
        integer, intent(out) :: info
        type(ideal_options_t) :: opt
        real(dp) :: ray, c, s, psi, pr, pz, prr, prz, pzz
        real(dp) :: depth, f, fx, p, px, pxx, gp(2), pp, det
        logical :: ok
        info = 1
        if (any(.not. ieee_is_finite([u, theta, rho, eta]))) return
        if (u <= 0 .or. u > 1 .or. rho <= 0 .or. eta < 0) return
        if (.not. allocated(g%coef) .or. .not. allocated(g%fcoef)) return
        depth = g%psi_edge - g%psi_axis
        if (.not. ieee_is_finite(depth) .or. abs(depth) <= tiny(1.0_dp)) return
        call flux_ray_point(g, u*u, theta, ray, ok)
        if (.not. ok) return
        c = cos(theta); s = sin(theta)
        jets%x = [g%r_axis + ray*c, g%z_axis + ray*s]
        if (jets%x(1) <= 0) return
        call psi_eval(g, jets%x(1), jets%x(2), psi, pr, pz, prr, prz, pzz)
        if (any(.not. ieee_is_finite([psi, pr, pz, prr, prz, pzz]))) return
        if (abs(pr*c + pz*s) <= tiny(1.0_dp)) return
        call profiles_of_psin(g, u*u, f, fx, p, px, pxx)
        if (any(.not. ieee_is_finite([f, fx, p, px, pxx]))) return
        if (p < 0) return
        jets%fluid%p = p; jets%fluid%rho = rho; jets%fluid%eta = eta
        call equilibrium_point(g, opt, jets%x(1), jets%x(2), &
                               jets%fluid%b0, jets%fluid%j0, gp, pp)
        call strip_field_jets(jets%x(1), depth, g%s_pol, pr, pz, prr, prz, &
                              pzz, f, fx, px, pxx, 1/mu0, &
                              jets%db(1, 1), jets%db(2, 1), jets%db(3, 1), &
                              jets%db(1, 2), jets%db(2, 2), jets%db(3, 2), &
                              jets%db(1, 3), jets%db(2, 3), jets%db(3, 3), &
                              jets%fluid%gp(1), jets%fluid%gp(2), jets%fluid%gp(3), &
                        jets%fluid%hp(1, 1), jets%fluid%hp(2, 1), jets%fluid%hp(3, 1), &
                        jets%fluid%hp(1, 2), jets%fluid%hp(2, 2), jets%fluid%hp(3, 2), &
                        jets%fluid%hp(1, 3), jets%fluid%hp(2, 3), jets%fluid%hp(3, 3), &
                       jets%curl_current(1), jets%curl_current(2), jets%curl_current(3))
        call strip_map_jets(u, depth, ray, c, s, pr, pz, prr, prz, pzz, &
                       jets%jac(1, 1), jets%jac(2, 1), jets%jac(1, 2), jets%jac(2, 2), &
                            jets%inverse(1, 1), jets%inverse(2, 1), &
                            jets%inverse(1, 2), jets%inverse(2, 2), &
                            jets%hmap(1, 1), jets%hmap(2, 1), jets%hmap(1, 2), &
                            jets%hmap(2, 2), jets%hmap(1, 3), jets%hmap(2, 3))
        det = jets%jac(1, 1)*jets%jac(2, 2) - jets%jac(1, 2)*jets%jac(2, 1)
        if (.not. ieee_is_finite(det) .or. abs(det) <= tiny(1.0_dp)) return
        jets%label_defect = abs((psi - g%psi_axis)/depth - u*u)
        jets%ampere_defect = norm2(jets%curl_current - jets%fluid%j0)
        info = 0
    end subroutine

    subroutine strip_equilibrium_grid(g, radial, nt, rho, eta, jets, info)
        type(psi_grid_t), intent(in) :: g
        real(dp), intent(in) :: radial(:), rho, eta
        integer, intent(in) :: nt
        type(strip_equilibrium_jets_t), intent(out) :: jets(size(radial), nt)
        integer, intent(out) :: info
        integer :: i, j, st, failures
        failures = 0
        if (nt < 5 .or. mod(nt, 2) == 0) then
            info = 1
            return
        end if
        !$omp parallel do collapse(2) private(st) reduction(+:failures) schedule(static)
        do j = 1, nt
            do i = 1, size(radial)
                call strip_equilibrium_point(g, radial(i), &
                                   2*acos(-1.0_dp)*(j - 1)/nt, rho, eta, jets(i, j), st)
                if (st /= 0) failures = failures + 1
            end do
        end do
        !$omp end parallel do
        info = failures
    end subroutine
end module
