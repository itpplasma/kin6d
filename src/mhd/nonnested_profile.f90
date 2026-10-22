module kin6d_nonnested_profile
    ! Closure B: fit once on a copied seed, then solve A+D with immutable
    ! spatial data. The weak source is S = -div(q_seed), not a nodal target
    ! penalty. No topology classification or fitting occurs inside Newton.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse, only: csc_t, csc_matvec
    use fortfem_meridian_fourier_curl, only: meridian_fourier_gradient
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_source_t
    use fortfem_meridian_feec_n, only: meridian_point_basis
    use kin6d_mhd_mesh, only: flux_mesh_locate
    use kin6d_equilibrium_nonnested, only: free_b_field_t
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_VALUES, &
                                       EQ_BR, EQ_BPHI, EQ_BZ, EQ_P
    use kin6d_nonnested_leaves, only: anisotropic_flux, coupled_current, &
                                      fitted_parallel_drive, profile_current_diagnostic
    implicit none
    private
    public :: nn_profile_t, nn_fitted_source_t, nn_fit_source, nn_profile_defect_t
    public :: nn_state_profile_t

    type, abstract :: nn_profile_t
        ! Reference-label/profile composition belongs to the caller. Values
        ! and gradients are physical, with vector order (R,phi,Z), lambda =
        ! J.B/|B|^2. A nonzero seed potential gradient must be supplied when
        ! fitting an existing A+D state. Implementations must own their data
        ! (no mutable pointers) and permit concurrent read-only evaluation.
    contains
        procedure(profile_i), deferred :: evaluate
    end type
    abstract interface
        subroutine profile_i(self, x, p, gradp, lambda, gradphi, status)
            import :: nn_profile_t, dp
            class(nn_profile_t), intent(in) :: self
            real(dp), intent(in) :: x(:, :) ! (R,Z,phi)
            real(dp), intent(out) :: p(:), gradp(:, :), lambda(:), gradphi(:, :)
            integer, intent(out) :: status
        end subroutine
    end interface

    type, extends(nn_profile_t) :: nn_state_profile_t
        ! Adapter for seed/fit/freeze/re-solve of an existing uniform-
        ! conductivity A+D state. All fields and potential coefficients are
        ! owned, with the V0 real-Fourier layout used by nc_solve.
        private
        type(free_b_field_t) :: field
        real(dp), allocatable :: pressure_gradient(:), potential_gradient(:)
        real(dp) :: voltage = 0, g_par = 1, g_perp = 1
    contains
        procedure :: init => state_profile_init
        procedure :: evaluate => state_profile_evaluate
    end type

    type :: nn_profile_defect_t
        real(dp) :: pressure_rms = 0, lambda_rms = 0, parallel_pressure_rms = 0
        real(dp) :: seed_current_fit_rms = 0
        character(32) :: claim = 'candidate'
        character(32) :: derivative_status = 'frozen-data-equation-level'
        character(32) :: bound_status = 'unavailable'
    end type

    type, extends(meridian_fourier_source_t) :: nn_fitted_source_t
        private
        class(equilibrium_field_t), allocatable :: seed
        class(nn_profile_t), allocatable :: profile
        real(dp) :: chi_par = 1, chi_perp = 1, g_par = 1, g_perp = 1, voltage = 0
    contains
        procedure :: source => fitted_source
        procedure :: data => fitted_data
        procedure :: ready => fitted_ready
        procedure :: defects => fitted_defects
    end type
contains
    subroutine state_profile_init(self, field, potential, voltage, g_par, &
                                  g_perp, status)
        class(nn_state_profile_t), intent(out) :: self
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: potential(:), voltage, g_par, g_perp
        integer, intent(out) :: status
        type(csc_t) :: gradient
        integer :: modes
        status = 1
        if (.not. allocated(field%p)) return
        modes = ubound(field%b, 2)
        if (size(potential) /= field%fs%n0*(2*modes + 1)) return
        if (size(field%p) /= size(potential)) return
        if (.not. all(ieee_is_finite(potential))) return
        if (.not. all(ieee_is_finite([voltage, g_par, g_perp]))) return
        if (g_perp <= 0 .or. g_par < g_perp) return
        call meridian_fourier_gradient(field%fs, modes, gradient, status)
        if (status /= 0) return
        self%field = field
        self%pressure_gradient = csc_matvec(gradient, field%p)
        self%potential_gradient = csc_matvec(gradient, potential)
        self%voltage = voltage; self%g_par = g_par; self%g_perp = g_perp
    end subroutine

    subroutine state_profile_evaluate(self, x, p, gradp, lambda, gradphi, status)
        class(nn_state_profile_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: p(:), gradp(:, :), lambda(:), gradphi(:, :)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2))
        real(dp) :: basis(3, self%field%fs%sp(1)%nloc)
        real(dp) :: xy(2), u, theta, xi, eta, b(3), j(3), angle, factor, parallel_p
        integer :: sts(size(x, 2)), i, t, n, n1, bk
        status = 1
        if (.not. allocated(self%pressure_gradient)) return
        if (size(x, 1) /= 3) return
        n1 = self%field%fs%sp(1)%ndof
        call self%field%sample(x, values, sts)
        if (any(sts /= 0)) return
        do i = 1, size(x, 2)
            call self%field%coordinates(self%field%centre, self%field%radial_scale, &
                                        x(:2, i), u, theta)
            call flux_mesh_locate(self%field%mesh, x(:2, i), u, theta, &
                                  t, xi, eta, status)
            if (status /= 0) return
            call meridian_point_basis(self%field%fs, self%field%pb, 1, t, &
                                      xi, eta, xy, basis, status)
            if (status /= 0) return
            p(i) = values(EQ_P, i); gradp(:, i) = 0; gradphi(:, i) = 0
            associate (dof => self%field%fs%sp(1)%dof(:, t), &
                       sign => self%field%fs%sp(1)%sgn(:, t))
                do bk = 0, size(self%pressure_gradient)/n1 - 1
                    factor = 1
                    n = (bk + 1)/2; angle = n*x(3, i)
                    if (bk > 0) then
                        if (mod(bk, 2) == 1) then
                            factor = cos(angle)
                        else
                            factor = sin(angle)
                        end if
                    end if
                    gradp(:, i) = gradp(:, i) + factor*matmul(basis, &
                                               sign*self%pressure_gradient(bk*n1 + dof))
                    gradphi(:, i) = gradphi(:, i) + factor*matmul(basis, &
                                              sign*self%potential_gradient(bk*n1 + dof))
                end do
            end associate
            b = [values(EQ_BR, i), values(EQ_BPHI, i), values(EQ_BZ, i)]
            status = 1
            if (dot_product(b, b) <= tiny(1.0_dp)) return
            call coupled_current(b(1), b(2), b(3), gradp(1, i), gradp(2, i), &
                             gradp(3, i), gradphi(1, i), gradphi(2, i), gradphi(3, i), &
                      self%voltage/(2*acos(-1.0_dp)*x(1, i)), self%g_par, self%g_perp, &
                                 j(1), j(2), j(3))
            call profile_current_diagnostic(b(1), b(2), b(3), j(1), j(2), j(3), &
                                gradp(1, i), gradp(2, i), gradp(3, i), 0.0_dp, 0.0_dp, &
                                            lambda(i), parallel_p)
        end do
        status = 0
    end subroutine

    subroutine nn_fit_source(seed, profile, chi_par, chi_perp, g_par, g_perp, &
                             voltage, fit, status)
        class(equilibrium_field_t), intent(in) :: seed
        class(nn_profile_t), intent(in) :: profile
        real(dp), intent(in) :: chi_par, chi_perp, g_par, g_perp, voltage
        type(nn_fitted_source_t), intent(out) :: fit
        integer, intent(out) :: status
        status = 1
        if (.not. all(ieee_is_finite([chi_par, chi_perp, g_par, g_perp, voltage]))) &
            return
        if (chi_perp <= 0 .or. chi_par < chi_perp) return
        if (g_perp <= 0 .or. g_par < g_perp) return
        allocate (fit%seed, source=seed)
        allocate (fit%profile, source=profile)
        fit%chi_par = chi_par; fit%chi_perp = chi_perp
        fit%g_par = g_par; fit%g_perp = g_perp; fit%voltage = voltage
        status = 0
    end subroutine

    logical function fitted_ready(self) result(ready)
        class(nn_fitted_source_t), intent(in) :: self
        ready = allocated(self%seed) .and. allocated(self%profile)
    end function

    subroutine fitted_data(self, x, flux, drive, status)
        class(nn_fitted_source_t), intent(in) :: self
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: flux(:, :), drive(:)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), p(size(x, 2))
        real(dp) :: gp(3, size(x, 2)), gf(3, size(x, 2)), lam(size(x, 2)), b(3)
        integer :: sts(size(x, 2)), i
        status = 1; flux = 0; drive = 0
        if (.not. self%ready()) return
        if (size(x, 1) /= 3) return
        if (.not. all(ieee_is_finite(x))) return
        if (any(x(1, :) <= 0)) return
        call self%seed%sample(x, values, sts)
        if (any(sts /= 0)) return
        call self%profile%evaluate(x, p, gp, lam, gf, status)
        if (status /= 0) return
        status = 1
        if (.not. all(ieee_is_finite(p))) return
        if (.not. all(ieee_is_finite(gp))) return
        if (.not. all(ieee_is_finite(gf))) return
        if (.not. all(ieee_is_finite(lam))) return
        do i = 1, size(x, 2)
            b = [values(EQ_BR, i), values(EQ_BPHI, i), values(EQ_BZ, i)]
            if (.not. all(ieee_is_finite(b))) return
            if (dot_product(b, b) <= tiny(1.0_dp)) return
            call anisotropic_flux(b(1), b(2), b(3), gp(1, i), gp(2, i), gp(3, i), &
                        self%chi_par, self%chi_perp, flux(1, i), flux(2, i), flux(3, i))
            call fitted_parallel_drive(b(1), b(2), b(3), gf(1, i), gf(2, i), &
                         gf(3, i), self%voltage/(2*acos(-1.0_dp)*x(1, i)), self%g_par, &
                                       lam(i), drive(i))
        end do
        if (.not. all(ieee_is_finite(flux))) return
        if (.not. all(ieee_is_finite(drive))) return
        status = 0
    end subroutine

    subroutine fitted_source(self, t, x, f0, f, status)
        class(nn_fitted_source_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: f0(:), f(:, :)
        integer, intent(out) :: status
        real(dp) :: drive(size(x, 2))
        f0 = 0; f = 0; status = 1
        if (t < 1) return
        call self%data(x, f, drive, status)
    end subroutine

    subroutine fitted_defects(self, field, x, weights, gradp, gradphi, voltage, &
                              report, status)
        ! Weighted physical RMS, with caller-declared positive volume weights.
        ! Report target mismatch separately from the seed current-law fit.
        ! These are diagnostics, not bounds or a current-solvability test.
        class(nn_fitted_source_t), intent(in) :: self
        class(equilibrium_field_t), intent(in) :: field
        real(dp), intent(in) :: x(:, :), weights(:), gradp(:, :), gradphi(:, :)
        real(dp), intent(in) :: voltage
        type(nn_profile_defect_t), intent(out) :: report
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), sv(EQ_VALUES, size(x, 2))
        real(dp) :: p(size(x, 2)), gp(3, size(x, 2)), gf(3, size(x, 2))
        real(dp) :: lam(size(x, 2)), drive(size(x, 2)), flux(3, size(x, 2))
        real(dp) :: b(3), bs(3), j(3), js(3), normb, dl, ds, par, ephi, seed_par
        integer :: sts(size(x, 2)), i
        status = 1
        if (size(x, 1) /= 3 .or. size(weights) /= size(x, 2)) return
        if (size(x, 2) == 0) return
        if (size(gradp, 1) /= 3 .or. size(gradp, 2) /= size(x, 2)) return
        if (size(gradphi, 1) /= 3 .or. size(gradphi, 2) /= size(x, 2)) return
        if (any(weights <= 0) .or. .not. all(ieee_is_finite(weights))) return
        if (.not. all(ieee_is_finite(gradp))) return
        if (.not. all(ieee_is_finite(gradphi))) return
        if (.not. ieee_is_finite(voltage)) return
        call self%data(x, flux, drive, status)
        if (status /= 0) return
        call self%profile%evaluate(x, p, gp, lam, gf, status)
        if (status /= 0) return
        call field%sample(x, values, sts)
        status = 1
        if (any(sts /= 0)) return
        call self%seed%sample(x, sv, sts)
        if (any(sts /= 0)) return
        do i = 1, size(x, 2)
            b = [values(EQ_BR, i), values(EQ_BPHI, i), values(EQ_BZ, i)]
            bs = [sv(EQ_BR, i), sv(EQ_BPHI, i), sv(EQ_BZ, i)]
            normb = dot_product(b, b)
            if (normb <= tiny(1.0_dp)) return
            ephi = voltage/(2*acos(-1.0_dp)*x(1, i))
            call coupled_current(b(1), b(2), b(3), gradp(1, i), gradp(2, i), &
                       gradp(3, i), gradphi(1, i), gradphi(2, i), gradphi(3, i), ephi, &
                                 self%g_par, self%g_perp, j(1), j(2), j(3))
            call coupled_current(bs(1), bs(2), bs(3), gp(1, i), gp(2, i), &
                                 gp(3, i), gf(1, i), gf(2, i), gf(3, i), &
                                 self%voltage/(2*acos(-1.0_dp)*x(1, i)), self%g_par, &
                                 self%g_perp, js(1), js(2), js(3))
            call profile_current_diagnostic(b(1), b(2), b(3), j(1), j(2), j(3), &
                       gradp(1, i), gradp(2, i), gradp(3, i), drive(i), lam(i), dl, par)
            call profile_current_diagnostic(bs(1), bs(2), bs(3), js(1), js(2), &
                    js(3), gp(1, i), gp(2, i), gp(3, i), drive(i), lam(i), ds, seed_par)
            report%pressure_rms = report%pressure_rms &
                                  + weights(i)*(values(EQ_P, i) - p(i))**2
            report%lambda_rms = report%lambda_rms + weights(i)*dl**2
            report%seed_current_fit_rms = report%seed_current_fit_rms + weights(i)*ds**2
         report%parallel_pressure_rms = report%parallel_pressure_rms + weights(i)*par**2
        end do
        report%pressure_rms = sqrt(report%pressure_rms/sum(weights))
        report%lambda_rms = sqrt(report%lambda_rms/sum(weights))
        report%seed_current_fit_rms = sqrt(report%seed_current_fit_rms/sum(weights))
        report%parallel_pressure_rms = sqrt(report%parallel_pressure_rms/sum(weights))
        if (.not. all(ieee_is_finite([report%pressure_rms, report%lambda_rms, &
                    report%seed_current_fit_rms, report%parallel_pressure_rms]))) return
        status = 0
    end subroutine
end module kin6d_nonnested_profile
