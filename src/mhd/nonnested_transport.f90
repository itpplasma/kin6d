module milca_nonnested_transport
    ! NN2 first increment: steady anisotropic pressure transport in a given
    ! magnetic field,
    !     -div(K grad p) = S in Omega,  p = p_b on the boundary,
    !     K = chi_par b b + chi_perp (I - b b),  b = B/|B|,
    ! with p in the Fourier-extended V0 space of the meridional complex.
    ! B is any equilibrium_field_t (NN1 free-B state, GS state or an analytic
    ! field); K comes from the generated anisotropic_conductivity leaf, whose
    ! JVP/VJP in B is the transport-field block of the later coupled Newton.
    ! Generic Fourier-V0 assembly lives in FortFEM; the sparse solve is the
    ! NN1 symmetric Dirichlet solve. See docs/equilibrium-nonnested-nn2.md.
    ! With estimate_error, the N = 0 block carries the equilibrated energy
    ! bound of milca_nonnested_error (nt_axisymmetric_error); a phi-dependent
    ! K with modes > 0 uses the full Fourier-Galerkin residual bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t, meridian_fourier_diffusion, &
        meridian_fourier_load, meridian_fourier_anisotropic_t, &
        meridian_fourier_anisotropic_apply
    use fortsparse, only: csc_t, csc_matvec
    use milca_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_VALUES
    use milca_equilibrium_nonnested, only: nn_factor_t, nn_factor, &
        nn_factor_solve_real, nn_factor_free
    use milca_nonnested_leaves, only: anisotropic_conductivity
    use milca_nonnested_error, only: nn_error_t, nt_axisymmetric_error
    implicit none
    private
    public :: nt_conductivity_t, nt_report_t, nt_solve

    ! Conductivity of a given field. The field must outlive the solve and be
    ! safe for concurrent sampling (element loops are OpenMP parallel).
    type, extends(meridian_fourier_anisotropic_t) :: nt_conductivity_t
        class(equilibrium_field_t), pointer :: field => null()
        real(dp) :: chi_par = 1, chi_perp = 1
    contains
        procedure :: tensor => nt_tensor
        procedure :: direction => nt_direction
    end type

    type :: nt_report_t
        real(dp) :: residual = 0, chi_ratio = 0
        integer :: ndof = 0, modes = 0, nphi = 0, refinements = 0
        character(32) :: claim = 'candidate'
        ! N = 0 energy bound (int R grad e . K_pol grad e dR dZ)^(1/2).
        type(nn_error_t) :: error
        integer :: error_status = -1
    end type
contains
    subroutine nt_tensor(self, t, x, k, status)
        ! K at points x(3, m) = (R, Z, phi) in physical (R, phi, Z) order.
        class(nt_conductivity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), b(3)
        integer :: statuses(size(x, 2)), p
        status = 1
        k = 0
        if (t < 1 .or. .not. associated(self%field)) return
        call self%field%sample(x, values, statuses)
        if (any(statuses /= 0)) return
        do p = 1, size(x, 2)
            b = [values(EQ_BR, p), values(EQ_BPHI, p), values(EQ_BZ, p)]
            if (.not. all(ieee_is_finite(b))) return
            if (dot_product(b, b) <= tiny(1.0_dp)) return
            call anisotropic_conductivity(b(1), b(2), b(3), self%chi_par, &
                self%chi_perp, k(1, 1, p), k(1, 2, p), k(1, 3, p), k(2, 2, p), &
                k(2, 3, p), k(3, 3, p))
            k(2, 1, p) = k(1, 2, p); k(3, 1, p) = k(1, 3, p)
            k(3, 2, p) = k(2, 3, p)
        end do
        status = 0
    end subroutine

    subroutine nt_direction(self, t, x, b, k_par, k_perp, status)
        class(nt_conductivity_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: b(:, :), k_par(:), k_perp(:)
        integer, intent(out) :: status
        real(dp) :: values(EQ_VALUES, size(x, 2)), bb
        integer :: statuses(size(x, 2)), p
        status = 1
        b = 0; k_par = self%chi_par; k_perp = self%chi_perp
        if (t < 1 .or. .not. associated(self%field)) return
        call self%field%sample(x, values, statuses)
        if (any(statuses /= 0)) return
        do p = 1, size(x, 2)
            b(:, p) = [values(EQ_BR, p), values(EQ_BPHI, p), values(EQ_BZ, p)]
            if (.not. all(ieee_is_finite(b(:, p)))) return
            bb = norm2(b(:, p))
            if (bb <= sqrt(tiny(1.0_dp))) return
            b(:, p) = b(:, p)/bb
        end do
        status = 0
    end subroutine

    subroutine nt_solve(fs, modes, nphi, conductivity, source, boundary, p, &
                        report, status, estimate_error)
        ! Dirichlet data are the boundary entries of boundary(:) in the
        ! FortFEM Fourier layout; p has the same layout. Derivatives with
        ! respect to S, p_b and chi are the same symmetric solve
        ! (nn_linear_solve) with differentiated loads; the B derivative uses
        ! anisotropic_conductivity_jvp/_vjp at the equation level.
        type(meridian_feec_t), intent(in) :: fs
        integer, intent(in) :: modes, nphi
        type(nt_conductivity_t), intent(in) :: conductivity
        class(meridian_fourier_source_t), intent(in) :: source
        real(dp), intent(in) :: boundary(:)
        real(dp), allocatable, intent(out) :: p(:)
        type(nt_report_t), intent(out) :: report
        integer, intent(out) :: status
        logical, intent(in), optional :: estimate_error
        type(csc_t) :: a, aabs
        real(dp), allocatable :: rhs(:), r(:)
        real(dp), allocatable :: dr(:), ap(:)
        logical, allocatable :: mask(:)
        type(nn_factor_t) :: factor
        integer, parameter :: max_refine = 30
        real(dp) :: step, last
        integer :: st, nb, b, ndof, it
        status = 1
        if (conductivity%chi_perp <= 0) return
        if (conductivity%chi_par < conductivity%chi_perp) return
        if (.not. all(ieee_is_finite([conductivity%chi_par, &
            conductivity%chi_perp]))) return
        nb = 2*modes + 1
        ndof = fs%n0*nb
        if (size(boundary) /= ndof) return
        if (.not. all(ieee_is_finite(boundary))) return
        call meridian_fourier_diffusion(fs, modes, nphi, conductivity, a, st)
        if (st /= 0) return
        call meridian_fourier_load(fs, modes, nphi, source, rhs, st)
        if (st /= 0) return
        allocate (mask(ndof), dr(ndof))
        do b = 0, nb - 1
            mask(b*fs%n0 + 1:(b + 1)*fs%n0) = fs%sp(0)%bnd
        end do
        p = merge(boundary, 0.0_dp, mask)
        r = rhs - csc_matvec(a, p)
        where (mask) r = 0
        call nn_factor(a, mask, factor, st)
        if (st /= 0) return
        call nn_factor_solve_real(factor, r, dr, st)
        if (st /= 0) then
            call nn_factor_free(factor); return
        end if
        p = p + dr
        ! At large chi_par/chi_perp, rounding of the assembled entries (size
        ! eps chi_par, unstructured) moves p by about the ratio times eps
        ! times the perpendicular condition, whatever the solver. The
        ! refinement residual is therefore the factored matrix-free form,
        ! whose rounding stays in the range of the parallel operator; the
        ! assembled factor only preconditions.
        aabs = a
        aabs%val = abs(aabs%val)
        last = huge(1.0_dp)
        do it = 1, max_refine
            call meridian_fourier_anisotropic_apply(fs, modes, nphi, &
                conductivity, p, ap, st)
            if (st == 0) then
                r = rhs - ap
                where (mask) r = 0
                call nn_factor_solve_real(factor, r, dr, st)
            end if
            if (st /= 0) then
                call nn_factor_free(factor); return
            end if
            step = maxval(abs(dr))/max(tiny(1.0_dp), maxval(abs(p)))
            ! A growing correction means the factor no longer contracts.
            if (step >= last) exit
            p = p + dr
            report%refinements = it
            if (step <= 4*epsilon(1.0_dp) .or. step >= last/2) exit
            last = step
        end do
        call nn_factor_free(factor)
        r = rhs - csc_matvec(a, p)
        where (mask) r = 0
        report%residual = norm2(r)/max(tiny(1.0_dp), &
            norm2(csc_matvec(aabs, abs(p))) + norm2(rhs))
        report%chi_ratio = conductivity%chi_par/conductivity%chi_perp
        report%ndof = ndof; report%modes = modes; report%nphi = nphi
        if (.not. all(ieee_is_finite(p))) return
        status = 0
        report%error%claim = 'unavailable-not-requested'
        if (.not. present(estimate_error)) return
        if (.not. estimate_error) return
        ! A failed estimate does not invalidate the solve; it is reported.
        call nt_axisymmetric_error(fs, modes, nphi, conductivity, &
            conductivity%chi_perp, source, p, report%error, report%error_status)
    end subroutine
end module milca_nonnested_transport
