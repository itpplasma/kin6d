module kin6d_nonnested_transport
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
    ! bound of kin6d_nonnested_error (nt_axisymmetric_error); a phi-dependent
    ! K with modes > 0 couples the modes and reports no bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_feec_t
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t, &
        meridian_fourier_source_t, meridian_fourier_diffusion, &
        meridian_fourier_load
    use fortsparse, only: csc_t, csc_matvec
    use kin6d_equilibrium_state, only: equilibrium_field_t, EQ_BR, EQ_BZ, &
        EQ_BPHI, EQ_VALUES
    use kin6d_equilibrium_nonnested, only: nn_linear_solve
    use kin6d_nonnested_leaves, only: anisotropic_conductivity
    use kin6d_nonnested_error, only: nn_error_t, nt_axisymmetric_error
    implicit none
    private
    public :: nt_conductivity_t, nt_report_t, nt_solve

    ! Conductivity of a given field. The field must outlive the solve and be
    ! safe for concurrent sampling (element loops are OpenMP parallel).
    type, extends(meridian_fourier_tensor_t) :: nt_conductivity_t
        class(equilibrium_field_t), pointer :: field => null()
        real(dp) :: chi_par = 1, chi_perp = 1
    contains
        procedure :: tensor => nt_tensor
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
        complex(dp), allocatable :: delta(:)
        logical, allocatable :: mask(:)
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
        allocate (mask(ndof), delta(ndof))
        do b = 0, nb - 1
            mask(b*fs%n0 + 1:(b + 1)*fs%n0) = fs%sp(0)%bnd
        end do
        p = merge(boundary, 0.0_dp, mask)
        r = rhs - csc_matvec(a, p)
        where (mask) r = 0
        call nn_linear_solve(a, mask, cmplx(r, 0.0_dp, dp), delta, st)
        if (st /= 0) return
        p = p + real(delta, dp)
        ! Large chi_par/chi_perp makes the system ill conditioned. The
        ! reported residual is the normwise backward error
        ! |r| / (| |A| |p| | + |rhs|), whose rounding floor does not grow
        ! with the ratio; up to two refinement steps restore it.
        aabs = a
        aabs%val = abs(aabs%val)
        do it = 0, 2
            r = rhs - csc_matvec(a, p)
            where (mask) r = 0
            report%residual = norm2(r)/max(tiny(1.0_dp), &
                norm2(csc_matvec(aabs, abs(p))) + norm2(rhs))
            report%refinements = it
            if (report%residual <= 1.0e-14_dp .or. it == 2) exit
            call nn_linear_solve(a, mask, cmplx(r, 0.0_dp, dp), delta, st)
            if (st /= 0) return
            p = p + real(delta, dp)
        end do
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
end module kin6d_nonnested_transport
