module kin6d_gs_equilibrated
    !! Equilibrated-flux energy bound for the fixed-boundary GS problem
    !! -div(grad psi / R) = f, psi = 0 on the boundary, on affine or
    !! isoparametric P1--P3 cells. FortFEM builds sigma ~ grad psi / R with
    !! RT_p equilibration (p = element order, A = 1/R), so the cost is local
    !! and linear in the mesh size, with no global solve and no estimated
    !! constant. Geometry Bernstein hulls supply A_T,min <= 1/max_T R
    !! and the Friedrichs constant of the
    !! bounding box W x H, C_F = sqrt(R_max)/(pi sqrt(W^-2 + H^-2)). The
    !! sampled integrals, geometry positivity between quadrature points and
    !! arithmetic remain approximate until an independent checker encloses them.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_triangle_basis, only: triangle_bernstein_controls
    use kin6d_gs_error_leaves, only: gs_error_coefficient, gs_error_flux, &
        gs_error_friedrichs
    use kin6d_gs_majorant, only: gs_point_source_fn
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    use kin6d_gs_majorant, only: gs_point_source_t
    implicit none
    private
    public :: gs_equilibrated_majorant, equilibrated_flux_estimate_t

    integer, parameter, public :: GS_EQ_OK = 0
    integer, parameter, public :: GS_EQ_BAD_INPUT = 1
    integer, parameter, public :: GS_EQ_CURVED = 2
    integer, parameter, public :: GS_EQ_CELL_FAILED = 3
    integer, parameter, public :: GS_EQ_EQUILIBRATION_FAILED = 4

    type, extends(gs_point_source_t) :: gs_point_source_adapter_t
        procedure(gs_point_source_fn), pointer, nopass :: callback => null()
    contains
        procedure :: evaluate => gs_point_source_adapter_t_evaluate
    end type

    interface gs_equilibrated_majorant
        module procedure gs_equilibrated_majorant_object, &
            gs_equilibrated_majorant_procedure
    end interface

contains

    subroutine gs_equilibrated_majorant_object(state, psi, source, estimate, status)
        !! estimate%bound >= |||psi* - psi|||, |||v|||^2 = int |grad v|^2 / R,
        !! for the solution psi* of -div(grad psi*/R) = source%evaluate(cell, R, Z,
        !! psi(R, Z)) frozen at the represented psi, which must vanish on the
        !! boundary. For a psi-dependent source this is the frozen-source
        !! estimate; dividing by the Jacobian's coercivity remains the caller's.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:)
        class(gs_point_source_t), target, intent(in) :: source
        type(equilibrated_flux_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), w(:), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:), geometry(:, :, :)
        real(dp) :: point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp) :: values(10), lo(2), hi(2), cf, ph, gradient(2), controls(2, 10)
        integer :: e, q, nd, nt, stat, failed, axis

        status = GS_EQ_BAD_INPUT
        if (.not. state%ready) return
        if (size(psi) /= size(state%points, 2)) return
        if (.not. all(ieee_is_finite(psi))) return
        if (state%element_order < 1 .or. state%element_order > 3) return
        if (any(abs(pack(psi, state%boundary)) > tiny(1.0_dp))) return
        nd = size(state%triangles, 1)
        nt = size(state%triangles, 2)
        ! Exact for the degree 2p + 2 local mass and divergence products; the
        ! extra degrees reduce the quadrature error of 1/R and the source.
        call triangle_duffy_quadrature(2*state%element_order + 6, xi, eta, w, stat)
        status = GS_EQ_CELL_FAILED
        if (stat /= 0) return
        allocate (coefficient(size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt), geometry(2, nd, nt))
        lo = huge(1.0_dp); hi = -huge(1.0_dp); failed = 0
        !$omp parallel do schedule(static) private(axis, controls) &
        !$omp reduction(min:lo) reduction(max:hi, failed)
        do e = 1, nt
            geometry(:, :, e) = state%points(:, state%triangles(:, e))
            do axis = 1, 2
                call triangle_bernstein_controls(geometry(axis, :, e), &
                    controls(axis, :nd))
            end do
            if (minval(controls(1, :nd)) <= 0) then
                failed = 1
                cycle
            end if
            call gs_error_coefficient(maxval(controls(1, :nd)), lower(e))
            lo = min(lo, minval(controls(:, :nd), dim=2))
            hi = max(hi, maxval(controls(:, :nd), dim=2))
        end do
        !$omp end parallel do
        if (failed /= 0 .or. any(hi <= lo)) return
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) &
        !$omp private(q, point, det, shape, grads, hess, values, ph, stat, gradient, axis)
        do e = 1, nt
            values(:nd) = psi(state%triangles(:, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [xi(q), eta(q)], point, det, shape(:nd), &
                    grads(:, :nd), hess(:, :, :nd), stat)
                if (stat /= 0) then
                    failed = 1
                    cycle
                end if
                ph = dot_product(shape(:nd), values(:nd))
                call gs_error_coefficient(point(1), coefficient(q, e))
                f(q, e) = source%evaluate(e, point(1), point(2), ph)
                do axis = 1, 2
                    gradient(axis) = dot_product(grads(axis, :nd), &
                        values(:nd) - values(1))
                end do
                call gs_error_flux(point(1), gradient(1), gradient(2), &
                    flux(1, q, e), flux(2, q, e))
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        call gs_error_friedrichs(hi(1), hi(1) - lo(1), hi(2) - lo(2), cf)
        call equilibrated_flux_estimate(state%points, state%triangles(1:3, :), &
            state%element_order, xi, eta, w, coefficient, f, flux, lower, cf, &
            estimate, stat, geometry_nodes=geometry)
        status = GS_EQ_EQUILIBRATION_FAILED
        if (stat /= EQUILIBRATION_OK) return
        status = GS_EQ_OK
    end subroutine gs_equilibrated_majorant_object


    subroutine gs_equilibrated_majorant_procedure(state, psi, source, estimate, status)
        !! estimate%bound >= |||psi* - psi|||, |||v|||^2 = int |grad v|^2 / R,
        !! for the solution psi* of -div(grad psi*/R) = source(cell, R, Z,
        !! psi(R, Z)) frozen at the represented psi, which must vanish on the
        !! boundary. For a psi-dependent source this is the frozen-source
        !! estimate; dividing by the Jacobian's coercivity remains the caller's.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: psi(:)
        procedure(gs_point_source_fn) :: source
        type(equilibrated_flux_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        type(gs_point_source_adapter_t), target :: source_adapter

        source_adapter%callback => source
        call gs_equilibrated_majorant_object(state, &
            psi, &
            source_adapter, &
            estimate, &
            status)
    end subroutine

    function gs_point_source_adapter_t_evaluate(self, cell, r, z, psi) result(f)
        class(gs_point_source_adapter_t), intent(in) :: self
            integer, intent(in) :: cell
            real(dp), intent(in) :: r, z, psi
            real(dp) :: f
        f = self%callback(cell, r, z, psi)
    end function

end module kin6d_gs_equilibrated
