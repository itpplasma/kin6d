module kin6d_gs_equilibrated
    !! Equilibrated-flux energy bound for the fixed-boundary GS problem
    !! -div(grad psi / R) = f, psi = 0 on the boundary, on a straight-sided
    !! mesh. The flux sigma ~ grad psi / R is built by FortFEM's patchwise
    !! RT_p equilibration (p = element order, A = 1/R), so the cost is local
    !! and linear in the mesh size, with no global solve and no estimated
    !! constant. The guaranteed constants are A_T,min = 1/max_T R (attained
    !! at a vertex of the affine cell) and the Friedrichs constant of the
    !! bounding box W x H, C_F = sqrt(R_max)/(pi sqrt(W^-2 + H^-2)). The
    !! remaining non-rigorous ingredient is the quadrature of the
    !! non-polynomial 1/R and source integrands.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov, only: gs_state_t, gs_cell_data
    use kin6d_gs_majorant, only: gs_point_source_fn
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    implicit none
    private
    public :: gs_equilibrated_majorant, equilibrated_flux_estimate_t

    integer, parameter, public :: GS_EQ_OK = 0
    integer, parameter, public :: GS_EQ_BAD_INPUT = 1
    integer, parameter, public :: GS_EQ_CURVED = 2
    integer, parameter, public :: GS_EQ_CELL_FAILED = 3
    integer, parameter, public :: GS_EQ_EQUILIBRATION_FAILED = 4

contains

    subroutine gs_equilibrated_majorant(state, psi, source, estimate, status)
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
        real(dp), allocatable :: xi(:), eta(:), w(:), coefficient(:, :), f(:, :)
        real(dp), allocatable :: flux(:, :, :), lower(:)
        real(dp) :: point(2), det, shape(10), grads(2, 10), hess(2, 2, 10)
        real(dp) :: values(10), lo(2), hi(2), cf, ph
        integer :: e, q, nd, nt, stat, failed

        status = GS_EQ_BAD_INPUT
        if (.not. state%ready .or. size(psi) /= size(state%points, 2)) return
        if (state%element_order < 1 .or. state%element_order > 3) return
        status = GS_EQ_CURVED
        if (allocated(state%curved_cell)) then
            if (any(state%curved_cell)) return
        end if
        nd = size(state%triangles, 1)
        nt = size(state%triangles, 2)
        ! Exact for the degree 2p + 2 local mass and divergence products; the
        ! extra degrees reduce the quadrature error of 1/R and the source.
        call triangle_duffy_quadrature(2*state%element_order + 6, xi, eta, w, stat)
        status = GS_EQ_CELL_FAILED
        if (stat /= 0) return
        allocate (coefficient(size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt))
        failed = 0
        !$omp parallel do schedule(static) reduction(max:failed) &
        !$omp private(q, point, det, shape, grads, hess, values, ph, stat)
        do e = 1, nt
            lower(e) = 1.0_dp/maxval(state%points(1, state%triangles(1:3, e)))
            values(:nd) = psi(state%triangles(:, e))
            do q = 1, size(w)
                call gs_cell_data(state, e, [xi(q), eta(q)], point, det, shape(:nd), &
                    grads(:, :nd), hess(:, :, :nd), stat)
                if (stat /= 0) failed = 1
                ph = dot_product(shape(:nd), values(:nd))
                coefficient(q, e) = 1.0_dp/point(1)
                f(q, e) = source(e, point(1), point(2), ph)
                flux(:, q, e) = matmul(grads(:, :nd), values(:nd) - values(1)) &
                    /point(1)
            end do
        end do
        !$omp end parallel do
        if (failed /= 0) return
        lo = minval(state%points, dim=2)
        hi = maxval(state%points, dim=2)
        cf = sqrt(hi(1))/(acos(-1.0_dp)*sqrt(sum(1.0_dp/(hi - lo)**2)))
        call equilibrated_flux_estimate(state%points, state%triangles(1:3, :), &
            state%element_order, xi, eta, w, coefficient, f, flux, lower, cf, &
            estimate, stat)
        status = GS_EQ_EQUILIBRATION_FAILED
        if (stat /= EQUILIBRATION_OK) return
        status = GS_EQ_OK
    end subroutine gs_equilibrated_majorant

end module kin6d_gs_equilibrated
