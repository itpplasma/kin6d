module kin6d_ideal_free_boundary_bem3d
    !! R2 adapter for the fixed coupling surface of the ideal-map shell.
    !! Plasma Gamma_p remains free and fully 3D; Gamma_c is fixed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
                                          bem_pair_report_t, build_p2_trace_connectivity
    use fortfem_generated_bem_frame, only: bem_frame
    use kin6d_equilibrium_chart, only: ideal_chart_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t, gamma_c_3d_build, &
        gamma_c_3d_build_periodic
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_attach_exterior_3d, &
                                   ideal_fb_full_turn_traces, IDEAL_FB_BAD_INPUT
    implicit none
    private
    public :: ideal_fb_use_bem3d, ideal_fb_surface, coupling_surface_t
    type, extends(bem_surface_t) :: coupling_surface_t
        class(ideal_chart_t), allocatable :: chart
        real(dp), allocatable :: edge(:, :, :), parameters(:, :, :)
        integer, allocatable :: panel_edge(:)
        integer :: order = 0
    contains
        procedure :: evaluate => evaluate_surface
        procedure,private :: assign_surface
        generic :: assignment(=)=>assign_surface
    end type
contains
    subroutine assign_surface(lhs, rhs)
        class(coupling_surface_t), intent(out) :: lhs
        type(coupling_surface_t), intent(in) :: rhs
        if (allocated(rhs%chart)) allocate (lhs%chart, source=rhs%chart)
        if (allocated(rhs%edge)) lhs%edge = rhs%edge
        if (allocated(rhs%parameters)) lhs%parameters = rhs%parameters
        if (allocated(rhs%panel_edge)) lhs%panel_edge = rhs%panel_edge
        lhs%order = rhs%order
    end subroutine
    pure subroutine edge_basis(order, t, w, dw)
        integer, intent(in) :: order
        real(dp), intent(in) :: t
        real(dp), intent(out) :: w(0:order), dw(0:order)
        real(dp) :: values(10), gradient(2, 10)
        integer :: index(4), nloc
        nloc = 6; index = [1, 4, 2, 0]
        if (order == 3) then
            nloc = 10; index = [1, 4, 5, 2]
        end if
        call triangle_basis_jet([t, 0.0_dp], values(:nloc), gradient(:, :nloc))
        w = values(index(:order + 1)); dw = gradient(1, index(:order + 1))
    end subroutine
    pure subroutine evaluate_surface(self, panel, uv, x, a, b, status)
        class(coupling_surface_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        real(dp) :: parameter(2), deriv(2, 2), q(3), qd(3, 2), w(0:3), dw(0:3)
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), jac, normal(3)
        integer :: e, k
        status = IDEAL_FB_BAD_INPUT; x = 0; a = 0; b = 0
        if (.not. allocated(self%panel_edge)) return
        if (panel < 1 .or. panel > size(self%panel_edge)) return
        if (.not. all(ieee_is_finite(uv))) return
        e = self%panel_edge(panel)
        parameter = self%parameters(:, 1, panel) + &
                 uv(1)*(self%parameters(:, 2, panel) - self%parameters(:, 1, panel)) + &
                    uv(2)*(self%parameters(:, 3, panel) - self%parameters(:, 1, panel))
        deriv(:, 1) = self%parameters(:, 2, panel) - self%parameters(:, 1, panel)
        deriv(:, 2) = self%parameters(:, 3, panel) - self%parameters(:, 1, panel)
        call edge_basis(self%order, parameter(1), w(:self%order), dw(:self%order))
        q(1:2) = matmul(self%edge(:, :, e), w(:self%order)); q(3) = parameter(2)
        do k = 1, 2
            qd(1:2, k) = matmul(self%edge(:, :, e), dw(:self%order))*deriv(1, k)
            qd(3, k) = deriv(2, k)
        end do
        call self%chart%jet(q, x, d1, d2, d3)
        a = matmul(d1, qd(:, 1)); b = matmul(d1, qd(:, 2))
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        call bem_frame(a(1), a(2), a(3), b(1), b(2), b(3), jac, &
                       normal(1), normal(2), normal(3))
        if (.not. ieee_is_finite(jac) .or. jac <= 0) return
        if (.not. all(ieee_is_finite(normal))) return
        status = 0
    end subroutine
    subroutine ideal_fb_surface(fb, nphi, surface, faces, restriction, status, subdivisions)
        !! Fixed polynomial geometry, independent P2 trace subdivision.
        !! Reference boundary edges are counter-clockwise, as in gamma_c_chain.
        !! The map orientation reverses face AND parameter order together.
        type(ideal_fb_t), intent(in) :: fb
        integer, intent(in) :: nphi
        type(coupling_surface_t), intent(out) :: surface
        integer, allocatable, intent(out) :: faces(:, :)
        real(dp), allocatable, intent(out) :: restriction(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: subdivisions
        integer, allocatable :: td(:, :), fd(:, :), number(:)
        real(dp) :: p0, p1, t0, t1, phi, t, w(0:3), dw(0:3), uv(2), trig
        real(dp) :: x(3), a(3), b(3)
        real(dp), parameter :: nodes(2, 6) = reshape( &
                               [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                                0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 6])
        integer :: e, j, p, i, m, k, next, after, ids(4), col, row, nf, sub, ntheta, segment, v
        integer :: refine, perm(3)
        status = IDEAL_FB_BAD_INPUT
        refine = 1
        if (present(subdivisions)) refine = subdivisions
        if (refine < 1) return
        if (fb%nfp < 1 .or. fb%ntor < 0) return
        if (fb%nmode /= 2*fb%ntor + 1) return
        if (nphi < 4 .or. nphi <= 2*fb%ntor*fb%nfp) return
        if (fb%order /= 2 .and. fb%order /= 3) return
        if (.not. allocated(fb%plasma%chart)) return
        if (abs(fb%plasma%sigma) /= 1) return
        if (fb%nedge < 3 .or. fb%ns < 1 .or. fb%nt < 1) return
        if (.not. allocated(fb%reference) .or. .not. allocated(fb%edge)) return
        if (.not. allocated(fb%trace)) return
        if (any(shape(fb%reference) /= [2, fb%ns])) return
        if (any(shape(fb%edge) /= [fb%order + 1, fb%nedge])) return
        if (lbound(fb%edge, 1) /= 0) return
        if (size(fb%trace) /= fb%nt) return
        if (minval(fb%edge) < 1 .or. maxval(fb%edge) > fb%ns) return
        if (minval(fb%trace) < 1 .or. maxval(fb%trace) > fb%ns) return
        if (.not. all(ieee_is_finite(fb%reference))) return
        do e = 1, fb%nedge
            next = mod(e, fb%nedge) + 1
            if (fb%edge(fb%order, e) /= fb%edge(0, next)) return
        end do
        ntheta = fb%nedge*refine
        nf = 2*ntheta*nphi
        perm = [1, 2, 3]
        if (fb%plasma%sigma > 0) perm = [1, 3, 2]
        surface%order = fb%order
        allocate (surface%chart, source=fb%plasma%chart)
        allocate (surface%edge(2, 0:fb%order, fb%nedge), &
                  surface%parameters(2, 3, nf), surface%panel_edge(nf), faces(3, nf))
        do e = 1, fb%nedge
            surface%edge(:, :, e) = fb%reference(:, fb%edge(:, e))
        end do
        p = 0
        do j = 1, nphi
            after = mod(j, nphi) + 1
            p0 = 2*acos(-1.0_dp)*real(j - 1, dp)/nphi
            p1 = 2*acos(-1.0_dp)*real(j, dp)/nphi
            do e = 1, fb%nedge
                do sub = 1, refine
                    segment = (e - 1)*refine + sub
                    next = mod(segment, ntheta) + 1
                    t0 = real(sub - 1, dp)/refine; t1 = real(sub, dp)/refine
                    ids = [segment + (j - 1)*ntheta, segment + (after - 1)*ntheta, &
                           next + (j - 1)*ntheta, next + (after - 1)*ntheta]
                    p = p + 1; faces(:, p) = ids([1, 2, 3]); surface%panel_edge(p) = e
                    surface%parameters(:, :, p) = reshape([t0, p0, t0, p1, t1, p0], [2, 3])
                    faces(:, p) = faces(perm, p)
                    surface%parameters(:, :, p) = surface%parameters(:, perm, p)
                    p = p + 1; faces(:, p) = ids([2, 4, 3]); surface%panel_edge(p) = e
                    surface%parameters(:, :, p) = reshape([t0, p1, t1, p1, t1, p0], [2, 3])
                    faces(:, p) = faces(perm, p)
                    surface%parameters(:, :, p) = surface%parameters(:, perm, p)
                end do
            end do
        end do
        do p = 1, nf
            call surface%evaluate(p, [1.0_dp/3, 1.0_dp/3], x, a, b, status)
            if (status /= 0) return
        end do
        call build_p2_trace_connectivity(faces, td, fd, status)
        if (status /= 0) return
        status = IDEAL_FB_BAD_INPUT
        allocate (restriction(maxval(td), fb%nt*fb%nmode), number(fb%ns))
        restriction = 0; number = 0
        do i = 1, fb%nt
            if (number(fb%trace(i)) /= 0) return
            number(fb%trace(i)) = i
        end do
        do e = 1, fb%nedge
            do v = 0, fb%order
                if (number(fb%edge(v, e)) == 0) return
            end do
        end do
        do p = 1, nf
            e = surface%panel_edge(p)
            do i = 1, 6
                uv = surface%parameters(:, 1, p) + nodes(1, i)* &
                     (surface%parameters(:, 2, p) - surface%parameters(:, 1, p)) + &
                 nodes(2, i)*(surface%parameters(:, 3, p) - surface%parameters(:, 1, p))
                t = uv(1); phi = uv(2); row = td(i, p)
                restriction(row, :) = 0
                call edge_basis(fb%order, t, w(:fb%order), dw(:fb%order))
                do m = 1, fb%nmode
                    trig = 1
                    if (m > 1) then
                        k = m/2; trig = cos(k*fb%nfp*phi)
                        if (mod(m, 2) == 1) trig = sin(k*fb%nfp*phi)
                    end if
                    do j = 0, fb%order
                        col = number(fb%edge(j, e)) + (m - 1)*fb%nt
                        restriction(row, col) = w(j)*trig
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine
    subroutine ideal_fb_use_bem3d(fb, nphi, options, report, status, subdivisions, &
                                  periodic)
        !! Existing Fourier drive uses its physical full-turn pairing.
        !! Set it before attaching. Gamma_c is fixed; no moving-BEM derivative.
        !! periodic: assemble only the target pairs of one field period
        !! (block-circulant, gamma_c_3d_build_periodic); needs nphi in nfp*Z
        !! and an nfp-symmetric chart. The same R^T S R up to the pair
        !! quadrature tolerance, at 1/nfp of the pair work and 1/nfp^3 of
        !! the dense DtN cost.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(in) :: nphi
        type(bem_pair_options_t), intent(in) :: options
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        integer, intent(in), optional :: subdivisions
        logical, intent(in), optional :: periodic
        type(coupling_surface_t) :: surface
        type(gamma_c_3d_t) :: boundary
        integer, allocatable :: faces(:, :)
        real(dp), allocatable :: restriction(:, :), drive(:), load(:)
        logical :: sectors
        status = IDEAL_FB_BAD_INPUT; report = bem_pair_report_t()
        sectors = .false.
        if (present(periodic)) sectors = periodic
        if (allocated(fb%exterior_3d)) return
        if (sectors .and. mod(nphi, max(fb%nfp, 1)) /= 0) return
        call ideal_fb_surface(fb, nphi, surface, faces, restriction, status, subdivisions)
        if (status /= 0) return
        if (sectors) then
            call gamma_c_3d_build_periodic(surface, faces, fb%nfp, restriction, &
                                           options, boundary, report, status)
        else
            call gamma_c_3d_build(surface, faces, restriction, options, boundary, &
                                  report, status)
        end if
        if (status /= 0) return
        allocate (drive(fb%nt*fb%nmode), load(fb%nt*fb%nmode))
        call ideal_fb_full_turn_traces(fb, fb%drive, fb%load, drive, load)
        call ideal_fb_attach_exterior_3d(fb, boundary, drive, load, status)
    end subroutine
end module
