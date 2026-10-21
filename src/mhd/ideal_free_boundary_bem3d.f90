module kin6d_ideal_free_boundary_bem3d
    !! R2 adapter for the fixed coupling surface of the ideal-map shell.
    !! Plasma Gamma_p remains free and fully 3D; Gamma_c is fixed.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
                                          bem_pair_report_t, build_p2_trace_connectivity
    use kin6d_equilibrium_chart, only: ideal_chart_t
    use kin6d_triangle_basis, only: triangle_basis_jet
    use kin6d_equilibrium_boundary_3d, only: gamma_c_3d_t, gamma_c_3d_build
    use kin6d_ideal_free_boundary, only: ideal_fb_t, ideal_fb_attach_exterior_3d, &
                                   ideal_fb_mode_weight, IDEAL_FB_OK, IDEAL_FB_BAD_INPUT
    implicit none
    private
    public :: ideal_fb_use_bem3d
    type, extends(bem_surface_t) :: coupling_surface_t
        class(ideal_chart_t), allocatable :: chart
        real(dp), allocatable :: edge(:, :, :), parameters(:, :, :)
        integer, allocatable :: panel_edge(:)
        integer :: order = 0
    contains
        procedure :: evaluate => evaluate_surface
    end type
contains
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
        real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
        integer :: e, k
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
        a = matmul(d1, qd(:, 1)); b = matmul(d1, qd(:, 2)); status = 0
    end subroutine
    subroutine ideal_fb_use_bem3d(fb, nphi, options, report, status)
        !! Existing Fourier drive is converted with its physical full-turn
        !! pairing. Set it BEFORE attaching. No scalar-potential coil fit here.
        type(ideal_fb_t), intent(inout) :: fb
        integer, intent(in) :: nphi
        type(bem_pair_options_t), intent(in) :: options
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(coupling_surface_t) :: surface
        type(gamma_c_3d_t) :: boundary
        integer, allocatable :: faces(:, :), td(:, :), fd(:, :), number(:)
        real(dp), allocatable :: restriction(:, :), drive(:), load(:)
        real(dp) :: p0, p1, phi, t, w(0:3), dw(0:3), uv(2), trig
        real(dp), parameter :: nodes(2, 6) = reshape( &
                               [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                                0.5_dp, 0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp], [2, 6])
        integer :: e, j, p, i, m, k, next, after, ids(4), col, row, nf
        status = IDEAL_FB_BAD_INPUT; report = bem_pair_report_t()
        if (nphi < 4 .or. nphi <= 2*fb%ntor*fb%nfp) return
        if (fb%order /= 2 .and. fb%order /= 3) return
        if (allocated(fb%exterior_3d)) return
        nf = 2*fb%nedge*nphi
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
                next = mod(e, fb%nedge) + 1
                ids = [e + (j - 1)*fb%nedge, e + (after - 1)*fb%nedge, &
                       next + (j - 1)*fb%nedge, next + (after - 1)*fb%nedge]
                p = p + 1; faces(:, p) = ids([1, 2, 3]); surface%panel_edge(p) = e
     surface%parameters(:, :, p) = reshape([0.0_dp, p0, 0.0_dp, p1, 1.0_dp, p0], [2, 3])
                p = p + 1; faces(:, p) = ids([2, 4, 3]); surface%panel_edge(p) = e
     surface%parameters(:, :, p) = reshape([0.0_dp, p1, 1.0_dp, p1, 1.0_dp, p0], [2, 3])
            end do
        end do
        call build_p2_trace_connectivity(faces, td, fd, status)
        if (status /= 0) return
        allocate (restriction(maxval(td), fb%nt*fb%nmode), number(fb%ns))
        restriction = 0; number = 0
        do i = 1, fb%nt
            number(fb%trace(i)) = i
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
   call gamma_c_3d_build(surface, faces, restriction, options, boundary, report, status)
        if (status /= 0) return
        drive = reshape(fb%drive, [fb%nt*fb%nmode])
        allocate (load(fb%nt*fb%nmode))
        do m = 1, fb%nmode
      load((m - 1)*fb%nt + 1:m*fb%nt) = fb%nfp*ideal_fb_mode_weight(fb, m)*fb%load(:, m)
        end do
        call ideal_fb_attach_exterior_3d(fb, boundary, drive, load, status)
    end subroutine
end module
