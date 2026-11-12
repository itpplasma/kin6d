module milca_vrmhd_torus
    !! Axisymmetric state of the circular-torus flux-pumping case: the
    !! cylinder case's p(psi_N) and FF'(psi_N) in a fixed-boundary
    !! Grad-Shafranov solve, its section mesh, sources and interpolant.
    !! See docs/vrmhd-torus.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    use milca_vrmhd_cylinder, only: vrmhd_cylinder_case_t, vrmhd_cylinder_profile_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_basis
    use milca_gs_symmetry, only: GS_ROTATION
    use milca_gs_profile_law, only: gs_profile_law_t, gs_law_primitives, gs_law_sources
    use milca_grad_shafranov, only: gs_state_t, initialize_gs, gs_axis_candidate, &
        gs_axis_safety_factor, gs_cell_data, gs_element_jet, gs_linear_response
    use milca_gs_nonlinear, only: solve_gs_cubic_profiles
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortfem_triangle_vector_interpolation, only: &
        interpolate_triangle_nedelec_second_kind
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    !> Axis of the solved state; shift = R_axis - R0.
    type, public :: vrmhd_torus_axis_t
        real(dp) :: position(2) = 0, flux = 0, shift = 0, q = 0
    end type

    ! The boundary callback of initialize_gs and FortFEM's interpolation
    ! callback take no context; the state is held here while they run.
    real(dp), save :: circle_major = 0
    type(gs_state_t), save :: held_flux, held_potential
    type(gs_profile_law_t), save :: held_law
    integer, save :: held_cell = 1

    public :: vrmhd_torus_law, vrmhd_torus_equilibrium, vrmhd_torus_mesh
    public :: vrmhd_torus_sources, vrmhd_torus_state
contains
    subroutine vrmhd_torus_law(c, prof, knots, law, status)
        !! Normalised law, piecewise linear in s = psi_N on `knots` uniform
        !! knots. psi = R0 int B_theta dr is the cylinder flux at R = R0, so
        !! dp/dpsi = T'/(R0 B_theta) and, from radial force balance and
        !! Ampere's law, FF' = -R0 (J_z + T'/B_theta).
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        integer, intent(in) :: knots
        type(gs_profile_law_t), intent(out) :: law
        integer, intent(out) :: status
        real(dp) :: pp(knots), ff(knots), s, w, dtr, btr, jz
        integer :: i, j, n
        status = 1
        n = size(prof%x)
        if (knots < 2 .or. n < 2) return
        j = 1
        do i = 1, knots
            s = real(i - 1, dp)/(knots - 1)
            do while (j < n - 1 .and. prof%psin(j + 1) < s)
                j = j + 1
            end do
            w = (s - prof%psin(j))/(prof%psin(j + 1) - prof%psin(j))
            dtr = (1 - w)*prof%dt_over_r(j) + w*prof%dt_over_r(j + 1)
            btr = (1 - w)*prof%btheta_over_r(j) + w*prof%btheta_over_r(j + 1)
            jz = (1 - w)*prof%jz(j) + w*prof%jz(j + 1)
            pp(i) = dtr/(c%r0*btr)
            ff(i) = -c%r0*(jz + dtr/btr)
        end do
        allocate (law%knots(knots), law%pp(4, knots - 1), law%ff(4, knots - 1))
        law%knots = [(real(i - 1, dp)/(knots - 1), i=1, knots)]
        law%pp = 0; law%ff = 0
        law%pp(1, :) = pp(:knots - 1); law%ff(1, :) = ff(:knots - 1)
        law%pp(2, :) = (pp(2:) - pp(:knots - 1))*(knots - 1)
        law%ff(2, :) = (ff(2:) - ff(:knots - 1))*(knots - 1)
        law%normalized = .true.
        law%p_edge = prof%t(n)
        law%f_edge = c%r0*prof%bz(n)
        law%axis = -c%r0*prof%psi(n)
        if (law%valid()) status = 0
    end subroutine

    subroutine vrmhd_torus_equilibrium(c, prof, edge, law, state, axis, status)
        !! Fixed-boundary solve in the circle of radius a = 1 about R0 with
        !! cubic elements of edge length about `edge`; mu0 = 1. Newton starts
        !! from the cylinder flux; from the solver's default start (psi = 0)
        !! it fails for these profiles.
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: edge
        type(gs_profile_law_t), intent(inout) :: law
        type(gs_state_t), intent(out) :: state
        type(vrmhd_torus_axis_t), intent(out) :: axis
        integer, intent(out) :: status
        real(dp), allocatable :: polygon(:, :)
        real(dp) :: hessian(2, 2), pressure, f, r, w
        integer :: nodes, i, j, cell, n
        status = 1
        if (edge <= 0 .or. edge > 1) return
        circle_major = c%r0
        nodes = max(8, nint(2*pi/edge))
        allocate (polygon(2, nodes))
        do i = 1, nodes
            call circle(real(i - 1, dp)/nodes, polygon(:, i))
        end do
        call initialize_gs(state, polygon, edge*edge/2, status, 3, curve=circle)
        if (status /= 0) return
        n = size(prof%x)
        do i = 1, size(state%psi)
            if (state%boundary(i)) then
                state%psi(i) = 0
                cycle
            end if
            r = min(norm2(state%points(:, i) - [c%r0, 0.0_dp]), 1.0_dp)*(n - 1)
            j = min(int(r) + 1, n - 1); w = r - (j - 1)
            state%psi(i) = c%r0*((1 - w)*prof%psi(j) + w*prof%psi(j + 1) - prof%psi(n))
        end do
        call solve_gs_cubic_profiles(state, 1.0_dp, law, status, 1.0e-11_dp, &
            quadrature_tolerance=1.0e-11_dp)
        if (status /= 0) return
        call gs_axis_candidate(state, axis%position, axis%flux, hessian, cell, status)
        if (status /= 0) return
        call gs_law_primitives(law, axis%flux, pressure, f, status)
        if (status /= 0) return
        call gs_axis_safety_factor(axis%position(1), f, hessian, axis%q, status)
        axis%shift = axis%position(1) - c%r0
    end subroutine

    subroutine vrmhd_torus_mesh(gs, mm, status)
        !! Straight section mesh on the vertices and triangles of the
        !! Grad-Shafranov mesh; triangle t is the same cell in both.
        type(gs_state_t), intent(in) :: gs
        type(meridian_mesh_t), intent(out) :: mm
        integer, intent(out) :: status
        type(mesh_2d_t) :: mesh
        mesh%n_vertices = maxval(gs%triangles(1:3, :))
        mesh%n_triangles = size(gs%triangles, 2); mesh%has_triangles = .true.
        mesh%vertices = gs%points(:, :mesh%n_vertices)
        mesh%triangles = gs%triangles(1:3, :)
        call meridian_mesh_straight(mesh, 1, mm, status)
    end subroutine

    subroutine vrmhd_torus_sources(op, c, prof, gs, law, source, voltage, status)
        !! Js = F' B_pol holds the poloidal current, as J_theta does in the
        !! cylinder; Qe = Qi = -(3/4) chi lap p keeps T = p(psi) under
        !! perpendicular conduction. The loop voltage V = 2 pi eta R J_phi is
        !! that of the axis, where eta = 1/S and R0 J_phi(R0) = -R0 J_z of the
        !! cylinder (vrmhd_torus_law); at p' = 0 it holds on every surface.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        type(gs_state_t), intent(in) :: gs
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(out) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: voltage
        integer, intent(out) :: status
        real(dp) :: basis(36, op%nlocal), x(2), weight, psi, grad(2), hessian(2, 2)
        real(dp) :: pp, ff, dpp, pressure, f, bpol(2), heat, jphi
        integer :: t, q, cell
        source = 0
        voltage = -2*pi*c%r0*prof%electric_field
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                if (status /= 0) return
                cell = t
                call flux_jet(gs, x, cell, psi, grad, hessian, status)
                if (status /= 0) return
                call gs_law_sources(law, psi, pp, ff, dpp=dpp)
                call gs_law_primitives(law, psi, pressure, f, status)
                if (status /= 0) return
                bpol = [-grad(2), grad(1)]/x(1)
                source(1:3, q, t) = ff/f*[bpol(1), 0.0_dp, bpol(2)]
                heat = -0.75_dp*c%chi_perpendicular*(dpp*sum(grad**2) &
                    + pp*(hessian(1, 1) + hessian(2, 2) + grad(1)/x(1)))
                jphi = x(1)*pp + ff/x(1)
                source(5, q, t) = heat - c%ohmic_fraction*voltage/(2*pi*x(1))*jphi
                source(6, q, t) = heat
            end do
        end do
        status = 0
    end subroutine

    subroutine vrmhd_torus_state(op, gs, law, u, status)
        !! Interpolant of the Grad-Shafranov state: rho = 1, v = 0,
        !! Te = Ti = p/2, A = (A_R, psi/R, A_Z) with
        !! A_pol = grad w x grad phi - F_edge ln R e_Z and Delta* w = F_edge - F,
        !! so curl A_pol = (F/R) e_phi, and J = F' B_pol + (R p' + FF'/R) e_phi
        !! (mu0 = 1). The vacuum part is exact; w carries only F - F_edge, so
        !! the second-derivative error of the cubic w stays small.
        type(vrmhd_section_t), intent(in) :: op
        type(gs_state_t), intent(in) :: gs
        type(gs_profile_law_t), intent(in) :: law
        real(dp), intent(out) :: u(op%ndof)
        integer, intent(out) :: status
        real(dp), allocatable :: xy(:, :), moments(:)
        real(dp) :: psi, grad(2), hessian(2, 2), pressure, f, pp, ff, vertices(2, 3)
        integer :: i, t, cell, field, nln
        status = 1
        u = 0
        if (op%symmetry%kind /= GS_ROTATION) return
        held_flux = gs; held_law = law
        call toroidal_potential(gs, law, held_potential, status)
        if (status /= 0) return
        cell = 1
        call meridian_lagrange_nodes(op%scalar_fs, xy)
        do i = 1, size(xy, 2)
            call flux_jet(gs, xy(:, i), cell, psi, grad, hessian, status)
            if (status /= 0) return
            call gs_law_primitives(law, psi, pressure, f, status)
            if (status /= 0) return
            u(op%offset(1) + i) = 1
            u(op%offset(5) + i) = pressure/2
            u(op%offset(6) + i) = pressure/2
        end do
        call meridian_lagrange_nodes(op%fs, xy)
        do i = 1, size(xy, 2)
            call flux_jet(gs, xy(:, i), cell, psi, grad, hessian, status)
            if (status /= 0) return
            call gs_law_sources(law, psi, pp, ff)
            u(op%offset(7) + op%fs%nn + i) = psi
            u(op%offset(8) + op%fs%nn + i) = xy(1, i)**2*pp + ff
        end do
        nln = op%fs%nln
        do field = 7, 8
            do t = 1, op%fs%mm%topology%n_triangles
                held_cell = t
                vertices = op%fs%mm%topology%vertices(:, &
                    op%fs%mm%topology%triangles(:, t))
                if (field == 7) then
                    call interpolate_triangle_nedelec_second_kind(vertices, &
                        op%fs%k, 2*op%fs%k + 4, section_potential, moments, status)
                else
                    call interpolate_triangle_nedelec_second_kind(vertices, &
                        op%fs%k, 2*op%fs%k + 4, section_current, moments, status)
                end if
                if (status /= 0) return
                status = 1
                if (.not. all(ieee_is_finite(moments))) return
                u(op%offset(field) + op%fs%sp(1)%dof(:nln, t)) = &
                    moments*op%fs%sp(1)%sgn(:nln, t)
            end do
        end do
        status = 0
    end subroutine

    subroutine toroidal_potential(gs, law, w, status)
        !! Delta* w = F_edge - F with w = 0 on the wall in the Grad-Shafranov
        !! element space: the weak form of psi with the source (F - F_edge)/R.
        type(gs_state_t), intent(in) :: gs
        type(gs_profile_law_t), intent(in) :: law
        type(gs_state_t), intent(out) :: w
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:), eta(:), qw(:), rhs(:), solution(:)
        real(dp) :: x(2), det, pressure, f
        real(dp) :: shape(size(gs%triangles, 1)), gradients(2, size(gs%triangles, 1))
        real(dp) :: hessians(2, 2, size(gs%triangles, 1))
        integer :: e, q
        call triangle_duffy_quadrature(10, xi, eta, qw, status)
        if (status /= 0) return
        allocate (rhs(size(gs%psi)))
        rhs = 0
        do e = 1, size(gs%triangles, 2)
            do q = 1, size(qw)
                call gs_cell_data(gs, e, [xi(q), eta(q)], x, det, shape, gradients, &
                    hessians, status)
                if (status /= 0) return
                call gs_law_primitives(law, dot_product(shape, &
                    gs%psi(gs%triangles(:, e))), pressure, f, status)
                if (status /= 0) return
                rhs(gs%triangles(:, e)) = rhs(gs%triangles(:, e)) &
                    + det*qw(q)*(f - law%f_edge)/x(1)*shape
            end do
        end do
        where (gs%boundary) rhs = 0
        call gs_linear_response(gs, rhs, solution, status, relative_tolerance=1.0e-12_dp)
        if (status /= 0) return
        w = gs
        w%psi = solution
    end subroutine

    subroutine flux_jet(state, x, cell, psi, grad, hessian, status)
        !! Jet of state%psi at x; `cell` is tried first and returns the cell found.
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: x(2)
        integer, intent(inout) :: cell
        real(dp), intent(out) :: psi, grad(2), hessian(2, 2)
        integer, intent(out) :: status
        integer :: e
        call gs_element_jet(state, cell, x(1), x(2), psi, grad, hessian, status)
        if (status == 0) return
        do e = 1, size(state%triangles, 2)
            call gs_element_jet(state, e, x(1), x(2), psi, grad, hessian, status)
            if (status /= 0) cycle
            cell = e
            return
        end do
    end subroutine

    subroutine section_potential(px, py, w)
        real(dp), intent(in) :: px, py
        real(dp), intent(out) :: w(2)
        real(dp) :: potential, grad(2), hessian(2, 2)
        integer :: cell, status
        cell = held_cell
        call flux_jet(held_potential, [px, py], cell, potential, grad, hessian, status)
        w = [-grad(2), grad(1)]/px - [0.0_dp, held_law%f_edge*log(px)]
        if (status /= 0) w = ieee_value(0.0_dp, ieee_quiet_nan)
    end subroutine

    subroutine section_current(px, py, w)
        real(dp), intent(in) :: px, py
        real(dp), intent(out) :: w(2)
        real(dp) :: psi, grad(2), hessian(2, 2), pp, ff, pressure, f
        integer :: cell, status
        cell = held_cell
        call flux_jet(held_flux, [px, py], cell, psi, grad, hessian, status)
        if (status == 0) call gs_law_primitives(held_law, psi, pressure, f, status)
        if (status /= 0) then
            w = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        call gs_law_sources(held_law, psi, pp, ff)
        w = ff/f*[-grad(2), grad(1)]/px
    end subroutine

    pure subroutine circle(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        point = [circle_major + cos(2*pi*t), sin(2*pi*t)]
    end subroutine
end module milca_vrmhd_torus
