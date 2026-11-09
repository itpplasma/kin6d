module milca_vrmhd_cylinder
    !! Helically symmetric straight periodic cylinder driven by a loop voltage:
    !! Ohmic profiles, polar section mesh, closures, sources and initial state.
    !! See docs/vrmhd-cylinder.md.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_gs_symmetry, only: gs_symmetry_t, GS_SCREW, GS_CHART_CARTESIAN
    use milca_vrmhd_parameters, only: vrmhd_parameters_t
    use milca_vrmhd_section, only: vrmhd_section_t, vrmhd_section_basis
    use milca_vrmhd_cylinder_leaves, only: vrmhd_cylinder_section, &
        vrmhd_cylinder_heat_shape, vrmhd_cylinder_source_shape
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use fortfem_triangle_vector_interpolation, only: &
        interpolate_triangle_nedelec_second_kind
    implicit none
    private
    real(dp), parameter :: pi = acos(-1.0_dp)

    !> Case of the JOREK cylinder recipe; units a = B0 = mu0 = rho = 1, time in
    !> Alfven times. chi_perpendicular is the toy-model diffusivity of
    !> d_t T = chi lap T; the section's internal energy 3 n T/2 rescales it.
    type, public :: vrmhd_cylinder_case_t
        real(dp) :: r0 = 10, q0 = 0.9_dp, beta_poloidal = 1, edge_temperature = 0.05_dp
        real(dp) :: heat_centre = 0.3_dp, heat_width = 0.1_dp
        real(dp) :: lundquist = 1.0e6_dp, prandtl = 10
        real(dp) :: chi_perpendicular = 2.4e-5_dp, chi_ratio = 1.0e6_dp
        !> S_j(0) in units of the source-free axis current; total current fixed.
        real(dp) :: source_amplitude = 0, source_width = 0.1_dp
        real(dp) :: ohmic_fraction = 0
        !> Isotropic particle diffusivity; JOREK's D_perp of these runs.
        real(dp) :: density_diffusion = 1.0e-6_dp
        integer :: grid = 20001
    end type

    !> Radial profiles on a uniform grid of x = r/a; quantities divided by r
    !> are stored with their axis limits.
    type, public :: vrmhd_cylinder_profile_t
        real(dp), allocatable :: x(:), psin(:), t(:), dt_over_r(:), heat(:)
        real(dp), allocatable :: jz(:), btheta_over_r(:), psi(:), bz(:)
        real(dp), allocatable :: atheta_over_r(:), jtheta_over_r(:), sj(:)
        !> Axial electric field E_z = eta (J_z - S_j), uniform.
        real(dp) :: electric_field = 0, temperature_axis = 0, heat_scale = 0
        real(dp) :: ohmic_factor = 1, q_axis = 0, axial_period = 0
    end type

    ! FortFEM's interpolation callback takes no context; a contained procedure
    ! would need an executable-stack trampoline, so the interpolated profile
    ! is held here for the duration of vrmhd_cylinder_state.
    type(vrmhd_cylinder_profile_t), save :: interpolated
    real(dp), save :: interpolated_pitch = 0

    integer, parameter, public :: CYL_PSIN = 1, CYL_T = 2, CYL_HEAT = 3, &
        CYL_JZ = 4, CYL_BTHETA_R = 5, CYL_PSI = 6, CYL_BZ = 7, CYL_ATHETA_R = 8, &
        CYL_JTHETA_R = 9, CYL_SJ = 10, CYL_DT_R = 11, CYL_NVALUE = 11

    public :: vrmhd_cylinder_profiles, vrmhd_cylinder_value
    public :: vrmhd_cylinder_mesh, vrmhd_cylinder_symmetry, vrmhd_cylinder_radii
    public :: vrmhd_cylinder_parameters, vrmhd_cylinder_sources
    public :: vrmhd_cylinder_state, vrmhd_cylinder_q, vrmhd_cylinder_diagnostics
contains
    subroutine vrmhd_cylinder_profiles(c, prof, status, current)
        !! Fixed point of make_case.py: uniform K, heat source shape(psi_N),
        !! j = c_ohm (T/T(0))^(3/2) + S_j, psi_N of the same current.
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(out) :: prof
        integer, intent(out) :: status
        !> Total current to keep with a source (int j r dr in shape units).
        real(dp), intent(in), optional :: current
        real(dp), allocatable :: x(:), psin(:), sh(:), g(:), gint(:), shape(:)
        real(dp), allocatable :: jshape(:), dpsi(:), psi(:), sg(:), work(:), bz2(:)
        real(dp) :: zj0, coefficient, dt, btheta_edge, pressure_average
        integer :: n, it, i
        status = 1
        n = c%grid
        if (n < 3 .or. c%r0 <= 0 .or. c%q0 <= 0 .or. c%lundquist <= 0) return
        if (c%source_amplitude > 0 .and. .not. present(current)) return
        allocate (x(n), psin(n), sh(n), g(n), gint(n), shape(n), jshape(n), &
            dpsi(n), psi(n), sg(n), work(n), bz2(n))
        x = [(real(i - 1, dp)/(n - 1), i=1, n)]
        zj0 = 2/c%q0
        psin = x**2
        sg = 0
        do it = 1, 200
            call vrmhd_cylinder_heat_shape(psin, c%heat_centre, c%heat_width, sh)
            call cumulative(sh*x, work)
            g = 0
            g(2:) = work(2:)/x(2:)
            call cumulative(g, gint)
            shape = c%edge_temperature + (gint(n) - gint)/gint(n)
            jshape = (shape/shape(1))**1.5_dp
            coefficient = 1
            if (c%source_amplitude > 0) then
                call vrmhd_cylinder_source_shape(psin, c%source_width, sg)
                sg = c%source_amplitude*sg
                call cumulative(sg*x, work)
                coefficient = current - work(n)
                call cumulative(jshape*x, work)
                coefficient = coefficient/work(n)
                jshape = coefficient*jshape + sg
            end if
            call cumulative(jshape*x, work)
            dpsi = 0
            dpsi(2:) = zj0*work(2:)/x(2:)
            call cumulative(dpsi, psi)
            work = psi/psi(n)
            if (maxval(abs(work - psin)) < 1.0e-12_dp) exit
            psin = work
        end do
        btheta_edge = dpsi(n)/c%r0
        call cumulative(shape*x, work)
        pressure_average = 2*work(n)
        dt = c%beta_poloidal*btheta_edge**2/2/pressure_average
        allocate (prof%x(n), source=x)
        allocate (prof%dt_over_r(n), prof%btheta_over_r(n), prof%atheta_over_r(n))
        prof%psin = psin
        prof%t = dt*shape
        ! T' = -dt g/gint(n); the axis limits use g ~ r sh(0)/2 and psi' ~ r j(0)/2.
        prof%dt_over_r = -dt*0.5_dp*sh(1)/gint(n)
        prof%dt_over_r(2:) = -dt*g(2:)/x(2:)/gint(n)
        prof%heat = dt*sh/gint(n)
        prof%jz = zj0*jshape/c%r0
        prof%btheta_over_r = zj0*jshape(1)/(2*c%r0)
        prof%btheta_over_r(2:) = dpsi(2:)/x(2:)/c%r0
        prof%psi = psi/c%r0
        prof%sj = zj0*sg/c%r0
        ! Radial force balance d(p + Bz^2/2)/dr = -B_theta J_z with p = T, Bz(a) = 1.
        work = prof%dt_over_r*x + prof%btheta_over_r*x*prof%jz
        call cumulative(work, bz2)
        prof%bz = sqrt(1 + 2*(bz2(n) - bz2))
        prof%jtheta_over_r = (prof%dt_over_r + prof%btheta_over_r*prof%jz)/prof%bz
        call cumulative(prof%bz*x, work)
        prof%atheta_over_r = prof%bz(1)/2
        prof%atheta_over_r(2:) = work(2:)/x(2:)**2
        prof%temperature_axis = prof%t(1)
        prof%ohmic_factor = coefficient
        prof%electric_field = zj0*coefficient/c%r0/c%lundquist
        prof%q_axis = prof%bz(1)/(c%r0*prof%btheta_over_r(1))
        prof%axial_period = 2*pi*c%r0
        status = 0
    contains
        subroutine cumulative(f, integral)
            real(dp), intent(in) :: f(:)
            real(dp), intent(out) :: integral(:)
            integer :: j
            integral(1) = 0
            do j = 2, size(f)
                integral(j) = integral(j - 1) &
                    + 0.5_dp*(f(j) + f(j - 1))*(x(j) - x(j - 1))
            end do
        end subroutine
    end subroutine

    pure subroutine vrmhd_cylinder_value(prof, r, values)
        !! Linear interpolation of every profile at radius r.
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(in) :: r
        real(dp), intent(out) :: values(CYL_NVALUE)
        real(dp) :: s, w
        integer :: i, n
        n = size(prof%x)
        s = min(max(r, 0.0_dp), 1.0_dp)*(n - 1)
        i = min(int(s) + 1, n - 1)
        w = s - (i - 1)
        values(CYL_PSIN) = lerp(prof%psin)
        values(CYL_T) = lerp(prof%t)
        values(CYL_HEAT) = lerp(prof%heat)
        values(CYL_JZ) = lerp(prof%jz)
        values(CYL_BTHETA_R) = lerp(prof%btheta_over_r)
        values(CYL_PSI) = lerp(prof%psi)
        values(CYL_BZ) = lerp(prof%bz)
        values(CYL_ATHETA_R) = lerp(prof%atheta_over_r)
        values(CYL_JTHETA_R) = lerp(prof%jtheta_over_r)
        values(CYL_SJ) = lerp(prof%sj)
        values(CYL_DT_R) = lerp(prof%dt_over_r)
    contains
        pure real(dp) function lerp(f)
            real(dp), intent(in) :: f(:)
            lerp = (1 - w)*f(i) + w*f(i + 1)
        end function
    end subroutine

    pure real(dp) function vrmhd_cylinder_q(prof, c, r) result(q)
        !! q = r Bz/(R0 B_theta) of the cylinder profiles.
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        type(vrmhd_cylinder_case_t), intent(in) :: c
        real(dp), intent(in) :: r
        real(dp) :: values(CYL_NVALUE)
        call vrmhd_cylinder_value(prof, r, values)
        q = values(CYL_BZ)/(c%r0*values(CYL_BTHETA_R))
    end function

    pure function vrmhd_cylinder_symmetry(c) result(symmetry)
        !! Fields depend on theta - z/R0: screw motion of pitch k = 1/R0.
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(gs_symmetry_t) :: symmetry
        symmetry = gs_symmetry_t(GS_SCREW, GS_CHART_CARTESIAN, 1/c%r0)
    end function

    subroutine vrmhd_cylinder_mesh(radii, inner, spacing, mm, status)
        !! Polar triangulation of the unit disk with rings at `radii`
        !! (increasing, last = 1). The ring point count starts at `inner` and
        !! doubles while the arc length exceeds `spacing` times the radial step,
        !! so cells stay anisotropic only up to that aspect ratio.
        real(dp), intent(in) :: radii(:), spacing
        integer, intent(in) :: inner
        type(meridian_mesh_t), intent(out) :: mm
        integer, intent(out) :: status
        type(mesh_2d_t) :: mesh
        integer, allocatable :: count(:), first(:), tri(:, :)
        real(dp), allocatable :: points(:, :)
        real(dp) :: step
        integer :: nring, i, j, nt, np, a, b, n0, n1
        status = 1
        nring = size(radii)
        if (nring < 1 .or. inner < 3 .or. spacing <= 0) return
        if (abs(radii(nring) - 1) > 1.0e-14_dp .or. radii(1) <= 0) return
        if (any(radii(2:) <= radii(:nring - 1))) return
        allocate (count(0:nring), first(0:nring))
        count(0) = 1; count(1) = inner
        do i = 2, nring
            step = radii(i) - radii(i - 1)
            count(i) = count(i - 1)
            if (2*pi*radii(i)/count(i - 1) > spacing*step) count(i) = 2*count(i - 1)
        end do
        first(0) = 1
        do i = 1, nring
            first(i) = first(i - 1) + count(i - 1)
        end do
        np = first(nring) + count(nring) - 1
        allocate (points(2, np), tri(3, 3*sum(count(1:))))
        points(:, 1) = 0
        do i = 1, nring
            do j = 0, count(i) - 1
                points(:, first(i) + j) = radii(i)*[cos(2*pi*j/count(i)), &
                    sin(2*pi*j/count(i))]
            end do
        end do
        nt = 0
        do j = 0, inner - 1
            call add(1, first(1) + j, first(1) + modulo(j + 1, inner))
        end do
        do i = 2, nring
            n0 = count(i - 1); n1 = count(i)
            do j = 0, n0 - 1
                a = first(i - 1) + j; b = first(i - 1) + modulo(j + 1, n0)
                if (n1 == n0) then
                    call add(a, first(i) + j, first(i) + modulo(j + 1, n1))
                    call add(a, first(i) + modulo(j + 1, n1), b)
                else
                    call add(a, first(i) + 2*j, first(i) + 2*j + 1)
                    call add(a, first(i) + 2*j + 1, b)
                    call add(b, first(i) + 2*j + 1, first(i) + modulo(2*j + 2, n1))
                end if
            end do
        end do
        mesh%n_vertices = np; mesh%vertices = points
        mesh%n_triangles = nt; mesh%triangles = tri(:, :nt)
        mesh%has_triangles = .true.
        call meridian_mesh_straight(mesh, 1, mm, status)
    contains
        subroutine add(p1, p2, p3)
            integer, intent(in) :: p1, p2, p3
            real(dp) :: e1(2), e2(2)
            nt = nt + 1
            tri(:, nt) = [p1, p2, p3]
            e1 = points(:, p2) - points(:, p1); e2 = points(:, p3) - points(:, p1)
            if (e1(1)*e2(2) - e1(2)*e2(1) < 0) tri(:, nt) = [p1, p3, p2]
        end subroutine
    end subroutine

    pure function vrmhd_cylinder_radii(n, centre, strength, uniform) result(radii)
        !! n ring radii in (0, 1], clustered at `centre` by
        !! r = (1-uniform)(centre + c sinh(strength (s - s0))) + uniform s,
        !! s in [0, 1], with c and s0 such that r(0) = 0 and r(1) = 1.
        integer, intent(in) :: n
        real(dp), intent(in) :: centre, strength, uniform
        real(dp) :: radii(n), s, s0, c, lo, hi, f
        integer :: i, it
        lo = 0; hi = 1
        ! sinh(b (1 - s0)) centre = sinh(b s0) (1 - centre) fixes s0 by bisection.
        do it = 1, 200
            s0 = (lo + hi)/2
            f = sinh(strength*(1 - s0))*centre - sinh(strength*s0)*(1 - centre)
            if (f > 0) then
                lo = s0
            else
                hi = s0
            end if
        end do
        c = (1 - centre)/sinh(strength*(1 - s0))
        do i = 1, n
            s = real(i, dp)/n
            radii(i) = (1 - uniform)*(centre + c*sinh(strength*(s - s0))) + uniform*s
        end do
        radii(n) = 1
    end function

    subroutine vrmhd_cylinder_parameters(c, prof, p)
        !! Closures of the JOREK cylinder runs: eta ~ Te^(-3/2) with
        !! eta(T(0)) = 1/S, constant viscosity Pm/S and constant conductivities.
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        type(vrmhd_parameters_t), intent(out) :: p
        p%mi = 1; p%mu0 = 1
        p%t0 = prof%temperature_axis/2
        p%eta0 = 1/c%lundquist; p%eta_exponent = -1.5_dp
        p%nu0 = c%prandtl/c%lundquist; p%nu_exponent = 0
        p%chi_perpendicular = 1.5_dp*c%chi_perpendicular
        p%chi_parallel0 = c%chi_ratio*p%chi_perpendicular
        p%chi_parallel_exponent = 0
        p%ohmic_fraction = c%ohmic_fraction
        p%d_perpendicular = c%density_diffusion; p%d_parallel = c%density_diffusion
    end subroutine

    subroutine vrmhd_cylinder_sources(op, c, prof, source, status)
        !! Js = J_theta e_theta + S_j e_z keeps the poloidal current of the
        !! force-balanced state and the central source; Qe, Qi keep T(r).
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_cylinder_case_t), intent(in) :: c
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(out) :: source(6, op%fs%nq, op%fs%mm%topology%n_triangles)
        integer, intent(out) :: status
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(CYL_NVALUE), heat
        integer :: t, q
        source = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q, x, weight, basis, status)
                if (status /= 0) return
                call vrmhd_cylinder_value(prof, norm2(x), v)
                source(1:3, q, t) = [-x(2)*v(CYL_JTHETA_R), x(1)*v(CYL_JTHETA_R), &
                    v(CYL_SJ)]
                heat = 1.5_dp*c%chi_perpendicular*v(CYL_HEAT)/2
                source(5, q, t) = heat - c%ohmic_fraction*prof%electric_field &
                    *v(CYL_JZ)
                source(6, q, t) = heat
            end do
        end do
        status = 0
    end subroutine

    subroutine vrmhd_cylinder_state(op, prof, u, status)
        !! Interpolant of the force-balanced Ohmic state: rho = 1, v = 0,
        !! Te = Ti = T/2, A and J from their defining moments and orbit nodes.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_cylinder_profile_t), intent(in) :: prof
        real(dp), intent(out) :: u(op%ndof)
        integer, intent(out) :: status
        real(dp), allocatable :: xy(:, :), moments(:)
        real(dp) :: v(CYL_NVALUE), k, u1, u2, a, vertices(2, 3)
        integer :: i, t, field, nln
        status = 1
        if (op%symmetry%kind /= GS_SCREW .or. .not. op%coordinate_coframe) return
        k = op%symmetry%k
        interpolated = prof; interpolated_pitch = k
        u = 0
        call meridian_lagrange_nodes(op%scalar_fs, xy)
        do i = 1, size(xy, 2)
            call vrmhd_cylinder_value(prof, norm2(xy(:, i)), v)
            u(op%offset(1) + i) = 1
            u(op%offset(5) + i) = v(CYL_T)/2
            u(op%offset(6) + i) = v(CYL_T)/2
        end do
        call meridian_lagrange_nodes(op%fs, xy)
        nln = op%fs%nln
        do field = 7, 8
            do i = 1, size(xy, 2)
                call vrmhd_cylinder_value(prof, norm2(xy(:, i)), v)
                if (field == 7) then
                    call vrmhd_cylinder_section(xy(1, i), xy(2, i), k, &
                        v(CYL_ATHETA_R), -v(CYL_PSI), u1, u2, a)
                else
                    call vrmhd_cylinder_section(xy(1, i), xy(2, i), k, &
                        v(CYL_JTHETA_R), v(CYL_JZ), u1, u2, a)
                end if
                u(op%offset(field) + op%fs%nn + i) = a
            end do
            do t = 1, op%fs%mm%topology%n_triangles
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
                u(op%offset(field) + op%fs%sp(1)%dof(:nln, t)) = &
                    moments*op%fs%sp(1)%sgn(:nln, t)
            end do
        end do
        status = 0
    end subroutine

    subroutine section_potential(px, py, w)
        real(dp), intent(in) :: px, py
        real(dp), intent(out) :: w(2)
        real(dp) :: pv(CYL_NVALUE), pa
        call vrmhd_cylinder_value(interpolated, hypot(px, py), pv)
        call vrmhd_cylinder_section(px, py, interpolated_pitch, pv(CYL_ATHETA_R), &
            -pv(CYL_PSI), w(1), w(2), pa)
    end subroutine

    subroutine section_current(px, py, w)
        real(dp), intent(in) :: px, py
        real(dp), intent(out) :: w(2)
        real(dp) :: pv(CYL_NVALUE), pa
        call vrmhd_cylinder_value(interpolated, hypot(px, py), pv)
        call vrmhd_cylinder_section(px, py, interpolated_pitch, pv(CYL_JTHETA_R), &
            pv(CYL_JZ), w(1), w(2), pa)
    end subroutine

    subroutine vrmhd_cylinder_diagnostics(op, r0, u, edges, q, energy, status)
        !! Annulus averages over the section quadrature: q = r Bz/(R0 B_theta)
        !! of the averaged field with Bz = 1 (the reduced-model convention) and
        !! with the averaged Bz in the annuli between `edges` (edges > 0), and
        !! the magnetic and kinetic energy per unit length of the m = 1 part,
        !! sum_annuli |integral f exp(-i theta) dA|^2/area over 400 equal annuli
        !! of the polar components (density taken as one in the kinetic part).
        type(vrmhd_section_t), intent(in) :: op
        real(dp), intent(in) :: r0, u(:), edges(:)
        real(dp), intent(out) :: q(2, size(edges) - 1), energy(2)
        integer, intent(out) :: status
        integer, parameter :: nfine = 400
        real(dp) :: basis(36, op%nlocal), x(2), weight, v(36), radius, polar(6)
        real(dp) :: area(size(edges) - 1), mean(7, size(edges) - 1)
        real(dp) :: fine_area(nfine)
        complex(dp) :: harmonic(6, nfine), phase
        integer :: t, q_point, bin, nbin
        nbin = size(edges) - 1
        area = 0; mean = 0; fine_area = 0; harmonic = 0
        do t = 1, op%fs%mm%topology%n_triangles
            do q_point = 1, op%fs%nq
                call vrmhd_section_basis(op, t, q_point, x, weight, basis, status)
                if (status /= 0) return
                v = matmul(basis, u(op%map(:, t)))
                call polar_components(x, v, polar)
                radius = norm2(x)
                phase = cmplx(x(1), -x(2), dp)/max(radius, tiny(1.0_dp))
                bin = min(nfine, int(radius*nfine) + 1)
                fine_area(bin) = fine_area(bin) + weight
                harmonic(:, bin) = harmonic(:, bin) + weight*polar*phase
                do bin = 1, nbin
                    if (radius < edges(bin) .or. radius >= edges(bin + 1)) cycle
                    area(bin) = area(bin) + weight
                    mean(:6, bin) = mean(:6, bin) + weight*polar
                    mean(7, bin) = mean(7, bin) + weight*polar(2)/radius
                end do
            end do
        end do
        status = 1
        if (any(area <= 0)) return
        energy = 0
        do bin = 1, nfine
            if (fine_area(bin) <= 0) cycle
            energy = energy + [sum(abs(harmonic(1:3, bin))**2), &
                sum(abs(harmonic(4:6, bin))**2)]/fine_area(bin)
        end do
        do bin = 1, nbin
            ! B_theta/r is regular on the axis, so its annulus mean has no
            ! bias from the area weighting of a linear B_theta.
            mean(:, bin) = mean(:, bin)/area(bin)
            q(1, bin) = 1/(r0*mean(7, bin))
            q(2, bin) = mean(3, bin)/(r0*mean(7, bin))
        end do
        status = 0
    contains
        pure subroutine polar_components(x, v, polar)
            real(dp), intent(in) :: x(2), v(36)
            real(dp), intent(out) :: polar(6)
            real(dp) :: rr, er(2), et(2)
            rr = max(norm2(x), tiny(1.0_dp))
            er = x/rr; et = [-er(2), er(1)]
            polar = [dot_product(er, v(25:26)), dot_product(et, v(25:26)), v(27), &
                dot_product(er, v(2:3)), dot_product(et, v(2:3)), v(4)]
        end subroutine
    end subroutine
end module milca_vrmhd_cylinder
