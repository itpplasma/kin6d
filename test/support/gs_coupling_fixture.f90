module gs_coupling_fixture
    !! Shared geometry, coils, current density and independent oracles for the
    !! FEM-BEM coupling tests (DESIGN section 5.1, T1 and T2).
    !!
    !! Gamma_c is an ellipse (centre R = 3, semi-axes 1.25 x 1.9); Gamma_p is
    !! the ellipse (centre 3.05, semi-axes 0.65 x 1.0) that bounds the compact
    !! prescribed current density j(rho) = j0 (1 - rho^2)^3, rho the elliptic
    !! radius about the centre of Gamma_p. j is C^2 across Gamma_p and carries
    !! a total current of 1 MA. Four rectangular coils (two circuits, one with
    !! a reversed winding) sit at least 0.3 m outside Gamma_c.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_axisymmetric_gs_kernel, only: gs_green, gs_mu0, gs_rect_coil_t
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, gs_element_jet
    use kin6d_gs_coils, only: gs_coil_t, gs_circuit_set_t, gs_circuit_flux_vectors
    implicit none
    private

    real(dp), parameter, public :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: gc_r0 = 3.0_dp, gc_a = 1.25_dp, gc_b = 1.9_dp
    real(dp), parameter, public :: gp_r0 = 3.05_dp, gp_a = 0.65_dp, gp_b = 1.0_dp
    real(dp), parameter, public :: plasma_current = 1.0e6_dp
    real(dp), parameter, public :: j0 = 4.0_dp * plasma_current / (pi * gp_a * gp_b)

    public :: outer_curve, build_state, build_coils, density, plasma_flux_ref
    public :: plasma_energy_ref, fe_value, elliptic_radius2

contains

    subroutine outer_curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)

        point = [gc_r0 + gc_a * cos(2.0_dp * pi * t), gc_b * sin(2.0_dp * pi * t)]
    end subroutine outer_curve

    subroutine build_state(state, n_outer, order, status)
        !! Two-loop mesh: Gamma_c with n_outer edges, the Gamma_p chain with
        !! n_outer/2 edges, curved elements of the given order.
        type(gs_state_t), intent(out) :: state
        integer, intent(in) :: n_outer, order
        integer, intent(out) :: status
        real(dp), allocatable :: contour(:, :), chain(:, :)
        real(dp) :: h, perimeter
        integer :: j, nc

        nc = n_outer / 2
        allocate (contour(2, n_outer), chain(2, nc))
        do j = 1, n_outer
            call outer_curve(real(j - 1, dp) / n_outer, contour(:, j))
        end do
        do j = 1, nc
            chain(:, j) = [gp_r0 + gp_a * cos(2.0_dp * pi * (j - 1) / nc), &
                           gp_b * sin(2.0_dp * pi * (j - 1) / nc)]
        end do
 perimeter = pi*(3.0_dp*(gc_a + gc_b) - sqrt((3.0_dp*gc_a + gc_b)*(gc_a + 3.0_dp*gc_b)))
        h = perimeter / n_outer
       call initialize_gs(state, contour, 0.4_dp * h * h, status, element_order=order, &
                           curved=.true., curve=outer_curve, interior_chain=chain)
    end subroutine build_state

    subroutine build_coils(set)
        !! Circuit 1: coils 1+2 in series; circuit 2: coil 3 minus coil 4.
        type(gs_circuit_set_t), intent(out) :: set
      real(dp), parameter :: centres(2, 4) = reshape([4.7_dp, 1.9_dp, 4.8_dp, -1.8_dp, &
                                               1.2_dp, 2.2_dp, 1.3_dp, -2.1_dp], [2, 4])
        integer :: c

        allocate (set % coils(4), set % incidence(4, 2))
        set % incidence = 0.0_dp
        set % incidence(1:2, 1) = 1.0_dp
        set % incidence(3, 2) = 1.0_dp
        set % incidence(4, 2) = -1.0_dp
        do c = 1, 4
          set % coils(c) % rect = gs_rect_coil_t(centres(1, c), centres(2, c), 0.3_dp, &
                                                 0.3_dp, merge(10.0_dp, 5.0_dp, c <= 2))
        end do
    end subroutine build_coils

    pure function elliptic_radius2(r, z) result(rho2)
        real(dp), intent(in) :: r, z
        real(dp) :: rho2

        rho2 = ((r - gp_r0) / gp_a)**2 + (z / gp_b)**2
    end function elliptic_radius2

    function density(r, z) result(j)
        !! Prescribed toroidal current density (A/m^2), independent of psi.
        real(dp), intent(in) :: r, z
        real(dp) :: j, rho2

        rho2 = elliptic_radius2(r, z)
        j = 0.0_dp
        if (rho2 < 1.0_dp) j = j0 * (1.0_dp - rho2)**3
    end function density

    function plasma_flux_ref(r, z, level) result(psi)
        !! psi_p(r, z) = mu0 int G j dA by tensor cubature. Inside the support
        !! the cubature is polar about the evaluation point (the logarithmic
        !! singularity of G is absorbed by the area element and a quadratic
        !! radial grading); outside, elliptic coordinates about Gamma_p with a
        !! periodic trapezoid rule in the angle. `level` scales the point counts.
        real(dp), intent(in) :: r, z
        integer, intent(in) :: level
        real(dp) :: psi
        real(dp), allocatable :: u(:), wu(:)
        real(dp) :: phi, c, s, a_, b_, c_, smax, sum_phi, sum_s, yr, yz, rho
        integer :: nu, nphi, i, k

        psi = 0.0_dp
        if (elliptic_radius2(r, z) < 1.0_dp) then
            nu = 40 * level
            nphi = 64 * level
            allocate (u(nu), wu(nu))
            call gauss_legendre(nu, u, wu)
            u = 0.5_dp * (u + 1.0_dp)
            wu = 0.5_dp * wu
            sum_phi = 0.0_dp
            do k = 1, nphi
                phi = 2.0_dp * pi * (k - 0.5_dp) / nphi
                c = cos(phi)
                s = sin(phi)
                a_ = (c / gp_a)**2 + (s / gp_b)**2
                b_ = (r - gp_r0) * c / gp_a**2 + z * s / gp_b**2
                c_ = elliptic_radius2(r, z) - 1.0_dp
                smax = (-b_ + sqrt(b_ * b_ - a_ * c_)) / a_
                sum_s = 0.0_dp
                do i = 1, nu
                    yr = r + smax * u(i)**2 * c
                    yz = z + smax * u(i)**2 * s
                    sum_s = sum_s + wu(i) * density(yr, yz) * gs_green(r, z, yr, yz) &
                            * smax * u(i)**2 * 2.0_dp * smax * u(i)
                end do
                sum_phi = sum_phi + sum_s
            end do
            psi = gs_mu0 * sum_phi * 2.0_dp * pi / nphi
        else
            nu = 24 * level
            nphi = 96 * level
            allocate (u(nu), wu(nu))
            call gauss_legendre(nu, u, wu)
            u = 0.5_dp * (u + 1.0_dp)
            wu = 0.5_dp * wu
            sum_phi = 0.0_dp
            do k = 1, nphi
                phi = 2.0_dp * pi * (k - 0.5_dp) / nphi
                do i = 1, nu
                    rho = sqrt(u(i))
                    yr = gp_r0 + gp_a * rho * cos(phi)
                    yz = gp_b * rho * sin(phi)
                    ! dA = a b rho drho dphi = a b du dphi / 2 with u = rho^2.
                    sum_phi = sum_phi + wu(i) * density(yr, yz) * gs_green(r, z, yr, yz)
                end do
            end do
            psi = gs_mu0 * sum_phi * 2.0_dp * pi / nphi * 0.5_dp * gp_a * gp_b
        end if
    end function plasma_flux_ref

    function plasma_energy_ref(level) result(energy)
        !! mu0 int j psi_p dA, the whole-plane energy of the plasma flux
        !! (Green's identity), on an elliptic Gauss grid with psi_p from the
        !! cubature above.
        integer, intent(in) :: level
        real(dp) :: energy
        real(dp), allocatable :: u(:), wu(:)
        real(dp) :: phi, rho, yr, yz
        integer :: nu, nphi, i, k

        nu = 12 * level
        nphi = 24 * level
        allocate (u(nu), wu(nu))
        call gauss_legendre(nu, u, wu)
        u = 0.5_dp * (u + 1.0_dp)
        wu = 0.5_dp * wu
        energy = 0.0_dp
        do k = 1, nphi
            phi = 2.0_dp * pi * (k - 0.5_dp) / nphi
            do i = 1, nu
                rho = sqrt(u(i))
                yr = gp_r0 + gp_a * rho * cos(phi)
                yz = gp_b * rho * sin(phi)
              energy = energy + wu(i) * density(yr, yz) * plasma_flux_ref(yr, yz, level)
            end do
        end do
        energy = gs_mu0 * energy * 2.0_dp * pi / nphi * 0.5_dp * gp_a * gp_b
    end function plasma_energy_ref

    subroutine fe_value(state, psi, r, z, value, found)
        !! Value of the finite element field `psi` at (r, z).
        type(gs_state_t), intent(inout) :: state
        real(dp), intent(in) :: psi(:), r, z
        real(dp), intent(out) :: value
        logical, intent(out) :: found
        real(dp) :: grad(2), hessian(2, 2)
        integer :: e, status

        state % psi = psi
        found = .false.
        value = 0.0_dp
        do e = 1, size(state % triangles, 2)
            call gs_element_jet(state, e, r, z, value, grad, hessian, status)
            if (status == 0) then
                found = .true.
                return
            end if
        end do
    end subroutine fe_value

end module gs_coupling_fixture
