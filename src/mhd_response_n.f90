module kin6d_mhd_response_n
    !! L4 (first stage): flux harmonics on the rings of a flux-aligned mesh,
    !! resonant sheet currents and ideal shielding of the per-n field.
    !!
    !! Angles.  On a ring (flux surface) with geometric angle theta_g, field
    !! line advance nu(theta_g), straight-field-line angle theta_s(theta_g)
    !! and safety factor q (conventions of kin6d_vacuum_surfaces), the
    !! straight toroidal angle is zeta = phi - nu + q theta_s and a per-n
    !! field X e^{i n phi} has the harmonics e^{i (n zeta - m theta_s)}.
    !!
    !! Flux harmonic.  For the field B e^{i n phi}
    !!
    !!     Phi_m = (2/2 pi) int (B . n) R ds e^{i n (nu - q theta_s)} e^{i m theta_s},
    !!
    !! the peak-amplitude flux coefficient phi_mn of kin6d_vacuum_surfaces
    !! (T m^2 per rad^2; (B . n) R ds dphi is the flux element).  It is the
    !! same for every straight-field-line angle when m = n q, so the resonant
    !! Phi_m is the coordinate-independent resonant flux.  On a ring edge,
    !! (B . n) R ds = q . n ds = qhat . nhat dtau with the reference RT values
    !! (contravariant Piola), so Phi_m of a FEEC field is exact up to the
    !! angle interpolation.
    !!
    !! Sheet current.  A surface current J = delta_s nhat x grad_n Lambda on
    !! the surface s (delta_s the surface delta, nhat the outward unit
    !! normal) is divergence free for every stream function Lambda.  The
    !! resonant sheet of ring s carries
    !!
    !!     Lambda = lambda_s e^{i n (q theta_s - nu)} e^{-i m theta_s}
    !!
    !! (lambda_s in A; the physical stream function is Lambda e^{i n phi} +
    !! c.c.).  Its load in the gauge a = 0 of kin6d_mhd_ampere_n is
    !!
    !!     (A', J) = i n int (A'_pol . t) Lambda dl,
    !!
    !! t = nhat x e_phi the counter-clockwise unit tangent: only the e_phi
    !! part i n Lambda/R of grad_n Lambda drives A_pol.  The tangential trace
    !! A' . t dl = Nhat . that dtau is geometry independent (covariant Piola).
    !!
    !! Ideal shielding.  The sheet amplitudes lambda_s are fixed by the
    !! ideal constraint that the total resonant flux Phi_{m_s} vanishes on
    !! every rational ring s.  With the (vacuum) response matrix
    !! M(r, s) = Phi_{m_r}(ring r) of the unit sheet s, the amplitudes solve
    !! M lambda = -Phi_drive.  M is the mutual inductance of the sheets: with
    !! q . n = i n A . t in the gauge a = 0, M = -(1/(pi mu0)) F^H K^{-1} F
    !! for the unit-sheet loads F and the Ampere operator K, so M is
    !! Hermitian negative definite (a discrete reciprocity identity).  The quantity GPEC reports as singflx (resonant flux of
    !! the shielding current through its own surface, per area) is
    !! -M(s, s) lambda_s / area_s with area_s = (1/2 pi) int R ds.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, &
        meridian_point_basis_t, meridian_point_basis
    use fortfem_meridian_magnetic_coupling, only: meridian_edge_point
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_ring_edge, &
        flux_mesh_geometry
    use kin6d_mhd_ampere_n, only: mu0
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    !> Angle data of one ring on a fine uniform theta_g grid.
    type, public :: ring_angles_t
        real(dp) :: q = 0.0_dp
        integer :: nfine = 0
        real(dp), allocatable :: ths(:)      !! theta_s - theta_g (periodic)
        real(dp), allocatable :: nup(:)      !! nu - q theta_g (periodic)
    end type ring_angles_t

    !> Gauss rule along ring edges.
    integer, parameter, public :: ring_gauss = 6

    !> Resonant sheets and the linear constraints that fix their amplitudes.
    !> Sheet s carries the harmonic m = sheet_m(s) on ring sheet_ring(s);
    !> constraint r reads sum_k w(r, k) Phi_{flux_m(k)}(flux_ring(k)) =
    !> target(r) for the total field.  L4 shielding: one sheet per rational
    !> ring with Phi_m = 0 there (w = identity); L6 (kin6d_mhd_closure_n):
    !> two layer-edge sheets per rational with vanishing large-solution
    !> coefficients on both sides.
    type, public :: sheet_closure_t
        integer :: ns = 0, nf = 0
        integer, allocatable :: sheet_ring(:), sheet_m(:)
        type(ring_angles_t), allocatable :: sheet_ra(:)
        integer, allocatable :: flux_ring(:), flux_m(:)
        type(ring_angles_t), allocatable :: flux_ra(:)
        real(dp), allocatable :: w(:, :)          !! (ns, nf)
        complex(dp), allocatable :: target(:)      !! (ns)
    end type sheet_closure_t

    public :: closure_shield_init
    public :: closure_fluxes
    public :: ring_angles_set
    public :: ring_angles_at
    public :: ring_points
    public :: ring_flux_harmonic
    public :: ring_flux_from_values
    public :: ring_sheet_load
    public :: ring_area

contains

    subroutine ring_angles_set(q, theta_s, nu, ra)
        !! From theta_s(theta_g) and nu(theta_g) sampled at theta_g = 2 pi (j -
        !! 1)/N (as kin6d_vacuum_surfaces surface_build provides them).
        real(dp), intent(in) :: q, theta_s(:), nu(:)
        type(ring_angles_t), intent(out) :: ra
        integer :: j, nf
        real(dp) :: th

        nf = size(theta_s)
        ra%q = q
        ra%nfine = nf
        allocate (ra%ths(nf), ra%nup(nf))
        do j = 1, nf
            th = 2.0_dp*pi*(j - 1)/nf
            ra%ths(j) = theta_s(j) - th
            ra%nup(j) = nu(j) - q*th
        end do
    end subroutine ring_angles_set

    pure subroutine ring_angles_at(ra, theta, ths, nu)
        !! theta_s and nu at the geometric angle theta (periodic linear
        !! interpolation of the periodic parts on the fine grid).
        type(ring_angles_t), intent(in) :: ra
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: ths, nu
        real(dp) :: s, w
        integer :: i0, i1

        s = modulo(theta, 2.0_dp*pi)/(2.0_dp*pi)*ra%nfine
        i0 = int(s)
        w = s - i0
        i1 = modulo(i0 + 1, ra%nfine) + 1
        i0 = modulo(i0, ra%nfine) + 1
        ths = theta + (1.0_dp - w)*ra%ths(i0) + w*ra%ths(i1)
        nu = ra%q*theta + (1.0_dp - w)*ra%nup(i0) + w*ra%nup(i1)
    end subroutine ring_angles_at

    subroutine ring_points(fm, j, x, normal_ds, theta, t, tau, w)
        !! Gauss points of the edges of ring j: physical points x(2, :),
        !! (dZ/dtau, -dR/dtau) = n ds/dtau (outward), geometric angles, the
        !! triangle, the edge parameter and the Gauss weight (in tau) of each
        !! point.  Points of edge i are (i - 1) ring_gauss + 1 .. i ring_gauss.
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j
        real(dp), allocatable, intent(out) :: x(:, :), normal_ds(:, :)
        real(dp), allocatable, intent(out) :: theta(:), tau(:), w(:)
        integer, allocatable, intent(out) :: t(:)
        real(dp) :: gx(ring_gauss), gw(ring_gauss), xi, eta, nh(2), jac(2, 2)
        real(dp) :: ta, tb
        integer :: i, g, k, ne, te

        ne = fm%nth(j)
        allocate (x(2, ne*ring_gauss), normal_ds(2, ne*ring_gauss), &
            theta(ne*ring_gauss), t(ne*ring_gauss), tau(ne*ring_gauss), &
            w(ne*ring_gauss))
        call ring_gauss_rule(gx, gw)
        do i = 1, ne
            call flux_mesh_ring_edge(fm, j, i, te, ta, tb)
            do g = 1, ring_gauss
                k = (i - 1)*ring_gauss + g
                call meridian_edge_point(1, gx(g), xi, eta, nh)
                call flux_mesh_geometry(fm, te, xi, eta, x(:, k), jac)
                normal_ds(:, k) = [jac(2, 1), -jac(1, 1)]
                theta(k) = ta + gx(g)*(tb - ta)
                t(k) = te
                tau(k) = gx(g)
                w(k) = gw(g)
            end do
        end do
    end subroutine ring_points

    subroutine ring_gauss_rule(x, w)
        real(dp), intent(out) :: x(ring_gauss), w(ring_gauss)
        real(dp) :: z, p1, p2, p3, pp, z1
        integer :: i, k, n

        n = ring_gauss
        do i = 1, (n + 1)/2
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do
                p1 = 1.0_dp
                p2 = 0.0_dp
                do k = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2*k - 1)*z*p2 - (k - 1)*p3)/k
                end do
                pp = n*(z*p1 - p2)/(z*z - 1.0_dp)
                z1 = z
                z = z1 - p1/pp
                if (abs(z - z1) < 1.0e-15_dp) exit
            end do
            x(i) = 0.5_dp*(1.0_dp - z)
            x(n + 1 - i) = 0.5_dp*(1.0_dp + z)
            w(i) = 1.0_dp/((1.0_dp - z*z)*pp*pp)
            w(n + 1 - i) = w(i)
        end do
    end subroutine ring_gauss_rule

    function ring_flux_harmonic(fs, pb, fm, j, ra, n, b2, m) result(phim)
        !! Phi_m of the FEEC field with V2 coefficients b2 on ring j.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j, n, m
        type(ring_angles_t), intent(in) :: ra
        complex(dp), intent(in) :: b2(:)
        complex(dp) :: phim
        real(dp) :: gx(ring_gauss), gw(ring_gauss), xi, eta, nh(2), x(2)
        real(dp) :: ta, tb, th, ths, nu
        real(dp), allocatable :: phi(:, :), ref(:, :)
        complex(dp) :: flux
        integer :: i, g, l, te, st

        call ring_gauss_rule(gx, gw)
        allocate (phi(3, fs%sp(2)%nloc), ref(2, fs%nln))
        phim = (0.0_dp, 0.0_dp)
        do i = 1, fm%nth(j)
            call flux_mesh_ring_edge(fm, j, i, te, ta, tb)
            do g = 1, ring_gauss
                call meridian_edge_point(1, gx(g), xi, eta, nh)
                call meridian_point_basis(fs, pb, 2, te, xi, eta, x, phi, st, ref)
                flux = (0.0_dp, 0.0_dp)
                do l = 1, fs%nln
                    flux = flux + real(fs%sp(2)%sgn(l, te), dp)* &
                        b2(fs%sp(2)%dof(l, te))*(ref(1, l)*nh(1) + ref(2, l)*nh(2))
                end do
                th = ta + gx(g)*(tb - ta)
                call ring_angles_at(ra, th, ths, nu)
                phim = phim + gw(g)*flux*exp(cmplx(0.0_dp, &
                    n*(nu - ra%q*ths) + m*ths, dp))
            end do
        end do
        phim = phim/pi
    end function ring_flux_harmonic

    function ring_flux_from_values(ra, n, m, theta, bn_rds, w) result(phim)
        !! Phi_m from samples bn_rds = (B . n) R ds/dtau at the ring points
        !! with angles theta and quadrature weights w (in tau).
        type(ring_angles_t), intent(in) :: ra
        integer, intent(in) :: n, m
        real(dp), intent(in) :: theta(:), w(:)
        complex(dp), intent(in) :: bn_rds(:)
        complex(dp) :: phim
        real(dp) :: ths, nu
        integer :: k

        phim = (0.0_dp, 0.0_dp)
        do k = 1, size(theta)
            call ring_angles_at(ra, theta(k), ths, nu)
            phim = phim + w(k)*bn_rds(k)*exp(cmplx(0.0_dp, &
                n*(nu - ra%q*ths) + m*ths, dp))
        end do
        phim = phim/pi
    end function ring_flux_from_values

    subroutine ring_sheet_load(fs, pb, fm, j, ra, n, m, lambda, f)
        !! f += mu0 (A', J) of the sheet with amplitude lambda on ring j.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j, n, m
        type(ring_angles_t), intent(in) :: ra
        complex(dp), intent(in) :: lambda
        complex(dp), intent(inout) :: f(:)
        real(dp) :: gx(ring_gauss), gw(ring_gauss), xi, eta, nh(2), x(2)
        real(dp) :: ta, tb, th, ths, nu
        real(dp), allocatable :: phi(:, :), ref(:, :)
        complex(dp) :: lam
        integer :: i, g, l, te, st

        call ring_gauss_rule(gx, gw)
        allocate (phi(3, fs%sp(1)%nloc), ref(2, fs%nln))
        do i = 1, fm%nth(j)
            call flux_mesh_ring_edge(fm, j, i, te, ta, tb)
            do g = 1, ring_gauss
                call meridian_edge_point(1, gx(g), xi, eta, nh)
                call meridian_point_basis(fs, pb, 1, te, xi, eta, x, phi, st, ref)
                th = ta + gx(g)*(tb - ta)
                call ring_angles_at(ra, th, ths, nu)
                lam = lambda*exp(cmplx(0.0_dp, n*(ra%q*ths - nu) - m*ths, dp))
                do l = 1, fs%nln
                    ! reference tangent of local edge 1 is (1, 0)
                    f(fs%sp(1)%dof(l, te)) = f(fs%sp(1)%dof(l, te)) &
                        + mu0*cmplx(0.0_dp, n, dp)*gw(g)*lam &
                        *real(fs%sp(1)%sgn(l, te), dp)*ref(1, l)
                end do
            end do
        end do
    end subroutine ring_sheet_load

    subroutine closure_shield_init(ring, ra, m, cl)
        !! L4 shielding: a sheet on each given ring, zero resonant flux there.
        integer, intent(in) :: ring(:), m(:)
        type(ring_angles_t), intent(in) :: ra(:)
        type(sheet_closure_t), intent(out) :: cl
        integer :: s

        cl%ns = size(ring)
        cl%nf = size(ring)
        cl%sheet_ring = ring
        cl%sheet_m = m
        cl%sheet_ra = ra
        cl%flux_ring = ring
        cl%flux_m = m
        cl%flux_ra = ra
        allocate (cl%w(cl%ns, cl%nf), cl%target(cl%ns))
        cl%w = 0.0_dp
        do s = 1, cl%ns
            cl%w(s, s) = 1.0_dp
        end do
        cl%target = (0.0_dp, 0.0_dp)
    end subroutine closure_shield_init

    function closure_fluxes(fs, pb, fm, cl, n, b2) result(phi)
        !! Resonant fluxes Phi_{flux_m(k)}(flux_ring(k)) of the FEEC field b2.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        type(sheet_closure_t), intent(in) :: cl
        integer, intent(in) :: n
        complex(dp), intent(in) :: b2(:)
        complex(dp) :: phi(cl%nf)
        integer :: k

        !$omp parallel do schedule(dynamic) if(cl%nf > 1)
        do k = 1, cl%nf
            phi(k) = ring_flux_harmonic(fs, pb, fm, cl%flux_ring(k), cl%flux_ra(k), &
                n, b2, cl%flux_m(k))
        end do
        !$omp end parallel do
    end function closure_fluxes

    function ring_area(fm, j) result(area)
        !! (1/2 pi) int R ds over ring j (the surface area per rad^2).
        type(flux_mesh_t), intent(in) :: fm
        integer, intent(in) :: j
        real(dp) :: area
        real(dp), allocatable :: x(:, :), nds(:, :), theta(:), tau(:), w(:)
        integer, allocatable :: t(:)
        integer :: k

        call ring_points(fm, j, x, nds, theta, t, tau, w)
        area = 0.0_dp
        do k = 1, size(theta)
            area = area + w(k)*x(1, k)*norm2(nds(:, k))
        end do
        area = area/(2.0_dp*pi)
    end function ring_area

end module kin6d_mhd_response_n
