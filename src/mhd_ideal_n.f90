module kin6d_mhd_ideal_n
    !! L4 (stage 3): volume current of the linear ideal-MHD perturbed
    !! equilibrium (static, omega = 0) in the form used by MEPHIT (Albert,
    !! Lainer and Biro, Comput. Phys. Commun. 277 (2022) 108401; Markl et al.,
    !! the MEPHIT papers), on the half rings of a flux-aligned mesh.
    !!
    !! Model.  For the total per-n field B (coil field plus response) the
    !! plasma current J follows from linear force balance
    !!     grad p_n = J x B0 + J0 x B,   div J = 0:
    !!   1. B0 . grad p_n = -p0' B . grad psi                 (pressure),
    !!   2. J_perp = B0 x (grad p_n - J0 x B)/B0^2,
    !!   3. B0 . grad h = -div J_perp,  J_par = h B0,  with
    !!      -div J_perp = -grad p_n . curl(B0/B0^2) + J0 . grad(B0 . B/B0^2)
    !!                    - B . grad sigma0,   sigma0 = J0 . B0/B0^2,
    !!      (div J0 = div B = 0; J0 . grad psi = 0, so J0 . grad needs only
    !!      derivatives within the surface).
    !! The magnetic differential equations are solved exactly on a surface
    !! sampled at uniform geometric angles: with nu = q theta + nu_p and
    !! B0 . grad f = (B_t/|x'|) f' + i n (F/R^2) f, the substitution
    !! f = e^{-i n nu_p} phi gives phi' + i n q phi = (|x'|/B_t) r e^{i n nu_p},
    !! diagonal in the Fourier index k of theta_g: phi_k = rho_k/(i (k + n q)),
    !! singular only on rational surfaces (k + n q = 0), which are mesh rings
    !! and never half rings.  The resonant homogeneous solutions are the sheet
    !! currents of kin6d_mhd_response_n, fixed by ideal shielding.
    !!
    !! Equilibrium (kin6d_vacuum_surfaces conventions): B0_pol = (s/R)
    !! (psi_Z, -psi_R), B0_phi = F/R, s = sgn(Ip) sgn(psi_edge - psi_axis);
    !! from mu0 curl(F/R e_phi) and radial force balance
    !!     mu0 J0_pol = (F'/R)(-psi_Z, psi_R),  J0_phi = -s (R p0' + F F'/(mu0 R)),
    !! with p0' = dp/dpsi and F' = dF/dpsi from the G-EQDSK profiles, so that
    !! J0 x B0 = grad p0 holds exactly at every point.
    !!
    !! Discretisation.  Half ring h (strip h between rings h-1 and h, h >= 2)
    !! is the flux surface u_h = (u_{h-1} + u_h)/2 sampled at np uniform
    !! geometric angles.  The current of strip h is that of its half ring,
    !! constant along u and interpolated in angle (first order in the ring
    !! spacing, as the strip-wise current of MEPHIT).  The fan around the axis
    !! carries no plasma current.  d/dpsi at fixed theta_g uses differences
    !! between half rings.  The volume current that leaves the plasma through
    !! the boundary ring continues in a field-aligned boundary sheet K =
    !! kappa B0, B0 . grad kappa = J . n, so the total current is divergence
    !! free; the boundary pressure jump of a cut equilibrium (p0_edge > 0
    !! with vacuum outside) is not modelled.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, &
        meridian_point_basis_t, meridian_point_basis, meridian_qpoint
    use fortfem_meridian_magnetic_coupling, only: meridian_edge_point
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_eval, profiles_of_psin, &
        surface_t, surface_build
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_locate, flux_mesh_ring_edge
    use kin6d_mhd_ampere_n, only: mu0
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c
    implicit none
    private

    real(dp), parameter :: pi = 3.14159265358979323846_dp

    !> Equilibrium and FEEC sampling data on one sample surface.
    type, public :: sample_ring_t
        integer :: strip = 0, np = 0
        real(dp) :: u = 0.0_dp, psi = 0.0_dp, q = 0.0_dp
        real(dp), allocatable :: x(:, :)        !! (2, np) points
        integer, allocatable :: tri(:)          !! triangle of each point
        real(dp), allocatable :: phi(:, :, :)   !! (3, nloc2, np) V2 basis
        real(dp), allocatable :: b0(:, :), j0(:, :), c1(:, :)  !! (3, np)
        real(dp), allocatable :: gpsi(:, :), gth(:, :), gsig(:, :)  !! (2, np)
        real(dp), allocatable :: pp(:), fac(:), nup(:), b0sq(:)
        complex(dp), allocatable :: bc(:, :)    !! (3, np) external field
    end type sample_ring_t

    !> Override of the equilibrium current for the force-free oracle:
    !> J0 = sigma B0/mu0 with constant sigma and p0' = 0.
    type, public :: ideal_options_t
        logical :: force_free = .false.
        real(dp) :: sigma = 0.0_dp
        logical :: boundary_sheet = .true.
        !> No volume response on surfaces with |q| < qlow (the core inside the
        !> integration start of DCON, qlow there; 0 keeps the whole plasma).
        real(dp) :: qlow = 0.0_dp
    end type ideal_options_t

    public :: equilibrium_point
    public :: sample_ring_build
    public :: sample_field
    public :: mde_solve
    public :: ideal_current
    public :: ideal_load
    public :: ideal_response

contains

    subroutine equilibrium_point(g, opt, r, z, b0, j0, gpsi, pp)
        !! B0, J0 (R, phi, Z), grad psi (R, Z) and p0' = dp/dpsi at (R, Z).
        type(psi_grid_t), intent(in) :: g
        type(ideal_options_t), intent(in) :: opt
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: b0(3), j0(3), gpsi(2), pp
        real(dp) :: psi, pr, pz, x, f, dfdx, p, dpdx, dpsi, fp

        call psi_eval(g, r, z, psi, pr, pz)
        dpsi = g%psi_edge - g%psi_axis
        x = (psi - g%psi_axis)/dpsi
        call profiles_of_psin(g, x, f, dfdx, p, dpdx)
        fp = dfdx/dpsi
        pp = dpdx/dpsi
        gpsi = [pr, pz]
        b0 = [g%s_pol*pz/r, f/r, -g%s_pol*pr/r]
        if (opt%force_free) then
            j0 = opt%sigma*b0/mu0
            pp = 0.0_dp
        else
            j0(1) = -fp*pz/(mu0*r)
            j0(3) = fp*pr/(mu0*r)
            j0(2) = -g%s_pol*(r*pp + f*fp/(mu0*r))
        end if
    end subroutine equilibrium_point

    subroutine sample_ring_build(g, opt, fm, fs, pb, u, strip, np, sr, status)
        !! Sample surface u (psin = u^2) at np uniform geometric angles,
        !! located in strip `strip` of the mesh, with the equilibrium data and
        !! the FEEC V2 basis at every point.
        type(psi_grid_t), intent(in) :: g
        type(ideal_options_t), intent(in) :: opt
        type(flux_mesh_t), intent(in) :: fm
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        real(dp), intent(in) :: u
        integer, intent(in) :: strip, np
        type(sample_ring_t), intent(out) :: sr
        integer, intent(out) :: status
        type(surface_t) :: s
        real(dp) :: xi, eta, x(2), th, h, bp(3), bm(3), jp(3), jm(3), gp(2), pq
        real(dp) :: vp(3), vm(3), vzp(3), vzm(3), rho2, ulab
        real(dp) :: b0(3), j0(3)
        integer :: k, st, nl, bad

        status = 1
        call surface_build(g, u*u, np, s, st)
        if (st /= 0) return
        sr%strip = strip
        sr%np = np
        sr%u = u
        sr%psi = g%psi_axis + u*u*(g%psi_edge - g%psi_axis)
        sr%q = s%q
        nl = fs%sp(2)%nloc
        allocate (sr%x(2, np), sr%tri(np), sr%phi(3, nl, np), sr%b0(3, np), &
            sr%j0(3, np), sr%c1(3, np), sr%gpsi(2, np), sr%gth(2, np), &
            sr%gsig(2, np), sr%pp(np), sr%fac(np), sr%nup(np), sr%b0sq(np))
        ulab = u
        if (strip > 1) ulab = 0.5_dp*(fm%u(strip - 1) + fm%u(strip))
        bad = 0
        !$omp parallel do schedule(dynamic, 16) reduction(max:bad) &
        !$omp private(k, th, st, xi, eta, x, b0, j0, gp, pq, rho2, h, vp, vm, &
        !$omp vzp, vzm, bp, bm, jp, jm)
        do k = 1, np
            th = s%theta_g(k)
            sr%x(:, k) = [s%r(k), s%z(k)]
            call flux_mesh_locate(fm, sr%x(:, k), ulab, th, sr%tri(k), xi, eta, st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            call meridian_point_basis(fs, pb, 2, sr%tri(k), xi, eta, x, &
                sr%phi(:, :, k), st)
            if (st /= 0) then
                bad = 1
                cycle
            end if
            call equilibrium_point(g, opt, s%r(k), s%z(k), b0, j0, gp, pq)
            sr%b0(:, k) = b0
            sr%j0(:, k) = j0
            sr%gpsi(:, k) = gp
            sr%pp(k) = pq
            sr%b0sq(k) = sum(b0**2)
            ! |x'|/B_t, B_t = B0_pol . t = -s sgn(dpsi) |grad psi|/R for the
            ! counter-clockwise tangent t = (-n_Z, n_R)
            sr%fac(k) = -s%dl(k)*s%r(k)/(g%s_pol*g%sgn_dpsi*norm2(gp))
            sr%nup(k) = s%nu(k) - s%q*th
            rho2 = (s%r(k) - g%r_axis)**2 + (s%z(k) - g%z_axis)**2
            sr%gth(:, k) = [-(s%z(k) - g%z_axis), s%r(k) - g%r_axis]/rho2
            ! curl(B0/B0^2) and grad sigma0 by central differences
            h = 1.0e-5_dp*s%r(k)
            call field_ratio(s%r(k) + h, s%z(k), vp, bp, jp)
            call field_ratio(s%r(k) - h, s%z(k), vm, bm, jm)
            call field_ratio(s%r(k), s%z(k) + h, vzp, bp, jp)
            call field_ratio(s%r(k), s%z(k) - h, vzm, bm, jm)
            sr%c1(1, k) = -(vzp(2) - vzm(2))/(2*h)
            sr%c1(2, k) = (vzp(1) - vzm(1))/(2*h) - (vp(3) - vm(3))/(2*h)
            sr%c1(3, k) = ((s%r(k) + h)*vp(2) - (s%r(k) - h)*vm(2))/(2*h*s%r(k))
            sr%gsig(1, k) = (sigma_at(s%r(k) + h, s%z(k)) &
                - sigma_at(s%r(k) - h, s%z(k)))/(2*h)
            sr%gsig(2, k) = (sigma_at(s%r(k), s%z(k) + h) &
                - sigma_at(s%r(k), s%z(k) - h))/(2*h)
        end do
        !$omp end parallel do
        if (bad /= 0) return
        allocate (sr%bc(3, np))
        sr%bc = (0.0_dp, 0.0_dp)
        status = 0
    contains
        subroutine field_ratio(rr, zz, v, bb, jj)
            real(dp), intent(in) :: rr, zz
            real(dp), intent(out) :: v(3), bb(3), jj(3)
            real(dp) :: gg(2), pq2
            call equilibrium_point(g, opt, rr, zz, bb, jj, gg, pq2)
            v = bb/sum(bb**2)
        end subroutine field_ratio
        real(dp) function sigma_at(rr, zz)
            real(dp), intent(in) :: rr, zz
            real(dp) :: bb(3), jj(3), gg(2), pq2
            call equilibrium_point(g, opt, rr, zz, bb, jj, gg, pq2)
            sigma_at = dot_product(jj, bb)/sum(bb**2)
        end function sigma_at
    end subroutine sample_ring_build

    subroutine sample_field(fs, sr, b2, b)
        !! FEEC field (R, phi, Z) at the sample points from V2 coefficients.
        type(meridian_feec_t), intent(in) :: fs
        type(sample_ring_t), intent(in) :: sr
        complex(dp), intent(in) :: b2(:)
        complex(dp), intent(out) :: b(:, :)
        integer :: k, l, t

        do k = 1, sr%np
            t = sr%tri(k)
            b(:, k) = (0.0_dp, 0.0_dp)
            do l = 1, fs%sp(2)%nloc
                b(:, k) = b(:, k) + real(fs%sp(2)%sgn(l, t), dp)* &
                    b2(fs%sp(2)%dof(l, t))*sr%phi(:, l, k)
            end do
        end do
    end subroutine sample_field

    pure integer function wavenumber(i, np)
        !! Fourier index of entry i of an unshifted length-np DFT.
        integer, intent(in) :: i, np
        wavenumber = i - 1
        if (2*(i - 1) >= np) wavenumber = i - 1 - np
    end function wavenumber

    subroutine dtheta(f, df, plan)
        !! Spectral derivative in theta_g of periodic samples (the Nyquist
        !! mode of an even length is dropped).
        complex(dp), intent(in) :: f(:)
        complex(dp), intent(out) :: df(:)
        type(fortnum_fft_plan_t), intent(in) :: plan
        integer :: np, i, k

        np = size(f)
        df = f
        call fft_c2c(df, -1, plan)
        do i = 1, np
            k = wavenumber(i, np)
            if (2*k == -np) then
                df(i) = (0.0_dp, 0.0_dp)
            else
                df(i) = cmplx(0.0_dp, k, dp)*df(i)/np
            end if
        end do
        call fft_c2c(df, 1, plan)
    end subroutine dtheta

    subroutine mde_solve(sr, n, r, f, plan)
        !! Solve B0 . grad (f e^{i n phi}) = r e^{i n phi} on the sample
        !! surface (see the module header).
        type(sample_ring_t), intent(in) :: sr
        integer, intent(in) :: n
        complex(dp), intent(in) :: r(:)
        complex(dp), intent(out) :: f(:)
        type(fortnum_fft_plan_t), intent(in), optional :: plan
        type(fortnum_fft_plan_t) :: own
        integer :: i, np

        np = sr%np
        f = r*sr%fac*exp(cmplx(0.0_dp, n*sr%nup, dp))
        if (present(plan)) then
            call fft_c2c(f, -1, plan)
        else
            call fft_c2c_plan_init(own, np)
            call fft_c2c(f, -1, own)
        end if
        do i = 1, np
            f(i) = f(i)/(cmplx(0.0_dp, wavenumber(i, np) + n*sr%q, dp)*np)
        end do
        if (present(plan)) then
            call fft_c2c(f, 1, plan)
        else
            call fft_c2c(f, 1, own)
        end if
        f = f*exp(cmplx(0.0_dp, -n*sr%nup, dp))
    end subroutine mde_solve

    subroutine ideal_current(rings, n, sgn_dpsi, qlow, b, jcur, kap, bnd)
        !! Plasma current jcur(3, np, ring) of the total field b(3, np, ring) on
        !! the sample rings (half rings 2 .. nring and, last, the boundary ring
        !! bnd), and the boundary-sheet coefficient kap(np) (K = kap B0).
        type(sample_ring_t), intent(in) :: rings(:)
        integer, intent(in) :: n, bnd
        real(dp), intent(in) :: sgn_dpsi, qlow
        complex(dp), intent(in) :: b(:, :, :)
        complex(dp), intent(out) :: jcur(:, :, :), kap(:)
        complex(dp), allocatable :: p(:, :), dpth(:), dpsi(:), g(:), dg(:)
        complex(dp), allocatable :: rhs(:), h(:), jn(:)
        complex(dp) :: gradp(3), v(3), jxb(3), bb(3)
        type(fortnum_fft_plan_t) :: plan
        integer :: nr, ir, k, np, lo, hi
        real(dp) :: ds

        nr = size(rings)
        np = rings(1)%np
        call fft_c2c_plan_init(plan, np)
        allocate (p(np, nr))
        ! 1. pressure on every sample ring
        !$omp parallel do schedule(dynamic) private(ir, rhs) if(nr > 1)
        do ir = 1, nr
            allocate (rhs(np))
            rhs = -rings(ir)%pp*(b(1, :, ir)*rings(ir)%gpsi(1, :) &
                + b(3, :, ir)*rings(ir)%gpsi(2, :))
            call mde_solve(rings(ir), n, rhs, p(:, ir), plan)
            deallocate (rhs)
        end do
        !$omp end parallel do
        ! 2. perpendicular and parallel current
        !$omp parallel do schedule(dynamic) if(nr > 1) &
        !$omp private(ir, k, lo, hi, ds, dpth, dpsi, g, dg, rhs, h, gradp, v, &
        !$omp jxb, bb)
        do ir = 1, nr
            allocate (dpth(np), dpsi(np), g(np), dg(np), rhs(np), h(np))
            call dtheta(p(:, ir), dpth, plan)
            ! d/dpsi at fixed theta_g between neighbouring rings of the
            ! sequence (the boundary ring follows the last half ring)
            lo = max(1, ir - 1)
            hi = min(nr, ir + 1)
            ds = rings(hi)%psi - rings(lo)%psi
            dpsi = (p(:, hi) - p(:, lo))/ds
            do k = 1, np
                bb = b(:, k, ir)
                gradp(1) = dpsi(k)*rings(ir)%gpsi(1, k) + dpth(k)*rings(ir)%gth(1, k)
                gradp(3) = dpsi(k)*rings(ir)%gpsi(2, k) + dpth(k)*rings(ir)%gth(2, k)
                gradp(2) = cmplx(0.0_dp, n, dp)*p(k, ir)/rings(ir)%x(1, k)
                jxb = cross_rc(rings(ir)%j0(:, k), bb)
                v = gradp - jxb
                jcur(:, k, ir) = cross_rc(rings(ir)%b0(:, k), v)/rings(ir)%b0sq(k)
                g(k) = dot_product(rings(ir)%b0(:, k), bb)/rings(ir)%b0sq(k)
            end do
            call dtheta(g, dg, plan)
            do k = 1, np
                bb = b(:, k, ir)
                gradp(1) = dpsi(k)*rings(ir)%gpsi(1, k) + dpth(k)*rings(ir)%gth(1, k)
                gradp(3) = dpsi(k)*rings(ir)%gpsi(2, k) + dpth(k)*rings(ir)%gth(2, k)
                gradp(2) = cmplx(0.0_dp, n, dp)*p(k, ir)/rings(ir)%x(1, k)
                rhs(k) = -sum(gradp*rings(ir)%c1(:, k)) &
                    + (rings(ir)%j0(1, k)*rings(ir)%gth(1, k) &
                    + rings(ir)%j0(3, k)*rings(ir)%gth(2, k))*dg(k) &
                    + rings(ir)%j0(2, k)/rings(ir)%x(1, k)*cmplx(0.0_dp, n, dp)*g(k) &
                    - (bb(1)*rings(ir)%gsig(1, k) + bb(3)*rings(ir)%gsig(2, k))
            end do
            call mde_solve(rings(ir), n, rhs, h, plan)
            do k = 1, np
                jcur(:, k, ir) = jcur(:, k, ir) + h(k)*rings(ir)%b0(:, k)
            end do
            if (abs(rings(ir)%q) < qlow) jcur(:, :, ir) = (0.0_dp, 0.0_dp)
            deallocate (dpth, dpsi, g, dg, rhs, h)
        end do
        !$omp end parallel do
        ! 3. boundary sheet: B0 . grad kap = J . n on the boundary ring
        kap = (0.0_dp, 0.0_dp)
        if (bnd > 0) then
            allocate (jn(np))
            do k = 1, np
                jn(k) = sgn_dpsi*(jcur(1, k, bnd)*rings(bnd)%gpsi(1, k) &
                    + jcur(3, k, bnd)*rings(bnd)%gpsi(2, k)) &
                    /norm2(rings(bnd)%gpsi(:, k))
            end do
            call mde_solve(rings(bnd), n, jn, kap, plan)
        end if
    end subroutine ideal_current

    pure function cross_rc(a, b) result(c)
        real(dp), intent(in) :: a(3)
        complex(dp), intent(in) :: b(3)
        complex(dp) :: c(3)
        ! (R, phi, Z) right handed: e_R x e_phi = e_Z
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_rc

    subroutine ideal_load(fs, pb, fm, g, jcur, kap, bnd, f)
        !! f = mu0 (A', J) for the Nedelec test functions: volume current of
        !! strip h from sample ring h - 1 (half ring of strip h), interpolated
        !! linearly in angle; boundary sheet K = kap B0 on the last ring with
        !! (A', K) = int (A' . t) kap B0_t R dl, A' . t dl = Nhat . (1, 0) dtau.
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        type(psi_grid_t), intent(in) :: g
        complex(dp), intent(in) :: jcur(:, :, :), kap(:)
        integer, intent(in) :: bnd
        complex(dp), intent(out) :: f(:)
        complex(dp), allocatable :: bl(:, :)
        real(dp), allocatable :: phi(:, :), ref(:, :)
        real(dp) :: x(2), wq, th, gx(6), gw(6), xi, eta, nh(2), ta, tb
        real(dp) :: psi, pr, pz, bt
        complex(dp) :: jj(3), kk
        integer :: nt, t, q, i, st, ir, nl, j, ig, te, jr

        nt = fs%mm%topology%n_triangles
        nl = fs%nln
        allocate (bl(nl, nt))
        bl = (0.0_dp, 0.0_dp)
        !$omp parallel default(shared) private(t, q, i, st, ir, jr, x, wq, th, jj, phi)
        allocate (phi(3, fs%sp(1)%nloc))
        !$omp do schedule(static)
        do t = fm%tlast(1) + 1, nt
            do jr = 2, fm%nring
                if (t <= fm%tlast(jr)) exit
            end do
            ir = jr - 1
            do q = 1, fs%nq
                call meridian_qpoint(fs, 1, t, q, x, wq, phi, st)
                th = atan2(x(2) - g%z_axis, x(1) - g%r_axis)
                jj = interp3(jcur(:, :, ir), th)
                do i = 1, nl
                    bl(i, t) = bl(i, t) + wq*(jj(1)*phi(1, i) + jj(2)*phi(2, i) &
                        + jj(3)*phi(3, i))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        f = (0.0_dp, 0.0_dp)
        do t = 1, nt
            do i = 1, nl
                f(fs%sp(1)%dof(i, t)) = f(fs%sp(1)%dof(i, t)) &
                    + mu0*real(fs%sp(1)%sgn(i, t), dp)*bl(i, t)
            end do
        end do
        if (bnd <= 0) return
        call gauss6(gx, gw)
        allocate (ref(2, nl), phi(3, fs%sp(1)%nloc))
        j = fm%nring
        do i = 1, fm%nth(j)
            call flux_mesh_ring_edge(fm, j, i, te, ta, tb)
            do ig = 1, 6
                call meridian_edge_point(1, gx(ig), xi, eta, nh)
                call meridian_point_basis(fs, pb, 1, te, xi, eta, x, phi, st, ref)
                th = ta + gx(ig)*(tb - ta)
                call psi_eval(g, x(1), x(2), psi, pr, pz)
                bt = -g%s_pol*g%sgn_dpsi*norm2([pr, pz])/x(1)
                kk = interp1(kap, th)*bt*x(1)
                do q = 1, nl
                    f(fs%sp(1)%dof(q, te)) = f(fs%sp(1)%dof(q, te)) &
                        + mu0*gw(ig)*kk*real(fs%sp(1)%sgn(q, te), dp)*ref(1, q)
                end do
            end do
        end do
    end subroutine ideal_load

    pure function interp3(arr, theta) result(v)
        !! Periodic linear interpolation of samples at uniform angles.
        complex(dp), intent(in) :: arr(:, :)
        real(dp), intent(in) :: theta
        complex(dp) :: v(size(arr, 1))
        real(dp) :: s, w
        integer :: i0, i1, np
        np = size(arr, 2)
        s = modulo(theta, 2.0_dp*pi)/(2.0_dp*pi)*np
        i0 = int(s)
        w = s - i0
        i1 = modulo(i0 + 1, np) + 1
        i0 = modulo(i0, np) + 1
        v = (1.0_dp - w)*arr(:, i0) + w*arr(:, i1)
    end function interp3

    pure function interp1(arr, theta) result(v)
        complex(dp), intent(in) :: arr(:)
        real(dp), intent(in) :: theta
        complex(dp) :: v
        real(dp) :: s, w
        integer :: i0, i1, np
        np = size(arr)
        s = modulo(theta, 2.0_dp*pi)/(2.0_dp*pi)*np
        i0 = int(s)
        w = s - i0
        i1 = modulo(i0 + 1, np) + 1
        i0 = modulo(i0, np) + 1
        v = (1.0_dp - w)*arr(i0) + w*arr(i1)
    end function interp1

    pure subroutine gauss6(xg, wg)
        !! 6-point Gauss-Legendre rule on [0, 1] (the ring rule).
        real(dp), intent(out) :: xg(6), wg(6)
        real(dp) :: z, p1, p2, p3, ppp, z1
        integer :: ii, jj
        do ii = 1, 3
            z = cos(pi*(ii - 0.25_dp)/6.5_dp)
            do
                p1 = 1
                p2 = 0
                do jj = 1, 6
                    p3 = p2
                    p2 = p1
                    p1 = ((2*jj - 1)*z*p2 - (jj - 1)*p3)/jj
                end do
                ppp = 6*(z*p1 - p2)/(z*z - 1)
                z1 = z
                z = z1 - p1/ppp
                if (abs(z - z1) < 1.0e-15_dp) exit
            end do
            xg(ii) = 0.5_dp*(1 - z)
            xg(7 - ii) = 0.5_dp*(1 + z)
            wg(ii) = 1/((1 - z*z)*ppp*ppp)
            wg(7 - ii) = wg(ii)
        end do
    end subroutine gauss6

    subroutine ideal_response(fs, pb, fm, g, amp, opt, n, rings, sheet_ring, &
            sheet_ra, sheet_m, msheet, b2sheet, phi_ext, a, lam, tol, maxit, &
            info, iters, resid)
        !! Linear ideal response: Nedelec coefficients a of the field of the
        !! volume (and boundary-sheet) plasma current, and the amplitudes lam
        !! of the resonant shielding sheets, for the external field stored in
        !! rings(:)%bc whose resonant fluxes on the sheet rings are phi_ext.
        !!
        !! The sheets are eliminated: for a volume field a the shielding
        !! amplitudes are lam(a) = -M^{-1} (phi_ext + Phi(curl a)) with the
        !! vacuum sheet matrix M (kin6d_mhd_response_n), so the total field
        !! B = B_ext + curl a + sum lam_s B_s has zero resonant flux.  The
        !! volume field solves the fixed point a = K^{-1} L(B(a)) (L the load
        !! of the plasma current), i.e. the linear system
        !!     a - K^{-1} L(B_hom(a)) = K^{-1} L(B_part),
        !! by restarted GMRES (fortnum).  rings(1 .. nr-1) are the half rings
        !! of strips 2 .. nring, rings(nr) the boundary ring.
        use fortnum_krylov, only: complex_gmres_operator
        use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_solve, ampere_n_curl
        use kin6d_mhd_response_n, only: ring_angles_t, ring_flux_harmonic
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        type(psi_grid_t), intent(in) :: g
        type(ampere_n_t), intent(inout) :: amp
        type(ideal_options_t), intent(in) :: opt
        integer, intent(in) :: n
        type(sample_ring_t), intent(in) :: rings(:)
        integer, intent(in) :: sheet_ring(:), sheet_m(:)
        type(ring_angles_t), intent(in) :: sheet_ra(:)
        complex(dp), intent(in) :: msheet(:, :), b2sheet(:, :), phi_ext(:)
        complex(dp), intent(inout) :: a(:)
        complex(dp), intent(out) :: lam(:)
        real(dp), intent(in) :: tol
        integer, intent(in) :: maxit
        integer, intent(out) :: info, iters
        real(dp), intent(out) :: resid
        complex(dp), allocatable :: bsheet(:, :, :, :), bext(:, :, :), rhs(:)
        complex(dp), allocatable :: minv(:, :), b(:, :, :)
        integer :: nr, ns, np, ir, s2, bnd, ipiv(max(1, size(sheet_m))), st
        logical :: debug
        character(len=8) :: env

        call get_environment_variable('KIN6D_L4_DEBUG', env)
        debug = len_trim(env) > 0
        nr = size(rings)
        ns = size(sheet_m)
        np = rings(1)%np
        bnd = 0
        if (opt%boundary_sheet) bnd = nr
        allocate (bsheet(3, np, nr, max(ns, 1)), bext(3, np, nr), b(3, np, nr))
        do s2 = 1, ns
            do ir = 1, nr
                call sample_field(fs, rings(ir), b2sheet(:, s2), &
                    bsheet(:, :, ir, s2))
            end do
        end do
        do ir = 1, nr
            bext(:, :, ir) = rings(ir)%bc
        end do
        allocate (minv(max(ns, 1), max(ns, 1)))
        if (ns > 0) then
            minv = (0.0_dp, 0.0_dp)
            do s2 = 1, ns
                minv(s2, s2) = (1.0_dp, 0.0_dp)
            end do
            block
                complex(dp) :: mcopy(ns, ns)
                mcopy = msheet
                call zgesv(ns, ns, mcopy, ns, ipiv, minv, ns, st)
            end block
            if (st /= 0) then
                info = -2
                return
            end if
        end if
        ! particular part: external field and its shielding
        if (ns > 0) lam = -matmul(minv, phi_ext)
        b = bext
        do s2 = 1, ns
            b = b + lam(s2)*bsheet(:, :, :, s2)
        end do
        allocate (rhs(amp%nn))
        call apply_k_inv_load(b, rhs)
        if (debug) print '(a,es12.4)', '  |rhs|', norm2(abs(rhs))
        call complex_gmres_operator(matvec, rhs, a, tol, maxit, 60, info, &
            iters, resid)
        ! total sheet amplitudes
        if (ns > 0) lam = lam + lam_hom(a)
    contains
        subroutine apply_k_inv_load(bf, y)
            complex(dp), intent(in) :: bf(:, :, :)
            complex(dp), intent(out) :: y(:)
            complex(dp), allocatable :: jc(:, :, :), kap(:), f(:)
            integer :: st2
            allocate (jc(3, np, nr), kap(np), f(amp%nn))
            call ideal_current(rings, n, g%sgn_dpsi, opt%qlow, bf, jc, kap, bnd)
            call ideal_load(fs, pb, fm, g, jc(:, :, 1:nr - 1), kap, bnd, f)
            call ampere_n_solve(amp, f, y, st2)
        end subroutine apply_k_inv_load
        function lam_hom(av) result(l)
            complex(dp), intent(in) :: av(:)
            complex(dp) :: l(ns), ph(ns)
            complex(dp), allocatable :: b2(:)
            integer :: r2
            call ampere_n_curl(amp, av, b2)
            do r2 = 1, ns
                ph(r2) = ring_flux_harmonic(fs, pb, fm, sheet_ring(r2), &
                    sheet_ra(r2), n, b2, sheet_m(r2))
            end do
            l = -matmul(minv, ph)
        end function lam_hom
        subroutine matvec(x, y)
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
            complex(dp), allocatable :: b2(:), bb(:, :, :), ky(:)
            complex(dp) :: lh(max(ns, 1))
            integer :: r2, s3
            call ampere_n_curl(amp, x, b2)
            allocate (bb(3, np, nr), ky(size(x)))
            do r2 = 1, nr
                call sample_field(fs, rings(r2), b2, bb(:, :, r2))
            end do
            if (ns > 0) then
                lh(1:ns) = lam_hom(x)
                do s3 = 1, ns
                    bb = bb + lh(s3)*bsheet(:, :, :, s3)
                end do
            end if
            call apply_k_inv_load(bb, ky)
            y = x - ky
            if (debug) print '(a,3es12.4)', '  matvec |x| |Kx| |y|', norm2(abs(x)), &
                norm2(abs(ky)), norm2(abs(y))
        end subroutine matvec
    end subroutine ideal_response

end module kin6d_mhd_ideal_n
