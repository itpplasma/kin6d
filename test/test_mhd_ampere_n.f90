program test_mhd_ampere_n
    !! L4 Ampere operator with the exact exterior (V3-type manufactured
    !! oracle with an independent field).
    !!
    !! Current J = curl_n W of a smooth vector bump W supported in a disk
    !! inside the circular domain |x - (R0, 0)| < a.  Then B = mu0 (W +
    !! grad_n chi) with chi(x) = int G3(x, y) div W(y) dV_y, G3 = 1/(4 pi |x -
    !! y|).  Outside the bump B = mu0 grad_n chi and the exterior potential is
    !! u = mu0 chi.  chi_n is evaluated by tensor quadrature over the bump
    !! and the trapezoidal rule in the toroidal angle (no toroidal functions,
    !! no FEM/BEM code), so the oracle is independent of the kernels used by
    !! the boundary operators.  Checked: B at points between the bump and
    !! Gamma, and the exterior potential trace on Gamma, converge at order k
    !! for elements of order k (geometry order k); the Hermitian operator is
    !! positive (energy check).
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_load, meridian_point_basis_t, meridian_point_basis_init
    use fortsparse, only: csc_z_t, csc_matvec, csc_conjugate_transpose, &
        fortsparse_status_t
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, flux_mesh_locate
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl, ampere_n_trace, ampere_n_field_point, mu0
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: r0 = 3.0_dp, amin = 1.0_dp
    real(dp), parameter :: xc(2) = [3.1_dp, 0.05_dp], rho = 0.4_dp
    complex(dp), parameter :: wvec(3) = [(1.0_dp, 0.5_dp), (-0.3_dp, 0.8_dp), &
        (0.6_dp, -0.2_dp)]
    integer, parameter :: nsample = 48
    real(dp), parameter :: usample(2) = [0.7_dp, 0.8_dp]
    logical :: all_passed
    integer :: nmode

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call convergence(1, 3)
    call convergence(2, 3)
    call convergence(2, 1)
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_ampere_n: all checks passed"

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        if (c) then
            print '(a,a)', "  PASS ", d
        else
            print '(a,a)', "  FAIL ", d
        end if
    end subroutine record

    subroutine circle_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + amin*u*cos(theta), amin*u*sin(theta)]
    end subroutine circle_map

    pure subroutine bump(x, b, br, bz)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b, br, bz
        real(dp) :: s, d(2)
        d = x - xc
        s = 1.0_dp - (d(1)**2 + d(2)**2)/rho**2
        if (s <= 0.0_dp) then
            b = 0.0_dp
            br = 0.0_dp
            bz = 0.0_dp
            return
        end if
        b = s**4
        br = -8.0_dp*s**3*d(1)/rho**2
        bz = -8.0_dp*s**3*d(2)/rho**2
    end subroutine bump

    subroutine wfun(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: b, br, bz
        call bump(x, b, br, bz)
        v = b*wvec/mu0
    end subroutine wfun

    subroutine chi_exact(x, chi, grad)
        !! chi_n and grad_n chi_n at x outside the bump.
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: chi, grad(3)
        integer, parameter :: nr = 24, na = 48, nd = 256
        real(dp) :: gr(nr), gw(nr), y(2), b, br, bz, rr, ang, wq, del, d, cr
        real(dp) :: g0, gdr, gdz, c
        complex(dp) :: divw
        integer :: i, j, k

        call gauss01(nr, gr, gw)
        chi = 0
        grad = 0
        do i = 1, nr
            rr = rho*gr(i)
            do j = 1, na
                ang = 2.0_dp*pi*(j - 1)/na
                y = xc + rr*[cos(ang), sin(ang)]
                call bump(y, b, br, bz)
                ! div_n W = W_R/R + d_R W_R + i n W_phi/R + d_Z W_Z
                divw = (b*wvec(1)/y(1) + br*wvec(1) &
                    + cmplx(0.0_dp, nmode, dp)*b*wvec(2)/y(1) + bz*wvec(3))/mu0
                wq = rho*gw(i)*rr*2.0_dp*pi/na*y(1)
                g0 = 0
                gdr = 0
                gdz = 0
                do k = 1, nd
                    del = 2.0_dp*pi*(k - 1)/nd
                    c = cos(nmode*del)
                    cr = cos(del)
                    d = sqrt(x(1)**2 + y(1)**2 - 2*x(1)*y(1)*cr &
                        + (x(2) - y(2))**2)
                    g0 = g0 + c/d
                    gdr = gdr - c*(x(1) - y(1)*cr)/d**3
                    gdz = gdz - c*(x(2) - y(2))/d**3
                end do
                g0 = g0*2.0_dp*pi/nd/(4.0_dp*pi)
                gdr = gdr*2.0_dp*pi/nd/(4.0_dp*pi)
                gdz = gdz*2.0_dp*pi/nd/(4.0_dp*pi)
                chi = chi + wq*divw*g0
                grad(1) = grad(1) + wq*divw*gdr
                grad(3) = grad(3) + wq*divw*gdz
            end do
        end do
        grad(2) = cmplx(0.0_dp, nmode, dp)*chi/x(1)
    end subroutine chi_exact

    subroutine gauss01(n, x, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)
        real(dp) :: z, p1, p2, p3, pp, z1
        integer :: i, j
        do i = 1, (n + 1)/2
            z = cos(pi*(i - 0.25_dp)/(n + 0.5_dp))
            do
                p1 = 1
                p2 = 0
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j
                end do
                pp = n*(z*p1 - p2)/(z*z - 1)
                z1 = z
                z = z1 - p1/pp
                if (abs(z - z1) < 1.0e-15_dp) exit
            end do
            x(i) = 0.5_dp*(1 - z)
            x(n + 1 - i) = 0.5_dp*(1 + z)
            w(i) = 1/((1 - z*z)*pp*pp)
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss01

    subroutine run_level(k, nrad, eb, eu, energy)
        integer, intent(in) :: k, nrad
        real(dp), intent(out) :: eb, eu, energy
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(ampere_n_t) :: amp
        type(meridian_point_basis_t) :: pb
        type(csc_z_t) :: ch
        type(fortsparse_status_t) :: fst
        complex(dp), allocatable :: g2(:), f(:), a(:), b2(:), u(:), uex(:)
        complex(dp) :: bh(3), bex(3), chi, grad(3)
        real(dp), allocatable :: ul(:)
        real(dp) :: x(2), xs(2), th, bmax, umax, xi, eta
        integer :: st, j, t, ic

        allocate (ul(nrad))
        do j = 1, nrad
            ul(j) = real(j, dp)/nrad
        end do
        call flux_mesh_build(ul, 8, 4096, k, circle_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, k, fs, st)
        if (st /= 0) error stop "feec"
        call ampere_n_build(fs, nmode, amp, st)
        if (st /= 0) error stop "ampere build"
        call meridian_load(fs, 2, wfun, g2, st)
        call csc_conjugate_transpose(amp%curl, ch, fst)
        f = mu0*csc_matvec(ch, g2)
        allocate (a(amp%nn))
        call ampere_n_solve(amp, f, a, st)
        if (st /= 0) error stop "solve"
        energy = real(dot_product(a, f), dp)
        call ampere_n_curl(amp, a, b2)
        call ampere_n_trace(amp, a, u)
        ! exterior potential trace
        allocate (uex(size(u)))
        !$omp parallel do private(j, chi, grad) schedule(dynamic)
        do j = 1, size(u)
            call chi_exact(amp%bd%curve%nodes(:, j), chi, grad)
            uex(j) = mu0*chi
        end do
        !$omp end parallel do
        umax = maxval(abs(uex))
        eu = maxval(abs(u - uex))/umax
        ! field at sample points
        call meridian_point_basis_init(fs, pb, st)
        eb = 0
        bmax = 0
        do ic = 1, 2
            do j = 1, nsample
                th = 2.0_dp*pi*(j - 0.5_dp)/nsample
                call circle_map(usample(ic), th, xs)
                call flux_mesh_locate(fm, xs, usample(ic), th, t, xi, eta, st)
                if (st /= 0) error stop "locate"
                call ampere_n_field_point(fs, pb, b2, t, xi, eta, x, bh, st)
                call chi_exact(xs, chi, grad)
                bex = mu0*grad
                eb = eb + sum(abs(bh - bex)**2)
                bmax = bmax + sum(abs(bex)**2)
            end do
        end do
        eb = sqrt(eb)
        bmax = sqrt(bmax)
        eb = eb/bmax
    end subroutine run_level

    subroutine convergence(k, n)
        integer, intent(in) :: k, n
        integer, parameter :: nlev = 3
        real(dp) :: eb(nlev), eu(nlev), en(nlev), rb, ru
        integer :: lev
        character(len=60) :: tag

        nmode = n
        do lev = 1, nlev
            call run_level(k, 4*2**lev, eb(lev), eu(lev), en(lev))
            print '(a,i0,a,i0,a,i3,a,2es11.3,a,es12.4)', "  k=", k, " n=", n, &
                " nrad=", 4*2**lev, " err B/u", eb(lev), eu(lev), &
                " energy", en(lev)
        end do
        ! mean rate over the ladder (pointwise samples make single-step
        ! rates oscillate with the position of the samples in their cells)
        rb = log(eb(1)/eb(nlev))/(log(2.0_dp)*(nlev - 1))
        ru = log(eu(1)/eu(nlev))/(log(2.0_dp)*(nlev - 1))
        write (tag, '(a,i0,a,i0)') "order ", k, ", n ", n
        print '(a,a,a,2f7.2)', "  ", trim(tag), " observed orders B/u", rb, ru
        call record(rb > k - 0.3_dp, "field order, "//trim(tag))
        call record(ru > k - 0.3_dp, "exterior trace order, "//trim(tag))
        call record(all(en > 0.0_dp), "positive energy, "//trim(tag))
        call record(eb(nlev) < 0.1_dp .and. eu(nlev) < 0.01_dp, &
            "fine-level accuracy, "//trim(tag))
    end subroutine convergence

end program test_mhd_ampere_n
