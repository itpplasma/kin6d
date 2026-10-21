program test_mhd_shield_torus
    !! L4 ideal shielding in a circular torus: mesh convergence of the
    !! resonant flux of the shielding sheet.
    !!
    !! Circular CHEASE equilibrium of aspect ratio 10 (q0 = 1.55, q = 2 at
    !! psin 0.37, q = 3 at psin 0.95 for n = 1) and the exact vacuum drive
    !! grad_n (R^n (1 + Z/2)).  The surfaces q = 2 and q = 3 are both
    !! resonant (the toroidal coupling gives the sideband m = 3 of the
    !! m = 2 drive); the m = 3 surface lies near the edge.  The shielding
    !! response (model 'ideal', closure 'sheet': MEPHIT-type volume current
    !! plus a sheet with zero resonant flux on every rational ring) is
    !! computed for several numbers of rings.  The strip current of the half rings has a pole
    !! 1/(u - u_s) at every resonant surface; a surface without sheet
    !! (resonant_surfaces_missing lists them all) makes the response of the
    !! other surfaces depend on where its pole falls between the half rings
    !! (observed in KIN6D response runs: |b_shield| of q = 2 at nrad 64 / 96 /
    !! 128 / 160 = 1.66 / 1.95 / 1.72 / 1.70e-4 T, converged and equal to
    !! GPEC's 1.633e-4 T with the q = 3 sheet).
    !! Checks (q = 2 sheet): |b_shield| differs by less than 2 % between all
    !! levels (0.5 % observed; 16 % with the q = 3 sheet omitted) and the
    !! successive changes do not grow.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_vacuum_surfaces, only: psi_grid_t, eqdsk_read, flux_ray_point, &
        surface_t, surface_build
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build, ring_labels_merge
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        ring_points, ring_flux_harmonic, ring_flux_from_values, ring_sheet_load, &
        ring_area, layer_closure_t, closure_shield_init, resonant_surfaces_missing
    use kin6d_mhd_ideal_n, only: sample_ring_t, sample_ring_build, &
        ideal_options_t, ideal_response
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: nlev = 4, nlist(nlev) = [48, 64, 80, 96]
    integer, parameter :: n = 1, np = 128
    type(psi_grid_t) :: g
    real(dp), allocatable :: psin_s(:)
    integer, allocatable :: m_s(:)
    real(dp) :: bsh(nlev), spread, dlast, dfirst
    logical :: all_passed
    integer :: l

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call make_equilibrium()
    do l = 1, nlev
        call shield_level(nlist(l), bsh(l))
        print '(a,i4,a,es14.6)', "  nrad", nlist(l), ": |b_shield| (T)", bsh(l)
    end do
    spread = (maxval(bsh) - minval(bsh))/bsh(nlev)
    dfirst = abs(bsh(2) - bsh(1))
    dlast = abs(bsh(nlev) - bsh(nlev - 1))
    print '(a,es10.2)', "  relative spread over the levels", spread
    call record(spread < 0.02_dp, "|b_shield| within 2 % over all levels")
    call record(dlast < 0.01_dp*bsh(nlev), "finest two levels agree to 1 %")
    call record(dlast <= max(dfirst, 1.0e-3_dp*bsh(nlev)), &
        "successive changes do not grow")
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_shield_torus: all checks passed"

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

    subroutine make_equilibrium()
        !! Coarse (65 x 65) copy of the circular CHEASE equilibrium of the
        !! KIN6D response study: A = 10, q0 = 1.55, q = 2 at psin 0.37 and
        !! q = 3 at psin 0.95 (n = 1).
        character(len=512) :: path
        real(dp), allocatable :: qs(:)
        integer :: st
        if (command_argument_count() /= 1) error stop "usage: test fixture.eqdsk"
        call get_command_argument(1, path)
        call eqdsk_read(trim(path), g, st)
        if (st /= 0) error stop "equilibrium"
        call resonant_surfaces_missing(g, n, 0.05_dp, 1.0_dp, 256, [integer ::], &
            psin_s, qs, m_s)
        if (size(psin_s) /= 2) error stop "expected the q = 2 and q = 3 surfaces"
        print '(a,2f12.8)', "  q = 2, 3 at psin", psin_s
    end subroutine make_equilibrium

    subroutine polar_map(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        real(dp) :: rho
        logical :: ok
        if (u <= 0.0_dp) then
            x = [g%r_axis, g%z_axis]
            return
        end if
        call flux_ray_point(g, u*u, theta, rho, ok)
        if (.not. ok) error stop "ray"
        x = [g%r_axis + rho*cos(theta), g%z_axis + rho*sin(theta)]
    end subroutine polar_map

    subroutine bext_fun(x, b)
        !! Vacuum per-n field grad_n (R^n (1 + Z/2)) (curl and divergence free).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: b(3)
        real(dp) :: u
        u = x(1)**n*(1.0_dp + 0.5_dp*x(2))
        b(1) = n*x(1)**(n - 1)*(1.0_dp + 0.5_dp*x(2))
        b(2) = cmplx(0.0_dp, real(n, dp), dp)*u/x(1)
        b(3) = 0.5_dp*x(1)**n
        b = b*1.0e-3_dp
    end subroutine bext_fun

    subroutine shield_level(nrad, bshield)
        !! |b_shield| = |M_ss lambda|/area of the first (lowest m) sheet.
        integer, intent(in) :: nrad
        real(dp), intent(out) :: bshield
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ideal_options_t) :: opt
        type(sample_ring_t), allocatable :: rings(:)
        type(layer_closure_t) :: cl
        type(ring_angles_t), allocatable :: ra(:)
        type(surface_t) :: surf
        real(dp), allocatable :: ul(:), xg(:, :), nds(:, :), thg(:), tau(:), wq(:)
        integer, allocatable :: tg(:), jr(:)
        complex(dp), allocatable :: f(:), a(:), bb(:), avol(:), bn(:), b2s(:, :)
        complex(dp), allocatable :: phivac(:), lam(:), phitot(:)
        complex(dp) :: mm, bv(3)
        real(dp) :: area, res
        integer :: j, k, nr, ns, s, st, info, its

        ns = size(psin_s)
        call ring_labels_merge(nrad, sqrt(psin_s), 0.3_dp, ul)
        nr = size(ul)
        allocate (jr(ns), ra(ns), b2s(1, ns), phivac(ns), lam(ns), phitot(ns))
        do s = 1, ns
            jr(s) = minloc(abs(ul - sqrt(psin_s(s))), 1)
        end do
        call flux_mesh_build(ul, 8, 128, 2, polar_map, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        call ampere_n_build(fs, n, amp, st)
        call meridian_point_basis_init(fs, pb, st)
        allocate (rings(nr))
        do j = 2, nr
            call sample_ring_build(g, opt, fm, fs, pb, 0.5_dp*(ul(j - 1) + ul(j)), &
                j, np, rings(j - 1), st)
            if (st /= 0) error stop "half ring"
        end do
        call sample_ring_build(g, opt, fm, fs, pb, ul(nr), nr, np, rings(nr), st)
        if (st /= 0) error stop "boundary ring"
        do j = 1, nr
            do k = 1, np
                call bext_fun(rings(j)%x(:, k), rings(j)%bc(:, k))
            end do
        end do
        deallocate (b2s)
        allocate (b2s(fs%sp(2)%ndof, ns), f(amp%nn), a(amp%nn))
        do s = 1, ns
            ! unit sheet on the rational ring, its self flux and the drive flux
            call surface_build(g, ul(jr(s))**2, 1024, surf, st)
            if (st /= 0) error stop "rational surface"
            call ring_angles_set(surf%q, surf%theta_s(:, 2), surf%nu, ra(s))
            f = (0.0_dp, 0.0_dp)
            call ring_sheet_load(fs, pb, fm, jr(s), ra(s), n, m_s(s), &
                (1.0_dp, 0.0_dp), f)
            call ampere_n_solve(amp, f, a, st)
            call ampere_n_curl(amp, a, bb)
            b2s(:, s) = bb
            if (s == 1) then
                mm = ring_flux_harmonic(fs, pb, fm, jr(s), ra(s), n, b2s(:, s), &
                    m_s(s))
                area = ring_area(fm, jr(s))
            end if
            call ring_points(fm, jr(s), xg, nds, thg, tg, tau, wq)
            if (allocated(bn)) deallocate (bn)
            allocate (bn(size(thg)))
            do k = 1, size(thg)
                call bext_fun(xg(:, k), bv)
                bn(k) = (bv(1)*nds(1, k) + bv(3)*nds(2, k))*xg(1, k)
            end do
            phivac(s) = ring_flux_from_values(ra(s), n, m_s(s), thg, bn, wq)
        end do
        call closure_shield_init(jr, ra, m_s, cl)
        allocate (avol(amp%nn))
        avol = (0.0_dp, 0.0_dp)
        call ideal_response(fs, pb, fm, g, amp, opt, n, rings, cl, b2s, phivac, &
            avol, lam, phitot, 1.0e-9_dp, 400, info, its, res)
        if (info /= 0) then
            print *, "GMRES info", info, its, res
            error stop "GMRES"
        end if
        bshield = abs(mm*lam(1))/area
        print '(a,i4,a,i4,a,i4)', "  nrad", nrad, ": rings", nr, ", GMRES", its
    end subroutine shield_level

end program test_mhd_shield_torus
