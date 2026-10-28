module kin6d_island_seed
    ! Energy-criterion island seeding on the NN1 FEEC Fourier complex
    ! (docs/island-seed.md; item 4 of docs/relaxation-comparators.md).
    !
    ! For every resonance iota(r_mn) = n/m of the axisymmetric (n = 0)
    ! background inside the represented annulus, two seed potentials
    !     dA = alpha((r - r_mn)/w) B0/|B0| {cos, sin}(m theta - n phi),
    ! |r - r_mn| < w = min(1/|m|, distance to the boundary) in units of the
    ! radial scale, are projected into V1 with zero boundary trace and moved
    ! to the R A_phi = 0 gauge. The magnetic energy W = int |B|^2/(2 mu0) is
    ! quadratic in their amplitudes a,
    !     W(B + sum a_i curl dA_i) = W(B) + a.g + a.G a/2,
    !     g_i = <J, dA_i> = (curl dA_i, B)/mu0,  G_ij = (curl dA_i, curl dA_j)/mu0,
    ! so the model is exact; the seed is a* = -G^{-1} g. The pendulum width
    ! of each chain is reported before and after seeding. The state, and the
    ! optional gauge potentials of the caller, are updated in place so the
    ! existing relaxations (nn_solve, relax_run) can continue from it.
    !
    ! Scope: the background n = 0 field must have nested surfaces close to
    ! circles about field%centre (circularity is reported); theta is the
    ! geometric chart angle. Only coprime (m, n) with 1 <= n <= modes and
    ! 1 <= |m| <= m_max are seeded. A symmetric background has g = 0 for
    ! n > 0 seeds by Fourier orthogonality, so seeding is non-trivial only
    ! for states with resonant n > 0 content (e.g. shielded nested states).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_meridian_feec_n, only: meridian_point_basis, meridian_op_complex, &
        meridian_qpoint, meridian_real_to_complex
    use fortsparse, only: csc_z_t, csc_matvec
    use kin6d_mhd_mesh, only: flux_mesh_locate
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_factor_t, nn_factor, &
        nn_factor_solve, nn_factor_free
    use kin6d_island_seed_leaves, only: circle_rotation, seed_potential, &
        pendulum_width
    implicit none
    private
    public :: seed_options_t, seed_resonance_t, seed_report_t
    public :: island_seed, island_resonances, seed_rotation, seed_resonant_br
    public :: seed_circle_quadrature

    real(dp), parameter :: pi = acos(-1.0_dp)

    type :: seed_options_t
        integer :: m_max = 6
        ! Highest toroidal mode seeded; -1 means every mode of the field.
        integer :: n_max = -1
        ! Angular scan samples per circle (at least 8 per ring vertex) for
        ! the element crossings of the circle quadrature, and radial scan
        ! intervals for the resonance roots.
        integer :: ntheta = 256, nscan = 64
        ! Half window of the least-squares d iota/du (units of the radial scale).
        real(dp) :: shear_step = 0.1_dp
        ! Radial margin (radial-scale units) kept from the boundaries/axis.
        real(dp) :: margin = 1.0e-3_dp
    end type

    type :: seed_resonance_t
        integer :: m = 0, n = 0
        ! Normalised and physical radius, iota, d iota/dr, support half width
        ! (normalised), averaged toroidal field and metric major radius.
        real(dp) :: u = 0, r = 0, iota = 0, shear = 0, support = 0
        real(dp) :: bphi = 0, rmajor = 0
        ! Background radial field relative to the poloidal field on r_mn.
        real(dp) :: circularity = 0
        ! Seed amplitudes of the cos and sin phase.
        real(dp) :: amplitude(2) = 0
        ! Resonant radial-field amplitude, pendulum width and O-point angle
        ! at phi = 0 (one of the m chain angles), before/after seeding.
        real(dp) :: dbr_before = 0, dbr_after = 0
        real(dp) :: width_before = 0, width_after = 0, theta_o = 0
    end type

    type :: seed_report_t
        integer :: status = 1, count = 0
        type(seed_resonance_t), allocatable :: res(:)
        real(dp), allocatable :: gradient(:), gram(:, :), amplitude(:)
        ! Energy change predicted by the model and recomputed from the
        ! updated V2 coefficients, both in the units of W.
        real(dp) :: energy_model = 0, energy_direct = 0
        real(dp) :: symmetry_defect = 0, min_eigenvalue = 0
        ! Largest |s_a| of the gauge-transformed seeds relative to the seed.
        real(dp) :: gauge_defect = 0
        character(64) :: message = 'not run'
    end type
contains
    subroutine island_seed(field, options, report, potentials)
        ! potentials(:, n), n = 1..modes: the caller's Nedelec potentials in
        ! the R A_phi = 0 gauge (nn_solve solutions); seeds are added.
        type(free_b_field_t), intent(inout) :: field
        type(seed_options_t), intent(in) :: options
        type(seed_report_t), intent(out) :: report
        complex(dp), intent(inout), optional :: potentials(:, :)
        type(seed_resonance_t), allocatable :: found(:)
        type(nn_factor_t) :: mass
        type(csc_z_t) :: curl, grad, m2
        complex(dp), allocatable :: sa(:, :), sb(:, :), load(:), x(:), chi(:)
        complex(dp), allocatable :: mb(:, :)
        real(dp), allocatable :: g(:), gm(:, :), work(:, :), eig(:)
        real(dp) :: u_lo, u_hi, scale, w_before, w_after, theta_unseeded
        integer, allocatable :: mode_of(:)
        integer :: nmax, k, i, j, nn, n1, st, info, n, ns
        nn = 0
        if (.not. allocated(field%b)) then
            report%message = 'field has no modes'; return
        end if
        nmax = ubound(field%b, 2)
        if (options%n_max >= 0) nmax = min(nmax, options%n_max)
        if (options%m_max < 1 .or. options%ntheta < 8 .or. options%nscan < 2 &
            .or. options%shear_step <= 0 .or. options%margin < 0) then
            report%message = 'invalid options'; return
        end if
        if (present(potentials)) then
            if (size(potentials, 1) /= field%fs%nn .or. &
                size(potentials, 2) < nmax) then
                report%message = 'potential shape'; return
            end if
        end if
        call radial_range(field, options%margin, u_lo, u_hi)
        call find_resonances(field, options, nmax, u_lo, u_hi, found, st)
        if (st /= 0) then
            report%message = 'resonance search failed'; return
        end if
        report%count = size(found)
        if (report%count == 0) then
            allocate (report%res(0), report%gradient(0), report%gram(0, 0), &
                report%amplitude(0))
            report%status = 0; report%message = 'no resonance in range'
            return
        end if
        ns = 2*report%count
        nn = field%fs%nn; n1 = field%fs%sp(1)%ndof
        allocate (sa(nn, ns), sb(field%fs%sp(2)%ndof, ns), mode_of(ns))
        call nn_factor(field%fs%mass(1), field%fs%sp(1)%bnd, mass, st)
        if (st /= 0) then
            report%message = 'V1 mass factor'; return
        end if
        do k = 1, report%count
            n = found(k)%n
            call seed_load(field, found(k), load, st)
            if (st /= 0) exit
            allocate (x(n1))
            call nn_factor_solve(mass, load, x, st)
            if (st /= 0) exit
            ! Gauge to R A_phi = 0: x - grad_n(a/(i n)); same curl exactly.
            call meridian_op_complex(field%fs%grad, n, grad, st)
            if (st /= 0) exit
            chi = x(nn + 1:)/cmplx(0.0_dp, real(n, dp), dp)
            x = x - csc_matvec(grad, chi)
            report%gauge_defect = max(report%gauge_defect, &
                maxval(abs(x(nn + 1:)))/max(maxval(abs(x(:nn))), tiny(1.0_dp)))
            x(nn + 1:) = 0
            call meridian_op_complex(field%fs%curl, n, curl, st)
            if (st /= 0) exit
            sa(:, 2*k - 1) = x(:nn)
            sb(:, 2*k - 1) = csc_matvec(curl, x)
            ! sin(m theta - n phi) has complex amplitude i e^{-i m theta}/2.
            sa(:, 2*k) = cmplx(0.0_dp, 1.0_dp, dp)*sa(:, 2*k - 1)
            sb(:, 2*k) = cmplx(0.0_dp, 1.0_dp, dp)*sb(:, 2*k - 1)
            mode_of(2*k - 1:2*k) = n
            deallocate (x)
        end do
        call nn_factor_free(mass)
        if (st /= 0) then
            report%message = 'seed projection'; return
        end if
        ! W = sum_n c_n b_n^H M2 b_n with c_0 = pi/mu0, c_n = 2 pi/mu0 (the
        ! real field 2 Re(b_n e^{i n phi}) and the dropped 2 pi of FortFEM).
        call meridian_real_to_complex(field%fs%mass(2), m2)
        allocate (mb(size(sb, 1), ns), g(ns), gm(ns, ns))
        do i = 1, ns
            mb(:, i) = csc_matvec(m2, sb(:, i))
        end do
        scale = 4*pi/field%mu0
        gm = 0
        do i = 1, ns
            g(i) = scale*real(dot_product(mb(:, i), field%b(:, mode_of(i))), dp)
            do j = 1, ns
                if (mode_of(j) /= mode_of(i)) cycle
                gm(i, j) = scale*real(dot_product(sb(:, i), mb(:, j)), dp)
            end do
        end do
        report%symmetry_defect = maxval(abs(gm - transpose(gm))) &
            /max(maxval(abs(gm)), tiny(1.0_dp))
        gm = (gm + transpose(gm))/2
        work = gm
        allocate (eig(ns))
        call symmetric_eigenvalues(work, eig, info)
        if (info /= 0) then
            report%message = 'Gram eigenvalues'; return
        end if
        report%min_eigenvalue = eig(1)/max(eig(ns), tiny(1.0_dp))
        work = gm
        report%amplitude = -g
        call cholesky_solve(work, report%amplitude, info)
        if (info /= 0) then
            report%message = 'Gram matrix not positive definite'; return
        end if
        report%gradient = g; report%gram = gm
        report%energy_model = dot_product(g, report%amplitude) &
            + dot_product(report%amplitude, matmul(gm, report%amplitude))/2
        ! Apply the seed; the direct energy change uses the updated modes.
        report%energy_direct = 0
        do n = 1, nmax
            if (.not. any(mode_of == n)) cycle
            report%energy_direct = report%energy_direct - mode_energy(n)
            do i = 1, ns
                if (mode_of(i) /= n) cycle
                field%b(:, n) = field%b(:, n) + report%amplitude(i)*sb(:, i)
                if (present(potentials)) potentials(:, n) = potentials(:, n) &
                    + report%amplitude(i)*sa(:, i)
            end do
            report%energy_direct = report%energy_direct + mode_energy(n)
        end do
        do k = 1, report%count
            found(k)%amplitude = report%amplitude(2*k - 1:2*k)
        end do
        ! Widths: dbr_before from the unseeded modes, recovered by linearity.
        do k = 1, report%count
            call resonant_width(field, found(k), options%ntheta, &
                found(k)%dbr_after, w_after, found(k)%theta_o, st)
            if (st /= 0) exit
            found(k)%width_after = w_after
            do i = 1, ns
                if (mode_of(i) == found(k)%n) field%b(:, found(k)%n) = &
                    field%b(:, found(k)%n) - report%amplitude(i)*sb(:, i)
            end do
            call resonant_width(field, found(k), options%ntheta, &
                found(k)%dbr_before, w_before, theta_unseeded, st)
            do i = 1, ns
                if (mode_of(i) == found(k)%n) field%b(:, found(k)%n) = &
                    field%b(:, found(k)%n) + report%amplitude(i)*sb(:, i)
            end do
            if (st /= 0) exit
            found(k)%width_before = w_before
        end do
        if (st /= 0) then
            report%message = 'width diagnostic'; return
        end if
        report%res = found
        report%status = 0; report%message = 'seeded'
    contains
        real(dp) function mode_energy(mode) result(e)
            integer, intent(in) :: mode
            e = scale/2*real(dot_product(field%b(:, mode), &
                csc_matvec(m2, field%b(:, mode))), dp)
        end function
    end subroutine

    subroutine seed_load(field, res, load, status)
        ! load_i = int dA . phi_i R dR dZ for the V1 basis, with the complex
        ! mode-n amplitude dA = alpha(xi) B0/|B0| e^{-i m theta}/2 of
        ! cos(m theta - n phi). B0 is evaluated in the same cell from its V2
        ! coefficients; no point location and no procedure argument.
        type(free_b_field_t), intent(in) :: field
        type(seed_resonance_t), intent(in) :: res
        complex(dp), allocatable, intent(out) :: load(:)
        integer, intent(out) :: status
        complex(dp), allocatable :: cell(:, :)
        real(dp), allocatable :: phi1(:, :), phi2(:, :)
        real(dp) :: x(2), wq, w2, u, theta, xi, b(3), a(3), radius
        complex(dp) :: v(3)
        integer :: nt, nl1, t, q, i, bad, st
        nt = field%fs%mm%topology%n_triangles
        nl1 = field%fs%sp(1)%nloc
        allocate (cell(nl1, nt), load(field%fs%sp(1)%ndof))
        bad = 0
        !$omp parallel default(shared) private(t, q, i, phi1, phi2, x, wq, w2, u, &
        !$omp& theta, xi, b, a, v, radius, st) reduction(max:bad)
        allocate (phi1(3, nl1), phi2(3, field%fs%sp(2)%nloc))
        !$omp do schedule(static)
        do t = 1, nt
            cell(:, t) = 0
            do q = 1, field%fs%nq
                call meridian_qpoint(field%fs, 1, t, q, x, wq, phi1, st)
                if (st /= 0) then
                    bad = 1; cycle
                end if
                call field%coordinates(field%centre, field%radial_scale, x, u, theta)
                xi = (u - res%u)/res%support
                if (abs(xi) >= 1) cycle
                call meridian_qpoint(field%fs, 2, t, q, x, w2, phi2, st)
                if (st /= 0) then
                    bad = 1; cycle
                end if
                b = matmul(phi2, real(field%b(field%fs%sp(2)%dof(:, t), 0), dp) &
                    *field%fs%sp(2)%sgn(:, t))
                radius = x(1)
                if (field%fs%cyl_r0 > 0) radius = field%fs%cyl_r0
                b(2) = b(2) + field%f_vac/radius
                call seed_potential(b(1), b(2), b(3), xi, a(1), a(2), a(3))
                v = a*exp(cmplx(0.0_dp, -res%m*theta, dp))/2
                do i = 1, nl1
                    cell(i, t) = cell(i, t) + wq*sum(v*phi1(:, i))
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        status = 1
        if (bad /= 0) return
        load = 0
        do t = 1, nt
            do i = 1, nl1
                load(field%fs%sp(1)%dof(i, t)) = load(field%fs%sp(1)%dof(i, t)) &
                    + field%fs%sp(1)%sgn(i, t)*cell(i, t)
            end do
        end do
        status = 0
    end subroutine

    subroutine island_resonances(field, options, found, status)
        ! The resonances island_seed would seed, without seeding.
        type(free_b_field_t), intent(in) :: field
        type(seed_options_t), intent(in) :: options
        type(seed_resonance_t), allocatable, intent(out) :: found(:)
        integer, intent(out) :: status
        real(dp) :: u_lo, u_hi
        integer :: nmax
        status = 1
        if (.not. allocated(field%b)) return
        nmax = ubound(field%b, 2)
        if (options%n_max >= 0) nmax = min(nmax, options%n_max)
        call radial_range(field, options%margin, u_lo, u_hi)
        call find_resonances(field, options, nmax, u_lo, u_hi, found, status)
    end subroutine

    subroutine radial_range(field, margin, u_lo, u_hi)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: margin
        real(dp), intent(out) :: u_lo, u_hi
        ! On a disc the axis cells do not resolve iota: start at half the
        ! first ring. An annulus starts at its inner boundary.
        u_lo = field%mesh%u(1)/2
        if (field%mesh%tlast(1) == 0) u_lo = field%mesh%u(1)
        u_hi = field%mesh%u(field%mesh%nring)
        u_lo = u_lo + margin*(u_hi - u_lo)
        u_hi = u_hi - margin*(u_hi - u_lo)
    end subroutine

    subroutine find_resonances(field, options, nmax, u_lo, u_hi, found, status)
        type(free_b_field_t), intent(in) :: field
        type(seed_options_t), intent(in) :: options
        integer, intent(in) :: nmax
        real(dp), intent(in) :: u_lo, u_hi
        type(seed_resonance_t), allocatable, intent(out) :: found(:)
        integer, intent(out) :: status
        type(seed_resonance_t) :: res
        real(dp), allocatable :: grid(:), iota(:)
        real(dp) :: a, b, fa, fm, mid, target, h, ip, im, inner, outer, circ
        integer, parameter :: nfit = 8
        integer :: i, m, n, it
        allocate (found(0))
        status = 1
        allocate (grid(0:options%nscan), iota(0:options%nscan))
        do i = 0, options%nscan
            grid(i) = u_lo + (u_hi - u_lo)*i/real(options%nscan, dp)
            call seed_rotation(field, grid(i), options%ntheta, iota(i), circ, &
                res%bphi, status)
            if (status /= 0) return
        end do
        inner = 0
        if (field%mesh%tlast(1) == 0) inner = field%mesh%u(1)
        outer = field%mesh%u(field%mesh%nring)
        do n = 1, nmax
            do m = -options%m_max, options%m_max
                if (m == 0) cycle
                if (gcd(abs(m), n) /= 1) cycle
                target = real(n, dp)/m
                do i = 1, options%nscan
                    if ((iota(i - 1) > target) .eqv. (iota(i) > target)) cycle
                    a = grid(i - 1); b = grid(i); fa = iota(i - 1) - target
                    do it = 1, 60
                        mid = (a + b)/2
                        call seed_rotation(field, mid, options%ntheta, fm, circ, &
                            res%bphi, status)
                        if (status /= 0) return
                        fm = fm - target
                        if (fa*fm <= 0) then
                            b = mid
                        else
                            a = mid; fa = fm
                        end if
                        if (b - a < 1.0e-13_dp*(u_hi - u_lo)) exit
                    end do
                    res%m = m; res%n = n; res%u = (a + b)/2
                    res%r = res%u*field%radial_scale
                    res%support = min(1.0_dp/abs(m), res%u - inner, outer - res%u)
                    if (res%support <= 0) cycle
                    ! The window spans at least two ring intervals: the
                    ! represented iota oscillates on the element scale.
                    h = min(max(options%shear_step, 2*ring_spacing(field, res%u)), &
                        0.5_dp*(res%u - u_lo), 0.5_dp*(u_hi - res%u))
                    if (h <= 0) cycle
                    ! Least-squares slope on 2 nfit + 1 points of [u - h, u + h]:
                    ! exact for quadratic iota, and it averages the element-scale
                    ! oscillation of the represented field.
                    ip = 0; im = 0
                    do it = -nfit, nfit
                        call seed_rotation(field, res%u + h*it/real(nfit, dp), &
                            options%ntheta, fm, circ, res%bphi, status)
                        if (status /= 0) return
                        ip = ip + it*fm; im = im + real(it, dp)**2
                    end do
                    res%shear = ip/(im*h/nfit)/field%radial_scale
                    call seed_rotation(field, res%u, options%ntheta, res%iota, &
                        res%circularity, res%bphi, status)
                    if (status /= 0) return
                    res%rmajor = field%centre(1)
                    if (field%fs%cyl_r0 > 0) res%rmajor = field%fs%cyl_r0
                    found = [found, res]
                end do
            end do
        end do
        status = 0
    end subroutine

    subroutine seed_circle_quadrature(field, r, nbase, theta, weight, status)
        ! Gauss-Legendre rule for the mean over theta on the circle of radius
        ! r about field%centre, split where the circle crosses element edges:
        ! represented fields are only piecewise smooth along the circle, and
        ! a uniform rule over their jumps converges at first order.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: r
        integer, intent(in) :: nbase
        real(dp), allocatable, intent(out) :: theta(:), weight(:)
        integer, intent(out) :: status
        integer, parameter :: ng = 6
        real(dp), parameter :: gx(ng) = [-0.932469514203152_dp, &
            -0.661209386466265_dp, -0.238619186083197_dp, 0.238619186083197_dp, &
            0.661209386466265_dp, 0.932469514203152_dp]
        real(dp), parameter :: gw(ng) = [0.171324492379170_dp, &
            0.360761573048139_dp, 0.467913934850384_dp, 0.467913934850384_dp, &
            0.360761573048139_dp, 0.171324492379170_dp]
        real(dp), allocatable :: cuts(:)
        real(dp) :: a, b, mid, half
        integer :: nb, j, ta, tb, tm, it, k, na
        nb = max(nbase, 8*maxval(field%mesh%nth))
        allocate (cuts(0))
        call cell_at(0.0_dp, ta, status)
        if (status /= 0) return
        do j = 1, nb
            b = 2*pi*j/nb
            call cell_at(b, tb, status)
            if (status /= 0) return
            if (tb /= ta) then
                a = 2*pi*(j - 1)/nb
                do it = 1, 60
                    mid = (a + b)/2
                    call cell_at(mid, tm, status)
                    if (status /= 0) return
                    if (tm == ta) then
                        a = mid
                    else
                        b = mid
                    end if
                    if (b - a < 1.0e-14_dp) exit
                end do
                cuts = [cuts, (a + b)/2]
                b = 2*pi*j/nb
            end if
            ta = tb
        end do
        if (size(cuts) == 0) cuts = [0.0_dp]
        na = size(cuts)
        allocate (theta(ng*na), weight(ng*na))
        do k = 1, na
            a = cuts(k)
            b = cuts(1) + 2*pi
            if (k < na) b = cuts(k + 1)
            half = (b - a)/2
            theta(ng*(k - 1) + 1:ng*k) = a + half*(1 + gx)
            weight(ng*(k - 1) + 1:ng*k) = half*gw/(2*pi)
        end do
        status = 0
    contains
        subroutine cell_at(angle, t, st)
            real(dp), intent(in) :: angle
            integer, intent(out) :: t, st
            real(dp) :: x(2), u, th, xi, eta
            x = field%centre + r*[cos(angle), sin(angle)]
            call field%coordinates(field%centre, field%radial_scale, x, u, th)
            call flux_mesh_locate(field%mesh, x, u, th, t, xi, eta, st)
        end subroutine
    end subroutine

    subroutine seed_rotation(field, u, ntheta, iota, circularity, bphi, status)
        ! iota = 1/<dphi/dtheta> on the circle of normalised radius u, the
        ! largest background radial-to-poloidal field ratio there, and <B_phi>.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: u
        integer, intent(in) :: ntheta
        real(dp), intent(out) :: iota, circularity, bphi
        integer, intent(out) :: status
        real(dp), allocatable :: theta(:), weight(:)
        real(dp) :: x(2), b(3), rot, radial, inverse, rmajor, rminor
        integer :: j
        iota = 0; circularity = 0; bphi = 0
        inverse = 0
        rminor = u*field%radial_scale
        call seed_circle_quadrature(field, rminor, ntheta, theta, weight, status)
        if (status /= 0) return
        do j = 1, size(theta)
            x = field%centre + rminor*[cos(theta(j)), sin(theta(j))]
            call background_point(field, x, b, status)
            if (status /= 0) return
            rmajor = x(1)
            if (field%fs%cyl_r0 > 0) rmajor = field%fs%cyl_r0
            call circle_rotation(rmajor, rminor, cos(theta(j)), sin(theta(j)), &
                b(1), b(3), b(2), rot, radial)
            status = 1
            if (abs(rot) <= tiny(1.0_dp)) return
            inverse = inverse + weight(j)/rot
            circularity = max(circularity, abs(radial)/max(abs(rot*rminor*b(2) &
                /rmajor), tiny(1.0_dp)))
            bphi = bphi + weight(j)*b(2)
        end do
        iota = 1/inverse
        status = 0
    end subroutine

    subroutine seed_resonant_br(field, res, ntheta, dbr, theta_o, status)
        ! Amplitude of the cos(m theta - n phi + phase) component of B_r on
        ! r_mn and the phi = 0 O-point angle of the pendulum model.
        type(free_b_field_t), intent(in) :: field
        type(seed_resonance_t), intent(in) :: res
        integer, intent(in) :: ntheta
        real(dp), intent(out) :: dbr, theta_o
        integer, intent(out) :: status
        real(dp), allocatable :: theta(:), weight(:)
        complex(dp) :: c, bn(3)
        real(dp) :: x(2), phase, chi
        integer :: j
        c = 0; dbr = 0; theta_o = 0
        call seed_circle_quadrature(field, res%r, ntheta, theta, weight, status)
        if (status /= 0) return
        do j = 1, size(theta)
            x = field%centre + res%r*[cos(theta(j)), sin(theta(j))]
            call mode_point(field, res%n, x, bn, status)
            if (status /= 0) return
            c = c + weight(j)*(cos(theta(j))*bn(1) + sin(theta(j))*bn(3)) &
                *exp(cmplx(0.0_dp, res%m*theta(j), dp))
        end do
        dbr = 2*abs(c)
        ! B_r = dbr cos(chi + phase), chi = m theta - n phi, phase = -arg c.
        ! dx/dphi = kappa cos(chi + phase), dchi/dphi = m iota' x: the
        ! O-point has sin(chi + phase) = sign(m iota').
        phase = -atan2(aimag(c), real(c, dp))
        chi = sign(pi/2, res%m*res%shear) - phase
        theta_o = modulo(chi/res%m, 2*pi/abs(res%m))
        status = 0
    end subroutine

    subroutine resonant_width(field, res, ntheta, dbr, width, theta_o, status)
        type(free_b_field_t), intent(in) :: field
        type(seed_resonance_t), intent(in) :: res
        integer, intent(in) :: ntheta
        real(dp), intent(out) :: dbr, width, theta_o
        integer, intent(out) :: status
        width = 0
        call seed_resonant_br(field, res, ntheta, dbr, theta_o, status)
        if (status /= 0) return
        status = 1
        if (abs(res%shear) <= tiny(1.0_dp) .or. abs(res%bphi) <= tiny(1.0_dp)) return
        call pendulum_width(res%rmajor, dbr, real(abs(res%m), dp), abs(res%shear), &
            abs(res%bphi), width)
        status = 0
    end subroutine

    subroutine background_point(field, x, b, status)
        ! Axisymmetric (n = 0) physical field (B_R, B_phi, B_Z) including the
        ! harmonic F/R e_phi.
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: b(3)
        integer, intent(out) :: status
        complex(dp) :: bn(3)
        real(dp) :: radius
        b = 0
        call mode_point(field, 0, x, bn, status)
        if (status /= 0) return
        radius = x(1)
        if (field%fs%cyl_r0 > 0) radius = field%fs%cyl_r0
        b = real(bn, dp)
        b(2) = b(2) + field%f_vac/radius
    end subroutine

    subroutine mode_point(field, n, x, bn, status)
        ! Complex amplitude of mode n of curl A at x, physical (R, phi, Z).
        type(free_b_field_t), intent(in) :: field
        integer, intent(in) :: n
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: bn(3)
        integer, intent(out) :: status
        real(dp) :: u, theta, xi, eta, y(2), v2(3, field%fs%sp(2)%nloc)
        integer :: t
        bn = 0; status = 1
        if (.not. all(ieee_is_finite(x))) return
        call field%coordinates(field%centre, field%radial_scale, x, u, theta)
        call flux_mesh_locate(field%mesh, x, u, theta, t, xi, eta, status)
        if (status /= 0) return
        call meridian_point_basis(field%fs, field%pb, 2, t, xi, eta, y, v2, status)
        if (status /= 0) return
        bn = matmul(v2, field%b(field%fs%sp(2)%dof(:, t), n) &
            *field%fs%sp(2)%sgn(:, t))
    end subroutine

    pure real(dp) function ring_spacing(field, u) result(h)
        type(free_b_field_t), intent(in) :: field
        real(dp), intent(in) :: u
        integer :: j
        h = field%mesh%u(1)
        do j = 2, field%mesh%nring
            if (field%mesh%u(j) >= u) then
                h = field%mesh%u(j) - field%mesh%u(j - 1)
                return
            end if
        end do
    end function

    pure recursive integer function gcd(a, b) result(d)
        integer, intent(in) :: a, b
        if (b == 0) then
            d = a
        else
            d = gcd(b, mod(a, b))
        end if
    end function

    subroutine symmetric_eigenvalues(a, w, info)
        real(dp), intent(inout) :: a(:, :)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: info
        interface
            subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
                import :: dp
                character, intent(in) :: jobz, uplo
                integer, intent(in) :: n, lda, lwork
                real(dp), intent(inout) :: a(lda, *)
                real(dp), intent(out) :: w(*), work(*)
                integer, intent(out) :: info
            end subroutine
        end interface
        real(dp) :: query(1)
        real(dp), allocatable :: work(:)
        integer :: n
        n = size(a, 1)
        call dsyev('N', 'U', n, a, n, w, query, -1, info)
        if (info /= 0) return
        allocate (work(max(1, int(query(1)))))
        call dsyev('N', 'U', n, a, n, w, work, size(work), info)
    end subroutine

    subroutine cholesky_solve(a, b, info)
        real(dp), intent(inout) :: a(:, :), b(:)
        integer, intent(out) :: info
        interface
            subroutine dpotrf(uplo, n, a, lda, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, lda
                real(dp), intent(inout) :: a(lda, *)
                integer, intent(out) :: info
            end subroutine
            subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
                import :: dp
                character, intent(in) :: uplo
                integer, intent(in) :: n, nrhs, lda, ldb
                real(dp), intent(in) :: a(lda, *)
                real(dp), intent(inout) :: b(ldb, *)
                integer, intent(out) :: info
            end subroutine
        end interface
        integer :: n
        n = size(a, 1)
        call dpotrf('L', n, a, n, info)
        if (info /= 0) return
        call dpotrs('L', n, 1, a, n, b, n, info)
    end subroutine
end module kin6d_island_seed
