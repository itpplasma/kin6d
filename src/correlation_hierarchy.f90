!> Classical correlation-hierarchy estimators of the periodic 3D OCP (H1a,
!> docs/correlation-hierarchy.md). Wigner-Seitz units: a = 1, n = 3/(4 pi),
!> beta v(r) = Gamma/r. Pair g2(r) and S(k) come from
!> kin6d_pair_correlation; this module adds
!>
!> - the triplet distribution g3(r12, r13, r23) on a cube of bins of width
!>   rmax/bins in each side (fully permutation-symmetric histogram), and its
!>   exact ideal-gas normalization: for one centre, the measure of ordered
!>   partner pairs is d3r2 d3r3 = 8 pi^2 r12 r13 r23 dr12 dr13 dr23 on the
!>   triangle domain, integrated per bin by a midpoint sub-grid in r12, r13
!>   (sub^2 points) and, in r23, 4-point Gauss-Legendre on the bin range
!>   clipped to the triangle inequality (the indicator is integrated
!>   exactly); the sub/2 grid gives the quadrature indicator;
!> - bin averages over the same measure of the Kirkwood product
!>   K = g12 g13 g23 and of the cumulant-truncation reference
!>   P = 1 + h12 + h13 + h23, from a fine g2 table; from these the connected
!>   triplet cumulant c3 = g3 - K (beyond-superposition part), the Ursell
!>   cumulant U3 = g3 - P, the reducible pair-product part K - P and the
!>   Kirkwood defect delta_K = g3/K - 1;
!> - the conditional mean force <F_1 . R-hat | r12 = R> (exact first YBG
!>   equation in the periodic box) and the YBG closure forces
!>   n int F(r13) X d3, X = g13 h23 (Kirkwood) or h23 (first order);
!> - pair-level observables of a g2 table: excess energy, screening moments
!>   and S(k) by radial Fourier transform.
module kin6d_correlation_hierarchy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use kin6d_correlation_hierarchy_leaves, only: ch_ybg_force_moment, &
        ch_excess_energy_density, ch_ocp_energy_fit, ch_ocp_energy_fit_jvp
    implicit none
    private
    public :: triplet_histogram_t, triplet_init, triplet_accumulate
    public :: triplet_distribution, triplet_references, triplet_norm
    public :: mean_force_accumulate, ybg_closure_force, ybg_potential
    public :: pair_excess_energy, pair_moments, pair_structure_factor, pair_value
    public :: CLOSURE_KIRKWOOD, CLOSURE_LINEAR, ocp_energy_literature

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: density = 3/(4*pi)
    integer, parameter :: CLOSURE_KIRKWOOD = 1, CLOSURE_LINEAR = 2
    !> Potekhin & Chabrier, Phys. Rev. E 62, 8554 (2000), fit of the
    !> Caillol (1999) and earlier MC data (A1, A2, B1, B2, B3, B4).
    real(dp), parameter :: pc2000(6) = [-0.907347_dp, 0.62849_dp, 4.50e-3_dp, &
                                        170.0_dp, -8.4e-5_dp, 3.70e-3_dp]

    type :: triplet_histogram_t
        integer :: bins = 0, sub = 0
        real(dp) :: rmax = 0
        real(dp), allocatable :: counts(:, :, :)
        !> 8 pi^2 int r1 r2 r3 over bin and triangle domain (sub and sub/2).
        real(dp), allocatable :: volume(:, :, :), volume_coarse(:, :, :)
        integer :: snapshots = 0
    end type
contains
    subroutine triplet_init(t, bins, rmax, sub)
        type(triplet_histogram_t), intent(out) :: t
        integer, intent(in) :: bins, sub
        real(dp), intent(in) :: rmax
        t%bins = bins; t%rmax = rmax; t%sub = max(2, 2*(sub/2))
        allocate (t%counts(bins, bins, bins), t%volume(bins, bins, bins), &
                  t%volume_coarse(bins, bins, bins))
        t%counts = 0
        call bin_measure(t, t%sub, t%volume)
        call bin_measure(t, t%sub/2, t%volume_coarse)
    end subroutine

    !> 8 pi^2 int_bin r1 r2 r3 [triangle] by a midpoint sub-grid.
    subroutine bin_measure(t, sub, volume)
        type(triplet_histogram_t), intent(in) :: t
        integer, intent(in) :: sub
        real(dp), intent(out) :: volume(:, :, :)
        integer :: a, b, c
        do c = 1, t%bins
            do b = 1, t%bins
                do a = 1, t%bins
                    volume(a, b, c) = 8*pi**2*bin_average(t, sub, a, b, c, one)
                end do
            end do
        end do
    contains
        real(dp) function one(r1, r2, r3)
            real(dp), intent(in) :: r1, r2, r3
            one = 1 + 0*(r1 + r2 + r3)
        end function
    end subroutine

    !> int_bin w f with w = r1 r2 r3 [triangle] (see module header).
    real(dp) function bin_average(t, sub, a, b, c, f) result(s)
        type(triplet_histogram_t), intent(in) :: t
        integer, intent(in) :: sub, a, b, c
        interface
            real(dp) function f(r1, r2, r3)
                import :: dp
                real(dp), intent(in) :: r1, r2, r3
            end function
        end interface
        real(dp), parameter :: x4(4) = [-0.8611363115940526_dp, -0.3399810435848563_dp, &
            0.3399810435848563_dp, 0.8611363115940526_dp]
        real(dp), parameter :: w4(4) = [0.3478548451374538_dp, 0.6521451548625461_dp, &
            0.6521451548625461_dp, 0.3478548451374538_dp]
        real(dp) :: h, d, r1, r2, r3, lo, hi
        integer :: i, j, k
        h = t%rmax/t%bins; d = h/sub; s = 0
        do j = 1, sub
            r2 = (b - 1)*h + (j - 0.5_dp)*d
            do i = 1, sub
                r1 = (a - 1)*h + (i - 0.5_dp)*d
                ! r3 over the bin clipped to the triangle range [|r1-r2|, r1+r2]:
                ! the indicator is integrated exactly, 4-point Gauss-Legendre
                ! for r3 f (exact for the measure).
                lo = max((c - 1)*h, abs(r1 - r2)); hi = min(c*h, r1 + r2)
                if (hi <= lo) cycle
                do k = 1, 4
                    r3 = (lo + hi)/2 + (hi - lo)/2*x4(k)
                    s = s + w4(k)*(hi - lo)/2/d*r1*r2*r3*f(r1, r2, r3)
                end do
            end do
        end do
        s = s*d**3
    end function

    !> Add all triplets of one snapshot x(3, n) with every side < rmax.
    subroutine triplet_accumulate(t, x, n, length)
        type(triplet_histogram_t), intent(inout) :: t
        real(dp), intent(in) :: x(:, :), length
        integer, intent(in) :: n
        real(dp) :: d(3, n), r(n), e(3), r23, scale
        integer :: idx(n), i, j, k, m, b1, b2, b3
        scale = t%bins/t%rmax
        do i = 1, n
            m = 0
            do j = 1, n
                if (j == i) cycle
                e = x(:, j) - x(:, i)
                e = e - length*anint(e/length)
                if (sum(e*e) >= t%rmax**2) cycle
                m = m + 1
                d(:, m) = e; r(m) = sqrt(sum(e*e)); idx(m) = 1 + int(r(m)*scale)
            end do
            do j = 1, m - 1
                b1 = min(t%bins, idx(j))
                do k = j + 1, m
                    e = d(:, j) - d(:, k)
                    r23 = sqrt(sum(e*e))
                    if (r23 >= t%rmax) cycle
                    b2 = min(t%bins, idx(k)); b3 = min(t%bins, 1 + int(r23*scale))
                    t%counts(b1, b2, b3) = t%counts(b1, b2, b3) + 1
                    t%counts(b2, b1, b3) = t%counts(b2, b1, b3) + 1
                end do
            end do
        end do
        t%snapshots = t%snapshots + 1
    end subroutine

    !> g3 per bin from counts: counts/(snapshots N (N-1)(N-2)/V^2 volume);
    !> 0 where the bin misses the triangle domain.
    subroutine triplet_distribution(t, counts, snapshots, n, length, g3)
        type(triplet_histogram_t), intent(in) :: t
        real(dp), intent(in) :: counts(:, :, :), length
        integer, intent(in) :: snapshots, n
        real(dp), intent(out) :: g3(:, :, :)
        real(dp) :: ideal
        ideal = real(snapshots, dp)*n*real(n - 1, dp)*(n - 2)/length**6
        where (t%volume > 0)
            g3 = counts/(ideal*t%volume)
        elsewhere
            g3 = 0
        end where
    end subroutine

    !> Bin averages (same measure as g3) of K = g12 g13 g23 and
    !> P = 1 + h12 + h13 + h23 from the g2 table (rg uniform bin centres),
    !> and K at the bin centre (discretization indicator).
    subroutine triplet_references(t, rg, g2, kirkwood, linear, centre)
        type(triplet_histogram_t), intent(in) :: t
        real(dp), intent(in) :: rg(:), g2(:)
        real(dp), intent(out) :: kirkwood(:, :, :), linear(:, :, :), centre(:, :, :)
        integer :: a, b, c
        real(dp) :: h, w
        h = t%rmax/t%bins
        do c = 1, t%bins
            do b = 1, t%bins
                do a = 1, t%bins
                    kirkwood(a, b, c) = 0; linear(a, b, c) = 0
                    centre(a, b, c) = pair_value(rg, g2, (a - 0.5_dp)*h)* &
                                      pair_value(rg, g2, (b - 0.5_dp)*h)* &
                                      pair_value(rg, g2, (c - 0.5_dp)*h)
                    if (t%volume(a, b, c) <= 0) cycle
                    w = t%volume(a, b, c)/(8*pi**2)
                    kirkwood(a, b, c) = bin_average(t, t%sub, a, b, c, product3)/w
                    linear(a, b, c) = bin_average(t, t%sub, a, b, c, sum3)/w
                end do
            end do
        end do
    contains
        real(dp) function product3(r1, r2, r3)
            real(dp), intent(in) :: r1, r2, r3
            product3 = pair_value(rg, g2, r1)*pair_value(rg, g2, r2)* &
                       pair_value(rg, g2, r3)
        end function
        real(dp) function sum3(r1, r2, r3)
            real(dp), intent(in) :: r1, r2, r3
            sum3 = pair_value(rg, g2, r1) + pair_value(rg, g2, r2) + &
                   pair_value(rg, g2, r3) - 2
        end function
    end subroutine

    !> Squared size n^2 int int X^2 d3r2 d3r3 over the triplet domain
    !> (all sides < rmax), X given per bin.
    real(dp) function triplet_norm(t, x) result(s)
        type(triplet_histogram_t), intent(in) :: t
        real(dp), intent(in) :: x(:, :, :)
        s = density**2*sum(x**2*t%volume)
    end function

    !> Linear interpolation of a pair table at r (1 beyond the table; the
    !> first value below the first centre).
    pure real(dp) function pair_value(rg, g, r) result(v)
        real(dp), intent(in) :: rg(:), g(:), r
        real(dp) :: dr, s
        integer :: k
        dr = rg(2) - rg(1)
        if (r >= rg(size(rg))) then
            v = 1
            return
        end if
        if (r <= rg(1)) then
            v = g(1)
            return
        end if
        s = (r - rg(1))/dr
        k = min(size(rg) - 1, 1 + int(s))
        s = s - (k - 1)
        v = (1 - s)*g(k) + s*g(k + 1)
    end function

    !> Add (F_i - F_j).R-hat/2 for every minimum-image pair with r < rmax
    !> into bins of width rmax/size(sums); R-hat = (x_i - x_j)/r.
    subroutine mean_force_accumulate(x, f, n, length, rmax, sums)
        real(dp), intent(in) :: x(:, :), f(:, :), length, rmax
        integer, intent(in) :: n
        real(dp), intent(inout) :: sums(:)
        integer :: i, j, k
        real(dp) :: d(3), r, scale
        scale = size(sums)/rmax
        do i = 1, n - 1
            do j = i + 1, n
                d = x(:, i) - x(:, j)
                d = d - length*anint(d/length)
                r = sqrt(sum(d*d))
                if (r >= rmax) cycle
                k = 1 + int(r*scale)
                sums(k) = sums(k) + sum((f(:, i) - f(:, j))*d)/(2*r)
            end do
        end do
    end subroutine

    !> YBG closure part I(R) (units k_B T/a) at the table centres rg:
    !> I(R) = 2 pi n int_0^smax ds int_-1^1 dmu [s^2 beta F.R-hat] X(s, t),
    !> t = |s - R|, X = g(s) h(t) (Kirkwood) or h(t) (first order); h = 0
    !> and g = 1 beyond the table. Midpoint rule in s (width ds), Gauss-
    !> Legendre in mu (nmu points).
    subroutine ybg_closure_force(rg, g2, gamma, closure, ds, nmu, force)
        real(dp), intent(in) :: rg(:), g2(:), gamma, ds
        integer, intent(in) :: closure, nmu
        real(dp), intent(out) :: force(:)
        real(dp) :: mu(nmu), wmu(nmu), s, t, big, moment, x, rmax
        integer :: k, i, j, ns
        rmax = rg(size(rg)) + (rg(2) - rg(1))/2
        call gauss_legendre(nmu, mu, wmu)
        do k = 1, size(rg)
            big = rg(k) + rmax
            ns = ceiling(big/ds)
            force(k) = 0
            do i = 1, ns
                s = (i - 0.5_dp)*big/ns
                x = 0
                do j = 1, nmu
                    t = sqrt(max(0.0_dp, s*s + rg(k)**2 - 2*s*rg(k)*mu(j)))
                    call ch_ybg_force_moment(gamma, s, mu(j), moment)
                    x = x + wmu(j)*moment*(pair_value(rg, g2, t) - 1)
                end do
                if (closure == CLOSURE_KIRKWOOD) x = x*pair_value(rg, g2, s)
                force(k) = force(k) + x*big/ns
            end do
            force(k) = 2*pi*density*force(k)
        end do
    end subroutine

    !> Potential of mean force w(r_k) = int_{r_k}^{r_last} force from the
    !> mean-force table on the centres (trapezoid), w(r_last) = 0; with
    !> `coulomb`, the bare part Gamma/r^2 is added analytically.
    subroutine ybg_potential(rg, force, gamma, coulomb, w)
        real(dp), intent(in) :: rg(:), force(:), gamma
        logical, intent(in) :: coulomb
        real(dp), intent(out) :: w(:)
        integer :: k, m
        m = size(rg)
        w(m) = 0
        do k = m - 1, 1, -1
            w(k) = w(k + 1) + (force(k) + force(k + 1))*(rg(k + 1) - rg(k))/2
        end do
        if (coulomb) w = w + gamma*(1/rg - 1/rg(m))
    end subroutine

    !> beta u_ex = int 2 pi n r^2 beta v h dr on the table (midpoint),
    !> without tail beyond the table.
    real(dp) function pair_excess_energy(rg, g2, gamma) result(u)
        real(dp), intent(in) :: rg(:), g2(:), gamma
        real(dp) :: e(size(rg))
        call ch_excess_energy_density(rg, gamma, g2 - 1, density, e)
        u = sum(e)*(rg(2) - rg(1))
    end function

    !> Perfect-screening moment m0 = n int h d3r (-> -1) and the
    !> Stillinger-Lovett moment m2 = n int r^2 h d3r kappa^2/6 (-> -1),
    !> kappa^2 = 3 Gamma, on the table.
    subroutine pair_moments(rg, g2, gamma, m0, m2)
        real(dp), intent(in) :: rg(:), g2(:), gamma
        real(dp), intent(out) :: m0, m2
        real(dp) :: dr
        dr = rg(2) - rg(1)
        m0 = density*sum(4*pi*rg**2*(g2 - 1))*dr
        m2 = density*sum(4*pi*rg**4*(g2 - 1))*dr*3*gamma/6
    end subroutine

    !> S(k) = 1 + n int h sin(kr)/(kr) 4 pi r^2 dr on the table.
    real(dp) function pair_structure_factor(rg, g2, k) result(s)
        real(dp), intent(in) :: rg(:), g2(:), k
        s = 1 + density*sum(4*pi*rg**2*(g2 - 1)*sin(k*rg)/(k*rg))*(rg(2) - rg(1))
    end function

    !> Literature OCP excess energy u = beta U_ex/N and du/dGamma (fit leaf
    !> and its generated JVP).
    elemental subroutine ocp_energy_literature(gamma, u, dudg)
        real(dp), intent(in) :: gamma
        real(dp), intent(out) :: u, dudg
        real(dp), parameter :: z = 0
        call ch_ocp_energy_fit(gamma, pc2000(1), pc2000(2), pc2000(3), pc2000(4), &
                               pc2000(5), pc2000(6), u)
        call ch_ocp_energy_fit_jvp(gamma, pc2000(1), pc2000(2), pc2000(3), pc2000(4), &
                                   pc2000(5), pc2000(6), 1.0_dp, z, z, z, z, z, z, dudg)
    end subroutine
end module kin6d_correlation_hierarchy
