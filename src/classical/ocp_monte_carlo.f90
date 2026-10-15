!> Canonical Metropolis Monte Carlo of N electrons in a periodic cube with a
!> uniform neutralizing background: the one classical periodic-cell sampler
!> of the correlation maps H1a (OCP, docs/correlation-hierarchy.md) and H1b
!> (finite-N classical parent, docs/quantum-classical-correlations.md).
!>
!> mc_setup_cell takes the cell edge L, beta = 1/kT and the Coulomb coupling
!> g = e^2/(4 pi eps0) in any consistent units, and a kin6d_pair_law
!> selector. mc_setup is the OCP form: Wigner-Seitz radius a = 1, g = 1, so
!> the reduced energy U is in e^2/(4 pi eps0 a) and beta = Gamma;
!> n = 3/(4 pi), L = (4 pi N/3)^(1/3). Gamma is checked against the shared
!> kin6d_species scales (one definition).
!>
!> Energy: Ewald sum of kin6d_coulomb_leaves with parameters from
!> kin6d_ewald (real space: minimum image for real_cutoff = 0, the OCP
!> choice, else the image cube |n|_inf <= real_cutoff of the minimum-image
!> separation, with the self-image constant; reciprocal cube
!> 0 < |n|_inf <= nk, tin-foil), self and background constants, plus the
!> declared short-range pair addition of kin6d_pair_law (Kelbg/Pauli) per
!> minimum-image pair. The Ewald truncation bound of kin6d_ewald holds per
!> pair in units g/L. A single-particle move updates the energy exactly for
!> this truncated Hamiltonian: O(N) real-space differences and an O(modes)
!> update of the structure factors rho_k = sum_j exp(i k.x_j). The chain
!> samples exp(-beta U) of that Hamiltonian; its distance from the
!> untruncated Ewald energy is the reported truncation bound.
!>
!> Moves: particle chosen uniformly, displacement uniform in [-d, d]^3,
!> wrapped into the cube; Metropolis acceptance. The step d (initially a/2)
!> may be tuned only during equilibration (mc_sweep tune=.true.); production
!> sweeps keep it fixed so that detailed balance holds.
module kin6d_ocp_monte_carlo
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use kin6d_coulomb_leaves, only: ewald_real_pair, ewald_reciprocal_weight, &
                                    ewald_self, ewald_background
    use kin6d_correlation_hierarchy_leaves, only: ch_ewald_real_potential
    use kin6d_ewald, only: ewald_parameters_t, ewald_select
    use kin6d_species, only: atomic_units, electron, plasma_scales, plasma_scales_t
    use kin6d_pair_law, only: pair_law_t, pair_law_excess, pair_coulomb
    implicit none
    private
    public :: ocp_mc_t, mc_setup, mc_setup_cell, mc_sweep, mc_energy, mc_forces, mc_acceptance
    public :: mc_wigner_seitz_length

    real(dp), parameter :: pi = 3.14159265358979323846_dp, two_pi = 2*pi

    type :: ocp_mc_t
        integer :: n = 0, nk = 0, modes = 0
        !> beta = 1/kT, Coulomb coupling g, Gamma = beta g/a.
        real(dp) :: gamma = 0, beta = 0, coupling = 1
        real(dp) :: length = 0, volume = 0, step = 0.5_dp
        type(ewald_parameters_t) :: ewald
        type(plasma_scales_t) :: scales
        type(pair_law_t) :: law
        !> Cached pair additions extra(i, j) of the pair law (zero for Coulomb).
        real(dp), allocatable :: extra(:, :), extra_new(:)
        !> Positions x(3, n) in [0, L)^3.
        real(dp), allocatable :: x(:, :)
        integer, allocatable :: wave(:, :)
        !> Half-space reciprocal weights (doubled for +/-k), coupling g.
        real(dp), allocatable :: weight(:), kvec(:, :)
        complex(dp), allocatable :: rho(:), e_old(:), e_new(:), axis(:, :)
        real(dp), allocatable :: r_old(:), r_new(:), u_old(:), u_new(:)
        !> Current energy (units of g/length) and its configuration constant.
        real(dp) :: energy = 0, constant = 0
        integer(int64) :: attempted = 0, accepted = 0
        type(rng_t) :: rng
    end type
contains
    !> Box edge of N particles at n = 3/(4 pi) (a = 1).
    pure real(dp) function mc_wigner_seitz_length(n) result(length)
        integer, intent(in) :: n
        length = (4*pi*n/3)**(1.0_dp/3)
    end function

    !> OCP chain of n electrons at coupling gamma (a = 1, g = 1, bare
    !> Coulomb), uniform random start.
    !> tolerance: relative Ewald truncation tolerance (kin6d_ewald units).
    subroutine mc_setup(mc, n, gamma, tolerance, seed, ok)
        type(ocp_mc_t), intent(out) :: mc
        integer, intent(in) :: n
        real(dp), intent(in) :: gamma, tolerance
        integer(int64), intent(in) :: seed
        logical, intent(out) :: ok
        type(plasma_scales_t) :: scales
        ok = .false.
        if (n < 2 .or. .not. (gamma > 0) .or. .not. ieee_is_finite(gamma)) return
        ! Gamma = e^2/(4 pi eps0 a k_B T) from the shared scale definitions.
        scales = plasma_scales(atomic_units(), electron(), 3/(4*pi), 1/gamma)
        if (abs(scales%gamma/gamma - 1) > 1.0e-12_dp .or. &
            abs(scales%wigner_seitz - 1) > 1.0e-12_dp) return
        call mc_setup_cell(mc, n, mc_wigner_seitz_length(n), gamma, 1.0_dp, &
                           tolerance, seed, ok, real_cutoff=0)
        ! Exact reduced values (a = 1 up to rounding in mc_setup_cell).
        mc%scales = scales; mc%gamma = gamma; mc%step = 0.5_dp
    end subroutine

    !> Chain of n electrons in the cube of edge length at beta = 1/kT with
    !> Coulomb coupling g; law (default bare Coulomb) adds the declared
    !> short-range pair term and must describe the same cell and coupling.
    !> real_cutoff as in ewald_select: 0 is the minimum-image real sum (large
    !> N); absent selects the cheapest real image cube 1..3 by the pair-count
    !> cost model (small N, where minimum image needs many reciprocal modes).
    subroutine mc_setup_cell(mc, n, length, beta, coupling, tolerance, seed, ok, law, &
                             real_cutoff)
        type(ocp_mc_t), intent(out) :: mc
        integer, intent(in) :: n
        real(dp), intent(in) :: length, beta, coupling, tolerance
        integer(int64), intent(in) :: seed
        logical, intent(out) :: ok
        type(pair_law_t), intent(in), optional :: law
        integer, intent(in), optional :: real_cutoff
        type(fortnum_status_t) :: status
        integer :: i, j, k, m, nk
        real(dp) :: k2, u
        ok = .false.
        if (n < 2 .or. .not. (length > 0 .and. beta > 0 .and. coupling > 0)) return
        if (.not. (ieee_is_finite(length) .and. ieee_is_finite(beta))) return
        if (present(law)) mc%law = law
        if (mc%law%model /= pair_coulomb .and. &
            (abs(mc%law%length/length - 1) > 1.0e-12_dp .or. &
             abs(mc%law%coupling/coupling - 1) > 1.0e-12_dp)) return
        mc%n = n; mc%beta = beta; mc%coupling = coupling
        mc%length = length; mc%volume = length**3
        mc%gamma = beta*coupling/(3*mc%volume/(4*pi*n))**(1.0_dp/3)
        mc%step = (3*mc%volume/(4*pi*n))**(1.0_dp/3)/2
        call ewald_select(mc%length, tolerance, mc%ewald, ok, particles=n, &
                          real_cutoff=real_cutoff)
        if (.not. ok) return
        ok = .false.
        nk = mc%ewald%reciprocal_cutoff; mc%nk = nk
        mc%modes = ((2*nk + 1)**3 - 1)/2
        allocate (mc%wave(3, mc%modes), mc%weight(mc%modes), mc%kvec(3, mc%modes))
        m = 0
        do i = 0, nk
            do j = -nk, nk
                do k = -nk, nk
                    if (i == 0 .and. (j < 0 .or. (j == 0 .and. k <= 0))) cycle
                    m = m + 1
                    mc%wave(:, m) = [i, j, k]
                    mc%kvec(:, m) = two_pi/mc%length*real([i, j, k], dp)
                    k2 = sum(mc%kvec(:, m)**2)
                    call ewald_reciprocal_weight(k2, mc%ewald%alpha, mc%volume, coupling, &
                                                 mc%weight(m))
                    mc%weight(m) = 2*mc%weight(m)
                end do
            end do
        end do
        allocate (mc%rho(mc%modes), mc%e_old(mc%modes), mc%e_new(mc%modes))
        allocate (mc%axis(-nk:nk, 3), mc%x(3, n))
        allocate (mc%r_old(n), mc%r_new(n), mc%u_old(n), mc%u_new(n))
        allocate (mc%extra(n, n), mc%extra_new(n))
        mc%extra = 0; mc%extra_new = 0
        call ewald_self(mc%ewald%alpha, coupling, u)
        mc%constant = n*u
        call ewald_background(mc%ewald%alpha, mc%volume, coupling*real(n, dp)**2, u)
        mc%constant = mc%constant + u/2
        ! Real-space interaction of every charge with its own images.
        mc%constant = mc%constant + n*coupling*image_sum(mc, [0.0_dp, 0.0_dp, 0.0_dp])/2
        call rng_seed(mc%rng, seed, status)
        if (status%code /= FORTNUM_OK) return
        do i = 1, n
            do j = 1, 3
                call rng_uniform(mc%rng, u)
                mc%x(j, i) = mc%length*u
            end do
        end do
        mc%energy = mc_energy(mc)
        ok = ieee_is_finite(mc%energy)
    end subroutine

    !> exp(i k.p) for every half-space mode.
    subroutine phases(mc, p, e)
        type(ocp_mc_t), intent(inout) :: mc
        real(dp), intent(in) :: p(3)
        complex(dp), intent(out) :: e(:)
        integer :: a, m
        real(dp) :: theta
        do a = 1, 3
            theta = two_pi*p(a)/mc%length
            do m = 0, mc%nk
                mc%axis(m, a) = cmplx(cos(m*theta), sin(m*theta), dp)
                mc%axis(-m, a) = conjg(mc%axis(m, a))
            end do
        end do
        do m = 1, mc%modes
            e(m) = mc%axis(mc%wave(1, m), 1)*mc%axis(mc%wave(2, m), 2)* &
                   mc%axis(mc%wave(3, m), 3)
        end do
    end subroutine

    !> Full reduced energy of the current configuration (also resets rho).
    real(dp) function mc_energy(mc) result(energy)
        type(ocp_mc_t), intent(inout) :: mc
        integer :: i, j
        real(dp) :: d(3)
        energy = mc%constant
        do i = 1, mc%n - 1
            do j = i + 1, mc%n
                d = mc%x(:, i) - mc%x(:, j)
                d = d - mc%length*anint(d/mc%length)
                mc%r_old(j) = sqrt(sum(d*d))
            end do
            if (any(mc%r_old(i + 1:mc%n) <= 0)) then
                energy = huge(1.0_dp)
                return
            end if
            call ch_ewald_real_potential(mc%r_old(i + 1:mc%n), mc%ewald%alpha, &
                                         mc%u_old(i + 1:mc%n))
            if (mc%ewald%real_cutoff > 0) then
                do j = i + 1, mc%n
                    d = mc%x(:, i) - mc%x(:, j)
                    mc%u_old(j) = mc%u_old(j) + image_sum(mc, d - mc%length*anint(d/mc%length))
                end do
            end if
            energy = energy + mc%coupling*sum(mc%u_old(i + 1:mc%n))
            if (mc%law%model == pair_coulomb) cycle
            do j = i + 1, mc%n
                d = mc%x(:, i) - mc%x(:, j)
                mc%extra(i, j) = pair_law_excess(mc%law, d - mc%length*anint(d/mc%length))
                mc%extra(j, i) = mc%extra(i, j)
            end do
            energy = energy + sum(mc%extra(i, i + 1:mc%n))
        end do
        mc%rho = 0
        do i = 1, mc%n
            call phases(mc, mc%x(:, i), mc%e_old)
            mc%rho = mc%rho + mc%e_old
        end do
        energy = energy + sum(mc%weight*(real(mc%rho)**2 + aimag(mc%rho)**2))/2
    end function

    !> One sweep of n single-particle Metropolis moves; with tune, the step
    !> is adapted towards 40% acceptance (equilibration only).
    subroutine mc_sweep(mc, tune)
        type(ocp_mc_t), intent(inout) :: mc
        logical, intent(in) :: tune
        integer :: move, i, a, accepted
        real(dp) :: u, trial(3), d(3), du, rate
        integer :: j
        accepted = 0
        do move = 1, mc%n
            call rng_uniform(mc%rng, u)
            i = min(mc%n, 1 + int(u*mc%n))
            do a = 1, 3
                call rng_uniform(mc%rng, u)
                trial(a) = mc%x(a, i) + mc%step*(2*u - 1)
                trial(a) = trial(a) - mc%length*floor(trial(a)/mc%length)
            end do
            do j = 1, mc%n
                d = mc%x(:, j) - trial
                d = d - mc%length*anint(d/mc%length)
                mc%r_new(j) = sqrt(sum(d*d))
                d = mc%x(:, j) - mc%x(:, i)
                d = d - mc%length*anint(d/mc%length)
                mc%r_old(j) = sqrt(sum(d*d))
            end do
            mc%r_new(i) = 1; mc%r_old(i) = 1
            mc%attempted = mc%attempted + 1
            ! A coincident trial position has infinite energy: reject.
            if (any(mc%r_new <= 0)) cycle
            call ch_ewald_real_potential(mc%r_new, mc%ewald%alpha, mc%u_new)
            call ch_ewald_real_potential(mc%r_old, mc%ewald%alpha, mc%u_old)
            if (mc%ewald%real_cutoff > 0) then
                do j = 1, mc%n
                    if (j == i) cycle
                    d = mc%x(:, j) - trial
                    mc%u_new(j) = mc%u_new(j) + image_sum(mc, d - mc%length*anint(d/mc%length))
                    d = mc%x(:, j) - mc%x(:, i)
                    mc%u_old(j) = mc%u_old(j) + image_sum(mc, d - mc%length*anint(d/mc%length))
                end do
            end if
            du = mc%coupling*(sum(mc%u_new) - sum(mc%u_old))
            if (mc%law%model /= pair_coulomb) then
                do j = 1, mc%n
                    if (j == i) cycle
                    d = trial - mc%x(:, j)
                    mc%extra_new(j) = pair_law_excess(mc%law, d - mc%length*anint(d/mc%length))
                    du = du + mc%extra_new(j) - mc%extra(i, j)
                end do
            end if
            call phases(mc, mc%x(:, i), mc%e_old)
            call phases(mc, trial, mc%e_new)
            mc%e_new = mc%e_new - mc%e_old
            du = du + sum(mc%weight*(real(conjg(mc%rho)*mc%e_new) + &
                      (real(mc%e_new)**2 + aimag(mc%e_new)**2)/2))
            call rng_uniform(mc%rng, u)
            if (du <= 0 .or. u < exp(-mc%beta*du)) then
                mc%x(:, i) = trial
                if (mc%law%model /= pair_coulomb) then
                    mc%extra_new(i) = 0
                    mc%extra(i, :) = mc%extra_new; mc%extra(:, i) = mc%extra_new
                end if
                mc%rho = mc%rho + mc%e_new
                mc%energy = mc%energy + du
                mc%accepted = mc%accepted + 1
                accepted = accepted + 1
            end if
        end do
        if (tune) then
            rate = real(accepted, dp)/mc%n
            mc%step = min(mc%length/2, mc%step*max(0.5_dp, min(2.0_dp, &
                          (rate + 0.05_dp)/0.45_dp)))
        end if
    end subroutine

    !> Unit-coupling real-space Ewald sum over the images 0 < |n|_inf <= nr
    !> of the minimum-image separation d (zero for the minimum-image sum).
    real(dp) function image_sum(mc, d) result(u)
        type(ocp_mc_t), intent(in) :: mc
        real(dp), intent(in) :: d(3)
        real(dp) :: r, v
        integer :: a, b, c, nr
        u = 0; nr = mc%ewald%real_cutoff
        do a = -nr, nr
            do b = -nr, nr
                do c = -nr, nr
                    if (a == 0 .and. b == 0 .and. c == 0) cycle
                    r = norm2(d + mc%length*[a, b, c])
                    call ch_ewald_real_potential(r, mc%ewald%alpha, v)
                    u = u + v
                end do
            end do
        end do
    end function

    !> Cumulative acceptance ratio.
    real(dp) function mc_acceptance(mc) result(rate)
        type(ocp_mc_t), intent(in) :: mc
        rate = real(mc%accepted, dp)/max(1_int64, mc%attempted)
    end function

    !> Total periodic force f(3, n) (units g/length^2) on every particle,
    !> the exact gradient of the sampled truncated Ewald Hamiltonian
    !> (real-space pairs plus reciprocal modes; background exerts none).
    !> Defined for the bare Coulomb pair law only (ok = .false. otherwise).
    subroutine mc_forces(mc, f, ok)
        type(ocp_mc_t), intent(inout) :: mc
        real(dp), intent(out) :: f(:, :)
        logical, intent(out) :: ok
        integer :: i, j, m, a, b, c, nr
        real(dp) :: d(3), g(3), u, t
        f(:, 1:mc%n) = 0
        nr = mc%ewald%real_cutoff
        ok = mc%law%model == pair_coulomb
        if (.not. ok) return
        do i = 1, mc%n - 1
            do j = i + 1, mc%n
                d = mc%x(:, i) - mc%x(:, j)
                d = d - mc%length*anint(d/mc%length)
                do a = -nr, nr
                    do b = -nr, nr
                        do c = -nr, nr
                            call ewald_real_pair(d(1) + mc%length*a, d(2) + mc%length*b, &
                                                 d(3) + mc%length*c, mc%ewald%alpha, &
                                                 mc%coupling, g(1), g(2), g(3), u)
                            f(:, i) = f(:, i) + g
                            f(:, j) = f(:, j) - g
                        end do
                    end do
                end do
            end do
        end do
        do i = 1, mc%n
            call phases(mc, mc%x(:, i), mc%e_old)
            do m = 1, mc%modes
                t = mc%weight(m)*aimag(conjg(mc%e_old(m))*mc%rho(m))
                f(:, i) = f(:, i) - t*mc%kvec(:, m)
            end do
        end do
    end subroutine
end module kin6d_ocp_monte_carlo
