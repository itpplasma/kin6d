!> Finite-N ideal spin-1/2 fermions (campaign case Q3): canonical (fixed-N)
!> partition sums, occupations, energy, pressure, one-body density matrix and
!> pair correlations for a periodic box or a harmonic trap, without Coulomb
!> interaction. Spin species are independent ensembles sharing a spectrum.
!> Canonical sums use the subtraction-free elementary-symmetric recursion
!> Z_k <- Z_k + y*Z_(k-1) per orbital, in binary64 for the candidate and in
!> outward FortNum intervals for the enclosure, including a rigorous bound on
!> the omitted orbitals above the spectral cutoff. The alternating power-sum
!> recursion is supplied only as a binary128 cross-check with its measured
!> cancellation factor. See docs/fermion-ensemble.md.
module kin6d_fermion_ensemble
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, hull, is_empty, mag, &
        operator(+), operator(-), operator(*), operator(/), log
    use fortnum_quadrature, only: gauss_legendre_ab
    use kin6d_fermi_gas, only: fermi_state_t, occupation_fermi_gas
    use kin6d_fermion_ensemble_leaves, only: fermion_weight_jet, &
        fermion_weight_jet_interval, fermion_box_fermi_energy, &
        fermion_box_fermi_energy_interval, fermion_inverse_temperature, &
        fermion_box_eos, fermion_box_eos_interval, fermion_box_tail, &
        fermion_box_tail_margin_interval, fermion_trap_tail, &
        fermion_pair_occupation, fermion_sinc, fermion_sinc_series, &
        fermion_lattice_argument, fermion_density_matrix_term, &
        fermion_grand_positive, fermion_grand_negative
    implicit none
    private
    public :: fermion_spectrum_t, canonical_state_t
    public :: box_spectrum, trap_spectrum, box_cutoff, box_parameters
    public :: solve_canonical, solve_boltzmann, solve_ground, box_eos
    public :: same_spin_pair, one_body_density, opposite_spin_pair
    public :: grand_canonical_density, power_sum_recursion
    integer, parameter, public :: box_kind = 1, trap_kind = 2

    !> Shell-ordered single-particle spectrum. A closed spectrum is the finite
    !> orbital model itself; otherwise the omitted orbitals above the cutoff
    !> belong to the model and enter every enclosure through a remainder bound.
    type :: fermion_spectrum_t
        integer :: kind = 0, cutoff = -1, shells = 0, orbitals = 0
        logical :: closed = .false.
        real(dp), allocatable :: level(:)
        integer, allocatable :: degeneracy(:), first(:)
        integer, allocatable :: vectors(:, :)
    end type fermion_spectrum_t

    !> One spin species with a fixed particle number. Energies in spectrum
    !> units; occupation(a) is the mean occupation of each orbital of shell a.
    type :: canonical_state_t
        integer :: particles = 0
        real(dp) :: b = 0.0_dp, shift = 0.0_dp
        real(dp) :: log_z = 0.0_dp, energy = 0.0_dp, variance = 0.0_dp
        real(dp), allocatable :: occupation(:), weight(:), same_shell_pair(:)
        type(interval_t) :: log_z_bound, energy_bound, variance_bound
        type(interval_t) :: omitted_bound
        type(interval_t), allocatable :: occupation_bound(:)
        logical :: valid = .false., certified = .false., boltzmann = .false.
        logical :: ground = .false.
        character(48) :: status = 'uninitialized'
    end type canonical_state_t

contains

    !> Periodic box orbitals n in Z**3 with |n|**2<=cutoff, grouped by |n|**2.
    subroutine box_spectrum(cutoff, spectrum, ok, closed)
        integer, intent(in) :: cutoff
        type(fermion_spectrum_t), intent(out) :: spectrum
        logical, intent(out) :: ok
        logical, intent(in), optional :: closed
        integer, allocatable :: count(:), position(:)
        integer :: m, i, j, k, q, total, shell
        ok = .false.
        if (cutoff < 0 .or. cutoff > 10**6) return
        m = int(sqrt(real(cutoff, dp))) + 1
        allocate (count(0:cutoff), source=0)
        do i = -m, m
            do j = -m, m
                do k = -m, m
                    q = i*i + j*j + k*k
                    if (q <= cutoff) count(q) = count(q) + 1
                end do
            end do
        end do
        spectrum%kind = box_kind
        spectrum%cutoff = cutoff
        if (present(closed)) spectrum%closed = closed
        spectrum%shells = count_nonzero(count)
        total = sum(count)
        spectrum%orbitals = total
        allocate (spectrum%level(spectrum%shells), spectrum%degeneracy(spectrum%shells), &
                  spectrum%first(spectrum%shells + 1), spectrum%vectors(3, total), &
                  position(0:cutoff))
        shell = 0
        total = 1
        do q = 0, cutoff
            position(q) = total
            if (count(q) == 0) cycle
            shell = shell + 1
            spectrum%level(shell) = real(q, dp)
            spectrum%degeneracy(shell) = count(q)
            spectrum%first(shell) = total
            total = total + count(q)
        end do
        spectrum%first(spectrum%shells + 1) = total
        do i = -m, m
            do j = -m, m
                do k = -m, m
                    q = i*i + j*j + k*k
                    if (q > cutoff) cycle
                    spectrum%vectors(:, position(q)) = [i, j, k]
                    position(q) = position(q) + 1
                end do
            end do
        end do
        ok = .true.
    end subroutine box_spectrum

    !> Isotropic 3D harmonic trap levels k=0..cutoff (units hbar*omega above
    !> the zero point), degeneracy (k+1)*(k+2)/2.
    subroutine trap_spectrum(cutoff, spectrum, ok, closed)
        integer, intent(in) :: cutoff
        type(fermion_spectrum_t), intent(out) :: spectrum
        logical, intent(out) :: ok
        logical, intent(in), optional :: closed
        integer :: k
        ok = .false.
        if (cutoff < 0 .or. cutoff > 10**4) return
        spectrum%kind = trap_kind
        spectrum%cutoff = cutoff
        if (present(closed)) spectrum%closed = closed
        spectrum%shells = cutoff + 1
        allocate (spectrum%level(cutoff + 1), spectrum%degeneracy(cutoff + 1), &
                  spectrum%first(cutoff + 2))
        spectrum%first(1) = 1
        do k = 0, cutoff
            spectrum%level(k + 1) = real(k, dp)
            spectrum%degeneracy(k + 1) = (k + 1)*(k + 2)/2
            spectrum%first(k + 2) = spectrum%first(k + 1) + spectrum%degeneracy(k + 1)
        end do
        spectrum%orbitals = spectrum%first(cutoff + 2) - 1
        ok = .true.
    end subroutine trap_spectrum

    pure integer function count_nonzero(histogram) result(n)
        integer, intent(in) :: histogram(0:)
        n = count(histogram > 0)
    end function count_nonzero

    !> Box cutoff with omitted weights below exp(-46) of the Fermi level and
    !> the remainder monotonicity margin b*a**2>=2 (a=sqrt(cutoff+1)-sqrt(3)).
    pure integer function box_cutoff(b, particles) result(cutoff)
        real(dp), intent(in) :: b
        integer, intent(in) :: particles
        real(dp) :: fermi_level, a
        fermi_level = (3.0_dp*real(particles, dp)/(4.0_dp*acos(-1.0_dp)))**(2.0_dp/3.0_dp)
        cutoff = ceiling(1.5_dp*fermi_level + 4.0_dp + 46.0_dp/b)
        do
            a = sqrt(real(cutoff + 1, dp)) - sqrt(3.0_dp)
            if (a > 0.0_dp .and. b*a*a > 2.5_dp) exit
            cutoff = cutoff + 1 + cutoff/4
        end do
    end function box_cutoff

    !> Thermodynamic-limit EF (spectrum units) of N electrons in the box and
    !> the inverse temperature b=1/(theta*EF).
    subroutine box_parameters(electrons, theta, ef, b)
        integer, intent(in) :: electrons
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: ef, b
        call fermion_box_fermi_energy(real(electrons, dp), ef)
        call fermion_inverse_temperature(theta, ef, b)
    end subroutine box_parameters

    !> Canonical fixed-N state of one spin species at inverse temperature b.
    !> certify (default true) adds outward enclosures of lnZ, E, Var(E) and
    !> the occupations, including the omitted-orbital remainder.
    subroutine solve_canonical(spectrum, particles, b, state, ok, certify)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: particles
        real(dp), intent(in) :: b
        type(canonical_state_t), intent(out) :: state
        logical, intent(out) :: ok
        logical, intent(in), optional :: certify
        real(dp), allocatable :: y(:, :), z(:, :), pre(:, :), suf(:, :), t(:)
        real(dp) :: d1, d2
        integer :: n, s, a, rep
        ok = .false.
        state%status = 'invalid-input'
        s = spectrum%shells
        if (s < 1 .or. .not. allocated(spectrum%level)) return
        if (particles < 1 .or. particles > spectrum%orbitals) return
        if (.not. ieee_is_finite(b) .or. b <= 0.0_dp) return
        n = particles
        state%particles = n
        state%b = b
        state%shift = scale_shift(spectrum, n, b)
        allocate (y(0:2, s), z(0:n, 0:2), pre(0:n, 0:s + 1), suf(0:n, 0:s + 1), &
                  t(0:n), state%occupation(s), state%weight(s), &
                  state%same_shell_pair(s))
        do a = 1, s
            call fermion_weight_jet(b, spectrum%level(a), state%shift, y(0, a), &
                                    y(1, a), y(2, a))
        end do
        state%weight = y(0, :)
        z = 0.0_dp
        z(0, 0) = 1.0_dp
        do a = 1, s
            do rep = 1, spectrum%degeneracy(a)
                call jet_multiply(z, y(:, a))
            end do
        end do
        state%status = 'candidate-range-failed'
        if (.not. all(ieee_is_finite(z)) .or. z(n, 0) <= 0.0_dp) return
        d1 = z(n, 1)/z(n, 0)
        d2 = z(n, 2)/z(n, 0)
        state%log_z = log(z(n, 0)) - b*real(n, dp)*state%shift
        state%energy = -d1 + real(n, dp)*state%shift
        state%variance = d2 - d1*d1
        call prefix_suffix(spectrum, y(0, :), pre, suf)
        do a = 1, s
            call reduced_shell(pre(:, a - 1), suf(:, a + 1), y(0, a), &
                               spectrum%degeneracy(a), n, t, state%same_shell_pair(a))
            state%occupation(a) = min(1.0_dp, y(0, a)*t(n - 1)/z(n, 0))
            state%same_shell_pair(a) = state%same_shell_pair(a)/z(n, 0)
        end do
        if (.not. all(ieee_is_finite(state%occupation))) return
        state%valid = .true.
        state%status = 'candidate'
        ok = .true.
        if (present(certify)) then
            if (.not. certify) return
        end if
        call certify_state(spectrum, state, ok)
    end subroutine solve_canonical

    !> Zero-temperature limit of the canonical ensemble: closed shells filled,
    !> m particles of an open shell of g orbitals spread uniformly over its
    !> C(g,m) determinants (n=m/g, same-shell pair m(m-1)/(g(g-1))). Exact
    !> rational occupations up to one rounding; b is reported as zero.
    subroutine solve_ground(spectrum, particles, state, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: particles
        type(canonical_state_t), intent(out) :: state
        logical, intent(out) :: ok
        integer :: a, left, take, g
        ok = .false.
        state%status = 'invalid-input'
        if (particles < 1 .or. particles > spectrum%orbitals) return
        allocate (state%occupation(spectrum%shells), state%weight(spectrum%shells), &
                  state%same_shell_pair(spectrum%shells))
        state%occupation = 0.0_dp; state%weight = 0.0_dp; state%same_shell_pair = 0.0_dp
        state%particles = particles
        left = particles
        do a = 1, spectrum%shells
            g = spectrum%degeneracy(a)
            take = min(left, g)
            if (take == 0) exit
            state%occupation(a) = real(take, dp)/real(g, dp)
            if (g >= 2) state%same_shell_pair(a) = real(take, dp)*real(take - 1, dp)/ &
                (real(g, dp)*real(g - 1, dp))
            state%energy = state%energy + real(take, dp)*spectrum%level(a)
            left = left - take
        end do
        state%omitted_bound = interval(0.0_dp)
        state%ground = .true.
        state%valid = .true.
        state%status = 'ground-state'
        ok = .true.
    end subroutine solve_ground

    !> Mean single-particle energy of the N-particle ground configuration;
    !> the exact binary64 shift keeps the scaled sums inside the range.
    pure real(dp) function ground_mean(spectrum, n) result(c)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: n
        integer :: a, left, take
        c = 0.0_dp
        left = n
        do a = 1, spectrum%shells
            take = min(left, spectrum%degeneracy(a))
            c = c + real(take, dp)*spectrum%level(a)
            left = left - take
            if (left == 0) exit
        end do
        c = c/real(n, dp)
    end function ground_mean

    !> Exact binary64 energy shift c of the weights exp(-b*(eps-c)): an
    !> estimate of the free energy per particle, F/N=-log(Z_N)/(b*N), from the
    !> grand-canonical (Chernoff) bound at the mu that holds N particles on the
    !> listed orbitals. It keeps Z_N near one; the remaining range of the
    !> lower-degree sums is about exp(N/e) hot and exp(0.11*N/theta) cold.
    !> Only scaling depends on it; every reported quantity is shift-free.
    real(dp) function scale_shift(spectrum, n, b) result(c)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: n
        real(dp), intent(in) :: b
        real(dp) :: lo, hi, mu, count, log_z
        integer :: iteration
        c = ground_mean(spectrum, n)
        lo = spectrum%level(1) - 1.0_dp/b
        do iteration = 1, 200
            call grand_sums(spectrum, b, lo, count, log_z)
            if (count < real(n, dp)) exit
            lo = lo - 2.0_dp*(abs(lo) + 1.0_dp/b)
        end do
        hi = spectrum%level(spectrum%shells) + 1.0_dp/b
        call grand_sums(spectrum, b, hi, count, log_z)
        if (count < real(n, dp)) return
        do iteration = 1, 100
            mu = 0.5_dp*(lo + hi)
            call grand_sums(spectrum, b, mu, count, log_z)
            if (count < real(n, dp)) then
                lo = mu
            else
                hi = mu
            end if
        end do
        call grand_sums(spectrum, b, mu, count, log_z)
        log_z = log_z - b*mu*real(n, dp)
        if (ieee_is_finite(log_z)) c = min(c, -log_z/(b*real(n, dp)))
    end function scale_shift

    pure subroutine grand_sums(spectrum, b, mu, count, log_z)
        type(fermion_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: b, mu
        real(dp), intent(out) :: count, log_z
        real(dp) :: x, occupation, factor
        integer :: a
        count = 0.0_dp; log_z = 0.0_dp
        do a = 1, spectrum%shells
            x = b*(mu - spectrum%level(a))
            if (x >= 0.0_dp) then
                call fermion_grand_positive(x, occupation, factor)
            else
                call fermion_grand_negative(x, occupation, factor)
            end if
            count = count + real(spectrum%degeneracy(a), dp)*occupation
            log_z = log_z + real(spectrum%degeneracy(a), dp)*factor
        end do
    end subroutine grand_sums

    !> Multiply the truncated generating polynomial (with b-jets) by (1+w*x).
    pure subroutine jet_multiply(z, w)
        real(dp), intent(inout) :: z(0:, 0:)
        real(dp), intent(in) :: w(0:2)
        integer :: k
        if (.not. (w(0) > 0.0_dp)) return
        do k = ubound(z, 1), 1, -1
            z(k, 2) = z(k, 2) + w(2)*z(k - 1, 0) + 2.0_dp*w(1)*z(k - 1, 1) + &
                      w(0)*z(k - 1, 2)
            z(k, 1) = z(k, 1) + w(1)*z(k - 1, 0) + w(0)*z(k - 1, 1)
            z(k, 0) = z(k, 0) + w(0)*z(k - 1, 0)
        end do
    end subroutine jet_multiply

    pure subroutine poly_multiply(p, w, times)
        real(dp), intent(inout) :: p(0:)
        real(dp), intent(in) :: w
        integer, intent(in) :: times
        integer :: k, rep
        if (w <= 0.0_dp) return
        do rep = 1, times
            do k = ubound(p, 1), 1, -1
                p(k) = p(k) + w*p(k - 1)
            end do
        end do
    end subroutine poly_multiply

    pure subroutine prefix_suffix(spectrum, y, pre, suf)
        type(fermion_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: pre(0:, 0:), suf(0:, 0:)
        integer :: a, s
        s = spectrum%shells
        pre = 0.0_dp; suf = 0.0_dp
        pre(0, 0) = 1.0_dp; suf(0, s + 1) = 1.0_dp
        do a = 1, s
            pre(:, a) = pre(:, a - 1)
            call poly_multiply(pre(:, a), y(a), spectrum%degeneracy(a))
            suf(:, s + 1 - a) = suf(:, s + 2 - a)
            call poly_multiply(suf(:, s + 1 - a), y(s + 1 - a), &
                               spectrum%degeneracy(s + 1 - a))
        end do
    end subroutine prefix_suffix

    !> Generating polynomial with one orbital of the shell removed (t), and the
    !> unnormalized same-shell pair weight y**2*[x**(n-2)] with two removed.
    pure subroutine reduced_shell(left, right, w, g, n, t, pair)
        real(dp), intent(in) :: left(0:), right(0:), w
        integer, intent(in) :: g, n
        real(dp), intent(out) :: t(0:), pair
        integer :: k
        do k = 0, n
            t(k) = sum(left(0:k)*right(k:0:-1))
        end do
        pair = 0.0_dp
        if (g >= 2) then
            call poly_multiply(t, w, g - 2)
            if (n >= 2) pair = w*w*t(n - 2)
            call poly_multiply(t, w, 1)
        end if
    end subroutine reduced_shell

    !> Outward enclosure of the same subtraction-free recursion plus the
    !> omitted-orbital remainder: with tau_m >= sum_omitted |eps-c|**m*y, the
    !> omitted elementary sums obey e_j<=tau0**j/j!, |e_j'|<=tau1*tau0**(j-1)/(j-1)!,
    !> |e_j''|<=tau2*tau0**(j-1)/(j-1)!+tau1**2*tau0**(j-2)/(j-2)!.
    subroutine certify_state(spectrum, state, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(canonical_state_t), intent(inout) :: state
        logical, intent(out) :: ok
        type(interval_t), allocatable :: y(:, :), z(:, :), pre(:, :), suf(:, :), t(:)
        type(interval_t) :: tau(0:2), full(0:2), below, rem(0:2), d1, nc, ib, gap
        integer :: n, s, a, rep, last
        ok = .false.
        state%status = 'certificate-failed'
        n = state%particles
        s = spectrum%shells
        allocate (y(0:2, s), z(0:n, 0:2), pre(0:n, 0:s + 1), suf(0:n, 0:s + 1), &
                  t(0:n), state%occupation_bound(s))
        ib = interval(state%b)
        do a = 1, s
            call fermion_weight_jet_interval(ib, interval(spectrum%level(a)), &
                interval(state%shift), y(0, a), y(1, a), y(2, a))
        end do
        call omitted_sums(spectrum, state%b, state%shift, tau, ok)
        if (.not. ok) return
        ok = .false.
        ! Listed shells whose candidate occupation is below 1e-18 join the
        ! omitted set: their exact weight sums are added to tau instead of
        ! being multiplied through. Any choice keeps the enclosure valid; this
        ! one only trades work for a negligible remainder width.
        last = s
        do while (last > 1)
            if (state%occupation(last) >= 1.0e-18_dp .or. &
                spectrum%first(last) <= n) exit
            last = last - 1
        end do
        do a = last + 1, s
            gap = interval(spectrum%level(a)) - interval(state%shift)
            gap = interval(0.0_dp, mag(gap))
            rem(0) = interval(real(spectrum%degeneracy(a), dp))*y(0, a)
            tau(0) = tau(0) + rem(0)
            tau(1) = tau(1) + rem(0)*gap
            tau(2) = tau(2) + rem(0)*gap*gap
        end do
        z = interval(0.0_dp)
        z(0, 0) = interval(1.0_dp)
        do a = 1, last
            do rep = 1, spectrum%degeneracy(a)
                call ijet_multiply(z, y(:, a))
            end do
        end do
        call remainder(z, n, tau, rem)
        full = [z(n, 0) + rem(0), z(n, 1) + rem(1), z(n, 2) + rem(2)]
        if (.not. all_finite(full) .or. full(0)%lo <= 0.0_dp) return
        call remainder(z, n - 1, tau, rem)
        below = z(n - 1, 0) + rem(0)
        nc = interval(real(n, dp))*interval(state%shift)
        d1 = full(1)/full(0)
        state%log_z_bound = hull(log(full(0)) - ib*nc, interval(state%log_z))
        state%energy_bound = hull(nc - d1, interval(state%energy))
        state%variance_bound = hull(full(2)/full(0) - d1*d1, interval(state%variance))
        state%omitted_bound = interval(0.0_dp, mag(tau(0)*below/full(0)))
        call iprefix_suffix(spectrum, y(0, :), last, pre, suf)
        do a = last + 1, s
            state%occupation_bound(a) = hull(interval(0.0_dp, &
                mag(y(0, a)*below/full(0))), interval(state%occupation(a)))
        end do
        do a = 1, last
            call ireduced_shell(pre(:, a - 1), suf(:, a + 1), y(0, a), &
                                spectrum%degeneracy(a), n, t)
            rem(0) = interval(0.0_dp, mag(tail_sum(t, n - 1, tau(0))))
            state%occupation_bound(a) = y(0, a)*(t(n - 1) + rem(0))/full(0)
            state%occupation_bound(a) = interval(max(0.0_dp, &
                state%occupation_bound(a)%lo), min(1.0_dp, state%occupation_bound(a)%hi))
            state%occupation_bound(a) = hull(state%occupation_bound(a), &
                                             interval(state%occupation(a)))
        end do
        if (.not. all_finite([state%log_z_bound, state%energy_bound, &
            state%variance_bound, state%omitted_bound])) return
        if (.not. all_finite(state%occupation_bound)) return
        state%status = 'certificate-budget-failed'
        if (.not. within_budget(state%log_z_bound, state%log_z)) return
        if (.not. within_budget(state%energy_bound, state%energy)) return
        if (.not. within_budget(state%variance_bound, state%variance)) return
        state%certified = .true.
        state%status = 'certified-model'
        ok = .true.
    end subroutine certify_state

    !> 1% relative output budget with a 1e-12 absolute floor (spectrum units).
    pure logical function within_budget(bound, center) result(ok)
        type(interval_t), intent(in) :: bound
        real(dp), intent(in) :: center
        ok = max(center - bound%lo, bound%hi - center) <= &
             max(0.01_dp*abs(center), 1.0e-12_dp)
    end function within_budget

    pure logical function all_finite(values) result(ok)
        type(interval_t), intent(in) :: values(:)
        integer :: j
        ok = .false.
        do j = 1, size(values)
            if (is_empty(values(j))) return
            if (.not. ieee_is_finite(values(j)%lo)) return
            if (.not. ieee_is_finite(values(j)%hi)) return
        end do
        ok = .true.
    end function all_finite

    !> Upper bounds tau_m on sum over omitted orbitals of |eps-c|**m*y.
    subroutine omitted_sums(spectrum, b, shift, tau, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: b, shift
        type(interval_t), intent(out) :: tau(0:2)
        logical, intent(out) :: ok
        type(interval_t) :: a, margin, ratio, deficit
        integer :: j
        ok = .false.
        tau = interval(0.0_dp)
        if (spectrum%closed) then
            ok = .true.
            return
        end if
        if (real(spectrum%cutoff + 1, dp) < shift) return
        deficit = interval(max(0.0_dp, -shift))
        if (spectrum%kind == box_kind) then
            call fermion_box_tail_margin_interval(interval(b), &
                interval(real(spectrum%cutoff, dp)), a, margin)
            if (.not. all_finite([a, margin])) return
            if (a%lo <= 0.0_dp .or. margin%lo < 0.0_dp) return
            call fermion_box_tail(interval(b), interval(shift), &
                interval(real(spectrum%cutoff, dp)), deficit, tau(0), tau(1), tau(2))
        else if (spectrum%kind == trap_kind) then
            call fermion_trap_tail(interval(b), interval(shift), &
                interval(real(spectrum%cutoff, dp)), deficit, tau(0), tau(1), tau(2), &
                ratio)
            if (.not. all_finite([ratio])) return
            if (ratio%hi >= 1.0_dp) return
        else
            return
        end if
        if (.not. all_finite(tau)) return
        do j = 0, 2
            tau(j) = interval(0.0_dp, tau(j)%hi)
        end do
        ok = .true.
    end subroutine omitted_sums

    !> Remainder jets of coefficient k from the omitted orbitals.
    pure subroutine remainder(z, k, tau, rem)
        type(interval_t), intent(in) :: z(0:, 0:), tau(0:2)
        integer, intent(in) :: k
        type(interval_t), intent(out) :: rem(0:2)
        type(interval_t) :: e0, e1, e2, previous, r0, r1, r2, t0, t1, t2
        integer :: j
        r0 = interval(0.0_dp); r1 = r0; r2 = r0
        t0 = interval(tau(0)%hi); t1 = interval(tau(1)%hi); t2 = interval(tau(2)%hi)
        e0 = interval(1.0_dp)
        previous = interval(0.0_dp)
        do j = 1, k
            ! previous=tau0**(j-2)/(j-2)!, e0 on entry tau0**(j-1)/(j-1)!.
            e1 = t1*e0
            e2 = t2*e0 + t1*t1*previous
            previous = e0
            e0 = e0*t0/interval(real(j, dp))
            r0 = r0 + interval(mag(z(k - j, 0)))*e0
            r1 = r1 + interval(mag(z(k - j, 1)))*e0 + interval(mag(z(k - j, 0)))*e1
            r2 = r2 + interval(mag(z(k - j, 2)))*e0 + &
                 interval(2.0_dp)*interval(mag(z(k - j, 1)))*e1 + &
                 interval(mag(z(k - j, 0)))*e2
        end do
        rem(0) = interval(0.0_dp, r0%hi)
        rem(1) = interval(-r1%hi, r1%hi)
        rem(2) = interval(-r2%hi, r2%hi)
    end subroutine remainder

    pure function tail_sum(t, k, tau0) result(r0)
        type(interval_t), intent(in) :: t(0:), tau0
        integer, intent(in) :: k
        type(interval_t) :: r0, e0
        integer :: j
        r0 = interval(0.0_dp)
        e0 = interval(1.0_dp)
        do j = 1, k
            e0 = e0*interval(tau0%hi)/interval(real(j, dp))
            r0 = r0 + interval(mag(t(k - j)))*e0
        end do
    end function tail_sum

    pure subroutine ijet_multiply(z, w)
        type(interval_t), intent(inout) :: z(0:, 0:)
        type(interval_t), intent(in) :: w(0:2)
        integer :: k
        if (.not. (w(0)%hi > 0.0_dp)) return
        do k = ubound(z, 1), 1, -1
            z(k, 2) = z(k, 2) + w(2)*z(k - 1, 0) + interval(2.0_dp)*w(1)*z(k - 1, 1) &
                      + w(0)*z(k - 1, 2)
            z(k, 1) = z(k, 1) + w(1)*z(k - 1, 0) + w(0)*z(k - 1, 1)
            z(k, 0) = z(k, 0) + w(0)*z(k - 1, 0)
        end do
    end subroutine ijet_multiply

    pure subroutine ipoly_multiply(p, w, times)
        type(interval_t), intent(inout) :: p(0:)
        type(interval_t), intent(in) :: w
        integer, intent(in) :: times
        integer :: k, rep
        if (w%hi <= 0.0_dp) return
        do rep = 1, times
            do k = ubound(p, 1), 1, -1
                p(k) = p(k) + w*p(k - 1)
            end do
        end do
    end subroutine ipoly_multiply

    pure subroutine iprefix_suffix(spectrum, y, s, pre, suf)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(interval_t), intent(in) :: y(:)
        integer, intent(in) :: s
        type(interval_t), intent(out) :: pre(0:, 0:), suf(0:, 0:)
        integer :: a
        pre = interval(0.0_dp); suf = interval(0.0_dp)
        pre(0, 0) = interval(1.0_dp); suf(0, s + 1) = interval(1.0_dp)
        do a = 1, s
            pre(:, a) = pre(:, a - 1)
            call ipoly_multiply(pre(:, a), y(a), spectrum%degeneracy(a))
            suf(:, s + 1 - a) = suf(:, s + 2 - a)
            call ipoly_multiply(suf(:, s + 1 - a), y(s + 1 - a), &
                                spectrum%degeneracy(s + 1 - a))
        end do
    end subroutine iprefix_suffix

    pure subroutine ireduced_shell(left, right, w, g, n, t)
        type(interval_t), intent(in) :: left(0:), right(0:), w
        integer, intent(in) :: g, n
        type(interval_t), intent(out) :: t(0:)
        integer :: k, j
        do k = 0, n
            t(k) = interval(0.0_dp)
            do j = 0, k
                t(k) = t(k) + left(j)*right(k - j)
            end do
        end do
        call ipoly_multiply(t, w, g - 1)
    end subroutine ireduced_shell

    !> Classical Maxwell-Boltzmann preparation of N particles in the same
    !> spectrum: Z=Z1**N/N!, no exclusion; occupations may exceed one.
    subroutine solve_boltzmann(spectrum, particles, b, state, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: particles
        real(dp), intent(in) :: b
        type(canonical_state_t), intent(out) :: state
        logical, intent(out) :: ok
        type(interval_t) :: factorial, count
        integer :: k
        ok = .false.
        if (particles < 1) return
        call solve_canonical(spectrum, 1, b, state, ok)
        if (.not. ok) return
        count = interval(real(particles, dp))
        factorial = interval(0.0_dp)
        do k = 2, particles
            factorial = factorial + log(interval(real(k, dp)))
        end do
        state%particles = particles
        state%boltzmann = .true.
        state%log_z_bound = count*state%log_z_bound - factorial
        state%log_z = real(particles, dp)*state%log_z - mid_value(factorial)
        state%log_z_bound = hull(state%log_z_bound, interval(state%log_z))
        state%energy = real(particles, dp)*state%energy
        state%energy_bound = hull(count*state%energy_bound, interval(state%energy))
        state%variance = real(particles, dp)*state%variance
        state%variance_bound = hull(count*state%variance_bound, interval(state%variance))
        state%occupation = real(particles, dp)*state%occupation
        do k = 1, size(state%occupation)
            state%occupation_bound(k) = hull(count*state%occupation_bound(k), &
                                             interval(state%occupation(k)))
        end do
        state%omitted_bound = count*state%omitted_bound
        state%same_shell_pair = 0.0_dp
        state%status = 'certified-model-boltzmann'
    end subroutine solve_boltzmann

    pure real(dp) function mid_value(x) result(m)
        type(interval_t), intent(in) :: x
        m = 0.5_dp*(x%lo + x%hi)
    end function mid_value

    !> Box EOS for two spin species at a common b, in units of the
    !> thermodynamic-limit EF of the same total density:
    !> values=[E/(N*EF), P/(n*EF), de/dtheta, dp/dtheta, B/(n*EF)].
    subroutine box_eos(up, down, values, bounds, ok, certified)
        type(canonical_state_t), intent(in) :: up, down
        real(dp), intent(out) :: values(5)
        type(interval_t), intent(out) :: bounds(5)
        logical, intent(out) :: ok, certified
        real(dp) :: ef, total
        type(interval_t) :: ief, itotal
        integer :: j
        ok = .false.; certified = .false.
        values = 0.0_dp; bounds = interval(0.0_dp)
        if (.not. (up%valid .and. down%valid)) return
        if (up%ground .or. down%ground .or. (up%boltzmann .neqv. down%boltzmann)) return
        if (up%b < down%b .or. up%b > down%b) return
        total = real(up%particles + down%particles, dp)
        call fermion_box_fermi_energy(total, ef)
        call fermion_box_eos(up%energy + down%energy, up%variance + down%variance, &
            up%b, total, ef, values(1), values(2), values(3), values(4), values(5))
        if (.not. all(ieee_is_finite(values))) return
        ok = .true.
        if (.not. (up%certified .and. down%certified)) return
        itotal = interval(total)
        call fermion_box_fermi_energy_interval(itotal, ief)
        call fermion_box_eos_interval(up%energy_bound + down%energy_bound, &
            up%variance_bound + down%variance_bound, interval(up%b), itotal, ief, &
            bounds(1), bounds(2), bounds(3), bounds(4), bounds(5))
        if (.not. all_finite(bounds)) return
        do j = 1, 5
            bounds(j) = hull(bounds(j), interval(values(j)))
            ! Declared normalized scale: 1% of max(1,|value|).
            if (max(values(j) - bounds(j)%lo, bounds(j)%hi - values(j)) > &
                0.01_dp*max(1.0_dp, abs(values(j)))) return
        end do
        certified = .true.
    end subroutine box_eos

    pure real(dp) function sphere_average(q, r) result(value)
        real(dp), intent(in) :: q, r
        real(dp) :: z
        call fermion_lattice_argument(q, r, z)
        if (abs(z) < 0.05_dp) then
            call fermion_sinc_series(z, value)
        else
            call fermion_sinc(z, value)
        end if
    end function sphere_average

    !> Sphere-averaged same-spin pair correlation g(r), r in units of L:
    !> g=[N(N-1)-sum_{i/=j}<n_i n_j>*sin(qr)/(qr)]/N**2, so g(0)=0 and
    !> g->1-1/N between exchange oscillations. Shells with occupation at or
    !> below tolerance are omitted; truncation is the resulting bound
    !> 2(N-1)*omitted occupation/N**2 on |delta g| in exact arithmetic.
    subroutine same_spin_pair(spectrum, state, radii, g, ok, tolerance, truncation)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(canonical_state_t), intent(in) :: state
        real(dp), intent(in) :: radii(:), tolerance
        real(dp), intent(out) :: g(:), truncation
        logical, intent(out) :: ok
        real(dp), allocatable :: c(:, :), h(:)
        integer, allocatable :: shell_of(:)
        real(dp) :: n, excluded, total
        integer :: last, a, bb, i, j, q, qmax, orbitals
        ok = .false.
        g = 0.0_dp; truncation = 0.0_dp
        if (spectrum%kind /= box_kind .or. .not. allocated(spectrum%vectors)) return
        if (.not. state%valid .or. state%boltzmann) return
        if (size(g) /= size(radii)) return
        n = real(state%particles, dp)
        last = 0
        do a = 1, spectrum%shells
            if (state%occupation(a) > tolerance) last = a
        end do
        if (last == 0) return
        excluded = state%omitted_bound%hi
        do a = last + 1, spectrum%shells
            excluded = excluded + real(spectrum%degeneracy(a), dp)*state%occupation(a)
        end do
        truncation = 2.0_dp*(n - 1.0_dp)*excluded/n**2
        if (state%particles < 2) then
            ok = .true.
            return
        end if
        allocate (c(last, last))
        do a = 1, last
            do bb = 1, last
                if (a == bb) then
                    c(a, a) = state%same_shell_pair(a)
                else if (state%ground) then
                    c(a, bb) = state%occupation(a)*state%occupation(bb)
                else
                    call fermion_pair_occupation(state%weight(a), state%weight(bb), &
                        state%occupation(a), state%occupation(bb), c(a, bb))
                end if
            end do
        end do
        orbitals = spectrum%first(last + 1) - 1
        allocate (shell_of(orbitals))
        do a = 1, last
            shell_of(spectrum%first(a):spectrum%first(a + 1) - 1) = a
        end do
        qmax = 4*int(spectrum%level(last)) + 1
        allocate (h(0:qmax), source=0.0_dp)
        do i = 1, orbitals - 1
            do j = i + 1, orbitals
                q = sum((spectrum%vectors(:, i) - spectrum%vectors(:, j))**2)
                h(q) = h(q) + c(shell_of(i), shell_of(j))
            end do
        end do
        do j = 1, size(radii)
            total = 0.0_dp
            do q = 1, qmax
                if (h(q) > 0.0_dp .or. h(q) < 0.0_dp) total = total + &
                    h(q)*sphere_average(real(q, dp), radii(j))
            end do
            g(j) = (n*(n - 1.0_dp) - 2.0_dp*total)/n**2
        end do
        ok = all(ieee_is_finite(g))
    end subroutine same_spin_pair

    !> Sphere-averaged one-body density matrix rho(r)/n_sigma, r in units of L.
    subroutine one_body_density(spectrum, state, radii, rho, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(canonical_state_t), intent(in) :: state
        real(dp), intent(in) :: radii(:)
        real(dp), intent(out) :: rho(:)
        logical, intent(out) :: ok
        integer :: a, j
        ok = .false.
        rho = 0.0_dp
        if (spectrum%kind /= box_kind .or. .not. state%valid) return
        if (size(rho) /= size(radii)) return
        do j = 1, size(radii)
            do a = 1, spectrum%shells
                rho(j) = rho(j) + real(spectrum%degeneracy(a), dp)*state%occupation(a)* &
                         sphere_average(spectrum%level(a), radii(j))
            end do
            rho(j) = rho(j)/real(state%particles, dp)
        end do
        ok = all(ieee_is_finite(rho))
    end subroutine one_body_density

    !> Opposite-spin pair correlation of independent species with uniform
    !> plane-wave densities: g=(sum n_up)*(sum n_down)/(N_up*N_down)=1.
    pure real(dp) function opposite_spin_pair(spectrum, up, down) result(g)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(canonical_state_t), intent(in) :: up, down
        g = sum(real(spectrum%degeneracy, dp)*up%occupation)* &
            sum(real(spectrum%degeneracy, dp)*down%occupation)/ &
            (real(up%particles, dp)*real(down%particles, dp))
    end function opposite_spin_pair

    !> Thermodynamic-limit (grand-canonical) same-density reference from the
    !> certified homogeneous Fermi gas: rho(x)/n=3*int u**2 f(u) sinc(u*x) du,
    !> x=kF*r, and the Wick hole g=1-rho**2. Gauss-Legendre candidate on
    !> panels split at the Fermi edge; a reference curve, not an enclosure.
    subroutine grand_canonical_density(fermi, x, rho, ok)
        type(fermi_state_t), intent(in) :: fermi
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: rho(:)
        logical, intent(out) :: ok
        integer, parameter :: nodes = 256
        real(dp) :: edges(4), u(nodes), w(nodes), f, lower, upper, term, s, z
        real(dp) :: edge, width
        integer :: panel, j, k, panels
        ok = .false.
        rho = 0.0_dp
        if (size(rho) /= size(x)) return
        edge = sqrt(max(fermi%mu, 0.0_dp))
        edges(4) = sqrt(max(fermi%mu, 0.0_dp) + 46.0_dp*fermi%theta)
        width = 30.0_dp*fermi%theta/max(edge, 1.0e-300_dp)
        if (edge > 0.0_dp .and. width < 0.5_dp*edge) then
            panels = 3
            edges(1:3) = [0.0_dp, edge - width, edge + width]
            edges(4) = max(edges(4), edge + 2.0_dp*width)
        else
            panels = 1
            edges(1:2) = [0.0_dp, edges(4)]
        end if
        do panel = 1, panels
            if (panels == 1) then
                call gauss_legendre_ab(nodes, edges(1), edges(2), u, w)
            else
                call gauss_legendre_ab(nodes, edges(panel), edges(panel + 1), u, w)
            end if
            do k = 1, nodes
                call occupation_fermi_gas(fermi, u(k), f, lower, upper, ok)
                if (.not. ok) return
                do j = 1, size(x)
                    z = u(k)*x(j)
                    if (abs(z) < 0.05_dp) then
                        call fermion_sinc_series(z, s)
                    else
                        call fermion_sinc(z, s)
                    end if
                    call fermion_density_matrix_term(u(k), f, s, term)
                    rho(j) = rho(j) + w(k)*term
                end do
            end do
        end do
        ok = all(ieee_is_finite(rho))
    end subroutine grand_canonical_density

    !> Cross-check only: the alternating Newton power-sum recursion
    !> Z_N=(1/N)sum_k (-1)**(k+1) Z_1(k*b) Z_(N-k) in binary128 on the listed
    !> orbitals (no omitted remainder), with the same exact shift. cancellation
    !> is max_k sum|terms|/|k*Z_k|; its product with 2**-112 estimates the
    !> relative rounding loss. Production sums never use this recursion.
    subroutine power_sum_recursion(spectrum, particles, b, log_z, cancellation, ok)
        type(fermion_spectrum_t), intent(in) :: spectrum
        integer, intent(in) :: particles
        real(dp), intent(in) :: b
        real(qp), intent(out) :: log_z, cancellation
        logical, intent(out) :: ok
        real(qp), allocatable :: power(:), e(:), y(:)
        real(qp) :: term, magnitude, c
        integer :: k, i
        ok = .false.
        log_z = 0.0_qp; cancellation = 0.0_qp
        if (particles < 1 .or. particles > spectrum%orbitals) return
        c = real(ground_mean(spectrum, particles), qp)
        allocate (power(particles), e(0:particles), y(spectrum%shells))
        y = exp(-real(b, qp)*(real(spectrum%level, qp) - c))
        do k = 1, particles
            power(k) = sum(real(spectrum%degeneracy, qp)*y**k)
        end do
        e(0) = 1.0_qp
        do k = 1, particles
            term = 0.0_qp
            magnitude = 0.0_qp
            do i = 1, k
                term = term + real((-1)**(i - 1), qp)*e(k - i)*power(i)
                magnitude = magnitude + abs(e(k - i)*power(i))
            end do
            e(k) = term/real(k, qp)
            if (abs(term) > 0.0_qp) cancellation = max(cancellation, magnitude/abs(term))
        end do
        if (.not. (e(particles) > 0.0_qp)) return
        log_z = log(e(particles)) - real(b, qp)*real(particles, qp)*c
        ok = .true.
    end subroutine power_sum_recursion

end module kin6d_fermion_ensemble
