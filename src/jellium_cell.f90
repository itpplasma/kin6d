!> Interacting finite-N periodic jellium cell (Q4): a consumer of the general
!> Fock engine. Spin orbitals are plane waves V^(-1/2) exp(i 2 pi n.x/L) with
!> spin; orbital 2a-1 is mode a up, 2a is mode a down. The Hamiltonian is the
!> engine's jellium Galerkin projection (bare Coulomb, q=0 removed by the
!> neutralizing background). Total momentum K and S_z are exact symmetries:
!> every solve uses a sector basis and the matrix-free coupled action, which
!> fails closed on leakage. Atomic units (hbar=m=g=1) come from kin6d_species.
!> Claims: finite-basis eigen/evolution numbers only; momentum-cutoff,
!> finite-size and thermodynamic-limit errors are indicators, not bounds.
module kin6d_jellium_cell
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, ipoint, operator(/), operator(-)
    use fortnum_krylov, only: real_gmres_operator, KRYLOV_OK
    use fortnum_symmetric_eigen, only: block_lanczos_lowest
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_ode_residual_certificate, only: certify_linear_trajectory
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_species, only: units_t, species_t, plasma_scales_t, atomic_units, &
        electron, plasma_scales, pair_coupling
    use kin6d_quantum_fock, only: fock_basis_t, fock_hamiltonian_t, &
        init_fock_sector_basis, prepare_jellium_hamiltonian, &
        prepare_fock_hamiltonian, get_fock_integrals, fock_matvec, &
        fock_matvec_between, fock_matvec_between_interval, fock_pair_indices, &
        fock_rdms, fock_transition_rdms, fock_rdm3_entries, fock_basis_info
    use kin6d_quantum_fock_leaves, only: qf_kinetic, qf_coulomb
    use kin6d_quantum_correlations, only: quantum_cumulant_norms
    use kin6d_jellium_cell_leaves, only: jc_lindhard_term, jc_rpa, jc_boltzmann, &
        jc_fermi, jc_shell_moment, jc_shell_moment_series, jc_shell_average, jc_pair_correlation, jc_susceptibility, &
        jc_half, jc_kick_response, jc_thermal_term, jc_truncation_interval
    use kin6d_sort_keys, only: sort_key_order
    implicit none
    private
    public :: jellium_cell_t, jellium_sector_t, jellium_spectrum_t
    public :: jellium_modes, init_jellium_cell, jellium_sector, jellium_sector_momenta
    public :: jellium_hamiltonian, jellium_interaction_model, jellium_density_model
    public :: jellium_lowest, jellium_build_spectrum, jellium_thermal_weights
    public :: jellium_thermal_rdms, jellium_pair_bins, jellium_bin_operator
    public :: jellium_cumulants, jellium_static_response, jellium_ideal_response
    public :: jellium_kick_evolution, jellium_coupling_derivative
    public :: jellium_coupling_vjp, jellium_orbital
    integer(i64), parameter :: model_budget = 4_i64*1024_i64**3

    type :: jellium_cell_t
        integer :: n = 0, m = 0, nmodes = 0
        integer, allocatable :: modes(:, :), momentum(:, :), spin(:)
        real(dp) :: length = 0, volume = 0, hbar = 1, mass = 1, coupling = 1
        real(dp) :: rs = 0, fermi = 0
    end type

    type :: jellium_sector_t
        integer :: nup = 0, d = 0
        integer :: k(3) = 0
        type(fock_basis_t) :: basis
    end type

    !> Low-lying spectrum over all (S_z, K) sectors; real vectors stored flat
    !> (state j occupies flat(offset(j)+1 : offset(j)+d of its sector)).
    type :: jellium_spectrum_t
        integer :: count = 0, nsectors = 0
        type(jellium_sector_t), allocatable :: sectors(:)
        integer, allocatable :: kept(:)
        real(dp), allocatable :: top(:)
        integer, allocatable :: sector_of(:)
        integer(i64), allocatable :: offset(:)
        real(dp), allocatable :: energy(:), residual(:), flat(:)
    end type
contains
    !> Integer mode vectors with |n|^2 <= cutoff ('spherical') or
    !> max|n_i| <= cutoff ('cubic'), ordered by |n|^2 then lexicographically.
    subroutine jellium_modes(kind, cutoff, modes)
        character(*), intent(in) :: kind
        integer, intent(in) :: cutoff
        integer, allocatable, intent(out) :: modes(:, :)
        integer :: buffer(3, 4096), count, a, b, c, r, s, i, j, tmp(3)
        logical :: take

        count = 0
        r = cutoff
        do a = -r, r
            do b = -r, r
                do c = -r, r
                    if (kind == 'cubic') then
                        take = .true.
                    else
                        take = a*a + b*b + c*c <= cutoff
                    end if
                    if (.not. take .or. count == size(buffer, 2)) cycle
                    count = count + 1; buffer(:, count) = [a, b, c]
                end do
            end do
        end do
        do i = 2, count
            tmp = buffer(:, i); s = sum(tmp**2); j = i - 1
            do while (j >= 1)
                if (sum(buffer(:, j)**2) <= s) exit
                buffer(:, j + 1) = buffer(:, j); j = j - 1
            end do
            buffer(:, j + 1) = tmp
        end do
        modes = buffer(:, 1:count)
    end subroutine

    !> Cell for N electrons in the given modes at density parameter rs.
    !> L follows from the shared rs definition (rs scales linearly with L);
    !> an optional coupling replaces g at fixed L (coupling scans).
    subroutine init_jellium_cell(cell, n, modes, rs, ok, coupling)
        type(jellium_cell_t), intent(out) :: cell
        integer, intent(in) :: n, modes(:, :)
        real(dp), intent(in) :: rs
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: coupling
        type(units_t) :: units
        type(species_t) :: e
        type(plasma_scales_t) :: unit_cell, scales
        integer :: a

        ok = .false.
        if (n < 1 .or. size(modes, 1) /= 3 .or. 2*size(modes, 2) > 62) return
        if (n > 2*size(modes, 2) .or. .not. (rs > 0)) return
        units = atomic_units(); e = electron()
        unit_cell = plasma_scales(units, e, real(n, dp), 1.0_dp)
        cell%length = rs/unit_cell%rs
        cell%volume = cell%length**3
        scales = plasma_scales(units, e, real(n, dp)/cell%volume, 1.0_dp)
        cell%rs = scales%rs; cell%fermi = scales%fermi
        cell%hbar = units%hbar; cell%mass = e%mass
        cell%coupling = pair_coupling(units, e, e)
        if (present(coupling)) cell%coupling = coupling
        cell%n = n; cell%nmodes = size(modes, 2); cell%m = 2*cell%nmodes
        cell%modes = modes
        allocate (cell%momentum(3, cell%m), cell%spin(cell%m))
        do a = 1, cell%nmodes
            cell%momentum(:, 2*a - 1) = modes(:, a); cell%spin(2*a - 1) = 1
            cell%momentum(:, 2*a) = modes(:, a); cell%spin(2*a) = -1
        end do
        ok = .true.
    end subroutine

    !> Spin orbital of mode vector k and spin s (+1 up, -1 down), 0 if absent.
    integer function jellium_orbital(cell, k, s) result(o)
        type(jellium_cell_t), intent(in) :: cell
        integer, intent(in) :: k(3), s
        integer :: a
        o = 0
        do a = 1, cell%nmodes
            if (all(cell%modes(:, a) == k)) then
                o = 2*a - 1
                if (s < 0) o = 2*a
                return
            end if
        end do
    end function

    !> Combinations of `count` modes as mode-bit masks (bit a-1 = mode a).
    subroutine mode_combinations(nmodes, count, masks, total)
        integer, intent(in) :: nmodes, count
        integer(i64), allocatable, intent(out) :: masks(:)
        integer, intent(out) :: total
        integer(i64) :: mask, low, ripple, last
        integer :: k
        total = 1
        do k = 1, count
            total = total*(nmodes - k + 1)/k
        end do
        allocate (masks(total))
        if (count == 0) then
            masks(1) = 0; return
        end if
        mask = ishft(1_i64, count) - 1; last = ishft(mask, nmodes - count)
        do k = 1, total
            masks(k) = mask
            if (mask == last) exit
            low = iand(mask, -mask); ripple = mask + low
            mask = ior(ripple, ishft(ieor(mask, ripple)/low, -2))
        end do
    end subroutine

    pure function mode_momentum(cell, mask) result(k)
        type(jellium_cell_t), intent(in) :: cell
        integer(i64), intent(in) :: mask
        integer :: k(3), a
        k = 0
        do a = 1, cell%nmodes
            if (btest(mask, a - 1)) k = k + cell%modes(:, a)
        end do
    end function

    pure integer(i64) function spin_mask(cell, up, down) result(mask)
        type(jellium_cell_t), intent(in) :: cell
        integer(i64), intent(in) :: up, down
        integer :: a
        mask = 0
        do a = 1, cell%nmodes
            if (btest(up, a - 1)) mask = ibset(mask, 2*a - 2)
            if (btest(down, a - 1)) mask = ibset(mask, 2*a - 1)
        end do
    end function

    !> Sector (or union of sectors) with nup up electrons and total momentum
    !> in the list ks(3, :); the basis is the sorted determinant list.
    subroutine jellium_sector(cell, nup, ks, sector, ok)
        type(jellium_cell_t), intent(in) :: cell
        integer, intent(in) :: nup, ks(:, :)
        type(jellium_sector_t), intent(out) :: sector
        logical, intent(out) :: ok
        integer(i64), allocatable :: up(:), down(:), masks(:), sorted(:)
        integer, allocatable :: kup(:, :), kdn(:, :), perm(:)
        integer :: nu, nd, i, j, l, count
        character(80) :: status

        ok = .false.
        if (nup < 0 .or. nup > cell%n .or. nup > cell%nmodes) return
        if (cell%n - nup > cell%nmodes) return
        call mode_combinations(cell%nmodes, nup, up, nu)
        call mode_combinations(cell%nmodes, cell%n - nup, down, nd)
        allocate (kup(3, nu), kdn(3, nd))
        do i = 1, nu
            kup(:, i) = mode_momentum(cell, up(i))
        end do
        do j = 1, nd
            kdn(:, j) = mode_momentum(cell, down(j))
        end do
        count = 0
        do l = 1, size(ks, 2)
            do i = 1, nu
                do j = 1, nd
                    if (all(kup(:, i) + kdn(:, j) == ks(:, l))) count = count + 1
                end do
            end do
        end do
        if (count == 0) return
        allocate (masks(count), perm(count), sorted(count))
        count = 0
        do l = 1, size(ks, 2)
            do i = 1, nu
                do j = 1, nd
                    if (.not. all(kup(:, i) + kdn(:, j) == ks(:, l))) cycle
                    count = count + 1; masks(count) = spin_mask(cell, up(i), down(j))
                end do
            end do
        end do
        call sort_key_order(masks, perm)
        sorted = masks(perm)
        call init_fock_sector_basis(sector%basis, cell%m, cell%n, sorted, &
                                    model_budget, ok, status)
        if (.not. ok) return
        sector%nup = nup; sector%d = count; sector%k = ks(:, 1)
    end subroutine

    !> Distinct total momenta (with nonempty sectors) for nup up electrons.
    subroutine jellium_sector_momenta(cell, nup, ks)
        type(jellium_cell_t), intent(in) :: cell
        integer, intent(in) :: nup
        integer, allocatable, intent(out) :: ks(:, :)
        integer(i64), allocatable :: up(:), down(:), keys(:)
        integer, allocatable :: perm(:), all_k(:, :)
        integer :: nu, nd, i, j, count, kk(3), b
        integer(i64), parameter :: span = 1000

        call mode_combinations(cell%nmodes, nup, up, nu)
        call mode_combinations(cell%nmodes, cell%n - nup, down, nd)
        allocate (keys(nu*nd), all_k(3, nu*nd), perm(nu*nd))
        count = 0
        do i = 1, nu
            do j = 1, nd
                count = count + 1
                kk = mode_momentum(cell, up(i)) + mode_momentum(cell, down(j))
                all_k(:, count) = kk
                keys(count) = ((kk(1) + span/2)*span + kk(2) + span/2)*span + kk(3) + span/2
            end do
        end do
        call sort_key_order(keys, perm)
        b = 0
        do i = 1, count
            if (i > 1) then
                if (keys(perm(i)) == keys(perm(i - 1))) cycle
            end if
            b = b + 1
        end do
        allocate (ks(3, b))
        b = 0
        do i = 1, count
            if (i > 1) then
                if (keys(perm(i)) == keys(perm(i - 1))) cycle
            end if
            b = b + 1; ks(:, b) = all_k(:, perm(i))
        end do
    end subroutine
    !> Jellium Hamiltonian at the cell coupling. The coefficient arrays do not
    !> depend on the sector, so one model serves every sector basis of the cell.
    subroutine jellium_hamiltonian(cell, basis, model, ok, coupling)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_basis_t), intent(in) :: basis
        type(fock_hamiltonian_t), intent(out) :: model
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: coupling
        character(120) :: status
        real(dp) :: g
        g = cell%coupling
        if (present(coupling)) g = coupling
        call prepare_jellium_hamiltonian(model, basis, cell%momentum, cell%spin, &
            cell%length, cell%hbar, cell%mass, g, model_budget, ok, status)
    end subroutine

    !> dH/dg = V1, the interaction at unit coupling (h = 0).
    subroutine jellium_interaction_model(cell, basis, model, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_basis_t), intent(in) :: basis
        type(fock_hamiltonian_t), intent(out) :: model
        logical, intent(out) :: ok
        type(fock_hamiltonian_t) :: unit
        real(dp), allocatable :: h(:, :), v(:, :)
        character(120) :: status
        integer :: p
        call jellium_hamiltonian(cell, basis, unit, ok, 1.0_dp)
        if (.not. ok) return
        p = cell%m*(cell%m - 1)/2
        allocate (h(cell%m, cell%m), v(p, p))
        call get_fock_integrals(unit, h, v, ok)
        if (.not. ok) return
        h = 0
        call prepare_fock_hamiltonian(model, basis, h, v, model_budget, ok, status)
    end subroutine

    !> One-body operator A = integral cos(q.x) n(x) d^3x: <k+-q,s|A|k,s> = 1/2.
    subroutine jellium_density_model(cell, basis, q, model, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: q(3)
        type(fock_hamiltonian_t), intent(out) :: model
        logical, intent(out) :: ok
        real(dp), allocatable :: h(:, :), v(:, :)
        character(120) :: status
        integer :: i, j, p
        p = cell%m*(cell%m - 1)/2
        allocate (h(cell%m, cell%m), v(p, p))
        h = 0; v = 0
        do j = 1, cell%m
            do i = 1, cell%m
                if (cell%spin(i) /= cell%spin(j)) cycle
                if (all(cell%momentum(:, i) - cell%momentum(:, j) == q) .or. &
                    all(cell%momentum(:, j) - cell%momentum(:, i) == q)) &
                    call jc_half(1.0_dp, h(i, j))
            end do
        end do
        call prepare_fock_hamiltonian(model, basis, h, v, model_budget, ok, status)
    end subroutine

    !> Lowest nev eigenpairs in one sector by FortNum block Lanczos on the
    !> matrix-free coupled action (block 8, full reorthogonalization).
    subroutine jellium_lowest(model, basis, nev, values, vectors, residuals, ok, &
                              iterations)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        integer, intent(in) :: nev
        real(dp), intent(out) :: values(:), vectors(:, :), residuals(:)
        logical, intent(out) :: ok
        integer, intent(out) :: iterations
        real(dp), allocatable :: start(:, :), zero(:), scratch(:)
        integer :: m, n, d, p, b, i, j, info, kmax, used
        integer(i64) :: bytes
        logical :: apply_ok

        ok = .false.; iterations = 0
        call fock_basis_info(basis, m, n, d, p, bytes, apply_ok)
        if (.not. apply_ok .or. nev < 1 .or. nev > d) return
        b = min(8, d)
        allocate (start(d, b), zero(d), scratch(d))
        zero = 0
        do j = 1, b
            do i = 1, d
                start(i, j) = cos(0.731_dp*i*j + 0.37_dp*j) + 0.01_dp*j
            end do
        end do
        apply_ok = .true.
        ! Without restarts, a nonconverged basis is retried from the same start
        ! with a larger (eventually complete) Krylov space.
        kmax = min(d, max(160, 24*nev))
        do
            call block_lanczos_lowest(apply, start, nev, kmax, 1.0e-10_dp, values, &
                                      vectors, residuals, info, used)
            iterations = iterations + used
            if (info == KRYLOV_OK .or. kmax == d) exit
            kmax = min(d, 3*kmax)
        end do
        ok = info == KRYLOV_OK .and. apply_ok
    contains
        subroutine apply(x, y)
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            logical :: step_ok
            call fock_matvec(model, basis, x, zero, y, scratch, step_ok)
            apply_ok = apply_ok .and. step_ok
        end subroutine
    end subroutine

    !> Lowest min(nev, d) states of every (nup, K) sector for nup = 0..N.
    subroutine jellium_build_spectrum(cell, model, nev, spectrum, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        integer, intent(in) :: nev
        type(jellium_spectrum_t), intent(out) :: spectrum
        logical, intent(out) :: ok
        integer, allocatable :: ks(:, :)
        real(dp), allocatable :: values(:), vectors(:, :), residuals(:), grow(:)
        integer(i64), allocatable :: grow_offset(:)
        integer, allocatable :: grow_int(:)
        real(dp), allocatable :: grow_real(:)
        integer :: nup, l, total, s, k, iterations, d
        integer(i64) :: used

        ok = .false.
        total = 0
        do nup = 0, cell%n
            if (nup > cell%nmodes .or. cell%n - nup > cell%nmodes) cycle
            call jellium_sector_momenta(cell, nup, ks)
            total = total + size(ks, 2)
        end do
        allocate (spectrum%sectors(total), spectrum%kept(total), spectrum%top(total))
        allocate (spectrum%sector_of(64), spectrum%offset(64), spectrum%energy(64), &
                  spectrum%residual(64), spectrum%flat(4096))
        s = 0; used = 0
        do nup = 0, cell%n
            if (nup > cell%nmodes .or. cell%n - nup > cell%nmodes) cycle
            call jellium_sector_momenta(cell, nup, ks)
            do l = 1, size(ks, 2)
                s = s + 1
                call jellium_sector(cell, nup, ks(:, l:l), spectrum%sectors(s), ok)
                if (.not. ok) return
                d = spectrum%sectors(s)%d
                k = min(nev, d)
                if (allocated(values)) deallocate (values, vectors, residuals)
                allocate (values(k), vectors(d, k), residuals(k))
                call jellium_lowest(model, spectrum%sectors(s)%basis, k, values, &
                                    vectors, residuals, ok, iterations)
                if (.not. ok) return
                spectrum%kept(s) = k; spectrum%top(s) = values(k)
                do while (spectrum%count + k > size(spectrum%energy))
                    call grow_states()
                end do
                do while (used + int(d, i64)*k > size(spectrum%flat, kind=i64))
                    allocate (grow(2*size(spectrum%flat)))
                    grow(1:used) = spectrum%flat(1:used)
                    call move_alloc(grow, spectrum%flat)
                end do
                do iterations = 1, k
                    spectrum%count = spectrum%count + 1
                    spectrum%sector_of(spectrum%count) = s
                    spectrum%offset(spectrum%count) = used
                    spectrum%energy(spectrum%count) = values(iterations)
                    spectrum%residual(spectrum%count) = residuals(iterations)
                    spectrum%flat(used + 1:used + d) = vectors(:, iterations)
                    used = used + d
                end do
            end do
        end do
        spectrum%nsectors = s
        ok = .true.
    contains
        subroutine grow_states()
            integer :: c
            c = size(spectrum%energy)
            allocate (grow_int(2*c)); grow_int(1:c) = spectrum%sector_of
            call move_alloc(grow_int, spectrum%sector_of)
            allocate (grow_offset(2*c)); grow_offset(1:c) = spectrum%offset
            call move_alloc(grow_offset, spectrum%offset)
            allocate (grow_real(2*c)); grow_real(1:c) = spectrum%energy
            call move_alloc(grow_real, spectrum%energy)
            allocate (grow_real(2*c)); grow_real(1:c) = spectrum%residual
            call move_alloc(grow_real, spectrum%residual)
        end subroutine
    end subroutine

    !> Canonical weights at Theta = T/T_F (normalized) and the omitted-weight
    !> estimate sum_s (d_s - kept_s) exp(-beta (top_s - E0))/Z, which bounds
    !> the omitted Boltzmann weight if every kept set is the lowest of its
    !> sector (candidate ordering claim). Theta <= 0: equal weights on the
    !> ground multiplet (energies within 1e-9 of E0).
    subroutine jellium_thermal_weights(cell, spectrum, theta, weights, truncation, &
                                       ground)
        type(jellium_cell_t), intent(in) :: cell
        type(jellium_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: weights(:), truncation, ground
        real(dp) :: beta, z, omitted
        type(interval_t) :: bound
        integer :: j, s

        ground = minval(spectrum%energy(1:spectrum%count))
        truncation = 0
        if (.not. (theta > 0)) then
            weights = 0
            where (spectrum%energy(1:spectrum%count) <= ground + &
                   1.0e-9_dp*max(1.0_dp, abs(ground))) weights = 1
            weights = weights/sum(weights)
            return
        end if
        beta = 1/(theta*cell%fermi)
        do j = 1, spectrum%count
            call jc_boltzmann(spectrum%energy(j), ground, beta, weights(j))
        end do
        z = sum(weights)
        omitted = 0
        do s = 1, spectrum%nsectors
            if (spectrum%kept(s) >= spectrum%sectors(s)%d) cycle
            call jc_truncation_interval(ipoint(real(spectrum%sectors(s)%d - spectrum%kept(s), &
                               dp)), ipoint(spectrum%top(s)), ipoint(ground), &
                               ipoint(beta), bound)
            omitted = omitted + bound%hi
        end do
        weights = weights/z
        truncation = omitted/z
    end subroutine
    !> Weighted real 1- and 2-RDMs of the ensemble (states with weight above
    !> threshold*max(weight)), via the engine's factorized RDM query.
    subroutine jellium_thermal_rdms(cell, spectrum, weights, gamma, gamma2, ok, &
                                    threshold)
        type(jellium_cell_t), intent(in) :: cell
        type(jellium_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: gamma(:, :), gamma2(:, :)
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: threshold
        real(dp), allocatable :: gr(:, :), gi(:, :), g2r(:, :), g2i(:, :), zero(:)
        real(dp) :: cut
        integer :: j, s, d, p
        integer(i64) :: o

        ok = .false.; gamma = 0; gamma2 = 0
        cut = 1.0e-12_dp
        if (present(threshold)) cut = threshold
        cut = cut*maxval(weights(1:spectrum%count))
        p = cell%m*(cell%m - 1)/2
        allocate (gr(cell%m, cell%m), gi(cell%m, cell%m), g2r(p, p), g2i(p, p))
        do j = 1, spectrum%count
            if (weights(j) <= cut) cycle
            s = spectrum%sector_of(j); d = spectrum%sectors(s)%d; o = spectrum%offset(j)
            if (allocated(zero)) deallocate (zero)
            allocate (zero(d)); zero = 0
            call fock_rdms(spectrum%sectors(s)%basis, spectrum%flat(o + 1:o + d), &
                           zero, gr, gi, g2r, g2i, ok)
            if (.not. ok) return
            gamma = gamma + weights(j)*gr
            gamma2 = gamma2 + weights(j)*g2r
        end do
        ok = .true.
    end subroutine

    !> Antisymmetrized pair-matrix coefficient v(p,q) of the radial-shell pair
    !> operator O_b = sum u_ijkl a_i^+ a_j^+ a_l a_k with
    !> u_ijkl = w_b(k_l-k_j) [s_i=s_k, s_j=s_l, channel] [k_i+k_j=k_k+k_l];
    !> channel 1 = same spin (both spins), 2 = opposite spin (both orders).
    !> w_b(Q) is the shell average of exp(iQ.r) over ra<|r|<rb (rb <= L/2).
    real(dp) function jellium_bin_operator(cell, ra, rb, channel, i, j, k, l) result(v)
        type(jellium_cell_t), intent(in) :: cell
        real(dp), intent(in) :: ra, rb
        integer, intent(in) :: channel, i, j, k, l
        v = term(i, j, k, l) - term(j, i, k, l) - term(i, j, l, k) + term(j, i, l, k)
    contains
        real(dp) function term(a, b, c, e) result(u)
            integer, intent(in) :: a, b, c, e
            integer :: transfer(3)
            real(dp) :: qn, fa, fb
            u = 0
            if (cell%spin(a) /= cell%spin(c) .or. cell%spin(b) /= cell%spin(e)) return
            if (channel == 1 .neqv. cell%spin(a) == cell%spin(b)) return
            if (any(cell%momentum(:, a) + cell%momentum(:, b) /= &
                    cell%momentum(:, c) + cell%momentum(:, e))) return
            transfer = cell%momentum(:, e) - cell%momentum(:, b)
            if (all(transfer == 0)) then
                u = 1; return
            end if
            qn = 2*acos(-1.0_dp)*sqrt(real(sum(transfer**2), dp))/cell%length
            call shell_moment(qn*ra, fa); call shell_moment(qn*rb, fb)
            call jc_shell_average(fa, fb, qn, ra, rb, u)
        end function
    end function

    !> Shell-averaged pair correlation g_same, g_opp from a (transition) 2-RDM:
    !> g = sum_pq v_b(p,q) Gamma(q,p) / (2 (N/2)^2) for both spin channels.
    subroutine jellium_pair_bins(cell, gamma2, edges, g_same, g_opp)
        type(jellium_cell_t), intent(in) :: cell
        real(dp), intent(in) :: gamma2(:, :), edges(:)
        real(dp), intent(out) :: g_same(:), g_opp(:)
        integer, allocatable :: pi(:), pj(:)
        integer :: p, q, b, np
        real(dp) :: half, value
        logical :: ok

        np = cell%m*(cell%m - 1)/2
        allocate (pi(np), pj(np))
        call pair_lists(cell%m, pi, pj)
        g_same = 0; g_opp = 0
        do q = 1, np
            do p = 1, np
                if (.not. abs(gamma2(q, p)) > 0) cycle
                do b = 1, size(edges) - 1
                    g_same(b) = g_same(b) + gamma2(q, p)*jellium_bin_operator(cell, &
                        edges(b), edges(b + 1), 1, pi(p), pj(p), pi(q), pj(q))
                    g_opp(b) = g_opp(b) + gamma2(q, p)*jellium_bin_operator(cell, &
                        edges(b), edges(b + 1), 2, pi(p), pj(p), pi(q), pj(q))
                end do
            end do
        end do
        half = 0.5_dp*cell%n
        do b = 1, size(edges) - 1
            call jc_pair_correlation(g_same(b), half, real(cell%n, dp), value)
            g_same(b) = value
            call jc_pair_correlation(g_opp(b), half, real(cell%n, dp), value)
            g_opp(b) = value
        end do
        ok = .true.
    end subroutine

    subroutine shell_moment(x, f)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: f
        if (x < 0.1_dp) then
            call jc_shell_moment_series(x, f)
        else
            call jc_shell_moment(x, f)
        end if
    end subroutine

    subroutine pair_lists(m, pi, pj)
        integer, intent(in) :: m
        integer, intent(out) :: pi(:), pj(:)
        integer :: i, j, p
        p = 0
        do i = 1, m - 1
            do j = i + 1, m
                p = p + 1; pi(p) = i; pj(p) = j
            end do
        end do
    end subroutine

    !> Pair and three-body cumulant norms of the ensemble (see
    !> quantum_cumulant_norms): norms = [|Gamma2|, |lambda2|, |Gamma3|, |lambda3|].
    subroutine jellium_cumulants(cell, spectrum, weights, norms, ok, threshold)
        type(jellium_cell_t), intent(in) :: cell
        type(jellium_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: norms(4)
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: threshold
        real(dp), allocatable :: gamma(:, :), gamma2(:, :), zero(:), vr(:), vi(:)
        real(dp), allocatable :: val(:), grow(:), g3(:)
        integer(i64), allocatable :: tk(:), uk(:), key(:), growk(:), t3(:), u3(:)
        integer, allocatable :: perm(:)
        integer(i64) :: m3, o
        real(dp) :: cut
        integer :: j, s, d, count, used, x, y, nred

        norms = 0; ok = .false.
        cut = 1.0e-10_dp
        if (present(threshold)) cut = threshold
        allocate (gamma(cell%m, cell%m), gamma2(cell%m*(cell%m - 1)/2, &
                  cell%m*(cell%m - 1)/2))
        call jellium_thermal_rdms(cell, spectrum, weights, gamma, gamma2, ok, cut)
        if (.not. ok) return
        m3 = int(cell%m, i64)**3
        allocate (key(1024), val(1024))
        used = 0
        cut = cut*maxval(weights(1:spectrum%count))
        do j = 1, spectrum%count
            if (weights(j) <= cut) cycle
            s = spectrum%sector_of(j); d = spectrum%sectors(s)%d; o = spectrum%offset(j)
            if (allocated(zero)) deallocate (zero)
            allocate (zero(d)); zero = 0
            call fock_rdm3_entries(spectrum%sectors(s)%basis, &
                spectrum%flat(o + 1:o + d), zero, tk, uk, vr, vi, count, ok)
            if (.not. ok) return
            do while (used + count > size(key))
                allocate (growk(2*size(key))); growk(1:used) = key(1:used)
                call move_alloc(growk, key)
                allocate (grow(2*size(val))); grow(1:used) = val(1:used)
                call move_alloc(grow, val)
            end do
            key(used + 1:used + count) = tk(1:count)*m3 + uk(1:count)
            val(used + 1:used + count) = weights(j)*vr(1:count)
            used = used + count
        end do
        allocate (perm(used), t3(used), u3(used), g3(used))
        call sort_key_order(key(1:used), perm)
        nred = 0; x = 1
        do while (x <= used)
            nred = nred + 1
            t3(nred) = key(perm(x))/m3; u3(nred) = mod(key(perm(x)), m3); g3(nred) = 0
            y = x
            do while (y <= used)
                if (key(perm(y)) /= key(perm(x))) exit
                g3(nred) = g3(nred) + val(perm(y)); y = y + 1
            end do
            x = y
        end do
        call quantum_cumulant_norms(cell%m, gamma, gamma2, t3(1:nred), u3(1:nred), &
            g3(1:nred), norms(1), norms(2), norms(3), norms(4), ok)
    end subroutine

    !> Solve P (H - shift) P x = P rhs by restarted GMRES on the coupled action
    !> (P projects out `project` when present). x starts at zero.
    subroutine shifted_solve(model, basis, shift, rhs, x, ok, iterations, project)
        type(fock_hamiltonian_t), intent(in) :: model
        type(fock_basis_t), intent(in) :: basis
        real(dp), intent(in) :: shift, rhs(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        integer, intent(out) :: iterations
        real(dp), intent(in), optional :: project(:)
        real(dp), allocatable :: zero(:), scratch(:), b(:), w(:)
        real(dp) :: residual
        integer :: info
        logical :: apply_ok

        allocate (zero(size(rhs)), scratch(size(rhs)), b(size(rhs)), w(size(rhs)))
        zero = 0; x = 0; b = rhs
        if (present(project)) b = b - dot_product(project, b)*project
        apply_ok = .true.
        call real_gmres_operator(apply, b, x, 1.0e-11_dp, 4000, 120, info, &
                                 iterations, residual)
        if (present(project)) x = x - dot_product(project, x)*project
        ok = info == KRYLOV_OK .and. apply_ok
    contains
        subroutine apply(v, y)
            real(dp), intent(in) :: v(:)
            real(dp), intent(out) :: y(:)
            logical :: step_ok
            w = v
            if (present(project)) w = w - dot_product(project, w)*project
            call fock_matvec(model, basis, w, zero, y, scratch, step_ok)
            apply_ok = apply_ok .and. step_ok
            y = y - shift*w
            if (present(project)) y = y - dot_product(project, y)*project
        end subroutine
    end subroutine

    !> Isothermal static density response per volume of the ensemble,
    !> chi = sum_n w_n chi_n, chi_n = -(4/V) <A psi_n|(H-E_n)^(-1)|A psi_n>.
    !> A psi_n lies in the K+-q sectors, orthogonal to psi_n, so <n|A|n> = 0 and
    !> the canonical -d2F/dlambda2 reduces to this weighted sum. Each solve is
    !> GMRES on the union K+-q basis; chi_state(j) is chi_n (0 if skipped).
    subroutine jellium_static_response(cell, model, spectrum, weights, q, chi, &
                                       chi_state, ok, threshold)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        type(jellium_spectrum_t), intent(in) :: spectrum
        real(dp), intent(in) :: weights(:)
        integer, intent(in) :: q(3)
        real(dp), intent(out) :: chi, chi_state(:)
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: threshold
        type(jellium_sector_t), allocatable :: unions(:)
        logical, allocatable :: built(:), empty(:)
        type(fock_hamiltonian_t) :: density
        real(dp), allocatable :: phi(:), phi_im(:), x(:), zero(:)
        real(dp) :: cut, value
        integer :: j, s, d, du, iterations, ks(3, 2)
        integer(i64) :: o

        ok = .false.; chi = 0; chi_state = 0
        cut = 1.0e-10_dp
        if (present(threshold)) cut = threshold
        cut = cut*maxval(weights(1:spectrum%count))
        call jellium_density_model(cell, spectrum%sectors(1)%basis, q, density, ok)
        if (.not. ok) return
        allocate (unions(spectrum%nsectors), built(spectrum%nsectors), &
                  empty(spectrum%nsectors))
        built = .false.; empty = .false.
        do j = 1, spectrum%count
            if (weights(j) <= cut) cycle
            s = spectrum%sector_of(j); d = spectrum%sectors(s)%d; o = spectrum%offset(j)
            if (.not. built(s)) then
                ks(:, 1) = spectrum%sectors(s)%k + q; ks(:, 2) = spectrum%sectors(s)%k - q
                call jellium_sector(cell, spectrum%sectors(s)%nup, ks, unions(s), ok)
                empty(s) = .not. ok
                built(s) = .true.
            end if
            if (empty(s)) cycle
            du = unions(s)%d
            if (allocated(phi)) deallocate (phi, phi_im, x, zero)
            allocate (phi(du), phi_im(du), x(du), zero(d))
            zero = 0
            call fock_matvec_between(density, spectrum%sectors(s)%basis, &
                unions(s)%basis, spectrum%flat(o + 1:o + d), zero, phi, phi_im, ok)
            if (.not. ok) return
            call shifted_solve(model, unions(s)%basis, spectrum%energy(j), phi, x, &
                               ok, iterations)
            if (.not. ok) return
            call jc_susceptibility(dot_product(phi, x), cell%volume, value)
            chi_state(j) = value
            chi = chi + weights(j)*value
        end do
        ok = .true.
    end subroutine

    !> Ideal (g = 0) response of the same finite cell and orbital set:
    !> canonical chi0 = sum_D p_D S_D/V over all N-determinants D (exact
    !> canonical ensemble; Theta <= 0 averages the ideal ground multiplet),
    !> grand-canonical Fermi-Dirac chi0 (Theta > 0) and the RPA forms
    !> chi0/(1 - c S0), c = 4 pi g/(V q^2).
    subroutine jellium_ideal_response(cell, theta, q, chi_can, chi_gc, rpa_can, &
                                      rpa_gc, ok)
        type(jellium_cell_t), intent(in) :: cell
        real(dp), intent(in) :: theta
        integer, intent(in) :: q(3)
        real(dp), intent(out) :: chi_can, chi_gc, rpa_can, rpa_gc
        logical, intent(out) :: ok
        real(dp) :: eps(62), occ(62), beta, ground, z, energy, cellsum, w, term
        real(dp) :: lo, hi, mu, coefficient, total
        integer :: partner(62), o, total_states, k, it
        integer(i64), allocatable :: masks(:)
        logical :: multiplet

        ok = .false.; chi_can = 0; chi_gc = 0; rpa_can = 0; rpa_gc = 0
        do o = 1, cell%m
            call qf_kinetic(cell%length, cell%hbar, cell%mass, &
                            real(sum(cell%momentum(:, o)**2), dp), eps(o))
            partner(o) = jellium_orbital(cell, cell%momentum(:, o) + q, cell%spin(o))
        end do
        call bit_combinations(cell%m, cell%n, masks, total_states)
        multiplet = .not. (theta > 0)
        beta = 0
        if (.not. multiplet) beta = 1/(theta*cell%fermi)
        ground = huge(1.0_dp)
        do k = 1, total_states
            ground = min(ground, determinant_energy(masks(k)))
        end do
        z = 0
        do k = 1, total_states
            energy = determinant_energy(masks(k))
            if (multiplet) then
                w = 0
                if (energy <= ground + 1.0e-12_dp*max(1.0_dp, abs(ground))) w = 1
            else
                call jc_boltzmann(energy, ground, beta, w)
            end if
            if (.not. w > 0) cycle
            occ = 0
            do o = 1, cell%m
                if (btest(masks(k), o - 1)) occ(o) = 1
            end do
            z = z + w
            chi_can = chi_can + w*lindhard_sum(occ)
        end do
        chi_can = chi_can/z
        call qf_coulomb(cell%length, cell%coupling, real(sum(q**2), dp), coefficient)
        call jc_rpa(chi_can/cell%volume, coefficient, chi_can, rpa_can)
        chi_can = chi_can/cell%volume
        if (.not. multiplet) then
            lo = minval(eps(1:cell%m)) - 50/beta; hi = maxval(eps(1:cell%m)) + 50/beta
            do it = 1, 200
                mu = 0.5_dp*(lo + hi)
                total = 0
                do o = 1, cell%m
                    call jc_fermi(eps(o), mu, beta, occ(o)); total = total + occ(o)
                end do
                if (total > cell%n) then
                    hi = mu
                else
                    lo = mu
                end if
            end do
            cellsum = lindhard_sum(occ)
            call jc_rpa(cellsum/cell%volume, coefficient, cellsum, rpa_gc)
            chi_gc = cellsum/cell%volume
        end if
        ok = ieee_is_finite(chi_can) .and. ieee_is_finite(rpa_can)
    contains
        real(dp) function determinant_energy(mask) result(e)
            integer(i64), intent(in) :: mask
            integer :: a
            e = 0
            do a = 1, cell%m
                if (btest(mask, a - 1)) e = e + eps(a)
            end do
        end function
        real(dp) function lindhard_sum(f) result(sum0)
            real(dp), intent(in) :: f(:)
            integer :: a
            sum0 = 0
            do a = 1, cell%m
                if (partner(a) == 0) cycle
                call jc_lindhard_term(f(a), f(partner(a)), eps(a), eps(partner(a)), term)
                sum0 = sum0 + term
            end do
        end function
    end subroutine

    subroutine bit_combinations(m, n, masks, total)
        integer, intent(in) :: m, n
        integer(i64), allocatable, intent(out) :: masks(:)
        integer, intent(out) :: total
        integer(i64) :: mask, low, ripple
        integer :: k
        total = 1
        do k = 1, min(n, m - n)
            total = int(int(total, i64)*(m - k + 1)/k)
        end do
        allocate (masks(total))
        mask = ishft(1_i64, n) - 1
        do k = 1, total
            masks(k) = mask
            if (k == total) exit
            low = iand(mask, -mask); ripple = mask + low
            mask = ior(ripple, ishft(ieor(mask, ripple)/low, -2))
        end do
    end subroutine

    !> Density-mode history after a weak kick kappa*A delta(t) of the state
    !> psi0 (energy e0, sector `sector`): to first order in kappa,
    !> delta<A>(t)/kappa = 2 Im(exp(i e0 t) C(t)), C(t) = <phi0|exp(-iHt) phi0>,
    !> phi0 = A psi0 (in the K+-q union basis). The untrusted candidate is
    !> fixed-step Cash-Karp RK5(4) on the matrix-free coupled action; FortNum's
    !> Hermite/Bernstein residual certificate with the strict interval action
    !> -i H x/hbar (skew-adjoint, mu = 0) bounds ||phi(t) - phi_h(t)|| on every
    !> closed cell. band(k) = 2 ||phi0|| radius(k) bounds the response error of
    !> the finite-basis model for the supplied psi0, e0 (certified-discrete;
    !> the Ritz errors of psi0/e0 and the first-order kick are separate).
    subroutine jellium_kick_evolution(cell, model, sector, psi0, e0, q, steps, dt, &
                                      times, response, band, ok, union_dim, &
                                      candidate_seconds, certificate_seconds)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_hamiltonian_t), intent(in) :: model
        type(jellium_sector_t), intent(in) :: sector
        real(dp), intent(in) :: psi0(:), e0, dt
        integer, intent(in) :: q(3), steps
        real(dp), intent(out) :: times(0:), response(0:), band(0:)
        logical, intent(out) :: ok
        integer, intent(out) :: union_dim
        real(dp), intent(out) :: candidate_seconds, certificate_seconds
        type(jellium_sector_t) :: union
        type(fock_hamiltonian_t) :: density
        type(fortnum_status_t) :: status
        real(dp), allocatable :: phi(:), phi_im(:), zero(:), y(:, :), f(:, :)
        real(dp), allocatable :: stage(:, :), residual(:), radius(:), hr(:), hi(:)
        type(interval_t), allocatable :: c(:, :), w(:, :), xr(:), xi(:), yr(:), yi(:)
        real(dp) :: t0, t1, norm0, cr, ci
        integer :: ks(3, 2), du, k, nfev
        logical :: apply_ok

        ok = .false.; union_dim = 0; candidate_seconds = 0; certificate_seconds = 0
        times = 0; response = 0; band = huge(1.0_dp)
        ks(:, 1) = sector%k + q; ks(:, 2) = sector%k - q
        call jellium_sector(cell, sector%nup, ks, union, ok)
        if (.not. ok) return
        du = union%d; union_dim = du
        call jellium_density_model(cell, sector%basis, q, density, ok)
        if (.not. ok) return
        allocate (phi(du), phi_im(du), zero(sector%d), hr(du), hi(du))
        zero = 0
        call fock_matvec_between(density, sector%basis, union%basis, psi0, zero, &
                                 phi, phi_im, ok)
        if (.not. ok) return
        norm0 = norm2(phi)
        allocate (y(2*du, 0:steps), f(2*du, 0:steps), stage(2*du, 9))
        allocate (residual(steps), radius(0:steps))
        allocate (c(2*du, 0:3), w(2*du, 0:3), xr(du), xi(du), yr(du), yi(du))
        y(1:du, 0) = phi; y(du + 1:, 0) = 0
        do k = 0, steps
            times(k) = k*dt
        end do
        apply_ok = .true.
        call cpu_time(t0)
        call rhs(0.0_dp, y(:, 0), f(:, 0))
        nfev = 0
        do k = 1, steps
            call cash_karp_step(rhs, times(k - 1), y(:, k - 1), dt, .false., &
                stage(:, 1), stage(:, 2), stage(:, 3), stage(:, 4), stage(:, 5), &
                stage(:, 6), stage(:, 7), stage(:, 8), stage(:, 9), nfev)
            y(:, k) = stage(:, 8)
            call rhs(times(k), y(:, k), f(:, k))
        end do
        call cpu_time(t1)
        candidate_seconds = t1 - t0
        if (.not. apply_ok) return
        call cpu_time(t0)
        call certify_linear_trajectory(action, times, y, f, ipoint(0.0_dp), 0.0_dp, &
                                       c, w, residual, radius, status)
        call cpu_time(t1)
        certificate_seconds = t1 - t0
        if (status%code /= FORTNUM_OK) return
        do k = 0, steps
            cr = dot_product(phi, y(1:du, k)); ci = dot_product(phi, y(du + 1:, k))
            call jc_kick_response(cr, ci, e0, times(k), response(k))
            band(k) = 2*norm0*radius(k)
        end do
        ok = .true.
    contains
        subroutine rhs(t, x, dxdt, ctx)
            real(dp), intent(in) :: t, x(:)
            real(dp), intent(out) :: dxdt(:)
            class(*), intent(in), optional :: ctx
            logical :: step_ok
            call fock_matvec(model, union%basis, x(1:du), x(du + 1:), hr, hi, step_ok)
            apply_ok = apply_ok .and. step_ok .and. ieee_is_finite(t)
            dxdt(1:du) = hi/cell%hbar; dxdt(du + 1:) = -hr/cell%hbar
        end subroutine
        subroutine action(x, ax, action_ok, ctx)
            type(interval_t), intent(in) :: x(:)
            type(interval_t), intent(out) :: ax(:)
            logical, intent(out) :: action_ok
            class(*), intent(in), optional :: ctx
            integer :: i
            xr = x(1:du); xi = x(du + 1:)
            call fock_matvec_between_interval(model, union%basis, union%basis, &
                                              xr, xi, yr, yi, action_ok)
            do i = 1, du
                ax(i) = yi(i)/ipoint(cell%hbar)
                ax(du + i) = (ipoint(0.0_dp) - yr(i))/ipoint(cell%hbar)
            end do
        end subroutine
    end subroutine

    !> Equation-level coupling derivative of a nondegenerate eigenstate:
    !> dE/dg = <psi|V1|psi> (Hellmann-Feynman), psi' = -(H-E)^+ (V1 - dE/dg) psi
    !> with psi' orthogonal to psi, and dg_b/dg = 2 <psi'|O_b|psi> from the
    !> transition 2-RDM. Returns the shell-averaged g as well.
    subroutine jellium_coupling_derivative(cell, basis, model, interaction, psi, &
                                           energy, edges, de, g_same, g_opp, &
                                           dg_same, dg_opp, ok, iterations)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_basis_t), intent(in) :: basis
        type(fock_hamiltonian_t), intent(in) :: model, interaction
        real(dp), intent(in) :: psi(:), energy, edges(:)
        real(dp), intent(out) :: de, g_same(:), g_opp(:), dg_same(:), dg_opp(:)
        logical, intent(out) :: ok
        integer, intent(out) :: iterations
        real(dp), allocatable :: v1(:), zero(:), scratch(:), x(:), rhs(:)
        real(dp), allocatable :: gr(:, :), gi(:, :), g2r(:, :), g2i(:, :)
        integer :: d, p

        d = size(psi); p = cell%m*(cell%m - 1)/2
        allocate (v1(d), zero(d), scratch(d), x(d), rhs(d))
        allocate (gr(cell%m, cell%m), gi(cell%m, cell%m), g2r(p, p), g2i(p, p))
        zero = 0
        call fock_matvec(interaction, basis, psi, zero, v1, scratch, ok)
        if (.not. ok) return
        de = dot_product(psi, v1)
        rhs = -(v1 - de*psi)
        call shifted_solve(model, basis, energy, rhs, x, ok, iterations, psi)
        if (.not. ok) return
        call fock_rdms(basis, psi, zero, gr, gi, g2r, g2i, ok)
        if (.not. ok) return
        call jellium_pair_bins(cell, g2r, edges, g_same, g_opp)
        call fock_transition_rdms(basis, x, zero, psi, zero, gr, gi, g2r, g2i, ok)
        if (.not. ok) return
        call jellium_pair_bins(cell, g2r, edges, dg_same, dg_opp)
        dg_same = 2*dg_same; dg_opp = 2*dg_opp
    end subroutine

    !> Adjoint (VJP) of (E, g_same, g_opp) with respect to g for a
    !> nondegenerate eigenstate: gbar = bar_e dE/dg + 2 <(V1-dE/dg) psi, lambda>,
    !> P (H-E) P lambda = -P O_bar psi, with O_bar = sum_b bar_b O_b / (N^2/2)
    !> applied matrix-free as a Fock pair operator.
    subroutine jellium_coupling_vjp(cell, basis, model, interaction, psi, energy, &
                                    edges, bar_e, bar_same, bar_opp, gbar, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(fock_basis_t), intent(in) :: basis
        type(fock_hamiltonian_t), intent(in) :: model, interaction
        real(dp), intent(in) :: psi(:), energy, edges(:), bar_e, bar_same(:), bar_opp(:)
        real(dp), intent(out) :: gbar
        logical, intent(out) :: ok
        type(fock_hamiltonian_t) :: observable
        real(dp), allocatable :: h(:, :), v(:, :), v1(:), zero(:), scratch(:)
        real(dp), allocatable :: opsi(:), lambda(:)
        integer, allocatable :: pi(:), pj(:)
        real(dp) :: de, scale
        integer :: d, p, a, b, k, iterations
        character(120) :: status

        d = size(psi); p = cell%m*(cell%m - 1)/2
        allocate (h(cell%m, cell%m), v(p, p), v1(d), zero(d), scratch(d), opsi(d), &
                  lambda(d), pi(p), pj(p))
        call pair_lists(cell%m, pi, pj)
        h = 0; v = 0; zero = 0
        scale = 2.0_dp/real(cell%n, dp)**2
        do b = 1, p
            do a = 1, p
                do k = 1, size(edges) - 1
                    v(a, b) = v(a, b) + scale*( &
                        bar_same(k)*jellium_bin_operator(cell, edges(k), edges(k + 1), &
                                                         1, pi(a), pj(a), pi(b), pj(b)) + &
                        bar_opp(k)*jellium_bin_operator(cell, edges(k), edges(k + 1), &
                                                        2, pi(a), pj(a), pi(b), pj(b)))
                end do
            end do
        end do
        v = 0.5_dp*(v + transpose(v))
        call prepare_fock_hamiltonian(observable, basis, h, v, model_budget, ok, status)
        if (.not. ok) return
        call fock_matvec(interaction, basis, psi, zero, v1, scratch, ok)
        if (.not. ok) return
        de = dot_product(psi, v1)
        call fock_matvec(observable, basis, psi, zero, opsi, scratch, ok)
        if (.not. ok) return
        call shifted_solve(model, basis, energy, -opsi, lambda, ok, iterations, psi)
        if (.not. ok) return
        gbar = bar_e*de + 2*dot_product(v1 - de*psi, lambda)
    end subroutine
end module kin6d_jellium_cell
