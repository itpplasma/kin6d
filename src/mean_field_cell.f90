!> Mean-field models of the Q4 jellium cell (campaign Q5): Hartree, restricted
!> and unrestricted Hartree-Fock (T = 0), restricted Mermin finite-temperature
!> HF and static TDHF (RPA with exchange) response, all in the SAME plane-wave
!> spin-orbital basis, N, spin sector, r_s and cutoff as the exact reference.
!> Orbitals are real combinations of the cell's plane-wave modes, one spatial
!> orbital set per spin (collinear UHF; no spin mixing, no complex orbitals).
!> The interaction is the generated jellium coefficient qf_coulomb (q = 0
!> removed), kinetic energies the generated qf_kinetic leaf, the density
!> perturbation the engine's jellium_density_model: the same definitions the
!> exact Fock-space Hamiltonian composes. A Slater state is expanded in the
!> engine's determinant basis to evaluate its full-parent residual.
!> Claims: finite-basis model quantities. A mean-field solution is a candidate
!> stationary point; the lowest of several starts is a candidate global
!> minimum. Temple bounds use a candidate (Lanczos-ordered) gap input.
module kin6d_mean_field_cell
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, ipoint
    use fortnum_symmetric_eigen, only: symmetric_eigen
    use kin6d_quantum_fock, only: fock_basis_t, fock_hamiltonian_t, &
        get_fock_integrals, fock_matvec, fock_basis_info, get_fock_states
    use kin6d_quantum_fock_leaves, only: qf_kinetic, qf_coulomb
    use kin6d_jellium_cell, only: jellium_cell_t, jellium_sector_t, &
        jellium_spectrum_t, jellium_sector, jellium_sector_momenta, &
        jellium_density_model, jellium_lowest
    use kin6d_jellium_cell_leaves, only: jc_fermi, jc_susceptibility
    use kin6d_mean_field_cell_leaves, only: mf_entropy, mf_free_energy, &
        mf_temple_interval
    implicit none
    private
    public :: mf_interaction_t, mf_state_t, mf_estimate_t
    public :: mf_interaction, mf_fock, mf_energy, mf_scf, mf_ground, mf_mermin
    public :: mf_slater_sector, mf_slater_vector, mf_as_spectrum, mf_residual
    public :: mf_static_response, mf_coupling_derivative, mf_start_orbitals
    public :: mf_frozen_energy

    !> Sparse momentum-conserving quartets (a, b, c, d) with k_a + k_b = k_c + k_d,
    !> k_a /= k_c and coefficient w = 4 pi g/(V |k_a - k_c|^2):
    !> V = (1/2) sum w c+_{a s} c+_{b s'} c_{d s'} c_{c s}.
    type :: mf_interaction_t
        integer :: nmodes = 0, count = 0
        integer, allocatable :: a(:), b(:), c(:), d(:)
        real(dp), allocatable :: w(:), kinetic(:)
        real(dp) :: coupling = 0
    end type

    !> Mean-field state: full orbital sets per spin (columns, real), orbital
    !> energies of the converged Fock operator and occupations (0/1 at T = 0,
    !> Fermi-Dirac for Mermin). energy = kinetic + hartree + exchange.
    type :: mf_state_t
        integer :: nup = 0, ndn = 0, iterations = 0
        logical :: restricted = .false., exchange = .true., converged = .false.
        logical :: aufbau = .false.
        real(dp), allocatable :: cup(:, :), cdn(:, :), eup(:), edn(:), fup(:), fdn(:)
        real(dp) :: energy = 0, kinetic = 0, hartree = 0, exchange_energy = 0
        real(dp) :: theta = 0, mu = 0, entropy = 0, free = 0, commutator = 0
    end type

    !> Full-parent error estimator of a Slater state Phi in the finite basis:
    !> rho = <Phi|H|Phi>, sigma = ||(H - rho) Phi||; e0/e1 the two lowest Ritz
    !> values (residuals r0/r1) of the symmetry space containing Phi;
    !> temple = rho - sigma^2/(beta - rho), beta = e1 - r1, when rho < beta.
    type :: mf_estimate_t
        real(dp) :: rho = 0, sigma = 0, norm = 0, e0 = 0, e1 = 0, r0 = 0, r1 = 0
        real(dp) :: beta = 0, temple = 0
        logical :: temple_valid = .false.
        integer :: dimension = 0, nk = 0
    end type
contains
    !> Quartet list and mode kinetic energies of the cell at coupling g
    !> (default: the cell coupling).
    subroutine mf_interaction(cell, inter, coupling)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_interaction_t), intent(out) :: inter
        real(dp), intent(in), optional :: coupling
        integer :: a, b, c, d, n, kd(3), q(3), pass
        real(dp) :: g

        g = cell%coupling
        if (present(coupling)) g = coupling
        n = cell%nmodes
        inter%nmodes = n; inter%coupling = g
        allocate (inter%kinetic(n))
        do a = 1, n
            call qf_kinetic(cell%length, cell%hbar, cell%mass, &
                            real(sum(cell%modes(:, a)**2), dp), inter%kinetic(a))
        end do
        do pass = 1, 2
            inter%count = 0
            do a = 1, n
                do c = 1, n
                    q = cell%modes(:, a) - cell%modes(:, c)
                    if (all(q == 0)) cycle
                    do b = 1, n
                        kd = cell%modes(:, b) + q
                        d = mode_index(cell, kd)
                        if (d == 0) cycle
                        inter%count = inter%count + 1
                        if (pass == 1) cycle
                        inter%a(inter%count) = a; inter%b(inter%count) = b
                        inter%c(inter%count) = c; inter%d(inter%count) = d
                        call qf_coulomb(cell%length, g, real(sum(q**2), dp), &
                                        inter%w(inter%count))
                    end do
                end do
            end do
            if (pass == 1) allocate (inter%a(inter%count), inter%b(inter%count), &
                inter%c(inter%count), inter%d(inter%count), inter%w(inter%count))
        end do
    end subroutine

    integer function mode_index(cell, k) result(index)
        type(jellium_cell_t), intent(in) :: cell
        integer, intent(in) :: k(3)
        integer :: a
        index = 0
        do a = 1, cell%nmodes
            if (all(cell%modes(:, a) == k)) then
                index = a; return
            end if
        end do
    end function

    !> Two-body mean-field matrices for real symmetric spin densities
    !> P_s(c, a) = <c+_{a s} c_{c s}>: Hartree G_H(a, c) = sum w n(d, b)
    !> (n = P_up + P_dn) and exchange G_X,s(a, d) = -sum w P_s(c, b).
    subroutine two_body(inter, pup, pdn, exchange, hartree, xup, xdn)
        type(mf_interaction_t), intent(in) :: inter
        real(dp), intent(in) :: pup(:, :), pdn(:, :)
        logical, intent(in) :: exchange
        real(dp), intent(out) :: hartree(:, :), xup(:, :), xdn(:, :)
        integer :: t, a, b, c, d
        real(dp) :: w
        hartree = 0; xup = 0; xdn = 0
        do t = 1, inter%count
            a = inter%a(t); b = inter%b(t); c = inter%c(t); d = inter%d(t)
            w = inter%w(t)
            hartree(a, c) = hartree(a, c) + w*(pup(d, b) + pdn(d, b))
            if (.not. exchange) cycle
            xup(a, d) = xup(a, d) - w*pup(c, b)
            xdn(a, d) = xdn(a, d) - w*pdn(c, b)
        end do
    end subroutine

    !> Fock matrices F_s = T + G_H + G_X,s and the energy components
    !> E = sum_s tr(T P_s) + (1/2) sum_s tr((G_H + G_X,s) P_s).
    subroutine mf_fock(inter, pup, pdn, exchange, fup, fdn, kinetic, hartree, &
                       exchange_energy)
        type(mf_interaction_t), intent(in) :: inter
        real(dp), intent(in) :: pup(:, :), pdn(:, :)
        logical, intent(in) :: exchange
        real(dp), intent(out) :: fup(:, :), fdn(:, :)
        real(dp), intent(out) :: kinetic, hartree, exchange_energy
        real(dp), allocatable :: gh(:, :), xu(:, :), xd(:, :)
        integer :: a, n
        n = inter%nmodes
        allocate (gh(n, n), xu(n, n), xd(n, n))
        call two_body(inter, pup, pdn, exchange, gh, xu, xd)
        kinetic = 0
        do a = 1, n
            kinetic = kinetic + inter%kinetic(a)*(pup(a, a) + pdn(a, a))
        end do
        hartree = 0.5_dp*sum(gh*(pup + pdn))
        exchange_energy = 0.5_dp*(sum(xu*pup) + sum(xd*pdn))
        fup = gh + xu; fdn = gh + xd
        do a = 1, n
            fup(a, a) = fup(a, a) + inter%kinetic(a)
            fdn(a, a) = fdn(a, a) + inter%kinetic(a)
        end do
    end subroutine

    !> Densities of a state: P_s = C_s diag(f_s) C_s^T.
    subroutine densities(state, pup, pdn)
        type(mf_state_t), intent(in) :: state
        real(dp), intent(out) :: pup(:, :), pdn(:, :)
        integer :: j
        pup = 0; pdn = 0
        do j = 1, size(state%fup)
            if (state%fup(j) > 0) pup = pup + state%fup(j)* &
                spread(state%cup(:, j), 2, size(pup, 1))*spread(state%cup(:, j), 1, size(pup, 1))
            if (state%fdn(j) > 0) pdn = pdn + state%fdn(j)* &
                spread(state%cdn(:, j), 2, size(pdn, 1))*spread(state%cdn(:, j), 1, size(pdn, 1))
        end do
    end subroutine

    !> Energy components of a state at the given interaction (recomputed).
    subroutine mf_energy(inter, state, kinetic, hartree, exchange_energy)
        type(mf_interaction_t), intent(in) :: inter
        type(mf_state_t), intent(in) :: state
        real(dp), intent(out) :: kinetic, hartree, exchange_energy
        real(dp), allocatable :: pup(:, :), pdn(:, :), fu(:, :), fd(:, :)
        integer :: n
        n = inter%nmodes
        allocate (pup(n, n), pdn(n, n), fu(n, n), fd(n, n))
        call densities(state, pup, pdn)
        call mf_fock(inter, pup, pdn, state%exchange, fu, fd, kinetic, hartree, &
                     exchange_energy)
    end subroutine

    !> Deterministic start orbitals per spin: eigenvectors of T + amplitude R,
    !> R a fixed pseudo-random symmetric matrix (seed 0 or amplitude 0: plane
    !> waves in kinetic order). Nonzero amplitudes break translation symmetry
    !> so that SCF can reach broken-symmetry minima.
    subroutine mf_start_orbitals(inter, seed, amplitude, c)
        type(mf_interaction_t), intent(in) :: inter
        integer, intent(in) :: seed
        real(dp), intent(in) :: amplitude
        real(dp), intent(out) :: c(:, :)
        real(dp), allocatable :: m(:, :), values(:)
        integer :: a, b, n, info
        n = inter%nmodes
        allocate (m(n, n), values(n))
        m = 0
        do a = 1, n
            m(a, a) = inter%kinetic(a) + 1.0e-9_dp*a*max(1.0_dp, inter%kinetic(n))
            if (seed == 0) cycle
            do b = 1, n
                m(a, b) = m(a, b) + amplitude*cos(0.917_dp*a*b + 1.71_dp*seed*(a + b) &
                                                  + 0.3_dp*seed)
            end do
        end do
        call symmetric_eigen(m, values, c, info)
        if (info /= 0) error stop 'mf_start_orbitals: eigen'
    end subroutine

    !> Self-consistent field at T = 0 with fixed (nup, N - nup) aufbau filling:
    !> Pulay DIIS (8 vectors) on the commutator F P - P F plus a constant
    !> virtual level shift (does not move stationary points). restricted
    !> forces C_up = C_dn (N even, nup = N/2); exchange = .false. gives
    !> Hartree. On exit the orbitals are canonical inside the occupied and
    !> virtual spaces, `aufbau` records max eps_occ < min eps_virt.
    subroutine mf_scf(cell, inter, nup, restricted, exchange, cup0, cdn0, state, ok, &
                      tolerance, max_iterations, shift)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_interaction_t), intent(in) :: inter
        integer, intent(in) :: nup
        logical, intent(in) :: restricted, exchange
        real(dp), intent(in) :: cup0(:, :), cdn0(:, :)
        type(mf_state_t), intent(out) :: state
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: tolerance, shift
        integer, intent(in), optional :: max_iterations
        integer, parameter :: depth = 8
        real(dp), allocatable :: pup(:, :), pdn(:, :), fu(:, :), fd(:, :), eu(:, :), &
            ed(:, :), hist_f(:, :, :, :), hist_e(:, :, :, :), bmat(:, :), coef(:), &
            gu(:, :), gd(:, :), values(:)
        real(dp) :: tol, lshift, ek, eh, ex, eold, err
        integer :: n, it, itmax, stored, slot, info, ndn
        ok = .false.
        n = inter%nmodes; ndn = cell%n - nup
        if (nup < 0 .or. ndn < 0 .or. nup > n .or. ndn > n) return
        if (restricted .and. nup /= ndn) return
        tol = 1.0e-10_dp; itmax = 400; lshift = 0.5_dp*cell%fermi
        if (present(tolerance)) tol = tolerance
        if (present(max_iterations)) itmax = max_iterations
        if (present(shift)) lshift = shift
        state%nup = nup; state%ndn = ndn; state%restricted = restricted
        state%exchange = exchange
        allocate (state%cup(n, n), state%cdn(n, n), state%eup(n), state%edn(n), &
                  state%fup(n), state%fdn(n))
        state%cup = cup0; state%cdn = cdn0
        if (restricted) state%cdn = cup0
        state%fup = 0; state%fup(1:nup) = 1; state%fdn = 0; state%fdn(1:ndn) = 1
        allocate (pup(n, n), pdn(n, n), fu(n, n), fd(n, n), eu(n, n), ed(n, n), &
                  gu(n, n), gd(n, n), values(n), hist_f(n, n, 2, depth), &
                  hist_e(n, n, 2, depth), bmat(depth + 1, depth + 1), coef(depth + 1))
        stored = 0; eold = huge(1.0_dp)
        do it = 1, itmax
            call densities(state, pup, pdn)
            call mf_fock(inter, pup, pdn, exchange, fu, fd, ek, eh, ex)
            if (restricted) then
                fu = 0.5_dp*(fu + fd); fd = fu
            end if
            eu = matmul(fu, pup) - matmul(pup, fu)
            ed = matmul(fd, pdn) - matmul(pdn, fd)
            err = sqrt(sum(eu**2) + sum(ed**2))
            state%energy = ek + eh + ex; state%kinetic = ek; state%hartree = eh
            state%exchange_energy = ex; state%iterations = it; state%commutator = err
            if (err < tol*max(1.0_dp, cell%fermi) .and. &
                abs(state%energy - eold) < 1.0e-12_dp*max(1.0_dp, abs(state%energy))) then
                state%converged = .true.; exit
            end if
            eold = state%energy
            slot = 1 + mod(it - 1, depth); stored = min(stored + 1, depth)
            hist_f(:, :, 1, slot) = fu; hist_f(:, :, 2, slot) = fd
            hist_e(:, :, 1, slot) = eu; hist_e(:, :, 2, slot) = ed
            if (stored >= 2) then
                call diis(hist_f(:, :, :, 1:stored), hist_e(:, :, :, 1:stored), fu, fd)
            end if
            gu = fu + lshift*(identity(n) - pup); gd = fd + lshift*(identity(n) - pdn)
            call symmetric_eigen(gu, values, state%cup, info)
            if (info /= 0) return
            if (restricted) then
                state%cdn = state%cup
            else
                call symmetric_eigen(gd, values, state%cdn, info)
                if (info /= 0) return
            end if
        end do
        call canonicalize(state%cup, fu, nup, state%eup)
        call canonicalize(state%cdn, fd, ndn, state%edn)
        state%aufbau = gap_ok(state%eup, nup) .and. gap_ok(state%edn, ndn)
        state%free = state%energy
        ok = state%converged .and. all(ieee_is_finite(state%eup)) .and. &
             all(ieee_is_finite(state%edn))
    contains
        subroutine diis(fh, eh_, fu_, fd_)
            real(dp), intent(in) :: fh(:, :, :, :), eh_(:, :, :, :)
            real(dp), intent(out) :: fu_(:, :), fd_(:, :)
            real(dp), allocatable :: b(:, :), vec(:, :), val(:), x(:), rhs(:)
            integer :: m, p, r, inf
            m = size(fh, 4)
            allocate (b(m + 1, m + 1), vec(m + 1, m + 1), val(m + 1), x(m + 1), &
                      rhs(m + 1))
            b = 0
            do p = 1, m
                do r = 1, m
                    b(p, r) = sum(eh_(:, :, :, p)*eh_(:, :, :, r))
                end do
                b(p, m + 1) = -1; b(m + 1, p) = -1
            end do
            rhs = 0; rhs(m + 1) = -1
            call symmetric_eigen(b, val, vec, inf)
            x = 0
            if (inf == 0) then
                do p = 1, m + 1
                    if (abs(val(p)) > 1.0e-14_dp*maxval(abs(val))) &
                        x = x + vec(:, p)*dot_product(vec(:, p), rhs)/val(p)
                end do
            end if
            if (inf /= 0 .or. .not. all(ieee_is_finite(x))) then
                x = 0; x(m) = 1
            end if
            fu_ = 0; fd_ = 0
            do p = 1, m
                fu_ = fu_ + x(p)*fh(:, :, 1, p); fd_ = fd_ + x(p)*fh(:, :, 2, p)
            end do
        end subroutine
    end subroutine

    pure function identity(n) result(e)
        integer, intent(in) :: n
        real(dp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function

    logical function gap_ok(eps, nocc) result(okay)
        real(dp), intent(in) :: eps(:)
        integer, intent(in) :: nocc
        okay = .true.
        if (nocc > 0 .and. nocc < size(eps)) okay = maxval(eps(1:nocc)) < &
            minval(eps(nocc + 1:)) + 1.0e-12_dp*max(1.0_dp, maxval(abs(eps)))
    end function

    !> Rotate the occupied (first nocc) and virtual columns of c separately to
    !> eigenvectors of the Fock matrix f; eps holds the corresponding values.
    subroutine canonicalize(c, f, nocc, eps)
        real(dp), intent(inout) :: c(:, :)
        real(dp), intent(in) :: f(:, :)
        integer, intent(in) :: nocc
        real(dp), intent(out) :: eps(:)
        call block(1, nocc)
        call block(nocc + 1, size(c, 2))
    contains
        subroutine block(lo, hi)
            integer, intent(in) :: lo, hi
            real(dp), allocatable :: sub(:, :), vec(:, :)
            integer :: info
            if (hi < lo) return
            sub = matmul(transpose(c(:, lo:hi)), matmul(f, c(:, lo:hi)))
            sub = 0.5_dp*(sub + transpose(sub))
            allocate (vec(hi - lo + 1, hi - lo + 1))
            call symmetric_eigen(sub, eps(lo:hi), vec, info)
            if (info /= 0) error stop 'canonicalize: eigen'
            c(:, lo:hi) = matmul(c(:, lo:hi), vec)
        end subroutine
    end subroutine

    !> Lowest converged T = 0 mean-field state over spin fillings and starts.
    !> UHF scans nup = ceil(N/2)..min(N, nmodes) (nup <-> ndn is symmetric);
    !> RHF uses nup = N/2 (N even). Start 0 is the plane-wave determinant,
    !> starts 1..nstarts-1 are symmetry-broken (mf_start_orbitals, amplitude
    !> E_F). nsolutions counts distinct converged energies (relative 1e-8);
    !> nconverged counts converged runs. Multiple solutions are expected for
    !> degenerate open shells and broken-symmetry minima at large r_s.
    subroutine mf_ground(cell, inter, restricted, exchange, nstarts, state, ok, &
                         nsolutions, nconverged)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_interaction_t), intent(in) :: inter
        logical, intent(in) :: restricted, exchange
        integer, intent(in) :: nstarts
        type(mf_state_t), intent(out) :: state
        logical, intent(out) :: ok
        integer, intent(out) :: nsolutions, nconverged
        type(mf_state_t) :: trial
        real(dp), allocatable :: cu(:, :), cd(:, :), found(:)
        integer :: nup, lo, hi, s, n
        logical :: run_ok
        ok = .false.; nsolutions = 0; nconverged = 0
        n = inter%nmodes
        allocate (cu(n, n), cd(n, n), found(0))
        if (restricted) then
            if (mod(cell%n, 2) /= 0) return
            lo = cell%n/2; hi = lo
        else
            lo = (cell%n + 1)/2; hi = min(cell%n, n)
        end if
        do nup = lo, hi
            do s = 0, nstarts - 1
                call mf_start_orbitals(inter, s, cell%fermi, cu)
                call mf_start_orbitals(inter, merge(0, s + 101, s == 0), cell%fermi, cd)
                call mf_scf(cell, inter, nup, restricted, exchange, cu, cd, trial, run_ok)
                if (.not. run_ok) cycle
                nconverged = nconverged + 1
                if (all(abs(found - trial%energy) > 1.0e-8_dp* &
                        max(1.0_dp, abs(trial%energy)))) found = [found, trial%energy]
                if (.not. ok) then
                    state = trial; ok = .true.
                else if (trial%energy < state%energy - 1.0e-12_dp*abs(state%energy)) then
                    state = trial
                end if
            end do
        end do
        nsolutions = size(found)
    end subroutine

    !> Restricted (paramagnetic) Mermin finite-temperature HF at Theta = T/T_F:
    !> stationary point of Omega = E[P] - T S[f] - mu N with Fermi occupations
    !> f = 1/(1 + exp((eps - mu)/T)) of the Fock eigenvalues and mu fixed by
    !> 2 sum f = N (grand-canonical orbital ensemble). Damped density mixing
    !> (alpha = 0.5) until ||P_new - P|| < tolerance. exchange = .false. gives
    !> finite-T Hartree. free = E - T S.
    subroutine mf_mermin(cell, inter, theta, exchange, state, ok, tolerance)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_interaction_t), intent(in) :: inter
        real(dp), intent(in) :: theta
        logical, intent(in) :: exchange
        type(mf_state_t), intent(out) :: state
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: tolerance
        real(dp), allocatable :: p(:, :), pnew(:, :), fu(:, :), fd(:, :), values(:)
        real(dp) :: tol, temperature, beta, ek, eh, ex, change, s, si
        integer :: n, it, j, info
        ok = .false.
        if (.not. (theta > 0)) return
        tol = 1.0e-11_dp
        if (present(tolerance)) tol = tolerance
        n = inter%nmodes
        temperature = theta*cell%fermi; beta = 1/temperature
        state%restricted = .true.; state%exchange = exchange; state%theta = theta
        allocate (state%cup(n, n), state%cdn(n, n), state%eup(n), state%edn(n), &
                  state%fup(n), state%fdn(n), p(n, n), pnew(n, n), fu(n, n), &
                  fd(n, n), values(n))
        call mf_start_orbitals(inter, 0, 0.0_dp, state%cup)
        call occupy(inter%kinetic)
        call density(p)
        change = huge(1.0_dp)
        do it = 1, 2000
            call mf_fock(inter, p, p, exchange, fu, fd, ek, eh, ex)
            call symmetric_eigen(fu, values, state%cup, info)
            if (info /= 0) return
            call occupy(values)
            call density(pnew)
            change = sqrt(sum((pnew - p)**2))
            state%iterations = it
            if (change < tol) then
                p = pnew; state%converged = .true.; exit
            end if
            p = 0.5_dp*(p + pnew)
        end do
        call mf_fock(inter, p, p, exchange, fu, fd, ek, eh, ex)
        state%eup = values; state%edn = values; state%cdn = state%cup
        state%fdn = state%fup
        state%kinetic = ek; state%hartree = eh; state%exchange_energy = ex
        state%energy = ek + eh + ex; state%commutator = change
        s = 0
        do j = 1, n
            if (state%fup(j) > 0 .and. state%fup(j) < 1) then
                call mf_entropy(state%fup(j), si); s = s + 2*si
            end if
        end do
        state%entropy = s
        call mf_free_energy(state%energy, s, temperature, state%free)
        ok = state%converged .and. ieee_is_finite(state%free)
    contains
        subroutine density(pm)
            real(dp), intent(out) :: pm(:, :)
            integer :: k
            pm = 0
            do k = 1, n
                pm = pm + state%fup(k)*spread(state%cup(:, k), 2, n)* &
                     spread(state%cup(:, k), 1, n)
            end do
        end subroutine
        subroutine occupy(eps)
            real(dp), intent(in) :: eps(:)
            real(dp) :: lo, hi, mid
            integer :: k, iter
            lo = minval(eps) - 60*temperature; hi = maxval(eps) + 60*temperature
            mid = lo
            do iter = 1, 200
                mid = 0.5_dp*(lo + hi)
                do k = 1, n
                    call jc_fermi(eps(k), mid, beta, state%fup(k))
                end do
                if (2*sum(state%fup) > cell%n) then
                    hi = mid
                else
                    lo = mid
                end if
                if (hi - lo < 1.0e-15_dp*max(1.0_dp, abs(mid))) exit
            end do
            state%mu = mid
        end subroutine
    end subroutine

    !> Energy of the N-1 determinant obtained by removing occupied orbital j of
    !> spin s (+1 up, -1 down) with all orbitals frozen (Koopmans: equals
    !> E_N - eps_j for a converged HF state).
    real(dp) function mf_frozen_energy(inter, state, j, s) result(e)
        type(mf_interaction_t), intent(in) :: inter
        type(mf_state_t), intent(in) :: state
        integer, intent(in) :: j, s
        type(mf_state_t) :: removed
        real(dp) :: ek, eh, ex
        removed = state
        if (s > 0) removed%fup(j) = 0
        if (s < 0) removed%fdn(j) = 0
        call mf_energy(inter, removed, ek, eh, ex)
        e = ek + eh + ex
    end function

    !> Amplitudes <S|Phi> = det(C_occ[S, :]) of the T = 0 Slater state in a
    !> determinant basis (spin orbital 2a-1 = mode a up, 2a = mode a down;
    !> occupied columns ordered up then down). norm = ||psi||^2 is the captured
    !> weight (1 when the basis holds every determinant of the spin filling).
    subroutine mf_slater_vector(cell, state, basis, psi, norm, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(in) :: state
        type(fock_basis_t), intent(in) :: basis
        real(dp), allocatable, intent(out) :: psi(:)
        real(dp), intent(out) :: norm
        logical, intent(out) :: ok
        integer(i64), allocatable :: masks(:)
        real(dp), allocatable :: cfull(:, :), mat(:, :)
        integer, allocatable :: occ(:), cols(:)
        integer :: m, n, d, p, k, r, i
        integer(i64) :: bytes
        ok = .false.; norm = 0
        call fock_basis_info(basis, m, n, d, p, bytes, ok)
        if (.not. ok .or. m /= cell%m .or. n /= cell%n) then
            ok = .false.; return
        end if
        allocate (masks(d), psi(d), cfull(m, n), mat(n, n), occ(n))
        call get_fock_states(basis, masks, ok)
        if (.not. ok) return
        cols = [pack([(k, k=1, size(state%fup))], state%fup > 0.5_dp)]
        cfull = 0
        do k = 1, state%nup
            cfull(1:m:2, k) = state%cup(:, cols(k))
        end do
        cols = [pack([(k, k=1, size(state%fdn))], state%fdn > 0.5_dp)]
        do k = 1, state%ndn
            cfull(2:m:2, state%nup + k) = state%cdn(:, cols(k))
        end do
        do k = 1, d
            r = 0
            do i = 1, m
                if (btest(masks(k), i - 1)) then
                    r = r + 1; occ(r) = i
                end if
            end do
            do r = 1, n
                mat(r, :) = cfull(occ(r), :)
            end do
            psi(k) = determinant(mat)
        end do
        norm = sum(psi**2)
        ok = all(ieee_is_finite(psi))
    end subroutine

    real(dp) function determinant(a) result(det)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: w(size(a, 1), size(a, 1)), row(size(a, 1))
        integer :: n, k, piv, i
        n = size(a, 1); w = a; det = 1
        do k = 1, n
            piv = k - 1 + maxloc(abs(w(k:n, k)), 1)
            if (.not. abs(w(piv, k)) > 0) then
                det = 0; return
            end if
            if (piv /= k) then
                row = w(k, :); w(k, :) = w(piv, :); w(piv, :) = row; det = -det
            end if
            det = det*w(k, k)
            do i = k + 1, n
                w(i, k:n) = w(i, k:n) - w(i, k)/w(k, k)*w(k, k:n)
            end do
        end do
    end function

    !> Smallest exact symmetry space holding the T = 0 Slater state: the
    !> (S_z, K) sectors of its spin filling in which it has weight above
    !> 1e-28. Returns that sector union, the amplitudes and the K count.
    subroutine mf_slater_sector(cell, state, sector, psi, norm, nk, ok)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(in) :: state
        type(jellium_sector_t), intent(out) :: sector
        real(dp), allocatable, intent(out) :: psi(:)
        real(dp), intent(out) :: norm
        integer, intent(out) :: nk
        logical, intent(out) :: ok
        type(jellium_sector_t) :: full
        integer, allocatable :: ks(:, :), support(:, :)
        integer(i64), allocatable :: masks(:)
        logical, allocatable :: used(:)
        integer :: k, l, i, kk(3)
        call jellium_sector_momenta(cell, state%nup, ks)
        call jellium_sector(cell, state%nup, ks, full, ok)
        if (.not. ok) return
        call mf_slater_vector(cell, state, full%basis, psi, norm, ok)
        if (.not. ok) return
        allocate (masks(full%d), used(size(ks, 2)))
        call get_fock_states(full%basis, masks, ok)
        used = .false.
        do k = 1, full%d
            if (.not. psi(k)**2 > 1.0e-28_dp) cycle
            kk = 0
            do i = 1, cell%m
                if (btest(masks(k), i - 1)) kk = kk + cell%momentum(:, i)
            end do
            do l = 1, size(ks, 2)
                if (all(ks(:, l) == kk)) used(l) = .true.
            end do
        end do
        nk = count(used)
        support = ks(:, pack([(l, l=1, size(ks, 2))], used))
        call jellium_sector(cell, state%nup, support, sector, ok)
        if (.not. ok) return
        call mf_slater_vector(cell, state, sector%basis, psi, norm, ok)
    end subroutine

    !> One-state spectrum record (weight 1) so the Q4 RDM, pair-bin and
    !> cumulant queries apply unchanged to a Slater state.
    subroutine mf_as_spectrum(sector, psi, energy, spectrum)
        type(jellium_sector_t), intent(in) :: sector
        real(dp), intent(in) :: psi(:), energy
        type(jellium_spectrum_t), intent(out) :: spectrum
        spectrum%count = 1; spectrum%nsectors = 1
        allocate (spectrum%sectors(1))
        spectrum%sectors(1) = sector
        spectrum%kept = [1]; spectrum%top = [energy]; spectrum%sector_of = [1]
        spectrum%offset = [0_i64]; spectrum%energy = [energy]
        spectrum%residual = [0.0_dp]; spectrum%flat = psi/sqrt(sum(psi**2))
    end subroutine

    !> Full-parent residual estimator of the Slater state (see mf_estimate_t)
    !> with the matrix-free Fock action on its symmetry space. e0/e1 come from
    !> FortNum block Lanczos (residual-checked, ordering candidate); Temple is
    !> evaluated in outward interval arithmetic (lower endpoint).
    subroutine mf_residual(model, sector, psi, estimate, ok)
        type(fock_hamiltonian_t), intent(in) :: model
        type(jellium_sector_t), intent(in) :: sector
        real(dp), intent(in) :: psi(:)
        type(mf_estimate_t), intent(out) :: estimate
        logical, intent(out) :: ok
        real(dp), allocatable :: x(:), hx(:), zero(:), scratch(:), vec(:, :)
        real(dp) :: values(2), residuals(2)
        type(interval_t) :: bound
        integer :: d, nev, iterations
        d = sector%d
        allocate (x(d), hx(d), zero(d), scratch(d))
        zero = 0
        estimate%norm = sqrt(sum(psi**2)); x = psi/estimate%norm
        call fock_matvec(model, sector%basis, x, zero, hx, scratch, ok)
        if (.not. ok) return
        estimate%rho = dot_product(x, hx)
        estimate%sigma = sqrt(sum((hx - estimate%rho*x)**2))
        estimate%dimension = d
        nev = min(2, d)
        allocate (vec(d, nev))
        call jellium_lowest(model, sector%basis, nev, values(1:nev), vec, &
                            residuals(1:nev), ok, iterations)
        if (.not. ok) return
        estimate%e0 = values(1); estimate%r0 = residuals(1)
        estimate%temple_valid = .false.
        if (nev < 2) return
        estimate%e1 = values(2); estimate%r1 = residuals(2)
        estimate%beta = values(2) - residuals(2)
        if (estimate%rho < estimate%beta) then
            call mf_temple_interval(ipoint(estimate%rho), ipoint(estimate%sigma), &
                                    ipoint(estimate%beta), bound)
            estimate%temple = bound%lo
            estimate%temple_valid = ieee_is_finite(bound%lo)
        end if
    end subroutine

    !> Static TDHF (coupled-perturbed HF, RPA with exchange) response per volume
    !> of a converged T = 0 state to lambda*A, A = integral cos(q.x) n(x):
    !> chi = -(4/V) a^T (A+B)^(-1) a with a_ai = <a|A|i> and the real orbital
    !> Hessian (A+B) x = (eps_a - eps_i) x_ai + [G(dP)]_ai, dP = sum x (|a><i| +
    !> |i><a|), G the linear two-body mean field of mf_fock. exchange = .false.
    !> in the state (Hartree reference) drops exchange in eps and kernel: RPA.
    !> hessian_min < 0 flags an unstable (saddle) mean-field solution; chi is
    !> then refused if a has weight on a nonpositive mode.
    subroutine mf_static_response(cell, inter, state, q, chi, hessian_min, ok, nexc)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_interaction_t), intent(in) :: inter
        type(mf_state_t), intent(in) :: state
        integer, intent(in) :: q(3)
        real(dp), intent(out) :: chi, hessian_min
        logical, intent(out) :: ok
        integer, intent(out) :: nexc
        real(dp), allocatable :: amode(:, :), hess(:, :), avec(:), values(:), &
            vec(:, :), du(:, :), dd(:, :), gh(:, :), xu(:, :), xd(:, :)
        integer, allocatable :: es(:), ea(:), ei(:)
        real(dp) :: overlap, proj, scale
        integer :: n, e, f, k, info
        ok = .false.; chi = 0; hessian_min = 0; nexc = 0
        n = inter%nmodes
        call density_operator(cell, state%nup, q, amode, ok)
        if (.not. ok) return
        ok = .false.
        call excitations(state%fup, 1); call excitations(state%fdn, -1)
        nexc = size(es)
        if (nexc == 0) return
        allocate (hess(nexc, nexc), avec(nexc), values(nexc), vec(nexc, nexc), &
                  du(n, n), dd(n, n), gh(n, n), xu(n, n), xd(n, n))
        do e = 1, nexc
            avec(e) = dot_product(orb(es(e), ea(e)), matmul(amode, orb(es(e), ei(e))))
            du = 0; dd = 0
            if (es(e) > 0) du = outer(ea(e), ei(e), 1)
            if (es(e) < 0) dd = outer(ea(e), ei(e), -1)
            call two_body(inter, du, dd, state%exchange, gh, xu, xd)
            xu = gh + xu; xd = gh + xd
            do f = 1, nexc
                if (es(f) > 0) hess(f, e) = dot_product(orb(1, ea(f)), &
                    matmul(xu, orb(1, ei(f))))
                if (es(f) < 0) hess(f, e) = dot_product(orb(-1, ea(f)), &
                    matmul(xd, orb(-1, ei(f))))
            end do
            hess(e, e) = hess(e, e) + eps(es(e), ea(e)) - eps(es(e), ei(e))
        end do
        hess = 0.5_dp*(hess + transpose(hess))
        call symmetric_eigen(hess, values, vec, info)
        if (info /= 0) return
        hessian_min = values(1)
        scale = max(1.0e-300_dp, maxval(abs(values)))
        overlap = 0
        do k = 1, nexc
            proj = dot_product(vec(:, k), avec)
            if (values(k) <= 1.0e-12_dp*scale) then
                if (abs(proj) > 1.0e-10_dp*max(1.0e-300_dp, sqrt(sum(avec**2)))) return
                cycle
            end if
            overlap = overlap + proj**2/values(k)
        end do
        call jc_susceptibility(overlap, cell%volume, chi)
        ok = ieee_is_finite(chi)
    contains
        subroutine excitations(occupation, spin)
            real(dp), intent(in) :: occupation(:)
            integer, intent(in) :: spin
            integer :: i, a
            if (.not. allocated(es)) allocate (es(0), ea(0), ei(0))
            do i = 1, n
                if (occupation(i) < 0.5_dp) cycle
                do a = 1, n
                    if (occupation(a) > 0.5_dp) cycle
                    es = [es, spin]; ea = [ea, a]; ei = [ei, i]
                end do
            end do
        end subroutine
        function orb(spin, j) result(c)
            integer, intent(in) :: spin, j
            real(dp) :: c(n)
            if (spin > 0) then
                c = state%cup(:, j)
            else
                c = state%cdn(:, j)
            end if
        end function
        real(dp) function eps(spin, j)
            integer, intent(in) :: spin, j
            eps = merge(state%eup(j), state%edn(j), spin > 0)
        end function
        function outer(a, i, spin) result(m)
            integer, intent(in) :: a, i, spin
            real(dp) :: m(n, n), ca(n), ci(n)
            ca = orb(spin, a); ci = orb(spin, i)
            m = spread(ca, 2, n)*spread(ci, 1, n) + spread(ci, 2, n)*spread(ca, 1, n)
        end function
    end subroutine

    !> Mode-space matrix of A = integral cos(q.x) n(x) (same for both spins)
    !> from the engine's one-body density operator of the cell.
    subroutine density_operator(cell, nup, q, amode, ok)
        type(jellium_cell_t), intent(in) :: cell
        integer, intent(in) :: nup, q(3)
        real(dp), allocatable, intent(out) :: amode(:, :)
        logical, intent(out) :: ok
        type(jellium_sector_t) :: sector
        type(fock_hamiltonian_t) :: model
        integer, allocatable :: ks(:, :)
        real(dp), allocatable :: h(:, :), v(:, :)
        integer :: p
        call jellium_sector_momenta(cell, nup, ks)
        call jellium_sector(cell, nup, ks(:, 1:1), sector, ok)
        if (.not. ok) return
        call jellium_density_model(cell, sector%basis, q, model, ok)
        if (.not. ok) return
        p = cell%m*(cell%m - 1)/2
        allocate (h(cell%m, cell%m), v(p, p))
        call get_fock_integrals(model, h, v, ok)
        if (.not. ok) return
        amode = h(1:cell%m:2, 1:cell%m:2)
        ok = .not. any(abs(amode - h(2:cell%m:2, 2:cell%m:2)) > 0)
    end subroutine

    !> Hellmann-Feynman coupling derivative of a converged (stationary)
    !> mean-field state: dE/dg = <Phi|V1|Phi> = E_int at unit coupling with
    !> frozen orbitals; orbital relaxation does not enter at stationarity.
    real(dp) function mf_coupling_derivative(cell, state) result(de)
        type(jellium_cell_t), intent(in) :: cell
        type(mf_state_t), intent(in) :: state
        type(mf_interaction_t) :: unit
        real(dp) :: ek, eh, ex
        call mf_interaction(cell, unit, 1.0_dp)
        call mf_energy(unit, state, ek, eh, ex)
        de = eh + ex
    end function
end module kin6d_mean_field_cell
