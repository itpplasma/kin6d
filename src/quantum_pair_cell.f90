!> H1b exact N = 2 quantum reference in the Q4 periodic jellium cell.
!>
!> Two electrons on the torus separate exactly into the centre of mass
!> (mass 2m, momenta K = 2 pi n_K/L) and the relative coordinate r = x1 - x2
!> (mass m/2) in the Ewald pair potential (Fourier coefficients 4 pi g/(V Q^2),
!> Q = 0 removed, the Q4 jellium coefficient). A lattice translation of one
!> particle shifts R by L/2 and r by L, so the relative wave function is
!> periodic when the component of n_K is even and antiperiodic when it is
!> odd: relative momenta k = pi kv/L with kv_i == n_K,i (mod 2). The relative
!> problem depends only on the twist class t = number of odd components
!> (cubic symmetry), and the centre-of-mass sums factorize per component.
!> Inversion r -> -r splits each twist block into even (spin singlet) and
!> odd (spin triplet, three spin states) spatial parity. All eigenstates of
!> the truncated relative basis |kv|^2 <= kmax2 are kept (dense LAPACK), so
!> the only discretization is that basis cutoff (indicator: a second cutoff).
!>
!> Observables use the common estimator contract of the classical cell
!> (docs/quantum-classical-correlations.md): spin-summed shell-averaged g,
!> spin channels in the Q4 normalization 2 (N/2)^2, the pair Fourier
!> moments f2(Q) on S = {1 <= |n|^2 <= 3}, and dg/dTheta by the exact
!> thermal covariance in the eigenbasis.
module kin6d_quantum_pair_cell
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_species, only: units_t, species_t, plasma_scales_t, atomic_units, &
        electron, plasma_scales, pair_coupling
    use kin6d_quantum_fock_leaves, only: qf_kinetic, qf_coulomb
    use kin6d_jellium_cell_leaves, only: jc_shell_moment, jc_shell_moment_series, &
        jc_shell_average
    implicit none
    private
    public :: quantum_pair_cell_t, init_quantum_pair_cell, quantum_pair_thermal
    public :: quantum_pair_relative_spectrum, quantum_pair_cutoff
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
        subroutine dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: ta, tb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
    end interface

    type :: pair_block_t
        integer :: d = 0
        integer, allocatable :: kv(:, :)
        real(dp), allocatable :: energy(:), vec(:, :)
    end type

    type :: quantum_pair_cell_t
        real(dp) :: length = 0, volume = 0, hbar = 1, mass = 1, coupling = 1
        real(dp) :: rs = 0, fermi = 0
        integer :: kmax2 = 0, kmax = 0
        !> block(t, 1) even, block(t, 2) odd parity, twist (1,..,1,0,..,0).
        type(pair_block_t) :: block(0:3, 2)
    end type
contains
    !> Smallest relative cutoff (units (pi/L)^2) whose kinetic energy is at
    !> least tau kT at Theta = theta_max, and at least floor2.
    integer function quantum_pair_cutoff(rs, theta_max, tau, floor2) result(kmax2)
        real(dp), intent(in) :: rs, theta_max, tau
        integer, intent(in) :: floor2
        type(plasma_scales_t) :: unit_cell, scales
        type(units_t) :: units
        type(species_t) :: e
        real(dp) :: length, kt, e1
        units = atomic_units(); e = electron()
        unit_cell = plasma_scales(atomic_units(), electron(), 2.0_dp, 1.0_dp)
        length = rs/unit_cell%rs
        scales = plasma_scales(atomic_units(), electron(), 2/length**3, 1.0_dp)
        kt = theta_max*scales%fermi
        ! Relative kinetic energy of kv: qf_kinetic(L, hbar, m/2, |kv|^2/4).
        call qf_kinetic(length, units%hbar, e%mass/2, 0.25_dp, e1)
        kmax2 = max(floor2, ceiling(tau*kt/e1))
    end function

    subroutine init_quantum_pair_cell(cell, rs, kmax2, ok, coupling)
        type(quantum_pair_cell_t), intent(out) :: cell
        real(dp), intent(in) :: rs
        integer, intent(in) :: kmax2
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: coupling
        type(units_t) :: units
        type(species_t) :: e
        type(plasma_scales_t) :: unit_cell, scales
        integer, allocatable :: reps(:, :)
        integer :: t, par
        ok = .false.
        if (.not. (rs > 0) .or. kmax2 < 3) return
        units = atomic_units(); e = electron()
        unit_cell = plasma_scales(units, e, 2.0_dp, 1.0_dp)
        cell%length = rs/unit_cell%rs; cell%volume = cell%length**3
        scales = plasma_scales(units, e, 2/cell%volume, 1.0_dp)
        cell%rs = scales%rs; cell%fermi = scales%fermi
        cell%hbar = units%hbar; cell%mass = e%mass
        cell%coupling = pair_coupling(units, e, e)
        if (present(coupling)) cell%coupling = coupling
        cell%kmax2 = kmax2; cell%kmax = int(sqrt(real(kmax2, dp)))
        do t = 0, 3
            call representatives(cell%kmax, kmax2, t, reps)
            do par = 1, 2
                call solve_block(cell, reps, par, cell%block(t, par), ok)
                if (.not. ok) return
            end do
        end do
    end subroutine

    !> Relative momenta kv (twist t: first t components odd) with kv = 0 or
    !> lexicographically positive, |kv|^2 <= kmax2.
    subroutine representatives(kmax, kmax2, t, reps)
        integer, intent(in) :: kmax, kmax2, t
        integer, allocatable, intent(out) :: reps(:, :)
        integer :: a, b, c, m, v(3), pass
        do pass = 1, 2
            m = 0
            do a = -kmax, kmax
                do b = -kmax, kmax
                    do c = -kmax, kmax
                        v = [a, b, c]
                        if (sum(v**2) > kmax2) cycle
                        if (modulo(a, 2) /= merge(1, 0, t >= 1)) cycle
                        if (modulo(b, 2) /= merge(1, 0, t >= 2)) cycle
                        if (modulo(c, 2) /= merge(1, 0, t >= 3)) cycle
                        if (.not. positive(v)) cycle
                        m = m + 1
                        if (pass == 2) reps(:, m) = v
                    end do
                end do
            end do
            if (pass == 1) allocate (reps(3, m))
        end do
    end subroutine

    pure logical function positive(v)
        integer, intent(in) :: v(3)
        positive = all(v == 0)
        if (v(1) /= 0) then
            positive = v(1) > 0
        else if (v(2) /= 0) then
            positive = v(2) > 0
        else if (v(3) /= 0) then
            positive = v(3) > 0
        end if
    end function

    !> Ewald pair Fourier coefficient for a relative momentum difference
    !> delta (units pi/L, even components); zero for delta = 0.
    real(dp) function pair_v(cell, delta) result(v)
        type(quantum_pair_cell_t), intent(in) :: cell
        integer, intent(in) :: delta(3)
        v = 0
        if (all(delta == 0)) return
        call qf_coulomb(cell%length, cell%coupling, real(sum(delta**2), dp)/4, v)
    end function

    !> Parity-adapted block: even e_k = (|k> + |-k>)/sqrt2 (|0> alone),
    !> odd o_k = (|k> - |-k>)/sqrt2; <e|V|e'> = v(k-k') + v(k+k').
    subroutine solve_block(cell, reps, par, block, ok)
        type(quantum_pair_cell_t), intent(in) :: cell
        integer, intent(in) :: reps(:, :), par
        type(pair_block_t), intent(out) :: block
        logical, intent(out) :: ok
        real(dp), allocatable :: h(:, :), work(:)
        real(dp) :: sgn, fa, fb, kin
        integer :: i, j, d, info
        logical :: zero
        ok = .false.
        d = 0
        allocate (block%kv(3, size(reps, 2)))
        do i = 1, size(reps, 2)
            zero = all(reps(:, i) == 0)
            if (par == 2 .and. zero) cycle
            d = d + 1; block%kv(:, d) = reps(:, i)
        end do
        block%kv = block%kv(:, 1:d)
        block%d = d
        allocate (block%energy(d))
        if (d == 0) then
            allocate (block%vec(0, 0)); ok = .true.; return
        end if
        sgn = merge(1.0_dp, -1.0_dp, par == 1)
        allocate (h(d, d))
        do j = 1, d
            fb = merge(sqrt(2.0_dp), 1.0_dp, all(block%kv(:, j) == 0))
            do i = 1, j
                fa = merge(sqrt(2.0_dp), 1.0_dp, all(block%kv(:, i) == 0))
                h(i, j) = (pair_v(cell, block%kv(:, i) - block%kv(:, j)) + &
                           sgn*pair_v(cell, block%kv(:, i) + block%kv(:, j)))
                ! The zero representative is a single plane wave.
                if (fa > 1 .and. fb > 1) then
                    h(i, j) = 0
                else if (fa > 1 .or. fb > 1) then
                    h(i, j) = h(i, j)/sqrt(2.0_dp)
                end if
                h(j, i) = h(i, j)
            end do
            call qf_kinetic(cell%length, cell%hbar, cell%mass/2, &
                            real(sum(block%kv(:, j)**2), dp)/4, kin)
            h(j, j) = h(j, j) + kin
        end do
        allocate (work(max(1, 3*d)*8))
        call dsyev('V', 'U', d, h, d, block%energy, work, size(work), info)
        if (info /= 0) return
        call move_alloc(h, block%vec)
        ok = all(ieee_is_finite(block%energy))
    end subroutine

    !> All eigenvalues of the relative Hamiltonian in an explicit plane-wave
    !> list kv (units pi/L, one twist parity), no symmetry reduction; the
    !> oracle hook for the Q4 two-particle sector equivalence.
    subroutine quantum_pair_relative_spectrum(cell, kv, energies, ok)
        type(quantum_pair_cell_t), intent(in) :: cell
        integer, intent(in) :: kv(:, :)
        real(dp), intent(out) :: energies(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: h(:, :), work(:)
        integer :: i, j, d, info
        d = size(kv, 2); ok = .false.
        if (size(energies) /= d) return
        allocate (h(d, d), work(8*max(1, 3*d)))
        do j = 1, d
            do i = 1, d
                h(i, j) = pair_v(cell, kv(:, i) - kv(:, j))
            end do
            call qf_kinetic(cell%length, cell%hbar, cell%mass/2, &
                            real(sum(kv(:, j)**2), dp)/4, work(1))
            h(j, j) = h(j, j) + work(1)
        end do
        call dsyev('N', 'U', d, h, d, energies, work, size(work), info)
        ok = info == 0
    end subroutine

    !> Canonical observables at Theta = kT/E_F > 0 (all spin states, all
    !> centre-of-mass momenta), or for theta <= 0 the equal mixture of the
    !> degenerate ground multiplet. g_sum is the spin-summed shell average
    !> V p(r); g_same/g_opp use the Q4 normalization 2 (N/2)^2 = 2; f2(Q) on
    !> the common set S (cell_q_set order). dg_dtheta: exact thermal
    !> covariance d<O>/dbeta = -(<O E> - <O><E>) times dbeta/dTheta.
    subroutine quantum_pair_thermal(cell, theta, edges, q, g_sum, g_same, g_opp, &
                                    f2, ok, dg_dtheta, energy)
        type(quantum_pair_cell_t), intent(in) :: cell
        real(dp), intent(in) :: theta, edges(:)
        integer, intent(in) :: q(:, :)
        real(dp), intent(out) :: g_sum(:), g_same(:), g_opp(:), f2(:)
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: dg_dtheta(:), energy
        real(dp), allocatable :: s(:, :, :, :), b1(:, :), b2(:, :), dm(:, :), &
            dme(:, :), w(:), we(:), shell(:, :), sq(:, :)
        integer, allocatable :: keep(:)
        real(dp) :: beta, zc(0:1), ec(0:1), eps, wcm(0:3), ecm(0:3), emin, wmax, z, &
            emean, ca, cb, sa, sb, qn, fa, fb, kern
        integer :: k, t, par, n, j, i, a, b, nk, sig, tau, d(3), m, mmax, nb
        logical :: ground

        ok = .false.
        k = cell%kmax; nb = size(edges) - 1
        ground = .not. (theta > 0)
        beta = 0
        if (.not. ground) beta = 1/(theta*cell%fermi)
        ! Centre-of-mass sums per Cartesian component and parity of n_K.
        zc = 0; ec = 0
        n = 0
        do
            call qf_kinetic(cell%length, cell%hbar, 2*cell%mass, real(n*n, dp), eps)
            if (ground) then
                if (n <= 1) zc(n) = zc(n) + merge(1.0_dp, 2.0_dp, n == 0)
                if (n <= 1) ec(n) = eps
                if (n >= 1) exit
            else
                if (n > 0 .and. beta*eps > 80) exit
                zc(mod(n, 2)) = zc(mod(n, 2)) + merge(1, 2, n == 0)*exp(-beta*eps)
                ec(mod(n, 2)) = ec(mod(n, 2)) + merge(1, 2, n == 0)*eps*exp(-beta*eps)
            end if
            n = n + 1
        end do
        if (.not. ground) ec = ec/zc
        do t = 0, 3
            wcm(t) = choose3(t)*zc(1)**t*zc(0)**(3 - t)
            ecm(t) = t*ec(1) + (3 - t)*ec(0)
        end do
        emin = huge(1.0_dp)
        do t = 0, 3
            do par = 1, 2
                if (cell%block(t, par)%d == 0) cycle
                if (ground) then
                    emin = min(emin, cell%block(t, par)%energy(1) + ecm(t))
                else
                    emin = min(emin, cell%block(t, par)%energy(1))
                end if
            end do
        end do
        wmax = 0
        do t = 0, 3
            do par = 1, 2
                if (cell%block(t, par)%d == 0) cycle
                wmax = max(wmax, weight(t, cell%block(t, par)%energy(1)))
            end do
        end do
        allocate (s(-k:k, -k:k, -k:k, 4))
        s = 0
        do t = 0, 3
            do par = 1, 2
                associate (bl => cell%block(t, par))
                    if (bl%d == 0) cycle
                    allocate (w(bl%d), we(bl%d))
                    do j = 1, bl%d
                        w(j) = weight(t, bl%energy(j))
                        we(j) = w(j)*(bl%energy(j) + ecm(t))
                    end do
                    keep = pack([(j, j=1, bl%d)], w > 1.0e-17_dp*wmax)
                    nk = size(keep)
                    if (nk == 0) then
                        deallocate (w, we); cycle
                    end if
                    allocate (b1(bl%d, nk), b2(bl%d, nk), dm(bl%d, bl%d), &
                              dme(bl%d, bl%d))
                    do j = 1, nk
                        b1(:, j) = bl%vec(:, keep(j))*w(keep(j))
                        b2(:, j) = bl%vec(:, keep(j))*we(keep(j))
                    end do
                    dm = 0; dme = 0
                    call dgemm('N', 'T', bl%d, bl%d, nk, 1.0_dp, b1, bl%d, &
                               bl%vec(:, keep), bl%d, 0.0_dp, dm, bl%d)
                    call dgemm('N', 'T', bl%d, bl%d, nk, 1.0_dp, b2, bl%d, &
                               bl%vec(:, keep), bl%d, 0.0_dp, dme, bl%d)
                    ! Expand the parity basis to plane waves: coefficient
                    ! c/sqrt2 at +k and +-c/sqrt2 at -k (zero: c at 0).
                    do b = 1, bl%d
                        cb = merge(1.0_dp, 1/sqrt(2.0_dp), all(bl%kv(:, b) == 0))
                        do a = 1, bl%d
                            ca = merge(1.0_dp, 1/sqrt(2.0_dp), all(bl%kv(:, a) == 0))
                            do sig = 1, -1, -2
                                if (sig < 0 .and. all(bl%kv(:, a) == 0)) cycle
                                sa = merge(1.0_dp, merge(1.0_dp, -1.0_dp, par == 1), &
                                           sig > 0)
                                do tau = 1, -1, -2
                                    if (tau < 0 .and. all(bl%kv(:, b) == 0)) cycle
                                    sb = merge(1.0_dp, merge(1.0_dp, -1.0_dp, par == 1), &
                                               tau > 0)
                                    d = (sig*bl%kv(:, a) - tau*bl%kv(:, b))/2
                                    s(d(1), d(2), d(3), par) = s(d(1), d(2), d(3), par) + &
                                        ca*cb*sa*sb*dm(a, b)
                                    s(d(1), d(2), d(3), par + 2) = &
                                        s(d(1), d(2), d(3), par + 2) + ca*cb*sa*sb*dme(a, b)
                                end do
                            end do
                        end do
                    end do
                    deallocate (w, we, b1, b2, dm, dme)
                end associate
            end do
        end do
        ! Spin multiplicities: even = singlet (1), odd = triplet (3).
        z = s(0, 0, 0, 1) + 3*s(0, 0, 0, 2)
        emean = (s(0, 0, 0, 3) + 3*s(0, 0, 0, 4))/z
        if (.not. (z > 0)) return
        ! Shell sums by |n|^2, then the radial shell kernel per bin.
        mmax = 3*k*k
        allocate (sq(0:mmax, 4), shell(0:mmax, nb))
        sq = 0
        do a = -k, k
            do b = -k, k
                do j = -k, k
                    m = a*a + b*b + j*j
                    sq(m, :) = sq(m, :) + s(a, b, j, :)
                end do
            end do
        end do
        do m = 0, mmax
            qn = 2*pi*sqrt(real(m, dp))/cell%length
            do i = 1, nb
                if (m == 0) then
                    shell(m, i) = 1
                else
                    call moment(qn*edges(i), fa); call moment(qn*edges(i + 1), fb)
                    call jc_shell_average(fa, fb, qn, edges(i), edges(i + 1), kern)
                    shell(m, i) = kern
                end if
            end do
        end do
        do i = 1, nb
            g_sum(i) = sum((sq(:, 1) + 3*sq(:, 2))*shell(:, i))/z
            g_same(i) = 2*sum(sq(:, 2)*shell(:, i))/z
            g_opp(i) = sum((sq(:, 1) + sq(:, 2))*shell(:, i))/z
            if (present(dg_dtheta)) then
                if (ground) then
                    dg_dtheta(i) = 0
                else
                    dg_dtheta(i) = -(sum((sq(:, 3) + 3*sq(:, 4))*shell(:, i))/z - &
                                     g_sum(i)*emean)*(-beta/theta)
                end if
            end if
        end do
        do j = 1, size(q, 2)
            f2(j) = cubic_average(s(:, :, :, 1) + 3*s(:, :, :, 2), k, q(:, j))/z
        end do
        if (present(energy)) energy = emean
        ok = all(ieee_is_finite(g_sum)) .and. all(ieee_is_finite(f2))
    contains
        real(dp) function weight(tt, e) result(wt)
            integer, intent(in) :: tt
            real(dp), intent(in) :: e
            if (ground) then
                wt = 0
                if (e + ecm(tt) <= emin + 1.0e-9_dp*max(1.0_dp, abs(emin))) &
                    wt = choose3(tt)*2.0_dp**tt
            else
                wt = wcm(tt)*exp(-beta*(e - emin))
            end if
        end function
        subroutine moment(x, f)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: f
            if (x < 0.1_dp) then
                call jc_shell_moment_series(x, f)
            else
                call jc_shell_moment(x, f)
            end if
        end subroutine
    end subroutine

    pure real(dp) function choose3(t) result(c)
        integer, intent(in) :: t
        c = merge(1.0_dp, 3.0_dp, t == 0 .or. t == 3)
    end function

    !> Average of f(R n) over the 48 signed permutations R of the cube group.
    real(dp) function cubic_average(f, k, n) result(avg)
        real(dp), intent(in) :: f(-k:, -k:, -k:)
        integer, intent(in) :: k, n(3)
        integer, parameter :: perm(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, &
                                                    2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        integer :: p, sx, sy, sz, v(3)
        avg = 0
        do p = 1, 6
            do sx = -1, 1, 2
                do sy = -1, 1, 2
                    do sz = -1, 1, 2
                        v = [sx*n(perm(1, p)), sy*n(perm(2, p)), sz*n(perm(3, p))]
                        if (maxval(abs(v)) <= k) avg = avg + f(v(1), v(2), v(3))
                    end do
                end do
            end do
        end do
        avg = avg/48
    end function
end module kin6d_quantum_pair_cell
