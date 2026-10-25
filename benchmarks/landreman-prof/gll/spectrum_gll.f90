program spectrum_gll
    ! Measurement-only probe (w3-landreman-profile): dense spectrum of the GL1
    ! damped unreduced map operator A = P H P + mu P + Q and of candidate
    ! preconditioned operators, plus dense PCG iteration counts at tol 1e-8.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_gll_spike
    use fortfem_tensor_operator, only: apply_tensor_operator
    use kin6d_landreman_benchmark, only: landreman_frozen_t
    implicit none
    real(dp), parameter :: delta = 1.0_dp/64
    type(landreman_frozen_t) :: field
    type(gll_map_t) :: map
    real(dp), allocatable :: x(:), base(:), g(:), e(:), y(:), a(:,:), k(:,:), m(:,:), &
        s(:,:), w(:), work(:), rhs(:), af(:,:), ev(:,:)
    integer, allocatable :: free(:), comp(:)
    real(dp) :: mu, r, z, sv, energy
    integer :: p, n, nd, nf, i, j, status, info, lw, nn, c, it
    character(32) :: arg
    p = 6; n = 1; mu = 1.e-2_dp
    call get_command_argument(1, arg); if (len_trim(arg) > 0) read(arg, *) p
    call get_command_argument(2, arg); if (len_trim(arg) > 0) read(arg, *) n
    call get_command_argument(3, arg); if (len_trim(arg) > 0) read(arg, *) mu
    field%delta = delta
    call build_gll_map(p, n, delta, field, map, x, status)
    if (status /= 0) error stop 'setup'
    nd = size(x); nn = map%nodes; base = x
    do i = 1, nn  ! same smooth perturbation as benchmark_gll.f90:120-127
        r = base(i); z = base(nn + i)
        sv = ((r*r - 1)**2/4 + z*z)/delta
        x(i) = x(i) + 1.e-4_dp*(1 - sv)*cos(17*z)
        x(nn + i) = x(nn + i) + 1.e-4_dp*(1 - sv)*sin(11*r)
        x(2*nn + i) = 1.e-3_dp*(1 - sv)*sin(13*r + 9*z)
    end do
    where (map%fixed) x = base
    allocate(g(nd), e(nd), y(nd))
    call evaluate_gll_map(map, x, energy, g, status, .true.)
    if (status /= 0) error stop 'linearize'
    free = pack([(i, i=1, nd)], .not. map%fixed); nf = size(free)
    comp = (free - 1)/nn + 1
    allocate(a(nf, nf), k(nf, nf), m(nf, nf))
    do j = 1, nf  ! dense free-free block of H via HVP columns
        e = 0; e(free(j)) = 1
        call apply_gll_hessian(map, e, y)
        a(:, j) = y(free)
    end do
    a = (a + transpose(a))/2
    ! Reference-metric scalar Laplacian per component (FortFEM diffusion), the
    ! operator approximated by GL1's scalar p-MG.
    k = 0
    do j = 1, nf
        c = comp(j)
        e = 0; e(free(j)) = 1
        call apply_tensor_operator(map%reference, e((c - 1)*nn + 1:c*nn), y(1:nn), .false.)
        do i = 1, nf
            if (comp(i) == c) k(i, j) = y(free(i) - (c - 1)*nn)
        end do
    end do
    k = (k + transpose(k))/2
    print '(a,3(i0,1x),es10.2)', 'p n free mu=', p, n, nf, mu
    call spectrum('H (undamped)', a, 0.0_dp)
    call spectrum('H+mu I', a, mu)
    rhs = -g(free)/norm2(g(free))
    ! Preconditioners: M symmetric positive definite, eigenvalues of M^-1 A.
    m = 0
    do i = 1, nf  ! (1) Jacobi
        m(i, i) = a(i, i) + mu
    end do
    call pspectrum('Jacobi', a, mu, m)
    m = 0  ! (2) exact field-split block diagonal: [H_RR H_RZ; H_ZR H_ZZ], H_ww
    do j = 1, nf
        do i = 1, nf
            if ((comp(i) == 3) .eqv. (comp(j) == 3)) m(i, j) = a(i, j)
        end do
        m(j, j) = m(j, j) + mu
    end do
    call pspectrum('fieldsplit (RZ | w) exact blocks', a, mu, m)
    m = k  ! (3) scalar reference Laplacian per component + mu (GL1 p-MG target)
    do i = 1, nf
        m(i, i) = m(i, i) + mu
    end do
    call pspectrum('scalar Laplacian+mu (exact MG limit)', a, mu, m)
    ! (4) undamped exact H (limit of a perfect Hessian-Galerkin p-MG at mu->0)
    m = a
    do i = 1, nf
        m(i, i) = m(i, i) + 1.e-2_dp*mu
    end do
    call pspectrum('H+mu/100 (shift-mismatch control)', a, mu, m)
contains
    subroutine spectrum(label, h, shift)
        character(*), intent(in) :: label
        real(dp), intent(in) :: h(:, :), shift
        real(dp), allocatable :: t(:, :), lam(:), wk(:)
        integer :: q, inf
        t = h
        do q = 1, size(t, 1)
            t(q, q) = t(q, q) + shift
        end do
        allocate(lam(size(t, 1)), wk(1))
        call dsyev('N', 'U', size(t, 1), t, size(t, 1), lam, wk, -1, inf)
        q = int(wk(1)); deallocate(wk); allocate(wk(q))
        call dsyev('N', 'U', size(t, 1), t, size(t, 1), lam, wk, q, inf)
        call report(label, lam)
    end subroutine
    subroutine pspectrum(label, h, shift, mm)
        character(*), intent(in) :: label
        real(dp), intent(in) :: h(:, :), shift, mm(:, :)
        real(dp), allocatable :: t(:, :), b(:, :), lam(:), wk(:)
        integer :: q, inf, its
        t = h
        do q = 1, size(t, 1)
            t(q, q) = t(q, q) + shift
        end do
        b = mm
        allocate(lam(size(t, 1)), wk(1))
        call dsygv(1, 'N', 'U', size(t, 1), t, size(t, 1), b, size(t, 1), lam, wk, -1, inf)
        q = int(wk(1)); deallocate(wk); allocate(wk(q))
        call dsygv(1, 'N', 'U', size(t, 1), t, size(t, 1), b, size(t, 1), lam, wk, q, inf)
        if (inf /= 0) then
            print '(a,a,i0)', label, ': dsygv info=', inf; return
        end if
        call report(label, lam)
        call dense_pcg(h, shift, mm, its)
        print '(a,a,i0)', label, ': dense PCG iterations (tol 1e-8, GL1 rhs)= ', its
    end subroutine
    subroutine report(label, lam)
        character(*), intent(in) :: label
        real(dp), intent(in) :: lam(:)
        real(dp) :: kap
        integer :: q
        kap = lam(size(lam))/lam(1)
        print '(a,a,es10.3,a,es10.3,a,es10.3,a,i0)', label, ': lmin=', lam(1), &
            ' lmax=', lam(size(lam)), ' kappa=', kap, ' CGbound=', &
            ceiling(0.5_dp*sqrt(kap)*log(2/1.e-8_dp))
        print '(a,9(es9.2,1x))', '   deciles: ', (lam(max(1, (q*size(lam))/10)), q=1, 9)
        print '(a,4(i0,1x))', '   count lam<1e1*lmin, <1e2*lmin, <1e3*lmin, <1e4*lmin: ', &
            count(lam < 1e1_dp*lam(1)), count(lam < 1e2_dp*lam(1)), &
            count(lam < 1e3_dp*lam(1)), count(lam < 1e4_dp*lam(1))
    end subroutine
    subroutine dense_pcg(h, shift, mm, its)
        real(dp), intent(in) :: h(:, :), shift, mm(:, :)
        integer, intent(out) :: its
        real(dp), allocatable :: l(:, :), xx(:), rr(:), zz(:), pp(:), qq(:)
        real(dp) :: rz, rzold, al, nb
        integer :: inf, q
        l = mm
        call dpotrf('L', size(l, 1), l, size(l, 1), inf)
        if (inf /= 0) then
            its = -1; return
        end if
        xx = 0*rhs; rr = rhs; nb = max(norm2(rhs), 1.0_dp)
        zz = rr; call dpotrs('L', size(l, 1), 1, l, size(l, 1), zz, size(l, 1), inf)
        pp = zz; rz = dot_product(rr, zz)
        do q = 1, 20000
            qq = matmul(h, pp) + shift*pp
            al = rz/dot_product(pp, qq)
            xx = xx + al*pp; rr = rr - al*qq
            if (norm2(rr) <= 1.e-8_dp*nb) exit
            zz = rr; call dpotrs('L', size(l, 1), 1, l, size(l, 1), zz, size(l, 1), inf)
            rzold = rz; rz = dot_product(rr, zz)
            pp = zz + (rz/rzold)*pp
        end do
        its = q
    end subroutine
end program
