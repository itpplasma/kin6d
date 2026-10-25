module kin6d_map_spectrum_probe
    ! Measurement only (perf/landreman-w3-20261025): dense spectra and PCG
    ! iteration counts of the reduced position-gauge map Hessian with candidate
    ! preconditioners. Not a solver; never part of a PR.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_ideal_map, only: ideal_map_t, ideal_map_hessian, ideal_map_energy, &
        ideal_map_reduce, ideal_map_rdof, ideal_map_solve, ideal_map_solve_report_t
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    private
    public :: map_spectrum_probe
contains
    subroutine map_spectrum_probe(map, xseed)
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: xseed(:, :, :)
        real(dp), allocatable :: hs(:, :), hf(:, :), g(:), gr(:), rhs(:), m(:, :)
        integer, allocatable :: free(:), nmode_of(:), node_of(:)
        integer :: nf, i, j, a, mm, jj, status
        real(dp) :: energy, mu
        free = pack([(i, i=1, map%nred)], .not. map%fixed); nf = size(free)
        allocate(nmode_of(map%nred), node_of(map%nred))
        nmode_of = 0; node_of = 0
        do a = 1, map%nodes
            do mm = 1, map%nmode
                do jj = 1, 2
                    i = ideal_map_rdof(map, jj, mm, a)
                    nmode_of(i) = mm/2  ! toroidal n: 1->0, 2,3->1, 4,5->2, ...
                    node_of(i) = a
                end do
            end do
        end do
        call dense_hessian(map, map%x, free, hf)
        call dense_hessian(map, xseed, free, hs)
        allocate(g(map%ndof), gr(map%nred))
        call ideal_map_energy(map, xseed, energy, status, g)
        call ideal_map_reduce(map, g, gr)
        rhs = -gr(free)/norm2(gr(free))
        print '(a,i0,a,i0)', 'probe free=', nf, ' nred=', map%nred
        call eig('H_seed', hs)
        call eig('H_final', hf)
        mu = 1.e-6_dp  ! initial LM damping, mu*|diag| as ideal_map.f90 Newton
        do i = 1, nf
            hs(i, i) = hs(i, i) + mu*abs(hs(i, i))
        end do
        allocate(m(nf, nf))
        call pcg('none', hs, rhs, identity=.true.)
        m = 0
        do i = 1, nf
            m(i, i) = hs(i, i)
        end do
        call pcg('Jacobi', hs, rhs, m)
        m = 0  ! node blocks: all modes and both components of one node
        do j = 1, nf
            do i = 1, nf
                if (node_of(free(i)) == node_of(free(j))) m(i, j) = hs(i, j)
            end do
        end do
        call pcg('node block (2*nmode)', hs, rhs, m)
        m = 0  ! toroidal-mode block diagonal (n = n'), all nodes coupled
        do j = 1, nf
            do i = 1, nf
                if (nmode_of(free(i)) == nmode_of(free(j))) m(i, j) = hs(i, j)
            end do
        end do
        call pcg('toroidal-n block diag', hs, rhs, m)
        call pcg('lagged exact H_final (+mu) for seed system', hs, rhs, hf + 0*hf)
        ! Newton-lag reuse: precondition the system at iterate k+1 with the exact
        ! Hessian of iterate k (states from ideal_map_solve truncated at k).
        block
            type(ideal_map_t) :: work
            type(ideal_map_solve_report_t) :: rep
            real(dp), allocatable :: hk(:, :), hk1(:, :), gk(:)
            integer :: k, st
            integer, parameter :: its(5) = [1, 2, 4, 8, 16]
            do k = 1, size(its)
                work = map; work%x = xseed
                call ideal_map_solve(work, 1.e-30_dp, its(k), rep, st)
                call dense_hessian(work, work%x, free, hk)
                work = map; work%x = xseed
                call ideal_map_solve(work, 1.e-30_dp, its(k) + 1, rep, st)
                call dense_hessian(work, work%x, free, hk1)
                call ideal_map_energy(work, work%x, energy, status, g)
                call ideal_map_reduce(work, g, gr)
                gk = -gr(free)/norm2(gr(free))
                do i = 1, nf
                    hk1(i, i) = hk1(i, i) + mu*abs(hk1(i, i))
                end do
                write (*, '(a,i0)', advance='no') 'lag k=', its(k)
                call pcg(' H_k for H_k+1', hk1, gk, hk)
                m = 0
                do j = 1, nf
                    do i = 1, nf
                        if (nmode_of(free(i)) == nmode_of(free(j))) m(i, j) = hk(i, j)
                    end do
                end do
                write (*, '(a,i0)', advance='no') 'lag k=', its(k)
                call pcg(' n-block(H_k) for H_k+1', hk1, gk, m)
            end do
        end block
    contains
        subroutine pcg(label, h, b, mm, identity)
            character(*), intent(in) :: label
            real(dp), intent(in) :: h(:, :), b(:)
            real(dp), intent(in), optional :: mm(:, :)
            logical, intent(in), optional :: identity
            real(dp), allocatable :: l(:, :), x(:), r(:), z(:), p(:), q(:), lam(:), &
                t(:, :), bb(:, :), wk(:)
            real(dp) :: rz, rzo, al
            integer :: k, inf, n, lw
            logical :: ident
            ident = present(identity)
            n = size(b)
            if (.not. ident) then
                l = mm
                do k = 1, n  ! same relative damping as the system
                    l(k, k) = l(k, k) + 1.e-6_dp*abs(l(k, k))
                end do
                t = h; bb = l
                allocate(lam(n), wk(1))
                call dsygv(1, 'N', 'U', n, t, n, bb, n, lam, wk, -1, inf)
                lw = int(wk(1)); deallocate(wk); allocate(wk(lw))
                call dsygv(1, 'N', 'U', n, t, n, bb, n, lam, wk, lw, inf)
                if (inf /= 0) then
                    print '(a,a,i0)', label, ': M not SPD or dsygv info=', inf; return
                end if
                print '(a,a,3(es10.3,1x),a,i0,1x,i0)', label, ': M^-1A lmin lmax kappa=', &
                    lam(1), lam(n), lam(n)/lam(1), ' n(lam<0.1) n(lam>10)=', &
                    count(lam < 0.1_dp), count(lam > 10)
                call dpotrf('L', n, l, n, inf)
                if (inf /= 0) then
                    print '(a,a)', label, ': dpotrf failed'; return
                end if
            end if
            x = 0*b; r = b; z = r
            if (.not. ident) call dpotrs('L', n, 1, l, n, z, n, inf)
            p = z; rz = dot_product(r, z)
            do k = 1, 20000
                q = matmul(h, p)
                al = rz/dot_product(p, q)
                if (.not. al > 0) then
                    print '(a,a,i0)', label, ': CG breakdown/negative curvature at it ', k
                    return
                end if
                x = x + al*p; r = r - al*q
                if (norm2(r) <= 1.e-8_dp*max(1.0_dp, norm2(b))) exit
                z = r
                if (.not. ident) call dpotrs('L', n, 1, l, n, z, n, inf)
                rzo = rz; rz = dot_product(r, z); p = z + (rz/rzo)*p
            end do
            print '(a,a,i0)', label, ': PCG iterations (tol 1e-8, seed Newton rhs)= ', k
        end subroutine
    end subroutine

    subroutine dense_hessian(map, x, free, h)
        type(ideal_map_t), intent(inout) :: map
        real(dp), intent(in) :: x(:, :, :)
        integer, intent(in) :: free(:)
        real(dp), allocatable, intent(out) :: h(:, :)
        type(csr_matrix_t) :: matrix
        integer, allocatable :: pos(:)
        integer :: i, k
        call ideal_map_hessian(map, x, matrix)
        allocate(pos(map%nred)); pos = 0
        do i = 1, size(free)
            pos(free(i)) = i
        end do
        allocate(h(size(free), size(free))); h = 0
        do i = 1, map%nred
            if (pos(i) == 0) cycle
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (pos(matrix%colidx(k)) > 0) h(pos(i), pos(matrix%colidx(k))) = matrix%val(k)
            end do
        end do
        h = (h + transpose(h))/2
    end subroutine

    subroutine eig(label, h)
        character(*), intent(in) :: label
        real(dp), intent(in) :: h(:, :)
        real(dp), allocatable :: t(:, :), lam(:), wk(:)
        integer :: n, inf, lw, q
        n = size(h, 1); t = h; allocate(lam(n), wk(1))
        call dsyev('N', 'U', n, t, n, lam, wk, -1, inf)
        lw = int(wk(1)); deallocate(wk); allocate(wk(lw))
        call dsyev('N', 'U', n, t, n, lam, wk, lw, inf)
        print '(a,a,3(es10.3,1x),a,i0)', label, ': lmin lmax |lmax/lmin|=', lam(1), lam(n), &
            abs(lam(n)/lam(1)), ' negatives=', count(lam < 0)
        print '(a,9(es9.2,1x))', '   deciles: ', (lam(max(1, (q*n)/10)), q=1, 9)
        print '(a,3(i0,1x))', '   count lam<1e-6*lmax, <1e-4*lmax, <1e-2*lmax: ', &
            count(lam < 1e-6_dp*lam(n)), count(lam < 1e-4_dp*lam(n)), count(lam < 1e-2_dp*lam(n))
    end subroutine
end module
