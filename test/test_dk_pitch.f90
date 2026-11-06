!> Lorentz pitch-angle operator on a uniform field: Legendre eigenfunctions,
!! symmetry, null space and the convergence of odd-l eigenvalues.
program test_dk_pitch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_pitch, only: pitch_operator_t, build_pitch_operator
    implicit none
    interface
        subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
            import :: dp
            integer, intent(in) :: itype, n, lda, ldb, lwork
            character, intent(in) :: jobz, uplo
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsygv
    end interface
    real(dp), parameter :: b = 1.3_dp, nu = 0.7_dp
    integer :: failures = 0

    call even_legendre_modes()
    call symmetry_and_null_space()
    call odd_mode_convergence()
    if (failures > 0) error stop 'FAIL'
    print '(a)', 'All tests passed!'
contains

    ! P_l(xi), xi = sigma sqrt(B) t, is a polynomial of degree l in t on each
    ! side, so for l <= p it lies in the discrete space and
    ! K v = nu l(l+1)/2 M v, for both parities.
    subroutine even_legendre_modes()
        type(pitch_operator_t) :: op
        real(dp), allocatable :: v(:), xi(:), r(:)
        integer :: l, p
        p = 3
        call build_pitch_operator(b, nu, 4, p, op)
        allocate (xi(op%n))
        xi = op%sigma*sqrt(max(0.0_dp, 1 - op%lambda*b))
        do l = 0, p
            v = legendre(l, xi)
            r = matmul(op%stiffness, v) - nu*l*(l + 1)/2.0_dp*matmul(op%mass, v)
            call check_small('P_'//digit(l)//' eigen-residual', &
                maxval(abs(r))/maxval(abs(matmul(op%mass, v))), 1e-12_dp)
        end do
    end subroutine even_legendre_modes

    subroutine symmetry_and_null_space()
        type(pitch_operator_t) :: op
        real(dp), allocatable :: mu(:)
        call build_pitch_operator(b, nu, 5, 2, op)
        call check_small('K symmetric', maxval(abs(op%stiffness - &
            transpose(op%stiffness)))/maxval(abs(op%stiffness)), 1e-14_dp)
        call check_small('K 1 = 0', maxval(abs(sum(op%stiffness, dim=2))) &
            /maxval(abs(op%stiffness)), 1e-13_dp)
        mu = eigenvalues(op)
        call check_small('lowest eigenvalue zero', abs(mu(1))/nu, 1e-12_dp)
        if (mu(2) < 0.5_dp*nu) then
            print '(a,es12.4)', 'FAIL null space larger than constants: ', mu(2)
            failures = failures + 1
        end if
    end subroutine symmetry_and_null_space

    ! The first mode outside the space (l = p + 1 = 3): check and report the rate.
    subroutine odd_mode_convergence()
        type(pitch_operator_t) :: op
        real(dp), allocatable :: mu(:)
        real(dp) :: err(4), rate
        integer :: k, n_el
        do k = 1, 4
            n_el = 2**k
            call build_pitch_operator(b, nu, n_el, 2, op)
            mu = eigenvalues(op)
            err(k) = abs(mu(4) - 6*nu)/(6*nu)
        end do
        rate = log(err(3)/err(4))/log(2.0_dp)
        print '(a,4es10.2,a,f5.2)', 'l=3 eigenvalue error, n_el 2..16, p=2:', err, &
            '  rate', rate
        ! l = p + 1 is not in the space; its eigenvalue converges like h^(2p).
        if (.not. all(err(2:4) < err(1:3)) .or. rate < 3.5_dp) then
            print '(a)', 'FAIL l=3 eigenvalue does not converge'
            failures = failures + 1
        end if
    end subroutine odd_mode_convergence

    function eigenvalues(op) result(mu)
        type(pitch_operator_t), intent(in) :: op
        real(dp), allocatable :: mu(:), a(:, :), m(:, :), work(:)
        integer :: info
        a = op%stiffness; m = op%mass
        allocate (mu(op%n), work(64*op%n))
        call dsygv(1, 'N', 'U', op%n, a, op%n, m, op%n, mu, work, size(work), info)
        if (info /= 0) error stop 'dsygv failed'
    end function eigenvalues

    elemental function legendre(l, x) result(pl)
        integer, intent(in) :: l
        real(dp), intent(in) :: x
        real(dp) :: pl, p0, p1
        integer :: k
        p0 = 1; p1 = x
        if (l == 0) then
            pl = 1
            return
        end if
        do k = 1, l - 1
            pl = ((2*k + 1)*x*p1 - k*p0)/(k + 1)
            p0 = p1; p1 = pl
        end do
        pl = p1
    end function legendre

    function digit(l) result(s)
        integer, intent(in) :: l
        character(:), allocatable :: s
        character(8) :: buf
        write (buf, '(i0)') l
        s = trim(buf)
    end function digit

    subroutine check_small(label, value, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tol
        if (.not. (value <= tol)) then
            print '(a,es12.4)', 'FAIL '//label//': ', value
            failures = failures + 1
        end if
    end subroutine check_small
end program test_dk_pitch
