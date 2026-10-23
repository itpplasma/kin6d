module kin6d_response_strip_preconditioner
    !! Fourier-diagonal radial-band approximation of the full strip action.
    !! Local normal/tangent components remove the rotating poloidal frame.
    !! Modes remain independent only in the preconditioner; FGMRES applies
    !! the unmodified coupled-harmonic operator. No global Jacobian is stored.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c, &
        fortnum_fft_workspace_t, fft_workspace_init
    use fortnum_krylov, only: complex_fgmres_operator, KRYLOV_OK, KRYLOV_STAGNATION
    use kin6d_response_toroidal_strip, only: toroidal_strip_t, &
        toroidal_strip_work_t, toroidal_strip_status_t, toroidal_strip_action, &
        toroidal_strip_vjp, toroidal_strip_defects
    use fortnum_krylov, only: complex_linear_operator_t, &
        complex_gmres, complex_fgmres
    implicit none
    private
    integer, parameter :: bw = 41, ld = 3*bw + 1, colors = 9
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    interface
        subroutine zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, kl, ku, ldab
            complex(dp), intent(inout) :: ab(ldab, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
        subroutine zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb, ipiv(*)
            complex(dp), intent(in) :: ab(ldab, *)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
    type, public :: strip_preconditioner_t
        private
        integer :: nr = 0, nt = 0
        logical :: ready = .false.
        complex(dp) :: gamma = 0
        complex(dp), allocatable :: lu(:, :, :), spectral(:, :, :), line(:, :)
        real(dp), allocatable :: normal(:, :, :), rows(:, :), cols(:, :)
        real(dp), allocatable :: scale(:, :, :)
        integer, allocatable :: piv(:, :), solve_info(:)
        type(fortnum_fft_workspace_t), allocatable :: fft_work(:, :)
        type(fortnum_fft_plan_t) :: fft
    end type
    public :: strip_preconditioner_init, strip_preconditioner_apply
    public :: strip_preconditioner_bytes, strip_preconditioned_solve
    type, extends(complex_linear_operator_t) :: strip_preconditioned_solve_matvec_t
        type(strip_preconditioner_t), pointer :: pc => null()
        logical, pointer :: reverse => null()
        type(toroidal_strip_t), pointer :: strip => null()
        complex(dp), pointer :: weighted(:) => null()
        type(toroidal_strip_work_t), pointer :: work => null()
    contains
        procedure :: apply => strip_preconditioned_solve_matvec_apply
    end type
    type, extends(complex_linear_operator_t) :: &
        strip_preconditioned_solve_precondition_t
        type(strip_preconditioner_t), pointer :: pc => null()
        integer, pointer :: pre_info => null()
        logical, pointer :: reverse => null()
    contains
        procedure :: apply => strip_preconditioned_solve_precondition_apply
    end type
contains
    subroutine strip_preconditioner_init(pc, strip, work, gamma, inverse, scale, info)
        type(strip_preconditioner_t), intent(out) :: pc
        type(toroidal_strip_t), intent(in) :: strip
        type(toroidal_strip_work_t), intent(inout) :: work
        complex(dp), intent(in) :: gamma
        real(dp), intent(in) :: inverse(:, :, :, :), scale(:, :, :)
        integer, intent(out) :: info
        complex(dp), allocatable :: x(:, :, :), y(:, :, :), mode(:)
        complex(dp), allocatable :: band(:, :), samples(:, :, :)
        real(dp) :: length, value
        integer :: nr, nt, nv, m, k, a, color, i, j, row, col, jj, mi
        integer :: lo, first, last, probes, p, modes(9)
        info = 1
        nr = size(scale, 2); nt = size(scale, 3); nv = 6*nr
        if (any(shape(scale) /= [6, nr, nt])) return
        if (any(shape(inverse) /= [2, 2, nr, nt])) return
        if (nr < 5 .or. nt < 5 .or. mod(nt, 2) == 0) return
        if (any(.not. ieee_is_finite(scale)) .or. any(scale <= 0)) return
        if (any(.not. ieee_is_finite(inverse))) return
        if (.not. ieee_is_finite(real(gamma, dp))) return
        if (.not. ieee_is_finite(aimag(gamma))) return
        if (abs(gamma) <= tiny(1.0_dp)) return
        probes = min(nt, 9)
        modes(1:3) = [0, 1, -1]
        if (nt == 5) then
            modes(4:5) = [2, -2]
        else
            m = (nt - 1)/2
            modes(4:9) = [m, -m, m - 1, 1 - m, m - 2, 2 - m]
        end if
        pc%nr = nr; pc%nt = nt; pc%gamma = gamma; pc%scale = scale
        allocate (pc%lu(ld, nv, nt), pc%piv(nv, nt), pc%normal(2, nr, nt), &
            pc%rows(nv, nt), pc%cols(nv, nt), pc%spectral(nt, 6, nr), &
            pc%solve_info(nt), pc%line(nv, nt), pc%fft_work(6, nr), &
            x(6, nr, nt), y(6, nr, nt), &
            mode(nt), band(ld, nv), samples(ld, nv, probes))
        call fft_c2c_plan_init(pc%fft, nt)
        do i = 1, nr
            do a = 1, 6
                call fft_workspace_init(pc%fft_work(a, i), pc%fft)
            end do
        end do
        do j = 1, nt
            do i = 1, nr
                length = norm2(inverse(1, :, i, j))
                if (length <= tiny(1.0_dp)) return
                pc%normal(:, i, j) = inverse(1, :, i, j)/length
            end do
        end do
        ! Stream colored columns. D has five local nodes; D D spans at most
        ! nine nodes (one-sided boundary rows reach six radii away).
        ! Thus nine colors recover every retained band entry without aliasing.
        do k = 1, probes
            m = modes(k)
            do j = 1, nt
                mode(j) = exp(iu*m*2*acos(-1.0_dp)*(j - 1)/nt)
            end do
            band = 0
            do a = 1, 6
                do color = 1, min(colors, nr)
                    x = 0
                    do i = color, nr, colors
                        x(a, i, :) = mode
                    end do
                    call rotate(pc, x, .false., .false.)
                    call toroidal_strip_action(strip, work, gamma, x, y, .true.)
                    call rotate(pc, y, .true., .true.)
                    do i = 1, nr
                        lo = min(max(i - 2, 1), nr - 4)
                        first = min(max(lo - 2, 1), nr - 4)
                        last = min(max(lo + 2, 1), nr - 4) + 4
                        do jj = first, last
                            if (mod(jj - color, colors) /= 0) cycle
                            col = 6*(jj - 1) + a
                            do mi = 1, 6
                                row = 6*(i - 1) + mi
                                if (abs(row - col) > bw) cycle
                                band(2*bw + 1 + row - col, col) = &
                                    sum(conjg(mode)*y(mi, i, :))/nt
                            end do
                        end do
                    end do
                end do
            end do
            samples(:, :, k) = band
        end do
        ! The continuous angular symbol is quadratic in m. Three probes retain
        ! its zeroth/first/second derivative terms, avoiding nt full probes.
        ! Six edge probes retain the dominant high-mode collocation aliasing.
        ! Remaining alias tails are omitted only from the preconditioner;
        ! the full operator always determines the solution.
        !$omp parallel do private(m,p,col,row,value) schedule(static)
        factors: do k = 1, nt
            pc%solve_info(k) = 1
            m = k - 1
            if (m > (nt - 1)/2) m = m - nt
            pc%lu(:, :, k) = (1.0_dp - real(m, dp)**2)*samples(:, :, 1) &
                + (real(m, dp)**2 + m)*samples(:, :, 2)/2 &
                + (real(m, dp)**2 - m)*samples(:, :, 3)/2
            do p = 1, probes
                if (m == modes(p)) pc%lu(:, :, k) = samples(:, :, p)
            end do
            ! Two-sided equilibration controls mixed SI units and gamma/eta
            ! magnitudes. No diagonal shift or change of the target equations.
            pc%rows(:, k) = 0
            do col = 1, nv
                do row = max(1, col - bw), min(nv, col + bw)
                    value = abs(pc%lu(2*bw + 1 + row - col, col, k))
                    pc%rows(row, k) = max(pc%rows(row, k), value)
                end do
            end do
            if (any(pc%rows(:, k) <= tiny(1.0_dp))) cycle factors
            do col = 1, nv
                do row = max(1, col - bw), min(nv, col + bw)
                    pc%lu(2*bw + 1 + row - col, col, k) = &
                        pc%lu(2*bw + 1 + row - col, col, k)/pc%rows(row, k)
                end do
                pc%cols(col, k) = maxval(abs(pc%lu(:, col, k)))
                if (pc%cols(col, k) <= tiny(1.0_dp)) cycle factors
                pc%lu(:, col, k) = pc%lu(:, col, k)/pc%cols(col, k)
            end do
            call zgbtrf(nv, nv, bw, bw, pc%lu(:, :, k), ld, pc%piv(:, k), &
                pc%solve_info(k))
        end do factors
        !$omp end parallel do
        info = maxval(abs(pc%solve_info))
        if (info /= 0) return
        pc%ready = .true.
        info = 0
    end subroutine
    subroutine rotate(pc, x, forward, equations)
        type(strip_preconditioner_t), intent(in) :: pc
        complex(dp), intent(inout) :: x(6, pc%nr, pc%nt)
        logical, intent(in) :: forward, equations
        complex(dp) :: r, z
        real(dp) :: c, s
        integer :: i, j
        !$omp parallel do collapse(2) private(c,s,r,z) schedule(static)
        do j = 1, pc%nt
            do i = 1, pc%nr
                if (equations .and. (i == 1 .or. i == pc%nr)) cycle
                c = pc%normal(1, i, j); s = pc%normal(2, i, j)
                if (.not. forward) s = -s
                r = x(1, i, j); z = x(2, i, j)
                x(1, i, j) = c*r + s*z; x(2, i, j) = -s*r + c*z
                r = x(3, i, j); z = x(5, i, j)
                x(3, i, j) = c*r + s*z; x(5, i, j) = -s*r + c*z
            end do
        end do
        !$omp end parallel do
    end subroutine
    subroutine strip_preconditioner_apply(pc, rhs, x, info, adjoint)
        !! Hermitian storage transpose when adjoint=.true. Scaling is included
        !! so this is the inverse approximation for A/scale, not for raw A.
        type(strip_preconditioner_t), intent(inout) :: pc
        complex(dp), intent(in) :: rhs(6, pc%nr, pc%nt)
        complex(dp), intent(out) :: x(6, pc%nr, pc%nt)
        integer, intent(out) :: info
        logical, intent(in), optional :: adjoint
        logical :: reverse
        character :: trans
        integer :: i, j, a, row
        info = 1
        if (.not. pc%ready) return
        reverse = .false.
        if (present(adjoint)) reverse = adjoint
        trans = 'N'
        if (reverse) trans = 'C'
        x = rhs
        if (.not. reverse) x = x*pc%scale
        call rotate(pc, x, .true.,.not. reverse)
        !$omp parallel do collapse(2) schedule(static)
        do i = 1, pc%nr
            do a = 1, 6
                pc%spectral(:, a, i) = x(a, i, :)
                call fft_c2c(pc%spectral(:, a, i), -1, pc%fft, pc%fft_work(a, i))
            end do
        end do
        !$omp end parallel do
        !$omp parallel do private(i,a,row) schedule(static)
        do j = 1, pc%nt
            do i = 1, pc%nr
                do a = 1, 6
                    row = 6*(i - 1) + a
                    if (reverse) then
                        pc%line(row, j) = pc%spectral(j, a, i)/pc%cols(row, j)
                    else
                        pc%line(row, j) = pc%spectral(j, a, i)/pc%rows(row, j)
                    end if
                end do
            end do
            call zgbtrs(trans, 6*pc%nr, bw, bw, 1, pc%lu(:, :, j), ld, &
                pc%piv(:, j), pc%line(:, j), 6*pc%nr, pc%solve_info(j))
            do i = 1, pc%nr
                do a = 1, 6
                    row = 6*(i - 1) + a
                    if (reverse) then
                        pc%spectral(j, a, i) = pc%line(row, j)/pc%rows(row, j)
                    else
                        pc%spectral(j, a, i) = pc%line(row, j)/pc%cols(row, j)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        info = maxval(abs(pc%solve_info))
        !$omp parallel do collapse(2) schedule(static)
        do i = 1, pc%nr
            do a = 1, 6
                call fft_c2c(pc%spectral(:, a, i), 1, pc%fft, pc%fft_work(a, i))
                x(a, i, :) = pc%spectral(:, a, i)/pc%nt
            end do
        end do
        !$omp end parallel do
        call rotate(pc, x, .false., reverse)
        if (reverse) x = x*pc%scale
    end subroutine
    integer(int64) function strip_preconditioner_bytes(pc) result(bytes)
        type(strip_preconditioner_t), intent(in) :: pc
        integer :: fft_length
        bytes = 0
        if (.not. allocated(pc%lu)) return
        fft_length = pc%nt
        if (pc%fft%bluestein) fft_length = pc%fft%m
        bytes = int(2*6*pc%nr*fft_length, int64)*16 &
            + int(size(pc%lu), int64)*16 &
            + int(size(pc%spectral) + size(pc%line), int64)*16 &
            + int(size(pc%normal) + size(pc%rows) + size(pc%cols) &
            + size(pc%scale), int64)*8 + int(size(pc%piv), int64)*4
    end function
    subroutine scale_vector(pc, values)
        !! Divide the canonical storage vector without materializing a reshaped
        !! copy of the full equation-scale field in every Krylov callback.
        type(strip_preconditioner_t), intent(in) :: pc
        complex(dp), intent(inout) :: values(:)
        integer :: a, i, j, k
        !$omp parallel do collapse(2) if(pc%nr*pc%nt >= 4096) &
        !$omp private(a,k) schedule(static)
        do j = 1, pc%nt
            do i = 1, pc%nr
                do a = 1, 6
                    k = 6*((j - 1)*pc%nr + i - 1) + a
                    values(k) = values(k)/pc%scale(a, i, j)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine
    subroutine strip_preconditioned_solve(pc, strip, work, rhs, x, tol, &
        maxit, restart, st, adjoint)
        type(strip_preconditioner_t), intent(inout), target :: pc
        type(toroidal_strip_t), intent(in), target :: strip
        type(toroidal_strip_work_t), intent(inout), target :: work
        complex(dp), intent(in) :: rhs(6, pc%nr, pc%nt)
        complex(dp), intent(inout) :: x(6, pc%nr, pc%nt)
        real(dp), intent(in) :: tol
        integer, intent(in) :: maxit, restart
        type(toroidal_strip_status_t), intent(out) :: st
        logical, intent(in), optional :: adjoint
        complex(dp), allocatable , target :: b(:), v(:), ax(:), weighted(:)
        logical , target :: reverse
        integer , target :: nn, pre_info
        real(dp) :: rhs_scale
        type(strip_preconditioned_solve_matvec_t) :: matvec_operator
        type(strip_preconditioned_solve_precondition_t) :: precondition_operator
        if (.not. pc%ready) return
        if (.not. ieee_is_finite(tol) .or. tol <= 0) return
        if (maxit < 1 .or. restart < 1) return
        reverse = .false.
        if (present(adjoint)) reverse = adjoint
        nn = 6*pc%nr*pc%nt
        allocate (b(nn), v(nn), ax(nn), weighted(nn))
        b = reshape(rhs, [nn]); v = reshape(x, [nn])
        if (.not. reverse) call scale_vector(pc, b)
        rhs_scale = norm2(abs(b))
        if (rhs_scale <= 0.0_dp) rhs_scale = 1.0_dp
        b = b/rhs_scale; v = v/rhs_scale
        pre_info = 0
        ! Reserve half the residual budget for reconstruction and the final
        ! physical action. Krylov stagnation at this tighter internal target
        ! is acceptable only if the independently recomputed residual meets tol.
        matvec_operator%pc => pc
        matvec_operator%reverse => reverse
        matvec_operator%strip => strip
        nullify(matvec_operator%weighted)
        if (allocated(weighted)) matvec_operator%weighted => weighted
        matvec_operator%work => work
        precondition_operator%pc => pc
        precondition_operator%pre_info => pre_info
        precondition_operator%reverse => reverse
        call complex_fgmres(matvec_operator, precondition_operator, b, v, 0.5_dp*tol, &
            maxit, restart, st%info, st%iterations, st%residual)
        v = v*rhs_scale
        b = reshape(rhs, [nn])
        if (.not. reverse) call scale_vector(pc, b)
        call matvec_operator%apply(v, ax)
        st%residual = norm2(abs(ax - b))/max(norm2(abs(b)), tiny(1.0_dp))
        if (st%info == KRYLOV_STAGNATION) then
            if (st%residual <= tol) st%info = KRYLOV_OK
        end if
        if (st%residual > tol .and. st%info == 0) st%info = 1
        if (pre_info /= 0) st%info = pre_info
        x = reshape(v, shape(x))
        if (.not. reverse) call toroidal_strip_defects(strip, work, st%div_b, st%div_j)
    end subroutine
    subroutine strip_preconditioned_solve_matvec_apply(self, input, output)
        class(strip_preconditioned_solve_matvec_t), intent(inout) :: self
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            if (self%reverse) then
                self%weighted = input
                call scale_vector(self%pc, self%weighted)
                call toroidal_strip_vjp(self%strip, self%work, self%pc%gamma, &
                    self%weighted, output, .true.)
            else
                call toroidal_strip_action(self%strip, self%work, self%pc%gamma, &
                    input, output, .true.)
                call scale_vector(self%pc, output)
            end if
    end subroutine strip_preconditioned_solve_matvec_apply
    subroutine strip_preconditioned_solve_precondition_apply(self, input, output)
        class(strip_preconditioned_solve_precondition_t), intent(inout) :: self
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            integer :: status
            call strip_preconditioner_apply(self%pc, input, output, status, &
                self%reverse)
            if (status /= 0) self%pre_info = status
    end subroutine strip_preconditioned_solve_precondition_apply
end module
