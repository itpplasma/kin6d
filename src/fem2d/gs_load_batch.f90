module kin6d_gs_load_batch
    !! Vectorised kernels of the fixed-rule smooth-cell load. All dummy arguments
    !! are explicit-shape and contiguous, so the loops over the quadrature
    !! points vectorise (`omp simd`). The quadrature index runs fastest in the
    !! transposed shape table shapes_t(np, ndof).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: batch_position_flux, batch_moments

contains

    subroutine batch_position_flux(np, ndof, points, local_psi, shapes_t, r, z, psi)
        !! Physical position and flux at all quadrature points of one cell.
        integer, intent(in) :: np, ndof
        real(dp), intent(in) :: points(2, ndof), local_psi(ndof), shapes_t(np, ndof)
        real(dp), intent(out) :: r(np), z(np), psi(np)
        integer :: a, q
        r = 0; z = 0; psi = local_psi(1)
        do a = 1, ndof
            !$omp simd
            do q = 1, np
                r(q) = r(q) + points(1, a)*shapes_t(q, a)
                z(q) = z(q) + points(2, a)*shapes_t(q, a)
                psi(q) = psi(q) + (local_psi(a) - local_psi(1))*shapes_t(q, a)
            end do
        end do
    end subroutine

    subroutine batch_moments(np, ndof, need_jacobian, c, f, df, af, shapes_t, &
                             integral)
        !! integral = [sum c f s_a, sum c df s_a s_b (b-major), sum c af s_a].
        !! c holds the point weights (rule weight times map determinant).
        integer, intent(in) :: np, ndof
        logical, intent(in) :: need_jacobian
        real(dp), intent(in) :: c(np), f(np), df(np), af(np), shapes_t(np, ndof)
        real(dp), intent(out) :: integral(:)
        real(dp) :: cf(np), cd(np), tmp(np), acc
        integer :: a, b, q
        !$omp simd
        do q = 1, np
            cf(q) = c(q)*f(q)
        end do
        do a = 1, ndof
            acc = 0
            !$omp simd reduction(+:acc)
            do q = 1, np
                acc = acc + cf(q)*shapes_t(q, a)
            end do
            integral(a) = acc
        end do
        if (.not. need_jacobian) return
        !$omp simd
        do q = 1, np
            cd(q) = c(q)*df(q)
            cf(q) = c(q)*af(q)
        end do
        do b = 1, ndof
            !$omp simd
            do q = 1, np
                tmp(q) = cd(q)*shapes_t(q, b)
            end do
            do a = 1, b
                acc = 0
                !$omp simd reduction(+:acc)
                do q = 1, np
                    acc = acc + tmp(q)*shapes_t(q, a)
                end do
                integral(ndof + (b - 1)*ndof + a) = acc
                integral(ndof + (a - 1)*ndof + b) = acc
            end do
        end do
        do a = 1, ndof
            acc = 0
            !$omp simd reduction(+:acc)
            do q = 1, np
                acc = acc + cf(q)*shapes_t(q, a)
            end do
            integral(ndof*(ndof + 1) + a) = acc
        end do
    end subroutine

end module kin6d_gs_load_batch
