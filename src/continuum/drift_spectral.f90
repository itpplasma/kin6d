!> Doubly periodic Fourier grid for the fluid drift-wave models of
!> docs/drift-turbulence.md.
!>
!> A real field f(x, y) on [0, lx) x [0, ly) is stored as its Fourier
!> coefficients fh(i, j), f = sum fh exp[i(kx(i) x + ky(j) y)], with
!> fh(-k) = conjg(fh(k)). Coefficients outside the retained band
!> |m| <= (n - 1)/3 per direction are zero, which removes the aliases of
!> quadratic products. The box average of a product is
!> <a b> = sum real(conjg(ah) bh) (drift_mean_product).
!>
!> The bracket {phi, f} = phi_x f_y - phi_y f_x is the average of its
!> advective and its two conservative forms. It is bilinear, so its tangent is
!> {dphi, f} + {phi, df}; under <.,.> its adjoint in f is g -> -{phi, g} and
!> its adjoint in phi is g -> {f, g}.
module milca_drift_spectral
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use fortnum_fft, only: fortnum_fft_plan_t, fortnum_fft_workspace_t, &
        fft_c2c_plan_init, fft_workspace_init, fft_c2c
    use milca_drift_fluid_leaves, only: bracket_split_products
    implicit none
    private
    public :: drift_grid_t, drift_work_t, drift_grid_init, drift_work_init
    public :: drift_to_spectral, drift_to_physical, drift_bracket
    public :: drift_poisson_inverse, drift_mean_product

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp

    type :: drift_grid_t
        integer :: nx = 0, ny = 0
        real(dp) :: lx = 0, ly = 0
        !> Signed wave numbers in transform order and kx^2 + ky^2.
        real(dp), allocatable :: kx(:), ky(:), k2(:, :)
        !> Retained band.
        logical, allocatable :: keep(:, :)
        !> Index of the coefficient at -k.
        integer, allocatable :: mx(:), my(:)
        type(fortnum_fft_plan_t) :: px, py
    end type

    !> Scratch of one concurrent caller.
    type :: drift_work_t
        type(fortnum_fft_workspace_t), allocatable :: wx(:), wy(:)
        complex(dp), allocatable :: line(:, :), c(:, :, :)
    end type
contains
    subroutine drift_grid_init(grid, nx, ny, lx, ly)
        type(drift_grid_t), intent(out) :: grid
        integer, intent(in) :: nx, ny
        real(dp), intent(in) :: lx, ly
        integer :: i, j

        if (nx < 4 .or. ny < 4) error stop 'drift_grid_init: need nx, ny >= 4'
        if (.not. (lx > 0 .and. ly > 0)) error stop 'drift_grid_init: need lx, ly > 0'
        grid%nx = nx; grid%ny = ny; grid%lx = lx; grid%ly = ly
        allocate (grid%kx(nx), grid%ky(ny), grid%k2(nx, ny), grid%keep(nx, ny))
        allocate (grid%mx(nx), grid%my(ny))
        do i = 1, nx
            grid%kx(i) = 2*pi/lx*signed_index(i, nx)
            grid%mx(i) = modulo(nx - i + 1, nx) + 1
        end do
        do j = 1, ny
            grid%ky(j) = 2*pi/ly*signed_index(j, ny)
            grid%my(j) = modulo(ny - j + 1, ny) + 1
        end do
        do j = 1, ny
            do i = 1, nx
                grid%k2(i, j) = grid%kx(i)**2 + grid%ky(j)**2
                grid%keep(i, j) = 3*abs(signed_index(i, nx)) <= nx - 1 .and. &
                                  3*abs(signed_index(j, ny)) <= ny - 1
            end do
        end do
        call fft_c2c_plan_init(grid%px, nx)
        call fft_c2c_plan_init(grid%py, ny)
    end subroutine

    pure integer function signed_index(i, n)
        integer, intent(in) :: i, n
        signed_index = i - 1
        if (signed_index > n/2) signed_index = signed_index - n
    end function

    subroutine drift_work_init(work, grid)
        type(drift_work_t), intent(out) :: work
        type(drift_grid_t), intent(in) :: grid
        integer :: threads, t

        threads = 1
        !$ threads = omp_get_max_threads()
        allocate (work%wx(threads), work%wy(threads), work%line(grid%ny, threads))
        allocate (work%c(grid%nx, grid%ny, 3))
        do t = 1, threads
            call fft_workspace_init(work%wx(t), grid%px)
            call fft_workspace_init(work%wy(t), grid%py)
        end do
    end subroutine

    !> In-place two-dimensional transform; sign = -1 forward without
    !> normalisation, +1 backward.
    subroutine transform(grid, wx, wy, line, a, sign)
        type(drift_grid_t), intent(in) :: grid
        type(fortnum_fft_workspace_t), intent(inout) :: wx(:), wy(:)
        complex(dp), intent(inout) :: line(:, :), a(:, :)
        integer, intent(in) :: sign
        integer :: i, j, t

        !$omp parallel private(i, j, t)
        t = 1
        !$ t = omp_get_thread_num() + 1
        !$omp do schedule(static)
        do j = 1, grid%ny
            call fft_c2c(a(:, j), sign, grid%px, wx(t))
        end do
        !$omp end do
        !$omp do schedule(static)
        do i = 1, grid%nx
            line(:, t) = a(i, :)
            call fft_c2c(line(:, t), sign, grid%py, wy(t))
            a(i, :) = line(:, t)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine

    !> Coefficients of a real field, restricted to the retained band.
    subroutine drift_to_spectral(grid, work, f, fh)
        type(drift_grid_t), intent(in) :: grid
        type(drift_work_t), intent(inout) :: work
        real(dp), intent(in) :: f(:, :)
        complex(dp), intent(out) :: fh(:, :)

        fh = f
        call transform(grid, work%wx, work%wy, work%line, fh, -1)
        where (grid%keep)
            fh = fh/(grid%nx*grid%ny)
        elsewhere
            fh = 0
        end where
    end subroutine

    subroutine drift_to_physical(grid, work, fh, f)
        type(drift_grid_t), intent(in) :: grid
        type(drift_work_t), intent(inout) :: work
        complex(dp), intent(in) :: fh(:, :)
        real(dp), intent(out) :: f(:, :)

        work%c(:, :, 1) = fh
        call transform(grid, work%wx, work%wy, work%line, work%c(:, :, 1), 1)
        f = real(work%c(:, :, 1), dp)
    end subroutine

    !> jh = coefficients of {phi, f}. phih and fh must lie in the retained band.
    !> Two real fields share one complex transform.
    subroutine drift_bracket(grid, work, phih, fh, jh)
        type(drift_grid_t), intent(in) :: grid
        type(drift_work_t), intent(inout) :: work
        complex(dp), intent(in) :: phih(:, :), fh(:, :)
        complex(dp), intent(out) :: jh(:, :)
        complex(dp), parameter :: iu = (0, 1)
        complex(dp) :: p12, p12m, p3
        real(dp) :: p1, p2, q3, scale
        integer :: i, j

        associate (c => work%c, wx => work%wx, wy => work%wy, line => work%line)
            !$omp parallel do private(i, j) schedule(static)
            do j = 1, grid%ny
                do i = 1, grid%nx
                    ! (phi_x, phi_y), (phi, f), (f_x, f_y) as real and imaginary parts.
                    c(i, j, 1) = iu*(grid%kx(i) + iu*grid%ky(j))*phih(i, j)
                    c(i, j, 2) = phih(i, j) + iu*fh(i, j)
                    c(i, j, 3) = iu*(grid%kx(i) + iu*grid%ky(j))*fh(i, j)
                end do
            end do
            !$omp end parallel do
            call transform(grid, wx, wy, line, c(:, :, 1), 1)
            call transform(grid, wx, wy, line, c(:, :, 2), 1)
            call transform(grid, wx, wy, line, c(:, :, 3), 1)
            !$omp parallel do private(i, j, p1, p2, q3) schedule(static)
            do j = 1, grid%ny
                do i = 1, grid%nx
                    call bracket_split_products(real(c(i, j, 1), dp), aimag(c(i, j, 1)), &
                        real(c(i, j, 2), dp), real(c(i, j, 3), dp), aimag(c(i, j, 3)), &
                        aimag(c(i, j, 2)), p1, p2, q3)
                    c(i, j, 1) = cmplx(p1, p2, dp)
                    c(i, j, 2) = q3
                end do
            end do
            !$omp end parallel do
            call transform(grid, wx, wy, line, c(:, :, 1), -1)
            call transform(grid, wx, wy, line, c(:, :, 2), -1)
            scale = 1.0_dp/(3*grid%nx*grid%ny)
            !$omp parallel do private(i, j, p12, p12m, p3) schedule(static)
            do j = 1, grid%ny
                do i = 1, grid%nx
                    if (.not. grid%keep(i, j)) then
                        jh(i, j) = 0
                        cycle
                    end if
                    p12 = c(i, j, 1)
                    p12m = conjg(c(grid%mx(i), grid%my(j), 1))
                    p3 = c(i, j, 2)
                    jh(i, j) = scale*((p12 + p12m)/2 + grid%kx(i)*(p12 - p12m)/2 &
                                      + iu*grid%ky(j)*p3)
                end do
            end do
            !$omp end parallel do
        end associate
    end subroutine

    !> phih of zero mean with laplace(phi) = zeta.
    subroutine drift_poisson_inverse(grid, zetah, phih)
        type(drift_grid_t), intent(in) :: grid
        complex(dp), intent(in) :: zetah(:, :)
        complex(dp), intent(out) :: phih(:, :)

        where (grid%k2 > 0)
            phih = -zetah/grid%k2
        elsewhere
            phih = 0
        end where
    end subroutine

    !> Box average <a b> of two real fields from their coefficients.
    function drift_mean_product(ah, bh) result(mean)
        complex(dp), intent(in) :: ah(:, :), bh(:, :)
        real(dp) :: mean
        integer :: i, j

        mean = 0
        !$omp parallel do private(i, j) reduction(+:mean) schedule(static)
        do j = 1, size(ah, 2)
            do i = 1, size(ah, 1)
                mean = mean + real(conjg(ah(i, j))*bh(i, j), dp)
            end do
        end do
        !$omp end parallel do
    end function
end module milca_drift_spectral
