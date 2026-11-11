!> Doubly periodic spectral grid: transforms, Poisson inverse and the split
!> bracket against a manufactured solution, its invariants and adjoints to
!> rounding, and spectral convergence. A collocated advective bracket without
!> truncation (test-local mutant) must violate the quadratic invariant.
program test_drift_spectral
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fft, only: fft_c2c
    use milca_drift_spectral
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: failures

    failures = 0
    call manufactured()
    call invariants()
    call convergence()
    if (failures > 0) error stop 'test_drift_spectral failed'
    write (*, '(a)') 'test_drift_spectral passed'
contains
    subroutine manufactured()
        integer, parameter :: nx = 30, ny = 40
        real(dp), parameter :: lx = 2*pi/0.15_dp, ly = 3*pi
        type(drift_grid_t) :: grid
        type(drift_work_t) :: work
        real(dp) :: x, y, a, b, c, d
        real(dp), dimension(nx, ny) :: phi, f, bracket, zeta, back
        complex(dp), dimension(nx, ny) :: phih, fh, jh, zetah, invh
        integer :: i, j

        call drift_grid_init(grid, nx, ny, lx, ly)
        call drift_work_init(work, grid)
        a = 2*pi/lx*3; b = 2*pi/ly*2; c = 2*pi/lx*2; d = -2*pi/ly*4
        do j = 1, ny
            do i = 1, nx
                x = lx*(i - 1)/nx; y = ly*(j - 1)/ny
                phi(i, j) = sin(a*x)*cos(b*y) + 0.3_dp*cos(c*x - b*y)
                f(i, j) = cos(c*x + d*y + 0.4_dp) + 0.7_dp
                zeta(i, j) = -(a**2 + b**2)*sin(a*x)*cos(b*y) &
                             - 0.3_dp*(c**2 + b**2)*cos(c*x - b*y)
                bracket(i, j) = (a*cos(a*x)*cos(b*y) - 0.3_dp*c*sin(c*x - b*y)) &
                                *(-d*sin(c*x + d*y + 0.4_dp)) &
                                - (-b*sin(a*x)*sin(b*y) + 0.3_dp*b*sin(c*x - b*y)) &
                                *(-c*sin(c*x + d*y + 0.4_dp))
            end do
        end do
        call drift_to_spectral(grid, work, phi, phih)
        call drift_to_spectral(grid, work, f, fh)
        call drift_to_physical(grid, work, phih, back)
        call check('round trip', maxval(abs(back - phi)), 1e-14_dp)
        call check('mean coefficient', abs(fh(1, 1) - 0.7_dp), 1e-15_dp)
        call drift_bracket(grid, work, phih, fh, jh)
        call drift_to_physical(grid, work, jh, back)
        call check('bracket vs manufactured solution', &
                   maxval(abs(back - bracket))/maxval(abs(bracket)), 1e-13_dp)
        call drift_to_spectral(grid, work, zeta, zetah)
        call drift_poisson_inverse(grid, zetah, invh)
        call check('Poisson inverse', maxval(abs(invh - phih)), 1e-14_dp)
        call check('mean product', abs(drift_mean_product(fh, fh) - 0.99_dp), 1e-14_dp)
    end subroutine

    subroutine invariants()
        integer, parameter :: nx = 48, ny = 36
        type(drift_grid_t) :: grid
        type(drift_work_t) :: work
        real(dp), dimension(nx, ny) :: phi, f, g, fx, fy, px, py
        complex(dp), dimension(nx, ny) :: phih, fh, gh, jh, kh
        real(dp) :: scale, pairing(3)

        call drift_grid_init(grid, nx, ny, 7.0_dp, 5.0_dp)
        call drift_work_init(work, grid)
        call noise(phi, 1); call noise(f, 2); call noise(g, 3)
        call drift_to_spectral(grid, work, phi, phih)
        call drift_to_spectral(grid, work, f, fh)
        call drift_to_spectral(grid, work, g, gh)
        if (count(abs(fh) > 0) /= 31*23) failures = failures + 1
        call drift_bracket(grid, work, phih, fh, jh)
        scale = sqrt(drift_mean_product(jh, jh))
        if (.not. scale > 0.1_dp) failures = failures + 1
        call check('<{phi,f}> = 0', abs(jh(1, 1))/scale, 1e-14_dp)
        call check('<f {phi,f}> = 0', abs(drift_mean_product(fh, jh)) &
                   /(scale*sqrt(drift_mean_product(fh, fh))), 1e-14_dp)
        call check('<phi {phi,f}> = 0', abs(drift_mean_product(phih, jh)) &
                   /(scale*sqrt(drift_mean_product(phih, phih))), 1e-14_dp)
        call drift_bracket(grid, work, fh, phih, kh)
        call check('antisymmetry', maxval(abs(jh + kh))/scale, 1e-14_dp)
        call drift_bracket(grid, work, fh, fh, kh)
        call check('{f,f} = 0', maxval(abs(kh))/scale, 1e-14_dp)
        ! Adjoints: <g {phi,f}> = -<f {phi,g}> = <phi {f,g}>.
        pairing(1) = drift_mean_product(gh, jh)
        call drift_bracket(grid, work, phih, gh, kh)
        pairing(2) = -drift_mean_product(fh, kh)
        call drift_bracket(grid, work, fh, gh, kh)
        pairing(3) = drift_mean_product(phih, kh)
        if (.not. abs(pairing(1)) > 1e-3_dp*scale) failures = failures + 1
        call check('adjoint in f and in phi', &
                   maxval(abs(pairing(2:3) - pairing(1)))/abs(pairing(1)), 1e-12_dp)

        ! Mutant: advective form at the collocation points, fields not truncated.
        call derivative(phi, 7.0_dp, 5.0_dp, px, py)
        call derivative(f, 7.0_dp, 5.0_dp, fx, fy)
        g = px*fy - py*fx
        call check_fails('advective bracket without truncation violates <f {phi,f}> = 0', &
                         abs(sum(f*g))/sqrt(sum(f**2)*sum(g**2)), 1e-4_dp)
    end subroutine

    subroutine convergence()
        integer, parameter :: sizes(3) = [12, 24, 48]
        type(drift_grid_t) :: grid
        type(drift_work_t) :: work
        real(dp), allocatable, dimension(:, :) :: phi, f, exact, back
        complex(dp), allocatable, dimension(:, :) :: phih, fh, jh
        real(dp) :: x, y, error(3)
        integer :: i, j, k, n

        do k = 1, 3
            n = sizes(k)
            call drift_grid_init(grid, n, n, 2*pi, 2*pi)
            call drift_work_init(work, grid)
            allocate (phi(n, n), f(n, n), exact(n, n), back(n, n))
            allocate (phih(n, n), fh(n, n), jh(n, n))
            do j = 1, n
                do i = 1, n
                    x = 2*pi*(i - 1)/n; y = 2*pi*(j - 1)/n
                    phi(i, j) = exp(sin(x) + cos(y))
                    f(i, j) = exp(cos(x - y))
                    exact(i, j) = phi(i, j)*f(i, j)*sin(x - y)*(cos(x) - sin(y))
                end do
            end do
            call drift_to_spectral(grid, work, phi, phih)
            call drift_to_spectral(grid, work, f, fh)
            call drift_bracket(grid, work, phih, fh, jh)
            call drift_to_physical(grid, work, jh, back)
            error(k) = maxval(abs(back - exact))/maxval(abs(exact))
            deallocate (phi, f, exact, back, phih, fh, jh)
        end do
        write (*, '(a,3es10.2)') 'bracket error at n = 12, 24, 48: ', error
        call check('spectral convergence 12 -> 24', error(2)/error(1), 1e-2_dp)
        call check('spectral convergence 24 -> 48', error(3)/error(2), 1e-4_dp)
        call check('bracket error at n = 48', error(3), 1e-9_dp)
    end subroutine

    !> Reproducible field of zero mean filling the whole collocation grid.
    subroutine noise(f, seed)
        real(dp), intent(out) :: f(:, :)
        integer, intent(in) :: seed
        integer :: i, j
        real(dp) :: u
        do j = 1, size(f, 2)
            do i = 1, size(f, 1)
                u = sin(12.9898_dp*i + 78.233_dp*j + 37.719_dp*seed)*43758.5453_dp
                f(i, j) = u - floor(u) - 0.5_dp
            end do
        end do
        f = f - sum(f)/size(f)
    end subroutine

    !> Collocation derivatives of the full discrete spectrum.
    subroutine derivative(f, lx, ly, fx, fy)
        real(dp), intent(in) :: f(:, :), lx, ly
        real(dp), intent(out) :: fx(:, :), fy(:, :)
        complex(dp) :: row(size(f, 1)), col(size(f, 2))
        integer :: i, j, m, nx, ny

        nx = size(f, 1); ny = size(f, 2)
        do j = 1, ny
            row = f(:, j)
            call fft_c2c(row, -1)
            do i = 1, nx
                m = i - 1
                if (m > nx/2) m = m - nx
                if (2*m == nx) m = 0
                row(i) = (0, 1)*2*pi/lx*m*row(i)/nx
            end do
            call fft_c2c(row, 1)
            fx(:, j) = real(row, dp)
        end do
        do i = 1, nx
            col = f(i, :)
            call fft_c2c(col, -1)
            do j = 1, ny
                m = j - 1
                if (m > ny/2) m = m - ny
                if (2*m == ny) m = 0
                col(j) = (0, 1)*2*pi/ly*m*col(j)/ny
            end do
            call fft_c2c(col, 1)
            fy(i, :) = real(col, dp)
        end do
    end subroutine

    subroutine check(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,es10.3)') 'FAIL ', name//': ', value, ' > ', tolerance
            failures = failures + 1
        end if
    end subroutine

    subroutine check_fails(name, value, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: value, tolerance
        if (value > tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3)') 'FAIL ', name//': ', value
            failures = failures + 1
        end if
    end subroutine
end program test_drift_spectral
