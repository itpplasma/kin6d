module kin6d_tiago_island_mse
    ! Fixed-wall synthetic inference; the physical closure stays in KIN6D.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_rng, only: rng_t, rng_normal
    use tiago_local_posterior, only: local_posterior_t
    use tiago_diagnostic_types, only: mse_channel_t
    use tiago_mse, only: mse_geometry, mse_angles
    use kin6d_island_provider, only: island_provider_t, island_provider_solve, &
        island_provider_sample, island_provider_jvp, IP_OK, IP_EPSILON, IP_NTHETA
    implicit none
    private
    public :: island_mse_posterior_t, island_mse_predict

    type, extends(local_posterior_t) :: island_mse_posterior_t
        type(island_provider_t) :: provider
        type(mse_channel_t), allocatable :: channels(:)
        real(dp), allocatable :: observed(:)
        real(dp) :: sigma = 1.0e-4_dp, prior_mean = 0, prior_sd = 10
    contains
        procedure :: residual
        procedure :: jacobian
        procedure :: n_residuals
        procedure :: log_normalization
        procedure :: proper_prior
        procedure :: draw_prior
    end type
contains
    subroutine island_mse_predict(provider, channels, epsilon, angles, ok, jac)
        type(island_provider_t), intent(inout) :: provider
        real(dp), intent(in) :: epsilon
        type(mse_channel_t), intent(in) :: channels(:)
        real(dp), allocatable, intent(out) :: angles(:)
        logical, intent(out) :: ok
        real(dp), allocatable, intent(out), optional :: jac(:)
        real(dp) :: theta(IP_NTHETA), direction(IP_NTHETA)
        real(dp) :: xyz(3, size(channels)), coefficients(6)
        real(dp) :: values(4, size(channels)), jv(4, size(channels))
        real(dp), allocatable :: field_gradient(:, :)
        character(:), allocatable :: message
        integer :: status, i
        integer(int32) :: ierr

        ok = .false.
        do i = 1, size(channels)
            call mse_geometry(channels(i), xyz(:, i), coefficients, ierr, message)
            if (ierr /= 0) return
        end do
        theta = provider%theta
        theta(IP_EPSILON) = epsilon
        call island_provider_solve(provider, theta, status)
        if (status /= IP_OK) return
        call island_provider_sample(provider, xyz, values, status)
        if (status /= IP_OK) return
        call mse_angles(channels, values(1:3, :), angles, ierr, message, &
            field_gradient)
        if (ierr /= 0) return
        if (present(jac)) then
            direction = 0
            direction(IP_EPSILON) = 1
            call island_provider_jvp(provider, xyz, direction, jv, status)
            if (status /= IP_OK) return
            jac = sum(field_gradient*jv(1:3, :), dim=1)
        end if
        ok = .true.
    end subroutine

    subroutine residual(self, theta, r, ok)
        class(island_mse_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: r(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: angles(:)
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: n

        r = ieee_value(0.0_dp, ieee_quiet_nan)
        ok = .false.
        if (size(theta) /= 1 .or. size(r) /= self%n_residuals()) return
        if (.not. self%proper_prior() .or. self%sigma <= 0) return
        call island_mse_predict(self%provider, self%channels, theta(1), &
            angles, ok)
        if (.not. ok) return
        n = size(self%observed)
        r(:n) = (modulo(angles - self%observed + pi/2, pi) - pi/2)/self%sigma
        r(n + 1) = (theta(1) - self%prior_mean)/self%prior_sd
    end subroutine

    subroutine jacobian(self, theta, jac)
        class(island_mse_posterior_t), intent(inout) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), allocatable, intent(out) :: jac(:, :)
        real(dp), allocatable :: angles(:), derivative(:)
        logical :: ok
        integer :: n

        n = size(self%observed)
        allocate (jac(n + 1, 1))
        jac = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(theta) /= 1) return
        if (.not. self%proper_prior() .or. self%sigma <= 0) return
        call island_mse_predict(self%provider, self%channels, theta(1), &
            angles, ok, derivative)
        if (.not. ok) return
        jac(:n, 1) = derivative/self%sigma
        jac(n + 1, 1) = 1/self%prior_sd
    end subroutine

    integer function n_residuals(self)
        class(island_mse_posterior_t), intent(in) :: self
        n_residuals = size(self%observed) + 1
    end function

    real(dp) function log_normalization(self)
        class(island_mse_posterior_t), intent(in) :: self
        log_normalization = -0.5_dp*self%n_residuals()*log(2*acos(-1.0_dp)) &
            - size(self%observed)*log(self%sigma) - log(self%prior_sd)
    end function

    logical function proper_prior(self)
        class(island_mse_posterior_t), intent(in) :: self
        proper_prior = self%prior_sd > 0
    end function

    subroutine draw_prior(self, rng, theta, ok)
        class(island_mse_posterior_t), intent(in) :: self
        type(rng_t), intent(inout) :: rng
        real(dp), intent(out) :: theta(:)
        logical, intent(out) :: ok
        ok = size(theta) == 1 .and. self%proper_prior()
        if (.not. ok) return
        call rng_normal(rng, theta(1))
        theta(1) = self%prior_mean + self%prior_sd*theta(1)
    end subroutine
end module
