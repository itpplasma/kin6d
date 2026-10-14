module kin6d_response_fluid
    !! Point residual of non-rotating, single-fluid R2-R4 with full pressure.
    !! Cartesian physical jets are supplied by the spatial discretisation.
    !! gamma = -i omega; eta is electrical resistivity in Ohm m (SI).
    !! This kernel does not discretise or solve a strip or an outer interface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_fluid_leaves, only: fluid_pressure, fluid_force
    implicit none
    private

    type, public :: response_fluid_equilibrium_t
        real(dp) :: b0(3) = 0.0_dp, j0(3) = 0.0_dp
        real(dp) :: p = 0.0_dp, gp(3) = 0.0_dp, hp(3, 3) = 0.0_dp
        real(dp) :: rho = 1.0_dp, eta = 0.0_dp, adiabatic = 5.0_dp/3.0_dp
    end type response_fluid_equilibrium_t

    type, public :: response_fluid_jet_t
        complex(dp) :: a(3) = (0.0_dp, 0.0_dp), xi(3) = (0.0_dp, 0.0_dp)
        ! dxi(a,k) = d xi_a / d x_k; gd = grad(div xi).
        complex(dp) :: dxi(3, 3) = (0.0_dp, 0.0_dp), gd(3) = (0.0_dp, 0.0_dp)
        complex(dp) :: grad_phi(3) = (0.0_dp, 0.0_dp)
        ! b = curl A; j = curl b / mu0, supplied by the magnetic operator.
        complex(dp) :: b(3) = (0.0_dp, 0.0_dp), j(3) = (0.0_dp, 0.0_dp)
    end type response_fluid_jet_t

    type, public :: response_fluid_result_t
        complex(dp) :: p1 = (0.0_dp, 0.0_dp), gp1(3) = (0.0_dp, 0.0_dp)
        complex(dp) :: momentum(3) = (0.0_dp, 0.0_dp)
        complex(dp) :: ohm(3) = (0.0_dp, 0.0_dp)
    end type response_fluid_result_t

    public :: response_fluid_residual, response_fluid_batch
contains
    pure subroutine response_fluid_residual(eq, jet, gamma, res)
        type(response_fluid_equilibrium_t), intent(in) :: eq
        type(response_fluid_jet_t), intent(in) :: jet
        complex(dp), intent(in) :: gamma
        type(response_fluid_result_t), intent(out) :: res
        real(dp) :: xp(3, 2), dxp(3, 3, 2), gdp(3, 2), pressure(4, 2)
        real(dp) :: ap(3, 2), ep(3, 2), vp(3, 2), jp(3, 2), bp(3, 2), force(6, 2)
        complex(dp) :: acc, electric, velocity
        integer :: a, k, part

        do a = 1, 3
            xp(a, 1) = real(jet%xi(a), dp)
            xp(a, 2) = aimag(jet%xi(a))
            gdp(a, 1) = real(jet%gd(a), dp)
            gdp(a, 2) = aimag(jet%gd(a))
            do k = 1, 3
                dxp(a, k, 1) = real(jet%dxi(a, k), dp)
                dxp(a, k, 2) = aimag(jet%dxi(a, k))
            end do
            acc = gamma**2*jet%xi(a)
            electric = -gamma*jet%a(a) - jet%grad_phi(a)
            velocity = gamma*jet%xi(a)
            ap(a, 1) = real(acc, dp)
            ap(a, 2) = aimag(acc)
            ep(a, 1) = real(electric, dp)
            ep(a, 2) = aimag(electric)
            vp(a, 1) = real(velocity, dp)
            vp(a, 2) = aimag(velocity)
            jp(a, 1) = real(jet%j(a), dp)
            jp(a, 2) = aimag(jet%j(a))
            bp(a, 1) = real(jet%b(a), dp)
            bp(a, 2) = aimag(jet%b(a))
        end do
        do part = 1, 2
            call pressure_part(eq, xp(:, part), dxp(:, :, part), gdp(:, part), &
                               pressure(:, part))
            call force_part(eq, ap(:, part), ep(:, part), vp(:, part), jp(:, part), &
                            bp(:, part), pressure(2:4, part), force(:, part))
        end do
        res%p1 = cmplx(pressure(1, 1), pressure(1, 2), dp)
        do a = 1, 3
            res%gp1(a) = cmplx(pressure(a + 1, 1), pressure(a + 1, 2), dp)
            res%momentum(a) = cmplx(force(a, 1), force(a, 2), dp)
            res%ohm(a) = cmplx(force(a + 3, 1), force(a + 3, 2), dp)
        end do
    end subroutine response_fluid_residual

    subroutine response_fluid_batch(eq, jet, gamma, res)
        type(response_fluid_equilibrium_t), intent(in) :: eq(:)
        type(response_fluid_jet_t), intent(in) :: jet(size(eq))
        complex(dp), intent(in) :: gamma(size(eq))
        type(response_fluid_result_t), intent(out) :: res(size(eq))
        integer :: k
        !$omp parallel do schedule(static)
        do k = 1, size(eq)
            call response_fluid_residual(eq(k), jet(k), gamma(k), res(k))
        end do
        !$omp end parallel do
    end subroutine response_fluid_batch

    pure subroutine pressure_part(eq, xi, dx, gd, out)
        type(response_fluid_equilibrium_t), intent(in) :: eq
        real(dp), intent(in) :: xi(3), dx(3, 3), gd(3)
        real(dp), intent(out) :: out(4)
        call fluid_pressure(eq%p, eq%adiabatic, eq%gp(1), xi(1), gd(1), &
                  eq%hp(1, 1), dx(1, 1), eq%hp(1, 2), dx(1, 2), eq%hp(1, 3), dx(1, 3), &
                            eq%gp(2), xi(2), gd(2), eq%hp(2, 1), dx(2, 1), &
                            eq%hp(2, 2), dx(2, 2), eq%hp(2, 3), dx(2, 3), &
                            eq%gp(3), xi(3), gd(3), eq%hp(3, 1), dx(3, 1), &
                            eq%hp(3, 2), dx(3, 2), eq%hp(3, 3), dx(3, 3), &
                            out(1), out(2), out(3), out(4))
    end subroutine pressure_part

    pure subroutine force_part(eq, acc, electric, velocity, current, magnetic, gp, out)
        type(response_fluid_equilibrium_t), intent(in) :: eq
        real(dp), intent(in) :: acc(3), electric(3), velocity(3)
        real(dp), intent(in) :: current(3), magnetic(3), gp(3)
        real(dp), intent(out) :: out(6)
        call fluid_force(eq%rho, eq%eta, eq%b0(1), eq%j0(1), acc(1), &
                         electric(1), velocity(1), current(1), magnetic(1), gp(1), &
                         eq%b0(2), eq%j0(2), acc(2), electric(2), velocity(2), &
                         current(2), magnetic(2), gp(2), eq%b0(3), eq%j0(3), acc(3), &
                         electric(3), velocity(3), current(3), magnetic(3), gp(3), &
                         out(1), out(2), out(3), out(4), out(5), out(6))
    end subroutine force_part
end module kin6d_response_fluid
