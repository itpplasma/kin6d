module kin6d_response_fluid
    !! Point residual of non-rotating, single-fluid R2-R4 with full pressure.
    !! Cartesian physical jets are supplied by the spatial discretisation.
    !! gamma = -i omega; eta is electrical resistivity in Ohm m (SI).
    !! This kernel does not discretise or solve a strip or an outer interface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_response_fluid_leaves, only: fluid_pressure, fluid_force, &
                                           fluid_pressure_vjp, fluid_force_vjp
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

    public :: response_fluid_residual, response_fluid_batch, response_fluid_adjoint
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

    pure subroutine response_fluid_adjoint(eq, gamma, seed, adj)
        !! Fixed-equilibrium complex Hermitian transpose of the point action.
        !! Parameter cotangents are not returned. The generated real leaves
        !! define both parts; complex time factors are transposed explicitly.
        type(response_fluid_equilibrium_t), intent(in) :: eq
        complex(dp), intent(in) :: gamma
        type(response_fluid_result_t), intent(in) :: seed
        type(response_fluid_jet_t), intent(out) :: adj
        real(dp) :: sm(3, 2), so(3, 2), sp(4, 2)
        real(dp) :: ac(3, 2), el(3, 2), vl(3, 2), jc(3, 2), bc(3, 2)
        real(dp) :: xp(3, 2), dx(3, 3, 2), gd(3, 2), gp(3, 2), unused(20)
        integer :: a, k, part

        sm(:, 1) = real(seed%momentum, dp); sm(:, 2) = aimag(seed%momentum)
        so(:, 1) = real(seed%ohm, dp); so(:, 2) = aimag(seed%ohm)
        sp(1, 1) = real(seed%p1, dp); sp(1, 2) = aimag(seed%p1)
        sp(2:4, 1) = real(seed%gp1, dp); sp(2:4, 2) = aimag(seed%gp1)
        do part = 1, 2
            call fluid_force_vjp( &
                eq%rho, eq%eta, eq%b0(1), eq%j0(1), &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, eq%b0(2), eq%j0(2), &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, eq%b0(3), eq%j0(3), &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, sm(1, part), sm(2, part), &
                sm(3, part), so(1, part), so(2, part), so(3, part), &
                unused(1), unused(2), unused(3), unused(4), ac(1, part), &
                el(1, part), vl(1, part), jc(1, part), bc(1, part), gp(1, part), &
                unused(5), unused(6), ac(2, part), el(2, part), vl(2, part), &
                jc(2, part), bc(2, part), gp(2, part), unused(7), unused(8), &
                ac(3, part), el(3, part), vl(3, part), jc(3, part), bc(3, part), &
                gp(3, part))
            sp(2:4, part) = sp(2:4, part) + gp(:, part)
            call fluid_pressure_vjp( &
                eq%p, eq%adiabatic, &
                eq%gp(1), 0.0_dp, 0.0_dp, eq%hp(1, 1), 0.0_dp, &
                eq%hp(1, 2), 0.0_dp, eq%hp(1, 3), 0.0_dp, &
                eq%gp(2), 0.0_dp, 0.0_dp, eq%hp(2, 1), 0.0_dp, &
                eq%hp(2, 2), 0.0_dp, eq%hp(2, 3), 0.0_dp, &
                eq%gp(3), 0.0_dp, 0.0_dp, eq%hp(3, 1), 0.0_dp, &
                eq%hp(3, 2), 0.0_dp, eq%hp(3, 3), 0.0_dp, &
                sp(1, part), sp(2, part), sp(3, part), sp(4, part), &
                unused(1), unused(2), unused(3), xp(1, part), gd(1, part), &
                unused(4), dx(1, 1, part), unused(5), dx(1, 2, part), &
                unused(6), dx(1, 3, part), unused(7), xp(2, part), gd(2, part), &
                unused(8), dx(2, 1, part), unused(9), dx(2, 2, part), &
                unused(10), dx(2, 3, part), unused(11), xp(3, part), gd(3, part), &
                unused(12), dx(3, 1, part), unused(13), dx(3, 2, part), &
                unused(14), dx(3, 3, part))
        end do
        do a = 1, 3
            adj%a(a) = -conjg(gamma)*cmplx(el(a, 1), el(a, 2), dp)
            adj%xi(a) = conjg(gamma**2)*cmplx(ac(a, 1), ac(a, 2), dp) &
                        + conjg(gamma)*cmplx(vl(a, 1), vl(a, 2), dp) &
                        + cmplx(xp(a, 1), xp(a, 2), dp)
            adj%grad_phi(a) = -cmplx(el(a, 1), el(a, 2), dp)
            adj%b(a) = cmplx(bc(a, 1), bc(a, 2), dp)
            adj%j(a) = cmplx(jc(a, 1), jc(a, 2), dp)
            adj%gd(a) = cmplx(gd(a, 1), gd(a, 2), dp)
            do k = 1, 3
                adj%dxi(a, k) = cmplx(dx(a, k, 1), dx(a, k, 2), dp)
            end do
        end do
    end subroutine response_fluid_adjoint

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
