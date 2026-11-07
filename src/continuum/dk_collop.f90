!> Linearised collision operator of the drift-kinetic solver on an energy
!! grid (#184): speed x = v/v_ta of each species on nodal Gauss-Lobatto-
!! Legendre (GLL) elements on [0, x_max], perturbation f_a = F_a h_a P_l(xi)
!! with F_a the Maxwellian and h_a in the element space of degree >= 2, so
!! that 1, x and x^2 (density, momentum, energy) lie in the space.
!!
!! The operator is stored per Legendre degree l as the Galerkin matrix
!!   a(i, j) = int phi_i P_l C_a[f] d^3v,   f_b = phi_j F_b P_l,
!! rows (species a, node i), columns (species b, node j), in m^-3 s^-1.
!! Its integrands are FortSym leaves (gen/gen_dk_collop.f90); the integrals
!! are composite Gauss-Legendre rules on panels split at element edges and
!! at the other species' edges mapped to this species' speed, converged to
!! rounding. The streaming terms of the kinetic equation are diagonal in the
!! nodes with the lumped GLL mass (collop_mass); M^-1 a gives nodal rates.
!!
!! Models: Lorentz (pitch-angle scattering only), test particle (pitch-angle
!! and energy scattering off Maxwellian partners, the thermal term retained
!! for unequal temperatures) and Krook (-nu_K f, a parity switch, never a
!! default). The test-particle operator is symmetric for equal temperatures
!! and has the Maxwellian (h = 1) as null vector; it exchanges momentum and
!! energy with its partners, which are reservoirs here.
!!
!! Tangents with respect to the species temperatures and densities
!! (collop_build_jvp) come from the leaf JVPs at fixed normalised grid.
module kin6d_dk_collop
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre, gauss_lobatto_legendre
    use fortnum_polynomial, only: lagrange_weights, lagrange_differentiation_matrix
    use kin6d_dk_collop_leaves, only: collop_tp, collop_tp_series, collop_tp_jvp, &
        collop_tp_series_jvp
    implicit none
    private
    public :: collop_species_t, collop_grid_t, collop_t
    public :: collop_make_grid, collop_build, collop_build_jvp, collop_mass
    public :: collop_basis, collop_nuhat, collop_vth
    integer, parameter, public :: COLLOP_LORENTZ = 1, COLLOP_TEST_PARTICLE = 2, &
        COLLOP_KROOK = 4
    integer, parameter, public :: COLLOP_OK = 0, COLLOP_BAD_INPUT = 1

    real(dp), parameter, public :: collop_qe = 1.602176634e-19_dp
    real(dp), parameter, public :: collop_eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    ! Gauss-Legendre points per panel; integrands are polynomials times
    ! Gaussians and erf on panels no wider than one thermal speed.
    integer, parameter :: nquad = 32

    type :: collop_species_t
        real(dp) :: m = 0 ! mass (kg)
        real(dp) :: z = 0 ! charge number
        real(dp) :: n = 0 ! density (m^-3)
        real(dp) :: t = 0 ! temperature (J)
    end type

    !> Continuous GLL elements on [0, x_max]: n = nelem*(order - 1) + 1 nodes.
    type :: collop_grid_t
        integer :: order = 0, nelem = 0, n = 0
        real(dp), allocatable :: edges(:), x(:), w(:), xi(:), d(:, :)
    end type

    type :: collop_t
        integer :: model = 0, ns = 0, lmax = -1
        real(dp) :: nu_krook = 0
        type(collop_grid_t) :: grid
        type(collop_species_t), allocatable :: sp(:)
        real(dp), allocatable :: lnlambda(:, :)
        ! Galerkin matrices (ns*n, ns*n, 0:lmax), species-major blocks.
        real(dp), allocatable :: a(:, :, :)
    end type

contains

    !> Grid on [0, x_max] with the element edges given, or nelem uniform
    !> elements, each with order GLL points (order >= 3: degree >= 2).
    subroutine collop_make_grid(x_max, nelem, order, grid, status, edges)
        real(dp), intent(in) :: x_max
        integer, intent(in) :: nelem, order
        type(collop_grid_t), intent(out) :: grid
        integer, intent(out) :: status
        real(dp), intent(in), optional :: edges(0:)
        real(dp), allocatable :: w(:)
        integer :: e, k, g

        status = COLLOP_BAD_INPUT
        if (order < 3 .or. nelem < 1 .or. x_max <= 0) return
        grid%order = order
        grid%nelem = nelem
        grid%n = nelem*(order - 1) + 1
        allocate (grid%edges(0:nelem), grid%xi(order), w(order), grid%d(order, order), &
            grid%x(grid%n), grid%w(grid%n))
        if (present(edges)) then
            if (size(edges) /= nelem + 1) return
            grid%edges = edges
            if (any(edges(1:) <= edges(:nelem - 1)) .or. edges(0) /= 0) return
        else
            grid%edges = [(x_max*e/nelem, e=0, nelem)]
        end if
        call gauss_lobatto_legendre(order, grid%xi, w)
        call lagrange_differentiation_matrix(order, grid%xi, grid%d)
        grid%w = 0
        do e = 1, nelem
            do k = 1, order
                g = (e - 1)*(order - 1) + k
                grid%x(g) = element_point(grid, e, grid%xi(k))
                grid%w(g) = grid%w(g) + w(k)*(grid%edges(e) - grid%edges(e - 1))/2
            end do
        end do
        grid%x(1) = 0
        status = COLLOP_OK
    end subroutine collop_make_grid

    pure real(dp) function element_point(grid, e, xi)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: e
        real(dp), intent(in) :: xi
        element_point = grid%edges(e - 1) + (xi + 1)*(grid%edges(e) - grid%edges(e - 1))/2
    end function element_point

    !> Values and x-derivatives of the order basis functions of the element
    !> holding x, and the global index of the first; phi = 0 beyond the grid.
    subroutine collop_basis(grid, x, first, phi, dphi)
        type(collop_grid_t), intent(in) :: grid
        real(dp), intent(in) :: x
        integer, intent(out) :: first
        real(dp), intent(out) :: phi(:), dphi(:)
        real(dp) :: xi, h
        integer :: e
        phi = 0
        dphi = 0
        first = 1
        if (x < 0 .or. x > grid%edges(grid%nelem)) return
        e = 1
        do while (e < grid%nelem)
            if (x <= grid%edges(e)) exit
            e = e + 1
        end do
        h = grid%edges(e) - grid%edges(e - 1)
        xi = 2*(x - grid%edges(e - 1))/h - 1
        call lagrange_weights(grid%order, xi, grid%xi, phi)
        dphi = matmul(phi, grid%d)*2/h
        first = (e - 1)*(grid%order - 1) + 1
    end subroutine collop_basis

    pure real(dp) function collop_vth(s)
        type(collop_species_t), intent(in) :: s
        collop_vth = sqrt(2*s%t/s%m)
    end function collop_vth

    !> nuhat_ab = n_b e_a^2 e_b^2 lnL/(4 pi eps0^2 m_a^2 v_ta^3) in 1/s.
    pure real(dp) function collop_nuhat(a, b, lnlambda)
        type(collop_species_t), intent(in) :: a, b
        real(dp), intent(in) :: lnlambda
        collop_nuhat = b%n*(a%z*b%z)**2*collop_qe**4*lnlambda &
            /(4*pi*collop_eps0**2*a%m**2*collop_vth(a)**3)
    end function collop_nuhat

    !> Lumped GLL mass int phi_k^2 F_a P_l^2 d^3v of every node and species
    !> (ns*n), in m^-3: diagonal of the streaming terms.
    function collop_mass(op, l) result(m)
        type(collop_t), intent(in) :: op
        integer, intent(in) :: l
        real(dp) :: m(op%ns*op%grid%n)
        integer :: a, n
        n = op%grid%n
        do a = 1, op%ns
            m((a - 1)*n + 1:a*n) = 4*pi/(2*l + 1)*op%sp(a)%n/pi**1.5_dp &
                *op%grid%w*op%grid%x**2*exp(-op%grid%x**2)
        end do
    end function collop_mass

    !> Galerkin matrices of degrees 0 ... lmax. lnlambda(a, b) is the Coulomb
    !> logarithm of each pair; nu_krook [1/s] is used by the Krook model only.
    subroutine collop_build(species, lnlambda, grid, model, lmax, op, status, nu_krook)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, lmax
        type(collop_t), intent(out) :: op
        integer, intent(out) :: status
        real(dp), intent(in), optional :: nu_krook
        real(dp), allocatable :: zero(:)

        status = COLLOP_BAD_INPUT
        if (.not. setup(species, lnlambda, grid, model, lmax, op)) return
        if (model == COLLOP_KROOK) then
            if (.not. present(nu_krook)) return
            if (nu_krook <= 0) return
            op%nu_krook = nu_krook
        end if
        allocate (zero(op%ns))
        zero = 0
        call assemble(op, zero, zero, op%a)
        status = COLLOP_OK
    end subroutine collop_build

    !> Tangent da of the matrices for temperature and density tangents dt
    !> [J] and dn [m^-3] of each species, at fixed normalised speeds x.
    subroutine collop_build_jvp(species, lnlambda, grid, model, lmax, dt, dn, op, da, &
            status, nu_krook)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :), dt(:), dn(:)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, lmax
        type(collop_t), intent(out) :: op
        real(dp), allocatable, intent(out) :: da(:, :, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: nu_krook

        call collop_build(species, lnlambda, grid, model, lmax, op, status, nu_krook)
        if (status /= COLLOP_OK) return
        status = COLLOP_BAD_INPUT
        if (size(dt) /= op%ns .or. size(dn) /= op%ns) return
        call assemble(op, dt, dn, da)
        status = COLLOP_OK
    end subroutine collop_build_jvp

    logical function setup(species, lnlambda, grid, model, lmax, op) result(ok)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, lmax
        type(collop_t), intent(out) :: op
        integer :: a
        ok = .false.
        if (size(species) < 1 .or. lmax < 0 .or. grid%n < 3) return
        if (size(lnlambda, 1) /= size(species) .or. size(lnlambda, 2) /= size(species)) return
        if (all(model /= [COLLOP_LORENTZ, COLLOP_TEST_PARTICLE, COLLOP_KROOK])) return
        do a = 1, size(species)
            if (species(a)%m <= 0 .or. species(a)%t <= 0 .or. species(a)%n <= 0) return
        end do
        op%model = model
        op%ns = size(species)
        op%lmax = lmax
        op%grid = grid
        op%sp = species
        op%lnlambda = lnlambda
        ok = .true.
    end function setup

    !> Matrices (dt = dn = 0) or their tangents along (dt, dn).
    subroutine assemble(op, dt, dn, a)
        type(collop_t), intent(in) :: op
        real(dp), intent(in) :: dt(:), dn(:)
        real(dp), allocatable, intent(out) :: a(:, :, :)
        integer :: n, sa, sb, l, ra, rb
        real(dp), allocatable :: block(:, :, :)
        n = op%grid%n
        allocate (a(op%ns*n, op%ns*n, 0:op%lmax))
        a = 0
        do sa = 1, op%ns
            ra = (sa - 1)*n
            if (op%model == COLLOP_KROOK) then
                call krook_block(op, sa, dn(sa), any(dt /= 0) .or. any(dn /= 0), block)
                a(ra + 1:ra + n, ra + 1:ra + n, :) = block
                cycle
            end if
            do sb = 1, op%ns
                rb = (sb - 1)*n
                call test_particle_block(op, sa, sb, dt, dn, block)
                do l = 0, op%lmax
                    a(ra + 1:ra + n, ra + 1:ra + n, l) = &
                        a(ra + 1:ra + n, ra + 1:ra + n, l) + block(:, :, l)
                end do
            end do
        end do
    end subroutine assemble

    !> Test-particle block of species a on Maxwellian partner b, all l. With
    !> pref = 4 pi/(2l+1) n_a nuhat_ab pi^(-3/2):
    !>   -pref [l(l+1)/2 int phi_i phi_j lorentz + int phi_i' (energy phi_j'
    !>   + thermal phi_j)] dx.
    subroutine test_particle_block(op, sa, sb, dt, dn, block)
        type(collop_t), intent(in) :: op
        integer, intent(in) :: sa, sb
        real(dp), intent(in) :: dt(:), dn(:)
        real(dp), allocatable, intent(out) :: block(:, :, :)
        real(dp), allocatable :: xq(:), wq(:)
        real(dp) :: rinv, tau, nu, drinv, dtau, dnu, lor, en, th, dlor, den, dth
        real(dp) :: phi(op%grid%order), dphi(op%grid%order), pref, scale
        real(dp), allocatable :: lorentz(:, :), energy(:, :)
        integer :: q, i, j, first, l, p, n
        logical :: tangent

        n = op%grid%n
        p = op%grid%order
        tangent = any(dt /= 0) .or. any(dn /= 0)
        rinv = collop_vth(op%sp(sa))/collop_vth(op%sp(sb))
        tau = op%sp(sa)%t/op%sp(sb)%t
        nu = op%sp(sa)%n*collop_nuhat(op%sp(sa), op%sp(sb), op%lnlambda(sa, sb))
        ! Relative tangents: rinv ~ sqrt(T_a/T_b), tau = T_a/T_b,
        ! n_a nuhat_ab ~ n_a n_b T_a^(-3/2).
        dtau = tau*(dt(sa)/op%sp(sa)%t - dt(sb)/op%sp(sb)%t)
        drinv = rinv*(dt(sa)/op%sp(sa)%t - dt(sb)/op%sp(sb)%t)/2
        dnu = nu*(dn(sa)/op%sp(sa)%n + dn(sb)/op%sp(sb)%n - 1.5_dp*dt(sa)/op%sp(sa)%t)

        call panel_rule(merge_breaks(op%grid%edges, op%grid%edges/rinv, &
            op%grid%edges(op%grid%nelem)), xq, wq)
        allocate (lorentz(n, n), energy(n, n))
        lorentz = 0
        energy = 0
        do q = 1, size(xq)
            if (xq(q)*rinv < 1) then
                call collop_tp_series(xq(q), rinv, tau, lor, en, th)
                call collop_tp_series_jvp(xq(q), rinv, tau, 0.0_dp, drinv, dtau, dlor, &
                    den, dth)
            else
                call collop_tp(xq(q), rinv, tau, lor, en, th)
                call collop_tp_jvp(xq(q), rinv, tau, 0.0_dp, drinv, dtau, dlor, den, dth)
            end if
            ! Product rule on nu times the weights.
            scale = nu
            if (tangent) then
                lor = dnu*lor + nu*dlor
                en = dnu*en + nu*den
                th = dnu*th + nu*dth
                scale = 1
            end if
            call collop_basis(op%grid, xq(q), first, phi, dphi)
            do j = 1, p
                do i = 1, p
                    lorentz(first + i - 1, first + j - 1) = lorentz(first + i - 1, first + j - 1) &
                        + wq(q)*scale*phi(i)*phi(j)*lor
                    if (op%model /= COLLOP_LORENTZ) &
                        energy(first + i - 1, first + j - 1) = &
                        energy(first + i - 1, first + j - 1) &
                        + wq(q)*scale*dphi(i)*(en*dphi(j) + th*phi(j))
                end do
            end do
        end do
        allocate (block(n, n, 0:op%lmax))
        do l = 0, op%lmax
            pref = 4*pi/(2*l + 1)/pi**1.5_dp
            block(:, :, l) = -pref*(l*(l + 1)/2.0_dp*lorentz + energy)
        end do
    end subroutine test_particle_block

    subroutine krook_block(op, sa, dn, tangent, block)
        type(collop_t), intent(in) :: op
        integer, intent(in) :: sa
        real(dp), intent(in) :: dn
        logical, intent(in) :: tangent
        real(dp), allocatable, intent(out) :: block(:, :, :)
        real(dp), allocatable :: xq(:), wq(:), m(:, :)
        real(dp) :: phi(op%grid%order), dphi(op%grid%order), dens
        integer :: q, i, j, first, l, p, n
        n = op%grid%n
        p = op%grid%order
        dens = op%sp(sa)%n
        if (tangent) dens = dn
        call panel_rule(op%grid%edges, xq, wq)
        allocate (m(n, n), block(n, n, 0:op%lmax))
        m = 0
        do q = 1, size(xq)
            call collop_basis(op%grid, xq(q), first, phi, dphi)
            do j = 1, p
                do i = 1, p
                    m(first + i - 1, first + j - 1) = m(first + i - 1, first + j - 1) &
                        + wq(q)*phi(i)*phi(j)*xq(q)**2*exp(-xq(q)**2)
                end do
            end do
        end do
        do l = 0, op%lmax
            block(:, :, l) = -op%nu_krook*4*pi/(2*l + 1)*dens/pi**1.5_dp*m
        end do
    end subroutine krook_block

    !> Sorted union of a and those points of b inside (0, top).
    function merge_breaks(a, b, top) result(c)
        real(dp), intent(in) :: a(:), b(:), top
        real(dp), allocatable :: c(:)
        real(dp), allocatable :: all(:)
        real(dp) :: t
        integer :: i, j, n
        all = [a, pack(b, b > 0 .and. b < top)]
        do i = 2, size(all)
            t = all(i)
            j = i - 1
            do while (j >= 1)
                if (all(j) <= t) exit
                all(j + 1) = all(j)
                j = j - 1
            end do
            all(j + 1) = t
        end do
        n = 1
        do i = 2, size(all)
            if (all(i) - all(n) > 1e-12_dp*top) then
                n = n + 1
                all(n) = all(i)
            end if
        end do
        c = all(:n)
    end function merge_breaks

    !> Composite Gauss-Legendre rule on the panels between sorted breaks.
    subroutine panel_rule(breaks, xq, wq)
        real(dp), intent(in) :: breaks(:)
        real(dp), allocatable, intent(out) :: xq(:), wq(:)
        real(dp) :: x0(nquad), w0(nquad), h
        integer :: k
        call gauss_legendre(nquad, x0, w0)
        allocate (xq(nquad*(size(breaks) - 1)), wq(nquad*(size(breaks) - 1)))
        do k = 1, size(breaks) - 1
            h = breaks(k + 1) - breaks(k)
            xq((k - 1)*nquad + 1:k*nquad) = breaks(k) + (x0 + 1)*h/2
            wq((k - 1)*nquad + 1:k*nquad) = w0*h/2
        end do
    end subroutine panel_rule

end module kin6d_dk_collop
