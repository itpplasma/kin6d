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
!! for unequal temperatures), full (test and field parts of every species
!! pair, conserving each species' density and total momentum and energy) and
!! Krook (-nu_K f, a parity switch, never a default). The test-particle operator
!! is symmetric for equal temperatures
!! and has the Maxwellian (h = 1) as null vector; it exchanges momentum and
!! energy with its partners, which are reservoirs here.
!!
!! Tangents with respect to the species temperatures and densities
!! (collop_build_jvp) come from the leaf JVPs at fixed normalised grid.
module milca_dk_collop
    use, intrinsic :: iso_fortran_env, only: dp => real64, real128
    use fortnum_quadrature, only: gauss_legendre, gauss_lobatto_legendre
    use fortnum_polynomial, only: lagrange_weights, lagrange_differentiation_matrix
    use milca_dk_collop_leaves, only: collop_tp, collop_tp_series, collop_tp_jvp, &
        collop_tp_series_jvp, collop_fp_lo, collop_fp_hi, collop_fp_lo_jvp, collop_fp_hi_jvp, &
        collop_fp_lo_moments, collop_fp_hi_moments
    implicit none
    private
    public :: collop_species_t, collop_grid_t, collop_t
    public :: collop_make_grid, collop_build, collop_build_jvp, collop_mass
    public :: collop_basis, collop_nuhat, collop_vth, collop_moments
    integer, parameter, public :: COLLOP_LORENTZ = 1, COLLOP_TEST_PARTICLE = 2, &
        COLLOP_FULL = 3, COLLOP_KROOK = 4
    integer, parameter, public :: COLLOP_OK = 0, COLLOP_BAD_INPUT = 1

    real(dp), parameter, public :: collop_qe = 1.602176634e-19_dp
    real(dp), parameter, public :: collop_eps0 = 8.8541878128e-12_dp
    real(dp), parameter :: pi = 4*atan(1.0_dp)
    ! Gauss-Legendre points per panel; integrands are polynomials times
    ! Gaussians and erf on panels no wider than one thermal speed.
    integer, parameter :: nquad = 32
    ! Unit panels integrating the polynomial tail beyond x_max: the
    ! integrals run over [0, oo), so that the collision invariants hold for
    ! any x_max, as in the continuous operator.
    integer, parameter :: tail = 6

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

    type :: potential_panel_t
        real(dp), allocatable :: y(:), weight(:), phi(:, :)
        integer, allocatable :: first(:)
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
    !> holding x, and the global index of the first. Beyond x_max the last
    !> element continues polynomially, so h is a polynomial on [x_max, oo).
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
        if (x < 0) return
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

    !> Radial moments over [0, oo) with the polynomial tail: the matrix
    !> int phi_i phi_j x^p exp(-x^2) dx and the vector int phi_i x^q
    !> exp(-x^2) dx, for the kinetic terms that a solver writes in the same
    !> Galerkin form as the collision operator.
    subroutine collop_moments(grid, p, q, matrix, vector)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: p, q
        real(dp), intent(out) :: matrix(grid%n, grid%n), vector(grid%n)
        real(dp), allocatable :: xq(:), wq(:)
        real(dp) :: phi(grid%order), dphi(grid%order), g
        integer :: k, i, j, first
        call panel_rule(tail_edges(grid), xq, wq)
        matrix = 0
        vector = 0
        do k = 1, size(xq)
            call collop_basis(grid, xq(k), first, phi, dphi)
            g = wq(k)*exp(-xq(k)**2)
            do j = 1, grid%order
                vector(first + j - 1) = vector(first + j - 1) + g*xq(k)**q*phi(j)
                do i = 1, grid%order
                    matrix(first + i - 1, first + j - 1) = matrix(first + i - 1, first + j - 1) &
                        + g*xq(k)**p*phi(i)*phi(j)
                end do
            end do
        end do
    end subroutine collop_moments

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
        if (all(abs(dt) <= 0) .and. all(abs(dn) <= 0)) then
            ! assemble uses zero directions to select the primal matrix.
            allocate (da(op%ns*grid%n, op%ns*grid%n, 0:lmax))
            da = 0
            status = COLLOP_OK
            return
        end if
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
        if (all(model /= [COLLOP_LORENTZ, COLLOP_TEST_PARTICLE, COLLOP_FULL, &
            COLLOP_KROOK])) return
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
                if (op%model /= COLLOP_FULL) cycle
                call field_particle_block(op, sa, sb, dt, dn, block)
                a(ra + 1:ra + n, rb + 1:rb + n, :) = a(ra + 1:ra + n, rb + 1:rb + n, :) + block
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
        real(dp), allocatable :: edges(:)
        ! Resolve self-collision cancellation before interspecies heat exchange.
        real(real128), allocatable :: lorentz(:, :), energy(:, :)
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

        edges = tail_edges(op%grid)
        call panel_rule(merge_breaks(edges, edges/rinv, edges(size(edges))), xq, wq)
        allocate (lorentz(n, n), energy(n, n))
        lorentz = 0
        energy = 0
        do q = 1, size(xq)
            if (xq(q)*rinv < 1) then
                call collop_tp_series(xq(q), rinv, tau, lor, en, th)
                if (tangent) call collop_tp_series_jvp(xq(q), rinv, tau, 0.0_dp, &
                    drinv, dtau, dlor, den, dth)
            else
                call collop_tp(xq(q), rinv, tau, lor, en, th)
                if (tangent) call collop_tp_jvp(xq(q), rinv, tau, 0.0_dp, &
                    drinv, dtau, dlor, den, dth)
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
            block(:, :, l) = real(-pref*(l*(l + 1)/2.0_dp*lorentz + energy), dp)
        end do
    end subroutine test_particle_block

    !> Field-particle block C_ab(F_a, f_b), f_b = phi_j F_b P_l, all l, in
    !> weak form. Its flux is F_a grad(Psi P_l), Psi from the Rosenbluth
    !> potentials of f_b, because (grad grad G).v = grad(v G' - G). In the
    !> speed x of a, with rho = v_tb/v_ta, mu = m_a/m_b and the test-particle
    !> prefactor pref,
    !>   pref (4 pi/(2l+1)) pi^(-3/2) int dx e^(-x^2) (x^2 phi_i' psi_j'
    !>        + l(l+1) phi_i psi_j),
    !>   psi_j(x) = int dy k1_l(x, rho y) phi_j(y) e^(-y^2) y^2.
    !> Density is conserved for every species identically (phi' = 0 for
    !> l = 0). The inner rule is split at y = x/rho, where the derivative of
    !> the kernel jumps, graded towards it, where the kernel (s/x)^l of high
    !> degree is concentrated, and geometrically beyond it, where the kernel
    !> decays like a power of y; below it the kernel is a polynomial in y.
    subroutine field_particle_block(op, sa, sb, dt, dn, block)
        type(collop_t), intent(in) :: op
        integer, intent(in) :: sa, sb
        real(dp), intent(in) :: dt(:), dn(:)
        real(dp), allocatable, intent(out) :: block(:, :, :)
        real(dp), allocatable :: xq(:), wq(:), yq(:), wy(:), ls(:), edges(:), origin(:)
        real(dp), allocatable :: k1(:), k1x(:), dk1(:), dk1x(:), psi(:, :), psix(:, :)
        real(dp) :: rho, mu, nu, drho, dnu, top, ystar, gauss, sig, scale, dsig
        real(dp) :: phi(op%grid%order), dphi(op%grid%order)
        real(dp) :: phj(op%grid%order), dphj(op%grid%order), c(0:op%lmax), cx(0:op%lmax)
        integer :: q, r, i, j, first, firstj, l, p, n, ngrade, k
        logical :: tangent

        n = op%grid%n
        p = op%grid%order
        edges = tail_edges(op%grid)
        top = edges(size(edges))
        tangent = any(dt /= 0) .or. any(dn /= 0)
        rho = collop_vth(op%sp(sb))/collop_vth(op%sp(sa))
        mu = op%sp(sa)%m/op%sp(sb)%m
        nu = op%sp(sa)%n*collop_nuhat(op%sp(sa), op%sp(sb), op%lnlambda(sa, sb))
        drho = rho*(dt(sb)/op%sp(sb)%t - dt(sa)/op%sp(sa)%t)/2
        dnu = nu*(dn(sa)/op%sp(sa)%n + dn(sb)/op%sp(sb)%n - 1.5_dp*dt(sa)/op%sp(sa)%t)
        if (abs(drho) <= 0) then
            call field_particle_sweep(op, rho, mu, merge(dnu, nu, tangent), block)
            return
        end if
        ls = [(real(l, dp), l=0, op%lmax)]
        ! For l >= 1 a basis function with h(0) /= 0 has potentials with
        ! x^l log(x) terms at the origin; geometric panels towards 0 keep
        ! the rule converging exponentially.
        origin = [(edges(2)*0.5_dp**k, k=1, 24)]
        ngrade = 2 + ceiling(log(real(op%lmax + 1, dp))/log(2.0_dp))
        allocate (block(n, n, 0:op%lmax), psi(n, 0:op%lmax), psix(n, 0:op%lmax), &
            k1(0:op%lmax), k1x(0:op%lmax), dk1(0:op%lmax), dk1x(0:op%lmax))
        block = 0

        call panel_rule(merge_breaks([edges, origin], [rho*edges, rho*origin], top), xq, wq)
        do q = 1, size(xq)
            ystar = xq(q)/rho
            call panel_rule(merge_breaks([edges, origin], [ystar, &
                (ystar*(1 - 0.5_dp**k), k=1, ngrade), (ystar*(1 + 0.5_dp**k), k=1, ngrade), &
                (ystar*2.0_dp**k, k=1, ceiling(log(top/ystar)/log(2.0_dp)))], top), yq, wy)
            psi = 0
            psix = 0
            do r = 1, size(yq)
                sig = rho*yq(r)
                dsig = yq(r)*drho
                if (sig < xq(q)) then
                    call collop_fp_lo(xq(q), sig/xq(q), ls, mu, k1, k1x)
                    call collop_fp_lo_jvp(xq(q), sig/xq(q), ls, mu, 0.0_dp, dsig/xq(q), &
                        0.0_dp, 0.0_dp, dk1, dk1x)
                else
                    call collop_fp_hi(sig, xq(q)/sig, ls, mu, k1, k1x)
                    call collop_fp_hi_jvp(sig, xq(q)/sig, ls, mu, dsig, -xq(q)*dsig/sig**2, &
                        0.0_dp, 0.0_dp, dk1, dk1x)
                end if
                ! Product rule on nu times the kernels for a tangent build.
                if (tangent) then
                    c = dnu*k1 + nu*dk1
                    cx = dnu*k1x + nu*dk1x
                else
                    c = nu*k1
                    cx = nu*k1x
                end if
                gauss = wy(r)*exp(-yq(r)**2)*yq(r)**2
                call collop_basis(op%grid, yq(r), firstj, phj, dphj)
                do j = 1, p
                    psi(firstj + j - 1, :) = psi(firstj + j - 1, :) + gauss*phj(j)*c
                    psix(firstj + j - 1, :) = psix(firstj + j - 1, :) + gauss*phj(j)*cx
                end do
            end do
            ! The rho-tangent of psi' has a boundary term: d_x k1 jumps where
            ! the split point y = x/rho moves.
            if (tangent .and. ystar < top) then
                call collop_fp_lo(xq(q), 1.0_dp, ls, mu, k1, k1x)
                call collop_fp_hi(xq(q), 1.0_dp, ls, mu, dk1, dk1x)
                call collop_basis(op%grid, ystar, firstj, phj, dphj)
                gauss = -nu*exp(-ystar**2)*ystar**2*ystar/rho*drho
                do j = 1, p
                    psix(firstj + j - 1, :) = psix(firstj + j - 1, :) &
                        + gauss*phj(j)*(k1x - dk1x)
                end do
            end if
            call collop_basis(op%grid, xq(q), first, phi, dphi)
            gauss = wq(q)*exp(-xq(q)**2)
            do l = 0, op%lmax
                scale = (4*pi/(2*l + 1))**2/pi**3
                do i = 1, p
                    block(first + i - 1, :, l) = block(first + i - 1, :, l) + scale*gauss &
                        *(xq(q)**2*dphi(i)*psix(:, l) + l*(l + 1)*phi(i)*psi(:, l))
                end do
            end do
        end do
    end subroutine field_particle_block

    !> Two scaled moment sweeps replace the repeated source integrals.
    !> Degrees are independent; the quadrature and basis tables are shared.
    subroutine field_particle_sweep(op, rho, mu, nu, block)
        type(collop_t), intent(in) :: op
        real(dp), intent(in) :: rho, mu, nu
        real(dp), allocatable, intent(out) :: block(:, :, :)
        type(potential_panel_t), allocatable :: panels(:)
        real(dp), allocatable :: edges(:), origin(:), xq(:), wq(:), phi(:, :), dphi(:, :)
        integer, allocatable :: first(:)
        real(dp) :: m0(op%grid%n), m2(op%grid%n), psi(op%grid%n), psix(op%grid%n)
        real(dp) :: oldx, oldy, y, ratio, power, scale, top, gauss
        integer :: n, p, nq, q, r, i, j, l, side, p0, p2, begin, finish, step, k, base
        logical :: lower

        n = op%grid%n
        p = op%grid%order
        edges = tail_edges(op%grid)
        top = edges(size(edges))
        origin = [(edges(2)*0.5_dp**k, k=1, 24)]
        call panel_rule(merge_breaks([edges, origin], [rho*edges, rho*origin], top), xq, wq)
        nq = size(xq)
        allocate (phi(p, nq), dphi(p, nq), first(nq), panels(nq), block(n, n, 0:op%lmax))
        do q = 1, nq
            call collop_basis(op%grid, xq(q), first(q), phi(:, q), dphi(:, q))
        end do
        block = 0
        do side = 1, 2
            lower = side == 1
            begin = merge(1, nq, lower)
            finish = merge(nq, 1, lower)
            step = merge(1, -1, lower)
            oldy = merge(0.0_dp, top, lower)
            do q = begin, finish, step
                y = min(xq(q)/rho, top)
                call potential_panel(op%grid, edges, min(oldy, y), max(oldy, y), &
                    y, lower, op%lmax, panels(q))
                oldy = y
            end do
            !$omp parallel do default(shared) schedule(static) &
            !$omp private(l,p0,p2,m0,m2,psi,psix,oldx,q,r,ratio,power,base,j,scale,gauss,i)
            do l = 0, op%lmax
                p0 = merge(l, l - 1, lower)
                p2 = p0 + 2
                m0 = 0
                m2 = 0
                oldx = xq(begin)
                scale = nu*(4*pi/(2*l + 1))**2/pi**3
                do q = begin, finish, step
                    ratio = merge(oldx/xq(q), xq(q)/oldx, lower)
                    m0 = m0*ratio**p0
                    m2 = m2*ratio**p2
                    do r = 1, size(panels(q)%y)
                        ratio = merge(rho*panels(q)%y(r)/xq(q), &
                            xq(q)/(rho*panels(q)%y(r)), lower)
                        power = ratio**p0
                        base = panels(q)%first(r) - 1
                        do j = 1, p
                            gauss = panels(q)%weight(r)*panels(q)%phi(j, r)*power
                            m0(base + j) = m0(base + j) + gauss
                            m2(base + j) = m2(base + j) + gauss*ratio**2
                        end do
                    end do
                    if (lower) then
                        call collop_fp_lo_moments(xq(q), real(l, dp), mu, m0, m2, psi, psix)
                    else
                        call collop_fp_hi_moments(xq(q), real(l, dp), mu, m0, m2, psi, psix)
                    end if
                    gauss = scale*wq(q)*exp(-xq(q)**2)
                    do i = 1, p
                        block(first(q) + i - 1, :, l) = block(first(q) + i - 1, :, l) &
                            + gauss*(xq(q)**2*dphi(i, q)*psix + l*(l + 1)*phi(i, q)*psi)
                    end do
                    oldx = xq(q)
                end do
            end do
            !$omp end parallel do
        end do
    end subroutine field_particle_sweep

    subroutine potential_panel(grid, edges, lo, hi, split, lower, lmax, panel)
        type(collop_grid_t), intent(in) :: grid
        real(dp), intent(in) :: edges(:), lo, hi, split
        logical, intent(in) :: lower
        integer, intent(in) :: lmax
        type(potential_panel_t), intent(out) :: panel
        real(dp), allocatable :: breaks(:), graded(:), w(:)
        real(dp) :: dphi(grid%order)
        integer :: k, ngrade
        ngrade = 2 + ceiling(log(real(lmax + 1, dp))/log(2.0_dp))
        graded = [(split*(1 + merge(-1.0_dp, 1.0_dp, lower)*0.5_dp**k), k=1, ngrade)]
        if (.not. lower) graded = [graded, (split*2.0_dp**k, &
            k=1, max(0, ceiling(log(hi/split)/log(2.0_dp))))]
        breaks = [pack(edges, edges > lo .and. edges < hi), &
            pack(graded, graded > lo .and. graded < hi)]
        call panel_rule(merge_breaks([lo, hi], breaks, hi), panel%y, w)
        allocate (panel%phi(grid%order, size(w)), panel%first(size(w)))
        panel%weight = w*exp(-panel%y**2)*panel%y**2
        do k = 1, size(w)
            call collop_basis(grid, panel%y(k), panel%first(k), panel%phi(:, k), dphi)
        end do
    end subroutine potential_panel

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
        call panel_rule(tail_edges(op%grid), xq, wq)
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

    !> Element edges continued by unit panels over the polynomial tail of
    !> the last element, to where exp(-x^2) has decayed below rounding.
    function tail_edges(grid) result(e)
        type(collop_grid_t), intent(in) :: grid
        real(dp), allocatable :: e(:)
        integer :: k
        e = [grid%edges, (grid%edges(grid%nelem) + k, k=1, tail)]
    end function tail_edges

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
    !> Panels wider than half a thermal speed are subdivided.
    subroutine panel_rule(breaks, xq, wq)
        real(dp), intent(in) :: breaks(:)
        real(dp), allocatable, intent(out) :: xq(:), wq(:)
        real(dp) :: x0(nquad), w0(nquad), h
        integer :: k, j, m, np, sub(size(breaks) - 1)
        call gauss_legendre(nquad, x0, w0)
        do k = 1, size(breaks) - 1
            sub(k) = max(1, ceiling(2*(breaks(k + 1) - breaks(k))))
        end do
        allocate (xq(nquad*sum(sub)), wq(nquad*sum(sub)))
        np = 0
        do k = 1, size(breaks) - 1
            h = (breaks(k + 1) - breaks(k))/sub(k)
            do j = 1, sub(k)
                m = np*nquad
                xq(m + 1:m + nquad) = breaks(k) + (j - 1)*h + (x0 + 1)*h/2
                wq(m + 1:m + nquad) = w0*h/2
                np = np + 1
            end do
        end do
    end subroutine panel_rule

end module milca_dk_collop
