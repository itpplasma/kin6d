!> Energy-coupled NTV kernel (procedures included in milca_dk_ntv) (#270): the operator of milca_dk_ntv on the
!! Galerkin speed basis of milca_dk_collop, with the linearised collision
!! operator in place of monoenergetic pitch-angle scattering.
!!
!! With f_a = F_a h_a and h_a = sum_k phi_k(x) h_ak(theta, xi), x = v/v_a,
!! each species obeys, tested with phi_k and the Legendre rows of dk_ntv,
!!   sum_j [v_a m3_kj S + m2_kj i n omega_E + (T_a/(Z_a e)) m4_kj P] h_aj
!!     - sum_b sum_j c_l(ak, bj) h_bj = (T_a/(Z_a e)) d_k (1 + xi^2) S_n/b,
!! m_p = int phi_k phi_j x^p e^-x^2 dx, S the streaming and mirror rows per
!! unit speed, P the precession rows per unit prefactor (both from dk_ntv at
!! v = 1), c_l = a_l (2l + 1)/(4 pi n_a pi^-3/2) from the Galerkin matrices
!! a_l of dk_collop, and d_k = int phi_k x^(2+2j) e^-x^2 dx for the drive of
!! A_j (j = 1, 2). Rows of node k are scaled by 1/(v_a m3_kk).
!!
!! Each species is factorised by pair-block elimination with its own
!! test-particle (or Lorentz) operator, blocks of 2 n_theta n_x unknowns.
!! For the full operator, whose field-particle part couples the species,
!! the coupled system is solved by GMRES right-preconditioned with these
!! factorisations.
!!
!! Gamma_psi of species a = -n_a sum_b sum_j k(a, b, j) A_j of species b,
!! k(a, b, j) = T_a/(Z_a e sqrt(pi) 2 J_0) sum_k d4_k R(h_ak) for the drive
!! of A_j of species b, R the flux functional of dk_ntv and J_0 = int
!! dtheta b^-2. The optional heat coefficients q(a, b, j) are the same with
!! d6_k in place of d4_k, so that the energy flux of species a is
!! -n_a T_a sum_b sum_j q(a, b, j) A_j of species b.

    subroutine ntv_energy_coefficients(sf, species, lnlambda, grid, model, omega_e, &
            n_theta, n_leg, k, iterations, status, q)
        type(ntv_surface_t), intent(in) :: sf
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :), omega_e
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_theta, n_leg
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out), optional :: q(:, :, :)
    type(energy_system_t) :: e, e0
    type(ntv_system_t) :: g, gf, t
    type(hybrid_t) :: lu, lu0
    complex(dp), allocatable :: angular(:, :), angular0(:, :), rhs(:, :, :, :), &
        x(:, :, :, :), h10(:, :, :, :)
    real(dp), allocatable :: b(:)
    real(dp) :: j0
    integer :: a, sb, jf, st, its, kk
    logical :: sources10

    status = NTV_BAD_INPUT
    iterations = 0
    if (size(k, 1) /= size(species) .or. size(k, 2) /= size(species) .or. size(k, 3) /= 2) &
        return
    if (present(q)) then
        if (any(shape(q) /= shape(k))) return
    end if
    sources10 = sf%source_dlf10
    call energy_build(sf, species, lnlambda, grid, model, omega_e, n_theta, n_leg, e, g, st)
    if (st /= NTV_OK) return
    status = NTV_SINGULAR
    call factor_hybrid(e, lu, st)
    if (st /= NTV_OK) return
    if (sources10) then
        ! The axisymmetric system: streaming and mirror only, bordered.
        e0 = e
        e0%gs%nq = 0
        e0%go%omega = 0
        e0%gp%precess = 0
        e0%border = .true.
        nullify (e0%cs, e0%co, e0%cp)
        call cache_blocks(e0)
        call factor_hybrid(e0, lu0, st)
        if (st /= NTV_OK) return
        ! f_10 drive (1 + xi^2) d_theta b/(q b): the drive of S = d_theta b/q.
        gf = g
        gf%dpar = .false.
        gf%dlf10 = .false.
        b = 2/(3*g%source)
        gf%s = log_b_slope(g)*b/sf%q
        call drive(gf, 1.0_dp, angular0)
    end if

    call drive(g, 1.0_dp, angular)
    j0 = 2*pi/g%n*sum((1.5_dp*g%source)**2)
    allocate (rhs(e%n*2, e%nx, 0:e%nb - 1, e%ns))
    do sb = 1, e%ns
        do jf = 1, 2
            rhs = 0
            do kk = 1, e%nx
                rhs(:, kk, :, sb) = e%row(kk, sb)*e%tz(sb) &
                    *merge(e%d4(kk), e%d6(kk), jf == 1)*angular
            end do
            if (sources10) then
                allocate (x, mold=rhs)
                x = 0
                do kk = 1, e%nx
                    x(:, kk, :, sb) = e%row(kk, sb)*e%tz(sb) &
                        *merge(e%d4(kk), e%d6(kk), jf == 1)*angular0
                    x(1:e%n, kk, 0, sb) = x(1:e%n, kk, 0, sb) &
                        - sum(x(1:e%n, kk, 0, sb))/e%n
                end do
                call solve_any(e0, lu0, x, h10, its, st)
                if (st /= NTV_OK) return
                iterations = iterations + its
                deallocate (x)
                t = g
                call add_f10_sources(e, t, h10, rhs)
            end if
            call solve_any(e, lu, rhs, x, its, st)
            if (st /= NTV_OK) return
            iterations = iterations + its
            do a = 1, e%ns
                k(a, sb, jf) = 0
                do kk = 1, e%nx
                    k(a, sb, jf) = k(a, sb, jf) + e%d4(kk)*flux(g, x(:, kk, :, a))
                end do
                k(a, sb, jf) = e%tz(a)/(sqrt(pi)*2*j0)*k(a, sb, jf)
                if (present(q)) then
                    q(a, sb, jf) = 0
                    do kk = 1, e%nx
                        q(a, sb, jf) = q(a, sb, jf) + e%d6(kk)*flux(g, x(:, kk, :, a))
                    end do
                    q(a, sb, jf) = e%tz(a)/(sqrt(pi)*2*j0)*q(a, sb, jf)
                end if
            end do
            deallocate (x)
        end do
    end do
    status = NTV_OK
    end subroutine ntv_energy_coefficients

    !> Axisymmetric transport matrix of the surface with the energy-coupled
    !> operator (n = 0: streaming, mirror and collisions, bordered by the
    !> density null space). Drives of species b, at each speed node:
    !>   j = 1, 2  the radial drift (1 + xi^2) S/b with S = b'/q (the f_10 drive)
    !>             times (T_b/(Z_b e)) d_k, d = d4, d6 for A_1, A_2;
    !>   j = 3     the parallel electric drive xi b times v_b d3_k.
    !> Responses of species a, i = 1, 2, 3: the radial particle and energy
    !> coefficients in the units of the NTV coefficients k and q,
    !>   T_a/(Z_a e sqrt(pi) 2 J_0) sum_k d4_k (d6_k) R(h_ak),
    !> and the parallel flow v_a/(sqrt(pi) 2 J_0) sum_k d3_k U(h_ak), U the
    !> functional of xi b with the weight of R, so that the matrix is
    !> symmetric under d(a, b, 1, 3) = -d(b, a, 3, 1) for a self-adjoint
    !> operator. U is proportional to <B u_par>: its weight is
    !> measure b^2 = 2/b, the Boozer surface weight 1/b^2 times b. kappa is
    !> b.grad theta/b; with pitch-angle scattering each node then solves the
    !> problem of dk_mono driven by (T/(Z e)) x (-2/(kappa q)) times its drift,
    !> so (3/4) (T/(Z e)) (-2/(kappa q)) relates the flow response to A_1 to
    !> the energy average of the monoenergetic D31. local(theta, a, b, j), on
    !> the n_theta collocation points theta_i = 2 pi (i - 1)/n_theta, is the
    !> parallel-flow response at each point in the same units, the quantity
    !> whose surface average with the Boozer weight 1/b^2 is d(a, b, 3, j):
    !> u_par B per unit drive.
    subroutine ntv_axisymmetric_coefficients(sf, species, lnlambda, grid, model, n_theta, &
            n_leg, d, iterations, status, kappa, local)
        type(ntv_surface_t), intent(in) :: sf
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_theta, n_leg
        real(dp), intent(out) :: d(:, :, :, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out), optional :: kappa
        real(dp), intent(out), optional :: local(:, :, :, :)
    type(energy_system_t) :: e
    type(ntv_system_t) :: g
    integer :: st

    status = NTV_BAD_INPUT
    iterations = 0
    if (any(shape(d) /= [size(species), size(species), 3, 3])) return
    if (present(local)) then
        if (any(shape(local) /= [n_theta, size(species), size(species), 3])) return
    end if
    call energy_build(sf, species, lnlambda, grid, model, 0.0_dp, n_theta, n_leg, e, g, st)
    if (st /= NTV_OK) return
    if (present(kappa)) kappa = 1.5_dp*g%stream(1)*g%source(1)
    call transport_matrix(e, g, &
        cmplx(log_b_slope(g)/(1.5_dp*g%source*sf%q), 0.0_dp, dp), e%nb, d, iterations, &
        status, local)
    end subroutine ntv_axisymmetric_coefficients

    !> The transport matrix of ntv_axisymmetric_coefficients on a 3D surface in
    !> Boozer angles with zero E x B: b(i, j) = B/B0 on n_theta x n_zeta (both
    !> odd) points of theta in [0, 2 pi) and zeta in one field period
    !> [0, 2 pi/nper), b.grad = kappa b (d_theta + q d_zeta) with kappa =
    !> iota B0/(G + iota I) [1/m], and the radial drift of the axisymmetric
    !> surface with d_theta b replaced by d_theta b - (I/G) d_zeta b
    !> (geometry of mono_coefficients_3d). Same drives, responses and units as
    !> the axisymmetric matrix, which an axisymmetric b reproduces; local is
    !> (n_theta n_zeta, ns, ns, 3), theta fastest. pairs (default all) is the
    !> number of lowest pair blocks factorised exactly, with the full speed and
    !> species coupling; with fewer, GMRES on the operator is right-
    !> preconditioned by them and, above them, by one direct block
    !> factorisation per collision mode (factor_modal), whose memory is that
    !> of the monoenergetic surface problem. pairs = 1 makes l <= 1, with
    !> density, momentum and energy conservation, exact.
    subroutine ntv_surface_coefficients(b, nper, q, i_over_g, kappa, species, &
            lnlambda, grid, model, n_leg, d, iterations, status, local, pairs)
        real(dp), intent(in) :: b(:, :), q, i_over_g, kappa
        integer, intent(in) :: nper
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_leg
        real(dp), intent(out) :: d(:, :, :, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out), optional :: local(:, :, :, :)
        integer, intent(in), optional :: pairs
    type(energy_system_t) :: e
    type(ntv_system_t) :: g
    real(dp), allocatable :: b_par(:), b_drift(:), radial(:)
    integer :: n, st, np

    status = NTV_BAD_INPUT
    iterations = 0
    n = size(b)
    if (any(shape(d) /= [size(species), size(species), 3, 3])) return
    if (mod(n_leg, 2) /= 0 .or. n_leg < 4 .or. abs(q) <= tiny(q) .or. &
        abs(kappa) <= tiny(kappa)) return
    if (present(local)) then
        if (any(shape(local) /= [n, size(species), size(species), 3])) return
    end if
    call mono_surface_operators(b, nper, q, i_over_g, g%d, b_par, b_drift, st)
    if (st /= MONO_OK) return
    g%n = n
    g%nblock = n_leg/2
    g%wsign = 1
    allocate (g%stream(n), g%mirror(n), g%source(n), g%measure(n), radial(n))
    call dk_ntv_surface_point(reshape(b, [n]), b_par, b_drift, kappa, q, g%stream, &
        g%mirror, radial, g%source, g%measure)
    allocate (g%drift(n), g%drift2(n), g%s(n))
    g%drift = 0
    g%drift2 = 0
    g%s = radial
    call energy_assemble(g, species, lnlambda, grid, model, 0.0_dp, e, st)
    if (st /= NTV_OK) return
    np = e%nb
    if (present(pairs)) np = max(0, min(pairs, e%nb))
    call transport_matrix(e, g, cmplx(radial, 0.0_dp, dp), np, d, iterations, status, &
        local)
    end subroutine ntv_surface_coefficients

    !> Shared n = 0 solve of the transport matrices: e and g assembled, radial
    !> the radial drive S on the surface points.
    subroutine transport_matrix(e, g, radial_s, pairs, d, iterations, status, local)
        type(energy_system_t), intent(inout) :: e
        type(ntv_system_t), intent(in) :: g
        complex(dp), intent(in) :: radial_s(:)
        integer, intent(in) :: pairs
        real(dp), intent(out) :: d(:, :, :, :)
        integer, intent(out) :: iterations, status
        real(dp), intent(out), optional :: local(:, :, :, :)
    type(ntv_system_t) :: gf
    type(hybrid_t) :: lu
    complex(dp), allocatable :: radial(:, :), parallel(:, :), rhs(:, :, :, :), &
        x(:, :, :, :)
    real(dp), allocatable :: b(:), d3(:)
    real(dp) :: j0
    integer :: a, sb, jf, st, kk, its

    status = NTV_SINGULAR
    iterations = 0
    e%gs%nq = 0
    e%go%omega = 0
    e%gp%precess = 0
    e%border = .true.
    nullify (e%cs, e%co, e%cp)
    call cache_blocks(e)
    ! Without rotation and precession only collisions damp the parallel-field
    ! response (~1/nu), and GMRES on the hybrid preconditioner stalls in the
    ! banana regime. By default every pair block therefore carries its full
    ! speed and species coupling, which makes the block factorisation the
    ! inverse.
    call factor_hybrid(e, lu, st, pairs, modal=pairs < e%nb)
    if (st /= NTV_OK) return
    b = 2/(3*g%source)
    gf = g
    gf%dpar = .false.
    gf%dlf10 = .false.
    gf%s = radial_s
    call drive(gf, 1.0_dp, radial)
    allocate (parallel, mold=radial)
    parallel = 0
    parallel(e%n + 1:2*e%n, 0) = b
    j0 = 2*pi/g%n*sum((1.5_dp*g%source)**2)
    d3 = diag(e%m3)
    allocate (rhs(e%n*2, e%nx, 0:e%nb - 1, e%ns))
    do sb = 1, e%ns
        do jf = 1, 3
            rhs = 0
            do kk = 1, e%nx
                if (jf < 3) then
                    rhs(:, kk, :, sb) = e%row(kk, sb)*e%tz(sb) &
                        *merge(e%d4(kk), e%d6(kk), jf == 1)*radial
                else
                    rhs(:, kk, :, sb) = e%row(kk, sb)*e%vt(sb)*d3(kk)*parallel
                end if
                rhs(1:e%n, kk, 0, sb) = rhs(1:e%n, kk, 0, sb) - sum(rhs(1:e%n, kk, 0, sb))/e%n
            end do
            if (pairs < e%nb) then
                call coupled_solve(e, lu, rhs, x, its, st)
                iterations = iterations + its
            else
                if (.not. allocated(x)) allocate (x, mold=rhs)
                call hybrid_solve(e, lu, rhs, x, st)
            end if
            if (st /= NTV_OK) then
                status = st
                return
            end if
            do a = 1, e%ns
                d(a, sb, :, jf) = 0
                do kk = 1, e%nx
                    d(a, sb, 1, jf) = d(a, sb, 1, jf) + e%d4(kk)*flux(gf, x(:, kk, :, a))
                    d(a, sb, 2, jf) = d(a, sb, 2, jf) + e%d6(kk)*flux(gf, x(:, kk, :, a))
                    d(a, sb, 3, jf) = d(a, sb, 3, jf) + d3(kk)*flow(x(:, kk, :, a))
                end do
                d(a, sb, 1:2, jf) = e%tz(a)/(sqrt(pi)*2*j0)*d(a, sb, 1:2, jf)
                d(a, sb, 3, jf) = e%vt(a)/(sqrt(pi)*2*j0)*d(a, sb, 3, jf)
                if (present(local)) then
                    local(:, a, sb, jf) = 0
                    do kk = 1, e%nx
                        local(:, a, sb, jf) = local(:, a, sb, jf) &
                            + d3(kk)*real(x(e%n + 1:2*e%n, kk, 0, a), dp)
                    end do
                    local(:, a, sb, jf) = e%vt(a)/sqrt(pi)*b*2/3.0_dp*local(:, a, sb, jf)
                end if
            end do
        end do
    end do
    status = NTV_OK
    contains
        !> 2 pi/n sum measure b^2 int xi h dxi, int xi P_1 dxi = 2/3.
        real(dp) function flow(h)
            complex(dp), intent(in) :: h(:, 0:)
            flow = 2*pi/g%n*sum(g%measure*b**2*2/3.0_dp*real(h(e%n + 1:2*e%n, 0), dp))
        end function flow

        pure function diag(m) result(v)
            real(dp), intent(in) :: m(:, :)
            real(dp) :: v(size(m, 1))
            integer :: i
            v = [(m(i, i), i=1, size(m, 1))]
        end function diag
    end subroutine transport_matrix

    !> The energy-coupled system of one surface, assembled but not factorised:
    !> streaming, rotation and precession blocks, the lumped kinetic masses,
    !> the drive moments d4 and d6, and the collision rows (test-particle
    !> preconditioner rows in pre, the full operator in full when coupled).
    subroutine energy_build(sf, species, lnlambda, grid, model, omega_e, n_theta, n_leg, &
            e, g, status)
        type(ntv_surface_t), intent(in) :: sf
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :), omega_e
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_theta, n_leg
        type(energy_system_t), intent(out) :: e
        type(ntv_system_t), intent(out) :: g
        integer, intent(out) :: status
    type(ntv_surface_t) :: base
    real(dp), allocatable :: zero10(:, :)
    integer :: st

    status = NTV_BAD_INPUT
    base = sf
    allocate (zero10(2*n_theta, n_leg/2))
    zero10 = 0
    call factor(base, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, n_theta, n_leg, g, st, &
        f10=zero10, assemble_only=.true.)
    if (st /= NTV_OK) return
    call energy_assemble(g, species, lnlambda, grid, model, omega_e, e, status)
    end subroutine energy_build

    !> The energy-coupled system of the assembled surface rows g.
    subroutine energy_assemble(g, species, lnlambda, grid, model, omega_e, e, status)
        type(ntv_system_t), intent(in) :: g
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :), omega_e
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model
        type(energy_system_t), intent(out) :: e
        integer, intent(out) :: status
    integer :: n_leg

    status = NTV_BAD_INPUT
    e%ns = size(species)
    ! Nodes 2..n: h = 0 at x = 0, exact for l >= 1 by regularity and of
    ! zero weight in every drive and flux for l = 0.
    e%nx = grid%n - 1
    if (e%ns < 1 .or. size(lnlambda, 1) /= e%ns .or. size(lnlambda, 2) /= e%ns) return
    n_leg = 2*g%nblock
    e%n = g%n
    e%nb = g%nblock
    e%nn = 2*g%n*e%nx
    e%gs = g
    e%gs%nu = 0; e%gs%omega = 0; e%gs%precess = 0
    e%go = g
    e%go%stream = 0; e%go%mirror = 0; e%go%nu = 0; e%go%precess = 0
    e%go%omega = omega_e
    e%gp = g
    e%gp%stream = 0; e%gp%mirror = 0; e%gp%nu = 0; e%gp%omega = 0
    call cache_blocks(e)
    call collision_part(species, lnlambda, grid, model, n_leg, e, status)
    end subroutine energy_assemble

    !> Speed masses, drive moments and the scaled collision rows of e (with
    !> e%ns and e%nx set): the part of the energy-coupled system that does
    !> not depend on the surface.
    subroutine collision_part(species, lnlambda, grid, model, n_leg, e, status)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_leg
        type(energy_system_t), intent(inout) :: e
        integer, intent(out) :: status
    type(collop_t) :: op
    real(dp), allocatable :: galerkin(:, :), moment(:), lorentz(:, :, :, :, :)
    integer :: a, st, kk

    allocate (e%d4(e%nx), e%d6(e%nx))
    ! Kinetic terms with the speed masses lumped at the GLL nodes (the
    ! spectral-element mass), so that only the collision rows couple nodes.
    allocate (galerkin(grid%n, grid%n), moment(grid%n))
    call collop_moments(grid, 3, 6, galerkin, moment)
    allocate (e%vt(e%ns), e%tz(e%ns), e%row(e%nx, e%ns), e%m2(e%nx, e%nx), &
        e%m3(e%nx, e%nx), e%m4(e%nx, e%nx))
    do a = 1, e%ns
        e%vt(a) = collop_vth(species(a))
        e%tz(a) = species(a)%t/(species(a)%z*charge)
        do kk = 1, e%nx
            e%row(kk, a) = 1/(e%vt(a)*galerkin(kk + 1, kk + 1))
        end do
    end do
    e%m2 = 0; e%m3 = 0; e%m4 = 0
    do kk = 1, e%nx
        associate (x => grid%x(kk + 1), w => grid%w(kk + 1))
            e%m2(kk, kk) = w*x**2*exp(-x**2)
            e%m3(kk, kk) = w*x**3*exp(-x**2)
            e%m4(kk, kk) = w*x**4*exp(-x**2)
        end associate
    end do
    e%d4 = grid%w(2:)*grid%x(2:)**4*exp(-grid%x(2:)**2)
    e%d6 = grid%w(2:)*grid%x(2:)**6*exp(-grid%x(2:)**2)

    status = NTV_SINGULAR
    e%coupled = model == COLLOP_FULL
    call collop_build(species, lnlambda, grid, merge(COLLOP_TEST_PARTICLE, model, &
        e%coupled), n_leg - 1, op, st)
    if (st /= COLLOP_OK) return
    call scaled_rows(e, op, species, e%pre)
    if (e%coupled) then
        call collop_build(species, lnlambda, grid, model, n_leg - 1, op, st)
        if (st /= COLLOP_OK) return
        call scaled_rows(e, op, species, e%full)
    end if
    ! Pitch-angle scattering lumped like the kinetic masses: nu_D at the
    ! nodes in place of its Galerkin matrix.
    call collop_build(species, lnlambda, grid, COLLOP_LORENTZ, n_leg - 1, op, st)
    if (st /= COLLOP_OK) return
    call scaled_rows(e, op, species, lorentz)
    call lump_lorentz(e, species, lnlambda, grid, lorentz)
    e%pre = e%pre - lorentz
    if (e%coupled) e%full = e%full - lorentz
    status = NTV_OK
    end subroutine collision_part

    !> Collision rows of the energy-coupled operator for solvers outside this
    !> module: rows(k, kk, a, b, l), node k of species a against node kk of
    !> species b in degree l, scaled by the row factor of (k, a), the operator's
    !> own (full or test-particle) rows; ratio(k, a) the streaming coefficient
    !> of node k after the same scaling (lumped over Galerkin speed mass); d3,
    !> d4, d6 the lumped moments w x^p exp(-x^2) of the nodes 2..n; vt, tz the
    !> thermal speed and T/(Z e) of each species.
    subroutine ntv_collision_rows(species, lnlambda, grid, model, n_leg, rows, ratio, &
            d3, d4, d6, vt, tz, status)
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: model, n_leg
        real(dp), allocatable, intent(out) :: rows(:, :, :, :, :), ratio(:, :), d3(:), &
            d4(:), d6(:), vt(:), tz(:)
        integer, intent(out) :: status
    type(energy_system_t) :: e
    integer :: a, b, k, l
    status = NTV_BAD_INPUT
    e%ns = size(species)
    e%nx = grid%n - 1
    if (e%ns < 1 .or. size(lnlambda, 1) /= e%ns .or. size(lnlambda, 2) /= e%ns) return
    if (n_leg < 2) return
    call collision_part(species, lnlambda, grid, model, n_leg, e, status)
    if (status /= NTV_OK) return
    allocate (rows(e%nx, e%nx, e%ns, e%ns, 0:n_leg - 1), ratio(e%nx, e%ns))
    do l = 0, n_leg - 1
        do b = 1, e%ns
            do a = 1, e%ns
                do k = 1, e%nx
                    rows(k, :, a, b, l) = collision_rows(e, k, a, b, l)
                end do
            end do
        end do
    end do
    do a = 1, e%ns
        do k = 1, e%nx
            ratio(k, a) = e%row(k, a)*e%vt(a)*e%m3(k, k)
        end do
    end do
    d3 = [(e%m3(k, k), k=1, e%nx)]
    d4 = e%d4
    d6 = e%d6
    vt = e%vt
    tz = e%tz
    end subroutine ntv_collision_rows

    !> The system of one speed node: the node-diagonal part of the operator,
    !> including the diagonal of the energy scattering.
    function node_system(e, k) result(f)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: k
        type(energy_system_t) :: f
        f%ns = e%ns; f%n = e%n; f%nb = e%nb; f%border = e%border
        f%vt = e%vt; f%tz = e%tz
        f%cs => e%cs; f%co => e%co; f%cp => e%cp
        f%nx = 1
        f%nn = 2*e%n
        f%m2 = e%m2(k:k, k:k)
        f%m3 = e%m3(k:k, k:k)
        f%m4 = e%m4(k:k, k:k)
        f%row = e%row(k:k, :)
        ! Explicit bounds: a section would start the degree index at 1, not 0.
        allocate (f%pre(1, 1, e%ns, e%ns, 0:ubound(e%pre, 5)))
        f%pre(1, 1, :, :, :) = e%pre(k, k, :, :, :)
        f%coupled = .false.
    end function node_system

    !> Preconditioner: one pair-block factorisation per species and node.
    subroutine factor_nodes(e, sys, lu, status)
        type(energy_system_t), intent(in) :: e
        type(energy_system_t), allocatable, intent(out) :: sys(:)
        type(species_lu_t), allocatable, intent(out) :: lu(:, :)
        integer, intent(out) :: status
        integer :: a, k
        allocate (sys(e%nx), lu(e%nx, e%ns))
        do k = 1, e%nx
            sys(k) = node_system(e, k)
        end do
        do a = 1, e%ns
            do k = 1, e%nx
                call factor_species(sys(k), a, lu(k, a), status)
                if (status /= NTV_OK) return
            end do
        end do
    end subroutine factor_nodes

    !> Dense block (j, jc) of all species and nodes with the operator's own
    !> collision rows, unknowns ordered (theta row, node, species).
    function full_block(e, j, jc) result(m)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: j, jc
        complex(dp) :: m(e%nn*e%ns, e%nn*e%ns)
        complex(dp) :: corr(e%nn*e%ns, e%nn*e%ns)
        integer :: a, b, k, kk, r, l, i, n2, row, col
        n2 = 2*e%n
        m = 0
        do a = 1, e%ns
            m((a - 1)*e%nn + 1:a*e%nn, (a - 1)*e%nn + 1:a*e%nn) = energy_block(e, a, j, jc)
        end do
        if (j /= jc) return
        corr = 0
        do a = 1, e%ns
            do b = 1, e%ns
                do r = 0, 1
                    l = 2*j + r
                    do k = 1, e%nx
                        do kk = 1, e%nx
                            do i = 1, e%n
                                row = (a - 1)*e%nn + (k - 1)*n2 + r*e%n + i
                                col = (b - 1)*e%nn + (kk - 1)*n2 + r*e%n + i
                                corr(row, col) = merge(e%pre(k, kk, a, a, l), 0.0_dp, &
                                    a == b) - collision_rows_entry(e, k, kk, a, b, l)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        if (e%border .and. j == 0) then
            do a = 1, e%ns
                do k = 1, e%nx
                    row = (a - 1)*e%nn + (k - 1)*n2
                    associate (rows => corr(row + 1:row + e%n, :))
                        rows = rows - spread(sum(rows, 1)/e%n, 1, e%n)
                    end associate
                end do
            end do
        end if
        m = m + corr
    end function full_block

    real(dp) function collision_rows_entry(e, k, kk, a, b, l) result(c)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: k, kk, a, b, l
        real(dp) :: rows(e%nx)
        rows = collision_rows(e, k, a, b, l)
        c = rows(kk)
    end function collision_rows_entry

    !> pairs: the number of lowest pair blocks factorised densely with the
    !> full speed coupling (default low_pairs); the others per node.
    subroutine factor_hybrid(e, h, status, pairs, modal)
        type(energy_system_t), intent(inout) :: e
        type(hybrid_t), intent(out) :: h
        integer, intent(out) :: status
        integer, intent(in), optional :: pairs
        ! Chain above the dense blocks per collision mode instead of per node.
        logical, intent(in), optional :: modal
        complex(dp), allocatable :: m(:, :)
        integer :: j, a, k, info, n2, off
        h%modal = .false.
        if (present(modal)) h%modal = modal
        if (h%modal) then
            call factor_modal(e, h, status)
        else
            call factor_nodes(e, h%sys, h%node, status)
        end if
        if (status /= NTV_OK) return
        h%j0 = min(low_pairs, e%nb)
        if (present(pairs)) h%j0 = min(pairs, e%nb)
        h%nf = e%nn*e%ns
        n2 = 2*e%n
        allocate (h%s_lu(h%nf, h%nf, 0:h%j0 - 1), h%y(h%nf, h%nf, h%j0 - 1), &
            h%piv(h%nf, 0:h%j0 - 1), m(h%nf, h%nf))
        status = NTV_SINGULAR
        do j = h%j0 - 1, 0, -1
            m = full_block(e, j, j)
            if (j == h%j0 - 1 .and. h%j0 < e%nb .and. h%modal) then
                call subtract_modal_coupling(e, h, j, m)
            else if (j == h%j0 - 1 .and. h%j0 < e%nb) then
                do a = 1, e%ns
                    do k = 1, e%nx
                        off = (a - 1)*e%nn + (k - 1)*n2
                        m(off + 1:off + n2, off + 1:off + n2) = &
                            m(off + 1:off + n2, off + 1:off + n2) &
                            - matmul(energy_block(h%sys(k), a, j, j + 1), &
                            h%node(k, a)%y(:, :, j + 1))
                    end do
                end do
            else if (j < h%j0 - 1) then
                call zgemm('N', 'N', h%nf, h%nf, h%nf, (-1.0_dp, 0.0_dp), &
                    full_block(e, j, j + 1), h%nf, h%y(:, :, j + 1), h%nf, &
                    (1.0_dp, 0.0_dp), m, h%nf)
            end if
            call zgetrf(h%nf, h%nf, m, h%nf, h%piv(:, j), info)
            if (info /= 0) return
            h%s_lu(:, :, j) = m
            if (j > 0) then
                h%y(:, :, j) = full_block(e, j, j - 1)
                call zgetrs('N', h%nf, h%nf, m, h%nf, h%piv(:, j), h%y(:, :, j), h%nf, info)
            end if
        end do
        status = NTV_OK
    end subroutine factor_hybrid

    !> x = P^-1 rhs with the hybrid preconditioner; arrays (2 n, nx, 0:nb-1, ns).
    subroutine hybrid_solve(e, h, rhs, x, status)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(in) :: h
        complex(dp), intent(in) :: rhs(:, :, 0:, :)
        complex(dp), intent(out) :: x(:, :, 0:, :)
        integer, intent(out) :: status
        complex(dp), allocatable :: z(:, :, :, :), t(:, :, :)
        integer :: a, k, j, info
        if (h%modal) then
            call modal_hybrid_solve(e, h, rhs, x, status)
            return
        end if
        allocate (z, mold=x)
        allocate (t(2*e%n, e%nx, e%ns))
        status = NTV_SINGULAR
        do a = 1, e%ns
            do k = 1, e%nx
                do j = e%nb - 1, h%j0, -1
                    z(:, k, j, a) = rhs(:, k, j, a)
                    if (j < e%nb - 1) z(:, k:k, j, a) = z(:, k:k, j, a) &
                        - offdiag_apply(h%sys(k), a, j, j + 1, z(:, k:k, j + 1, a))
                    call zgetrs('N', 2*e%n, 1, h%node(k, a)%s_lu(:, :, j), 2*e%n, &
                        h%node(k, a)%piv(:, j), z(:, k, j, a), 2*e%n, info)
                    if (info /= 0) return
                end do
            end do
        end do
        do j = h%j0 - 1, 0, -1
            t = rhs(:, :, j, :)
            if (j == h%j0 - 1 .and. h%j0 < e%nb) then
                do a = 1, e%ns
                    do k = 1, e%nx
                        t(:, k:k, a) = t(:, k:k, a) &
                            - offdiag_apply(h%sys(k), a, j, j + 1, z(:, k:k, j + 1, a))
                    end do
                end do
            else if (j < h%j0 - 1) then
                do a = 1, e%ns
                    t(:, :, a) = t(:, :, a) - offdiag_apply(e, a, j, j + 1, z(:, :, j + 1, a))
                end do
            end if
            call zgetrs('N', h%nf, 1, h%s_lu(:, :, j), h%nf, h%piv(:, j), t, h%nf, info)
            if (info /= 0) return
            z(:, :, j, :) = t
        end do
        x(:, :, 0, :) = z(:, :, 0, :)
        do j = 1, e%nb - 1
            if (j < h%j0) then
                x(:, :, j, :) = z(:, :, j, :) - reshape(matmul(h%y(:, :, j), &
                    reshape(x(:, :, j - 1, :), [h%nf])), [2*e%n, e%nx, e%ns])
            else
                do a = 1, e%ns
                    do k = 1, e%nx
                        x(:, k, j, a) = z(:, k, j, a) - matmul(h%node(k, a)%y(:, :, j), &
                            x(:, k, j - 1, a))
                    end do
                end do
            end if
        end do
        status = NTV_OK
    end subroutine hybrid_solve

    !> hybrid_solve with the chain above the dense blocks per collision mode:
    !> modal coordinates V^-1 R^-1 rhs for pair blocks j >= j0, the dense
    !> blocks j < j0 in speed and species, joined through R V at the boundary.
    subroutine modal_hybrid_solve(e, h, rhs, x, status)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(in) :: h
        complex(dp), intent(in) :: rhs(:, :, 0:, :)
        complex(dp), intent(out) :: x(:, :, 0:, :)
        integer, intent(out) :: status
        complex(dp), allocatable :: zm(:, :, :), ym(:, :, :), u(:, :), t(:, :), xd(:, :, :)
        integer :: m, j, info, nm, n2, bad, j0
        nm = e%nx*e%ns
        n2 = 2*e%n
        j0 = h%j0
        allocate (zm(n2, j0:e%nb - 1, nm), ym(n2, j0:e%nb - 1, nm), u(n2, nm), &
            t(n2, nm), xd(n2, nm, 0:max(j0 - 1, 0)))
        do j = j0, e%nb - 1
            zm(:, j, :) = matmul(reshape(rhs(:, :, j, :), [n2, nm]), transpose(h%vinv))
        end do
        bad = 0
        !$omp parallel do schedule(dynamic, 1) default(shared) private(m, j, info) &
        !$omp reduction(+:bad)
        do m = 1, nm
            associate (f => h%node(m, 1), s => h%sys(m))
                do j = e%nb - 1, j0, -1
                    if (j < e%nb - 1) zm(:, j:j, m) = zm(:, j:j, m) &
                        - offdiag_apply(s, 1, j, j + 1, zm(:, j + 1:j + 1, m))
                    call zgetrs('N', n2, 1, f%s_lu(:, :, j), n2, f%piv(:, j), &
                        zm(:, j, m), n2, info)
                    if (info /= 0) bad = bad + 1
                end do
            end associate
        end do
        !$omp end parallel do
        status = NTV_SINGULAR
        if (bad > 0) return
        ! Dense blocks j < j0, the boundary coupled through (I x R V).
        do j = j0 - 1, 0, -1
            t = reshape(rhs(:, :, j, :), [n2, nm])
            if (j == j0 - 1 .and. j0 < e%nb) then
                do m = 1, nm
                    u(:, m:m) = offdiag_apply(h%sys(m), 1, j, j + 1, zm(:, j + 1:j + 1, m))
                end do
                t = t - matmul(u, transpose(h%rv))
            else if (j < j0 - 1) then
                do m = 1, e%ns
                    t(:, (m - 1)*e%nx + 1:m*e%nx) = t(:, (m - 1)*e%nx + 1:m*e%nx) &
                        - offdiag_apply(e, m, j, j + 1, xd(:, (m - 1)*e%nx + 1:m*e%nx, j + 1))
                end do
            end if
            call zgetrs('N', h%nf, 1, h%s_lu(:, :, j), h%nf, h%piv(:, j), t, h%nf, info)
            if (info /= 0) return
            xd(:, :, j) = t
        end do
        do j = 1, j0 - 1
            xd(:, :, j) = xd(:, :, j) - reshape(matmul(h%y(:, :, j), &
                reshape(xd(:, :, j - 1), [h%nf])), [n2, nm])
        end do
        ! Modal back substitution; at j0 the dense x_{j0-1} in modal coordinates.
        if (j0 > 0) then
            u = matmul(xd(:, :, j0 - 1), transpose(h%vu))
            do m = 1, nm
                ym(:, j0, m) = zm(:, j0, m) - matmul(h%node(m, 1)%y(:, :, j0), u(:, m))
            end do
        else
            ym(:, 0, :) = zm(:, 0, :)
        end if
        !$omp parallel do schedule(dynamic, 1) default(shared) private(m, j)
        do m = 1, nm
            do j = j0 + 1, e%nb - 1
                ym(:, j, m) = zm(:, j, m) - matmul(h%node(m, 1)%y(:, :, j), ym(:, j - 1, m))
            end do
        end do
        !$omp end parallel do
        do j = 0, j0 - 1
            x(:, :, j, :) = reshape(xd(:, :, j), [n2, e%nx, e%ns])
        end do
        do j = j0, e%nb - 1
            x(:, :, j, :) = reshape(matmul(ym(:, j, :), transpose(h%v)), [n2, e%nx, e%ns])
        end do
        status = NTV_OK
    end subroutine modal_hybrid_solve

    !> GMRES on the operator, right-preconditioned node by node.
    subroutine solve_any(e, lu, rhs, x, iterations, status)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(in) :: lu
        complex(dp), intent(in) :: rhs(:, :, 0:, :)
        complex(dp), allocatable, intent(out) :: x(:, :, :, :)
        integer, intent(out) :: iterations, status
        call coupled_solve(e, lu, rhs, x, iterations, status)
    end subroutine solve_any

    !> Collision-mode chain of the n = 0 preconditioner. After the row scaling
    !> the operator is S x R - I x C_l: streaming S per unit speed, R the
    !> diagonal ratio of lumped to Galerkin speed mass of each node, C_l the
    !> speed-species collision rows of degree l. In the eigenvectors V of
    !> R^-1 C_1 the l = 1 rows, which carry momentum conservation, are
    !> diagonal, so each mode m is a monoenergetic surface problem with the
    !> rates diag(V^-1 R^-1 C_l V)_m, factorised directly over all pair
    !> blocks. Requires a real spectrum of R^-1 C_1, which holds where C is
    !> self-adjoint in the Maxwellian-weighted product. vinv holds V^-1 R^-1,
    !> vu V^-1 and rv R V.
    subroutine factor_modal(e, h, status)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(inout) :: h
        integer, intent(out) :: status
        real(dp), allocatable :: c(:, :), w(:, :), wr(:), wi(:), work(:), vl(:, :), &
            chat(:, :)
        integer, allocatable :: piv(:), mode_status(:)
        integer :: nm, l, info, m, a, k, q
        nm = e%nx*e%ns
        status = NTV_SINGULAR
        allocate (c(nm, nm), w(nm, nm), wr(nm), wi(nm), work(8*nm), vl(1, 1), &
            h%v(nm, nm), h%vinv(nm, nm), h%vu(nm, nm), h%rv(nm, nm), piv(nm), &
            chat(nm, 0:2*e%nb - 1))
        c = collision_matrix(e, 1)
        w = c
        call dgeev('N', 'V', nm, w, nm, wr, wi, vl, 1, h%v, nm, work, size(work), info)
        if (info /= 0) return
        if (any(abs(wi) > 1e-10_dp*maxval(abs(wr)))) return
        w = h%v
        h%vu = 0
        do m = 1, nm
            h%vu(m, m) = 1
        end do
        call dgesv(nm, nm, w, nm, piv, h%vu, nm, info)
        if (info /= 0) return
        do l = 0, 2*e%nb - 1
            w = matmul(h%vu, matmul(collision_matrix(e, l), h%v))
            do m = 1, nm
                chat(m, l) = w(m, m)
            end do
        end do
        do a = 1, e%ns
            do k = 1, e%nx
                q = (a - 1)*e%nx + k
                h%vinv(:, q) = h%vu(:, q)/stream_ratio(e, k, a)
                h%rv(q, :) = stream_ratio(e, k, a)*h%v(q, :)
            end do
        end do
        h%modal = .true.
        allocate (h%sys(nm), h%node(nm, 1), mode_status(nm))
        !$omp parallel do schedule(dynamic, 1) default(shared) private(m)
        do m = 1, nm
            h%sys(m) = mode_system(e, chat(m, :))
            call factor_species(h%sys(m), 1, h%node(m, 1), mode_status(m))
        end do
        !$omp end parallel do
        status = merge(NTV_OK, NTV_SINGULAR, all(mode_status == NTV_OK))
    end subroutine factor_modal

    !> m = m - (I x R V) [U_m Y_m] (I x V^-1) at the boundary pair block j of
    !> the dense part, U_m Y_m the coupling of mode m to the chain above.
    subroutine subtract_modal_coupling(e, h, j, m)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(in) :: h
        integer, intent(in) :: j
        complex(dp), intent(inout) :: m(:, :)
        complex(dp), allocatable :: w(:, :, :)
        integer :: nm, n2, mm, q, qc
        nm = e%nx*e%ns
        n2 = 2*e%n
        allocate (w(n2, n2, nm))
        !$omp parallel do schedule(dynamic, 1) default(shared) private(mm)
        do mm = 1, nm
            w(:, :, mm) = matmul(energy_block(h%sys(mm), 1, j, j + 1), &
                h%node(mm, 1)%y(:, :, j + 1))
        end do
        !$omp end parallel do
        !$omp parallel do collapse(2) schedule(static) default(shared) private(q, qc, mm)
        do qc = 1, nm
            do q = 1, nm
                do mm = 1, nm
                    m((q - 1)*n2 + 1:q*n2, (qc - 1)*n2 + 1:qc*n2) = &
                        m((q - 1)*n2 + 1:q*n2, (qc - 1)*n2 + 1:qc*n2) &
                        - h%rv(q, mm)*h%vu(mm, qc)*w(:, :, mm)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine subtract_modal_coupling

    !> Scaled collision rows of degree l over (node, species), node fastest.
    function collision_matrix(e, l) result(c)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: l
        real(dp) :: c(e%nx*e%ns, e%nx*e%ns)
        integer :: a, b, k
        do b = 1, e%ns
            do a = 1, e%ns
                do k = 1, e%nx
                    c((a - 1)*e%nx + k, (b - 1)*e%nx + 1:b*e%nx) = &
                        collision_rows(e, k, a, b, l)/stream_ratio(e, k, a)
                end do
            end do
        end do
    end function collision_matrix

    !> Streaming coefficient of node k of species a after the row scaling.
    pure real(dp) function stream_ratio(e, k, a)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: k, a
        stream_ratio = e%row(k, a)*e%vt(a)*e%m3(k, k)
    end function stream_ratio

    !> Unit-speed surface system of one collision mode with rates rate(l).
    function mode_system(e, rate) result(f)
        type(energy_system_t), intent(in) :: e
        real(dp), intent(in) :: rate(0:)
        type(energy_system_t) :: f
        f%ns = 1; f%n = e%n; f%nb = e%nb; f%border = e%border
        f%cs => e%cs; f%co => e%co; f%cp => e%cp
        f%nx = 1
        f%nn = 2*e%n
        f%vt = [1.0_dp]
        f%tz = [1.0_dp]
        f%row = reshape([1.0_dp], [1, 1])
        f%m2 = reshape([1.0_dp], [1, 1])
        f%m3 = reshape([1.0_dp], [1, 1])
        f%m4 = reshape([1.0_dp], [1, 1])
        allocate (f%pre(1, 1, 1, 1, 0:ubound(rate, 1)))
        f%pre(1, 1, 1, 1, :) = rate
        f%coupled = .false.
    end function mode_system

    !> The perturbed-operator sources of f_10 (module header of dk_ntv) per
    !> node: the mirror force on the streaming mass, 1.5 source stream sdb
    !> (1 - xi^2)/2 d_xi f_10, and the collision source -1.5 source db C[f_10]
    !> with the operator's collision rows.
    subroutine add_f10_sources(e, t, h10, rhs)
        type(energy_system_t), intent(in) :: e
        type(ntv_system_t), intent(inout) :: t
        complex(dp), intent(in) :: h10(:, :, 0:, :)
        complex(dp), intent(inout) :: rhs(:, :, 0:, :)
        complex(dp), allocatable :: fold(:, :, :)
        complex(dp) :: wm(e%n), wc(e%n)
        integer :: a, b, kk, k, j, r, l
        wm = 1.5_dp*t%source*t%stream*t%sdb
        wc = 1.5_dp*t%source*t%db
        allocate (fold(2*e%n, e%nx, 0:e%nb - 1))
        do a = 1, e%ns
            do kk = 1, e%nx
                t%f10 = real(h10(:, kk, :, a), dp)
                call f10_foldings(t)
                fold(:, kk, :) = t%f10f
            end do
            do j = 0, e%nb - 1
                do k = 1, e%nx
                    do r = 0, 1
                        l = 2*j + r
                        rhs(r*e%n + 1:(r + 1)*e%n, k, j, a) = &
                            rhs(r*e%n + 1:(r + 1)*e%n, k, j, a) + e%row(k, a)*e%vt(a)*wm &
                            *matmul(fold(r*e%n + 1:(r + 1)*e%n, :, j), e%m3(k, :))
                        do b = 1, e%ns
                            rhs(r*e%n + 1:(r + 1)*e%n, k, j, a) = &
                                rhs(r*e%n + 1:(r + 1)*e%n, k, j, a) - wc &
                                *matmul(real(h10(r*e%n + 1:(r + 1)*e%n, :, j, b), dp), &
                                collision_rows(e, k, a, b, l))
                        end do
                    end do
                end do
            end do
        end do
    end subroutine add_f10_sources

    !> lorentz becomes the Galerkin Lorentz rows minus their lumped form
    !> -l (l + 1)/2 nu_D(x_k) m2_kk (scaled), so that subtracting it from a
    !> model's rows replaces its Galerkin pitch-angle part by the lumped one.
    subroutine lump_lorentz(e, species, lnlambda, grid, lorentz)
        type(energy_system_t), intent(in) :: e
        type(collop_species_t), intent(in) :: species(:)
        real(dp), intent(in) :: lnlambda(:, :)
        type(collop_grid_t), intent(in) :: grid
        real(dp), intent(inout) :: lorentz(:, :, :, :, 0:)
        real(dp) :: nu, x, y, g
        integer :: a, b, k, l
        do a = 1, e%ns
            do k = 1, e%nx
                x = grid%x(k + 1)
                nu = 0
                do b = 1, e%ns
                    y = x*e%vt(a)/e%vt(b)
                    g = (erf(y) - 2*y*exp(-y*y)/sqrt(pi))/(2*y*y)
                    nu = nu + collop_nuhat(species(a), species(b), lnlambda(a, b)) &
                        *(erf(y) - g)/x**3
                end do
                do l = 0, ubound(lorentz, 5)
                    lorentz(k, k, a, a, l) = lorentz(k, k, a, a, l) &
                        + e%row(k, a)*l*(l + 1)/2.0_dp*nu*e%m2(k, k)
                end do
            end do
        end do
    end subroutine lump_lorentz

    !> Fill the theta-block cache of e from its three parts.
    subroutine cache_blocks(e)
        type(energy_system_t), intent(inout) :: e
        integer :: j, d
        allocate (e%cs(2*e%n, 2*e%n, -1:1, 0:e%nb - 1), e%co(2*e%n, 2*e%n, -1:1, 0:e%nb - 1), &
            e%cp(2*e%n, 2*e%n, -1:1, 0:e%nb - 1))
        e%cs = 0; e%co = 0; e%cp = 0
        do j = 0, e%nb - 1
            do d = -1, 1
                if (j + d < 0 .or. j + d > e%nb - 1) cycle
                e%cs(:, :, d, j) = block(e%gs, j, j + d)
                e%co(:, :, d, j) = block(e%go, j, j + d)
                e%cp(:, :, d, j) = block(e%gp, j, j + d)
            end do
        end do
    end subroutine cache_blocks

    !> c_l(ak, bj) = a_l(ak, bj) (2l + 1)/(4 pi n_a pi^-3/2), times the row
    !> scale of node k of species a.
    subroutine scaled_rows(e, op, species, c)
        type(energy_system_t), intent(in) :: e
        type(collop_t), intent(in) :: op
        type(collop_species_t), intent(in) :: species(:)
        real(dp), allocatable, intent(out) :: c(:, :, :, :, :)
        integer :: a, b, l, kk, nx, ng
        nx = e%nx
        ng = op%grid%n
        allocate (c(nx, nx, e%ns, e%ns, 0:op%lmax))
        do l = 0, op%lmax
            do b = 1, e%ns
                do a = 1, e%ns
                    c(:, :, a, b, l) = op%a((a - 1)*ng + 2:a*ng, (b - 1)*ng + 2:b*ng, l) &
                        *(2*l + 1)*pi**1.5_dp/(4*pi*species(a)%n)
                    do kk = 1, nx
                        c(kk, :, a, b, l) = e%row(kk, a)*c(kk, :, a, b, l)
                    end do
                end do
            end do
        end do
    end subroutine scaled_rows

    !> Dense block (j, jc) of species a with its preconditioner collisions,
    !> unknowns ordered (theta row of the pair block, node).
    function energy_block(e, a, j, jc) result(m)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: a, j, jc
        complex(dp) :: m(e%nn, e%nn)
        complex(dp), dimension(2*e%n, 2*e%n) :: bs, bo, bp
        integer :: k, kk, r, l, i, n2
        n2 = 2*e%n
        bs = e%cs(:, :, jc - j, j)
        bo = e%co(:, :, jc - j, j)
        bp = e%cp(:, :, jc - j, j)
        do kk = 1, e%nx
            do k = 1, e%nx
                m((k - 1)*n2 + 1:k*n2, (kk - 1)*n2 + 1:kk*n2) = e%row(k, a) &
                    *(e%vt(a)*e%m3(k, kk)*bs + e%m2(k, kk)*bo + e%tz(a)*e%m4(k, kk)*bp)
            end do
        end do
        if (e%border .and. j == 0) then
            do k = 1, e%nx
                associate (rows => m((k - 1)*n2 + 1:(k - 1)*n2 + e%n, :))
                    rows = rows - spread(sum(rows, 1)/e%n, 1, e%n)
                end associate
            end do
        end if
        if (j /= jc) return
        do kk = 1, e%nx
            do k = 1, e%nx
                do r = 0, 1
                    l = 2*j + r
                    do i = 1, e%n
                        m((k - 1)*n2 + r*e%n + i, (kk - 1)*n2 + r*e%n + i) = &
                            m((k - 1)*n2 + r*e%n + i, (kk - 1)*n2 + r*e%n + i) &
                            - e%pre(k, kk, a, a, l)
                    end do
                end do
            end do
        end do
        if (e%border .and. j == 0) then
            do k = 1, e%nx
                associate (rows => m((k - 1)*n2 + 1:(k - 1)*n2 + e%n, &
                        (k - 1)*n2 + 1:(k - 1)*n2 + e%n))
                    rows = rows - spread(sum(rows, 1)/e%n, 1, e%n) + 1.0_dp/e%n
                end associate
            end do
        end if
    end function energy_block

    !> Block (j, jc) of species a without collisions applied to x.
    function offdiag_apply(e, a, j, jc, x) result(y)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: a, j, jc
        complex(dp), intent(in) :: x(:, :)
        complex(dp) :: y(2*e%n, e%nx)
        complex(dp), dimension(2*e%n, e%nx) :: ts, tp
        integer :: k
        ts = matmul(e%cs(:, :, jc - j, j), x)
        tp = matmul(e%cp(:, :, jc - j, j), x)
        y = matmul(ts, transpose(e%vt(a)*e%m3)) + matmul(tp, transpose(e%tz(a)*e%m4))
        if (j == jc) y = y + matmul(matmul(e%co(:, :, 0, j), x), transpose(e%m2))
        do k = 1, e%nx
            y(:, k) = e%row(k, a)*y(:, k)
            if (e%border .and. j == 0 .and. jc /= 0) &
                y(1:e%n, k) = y(1:e%n, k) - sum(y(1:e%n, k))/e%n
        end do
    end function offdiag_apply

    subroutine factor_species(e, a, f, status)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: a
        type(species_lu_t), intent(out) :: f
        integer, intent(out) :: status
        complex(dp), allocatable :: m(:, :)
        integer :: j, info
        allocate (f%s_lu(e%nn, e%nn, 0:e%nb - 1), f%y(e%nn, e%nn, e%nb - 1), &
            f%piv(e%nn, 0:e%nb - 1), m(e%nn, e%nn))
        status = NTV_SINGULAR
        do j = e%nb - 1, 0, -1
            m = energy_block(e, a, j, j)
            if (j < e%nb - 1) call zgemm('N', 'N', e%nn, e%nn, e%nn, (-1.0_dp, 0.0_dp), &
                energy_block(e, a, j, j + 1), e%nn, f%y(:, :, j + 1), e%nn, &
                (1.0_dp, 0.0_dp), m, e%nn)
            call zgetrf(e%nn, e%nn, m, e%nn, f%piv(:, j), info)
            if (info /= 0) return
            f%s_lu(:, :, j) = m
            if (j > 0) then
                f%y(:, :, j) = energy_block(e, a, j, j - 1)
                call zgetrs('N', e%nn, e%nn, m, e%nn, f%piv(:, j), f%y(:, :, j), e%nn, info)
            end if
        end do
        status = NTV_OK
    end subroutine factor_species

    !> GMRES on the species-coupled operator A = M + (C_pre - C_full),
    !> right-preconditioned with the species factorisations M.
    subroutine coupled_solve(e, lu, rhs, x, iterations, status)
        type(energy_system_t), intent(in) :: e
        type(hybrid_t), intent(in) :: lu
        complex(dp), intent(in) :: rhs(:, :, 0:, :)
        complex(dp), allocatable, intent(out) :: x(:, :, :, :)
        integer, intent(out) :: iterations, status
        complex(dp), allocatable :: b(:), u(:), w(:)
        real(dp) :: scale, residual
        integer :: info, ntot
        ntot = size(rhs)
        allocate (x, mold=rhs)
        b = reshape(rhs, [ntot])
        scale = norm2(abs(b))
        iterations = 0
        x = 0
        status = NTV_OK
        if (scale == 0) return
        b = b/scale
        allocate (u(ntot), w(ntot))
        u = 0
        call complex_gmres_operator(operator, b, u, energy_tolerance, 1000, 200, info, &
            iterations, residual)
        if (status /= NTV_OK) return
        status = NTV_NO_CONVERGENCE
        if (info /= KRYLOV_OK) return
        status = NTV_OK
        call precondition(u, w)
        if (status /= NTV_OK) return
        x = scale*reshape(w, shape(rhs))
    contains
        subroutine precondition(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            complex(dp), allocatable :: v(:, :, :, :), y(:, :, :, :)
            integer :: st
            v = reshape(input, shape(rhs))
            allocate (y, mold=v)
            call hybrid_solve(e, lu, v, y, st)
            if (st /= NTV_OK) status = st
            output = reshape(y, [ntot])
        end subroutine precondition
        subroutine operator(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            complex(dp), allocatable :: t(:)
            allocate (t(ntot))
            call precondition(input, t)
            output = reshape(apply_full(e, reshape(t, shape(rhs))), [ntot])
        end subroutine operator
    end subroutine coupled_solve

    !> Row k of the scaled collision block (a, b) of degree l of the operator.
    function collision_rows(e, k, a, b, l) result(c)
        type(energy_system_t), intent(in) :: e
        integer, intent(in) :: k, a, b, l
        real(dp) :: c(e%nx)
        if (e%coupled) then
            c = e%full(k, :, a, b, l)
        else
            c = merge(e%pre(k, :, a, b, l), 0.0_dp, a == b)
        end if
    end function collision_rows

    !> The species-coupled operator with the full collision rows.
    function apply_full(e, x) result(y)
        type(energy_system_t), intent(in) :: e
        complex(dp), intent(in) :: x(:, :, 0:, :)
        complex(dp) :: y(size(x, 1), size(x, 2), 0:size(x, 3) - 1, size(x, 4))
        integer :: a, b, j, r, l, k
        real(dp) :: corr(e%nx, e%nx)
        y = 0
        !$omp parallel do collapse(2) schedule(static) default(shared) &
        !$omp private(a, b, j, r, l, k, corr)
        do a = 1, e%ns
            do j = 0, e%nb - 1
                y(:, :, j, a) = offdiag_apply(e, a, j, j, x(:, :, j, a))
                do r = 0, 1
                    l = 2*j + r
                    y(r*e%n + 1:(r + 1)*e%n, :, j, a) = y(r*e%n + 1:(r + 1)*e%n, :, j, a) &
                        - matmul(x(r*e%n + 1:(r + 1)*e%n, :, j, a), &
                        transpose(e%pre(:, :, a, a, l)))
                end do
                if (j > 0) y(:, :, j, a) = y(:, :, j, a) &
                    + offdiag_apply(e, a, j, j - 1, x(:, :, j - 1, a))
                if (j < e%nb - 1) y(:, :, j, a) = y(:, :, j, a) &
                    + offdiag_apply(e, a, j, j + 1, x(:, :, j + 1, a))
                do b = 1, e%ns
                    do r = 0, 1
                        l = 2*j + r
                        if (e%coupled) then
                            corr = -e%full(:, :, a, b, l)
                            if (a == b) corr = corr + e%pre(:, :, a, a, l)
                        else
                            cycle
                        end if
                        y(r*e%n + 1:(r + 1)*e%n, :, j, a) = y(r*e%n + 1:(r + 1)*e%n, :, j, a) &
                            + matmul(x(r*e%n + 1:(r + 1)*e%n, :, j, b), transpose(corr))
                    end do
                end do
                if (e%border .and. j == 0) then
                    do k = 1, e%nx
                        y(1:e%n, k, 0, a) = y(1:e%n, k, 0, a) - sum(y(1:e%n, k, 0, a))/e%n &
                            + sum(x(1:e%n, k, 0, a))/e%n
                    end do
                end if
            end do
        end do
        !$omp end parallel do
    end function apply_full
