!> Axisymmetric (n = 0) drift-kinetic solve on the energy grid with the
!! linearised collision operator of dk_collop (#216, #184): the parallel
!! flows driven by the radial gradients for one species, the coefficients
!! D31 and D32 and the poloidal-rotation coefficient k = 5/2 - D32/D31.
!!
!! The perturbation is h(x, theta, xi) = sum_k phi_k(x) h_k(theta, xi) with
!! the GLL basis phi_k of dk_collop and h_k expanded as in dk_mono (Legendre
!! P_l in xi, Fourier collocation in theta). The energy direction is
!! Galerkin throughout, with every integral over [0, oo) including the
!! polynomial tail of the last element: per degree l and test function
!! phi_k,
!!   sum_j m_kj S h_j - c_kj h_j = d_k s_1,
!! m_kj = int phi_k phi_j x^3 exp(-x^2) dx (streaming carries the speed),
!! c the collision matrix of dk_collop per unit 4 pi/(2l+1) n pi^(-3/2) v_t,
!! S the streaming and mirror operator per unit length (dk_mono), s_1 the
!! radial drift per unit speed with I/Omega = 1, in units u v_t = 1, and
!! d_k = int phi_k x^4 W exp(-x^2) dx for the energy weight W = 1 (A_1) or
!! x^2 (A_2). Each node's rows are divided by m_kk. A lumped (row-sum) mass
!! is not usable: for degree >= 2 elements with the Maxwellian weight some
!! row sums are negative, which flips the sign of the collision rows there.
!!
!! S couples l to l -+ 1 and the mass and collision matrices couple the
!! nodes within a degree, so the system is block tridiagonal in l with
!! dense blocks of n_theta n_x unknowns. It is factorised once by band LU
!! with partial pivoting and reused for every right-hand side (drives and
!! tangents), followed by one step of iterative refinement. Two cheaper
!! routes were measured and rejected: a Krylov iteration preconditioned by
!! one monoenergetic factorisation per node did not converge in 400 steps
!! outside a band of about two decades in collisionality (energy scattering
!! dominates where that preconditioner omits it), and block elimination in
!! l without pivoting across degrees left backward errors up to 1e-2 at low
!! collisionality.
!!
!! The theta-constant degree-0 component of each node is bordered by a
!! source and the constraint that its flux-surface average vanishes, as in
!! dk_mono; it enters only the theta averages of the degree-0 rows and not
!! the flows.
!!
!! Flows: D3j = (4/sqrt(pi)) sum_k int phi_k x^3 exp(-x^2) dx <b g_k,1>/3
!! with <.> the flux-surface average of dk_mono. With a momentum-conserving
!! operator the D31 drive has the exact solution h = x xi/b, a rigid
!! parallel shift, so D31 = (4/sqrt(pi)) int x^4 exp(-x^2) dx/3 = 1/2 up to
!! the Fourier truncation of 1/b.
module kin6d_dk_energy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_dk_energy_leaves, only: dk_energy_rotation, dk_energy_rotation_jvp
    use kin6d_dk_mono, only: mono_factor_t, mono_factor, mono_coupling, mono_drives, &
        MONO_OK, trig_resample
    use kin6d_dk_collop, only: collop_t, collop_species_t, collop_grid_t, collop_build, &
        collop_vth, collop_nuhat, collop_moments, collop_make_grid, COLLOP_OK, COLLOP_LORENTZ, &
        COLLOP_TEST_PARTICLE, COLLOP_FULL, COLLOP_BAD_INPUT
    implicit none
    private
    public :: energy_axi_t, energy_axi_coefficients, energy_axi_batch
    public :: energy_estimate_t, energy_axi_estimate
    integer, parameter, public :: ENERGY_OK = 0, ENERGY_BAD_INPUT = 1, ENERGY_SINGULAR = 2
    real(dp), parameter :: pi = 4*atan(1.0_dp)

    type :: energy_axi_t
        real(dp) :: d31 = 0, d32 = 0, k = 0
        ! D31 of the rigid shift h = x xi/b in this discretisation, the
        ! exact value for a momentum-conserving operator.
        real(dp) :: d31_rigid = 0
        ! Componentwise backward error of the bordered solve, both drives.
        real(dp) :: residual(2) = 0
        ! Exact tangents with respect to s = ln(n/T^2): with lnLambda fixed the
        ! collision rows scale like nuhat/v_t ~ n T^(-2), and nothing else in
        ! the scaled system depends on n or T, so dk/dn = dk_ds/n and
        ! dk/dT = -2 dk_ds/T. k does not depend on the gradients.
        real(dp) :: dd31_ds = 0, dd32_ds = 0, dk_ds = 0
        character(len=11) :: certification = 'unavailable'
    end type

    !> Signed refined-minus-coarse differences, not continuous error bounds.
    type :: energy_estimate_t
        real(dp) :: d31 = 0, d32 = 0, k = 0
    end type

    !> Band LU of the bordered system and the operator data.
    type :: system_t
        integer :: nth = 0, lmax = 0, nx = 0, nb = 0, nnull = 0
        type(collop_t) :: op
        type(mono_factor_t) :: f
        ! Streaming mass m_kj, row scaling 1/m_kk, collision c_kj per degree.
        real(dp), allocatable :: m(:, :), row(:), c3(:, :, :), weight(:)
        real(dp), allocatable :: lo(:, :, :), up(:, :, :)
        integer :: kband = 0, nall = 0
        real(dp), allocatable :: band(:, :)
        integer, allocatable :: piv(:)
    end type

    interface
        subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, kl, ku, ldab
            real(dp), intent(inout) :: ab(ldab, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgbtrf
        subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb, ipiv(*)
            real(dp), intent(in) :: ab(ldab, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgbtrs
    end interface

contains

    !> Double pitch resolution and bisect each speed element; resample the
    !> same trigonometric field at 2*n_theta+1 points. Missing field harmonics
    !> and the polynomial speed tail are not qualified by this indicator.
    subroutine energy_axi_estimate(b, kappa, species, lnlambda, grid, n_leg, model, &
            res, estimate, status)
        real(dp), intent(in) :: b(:), kappa, lnlambda
        type(collop_species_t), intent(in) :: species
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: n_leg, model
        type(energy_axi_t), intent(out) :: res
        type(energy_estimate_t), intent(out) :: estimate
        integer, intent(out) :: status
        type(collop_grid_t) :: refined_grid
        type(energy_axi_t) :: refined
        integer :: st

        call energy_axi_coefficients(b, kappa, species, lnlambda, grid, n_leg, model, res, status)
        if (status /= ENERGY_OK) return
        call refinement_grid(grid, refined_grid, st)
        status = ENERGY_BAD_INPUT
        if (st /= COLLOP_OK) return
        call energy_axi_coefficients(trig_resample(b, 2*size(b) + 1), kappa, species, &
            lnlambda, refined_grid, 2*n_leg, model, refined, status)
        if (status /= ENERGY_OK) return
        estimate%d31 = refined%d31 - res%d31
        estimate%d32 = refined%d32 - res%d32
        estimate%k = refined%k - res%k
    end subroutine energy_axi_estimate

    subroutine refinement_grid(grid, refined, status)
        type(collop_grid_t), intent(in) :: grid
        type(collop_grid_t), intent(out) :: refined
        integer, intent(out) :: status
        real(dp), allocatable :: edges(:)
        integer :: e
        status = COLLOP_BAD_INPUT
        if (.not. same_grid(grid, grid)) return
        allocate (edges(0:2*grid%nelem))
        do e = 0, grid%nelem
            edges(2*e) = grid%edges(e)
        end do
        do e = 1, grid%nelem
            edges(2*e - 1) = (grid%edges(e - 1) + grid%edges(e))/2
        end do
        call collop_make_grid(edges(2*grid%nelem), 2*grid%nelem, grid%order, &
            refined, status, edges)
    end subroutine refinement_grid

    !> Shared energy grid and collision model; one status per surface so a
    !> rejected input does not discard the remaining profile. Optional
    !> estimates use the same refinement as energy_axi_estimate.
    recursive subroutine energy_axi_batch(b, kappa, species, lnlambda, grid, n_leg, model, &
            res, status, estimate)
        real(dp), intent(in) :: b(:, :), kappa(:), lnlambda(:)
        type(collop_species_t), intent(in) :: species(:)
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: n_leg, model
        type(energy_axi_t), intent(out) :: res(:)
        integer, intent(out) :: status(:)
        type(energy_estimate_t), intent(out), optional :: estimate(:)
        type(collop_t) :: op
        type(collop_grid_t) :: refined_grid
        type(energy_axi_t), allocatable :: refined(:)
        real(dp), allocatable :: fine_b(:, :)
        integer, allocatable :: fine_status(:)
        real(dp) :: lnl(1, 1)
        integer :: j, st, ns

        res = energy_axi_t()
        status = ENERGY_BAD_INPUT
        if (present(estimate)) estimate = energy_estimate_t()
        ns = size(species)
        if (any([size(b, 2), size(kappa), size(lnlambda), size(res), size(status)] /= ns)) return
        if (present(estimate)) then
            if (size(estimate) /= ns) return
        end if
        if (n_leg < 4 .or. .not. same_grid(grid, grid)) return
        do j = 1, ns
            if (.not. all(ieee_is_finite([species(j)%m, species(j)%z, species(j)%n, &
                species(j)%t, lnlambda(j)]))) cycle
            if (min(species(j)%m, species(j)%n, species(j)%t, lnlambda(j)) <= 0) cycle
            if (abs(species(j)%z) <= tiny(1.0_dp)) cycle
            lnl = lnlambda(j)
            call collop_build(species(j:j), lnl, grid, model, n_leg - 1, op, st)
            if (st == COLLOP_OK) exit
        end do
        if (j > ns) return
        do j = 1, ns
            call energy_axi_coefficients(b(:, j), kappa(j), species(j), lnlambda(j), &
                grid, n_leg, model, res(j), status(j), operator=op)
        end do
        if (.not. present(estimate)) return
        call refinement_grid(grid, refined_grid, st)
        if (st /= COLLOP_OK) then
            where (status == ENERGY_OK) status = ENERGY_BAD_INPUT
            return
        end if
        allocate (fine_b(2*size(b, 1) + 1, ns), refined(ns), fine_status(ns))
        fine_b = 0
        do j = 1, ns
            if (status(j) /= ENERGY_OK) cycle
            fine_b(:, j) = trig_resample(b(:, j), size(fine_b, 1))
        end do
        call energy_axi_batch(fine_b, kappa, species, lnlambda, refined_grid, &
            2*n_leg, model, refined, fine_status)
        do j = 1, ns
            if (status(j) /= ENERGY_OK) cycle
            status(j) = fine_status(j)
            if (status(j) /= ENERGY_OK) cycle
            estimate(j)%d31 = refined(j)%d31 - res(j)%d31
            estimate(j)%d32 = refined(j)%d32 - res(j)%d32
            estimate(j)%k = refined(j)%k - res(j)%k
        end do
    end subroutine energy_axi_batch

    !> D31, D32 and k of one species with self-collisions of the selected
    !> collision model on the surface b(theta) (n_theta odd, uniform theta),
    !> kappa = iota B0/(G + iota I) [1/m], n_leg Legendre polynomials and the
    !> energy grid. For one species the self-collision operator is n nuhat
    !> times a matrix that depends on nothing else: a surface scan passes the
    !> operator built once (any reference species, same grid and model, at
    !> least n_leg - 1 degrees) as `operator`, and it is rescaled.
    subroutine energy_axi_coefficients(b, kappa, species, lnlambda, grid, n_leg, model, &
            res, status, operator)
        real(dp), intent(in) :: b(:), kappa, lnlambda
        type(collop_species_t), intent(in) :: species
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: n_leg, model
        type(energy_axi_t), intent(out) :: res
        integer, intent(out) :: status
        type(collop_t), intent(in), optional :: operator
        type(system_t) :: sys
        real(dp), allocatable :: g(:, :, :, :), rhs(:, :, :, :), s(:, :, :), r(:, :, :), mu(:, :)
        real(dp), allocatable :: s3(:, :, :), corr(:, :, :, :), dg(:, :, :, :), dmu(:, :)
        real(dp) :: flow(2), dflow(2)
        real(dp), allocatable :: mom(:, :), d(:)
        integer :: j, k, it

        call setup(b, kappa, species, lnlambda, grid, n_leg, model, sys, status, operator)
        if (status /= ENERGY_OK) return
        allocate (mom(sys%nx, sys%nx), d(sys%nx))
        s = mono_drives(sys%f)
        allocate (rhs(sys%nth, 0:sys%lmax, sys%nx, 2), r(sys%nth, 0:sys%lmax, sys%nx), &
            s3(sys%nth, 0:sys%lmax, sys%nx), corr(sys%nth, 0:sys%lmax, sys%nx, 2))
        do j = 1, 2
            call collop_moments(grid, 0, 2 + 2*j, mom, d)
            do k = 1, sys%nx
                rhs(:, :, k, j) = sys%row(k)*d(k)*s(:, :, 1)
            end do
        end do
        call solve(sys, rhs, g, mu, status)
        if (status /= ENERGY_OK) return
        ! One step of iterative refinement with the operator.
        do it = 1, 1
            do j = 1, 2
                call apply(sys, g(:, :, :, j), mu(:, j), r)
                corr(:, :, :, j) = rhs(:, :, :, j) - r
            end do
            call solve(sys, corr, dg, dmu, status)
            if (status /= ENERGY_OK) return
            g = g + dg
            mu = mu + dmu
        end do
        ! Tangents: A dg = C g ds with C the scaled collision rows (A = S - C).
        do j = 1, 2
            call apply_collision(sys, g(:, :, :, j), corr(:, :, :, j))
        end do
        call solve(sys, corr, dg, dmu, status)
        if (status /= ENERGY_OK) return
        call collop_moments(grid, 0, 3, mom, d)
        d = 4/sqrt(pi)*d/3
        flow = 0
        dflow = 0
        do k = 1, sys%nx
            res%d31_rigid = res%d31_rigid + d(k)*grid%x(k)
            do j = 1, 2
                flow(j) = flow(j) + d(k)*sum(sys%weight*b*g(:, 1, k, j))
                dflow(j) = dflow(j) + d(k)*sum(sys%weight*b*dg(:, 1, k, j))
            end do
        end do
        status = ENERGY_SINGULAR
        if (.not. all(ieee_is_finite(flow)) .or. .not. all(ieee_is_finite(dflow))) return
        if (abs(flow(1)) <= tiny(1.0_dp)) return
        status = ENERGY_OK
        res%dd31_ds = dflow(1)
        res%dd32_ds = dflow(2)
        call dk_energy_rotation_jvp(flow(1), flow(2), dflow(1), dflow(2), res%dk_ds)
        ! Componentwise backward error |A h - rhs|/(|A| |h| + |rhs|).
        do j = 1, 2
            call apply(sys, g(:, :, :, j), mu(:, j), r)
            r = abs(r - rhs(:, :, :, j))
            call apply(sys, g(:, :, :, j), mu(:, j), s3, absolute=.true.)
            res%residual(j) = maxval(r/(s3 + abs(rhs(:, :, :, j)) + tiny(1.0_dp)))
        end do
        res%d31 = flow(1)
        res%d32 = flow(2)
        call dk_energy_rotation(flow(1), flow(2), res%k)
        if (.not. all(ieee_is_finite([res%k, res%dk_ds, res%residual]))) &
            status = ENERGY_SINGULAR
    end subroutine energy_axi_coefficients

    subroutine setup(b, kappa, species, lnlambda, grid, n_leg, model, sys, status, operator)
        real(dp), intent(in) :: b(:), kappa, lnlambda
        type(collop_species_t), intent(in) :: species
        type(collop_grid_t), intent(in) :: grid
        integer, intent(in) :: n_leg, model
        type(system_t), intent(out) :: sys
        integer, intent(out) :: status
        type(collop_t), intent(in), optional :: operator
        real(dp) :: lnl(1, 1), rate, reference_rate
        integer :: l, st

        status = ENERGY_BAD_INPUT
        if (n_leg < 4 .or. grid%n < 3) return
        if (.not. all(ieee_is_finite([species%m, species%z, species%n, species%t, &
            lnlambda, kappa]))) return
        if (min(species%m, species%n, species%t, lnlambda) <= 0) return
        if (abs(species%z) <= tiny(1.0_dp) .or. abs(kappa) <= tiny(1.0_dp)) return
        if (size(b) < 3 .or. mod(size(b), 2) == 0) return
        if (.not. all(ieee_is_finite(b)) .or. any(b <= 0)) return
        if (all(model /= [COLLOP_LORENTZ, COLLOP_TEST_PARTICLE, COLLOP_FULL])) return
        if (.not. same_grid(grid, grid)) return
        sys%nth = size(b)
        sys%lmax = n_leg - 1
        sys%nx = grid%n
        sys%nb = sys%nth*sys%nx
        lnl = lnlambda
        if (present(operator)) then
            if (operator%ns /= 1 .or. operator%lmax < sys%lmax .or. &
                operator%model /= model .or. operator%grid%n /= grid%n) return
            if (.not. same_grid(grid, operator%grid)) return
            if (.not. allocated(operator%sp) .or. .not. allocated(operator%lnlambda) .or. &
                .not. allocated(operator%a)) return
            if (size(operator%sp) /= 1 .or. any(shape(operator%lnlambda) /= [1, 1])) return
            if (size(operator%a, 1) /= grid%n .or. size(operator%a, 2) /= grid%n) return
            if (lbound(operator%a, 3) > 0 .or. ubound(operator%a, 3) < sys%lmax) return
            if (.not. all(ieee_is_finite(operator%a))) return
            if (.not. all(ieee_is_finite([operator%sp(1)%m, operator%sp(1)%z, &
                operator%sp(1)%n, operator%sp(1)%t, operator%lnlambda(1, 1)]))) return
            if (min(operator%sp(1)%m, operator%sp(1)%n, operator%sp(1)%t, &
                operator%lnlambda(1, 1)) <= 0) return
            if (abs(operator%sp(1)%z) <= tiny(1.0_dp)) return
            reference_rate = operator%sp(1)%n*collop_nuhat(operator%sp(1), &
                operator%sp(1), operator%lnlambda(1, 1))
            rate = species%n*collop_nuhat(species, species, lnlambda)
            if (.not. all(ieee_is_finite([rate, reference_rate]))) return
            if (min(rate, reference_rate) <= 0) return
            sys%op = operator
            sys%op%a = operator%a(:, :, 0:sys%lmax)*(rate/reference_rate)
            sys%op%sp = [species]
            sys%op%lnlambda = lnl
            sys%op%lmax = sys%lmax
        else
            call collop_build([species], lnl, grid, model, sys%lmax, sys%op, st)
            if (st /= COLLOP_OK) return
        end if
        allocate (sys%m(sys%nx, sys%nx), sys%row(sys%nx), sys%c3(sys%nx, sys%nx, 0:sys%lmax))
        call collop_moments(grid, 3, 0, sys%m, sys%row)
        sys%row = [(1/sys%m(l, l), l=1, sys%nx)]
        do l = 0, sys%lmax
            sys%c3(:, :, l) = sys%op%a(:, :, l)*(2*l + 1)/(4*pi*species%n/pi**1.5_dp &
                *collop_vth(species))
        end do
        ! The streaming blocks come from dk_mono; its deflection frequency
        ! is not used.
        status = ENERGY_SINGULAR
        call mono_factor(b, kappa, 1.0_dp, n_leg, sys%f, st)
        if (st /= MONO_OK) return
        sys%weight = sys%f%p%weight
        allocate (sys%lo(sys%nth, sys%nth, 0:sys%lmax), sys%up(sys%nth, sys%nth, 0:sys%lmax))
        do l = 0, sys%lmax
            call mono_coupling(sys%f, l, sys%lo(:, :, l), sys%up(:, :, l))
        end do
        call null_space(sys)
        call factor(sys, st)
        if (st /= 0) return
        status = ENERGY_OK
    end subroutine setup

    logical function same_grid(a, b) result(ok)
        type(collop_grid_t), intent(in) :: a, b
        ok = .false.
        if (a%n < 3 .or. a%order < 3 .or. a%nelem < 1) return
        if (a%n /= a%nelem*(a%order - 1) + 1) return
        if (a%n /= b%n .or. a%order /= b%order .or. a%nelem /= b%nelem) return
        if (.not. allocated(a%x) .or. .not. allocated(a%edges) .or. &
            .not. allocated(a%xi) .or. .not. allocated(a%d)) return
        if (.not. allocated(b%x) .or. .not. allocated(b%edges) .or. &
            .not. allocated(b%xi) .or. .not. allocated(b%d)) return
        if (size(a%x) /= a%n .or. size(b%x) /= b%n) return
        if (size(a%edges) /= a%nelem + 1 .or. size(b%edges) /= b%nelem + 1) return
        if (lbound(a%edges, 1) /= 0 .or. lbound(b%edges, 1) /= 0) return
        if (size(a%xi) /= a%order .or. size(b%xi) /= b%order) return
        if (any(shape(a%d) /= [a%order, a%order]) .or. &
            any(shape(b%d) /= [b%order, b%order])) return
        if (.not. all(ieee_is_finite(a%x)) .or. .not. all(ieee_is_finite(a%edges)) .or. &
            .not. all(ieee_is_finite(a%xi)) .or. .not. all(ieee_is_finite(a%d))) return
        if (abs(a%edges(0)) > 0 .or. any(a%edges(1:) <= a%edges(:a%nelem - 1))) return
        ok = all(abs(a%x - b%x) <= 0) .and. all(abs(a%edges - b%edges) <= 0) &
            .and. all(abs(a%xi - b%xi) <= 0) .and. all(abs(a%d - b%d) <= 0)
    end function same_grid

    !> Bordering of the theta-constant degree-0 component of each node by a
    !> source on that node and the constraint <h_k,0> = 0. That component
    !> appears only in the theta averages of the degree-0 rows (streaming
    !> annihilates it), so the bordering leaves every other component exact.
    !> Determining it from the weak degree-0 collision coupling instead
    !> amplifies the discrete flux defect of the Fourier streaming by
    !> 1/nu and loses the solution at low collisionality.
    subroutine null_space(sys)
        type(system_t), intent(inout) :: sys
        sys%nnull = sys%nx
    end subroutine null_space

    !> Diagonal block of degree l: scaled collision rows, theta-diagonal.
    function diagonal(sys, l) result(m)
        type(system_t), intent(in) :: sys
        integer, intent(in) :: l
        real(dp) :: m(sys%nb, sys%nb)
        integer :: k, j, i
        m = 0
        do j = 1, sys%nx
            do k = 1, sys%nx
                do i = 1, sys%nth
                    m((k - 1)*sys%nth + i, (j - 1)*sys%nth + i) = -sys%row(k)*sys%c3(k, j, l)
                end do
            end do
        end do
    end function diagonal

    !> Materialize the streaming block directly instead of multiplying by
    !> an identity through two dense matrix products.
    function couple(sys, blk) result(y)
        type(system_t), intent(in) :: sys
        real(dp), intent(in) :: blk(:, :)
        real(dp) :: y(sys%nb, sys%nb)
        integer :: k, j
        do j = 1, sys%nx
            do k = 1, sys%nx
                y((k - 1)*sys%nth + 1:k*sys%nth, (j - 1)*sys%nth + 1:j*sys%nth) = &
                    (sys%row(k)*sys%m(k, j))*blk
            end do
        end do
    end function couple

    !> Band LU with partial pivoting (LAPACK dgbtrf) of the bordered system.
    !> Unknowns: block 0 = (h_0 on all nodes, the sources mu), then h_l,
    !> l >= 1. Inspect the assembled nonzeros before allocating the band so
    !> sparse speed-mass blocks do not reserve unused diagonals. Block
    !> elimination in l without pivoting across degrees was tried first and
    !> left backward errors up to 1e-2 at low collisionality, which
    !> refinement did not remove.
    subroutine factor(sys, info)
        type(system_t), intent(inout) :: sys
        integer, intent(out) :: info
        integer :: l, nb, n0, k, pass
        real(dp), allocatable :: blk(:, :)
        nb = sys%nb
        n0 = nb + sys%nx
        sys%nall = n0 + sys%lmax*nb
        sys%kband = 0
        do pass = 1, 2
            if (pass == 2) then
                allocate (sys%band(3*sys%kband + 1, sys%nall), sys%piv(sys%nall))
                sys%band = 0
            end if
            do l = 0, sys%lmax
                call put(start(l), start(l), diagonal(sys, l))
                if (l > 0) call put(start(l), start(l - 1), couple(sys, sys%lo(:, :, l)))
                if (l < sys%lmax) call put(start(l), start(l + 1), &
                    couple(sys, sys%up(:, :, l)))
            end do
            allocate (blk(nb, sys%nx))
            blk = 0
            do k = 1, sys%nx
                blk((k - 1)*sys%nth + 1:k*sys%nth, k) = 1
            end do
            call put(1, nb + 1, blk)
            blk = 0
            do k = 1, sys%nx
                blk((k - 1)*sys%nth + 1:k*sys%nth, k) = sys%weight
            end do
            call put(nb + 1, 1, transpose(blk))
            deallocate (blk)
        end do
        call dgbtrf(sys%nall, sys%nall, sys%kband, sys%kband, sys%band, 3*sys%kband + 1, &
            sys%piv, info)
    contains
        integer function start(l)
            integer, intent(in) :: l
            start = 1
            if (l > 0) start = n0 + (l - 1)*nb + 1
        end function start
        subroutine put(i0, j0, m)
            integer, intent(in) :: i0, j0
            real(dp), intent(in) :: m(:, :)
            integer :: i, j
            do j = 1, size(m, 2)
                do i = 1, size(m, 1)
                    if (m(i, j) == 0) cycle
                    if (pass == 1) then
                        sys%kband = max(sys%kband, abs(i0 + i - j0 - j))
                    else
                        sys%band(2*sys%kband + 1 + i0 + i - j0 - j, j0 + j - 1) = m(i, j)
                    end if
                end do
            end do
        end subroutine put
    end subroutine factor

    !> Solve for the right-hand sides rhs(theta, l, node, :); mu are the
    !> source amplitudes of the border.
    subroutine solve(sys, rhs, x, mu, status)
        type(system_t), intent(in) :: sys
        real(dp), intent(in) :: rhs(:, 0:, :, :)
        real(dp), allocatable, intent(out) :: x(:, :, :, :), mu(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: z(:, :)
        integer :: l, nb, nr, info, n0, i0
        nb = sys%nb
        n0 = nb + sys%nx
        nr = size(rhs, 4)
        allocate (z(sys%nall, nr), x(sys%nth, 0:sys%lmax, sys%nx, nr))
        z = 0
        z(1:nb, :) = reshape(rhs(:, 0, :, :), [nb, nr])
        do l = 1, sys%lmax
            i0 = n0 + (l - 1)*nb
            z(i0 + 1:i0 + nb, :) = reshape(rhs(:, l, :, :), [nb, nr])
        end do
        status = ENERGY_SINGULAR
        call dgbtrs('N', sys%nall, sys%kband, sys%kband, nr, sys%band, 3*sys%kband + 1, &
            sys%piv, z, sys%nall, info)
        if (info /= 0) return
        x(:, 0, :, :) = reshape(z(1:nb, :), [sys%nth, sys%nx, nr])
        mu = z(nb + 1:n0, :)
        do l = 1, sys%lmax
            i0 = n0 + (l - 1)*nb
            x(:, l, :, :) = reshape(z(i0 + 1:i0 + nb, :), [sys%nth, sys%nx, nr])
        end do
        status = ENERGY_OK
    end subroutine solve

    !> Bordered operator: the scaled rows plus the sources mu per node. With
    !> absolute, every matrix and h enter by modulus: |A| |h| for the
    !> backward error of a solution.
    subroutine apply(sys, h, mu, y, absolute)
        type(system_t), intent(in) :: sys
        real(dp), intent(in) :: h(:, 0:, :), mu(:)
        real(dp), intent(out) :: y(:, 0:, :)
        logical, intent(in), optional :: absolute
        real(dp) :: t(sys%nth, sys%nx), scaled(sys%nx, sys%nx), coll(sys%nx, sys%nx)
        real(dp) :: lo(sys%nth, sys%nth), up(sys%nth, sys%nth), hh(sys%nth, 0:sys%lmax, sys%nx)
        real(dp) :: sgn
        logical :: ab
        integer :: k, l
        ab = .false.
        if (present(absolute)) ab = absolute
        sgn = -1
        if (ab) sgn = 1
        scaled = sys%m*spread(sys%row, 2, sys%nx)
        hh = h
        if (ab) then
            scaled = abs(scaled)
            hh = abs(h)
        end if
        do l = 0, sys%lmax
            lo = sys%lo(:, :, l)
            up = sys%up(:, :, l)
            coll = sys%c3(:, :, l)*spread(sys%row, 2, sys%nx)
            if (ab) then
                lo = abs(lo)
                up = abs(up)
                coll = abs(coll)
            end if
            t = 0
            if (l > 0) t = matmul(lo, hh(:, l - 1, :))
            if (l < sys%lmax) t = t + matmul(up, hh(:, l + 1, :))
            y(:, l, :) = matmul(t, transpose(scaled)) + sgn*matmul(hh(:, l, :), transpose(coll))
        end do
        do k = 1, sys%nx
            if (ab) then
                y(:, 0, k) = y(:, 0, k) + abs(mu(k))
            else
                y(:, 0, k) = y(:, 0, k) + mu(k)
            end if
        end do
    end subroutine apply

    !> Scaled collision rows C h (the operator is A = S - C + border).
    subroutine apply_collision(sys, h, y)
        type(system_t), intent(in) :: sys
        real(dp), intent(in) :: h(:, 0:, :)
        real(dp), intent(out) :: y(:, 0:, :)
        integer :: l
        do l = 0, sys%lmax
            y(:, l, :) = matmul(h(:, l, :), transpose(sys%c3(:, :, l)*spread(sys%row, 2, sys%nx)))
        end do
    end subroutine apply_collision

end module kin6d_dk_energy
