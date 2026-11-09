!> Monoenergetic NTV flux kernel of one toroidal harmonic n with Lorentz
!! pitch-angle scattering (#183), the simpler reference discretisation of
!! drift-kinetics.md: Legendre polynomials P_l(xi) and Fourier collocation in
!! theta, as for the axisymmetric problem (dk_mono).
!!
!! With f = F(theta, xi) exp(-i n zeta), the local equation per unit speed is
!!   xi k (d_theta - i n q) F - (1 - xi^2)/2 k d_theta(ln b) d_xi F
!!     + i (n/v) Omega F - (nu/2) d_xi((1 - xi^2) d_xi F) = lambda S/v,
!! b = |B|/B0, lambda = (1 - xi^2)/b, nu = nu_D/v,
!! Omega = omega_E + pre lambda (cos(theta) - offset) (grad-B precession,
!! bounce average -(mu B0/(e chi')) eps' (2E/K - 1 - offset)) and
!! S = i n sum_m c_m exp(i m theta) for B1 = Re sum_m c_m exp(i(m theta - n zeta)).
!! Its bounce average over a trapped orbit is the equation of the stage-1
!! bounce reference, i n Omega I_th f - 2 nu_D d_lambda(lambda I_xi d_lambda f)
!! = lambda S_hat, and the flux kernel is that reference's
!!   R = int dlambda oint dtheta/|xi| lambda Re(S^* f)
!!     = int dtheta int dxi (2/b) lambda Re(S^* F),
!! here over trapped and passing particles. (1 - xi^2) couples l to l +- 2,
!! so the rows l = 2j, 2j + 1 form pair blocks of size 2 n_theta, block
!! tridiagonal in j, factorised once by block elimination from the highest j.
!! Tangents with respect to nu and omega_E, the perturbation amplitudes c_m
!! and a direction in (eps, k, offset, c_J) reuse the factorisation. With
!! A F = H and R = <H, F>_W, one coefficient-space adjoint solves
!! A^H lambda = W H; state derivatives contract lambda with dH - dA F.
!! Direct goal derivatives retain the physical weights and FortSym JVPs.
!!
!! Three model switches, each with a measurable defect, select the equation of
!! NEO-2-QL and NEO-RT instead of the stage-1 bounce reference:
!!   weight  NTV_WEIGHT_MU (default): drift and drive carry v_perp^2/(2 v^2 b)
!!           = lambda, the grad-B part of the reference; NTV_WEIGHT_FULL: they
!!           carry (v_par^2 + v_perp^2/2)/(v^2 b) = (1 + xi^2)/b, grad-B plus
!!           curvature at low beta, equal to lambda at xi = 0;
!!   current c_J  adds the parallel-current part of the curvature drift,
!!           pre c_J xi^2, c_J = 2 B0 (I' + q G')/((G + iota I) q eps') with the
!!           covariant Boozer components G, I and ' = d/ds (zero by default);
!!   boozer_stream  b.grad theta = k b (Boozer angles, Jacobian ~ 1/B^2)
!!           instead of the constant k, with the surface measure
!!           dtheta/(b.grad theta B) ~ 1/b^2 in the flux.
!!
!! Two further source terms complete the right-hand side of Kasilov et al.,
!! Phys. Plasmas 21, 092506 (2014), Eq. (61), which NEO-2-QL solves (#192;
!! the circular covariant factors reduce to 1/q at leading order in eps).
!! The rest of the perturbed radial drift beyond S,
!! (1 + wsign xi^2)/(b v q) sum_m (i (m - n q) - 2 d_theta ln b) c_m e^{i m theta}:
!! its parallel-derivative part and the change -2 (dB/B) d_theta ln B of the
!! drift factor. And the perturbed operator acting on the axisymmetric
!! solution f_10,
!! stream/b sum_m (i (m - n q) - d_theta ln b) c_m e^{i m theta}
!! (1 - xi^2)/2 d_xi f_10 - sum_m c_m e^{i m theta}/b C[f_10]: the mirror
!! force of the parallel gradient b (d_theta - i n q)(dB/B), along the line
!! in this frame as in Eq. (61) at fixed field-line label, and the collision
!! term, with C the Lorentz operator of the
!! equation above. They are switched individually by source_dpar and
!! source_dlf10; the second needs f_10 on this kernel's grid through the
!! optional argument of the same name, one theta point per row (rows 1:n
!! even degrees 2j, rows n+1:2n odd 2j+1, one column per pair block) per
!! unit speed in the kernel's normalisation. With source_dlf10 on, the
!! directional surface tangent is unavailable (its geometry tangents of the
!! f_10 sources are not implemented) and the two-resolution refinement
!! would need f_10 on the refined grid, so indicators must be passed as
!! false; source_dpar needs q /= 0 (q = 1/iota keeps its sign). The
!! amplitude sensitivities keep the adjoint form with all sources of the
!! response contracted into it, and
!! the torque functional still pairs only with S.
!! A general field line enters through an ntv_line_t provider: its sample
!! returns, at the theta points of the solver, the coefficients the circular
!! leaf returns otherwise (b.grad theta, the mirror factor
!! b.grad theta d_theta(ln b)/2, the drift multiplying (1 + wsign xi^2), the
!! pure xi^2 drift, the source factor 2/(3 b) and the flux measure), with b
!! = |B|/B_ref. The perturbation spectrum (any m, n) and q stay in the
!! surface. Tangents along (eps, k, offset, c_J) exist only for the circular
!! surface.
module milca_dk_ntv
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use fortnum_fourier, only: fourier_derivative
    use milca_dk_ntv_leaves, only: dk_ntv_point, dk_ntv_point_jvp
    use fortnum_krylov, only: complex_gmres_operator, KRYLOV_OK
    use milca_dk_collop, only: collop_species_t, collop_grid_t, collop_t, collop_build, &
        collop_moments, collop_vth, collop_nuhat, COLLOP_OK, COLLOP_FULL, &
        COLLOP_TEST_PARTICLE, COLLOP_LORENTZ
    use milca_dk_mono, only: mono_factor_t, mono_factor, mono_solve, mono_drives, MONO_OK
    implicit none
    private
    public :: ntv_surface_t, ntv_kernel_t, ntv_kernel, ntv_line_t, ntv_circular_line_t
    public :: ntv_f10, ntv_factor_t, ntv_energy_coefficients
    integer, parameter, public :: NTV_OK = 0, NTV_BAD_INPUT = 1, NTV_SINGULAR = 2
    integer, parameter, public :: NTV_WEIGHT_MU = 0, NTV_WEIGHT_FULL = 1

    !> Field-line coefficients of a general surface (see the module header).
    type, abstract :: ntv_line_t
    contains
        procedure(line_sample), deferred :: sample
    end type

    abstract interface
        subroutine line_sample(self, theta, stream, mirror, drift, drift2, source, &
                measure)
            import :: ntv_line_t, dp
            class(ntv_line_t), intent(in) :: self
            real(dp), intent(in) :: theta(:)
            real(dp), intent(out) :: stream(:), mirror(:), drift(:), drift2(:), &
                source(:), measure(:)
        end subroutine line_sample
    end interface

    !> The circular stage-1 field line as a provider: the same leaf and
    !> arguments as a surface without line.
    type, extends(ntv_line_t) :: ntv_circular_line_t
        real(dp) :: eps = 0, k = 0, offset = 0, current = 0
        logical :: boozer_stream = .false.
    contains
        procedure :: sample => circular_sample
    end type

    !> Surface of the frozen-field benchmark and the perturbation spectrum;
    !> with line allocated, eps, k, offset, current and boozer_stream are
    !> replaced by the provider.
    type :: ntv_surface_t
        real(dp) :: eps = 0, k = 0, q = 0, offset = 0
        integer :: n = 0
        integer :: weight = NTV_WEIGHT_MU
        real(dp) :: current = 0
        logical :: boozer_stream = .false.
        ! The two #192 quasilinear source switches of the module header.
        logical :: source_dpar = .false.
        logical :: source_dlf10 = .false.
        class(ntv_line_t), allocatable :: line
        integer, allocatable :: m(:)
        complex(dp), allocatable :: c(:)
    end type

    !> Flux kernel, its tangents and its two-resolution indicator. dr_dc(j)
    !> = dR/dRe(c_j) + i dR/dIm(c_j); dr_dsurface is the derivative along
    !> the (eps, k, offset, current) components of the optional direction.
    type :: ntv_kernel_t
        real(dp) :: r = 0, dr_dnu = 0, dr_domega = 0, dr_dsurface = 0
        complex(dp), allocatable :: dr_dc(:)
        ! Legendre spectrum of the solution: max over theta of |F_l|, l = 0..L-1.
        real(dp), allocatable :: spectrum(:)
        ! Directional two-resolution indicators: theta points 2 n + 1, and
        ! Legendre polynomials 2 L; their sum indicates the error of r.
        real(dp) :: indicator_theta = 0, indicator_leg = 0
        ! GMRES iterations of the forward and adjoint solves preconditioned
        ! by a kept factorisation; zero for a direct solve.
        integer :: iterations = 0
    end type

    type :: ntv_system_t
        integer :: n = 0, nblock = 0
        ! Legendre factors of the drift and drive weight in xi^2: (1 - xi^2)
        ! has sign -1, (1 + xi^2) sign +1.
        real(dp) :: wsign = -1
        real(dp) :: nu = 0, nq = 0
        complex(dp) :: shift = 0, omega = 0, precess = 0
        real(dp), allocatable :: d(:, :), stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:)
        complex(dp), allocatable :: s(:), s_lu(:, :, :), y(:, :, :)
        ! The parallel-derivative drive sd, the raw perturbation profiles
        ! sdb = sum i m c_m e^{i m theta} and db = sum c_m e^{i m theta}, and
        ! the axisymmetric solution f_10 with its two operator foldings:
        ! f10f = (1 - xi^2)/2 d_xi f_10 per degree, f10c = f_10 scaled by
        ! nu l(l + 1)/2, the Lorentz rows of the operator.
        logical :: dpar = .false., dlf10 = .false.
        complex(dp), allocatable :: sd(:), sdb(:), db(:)
        real(dp), allocatable :: f10(:, :), f10f(:, :), f10c(:, :)
        integer, allocatable :: piv(:, :)
    end type

    !> The factorisation of one kernel, kept to precondition kernels of
    !> nearby surfaces, energies or parameters on the same grid.
    type :: ntv_factor_t
        private
        type(ntv_system_t), allocatable :: sys
    end type

    real(dp), parameter :: pi = 4*atan(1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    ! Relative residual of the preconditioned solves, well below the
    ! discretisation error of the torque.
    real(dp), parameter :: gmres_tolerance = 1e-11_dp

    interface
        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            complex(dp), intent(inout) :: c(ldc, *)
        end subroutine zgemm
        subroutine zgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in) :: m, n, lda
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgetrf
        subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
            complex(dp), intent(in) :: a(lda, *)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgetrs
    end interface


    ! Energy-coupled kernel (dk_ntv_energy.inc).
    real(dp), parameter :: charge = 1.602176634e-19_dp
    ! Relative residual of the energy-coupled GMRES solves.
    real(dp), parameter :: energy_tolerance = 1e-7_dp

    !> Factorised test-particle operator of one species.
    type :: species_lu_t
        complex(dp), allocatable :: s_lu(:, :, :), y(:, :, :)
        integer, allocatable :: piv(:, :)
    end type

    !> Hybrid preconditioner: the node factors for pair blocks j >= j0 and
    !> one dense factorisation of all nodes and species with the operator's
    !> own collision rows for j < j0, where the field-particle part and the
    !> energy scattering are strongest.
    type :: hybrid_t
        integer :: j0 = 0, nf = 0
        type(species_lu_t), allocatable :: node(:, :)
        ! The single-node systems, sharing the theta-block caches.
        type(energy_system_t), allocatable :: sys(:)
        complex(dp), allocatable :: s_lu(:, :, :), y(:, :, :)
        integer, allocatable :: piv(:, :)
    end type

    ! Pair blocks factorised densely: degrees l < 2 low_pairs.
    integer, parameter :: low_pairs = 2

    !> Everything the operator applies need.
    type :: energy_system_t
        integer :: ns = 0, nx = 0, nn = 0, nb = 0, n = 0
        type(ntv_system_t) :: gs, go, gp
        real(dp), allocatable :: m2(:, :), m3(:, :), m4(:, :), d4(:), d6(:)
        real(dp), allocatable :: vt(:), tz(:), row(:, :)
        ! Scaled collision rows c_l per species pair: (nx, nx, ns, ns, 0:lmax)
        ! for the preconditioner (pre) and the full operator (full).
        real(dp), allocatable :: pre(:, :, :, :, :), full(:, :, :, :, :)
        logical :: coupled = .false.
        ! n = 0 system: each node's theta-averaged degree-0 row is replaced by
        ! the constraint that the theta average of its degree-0 part vanish.
        logical :: border = .false.
        ! theta blocks (j, j + d), d = -1, 0, 1, of the streaming, rotation and
        ! precession parts, shared by every species and node.
        complex(dp), pointer, contiguous :: cs(:, :, :, :) => null(), &
            co(:, :, :, :) => null(), cp(:, :, :, :) => null()
    end type

contains

    !> The axisymmetric solution f_10 of source_dlf10 on this kernel's grid
    !> and layout: the n = 0 equation driven by the axisymmetric radial drift
    !> (1 + xi^2)/(b v q) d_theta b, of which the perturbed drift
    !> (1 + xi^2)/(b v) (i m/q) c_m e^{i m theta} of the two drives above is
    !> the first-order part. With b.grad theta = kappa b this is the dk_mono
    !> radial-drive solution, drive -kappa (d_theta b/b)(1 + xi^2)/2, times
    !> -2/(v q kappa), so the surface needs Boozer streaming and the weight
    !> (1 + xi^2)/b; the density null space is fixed by <f_10> = 0, which
    !> both f_10 sources annihilate.
    subroutine ntv_f10(sf, v, nu_d, n_theta, n_leg, f10, status)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d
        integer, intent(in) :: n_theta, n_leg
        real(dp), allocatable, intent(out) :: f10(:, :)
        integer, intent(out) :: status
        type(mono_factor_t) :: f
        real(dp), allocatable :: theta(:), b(:), stream(:), mirror(:), drift(:), &
            drift2(:), source(:), measure(:), stream_b(:), mirror_b(:), &
            measure_b(:), s(:, :, :), g(:, :, :)
        real(dp) :: kappa, mu(1)
        integer :: n, i, j

        status = NTV_BAD_INPUT
        if (mod(n_theta, 2) == 0 .or. n_theta < 3 .or. mod(n_leg, 2) /= 0 .or. &
            n_leg < 4 .or. v <= 0 .or. nu_d <= 0 .or. sf%q == 0 .or. &
            sf%weight /= NTV_WEIGHT_FULL) return
        if (.not. allocated(sf%line) .and. .not. sf%boozer_stream) return
        n = n_theta
        allocate (theta(n), b(n), stream(n), mirror(n), drift(n), drift2(n), &
            source(n), measure(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        if (allocated(sf%line)) then
            call sf%line%sample(theta, stream, mirror, drift, drift2, source, measure)
            b = 2/(3*source)
        else
            allocate (stream_b(n), mirror_b(n), measure_b(n))
            call dk_ntv_point(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, b, &
                mirror, drift, source, measure, stream_b, mirror_b, measure_b)
            stream = stream_b
        end if
        kappa = stream(1)/b(1)
        if (any(abs(stream/b - kappa) > 1e-12_dp*abs(kappa))) return
        status = NTV_SINGULAR
        call mono_factor(b, kappa, nu_d/v, n_leg, f, i)
        if (i /= MONO_OK) return
        s = mono_drives(f)
        call mono_solve(f, s(:, :, 1:1), [0.0_dp], g, mu, i)
        if (i /= MONO_OK) return
        allocate (f10(2*n, 0:n_leg/2 - 1))
        do j = 0, n_leg/2 - 1
            f10(1:n, j) = g(:, 2*j, 1)
            f10(n + 1:, j) = g(:, 2*j + 1, 1)
        end do
        f10 = -2/(v*sf%q*kappa)*f10
        status = NTV_OK
    end subroutine ntv_f10

    !> Kernel R at speed v [m/s], deflection frequency nu_d [1/s], omega_E
    !> [rad/s] and precession prefactor pre [rad/s], on n_theta (odd) points
    !> and n_leg (even) Legendre polynomials, with the directional
    !> indicators at 2 n_theta + 1 points and at 2 n_leg polynomials.
    subroutine ntv_kernel(sf, v, nu_d, omega_e, pre, n_theta, n_leg, kern, status, &
            indicators, direction, f10, keep, precondition)
        type(ntv_surface_t), intent(in) :: sf
        ! Tangent direction: its eps, k, offset and current components.
        type(ntv_surface_t), intent(in), optional :: direction
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        type(ntv_kernel_t), intent(out) :: kern
        integer, intent(out) :: status
        ! A resolution ladder supplies its own two-resolution estimate and
        ! can skip the refined solves, which dominate memory.
        logical, intent(in), optional :: indicators
        ! The axisymmetric solution for the perturbed-operator source
        ! (source_dlf10): one theta point per row, even degrees above odd.
        real(dp), intent(in), optional, dimension(:, :) :: f10
        ! keep receives this kernel's factorisation; with precondition the
        ! kernel is not factorised but solved by GMRES, right-preconditioned
        ! by that factorisation of a nearby kernel on the same grid.
        type(ntv_factor_t), intent(out), optional :: keep
        type(ntv_factor_t), intent(inout), optional :: precondition
        type(ntv_system_t), allocatable :: a
        type(ntv_system_t) :: da
        logical :: krylov
        complex(dp), allocatable :: x(:, :), adj(:, :), rhs(:, :), ds(:), gradient(:), &
            xf(:), xw(:)
        real(dp), allocatable :: prof(:, :)
        real(dp) :: fine
        integer :: j, i

        if (present(direction)) then
            status = NTV_BAD_INPUT
            if (allocated(sf%line)) return
        end if
        if (sf%source_dlf10) then
            status = NTV_BAD_INPUT
            if (.not. present(f10) .or. present(direction)) return
            if (.not. present(indicators) .or. indicators) return
        end if
        krylov = present(precondition)
        if (krylov) then
            status = NTV_BAD_INPUT
            if (.not. allocated(precondition%sys)) return
            if (precondition%sys%n /= n_theta .or. 2*precondition%sys%nblock /= n_leg) return
        end if
        allocate (a)
        call factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status, f10, krylov)
        if (status /= NTV_OK) return
        call drive(a, v, rhs)
        if (krylov) then
            call gmres_solve(a, precondition%sys, rhs, x, status, .false., kern%iterations)
        else
            call solve(a, rhs, x, status)
        end if
        if (status /= NTV_OK) return
        kern%r = flux(a, x)
        allocate (kern%spectrum(0:2*a%nblock - 1))
        do j = 0, a%nblock - 1
            kern%spectrum(2*j) = maxval(abs(x(1:a%n, j)))
            kern%spectrum(2*j + 1) = maxval(abs(x(a%n + 1:, j)))
        end do

        ! One torque adjoint replaces one tangent solve per parameter.
        call goal(a, rhs)
        if (krylov) then
            call gmres_solve(a, precondition%sys, rhs, adj, status, .true., i)
            kern%iterations = kern%iterations + i
        else
            call solve(a, rhs, adj, status, adjoint=.true.)
        end if
        if (status /= NTV_OK) return
        kern%dr_dnu = -real(sum(conjg(adj)*lorentz(a, x)), dp)/v
        kern%dr_domega = -real(sum(conjg(adj)*a%shift/v*x), dp)

        allocate (kern%dr_dc(size(sf%m)), ds(a%n))
        da = coefficients(a)
        da%s = 1
        ! The unit rhs carries the folding shared by S and the parallel-
        ! derivative drive; the f_10 sources enter per amplitude below.
        da%dpar = .false.
        da%dlf10 = .false.
        call drive(da, v, rhs)
        ! Only S pairs with F in the torque, every drive with the adjoint.
        xw = conjg(conjg(adj(1:a%n, 0))*rhs(1:a%n, 0) &
            + conjg(adj(1:a%n, 1))*rhs(1:a%n, 1))
        gradient = 2*pi/a%n*a%measure*moment(a, x) + xw
        do j = 1, size(sf%m)
            ds = [(exp(iu*sf%m(j)*2*pi*(i - 1)/a%n), i=1, a%n)]
            kern%dr_dc(j) = sum(conjg(iu*sf%n*ds)*gradient)
            if (a%dpar) kern%dr_dc(j) = kern%dr_dc(j) &
                + sum(conjg((iu*(sf%m(j) - sf%n*sf%q) - 2*log_b_slope(a))/sf%q*ds)*xw)
        end do
        if (a%dlf10) then
            ! dH/dc_j carries e^{i m_j theta} ((i (m_j - n q) - d_theta ln b) stream
            ! f10f + f10c)/b; contract the adjoint against the two profiles.
            allocate (prof(2*a%n, 0:a%nblock - 1))
            do j = 0, a%nblock - 1
                do i = 1, a%n
                    prof(i, j) = 1.5_dp*a%source(i)*a%stream(i)*a%f10f(i, j)
                    prof(a%n + i, j) = 1.5_dp*a%source(i)*a%stream(i)*a%f10f(a%n + i, j)
                end do
            end do
            xf = f10_collapse(a, adj, prof)
            do j = 0, a%nblock - 1
                do i = 1, a%n
                    prof(i, j) = 1.5_dp*a%source(i)*a%f10c(i, j)
                    prof(a%n + i, j) = 1.5_dp*a%source(i)*a%f10c(a%n + i, j)
                end do
            end do
            xw = f10_collapse(a, adj, prof)
            do j = 1, size(sf%m)
                ds = [(exp(iu*sf%m(j)*2*pi*(i - 1)/a%n), i=1, a%n)]
                kern%dr_dc(j) = kern%dr_dc(j) &
                    + sum(conjg((iu*(sf%m(j) - sf%n*sf%q) - log_b_slope(a))*ds*xf)) &
                    + sum(conjg(ds*xw))
            end do
        end if

        if (present(direction)) then
            call tangent_system(sf, direction, a, da)
            call drive(da, v, rhs)
            rhs = rhs - apply(da, x)
            kern%dr_dsurface = real(sum(conjg(adj)*rhs), dp) &
                + 2*pi/a%n*sum(da%measure*real(conjg(a%s)*moment(a, x), dp))
        end if
        if (present(keep) .and. .not. krylov) call move_alloc(a, keep%sys)
        if (present(indicators)) then
            if (.not. indicators) return
        end if

        call refined(sf, v, nu_d, omega_e, pre, 2*n_theta + 1, n_leg, fine, status)
        if (status /= NTV_OK) return
        kern%indicator_theta = fine - kern%r
        call refined(sf, v, nu_d, omega_e, pre, n_theta, 2*n_leg, fine, status)
        if (status /= NTV_OK) return
        kern%indicator_leg = fine - kern%r
    end subroutine ntv_kernel

    subroutine refined(sf, v, nu_d, omega_e, pre, n_theta, n_leg, r, status)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        real(dp), intent(out) :: r
        integer, intent(out) :: status
        type(ntv_system_t) :: a
        complex(dp), allocatable :: x(:, :), rhs(:, :)
        call factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status)
        if (status /= NTV_OK) return
        call drive(a, v, rhs)
        call solve(a, rhs, x, status)
        if (status /= NTV_OK) return
        r = flux(a, x)
    end subroutine refined

    subroutine factor(sf, v, nu_d, omega_e, pre, n_theta, n_leg, a, status, f10, assemble_only)
        type(ntv_surface_t), intent(in) :: sf
        real(dp), intent(in) :: v, nu_d, omega_e, pre
        integer, intent(in) :: n_theta, n_leg
        type(ntv_system_t), intent(out) :: a
        integer, intent(out) :: status
        real(dp), intent(in), optional, dimension(:, :) :: f10
        logical, intent(in), optional :: assemble_only
        real(dp), allocatable :: theta(:), b(:), stream_b(:), mirror_b(:), measure_b(:)
        complex(dp), allocatable :: m(:, :), upper(:, :)
        integer :: n, j, i, info

        status = NTV_BAD_INPUT
        if (mod(n_theta, 2) == 0 .or. n_theta < 3 .or. mod(n_leg, 2) /= 0 .or. &
            n_leg < 4 .or. v <= 0 .or. nu_d <= 0 .or. sf%n == 0 .or. &
            .not. allocated(sf%m) .or. &
            (sf%weight /= NTV_WEIGHT_MU .and. sf%weight /= NTV_WEIGHT_FULL) .or. &
            (sf%source_dpar .and. sf%q == 0)) return
        ! Nested: Fortran need not short-circuit, and f10 may be absent.
        if (sf%source_dlf10) then
            if (.not. present(f10)) return
            if (size(f10, 1) /= 2*n_theta .or. size(f10, 2) /= n_leg/2) return
            if (any(ieee_is_nan(f10))) return
        end if
        n = n_theta
        a%n = n
        a%nblock = n_leg/2
        a%nu = nu_d/v
        a%shift = iu*sf%n/v
        a%nq = sf%n*sf%q
        a%omega = omega_e
        a%precess = pre
        allocate (theta(n), b(n), a%stream(n), a%mirror(n), a%drift(n), &
            a%drift2(n), a%source(n), a%measure(n), a%s(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        allocate (a%d(n, n))
        call fourier_derivative(n, a%d)
        if (allocated(sf%line)) then
            call sf%line%sample(theta, a%stream, a%mirror, a%drift, a%drift2, &
                a%source, a%measure)
        else
            ! Streaming along the field line in the frame exp(-i n zeta):
            ! d_theta - i n q.
            a%stream = sf%k
            allocate (stream_b(n), mirror_b(n), measure_b(n))
            call dk_ntv_point(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, b, &
                a%mirror, a%drift, a%source, a%measure, stream_b, mirror_b, measure_b)
            if (sf%boozer_stream) then
                a%stream = stream_b
                a%mirror = mirror_b
                a%measure = measure_b
            end if
            a%drift2 = sf%current
        end if
        if (sf%weight == NTV_WEIGHT_FULL) a%wsign = 1
        a%s = 0
        do j = 1, size(sf%m)
            a%s = a%s + iu*sf%n*sf%c(j)*exp(iu*sf%m(j)*theta)
        end do
        if (sf%source_dpar) then
            a%dpar = .true.
            allocate (a%sd(n))
            a%sd = 0
            do j = 1, size(sf%m)
                a%sd = a%sd + (iu*(sf%m(j) - sf%n*sf%q) - 2*log_b_slope(a))/sf%q &
                    *sf%c(j)*exp(iu*sf%m(j)*theta)
            end do
        end if
        if (sf%source_dlf10) then
            a%dlf10 = .true.
            allocate (a%sdb(n), a%db(n), a%f10(2*n, 0:a%nblock - 1), &
                a%f10f(2*n, 0:a%nblock - 1), a%f10c(2*n, 0:a%nblock - 1))
            a%sdb = 0
            a%db = 0
            do j = 1, size(sf%m)
                a%sdb = a%sdb + (iu*(sf%m(j) - sf%n*sf%q) - log_b_slope(a))*sf%c(j) &
                    *exp(iu*sf%m(j)*theta)
                a%db = a%db + sf%c(j)*exp(iu*sf%m(j)*theta)
            end do
            a%f10 = f10
            call f10_foldings(a)
        end if
        if (present(assemble_only)) then
            if (assemble_only) then
                status = NTV_OK
                return
            end if
        end if

        ! S_J = Di_J; S_j = Di_j - Up_j Y_{j+1}, Y_j = S_j^-1 Lo_j.
        allocate (a%s_lu(2*n, 2*n, 0:a%nblock - 1), a%y(2*n, 2*n, a%nblock - 1), &
            a%piv(2*n, 0:a%nblock - 1), m(2*n, 2*n), upper(2*n, 2*n))
        status = NTV_SINGULAR
        do j = a%nblock - 1, 0, -1
            m = block(a, j, j)
            if (j < a%nblock - 1) then
                upper = block(a, j, j + 1)
                call zgemm('N', 'N', 2*n, 2*n, 2*n, (-1.0_dp, 0.0_dp), upper, 2*n, &
                    a%y(:, :, j + 1), 2*n, (1.0_dp, 0.0_dp), m, 2*n)
            end if
            call zgetrf(2*n, 2*n, m, 2*n, a%piv(:, j), info)
            if (info /= 0) return
            a%s_lu(:, :, j) = m
            if (j > 0) then
                a%y(:, :, j) = block(a, j, j - 1)
                call zgetrs('N', 2*n, 2*n, m, 2*n, a%piv(:, j), a%y(:, :, j), 2*n, info)
            end if
        end do
        status = NTV_OK
    end subroutine factor

    !> d_theta ln b on the theta points from the mirror factor
    !> stream d_theta(ln b)/2 of the line contract.
    pure function log_b_slope(a) result(g)
        type(ntv_system_t), intent(in) :: a
        real(dp) :: g(a%n)
        g = 2*a%mirror/a%stream
    end function log_b_slope

    !> The two operator foldings of the axisymmetric solution: f10f is
    !> (1 - xi^2)/2 d_xi f_10 per Legendre degree from
    !> (1 - xi^2) P'_l = l(l + 1)/(2l + 1)(P_{l - 1} - P_{l + 1}), f10c is
    !> f_10 scaled by nu l(l + 1)/2, the Lorentz rows of the operator.
    subroutine f10_foldings(a)
        type(ntv_system_t), intent(inout) :: a
        integer :: n, j, i, l
        n = a%n
        a%f10f = 0
        a%f10c = 0
        do j = 0, a%nblock - 1
            l = 2*j
            do i = 1, n
                a%f10f(i, j) = (l + 1)*(l + 2)/real(2*(2*l + 3), dp)*a%f10(n + i, j)
                if (j > 0) a%f10f(i, j) = a%f10f(i, j) &
                    - (l - 1)*l/real(2*(2*l - 1), dp)*a%f10(n + i, j - 1)
                if (j < a%nblock - 1) a%f10f(n + i, j) = a%f10f(n + i, j) &
                    + (l + 2)*(l + 3)/real(2*(2*l + 5), dp)*a%f10(i, j + 1)
                a%f10f(n + i, j) = a%f10f(n + i, j) &
                    - l*(l + 1)/real(2*(2*l + 1), dp)*a%f10(i, j)
                a%f10c(i, j) = a%nu*l*(l + 1)/2.0_dp*a%f10(i, j)
                a%f10c(n + i, j) = a%nu*(l + 1)*(l + 2)/2.0_dp*a%f10(n + i, j)
            end do
        end do
    end subroutine f10_foldings

    !> Per theta point i, sum_j conjg(adj(rows, j)) profile(rows, j) over
    !> both parities of every block: the adjoint contraction of a response
    !> profile that occupies all degrees.
    function f10_collapse(a, adj, profile) result(x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: adj(:, 0:)
        real(dp), intent(in) :: profile(:, 0:)
        complex(dp) :: x(a%n)
        complex(dp) :: acc
        integer :: i, j
        do i = 1, a%n
            acc = 0
            do j = 0, a%nblock - 1
                acc = acc + conjg(adj(i, j))*profile(i, j) &
                    + conjg(adj(a%n + i, j))*profile(a%n + i, j)
            end do
            x(i) = acc
        end do
    end function f10_collapse

    !> Solve A x = rhs, pair block j holding l = 2j (rows 1:n) and 2j + 1.
    subroutine solve(a, rhs, x, status, adjoint)
        type(ntv_system_t), intent(inout) :: a
        complex(dp), intent(in) :: rhs(:, 0:)
        complex(dp), allocatable, intent(out) :: x(:, :)
        integer, intent(out) :: status
        complex(dp), allocatable :: z(:, :)
        integer :: j, info, nn
        logical, intent(in), optional :: adjoint
        logical :: reverse
        reverse = .false.
        if (present(adjoint)) reverse = adjoint
        nn = 2*a%n
        allocate (x(nn, 0:a%nblock - 1))
        allocate (z(nn, 0:a%nblock - 1))
        status = NTV_SINGULAR
        if (reverse) then
            ! A = U L for the retained descending block elimination.
            z = rhs
            do j = a%nblock - 2, 0, -1
                z(:, j) = z(:, j) - conjg(matmul(conjg(z(:, j + 1)), a%y(:, :, j + 1)))
            end do
            do j = 0, a%nblock - 1
                x(:, j) = z(:, j)
                if (j > 0) x(:, j) = x(:, j) - block_apply(a, j - 1, j, x(:, j - 1), .true.)
                call zgetrs('C', nn, 1, a%s_lu(:, :, j), nn, a%piv(:, j), x(:, j), nn, info)
                if (info /= 0) return
            end do
            status = NTV_OK
            return
        end if
        do j = a%nblock - 1, 0, -1
            z(:, j) = rhs(:, j)
            if (j < a%nblock - 1) z(:, j) = z(:, j) &
                - block_apply(a, j, j + 1, z(:, j + 1), .false.)
            call zgetrs('N', nn, 1, a%s_lu(:, :, j), nn, a%piv(:, j), z(:, j), nn, info)
            if (info /= 0) return
        end do
        x(:, 0) = z(:, 0)
        do j = 1, a%nblock - 1
            x(:, j) = z(:, j) - matmul(a%y(:, :, j), x(:, j - 1))
        end do
        status = NTV_OK
    end subroutine solve

    !> Solve A x = rhs, or A^H x = rhs, by restarted GMRES on A M^-1 with
    !> M the factorised operator of a nearby kernel (M^-H for the adjoint).
    subroutine gmres_solve(a, m, rhs, x, status, reverse, iterations)
        type(ntv_system_t), intent(inout) :: a, m
        complex(dp), intent(in) :: rhs(:, 0:)
        complex(dp), allocatable, intent(out) :: x(:, :)
        integer, intent(out) :: status, iterations
        logical, intent(in) :: reverse
        complex(dp), allocatable :: b(:), u(:), w(:)
        real(dp) :: residual, scale
        integer :: info, nn, nb
        nn = 2*a%n
        nb = a%nblock
        b = reshape(rhs, [nn*nb])
        scale = norm2(abs(b))
        allocate (u(nn*nb), w(nn*nb), x(nn, 0:nb - 1))
        x = 0
        iterations = 0
        status = NTV_OK
        if (scale == 0) return
        ! Unit right-hand side: the GMRES target is absolute below norm one.
        b = b/scale
        u = 0
        call complex_gmres_operator(operator, b, u, gmres_tolerance, 400, 100, info, &
            iterations, residual)
        if (status /= NTV_OK) return
        status = NTV_SINGULAR
        if (info /= KRYLOV_OK) return
        status = NTV_OK
        call precondition(u, w)
        if (status /= NTV_OK) return
        x = scale*reshape(w, [nn, nb])
    contains
        subroutine precondition(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            complex(dp), allocatable :: y(:, :)
            integer :: st
            call solve(m, reshape(input, [nn, nb]), y, st, adjoint=reverse)
            if (st /= NTV_OK) status = st
            output = reshape(y, [nn*nb])
        end subroutine precondition
        subroutine operator(input, output)
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
            complex(dp), allocatable :: t(:)
            allocate (t(size(input)))
            call precondition(input, t)
            if (reverse) then
                output = reshape(apply_adjoint(a, reshape(t, [nn, nb])), [nn*nb])
            else
                output = reshape(apply(a, reshape(t, [nn, nb])), [nn*nb])
            end if
        end subroutine operator
    end subroutine gmres_solve

    !> Rows of pair block j, columns of pair block jc: the couplings of
    !> rows l to columns l' with |l - l'| <= 2.
    function block(a, j, jc) result(m)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: j, jc
        complex(dp) :: m(2*a%n, 2*a%n)
        integer :: r, c, l, lc
        m = 0
        do r = 0, 1
            do c = 0, 1
                l = 2*j + r
                lc = 2*jc + c
                m(r*a%n + 1:(r + 1)*a%n, c*a%n + 1:(c + 1)*a%n) = entry(a, l, lc)
            end do
        end do
    end function block

    !> Coefficient matrix of x_lc in row l.
    function entry(a, l, lc) result(m)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: l, lc
        complex(dp) :: m(a%n, a%n)
        integer :: i
        complex(dp) :: fd, diag
        m = 0
        if (abs(l - lc) > 2) return
        do i = 1, a%n
            call coupling(a, l, lc, i, fd, diag)
            if (fd /= 0) m(i, :) = fd*a%d(i, :)
            m(i, i) = m(i, i) + diag
        end do
    end function entry

    !> Row (l, i), column block lc: the factor of the theta derivative row
    !> d(i, :) and the diagonal coefficient.
    pure subroutine coupling(a, l, lc, i, fd, diag)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: l, lc, i
        complex(dp), intent(out) :: fd, diag
        real(dp) :: f
        complex(dp) :: p, pj
        fd = 0
        diag = 0
        p = a%shift*a%precess*a%drift(i)
        pj = a%shift*a%precess*a%drift2(i)
        select case (lc - l)
        case (-1)
            f = l/real(2*l - 1, dp)
            fd = a%stream(i)*f
            diag = -a%stream(i)*f*iu*a%nq + a%mirror(i)*(l - 1)*l/real(2*l - 1, dp)
        case (1)
            f = (l + 1)/real(2*l + 3, dp)
            fd = a%stream(i)*f
            diag = -a%stream(i)*f*iu*a%nq - a%mirror(i)*(l + 1)*(l + 2)/real(2*l + 3, dp)
        case (0)
            diag = a%nu*l*(l + 1)/2.0_dp + a%shift*a%omega &
                + p*(1 + a%wsign*c0(l)) + pj*c0(l)
        case (-2)
            diag = (a%wsign*p + pj)*cp(lc)
        case (2)
            diag = (a%wsign*p + pj)*cm(lc)
        end select
    end subroutine coupling

    ! xi^2 P_l = cp P_{l+2} + c0 P_l + cm P_{l-2}, so (1 + wsign xi^2) P_l
    ! has wsign cp, 1 + wsign c0 and wsign cm.
    pure real(dp) function cp(l)
        integer, intent(in) :: l
        cp = (l + 1)*(l + 2)/real((2*l + 1)*(2*l + 3), dp)
    end function cp

    pure real(dp) function c0(l)
        integer, intent(in) :: l
        c0 = (l + 1)**2/real((2*l + 1)*(2*l + 3), dp)
        if (l > 0) c0 = c0 + l**2/real((2*l + 1)*(2*l - 1), dp)
    end function c0

    pure real(dp) function cm(l)
        integer, intent(in) :: l
        cm = 0
        if (l > 1) cm = l*(l - 1)/real((2*l + 1)*(2*l - 1), dp)
    end function cm

    function lorentz(a, x) result(y)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: j, l
        do j = 0, ubound(x, 2)
            l = 2*j
            y(1:a%n, j) = l*(l + 1)/2.0_dp*x(1:a%n, j)
            y(a%n + 1:, j) = (l + 1)*(l + 2)/2.0_dp*x(a%n + 1:, j)
        end do
    end function lorentz

    subroutine circular_sample(self, theta, stream, mirror, drift, drift2, source, &
            measure)
        class(ntv_circular_line_t), intent(in) :: self
        real(dp), intent(in) :: theta(:)
        real(dp), intent(out) :: stream(:), mirror(:), drift(:), drift2(:), &
            source(:), measure(:)
        real(dp) :: b(size(theta)), stream_b(size(theta)), mirror_b(size(theta)), &
            measure_b(size(theta))
        stream = self%k
        call dk_ntv_point(cos(theta), sin(theta), self%eps, self%k, self%offset, b, &
            mirror, drift, source, measure, stream_b, mirror_b, measure_b)
        if (self%boozer_stream) then
            stream = stream_b
            mirror = mirror_b
            measure = measure_b
        end if
        drift2 = self%current
    end subroutine circular_sample

    !> <ds, x>_W: the flux functional with the conjugated drive replaced.
    real(dp) function pair(a, ds, x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: ds(:), x(:, 0:)
        pair = 2*pi/a%n*sum(a%measure*real(conjg(ds)*moment(a, x), dp))
    end function pair

    pure function moment_weights(a) result(w)
        type(ntv_system_t), intent(in) :: a
        real(dp) :: w(2)
        w = [2 + 2*a%wsign/3, 4*a%wsign/15]
    end function moment_weights

    !> int (1 + wsign xi^2) F dxi per theta point.
    function moment(a, x) result(mom)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: mom(a%n)
        real(dp) :: w(2)
        w = moment_weights(a)
        mom = w(1)*x(1:a%n, 0)
        if (ubound(x, 2) >= 1) mom = mom + w(2)*x(1:a%n, 1)
    end function moment

    subroutine goal(a, rhs)
        type(ntv_system_t), intent(in) :: a
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        real(dp) :: w(2)
        w = moment_weights(a)
        allocate (rhs(2*a%n, 0:a%nblock - 1))
        rhs = 0
        rhs(1:a%n, 0) = 2*pi/a%n*a%measure*a%s*w(1)
        rhs(1:a%n, 1) = 2*pi/a%n*a%measure*a%s*w(2)
    end subroutine goal

    !> Tangent of the operator and the pointwise leaves along a direction in
    !> (eps, k, offset, current): da holds dA (no collision or E x B part,
    !> which do not depend on these), dsource and dmeasure.
    subroutine tangent_system(sf, dir, a, da)
        type(ntv_surface_t), intent(in) :: sf, dir
        type(ntv_system_t), intent(in) :: a
        type(ntv_system_t), intent(out) :: da
        real(dp), allocatable :: theta(:), db(:), dstream_b(:), dmirror_b(:), &
            dmeasure_b(:)
        integer :: i, n
        n = a%n
        da = coefficients(a)
        da%nu = 0
        da%omega = 0
        da%drift2 = dir%current
        allocate (theta(n), db(n), dstream_b(n), dmirror_b(n), dmeasure_b(n))
        theta = [(2*pi*(i - 1)/n, i=1, n)]
        call dk_ntv_point_jvp(cos(theta), sin(theta), sf%eps, sf%k, sf%offset, &
            0*theta, 0*theta, dir%eps + 0*theta, dir%k + 0*theta, dir%offset + 0*theta, &
            db, da%mirror, da%drift, da%source, da%measure, dstream_b, dmirror_b, &
            dmeasure_b)
        da%stream = dir%k
        if (sf%boozer_stream) then
            da%stream = dstream_b
            da%mirror = dmirror_b
            da%measure = dmeasure_b
        end if
    end subroutine tangent_system

    !> The operator coefficients of a without its factors.
    function coefficients(a) result(c)
        type(ntv_system_t), intent(in) :: a
        type(ntv_system_t) :: c
        c%n = a%n; c%nblock = a%nblock; c%wsign = a%wsign; c%nu = a%nu; c%nq = a%nq
        c%drift2 = a%drift2; c%shift = a%shift; c%omega = a%omega
        c%precess = a%precess
        c%d = a%d; c%stream = a%stream; c%mirror = a%mirror; c%drift = a%drift
        c%source = a%source; c%measure = a%measure; c%s = a%s
        c%dpar = a%dpar; c%dlf10 = a%dlf10
        if (a%dpar) c%sd = a%sd
        if (a%dlf10) then
            c%sdb = a%sdb; c%db = a%db
            c%f10 = a%f10; c%f10f = a%f10f; c%f10c = a%f10c
        end if
    end function coefficients

    !> y = A x for the pair-block tridiagonal operator.
    function apply(a, x) result(y)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: j
        do j = 0, a%nblock - 1
            y(:, j) = block_apply(a, j, j, x(:, j), .false.)
            if (j > 0) y(:, j) = y(:, j) + block_apply(a, j, j - 1, x(:, j - 1), .false.)
            if (j < a%nblock - 1) y(:, j) = y(:, j) &
                + block_apply(a, j, j + 1, x(:, j + 1), .false.)
        end do
    end function apply

    !> y = A^H x.
    function apply_adjoint(a, x) result(y)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        complex(dp) :: y(size(x, 1), 0:ubound(x, 2))
        integer :: j
        do j = 0, a%nblock - 1
            y(:, j) = block_apply(a, j, j, x(:, j), .true.)
            if (j > 0) y(:, j) = y(:, j) + block_apply(a, j - 1, j, x(:, j - 1), .true.)
            if (j < a%nblock - 1) y(:, j) = y(:, j) &
                + block_apply(a, j + 1, j, x(:, j + 1), .true.)
        end do
    end function apply_adjoint

    !> Block (j, jc) times x, or its adjoint times x (x then in block j and
    !> the result in block jc), from the couplings without forming the
    !> block: each quadrant is diag(fd) D + diag(diag).
    function block_apply(a, j, jc, x, adjoint) result(y)
        type(ntv_system_t), intent(in) :: a
        integer, intent(in) :: j, jc
        complex(dp), intent(in) :: x(:)
        logical, intent(in) :: adjoint
        complex(dp) :: y(2*a%n)
        complex(dp) :: fd(a%n), dg(a%n), xs(a%n)
        integer :: r, c, l, lc, i, n
        n = a%n
        y = 0
        do c = 0, 1
            do r = 0, 1
                l = 2*j + r
                lc = 2*jc + c
                if (abs(l - lc) > 2) cycle
                do i = 1, n
                    call coupling(a, l, lc, i, fd(i), dg(i))
                end do
                if (adjoint) then
                    xs = x(r*n + 1:(r + 1)*n)
                    y(c*n + 1:(c + 1)*n) = y(c*n + 1:(c + 1)*n) + conjg(dg)*xs
                    if (any(fd /= 0)) y(c*n + 1:(c + 1)*n) = y(c*n + 1:(c + 1)*n) &
                        + matmul(conjg(fd)*xs, a%d)
                else
                    xs = x(c*n + 1:(c + 1)*n)
                    y(r*n + 1:(r + 1)*n) = y(r*n + 1:(r + 1)*n) + dg*xs
                    if (any(fd /= 0)) y(r*n + 1:(r + 1)*n) = y(r*n + 1:(r + 1)*n) &
                        + fd*matmul(a%d, xs)
                end if
            end do
        end do
    end function block_apply

    !> R = (2 pi/n) sum (2/b^2) Re(S^* int (1 + wsign xi^2) F dxi), with
    !> int (1 + wsign xi^2) F dxi = (2 + 2 wsign/3) a_0 + (4 wsign/15) a_2.
    real(dp) function flux(a, x)
        type(ntv_system_t), intent(in) :: a
        complex(dp), intent(in) :: x(:, 0:)
        flux = pair(a, a%s, x)
    end function flux

    !> Drive (1 + wsign xi^2) S/(b v): source = 2/(3 b) per unit of the
    !> Legendre components 1 + wsign/3 (P_0) and 2 wsign/3 (P_2).
    subroutine drive(a, v, rhs)
        type(ntv_system_t), intent(in) :: a
        real(dp), intent(in) :: v
        complex(dp), allocatable, intent(out) :: rhs(:, :)
        integer :: i, j
        allocate (rhs(2*a%n, 0:a%nblock - 1))
        rhs = 0
        rhs(1:a%n, 0) = (1.5_dp + a%wsign/2)*a%source*a%s/v
        rhs(1:a%n, 1) = a%wsign*a%source*a%s/v
        if (a%dpar) then
            rhs(1:a%n, 0) = rhs(1:a%n, 0) + (1.5_dp + a%wsign/2)*a%source*a%sd/v
            rhs(1:a%n, 1) = rhs(1:a%n, 1) + a%wsign*a%source*a%sd/v
        end if
        if (a%dlf10) then
            ! 1/b = 3 source/2 by the provider contract of source = 2/(3 b);
            ! f10c holds -C[f_10], so -db/b C[f_10] enters with a plus.
            do j = 0, a%nblock - 1
                do i = 1, a%n
                    rhs(i, j) = rhs(i, j) + 1.5_dp*a%source(i)*a%stream(i)*a%sdb(i) &
                        *a%f10f(i, j) + 1.5_dp*a%source(i)*a%db(i)*a%f10c(i, j)
                    rhs(a%n + i, j) = rhs(a%n + i, j) &
                        + 1.5_dp*a%source(i)*a%stream(i)*a%sdb(i)*a%f10f(a%n + i, j) &
                        + 1.5_dp*a%source(i)*a%db(i)*a%f10c(a%n + i, j)
                end do
            end do
        end if
    end subroutine drive

    include 'dk_ntv_energy.inc'

end module milca_dk_ntv
