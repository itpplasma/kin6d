! Two-body Reid68 impulse-approximation pp fusion, MeV/fm/seconds.
module kin6d_pp_fusion
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_radial_dvr, only: radial_mesh_t, build_radial_mesh, &
        radial_band_hamiltonian, radial_coupled_hamiltonian, &
        radial_band_eigen_window, radial_band_matvec, radial_bound_response
    use kin6d_pp_fusion_leaves, only: pp_reid, pp_sfactor, &
        pp_lambda_integrand, pp_gv, pp_phase_quad, pp_phase_quad_jvp, &
        pp_match, pp_match_jvp, pp_match_norm, pp_match_norm_jvp, &
        pp_effective_range, pp_effective_range_jvp, pp_overlap, pp_overlap_jvp, &
        pp_sfactor_jvp
    use fortnum_special_coulomb, only: coulomb_fg0, coulomb_c0_squared, coulomb_h
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    private
    ! CODATA 2022 (Rev.Mod.Phys.97,025002); weak inputs frozen to Solar Fusion II.
    real(dp), parameter, public :: pp_hc = 197.3269804_dp
    real(dp), parameter, public :: pp_hbar = 6.582119569e-22_dp
    real(dp), parameter, public :: pp_alpha = 7.2973525643e-3_dp
    real(dp), parameter, public :: pp_mp = 938.27208943_dp
    real(dp), parameter, public :: pp_mn = 939.56542194_dp
    real(dp), parameter, public :: pp_me = 0.51099895069_dp
    real(dp), parameter, public :: axial = 1.2695_dp, ft0 = 3071.4_dp
    real(dp), parameter :: outer_radiative = 1.0162_dp
    ! Reid's fitted kinetic coefficient. Using modern masses instead shifts B.
    real(dp), parameter, public :: reid_kinetic = 41.47_dp
    type, public :: pp_fusion_config_t
        real(dp) :: rmax = 100, hmax = 0.8_dp, h0 = 0.05_dp
        real(dp) :: energy = 0.001_dp ! representative 1 keV c.m.
        real(dp) :: strength(2) = [1.0_dp, 1.0_dp] ! singlet/triplet nuclear scales
        integer :: order = 8
    end type
    type, public :: pp_fusion_result_t
        type(radial_mesh_t) :: mesh
        real(dp) :: binding = 0, d_probability = 0, a_pp = 0
        real(dp) :: lambda2 = 0, sfactor = 0, slope = 0, phase_space = 0
        real(dp) :: phase_shift = 0, energy = 0, eigen_residual = 0
        real(dp) :: scattering_residual = 0, gamma = 0
        real(dp) :: energy_extrapolation = 0
        real(dp), allocatable :: us(:), ud(:), upp(:), integrand(:)
        real(dp), allocatable :: potentials(:, :)
        character(32) :: status = 'candidate'
        character(48) :: numerical_bound = 'unavailable; convergence indicators only'
        character(48) :: model_bound = 'unavailable; one-body Reid68 model'
    end type
    public :: pp_fusion_solve, pp_scattering, pp_phase_space, pp_vector_coupling
    public :: pp_fusion_jvp, pp_fusion_vjp
    interface
        subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
            real(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine
    end interface
contains
    ! Weak vector coupling inferred from the corrected superallowed ft value.
    ! It includes the inner radiative convention, avoiding its double counting.
    function pp_vector_coupling() result(gv)
        real(dp) :: gv
        call pp_gv(ft0, pp_me, pp_hbar, gv)
    end function

    function pp_phase_space(energy, binding, dbinding) result(f)
        real(dp), intent(in) :: energy, binding
        real(dp), intent(out), optional :: dbinding
        real(dp) :: f, x(128), w(128), endpoint, term, dterm, df
        integer :: i
        endpoint = (pp_mp - pp_mn + binding + energy)/pp_me
        f = 0; df = 0
        if (present(dbinding)) dbinding = 0
        if (endpoint <= 1) return
        call gauss_legendre(size(x), x, w)
        do i = 1, size(x)
            call pp_phase_quad(x(i), endpoint, pp_alpha, term)
            f = f + w(i)*term
            if (.not. present(dbinding)) cycle
            call pp_phase_quad_jvp(x(i), endpoint, pp_alpha, &
                0.0_dp, 1/pp_me, 0.0_dp, dterm)
            df = df + w(i)*dterm
        end do
        f = f*outer_radiative
        if (present(dbinding)) dbinding = df*outer_radiative
    end function

    subroutine pp_fusion_solve(cfg, result, ok, jacobian)
        type(pp_fusion_config_t), intent(in) :: cfg
        type(pp_fusion_result_t), intent(out) :: result
        logical, intent(out) :: ok
        real(dp), intent(out), optional :: jacobian(5, 2)
        type(radial_mesh_t) :: bound_mesh
        real(dp), allocatable :: pot(:, :, :), hb(:, :), e(:), vec(:, :), res(:)
        real(dp), allocatable :: dhb(:, :), dv(:), dus(:, :), dwave(:, :)
        real(dp) :: de, db(2), dgam(2), dlam(3, 2), dkcm(3, 2), dphase, df, ov, dov
        real(dp) :: dl0(2), dummy
        integer :: axis
        real(dp), allocatable :: wave(:, :), phases(:), er(:), effective(:)
        real(dp) :: orth, energies(3), l2(3), sfe(3), vs, vss, vsd, vdd, sign_u
        real(dp) :: kval, eta, c2, kcm, kinetic_pp
        integer :: i, j, n
        logical :: good(3)
        ok = .false.
        if (.not. all(ieee_is_finite([cfg%rmax, cfg%hmax, cfg%h0, &
            cfg%energy, cfg%strength]))) return
        if (any(cfg%strength < 0.95_dp) .or. any(cfg%strength > 1.05_dp)) return
        if (cfg%energy <= 0 .or. cfg%energy > 0.01_dp) return
        if (cfg%rmax < 30 .or. cfg%rmax > 100) return
        call build_radial_mesh(bound_mesh, cfg%rmax, cfg%order, cfg%hmax, &
            cfg%h0, 1.15_dp, ok)
        if (.not. ok) return
        ok = .false.
        n = bound_mesh%n
        allocate (pot(2, 2, n), hb(2*cfg%order - 1, 2*n))
        !$omp parallel do private(vs,vss,vsd,vdd) schedule(static)
        do i = 1, n
            call pp_reid(bound_mesh%r(i), reid_kinetic, vs, vss, vsd, vdd)
            pot(1, 1, i) = cfg%strength(2)*vss
            pot(2, 2, i) = cfg%strength(2)*(vdd - 6*reid_kinetic/bound_mesh%r(i)**2) &
                + 6*reid_kinetic/bound_mesh%r(i)**2
            pot(1, 2, i) = cfg%strength(2)*vsd; pot(2, 1, i) = cfg%strength(2)*vsd
        end do
        !$omp end parallel do
        call radial_coupled_hamiltonian(bound_mesh, reid_kinetic, pot, hb)
        call radial_band_eigen_window(hb, 0.0_dp, e, vec, res, orth, ok, emin=-100.0_dp)
        if (.not. ok) return
        ok = .false.
        if (size(e) /= 1) return
        result%binding = -e(1); result%eigen_residual = res(1)
        result%d_probability = sum(vec(2::2, 1)**2)
        result%gamma = sqrt(result%binding/reid_kinetic)
        sign_u = sign(1.0_dp, vec(2*minloc(abs(bound_mesh%r - 3), dim=1) - 1, 1))
        allocate (dus(n + 1, 2))
        dus = 0; db = 0; dgam = 0
        if (present(jacobian)) then
            jacobian = 0
            allocate (dhb(size(hb, 1), size(hb, 2)), dv(2*n))
            ! dH/d strength_triplet contains nuclear potential only.
            do i = 1, n
                pot(2, 2, i) = pot(2, 2, i) - 6*reid_kinetic/bound_mesh%r(i)**2
            end do
            pot = pot/cfg%strength(2)
            call radial_coupled_hamiltonian(bound_mesh, 0.0_dp, pot, dhb)
            call radial_bound_response(hb, e(1), vec(:, 1), dhb, de, dv, ok)
            if (.not. ok) return
            ok = .false.
            db(2) = -de; dgam(2) = db(2)/(2*reid_kinetic*result%gamma)
            dus = 0
            dus(:n, 2) = sign_u*dv(1::2)/sqrt(bound_mesh%w)
            jacobian(1, 2) = db(2)
            jacobian(2, 2) = 2*dot_product(vec(2::2, 1), dv(2::2))
        end if
        ! The scattering mesh retains r=R. Its other nodes/weights coincide
        ! with the bound mesh; u_d(R)=0 is its declared Dirichlet condition.
        call build_radial_mesh(result%mesh, cfg%rmax, cfg%order, cfg%hmax, &
            cfg%h0, 1.15_dp, ok, retain_outer=.true.)
        if (.not. ok) return
        ok = .false.
        n = result%mesh%n
        allocate (result%us(n), result%ud(n), result%potentials(4, n))
        result%us = 0; result%ud = 0
        result%us(:n - 1) = sign_u*vec(1::2, 1)/sqrt(bound_mesh%w)
        result%ud(:n - 1) = sign_u*vec(2::2, 1)/sqrt(bound_mesh%w)
        !$omp parallel do private(vs,vss,vsd,vdd) schedule(static)
        do i = 1, n
            call pp_reid(result%mesh%r(i), reid_kinetic, vs, vss, vsd, vdd)
            result%potentials(:, i) = [cfg%strength(1)*vs + &
                pp_alpha*pp_hc/result%mesh%r(i), cfg%strength(2)*vss, &
                cfg%strength(2)*vsd, &
                cfg%strength(2)*(vdd - 6*reid_kinetic/result%mesh%r(i)**2)]
        end do
        !$omp end parallel do
        ! Quadratic E->0 extrapolation, with the two-point change retained.
        energies = [0.0001_dp, 0.0002_dp, 0.0004_dp]
        allocate (wave(n, 3), phases(3), er(3), effective(3))
        kinetic_pp = pp_hc**2/pp_mp
        !$omp parallel do private(j,kval,eta,c2,kcm) schedule(static)
        do j = 1, 3
            call pp_scattering(result%mesh, energies(j), wave(:, j), phases(j), &
                effective(j), er(j), good(j), strength=cfg%strength(1))
            if (.not. good(j)) cycle
            kval = sqrt(energies(j)/kinetic_pp)
            eta = pp_alpha*pp_hc/(2*kinetic_pp*kval)
            c2 = coulomb_c0_squared(eta)
            call pp_overlap(result%gamma, &
                sum(result%mesh%w*result%us*wave(:, j)), c2, kval, l2(j))
            call pp_sfactor(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, result%gamma, &
                l2(j), pp_phase_space(energies(j), result%binding), ft0, sfe(j))
        end do
        !$omp end parallel do
        if (.not. all(good)) return
        if (present(jacobian)) then
            allocate (dwave(n, 1))
            dlam = 0; dkcm = 0
            do j = 1, 3
                call pp_scattering(result%mesh, energies(j), wave(:, j), phases(j), &
                    effective(j), er(j), good(j), cfg%strength(1), dwave(:, 1), &
                    dphase, dkcm(j, 1))
                if (.not. good(j)) return
                kval = sqrt(energies(j)/kinetic_pp)
                eta = pp_alpha*pp_hc/(2*kinetic_pp*kval)
                c2 = coulomb_c0_squared(eta)
                ov = sum(result%mesh%w*result%us*wave(:, j))
                do axis = 1, 2
                    dov = sum(result%mesh%w*dus(:, axis)*wave(:, j))
                    if (axis == 1) dov = dov + sum(result%mesh%w*result%us*dwave(:, 1))
                    call pp_overlap_jvp(result%gamma, ov, c2, kval, &
                        dgam(axis), dov, 0.0_dp, 0.0_dp, dlam(j, axis))
                end do
            end do
        end if
        result%lambda2 = (8*l2(1) - 6*l2(2) + l2(3))/3
        result%energy_extrapolation = abs(result%lambda2 - (2*l2(1) - l2(2)))
        kcm = (8*effective(1) - 6*effective(2) + effective(3))/3
        result%a_pp = -1/kcm
        result%phase_space = pp_phase_space(0.0_dp, result%binding)
        call pp_sfactor(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, result%gamma, &
            result%lambda2, result%phase_space, ft0, result%sfactor)
        if (present(jacobian)) then
            dl0 = (8*dlam(1, :) - 6*dlam(2, :) + dlam(3, :))/3
            jacobian(3, :) = (8*dkcm(1, :) - 6*dkcm(2, :) + dkcm(3, :))/ &
                (3*kcm**2)
            jacobian(4, :) = dl0
            dummy = pp_phase_space(0.0_dp, result%binding, df)
            do axis = 1, 2
                call pp_sfactor_jvp(pp_mp, pp_alpha, pp_hbar, pp_hc, axial, &
                    result%gamma, result%lambda2, result%phase_space, ft0, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, dgam(axis), &
                    dl0(axis), df*db(axis), 0.0_dp, jacobian(5, axis))
            end do
        end if
        result%slope = (-2*sfe(1) + 2.5_dp*sfe(2) - 0.5_dp*sfe(3))/ &
            (energies(1)*result%sfactor)
        allocate (result%upp(n), result%integrand(n))
        call pp_scattering(result%mesh, cfg%energy, result%upp, result%phase_shift, &
            kcm, result%scattering_residual, ok, strength=cfg%strength(1))
        if (.not. ok) return
        result%energy = cfg%energy
        kval = sqrt(cfg%energy/kinetic_pp)
        eta = pp_alpha*pp_hc/(2*kinetic_pp*kval)
        c2 = coulomb_c0_squared(eta)
        do i = 1, n
            call pp_lambda_integrand(result%gamma, result%us(i), &
                result%upp(i)/(sqrt(c2)*kval), result%integrand(i))
        end do
        ok = all(ieee_is_finite([result%binding, result%d_probability, result%a_pp, &
            result%lambda2, result%sfactor, result%slope]))
    end subroutine

    ! Fixed-energy regular solution. A prescribed u(R)=1 supplies a linear
    ! Dirichlet solve; its weak boundary residual gives the physical u'(R).
    ! Matching then gives chi=cos(delta)F+sin(delta)G with unit amplitude.
    subroutine pp_scattering(mesh, energy, u, delta, effective, residual, ok, &
        strength, &
        du, ddelta, deffective)
        type(radial_mesh_t), intent(in) :: mesh
        real(dp), intent(in) :: energy
        real(dp), intent(out) :: u(:), delta, effective, residual
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: strength
        real(dp), intent(out), optional :: du(:), ddelta, deffective
        real(dp), allocatable :: hb(:, :), ab(:, :), rhs(:, :), coeff(:, :), hv(:, :)
        real(dp) :: kinetic, k, eta, c2, vs, vss, vsd, vdd, scale, deriv, f, g, fp, &
            gp, aa, bb
        real(dp) :: dderiv, da, db, ds
        real(dp), allocatable :: dc(:, :), dvpot(:)
        integer, allocatable :: piv(:)
        integer :: n, kd, i, j, info
        ok = .false.; u = 0; delta = 0; effective = 0; residual = huge(1.0_dp)
        if (.not. ieee_is_finite(energy)) return
        if (energy <= 0) return
        n = mesh%n; kd = mesh%order - 1
        if (size(u) /= n) return
        if (abs(mesh%r(n) - mesh%rmax) > epsilon(mesh%rmax)*mesh%rmax) return
        kinetic = pp_hc**2/pp_mp; k = sqrt(energy/kinetic)
        eta = pp_alpha*pp_hc/(2*kinetic*k)
        call coulomb_fg0(eta, k*mesh%rmax, f, g, fp, gp, ok)
        if (.not. ok) return
        ok = .false.
        allocate (hb(kd + 1, n), ab(3*kd + 1, n - 1), rhs(n - 1, 1), piv(n - 1), &
            coeff(n, 1), hv(n, 1), dvpot(n))
        scale = 1
        if (present(strength)) scale = strength
        do i = 1, n
            call pp_reid(mesh%r(i), reid_kinetic, vs, vss, vsd, vdd)
            dvpot(i) = vs
        end do
        call radial_band_hamiltonian(mesh, 0, pp_alpha*pp_hc, &
            pp_mp/2, pp_hc, hb, potential=scale*dvpot)
        hb(1,:) = hb(1,:) - energy
        ab = 0; rhs = 0
        do j = 1, n - 1
            do i = max(1, j - kd), min(n - 1, j + kd)
                ab(2*kd + 1 + i - j, j) = hb(1 + abs(i - j), min(i, j))
            end do
            if (n - j <= kd) rhs(j, 1) = -hb(1 + n - j, j)*sqrt(mesh%w(n))
        end do
        call dgbsv(n - 1, kd, kd, 1, ab, 3*kd + 1, piv, rhs, n - 1, info)
        if (info /= 0) return
        coeff(:n - 1, 1) = rhs(:, 1); coeff(n, 1) = sqrt(mesh%w(n))
        call radial_band_matvec(hb, coeff, hv)
        residual = norm2(hv(:n - 1, 1))/norm2(coeff(:, 1))
        deriv = hv(n, 1)*sqrt(mesh%w(n))/kinetic
        call pp_match(f, g, fp, gp, deriv, k, aa, bb)
        call pp_match_norm(aa, bb, scale)
        if (scale <= tiny(scale)) return
        u = coeff(:, 1)/(sqrt(mesh%w)*scale)
        delta = atan2(bb, aa)
        c2 = coulomb_c0_squared(eta)
        if (abs(bb) <= tiny(bb)) then
            effective = huge(1.0_dp)
        else
            call pp_effective_range(aa, bb, c2, k, eta, coulomb_h(eta), effective)
        end if
        if (present(du)) then
            if (.not. present(ddelta)) return
            if (.not. present(deffective)) return
            allocate (dc(n, 1))
            rhs(:, 1) = -dvpot(:n - 1)*coeff(:n - 1, 1)
            ab = 0
            do j = 1, n - 1
                do i = max(1, j - kd), min(n - 1, j + kd)
                    ab(2*kd + 1 + i - j, j) = hb(1 + abs(i - j), min(i, j))
                end do
            end do
            call dgbsv(n - 1, kd, kd, 1, ab, 3*kd + 1, piv, rhs, n - 1, info)
            if (info /= 0) return
            dc(:n - 1, 1) = rhs(:, 1); dc(n, 1) = 0
            call radial_band_matvec(hb, dc, hv)
            dderiv = (hv(n, 1) + dvpot(n)*coeff(n, 1))*sqrt(mesh%w(n))/kinetic
            call pp_match_jvp(f, g, fp, gp, deriv, k, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, dderiv, 0.0_dp, da, db)
            call pp_match_norm_jvp(aa, bb, da, db, ds)
            du = dc(:, 1)/(sqrt(mesh%w)*scale) - u*ds/scale
            ddelta = (aa*db - bb*da)/scale**2
            call pp_effective_range_jvp(aa, bb, c2, k, eta, coulomb_h(eta), &
                da, db, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, deffective)
        end if
        ok = all(ieee_is_finite(u)) .and. ieee_is_finite(effective)
    end subroutine
    ! Parameter products for the explicit two nuclear-strength controls.
    ! Output order: B_d, P_D, a_pp, Lambda^2, S(0). Euclidean scalar pairing.
    subroutine pp_fusion_jvp(cfg, direction, result, product, ok)
        type(pp_fusion_config_t), intent(in) :: cfg
        real(dp), intent(in) :: direction(2)
        type(pp_fusion_result_t), intent(out) :: result
        real(dp), intent(out) :: product(5)
        logical, intent(out) :: ok
        real(dp) :: jac(5, 2)
        call pp_fusion_solve(cfg, result, ok, jac)
        product = 0
        if (ok) product = matmul(jac, direction)
    end subroutine

    subroutine pp_fusion_vjp(cfg, weights, result, product, ok)
        type(pp_fusion_config_t), intent(in) :: cfg
        real(dp), intent(in) :: weights(5)
        type(pp_fusion_result_t), intent(out) :: result
        real(dp), intent(out) :: product(2)
        logical, intent(out) :: ok
        real(dp) :: jac(5, 2)
        call pp_fusion_solve(cfg, result, ok, jac)
        product = 0
        if (ok) product = matmul(transpose(jac), weights)
    end subroutine
end module
