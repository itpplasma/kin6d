program test_pp_fusion_eft
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use kin6d_pp_fusion_eft
    use kin6d_pp_fusion, only: pp_mp, pp_alpha, pp_fusion_config_t, &
        pp_fusion_result_t, pp_fusion_solve
    use kin6d_pp_fusion_eft_leaves, only: eft_orders, eft_coulomb_slope
    implicit none
    type(pp_eft_config_t) :: cfg, other
    type(pp_eft_result_t) :: answer, plus, minus, tmp
    type(pp_fusion_config_t) :: reid_cfg
    type(pp_fusion_result_t) :: reid
    real(dp), parameter :: chi_ref(5) = [.01_dp, .05_dp, .15_dp, .3_dp, 1._dp]
    real(dp), parameter :: j_ref(5) = [ &
        .95921488556543574153389580285689_dp, &
        .87027848251196933391839669931476_dp, &
        .74478074260397175047588530308308_dp, &
        .63323931847582433304191017240908_dp, &
        .40365263767680592565892150063072_dp]
    real(dp) :: jc, indicator, slope, gamma, chi, cc(3), res, lb, ref(16), vals(16)
    real(dp) :: scaled(3), rs, lbs, unit_scale, jac(8, 5), v(5), weights(8)
    real(dp) :: prod(8), adj(5), fd(8), h, scale(8), mu, screen, delta(5)
    integer :: i, j, unit
    logical :: ok
    character(1024) :: path
    call get_command_argument(1, path)
    ! Values of E1 independently evaluated at 70 decimal digits, not a shared leaf.
    do i = 1, size(chi_ref)
        call pp_eft_coulomb_integral(chi_ref(i), jc, indicator, ok)
        print *, 'Coulomb chi,J,indicator,ok:', chi_ref(i), jc, indicator, ok
        call require(ok, 'Coulomb integral status')
        call require(abs(jc - j_ref(i)) < 3e-14_dp, 'independent E1 oracle')
        call eft_coulomb_slope(chi_ref(i), jc, slope)
        h = 1e-5_dp
        if (i == 1 .or. i == 5) cycle ! perturb both sides only within domain
        call pp_eft_coulomb_integral(chi_ref(i) + h, res, indicator, ok)
        call require(ok, 'positive integral perturbation')
        call pp_eft_coulomb_integral(chi_ref(i) - h, lb, indicator, ok)
        call require(ok, 'negative integral perturbation')
        call require(abs((res - lb)/(2*h) - slope) < 2e-8_dp, 'integral derivative')
    end do
    call pp_eft_coulomb_integral(0._dp, jc, indicator, ok)
    call require(ok .and. abs(jc - 1) < 1e-15_dp, 'neutral integral limit')
    ! No Coulomb and zero effective ranges: contact/overlap identity 1-gamma*a.
    call eft_orders(.23_dp, -7.8_dp, 0._dp, 0._dp, 0._dp, 1.26_dp, &
        .7_dp, 0._dp, jc, 0._dp, cc(1), cc(2), cc(3), res, lb)
    call require(abs(cc(1) - (1 + .23_dp*7.8_dp)) < 1e-14_dp, &
        'neutral analytic overlap')
    call require(maxval(abs(cc(2:))) < 1e-14_dp, 'neutral zero-range orders')
    call require(abs(res - cc(1)) < 1e-14_dp, 'zero-range resummation')
    ! Published KR (2001) (33),(55): their stated low-energy parameters.
    gamma = 45.7_dp/197.327_dp; chi = 1._dp/137*938.9_dp/(197.327_dp*gamma)
    call pp_eft_coulomb_integral(chi, jc, indicator, ok)
    call require(ok, 'KR Coulomb integral')
    call eft_orders(gamma, -7.82_dp, 2.79_dp, 1.76_dp, 0._dp, 1.26_dp, &
        1e8_dp, chi, jc, 0._dp, cc(1), cc(2), cc(3), res, lb)
    call require(abs(cc(1) - 2.51_dp) < .005_dp, 'KR LO amplitude 2.51')
    call require(abs(cc(1)**2 - 6.30_dp) < .02_dp, 'KR LO squared 6.30')
    ! KR's rounded prose 6.45 does not equal (55) with their stated inputs.
    ! Test the equation against a separate high-precision evaluation instead.
    call require(abs(sum(cc(:2))**2 - 6.478515339676822_dp) < 3e-9_dp, &
        'KR (55) independent evaluation, rho_d=1.76 fm')
    ! BC (13): c0*(1+0.039+0.029+...) for L1A=6, original inputs.
    mu = 140._dp/197.327_dp
    gamma = 45.69_dp/197.327_dp
    chi = 938.9_dp/(137*197.327_dp*gamma)
    screen = 938.9_dp/(137*197.327_dp)*(log(mu*sqrt(acos(-1._dp))* &
        137*197.327_dp/938.9_dp) + 1 - 1.5_dp*.577_dp)
    call pp_eft_coulomb_integral(chi, jc, indicator, ok)
    call require(ok, 'BC Coulomb integral')
    call eft_orders(gamma, -7.82_dp, 2.79_dp, 1.764_dp, 6._dp, 1.26_dp, &
        mu, chi, jc, screen, cc(1), cc(2), cc(3), res, lb)
    call require(abs(sum(cc(:2)) - 2.51_dp*(1 + .039_dp)) < .005_dp, &
        'BC published cumulative NLO amplitude at quoted input precision')
    call require(abs(cc(2)/cc(1) - .0378990249535741_dp) < 3e-14_dp, &
        'BC (7),(9),(10) independent NLO coefficient')
    call require(abs(cc(3)/cc(1) - .029_dp) < .001_dp, 'BC published NNLO coefficient')
    ! Ando (29): LO at mu->infinity and L1A=0 is the effective-range result.
    gamma = 45.70_dp/197.327_dp
    chi = pp_alpha*pp_mp/(197.327_dp*gamma)
    call pp_eft_coulomb_integral(chi, jc, indicator, ok)
    call require(ok, 'Ando integral')
    call eft_orders(gamma, -7.8063_dp, 2.794_dp, 1.760_dp, 0._dp, 1.2695_dp, &
        1e8_dp, chi, jc, 0._dp, cc(1), cc(2), cc(3), res, lb)
    call require(abs(res - 2.641_dp) < .002_dp, 'Ando resummed LO 2.641')
    ! Simultaneous change of length unit: gamma,mu,screen inverse length;
    ! a,r,rho lengths; L1A volume. All amplitudes must remain dimensionless.
    unit_scale = 10
    call eft_orders(gamma/unit_scale, -7.8063_dp*unit_scale, &
        2.794_dp*unit_scale, 1.760_dp*unit_scale, 0._dp, 1.2695_dp, &
        1e8_dp/unit_scale, chi, jc, 0._dp, scaled(1), scaled(2), scaled(3), rs, lbs)
    call require(maxval(abs(scaled - cc)) < 1e-13_dp, 'length dimensions')
    call require(abs(rs - res) < 1e-13_dp, 'resummed dimensions')
    open (newunit=unit, file=trim(path), status='old')
    do i = 1, 2
        cfg%l1a = 0
        if (i == 2) cfg%l1a = 3.9_dp
        call pp_eft_evaluate(cfg, answer, ok, jac)
        call require(ok, 'physical EFT evaluation')
        read (unit, *) ref
        vals(:4) = answer%lambda; vals(5:8) = answer%lambda2
        vals(9:12) = answer%sfactor; vals(13:) = answer%slope
        ! Phase integral's endpoint singularity limits the inherited 128-node rule.
        call require(maxval(abs(vals(:8) - ref(:8))) < 3e-13_dp, &
            'independent high-precision contact/E1 amplitude and rate')
        call require(maxval(abs(vals(9:12)/ref(9:12) - 1)) < 2e-7_dp, &
            'independent adaptive phase space and MeV barn normalization')
        call require(maxval(abs(vals(13:) - ref(13:))) < 2e-5_dp, &
            'independent finite-energy spectral integral and phase slope')
    end do
    close (unit)
    print *, 'Lambda2:', answer%lambda2
    print *, 'S [MeV b]:', answer%sfactor
    print *, 'S prime/S [MeV^-1]:', answer%slope
    call require(abs(answer%slope(3) - 11.3_dp) < .1_dp, 'Chen-Liu-Yu published slope')
    print *, 'truncation estimates:', answer%truncation_s
    ! Directional finite differences in every admitted low-energy input.
    scale(:4) = 1; scale(5:) = 1e-25_dp
    do j = 1, 5
        v = 0; v(j) = 1; h = 2e-5_dp
        other = cfg
        delta = h*v
        call shift(other, delta)
        call pp_eft_evaluate(other, plus, ok)
        call require(ok, 'positive input perturbation')
        other = cfg
        delta = -h*v
        call shift(other, delta)
        call pp_eft_evaluate(other, minus, ok)
        call require(ok, 'negative input perturbation')
        fd(:4) = (plus%lambda2 - minus%lambda2)/(2*h)
        fd(5:) = (plus%sfactor - minus%sfactor)/(2*h)
        call require(maxval(abs(fd - jac(:, j))/scale) < 2e-7_dp, &
            'full generated input derivative including Coulomb and phase space')
    end do
    v = [.03_dp, -.2_dp, .4_dp, -.1_dp, .7_dp]
    weights = [.3_dp, -.2_dp, .1_dp, .4_dp, 1e24_dp, -2e24_dp, 3e24_dp, 4e24_dp]
    call pp_eft_jvp(cfg, v, tmp, prod, ok)
    call require(ok, 'EFT JVP')
    call pp_eft_vjp(cfg, weights, tmp, adj, ok)
    call require(ok, 'EFT VJP')
    call require(abs(dot_product(weights, prod) - dot_product(v, adj)) < 1e-12_dp, &
        'physical-input JVP/VJP pairing')
    call pp_fusion_solve(reid_cfg, reid, ok)
    call require(ok, 'Reid68 comparison solve')
    call require(abs(answer%sfactor(3) - reid%sfactor) < answer%truncation_s(3), &
        'Reid68 lies within the NNLO EFT truncation estimate')
    call require(abs(answer%sfactor(3) - 4.01e-25_dp) < answer%truncation_s(3), &
        'Adelberger lies within NNLO EFT truncation estimate')
    other = cfg; other%rho_d = 5
    call pp_eft_evaluate(other, tmp, ok)
    call require(.not. ok, 'deuteron pole with negative residue rejected')
    other = cfg; other%binding = -1
    call pp_eft_evaluate(other, tmp, ok)
    call require(.not. ok, 'negative binding rejected')
    other = cfg; other%r_pp = ieee_value(0._dp, ieee_quiet_nan)
    call pp_eft_evaluate(other, tmp, ok)
    call require(.not. ok, 'nonfinite input rejected')
    call pp_eft_coulomb_integral(-.1_dp, jc, indicator, ok)
    call require(.not. ok, 'attractive Coulomb outside model rejected')
    print *, 'pionless EFT independent oracles passed'
contains
    subroutine shift(c, delta)
        type(pp_eft_config_t), intent(inout) :: c
        real(dp), intent(in) :: delta(5)
        c%binding = c%binding + delta(1); c%a_pp = c%a_pp + delta(2)
        c%r_pp = c%r_pp + delta(3); c%rho_d = c%rho_d + delta(4)
        c%l1a = c%l1a + delta(5)
    end subroutine
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine
end program
