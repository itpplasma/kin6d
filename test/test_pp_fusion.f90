program test_pp_fusion
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use kin6d_pp_fusion
    use kin6d_radial_dvr, only: radial_locate, radial_band_hamiltonian
    use kin6d_pp_fusion_leaves, only: pp_reid, pp_reid_jvp, pp_reid_vjp
    use fortnum_special_coulomb, only: coulomb_fg0, coulomb_c0_squared
    implicit none
    type(pp_fusion_config_t) :: cfg, other
    type(pp_fusion_result_t) :: result, fine, coarse, plus, minus, tmp
    real(dp) :: ref(7), sample(4), vals(5), base(5), jac(5, 2), fd(5), &
        fd_coarse(5), prod(5), adj(2)
    real(dp) :: direction(2), weights(5), scale(5), error, step, f, g, fp, gp, k, eta
    real(dp) :: potential(4), dv(4), vp(4), vm(4), grad(2), val(40)
    real(dp), allocatable :: u(:), hb(:, :), hb0(:, :), zero(:)
    real(dp) :: phase, effective, resid, neutral_error
    integer :: unit, ios, i, j, idx(40)
    logical :: ok
    character(1024) :: path
    call get_command_argument(1, path)
    open (newunit=unit, file=trim(path), status='old')
    read (unit, *) ref
    call pp_fusion_solve(cfg, result, ok, jac)
    call require(ok, 'default solve')
    base = observables(result)
    print *, 'B,PD,a,lambda2,S,slope,f:', base, result%slope, result%phase_space
    print *, 'phase,residual:', result%phase_shift, result%eigen_residual, &
        result%scattering_residual, result%energy_extrapolation
    call require(abs(result%binding - 2.2246_dp) < 2e-5_dp, 'Reid binding energy')
    call require(abs(result%d_probability - .064696_dp) < 3e-6_dp, 'Reid D state')
    call require(abs(result%a_pp + 7.8_dp) < 0.1_dp, &
        'Coulomb-modified scattering length')
    call require(abs(result%lambda2 - 7.0_dp) < .15_dp, &
        'one-body overlap literature range')
    call require(abs(result%sfactor/4.01e-25_dp - 1) < .025_dp, &
        'Solar Fusion II few percent')
    call require(abs(result%slope - 11.2_dp) < .2_dp, 'Bahcall-May energy slope')
    call require(abs(result%phase_space - .144_dp) < .0006_dp, 'computed phase space')
    call require(result%eigen_residual < 2e-8_dp, 'bound residual')
    call require(result%scattering_residual < 2e-9_dp, 'scattering residual')
    call require(result%energy_extrapolation < 1e-5_dp, 'zero-energy extrapolation')
    call require(abs(sum(result%mesh%w*(result%us**2 + result%ud**2)) - 1) < 2e-12_dp, &
        'coupled normalization')
    call require(maxval(abs(base - ref(:5))/max(abs(ref(:5)), 1e-30_dp)) < 3e-7_dp, &
        'independent BVP/ODE scalar oracle')
    call require(abs(result%phase_shift/ref(6) - 1) < 3e-7_dp, &
        'independent phase shift')
    do
        read (unit, *, iostat=ios) sample
        if (ios < 0) exit
        call require(ios == 0, 'reference table read')
        call radial_locate(result%mesh, sample(1), idx(:cfg%order), val(:cfg%order))
        vals(:3) = [sum(val(:cfg%order)*sqrt(result%mesh%w(idx(:cfg%order)))* &
            result%us(idx(:cfg%order))), &
            sum(val(:cfg%order)*sqrt(result%mesh%w(idx(:cfg%order)))* &
            result%ud(idx(:cfg%order))), &
            sum(val(:cfg%order)*sqrt(result%mesh%w(idx(:cfg%order)))* &
            result%upp(idx(:cfg%order)))]
        call require(maxval(abs(vals(:3) - sample(2:))/max(abs(sample(2:)), 1e-20_dp)) &
            < 3e-4_dp, 'independent radial functions')
    end do
    close (unit)
    ! Genuine mesh/order and domain changes, compared at matched observables.
    other = cfg; other%order = 6; other%hmax = 1.2_dp
    call pp_fusion_solve(other, coarse, ok)
    call require(ok, 'coarse mesh')
    other = cfg; other%order = 10; other%hmax = .6_dp
    call pp_fusion_solve(other, fine, ok)
    call require(ok, 'fine mesh')
    scale = max(abs(base), [1.0_dp, .01_dp, 1.0_dp, 1.0_dp, 1e-25_dp])
    call require(maxval(abs(observables(fine) - base)/scale) < 2e-7_dp, &
        'mesh convergence')
    call require(maxval(abs(observables(coarse) - base)/scale) < 3e-5_dp, &
        'coarse convergence')
    other = cfg; other%rmax = 80
    call pp_fusion_solve(other, fine, ok)
    call require(ok, 'larger domain')
    print *, 'domain difference:', observables(fine) - base
    call require(maxval(abs(observables(fine) - base)/scale) < 2e-6_dp, &
        'tail convergence')
    ! Unit-amplitude Coulomb state with nuclear strength exactly zero.
    allocate (u(result%mesh%n), hb(cfg%order, result%mesh%n), &
        hb0(cfg%order, result%mesh%n), zero(result%mesh%n))
    call pp_scattering(result%mesh, cfg%energy, u, phase, effective, resid, ok, 0.0_dp)
    call require(ok, 'pure Coulomb solve')
    k = sqrt(cfg%energy*pp_mp/pp_hc**2)
    eta = pp_alpha*pp_hc/(2*(pp_hc**2/pp_mp)*k)
    neutral_error = 0
    do i = 1, result%mesh%n
        call coulomb_fg0(eta, k*result%mesh%r(i), f, g, fp, gp, ok)
        call require(ok, 'pure Coulomb reference')
        neutral_error = max(neutral_error, abs(u(i) - f)/sqrt(coulomb_c0_squared(eta)))
    end do
    call require(neutral_error < 2e-9_dp, 'F0 normalization without nuclear potential')
    zero = 0
    call radial_band_hamiltonian(result%mesh, 0, pp_alpha*pp_hc, pp_mp/2, pp_hc, hb)
    call radial_band_hamiltonian(result%mesh, 0, pp_alpha*pp_hc, pp_mp/2, pp_hc, &
        hb0, zero)
    call require(all(transfer(hb, [0_int64], size(hb)) == &
        transfer(hb0, [0_int64], size(hb0))), 'Coulomb assembly bit compatibility')
    ! Implicit eigen/linear-solve products against independently perturbed solves.
    step = 5e-4_dp
    do j = 1, 2
        other = cfg; other%strength(j) = 1 + step
        call pp_fusion_solve(other, plus, ok)
        call require(ok, 'positive strength perturbation')
        other%strength(j) = 1 - step
        call pp_fusion_solve(other, minus, ok)
        call require(ok, 'negative strength perturbation')
        fd_coarse = (observables(plus) - observables(minus))/(2*step)
        other = cfg; other%strength(j) = 1 + step/2
        call pp_fusion_solve(other, plus, ok)
        call require(ok, 'half-step positive perturbation')
        other%strength(j) = 1 - step/2
        call pp_fusion_solve(other, minus, ok)
        call require(ok, 'half-step negative perturbation')
        fd = (4*(observables(plus) - observables(minus))/step - fd_coarse)/3
        error = maxval(abs(jac(:, j) - fd)/max(abs(fd), scale))
        print *, 'implicit derivative axis,error:', j, error
        print *, 'jac, FD:', jac(:, j), fd
        call require(error < 3e-6_dp, 'implicit parameter response')
    end do
    direction = [.3_dp, -.7_dp]
    weights = [.2_dp, 2.0_dp, -.03_dp, .4_dp, 1e24_dp]
    call pp_fusion_jvp(cfg, direction, tmp, prod, ok)
    call require(ok, 'JVP')
    call pp_fusion_vjp(cfg, weights, tmp, adj, ok)
    call require(ok, 'VJP')
    call require(abs(dot_product(weights, prod) - dot_product(direction, adj)) &
        < 1e-11_dp, 'JVP/VJP pairing')
    ! Generated pointwise derivatives (radial and kinetic directions).
    step = 2e-5_dp
    call pp_reid(1.1_dp, reid_kinetic, potential(1), potential(2), potential(3), &
        potential(4))
    call pp_reid_jvp(1.1_dp, reid_kinetic, .3_dp, .7_dp, dv(1), dv(2), dv(3), dv(4))
    call pp_reid(1.1_dp + step*.3_dp, reid_kinetic + step*.7_dp, vp(1), vp(2), &
        vp(3), vp(4))
    call pp_reid(1.1_dp - step*.3_dp, reid_kinetic - step*.7_dp, vm(1), vm(2), &
        vm(3), vm(4))
    call require(maxval(abs(dv - (vp - vm)/(2*step))) < 2e-6_dp, &
        'generated potential JVP')
    call pp_reid_vjp(1.1_dp, reid_kinetic, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, grad(1), &
        grad(2))
    call require(abs(dot_product(dv, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]) - &
        dot_product(grad, [.3_dp, .7_dp])) < 2e-11_dp, 'generated potential transpose')
    other = cfg; other%energy = -1
    call pp_fusion_solve(other, tmp, ok)
    call require(.not. ok, 'invalid energy rejected')
    print *, 'pp fusion behavioral oracles passed'
contains
    function observables(r) result(v)
        type(pp_fusion_result_t), intent(in) :: r
        real(dp) :: v(5)
        v = [r%binding, r%d_probability, r%a_pp, r%lambda2, r%sfactor]
    end function
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine
end program
