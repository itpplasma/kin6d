program test_quantum_reference_gallery
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, imul
    use kin6d_quantum_reference_gallery, only: quantum_observable_packet_t, &
        prepare_quantum_observables, sample_quantum_observables, &
        quantum_initial_norm, quantum_matrix_defect, quantum_physical_radius, &
        prepare_quantum_pair_sector
    use kin6d_quantum_reference_gallery_leaves, only: qgallery_fourier, &
        qgallery_fourier_jvp, qgallery_fourier_vjp
    implicit none
    type(quantum_observable_packet_t) :: packet
    integer :: sector_momentum(3, 2)
    type(interval_t) :: gr(4, 4), gi(4, 4), pr(6, 6), pi(6, 6), n0
    type(interval_t) :: hrb(2, 2), hib(2, 2)
    integer :: momentum(3, 4), spin(4), left(6), right(6), j
    real(dp) :: value(6), lower(6), upper(6)
    real(dp) :: ll(6), lh(6), hl(6), hh(6), x, norm, radius, defect
    real(dp) :: hr(2, 2), hi(2, 2), derivative, ga, gb, gn, gx, fd, p, m
    real(qp) :: truth(6), a, b, aq, bq, z
    logical :: ok
    momentum = 0; momentum(1, :) = [0, 0, 1, -1]
    spin = [1, -1, 1, -1]
    left = [1, 1, 1, 2, 2, 3]; right = [2, 3, 4, 3, 4, 4]
    a = real(0.6_dp, qp); b = real(0.8_dp, qp)
    call quantum_initial_norm([0.6_dp, 0.0_dp], [0.0_dp, 0.8_dp], n0, ok)
    call require(ok, 'initial raw-state norm')
    gr = interval(0.0_dp); gi = interval(0.0_dp)
    pr = interval(0.0_dp); pi = interval(0.0_dp)
    gr(1, 1) = imul(interval(0.6_dp), interval(0.6_dp)); gr(2, 2) = gr(1, 1)
    gr(3, 3) = imul(interval(0.8_dp), interval(0.8_dp)); gr(4, 4) = gr(3, 3)
    pr(1, 1) = gr(1, 1); pr(6, 6) = gr(3, 3)
    pi(1, 6) = imul(interval(-0.6_dp), interval(0.8_dp))
    pi(6, 1) = imul(interval(0.6_dp), interval(0.8_dp))
    do j = 0, 20
        x = -0.5_dp + real(j, dp)/20.0_dp
        call prepare_quantum_observables(packet, momentum, spin, left, right, &
            gr, gi, pr, pi, n0, n0, 0.0_dp, 0.0_dp, ok)
        call require(ok, 'point packet')
        call sample_quantum_observables(packet, x, value, lower, upper, ok)
        call require(ok, 'point sampling')
        call reference(real(x, qp), a, b, a*a + b*b, truth)
        call require(all(real(lower, qp) <= truth), 'complex-state lower containment')
        call require(all(real(upper, qp) >= truth), 'complex-state upper containment')
        call require(maxval(abs(value - real(truth, dp))) < 1e-12_dp, &
            'independent plane-wave pair amplitudes and Wick oracle')
    end do
    call prepare_quantum_observables(packet, momentum, spin, left, right, &
        gr, gi, pr, pi, n0, n0, 0.002_dp, 0.125_dp, ok)
    call require(ok, 'norm-to-observable and complete-cell packet')
    aq = a + real(0.001_dp, qp); bq = b - real(0.0005_dp, qp)
    call sample_quantum_observables(packet, 0.137_dp, value, ll, lh, ok)
    call require(ok, 'left node')
    call sample_quantum_observables(packet, 0.262_dp, value, hl, hh, ok)
    call require(ok, 'right node')
    call reference((real(0.137_dp, qp) + real(0.262_dp, qp))/2, &
        aq, bq, a*a + b*b, truth)
    call require(all((real(ll, qp) + real(hl, qp))/2 <= truth), &
        'connected lower chord contains perturbed true point')
    call require(all((real(lh, qp) + real(hh, qp))/2 >= truth), &
        'connected upper chord contains perturbed true point')
    call prepare_quantum_observables(packet, momentum, spin, left, right, &
        gr, gi, pr, pi, n0, n0, 0.0_dp, 0.0_dp, ok, [0.0625_dp, 1.0_dp, 1.0_dp])
    call require(ok, 'matched full-transverse slab packet')
    x = 0.173_dp
    call sample_quantum_observables(packet, x, value, lower, upper, ok)
    call reference(real(x, qp), a, b, a*a + b*b, truth)
    z = acos(-1.0_qp)/16
    truth(1) = 1 - 2*a*b*sin(2*acos(-1.0_qp)*real(x, qp))*sin(z)/z/(a*a + b*b)
    truth(2) = truth(1) - 1
    truth(4) = -2*a*a*b*b*(1 - cos(2*acos(-1.0_qp)*real(x, qp))*sin(z)/z) &
        /(a*a + b*b)**2
    call require(all(real(lower, qp) <= truth), 'analytic 3D slab lower containment')
    call require(all(real(upper, qp) >= truth), 'analytic 3D slab upper containment')
    sector_momentum = 0; sector_momentum(1, :) = [0, 1]
    call prepare_quantum_pair_sector(packet, sector_momentum, [0.6_dp, 0.0_dp], &
        [0.0_dp, 0.8_dp], n0, 0.0_dp, 0.0_dp, [0.0625_dp, 1.0_dp, 1.0_dp], ok)
    call require(ok, 'exact-sector observable packet')
    call sample_quantum_observables(packet, x, value, lower, upper, ok)
    call require(ok, 'exact-sector spatial sample')
    call require(all(real(lower, qp) <= truth), &
        'sector analytic slab lower containment')
    call require(all(real(upper, qp) >= truth), &
        'sector analytic slab upper containment')
    call require(maxval(abs(value - real(truth, dp))) < 1e-12_dp, &
        'sector and independent full-RDM formulas agree')
    momentum(2, 3) = 1
    call prepare_quantum_observables(packet, momentum, spin, left, right, &
        gr, gi, pr, pi, n0, n0, 0.0_dp, 0.0_dp, ok, [0.0625_dp, 1.0_dp, 1.0_dp])
    call sample_quantum_observables(packet, x, value, lower, upper, ok)
    call require(ok, '3D transverse mode')
    call require(lower(1) <= 1.0_dp .and. upper(1) >= 1.0_dp, &
        'full transverse average removes nonzero transverse modes')
    call prepare_quantum_observables(packet, momentum, spin, left, right, &
        gr, gi, pr, pi, n0, n0, ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, ok)
    call require(.not. ok, 'NaN radius rejected')
    call sample_quantum_observables(packet, x, value, lower, upper, ok)
    call require(.not. ok, 'failed preparation invalidates old packet')
    hr = 0.0_dp; hi = 0.0_dp; hrb = interval(0.0_dp); hib = interval(0.0_dp)
    hrb(1, 1) = interval(-0.001_dp, 0.001_dp)
    call quantum_matrix_defect(hr, hi, hrb, hib, defect, ok)
    call require(ok, 'Hamiltonian rounding-defect norm')
    call require(defect >= 0.001_dp .and. defect < 0.0010001_dp, 'Frobenius oracle')
    call quantum_physical_radius(n0, 0.001_dp, defect, 2.0_dp, 1.0_dp, radius, ok)
    call require(ok, 'physical coefficient propagation')
    norm = sqrt(real(a*a + b*b, dp))
    call require(radius >= 0.001_dp + 2.0_dp*norm*0.001_dp, &
        'unitary perturbation bound')
    call qgallery_fourier_jvp(0.4_dp, -0.2_dp, 2.0_dp, 0.13_dp, &
        0.1_dp, -0.3_dp, 0.0_dp, 0.07_dp, derivative)
    call qgallery_fourier(0.4_dp + 1e-6_dp*0.1_dp, -0.2_dp - 1e-6_dp*0.3_dp, &
        2.0_dp, 0.13_dp + 1e-6_dp*0.07_dp, p)
    call qgallery_fourier(0.4_dp - 1e-6_dp*0.1_dp, -0.2_dp + 1e-6_dp*0.3_dp, &
        2.0_dp, 0.13_dp - 1e-6_dp*0.07_dp, m)
    fd = (p - m)/2e-6_dp
    call require(abs(fd - derivative) < 1e-9_dp, 'independent Fourier JVP difference')
    call qgallery_fourier_vjp(0.4_dp, -0.2_dp, 2.0_dp, 0.13_dp, 0.6_dp, &
        ga, gb, gn, gx)
    call require(abs(0.6_dp*derivative - (0.1_dp*ga - 0.3_dp*gb + 0.07_dp*gx)) &
        < 1e-13_dp, 'Fourier JVP/VJP real pairing')
    call pair_basis_guards()
    print *, 'PASS: quantum spatial observables, slab bands and Hilbert bridge'
contains
    subroutine pair_basis_guards()
        integer :: momenta(3, 3), spins(3), first(3), second(3)
        type(interval_t) :: gamma(3, 3), gamma_i(3, 3)
        type(interval_t) :: pairs(3, 3), pairs_i(3, 3)
        real(dp) :: result(6), lo(6), hi(6)
        logical :: valid
        momenta = 0; momenta(1, 3) = 1; spins = [1, -1, 1]
        gamma = interval(0.0_dp); gamma_i = interval(0.0_dp)
        pairs = interval(0.0_dp); pairs_i = interval(0.0_dp)
        gamma(1, 1) = interval(1.0_dp); gamma(2, 2) = interval(1.0_dp)
        pairs(1, 1) = interval(1.0_dp)
        first = [1, 1, 2]; second = [2, 3, 3]
        call prepare_quantum_observables(packet, momenta, spins, first, second, &
            gamma, gamma_i, pairs, pairs_i, interval(1.0_dp), interval(1.0_dp), &
            0.0_dp, 0.0_dp, valid)
        call require(valid, 'canonical three-orbital Slater pair basis')
        call sample_quantum_observables(packet, 0.0_dp, result, lo, hi, valid)
        call require(valid, 'Slater pair observable')
        call require(lo(2) <= 0.0_dp .and. hi(2) >= 0.0_dp, &
            'independent Slater Wick cumulant vanishes')
        call require(abs(result(2)) < 1e-12_dp, 'Slater candidate cumulant')
        first = [1, 1, 2]; second = [2, 2, 3]
        call prepare_quantum_observables(packet, momenta, spins, first, second, &
            gamma, gamma_i, pairs, pairs_i, interval(1.0_dp), interval(1.0_dp), &
            0.0_dp, 0.0_dp, valid)
        call require(.not. valid, 'repeated ordered pair is not a pair basis')
        call sample_quantum_observables(packet, 0.0_dp, result, lo, hi, valid)
        call require(.not. valid, 'failed preparation invalidates old packet')
        first = [1, 1, 2]; second = [2, 3, 3]; momenta(1, 3) = 0
        call prepare_quantum_observables(packet, momenta, spins, first, second, &
            gamma, gamma_i, pairs, pairs_i, interval(1.0_dp), interval(1.0_dp), &
            0.0_dp, 0.0_dp, valid)
        call require(.not. valid, 'repeated spin orbital is not orthonormal')
    end subroutine
    subroutine reference(x, a, b, norm0, result)
        real(qp), intent(in) :: x, a, b, norm0
        real(qp), intent(out) :: result(6)
        real(qp) :: phase, n
        phase = 2*acos(-1.0_qp)*x; n = (a*a + b*b)/norm0
        result(1) = (a*a + b*b - 2*a*b*sin(phase))/norm0
        result(2) = result(1) - n*n
        result(3) = 0
        result(4) = -2*a*a*b*b*(1 - cos(phase))/norm0**2
        result(5:6) = n
    end subroutine
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program
