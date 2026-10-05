program test_quantum_evolution
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t
    use kin6d_quantum_evolution, only: quantum_evolution_t, &
        initialize_quantum_evolution, solve_quantum_evolution, &
        sample_quantum_evolution, sample_quantum_evolution_box, &
        quantum_evolution_accepted, solve_quantum_evolution_tangent, &
        sample_quantum_evolution_tangent, quantum_evolution_vjp, &
        quantum_hamiltonian_action
    use kin6d_quantum_evolution_leaves, only: qe_rhs_term, qe_rhs_term_jvp, &
                                              qe_rhs_term_vjp
    implicit none
    type(quantum_evolution_t) :: state, plus, minus, coarse, invalid
    type(interval_t) :: ri(2), ii(2)
    real(dp) :: hr(2, 2), hi(2, 2), pr(2), pi(2), re(2), im(2), radius
    real(dp) :: dhr(2, 2), dhi(2, 2), dr(2), di(2), tr(2), ti(2), tradius
    real(dp) :: rp(2), ip(2), rm(2), imm(2), wr(2), wi(2), gr(2), gi(2)
    real(dp) :: ghr(2, 2), ghi(2, 2), grad_radius, left, right, error(2)
    real(dp) :: eps, h, ar, ai, x(5), dx(5), weights(2), jvp(2), vjp(5)
    complex(qp) :: exact(2), tangent_exact(2), initial(2), exact_p(2), exact_m(2)
    real(qp) :: t, norm_error
    integer :: cell, level
    logical :: ok, bounded, gb, hgb

    hr = reshape([.75_dp, .25_dp, .25_dp, -.5_dp], [2, 2])
    hi = reshape([0._dp, .375_dp, -.375_dp, 0._dp], [2, 2])
    pr = [-.5_dp, .375_dp]; pi = [.25_dp, -.125_dp]
    initial = cmplx(real(pr, qp), real(pi, qp), qp)
    h = .025_dp
    call initialize_quantum_evolution(state, hr, hi, pr, pi, 80, h, 1.25_dp, ok)
    call check(ok, 'arbitrary complex Hermitian input and signed amplitudes')
    call solve_quantum_evolution(state, ok)
    call check(ok, 'Cash-Karp candidate')
    call check(quantum_evolution_accepted(state), 'default 1% norm-scale target')
    do cell = 1, 80
        t = (real(cell - 1, qp) + .375_qp)*real(h, qp)
        call two_level(hr, hi, initial, 1.25_qp, t, exact)
        call sample_quantum_evolution(state, cell, .375_dp, re, im, radius, &
                                      ok, bounded)
        call check(ok .and. bounded, 'interior whole-time certificate available')
        norm_error = sqrt(sum(abs(cmplx(real(re, qp), real(im, qp), qp) - exact)**2))
        call check(norm_error <= real(radius, qp), 'independent interior Hilbert bound')
        call sample_quantum_evolution_box(state, cell, .375_dp, ri, ii, radius, &
                                          ok, bounded)
        call check(ok .and. bounded, 'strict component query')
        call check(all(real(ri%lo, qp) <= real(exact, qp)) .and. &
               all(real(exact, qp) <= real(ri%hi, qp)), 'true real components enclosed')
        call check(all(real(ii%lo, qp) <= aimag(exact)) .and. &
             all(aimag(exact) <= real(ii%hi, qp)), 'true imaginary components enclosed')
    end do
    do level = 1, 2
        call initialize_quantum_evolution(coarse, hr, hi, pr, pi, 8*2**level, &
                                          2._dp/real(8*2**level, dp), 1.25_dp, ok)
        call solve_quantum_evolution(coarse, ok)
        call sample_quantum_evolution(coarse, 8*2**level, 1._dp, re, im, &
                                      radius, ok, bounded)
        call check(ok .and. bounded, 'refinement certificate')
        call two_level(hr, hi, initial, 1.25_qp, 2._qp, exact)
        error(level) = real(sqrt(sum(abs(cmplx(real(re, qp), real(im, qp), qp) - &
                                         exact)**2)), dp)
    end do
    call check(error(1)/error(2) > 20._dp, 'fifth order candidate refinement')
    dhr = reshape([.125_dp, -.25_dp, -.25_dp, .5_dp], [2, 2])
    dhi = reshape([0._dp, -.125_dp, .125_dp, 0._dp], [2, 2])
    dr = [.125_dp, -.25_dp]; di = [-.125_dp, .375_dp]
    call solve_quantum_evolution_tangent(state, dhr, dhi, dr, di, ok)
    call check(ok, 'equation-level directional solve')
    call sample_quantum_evolution_tangent(state, 80, .625_dp, tr, ti, tradius, &
                                          ok, bounded)
    call check(ok .and. bounded, 'continuous tangent certificate')
    eps = 2._dp**(-18)
    t = (79._qp + .625_qp)*real(h, qp)
    call two_level_qp(cmplx(real(hr, qp), real(hi, qp), qp) + &
                      real(eps, qp)*cmplx(real(dhr, qp), real(dhi, qp), qp), &
                      initial + real(eps, qp)*cmplx(real(dr, qp), real(di, qp), qp), &
                      1.25_qp, t, exact_p)
    call two_level_qp(cmplx(real(hr, qp), real(hi, qp), qp) - &
                      real(eps, qp)*cmplx(real(dhr, qp), real(dhi, qp), qp), &
                      initial - real(eps, qp)*cmplx(real(dr, qp), real(di, qp), qp), &
                      1.25_qp, t, exact_m)
    tangent_exact = (exact_p - exact_m)/(2*real(eps, qp))
   norm_error = sqrt(sum(abs(cmplx(real(tr, qp), real(ti, qp), qp) - tangent_exact)**2))
   call check(norm_error < real(tradius, qp), 'independent quantum tangent containment')
    call initialize_quantum_evolution(plus, hr + eps*dhr, hi + eps*dhi, &
                                      pr + eps*dr, pi + eps*di, 80, h, 1.25_dp, ok)
    call initialize_quantum_evolution(minus, hr - eps*dhr, hi - eps*dhi, &
                                      pr - eps*dr, pi - eps*di, 80, h, 1.25_dp, ok)
    call solve_quantum_evolution(plus, ok, certify=.false.)
    call solve_quantum_evolution(minus, ok, certify=.false.)
    call sample_quantum_evolution(plus, 80, 1._dp, rp, ip, radius, ok, bounded)
    call check(ok .and. .not. bounded, 'candidate-only flag fails closed')
    call sample_quantum_evolution(minus, 80, 1._dp, rm, imm, radius, ok, bounded)
    call sample_quantum_evolution_tangent(state, 80, 1._dp, tr, ti, tradius, &
                                          ok, bounded)
    call check(maxval(abs((rp - rm)/(2*eps) - tr)) < 2.e-7_dp .and. &
 maxval(abs((ip - imm)/(2*eps) - ti)) < 2.e-7_dp, 'independent solve finite difference')
    wr = [.5_dp, -.25_dp]; wi = [.125_dp, .375_dp]
    call quantum_evolution_vjp(state, wr, wi, gr, gi, ghr, ghi, grad_radius, &
                               ok, gb, hgb)
    call check(ok .and. gb .and. .not. hgb, &
               'initial VJP bound; H gradient bound unavailable')
    call two_level_qp(-cmplx(real(hr, qp), real(hi, qp), qp), &
                      cmplx(real(wr, qp), real(wi, qp), qp), 1.25_qp, &
                      80*real(h, qp), exact)
    norm_error = sqrt(sum(abs(cmplx(real(gr, qp), real(gi, qp), qp) - exact)**2))
call check(norm_error <= real(grad_radius, qp), 'independent certified initial adjoint')
    left = dot_product(wr, tr) + dot_product(wi, ti)
    right = dot_product(gr, dr) + dot_product(gi, di) + sum(ghr*dhr) + sum(ghi*dhi)
    call check(abs(left - right) < 2.e-7_dp, &
               'real complex Hilbert terminal adjoint pairing')
    x = [.5_dp, -.25_dp, -.75_dp, .125_dp, 1.25_dp]
    dx = [.125_dp, -.375_dp, .5_dp, -.25_dp, .125_dp]
    weights = [.25_dp, -.5_dp]
    call qe_rhs_term_jvp(x(1), x(2), x(3), x(4), x(5), dx(1), dx(2), dx(3), &
                         dx(4), dx(5), jvp(1), jvp(2))
    call qe_rhs_term_vjp(x(1), x(2), x(3), x(4), x(5), weights(1), weights(2), &
                         vjp(1), vjp(2), vjp(3), vjp(4), vjp(5))
    call check(abs(dot_product(jvp, weights) - dot_product(dx, vjp)) < 1.e-14_dp, &
               'local RHS generated JVP/VJP pairing')
    call qe_rhs_term(x(1) + eps*dx(1), x(2) + eps*dx(2), x(3) + eps*dx(3), &
                     x(4) + eps*dx(4), x(5) + eps*dx(5), ar, ai)
    call qe_rhs_term(x(1) - eps*dx(1), x(2) - eps*dx(2), x(3) - eps*dx(3), &
                     x(4) - eps*dx(4), x(5) - eps*dx(5), left, right)
    call check(maxval(abs(([ar, ai] - [left, right])/(2*eps) - jvp)) < 2.e-10_dp, &
               'local RHS independent directional finite difference')
    call other_dimensions()
    call failures()
    if (command_argument_count() > 0) call benchmark()
    write (*, '(a,2es12.3)') 'Two-level candidate refinement errors: ', error
    write (*, '(a)') 'PASS: general Hermitian IVP, full-cell bounds, '// &
        'tangents, adjoint and failures'
contains
    subroutine two_level(hr, hi, psi, hbar, time, out)
        real(dp), intent(in) :: hr(2, 2), hi(2, 2)
        complex(qp), intent(in) :: psi(2)
        real(qp), intent(in) :: hbar, time
        complex(qp), intent(out) :: out(2)
        call two_level_qp(cmplx(real(hr, qp), real(hi, qp), qp), &
                          psi, hbar, time, out)
    end subroutine

    subroutine two_level_qp(matrix, psi, hbar, time, out)
        complex(qp), intent(in) :: matrix(2, 2), psi(2)
        real(qp), intent(in) :: hbar, time
        complex(qp), intent(out) :: out(2)
        complex(qp) :: traceless(2, 2), phase
        real(qp) :: mean, omega
        mean = real(matrix(1, 1) + matrix(2, 2), qp)/2
        traceless = matrix
        traceless(1, 1) = traceless(1, 1) - mean
        traceless(2, 2) = traceless(2, 2) - mean
        omega = sqrt(real(traceless(1, 1), qp)**2 + abs(matrix(1, 2))**2)
        phase = exp(cmplx(0._qp, -mean*time/hbar, qp))
        out = phase*(cos(omega*time/hbar)*psi - cmplx(0._qp, 1._qp, qp)* &
                     sin(omega*time/hbar)/omega*matmul(traceless, psi))
    end subroutine

    subroutine other_dimensions()
        type(quantum_evolution_t) :: packet
        real(dp) :: hr(5, 5), hi(5, 5), pr(5), pi(5), re(5), im(5), radius
        real(dp) :: out_re(5), out_im(5)
        complex(qp) :: exact(5), psi(5), action(5)
        integer :: j
        logical :: ok, bounded
        hr = 0._dp; hi = 0._dp
        do j = 1, 5
            hr(j, j) = real(j - 3, dp)/2
        end do
        pr = [-.5_dp, .125_dp, 0._dp, .75_dp, -.375_dp]
        pi = [.25_dp, -.75_dp, .5_dp, -.125_dp, .625_dp]
        psi = cmplx(real(pr, qp), real(pi, qp), qp)
        call initialize_quantum_evolution(packet, hr, hi, pr, pi, 100, &
                                          .03125_dp, 2._dp, ok)
        call solve_quantum_evolution(packet, ok)
        call sample_quantum_evolution(packet, 100, .5_dp, re, im, radius, ok, bounded)
        call check(ok .and. bounded, 'general five-dimensional Hermitian model')
        do j = 1, 5
            exact(j) = exp(cmplx(0._qp, -real(hr(j, j), qp)* &
                                 99.5_qp*.03125_qp/2, qp))*psi(j)
        end do
        call check(sqrt(sum(abs(cmplx(real(re, qp), real(im, qp), qp) - exact)**2)) &
                   <= real(radius, qp), 'signed diagonal unitary independent oracle')
        call quantum_hamiltonian_action(hr, hi, pr, pi, out_re, out_im, ok)
        action = matmul(cmplx(real(hr, qp), real(hi, qp), qp), psi)
        call check(ok .and. maxval(abs(cmplx(real(out_re, qp), real(out_im, qp), &
                     qp) - action)) < 1.e-14_qp, 'assembled general Hamiltonian action')
        hr = 0._dp
        call initialize_quantum_evolution(packet, hr, hi, pr, pi, 8, &
                                          .125_dp, .5_dp, ok)
        call solve_quantum_evolution(packet, ok)
        call sample_quantum_evolution(packet, 5, .375_dp, re, im, radius, ok, bounded)
        call check(ok .and. bounded, 'zero-generator free unitary certificate')
        call check(maxval(abs(re - pr)) + maxval(abs(im - pi)) <= radius, &
                   'arbitrary complex free-state reference')
    end subroutine

    subroutine failures()
        real(dp) :: bad(2, 2), badim(2, 2), nan, re(2), im(2), radius
        logical :: ok, bounded
        nan = ieee_value(0._dp, ieee_quiet_nan)
        bad = hr; bad(1, 2) = .5_dp
        call initialize_quantum_evolution(invalid, bad, hi, pr, pi, 2, .5_dp, &
                                          1._dp, ok)
        call check(.not. ok, 'non-Hermitian real input rejected')
        badim = hi; badim(1, 1) = .25_dp
        call initialize_quantum_evolution(invalid, hr, badim, pr, pi, 2, .5_dp, &
                                          1._dp, ok)
        call check(.not. ok, 'imaginary diagonal rejected')
        call initialize_quantum_evolution(invalid, hr, hi, pr, pi, 2, .5_dp, &
                                          0._dp, ok)
        call check(.not. ok, 'zero hbar rejected')
        call initialize_quantum_evolution(invalid, hr, hi, pr, pi, 2, nan, &
                                          1._dp, ok)
        call check(.not. ok, 'nonfinite time input rejected')
        call initialize_quantum_evolution(invalid, hr(:, :1), hi, pr, pi, 2, &
                                          .5_dp, 1._dp, ok)
        call check(.not. ok, 'matrix shape rejected before indexing')
        call solve_quantum_evolution(invalid, ok)
        call check(.not. ok, 'uninitialized solve rejected')
        call sample_quantum_evolution(state, 0, .5_dp, re, im, radius, ok, bounded)
        call check(.not. ok .and. .not. bounded, 'invalid time cell rejected')
        call sample_quantum_evolution(state, 1, nan, re, im, radius, ok, bounded)
        call check(.not. ok .and. .not. bounded, 'NaN fraction rejected')
        call sample_quantum_evolution(state, 1, .5_dp, re(:1), im, radius, &
                                      ok, bounded)
        call check(.not. ok .and. .not. bounded, 'sample output shape rejected')
        call initialize_quantum_evolution(invalid, hr, hi, pr, pi, 2, .5_dp, &
                                          1._dp, ok, budget=0._dp)
        call solve_quantum_evolution(invalid, ok)
        call check(ok .and. .not. quantum_evolution_accepted(invalid), &
                   'unmet accuracy budget fails acceptance')
        call sample_quantum_evolution(invalid, 2, 1._dp, re, im, radius, ok, bounded)
        call check(ok .and. bounded, 'valid wider enclosure retains honest status')
        call solve_quantum_evolution(invalid, ok, certify=.false.)
        call sample_quantum_evolution_box(invalid, 2, 1._dp, ri, ii, radius, &
                                          ok, bounded)
   call check(.not. ok .and. .not. bounded, 'candidate-only rebuild invalidates packet')
        bad = hr; bad(1, 1) = huge(0._dp)
        call initialize_quantum_evolution(invalid, bad, hi, pr, pi, 4, 10._dp, &
                                          1._dp, ok)
        if (ok) then
            call solve_quantum_evolution(invalid, ok)
            call check(.not. ok, 'overflow candidate fails closed')
        end if
    end subroutine

    subroutine benchmark()
        type(quantum_evolution_t) :: packet
        real(dp) :: t0, t1, candidate, complete, tangent_candidate, tangent_complete
        real(dp) :: radius, tangent_radius, re(2), im(2)
        integer :: trial, repeat
        integer, parameter :: repeats = 100, cells = 4
        logical :: ok, bounded
        call initialize_quantum_evolution(packet, hr, hi, pr, pi, cells, &
                                          .5_dp, 1.25_dp, ok)
        call check(ok, 'benchmark initialization')
        call solve_quantum_evolution(packet, ok)
        call solve_quantum_evolution_tangent(packet, dhr, dhi, dr, di, ok)
        call check(ok, 'warm primal and tangent setup')
        write (*, '(a)') 'trial,candidate_s,complete_s,extra_s,factor,'// &
            'primal_tangent_s,complete_tangent_s,extra_tangent_s,tangent_factor,'// &
            'state_radius,tangent_radius,accepted'
        do trial = 1, 20
            call cpu_time(t0)
            do repeat = 1, repeats
                call solve_quantum_evolution(packet, ok, certify=.false.)
            end do
            call cpu_time(t1)
            candidate = (t1 - t0)/repeats
            call check(ok, 'benchmark primal candidate')
            call cpu_time(t0)
            do repeat = 1, repeats
                call solve_quantum_evolution(packet, ok)
            end do
            call cpu_time(t1)
            complete = (t1 - t0)/repeats
            call check(ok, 'benchmark complete primal')
            call sample_quantum_evolution(packet, cells, 1._dp, re, im, &
                                          radius, ok, bounded)
            call check(ok .and. bounded, 'benchmark primal validation')
            call cpu_time(t0)
            do repeat = 1, repeats
                call solve_quantum_evolution(packet, ok, certify=.false.)
                call solve_quantum_evolution_tangent(packet, dhr, dhi, dr, di, ok, &
                                                     certify=.false.)
            end do
            call cpu_time(t1)
            tangent_candidate = (t1 - t0)/repeats
            call check(ok, 'benchmark primal plus tangent candidate')
            call cpu_time(t0)
            do repeat = 1, repeats
                call solve_quantum_evolution(packet, ok, certify=.false.)
                call solve_quantum_evolution_tangent(packet, dhr, dhi, dr, di, ok)
            end do
            call cpu_time(t1)
            tangent_complete = (t1 - t0)/repeats
            call check(ok, 'benchmark complete primal and tangent')
            call sample_quantum_evolution_tangent(packet, cells, 1._dp, re, im, &
                                                  tangent_radius, ok, bounded)
            call check(ok .and. bounded, 'benchmark tangent validation')
            call check(quantum_evolution_accepted(packet) .and. &
                       tangent_radius <= .01_dp, 'matched 1% state/tangent targets')
            write (*, '(i0,10(",",es24.16e3),",",l1)') trial, candidate, complete, &
                complete - candidate, complete/candidate, tangent_candidate, &
                tangent_complete, tangent_complete - tangent_candidate, &
                tangent_complete/tangent_candidate, radius, tangent_radius, &
                quantum_evolution_accepted(packet)
        end do
    end subroutine

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            write (*, '(a)') 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program test_quantum_evolution
