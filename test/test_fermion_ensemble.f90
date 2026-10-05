program test_fermion_ensemble
    ! Q3 finite-N ideal fermions. Oracles written here, independent of the
    ! production recursion: brute-force enumeration of all N-orbital Slater
    ! determinants (binary128), the zero-temperature thermodynamic-limit
    ! exchange hole, the certified homogeneous Fermi gas, Boltzmann limits,
    ! ground-state asymptotics and finite differences of independent solves.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, interval, operator(+), operator(*)
    use kin6d_fermion_ensemble
    use kin6d_fermi_gas, only: fermi_workspace_t, fermi_state_t, &
        initialize_fermi_workspace, solve_fermi_gas, occupation_fermi_gas
    use kin6d_fermion_ensemble_leaves, only: fermion_box_eos, fermion_box_eos_jvp
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), scan(3) = [0.05_dp, 1.0_dp, 10.0_dp]
    type(fermion_spectrum_t) :: spec, wide
    type(canonical_state_t) :: st, other, plus, minus
    type(fermi_workspace_t) :: work
    type(fermi_state_t) :: fermi
    real(dp) :: radii(6), g(6), trunc, ef, b, values(5), vp(5), vm(5), h, fd
    real(dp) :: rho(1), dev(3), occupation_error, lower, upper, f, tangent(5)
    type(interval_t) :: bounds(5)
    real(qp) :: log_z, cancellation
    logical :: ok, cert
    integer :: n, k, a

    radii = [0.0_dp, 0.05_dp, 0.13_dp, 0.21_dp, 0.34_dp, 0.47_dp]
    ! A. Brute force on the closed 19-orbital box set |n|**2<=2.
    call box_spectrum(2, spec, ok, closed=.true.)
    call require(ok .and. spec%orbitals == 19 .and. spec%shells == 3, 'box spectrum')
    do n = 2, 4
        call brute_box(n, 0.4_dp)
        call brute_box(n, 2.0_dp)
    end do
    ! B. Brute force on the closed 10-orbital trap set k<=2.
    call trap_spectrum(2, spec, ok, closed=.true.)
    call require(ok .and. spec%orbitals == 10, 'trap spectrum')
    call brute_trap(3, 0.7_dp)
    ! C. Open trap (rigorous remainder) versus a far larger closed set.
    call trap_spectrum(60, spec, ok)
    call solve_canonical(spec, 20, 0.5_dp, st, ok)
    call require(ok .and. st%certified, 'open trap certificate')
    call trap_spectrum(200, wide, ok, closed=.true.)
    call solve_canonical(wide, 20, 0.5_dp, other, ok)
    call require(ok .and. inside(st%energy_bound, real(other%energy, qp)) .and. &
                 inside(st%log_z_bound, real(other%log_z, qp)), 'trap remainder')
    ! D. Alternating power-sum recursion: agreement when well conditioned,
    !    catastrophic cancellation cold; cold DP against ground asymptotics.
    call box_spectrum(12, spec, ok, closed=.true.)
    call solve_canonical(spec, 7, 0.4_dp, st, ok)
    call power_sum_recursion(spec, 7, 0.4_dp, log_z, cancellation, ok)
    call require(ok .and. abs(real(log_z, dp) - st%log_z) < 1.0e-12_dp* &
                 abs(st%log_z) + 1.0e-12_dp, 'power-sum recursion agreement')
    call solve_canonical(spec, 19, 20.0_dp, st, ok)
    call require(ok .and. st%certified, 'cold certificate')
    ! Closed shells 1+6+12, E0=6+24; lowest excitations 12*8 at gap one.
    call require(abs(st%log_z + 20.0_dp*30.0_dp - 96.0_dp*exp(-20.0_dp)) < &
                 1.0e-12_dp, 'ground-state asymptotics of lnZ')
    call power_sum_recursion(spec, 19, 20.0_dp, log_z, cancellation, ok)
    call require(cancellation > 1.0e30_qp, 'cold power-sum cancellation exposed')
    ! E. Pauli hole, sum rules and opposite-spin independence at theta=1.
    call box_parameters(66, 1.0_dp, ef, b)
    call box_spectrum(box_cutoff(b, 33), spec, ok)
    call solve_canonical(spec, 33, b, st, ok)
    call require(ok .and. st%certified, 'box certificate theta=1')
    call require(inside(sum_bound(spec, st), 33.0_qp), 'particle number')
    call same_spin_pair(spec, st, radii(1:1), g(1:1), ok, 1.0e-12_dp, trunc)
    call require(ok .and. abs(g(1)) < 1.0e-11_dp + trunc, 'Pauli g(0)=0 sum rule')
    call require(abs(opposite_spin_pair(spec, st, st) - 1.0_dp) < 1.0e-12_dp, &
                 'opposite-spin g=1')
    call one_body_density(spec, st, radii(1:1), rho, ok)
    call require(ok .and. abs(rho(1) - 1.0_dp) < 1.0e-12_dp, 'rho(0)=n')
    ! F. Thermodynamic limits (convergence indicators, not equalities).
    call box_spectrum(60, wide, ok)
    call cold_hole(33, dev(1)); call cold_hole(251, dev(2)); call cold_hole(1021, dev(3))
    call require(dev(3) < dev(2) .and. dev(2) < dev(1) .and. dev(3) < 2.0e-4_dp, &
                 'cold exchange hole converges to 1-9(j1(x)/x)**2')
    call initialize_fermi_workspace(work, 16384, ok)
    call solve_fermi_gas(work, 1.0_dp, fermi, ok)
    call require(ok, 'thermodynamic-limit reference')
    call box_eos(st, st, values, bounds, ok, cert)
    call require(ok .and. cert, 'box EOS certificate')
    dev(1) = abs(values(1) - fermi%energy)
    call box_parameters(246, 1.0_dp, ef, b)
    call box_spectrum(box_cutoff(b, 123), spec, ok)
    call solve_canonical(spec, 123, b, other, ok)
    call box_eos(other, other, vp, bounds, ok, cert)
    dev(2) = abs(vp(1) - fermi%energy)
    call require(dev(2) < dev(1) .and. dev(2) < 2.0e-3_dp*fermi%energy, &
                 'finite-N energy approaches grand-canonical limit')
    occupation_error = 0.0_dp
    do a = 1, spec%shells
        call occupation_fermi_gas(fermi, sqrt(spec%level(a)/ef), f, lower, upper, ok)
        occupation_error = max(occupation_error, abs(other%occupation(a) - f))
    end do
    call require(occupation_error < 0.02_dp, 'canonical versus Fermi-Dirac occupations')
    print '(a,3es11.3)', 'theta=1 energy deviation N=33,123; occupation:', dev(1:2), &
        occupation_error
    ! G. Classical limits: fermions at theta=10 and the Boltzmann preparation.
    call box_parameters(66, 10.0_dp, ef, b)
    call box_spectrum(box_cutoff(b, 33), spec, ok)
    call solve_canonical(spec, 33, b, st, ok)
    call box_eos(st, st, values, bounds, ok, cert)
    call require(ok .and. cert .and. abs(values(1)/15.0_dp - 1.0_dp) < 0.01_dp, &
                 'fermion classical limit e=3theta/2')
    call solve_boltzmann(spec, 33, b, other, ok)
    call box_eos(other, other, vp, bounds, ok, cert)
    call require(ok .and. abs(vp(1)/15.0_dp - 1.0_dp) < 1.0e-9_dp .and. &
                 abs(vp(2)/10.0_dp - 1.0_dp) < 1.0e-9_dp, 'Boltzmann box e, p')
    call require(values(1) > vp(1), 'Fermi pressure exceeds classical')
    ! H. Derivatives: de/dtheta, P=2E/(3V) via lnZ and the bulk modulus.
    call box_parameters(66, 0.5_dp, ef, b)
    call box_spectrum(box_cutoff(b, 33), spec, ok)
    call solve_canonical(spec, 33, b, st, ok)
    call box_eos(st, st, values, bounds, ok, cert)
    h = 1.0e-4_dp
    call solve_canonical(spec, 33, 1.0_dp/((0.5_dp + h)*ef), plus, ok)
    call solve_canonical(spec, 33, 1.0_dp/((0.5_dp - h)*ef), minus, ok)
    call box_eos(plus, plus, vp, bounds, ok, cert)
    call box_eos(minus, minus, vm, bounds, ok, cert)
    fd = (vp(1) - vm(1))/(2.0_dp*h)
    call require(abs(fd - values(3)) < 1.0e-6_dp*abs(fd), 'de/dtheta finite difference')
    call solve_canonical(spec, 33, b*(1.0_dp + h), plus, ok)
    call solve_canonical(spec, 33, b*(1.0_dp - h), minus, ok)
    fd = -(2.0_dp*plus%log_z - 2.0_dp*minus%log_z)/(2.0_dp*h*b)/(3.0_dp*33.0_dp*ef)
    call require(abs(fd - values(2)) < 1.0e-6_dp*values(2), 'P=2E/(3V) from dlnZ/db')
    ! Box scaled by lambda at fixed kB*T: b->b/lambda**2, P->lambda**-5*p.
    call solve_canonical(spec, 33, b/(1.0_dp + h)**2, plus, ok)
    call solve_canonical(spec, 33, b/(1.0_dp - h)**2, minus, ok)
    call box_eos(plus, plus, vp, bounds, ok, cert)
    call box_eos(minus, minus, vm, bounds, ok, cert)
    fd = -((1.0_dp + h)**(-5)*vp(2) - (1.0_dp - h)**(-5)*vm(2))/(2.0_dp*h)/3.0_dp
    call require(abs(fd - values(5)) < 1.0e-6_dp*values(5), 'bulk modulus n dP/dn')
    call fermion_box_eos_jvp(st%energy*2, st%variance*2, b, 66.0_dp, ef, 1.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, tangent(1), tangent(2), tangent(3), &
        tangent(4), tangent(5))
    call require(abs(tangent(1)*st%energy*2 - values(1)) < 1.0e-14_dp*values(1), &
                 'generated EOS tangent')
    print '(a,5f10.6)', 'theta=0.5 N=33 [e,p,de/dtheta,dp/dtheta,B]:', values
    ! I. Certified widths across the scan.
    do k = 1, 3
        call box_parameters(66, scan(k), ef, b)
        call box_spectrum(box_cutoff(b, 33), spec, ok)
        call solve_canonical(spec, 33, b, st, ok)
        call box_eos(st, st, values, bounds, ok, cert)
        print '(a,f6.2,l2,5es10.2)', 'theta, certified, EOS widths:', scan(k), cert, bounds%hi - bounds%lo
        call require(ok .and. cert .and. maxval(bounds(1:2)%hi - bounds(1:2)%lo) < &
            1.0e-8_dp*values(1) .and. maxval(bounds(3:5)%hi - bounds(3:5)%lo) < &
            1.0e-4_dp*maxval(abs(values(3:5))), 'scan certificate widths')
    end do
    print '(a)', 'finite-N ideal fermion oracles passed'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    logical function inside(bound, value)
        type(interval_t), intent(in) :: bound
        real(qp), intent(in) :: value
        inside = real(bound%lo, qp) <= value .and. value <= real(bound%hi, qp)
    end function inside

    function sum_bound(spectrum, state) result(total)
        type(fermion_spectrum_t), intent(in) :: spectrum
        type(canonical_state_t), intent(in) :: state
        type(interval_t) :: total
        integer :: j
        total = state%omitted_bound
        do j = 1, spectrum%shells
            total = total + interval(real(spectrum%degeneracy(j), dp))* &
                    state%occupation_bound(j)
        end do
    end function sum_bound

    real(qp) function sinc(z)
        real(qp), intent(in) :: z
        sinc = 1.0_qp
        if (abs(z) > 0.0_qp) sinc = sin(z)/z
    end function sinc

    !> All C(19,n) determinants of the |n|**2<=2 box orbitals.
    subroutine brute_box(n, beta)
        integer, intent(in) :: n
        real(dp), intent(in) :: beta
        integer :: vec(3, 19), level(19), m, i, j, l, mask, members(19), count
        real(qp) :: z, e, e2, w, occ(19), pair(6), energy, mean, shell_occ
        real(dp) :: production(6)
        m = 0
        do i = -1, 1
            do j = -1, 1
                do l = -1, 1
                    if (i*i + j*j + l*l > 2) cycle
                    m = m + 1
                    vec(:, m) = [i, j, l]
                    level(m) = i*i + j*j + l*l
                end do
            end do
        end do
        z = 0; e = 0; e2 = 0; occ = 0; pair = 0
        do mask = 0, 2**19 - 1
            if (popcnt(mask) /= n) cycle
            count = 0
            do i = 1, 19
                if (btest(mask, i - 1)) then
                    count = count + 1
                    members(count) = i
                end if
            end do
            energy = real(sum(level(members(1:n))), qp)
            w = exp(-real(beta, qp)*energy)
            z = z + w; e = e + w*energy; e2 = e2 + w*energy**2
            occ(members(1:n)) = occ(members(1:n)) + w
            do i = 1, n
                do j = 1, n
                    if (i == j) cycle
                    do l = 1, 6
                        pair(l) = pair(l) + w*(1.0_qp - sinc(2.0_qp*acos(-1.0_qp)* &
                            sqrt(real(sum((vec(:, members(i)) - vec(:, members(j)))**2), &
                            qp))*real(radii(l), qp)))
                    end do
                end do
            end do
        end do
        mean = e/z
        call solve_canonical(spec, n, beta, st, ok)
        call require(ok .and. st%certified, 'brute-force box certificate')
        call require(abs(st%log_z - real(log(z), dp)) < 1.0e-12_dp*max(1.0_dp, &
            abs(st%log_z)), 'brute-force lnZ')
        call require(abs(st%energy - real(mean, dp)) < 1.0e-12_dp*mean, 'brute-force E')
        call require(abs(st%variance - real(e2/z - mean**2, dp)) < 1.0e-11_dp* &
            max(1.0_dp, st%variance), 'brute-force Var(E)')
        call require(inside(st%log_z_bound, log(z)) .and. inside(st%energy_bound, &
            mean) .and. inside(st%variance_bound, e2/z - mean**2), &
            'brute-force containment')
        do i = 1, 3
            shell_occ = 0
            count = 0
            do j = 1, 19
                if (level(j) == nint(spec%level(i))) then
                    shell_occ = shell_occ + occ(j)/z
                    count = count + 1
                end if
            end do
            shell_occ = shell_occ/real(count, qp)
            call require(abs(st%occupation(i) - real(shell_occ, dp)) < 1.0e-13_dp .and. &
                inside(st%occupation_bound(i), shell_occ), 'brute-force occupation')
        end do
        call same_spin_pair(spec, st, radii, g, ok, 0.0_dp, trunc)
        production = g
        call require(ok .and. maxval(abs(production - real(pair/(z*n*n), dp))) < &
                     1.0e-12_dp, 'brute-force pair correlation')
    end subroutine brute_box

    !> All C(10,n) determinants of trap levels nx+ny+nz<=2.
    subroutine brute_trap(n, beta)
        integer, intent(in) :: n
        real(dp), intent(in) :: beta
        integer :: level(10), m, i, j, l, mask
        real(qp) :: z, e, w, occ(3), energy
        m = 0
        do i = 0, 2
            do j = 0, 2 - i
                do l = 0, 2 - i - j
                    m = m + 1
                    level(m) = i + j + l
                end do
            end do
        end do
        z = 0; e = 0; occ = 0
        do mask = 0, 2**10 - 1
            if (popcnt(mask) /= n) cycle
            energy = 0
            do i = 1, 10
                if (btest(mask, i - 1)) energy = energy + real(level(i), qp)
            end do
            w = exp(-real(beta, qp)*energy)
            z = z + w; e = e + w*energy
            do i = 1, 10
                if (btest(mask, i - 1)) occ(level(i) + 1) = occ(level(i) + 1) + w
            end do
        end do
        occ = occ/(z*real([1, 3, 6], qp))
        call solve_canonical(spec, n, beta, st, ok)
        call require(ok .and. st%certified, 'brute-force trap certificate')
        call require(abs(st%log_z - real(log(z), dp)) < 1.0e-12_dp .and. &
            abs(st%energy - real(e/z, dp)) < 1.0e-12_dp .and. &
            maxval(abs(st%occupation - real(occ, dp))) < 1.0e-13_dp, 'brute-force trap')
        call require(inside(st%energy_bound, e/z) .and. all([(inside( &
            st%occupation_bound(i), occ(i)), i=1, 3)]), 'brute-force trap containment')
    end subroutine brute_trap

    !> Max deviation of the ground-state hole from 1-9*(j1(x)/x)**2, x=kF*r.
    subroutine cold_hole(n, deviation)
        integer, intent(in) :: n
        real(dp), intent(out) :: deviation
        real(dp) :: r(101), gr(101)
        real(qp) :: x, j1
        type(canonical_state_t) :: ground
        integer :: i
        do i = 1, 101
            r(i) = 0.5_dp*real(i - 1, dp)/100.0_dp
        end do
        call solve_ground(wide, n, ground, ok)
        call same_spin_pair(wide, ground, r, gr, ok, 0.0_dp, trunc)
        call require(ok, 'ground-state pair correlation')
        deviation = abs(gr(1))
        do i = 2, 101
            x = (6.0_qp*acos(-1.0_qp)**2*real(n, qp))**(1.0_qp/3.0_qp)*real(r(i), qp)
            j1 = 3.0_qp*(sin(x) - x*cos(x))/x**3
            deviation = max(deviation, abs(gr(i) - real(1.0_qp - j1**2, dp)))
        end do
    end subroutine cold_hole
end program test_fermion_ensemble
