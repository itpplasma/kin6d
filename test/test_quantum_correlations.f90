program test_quantum_correlations
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_quantum_correlations, only: quantum_wick_pair, quantum_wick_pair_jvp, &
        quantum_wick_pair_vjp, quantum_wick_pair_interval, quantum_cumulant_pair, &
        quantum_cumulant_pair_jvp, quantum_cumulant_pair_vjp, &
        quantum_cumulant_pair_interval, quantum_pauli, &
        pauli_certified, pauli_violated, pauli_undetermined
    implicit none
    integer, parameter :: m = 4, d = 6, masks(d) = [3,5,6,9,10,12]
    complex(dp) :: rho(d,d), gamma(m,m), two(m,m,m,m), wick, cumulant
    complex(dp) :: u(4), du(4), expected, tangent, g, dg
    complex(qp) :: uq(4), expected_q
    real(dp) :: wr, wi, cr, ci, jr, ji, v(8), dv(8), grad(10), lhs, rhs
    real(dp) :: gr(m,m), gi(m,m), trace_two, trace_wick, partial, eps
    real(dp) :: plus_r, plus_i, minus_r, minus_i
    type(interval_t) :: args(10), ir, ii
    logical :: ok
    integer :: i,j,k,l,a,status,corner,p

    ! Independent CAR matrix oracle, not the Wick formula under test.
    rho = 0; rho(1,1) = 1
    call form_rdms(rho,gamma,two)
    call quantum_pauli(real(gamma,dp),aimag(gamma),status,ok)
    call check(ok .and. status == pauli_certified,'diagonal Slater Pauli endpoints')
    call check(abs(sum([(gamma(i,i),i=1,m)])-2) < 1e-15_dp,'Slater Tr gamma=N')
    do l = 1, m
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    call wick_entry(gamma,i,j,k,l,wick)
                    call check(abs(wick-two(i,j,k,l)) < 1e-15_dp, &
                        'Slater CAR equals Wick')
                end do
            end do
        end do
    end do
    ! Complex single Slater determinant: orbital1 wedge (orbital2+i orbital3)/sqrt2.
    rho = 0; rho(1,1) = 0.5_dp; rho(2,2) = 0.5_dp
    rho(1,2) = cmplx(0.0_dp,-0.5_dp,dp); rho(2,1) = conjg(rho(1,2))
    call form_rdms(rho,gamma,two)
    call quantum_pauli(real(gamma,dp),aimag(gamma),status,ok)
    call check(ok .and. status == pauli_undetermined, &
        'physical rotated Slater may fail sufficient Gershgorin after rounding')
    do l = 1, m
        do k = 1, m
            do j = 1, m
                do i = 1, m
                    call wick_entry(gamma,i,j,k,l,wick)
                    call check(abs(wick-two(i,j,k,l)) < 1e-15_dp, &
                        'complex Slater CAR equals Wick')
                end do
            end do
        end do
    end do
    ! Non-Gaussian fixed-N coherent state, with fractional occupations.
    rho = 0; rho(1,1) = 0.5_dp; rho(6,6) = 0.5_dp
    rho(1,6) = 0.5_dp; rho(6,1) = 0.5_dp
    call form_rdms(rho,gamma,two)
    call quantum_pauli(real(gamma,dp),aimag(gamma),status,ok)
    call check(ok .and. status == pauli_certified,'correlated fractional Pauli')
    call cumulant_entry(gamma,two,1,2,3,4,cumulant)
    call check(abs(cumulant-0.5_dp) < 1e-15_dp,'beyond-Slater pair coherence')
    trace_two = 0; trace_wick = 0
    do i = 1, m
        partial = 0
        do a = 1, m
            partial = partial+real(two(i,a,i,a),dp)
            trace_two = trace_two+real(two(i,a,i,a),dp)
            call wick_entry(gamma,i,a,i,a,wick)
            trace_wick = trace_wick+real(wick,dp)
        end do
        call check(abs(partial-real(gamma(i,i),dp)) < 1e-15_dp, &
            'fixed-N contraction')
    end do
    call check(abs(trace_two-2) < 1e-15_dp,'Tr Gamma=N(N-1)')
    call check(abs(trace_wick-3) < 1e-15_dp, &
        'fractional Wick grandcanonical number variance')
    partial = 0
    do a = 1, m
        call wick_entry(gamma,1,a,1,a,wick)
        partial = partial+real(wick,dp)
    end do
    call check(abs(partial-0.75_dp) < 1e-15_dp, &
        'fractional Wick contraction differs fixed N')
    rho(1,6) = 0; rho(6,1) = 0
    call form_rdms(rho,gamma,two)
    call cumulant_entry(gamma,two,1,2,3,4,cumulant)
    call check(abs(cumulant) < 1e-15_dp,'incoherent mixture loses pair coherence')

    ! Independent native complex algebra and analytic directional product derivative.
    u = [cmplx(0.25_dp,0.125_dp,dp),cmplx(-0.5_dp,0.25_dp,dp), &
        cmplx(0.125_dp,-0.25_dp,dp),cmplx(0.5_dp,0.0625_dp,dp)]
    du = [cmplx(-0.0625_dp,0.125_dp,dp),cmplx(0.03125_dp,-0.0625_dp,dp), &
        cmplx(0.0625_dp,0.03125_dp,dp),cmplx(-0.125_dp,0.0625_dp,dp)]
    v(1:8:2) = real(u,dp); v(2:8:2) = aimag(u)
    dv(1:8:2) = real(du,dp); dv(2:8:2) = aimag(du)
    g = cmplx(0.375_dp,-0.125_dp,dp); dg = cmplx(0.0625_dp,0.03125_dp,dp)
    expected = u(1)*u(2)-u(3)*u(4)
    tangent = du(1)*u(2)+u(1)*du(2)-du(3)*u(4)-u(3)*du(4)
    call quantum_wick_pair(v(1),v(2),v(3),v(4),v(5),v(6),v(7),v(8),wr,wi)
    call check(abs(cmplx(wr,wi,dp)-expected) < 1e-15_dp,'complex Wick oracle')
    call quantum_wick_pair_jvp(v(1),v(2),v(3),v(4),v(5),v(6),v(7),v(8), &
        dv(1),dv(2),dv(3),dv(4),dv(5),dv(6),dv(7),dv(8),jr,ji)
    call check(abs(cmplx(jr,ji,dp)-tangent) < 1e-15_dp,'complex derivative oracle')
    call quantum_wick_pair_vjp(v(1),v(2),v(3),v(4),v(5),v(6),v(7),v(8), &
        0.25_dp,-0.125_dp,grad(1),grad(2),grad(3),grad(4),grad(5),grad(6), &
        grad(7),grad(8))
    call check(abs(sum(grad(1:8)*dv)-(0.25_dp*jr-0.125_dp*ji)) < 1e-15_dp, &
        'Wick real pairing')
    call quantum_cumulant_pair(real(g,dp),aimag(g),v(1),v(2),v(3),v(4), &
        v(5),v(6),v(7),v(8),cr,ci)
    call check(abs(cmplx(cr,ci,dp)-(g-expected)) < 1e-15_dp, &
        'cumulant complex oracle')
    call quantum_cumulant_pair_jvp(real(g,dp),aimag(g),v(1),v(2),v(3),v(4), &
        v(5),v(6),v(7),v(8),real(dg,dp),aimag(dg),dv(1),dv(2),dv(3),dv(4), &
        dv(5),dv(6),dv(7),dv(8),jr,ji)
    call check(abs(cmplx(jr,ji,dp)-(dg-tangent)) < 1e-15_dp,'cumulant tangent oracle')
    call quantum_cumulant_pair_vjp(real(g,dp),aimag(g),v(1),v(2),v(3),v(4), &
        v(5),v(6),v(7),v(8),0.25_dp,-0.125_dp,grad(1),grad(2),grad(3),grad(4), &
        grad(5),grad(6),grad(7),grad(8),grad(9),grad(10))
    lhs = 0.25_dp*jr-0.125_dp*ji
    rhs = grad(1)*real(dg,dp)+grad(2)*aimag(dg)+sum(grad(3:10)*dv)
    call check(abs(lhs-rhs) < 1e-15_dp,'cumulant real pairing')
    eps = 2.0_dp**(-16)
    call quantum_cumulant_pair(real(g+eps*dg,dp),aimag(g+eps*dg), &
        v(1)+eps*dv(1),v(2)+eps*dv(2),v(3)+eps*dv(3),v(4)+eps*dv(4), &
        v(5)+eps*dv(5),v(6)+eps*dv(6),v(7)+eps*dv(7),v(8)+eps*dv(8),plus_r,plus_i)
    call quantum_cumulant_pair(real(g-eps*dg,dp),aimag(g-eps*dg), &
        v(1)-eps*dv(1),v(2)-eps*dv(2),v(3)-eps*dv(3),v(4)-eps*dv(4), &
        v(5)-eps*dv(5),v(6)-eps*dv(6),v(7)-eps*dv(7),v(8)-eps*dv(8),minus_r,minus_i)
    call check(abs((plus_r-minus_r)/(2*eps)-jr) < 1e-13_dp &
        .and. abs((plus_i-minus_i)/(2*eps)-ji) < 1e-13_dp,'cumulant finite difference')

    ! All 256 dyadic zero-crossing complex box corners; binary128 oracle.
    args(1:8) = interval(-0.25_dp,0.375_dp)
    call quantum_wick_pair_interval(args(1),args(2),args(3),args(4), &
        args(5),args(6),args(7),args(8),ir,ii)
    do corner = 0, 255
        do p = 1, 8
            v(p) = args(p)%lo
            if (btest(corner,p-1)) v(p) = args(p)%hi
        end do
        uq = cmplx(real(v(1:8:2),qp),real(v(2:8:2),qp),qp)
        expected_q = uq(1)*uq(2)-uq(3)*uq(4)
        call check(real(ir%lo,qp) <= real(expected_q,qp) &
            .and. real(ir%hi,qp) >= real(expected_q,qp),'strict Wick real corners')
        call check(real(ii%lo,qp) <= aimag(expected_q) &
            .and. real(ii%hi,qp) >= aimag(expected_q),'strict Wick imag corners')
    end do
    args(9) = interval(real(g,dp)); args(10) = interval(aimag(g))
    call quantum_cumulant_pair_interval(args(9),args(10),args(1),args(2),args(3), &
        args(4),args(5),args(6),args(7),args(8),ir,ii)
    expected_q = cmplx(real(g,qp),real(aimag(g),qp),qp)-expected_q
    call check(real(ir%lo,qp) <= real(expected_q,qp) &
        .and. real(ir%hi,qp) >= real(expected_q,qp),'strict cumulant real')
    call check(real(ii%lo,qp) <= aimag(expected_q) &
        .and. real(ii%hi,qp) >= aimag(expected_q),'strict cumulant imag')

    gr = 0; gi = 0
    do i = 1, m
        gr(i,i) = 0.5_dp
    end do
    gr(1,2) = 0.0625_dp; gr(2,1) = gr(1,2)
    gi(1,3) = 0.03125_dp; gi(3,1) = -gi(1,3)
    call quantum_pauli(gr,gi,status,ok)
    call check(ok .and. status == pauli_certified,'complex Gershgorin sufficient')
    gr(1,1) = 1.125_dp
    call quantum_pauli(gr,gi,status,ok)
    call check(ok .and. status == pauli_violated,'Pauli diagonal violation')
    gr(1,1) = 0.5_dp; gr(1,2) = 0.75_dp; gr(2,1) = 0.75_dp
    call quantum_pauli(gr,gi,status,ok)
    call check(ok .and. status == pauli_undetermined,'insufficient Gershgorin')
    gi(1,1) = 0.125_dp
    call quantum_pauli(gr,gi,status,ok)
    call check(ok .and. status == pauli_violated,'non-Hermitian Pauli')
    gi(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call quantum_pauli(gr,gi,status,ok); call check(.not. ok,'NaN rejection')
    call quantum_pauli(gr(:,1:3),gi,status,ok)
    call check(.not. ok,'shape rejection')
    gr = 0.5_dp; gi = 0
    gr(1,2) = huge(1.0_dp); gr(2,1) = huge(1.0_dp)
    gi(1,2) = huge(1.0_dp); gi(2,1) = -huge(1.0_dp)
    call quantum_pauli(gr,gi,status,ok); call check(.not. ok,'strict radius overflow')
    print '(a)', &
        'PASS: quantum Wick/cumulant CAR, fixed-N distinction, AD, intervals, Pauli'
contains
    subroutine check(condition,label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        print '(a)', 'FAIL: '//label
        error stop 1
    end subroutine
    function word_matrix(orbitals,creation) result(matrix)
        integer, intent(in) :: orbitals(:)
        logical, intent(in) :: creation(:)
        real(dp) :: matrix(d,d)
        integer :: ket,bra,p,bits,signum,bit
        logical :: live
        matrix = 0
        do ket = 1, d
            bits = masks(ket); signum = 1; live = .true.
            do p = size(orbitals), 1, -1
                bit = orbitals(p)-1
                if (btest(bits,bit) .eqv. creation(p)) then
                    live = .false.
                    exit
                end if
                if (mod(popcnt(iand(bits,2**bit-1)),2) /= 0) signum = -signum
                if (creation(p)) then
                    bits = ibset(bits,bit)
                else
                    bits = ibclr(bits,bit)
                end if
            end do
            if (.not. live) cycle
            do bra = 1, d
                if (masks(bra) == bits) matrix(bra,ket) = real(signum,dp)
            end do
        end do
    end function
    subroutine form_rdms(density,gamma,two)
        complex(dp), intent(in) :: density(d,d)
        complex(dp), intent(out) :: gamma(m,m), two(m,m,m,m)
        integer :: i,j,k,l
        do j = 1, m
            do i = 1, m
                gamma(i,j) = sum(density*transpose(word_matrix([j,i],[.true.,.false.])))
                do l = 1, m
                    do k = 1, m
                        two(i,j,k,l) = sum(density*transpose(word_matrix([k,l,j,i], &
                            [.true.,.true.,.false.,.false.])))
                    end do
                end do
            end do
        end do
    end subroutine
    subroutine wick_entry(gamma,i,j,k,l,value)
        complex(dp), intent(in) :: gamma(m,m)
        integer, intent(in) :: i,j,k,l
        complex(dp), intent(out) :: value
        real(dp) :: r,s
        call quantum_wick_pair(real(gamma(i,k),dp),aimag(gamma(i,k)), &
            real(gamma(j,l),dp),aimag(gamma(j,l)),real(gamma(i,l),dp), &
            aimag(gamma(i,l)),real(gamma(j,k),dp),aimag(gamma(j,k)),r,s)
        value = cmplx(r,s,dp)
    end subroutine
    subroutine cumulant_entry(gamma,two,i,j,k,l,value)
        complex(dp), intent(in) :: gamma(m,m), two(m,m,m,m)
        integer, intent(in) :: i,j,k,l
        complex(dp), intent(out) :: value
        real(dp) :: r,s
        call quantum_cumulant_pair(real(two(i,j,k,l),dp),aimag(two(i,j,k,l)), &
            real(gamma(i,k),dp),aimag(gamma(i,k)),real(gamma(j,l),dp), &
            aimag(gamma(j,l)),real(gamma(i,l),dp),aimag(gamma(i,l)), &
            real(gamma(j,k),dp),aimag(gamma(j,k)),r,s)
        value = cmplx(r,s,dp)
    end subroutine
end program
