module vrmhd_structural_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_jacobian
    use rm_mhd_oracles, only: uniform_waves, oracle_ok
    use fortsparse, only: csc_t, csc_matvec
    implicit none
    private
    public :: structural_spectrum, match_folded
    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine
        subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
            import dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            complex(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), rwork(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, &
                work, lwork, rwork, info)
            import dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(a)', message
        error stop 'RM2 structural spectrum gate'
    end subroutine

    subroutine structural_spectrum(op, p, u, source, transform, stiffness, &
            mass, ns, na, n, repeat, degree, b0, exact, frequencies)
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:), source(:, :, :), b0(3), exact(3)
        complex(dp), intent(in) :: transform(:, :), stiffness(:, :), mass(:, :)
        integer, intent(in) :: ns, na, n, repeat, degree
        real(dp), intent(out) :: frequencies(3)
        complex(dp), allocatable :: f(:, :), e(:, :), ef(:, :), q(:, :), w(:, :)
        complex(dp), allocatable :: h(:, :), modes(:, :), work(:), a(:, :), values(:)
        complex(dp), allocatable :: curls(:, :), weighted(:, :), full(:), image(:)
        complex(dp), allocatable :: kd(:, :), fd(:, :), l(:, :), lift(:, :)
        complex(dp) :: dummy(1, 1), dummy_right(1, 1)
        real(dp), allocatable :: energies(:), omega(:), rw(:)
        real(dp) :: pressure(3), skew, floor, f_norm, defect, maximum_defect
        real(dp) :: scale, reality, diffusivity, cell_wave, slowest, growth
        integer :: d, nb, i, j, row, entry, a0, j0, info, ne, nz, static_count
        integer :: offsets(3), status
        type(vrmhd_parameters_t) :: dissipative
        type(csc_t) :: jac
        d = 6*ns + na; nb = size(stiffness, 1)
        allocate(f(d,d), e(d,d), ef(d,d), q(d,d), energies(d), rw(8*d), work(32*d))
        allocate(l(nb,nb),lift(nb,d))
        call evolution(stiffness, mass, ns, na, p%mu0, f, l, lift)
        e = 0
        e(ns+1:4*ns,ns+1:4*ns) = mass(ns+1:4*ns,ns+1:4*ns)
        ! Independent quadratic wave energy: |v|^2 + |curl A|^2/mu0 + p'^2/(gamma p0).
        pressure = [0.6_dp, 1.0_dp, 1.0_dp]
        offsets = [0, 4*ns, 5*ns]
        do i = 1, 3
            do j = 1, 3
                e(offsets(i)+1:offsets(i)+ns,offsets(j)+1:offsets(j)+ns) = &
                    pressure(i)*pressure(j)*mass(:ns,:ns)
            end do
        end do
        a0 = 6*op%scalar_fs%n0; j0 = a0 + op%fs%sp(1)%ndof
        allocate(curls(op%fs%sp(2)%ndof,na), weighted(op%fs%sp(2)%ndof,na))
        do j = 1, na
            full = transform(a0+1:j0,6*ns+j)
            curls(:,j) = cmplx(csc_matvec(op%fs%curl%a,real(full,dp)), &
                csc_matvec(op%fs%curl%a,aimag(full)),dp)
            weighted(:,j) = cmplx(csc_matvec(op%fs%mass(2),real(curls(:,j),dp)), &
                csc_matvec(op%fs%mass(2),aimag(curls(:,j))),dp)
        end do
        ! A is measured in sqrt(mu0) units in evolution(), cancelling 1/mu0 here.
        e(6*ns+1:,6*ns+1:) = matmul(transpose(conjg(curls)),weighted)
        ef = matmul(e,f)
        skew = sqrt(sum(abs(ef+transpose(conjg(ef)))**2)/sum(abs(ef)**2))
        print '(a,es20.12)', '||E F + (E F)^*||_F / ||E F||_F: ', skew
        call check(skew <= 1.0e-12_dp, 'ideal energy skew-adjointness exceeds 1e-12')
        q = e
        call zheev('V','U',d,q,d,energies,work,size(work),rw,info)
        call check(info == 0, 'quadratic wave energy eigendecomposition')
        floor = 1024*epsilon(1.0_dp)*maxval(energies)
        call check(minval(energies) >= -floor, 'nonnegative conserved quadratic energy')
        ne = count(energies > floor); nz = d-ne
        call check(nz == 3*ns, 'two entropy fields and one gauge field have zero energy')
        f_norm = sqrt(sum(abs(f)**2))
        maximum_defect = 0
        do j = 1, nz
            defect = sqrt(sum(abs(matmul(f,q(:,j)))**2))/f_norm
            call check(defect <= 1.0e-12_dp, 'energy-free entropy/gauge kernel is static')
            defect = static_defect(q(:,j))
            maximum_defect = max(maximum_defect,defect)
            call check(defect <= 1.0e-12_dp, 'entropy/gauge null satisfies the full spatial operator')
        end do
        allocate(w(d,ne), h(ne,ne), omega(ne), modes(d,ne), values(d), a(d,d))
        do j = 1, ne
            w(:,j) = q(:,nz+j)/sqrt(energies(nz+j))
        end do
        h = cmplx(0.0_dp,1.0_dp,dp)*matmul(transpose(conjg(w)),matmul(ef,w))
        ! The roundoff Hermitian representative is used only after the structural gate.
        h = (h+transpose(conjg(h)))/2
        call zheev('V','U',ne,h,ne,omega,work,size(work),rw,info)
        call check(info == 0, 'energy-quotient ideal frequencies')
        modes = matmul(w,h)
        a = f
        call zgeev('N','N',d,a,d,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'diagnostic nonsymmetric ideal eigensolve')
        scale = maxval(abs(values**2))
        reality = maxval(abs(aimag(-values**2)))/scale
        print '(a,es20.12)', 'diagnostic max |Im(omega^2)| / max |omega^2|: ', reality
        print '(a,*(es20.12))', 'energy-quotient signed angular frequencies [1/s]: ', omega
        ! The 4..32 ladder is resolved for branch isolation from 16 cells onward.
        call match_folded(omega, n, repeat, degree, b0, exact, frequencies, &
            resolved=n >= 16)
        static_count = 0
        do j = 1, ne
            if (abs(omega(j)) > 1.0e-12_dp*maxval(abs(omega))) cycle
            static_count = static_count+1
            defect = static_defect(modes(:,j))
            maximum_defect = max(maximum_defect,defect)
            call check(defect <= 1.0e-12_dp, 'energetic null mode is exactly static')
        end do
        print '(a,3i6)', 'static total / energetic / cells: ', nz+static_count, &
            static_count, repeat**2
        print '(a,f12.6)', 'energetic static modes per cell: ', &
            real(static_count,dp)/repeat**2
        print '(a,es20.12)', 'maximum ||L x|| / (||L||_F ||x||): ', maximum_defect
        call check(static_count == degree*repeat**2, 'static count agrees with pairing study')
        diffusivity = 1.0e-3_dp
        dissipative = p
        dissipative%eta0 = p%mu0*diffusivity*0.3_dp**1.5_dp
        dissipative%nu0 = diffusivity*0.3_dp**1.5_dp
        call vrmhd_section_jacobian(op,dissipative,u,0*u,source,0.0_dp,jac,status)
        call check(status == 0, 'production resistive and viscous linearization')
        allocate(kd(nb,nb), image(op%ndof), fd(d,d))
        do j = 1, nb
            image = 0
            do i = 1, op%ndof
                do entry = jac%col_ptr(i), jac%col_ptr(i+1)-1
                    row = jac%row_idx(entry)
                    image(row) = image(row)+jac%val(entry)*transform(i,j)
                end do
            end do
            kd(:,j) = matmul(transpose(conjg(transform)),image)
        end do
        call evolution(kd,mass,ns,na,p%mu0,fd)
        a = fd
        call zgeev('N','N',d,a,d,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'dissipative full spectrum')
        growth = maxval(real(values,dp))
        print '(a,es20.12)', 'maximum dissipative Re(lambda) [1/s]: ', growth
        call check(growth <= 1.0e-12_dp*sqrt(sum(abs(fd)**2)), 'dissipative spectrum has no growth')
        h = matmul(transpose(conjg(w)),matmul(e,matmul(fd,w)))
        call zgeev('N','N',ne,h,ne,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'dissipative energetic quotient')
        slowest = -maxval(real(values(:ne),dp))
        cell_wave = real(n,dp)
        print '(a,2es20.12)', 'slowest energetic decay [1/s] / (eta/mu0) k_cell^2: ', &
            slowest, slowest/(diffusivity*cell_wave**2)
        call check(slowest > 1.0e-12_dp*sqrt(sum(abs(fd)**2)), &
            'every energetic mode, including the ideal static kernel, decays')
    contains
        real(dp) function static_defect(x) result(value)
            complex(dp), intent(in) :: x(:)
            complex(dp) :: full(size(l,1))
            full = matmul(lift,x)
            value = sqrt(sum(abs(matmul(l,full))**2)) &
                /sqrt(sum(abs(l)**2)*sum(abs(full)**2))
        end function
    end subroutine

    subroutine evolution(stiffness, mass, ns, na, mu0, f, l, lift)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :)
        integer, intent(in) :: ns, na
        real(dp), intent(in) :: mu0
        complex(dp), intent(out) :: f(:, :)
        complex(dp), intent(out), optional :: l(:, :), lift(:, :)
        complex(dp) :: k(size(stiffness,1),size(stiffness,2))
        complex(dp) :: m(size(mass,1),size(mass,2)), jj(na,na)
        complex(dp) :: current(na,size(f,1)), dynamic(size(f,1),size(f,1))
        real(dp) :: columns(size(k,1)), rows(size(k,1))
        integer :: j, d, info, pivot(size(k,1))
        d = size(f,1)
        columns = 1; rows = 1
        columns(6*ns+1:d) = sqrt(mu0); columns(d+1:) = 1/sqrt(mu0)
        rows(6*ns+1:) = 1/sqrt(mu0)
        do j = 1, size(k,1)
            k(:,j) = rows*stiffness(:,j)*columns(j)
            m(:,j) = rows*mass(:,j)*columns(j)
        end do
        jj = k(d+1:,d+1:); current = -k(d+1:,:d)
        call zgesv(na,d,jj,na,pivot,current,na,info)
        call check(info == 0, 'eliminate algebraic current through its mass')
        if (present(l)) l = -k
        if (present(lift)) then
            lift = 0
            do j = 1, d
                lift(j,j) = 1
            end do
            lift(d+1:,:) = current
        end if
        f = -k(:d,:d)-matmul(k(:d,d+1:),current)
        dynamic = m(:d,:d)
        call zgesv(d,d,dynamic,d,pivot,f,d,info)
        call check(info == 0, 'invert physical time mass')
    end subroutine

    subroutine match_folded(omega, n, repeat, degree, b0, exact, frequencies, resolved)
        real(dp), intent(in) :: omega(:), b0(3), exact(3)
        integer, intent(in) :: n, repeat, degree
        real(dp), intent(out) :: frequencies(3)
        logical, intent(in), optional :: resolved
        real(dp), allocatable :: oracle(:)
        real(dp) :: wave(3), squared(3)
        real(dp) :: cutoff, separation, tolerance
        logical, allocatable :: candidates(:,:), visited(:)
        integer :: jx, jy, j, i, count_oracle, status, direction, matched, base(3)
        integer, allocatable :: owner(:,:)
        integer :: expected, computed, nearest
        ! Avoid gfortran's extent-order bug: test/reproducers/gfortran_host_extent.f90.
        allocate(oracle(3*(repeat*degree)**2))
        allocate(candidates(size(omega),size(oracle)),visited(size(oracle)), &
            owner(size(oracle),2))
        cutoff = maxval(exact); tolerance = 1.0e-12_dp*maxval(abs(omega))
        count_oracle = 0; base = 0
        ! The resolved band has repeat*degree samples in each periodic direction.
        do jy = -(repeat*degree)/2, (repeat*degree-1)/2
            do jx = -(repeat*degree)/2, (repeat*degree-1)/2
                wave = [1.0_dp+real(n*jx,dp)/repeat, real(n*jy,dp)/repeat, 0.0_dp]
                call uniform_waves(wave,b0,1.0_dp,0.6_dp,5.0_dp/3,squared,status)
                call check(status == oracle_ok, 'Bloch-folded uniform-wave oracle')
                do i = 1, 3
                    count_oracle = count_oracle+1
                    oracle(count_oracle) = sqrt(squared(i))
                    if (jx == 0 .and. jy == 0) base(i) = count_oracle
                end do
            end do
        end do
        candidates = .false.; owner = 0
        do i = 1, count_oracle
            if (oracle(i) > cutoff+tolerance) cycle
            separation = huge(1.0_dp)
            do j = 1, count_oracle
                if (abs(oracle(j)-oracle(i)) <= tolerance) cycle
                separation = min(separation,abs(oracle(j)-oracle(i))/2)
            end do
            candidates(:,i) = abs(abs(omega)-oracle(i)) < separation
        end do
        expected = 2*count(oracle <= cutoff+tolerance)
        computed = count(abs(omega) > tolerance .and. abs(omega) <= cutoff+tolerance)
        print '(a,2i6)', 'computed / expected signed waves below cutoff: ', computed,expected
        call check(computed == expected, 'folded oracle multiplicity below cutoff')
        if (present(resolved)) then
            if (.not. resolved) then
                ! Coarse estimates enter the ladder without asserting branch identity.
                do i = 1, 3
                    nearest = minloc(abs(omega-exact(i)),dim=1,mask=omega > tolerance)
                    call check(nearest > 0, 'positive coarse-mesh frequency exists')
                    frequencies(i) = omega(nearest)
                end do
                return
            end if
        end if
        matched = 0
        do j = 1, size(omega)
            if (abs(omega(j)) <= tolerance .or. abs(omega(j)) > cutoff+tolerance) cycle
            direction = 1
            if (omega(j) > 0) direction = 2
            visited = .false.
            if (assign(j,direction)) matched = matched+1
            print '(a,es20.12,i6)', 'below-cutoff omega [1/s] / candidate aliases: ', &
                omega(j), count(candidates(j,:))
        end do
        print '(a,i6)', 'signed waves with a one-to-one folded branch match: ', matched
        call check(matched == computed, 'every nonzero frequency matches a folded oracle branch')
        call check(matched == expected, 'folded oracle multiplicity below cutoff')
        do i = 1, 3
            call check(owner(base(i),2) > 0, 'positive base-wavenumber branch found')
            frequencies(i) = omega(owner(base(i),2))
        end do
    contains
        recursive logical function assign(mode, sign_index) result(found)
            integer, value, intent(in) :: mode, sign_index
            integer :: alias
            found = .false.
            do alias = 1, count_oracle
                if (.not. candidates(mode,alias) .or. visited(alias)) cycle
                visited(alias) = .true.
                if (owner(alias,sign_index) /= 0) then
                    if (.not. assign(owner(alias,sign_index),sign_index)) cycle
                end if
                owner(alias,sign_index) = mode
                found = .true.
                return
            end do
        end function
    end subroutine
end module
