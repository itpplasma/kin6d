module vrmhd_structural_spectrum
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_jacobian, &
        vrmhd_section_basis, vrmhd_section_residual
    use rm_mhd_oracles, only: uniform_waves, oracle_ok
    use fortsparse, only: csc_t, csc_matvec
    implicit none
    private
    public :: structural_spectrum, match_folded, bloch_projection, fourier_content
    public :: bloch_rows, bloch_sandwich, bloch_adjoint, check_bloch_action
    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*), info
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
        end subroutine
        subroutine zgecon(norm, n, a, lda, anorm, rcond, work, rwork, info)
            import dp
            character, intent(in) :: norm
            integer, intent(in) :: n, lda
            complex(dp), intent(in) :: a(lda, *)
            real(dp), intent(in) :: anorm
            real(dp), intent(out) :: rcond, rwork(*)
            complex(dp), intent(out) :: work(*)
            integer, intent(out) :: info
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
        subroutine zpotrf(uplo, n, a, lda, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine zgeqrf(m, n, a, lda, tau, work, lwork, info)
            import dp
            integer, intent(in) :: m, n, lda, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(out) :: tau(*), work(*)
            integer, intent(out) :: info
        end subroutine
        subroutine zungqr(m, n, k, a, lda, tau, work, lwork, info)
            import dp
            integer, intent(in) :: m, n, k, lda, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(in) :: tau(*)
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

    function bloch_projection(transform, image) result(projected)
        complex(dp), intent(in) :: transform(:, :)
        complex(qp), intent(in) :: image(:)
        complex(dp) :: projected(size(transform,2))
        integer, allocatable :: rows(:)
        integer :: i, j
        ! Exact zero support avoids dense quadruple-precision work on Bloch columns.
        do j = 1, size(transform,2)
            rows = pack([(i,i=1,size(transform,1))],abs(transform(:,j)) > 0)
            projected(j) = cmplx(dot_product(cmplx(transform(rows,j),kind=qp), &
                image(rows)),kind=dp)
        end do
    end function

    subroutine bloch_rows(transform, column, phase, norm)
        !! The single Bloch column and phase of each coefficient row.
        complex(dp), intent(in) :: transform(:, :)
        integer, allocatable, intent(out) :: column(:)
        complex(dp), allocatable, intent(out) :: phase(:)
        complex(qp), allocatable, intent(out) :: norm(:)
        integer :: i, j
        allocate(column(size(transform,1)), phase(size(transform,1)), &
            norm(size(transform,2)))
        column = 0; phase = 0; norm = 0
        do j = 1, size(transform,2)
            do i = 1, size(transform,1)
                if (abs(transform(i,j)) > 0) then
                    call check(column(i) == 0, 'one Bloch column per coefficient')
                    column(i) = j; phase(i) = transform(i,j)
                    norm(j) = norm(j) + abs(cmplx(transform(i,j),kind=qp))**2
                end if
            end do
        end do
        call check(all(column > 0), 'every coefficient in a Bloch column')
    end subroutine

    function bloch_adjoint(column, phase, image, nb) result(projected)
        !! T^H image for a transform with one phase per row.
        integer, intent(in) :: column(:), nb
        complex(dp), intent(in) :: phase(:)
        complex(qp), intent(in) :: image(:)
        complex(qp) :: projected(nb)
        integer :: i
        projected = 0
        do i = 1, size(column)
            projected(column(i)) = projected(column(i)) &
                + conjg(cmplx(phase(i),kind=qp))*image(i)
        end do
    end function

    function bloch_sandwich(a, column, phase, nb) result(block)
        !! T^H A T in one pass over A, with the quadruple-precision sums of
        !! the column-by-column projection it replaces.
        type(csc_t), intent(in) :: a
        integer, intent(in) :: column(:), nb
        complex(dp), intent(in) :: phase(:)
        complex(dp) :: block(nb, nb)
        complex(qp), allocatable :: total(:, :)
        complex(qp) :: right
        integer :: i, entry, row
        allocate(total(nb, nb))
        total = 0
        do i = 1, a%ncol
            right = cmplx(phase(i),kind=qp)
            do entry = a%col_ptr(i), a%col_ptr(i+1)-1
                row = a%row_idx(entry)
                total(column(row),column(i)) = total(column(row),column(i)) &
                    + conjg(cmplx(phase(row),kind=qp))*(real(a%val(entry),qp)*right)
            end do
        end do
        block = cmplx(total, kind=dp)
    end function

    subroutine check_bloch_action(op, p, u, column, phase, block, time_rate, label)
        !! Applies the production residual derivative at the section's own
        !! quadrature to the sum of all Bloch columns, whose coefficients are
        !! the row phases, and compares it with the assembled Bloch block.
        type(vrmhd_section_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:)
        integer, intent(in) :: column(:)
        complex(dp), intent(in) :: phase(:), block(:, :)
        logical, intent(in) :: time_rate
        character(*), intent(in) :: label
        real(dp), allocatable :: source(:, :, :)
        real(dp) :: zero(op%ndof), r(op%ndof), real_jv(op%ndof), imag_jv(op%ndof)
        complex(dp) :: difference(size(block, 1))
        real(dp) :: defect
        integer :: status
        allocate(source(6, op%fs%nq, op%fs%mm%topology%n_triangles))
        source = 0; zero = 0
        if (time_rate) then
            call vrmhd_section_residual(op, p, u, zero, zero, real(phase, dp), &
                source, r, real_jv, status)
            call check(status == 0, 'matrix-free real time-mass action')
            call vrmhd_section_residual(op, p, u, zero, zero, aimag(phase), &
                source, r, imag_jv, status)
        else
            call vrmhd_section_residual(op, p, u, zero, real(phase, dp), zero, &
                source, r, real_jv, status)
            call check(status == 0, 'matrix-free real Jacobian action')
            call vrmhd_section_residual(op, p, u, zero, aimag(phase), zero, &
                source, r, imag_jv, status)
        end if
        call check(status == 0, 'matrix-free imaginary action')
        difference = cmplx(bloch_adjoint(column, phase, cmplx(real_jv, imag_jv, qp), &
            size(block, 1)), kind=dp) - sum(block, dim=2)
        defect = sqrt(sum(abs(difference)**2)/(size(block, 2)*sum(abs(block)**2)))
        print '(3a,es20.12)', 'assembled / matrix-free ', label, ' defect: ', defect
        call check(defect <= 1.0e-12_dp, &
            'assembled Bloch block equals the production matrix-free action')
    end subroutine

    subroutine structural_spectrum(op, assembly, p, u, source, transform, column, &
            row_phase, stiffness, mass, ns, na, n, repeat, degree, b0, exact, &
            frequencies, mass_condition)
        type(vrmhd_section_t), intent(in) :: op, assembly
        integer, intent(in) :: column(:)
        complex(dp), intent(in) :: row_phase(:)
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u(:), source(:, :, :), b0(3), exact(3)
        complex(dp), intent(in) :: transform(:, :), stiffness(:, :), mass(:, :)
        integer, intent(in) :: ns, na, n, repeat, degree
        real(dp), intent(out) :: frequencies(3)
        real(dp), intent(out), optional :: mass_condition
        complex(dp), allocatable :: f(:, :), e(:, :), ef(:, :), q(:, :), w(:, :)
        complex(dp), allocatable :: h(:, :), modes(:, :), work(:), a(:, :), values(:)
        complex(dp), allocatable :: curls(:, :), weighted(:, :), full(:)
        complex(dp), allocatable :: kd(:, :), fd(:, :), l(:, :), lift(:, :)
        complex(dp), allocatable :: scalar_gram(:, :), scalar_weighted(:, :)
        complex(dp) :: dummy(1, 1), dummy_right(1, 1)
        real(dp), allocatable :: energies(:), omega(:), rw(:), shares(:,:), captured(:)
        real(dp) :: pressure(3), skew, f_norm, defect, maximum_defect
        real(dp) :: scale, reality, diffusivity, cell_wave, slowest, growth
        real(dp) :: condition, skew_bound, coupling_defect
        integer :: d, nb, i, j, a0, j0, info, ne, nz, static_count
        integer :: offsets(3), status
        type(vrmhd_parameters_t) :: dissipative
        type(csc_t) :: jac
        d = 6*ns + na; nb = size(stiffness, 1)
        allocate(f(d,d), e(d,d), ef(d,d), q(d,d), energies(d), rw(8*d), work(32*d))
        allocate(l(nb,nb),lift(nb,d))
        call evolution(stiffness, mass, ns, na, p%mu0, f, l, lift, &
            condition, coupling_defect)
        e = 0
        allocate(scalar_weighted(op%scalar_fs%n0,ns))
        do j = 1, ns
            full = transform(:op%scalar_fs%n0,j)
            scalar_weighted(:,j) = cmplx( &
                csc_matvec(op%scalar_fs%mass(0),real(full,dp)), &
                csc_matvec(op%scalar_fs%mass(0),aimag(full)),dp)
        end do
        allocate(scalar_gram(ns,ns))
        do j = 1, ns
            scalar_gram(:,j) = bloch_projection( &
                transform(:op%scalar_fs%n0,:ns),cmplx(scalar_weighted(:,j),kind=qp))
        end do
        do i = 1, 3
            e(i*ns+1:(i+1)*ns,i*ns+1:(i+1)*ns) = scalar_gram
        end do
        ! Independent quadratic wave energy: |v|^2 + |curl A|^2/mu0 + p'^2/(gamma p0).
        pressure = [0.6_dp, 1.0_dp, 1.0_dp]
        offsets = [0, 4*ns, 5*ns]
        do i = 1, 3
            do j = 1, 3
                e(offsets(i)+1:offsets(i)+ns,offsets(j)+1:offsets(j)+ns) = &
                    pressure(i)*pressure(j)*scalar_gram
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
        do j = 1, na
            e(6*ns+1:,6*ns+j) = bloch_projection(curls,cmplx(weighted(:,j),kind=qp))
        end do
        ef = cmplx(matmul(cmplx(e,kind=qp),cmplx(f,kind=qp)),kind=dp)
        skew = real(sqrt(sum(abs(cmplx(ef,kind=qp) &
            + transpose(conjg(cmplx(ef,kind=qp))))**2) &
            /sum(abs(cmplx(ef,kind=qp))**2)),dp)
        skew_bound = 100*epsilon(1.0_dp)*condition
        print '(a,es20.12)', 'physical time mass M_r 1-norm condition estimate: ', condition
        if (present(mass_condition)) mass_condition = condition
        print '(a,es20.12)', 'diagnostic 100 epsilon kappa_1(M_r): ', skew_bound
        print '(a,es20.12)', '||E F + (E F)^*||_F / ||E F||_F: ', skew
        print '(a,es20.12)', 'inverse-free ||L_vJ + L_Av^*||_F / ||L_vJ||_F: ', coupling_defect
        call inverse_free_structure(l,mass(:d,:d),e,scalar_gram,ns,na, &
            coupling_defect)
        nz = 2*ns + (repeat*op%fs%k0)**2; ne = d-nz
        call energy_quotient(op,transform,column,row_phase,e,ns,na,repeat,degree, &
            q,energies)
        f_norm = real(sqrt(sum(abs(cmplx(f,kind=qp))**2)),dp)
        maximum_defect = 0
        do j = 1, nz
            defect = real(sqrt(sum(abs(matmul(cmplx(f,kind=qp), &
                cmplx(q(:,j),kind=qp)))**2))/f_norm,dp)
            call check(defect <= 1.0e-12_dp, 'energy-free entropy/gauge kernel is static')
            defect = static_defect(q(:,j))
            maximum_defect = max(maximum_defect,defect)
            call check(defect <= 1.0e-12_dp, 'entropy/gauge null satisfies the full spatial operator')
        end do
        allocate(w(d,ne), h(ne,ne), omega(ne), modes(d,ne), values(d), a(d,d))
        do j = 1, ne
            w(:,j) = q(:,nz+j)/sqrt(energies(nz+j))
        end do
        h = cmplx(cmplx(0.0_qp,1.0_qp,qp)*matmul(transpose(conjg(cmplx(w,kind=qp))), &
            matmul(cmplx(ef,kind=qp),cmplx(w,kind=qp))),kind=dp)
        ! The roundoff Hermitian representative is used only after the structural gate.
        h = (h+transpose(conjg(h)))/2
        call zheev('V','U',ne,h,ne,omega,work,size(work),rw,info)
        call check(info == 0, 'energy-quotient ideal frequencies')
        modes = cmplx(matmul(cmplx(w,kind=qp),cmplx(h,kind=qp)),kind=dp)
        a = f
        call zgeev('N','N',d,a,d,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'diagnostic nonsymmetric ideal eigensolve')
        scale = maxval(abs(values**2))
        reality = maxval(abs(aimag(-values**2)))/scale
        print '(a,es20.12)', 'diagnostic max |Im(omega^2)| / max |omega^2|: ', reality
        print '(a,*(es20.12))', 'energy-quotient signed angular frequencies [1/s]: ', omega
        call fourier_content(op,transform,modes,n,repeat,degree,shares,captured, &
            column,row_phase)
        ! Branch isolation requires at least 16 cells per period.
        call match_folded(omega, n, repeat, degree, b0, exact, frequencies, &
            resolved=n >= 16, shares=shares, captured=captured)
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
        ! Per cell, A modulo gradients has na/repeat**2 - k0**2 directions; the
        ! Lorentz force tests them only with the 2 degree**2 velocity
        ! components normal to B0, so the rest is magnetically static.
        call check(static_count == na - (repeat*op%fs%k0)**2 - 2*ns, &
            'static count agrees with the magnetic quotient minus Lorentz rank')
        diffusivity = 1.0e-3_dp
        dissipative = p
        dissipative%eta0 = p%mu0*diffusivity*0.3_dp**1.5_dp
        dissipative%nu0 = diffusivity*0.3_dp**1.5_dp
        call vrmhd_section_jacobian(assembly,dissipative,u,0*u,source,0.0_dp,jac,status)
        call check(status == 0, 'production resistive and viscous linearization')
        allocate(kd(nb,nb), fd(d,d))
        kd = bloch_sandwich(jac,column,row_phase,nb)
        call check_bloch_action(op,dissipative,u,column,row_phase,kd,.false., &
            'dissipative Jacobian')
        call evolution(kd,mass,ns,na,p%mu0,fd)
        a = fd
        call zgeev('N','N',d,a,d,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'dissipative full spectrum')
        growth = maxval(real(values,dp))
        print '(a,es20.12)', 'maximum dissipative Re(lambda) [1/s]: ', growth
        call check(growth <= 1.0e-12_dp*sqrt(sum(abs(cmplx(fd,kind=qp))**2)), &
            'dissipative spectrum has no growth')
        h = cmplx(matmul(transpose(conjg(cmplx(w,kind=qp))), &
            matmul(cmplx(e,kind=qp), &
            matmul(cmplx(fd,kind=qp),cmplx(w,kind=qp)))),kind=dp)
        call zgeev('N','N',ne,h,ne,values,dummy,1,dummy_right,1,work,size(work),rw,info)
        call check(info == 0, 'dissipative energetic quotient')
        slowest = -maxval(real(values(:ne),dp))
        cell_wave = real(n,dp)
        print '(a,2es20.12)', 'slowest energetic decay [1/s] / (eta/mu0) k_cell^2: ', &
            slowest, slowest/(diffusivity*cell_wave**2)
        call check(slowest > 1.0e-12_dp*sqrt(sum(abs(cmplx(fd,kind=qp))**2)), &
            'every energetic mode, including the ideal static kernel, decays')
    contains
        real(dp) function static_defect(x) result(value)
            complex(dp), intent(in) :: x(:)
            complex(qp) :: full(size(l,1))
            full = matmul(cmplx(lift,kind=qp),cmplx(x,kind=qp))
            value = real(sqrt(sum(abs(matmul(cmplx(l,kind=qp),full))**2)) &
                /sqrt(sum(abs(cmplx(l,kind=qp))**2)*sum(abs(full)**2)),dp)
        end function
    end subroutine

    subroutine fourier_content(op, transform, modes, n, repeat, degree, &
            shares, captured, column, row_phase)
        type(vrmhd_section_t), intent(in) :: op
        complex(dp), intent(in) :: transform(:,:), modes(:,:)
        integer, intent(in) :: n, repeat, degree
        real(dp), allocatable, intent(out) :: shares(:,:), captured(:)
        integer, intent(in), optional :: column(:)
        complex(dp), intent(in), optional :: row_phase(:)
        complex(dp), allocatable :: coefficients(:,:), amplitude(:,:,:)
        complex(dp) :: field(7,size(modes,2)), phase
        real(dp) :: basis(36,op%nlocal), readback(7,op%nlocal), x(2), weight
        real(dp) :: wave(2,(repeat*degree)**2), volume, energy(size(modes,2))
        integer :: t, iq, status, jx, jy, g, j
        allocate(amplitude(7,(repeat*degree)**2,size(modes,2)))
        allocate(shares((repeat*degree)**2,size(modes,2)),captured(size(modes,2)))
        if (present(column)) then
            ! Each row has one Bloch term, so this equals the dense product.
            allocate(coefficients(size(transform,1),size(modes,2)))
            coefficients = 0
            do t = 1, size(column)
                if (column(t) <= size(modes,1)) &
                    coefficients(t,:) = row_phase(t)*modes(column(t),:)
            end do
        else
            coefficients = matmul(transform(:,:size(modes,1)),modes)
        end if
        g = 0
        do jy = -(repeat*degree)/2, (repeat*degree-1)/2
            do jx = -(repeat*degree)/2, (repeat*degree-1)/2
                g = g+1
                wave(:,g) = [1.0_dp+real(n*jx,dp)/repeat,real(n*jy,dp)/repeat]
            end do
        end do
        amplitude = 0; volume = 0; energy = 0
        ! Physical energy fields avoid arbitrary entropy and gauge content.
        do t = 1, op%fs%mm%topology%n_triangles
            do iq = 1, op%fs%nq
                call vrmhd_section_basis(op,t,iq,x,weight,basis,status)
                call check(status == 0, 'Fourier physical field readback')
                readback(:3,:) = basis(2:4,:)
                readback(4:6,:) = basis(25:27,:)
                readback(7,:) = 0.6_dp*basis(1,:)+basis(11,:)+basis(12,:)
                field = matmul(readback,coefficients(op%map(:,t),:))
                volume = volume+weight
                energy = energy+weight*sum(abs(field)**2,dim=1)
                do g = 1, size(wave,2)
                    phase = exp(cmplx(0.0_dp,-dot_product(wave(:,g),x),dp))
                    amplitude(:,g,:) = amplitude(:,g,:)+weight*phase*field
                end do
            end do
        end do
        shares = sum(abs(amplitude)**2,dim=1)
        do j = 1, size(modes,2)
            captured(j) = sum(shares(:,j))/(volume*energy(j))
            call check(sum(shares(:,j)) > tiny(1.0_dp), &
                'nonzero retained Fourier power')
            shares(:,j) = shares(:,j)/sum(shares(:,j))
        end do
    end subroutine

    subroutine inverse_free_structure(l_dp, mr_dp, e_dp, scalar_gram_dp, &
            ns, na, coupling_defect)
        complex(dp), intent(in) :: l_dp(:, :), mr_dp(:, :), e_dp(:, :)
        complex(dp), intent(in) :: scalar_gram_dp(:, :)
        integer, intent(in) :: ns, na
        real(dp), intent(in) :: coupling_defect
        complex(qp) :: l(size(l_dp,1),size(l_dp,2))
        complex(qp) :: mr(size(mr_dp,1),size(mr_dp,2))
        complex(qp) :: e(size(e_dp,1),size(e_dp,2))
        complex(qp) :: scalar_gram(size(scalar_gram_dp,1),size(scalar_gram_dp,2))
        complex(qp) :: compression(ns,3*ns), expected(size(mr,1),size(mr,2))
        complex(qp) :: remaining(size(l,1),size(l,2))
        real(qp), parameter :: pressure(3) = real([0.6_dp,1.0_dp,1.0_dp],qp)
        real(dp) :: defect
        integer :: d, i, j, offset(3)
        l = cmplx(l_dp,kind=qp); mr = cmplx(mr_dp,kind=qp)
        e = cmplx(e_dp,kind=qp); scalar_gram = cmplx(scalar_gram_dp,kind=qp)
        d = 6*ns+na; offset = [0,4*ns,5*ns]
        call check(coupling_defect <= 1.0e-14_dp, &
            'inverse-free Lorentz/induction negative adjoints')
        ! Total internal energy gives p_dot without inverting a spatial mass.
        compression = real(2.0_dp/3,qp)*(l(4*ns+1:5*ns,ns+1:4*ns) &
            + l(5*ns+1:6*ns,ns+1:4*ns))
        do i = 1, 3
            defect = relative_norm(l(ns+1:4*ns,offset(i)+1:offset(i)+ns) &
                + pressure(i)*transpose(conjg(compression)), &
                l(ns+1:4*ns,offset(i)+1:offset(i)+ns))
            print '(a,i2,es20.12)', 'pressure/compression adjoint defect, field ', i, defect
            call check(defect <= 1.0e-14_dp, &
                'inverse-free pressure/compression negative adjoints')
        end do
        do i = 2, 3
            defect = relative_norm(l(offset(i)+1:offset(i)+ns,ns+1:4*ns) &
                - real(0.75_dp,qp)*l(:ns,ns+1:4*ns), &
                l(offset(i)+1:offset(i)+ns,ns+1:4*ns))
            print '(a,i2,es20.12)', 'species advection/compression consistency, field ', i, defect
            call check(defect <= 1.0e-14_dp, 'species compression agrees with continuity')
        end do
        defect = relative_norm(compression-l(:ns,ns+1:4*ns),compression)
        print '(a,es20.12)', 'density advection / pressure compression defect: ', defect
        call check(defect <= 1.0e-14_dp, 'continuity and pressure compression agree')
        remaining = l
        remaining(:ns,ns+1:4*ns) = 0
        remaining(ns+1:4*ns,:ns) = 0
        remaining(ns+1:4*ns,4*ns+1:6*ns) = 0
        remaining(ns+1:4*ns,d+1:) = 0
        remaining(4*ns+1:6*ns,ns+1:4*ns) = 0
        remaining(6*ns+1:d,ns+1:4*ns) = 0
        remaining(d+1:,6*ns+1:) = 0
        defect = relative_norm(remaining,l(ns+1:4*ns,d+1:))
        print '(a,es20.12)', 'remaining uniform-state ideal blocks / Lorentz norm: ', defect
        call check(defect <= 1.0e-14_dp, 'no unpaired uniform-state ideal blocks')

        defect = relative_norm(e(6*ns+1:,6*ns+1:)-l(d+1:,6*ns+1:d), &
            l(d+1:,6*ns+1:d))
        print '(a,es20.12)', 'magnetic energy / Ampere curl-curl Gram defect: ', defect
        call check(defect <= 1.0e-14_dp, 'magnetic energy equals Ampere curl-curl stiffness')
        expected = 0
        expected(:ns,:ns) = scalar_gram
        do i = 1, 3
            expected(i*ns+1:(i+1)*ns,i*ns+1:(i+1)*ns) = &
                e(i*ns+1:(i+1)*ns,i*ns+1:(i+1)*ns)
        end do
        ! Primitive temperature rates include the equilibrium density-rate term.
        do i = 2, 3
            expected(offset(i)+1:offset(i)+ns,:ns) = real(0.45_dp,qp)*scalar_gram
            expected(offset(i)+1:offset(i)+ns,offset(i)+1:offset(i)+ns) = &
                real(1.5_dp,qp)*scalar_gram
        end do
        expected(6*ns+1:,6*ns+1:) = -l(d+1:,d+1:)
        defect = relative_norm(mr-expected,expected)
        print '(a,es20.12)', 'kinetic / thermal / induction time-mass Gram defect: ', defect
        call check(defect <= 1.0e-14_dp, 'time masses equal independent finite-element Grams')
        do i = 1, 3
            do j = 1, 3
                defect = relative_norm(e(offset(i)+1:offset(i)+ns,offset(j)+1:offset(j)+ns) &
                    - pressure(i)*pressure(j)*mr(:ns,:ns), &
                    e(offset(i)+1:offset(i)+ns,offset(j)+1:offset(j)+ns))
                print '(a,2i2,es20.12)', 'thermal pressure-energy / scalar mass defect: ', i,j,defect
                call check(defect <= 1.0e-14_dp, 'pressure energy uses the scalar mass Gram')
            end do
        end do
        defect = relative_norm(e-transpose(conjg(e)),e)
        print '(a,es20.12)', 'energy Gram Hermitian defect: ', defect
        call check(defect <= 1.0e-14_dp, 'energy Gram is Hermitian')
    end subroutine

    subroutine energy_quotient(op, transform, column, row_phase, e, ns, na, repeat, &
            degree, q, energies)
        type(vrmhd_section_t), intent(in) :: op
        complex(dp), intent(in) :: transform(:, :), e(:, :), row_phase(:)
        integer, intent(in) :: column(:)
        integer, intent(in) :: ns, na, repeat, degree
        complex(dp), intent(out) :: q(:, :)
        real(dp), intent(out) :: energies(:)
        complex(dp), allocatable :: z(:,:), tau(:), quotient(:,:), factor(:,:)
        complex(dp) :: work(32*size(e,1))
        complex(dp) :: gradient(op%fs%sp(1)%ndof), scalar(op%fs%n0)
        real(dp) :: rw(8*size(e,1)), defect, maximum_defect, norm_e
        complex(qp) :: projected(size(e,1)), norm(size(e,1))
        integer :: d, nz, ne, ng, a0, j0, i, j, info, orbit
        d = size(e,1); ng = (repeat*op%fs%k0)**2; nz = 2*ns+ng; ne = d-nz
        allocate(z(d,nz), tau(nz), quotient(ne,ne), factor(ne,ne))
        a0 = 6*op%scalar_fs%n0; j0 = a0+op%fs%sp(1)%ndof
        call check(ns == (repeat*degree)**2, 'analytic Bloch scalar dimension')
        call check(all(column(a0+1:j0) > 6*ns .and. column(a0+1:j0) <= 6*ns+na), &
            'A coefficients lie in A Bloch columns')
        norm = 0
        do i = a0+1, j0
            norm(column(i)) = norm(column(i)) + abs(cmplx(row_phase(i),kind=qp))**2
        end do
        ! The gauge scalars are the edge complex's Lagrange space; its Bloch
        ! basis is the one of the orbit components, after the edge columns.
        orbit = 6*ns + na - ng
        z = 0
        do j = 1, ns
            z(j,j) = 1; z(4*ns+j,j) = -0.6_dp
            z(4*ns+j,ns+j) = 1; z(5*ns+j,ns+j) = -1
        end do
        do j = 1, ng
            scalar = transform(a0+op%fs%nn+1:j0,orbit+j)
            call check(any(abs(scalar) > 0), 'gauge Bloch scalar column')
            gradient = cmplx(csc_matvec(op%fs%grad%a,real(scalar,dp)), &
                csc_matvec(op%fs%grad%a,aimag(scalar)),dp)
            ! The quotient sums run over the A rows only, as before.
            projected = bloch_adjoint(column(a0+1:j0), row_phase(a0+1:j0), &
                cmplx(gradient,kind=qp), size(e,1))
            z(6*ns+1:6*ns+na,2*ns+j) = cmplx(projected(6*ns+1:6*ns+na) &
                /norm(6*ns+1:6*ns+na),kind=dp)
        end do
        print '(a,4i6)', 'analytic energy nulls: entropy / gauge / total / cells: ', &
            2*ns,ng,nz,repeat**2
        norm_e = real(sqrt(sum(abs(cmplx(e,kind=qp))**2)),dp)
        maximum_defect = 0
        do j = 1, nz
            defect = real(sqrt(sum(abs(matmul(cmplx(e,kind=qp), &
                cmplx(z(:,j),kind=qp)))**2)) &
                /(norm_e*sqrt(sum(abs(cmplx(z(:,j),kind=qp))**2))),dp)
            maximum_defect = max(maximum_defect,defect)
            call check(defect <= 1.0e-14_dp, 'analytic entropy/gauge basis has zero energy')
        end do
        print '(a,es20.12)', 'maximum ||E Z_j|| / (||E||_F ||Z_j||): ', maximum_defect
        q = 0; q(:,:nz) = z
        call zgeqrf(d,nz,q,d,tau,work,size(work),info)
        call check(info == 0, 'analytic null basis QR')
        do j = 1, nz
            call check(abs(q(j,j)) &
                > 1.0e-14_dp*sqrt(sum(abs(cmplx(z(:,j),kind=qp))**2)), &
                'analytic null basis has full rank')
        end do
        call zungqr(d,d,nz,q,d,tau,work,size(work),info)
        call check(info == 0, 'orthonormal representatives of the energy quotient')
        ! E annihilates Z, so its Euclidean complement represents the E-quotient.
        quotient = cmplx(matmul(transpose(conjg(cmplx(q(:,nz+1:),kind=qp))), &
            matmul(cmplx(e,kind=qp),cmplx(q(:,nz+1:),kind=qp))),kind=dp)
        factor = quotient
        call zpotrf('U',ne,factor,ne,info)
        call check(info == 0, 'energy is positive definite modulo analytic physical nulls')
        energies = 0
        call zheev('V','U',ne,quotient,ne,energies(nz+1:),work,size(work),rw,info)
        call check(info == 0, 'energy quotient eigendecomposition')
        call check(energies(nz+1) > 0, 'energy quotient has strictly positive eigenvalues')
        print '(a,es20.12)', 'quotient energy smallest / largest eigenvalue: ', &
            energies(nz+1)/energies(d)
        q(:,nz+1:) = cmplx(matmul(cmplx(q(:,nz+1:),kind=qp), &
            cmplx(quotient,kind=qp)),kind=dp)
    end subroutine

    real(dp) function relative_norm(difference, reference) result(value)
        complex(qp), intent(in) :: difference(:, :), reference(:, :)
        value = real(sqrt(sum(abs(difference)**2)/sum(abs(reference)**2)),dp)
    end function

    subroutine evolution(stiffness, mass, ns, na, mu0, f, l, lift, condition, coupling_defect)
        complex(dp), intent(in) :: stiffness(:, :), mass(:, :)
        integer, intent(in) :: ns, na
        real(dp), intent(in) :: mu0
        complex(dp), intent(out) :: f(:, :)
        complex(dp), intent(out), optional :: l(:, :), lift(:, :)
        real(dp), intent(out), optional :: condition, coupling_defect
        complex(dp) :: k(size(stiffness,1),size(stiffness,2))
        complex(dp) :: m(size(mass,1),size(mass,2)), jj(na,na)
        complex(dp) :: current(na,size(f,1)), dynamic(size(f,1),size(f,1))
        complex(dp) :: condition_work(2*size(f,1))
        real(dp) :: columns(size(k,1)), rows(size(k,1))
        real(dp) :: mass_norm, reciprocal_condition, condition_rwork(2*size(f,1))
        integer :: j, d, info, pivot(size(k,1))
        d = size(f,1)
        columns = 1; rows = 1
        columns(6*ns+1:d) = sqrt(mu0); columns(d+1:) = 1/sqrt(mu0)
        rows(6*ns+1:) = 1/sqrt(mu0)
        do j = 1, size(k,1)
            k(:,j) = rows*stiffness(:,j)*columns(j)
            m(:,j) = rows*mass(:,j)*columns(j)
        end do
        if (present(coupling_defect)) then
            ! Magnetic energy tests induction with J, paired with the velocity test.
            coupling_defect = real(sqrt(sum(abs(cmplx(k(ns+1:4*ns,d+1:),kind=qp) &
                + transpose(conjg(cmplx(k(6*ns+1:d,ns+1:4*ns),kind=qp))))**2) &
                /sum(abs(cmplx(k(ns+1:4*ns,d+1:),kind=qp))**2)),dp)
        end if
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
        f = cmplx(-cmplx(k(:d,:d),kind=qp) &
            - matmul(cmplx(k(:d,d+1:),kind=qp),cmplx(current,kind=qp)),kind=dp)
        dynamic = m(:d,:d)
        mass_norm = real(maxval(sum(abs(cmplx(dynamic,kind=qp)),dim=1)),dp)
        call zgesv(d,d,dynamic,d,pivot,f,d,info)
        call check(info == 0, 'invert physical time mass')
        if (present(condition)) then
            call zgecon('1',d,dynamic,d,mass_norm,reciprocal_condition, &
                condition_work,condition_rwork,info)
            call check(info == 0, 'estimate physical time mass condition from its LU factor')
            call check(reciprocal_condition > 0, 'physical time mass has finite condition estimate')
            condition = 1/reciprocal_condition
        end if
    end subroutine

    subroutine match_folded(omega, n, repeat, degree, b0, exact, frequencies, resolved, &
            shares, captured)
        real(dp), intent(in) :: omega(:), b0(3), exact(3)
        integer, intent(in) :: n, repeat, degree
        real(dp), intent(out) :: frequencies(3)
        logical, intent(in), optional :: resolved
        real(dp), intent(in), optional :: shares(:,:), captured(:)
        real(dp), allocatable :: oracle(:), resolved_power(:)
        real(dp) :: wave(3), squared(3)
        real(dp) :: cutoff, separation, tolerance, last_counted, next_folded
        logical, allocatable :: candidates(:,:), visited(:), selected(:), retained(:)
        integer :: jx, jy, j, i, count_oracle, status, direction, matched, base(3)
        integer, allocatable :: owner(:,:)
        integer :: expected, computed, nearest, base_alias
        ! Avoid gfortran's extent-order bug: test/reproducers/gfortran_host_extent.f90.
        allocate(oracle(3*(repeat*degree)**2))
        allocate(candidates(size(omega),size(oracle)),visited(size(oracle)), &
            owner(size(oracle),2))
        allocate(selected(size(omega)),retained(size(oracle)),resolved_power(size(omega)))
        tolerance = 1.0e-12_dp*maxval(abs(omega))
        count_oracle = 0; base = 0; base_alias = 0
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
                if (jx == 0 .and. jy == 0) base_alias = count_oracle/3
            end do
        end do
        last_counted = maxval(oracle,mask=oracle <= maxval(exact)+tolerance)
        next_folded = minval(oracle,mask=oracle > last_counted+tolerance)
        if (next_folded >= huge(1.0_dp)) then
            ! The one-sample band needs its adjacent aliases to bracket the cutoff.
            do jy = -1, 1
                do jx = -1, 1
                    wave = [1.0_dp+real(n*jx,dp)/repeat, real(n*jy,dp)/repeat, 0.0_dp]
                    call uniform_waves(wave,b0,1.0_dp,0.6_dp,5.0_dp/3,squared,status)
                    call check(status == oracle_ok, 'cutoff guard-alias oracle')
                    next_folded = min(next_folded, &
                        minval(sqrt(squared),mask=sqrt(squared) > last_counted+tolerance))
                end do
            end do
        end if
        call check(next_folded < huge(1.0_dp), 'folded oracle gap exists above counted band')
        cutoff = (last_counted+next_folded)/2
        print '(a,3es20.12)', 'last counted / next folded / cutoff [1/s]: ', &
            last_counted, next_folded, cutoff
        selected = abs(omega) > tolerance .and. abs(omega) <= cutoff+tolerance
        retained = .true.
        if (n >= 16) then
            call check(present(shares) .and. present(captured), &
                'resolved-mesh matching requires physical Fourier power')
            resolved_power = captured*shares(base_alias,:)
            retained = .false.; retained(base) = .true.
        end if
        if (present(shares)) then
            print '(a,3i6)', 'Fourier case cells / repeat / degree: ', n,repeat,degree
            j = 0
            do jy = -(repeat*degree)/2, (repeat*degree-1)/2
                do jx = -(repeat*degree)/2, (repeat*degree-1)/2
                    j = j+1
                    print '(a,i4,2es20.12)', 'Fourier alias index / G [1/m]: ', &
                        j,real(n*[jx,jy],dp)/repeat
                end do
            end do
            do j = 1, size(omega)
                if (.not. selected(j)) cycle
                print '(a,es20.12)', 'Fourier mode omega [rad/s]: ',omega(j)
                if (present(captured)) print '(a,es20.12)', &
                    'captured L2 fraction: ',captured(j)
                print '(a,*(es20.12))', 'L2 shares per retained alias: ',shares(:,j)
                if (n >= 16) print '(a,es20.12)', &
                    'resolved physical L2 fraction r: ',resolved_power(j)
            end do
        end if
        if (n >= 16) then
            call check(all(pack(resolved_power <= 0.1_dp .or. resolved_power >= 0.9_dp, &
                selected)), 'low-mode resolved power lies outside the ambiguous interval')
            selected = selected .and. resolved_power >= 0.9_dp
        end if
        candidates = .false.; owner = 0
        do i = 1, count_oracle
            if (.not. retained(i)) cycle
            if (oracle(i) > cutoff+tolerance) cycle
            separation = huge(1.0_dp)
            do j = 1, count_oracle
                if (abs(oracle(j)-oracle(i)) <= tolerance) cycle
                separation = min(separation,abs(oracle(j)-oracle(i))/2)
            end do
            candidates(:,i) = abs(abs(omega)-oracle(i)) < separation
        end do
        expected = 2*count(retained .and. oracle <= cutoff+tolerance)
        computed = count(selected)
        print '(a,2i6)', 'selected / expected signed waves below cutoff: ', computed,expected
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
        call check(computed == expected, 'selected oracle multiplicity below cutoff')
        matched = 0
        do j = 1, size(omega)
            if (.not. selected(j)) cycle
            direction = 1
            if (omega(j) > 0) direction = 2
            visited = .false.
            if (assign(j,direction)) matched = matched+1
            print '(a,es20.12,i6)', 'below-cutoff omega [1/s] / candidate aliases: ', &
                omega(j), count(candidates(j,:))
        end do
        print '(a,i6)', 'signed waves with a one-to-one folded branch match: ', matched
        call check(matched == computed, &
            'every selected frequency matches a retained oracle branch')
        call check(matched == expected, 'selected oracle multiplicity below cutoff')
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
