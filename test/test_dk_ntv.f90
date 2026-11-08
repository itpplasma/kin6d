!> NTV flux kernel of the Legendre reference (#183).
!!   jvp          tangents in nu and omega_E against centred differences;
!!   tangents     tangents in the amplitudes c_m and along (eps, k, offset,
!!                c_J) against centred differences, for both weights and
!!                both streaming models;
!!   line         the circular field through the general field-line provider
!!                reproduces the surface fields bit for bit;
!!   assembly     both drift/drive weights and both streaming models against a
!!                dense Legendre-Galerkin x Fourier solve assembled by Gauss
!!                quadrature and the DFT derivative on a small grid;
!!   kernel FILE  kernels of a surface at the speeds listed in FILE (the
!!                stage-1 gate driver feeds this mode).
program test_dk_ntv
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_dk_ntv, only: ntv_surface_t, ntv_kernel_t, ntv_kernel, ntv_f10, NTV_OK, &
        NTV_WEIGHT_MU, NTV_WEIGHT_FULL, ntv_circular_line_t, NTV_BAD_INPUT
    use fortnum_quadrature, only: gauss_legendre
    implicit none
    character(256) :: check, path
    logical :: ok

    call get_command_argument(1, check)
    select case (trim(check))
    case ('jvp'); ok = jvp()
    case ('tangents'); ok = tangents()
    case ('line'); ok = line()
    case ('ql_dpar'); ok = ql_dpar()
    case ('ql_f10'); ok = ql_f10()
    case ('assembly'); ok = assembly()
    case ('f10_provider'); ok = f10_provider()
    case ('sparse_assembly'); ok = assembly(2)
    case ('sparse_jvp'); ok = jvp(6)
    case ('sparse_tangents'); ok = tangents(6)
    case ('sparse_refinement'); ok = sparse_refinement()
    case ('kernel')
        call get_command_argument(2, path)
        call kernel_table(trim(path))
        ok = .true.
    case default
        print '(a)', 'unknown check '//trim(check)
        ok = .false.
    end select
    if (.not. ok) error stop 1

contains

    !> One (m, n) = (6, 3) harmonic on eps = 0.05, q = 2.2 with collisions
    !> and precession of the stage-1 order.
    function sample() result(sf)
        type(ntv_surface_t) :: sf
        sf%eps = 0.05_dp
        sf%q = 2.2_dp
        sf%k = 1/(sf%q*6.2_dp)
        sf%offset = 0.01_dp
        sf%n = 3
        sf%m = [6]
        sf%c = [(1e-3_dp, 2e-4_dp)]
    end function sample

    logical function jvp(theta_order) result(ok)
        integer, intent(in), optional :: theta_order
        integer :: order, degree
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        real(dp) :: h, fd
        type(ntv_kernel_t) :: k0, up, down
        integer :: status
        order = 0
        if (present(theta_order)) order = theta_order
        degree = merge(12,80,order > 0)
        call ntv_kernel(sample(), v, nu, omega, pre, 33, degree, k0, status, &
            theta_order=order)
        if (status /= NTV_OK) error stop 'kernel failed'
        h = 1e-6_dp*nu
        call ntv_kernel(sample(), v, nu + h, omega, pre, 33, degree, up, status, &
            theta_order=order)
        call ntv_kernel(sample(), v, nu - h, omega, pre, 33, degree, down, status, &
            theta_order=order)
        fd = (up%r - down%r)/(2*h)
        print '(a,2es14.6,es9.1)', 'dR/dnu      ', k0%dr_dnu, fd, abs(k0%dr_dnu/fd - 1)
        ok = abs(k0%dr_dnu/fd - 1) < 1e-6_dp
        h = 1e-6_dp*omega
        call ntv_kernel(sample(), v, nu, omega + h, pre, 33, degree, up, status, &
            theta_order=order)
        call ntv_kernel(sample(), v, nu, omega - h, pre, 33, degree, down, status, &
            theta_order=order)
        fd = (up%r - down%r)/(2*h)
        print '(a,2es14.6,es9.1)', 'dR/domega_E ', k0%dr_domega, fd, abs(k0%dr_domega/fd - 1)
        ok = ok .and. abs(k0%dr_domega/fd - 1) < 1e-6_dp
    end function jvp

    logical function tangents(theta_order) result(ok)
        integer, intent(in), optional :: theta_order
        integer :: order, degree
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        type(ntv_surface_t) :: sf, dir, up, down
        type(ntv_kernel_t) :: k0, kp, km
        real(dp) :: h, fd
        complex(dp) :: dc
        integer :: status, j, ib, iw, part
        order = 0
        if (present(theta_order)) order = theta_order
        degree = merge(12,80,order > 0)
        ok = .true.
        dir%eps = 1; dir%k = 0.3_dp*sample_k(); dir%offset = 0.5_dp; dir%current = 0.7_dp
        do ib = 0, 1
            do iw = 1, 2
                sf = sample()
                sf%m = [5, 6, 7]
                sf%c = [(2e-4_dp, -1e-4_dp), (1e-3_dp, 2e-4_dp), (-3e-4_dp, 5e-5_dp)]
                sf%boozer_stream = ib == 1
                sf%weight = merge(NTV_WEIGHT_MU, NTV_WEIGHT_FULL, iw == 1)
                sf%current = merge(0.0_dp, 0.37_dp, iw == 1)
                call ntv_kernel(sf, v, nu, omega, pre, 33, degree, k0, status, &
                    indicators=.false., direction=dir, &
                    theta_order=order)
                if (status /= NTV_OK) error stop 'kernel failed'
                do j = 1, size(sf%m)
                    do part = 1, 2
                        dc = merge((1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), part == 1)
                        h = 1e-4_dp*abs(sf%c(j))
                        up = sf; up%c(j) = sf%c(j) + h*dc
                        down = sf; down%c(j) = sf%c(j) - h*dc
                        call ntv_kernel(up, v, nu, omega, pre, 33, degree, kp, status, &
                            indicators=.false., theta_order=order)
                        call ntv_kernel(down, v, nu, omega, pre, 33, degree, km, status, &
                            indicators=.false., theta_order=order)
                        fd = (kp%r - km%r)/(2*h)
                        ok = ok .and. near(merge(real(k0%dr_dc(j)), aimag(k0%dr_dc(j)), &
                            part == 1), fd, 'dR/dc')
                    end do
                end do
                h = 1e-5_dp
                up = sf; up%eps = sf%eps + h*dir%eps; up%k = sf%k + h*dir%k
                up%offset = sf%offset + h*dir%offset; up%current = sf%current + h*dir%current
                down = sf; down%eps = sf%eps - h*dir%eps; down%k = sf%k - h*dir%k
                down%offset = sf%offset - h*dir%offset; down%current = sf%current - h*dir%current
                call ntv_kernel(up, v, nu, omega, pre, 33, degree, kp, status, indicators=.false., &
                    theta_order=order)
                call ntv_kernel(down, v, nu, omega, pre, 33, degree, km, status, indicators=.false., &
                    theta_order=order)
                fd = (kp%r - km%r)/(2*h)
                ok = ok .and. near(k0%dr_dsurface, fd, 'dR/dsurface')
            end do
        end do
    end function tangents

    logical function near(got, want, label)
        real(dp), intent(in) :: got, want
        character(*), intent(in) :: label
        print '(a,2es14.6,es9.1)', label//' ', got, want, abs(got/want - 1)
        near = abs(got/want - 1) < 1e-6_dp
    end function near

    logical function line() result(ok)
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        type(ntv_surface_t) :: sf, sl
        type(ntv_kernel_t) :: kf, kl
        integer :: status, ib, iw
        ok = .true.
        do ib = 0, 1
            do iw = 1, 2
                sf = sample()
                sf%m = [5, 6, 7]
                sf%c = [(2e-4_dp, -1e-4_dp), (1e-3_dp, 2e-4_dp), (-3e-4_dp, 5e-5_dp)]
                sf%boozer_stream = ib == 1
                sf%weight = merge(NTV_WEIGHT_MU, NTV_WEIGHT_FULL, iw == 1)
                sf%current = merge(0.0_dp, 0.37_dp, iw == 1)
                sl = sf
                sl%line = ntv_circular_line_t(eps=sf%eps, k=sf%k, offset=sf%offset, &
                    current=sf%current, boozer_stream=sf%boozer_stream)
                sl%eps = 0; sl%k = 0; sl%offset = 0; sl%current = 0
                sl%boozer_stream = .not. sf%boozer_stream
                call ntv_kernel(sf, v, nu, omega, pre, 33, 80, kf, status)
                if (status /= NTV_OK) error stop 'kernel failed'
                call ntv_kernel(sl, v, nu, omega, pre, 33, 80, kl, status)
                if (status /= NTV_OK) error stop 'line kernel failed'
                print '(a,2es24.16)', 'surface/line R ', kf%r, kl%r
                ok = ok .and. kf%r == kl%r .and. kf%dr_dnu == kl%dr_dnu &
                    .and. kf%dr_domega == kl%dr_domega .and. all(kf%dr_dc == kl%dr_dc) &
                    .and. kf%indicator_theta == kl%indicator_theta &
                    .and. kf%indicator_leg == kl%indicator_leg
                call ntv_kernel(sl, v, nu, omega, pre, 33, 80, kl, status, direction=sf)
                ok = ok .and. status == NTV_BAD_INPUT
            end do
        end do
    end function line

    logical function sparse_refinement() result(ok)
        type(ntv_kernel_t) :: reference, local
        real(dp) :: errors(3)
        integer, parameter :: points(3) = [33,65,129]
        integer :: status, j
        call ntv_kernel(sample(),1.0_dp,0.02_dp,0.005_dp,0.02_dp,129,12, &
            reference,status,indicators=.false.)
        if (status /= NTV_OK) error stop 'reference failed'
        do j = 1,3
            call ntv_kernel(sample(),1.0_dp,0.02_dp,0.005_dp,0.02_dp,points(j),12, &
                local,status,indicators=.false.,theta_order=6)
            if (status /= NTV_OK) error stop 'local refinement failed'
            errors(j) = abs(local%r-reference%r)/abs(reference%r)
        end do
        print '(a,3es15.6)', 'relative theta errors ',errors
        ok = errors(2) < errors(1)/16 .and. errors(3) < errors(2)/16 &
            .and. errors(3) < 0.01_dp
    end function sparse_refinement

    real(dp) function sample_k()
        type(ntv_surface_t) :: sf
        sf = sample()
        sample_k = sf%k
    end function sample_k

    !> Dense oracle on a small grid for both weights and both streaming
    !> models: the equation
    !>   k_s xi (d_theta - i n q) F - (1 - xi^2)/2 k_s d_theta(ln b) d_xi F
    !>   + i (n/v)(omega + pre (d (1 + w xi^2) + c_J xi^2)) F
    !>   - (nu/2) d_xi((1 - xi^2) d_xi F)
    !>   = (1 + w xi^2) S/(b v),
    !> k_s = k or k b, is projected on P_i by Gauss quadrature (not the
    !> recurrences of the kernel), d_theta is the DFT spectral derivative,
    !> and R = (2 pi/n_theta) sum m_s Re(S^* int (1 + w xi^2) F), m_s = 2/b^2
    !> or 2/b^3.
    logical function assembly(theta_order) result(ok)
        integer, intent(in), optional :: theta_order
        integer :: order
        integer, parameter :: nt = 7, nl = 12, nq = 40, nn = nt*nl
        real(dp), parameter :: v = 1.3_dp, nu = 0.3_dp, omega = 0.2_dp, pre = 0.5_dp
        real(dp), parameter :: pi = 4*atan(1.0_dp)
        complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: kern
        real(dp) :: xq(nq), wq(nq), p(nq, nl), dp_(nq, nl), w, theta(nt), b(nt), d(nt)
        real(dp) :: ks(nt), dlnb(nt), meas(nt), wgt(nq), r
        complex(dp) :: a(nn, nn), rhs(nn), s(nt), dft(nt, nt), moment
        integer :: it, jt, i, j, status, info, piv(nn), iw, ib, mode, row, col
        interface
            subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
                integer, intent(out) :: ipiv(*), info
            end subroutine zgesv
        end interface
        order = 0
        if (present(theta_order)) order = theta_order
        call gauss_legendre(nq, xq, wq)
        p(:, 1) = 1; p(:, 2) = xq; dp_(:, 1) = 0; dp_(:, 2) = 1
        do j = 3, nl
            p(:, j) = ((2*j - 3)*xq*p(:, j - 1) - (j - 2)*p(:, j - 2))/(j - 1)
            dp_(:, j) = dp_(:, j - 2) + (2*j - 3)*p(:, j - 1)
        end do
        theta = [(2*pi*(it - 1)/nt, it=1, nt)]
        dft = 0
        do jt = 1, nt
            do it = 1, nt
                do mode = -(nt - 1)/2, (nt - 1)/2
                    dft(it, jt) = dft(it, jt) + iu*mode*exp(iu*mode*(theta(it) - theta(jt)))/nt
                end do
            end do
        end do
        if (order == 2) then
            dft = 0
            do it = 1, nt
                dft(it,modulo(it,nt)+1) = nt/(4*pi)
                dft(it,modulo(it-2,nt)+1) = -nt/(4*pi)
            end do
        end if
        ok = .true.
        do ib = 0, 1
            do iw = 1, 2
                sf = sample()
                sf%boozer_stream = ib == 1
                sf%weight = merge(NTV_WEIGHT_MU, NTV_WEIGHT_FULL, iw == 1)
                sf%current = merge(0.0_dp, 0.37_dp, iw == 1)
                w = merge(-1.0_dp, 1.0_dp, iw == 1)
                wgt = 1 + w*xq**2
                b = 1 - sf%eps*cos(theta)
                d = (cos(theta) - sf%offset)/b
                dlnb = sf%eps*sin(theta)/b
                ks = sf%k*merge(b, b**0, ib == 1)
                meas = 2/b**(2 + ib)
                s = iu*sf%n*sf%c(1)*exp(iu*sf%m(1)*theta)
                call ntv_kernel(sf, v, nu*v, omega, pre, nt, nl, kern, status, &
                    indicators=.false., &
                    theta_order=order)
                if (status /= NTV_OK) error stop 'kernel failed'
                a = 0
                do it = 1, nt
                    do i = 1, nl
                        row = (i - 1)*nt + it
                        rhs(row) = sum(wq*p(:, i)*wgt)*s(it)/(b(it)*v)
                        do j = 1, nl
                            col = (j - 1)*nt + it
                            a(row, col) = sum(wq*(nu/2*(1 - xq**2)*dp_(:, i)*dp_(:, j) &
                                + iu*sf%n/v*(omega + pre*(d(it)*wgt + sf%current*xq**2)) &
                                *p(:, i)*p(:, j) &
                                - ks(it)*iu*sf%n*sf%q*xq*p(:, i)*p(:, j) &
                                - ks(it)*dlnb(it)/2*(1 - xq**2)*p(:, i)*dp_(:, j)))
                            do jt = 1, nt
                                a(row, (j - 1)*nt + jt) = a(row, (j - 1)*nt + jt) &
                                    + ks(it)*dft(it, jt)*sum(wq*xq*p(:, i)*p(:, j))
                            end do
                        end do
                    end do
                end do
                call zgesv(nn, 1, a, nn, piv, rhs, nn, info)
                if (info /= 0) error stop 'oracle solve failed'
                r = 0
                do it = 1, nt
                    moment = 0
                    do j = 1, nl
                        moment = moment + rhs((j - 1)*nt + it)*sum(wq*wgt*p(:, j))
                    end do
                    r = r + meas(it)*real(conjg(s(it))*moment, dp)
                end do
                r = 2*pi/nt*r
                print '(a,i2,a,l2,2es22.14,es9.1)', 'weight', sf%weight, ' boozer', &
                    sf%boozer_stream, kern%r, r, abs(kern%r/r - 1)
                ok = ok .and. abs(kern%r/r - 1) < 1e-10_dp
            end do
        end do
    end function assembly

    !> FILE: eps k q offset n nm [weight boozer_stream(0|1) [c_J]] /
    !> nm lines m Re(c) Im(c) / n_theta n_leg [indicators 0|1 [theta order]] /
    !> nx / nx lines v nu_d omega_E pre. Prints per node R, the theta and
    !> Legendre indicators, dR/dnu_D and dR/domega_E.
    subroutine kernel_table(file)
        character(*), intent(in) :: file
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: k
        integer :: u, nm, j, n_theta, n_leg, nx, status, ios, with_indicators, weight, stream, &
            order
        real(dp) :: re, im, v, nu_d, omega, pre
        character(256) :: line, spectrum_file
        integer :: us
        ! Optional third argument: file for the Legendre spectrum per node.
        call get_command_argument(3, spectrum_file)
        if (len_trim(spectrum_file) > 0) open (newunit=us, file=trim(spectrum_file), &
            status='replace', action='write')
        open (newunit=u, file=file, status='old', action='read')
        read (u, '(a)') line
        read (line, *, iostat=ios) sf%eps, sf%k, sf%q, sf%offset, sf%n, nm, weight, stream, &
            sf%current
        if (ios /= 0) sf%current = 0
        if (ios /= 0) read (line, *, iostat=ios) sf%eps, sf%k, sf%q, sf%offset, sf%n, nm, &
            weight, stream
        if (ios /= 0) then
            read (line, *) sf%eps, sf%k, sf%q, sf%offset, sf%n, nm
        else
            sf%weight = weight
            sf%boozer_stream = stream /= 0
        end if
        allocate (sf%m(nm), sf%c(nm))
        do j = 1, nm
            read (u, *) sf%m(j), re, im
            sf%c(j) = cmplx(re, im, dp)
        end do
        read (u, '(a)') line
        read (line, *, iostat=ios) n_theta, n_leg, with_indicators, order
        if (ios /= 0) then
            order = 0
            read (line, *, iostat=ios) n_theta, n_leg, with_indicators
        end if
        if (ios /= 0) then
            read (line, *) n_theta, n_leg
            with_indicators = 1
        end if
        read (u, *) nx
        do j = 1, nx
            read (u, *) v, nu_d, omega, pre
            call ntv_kernel(sf, v, nu_d, omega, pre, n_theta, n_leg, k, status, &
                indicators=with_indicators /= 0, theta_order=order)
            if (status /= NTV_OK) error stop 'kernel failed'
            print '(5es24.15)', k%r, k%indicator_theta, k%indicator_leg, k%dr_dnu, k%dr_domega
            if (len_trim(spectrum_file) > 0) write (us, '(*(es12.4))') k%spectrum
        end do
        close (u)
        if (len_trim(spectrum_file) > 0) close (us)
    end subroutine kernel_table

    !> The f_10 provider against an independent n = 0 oracle: the equation
    !>   k b xi d_theta f - (1 - xi^2)/2 k b d_theta(ln b) d_xi f
    !>   - (nu/2) d_xi((1 - xi^2) d_xi f) = (1 + xi^2) d_theta b/(b v q),
    !> projected on P_i by Gauss quadrature with the DFT derivative, holds for
    !> the provided f_10 up to the density source, a theta-constant residual
    !> in the P_0 rows; the provided f_10 then drives the kernel.
    logical function f10_provider() result(ok)
        integer, parameter :: nt = 7, nl = 12, nq = 40, nn = nt*nl
        real(dp), parameter :: v = 1.3_dp, nu = 0.3_dp, pi = 4*atan(1.0_dp)
        type(ntv_surface_t) :: sf
        type(ntv_kernel_t) :: k0, k1
        real(dp), allocatable :: f10(:, :)
        real(dp) :: xq(nq), wq(nq), p(nq, nl), dp_(nq, nl), theta(nt), b(nt), &
            ks(nt), dlnb(nt), dft(nt, nt), a(nn, nn), rhs(nn), x(nn), r(nn), scale
        integer :: it, jt, i, j, mode, row, status

        sf = sample()
        sf%boozer_stream = .true.
        sf%weight = NTV_WEIGHT_FULL
        call ntv_f10(sf, v, nu*v, nt, nl, f10, status)
        if (status /= NTV_OK) error stop 'f10 provider failed'
        call gauss_legendre(nq, xq, wq)
        p(:, 1) = 1; p(:, 2) = xq; dp_(:, 1) = 0; dp_(:, 2) = 1
        do j = 3, nl
            p(:, j) = ((2*j - 3)*xq*p(:, j - 1) - (j - 2)*p(:, j - 2))/(j - 1)
            dp_(:, j) = dp_(:, j - 2) + (2*j - 3)*p(:, j - 1)
        end do
        theta = [(2*pi*(it - 1)/nt, it=1, nt)]
        dft = 0
        do jt = 1, nt
            do it = 1, nt
                do mode = 1, (nt - 1)/2
                    dft(it, jt) = dft(it, jt) - 2*mode*sin(mode*(theta(it) - theta(jt)))/nt
                end do
            end do
        end do
        b = 1 - sf%eps*cos(theta)
        dlnb = sf%eps*sin(theta)/b
        ks = sf%k*b
        a = 0
        do it = 1, nt
            do i = 1, nl
                row = (i - 1)*nt + it
                rhs(row) = sum(wq*p(:, i)*(1 + xq**2))*sf%eps*sin(theta(it)) &
                    /(b(it)*v*sf%q)
                do j = 1, nl
                    a(row, (j - 1)*nt + it) = sum(wq*(nu/2*(1 - xq**2)*dp_(:, i)*dp_(:, j) &
                        - ks(it)*dlnb(it)/2*(1 - xq**2)*p(:, i)*dp_(:, j)))
                    do jt = 1, nt
                        a(row, (j - 1)*nt + jt) = a(row, (j - 1)*nt + jt) &
                            + ks(it)*dft(it, jt)*sum(wq*xq*p(:, i)*p(:, j))
                    end do
                end do
            end do
        end do
        do j = 1, nl
            do it = 1, nt
                x((j - 1)*nt + it) = f10(merge(0, nt, mod(j - 1, 2) == 0) + it, (j - 1)/2)
            end do
        end do
        r = matmul(a, x) - rhs
        scale = maxval(abs(rhs))
        print '(a,2es9.1)', 'f_10 oracle residual (P_0 spread, P_l>0): ', &
            (maxval(r(1:nt)) - minval(r(1:nt)))/scale, maxval(abs(r(nt + 1:)))/scale
        ok = maxval(r(1:nt)) - minval(r(1:nt)) < 1e-12_dp*scale .and. &
            maxval(abs(r(nt + 1:))) < 1e-12_dp*scale

        ! The provided f_10 drives the kernel through source_dlf10.
        call ntv_kernel(sf, v, nu*v, 0.2_dp, 0.5_dp, nt, nl, k0, status, indicators=.false.)
        if (status /= NTV_OK) error stop 'kernel failed'
        sf%source_dlf10 = .true.
        call ntv_kernel(sf, v, nu*v, 0.2_dp, 0.5_dp, nt, nl, k1, status, &
            indicators=.false., f10=f10)
        if (status /= NTV_OK) error stop 'kernel failed'
        print '(a,es9.1)', 'response change from the provided f_10: ', abs(k1%r/k0%r - 1)
        ok = ok .and. abs(k1%r/k0%r - 1) > 1e-6_dp
    end function f10_provider

    !> The parallel-derivative quasilinear drive (#192): the two pieces of
    !> the perturbed radial drift together form the i m/q drive, it
    !> vanishes on the resonant m = n q, and the amplitude sensitivities
    !> keep their adjoint form with the switch on.
    logical function ql_dpar() result(ok)
        type(ntv_surface_t) :: sf, up, down
        type(ntv_kernel_t) :: k0, kp, km
        integer :: status, part
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        real(dp) :: h, fd
        complex(dp) :: dc

        sf = sample()
        ! With i n + i (m - n q)/q = i m/q the switch scales the drive, and
        ! so the response, by m/(n q), while the torque still pairs with S.
        call ntv_kernel(sf, v, nu, omega, pre, 33, 48, k0, status, indicators=.false.)
        if (status /= NTV_OK) error stop 'kernel failed'
        up = sf
        up%source_dpar = .true.
        call ntv_kernel(up, v, nu, omega, pre, 33, 48, kp, status, indicators=.false.)
        if (status /= NTV_OK) error stop 'kernel failed'
        fd = real(sf%m(1), dp)/(sf%n*sf%q)
        print '(a,es9.1)', 'dpar on the i m/q identity: ', abs(kp%r/(fd*k0%r) - 1)
        ok = abs(kp%r/(fd*k0%r) - 1) < 1e-12_dp

        ! On the resonant m = n q the parallel derivative of the
        ! perturbation along the line vanishes.
        up = sample()
        up%q = 2.0_dp
        up%k = 1/(up%q*6.2_dp)
        call ntv_kernel(up, v, nu, omega, pre, 33, 48, k0, status, indicators=.false.)
        up%source_dpar = .true.
        call ntv_kernel(up, v, nu, omega, pre, 33, 48, kp, status, indicators=.false.)
        print '(a,es9.1)', 'dpar on the resonant m = n q: ', abs(kp%r - k0%r)
        ok = ok .and. abs(kp%r - k0%r) <= 1e-13_dp*abs(k0%r)

        ! Amplitude sensitivity with the switch on against differences.
        sf = sample()
        sf%source_dpar = .true.
        call ntv_kernel(sf, v, nu, omega, pre, 33, 48, k0, status, indicators=.false.)
        do part = 1, 2
            dc = merge((1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), part == 1)
            h = 1e-4_dp*abs(sf%c(1))
            up = sf
            up%c(1) = sf%c(1) + h*dc
            down = sf
            down%c(1) = sf%c(1) - h*dc
            call ntv_kernel(up, v, nu, omega, pre, 33, 48, kp, status, indicators=.false.)
            call ntv_kernel(down, v, nu, omega, pre, 33, 48, km, status, indicators=.false.)
            fd = (kp%r - km%r)/(2*h)
            ok = ok .and. near(merge(real(k0%dr_dc(1)), aimag(k0%dr_dc(1)), &
                part == 1), fd, 'dR/dc dpar')
        end do
    end function ql_dpar

    !> The perturbed-operator sources on the axisymmetric solution f_10
    !> (#192): a pitch- and theta-constant f_10 is annihilated by both
    !> foldings, and with a nontrivial f_10 the response and its amplitude
    !> sensitivities with both source switches on agree with differences.
    logical function ql_f10() result(ok)
        type(ntv_surface_t) :: sf, up, down
        type(ntv_kernel_t) :: k0, kp, km
        integer :: status, i, j, part
        integer, parameter :: nth = 33, nleg = 48, nb = nleg/2
        real(dp), parameter :: v = 1.0_dp, nu = 2e-3_dp, omega = 5e-3_dp, pre = 2e-2_dp
        real(dp) :: f10(2*nth, 0:nb - 1), h, fd, theta
        real(dp), parameter :: pi = 4*atan(1.0_dp)
        complex(dp) :: dc

        ! Lorentz and mirror foldings annihilate the constant.
        sf = sample()
        sf%source_dlf10 = .true.
        f10 = 0
        f10(1:nth, 0) = 1.0_dp
        call ntv_kernel(sf, v, nu, omega, pre, nth, nleg, k0, status, &
            indicators=.false., f10=f10)
        if (status /= NTV_OK) error stop 'kernel failed'
        up = sf
        up%source_dlf10 = .false.
        call ntv_kernel(up, v, nu, omega, pre, nth, nleg, kp, status, indicators=.false.)
        if (status /= NTV_OK) error stop 'kernel failed'
        print '(a,es9.1)', 'dlf10 annihilates the constant f_10: ', abs(k0%r - kp%r)
        ok = abs(k0%r - kp%r) <= 1e-14_dp*abs(kp%r)

        ! A smooth nontrivial f_10 changes the response; both switches on.
        f10 = 0
        do i = 1, nth
            theta = 2*pi*(i - 1)/nth
            f10(i, 0) = 1 + 0.5_dp*cos(theta)
            f10(i, 1) = 0.3_dp*cos(theta)
            f10(nth + i, 0) = 0.2_dp*sin(theta)
        end do
        sf%source_dpar = .true.
        call ntv_kernel(sf, v, nu, omega, pre, nth, nleg, k0, status, &
            indicators=.false., f10=f10)
        if (status /= NTV_OK) error stop 'kernel failed'
        up = sf
        up%source_dlf10 = .false.
        call ntv_kernel(up, v, nu, omega, pre, nth, nleg, kp, status, indicators=.false.)
        if (status /= NTV_OK) error stop 'kernel failed'
        print '(a,es9.1)', 'response change from the f_10 sources: ', abs(k0%r/kp%r - 1)
        ok = ok .and. abs(k0%r/kp%r - 1) > 1e-6_dp
        do part = 1, 2
            dc = merge((1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), part == 1)
            h = 1e-4_dp*abs(sf%c(1))
            up = sf; up%c(1) = sf%c(1) + h*dc
            down = sf; down%c(1) = sf%c(1) - h*dc
            call ntv_kernel(up, v, nu, omega, pre, nth, nleg, kp, status, &
                indicators=.false., f10=f10)
            call ntv_kernel(down, v, nu, omega, pre, nth, nleg, km, status, &
                indicators=.false., f10=f10)
            fd = (kp%r - km%r)/(2*h)
            ok = ok .and. near(merge(real(k0%dr_dc(1)), aimag(k0%dr_dc(1)), &
                part == 1), fd, 'dR/dc f10')
        end do
    end function ql_f10

end program test_dk_ntv
