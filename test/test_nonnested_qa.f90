module nn_qa_oracle
    ! Landreman-Paul precise-QA vacuum benchmark for mapped NN1
    ! (docs/nn1-qa-benchmark.md). Oracle: Biot-Savart field and vector potential
    ! of the published 16-coil set (polygonal filaments, exact per segment).
    ! Wall: the published VMEC boundary R, Z = sum RBC cos / ZBS sin (m t - nfp n p)
    ! as the image of a reference disk, with the polynomial (harmonic) extension
    !   (R, Z) = sum RBC Re(w) - b (1 - |z|^2) m ZBS Re(w), ... Im(w) with
    !   w = z^m exp(-i nfp n phi), z = ((R - 1) + i Z) / a_ref,
    ! where the b-term vanishes on the boundary and removes the fold of the
    ! plain harmonic extension at the inboard dent (det DF > 0 for b = 0.1).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: qa_load, qa_map, qa_potential, qa_field, qa_flux, qa_bs
    public :: qa_circle, qa_nfp, qa_aref, qa_mapped_modes, qa_surface_residual
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4.0e-7_dp*pi
    real(dp), parameter :: qa_aref = 0.2_dp, qa_fold = 0.1_dp
    integer :: qa_nfp = 0, qa_mapped_modes = 0
    integer :: ncoil = 0, ntab = 0
    integer, allocatable :: tn(:), tm(:)
    real(dp), allocatable :: rbc(:), zbs(:)
    real(dp), allocatable :: p0(:, :), dl(:, :), seglen(:), cur(:)
contains
    subroutine qa_load(dir)
        character(*), intent(in) :: dir
        integer :: u, i, c, k, npt, nseg
        real(dp), allocatable :: pts(:, :)
        open (newunit=u, file=trim(dir)//'/qa_boundary.dat', status='old', action='read')
        read (u, *) qa_nfp, ntab
        allocate (tn(ntab), tm(ntab), rbc(ntab), zbs(ntab))
        do i = 1, ntab
            read (u, *) tn(i), tm(i), rbc(i), zbs(i)
        end do
        close (u)
        qa_mapped_modes = qa_nfp*maxval(abs(tn))
        open (newunit=u, file=trim(dir)//'/qa_coils.dat', status='old', action='read')
        read (u, *) ncoil, npt
        nseg = ncoil*npt
        allocate (p0(3, nseg), dl(3, nseg), seglen(nseg), cur(nseg), pts(3, npt))
        do c = 1, ncoil
            read (u, *) cur((c - 1)*npt + 1)
            cur((c - 1)*npt + 1:c*npt) = cur((c - 1)*npt + 1)
            do k = 1, npt
                read (u, *) pts(:, k)
            end do
            do k = 1, npt
                p0(:, (c - 1)*npt + k) = pts(:, k)
                dl(:, (c - 1)*npt + k) = pts(:, modulo(k, npt) + 1) - pts(:, k)
            end do
        end do
        close (u)
        seglen = sqrt(sum(dl**2, 1))
    end subroutine

    pure subroutine qa_bs(x, a, b)
        ! Cartesian vector potential a and field b of the coil set at x [T, T m].
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3), b(3)
        real(dp) :: r1(3), r2(3), n1, n2, f, cr(3)
        integer :: s
        a = 0; b = 0
        do s = 1, size(seglen)
            r1 = x - p0(:, s)
            r2 = r1 - dl(:, s)
            n1 = sqrt(sum(r1**2)); n2 = sqrt(sum(r2**2))
            f = cur(s)*1.0e-7_dp
            a = a + f*dl(:, s)/seglen(s)*log((n1 + n2 + seglen(s))/(n1 + n2 - seglen(s)))
            cr = [dl(2, s)*r1(3) - dl(3, s)*r1(2), dl(3, s)*r1(1) - dl(1, s)*r1(3), &
                  dl(1, s)*r1(2) - dl(2, s)*r1(1)]
            b = b + f*(n1 + n2)/(n1*n2*(n1*n2 + sum(r1*r2)))*cr
        end do
    end subroutine

    subroutine qa_potential(x, v)
        ! Physical (R, phi, Z) components of A at x = (R, Z, phi).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: a(3), b(3), c, s
        call qa_bs([x(1)*cos(x(3)), x(1)*sin(x(3)), x(2)], a, b)
        c = cos(x(3)); s = sin(x(3))
        v = [a(1)*c + a(2)*s, -a(1)*s + a(2)*c, a(3)]
    end subroutine

    subroutine qa_field(x, v)
        ! Physical (R, phi, Z) components of B at x = (R, Z, phi).
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: v(3)
        real(dp) :: a(3), b(3), c, s
        call qa_bs([x(1)*cos(x(3)), x(1)*sin(x(3)), x(2)], a, b)
        c = cos(x(3)); s = sin(x(3))
        v = [b(1)*c + b(2)*s, -b(1)*s + b(2)*c, b(3)]
    end subroutine

    function qa_flux(r, z) result(f)
        ! F = (1/2 pi) int R B_phi dphi on the circle (R, Z).
        real(dp), intent(in) :: r, z
        real(dp) :: f, v(3)
        integer :: j
        integer, parameter :: n = 128
        f = 0
        do j = 1, n
            call qa_field([r, z, 2*pi*(j - 1)/n], v)
            f = f + r*v(2)/n
        end do
    end function

    subroutine qa_circle(u, theta, x)
        ! Reference disk: centre (1, 0), radius qa_aref.
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [1 + qa_aref*u*cos(theta), qa_aref*u*sin(theta)]
    end subroutine

    subroutine qa_map(x, phi, y)
        real(dp), intent(in) :: x(2), phi
        real(dp), intent(out) :: y(2)
        complex(dp) :: z, w
        real(dp) :: damp
        integer :: i
        z = cmplx(x(1) - 1, x(2), dp)/qa_aref
        damp = 1 - abs(z)**2
        y = 0
        do i = 1, ntab
            w = z**tm(i)*exp(cmplx(0.0_dp, -real(qa_nfp*tn(i), dp)*phi, dp))
            y(1) = y(1) + (rbc(i) - qa_fold*damp*tm(i)*zbs(i))*real(w, dp)
            y(2) = y(2) + (zbs(i) - qa_fold*damp*tm(i)*rbc(i))*aimag(w)
        end do
    end subroutine

    subroutine qa_surface_residual(mean_rel, max_rel, bmean)
        ! B.n / |B| of the coil field on the published surface (area weighted).
        real(dp), intent(out) :: mean_rel, max_rel, bmean
        integer, parameter :: nt = 96, np = 48
        real(dp) :: th, ph, r, dr_t, dz_t, dr_p, dz_p, a, tt(3), tp(3), nv(3), bb(3), aa(3)
        real(dp) :: x(3), area, sa, sbn, sb
        integer :: i, j, k
        sa = 0; sbn = 0; sb = 0; max_rel = 0
        do j = 1, np
            ph = 2*pi*(j - 0.5_dp)/(qa_nfp*np)
            do i = 1, nt
                th = 2*pi*(i - 0.5_dp)/nt
                r = 0; dr_t = 0; dz_t = 0; dr_p = 0; dz_p = 0; x(3) = 0
                do k = 1, ntab
                    a = tm(k)*th - qa_nfp*tn(k)*ph
                    r = r + rbc(k)*cos(a)
                    x(3) = x(3) + zbs(k)*sin(a)
                    dr_t = dr_t - rbc(k)*tm(k)*sin(a)
                    dz_t = dz_t + zbs(k)*tm(k)*cos(a)
                    dr_p = dr_p + rbc(k)*qa_nfp*tn(k)*sin(a)
                    dz_p = dz_p - zbs(k)*qa_nfp*tn(k)*cos(a)
                end do
                tt = [dr_t*cos(ph), dr_t*sin(ph), dz_t]
                tp = [dr_p*cos(ph) - r*sin(ph), dr_p*sin(ph) + r*cos(ph), dz_p]
                nv = [tt(2)*tp(3) - tt(3)*tp(2), tt(3)*tp(1) - tt(1)*tp(3), &
                      tt(1)*tp(2) - tt(2)*tp(1)]
                area = sqrt(sum(nv**2)); nv = nv/area
                call qa_bs([r*cos(ph), r*sin(ph), x(3)], aa, bb)
                sa = sa + area; sb = sb + area*sqrt(sum(bb**2))
                sbn = sbn + area*abs(sum(bb*nv))
                max_rel = max(max_rel, abs(sum(bb*nv))/sqrt(sum(bb**2)))
            end do
        end do
        mean_rel = sbn/sb; bmean = sb/sa
    end subroutine
end module nn_qa_oracle

program test_nonnested_qa
    ! Usage: DATADIR check | case nr N [nphi] | hladder N nr... | modes nr N...
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_qa_oracle
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_build
    use kin6d_nonnested_mapped, only: nn_map_t, nn_map_build, nn_mapped_data, &
        nn_mapped_solve, nn_mapped_error, nn_mapped_report_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Pre-declared rule (docs/nn1-qa-benchmark.md, section Rule).
    real(dp), parameter :: gate_order = 1.2_dp, gate_final = 0.05_dp
    real(dp), parameter :: gate_mode_gain = 0.5_dp, gate_mode_sat = 0.2_dp
    character(256) :: dir, mode, arg
    integer, allocatable :: list(:)
    integer :: i, fixed_arg, nphi_arg
    real(dp), allocatable :: errs(:)
    real(dp) :: ord
    logical :: ok
    ok = .true.
    call get_command_argument(1, dir)
    call get_command_argument(2, mode)
    call qa_load(trim(dir))
    allocate (list(max(0, command_argument_count() - 3)))
    fixed_arg = 0; nphi_arg = 0
    if (command_argument_count() > 2) then
        call get_command_argument(3, arg)
        read (arg, *) fixed_arg
    end if
    do i = 1, size(list)
        call get_command_argument(i + 3, arg)
        read (arg, *) list(i)
    end do
    select case (trim(mode))
    case ('check')
        call check()
    case ('case')
        if (size(list) > 1) nphi_arg = list(2)
        allocate (errs(1))
        call qa_case(fixed_arg, list(1), nphi_arg, errs(1))
    case ('hladder')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call qa_case(list(i), fixed_arg, 0, errs(i))
            if (i > 1) then
                ord = log(errs(i - 1)/errs(i))/log(real(list(i), dp)/list(i - 1))
                write (*, '(a,i0,a,f6.3)') 'order nr=', list(i), ': ', ord
                if (ord < gate_order) then
                    ok = .false.
                    write (*, '(a,i0)') 'GATE order violated at nr=', list(i)
                end if
            end if
        end do
        if (errs(size(errs)) > gate_final) then
            ok = .false.
            write (*, '(a,es9.2)') 'GATE final error violated: ', errs(size(errs))
        end if
        if (.not. ok) error stop 'QA h ladder violates the declared rule'
    case ('modes')
        allocate (errs(size(list)))
        do i = 1, size(list)
            call qa_case(fixed_arg, list(i), 0, errs(i))
        end do
        if (errs(size(errs)) > gate_mode_gain*errs(1)) then
            ok = .false.; write (*, '(a)') 'GATE mode gain violated'
        end if
        if (size(errs) > 1) then
            if (abs(errs(size(errs)) - errs(size(errs) - 1)) &
                > gate_mode_sat*errs(size(errs) - 1)) then
                ok = .false.; write (*, '(a)') 'GATE mode saturation violated'
            end if
        end if
        if (.not. ok) error stop 'QA mode ladder violates the declared rule'
    case default
        error stop 'usage: DATADIR check | case nr N [nphi] | hladder N nr... | modes nr N...'
    end select
contains
    subroutine check()
        ! Oracle consistency: F constant, curl A = B (finite differences), vacuum
        ! residual of the published surface, field values for cross-checks.
        real(dp) :: x(3), a1(3), a2(3), b(3), curl(3), h, mr, xr, bm, v(3), f1, f2
        real(dp) :: pts(3, 3)
        integer :: k, c
        f1 = qa_flux(1.0_dp, 0.0_dp); f2 = qa_flux(1.08_dp, 0.12_dp)
        write (*, '(a,es20.12,a,es20.12)') 'F(1,0)=', f1, ' F(1.08,0.12)=', f2
        pts = reshape([1.05_dp, 0.02_dp, 0.3_dp, 0.97_dp, -0.05_dp, 1.9_dp, &
                       1.02_dp, 0.10_dp, 4.4_dp], [3, 3])
        h = 1.0e-5_dp
        do k = 1, 3
            call qa_field(pts(:, k), v)
            write (*, '(a,3f12.8,a,3es22.14)') 'B(R,Z,phi)=', pts(:, k), ' (R,phi,Z)=', v
            ! curl A in cylindrical components by central differences
            x = pts(:, k)
            curl = 0
            call diff_curl(x, h, curl)
            write (*, '(a,es10.3)') '  |curl A - B|/|B| = ', &
                sqrt(sum((curl - v)**2))/sqrt(sum(v**2))
        end do
        call qa_surface_residual(mr, xr, bm)
        write (*, '(a,es10.3,a,es10.3,a,es10.3)') 'surface: <|B.n|>/<|B|>=', mr, &
            ' max|B.n|/|B|=', xr, ' <|B|>=', bm
        c = qa_nfp
        write (*, '(a,i0,a,i0)') 'nfp=', c, ' mapped modes=', qa_mapped_modes
        if (abs(f1 - f2) > 1.0e-9_dp*abs(f1)) error stop 'flux not constant'
        if (mr > 1.0e-3_dp) error stop 'surface residual above 1e-3'
    end subroutine

    subroutine diff_curl(x, h, curl)
        real(dp), intent(in) :: x(3), h
        real(dp), intent(out) :: curl(3)
        real(dp) :: ap(3), am(3), dar(3), daz(3), daph(3), r
        r = x(1)
        call qa_potential(x + [h, 0.0_dp, 0.0_dp], ap); call qa_potential(x - [h, 0.0_dp, 0.0_dp], am)
        dar = (ap - am)/(2*h)
        call qa_potential(x + [0.0_dp, h, 0.0_dp], ap); call qa_potential(x - [0.0_dp, h, 0.0_dp], am)
        daz = (ap - am)/(2*h)
        call qa_potential(x + [0.0_dp, 0.0_dp, h], ap); call qa_potential(x - [0.0_dp, 0.0_dp, h], am)
        daph = (ap - am)/(2*h)
        call qa_potential(x, ap)
        ! components (R, phi, Z) of A; curl in cylindrical coordinates
        curl(1) = daph(3)/r - daz(2)
        curl(2) = daz(1) - dar(3)
        curl(3) = (ap(2) + r*dar(2) - daph(1))/r
    end subroutine

    subroutine qa_case(nr, modes, nphi_in, rel)
        integer, intent(in) :: nr, modes, nphi_in
        real(dp), intent(out) :: rel
        type(free_b_field_t), target :: field
        type(flux_mesh_t) :: fm
        type(nn_map_t), target :: map
        type(nn_mapped_report_t) :: report
        real(dp), allocatable :: data(:)
        real(dp) :: rings(nr), f, err, dev, eh, ee
        integer :: j, st, nphi, c0, c1, rate, tsetup
        nphi = merge(nphi_in, 4*modes + 4, nphi_in > 0)
        call system_clock(c0, rate)
        do j = 1, nr
            rings(j) = real(j, dp)/nr
        end do
        call flux_mesh_build(rings, 8, 4096, 2, qa_circle, fm, st)
        if (st /= 0) error stop 'mesh'
        f = qa_flux(1.0_dp, 0.0_dp)
        call nn_build(fm, 2, modes, [1.0_dp, 0.0_dp], qa_aref, f, 1.0_dp, field, st)
        if (st /= 0) error stop 'build'
        call nn_map_build(field%fs, qa_mapped_modes, 64, qa_map, map, st)
        if (st /= 0) error stop 'map'
        call nn_mapped_data(field, map, nphi, qa_potential, data, st)
        if (st /= 0) error stop 'data'
        call system_clock(c1); tsetup = c1 - c0
        call nn_mapped_solve(field, map, nphi, data, 1.0e-12_dp, report, st)
        if (st /= 0) error stop 'solve'
        call nn_mapped_error(field, map, 2*nphi + 1, qa_field, err, dev, eh, ee, st)
        if (st /= 0) error stop 'error (map not invertible?)'
        call system_clock(c1)
        rel = err/dev
        write (*, '(a,i0,a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,es10.3,a,i0,a,es9.2,a,es9.2,a,f8.1,a,f8.1)') &
            'QA nr=', nr, ' N=', modes, ' nphi=', nphi, ' ndof=', report%ndof, &
            ' err/|B-Bphi|=', rel, ' err/|B|=', err/sqrt(2*ee), ' energy rel=', &
            abs(report%energy - ee)/ee, ' fgmres=', report%iterations, &
            ' residual=', report%residual, ' div=', report%divergence, &
            ' s_setup=', real(tsetup, dp)/rate, ' s_total=', real(c1 - c0, dp)/rate
        if (report%divergence > 1.0e-12_dp) error stop 'divergence'
        if (report%residual > 1.0e-9_dp) error stop 'residual'
    end subroutine
end program test_nonnested_qa
