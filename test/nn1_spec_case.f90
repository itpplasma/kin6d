! NN1 (vacuum, prescribed traces) on the common KIN6D/SPECTRE case: a circular annulus
! a_in < rho < a_out about (R0, 0); the traces and the harmonic F come from N helical closed
! polygonal filaments outside the domain (direct Biot-Savart is the shared oracle).
! Usage: nn1_spec_case case.nml grid.txt out_prefix nr nmodes [nphi]
!   grid.txt: lines "R Z phi" (m, m, rad).  Writes out_prefix.B.txt: R Z phi BR BZ BPHI | direct BR BZ BPHI
module nn1_spec_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d, only: coil_set_t, coil_set_from_polylines, coil_set_b_cart, &
        coil_set_modes, coil_mu0
    use kin6d_nonnested_leaves, only: segment_potential, vacuum_gauge
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    type(coil_set_t) :: coils
    real(dp) :: r0 = 3, ftor = 10, rc = 1.6_dp, current = 2.0e5_dp, a_in = 0.3_dp, a_out = 0.9_dp
    integer :: nfil = 2, kwind = 1, nseg = 256, nfp = 2, nphi = 64, active_mode = 1
contains
    subroutine read_case(file)
        character(*), intent(in) :: file
        integer :: u
        namelist /case/ r0, ftor, rc, current, a_in, a_out, nfil, kwind, nseg, nfp
        open (newunit=u, file=file, status='old', action='read')
        read (u, nml=case)
        close (u)
    end subroutine
    subroutine make_coil()
        real(dp), allocatable :: pts(:, :)
        integer, allocatable :: counts(:)
        real(dp) :: phi, th, r
        integer :: j, i, st
        allocate (pts(3, nfil*nseg), counts(nfil))
        counts = nseg
        do j = 1, nfil
            do i = 1, nseg
                phi = 2*pi*(i - 1)/nseg
                th = kwind*phi + 2*pi*(j - 1)/nfil
                r = r0 + rc*cos(th)
                pts(:, (j - 1)*nseg + i) = [r*cos(phi), r*sin(phi), rc*sin(th)]
            end do
        end do
        call coil_set_from_polylines(coils, counts, pts, status=st)
        if (st /= 0) error stop 'coil'
    end subroutine
    subroutine circle(u, theta, x)
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [r0 + u*cos(theta), u*sin(theta)]
    end subroutine
    subroutine gauge_potential(x, a)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: a(3)
        real(dp) :: rz(2, 1), ar, az, br, bz
        complex(dp) :: bn(3, 1, 1)
        real(dp), allocatable :: w(:)
        integer :: st, part
        rz(:, 1) = x
        allocate (w(nfil)); w = current
        call coil_set_modes(coils, w, rz, nphi, [active_mode], bn, st)
        if (st /= 0) error stop 'coil Fourier'
        a = 0
        do part = 1, 2
            if (part == 1) then
                br = real(bn(1, 1, 1), dp); bz = real(bn(3, 1, 1), dp)
            else
                br = aimag(bn(1, 1, 1)); bz = aimag(bn(3, 1, 1))
            end if
            call vacuum_gauge(x(1), br, bz, real(active_mode, dp), ar, az)
            if (part == 1) then
                a = cmplx(0.0_dp, [ar, 0.0_dp, az], dp)
            else
                a = a - [ar, 0.0_dp, az]
            end if
        end do
    end subroutine
    real(dp) function mean_psi(x) result(psi)
        real(dp), intent(in) :: x(2)
        real(dp) :: phi, p(3), a(3), coefficient, len, r1, r2
        integer :: j, s
        psi = 0
        do j = 1, nphi
            phi = 2*pi*(j - 1)/nphi
            p = [x(1)*cos(phi), x(1)*sin(phi), x(2)]
            a = 0
            do s = 1, coils%nseg
                len = norm2(coils%l(:, s))
                r1 = norm2(p - coils%a(:, s))
                r2 = norm2(p - coils%a(:, s) - coils%l(:, s))
                call segment_potential(r1, r2, len, coefficient)
                a = a + coefficient*coils%l(:, s)/len
            end do
            psi = psi + x(1)*(-a(1)*sin(phi) + a(2)*cos(phi))/nphi
        end do
        psi = psi*coil_mu0*current/(4*pi)
    end function
    subroutine direct_field(r, z, phi, b)
        real(dp), intent(in) :: r, z, phi
        real(dp), intent(out) :: b(3) ! (R, Z, phi), total with F/R
        real(dp) :: cart(3), p(3)
        real(dp), allocatable :: w(:)
        allocate (w(nfil)); w = current
        p = [r*cos(phi), r*sin(phi), z]
        call coil_set_b_cart(coils, w, p, cart)
        b = [cart(1)*cos(phi) + cart(2)*sin(phi), cart(3), &
            -cart(1)*sin(phi) + cart(2)*cos(phi) + ftor/r]
    end subroutine
end module

program nn1_spec_case
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn1_spec_fixture
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_feec_n, only: meridian_lagrange_nodes
    use kin6d_equilibrium_nonnested, only: free_b_field_t, nn_report_t, nn_build, &
        nn_solve, nn_vacuum_mode, nn_point
    use kin6d_equilibrium_state, only: EQ_VALUES, EQ_BR, EQ_BZ, EQ_BPHI
    implicit none
    type(flux_mesh_t) :: fm
    type(mesh_2d_t) :: topo
    type(nn_report_t) :: report
    type(free_b_field_t) :: field
    real(dp), allocatable :: rings(:), nodes(:, :), grid(:, :)
    complex(dp), allocatable :: bc(:), seed(:)
    real(dp) :: val(EQ_VALUES), b(3), t0, t1, worst_div, feff, fcheck, pr, pz
    character(256) :: arg(6)
    integer :: nr, nmodes, fan, nt, i, n, st, u, np, ios
    do i = 1, 5
        call get_command_argument(i, arg(i))
    end do
    read (arg(4), *) nr
    read (arg(5), *) nmodes
    call get_command_argument(6, arg(6))
    call read_case(trim(arg(1)))
    if (len_trim(arg(6)) > 0) read (arg(6), *) nphi
    call make_coil()
    ! The filaments carry a net link current through every toroidal loop in the domain, so
    ! their n = 0 B_phi is a harmonic mu0 I_link/(2 pi R): F_eff = F + R <B_phi,filaments>_phi.
    feff = 0; fcheck = 0
    do i = 1, 64
        call direct_field(r0, 0.0_dp, 2*pi*(i - 0.5_dp)/64, b)
        feff = feff + (b(3) - ftor/r0)*r0/64
        pr = r0 + 0.5_dp*a_out; pz = 0.3_dp*a_out
        call direct_field(pr, pz, 2*pi*(i - 0.5_dp)/64, b)
        fcheck = fcheck + (b(3) - ftor/pr)*pr/64
    end do
    print '(a,2es20.12)', 'F_eff - F at two points (must agree) ', feff, fcheck
    feff = ftor + feff
    open (newunit=u, file=trim(arg(2)), status='old', action='read')
    np = 0
    do
        read (u, *, iostat=ios) b
        if (ios /= 0) exit
        np = np + 1
    end do
    rewind (u)
    allocate (grid(3, np))
    do i = 1, np
        read (u, *) grid(:, i)
    end do
    close (u)
    call cpu_time(t0)
    allocate (rings(nr + 1))
    do i = 1, nr + 1
        rings(i) = a_in + (a_out - a_in)*(i - 1)/nr
    end do
    call flux_mesh_build(rings, 8*nr, 8*nr, 2, circle, fm, st)
    if (st /= 0) error stop 'mesh'
    fan = fm%tlast(1)
    nt = fm%mm%topology%n_triangles - fan
    topo = mesh_2d_t()
    topo%n_triangles = nt; topo%n_vertices = fm%mm%topology%n_vertices - 1
    topo%has_triangles = .true.
    topo%vertices = fm%mm%topology%vertices(:, 2:)
    topo%triangles = fm%mm%topology%triangles(:, fan + 1:) - 1
    fm%mm%topology = topo
    fm%mm%nodes = fm%mm%nodes(:, :, fan + 1:)
    fm%tlast = fm%tlast - fan; fm%vfirst = fm%vfirst - 1
    fm%ring_tri = fm%ring_tri - fan
    call nn_build(fm, 2, nmodes, [r0, 0.0_dp], 1.0_dp, feff, coil_mu0, field, st)
    if (st /= 0) error stop 'FEEC build'
    call meridian_lagrange_nodes(field%fs, nodes)
    allocate (bc(field%fs%n0), seed(field%fs%n0))
    bc = 0; seed = 0
    !$omp parallel do default(shared) private(i) schedule(static)
    do i = 1, size(bc)
        if (field%fs%sp(0)%bnd(i)) bc(i) = mean_psi(nodes(:, i))
    end do
    !$omp end parallel do
    call nn_solve(field, 0, bc, seed, report, st)
    if (st /= 0) error stop 'n=0 solve'
    print '(a,es12.4)', 'n=0 divergence ', report%divergence
    worst_div = report%divergence
    do n = 1, nmodes
        if (mod(n, nfp) /= 0) cycle
        active_mode = n
        call nn_vacuum_mode(field, n, gauge_potential, report, st)
        if (st /= 0) error stop 'mode solve'
        print '(a,i0,a,3es12.4)', 'n=', n, ' div, residual, energy_after ', &
            report%divergence, report%residual, report%energy_after
        worst_div = max(worst_div, report%divergence)
    end do
    call cpu_time(t1)
    print '(a,es12.4,a,f10.2)', 'worst divergence ', worst_div, '  solve cpu s ', t1 - t0
    open (newunit=u, file=trim(arg(3))//'.B.txt', status='replace', action='write')
    write (u, '(a)') '# R Z phi BR BZ BPHI | direct BR BZ BPHI (total B, T)'
    do i = 1, np
        call nn_point(field, grid(:, i), val, st)
        if (st /= 0) then
            print '(a,3es14.6)', 'point outside mesh ', grid(:, i)
            error stop 'sampling'
        end if
        call direct_field(grid(1, i), grid(2, i), grid(3, i), b)
        write (u, '(9es24.16)') grid(:, i), val(EQ_BR), val(EQ_BZ), val(EQ_BPHI), b
    end do
    close (u)
    print '(a,i0,a)', 'wrote ', np, ' points'
end program
