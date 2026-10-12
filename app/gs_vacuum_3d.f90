program gs_vacuum_3d
    !! L3: per-n vacuum field of 3D coils on the flux surfaces of a G-EQDSK
    !! equilibrium, projected onto PEST and Hamada harmonics b_mn.
    !! Usage: gs_vacuum_3d input.nml   (namelist /l3/, see below)
    !! Conventions are fixed in module kin6d_vacuum_surfaces and kin6d_coils_3d.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_coils_3d
    use kin6d_vacuum_surfaces
    implicit none
    integer, parameter :: maxs = 64, maxn = 16
    character(len=512) :: eqdsk = '', coils = '', currents = '', prefix = 'l3', nml
    integer :: nphi = 512, ntheta = 512, mmax = 48, nmode = 0, nrat = 0
    integer :: n_list(maxn) = 0, rat_n(maxn) = 0, rat_num(maxn) = 0, rat_den(maxn) = 0
    real(dp) :: psin(maxs) = 0, psin_rat_max = 0.997_dp, psin_rat_min = 0.2_dp
    namelist /l3/ eqdsk, coils, currents, prefix, nphi, ntheta, mmax, n_list, psin, &
        rat_n, rat_num, rat_den, psin_rat_max, psin_rat_min

    type(psi_grid_t) :: g
    type(coil_set_t) :: set
    type(surface_t), allocatable :: surf(:)
    real(dp), allocatable :: w(:), rz(:, :), psis(:)
    complex(dp), allocatable :: bn(:, :, :), bmn(:), phimn(:), bb(:, :)
    integer :: unit, uq, ub, ur, st, i, j, ns, ir, ia, in, m, off
    real(dp) :: qref, target, root, t0, t1
    type(surface_t) :: sref
    logical :: found
    character(len=8), parameter :: aname(2) = ['pest   ', 'hamada ']

    if (command_argument_count() /= 1) error stop 'usage: gs_vacuum_3d input.nml'
    call get_command_argument(1, nml)
    open (newunit=unit, file=trim(nml), status='old', action='read')
    read (unit, nml=l3)
    close (unit)
    nmode = count(n_list > 0)
    nrat = count(rat_n > 0)
    if (nmode == 0) error stop 'n_list is empty'

    call eqdsk_read(trim(eqdsk), g, st)
    if (st /= 0) error stop 'cannot read EQDSK'
    call coil_set_read(set, trim(coils), st)
    if (st /= 0) error stop 'cannot read coil table'
    allocate (w(set%ncoil))
    call coil_read_currents(trim(currents), set%ncoil, w, st)
    if (st /= 0) error stop 'cannot read currents'

    ! surfaces: explicit list followed by resonant surfaces
    allocate (psis(maxs + maxn))
    ns = 0
    do i = 1, maxs
        if (psin(i) > 0) then
            ns = ns + 1; psis(ns) = psin(i)
        end if
    end do
    call surface_build(g, 0.5_dp, ntheta, sref, st)
    qref = sref%q
    do ir = 1, nrat
        target = sign(real(rat_num(ir), dp)/rat_den(ir), qref)
found = find_surface_q(g, target, psin_rat_min, psin_rat_max, 256, 1.0e-13_dp, root, st)
        if (found) then
            ns = ns + 1; psis(ns) = root
      write (*, '(a,i0,a,i0,a,f12.9)') 'rational q = ', rat_num(ir), '/', rat_den(ir), &
                ' at psin = ', root
        else
            write (*, '(a,i0,a,i0,a)') 'rational q = ', rat_num(ir), '/', rat_den(ir), &
                ' not found in range'
        end if
    end do

    allocate (surf(ns), rz(2, ns*ntheta), bn(3, nmode, ns*ntheta))
    do i = 1, ns
        call surface_build(g, psis(i), ntheta, surf(i), st)
        if (st /= 0) error stop 'surface construction failed'
        do j = 1, ntheta
            rz(:, (i - 1)*ntheta + j) = [surf(i)%r(j), surf(i)%z(j)]
        end do
    end do
    call cpu_time(t0)
    call coil_set_modes(set, w, rz, nphi, pack(n_list, n_list > 0), bn, st)
    call cpu_time(t1)
    if (st /= 0) error stop 'coil projection failed'
    write (*, '(a,f10.2,a)') 'coil projection cpu time: ', t1 - t0, ' s'

    open (newunit=uq, file=trim(prefix)//'_q.dat', status='replace')
    write (uq, '(a)') '# psin qs (signed; qs<0: field lines run clockwise in (R,Z)'// &
        ' as phi increases)'
    open (newunit=ub, file=trim(prefix)//'_bmn.dat', status='replace')
    write (ub, '(a)') '# angle n psin qs m Re(b_mn) Im(b_mn) [T, peak amplitude] '// &
        'Re(phi_mn) Im(phi_mn) [T m^2 per rad^2, peak]'
    open (newunit=ur, file=trim(prefix)//'_real.dat', status='replace')
    write (ur, '(a)') '# psin n theta_g R Z nR nZ Re(2 B_n.n) Im(2 B_n.n) '// &
        'theta_pest theta_hamada nu'
    allocate (bmn(-mmax:mmax), phimn(-mmax:mmax), bb(2, ntheta))
    do i = 1, ns
        write (uq, '(2es24.15e3)') psis(i), surf(i)%q
        off = (i - 1)*ntheta
        do in = 1, nmode
            do j = 1, ntheta
                bb(1, j) = bn(1, in, off + j)
                bb(2, j) = bn(3, in, off + j)
       write (ur, '(es24.14e3,i4,es24.14e3,2f16.10,2es24.14e3,2es24.15e3,3es24.14e3)') &
                  psis(i), n_list(in), surf(i)%theta_g(j), surf(i)%r(j), surf(i)%z(j), &
                    surf(i)%nr(j), surf(i)%nz(j), &
                    2*(bb(1, j)*surf(i)%nr(j) + bb(2, j)*surf(i)%nz(j)), &
                    surf(i)%theta_s(j, 1), surf(i)%theta_s(j, 2), surf(i)%nu(j)
            end do
            do ia = 1, 2
            call surface_harmonics(surf(i), n_list(in), bb, ia, -mmax, mmax, bmn, phimn)
                do m = -mmax, mmax
                    write (ub, '(a,i4,es22.14e3,es24.15e3,i5,4es24.15e3)') aname(ia), &
                        n_list(in), psis(i), surf(i)%q, m, bmn(m), phimn(m)
                end do
            end do
        end do
    end do
    close (uq); close (ub); close (ur)
end program
