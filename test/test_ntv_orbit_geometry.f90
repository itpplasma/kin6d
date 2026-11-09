!> Finite-orbit geometry of the orbit-class spectral elements
!! (milca_ntv_orbit_geometry) in its thin-orbit limit. On a synthetic Boozer
!! file with B = 2.3 (1 - 0.1 sqrt(s) cos(theta)), constant G, I = 0 and no
!! magnetic shear (with shear the orbit also precesses through -q' <theta
!! psi_dot>, which the local thin field lacks) the
!! orbit has theta_dot = v_par iota B/G + H_psi and zeta_dot = v_par B/G, so in
!! the frame exp(-i n zeta) it precesses at n q H_psi with
!!   H_psi = (m v^2/(2 Z e B)) (1 + xi^2) dB/dpsi.
!! The thin field (hp_field_t) with Boozer streaming, the full drift weight
!! and the precession prefactor
!!   pre = q m v^2 dB/ds(cos)/(2 Z e B_ref psi'),
!! (offset 0) solves the same problem to first order in the orbit width, so
!! the NTV kernels of the two geometries agree when the orbit is thin.
program test_ntv_orbit_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_boozer_file, only: boozer_file_t, load_boozer_file, BOOZER_OK
    use milca_boozer_orbit, only: orbit_particle_t, orbit_state
    use milca_dk_hp, only: hp_field_t, hp_field, hp_mesh_t, hp_default_mesh, hp_operator_t, &
        hp_ntv_kernel, HP_OK
    use milca_ntv_orbit_geometry, only: hp_orbit_geometry_t, hp_orbit_geometry
    implicit none
    real(dp), parameter :: pi = 4*atan(1.0_dp), e = 1.602176634e-19_dp
    real(dp), parameter :: md = 3.3435837768e-27_dp, s_star = 0.45_dp, b_ref = 2.3_dp
    character(*), parameter :: file = 'ntv_orbit_geometry.bc'
    type(boozer_file_t), target :: bf
    type(orbit_particle_t) :: p
    type(hp_orbit_geometry_t) :: g
    type(hp_field_t) :: f
    type(hp_mesh_t) :: mesh_g, mesh_f
    type(hp_operator_t) :: op_g, op_f
    real(dp) :: b, b_s, gg, i, g_s, i_s, iota, psi_prime, v, r_g, r_f, samples(33), &
        omega_d, factor
    integer :: status, j
    integer(8) :: t0, t1, rate

    call write_file(file)
    call load_boozer_file(file, bf, status)
    if (status /= BOOZER_OK) error stop 'cannot read the test file'
    factor = 1e-3_dp
    if (command_argument_count() > 0) then
        block
            character(32) :: arg
            call get_command_argument(1, arg)
            read (arg, *) factor
        end block
    end if
    p%mass = md*factor
    p%charge = e
    v = sqrt(2*10e3_dp*e/md)
    p%energy = p%mass*v**2/2
    p%s_ref = s_star
    call hp_orbit_geometry(bf, p, s_star, b_ref, g, status)
    if (status /= 0) error stop 'orbit geometry failed'
    call orbit_state(bf, s_star, 0.0_dp, b, b_s, gg, i, g_s, i_s, iota, psi_prime, status)
    do j = 1, 33
        call orbit_state(bf, s_star, 2*pi*(j - 1)/33, b, b_s, gg, i, g_s, i_s, iota, &
            psi_prime, status)
        samples(j) = b/b_ref
    end do
    call hp_field(samples, iota*b_ref/gg, f, status)
    if (status /= HP_OK) error stop 'hp_field failed'

    ! dB/ds = -0.23/(2 sqrt(s)) cos(theta) for this file.
    op_f%pre = (1/iota)*p%mass*v**2*(-0.23_dp/(2*sqrt(s_star)))/(2*e*b_ref*psi_prime)
    omega_d = abs(op_f%pre)/v
    op_f%k = iota*b_ref/gg
    op_f%k_power = 1
    op_f%q = 1/iota
    op_f%n = 3
    op_f%m = [3]
    op_f%c = [(1e-3_dp, 0.0_dp)]
    op_f%v = v
    op_f%wsign = 1
    op_f%offset = 0
    op_f%nu = 0.1_dp*op_f%n*omega_d
    op_g = op_f
    op_g%pre = 0

    call hp_default_mesh(f, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_f)
    call hp_ntv_kernel(f, mesh_f, op_f, r_f, status)
    if (status /= HP_OK) error stop 'thin kernel failed'
    call system_clock(t0, rate)
    call hp_default_mesh(g, 4, 16, 6, 1e-3_dp, 0.3_dp, mesh_g)
    call hp_ntv_kernel(g, mesh_g, op_g, r_g, status)
    call system_clock(t1)
    if (status /= HP_OK) error stop 'orbit kernel failed'
    print '(a,2es20.11,a,es9.2,a,f8.2)', 'R thin field, orbit geometry ', r_f, r_g, ' rel ', &
        abs(r_g/r_f - 1), ' orbit s ', real(t1 - t0, dp)/rate
    if (abs(r_g/r_f - 1) > 1e-3_dp) error stop 1

contains

    subroutine write_file(name)
        character(*), intent(in) :: name
        integer :: unit, k
        real(dp) :: s
        open (newunit=unit, file=name, status='replace', action='write')
        write (unit, '(a)') 'CC synthetic test file'
        write (unit, '(a)') ' m0b n0b nsurf nper flux/[Tm^2] a/[m] R/[m]'
        write (unit, '(a)') ' 1 0 19 1 -1.3 0.5 1.7'
        do k = 1, 19
            s = 0.05_dp*k
            write (unit, '(a)') ' s iota curr_pol/nper curr_tor pprime sqrt g(0,0)'
            write (unit, '(a)') ' [A] [A] dp/ds,[Pa] (dV/ds)/nper'
            write (unit, '(6es24.16)') s, 0.765_dp, -1.2e7_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') ' m n r/[m] z/[m] (phib-phi)*nper/twopi bmn/[T]'
            write (unit, '(2i5,4es24.16)') 0, 0, 1.7_dp, 0.0_dp, 0.0_dp, 2.3_dp
            write (unit, '(2i5,4es24.16)') 1, 0, 0.5_dp*sqrt(s), 0.5_dp*sqrt(s), 0.0_dp, &
                -0.23_dp*sqrt(s)
        end do
        close (unit)
    end subroutine write_file
end program test_ntv_orbit_geometry
