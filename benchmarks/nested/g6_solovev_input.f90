! Matched fixed-boundary VMEC input for the NN2 Solov'ev regression.
! This is independent oracle/adapter code, not a production closure fit.
program g6_solovev_input
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nn_solovev_fixture, only: r0, minor, pp, c, mapping
    use g6_solovev_oracle, only: exact_flux, exact_iota
    implicit none
    integer, parameter :: np = 65, nm = 16, nt = 512
    real(dp), parameter :: pi = acos(-1.0_dp), mu0 = 4*pi*1.0e-7_dp
    real(dp) :: rbc(0:0, 0:nm), zbs(0:0, 0:nm), raxis_cc(0:0)
    real(dp) :: am_aux_s(np), am_aux_f(np), ai_aux_s(np), ai_aux_f(np)
    real(dp) :: phiedge, theta, t, x(2), reconstructed(2), defect
    real(dp) :: ftol_array(3) = 1.0e-12_dp
    integer :: ns_array(3) = [33, 65, 129], niter_array(3) = 10000
    integer :: nfp = 1, mpol = nm + 1, ntor = 0, ncurr = 0
    logical :: lasym = .false., lfreeb = .false.
    character(32) :: pmass_type = 'cubic_spline', piota_type = 'cubic_spline'
    character(1024) :: path
    integer :: i, m, unit
    namelist /indata/ nfp, lasym, lfreeb, mpol, ntor, ncurr, ns_array, &
        niter_array, ftol_array, phiedge, pmass_type, piota_type, &
        am_aux_s, am_aux_f, ai_aux_s, ai_aux_f, rbc, zbs, raxis_cc
    call get_command_argument(1, path)
    if (len_trim(path) == 0) error stop 'usage: benchmark_g6_solovev_input OUTPUT'
    rbc = 0; zbs = 0
    do i = 1, nt
        theta = 2*pi*(i - 1)/nt
        call mapping(1.0_dp, theta, x)
        do m = 0, nm
            rbc(0, m) = rbc(0, m) + x(1)*cos(m*theta)/nt
            zbs(0, m) = zbs(0, m) + x(2)*sin(m*theta)/nt
        end do
    end do
    rbc(0, 1:) = 2*rbc(0, 1:); zbs(0, 1:) = 2*zbs(0, 1:)
    defect = 0
    do i = 1, nt
        theta = 2*pi*(i - 0.37_dp)/nt
        call mapping(1.0_dp, theta, x)
        reconstructed = 0
        do m = 0, nm
            reconstructed = reconstructed &
                + [rbc(0, m)*cos(m*theta), zbs(0, m)*sin(m*theta)]
        end do
        defect = max(defect, norm2(reconstructed - x)/minor)
    end do
    if (defect > 1.0e-10_dp) error stop 'Solovev boundary truncation'
    phiedge = exact_flux(1.0_dp)
    raxis_cc = r0
    do i = 1, np
        t = real(i - 1, dp)/(np - 1)
        am_aux_s(i) = exact_flux(sqrt(t))/phiedge
        am_aux_f(i) = pp*4*c*r0*r0*minor*minor*(t - 1)/mu0
        ai_aux_f(i) = exact_iota(sqrt(t))
    end do
    ai_aux_s = am_aux_s
    open (newunit=unit, file=trim(path), status='replace', action='write')
    write (unit, nml=indata)
    close (unit)
    print '(a,3es16.7)', '[V] input flux/axis iota/boundary defect ', &
        phiedge, ai_aux_f(1), defect
end program
