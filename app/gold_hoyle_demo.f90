program gold_hoyle_demo
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_grad_shafranov_leaves, only: gs_gold_hoyle
    use kin6d_result_io, only: run_t, result_file_t
    implicit none
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp) :: r0 = 6.2_dp, a = 0.62_dp, b0 = 5.3_dp, q = 1.5_dp, k
    real(dp), allocatable :: radius(:), bz(:), bt(:), psi(:)
    integer :: samples = 257, unit, status, i
    character(256) :: message
    namelist /gold_hoyle/ r0, a, b0, q, samples
    call run%start('gold_hoyle')
    status = 0; message = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=gold_hoyle, iostat=status, iomsg=message)
    call run%close_config(unit, 'gold_hoyle', status, message)
    if (r0 <= 0 .or. a <= 0 .or. q == 0 .or. samples < 2) &
        error stop 'invalid Gold-Hoyle cylindrical model'
    call run%open_resolved(unit)
    write (unit, nml=gold_hoyle)
    call run%close_resolved(unit)
    allocate (radius(samples), bz(samples), bt(samples), psi(samples))
    k = 1/(q*r0)
    do i = 1, samples
        radius(i) = a*real(i - 1, dp)/(samples - 1)
        call gs_gold_hoyle(radius(i), r0, b0, k, bz(i), bt(i), psi(i))
    end do
    call result%create(run)
    call result%attribute('model', 'straight cylindrical Gold-Hoyle; not finite-R GS')
    call result%attribute('claim', 'analytical model; sampled floating-point evaluation')
    call result%write('r', radius, 'm', 'distance from cylinder axis')
    call result%write('B_z', bz, 'T', 'axial magnetic field')
    call result%write('B_theta', bt, 'T', 'cylindrical poloidal field')
    call result%write('psi', psi, 'Wb/rad', 'cylindrical analogue poloidal flux')
    call result%write('q', q, '1', 'constant safety factor')
    call result%close()
end program gold_hoyle_demo
