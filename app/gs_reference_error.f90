program gs_reference_error
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, interval_pi, sqrt, &
        operator(+), operator(*), operator(/)
    use kin6d_gs_reference_error, only: gs_reference_parameters_t, &
        enclose_gs_reference_cell
    implicit none
    type(gs_reference_parameters_t) :: parameters
    type(interval_t) :: total(3), cell(3), area, volume, a, v, maximum, peak
    type(interval_t) :: rms, energy, current, force
    real(dp), allocatable :: points(:,:), psi(:)
    integer, allocatable :: triangles(:,:)
    real(dp) :: r0, c, edge, f, mu0, pp
    integer :: depth, unit, n, m, i, marker, stat
    character(2048) :: path, param_path, header
    logical :: ok
    namelist /gs_reference/ r0,c,edge,f,mu0,pp,depth
    if (command_argument_count() /= 2) &
        error stop 'usage: gs_reference_error mesh params'
    r0=ieee_value(0.0_dp,ieee_quiet_nan)
    c=r0;edge=r0;f=r0;mu0=r0;pp=r0;depth=-1
    call get_command_argument(1,path); call get_command_argument(2,param_path)
    open(newunit=unit,file=trim(param_path),status='old',action='read')
    read(unit,nml=gs_reference,iostat=stat); close(unit)
    if (stat /= 0) error stop 'invalid reference configuration'
    parameters=gs_reference_parameters_t(r0,c,edge,f,mu0,pp)
    open(newunit=unit,file=trim(path),status='old',action='read')
    read(unit,'(a)') header
    read(unit,*) n,m
    if (n<3 .or. m<1) error stop 'empty mesh'
    allocate(points(2,n),psi(n),triangles(3,m))
    do i=1,n
        read(unit,*) points(:,i),psi(i),marker
    end do
    read(unit,'(a)') header
    do i=1,m
        read(unit,*) triangles(:,i)
    end do
    close(unit)
    if (minval(triangles)<1 .or. maxval(triangles)>n) error stop 'bad connectivity'
    total=interval(0); area=interval(0); volume=interval(0); peak=interval(0)
    do i=1,m
        call enclose_gs_reference_cell(points(:,triangles(:,i)), &
            psi(triangles(:,i)),parameters,depth,cell,a,v,maximum,ok)
        if (.not. ok) error stop 'reference enclosure failed'
        total=total+cell; area=area+a; volume=volume+v
        peak%hi=max(peak%hi,maximum%hi)
    end do
    rms=sqrt(total(1)/area)
    energy=sqrt(total(2))
    current=energy/interval(mu0)
    force=sqrt(2*interval_pi()*total(3))
    current%lo=0; force%lo=0
    print '(a)', '{'
    print '(a,i0,a)', '  "nodes":',n,','
    print '(a,i0,a)', '  "triangles":',m,','
    print '(a,i0,a)', '  "cover_depth":',depth,','
    call emit('area_m2',area,.true.)
    call emit('toroidal_volume_m3',2*interval_pi()*volume,.true.)
    call emit('Bpol_area_RMS_error_T',rms,.true.)
    call emit('Bpol_whole_cell_max_error_T',sqrt(peak),.true.)
    call emit('weighted_gradient_error_Wb_per_rad_per_sqrt_m',energy,.true.)
    call emit('current_dual_norm_upper_A_sqrt_m',current,.true.)
    call emit('force_dual_norm_upper_N_per_sqrt_m',force,.false.)
    print '(a)', '}'
contains
    subroutine emit(name,value,comma)
        character(*),intent(in)::name
        type(interval_t),intent(in)::value
        logical,intent(in)::comma
        character(1)::suffix
        suffix=' '
        if(comma)suffix=','
        print '(a,es25.16e3,a,es25.16e3,a)', &
            '  "'//name//'": [',value%lo,',',value%hi,']'//suffix
    end subroutine emit
end program gs_reference_error
