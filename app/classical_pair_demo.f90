program classical_pair_demo
    use, intrinsic :: iso_fortran_env, only: dp=>real64,int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_classical_pair, only: classical_pair_model_t, initialize_classical_pair, &
        classical_pair_sample, classical_pair_wrap, classical_pair_state, &
        classical_pair_bin_center
    use kin6d_classical_pair_leaves, only: cp_statistic
    use kin6d_nbody, only: nbody_step, nbody_energy, nbody_encounter_step
    use fortplot, only: figure,figure_t,get_global_figure,subplot,scatter,plot, &
        errorbar,title,xlabel,ylabel,xlim,ylim,legend,tight_layout,savefig,suptitle
    use kin6d_demo_style, only: vermillion,black,movie_fps,hold_frames,num
    use fortplot_animation, only: animation_t,animation_func=>FuncAnimation,save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    use kin6d_result_io, only: run_t,result_file_t
    use kin6d_periodic_pair_config, only: periodic_pair_t,read_periodic_pair, &
        write_periodic_pair,declared_periodic_pair,periodic_pair_bins
    implicit none
    ! The 61 shared periodic-pair frames are each shown three times at 24 fps
    ! (8 physical frames per second, 7.6 s) plus a 1 s hold; no interpolation.
    integer,parameter :: bins=periodic_pair_bins,visible=12,repeat_frames=3
    ! Width/height of the projection axes box in pixels (measured): equal aspect.
    real(dp),parameter :: box_aspect=1.08_dp
    integer :: samples=4096
    ! Encounter-resolved step: eta times the Kepler/crossing time, capped.
    real(dp) :: encounter_eta=0.02_dp,max_step=0.02_dp
    ! Hoeffding exponent h: 2*bins*frames*exp(-h) < 0.05 simultaneous coverage.
    real(dp) :: confidence=11.0_dp
    real(dp),parameter :: blue(3)=[0.0_dp,0.4471_dp,0.6980_dp]
    type(classical_pair_model_t) :: model
    type(fortnum_status_t) :: status
    type(figure_t),pointer :: canvas
    type(animation_t) :: animation
    type(run_t) :: run
    type(result_file_t) :: result
    type(periodic_pair_t) :: pair,declared
    real(dp),allocatable :: state(:,:),energy0(:),time(:),mean(:,:),sampling_radius(:,:)
    real(dp),allocatable :: particle(:,:,:,:),max_defect(:),solve_cpu(:),extra_cpu(:)
    real(dp) :: center(3),relative(3),energy,wrapped(3)
    real(dp) :: xbin(bins),xp(visible),yp(visible),view,now,dt,start,finish,limit
    real(dp) :: defect
    integer,allocatable :: counts(:,:)
    integer :: frames,frame,sample,j,b,unit,snapshot_frame,encode_status,steps,ios
    integer(int64) :: seed=104729_int64
    logical :: ok,snapshot_only
    character(64) :: argument,factor_text
    character(150) :: heading
    character(256) :: msg
    namelist /classical_pair/ samples,encounter_eta,max_step,confidence,seed
    call run%start('classical_pair_demo')
    call read_periodic_pair(run,pair)
    ios=0;msg=''
    call run%open_config(unit)
    if(run%has_config()) read(unit,nml=classical_pair,iostat=ios,iomsg=msg)
    call run%close_config(unit,'classical_pair',ios,msg)
    call run%open_resolved(unit)
    call write_periodic_pair(unit,pair)
    write(unit,nml=classical_pair)
    call run%close_resolved(unit)
    ! The classical force law is fixed by the declared electron pair; the
    ! shared &periodic_pair may not silently describe a different system.
    declared=declared_periodic_pair()
    if(pair%mass/=declared%mass.or.pair%coupling/=declared%coupling) error stop &
        '&periodic_pair mass/coupling differ from the declared electron pair; '// &
        'classical_pair_demo realizes only that pair'
    if(.not.(pair%length>0.0_dp.and.pair%time_max>0.0_dp)) error stop &
        '&periodic_pair needs positive length and time_max'
    frames=pair%frames
    if(frames<2.or.frames>1001) error stop '&periodic_pair frames outside 2..1001'
    if(samples<1.or.visible>samples) error stop 'samples must cover the visible pairs'
    if(.not.(encounter_eta>0.0_dp.and.max_step>0.0_dp)) error stop 'invalid step control'
    if(seed<1_int64) error stop 'seed must be positive'
    if(2*bins*frames*exp(-confidence)>=0.05_dp) error stop &
        'confidence exponent too small for simultaneous 95% coverage'
    allocate(state(12,samples),energy0(samples),time(frames),mean(bins,frames), &
        sampling_radius(bins,frames),particle(3,2,visible,frames),max_defect(frames), &
        solve_cpu(frames),extra_cpu(frames),counts(bins,frames))
    snapshot_only=.false.;snapshot_frame=(frames+1)/2
    if(run%argument_count()>2) error stop &
        'usage: classical_pair_demo [--config FILE] [--snapshot [frame]]'
    if(run%argument_count()>=1) then
        call run%argument(1,argument)
        if(trim(argument)/='--snapshot') error stop 'expected --snapshot'
        snapshot_only=.true.
    end if
    if(run%argument_count()==2) then
        call run%argument(2,argument)
        read(argument,*,iostat=encode_status) snapshot_frame
        if(encode_status/=0) error stop 'invalid snapshot frame'
    end if
    if(snapshot_frame<1.or.snapshot_frame>frames) error stop 'snapshot index outside frames'
    model%length=pair%length
    call initialize_classical_pair(model,ok)
    if(.not.ok) error stop 'Ewald setup failed'
    call create_directory_runtime('output',ok)
    if(.not.ok) error stop 'output directory failed'
    solve_cpu=0.0_dp;extra_cpu=0.0_dp;max_defect=0.0_dp;counts=0
    do sample=1,samples
        call classical_pair_sample(seed,model%length,center,relative)
        call classical_pair_wrap(relative,model%length,wrapped)
        call classical_pair_state(center,wrapped,state(:,sample))
        call nbody_energy(model%system,state(:,sample),energy0(sample),status)
        if(status%code/=FORTNUM_OK) call failed(sample,1)
    end do
    do frame=1,frames
        time(frame)=pair%time_max*real(frame-1,dp)/real(frames-1,dp)
        if(frame>1) then
            call cpu_time(start)
            do sample=1,samples
                now=time(frame-1);steps=0
                do while(now<time(frame))
                    call nbody_encounter_step(model%system,state(:,sample),encounter_eta, &
                                              limit,status)
                    if(status%code/=FORTNUM_OK) call failed(sample,frame)
                    dt=min(limit,max_step,time(frame)-now)
                    call nbody_step(model%system,state(:,sample),dt,status)
                    if(status%code/=FORTNUM_OK) call failed(sample,frame)
                    now=now+dt
                    steps=steps+1
                    if(steps>1000000) call failed(sample,frame)
                end do
            end do
            call cpu_time(finish);solve_cpu(frame)=finish-start
        end if
        call cpu_time(start)
        do sample=1,samples
            call classical_pair_wrap(state(:3,sample)-state(4:6,sample),model%length,wrapped)
            b=min(bins,1+floor(bins*(wrapped(1)/model%length+0.5_dp)))
            if(b<1) call failed(sample,frame)
            counts(b,frame)=counts(b,frame)+1
            call nbody_energy(model%system,state(:,sample),energy,status)
            if(status%code/=FORTNUM_OK) call failed(sample,frame)
            defect=abs(energy-energy0(sample))/(1+abs(energy0(sample)))
            max_defect(frame)=max(max_defect(frame),defect)
            if(sample<=visible) then
                do j=1,2
                    call classical_pair_wrap(state(3*j-2:3*j,sample),model%length, &
                                             particle(:,j,sample,frame))
                end do
            end if
        end do
        do b=1,bins
            call cp_statistic(real(counts(b,frame),dp),real(samples,dp),confidence, &
                              mean(b,frame),sampling_radius(b,frame))
        end do
        call cpu_time(finish);extra_cpu(frame)=finish-start
    end do
    open(newunit=unit,file='output/classical-pair-bins.csv',status='replace')
    write(unit,'(a,es24.16e3,a)') '# periodic3D; L=',model%length, &
        ',m=1,g=2 (periodic_pair_units); opposite-spin labels do not affect classical force'
    write(unit,'(a)') '# numerical trajectory/truncation/model certification unavailable; no rejected draws'
    write(unit,'(a)') '# simultaneous >=95% Hoeffding sampling intervals for numerical-map bin law'
    write(unit,'(a,es10.3,a,es10.3,a,i0,a,i0,a,es10.3,a)') '# bins widths/L=(1/16,1,1); Ewald tol=', &
        model%tolerance,',alpha=',model%ewald%alpha,',real=',model%ewald%real_cutoff, &
        ',reciprocal=',model%ewald%reciprocal_cutoff,',bound=', &
        model%ewald%potential_bound*model%length, &
        ' g/L;k=0 omitted'
    write(unit,'(a)') 'frame,t,bin,x_over_L,count,samples,g,sampling_radius,status'
    do frame=1,frames
        do b=1,bins
            xbin(b)=classical_pair_bin_center(b)
            write(unit,'(i0,a,es24.16e3,a,i0,a,es24.16e3,2(a,i0),2(a,es24.16e3),a)') &
                frame,',',time(frame),',',b,',',xbin(b),',',counts(b,frame),',',samples, &
                ',',mean(b,frame),',',sampling_radius(b,frame),',candidate+sampling-CI'
        end do
    end do
    close(unit)
    open(newunit=unit,file='output/classical-pair-timing.csv',status='replace')
    write(unit,'(a)') 'frame,t,solve_cpu_s,diagnostic_sampling_extra_cpu_s,max_energy_defect,extra_multiplier,failed_samples'
    do frame=1,frames
        factor_text='unavailable'
        if(solve_cpu(frame)>0.0_dp) then
            write(factor_text,'(es24.16e3)') 1+extra_cpu(frame)/solve_cpu(frame)
        end if
        write(unit,'(i0,4(a,es24.16e3),a)') frame,',',time(frame),',',solve_cpu(frame), &
            ',',extra_cpu(frame),',',max_defect(frame),','//trim(adjustl(factor_text))//',0'
    end do
    close(unit)
    open(newunit=unit,file='output/classical-pair-particles.csv',status='replace')
    write(unit,'(a)') 'frame,t,ensemble,particle,x,y,z,status'
    do frame=1,frames
        do sample=1,visible
            do j=1,2
                write(unit,'(i0,a,es24.16e3,2(a,i0),3(a,es24.16e3),a)') frame,',',time(frame), &
                    ',',sample,',',j,',',particle(1,j,sample,frame),',', &
                    particle(2,j,sample,frame),',',particle(3,j,sample,frame),',candidate'
            end do
        end do
    end do
    close(unit)
    call result%create(run)
    call result%attribute('certificate_scope','simultaneous >=95% Hoeffding sampling '// &
        'intervals for the numerical-map bin law; trajectory/truncation certification unavailable')
    call result%write('time',time,'periodic-pair time unit','frame time',dims='frame')
    call result%write('bin_center',xbin,'1','relative displacement bin center x/L',dims='bin')
    call result%write('pair_correlation',mean,'1','slab pair correlation g = bins count/samples', &
        dims='bin,frame',status='candidate+sampling-CI')
    call result%write('sampling_radius',sampling_radius,'1', &
        'simultaneous 95% Hoeffding sampling radius of g',dims='bin,frame', &
        status='probabilistically-certified sampling interval')
    call result%write('max_energy_defect',max_defect,'1', &
        'max relative energy defect over samples',dims='frame',status='diagnostic')
    call result%write('solve_cpu',solve_cpu,'s','trajectory solve CPU time',dims='frame')
    call result%close()
    view=maxval(mean)+maxval(sampling_radius)+1
    call figure(figsize=[12.8_dp,6.4_dp]);canvas=>get_global_figure()
    call draw_frame(snapshot_frame)
    call savefig('output/classical-coulomb-pair.png')
    if(.not.snapshot_only) then
        animation=animation_func(update_frame,frames=repeat_frames*frames+hold_frames, &
            interval=42,fig=canvas)
        call save_animation(animation,'output/classical-coulomb-pair.mp4',fps=movie_fps, &
            status=encode_status)
        if(encode_status/=0) error stop 'movie encoding failed'
    end if
contains
    subroutine failed(sample,index)
        integer,intent(in) :: sample,index
        integer :: failure_unit
        open(newunit=failure_unit,file='output/classical-pair-failed.txt',status='replace')
        write(failure_unit,'(a,2(i0,1x))') 'No sample replaced; unresolved trajectory sample,frame: ',sample,index
        close(failure_unit)
        error stop 'unresolved classical trajectory; no movie promoted'
    end subroutine
    subroutine update_frame(number)
        integer,intent(in) :: number
        call draw_frame(min((number-1)/repeat_frames+1,frames))
    end subroutine
    subroutine draw_frame(index)
        integer,intent(in) :: index
        integer :: p,q
        if(index<1.or.index>frames) error stop 'animation index outside frames'
        call canvas%clear()
        call suptitle('Classical Coulomb repulsion empties small separations: g dips '// &
            'at x = 0 (bars: sampling CIs, not numerical bounds)')
        call subplot(1,2,1)
        do q=1,2
            do p=1,visible
                xp(p)=(particle(1,q,p,index)+particle(3,q,p,index)/4)/model%length
                yp(p)=(particle(2,q,p,index)+particle(3,q,p,index)/6)/model%length
            end do
            if(q==1) then
                call scatter(xp,yp,color=blue,marker='o',s=36.0_dp)
            else
                call scatter(xp,yp,color=vermillion,marker='s',s=36.0_dp)
            end if
        end do
        call title('12 independent 3D pairs (o: electron 1, s: electron 2)  |  t = '// &
            num(time(index),2))
        call xlabel('(x + z/4) / L');call ylabel('(y + z/6) / L')
        call xlim(-0.775_dp*box_aspect,0.775_dp*box_aspect);call ylim(-0.65_dp,0.9_dp)
        call subplot(1,2,2)
        write(heading,'(i0)') samples
        call errorbar(xbin,mean(:,index),yerr=sampling_radius(:,index),fmt='o', &
                      color=blue,capsize=3.0_dp,label=trim(heading)//' pairs ± simultaneous 95% CI')
        call plot([-0.5_dp,0.5_dp],[1.0_dp,1.0_dp],linestyle='--',color=black, &
                  linewidth=1.2_dp,label='uniform initial law')
        call title('Relative-position slab histogram g(x)')
        call xlabel('relative displacement x / L   (slab 1/16 x 1 x 1)')
        write(heading,'(a,i0,a,i0)') 'g = ',bins,' count / ',samples
        call ylabel(trim(heading));call xlim(-0.5_dp,0.5_dp)
        call ylim(-0.7_dp,view);call legend(loc='upper center');call tight_layout()
        call canvas%set_rendered(.false.)
    end subroutine
end program
