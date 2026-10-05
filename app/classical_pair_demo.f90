program classical_pair_demo
    use, intrinsic :: iso_fortran_env, only: dp=>real64,int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use kin6d_classical_pair, only: classical_pair_model_t, initialize_classical_pair, &
        classical_pair_sample, classical_pair_wrap, classical_pair_state, &
        classical_pair_bin_center
    use kin6d_classical_pair_leaves, only: cp_statistic
    use kin6d_nbody, only: nbody_step, nbody_energy, nbody_encounter_step
    use fortplot, only: figure,figure_t,get_global_figure,subplot,scatter,plot, &
        errorbar,title,xlabel,ylabel,xlim,ylim,legend,tight_layout,savefig
    use fortplot_animation, only: animation_t,animation_func=>FuncAnimation,save_animation
    use fortplot_system_runtime, only: create_directory_runtime
    implicit none
    integer,parameter :: samples=4096,frames=61,bins=16,visible=12,fps=12
    ! Encounter-resolved step: eta times the Kepler/crossing time, capped.
    real(dp),parameter :: encounter_eta=0.02_dp,max_step=0.02_dp
    ! Hoeffding exponent h: 2*16*61*exp(-11) < 0.05 simultaneous coverage.
    real(dp),parameter :: confidence=11.0_dp
    real(dp),parameter :: blue(3)=[0.0_dp,0.4471_dp,0.6980_dp]
    real(dp),parameter :: orange(3)=[0.8353_dp,0.3686_dp,0.0_dp]
    type(classical_pair_model_t) :: model
    type(fortnum_status_t) :: status
    type(figure_t),pointer :: canvas
    type(animation_t) :: animation
    real(dp) :: state(12,samples),center(3),relative(3),energy0(samples),energy,wrapped(3)
    real(dp) :: time(frames),mean(bins,frames),sampling_radius(bins,frames)
    real(dp) :: particle(3,2,visible,frames),max_defect(frames),solve_cpu(frames),extra_cpu(frames)
    real(dp) :: xbin(bins),xp(visible),yp(visible),view,now,dt,start,finish,limit
    real(dp) :: defect
    integer :: counts(bins,frames),frame,sample,j,b,unit,snapshot_frame,encode_status,steps
    integer(int64) :: seed=104729_int64
    logical :: ok,snapshot_only
    character(64) :: argument,factor_text
    character(150) :: heading
    snapshot_only=.false.;snapshot_frame=31
    if(command_argument_count()>2) error stop 'usage: classical_pair_demo [--snapshot [frame]]'
    if(command_argument_count()>=1) then
        call get_command_argument(1,argument)
        if(trim(argument)/='--snapshot') error stop 'expected --snapshot'
        snapshot_only=.true.
    end if
    if(command_argument_count()==2) then
        call get_command_argument(2,argument)
        read(argument,*,iostat=encode_status) snapshot_frame
        if(encode_status/=0) error stop 'invalid snapshot frame'
    end if
    if(snapshot_frame<1.or.snapshot_frame>frames) error stop 'snapshot index outside 1..61'
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
        time(frame)=4*real(frame-1,dp)/15
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
            b=min(bins,1+floor(16*(wrapped(1)/model%length+0.5_dp)))
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
    write(unit,'(a)') '# periodic3D; L=2pi,m=1,g=2 (periodic_pair_units); opposite-spin labels do not affect classical force'
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
    view=maxval(mean)+maxval(sampling_radius)+1
    call figure(figsize=[12.0_dp,5.5_dp]);canvas=>get_global_figure()
    call draw_frame(snapshot_frame)
    call savefig('output/classical-coulomb-pair.png')
    if(.not.snapshot_only) then
        animation=animation_func(draw_frame,frames=frames,interval=83,fig=canvas)
        call save_animation(animation,'output/classical-coulomb-pair.mp4',fps=fps,status=encode_status)
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
    subroutine draw_frame(index)
        integer,intent(in) :: index
        integer :: p,q
        if(index<1.or.index>frames) error stop 'animation index outside 1..61'
        call canvas%clear()
        call subplot(1,2,1)
        do q=1,2
            do p=1,visible
                xp(p)=(particle(1,q,p,index)+particle(3,q,p,index)/4)/model%length
                yp(p)=(particle(2,q,p,index)+particle(3,q,p,index)/6)/model%length
            end do
            if(q==1) then
                call scatter(xp,yp,color=blue,marker='o',s=36.0_dp,label='particle 1')
            else
                call scatter(xp,yp,color=orange,marker='s',s=36.0_dp,label='particle 2')
            end if
        end do
        write(heading,'(a,f5.2)') '12 independent classical 3D pairs; t=',time(index)
        call title(trim(heading));call xlabel('(x+z/4)/L');call ylabel('(y+z/6)/L')
        call xlim(-0.7_dp,0.7_dp);call ylim(-0.65_dp,0.9_dp);call legend()
        call subplot(1,2,2)
        call errorbar(xbin,mean(:,index),yerr=sampling_radius(:,index),fmt='o', &
                      color=blue,capsize=3.0_dp,label='4096 pairs: simultaneous 95% sampling CI')
        call plot([-0.5_dp,0.5_dp],[1.0_dp,1.0_dp],linestyle='--',color=orange,label='uniform initial law')
        call title('3D slab pair correlation; numerical bounds unavailable')
        call xlabel('relative displacement x/L; slab widths/L=(1/16,1,1)')
        call ylabel('g = 16 count / 4096');call xlim(-0.5_dp,0.5_dp)
        call ylim(-0.7_dp,view);call legend();call tight_layout()
        call canvas%set_rendered(.false.)
    end subroutine
end program
