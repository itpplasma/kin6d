program test_streaming_particle_density
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval
    use kin6d_streaming_comparison, only: streaming_state_t, streaming_particles_t, &
        init_streaming, init_streaming_particles
    use kin6d_streaming_particle_density
    implicit none

    call uniform_oracle()
    call arbitrary_p1_oracle()
    call rejection_oracles()
contains
    subroutine uniform_oracle()
        type(streaming_state_t) :: initial
        type(streaming_particles_t) :: particles
        type(particle_density_workspace_t) :: work
        type(particle_density_certificate_t) :: certificate
        real(dp) :: grid(2), a(0:1,2), b(0:1,2), c(0:1), s(0:1), t, x
        type(interval_t) :: ci(0:1), si(0:1), lower, upper, rounding, quad, padding
        type(interval_t) :: left_lower,left_upper,right_lower,right_upper
        real(qp) :: exact, amplitude, tau, chord_lo, chord_hi
        logical :: ok
        integer :: j, frame

        grid = [-1.0_dp,1.0_dp]
        a(0,:) = 1.0_dp
        a(1,:) = 0.5_dp
        b = 0.0_dp
        call init_streaming(initial,grid,a,b,ok)
        call require(ok,'uniform source')
        call init_streaming_particles(initial,32,24,particles,ok)
        call require(ok,'uniform tensor midpoint loader')
        call prepare_particle_density_workspace(initial,particles,work,ok)
        call require(ok,'preallocated particle workspace')
        do frame = 0, 3
            t = real(frame,dp)*0.75_dp
            call project_particle_density(particles,t,c,s,ok,work)
            call require(ok,'candidate weak Fourier projection')
            call certify_particle_density(initial,particles,t,certificate,ok,work)
            call require(ok,'uniform weak density certification')
            call particle_density_coefficients(initial,particles,t,certificate,ci,si,ok)
            call require(ok,'coefficient interval getter')
            amplitude = 0.5_qp
            if (t /= 0.0_dp) amplitude = 0.5_qp*sin(real(t,qp))/real(t,qp)
            call contains(ci(0),1.0_qp,'analytic mean')
            call contains(ci(1),amplitude,'analytic sinc amplitude')
            call contains(si(1),0.0_qp,'analytic sine projection')
            call contains(ci(0),real(c(0),qp),'stored candidate mean')
            call contains(ci(1),real(c(1),qp),'stored candidate harmonic')
            call particle_density_budgets(certificate,rounding,quad,ok)
            call require(ok .and. rounding%hi > 0.0_dp,'rounding budget retained')
            if (t == 0.0_dp) call require(quad%hi < 1.0e-300_dp,'P1 mean/time-zero exactness')
            do j = 0, 17
                x = real(j,dp)*0.31_dp
                call sample_particle_density_band(initial,particles,t,certificate, &
                    interval(x),lower,upper,ok)
                call require(ok,'pointwise finite-Fourier density enclosure')
                exact = 1.0_qp+amplitude*cos(real(x,qp))
                call contains(interval(lower%lo,upper%hi),exact,'uniform analytic density')
            end do
            call particle_density_chord(initial,particles,t,certificate,interval(0.0_dp), &
                interval(0.7_dp),padding,ok)
            call require(ok,'density chord remainder')
            call require(real(padding%hi,qp) >= abs(amplitude)*0.7_qp**2/8, &
                'independent chord curvature bound')
            call sample_particle_density_band(initial,particles,t,certificate, &
                interval(0.0_dp),left_lower,left_upper,ok)
            call require(ok,'left density chord endpoint')
            call sample_particle_density_band(initial,particles,t,certificate, &
                interval(0.7_dp),right_lower,right_upper,ok)
            call require(ok,'right density chord endpoint')
            do j = 1,15
                tau = real(j,qp)/16
                exact = 1.0_qp+amplitude*cos(tau*real(0.7_dp,qp))
                chord_lo = (1-tau)*real(left_lower%lo,qp)+ &
                    tau*real(right_lower%lo,qp)-real(padding%hi,qp)
                chord_hi = (1-tau)*real(left_upper%hi,qp)+ &
                    tau*real(right_upper%hi,qp)+real(padding%hi,qp)
                call require(chord_lo <= exact .and. exact <= chord_hi, &
                    'analytic density enclosed between connected vertices')
            end do
        end do
    end subroutine

    subroutine arbitrary_p1_oracle()
        type(streaming_state_t) :: initial
        type(streaming_particles_t) :: particles
        type(particle_density_workspace_t) :: work
        type(particle_density_certificate_t) :: certificate
        real(dp) :: grid(4), a(0:2,4), b(0:2,4), t, x, c(0:2), s(0:2)
        type(interval_t) :: ci(0:2), si(0:2), lower, upper, rounding, quad
        real(qp) :: cq(0:2),sq(0:2),exact
        real(dp) :: coarse
        logical :: ok
        integer :: k,j,nv,frame

        grid = [-1.0_dp,-0.625_dp,0.25_dp,1.0_dp]
        a(0,:) = [1.0_dp,1.25_dp,0.75_dp,1.125_dp]
        a(1,:) = [0.5_dp,-0.125_dp,0.25_dp,-0.25_dp]
        a(2,:) = [-0.125_dp,0.375_dp,0.0625_dp,-0.125_dp]
        b(0,:) = 0.0_dp
        b(1,:) = [-0.25_dp,0.0625_dp,0.375_dp,0.125_dp]
        b(2,:) = [0.25_dp,0.125_dp,-0.25_dp,0.0625_dp]
        call init_streaming(initial,grid,a,b,ok)
        call require(ok,'nonuniform arbitrary P1 source')
        do nv = 3, 6, 3
            call init_streaming_particles(initial,12,nv,particles,ok)
            call require(ok,'midpoints subdivide each original P1 cell')
            call prepare_particle_density_workspace(initial,particles,work,ok)
            call require(ok,'reuse workspace after changed particle count')
            do frame = 0, 2
                t = -1.25_dp+real(frame,dp)*1.25_dp
                call exact_coefficients(grid,a,b,t,cq,sq)
                call project_particle_density(particles,t,c,s,ok,work)
                call require(ok,'arbitrary candidate projection')
                call certify_particle_density(initial,particles,t,certificate,ok,work)
                call require(ok,'arbitrary P1 weak certificate')
                call particle_density_coefficients(initial,particles,t,certificate,ci,si,ok)
                call require(ok,'arbitrary coefficient enclosure')
                do k = 0, 2
                    call contains(ci(k),cq(k),'independent analytic P1 cosine integral')
                    call contains(si(k),sq(k),'independent analytic P1 sine integral')
                    call contains(ci(k),real(c(k),qp),'candidate cosine inside band')
                    call contains(si(k),real(s(k),qp),'candidate sine inside band')
                end do
                do j = 0, 19
                    x = real(j,dp)*0.23_dp
                    exact = cq(0)
                    do k = 1, 2
                        exact = exact+cq(k)*cos(real(k,qp)*real(x,qp))+ &
                            sq(k)*sin(real(k,qp)*real(x,qp))
                    end do
                    call sample_particle_density_band(initial,particles,t,certificate, &
                        interval(x),lower,upper,ok)
                    call require(ok,'arbitrary point density band')
                    call contains(interval(lower%lo,upper%hi),exact,'analytic arbitrary density')
                end do
                call particle_density_budgets(certificate,rounding,quad,ok)
                call require(ok,'separate quadrature and rounding budgets')
                if (frame == 0) then
                    if (nv == 3) coarse = quad%hi
                    if (nv == 6) call require(abs(4.0_dp*quad%hi/coarse-1.0_dp) < &
                        1.0e-12_dp,'independent quadratic subdivision scaling')
                end if
            end do
        end do
    end subroutine

    subroutine rejection_oracles()
        type(streaming_state_t) :: initial,other
        type(streaming_particles_t) :: particles,aliased,changed,empty
        type(particle_density_workspace_t) :: work,badwork
        type(particle_density_certificate_t) :: certificate
        real(dp) :: grid(3),a(0:2,3),b(0:2,3),c(0:2),s(0:2),nan
        type(interval_t) :: lower,upper,padding
        logical :: ok

        grid = [-1.0_dp,0.0_dp,1.0_dp]
        a = 0.25_dp
        a(0,:) = 1.0_dp
        b = 0.0_dp
        call init_streaming(initial,grid,a,b,ok)
        call require(ok,'rejection source')
        call init_streaming_particles(initial,8,4,particles,ok)
        call require(ok,'rejection particle loader')
        call prepare_particle_density_workspace(initial,particles,work,ok)
        call require(ok,'rejection workspace')
        call certify_particle_density(initial,particles,1.0_dp,certificate,ok,work)
        call require(ok,'rejection baseline certificate')
        call sample_particle_density_band(initial,particles,1.25_dp,certificate, &
            interval(0.0_dp),lower,upper,ok)
        call require(.not. ok,'reject stale time')
        call init_streaming_particles(initial,12,4,changed,ok)
        call require(ok,'different particle grid')
        call sample_particle_density_band(initial,changed,1.0_dp,certificate, &
            interval(0.0_dp),lower,upper,ok)
        call require(.not. ok,'reject stale particle geometry')
        a(1,2) = 0.375_dp
        call init_streaming(other,grid,a,b,ok)
        call require(ok,'different P1 source')
        call sample_particle_density_band(other,particles,1.0_dp,certificate, &
            interval(0.0_dp),lower,upper,ok)
        call require(.not. ok,'reject stale source snapshot')
        call init_streaming_particles(initial,4,4,aliased,ok)
        call require(ok,'loader permits candidate aliased schedule')
        call prepare_particle_density_workspace(initial,aliased,badwork,ok)
        call require(.not. ok,'reject insufficient exact Fourier quadrature degree')
        call project_particle_density(aliased,1.0_dp,c,s,ok)
        call require(.not. ok,'candidate rejects insufficient projection degree')
        call project_particle_density(empty,1.0_dp,c,s,ok)
        call require(.not. ok,'reject uninitialized particles')
        call particle_density_chord(initial,particles,1.0_dp,certificate, &
            interval(1.0_dp),interval(0.0_dp),padding,ok)
        call require(.not. ok,'reject reversed chord geometry')
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        call sample_particle_density_band(initial,particles,1.0_dp,certificate, &
            interval(nan),lower,upper,ok)
        call require(.not. ok,'reject nonfinite position')
        call certify_particle_density(initial,particles,nan,certificate,ok,work)
        call require(.not. ok,'reject nonfinite time')
        call sample_particle_density_band(initial,particles,1.0_dp,certificate, &
            interval(0.0_dp),lower,upper,ok)
        call require(.not. ok,'failed recertification invalidates stale certificate')
    end subroutine

    ! Independent real128 elementary antiderivatives, not generated expressions
    ! or a residual of the tested code. Coefficients are integrated cell by cell.
    subroutine exact_coefficients(grid,a,b,t,c,s)
        real(dp), intent(in) :: grid(:),a(0:,:),b(0:,:),t
        real(qp), intent(out) :: c(0:),s(0:)
        real(qp) :: vl,vr,al,ar,bl,br,q,sa,sb,ia,ib,jc,js,ic,is
        integer :: k,j

        c = 0.0_qp
        s = 0.0_qp
        do k = 0,size(a,1)-1
            q = real(k,qp)*real(t,qp)
            do j = 1,size(grid)-1
                vl = real(grid(j),qp); vr = real(grid(j+1),qp)
                al = real(a(k,j),qp); ar = real(a(k,j+1),qp)
                bl = real(b(k,j),qp); br = real(b(k,j+1),qp)
                if (q == 0.0_qp) then
                    c(k) = c(k)+(vr-vl)*(al+ar)/4.0_qp
                    s(k) = s(k)+(vr-vl)*(bl+br)/4.0_qp
                else
                    sa = (ar-al)/(vr-vl); sb = (br-bl)/(vr-vl)
                    ia = al-sa*vl; ib = bl-sb*vl
                    ic = (sin(q*vr)-sin(q*vl))/q
                    is = (cos(q*vl)-cos(q*vr))/q
                    jc = (vr*sin(q*vr)-vl*sin(q*vl))/q+ &
                        (cos(q*vr)-cos(q*vl))/q**2
                    js = (vl*cos(q*vl)-vr*cos(q*vr))/q+ &
                        (sin(q*vr)-sin(q*vl))/q**2
                    c(k) = c(k)+(ia*ic+sa*jc-ib*is-sb*js)/2.0_qp
                    s(k) = s(k)+(ia*is+sa*js+ib*ic+sb*jc)/2.0_qp
                end if
            end do
        end do
    end subroutine

    subroutine contains(box,value,label)
        type(interval_t), intent(in) :: box
        real(qp), intent(in) :: value
        character(*), intent(in) :: label

        if (real(box%lo,qp) <= value .and. value <= real(box%hi,qp)) return
        print '(a,3es25.16)',trim(label)//': ',box%lo,real(value,dp),box%hi
        error stop 1
    end subroutine

    subroutine require(condition,label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        print '(a)',trim(label)
        error stop 1
    end subroutine
end program test_streaming_particle_density
