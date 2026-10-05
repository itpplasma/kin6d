module kin6d_classical_pair
    use, intrinsic :: iso_fortran_env, only: dp=>real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_classical_pair_leaves, only: &
        cp_axpy, cp_background, cp_chart, &
        cp_difference, cp_energy, cp_fraction, &
        cp_phase_product, cp_phasor, cp_real, &
        cp_real_jvp, cp_recip_jvp, cp_recip_phase, &
        cp_relative_kick, cp_scale, cp_trig, &
        cp_weight
    implicit none
    private
    public :: classical_pair_model_t, initialize_classical_pair, classical_pair_force
    public :: classical_pair_step, classical_pair_sample, classical_pair_wrap
    public :: classical_pair_force_jvp
    type :: classical_pair_model_t
        real(dp) :: length=0.0_dp, alpha=0.4_dp, coupling=2.0_dp
        integer :: real_cutoff=2, reciprocal_cutoff=4, modes=0
        integer :: wave(3,364)=0
        real(dp) :: weight(364)=0.0_dp
        logical :: ready=.false.
    end type
contains
    subroutine initialize_classical_pair(model,ok)
        type(classical_pair_model_t), intent(inout) :: model
        logical, intent(out) :: ok
        integer :: i,j,k,m
        real(dp) :: next
        ok=.false.; model%ready=.false.; model%modes=0
        if (.not.all(ieee_is_finite([model%alpha,model%coupling]))) return
        if (model%alpha<=0.0_dp.or.model%coupling<=0.0_dp) return
        if (model%real_cutoff<1.or.model%reciprocal_cutoff<1) return
        if (model%reciprocal_cutoff>4) return
        call cp_scale(1.0_dp,model%length)
        do i=-model%reciprocal_cutoff,model%reciprocal_cutoff
            do j=-model%reciprocal_cutoff,model%reciprocal_cutoff
                do k=-model%reciprocal_cutoff,model%reciprocal_cutoff
                    ! One of each +/-k pair; doubled weight below, k=0 omitted.
                    if (i<0) cycle
                    if (i==0.and.j<0) cycle
                    if (i==0.and.j==0.and.k<=0) cycle
                    m=model%modes+1; model%modes=m
                    model%wave(:,m)=[i,j,k]
                    call cp_weight(real(i,dp),real(j,dp),real(k,dp), &
                                   model%alpha,model%coupling,model%weight(m))
                    call cp_axpy(model%weight(m),model%weight(m),1.0_dp,next)
                    model%weight(m)=next
                end do
            end do
        end do
        model%ready=.true.; ok=.true.
    end subroutine

    subroutine classical_pair_wrap(position,length,wrapped)
        real(dp), intent(in) :: position(3),length
        real(dp), intent(out) :: wrapped(3)
        integer :: j,chart
        real(dp) :: fraction,offset
        do j=1,3
            call cp_chart(position(j),length,fraction)
            chart=floor(fraction)
            call cp_scale(real(chart,dp),offset)
            call cp_difference(position(j),offset,wrapped(j))
        end do
    end subroutine

    subroutine classical_pair_force(model,position,force,potential,ok)
        type(classical_pair_model_t), intent(in) :: model
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: force(3),potential
        logical, intent(out) :: ok
        real(dp) :: r(3),image(3),f(3),u,c(-4:4,3),s(-4:4,3),cc,ss,next(3),nu
        integer :: i,j,k,m,a
        ok=.false.; force=0.0_dp; potential=0.0_dp
        if (.not.model%ready.or..not.all(ieee_is_finite(position))) return
        call classical_pair_wrap(position,model%length,r)
        if (maxval(abs(r))<1.0e-12_dp) return
        call cp_background(model%alpha,model%coupling,potential)
        do i=-model%real_cutoff,model%real_cutoff
            do j=-model%real_cutoff,model%real_cutoff
                do k=-model%real_cutoff,model%real_cutoff
                    call cp_axpy(r,model%length,real([i,j,k],dp),image)
                    call cp_real(image(1),image(2),image(3),model%alpha, &
                                 model%coupling,f(1),f(2),f(3),u)
                    call cp_axpy(force,f,1.0_dp,next); force=next
                    call cp_axpy(potential,u,1.0_dp,nu); potential=nu
                end do
            end do
        end do
        c(0,:)=1.0_dp; s(0,:)=0.0_dp
        do a=1,3
            call cp_trig(r(a),c(1,a),s(1,a))
            do j=2,model%reciprocal_cutoff
                call cp_phasor(c(j-1,a),s(j-1,a),c(1,a),s(1,a),c(j,a),s(j,a))
            end do
            do j=1,model%reciprocal_cutoff
                c(-j,a)=c(j,a); s(-j,a)=-s(j,a)
            end do
        end do
        do m=1,model%modes
            i=model%wave(1,m); j=model%wave(2,m); k=model%wave(3,m)
            call cp_phase_product(c(i,1),s(i,1),c(j,2),s(j,2),c(k,3),s(k,3),cc,ss)
            call cp_recip_phase(real(i,dp),real(j,dp),real(k,dp),model%weight(m), &
                                cc,ss,f(1),f(2),f(3),u)
            call cp_axpy(force,f,1.0_dp,next); force=next
            call cp_axpy(potential,u,1.0_dp,nu); potential=nu
        end do
        ok=all(ieee_is_finite(force)).and.ieee_is_finite(potential)
    end subroutine

    subroutine classical_pair_force_jvp(model,position,dposition,dforce,du,ok)
        type(classical_pair_model_t), intent(in) :: model
        real(dp), intent(in) :: position(3),dposition(3)
        real(dp), intent(out) :: dforce(3),du
        logical, intent(out) :: ok
        real(dp) :: r(3),image(3),f(3),u,next(3),nu
        integer :: i,j,k,m
        ok=.false.; dforce=0.0_dp; du=0.0_dp
        if (.not.model%ready) return
        call classical_pair_wrap(position,model%length,r)
        if (maxval(abs(r))<1.0e-12_dp) return
        do i=-model%real_cutoff,model%real_cutoff
            do j=-model%real_cutoff,model%real_cutoff
                do k=-model%real_cutoff,model%real_cutoff
                    call cp_axpy(r,model%length,real([i,j,k],dp),image)
                    call cp_real_jvp(image(1),image(2),image(3),model%alpha, &
                        model%coupling,dposition(1),dposition(2),dposition(3), &
                        0.0_dp,0.0_dp,f(1),f(2),f(3),u)
                    call cp_axpy(dforce,f,1.0_dp,next); dforce=next
                    call cp_axpy(du,u,1.0_dp,nu); du=nu
                end do
            end do
        end do
        do m=1,model%modes
            call cp_recip_jvp(r(1),r(2),r(3),real(model%wave(1,m),dp), &
                real(model%wave(2,m),dp),real(model%wave(3,m),dp),model%weight(m), &
                dposition(1),dposition(2),dposition(3),0.0_dp,0.0_dp,0.0_dp, &
                0.0_dp,f(1),f(2),f(3),u)
            call cp_axpy(dforce,f,1.0_dp,next); dforce=next
            call cp_axpy(du,u,1.0_dp,nu); du=nu
        end do
        ok=all(ieee_is_finite(dforce)).and.ieee_is_finite(du)
    end subroutine

    subroutine classical_pair_step(model,position,velocity,force,dt,energy,ok)
        type(classical_pair_model_t), intent(in) :: model
        real(dp), intent(inout) :: position(3),velocity(3),force(3)
        real(dp), intent(in) :: dt
        real(dp), intent(out) :: energy
        logical, intent(out) :: ok
        real(dp) :: next(3),u,half
        call cp_fraction(dt,2.0_dp,half)
        call cp_relative_kick(velocity,force,half,next); velocity=next
        call cp_axpy(position,velocity,dt,next); position=next
        call classical_pair_force(model,position,force,u,ok)
        if (.not.ok) return
        call cp_relative_kick(velocity,force,half,next); velocity=next
        call cp_energy(velocity(1),velocity(2),velocity(3),u,energy)
    end subroutine

    subroutine classical_pair_sample(seed,center,relative)
        integer(int64), intent(inout) :: seed
        real(dp), intent(out) :: center(3),relative(3)
        real(dp) :: draws(6),u
        integer :: j
        ! Park--Miller; exact integer arithmetic, independent draws by stream.
        do j=1,6
            seed=modulo(16807_int64*seed,2147483647_int64)
            call cp_fraction(real(seed,dp),2147483647.0_dp,u)
            call cp_scale(u,draws(j))
        end do
        center=draws(:3); relative=draws(4:)
        ! Uniform torus center and independent uniform relative displacement
        ! imply uniform independent particle positions, each p=0.
    end subroutine
end module
