program test_classical_pair
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use kin6d_classical_pair, only: classical_pair_model_t, &
        initialize_classical_pair, classical_pair_force, classical_pair_step, &
        classical_pair_force_jvp
    use kin6d_classical_pair_leaves, only: cp_real,cp_real_jvp,cp_real_vjp,cp_statistic,cp_bin_center
    implicit none
    type(classical_pair_model_t) :: a,b
    real(dp) :: r(3),f(3),fm(3),df(3),u,um,du,fp(3),up,h,weights(4),g(5)
    real(dp) :: position(3),v(3),energy,e0,dot1,dot2
    integer :: j
    logical :: ok
    call initialize_classical_pair(a,ok)
    if (.not.ok) error stop 'model initialization'
    b%alpha=0.35_dp; b%real_cutoff=3
    call initialize_classical_pair(b,ok)
    if (.not.ok) error stop 'second model initialization'
    r=[0.7_dp,-1.1_dp,0.3_dp]
    call classical_pair_force(a,r,f,u,ok)
    if (.not.ok) error stop 'force'
    call classical_pair_force(a,-r,fm,um,ok)
    if (maxval(abs(f+fm))>1e-12_dp) error stop 'pair antisymmetry'
    if (abs(u-um)>1e-12_dp) error stop 'pair parity'
    call classical_pair_force(b,r,fm,um,ok)
    if (maxval(abs(f-fm))>1e-8_dp) error stop 'independent Ewald splitting'
    if (abs(u-um)>1e-8_dp) error stop 'background energy splitting'
    h=1e-5_dp
    do j=1,3
        position=r; position(j)=r(j)+h
        call classical_pair_force(a,position,fp,up,ok)
        position(j)=r(j)-h
        call classical_pair_force(a,position,fm,um,ok)
        if (abs((up-um)/(2*h)+f(j))>1e-8_dp) error stop 'energy gradient'
        position=0.0_dp; position(j)=1.0_dp
        call classical_pair_force_jvp(a,r,position,df,du,ok)
        if (maxval(abs(df-(fp-fm)/(2*h)))>1e-8_dp) error stop 'force JVP'
    end do
    weights=[0.2_dp,-0.7_dp,0.4_dp,0.3_dp]
    call cp_real_jvp(r(1),r(2),r(3),a%alpha,a%coupling, &
        0.2_dp,-0.1_dp,0.3_dp,0.4_dp,-0.2_dp,df(1),df(2),df(3),du)
    call cp_real_vjp(r(1),r(2),r(3),a%alpha,a%coupling, &
        weights(1),weights(2),weights(3),weights(4),g(1),g(2),g(3),g(4),g(5))
    dot1=dot_product(weights,[df,du]);dot2=dot_product(g,[0.2_dp,-0.1_dp,0.3_dp,0.4_dp,-0.2_dp])
    if(abs(dot1-dot2)>1e-12_dp) error stop 'leaf JVP/VJP pairing'
    position=r; v=0.0_dp
    call classical_pair_force(a,position,f,e0,ok)
    do j=1,100
        call classical_pair_step(a,position,v,f,0.001_dp,energy,ok)
        if(.not.ok) error stop 'step'
    end do
    if(abs(energy-e0)>1e-7_dp) error stop 'short-run Hamiltonian drift'
    call cp_bin_center(1.0_dp,u)
    if(abs(u+15.0_dp/32)>1e-14_dp) error stop 'first bin center'
    call cp_bin_center(16.0_dp,u)
    if(abs(u-15.0_dp/32)>1e-14_dp) error stop 'last bin center'
    call cp_statistic(256.0_dp,4096.0_dp,11.0_dp,u,um)
    if(abs(u-1)>1e-14_dp) error stop 'histogram normalization'
    if(abs(um-16*sqrt(11.0_dp/(2*4096)))>1e-14_dp) error stop 'sampling radius'
    if(2*16*61*exp(-11.0_dp)>=.05_dp) error stop 'simultaneous coverage allocation'
    print '(a)', 'PASS: Ewald splitting, parity, energy gradient, AD and Hamiltonian'
end program
