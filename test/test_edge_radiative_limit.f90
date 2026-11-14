program test_edge_radiative_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_edge_radiative_limit
    implicit none
    type(edge_model_t) :: p
    real(dp) :: n, ft, nlimit, t_fold, ngw, n2, t2, nl, nr, sl, sr
    integer :: rc
    p%j_parallel = 1e6_dp/acos(-1.0_dp)
    p%p_external = 3e5_dp
    call greenwald_reference_density(1e6_dp,1.0_dp,ngw,rc)
    if (rc /= EDGE_OK .or. abs(ngw/(1e20_dp/acos(-1.0_dp))-1)>1e-13_dp) &
        error stop 'Greenwald empirical reference'
    ! The algebraic equation must actually balance heating, transport, radiation.
    call edge_equilibrium_density(p,100.0_dp,n,ft,rc)
    if (rc /= EDGE_OK .or. n < 4e19_dp .or. n > 6e19_dp) &
        error stop 'hot branch equilibrium'
    call edge_density_limit(p,5.0_dp,800.0_dp,nlimit,t_fold,rc)
    if (rc /= EDGE_OK) error stop 'thermal fold missed'
    if (t_fold < 100 .or. t_fold > 140) error stop 'fold temperature'
    if (nlimit < 4e19_dp .or. nlimit > 6e19_dp) error stop 'fold density'
    ! Independent fold oracle: n_eq(T) has a local maximum, and F_T flips sign.
    call edge_equilibrium_density(p,t_fold*0.99_dp,nl,sl,rc)
    if (rc /= EDGE_OK) error stop 'low temperature'
    call edge_equilibrium_density(p,t_fold*1.01_dp,nr,sr,rc)
    if (rc /= EDGE_OK) error stop 'high temperature'
    if (nl >= nlimit .or. nr >= nlimit) error stop 'no genuine fold'
    if (sl <= 0 .or. sr >= 0) error stop 'wrong stability branch'
    ! Adding power raises critical density for this fixed transport/radiation law.
    p%p_external = 5e5_dp
    call edge_density_limit(p,5.0_dp,800.0_dp,n2,t2,rc)
    if (rc /= EDGE_OK .or. n2 <= nlimit) error stop 'power trend'
    ! Reject nonphysical configuration, report missing fold rather than an
    ! invented density limit or returning the Greenwald reference as theory.
    p%impurity_fraction = 0
    call edge_density_limit(p,5.0_dp,800.0_dp,n2,t2,rc)
    if (rc /= EDGE_INVALID) error stop 'missing impurity emission accepted'
    p%impurity_fraction = 0.005_dp
    p%p_external = 0
    call edge_density_limit(p,5.0_dp,800.0_dp,n2,t2,rc)
    if (rc /= EDGE_NO_FOLD) error stop 'no-fold case misidentified'
    print '(a,3es14.5)', 'edge fold n,T and empirical nGW ', nlimit,t_fold,ngw
    print '(a)', 'PASS edge radiative fold model'
end program
