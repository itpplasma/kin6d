program benchmark_error_continuum
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval
    use kin6d_laboratory, only: laboratory_t, initialize, advance, certificate_radius
    use kin6d_fp_laboratory, only: fp_state_t, initialize_fp, solve_fp, &
        fp_certificate_radius
    use poisson_laboratory, only: coupled_steady_t, solve_coupled_steady, &
        coupled_certificate_radius, poisson_mode_t, solve_poisson_mode
    use kin6d_poisson_leaves, only: poisson_coefficients_interval
    use kin6d_streaming_comparison, only: streaming_state_t, init_streaming, &
        advance_streaming, streaming_radius
    use kin6d_transport_laboratory, only: transport_state_t, initialize_transport, &
        solve_transport
    use kin6d_transport_certificate, only: transport_certificate_t, certify_transport
    use kin6d_slab_oscillation, only: slab_oscillation_state_t, &
        initialize_slab_oscillation, solve_slab_oscillation, slab_oscillation_radius
    use kin6d_shielding_dynamics, only: shielding_state_t, init_shielding, &
        advance_shielding, sample_shielding
    use kin6d_debye_sphere_reference, only: debye_sphere_values, debye_sphere_enclosure
    use kin6d_interaction_laboratory, only: interaction_state_t, init_interaction, &
        step_interaction, sample_interaction
    use kin6d_species, only: atomic_units, electron
    use error_benchmark_support
    implicit none
    type(laboratory_t) :: stream
    type(fp_state_t) :: fp
    type(coupled_steady_t) :: coupled
    type(poisson_mode_t) :: poisson
    type(streaming_state_t) :: comparison
    type(transport_state_t) :: transport
    type(transport_certificate_t) :: cert
    type(slab_oscillation_state_t) :: slab
    type(shielding_state_t) :: shield
    type(interaction_state_t) :: encounter
    type(interval_t) :: energy, pbox, bbox, tbox, mbox, mbbox, mtbox, rb(3), ub(3)
    type(interval_t) :: pc,ps,ec,es
    real(dp) :: a, b, c, d, solve, check, radius, residual, p, bare, tangent
    real(dp) :: lo(3), hi(3), values(3), r(3), u(3)
    real(dp) :: velocity(129),cosine(0:1,129),sine(0:1,129)
    integer :: trial, j, mode
    logical :: ok, bounded
    call header()
    do trial = 1, 5
        call initialize(stream, 256)
        call cpu_time(a)
        do j = 1, 64
            call advance(stream, 1._dp/64)
        end do
        call cpu_time(b)
        radius = certificate_radius(stream)
        call cpu_time(c)
        call row('streaming', 'S', trial, b-a, 0._dp, c-b, radius, 1._dp, &
            'certified-continuous', 'cells=256;steps=64;t=1', 'direct')
        velocity=[(-1._dp+real(j,dp)/64,j=0,128)]
        cosine(0,:)=1; cosine(1,:)=.5_dp; sine=0
        call init_streaming(comparison,velocity,cosine,sine,ok)
        call require(ok)
        call cpu_time(a)
        call advance_streaming(comparison,1._dp,ok,certify=.false.)
        call cpu_time(b)
        call require(ok); solve=b-a
        call cpu_time(a)
        call advance_streaming(comparison,1._dp,ok)
        call cpu_time(b)
        call require(ok)
        call streaming_radius(comparison,radius,ok)
        call require(ok)
        call row('streaming-comparison','S',trial,solve,0._dp,max(0._dp,b-a-solve), &
            radius,1._dp,'certified-continuous','cells=128;modes=1;t=1', &
            'paired-difference;prepared-initial-data')
        call initialize_fp(fp, 256)
        call cpu_time(a)
        do j = 1, 100
            call solve_fp(fp, residual)
        end do
        call cpu_time(b)
        solve = (b-a)/100
        call cpu_time(a)
        radius = fp_certificate_radius(fp)
        call cpu_time(b)
        call row('FP', 'S', trial, solve, 0._dp, b-a, radius, 1._dp, &
            'certified-continuous', 'cells=256;theta=1;vmax=2', 'direct')
        call cpu_time(a)
        do j = 1, 1000
            call solve_coupled_steady(.1_dp, 1._dp, 1, coupled, ok, &
                fp_rate=1._dp, temperature=1._dp)
        end do
        call cpu_time(b)
        call require(ok)
        solve = (b-a)/1000
        call cpu_time(a)
        do j = 1, 100
            call coupled_certificate_radius(coupled, radius, p, bare, energy, ok)
        end do
        call cpu_time(b)
        call require(ok)
        call row('Poisson-VPFP', 'S', trial, solve, 0._dp, (b-a)/100, &
            max(radius,p,bare), 1._dp, 'certified-continuous', &
            'k=1;nu=1;fp=1;drive=.1', 'direct')
        call cpu_time(a)
        do j=1,1000
            call solve_poisson_mode(2,.75_dp,0._dp,.1_dp,-.2_dp,poisson,ok)
        end do
        call cpu_time(b)
        call require(ok); solve=(b-a)/1000
        call cpu_time(a)
        do j=1,100
            call poisson_coefficients_interval(interval(2._dp),interval(.75_dp), &
                interval(.1_dp),interval(-.2_dp),pc,ps,ec,es)
        end do
        call cpu_time(b)
        radius=max(poisson%phi_cos-pc%lo,pc%hi-poisson%phi_cos, &
            poisson%phi_sin-ps%lo,ps%hi-poisson%phi_sin, &
            poisson%e_cos-ec%lo,ec%hi-poisson%e_cos, &
            poisson%e_sin-es%lo,es%hi-poisson%e_sin)
        call row('Poisson-mode','S',trial,solve,0._dp,(b-a)/100,radius,1._dp, &
            'certified-discrete coefficients','k=2;eps=.75;rho_cos=.1;rho_sin=-.2', &
            'direct')
        call initialize_transport(transport, 64, 2, 1._dp, ok)
        call require(ok)
        transport%source_cos = .1_dp; transport%source_sin = .05_dp
        transport%source_sin(0,:) = 0
        call cpu_time(a)
        do j = 1, 100
            call solve_transport(transport, ok)
        end do
        call cpu_time(b)
        call require(ok)
        solve = (b-a)/100
        call cpu_time(a)
        call certify_transport(transport, cert, ok)
        call cpu_time(b)
        call require(ok)
        call row('transport', 'D', trial, solve, 0._dp, b-a, cert%state_radius, &
            1._dp, 'certified-continuous', 'cells=64;modes=2;nu=1', 'direct')
        call initialize_slab_oscillation(slab, 128, 4._dp, 1._dp, 0._dp, ok)
        call cpu_time(a)
        call solve_slab_oscillation(slab, ok, certify=.false.)
        call cpu_time(b)
        call require(ok)
        solve = b-a
        call cpu_time(a)
        call solve_slab_oscillation(slab, ok)
        call cpu_time(b)
        call require(ok)
        call slab_oscillation_radius(slab, radius, ok)
        call require(ok)
        call row('slab', 'S', trial, solve, 0._dp, max(0._dp,b-a-solve), &
            radius, 1._dp, 'certified-continuous', 'steps=128;t=4', 'paired-difference')
        do mode = 1, 2
            call init_shielding(shield, [1._dp,.2_dp,.1_dp], ok)
            call cpu_time(a)
            call advance_shielding(shield, 1._dp, .001_dp, ok, certify=mode==2)
            call cpu_time(b)
            call require(ok)
            if (mode==1) solve = b-a
        end do
        call sample_shielding(shield, .25_dp, values, lo, hi, ok)
        call require(ok)
        call row('shielding', 'S', trial, solve, 0._dp, max(0._dp,b-a-solve), &
            maxval(max(values-lo,hi-values)), 1._dp, 'certified-continuous', &
            'modes=3;dt=.001;t=1;x=.25', 'paired-difference')
        call cpu_time(a)
        do j = 1, 1000
            call debye_sphere_values(.75_dp, 1.75_dp, p, bare, tangent, ok)
        end do
        call cpu_time(b)
        solve = (b-a)/1000
        call cpu_time(a)
        do j = 1, 100
            call debye_sphere_enclosure(interval(.75_dp), interval(1.75_dp), &
                pbox, bbox, tbox, mbox, mbbox, mtbox, ok)
        end do
        call cpu_time(b)
        call require(ok)
        radius = max(p-pbox%lo,pbox%hi-p)
        call row('Debye-sphere', 'S', trial, solve, 0._dp, (b-a)/100, radius, &
            1._dp, 'certified-continuous', 'r=.75;kappa=1.75', 'direct')
        do mode = 1, 2
            call init_interaction(encounter, [-6._dp,.6_dp,0._dp], &
                [2._dp,0._dp,0._dp], atomic_units(), electron(), electron(), ok)
            call cpu_time(a)
            do j = 1, 200
                call step_interaction(encounter, .0005_dp, ok, certify=mode==2)
                call require(ok)
            end do
            call cpu_time(b)
            if (mode==1) solve = b-a
        end do
        call sample_interaction(encounter, r, u, rb, ub, ok, bounded)
        call require(ok .and. bounded)
        radius = maxval(max(r-rb%lo,rb%hi-r))
        call row('C0-flow-tube', 'S', trial, solve, 0._dp, max(0._dp,b-a-solve), &
            radius, 6._dp, 'certified-continuous', 'steps=200;dt=.0005', &
            'paired-difference')
    end do
end program
