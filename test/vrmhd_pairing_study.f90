program vrmhd_pairing_study
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use vrmhd_spectrum_support, only: spectrum, require, spectrum_result_t
    use rm_mhd_oracles, only: uniform_waves, oracle_mu0, oracle_ok
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build
    use kin6d_gs_symmetry, only: gs_symmetry_t, GS_TRANSLATION, GS_CHART_CARTESIAN
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, meridian_mesh_straight
    implicit none
    type(spectrum_result_t) :: result
    type(mesh_2d_t) :: mesh
    type(meridian_mesh_t) :: mm
    type(vrmhd_section_t) :: op
    real(dp) :: exact(3), b0(3), frequencies(3), optical
    integer :: ks, ka, kind, n, level, nb, unmatched, status, unit, repeat
    character(1024) :: argument
    logical :: smoke
    call get_command_argument(1, argument)
    smoke = trim(argument) == 'smoke'
    b0 = sqrt(oracle_mu0)*[0.6_dp, 0.0_dp, 0.8_dp]
    call uniform_waves([1.0_dp, 0.0_dp, 0.0_dp], b0, 1.0_dp, 0.6_dp, &
        5.0_dp/3, exact, status)
    call require(status == oracle_ok, 'independent uniform-wave oracle')
    exact = sqrt(exact)
    if (smoke) then
        call mesh%create_rectangular(2, 2, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
        call meridian_mesh_straight(mesh, 1, mm, status)
        call require(status == 0, 'admission check mesh')
        call vrmhd_section_build(mm, &
            gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
            1, 1, op, status)
        call require(status == 0, 'default equal-degree constructor')
        call require(op%scalar_fs%n0 == 4, 'default P1 has one scalar dof per vertex')
        call vrmhd_section_build(mm, &
            gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
            1, 1, op, status, scalar_degree=1)
        call require(status == 0, 'ordinary constructor admits equal order')
        call vrmhd_section_build(mm, &
            gs_symmetry_t(GS_TRANSLATION, GS_CHART_CARTESIAN, 0.0_dp), &
            1, 1, op, status, scalar_degree=2)
        call require(status /= 0, 'unequal degrees require the study option')
        call spectrum(4, 1, 1, 14, b0, exact, .true., frequencies, optical, &
            unmatched, result)
        call require(result%magnetic >= 1, 'equal-order magnetic checkerboard detected')
        call require(result%real_ok .and. result%count_ok, 'baseline QZ diagnostics')
        call require(result%backward_error < 1.0e-12_dp, 'baseline pencil residual')
        call spectrum(4, 1, 2, 32, b0, exact, .true., frequencies, optical, &
            unmatched, result)
        call require(result%shear_nulls == 2, 'P2 velocity / Nedelec 1 blind shear')
        call require(result%static_velocity >= 2, 'null quotient detects kinetic energy')
        call spectrum(4, 1, 1, 56, b0, exact, .true., frequencies, optical, &
            unmatched, result, perturbed=.true.)
        call require(result%backward_error < 1.0e-12_dp, 'perturbed pencil residual')
        print '(a)', 'pairing measurement smoke: PASS'
        stop
    end if
    if (len_trim(argument) == 0) argument = 'vrmhd-pairings.csv'
    open(newunit=unit, file=trim(argument), status='replace', action='write')
    write(unit,'(a)') 'mesh,cells_per_period,ks,ka,block_cells,dofs_per_cell,block_dofs,' &
        //'finite,zeros,positive,unmatched,magnetic,velocity,mixed,svd_zeros,' &
        //'energetic_nulls,static_magnetic,static_velocity,static_mixed,vz_blind,' &
        //'alfven_error,fast_error,slow_error,' &
        //'omega2_imag_s_2,backward_error,equilibrium_defect,analytic_null_defect,' &
        //'real_ok,count_ok,wave_multiplicity_ok,null_split_ok,seconds'
    do kind = 0, 1
        repeat = kind+1
        do level = 1, 2
            n = 2**(level+1)
            do ks = 1, 4
                do ka = 1, min(4,ks+1)
                    nb = repeat**2*(6*ks**2+2*(3*ka**2+ka))
                    print '(a,4i5)', 'study mesh / n / ks / ka: ',kind,n,ks,ka
                    call spectrum(n,ka,ks,nb,b0,exact,.true.,frequencies,optical, &
                        unmatched,result,perturbed=kind == 1)
                    call require(result%magnetic+result%velocity+result%mixed == unmatched, &
                        'energy classification accounts for every unmatched eigenvalue')
                    write(unit,'(a,19(",",i0),7(",",es24.16),4(",",l1),",",f12.6)') &
                        trim(merge('perturbed ', 'structured',kind == 1)), &
                        n,ks,ka,repeat**2,nb/repeat**2,nb,result%finite,result%zeros, &
                        result%positive,unmatched,result%magnetic,result%velocity, &
                        result%mixed,result%svd_zeros,result%energetic_nulls,result%static_magnetic, &
                        result%static_velocity,result%static_mixed,result%shear_nulls, &
                        abs(frequencies/exact-1),result%reality,result%backward_error, &
                        result%equilibrium,result%analytic_null_defect, &
                        result%real_ok,result%count_ok,all(result%signed_counts == 1), &
                        result%null_split_ok,result%seconds
                    flush(unit)
                end do
            end do
        end do
    end do
    close(unit)
end program
