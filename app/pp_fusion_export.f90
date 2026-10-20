! Teaching export; canonical scientific result also written through FortIO.
program pp_fusion_export
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version, &
        compiler_options
    use kin6d_pp_fusion, only: pp_fusion_config_t, pp_fusion_result_t, &
        pp_fusion_solve, &
        pp_hc, pp_mp, pp_alpha, pp_vector_coupling
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    use fortnum_special_coulomb, only: coulomb_c0_squared
    implicit none
    include 'kin6d_build_info.inc'
    type(pp_fusion_config_t) :: cfg
    type(pp_fusion_result_t) :: answer
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp) :: rmax, hmax, h0, energy, strength(2), k, eta, c2, jac(5, 2)
    integer :: order, unit, ios, i, status, json_unit
    logical :: ok
    character(4096) :: arg, json_path, message
    character(:), allocatable :: commit, tree
    namelist /pp_fusion/ rmax, hmax, h0, energy, order, strength
    rmax = cfg%rmax; hmax = cfg%hmax; h0 = cfg%h0; energy = cfg%energy
    order = cfg%order; strength = cfg%strength
    json_path = 'output/pp_fusion_kin6d.json'
    call run%start('pp_fusion_export')
    call run%open_config(unit)
    ios = 0; message = ''
    if (run%has_config()) read (unit, nml=pp_fusion, iostat=ios, iomsg=message)
    call run%close_config(unit, 'pp_fusion', ios, message)
    i = 1
    do while (i <= run%argument_count())
        call run%argument(i, arg)
        if (trim(arg) /= '--json') &
            error stop 'usage: pp_fusion_export [--config FILE] [--json FILE]'
        i = i + 1
        if (i > run%argument_count()) error stop '--json needs a path'
        call run%argument(i, json_path)
        i = i + 1
    end do
    cfg%rmax = rmax; cfg%hmax = hmax; cfg%h0 = h0; cfg%energy = energy
    cfg%order = order; cfg%strength = strength
    call run%open_resolved(unit)
    write (unit, nml=pp_fusion)
    call run%close_resolved(unit)
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output/'
    call pp_fusion_solve(cfg, answer, ok, jac)
    if (.not. ok) error stop 'pp fusion model or numerical solve failed'
    call result%create(run)
    call result%attribute('physical_model', 'Reid68 + point Coulomb; one-body GT')
    call result%attribute('potential_doi', '10.1016/0003-4916(68)90126-7')
    call result%attribute('normalization_doi', '10.1103/PhysRevC.58.1263')
    call result%attribute('numerical_bound', answer%numerical_bound)
    call result%attribute('model_bound', answer%model_bound)
    call result%write('r', answer%mesh%r, 'fm', 'radial GLL nodes', status='candidate')
    call result%write('u_d_S', answer%us, 'fm^-1/2', 'deuteron S wave', &
        status='candidate')
    call result%write('u_d_D', answer%ud, 'fm^-1/2', 'deuteron D wave', &
        status='candidate')
    call result%write('u_pp', answer%upp, '1', &
        'unit Coulomb amplitude scattering wave', status='candidate')
    call result%write('lambda_integrand', answer%integrand, 'fm^-1', &
        'Lambda integrand', status='candidate')
    call result%write('binding', answer%binding, 'MeV', 'deuteron binding energy', &
        'candidate')
    call result%write('D_probability', answer%d_probability, '1', &
        'D state probability', status='candidate')
    call result%write('a_pp', answer%a_pp, 'fm', &
        'Coulomb-modified scattering length', status='candidate')
    call result%write('lambda2', answer%lambda2, '1', 'one-body overlap squared', &
        'candidate')
    call result%write('S_pp_0', answer%sfactor, 'MeV b', &
        'one-body S factor at zero energy', status='candidate')
    call result%write('S_prime_over_S', answer%slope, 'MeV^-1', &
        'zero energy slope', status='candidate')
    call result%write('phase_space', answer%phase_space, '1', &
        'computed phase space', status='candidate')
    call result%write('potentials', answer%potentials, 'MeV', &
        'pp total; SS; SD; DD nuclear', status='candidate')
    call result%write('strength_jacobian', jac, 'mixed', &
        'rows B,PD,a,lambda2,S; cols singlet,triplet', status='candidate')
    call result%close()
    call runtime_git(commit, tree, 'pp_fusion_export')
    open (newunit=json_unit, file=trim(json_path), status='replace', iostat=ios)
    if (ios /= 0) &
        error stop 'cannot open requested JSON path (create its directory first)'
    write (json_unit, '(a)') '{'
    call string('schema', 'kin6d-pp-fusion-v1')
    call string('claim_status', trim(answer%status))
    call string('numerical_bound', trim(answer%numerical_bound))
    call string('model_bound', trim(answer%model_bound))
    call string('pointwise_bound', 'unavailable')
    call string('physical_model', 'Reid68 + point Coulomb; one-body Gamow-Teller')
    call string('potential_doi', '10.1016/0003-4916(68)90126-7')
    call string('overlap_doi', '10.1103/PhysRevC.58.1263')
    call string('weak_normalization_doi', '10.1103/RevModPhys.83.195')
    call string('constants_doi', '10.1103/RevModPhys.97.025002')
    call string('git_commit', commit)
    call string('git_tree', tree)
    call string('fortnum_revision', kin6d_fortnum_revision)
    call string('fortsym_revision', kin6d_fortsym_revision)
    call string('compiler', compiler_version())
    call string('compiler_options', compiler_options())
    call string('command', run%command_line)
    call string('literal_configuration', run%literal)
    call string('resolved_configuration', run%resolved)
    call scalar('B_d_MeV', answer%binding)
    call scalar('P_D', answer%d_probability)
    call scalar('a_pp_fm', answer%a_pp)
    call scalar('Lambda2_0', answer%lambda2)
    call scalar('S_pp_0_MeV_b', answer%sfactor)
    call scalar('S_prime_over_S_per_MeV', answer%slope)
    call scalar('f_pp_R', answer%phase_space)
    call scalar('energy_MeV', answer%energy)
    call scalar('phase_shift_radians', answer%phase_shift)
    call scalar('eigen_residual_MeV', answer%eigen_residual)
    call scalar('scattering_residual_MeV', answer%scattering_residual)
    call scalar('Lambda2_energy_extrapolation_indicator', answer%energy_extrapolation)
    k = sqrt(energy*pp_mp/pp_hc**2)
    eta = pp_alpha*pp_hc/(2*(pp_hc**2/pp_mp)*k)
    c2 = coulomb_c0_squared(eta)
    call scalar('k_per_fm', k)
    call scalar('G_V_MeV^-2', pp_vector_coupling())
    call scalar('eta', eta)
    call scalar('C0', sqrt(c2))
    write (json_unit, '(a)') '"radial_columns": ["r_fm", "u_d_S_fm^-1/2", '// &
        '"u_d_D_fm^-1/2", "u_pp", "u_pp_over_Ck_fm", "V_pp_MeV", '// &
        '"V_d_SS_MeV", "V_d_SD_MeV", "V_d_DD_nuclear_MeV", '// &
        '"Lambda_integrand_per_fm", "quadrature_weight_fm"],'
    write (json_unit, '(a)') '"radial": ['
    do i = 1, answer%mesh%n
        write (json_unit, '(a,10(es25.16e3,a),es25.16e3,a)', advance='no') &
            '[', answer%mesh%r(i), ',', answer%us(i), ',', answer%ud(i), ',', &
            answer%upp(i), ',', answer%upp(i)/(sqrt(c2)*k), ',', &
            answer%potentials(1, i), ',', answer%potentials(2, i), ',', &
            answer%potentials(3, i), ',', answer%potentials(4, i), ',', &
            answer%integrand(i), ',', answer%mesh%w(i), ']'
        if (i < answer%mesh%n) then
            write (json_unit, '(a)') ','
        else
            write (json_unit, '(a)') ''
        end if
    end do
    write (json_unit, '(a)') ']}'
    close (json_unit)
    print '(a,es16.8,a)', 'S_pp(0) = ', answer%sfactor, ' MeV b (one-body candidate)'
    print '(a)', trim(json_path)
contains
    subroutine scalar(key, value)
        character(*), intent(in) :: key
        real(dp), intent(in) :: value
        write (json_unit, '(a,es25.16e3,a)') '"'//key//'":', value, ','
    end subroutine
    subroutine string(key, value)
        character(*), intent(in) :: key, value
        integer :: j, c
        write (json_unit, '(a)', advance='no') '"'//key//'":"'
        do j = 1, len(value)
            c = iachar(value(j:j))
            select case (c)
            case (10)
                write (json_unit, '(a)', advance='no') '\n'
            case (13)
                write (json_unit, '(a)', advance='no') '\r'
            case (9)
                write (json_unit, '(a)', advance='no') '\t'
            case (34, 92)
                write (json_unit, '(a)', advance='no') achar(92)//value(j:j)
            case default
                if (c < 32) error stop 'unsupported JSON control character'
                write (json_unit, '(a)', advance='no') value(j:j)
            end select
        end do
        write (json_unit, '(a)') '",'
    end subroutine
end program
