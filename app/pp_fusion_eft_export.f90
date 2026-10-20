! Minimal-input pionless EFT; JSON interchange plus native FortIO provenance.
program pp_fusion_eft_export
    use, intrinsic :: iso_fortran_env, only: dp => real64, compiler_version, &
        compiler_options
    use kin6d_pp_fusion_eft, only: pp_eft_config_t, pp_eft_result_t, pp_eft_evaluate, &
        pp_eft_lattice_errors, pp_eft_input_errors
    use kin6d_pp_fusion, only: pp_fusion_config_t, pp_fusion_result_t, &
        pp_fusion_solve, &
        pp_vector_coupling, pp_hc, pp_mp, pp_mn, pp_alpha, axial, ft0
    use kin6d_result_io, only: run_t, result_file_t, runtime_git
    implicit none
    include 'kin6d_build_info.inc'
    type(pp_eft_config_t) :: cfg, zero
    type(pp_eft_result_t) :: answer(2)
    type(pp_fusion_config_t) :: reid_cfg
    type(pp_fusion_result_t) :: reid
    type(run_t) :: run
    type(result_file_t) :: result
    real(dp) :: binding, a_pp, r_pp, rho_d, l1a, mu_mev, jac(8, 5, 2)
    real(dp) :: input_values(6), error_values(4)
    character(24), parameter :: order_names(4) = &
        [character(24) :: 'LO', 'NLO', 'NNLO', 'range_resummed_NNLO']
    character(24), parameter :: input_names(6) = &
        [character(24) :: 'B_d', 'a_pp', 'r_pp', 'rho_d', 'L1A', 'mu']
    character(8), parameter :: input_units(6) = &
        [character(8) :: 'MeV', 'fm', 'fm', 'fm', 'fm^3', 'MeV']
    character(64), parameter :: input_sources(6) = [character(64) :: &
        '10.1103/PhysRevC.51.38 Table X', &
        '10.1103/PhysRevC.38.15; 10.1103/PhysRevC.51.38 Table VIII', &
        '10.1103/PhysRevC.38.15; 10.1103/PhysRevC.51.38 Table VIII', &
        '10.1016/S0370-2693(01)01152-2 (7); 1610.04545v2 supplement', &
        '10.1103/PhysRevLett.119.062002 (14), arXiv v2', &
        '10.1016/j.physletb.2013.02.019 (11)']
    character(48), parameter :: input_kinds(6) = [character(48) :: &
        'measured binding energy', 'fit to Coulomb-modified scattering', &
        'fit to Coulomb-modified scattering', 'fit to deuteron pole residue', &
        'lattice matching with quark-mass assumption', &
        'PDS convention; also EFT breakdown scale']
    character(32), parameter :: lerrors(4) = [character(32) :: &
        'statistical', 'fitting_and_analysis', 'mass_extrapolation', 'EFT_matching']
    character(4096) :: arg, json_path, message
    character(:), allocatable :: commit, tree
    integer :: unit, json_unit, ios, i, j, k, status
    logical :: ok
    namelist /pp_fusion_eft/ binding, a_pp, r_pp, rho_d, l1a, mu_mev
    binding = cfg%binding; a_pp = cfg%a_pp; r_pp = cfg%r_pp
    rho_d = cfg%rho_d; l1a = cfg%l1a; mu_mev = cfg%mu_mev
    json_path = 'output/pp_fusion_eft.json'
    call run%start('pp_fusion_eft_export')
    call run%open_config(unit)
    ios = 0; message = ''
    if (run%has_config()) read (unit, nml=pp_fusion_eft, iostat=ios, iomsg=message)
    call run%close_config(unit, 'pp_fusion_eft', ios, message)
    i = 1
    do while (i <= run%argument_count())
        call run%argument(i, arg)
        if (trim(arg) /= '--json') &
            error stop 'usage: pp_fusion_eft_export [--config FILE] [--json FILE]'
        i = i + 1
        if (i > run%argument_count()) error stop '--json needs a path'
        call run%argument(i, json_path)
        i = i + 1
    end do
    cfg%binding = binding; cfg%a_pp = a_pp; cfg%r_pp = r_pp
    cfg%rho_d = rho_d; cfg%l1a = l1a; cfg%mu_mev = mu_mev
    call run%open_resolved(unit)
    write (unit, nml=pp_fusion_eft)
    call run%close_resolved(unit)
    zero = cfg; zero%l1a = 0
    call pp_eft_evaluate(zero, answer(1), ok, jac(:, :, 1))
    if (.not. ok) error stop 'invalid EFT input or failed Coulomb quadrature'
    call pp_eft_evaluate(cfg, answer(2), ok, jac(:, :, 2))
    if (.not. ok) error stop 'invalid EFT input or failed Coulomb quadrature'
    call pp_fusion_solve(reid_cfg, reid, ok)
    if (.not. ok) error stop 'Reid68 comparison failed'
    call execute_command_line('mkdir -p output', exitstat=status)
    if (status /= 0) error stop 'cannot create output/'
    call result%create(run)
    call result%attribute('physical_model', 'pionless EFT; point Coulomb; threshold GT')
    call result%attribute('amplitude_source', &
        'Butler-Chen (7),(9),(10); Chen-Liu-Yu (12)')
    call result%attribute('numerical_bound', answer(2)%numerical_bound)
    call result%attribute('model_bound', answer(2)%model_bound)
    call result%attribute('lattice_doi', '10.1103/PhysRevLett.119.062002')
    call result%write('lambda2_no_L1A', answer(1)%lambda2, '1', &
        'LO,NLO,NNLO,range-resummed diagnostic; L1A(mu)=0', status='candidate')
    call result%write('lambda2_with_L1A', answer(2)%lambda2, '1', &
        'LO,NLO,NNLO,range-resummed diagnostic', status='candidate')
    call result%write('S_no_L1A', answer(1)%sfactor, 'MeV b', &
        'zero-energy S factor; L1A(mu)=0', status='candidate')
    call result%write('S_with_L1A', answer(2)%sfactor, 'MeV b', &
        'zero-energy S factor', status='candidate')
    call result%write('S_prime_over_S', answer(2)%slope, 'MeV^-1', &
        'threshold slope; uncertainty unavailable', status='candidate')
    call result%write('S_truncation_estimate', answer(2)%truncation_s, 'MeV b', &
        'power-counting estimate; neither probability nor certificate', &
        status='indicator')
    call result%write('input_jacobian', jac(:, :, 2), 'mixed', &
        'rows lambda2 then S; columns B,a,r,rho,L1A', status='candidate')
    call result%write('Reid68_S_0', reid%sfactor, 'MeV b', &
        'independently solved Reid68 one-body comparison', status='candidate')
    call result%close()
    call runtime_git(commit, tree, 'pp_fusion_eft_export')
    open (newunit=json_unit, file=trim(json_path), status='replace', iostat=ios)
    if (ios /= 0) error stop 'cannot open JSON path; create its directory first'
    write (json_unit, '(a)') '{'
    call string('schema', 'kin6d-pp-fusion-eft-v1')
    call string('claim_status', 'candidate')
    call string('rate_expansion', &
        'Lambda is the truncated amplitude; Lambda2 is the consistently expanded rate')
    call string('verification', &
        'verified implementation; model assumptions remain hypotheses')
    call string('numerical_bound', answer(2)%numerical_bound)
    call string('model_bound', answer(2)%model_bound)
    call string('input_uncertainty', &
        'linear sensitivities; separate source components, no joint law')
    call string('truncation_semantics', &
        'power-counting estimate, not a rigorous bound or confidence interval')
    call string('counterterm_convention', &
        'Butler-Chen PDS L1A in fm^3; gA divided in invariant Lbar')
    call string('zero_counterterm', &
        'L1A(mu)=0; invariant Lbar is nonzero; scale-dependent diagnostic')
    call string('lattice_assumption', &
        'm_pi~806 MeV; mild mass dependence; final arXiv v2 uncertainties')
    call string('resummation_status', &
        'diagnostic, not a new EFT order; omits K1A and shape terms')
    call string('source_discrepancy', &
        'resummed BC matching differs from Savage rounded 2.62+0.0105 L1A; see'// &
        ' documentation')
    call string('S_prime_over_S_status', &
        'candidate; Chen-Liu-Yu (7)-(11) threshold derivative; uncertainty unavailable')
    call string('slope_input_derivatives', &
        'unavailable in the high-level input Jacobian; generated leaf products'// &
        ' retained')
    call string('weak_normalization_doi', '10.1103/RevModPhys.83.195 Eqs (20)-(22)')
    call string('constants_doi', '10.1103/RevModPhys.97.025002')
    call string('amplitude_doi', '10.1016/S0370-2693(01)01152-2 Eqs (7),(9),(10)')
    call string('rate_doi', '10.1016/j.physletb.2013.02.019 Eq (12)')
    call string('git_commit', commit)
    call string('git_tree', tree)
    call string('fortnum_revision', kin6d_fortnum_revision)
    call string('fortsym_revision', kin6d_fortsym_revision)
    call string('compiler', compiler_version())
    call string('compiler_options', compiler_options())
    call string('command', run%command_line)
    call string('literal_configuration', run%literal)
    call string('resolved_configuration', run%resolved)
    call scalar('G_V_MeV^-2', pp_vector_coupling())
    call scalar('g_A', axial)
    call scalar('ft_s', ft0)
    call scalar('alpha', pp_alpha)
    call scalar('hbar_c_MeV_fm', pp_hc)
    call scalar('m_p_MeV', pp_mp)
    call scalar('m_n_MeV', pp_mn)
    call scalar('gamma_per_fm', answer(2)%gamma)
    call scalar('chi', answer(2)%chi)
    call scalar('chi_I_chi', answer(2)%coulomb_integral)
    call scalar('f_pp_R', answer(2)%phase_space)
    call scalar('Q', answer(2)%expansion_q)
    call scalar('coulomb_quadrature_indicator', answer(2)%quadrature_indicator)
    call scalar('Reid68_S_prime_over_S_per_MeV', reid%slope)
    call scalar('Reid68_Lambda2', reid%lambda2)
    call scalar('Reid68_S_pp_0_MeV_b', reid%sfactor)
    call string('Reid68_uncertainty', 'unavailable; one-body potential-model candidate')
    call scalar('Adelberger2011_S_pp_0_MeV_b', 4.01e-25_dp)
    call scalar('Adelberger2011_uncertainty_MeV_b', .04e-25_dp)
    input_values = [binding, a_pp, r_pp, rho_d, l1a, mu_mev]
    error_values = pp_eft_input_errors
    write (json_unit, '(a)') '"inputs":['
    do i = 1, 6
        write (json_unit, '(a)') '{'
        call string('name', trim(input_names(i)))
        call string('unit', trim(input_units(i)))
        call string('source', trim(input_sources(i)))
        call string('kind', trim(input_kinds(i)))
        call scalar('value', input_values(i))
        if (i == 5) then
            do k = 1, 4
                call scalar(trim(lerrors(k))//'_error', pp_eft_lattice_errors(k))
            end do
        end if
        if (i <= 4) then
            call scalar('quoted_error', error_values(i))
        else
            write (json_unit, '(a)') '"quoted_error":null,'
        end if
        write (json_unit, &
            '(a)') '"error_status":"L1A components above; mu fixed; other erro'// &
            'rs quoted source values"}'
        if (i < 6) write (json_unit, '(a)') ','
    end do
    write (json_unit, '(a)') '],"orders":['
    do j = 1, 2
        do i = 1, 4
            write (json_unit, '(a)') '{'
            call string('order', trim(order_names(i)))
            if (j == 1) then
                call string('scenario', 'L1A_zero')
                call scalar('L1A_fm3', 0._dp)
            else
                call string('scenario', 'supplied_L1A')
                call scalar('L1A_fm3', l1a)
            end if
            call scalar('Lambda', answer(j)%lambda(i))
            call scalar('Lambda2', answer(j)%lambda2(i))
            call scalar('S_pp_0_MeV_b', answer(j)%sfactor(i))
            call scalar('S_prime_over_S_per_MeV', answer(j)%slope(i))
            call scalar('truncation_Lambda_estimate', answer(j)%truncation_lambda(i))
            call scalar('truncation_Lambda2_estimate', answer(j)%truncation_lambda2(i))
            call scalar('truncation_S_estimate_MeV_b', answer(j)%truncation_s(i))
            call scalar('Lbar_fm', answer(j)%lbar)
            write (json_unit, '(a)') '"input_error_components":{'
            do k = 1, 4
                call scalar(trim(input_names(k))//'_Lambda2', abs(jac(i, k, &
                    j))*pp_eft_input_errors(k))
                call scalar(trim(input_names(k))//'_S_MeV_b', abs(jac(i + 4, k, &
                    j))*pp_eft_input_errors(k))
            end do
            if (j == 2) then
                do k = 1, 4
                    call scalar('L1A_'//trim(lerrors(k))//'_Lambda2', &
                        abs(jac(i, 5, j))*pp_eft_lattice_errors(k))
                    call scalar('L1A_'//trim(lerrors(k))//'_S_MeV_b', &
                        abs(jac(i + 4, 5, j))*pp_eft_lattice_errors(k))
                end do
            end if
            write (json_unit, &
                '(a)') '"combination_status":"not combined; matching and trunc'// &
                'ation may overlap"}}'
            if (i < 4 .or. j < 2) write (json_unit, '(a)') ','
        end do
    end do
    write (json_unit, '(a)') ']}'
    close (json_unit)
    print '(a)', &
        'Order                         Lambda2        S / (1e-25 MeV b)  trunc'// &
        'ation estimate'
    do i = 1, 4
        print '(a24,3f17.8)', order_names(i), answer(2)%lambda2(i), &
            answer(2)%sfactor(i)/1e-25_dp, answer(2)%truncation_s(i)/1e-25_dp
    end do
    print '(a)', trim(json_path)
contains
    subroutine scalar(key, value)
        character(*), intent(in) :: key
        real(dp), intent(in) :: value
        write (json_unit, '(a,es25.16e3,a)') '"'//key//'":', value, ','
    end subroutine
    subroutine string(key, value)
        character(*), intent(in) :: key, value
        integer :: n, c
        write (json_unit, '(a)', advance='no') '"'//key//'":"'
        do n = 1, len(value)
            c = iachar(value(n:n))
            select case (c)
            case (10)
                write (json_unit, '(a)', advance='no') '\n'
            case (13)
                write (json_unit, '(a)', advance='no') '\r'
            case (9)
                write (json_unit, '(a)', advance='no') '\t'
            case (34, 92)
                write (json_unit, '(a)', advance='no') achar(92)//value(n:n)
            case default
                if (c < 32) error stop 'unsupported JSON control character'
                write (json_unit, '(a)', advance='no') value(n:n)
            end select
        end do
        write (json_unit, '(a)') '",'
    end subroutine
end program
