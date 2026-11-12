program test_vrmhd_torus
    !! Circular-torus flux-pumping case (docs/vrmhd-torus.md): axisymmetric
    !! Grad-Shafranov state of the cylinder profiles.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_vrmhd_cylinder, only: vrmhd_cylinder_case_t, vrmhd_cylinder_profile_t, &
        vrmhd_cylinder_profiles
    use milca_vrmhd_torus, only: vrmhd_torus_law, vrmhd_torus_equilibrium, &
        vrmhd_torus_axis_t
    use milca_gs_profile_law, only: gs_profile_law_t
    use milca_grad_shafranov, only: gs_state_t
    use rm_mhd_oracles, only: shafranov_formula
    implicit none
    character(64) :: gate, arg
    real(dp) :: edge = 0.1_dp
    integer :: knots = 129, k, e, failures = 0
    call get_command_argument(1, gate)
    do k = 2, command_argument_count()
        call get_command_argument(k, arg)
        e = index(arg, '=')
        if (e == 0) error stop 'arguments are key=value'
        select case (arg(:e - 1))
        case ('edge'); read (arg(e + 1:), *) edge
        case ('knots'); read (arg(e + 1:), *) knots
        case default; error stop 'unknown argument'
        end select
    end do
    select case (trim(gate))
    case ('shafranov')
        call shafranov()
    case default
        error stop 'expected shafranov'
    end select
    if (failures /= 0) error stop 'torus checks failed'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        if (condition) then
            print '(2a)', 'PASS ', description
        else
            print '(2a)', 'FAIL ', description
            failures = failures + 1
        end if
    end subroutine

    ! Axis shift against the large-aspect-ratio formula where the shifted
    ! circles stay well apart (max |dDelta/dr| <= 0.3), and elsewhere against
    ! the polar finite-difference Grad-Shafranov solve of docs/vrmhd-torus.md.
    ! The unstructured mesh is not up-down symmetric; its axis height is a
    ! discretisation error of about 1e-4 at edge = 0.2.
    subroutine shafranov()
        type(vrmhd_cylinder_case_t) :: c
        type(vrmhd_cylinder_profile_t) :: prof
        type(gs_profile_law_t) :: law
        type(gs_state_t) :: state
        type(vrmhd_torus_axis_t) :: axis
        real(dp) :: formula, slope, reference, ratio
        real(dp) :: cases(3, 4)
        integer :: i, status
        ! R0/a, beta_p, Delta(0) of the finite-difference solve
        cases = reshape([3.0_dp, 0.02_dp, 0.0755_dp, 3.0_dp, 1.0_dp, 0.157_dp, &
            3.0_dp, 3.0_dp, 0.289_dp, 10.0_dp, 1.0_dp, 0.0476_dp], [3, 4])
        do i = 1, size(cases, 2)
            c%r0 = cases(1, i); c%beta_poloidal = cases(2, i)
            call vrmhd_cylinder_profiles(c, prof, status)
            call check(status == 0, 'cylinder profiles')
            call vrmhd_torus_law(c, prof, knots, law, status)
            call check(status == 0, 'torus profile law')
            call vrmhd_torus_equilibrium(c, prof, edge, law, state, axis, status)
            call check(status == 0, 'torus Grad-Shafranov solve')
            if (status /= 0) cycle
            call shafranov_formula(prof%x, prof%t, prof%x*prof%btheta_over_r, &
                c%r0, formula, slope)
            reference = cases(3, i)
            if (slope <= 0.3_dp) reference = formula
            ratio = axis%shift/reference
            print '(a,f5.1,a,f5.2,a,4f9.5,a,f8.5,a,es10.2)', 'R0/a', c%r0, ' beta_p', &
                c%beta_poloidal, ' shift, formula, FD, ratio', axis%shift, &
                formula, cases(3, i), ratio, ' q0', axis%q, ' Z', axis%position(2)
            call check(abs(axis%position(2)) < 1.0e-3_dp, 'axis on the midplane')
            if (slope <= 0.3_dp) then
                call check(abs(ratio - 1) < 0.03_dp, 'axis shift within 3 % of the formula')
            else
                call check(abs(ratio - 1) < 0.01_dp, &
                    'axis shift within 1 % of the finite-difference solve')
            end if
        end do
    end subroutine
end program test_vrmhd_torus
