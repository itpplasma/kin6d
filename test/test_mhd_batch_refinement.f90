program test_mhd_batch_refinement
    ! Transverse uniform-field slab Fourier mode: L xi = B^2 k^2 xi/mu.
    ! For this positive scalar operator the relative solution infinity error
    ! equals the relative residual; the regression gate follows the stop test.
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32
    use kin6d_ideal_xi, only: xi_bilinear_batch
    use fortnum_precision, only: precision_policy_t, precision_pair_t, &
            precision_preset, FP64, FP32, FP16, BF16, STAGE_PRECONDITION, STAGE_RESIDUAL
    use fortnum_precision_vector, only: precision_vector_t
    use fortnum_refinement, only: mixed_refine
    implicit none
    integer, parameter :: n = 257, formats(4) = [FP64, FP32, FP16, BF16]
  real(dp), parameter :: field = 1.1_dp, wave = 1.3_dp, mu = 1.0_dp, tolerance = 1e-6_dp
    real(dp), parameter :: lambda = field**2*wave**2/mu
    real(dp) :: inputs(n, 47), outputs(n, 25), rhs(n), exact(n), x(n), ax(n)
    real(dp) :: low_coefficient(n), wide_coefficient(n), residual, defect, rounding
    real(sp) :: single_inputs(n, 47), single_outputs(n, 25)
    type(precision_policy_t) :: policy
    type(precision_pair_t) :: pair
    type(precision_vector_t) :: storage
    integer :: i, f, info, steps
    logical :: ok, accepted
    inputs = 0; inputs(:, 3) = field; inputs(:, 5) = mu
    inputs(:, 12) = wave; inputs(:, 21) = wave; inputs(:, 47) = 1
    call xi_bilinear_batch(n, inputs, outputs)
    if (maxval(abs(outputs(:, 1) - lambda)) > 2e-12_dp*lambda) &
        error stop 'shared Hessian vs analytical Alfven coefficient'
    single_inputs = real(inputs, sp)
    call xi_bilinear_batch(n, single_inputs, single_outputs)
    do i = 1, n
        exact(i) = sin(acos(-1.0_dp)*real(i, dp)/(n + 1))
    end do
    rhs = lambda*exact
print '(a)', '# storage,steps,relative_residual,relative_state_error,coefficient_defect'
    do f = 1, size(formats)
        policy = precision_preset('cpu-mixed', ok)
        if (.not. ok) error stop 'precision preset'
        call policy%set(STAGE_PRECONDITION, formats(f), &
                        merge(FP64, FP32, formats(f) == FP64), ok)
        if (.not. ok) error stop 'precondition policy'
        wide_coefficient = real(single_outputs(:, 1), dp)
        if (formats(f) == FP64) wide_coefficient = outputs(:, 1)
        call storage%store(wide_coefficient, formats(f), ok)
        if (.not. ok) error stop 'coefficient storage'
        call storage%load(low_coefficient, ok)
        if (.not. ok) error stop 'coefficient load'
        rounding = maxval(abs(low_coefficient/lambda - 1))
        ! The exact scalar conditioning is one; the allowed physical error
        ! is 1%, with the policy's 10% fraction for candidate rounding.
        call policy%assess(STAGE_PRECONDITION, rounding, 1e-2_dp, accepted)
        if (.not. accepted) error stop 'storage requires precision escalation'
        pair = policy%get(STAGE_RESIDUAL)
        if (pair%compute /= FP64 .or. pair%storage /= FP64) error stop 'residual lock'
        x = 0
        call mixed_refine(exact_action, correction, rhs, x, policy, tolerance, &
                          8, info, steps, residual)
        if (info /= 0) error stop 'MHD mixed refinement'
        call exact_action(x, ax)
        residual = maxval(abs(rhs - ax))/maxval(abs(rhs))
        defect = maxval(abs(x - exact))/maxval(abs(exact))
        if (residual > tolerance) error stop 'independent FP64 residual'
   if (defect > tolerance + 2e-12_dp) error stop 'condition-one residual/error identity'
        print '(i0,a,i0,3(a,es24.16))', formats(f), ',', steps, ',', residual, &
            ',', defect, ',', rounding
    end do
contains
    subroutine exact_action(state, action)
        real(dp), intent(in) :: state(:)
        real(dp), intent(out) :: action(:)
        inputs(:, 12) = wave*state
        call xi_bilinear_batch(n, inputs, outputs)
        action = outputs(:, 1)
    end subroutine exact_action
    subroutine correction(defect, delta, precision, status)
        real(dp), intent(in) :: defect(:)
        real(dp), intent(out) :: delta(:)
        type(precision_pair_t), intent(in) :: precision
        integer, intent(out) :: status
        status = 0
        if (precision%compute == FP64) then
            delta = defect/low_coefficient
        else
            delta = real(real(defect, sp)/real(low_coefficient, sp), dp)
        end if
    end subroutine correction
end program test_mhd_batch_refinement
