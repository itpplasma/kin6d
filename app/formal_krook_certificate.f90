program formal_krook_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_support_datatype
    use kin6d_poisson_leaves, only: coupled_amplitude
    implicit none
    real(dp) :: amplitude
    integer(int64) :: bits, significand
    integer :: biased_exponent, power, unit
    character(512) :: path
    character(64) :: numerator, denominator
    character(:), allocatable :: rational

    ! Prescribed E=0, nu=1, exact source amplitude=1/2. This scalar direct
    ! solve supplies the coefficient of the continuous polynomial/cosine
    ! reconstruction defined and differentiated by the Lean theorem.
    call coupled_amplitude(0.5_dp, 1.0_dp, 0.0_dp, 1.0_dp, amplitude)
    if (.not. ieee_is_finite(amplitude)) error stop 'nonfinite formal candidate'
    if (amplitude < 0.0_dp) error stop 'unsupported negative certificate coefficient'
    if (storage_size(amplitude) /= 64) error stop 'certificate requires IEEE binary64'
    if (.not. ieee_support_datatype(amplitude)) error stop 'IEEE data unsupported'
    if (radix(amplitude) /= 2) error stop 'certificate requires binary radix'
    if (digits(amplitude) /= 53) error stop 'certificate requires binary64 precision'
    if (maxexponent(amplitude) /= 1024) error stop 'unsupported exponent range'
    if (minexponent(amplitude) /= -1021) error stop 'unsupported exponent range'
    if (bit_size(0_int64) /= 64) error stop 'certificate requires 64-bit integers'
    if (transfer(1.0_dp, 0_int64) /= int(z'3ff0000000000000', int64)) &
        error stop 'unsupported binary64 bit interpretation'
    ! Serialization interprets the actual candidate bytes as an exact dyadic.
    ! These are integer representation operations, not numerical physics.
    bits = transfer(amplitude, 0_int64)
    biased_exponent = int(ibits(bits, 52, 11))
    significand = ibits(bits, 0, 52)
    if (biased_exponent == 0) then
        power = -1074
    else
        significand = ibset(significand, 52)
        power = biased_exponent - 1023 - 52
    end if
    write (numerator, '(i0)') significand
    write (denominator, '(i0)') abs(power)
    if (power < 0) then
        rational = '((' // trim(numerator) // ' : Rat) / 2 ^ ' // &
            trim(denominator) // ')'
    else
        rational = '((' // trim(numerator) // ' : Rat) * 2 ^ ' // &
            trim(denominator) // ')'
    end if
    path = 'output/ComputedKrook.lean'
    if (command_argument_count() > 0) call get_command_argument(1, path)
    open (newunit=unit, file=trim(path), status='replace')
    write (unit, '(a)') 'import Kin6d'
    write (unit, '(a)') '-- Exact dyadic interpretation of the Fortran candidate bytes.'
    write (unit, '(a)') 'namespace Kin6d.Computed'
    write (unit, '(a)') 'def candidate : Certificate := '
    write (unit, '(a)') '  { exactAmplitude := 1/2, candidateAmplitude := ' // &
        rational // ', residualBound := 0, radius := 0 }'
    write (unit, '(a)') 'theorem accepted : candidate.Valid := by'
    write (unit, '(a)') '  norm_num [candidate, Certificate.Valid]'
    write (unit, '(a)') 'theorem continuous : '
    write (unit, '(a)') '  ∀ v ∈ Set.Icc (-1 : Real) 1, ∀ x : Real,'
    write (unit, '(a)') '    abs (profile (candidate.candidateAmplitude : Real) v x -'
    write (unit, '(a)') '      profile (candidate.exactAmplitude : Real) v x) ≤'
    write (unit, '(a)') '      (candidate.radius : Real) := candidate.sound accepted'
    write (unit, '(a)') '#print axioms continuous'
    write (unit, '(a)') 'end Kin6d.Computed'
    close (unit)
    print '(a,es24.16)', 'Formal Krook candidate amplitude=', amplitude
    print '(a)', 'Wrote exact dyadic certificate: ' // trim(path)
end program
