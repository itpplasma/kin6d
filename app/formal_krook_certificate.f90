program formal_krook_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_support_datatype
    use kin6d_poisson_leaves, only: coupled_amplitude
    use kin6d_result_io, only: run_t, result_file_t, result_directory
    implicit none
    ! Inputs of the scalar direct solve. They are fixed by the Lean model
    ! (exactAmplitude := 1/2); the namelist records them and rejects changes.
    real(dp), parameter :: lean_source = 0.5_dp, lean_nu = 1.0_dp
    real(dp), parameter :: lean_field = 0.0_dp, lean_temperature = 1.0_dp
    real(dp) :: source_amplitude = lean_source, nu = lean_nu
    real(dp) :: electric_field = lean_field, temperature = lean_temperature
    namelist /formal_krook/ source_amplitude, nu, electric_field, temperature
    type(run_t) :: run
    type(result_file_t) :: result
    integer :: ios
    character(256) :: msg
    real(dp) :: amplitude
    integer(int64) :: bits, significand
    integer :: biased_exponent, power, unit
    character(512) :: path
    character(64) :: numerator, denominator
    character(:), allocatable :: rational

    call run%start('formal_krook_certificate')
    if (run%argument_count() > 1) error stop &
        'usage: formal_krook_certificate [--config FILE] [lean_path]'
    ios = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=formal_krook, iostat=ios, iomsg=msg)
    call run%close_config(unit, 'formal_krook', ios, msg)
    if (source_amplitude /= lean_source .or. nu /= lean_nu .or. &
        electric_field /= lean_field .or. temperature /= lean_temperature) &
        error stop '&formal_krook values are fixed by the Lean model '// &
        '(source 1/2, nu 1, E 0, temperature 1)'
    call run%open_resolved(unit)
    write (unit, nml=formal_krook)
    call run%close_resolved(unit)
    ! Prescribed E=0, nu=1, exact source amplitude=1/2. This scalar direct
    ! solve supplies the coefficient of the continuous polynomial/cosine
    ! reconstruction defined and differentiated by the Lean theorem.
    call coupled_amplitude(source_amplitude, nu, electric_field, temperature, &
        amplitude)
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
    if (run%argument_count() > 0) call run%argument(1, path)
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
    call execute_command_line('mkdir -p '//result_directory, exitstat=ios)
    if (ios /= 0) error stop 'cannot create result directory'
    call result%create(run)
    call result%attribute('lean_certificate_path', trim(path))
    call result%write('candidate_amplitude', amplitude, '1', &
        'binary64 candidate amplitude serialized as an exact dyadic', &
        status='candidate; formally-certified-continuous only after the Lean '// &
        'check (proof/check.sh) accepts '//trim(path))
    call result%write('exact_amplitude', 0.5_dp, '1', &
        'exact amplitude of the Lean model', status='exact model value')
    call result%close()
    print '(a,es24.16)', 'Formal Krook candidate amplitude=', amplitude
    print '(a)', 'Wrote exact dyadic certificate: ' // trim(path)
end program
