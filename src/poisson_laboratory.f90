module poisson_laboratory
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, mag
    use kin6d_poisson_leaves, only: poisson_coefficients, poisson_reconstruct, &
                                   poisson_coefficients_interval, &
                                   poisson_energy_interval, coupled_denominator, &
                                   coupled_amplitude, coupled_amplitude_interval, &
                                   coupled_density, coupled_distribution, &
                                   coupled_reference_value, coupled_source, &
                                   coupled_residual, coupled_radius_interval, &
                                   difference_interval, sum_interval
    implicit none
    private
    public :: poisson_mode_t, coupled_steady_t, solve_poisson_mode, poisson_sample, &
              solve_coupled_steady, coupled_sample, coupled_reference, &
              coupled_certificate_radius, coupled_certificate_accepts

    type :: poisson_mode_t
        integer :: k = 0
        real(dp) :: epsilon0 = 1.0_dp
        real(dp) :: rho_mean = 0.0_dp, rho_cos = 0.0_dp, rho_sin = 0.0_dp
        real(dp) :: phi_mean = 0.0_dp, phi_cos = 0.0_dp, phi_sin = 0.0_dp
        real(dp) :: e_cos = 0.0_dp, e_sin = 0.0_dp
    end type

    type :: coupled_steady_t
        integer :: k = 0
        real(dp) :: amplitude = 0.0_dp, source_drive = 0.0_dp, nu = 0.0_dp
        real(dp) :: fp_rate = 0.0_dp, temperature = 1.0_dp
        type(poisson_mode_t) :: field
    end type
contains
    subroutine solve_poisson_mode(k, epsilon0, rho_mean, rho_cos, rho_sin, field, ok)
        integer, intent(in) :: k
        real(dp), intent(in) :: epsilon0, rho_mean, rho_cos, rho_sin
        type(poisson_mode_t), intent(out) :: field
        logical, intent(out) :: ok
        ok = .false.
        if (k <= 0) return
        if (any(.not. ieee_is_finite([epsilon0, rho_mean, rho_cos, rho_sin]))) return
        if (epsilon0 <= 0.0_dp) return
        ! Exact compatibility: never project away incompatible net charge.
        if (abs(rho_mean) > 0.0_dp) return
        field%k = k; field%epsilon0 = epsilon0
        field%rho_mean = rho_mean; field%rho_cos = rho_cos; field%rho_sin = rho_sin
        field%phi_mean = 0.0_dp
        call poisson_coefficients(real(k, dp), epsilon0, rho_cos, rho_sin, &
                                  field%phi_cos, field%phi_sin, &
                                  field%e_cos, field%e_sin)
        ok = valid_field(field)
    end subroutine

    subroutine coupled_certificate_radius(state, f_radius, phi_radius, &
                                          e_radius, energy, ok)
        ! Known continuous source-balanced solution gives existence directly.
        ! This certifies proximity to that solution, not general PDE uniqueness.
        ! Uniform pointwise radii cover the entire periodic slab and [-1,1].
        type(coupled_steady_t), intent(in) :: state
        real(dp), intent(out) :: f_radius, phi_radius, e_radius
        type(interval_t), intent(out) :: energy
        logical, intent(out) :: ok
        type(interval_t) :: exact_a, pc, ps, ec, es, delta, pa, pb, ea, eb
        type(interval_t) :: aradius, pradius, eradius, fr, pr, er
        f_radius = 0.0_dp; phi_radius = 0.0_dp; e_radius = 0.0_dp
        energy = interval(0.0_dp); ok = .false.
        if (.not. valid_state(state)) return
        call coupled_amplitude_interval(interval(state%source_drive), &
                                        interval(state%nu), interval(state%fp_rate), &
                                        interval(state%temperature), exact_a)
        if (.not. finite_interval(exact_a)) return
        if (exact_a%lo <= -1.0_dp .or. exact_a%hi >= 1.0_dp) return
        call difference_interval(interval(state%amplitude), exact_a, delta)
        aradius = interval(mag(delta))
        call poisson_coefficients_interval(interval(real(state%k, dp)), &
                                           interval(1.0_dp), exact_a, &
                                           interval(0.0_dp), pc, ps, ec, es)
        call difference_interval(interval(state%field%phi_cos), pc, pa)
        call difference_interval(interval(state%field%phi_sin), ps, pb)
        call difference_interval(interval(state%field%e_cos), ec, ea)
        call difference_interval(interval(state%field%e_sin), es, eb)
        call sum_interval(interval(mag(pa)), interval(mag(pb)), pradius)
        call sum_interval(interval(mag(ea)), interval(mag(eb)), eradius)
        call coupled_radius_interval(aradius, pradius, eradius, fr, pr, er)
        call poisson_energy_interval(interval(real(state%k, dp)), &
                                     interval(1.0_dp), exact_a, interval(0.0_dp), &
                                     energy)
        if (.not. finite_interval(fr)) return
        if (.not. finite_interval(pr)) return
        if (.not. finite_interval(er)) return
        if (.not. finite_interval(energy)) return
        f_radius = fr%hi; phi_radius = pr%hi; e_radius = er%hi
        ok = .true.
    end subroutine

    logical function coupled_certificate_accepts(state, f_radius, phi_radius, &
                                                 e_radius) result(accepts)
        type(coupled_steady_t), intent(in) :: state
        real(dp), intent(in) :: f_radius, phi_radius, e_radius
        real(dp) :: needed_f, needed_phi, needed_e
        type(interval_t) :: energy
        logical :: ok
        accepts = .false.
        if (any(.not. ieee_is_finite([f_radius, phi_radius, e_radius]))) return
        if (min(f_radius, phi_radius, e_radius) < 0.0_dp) return
        call coupled_certificate_radius(state, needed_f, needed_phi, needed_e, &
                                          energy, ok)
        if (.not. ok) return
        accepts = f_radius >= needed_f .and. phi_radius >= needed_phi &
                  .and. e_radius >= needed_e
    end function

    logical function finite_interval(value) result(valid)
        type(interval_t), intent(in) :: value
        valid = all(ieee_is_finite([value%lo, value%hi]))
        if (.not. valid) return
        valid = value%lo <= value%hi
    end function

    logical function valid_field(field) result(valid)
        type(poisson_mode_t), intent(in) :: field
        valid = .false.
        if (field%k <= 0) return
        if (any(.not. ieee_is_finite([field%epsilon0, field%rho_mean, &
                                     field%rho_cos, field%rho_sin, field%phi_mean, &
                                     field%phi_cos, field%phi_sin, &
                                     field%e_cos, field%e_sin]))) return
        if (field%epsilon0 <= 0.0_dp) return
        if (abs(field%rho_mean) > 0.0_dp) return
        if (abs(field%phi_mean) > 0.0_dp) return
        valid = .true.
    end function

    logical function valid_parameters(state) result(valid)
        type(coupled_steady_t), intent(in) :: state
        real(dp) :: denominator, exact_amplitude
        valid = .false.
        if (state%k <= 0) return
        if (any(.not. ieee_is_finite([state%source_drive, state%nu, state%fp_rate, &
                                     state%temperature]))) return
        if (state%nu < 0.0_dp .or. state%fp_rate < 0.0_dp) return
        if (state%temperature <= 0.0_dp) return
        call coupled_denominator(state%nu, state%fp_rate, state%temperature, &
                                  denominator)
        if (.not. ieee_is_finite(denominator)) return
        if (denominator <= 0.0_dp) return
        call coupled_amplitude(state%source_drive, state%nu, state%fp_rate, &
                               state%temperature, exact_amplitude)
        if (.not. ieee_is_finite(exact_amplitude)) return
        if (abs(exact_amplitude) >= 1.0_dp) return
        valid = .true.
    end function

    logical function valid_state(state) result(valid)
        type(coupled_steady_t), intent(in) :: state
        valid = .false.
        if (.not. valid_parameters(state)) return
        if (.not. valid_field(state%field)) return
        if (state%field%k /= state%k) return
        if (state%field%epsilon0 < 1.0_dp .or. state%field%epsilon0 > 1.0_dp) return
        if (.not. ieee_is_finite(state%amplitude)) return
        if (abs(state%amplitude) >= 1.0_dp) return
        if (state%field%rho_cos < state%amplitude) return
        if (state%field%rho_cos > state%amplitude) return
        if (abs(state%field%rho_sin) > 0.0_dp) return
        valid = .true.
    end function
    subroutine poisson_sample(field, x, potential, electric, ok)
        type(poisson_mode_t), intent(in) :: field
        real(dp), intent(in) :: x
        real(dp), intent(out) :: potential, electric
        logical, intent(out) :: ok
        potential = 0.0_dp; electric = 0.0_dp
        ok = .false.
        if (.not. valid_field(field)) return
        if (.not. ieee_is_finite(x)) return
        call poisson_reconstruct(real(field%k, dp), field%phi_cos, field%phi_sin, &
                                 field%e_cos, field%e_sin, x, potential, electric)
        ok = all(ieee_is_finite([potential, electric]))
    end subroutine

    subroutine solve_coupled_steady(source_drive, nu, k, state, ok, &
                                    fp_rate, temperature)
        ! Direct solve in the finite polynomial/Fourier manufactured family.
        ! Its density moment is solved first; Poisson supplies the actual field.
        ! Sources and zero velocity Dirichlet data define an open system.
        real(dp), intent(in) :: source_drive, nu
        integer, intent(in) :: k
        type(coupled_steady_t), intent(out) :: state
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: fp_rate, temperature
        real(dp) :: density_mean, density_cos
        ok = .false.
        state%k = k; state%source_drive = source_drive; state%nu = nu
        if (present(fp_rate)) state%fp_rate = fp_rate
        if (present(temperature)) state%temperature = temperature
        if (.not. valid_parameters(state)) return
        call coupled_amplitude(source_drive, nu, state%fp_rate, &
                               state%temperature, state%amplitude)
        if (.not. ieee_is_finite(state%amplitude)) return
        if (abs(state%amplitude) >= 1.0_dp) return
        call coupled_density(state%amplitude, density_mean, density_cos)
        ! Generated exact mass is one; fixed ions carry charge -1 exactly.
        if (density_mean < 1.0_dp .or. density_mean > 1.0_dp) return
        call solve_poisson_mode(k, 1.0_dp, 0.0_dp, density_cos, 0.0_dp, &
                                state%field, ok)
    end subroutine

    subroutine coupled_sample(state, x, v, f, electric, potential, source, &
                              residual, ok)
        type(coupled_steady_t), intent(in) :: state
        real(dp), intent(in) :: x, v
        real(dp), intent(out) :: f, electric, potential, source, residual
        logical, intent(out) :: ok
        f = 0.0_dp; electric = 0.0_dp; potential = 0.0_dp
        source = 0.0_dp; residual = 0.0_dp; ok = .false.
        if (.not. valid_state(state)) return
        if (any(.not. ieee_is_finite([x, v]))) return
        if (abs(v) > 1.0_dp) return
        call coupled_distribution(state%amplitude, real(state%k, dp), x, v, f)
        call poisson_sample(state%field, x, potential, electric, ok)
        if (.not. ok) return
        call coupled_source(state%source_drive, state%nu, real(state%k, dp), x, v, &
                            state%fp_rate, state%temperature, source)
        call coupled_residual(state%amplitude, real(state%k, dp), state%nu, x, v, &
                              state%field%e_cos, state%field%e_sin, source, &
                              state%fp_rate, state%temperature, residual)
        ok = all(ieee_is_finite([f, electric, potential, source, residual]))
    end subroutine

    subroutine coupled_reference(state, x, v, f, electric, potential, ok)
        type(coupled_steady_t), intent(in) :: state
        real(dp), intent(in) :: x, v
        real(dp), intent(out) :: f, electric, potential
        logical, intent(out) :: ok
        f = 0.0_dp; electric = 0.0_dp; potential = 0.0_dp; ok = .false.
        if (.not. valid_parameters(state)) return
        if (any(.not. ieee_is_finite([x, v]))) return
        if (abs(v) > 1.0_dp) return
        call coupled_reference_value(state%source_drive, state%nu, state%fp_rate, &
                                     state%temperature, real(state%k, dp), x, v, &
                                     f, electric, potential)
        ok = all(ieee_is_finite([f, electric, potential]))
    end subroutine
end module
