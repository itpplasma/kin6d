module kin6d_microscopic_cloud
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
                              FORTNUM_DOMAIN_ERROR, status_set
    use fortnum_ode_geometric, only: geometric_strang_step
    use kin6d_interaction_leaves, only: interaction_pair_force
    use kin6d_cloud_leaves, only: cloud_background_inside, &
                       cloud_background_outside, cloud_axpy, cloud_norm, cloud_energy, &
                        cloud_initial_position, cloud_initial_velocity, cloud_difference
    implicit none
    private
   public :: cloud_model_t, cloud_initialize, cloud_force, cloud_step, cloud_hamiltonian
    type :: cloud_model_t
        real(dp) :: sphere_charge = 32.0_dp
        real(dp) :: sphere_radius = 3.0_dp
        real(dp) :: impurity_coupling = 0.0_dp
    end type
contains
    subroutine cloud_initialize(state)
        real(dp), intent(out) :: state(192)
        integer :: p, j
        do p = 1, 32
            j = 3*p - 2
          call cloud_initial_position(real(p, dp), state(j), state(j + 1), state(j + 2))
   call cloud_initial_velocity(real(p, dp), state(96 + j), state(97 + j), state(98 + j))
        end do
    end subroutine

    subroutine cloud_force(state, model, force, potential, status)
        real(dp), intent(in) :: state(:)
        type(cloud_model_t), intent(in) :: model
        real(dp), intent(out) :: force(:), potential
        type(fortnum_status_t), intent(out) :: status
        integer :: n, p, q, j, k
        real(dp) :: radius, f(3), energy, dr(3), next, next_force(3)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid cloud state/model')
        potential = 0.0_dp; force = 0.0_dp
        if (mod(size(state), 6) /= 0) return
        n = size(state)/6
        if (size(force) /= 3*n .or. n < 1) return
        if (.not. all(ieee_is_finite(state))) return
        if (.not. ieee_is_finite(model%sphere_radius)) return
        if (model%sphere_radius <= 0.0_dp) return
        if (.not. ieee_is_finite(model%sphere_charge)) return
        if (model%sphere_charge < 0.0_dp) return
        if (.not. ieee_is_finite(model%impurity_coupling)) return
        if (model%impurity_coupling < 0.0_dp) return
        do p = 1, n
            j = 3*p - 2
            call cloud_norm(state(j), state(j + 1), state(j + 2), radius)
            if (.not. ieee_is_finite(radius)) return
            if (radius <= model%sphere_radius) then
                call cloud_background_inside(state(j), state(j + 1), &
                               state(j + 2), model%sphere_charge, model%sphere_radius, &
                                             f(1), f(2), f(3), energy)
            else
                call cloud_background_outside(state(j), state(j + 1), &
                               state(j + 2), model%sphere_charge, model%sphere_radius, &
                                              f(1), f(2), f(3), energy)
            end if
            call cloud_axpy(force(j:j + 2), f, 1.0_dp, next_force)
            force(j:j + 2) = next_force
            call cloud_axpy(potential, energy, 1.0_dp, next)
            potential = next
            if (abs(model%impurity_coupling) > 0.0_dp) then
                if (radius <= 1.0e-10_dp) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                                  'impurity encounter unresolved; no softening applied')
                    return
                end if
                call interaction_pair_force(state(j), state(j + 1), state(j + 2), &
                                      model%impurity_coupling, f(1), f(2), f(3), energy)
                call cloud_axpy(force(j:j + 2), f, 1.0_dp, next_force)
                force(j:j + 2) = next_force
                call cloud_axpy(potential, energy, 1.0_dp, next)
                potential = next
            end if
            do q = p + 1, n
                k = 3*q - 2
                call cloud_difference(state(j:j + 2), state(k:k + 2), dr)
                call cloud_norm(dr(1), dr(2), dr(3), radius)
                if (.not. ieee_is_finite(radius)) return
                if (radius <= 1.0e-10_dp) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                                    'pair encounter unresolved; no softening applied')
                    return
                end if
                call interaction_pair_force(dr(1), dr(2), dr(3), 1.0_dp, &
                                            f(1), f(2), f(3), energy)
                call cloud_axpy(force(j:j + 2), f, 1.0_dp, next_force)
                force(j:j + 2) = next_force
                call cloud_axpy(force(k:k + 2), f, -1.0_dp, next_force)
                force(k:k + 2) = next_force
                call cloud_axpy(potential, energy, 1.0_dp, next)
                potential = next
            end do
        end do
        if (.not. all(ieee_is_finite(force))) return
        if (.not. ieee_is_finite(potential)) return
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine cloud_step(state, dt, model, status)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(cloud_model_t), intent(in) :: model
        type(fortnum_status_t), intent(out) :: status
        call geometric_strang_step(state, dt, kick, drift, status, model)
    end subroutine

    subroutine kick(state, dt, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context
        real(dp) :: force(96), next_state(96), energy
        integer :: n3
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'cloud context required')
        if (.not. present(context)) return
        select type (context)
        type is (cloud_model_t)
            n3 = size(state)/2
            if (n3 > 96) return
            call cloud_force(state, context, force(:n3), energy, status)
            if (status%code /= FORTNUM_OK) return
            n3 = size(state)/2
            call cloud_axpy(state(n3 + 1:), force(:n3), dt, next_state(:n3))
            state(n3 + 1:) = next_state(:n3)
        end select
    end subroutine

    subroutine drift(state, dt, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: dt
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context
        integer :: n3
        real(dp) :: next_state(96)
        n3 = size(state)/2
        if (n3 > 96) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'cloud exceeds 32 particles')
            return
        end if
        call cloud_axpy(state(:n3), state(n3 + 1:), dt, next_state(:n3))
        state(:n3) = next_state(:n3)
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine cloud_hamiltonian(state, model, energy, status)
        real(dp), intent(in) :: state(:)
        type(cloud_model_t), intent(in) :: model
        real(dp), intent(out) :: energy
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: force(96), next
        integer :: j, n3
        n3 = size(state)/2
        if (n3 > 96) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'cloud exceeds 32 particles')
            return
        end if
        call cloud_force(state, model, force(:n3), energy, status)
        if (status%code /= FORTNUM_OK) return
        do j = size(state)/2 + 1, size(state), 3
            call cloud_energy(energy, state(j), state(j + 1), state(j + 2), next)
            energy = next
        end do
    end subroutine
end module
