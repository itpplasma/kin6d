module kin6d_streaming_particle_density
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, rad, hull, &
                                linear_interpolation_remainder
    use kin6d_streaming_comparison, only: streaming_state_t, streaming_particles_t, &
                      streaming_shape, get_streaming_initial, get_streaming_particles, &
                             streaming_particle_info, streaming_particles_source_matches
    use kin6d_streaming_leaves, only: streaming_characteristic, &
                           streaming_characteristic_interval, streaming_particle_mean, &
                            streaming_particle_mean_interval, streaming_particle_mode, &
                  streaming_particle_mode_interval, streaming_cell_curvature_interval, &
                       streaming_particle_cell_error_interval, streaming_add_interval, &
                                      streaming_particle_coefficient_band_interval, &
                      streaming_mode_accumulate_interval, streaming_x_curvature_interval
    implicit none
    private
    public :: particle_density_workspace_t, prepare_particle_density_workspace
    public :: particle_density_certificate_t, project_particle_density
    public :: certify_particle_density, sample_particle_density_band
    public :: particle_density_chord, particle_density_budgets, &
        particle_density_coefficients

    ! The observable is rho(x,t)=1/2 integral[-1,1] f(x,v,t)dv, reconstructed
    ! in the admitted finite Fourier space. This is not a density of Dirac masses.
    ! Sensitivities of this weak certificate are not currently supported.
    type :: particle_density_workspace_t
        private
        type(streaming_state_t) :: source
        real(dp), allocatable :: x(:), v(:), w(:), grid(:), c0(:, :), s0(:, :)
        type(interval_t), allocatable :: xi(:), vi(:), wi(:), ac(:), bc(:)
        integer :: nx = 0, nv = 0, cells = 0, modes = -1
        logical :: global_velocity = .false., ready = .false.
    end type

    type :: particle_density_certificate_t
        private
        type(streaming_state_t) :: source
        type(interval_t), allocatable :: cosine(:), sine(:)
        type(interval_t) :: rounding, quadrature
        real(dp) :: t = 0.0_dp
        integer :: nx = 0, nv = 0, cells = 0, modes = -1
        logical :: global_velocity = .false., valid = .false.
    end type
contains
    subroutine prepare_particle_density_workspace(initial, particles, workspace, ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        type(particle_density_workspace_t), intent(inout) :: workspace
        logical, intent(out) :: ok
        integer :: nodes, count

        ok = .false.
        workspace%ready = .false.
        if (.not. streaming_particles_source_matches(initial, particles)) return
        call streaming_shape(initial, workspace%modes, nodes, ok)
        if (.not. ok) return
        call streaming_particle_info(particles, workspace%nx, workspace%nv, &
                                     workspace%cells, workspace%global_velocity, ok)
        if (.not. ok) return
        ok = .false.
        if (workspace%nx <= 2*workspace%modes) return
        count = workspace%nx*workspace%nv*workspace%cells
        if (allocated(workspace%x)) then
            if (size(workspace%x) /= count .or. size(workspace%grid) /= nodes .or. &
                size(workspace%ac) /= workspace%modes + 1) then
                deallocate (workspace%x, workspace%v, workspace%w, workspace%xi, &
                            workspace%vi, workspace%wi, workspace%grid, workspace%c0, &
                            workspace%s0, workspace%ac, workspace%bc)
            end if
        end if
        if (.not. allocated(workspace%x)) then
            allocate (workspace%x(count), workspace%v(count), workspace%w(count), &
                      workspace%xi(count), workspace%vi(count), workspace%wi(count), &
                      workspace%grid(nodes), workspace%c0(0:workspace%modes, nodes), &
                      workspace%s0(0:workspace%modes, nodes), &
                      workspace%ac(0:workspace%modes), workspace%bc(0:workspace%modes))
        end if
        call get_streaming_particles(particles, workspace%x, workspace%v, &
                              workspace%w, ok, workspace%xi, workspace%vi, workspace%wi)
        if (.not. ok) return
        call get_streaming_initial(initial, workspace%grid, workspace%c0, &
                                   workspace%s0, ok)
        if (.not. ok) return
        workspace%source = initial
        workspace%ready = .true.
    end subroutine

    subroutine project_particle_density(particles, t, cosine, sine, ok, workspace)
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        real(dp), intent(out) :: cosine(0:), sine(0:)
        logical, intent(out) :: ok
        type(particle_density_workspace_t), intent(in), optional :: workspace
        real(dp), allocatable :: x(:), v(:), w(:)
        integer :: nx, nv, cells, modes
        logical :: global_velocity

        ok = .false.
        if (.not. ieee_is_finite(t)) return
        if (size(cosine) /= size(sine) .or. size(cosine) < 1) return
        if (present(workspace)) then
            if (.not. workspace_matches(workspace, particles)) return
            if (size(cosine) /= workspace%modes + 1) return
            call float_projection(workspace%x, workspace%v, workspace%w, t, &
                                  cosine, sine, ok)
        else
      call streaming_particle_info(particles, nx, nv, cells, global_velocity, ok, modes)
            if (.not. ok) return
            ok = .false.
            if (size(cosine) /= modes + 1) return
            if (nx <= 2*modes) return
            allocate (x(nx*nv*cells), v(nx*nv*cells), w(nx*nv*cells))
            call get_streaming_particles(particles, x, v, w, ok)
            if (.not. ok) return
            call float_projection(x, v, w, t, cosine, sine, ok)
        end if
    end subroutine

    subroutine float_projection(x0, v, w, t, cosine, sine, ok)
        real(dp), intent(in) :: x0(:), v(:), w(:), t
        real(dp), intent(out) :: cosine(0:), sine(0:)
        logical, intent(out) :: ok
        real(dp) :: x, a, b
        integer :: p, k

        cosine = 0.0_dp
        sine = 0.0_dp
        do p = 1, size(x0)
            call streaming_characteristic(x0(p), v(p), t, x)
            call streaming_particle_mean(w(p), cosine(0), a)
            cosine(0) = a
            do k = 1, size(cosine) - 1
                call streaming_particle_mode(x, real(k, dp), w(p), cosine(k), &
                                             sine(k), a, b)
                cosine(k) = a
                sine(k) = b
            end do
        end do
        ok = all(ieee_is_finite(cosine)) .and. all(ieee_is_finite(sine))
    end subroutine

  subroutine certify_particle_density(initial, particles, t, certificate, ok, workspace)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(inout) :: certificate
        logical, intent(out) :: ok
        type(particle_density_workspace_t), intent(inout), optional :: workspace
        type(particle_density_workspace_t) :: local

        ok = .false.
        certificate%valid = .false.
        if (.not. ieee_is_finite(t)) return
        if (.not. streaming_particles_source_matches(initial, particles)) return
        if (present(workspace)) then
            if (.not. workspace_matches(workspace, particles)) return
            call interval_projection(workspace, t, certificate, ok)
        else
            call prepare_particle_density_workspace(initial, particles, local, ok)
            if (.not. ok) return
            call interval_projection(local, t, certificate, ok)
        end if
    end subroutine

    subroutine interval_projection(work, t, certificate, ok)
        type(particle_density_workspace_t), intent(inout) :: work
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(inout) :: certificate
        logical, intent(out) :: ok
        type(interval_t) :: x, a, b, next, bound, error, lower, upper
        type(interval_t) :: rounding, quadrature, coeff_error
        integer :: p, k, j, last

        ok = .false.
        work%ac = interval(0.0_dp)
        work%bc = interval(0.0_dp)
        do p = 1, size(work%x)
            ! Hulls cover both ideal mathematical geometry and stored float markers.
         call streaming_characteristic_interval(hull(work%xi(p), interval(work%x(p))), &
                                  hull(work%vi(p), interval(work%v(p))), interval(t), x)
          call streaming_particle_mean_interval(hull(work%wi(p), interval(work%w(p))), &
                                                  work%ac(0), a)
            work%ac(0) = a
            do k = 1, work%modes
                call streaming_particle_mode_interval(x, interval(real(k, dp)), &
                    hull(work%wi(p), interval(work%w(p))), work%ac(k), work%bc(k), a, b)
                work%ac(k) = a
                work%bc(k) = b
            end do
        end do
        if (.not. all(finite_interval(work%ac))) return
        if (.not. all(finite_interval(work%bc))) return
        if (allocated(certificate%cosine)) then
            if (size(certificate%cosine) /= work%modes + 1) &
                deallocate (certificate%cosine, certificate%sine)
        end if
        if (.not. allocated(certificate%cosine)) &
            allocate (certificate%cosine(0:work%modes), certificate%sine(0:work%modes))
        certificate%cosine = work%ac
        certificate%sine = work%bc
        call streaming_add_interval(interval(rad(work%ac(0))), &
                                    interval(rad(work%ac(0))), rounding)
        quadrature = interval(0.0_dp)
        do k = 1, work%modes
            ! Twice the interval radius bounds displacement about any contained
            ! float candidate, including the rounded interval midpoint offset.
            call add_candidate_rounding(rounding, work%ac(k), next)
            rounding = next
            call add_candidate_rounding(rounding, work%bc(k), next)
            rounding = next
            coeff_error = interval(0.0_dp)
            last = size(work%grid) - 1
            if (work%global_velocity) last = 1
            do j = 1, last
                if (work%global_velocity) then
            call streaming_cell_curvature_interval(interval(t), interval(real(k, dp)), &
                         interval(-1.0_dp), interval(1.0_dp), interval(work%c0(k, 1)), &
                                     interval(work%c0(k, 1)), interval(work%s0(k, 1)), &
                                                         interval(work%s0(k, 1)), bound)
                    call streaming_particle_cell_error_interval(interval(-1.0_dp), &
                            interval(1.0_dp), interval(real(work%nv, dp)), bound, error)
                else
            call streaming_cell_curvature_interval(interval(t), interval(real(k, dp)), &
                                   interval(work%grid(j)), interval(work%grid(j + 1)), &
                                 interval(work%c0(k, j)), interval(work%c0(k, j + 1)), &
                            interval(work%s0(k, j)), interval(work%s0(k, j + 1)), bound)
                   call streaming_particle_cell_error_interval(interval(work%grid(j)), &
                  interval(work%grid(j + 1)), interval(real(work%nv, dp)), bound, error)
                end if
                if (.not. finite_interval(error)) return
                error = interval(0.0_dp, max(0.0_dp, error%hi))
                call streaming_add_interval(coeff_error, error, next)
                coeff_error = next
            end do
            ! Individual C and S remainders each use M; their L1 sum uses two.
call streaming_particle_coefficient_band_interval(work%ac(k), coeff_error, lower, upper)
            certificate%cosine(k) = interval(lower%lo, upper%hi)
call streaming_particle_coefficient_band_interval(work%bc(k), coeff_error, lower, upper)
            certificate%sine(k) = interval(lower%lo, upper%hi)
            call streaming_add_interval(quadrature, coeff_error, next)
            quadrature = next
            call streaming_add_interval(quadrature, coeff_error, next)
            quadrature = next
        end do
        if (.not. all(finite_interval(certificate%cosine))) return
        if (.not. all(finite_interval(certificate%sine))) return
      if (.not. finite_interval(rounding) .or. .not. finite_interval(quadrature)) return
        certificate%rounding = rounding
        certificate%quadrature = quadrature
        certificate%source = work%source
        certificate%t = t
        certificate%nx = work%nx
        certificate%nv = work%nv
        certificate%cells = work%cells
        certificate%modes = work%modes
        certificate%global_velocity = work%global_velocity
        certificate%valid = .true.
        ok = .true.
    end subroutine

    subroutine sample_particle_density_band(initial, particles, t, certificate, x, &
        lower, upper, ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(in) :: certificate
        type(interval_t), intent(in) :: x
        type(interval_t), intent(out) :: lower, upper
        logical, intent(out) :: ok
        type(interval_t) :: acc, next
        integer :: k

        ok = .false.
        if (.not. certificate_matches(initial, particles, t, certificate)) return
        if (.not. finite_interval(x)) return
        acc = certificate%cosine(0)
        do k = 1, certificate%modes
            call streaming_mode_accumulate_interval(x, interval(real(k, dp)), &
                                  certificate%cosine(k), certificate%sine(k), acc, next)
            acc = next
        end do
        if (.not. finite_interval(acc)) return
        lower = interval(acc%lo)
        upper = interval(acc%hi)
        ok = .true.
    end subroutine

    subroutine particle_density_chord(initial,particles,t,certificate,xl,xr,padding,ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(in) :: certificate
        type(interval_t), intent(in) :: xl, xr
        type(interval_t), intent(out) :: padding
        logical, intent(out) :: ok
        type(interval_t) :: curvature, next
        integer :: k

        ok = .false.
        if (.not. certificate_matches(initial, particles, t, certificate)) return
        if (.not. finite_interval(xl) .or. .not. finite_interval(xr)) return
        curvature = interval(0.0_dp)
        do k = 1, certificate%modes
            call streaming_x_curvature_interval(interval(real(k, dp)), &
                            certificate%cosine(k), certificate%sine(k), curvature, next)
            curvature = next
        end do
        if (.not. finite_interval(curvature)) return
        curvature%lo = max(0.0_dp, curvature%lo)
        padding = linear_interpolation_remainder(xl, xr, curvature)
        ok = finite_interval(padding)
    end subroutine

    subroutine particle_density_budgets(certificate, rounding, quadrature, ok)
        type(particle_density_certificate_t), intent(in) :: certificate
        type(interval_t), intent(out) :: rounding, quadrature
        logical, intent(out) :: ok

        ok = certificate%valid
        if (.not. ok) return
        rounding = certificate%rounding
        quadrature = certificate%quadrature
    end subroutine

    subroutine add_candidate_rounding(acc, coefficient, next)
        type(interval_t), intent(in) :: acc, coefficient
        type(interval_t), intent(out) :: next
        type(interval_t) :: diameter

        call streaming_add_interval(interval(rad(coefficient)), &
                                    interval(rad(coefficient)), diameter)
        call streaming_add_interval(acc, diameter, next)
    end subroutine

subroutine particle_density_coefficients(initial,particles,t,certificate,cosine,sine,ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(in) :: certificate
        type(interval_t), intent(out) :: cosine(0:), sine(0:)
        logical, intent(out) :: ok

        ok = .false.
        if (.not. certificate_matches(initial, particles, t, certificate)) return
       if (size(cosine) /= certificate%modes + 1 .or. size(sine) /= size(cosine)) return
        cosine = certificate%cosine
        sine = certificate%sine
        ok = .true.
    end subroutine

    logical function workspace_matches(work, particles) result(ok)
        type(particle_density_workspace_t), intent(in) :: work
        type(streaming_particles_t), intent(in) :: particles
        integer :: nx, nv, cells
        logical :: global_velocity, prepared

        ok = .false.
        if (.not. work%ready) return
        if (.not. streaming_particles_source_matches(work%source, particles)) return
       call streaming_particle_info(particles, nx, nv, cells, global_velocity, prepared)
        if (.not. prepared) return
        ok = nx == work%nx .and. nv == work%nv .and. cells == work%cells .and. &
             (global_velocity .eqv. work%global_velocity)
    end function

    logical function certificate_matches(initial, particles, t, certificate) result(ok)
        type(streaming_state_t), intent(in) :: initial
        type(streaming_particles_t), intent(in) :: particles
        real(dp), intent(in) :: t
        type(particle_density_certificate_t), intent(in) :: certificate
        integer :: nx, nv, cells
        logical :: global_velocity, prepared

        ok = .false.
        if (.not. certificate%valid .or. .not. ieee_is_finite(t)) return
        if (t < certificate%t .or. t > certificate%t) return
        if (.not. streaming_particles_source_matches(initial, particles)) return
     if (.not. streaming_particles_source_matches(certificate%source, particles)) return
       call streaming_particle_info(particles, nx, nv, cells, global_velocity, prepared)
        if (.not. prepared) return
        ok = nx == certificate%nx .and. nv == certificate%nv .and. &
             cells == certificate%cells .and. &
             (global_velocity .eqv. certificate%global_velocity)
    end function

    pure elemental logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value

 ok = ieee_is_finite(value%lo) .and. ieee_is_finite(value%hi) .and. value%lo <= value%hi
    end function
end module kin6d_streaming_particle_density
