module kin6d_geqdsk_adapter
    ! Outer interchange boundary: SI / COCOS 3 inside, libneo CGS outside.
    ! Files and classification are discrete operations, without derivatives or
    ! continuous error enclosures. Imported and extended fields are candidates.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use geqdsk_tools, only: geqdsk_t, geqdsk_read, geqdsk_write, &
                            geqdsk_classify, geqdsk_standardise
    use math_constants, only: length_si_to_cgs, magneticflux_si_to_cgs, &
                          magneticfield_si_to_cgs, current_si_to_cgs, pressure_si_to_cgs
    use kin6d_real_identity, only: same_real_value
    use kin6d_equilibrium_state, only: equilibrium_state_t, equilibrium_field_t, &
                                       EQ_PSI, EQ_VALUES
    use kin6d_equilibrium_readback, only: equilibrium_readback_t, read_equilibrium
    use kin6d_gs_equilibrium, only: gs_field_t, gs_equilibrium
    use kin6d_gs_profile_law, only: gs_profile_law_t, gs_law_sources
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use kin6d_gs_nonlinear, only: solve_gs_cubic_profiles
    use kin6d_vacuum_surfaces, only: psi_grid_t, psi_grid_init
    implicit none
    private
    public :: geqdsk_data_t, geqdsk_import_t, export_geqdsk, import_geqdsk
    public :: geqdsk_grid, resolve_geqdsk, native_cocos
    integer, parameter :: native_cocos = 3

    type :: geqdsk_data_t
        ! Uniform R,Z grid and profiles on uniform axis-to-edge flux, in SI.
        real(dp), allocatable :: r(:), z(:), psi(:, :), flux(:)
        real(dp), allocatable :: f(:), p(:), pp(:), ff(:), q(:)
        real(dp), allocatable :: lcfs(:, :), limiter(:, :)
        real(dp) :: axis(2) = 0, psi_axis = 0, psi_edge = 0
        real(dp) :: rcentr = 0, bcentr = 0, ip = 0
        integer :: cocos = native_cocos
        character(48) :: header = ''
        character(32) :: claim = 'candidate', bound_status = 'unavailable'
        character(80) :: exterior = 'source grid; vacuum physics unspecified'
    end type

    type :: geqdsk_import_t
        type(geqdsk_data_t) :: original, data
        integer :: source_cocos = 0, classified_cocos = 0
        logical :: profiles_repaired = .false.
        character(1024) :: source = ''
    end type
contains
    subroutine import_geqdsk(path, source_cocos, imported, status)
        character(*), intent(in) :: path
        integer, intent(in) :: source_cocos
        type(geqdsk_import_t), intent(out) :: imported
        integer, intent(out) :: status
        type(geqdsk_t) :: g
        logical :: exists
        status = 1
        ! libneo's gfile classifier assumes right-handed cylindrical axes.
        if (.not. any(source_cocos == [1, 3, 5, 7, 11, 13, 15, 17])) return
        inquire (file=path, exist=exists)
        if (.not. exists) return
        imported%source = path
        imported%source_cocos = source_cocos
        call geqdsk_read(g, path)
        call from_libneo(g, imported%original)
        imported%original%cocos = source_cocos
        call geqdsk_classify(g)
        imported%classified_cocos = g%cocos%index
        imported%profiles_repaired = &
            any(.not. same_real_value(g%pprime/ &
                   (pressure_si_to_cgs/magneticflux_si_to_cgs), imported%original%pp)) &
            .or. any(.not. same_real_value(g%ffprim/magneticfield_si_to_cgs, &
                                           imported%original%ff)) &
            .or. any(.not. same_real_value(g%fpol/ &
                       (magneticfield_si_to_cgs*length_si_to_cgs), imported%original%f))
        status = 2
        if (g%cocos%index /= source_cocos) return
        call geqdsk_standardise(g)
        call from_libneo(g, imported%data)
        imported%data%cocos = native_cocos
        status = 0
    end subroutine

    subroutine export_geqdsk(path, eq, rectangle, nr, nz, lcfs, limiter, ip, &
                             rcentr, cocos, data, status, continuation)
        character(*), intent(in) :: path
        type(equilibrium_state_t), intent(in) :: eq
        real(dp), intent(in) :: rectangle(4), lcfs(:, :), limiter(:, :), ip, rcentr
        integer, intent(in) :: nr, nz, cocos
        type(geqdsk_data_t), intent(out) :: data
        integer, intent(out) :: status
        ! Explicit candidate continuation for points outside the native domain.
        class(equilibrium_field_t), intent(in), optional :: continuation
        type(equilibrium_readback_t) :: readback
        type(geqdsk_t) :: g
        real(dp), allocatable :: points(:, :), values(:, :), rho(:)
        integer, allocatable :: statuses(:)
        real(dp) :: pp, ff, dpp, dff, app, aff
        integer :: i, j, k
        status = 1
        if (cocos /= native_cocos) return
        if (nr < 11 .or. nz < 11) return
        if (.not. allocated(eq%field)) return
        if (.not. all(ieee_is_finite(rectangle))) return
        if (rectangle(1) <= 0 .or. rectangle(2) <= rectangle(1)) return
        if (rectangle(4) <= rectangle(3) .or. rcentr <= 0) return
        if (size(lcfs, 1) /= 2 .or. size(lcfs, 2) < 3) return
        if (size(limiter, 1) /= 2) return
        allocate (rho(nr))
        rho = [(sqrt(real(i - 1, dp)/(nr - 1)), i=1, nr)]
        call read_equilibrium(eq, nr, 96, readback, status, rho_grid=rho)
        if (status /= 0) return
        status = 1
        allocate (data%r(nr), data%z(nz), data%psi(nr, nz), data%f(nr), &
                  data%p(nr), data%pp(nr), data%ff(nr), points(2, nr*nz), &
                  values(EQ_VALUES, nr*nz), statuses(nr*nz))
        data%r = rectangle(1) + [(real(i - 1, dp)/(nr - 1), i=1, nr)] &
                 *(rectangle(2) - rectangle(1))
        data%z = rectangle(3) + [(real(i - 1, dp)/(nz - 1), i=1, nz)] &
                 *(rectangle(4) - rectangle(3))
        do j = 1, nz
            do i = 1, nr
                k = i + (j - 1)*nr
                points(:, k) = [data%r(i), data%z(j)]
            end do
        end do
        call eq%field%sample(points, values, statuses)
        data%exterior = 'represented by source D7 evaluator'
        if (any(statuses /= 0)) then
            if (.not. present(continuation)) return
            do k = 1, nr*nz
                if (statuses(k) == 0) cycle
                call continuation%sample(points(:, k:k), values(:, k:k), &
                                         statuses(k:k))
            end do
            if (any(statuses /= 0)) return
            data%exterior = 'explicit caller continuation; candidate, no vacuum claim'
        end if
        data%psi = reshape(values(EQ_PSI, :), [nr, nz])
        data%flux = readback%psi
        data%q = readback%q
        data%axis = readback%axis
        data%psi_axis = readback%psi_axis; data%psi_edge = eq%psi_bound
        data%ip = ip; data%rcentr = rcentr; data%bcentr = eq%f_vac/rcentr
        data%lcfs = lcfs; data%limiter = limiter
        data%header = 'KIN6D COCOS 3 SI; candidate equilibrium'
        select type (field => eq%field)
        class is (gs_field_t)
            do i = 1, nr
                call field%profiles(data%flux(i) - eq%psi_bound, &
                                    data%f(i), data%p(i))
                data%pp(i) = field%pprime; data%ff(i) = field%ffprime
                if (allocated(field%law)) then
                    call gs_law_sources(field%law, data%flux(i) - eq%psi_bound, &
                                        pp, ff, dpp, dff, app, aff)
                    data%pp(i) = pp; data%ff(i) = ff
                end if
            end do
        class default
            return
        end select
        if (.not. all(ieee_is_finite(values))) return
        call to_libneo(data, g)
        call geqdsk_classify(g)
        status = 2
        if (g%cocos%index /= cocos) return
        ! Qualification may repair inconsistent profiles. Never export a
        ! repaired state while returning the caller's unrepaired SI samples.
        status = 3
        if (any(.not. same_real_value(g%fpol, data%f* &
                                    (magneticfield_si_to_cgs*length_si_to_cgs)))) return
        if (any(.not. same_real_value(g%pprime, data%pp* &
                                   (pressure_si_to_cgs/magneticflux_si_to_cgs)))) return
        if (any(.not. same_real_value(g%ffprim, &
                                      data%ff*magneticfield_si_to_cgs))) return
        call geqdsk_write(g, path)
        status = 0
    end subroutine

    subroutine geqdsk_grid(data, grid, status)
        type(geqdsk_data_t), intent(in) :: data
        type(psi_grid_t), intent(out) :: grid
        integer, intent(out) :: status
        integer :: nr, nz
        status = 1
        if (data%cocos /= native_cocos) return
        nr = size(data%r); nz = size(data%z)
        call psi_grid_init(grid, data%psi, data%r(1), data%r(nr) - data%r(1), &
                           (data%z(1) + data%z(nz))/2, data%z(nz) - data%z(1), &
                           data%axis(1), data%axis(2), data%psi_edge, data%ip, data%f, &
                           status, pressure=data%p)
    end subroutine

    subroutine resolve_geqdsk(data, max_area, state, eq, status, edge_pressure, &
                              curved_lcfs)
        type(geqdsk_data_t), intent(in) :: data
        real(dp), intent(in) :: max_area
        type(gs_state_t), intent(out) :: state
        type(equilibrium_state_t), intent(out) :: eq
        integer, intent(out) :: status
        real(dp), intent(in), optional :: edge_pressure
        logical, intent(in), optional :: curved_lcfs
        type(gs_profile_law_t) :: law
        real(dp), allocatable :: flux(:), pp(:), ff(:)
        integer :: n, i, nb
        status = 1
        if (data%cocos /= native_cocos) return
        nb = size(data%lcfs, 2)
        if (nb < 3) return
        ! Gfiles may repeat the closing point; the native mesher closes it.
        if (norm2(data%lcfs(:, nb) - data%lcfs(:, 1)) < 1.0e-10_dp) nb = nb - 1
        ! Explicit smooth-contour reconstruction. Without this choice the
        ! native polygon is retained, including any corners. Interior refinement
        ! alone cannot remove that polygon's boundary reconstruction error.
        call initialize_gs(state, data%lcfs(:, :nb), max_area, status, &
                           element_order=2, curved=curved_lcfs)
        if (status /= 0) return
        flux = data%flux - data%psi_edge
        pp = data%pp; ff = data%ff
        n = size(flux)
        if (flux(n) < flux(1)) then
            flux = flux(n:1:-1); pp = pp(n:1:-1); ff = ff(n:1:-1)
        end if
        allocate (law%knots(n), law%pp(4, n - 1), law%ff(4, n - 1))
        law%knots = flux; law%pp = 0; law%ff = 0
        do i = 1, n - 1
            law%pp(1, i) = pp(i); law%ff(1, i) = ff(i)
            law%pp(2, i) = (pp(i + 1) - pp(i))/(flux(i + 1) - flux(i))
            law%ff(2, i) = (ff(i + 1) - ff(i))/(flux(i + 1) - flux(i))
        end do
        law%p_edge = data%p(n); law%f_edge = data%f(n)
        if (present(edge_pressure)) law%p_edge = edge_pressure
        law%axis = data%psi_axis - data%psi_edge
        ! Cold native Newton solve; no interior imported psi is copied.
        call solve_gs_cubic_profiles(state, eq%mu0, law, status, &
                          relative_tolerance=1.0e-9_dp, quadrature_tolerance=1.0e-11_dp)
        if (status /= 0) return
        call gs_equilibrium(state, eq, status, data%f(n), 0.0_dp, law=law)
        if (status /= 0) return
        write (eq%seed, '(a,l1,a,es24.16)') &
            'G-EQDSK LCFS/profiles; cold solve; curved LCFS ', state%mapped, &
            '; edge pressure offset [Pa] ', &
            law%p_edge - data%p(n)
    end subroutine

    subroutine from_libneo(g, data)
        type(geqdsk_t), intent(in) :: g
        type(geqdsk_data_t), intent(out) :: data
        data%r = g%R_eqd/length_si_to_cgs; data%z = g%Z_eqd/length_si_to_cgs
        data%psi = g%psirz/magneticflux_si_to_cgs
        data%flux = g%psi_eqd/magneticflux_si_to_cgs
        data%f = g%fpol/(magneticfield_si_to_cgs*length_si_to_cgs)
        data%p = g%pres/pressure_si_to_cgs
        data%pp = g%pprime/(pressure_si_to_cgs/magneticflux_si_to_cgs)
        data%ff = g%ffprim/magneticfield_si_to_cgs; data%q = g%qpsi
        data%axis = [g%rmaxis, g%zmaxis]/length_si_to_cgs
        data%psi_axis = g%simag/magneticflux_si_to_cgs
        data%psi_edge = g%sibry/magneticflux_si_to_cgs
        data%rcentr = g%rcentr/length_si_to_cgs
        data%bcentr = g%bcentr/magneticfield_si_to_cgs
        data%ip = g%current/current_si_to_cgs
        allocate (data%lcfs(2, g%nbbbs), data%limiter(2, g%limitr))
        data%lcfs(1, :) = g%rbbbs/length_si_to_cgs
        data%lcfs(2, :) = g%zbbbs/length_si_to_cgs
        data%limiter(1, :) = g%rlim/length_si_to_cgs
        data%limiter(2, :) = g%zlim/length_si_to_cgs
        data%header = g%header
    end subroutine

    subroutine to_libneo(data, g)
        type(geqdsk_data_t), intent(in) :: data
        type(geqdsk_t), intent(out) :: g
        integer :: nr, nz
        nr = size(data%r); nz = size(data%z)
        g%nw = nr; g%nh = nz; g%nbbbs = size(data%lcfs, 2)
        g%limitr = size(data%limiter, 2); g%header = data%header
        g%rdim = (data%r(nr) - data%r(1))*length_si_to_cgs
        g%rleft = data%r(1)*length_si_to_cgs
        g%zdim = (data%z(nz) - data%z(1))*length_si_to_cgs
        g%zmid = (data%z(nz) + data%z(1))/2*length_si_to_cgs
        g%rcentr = data%rcentr*length_si_to_cgs
        g%rmaxis = data%axis(1)*length_si_to_cgs
        g%zmaxis = data%axis(2)*length_si_to_cgs
        g%simag = data%psi_axis*magneticflux_si_to_cgs
        g%sibry = data%psi_edge*magneticflux_si_to_cgs
        g%bcentr = data%bcentr*magneticfield_si_to_cgs
        g%current = data%ip*current_si_to_cgs
        g%fpol = data%f*(magneticfield_si_to_cgs*length_si_to_cgs)
        g%pres = data%p*pressure_si_to_cgs
        g%pprime = data%pp*(pressure_si_to_cgs/magneticflux_si_to_cgs)
        g%ffprim = data%ff*magneticfield_si_to_cgs
        g%qpsi = data%q; g%psirz = data%psi*magneticflux_si_to_cgs
        g%psi_eqd = data%flux*magneticflux_si_to_cgs
        g%R_eqd = data%r*length_si_to_cgs; g%Z_eqd = data%z*length_si_to_cgs
        g%rbbbs = data%lcfs(1, :)*length_si_to_cgs
        g%zbbbs = data%lcfs(2, :)*length_si_to_cgs
        g%rlim = data%limiter(1, :)*length_si_to_cgs
        g%zlim = data%limiter(2, :)*length_si_to_cgs
        g%fprime = g%ffprim/g%fpol
    end subroutine
end module
