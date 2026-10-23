module test_mhd_shielding_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_point_basis_t, meridian_point_basis_init
    use kin6d_mhd_mesh, only: flux_mesh_t, flux_mesh_build
    use kin6d_mhd_ampere_n, only: ampere_n_t, ampere_n_build, ampere_n_solve, &
        ampere_n_curl
    use kin6d_mhd_response_n, only: ring_angles_t, ring_angles_set, &
        ring_flux_harmonic, ring_sheet_load

    use kin6d_mhd_mesh, only: polar_map_t
    implicit none
    private
    public :: run_test_mhd_shielding
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    integer, parameter :: nrad = 32, jc = 24, js = 16
    type :: run_test_mhd_shielding_context_t
        real(dp), pointer :: r0 => null()
    end type

    type, extends(polar_map_t) :: run_test_mhd_shielding_circle_map_t
        type(run_test_mhd_shielding_context_t), pointer :: data => null()
    contains
        procedure :: evaluate => run_test_mhd_shielding_circle_map_binding
    end type

contains
    subroutine run_test_mhd_shielding()
    !! L4 resonant sheet currents and ideal shielding.
    !!
    !! 1. Screened response, cylinder limit (asymptotic oracle).  In a
    !!    large-aspect-ratio torus (circular rings about (R0, 0), angles of a
    !!    straight cylinder: theta_s = theta_g, nu = q theta_g) a drive sheet
    !!    with stream function e^{-i m theta} on r_c is shielded by a sheet on
    !!    r_s < r_c that cancels the resonant flux Phi_m(r_s).  In the periodic
    !!    cylinder (k = n/R0) the field off the sheets is mu0 grad phi with
    !!    [phi] = Lambda across a sheet, which gives
    !!        lambda_s/lambda_c = -x_c K_m'(x_c) / (x_s K_m'(x_s)),  x = k r.
    !!    K_m' is evaluated from its integral representation (independent of
    !!    the FEEC/BEM code).  The torus result must approach this limit as
    !!    the inverse aspect ratio a/R0 -> 0 (observed: second order).
    !! 2. Reciprocity: the sheet response matrix M(r, s) = Phi_{m_r}(ring r)
    !!    of the unit sheet s is Hermitian (M = -(1/(pi mu0)) F^H K^{-1} F).
    !! 3. Shielding: with the amplitudes from M lambda = -Phi_drive the total
    !!    resonant flux on every shielded ring vanishes to roundoff.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    real(dp) , target :: r0
    logical :: all_passed
        type(run_test_mhd_shielding_context_t), target :: ctx
        type(run_test_mhd_shielding_circle_map_t), target :: circle_map_operator

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call cylinder_limit(3, 2)
    call cylinder_limit(2, 1)
    call reciprocity_and_shielding()
    if (.not. all_passed) error stop 1
    print '(a)', "test_mhd_shielding: all checks passed"

contains
    subroutine bind_run_test_mhd_shielding_context()
        ctx%r0 => r0
        circle_map_operator%data => ctx
    end subroutine bind_run_test_mhd_shielding_context


    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        if (c) then
            print '(a,a)', "  PASS ", d
        else
            print '(a,a)', "  FAIL ", d
        end if
    end subroutine record


    subroutine cylinder_angles(q, ra)
        real(dp), intent(in) :: q
        type(ring_angles_t), intent(out) :: ra
        integer, parameter :: nf = 256
        real(dp) :: th(nf)
        integer :: i
        do i = 1, nf
            th(i) = 2.0_dp*pi*(i - 1)/nf
        end do
        call ring_angles_set(q, th, q*th, ra)
    end subroutine cylinder_angles

    real(dp) function dbesk(m, x)
        !! K_m'(x) = -int_0^inf cosh t exp(-x cosh t) cosh(m t) dt.
        integer, intent(in) :: m
        real(dp), intent(in) :: x
        real(dp) :: t, h, tmax, s
        integer :: i, nstep
        h = 2.0e-3_dp
        tmax = acosh(80.0_dp/x)
        nstep = int(tmax/h) + 1
        s = 0.5_dp*exp(-x)
        do i = 1, nstep
            t = i*h
            s = s + cosh(t)*exp(-x*cosh(t))*cosh(m*t)
        end do
        dbesk = -s*h
    end function dbesk

    subroutine setup(fm, fs, pb, amp, n)
        type(flux_mesh_t), intent(out) :: fm
        type(meridian_feec_t), intent(out) :: fs
        type(meridian_point_basis_t), intent(out) :: pb
        type(ampere_n_t), intent(inout) :: amp
        integer, intent(in) :: n
        real(dp) :: ul(nrad)
        integer :: j, st
        do j = 1, nrad
            ul(j) = real(j, dp)/nrad
        end do
        call bind_run_test_mhd_shielding_context()
        call flux_mesh_build(ul, 8, 4096, 2, circle_map_operator, fm, st)
        if (st /= 0) error stop "mesh"
        call meridian_feec_build(fm%mm, 2, fs, st)
        if (st /= 0) error stop "feec"
        call ampere_n_build(fs, n, amp, st)
        if (st /= 0) error stop "ampere"
        call meridian_point_basis_init(fs, pb, st)
    end subroutine setup

    function sheet_field(fs, pb, fm, amp, j, ra, n, m) result(b2)
        type(meridian_feec_t), intent(in) :: fs
        type(meridian_point_basis_t), intent(in) :: pb
        type(flux_mesh_t), intent(in) :: fm
        type(ampere_n_t), intent(inout) :: amp
        integer, intent(in) :: j, n, m
        type(ring_angles_t), intent(in) :: ra
        complex(dp), allocatable :: b2(:), f(:), a(:)
        integer :: st
        allocate (f(amp%nn), a(amp%nn))
        f = (0.0_dp, 0.0_dp)
        call ring_sheet_load(fs, pb, fm, j, ra, n, m, (1.0_dp, 0.0_dp), f)
        call ampere_n_solve(amp, f, a, st)
        if (st /= 0) error stop "solve"
        call ampere_n_curl(amp, a, b2)
    end function sheet_field

    subroutine cylinder_limit(m, n)
        integer, intent(in) :: m, n
        real(dp), parameter :: r0s(3) = [10.0_dp, 20.0_dp, 40.0_dp]
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ring_angles_t) :: ra_c, ra_s
        complex(dp), allocatable :: bc(:), bs(:)
        complex(dp) :: phic, phis, lam
        real(dp) :: err(3), xc, xs, ref, rate
        integer :: i
        character(len=40) :: tag

        do i = 1, 3
            r0 = r0s(i)
            call setup(fm, fs, pb, amp, n)
            call cylinder_angles(1.7_dp, ra_c)
            call cylinder_angles(real(m, dp)/n, ra_s)
            bc = sheet_field(fs, pb, fm, amp, jc, ra_c, n, m)
            bs = sheet_field(fs, pb, fm, amp, js, ra_s, n, m)
            phic = ring_flux_harmonic(fs, pb, fm, js, ra_s, n, bc, m)
            phis = ring_flux_harmonic(fs, pb, fm, js, ra_s, n, bs, m)
            lam = -phic/phis
            xc = n*fm%u(jc)/r0
            xs = n*fm%u(js)/r0
            ref = -xc*dbesk(m, xc)/(xs*dbesk(m, xs))
            err(i) = abs(lam - ref)/abs(ref)
            print '(a,i0,a,i0,a,f5.1,a,2es13.5,a,es10.2)', "  m=", m, " n=", n, &
                " R0=", r0, " lambda_s", lam, " rel. error", err(i)
        end do
        rate = log(err(1)/err(3))/log(4.0_dp)
        write (tag, '(a,i0,a,i0)') "m ", m, ", n ", n
        print '(a,a,f6.2)', "  ", trim(tag)//": order in a/R0", rate
        call record(rate > 1.7_dp, "cylinder limit at second order in a/R0, " &
            //trim(tag))
        call record(err(3) < 1.0e-4_dp, "cylinder limit within 1e-4 at R0 = 40, " &
            //trim(tag))
    end subroutine cylinder_limit

    subroutine reciprocity_and_shielding()
        integer, parameter :: n = 3, nsh = 3
        integer, parameter :: ring(nsh) = [10, 16, 26], mm(nsh) = [4, 5, 7]
        type(flux_mesh_t) :: fm
        type(meridian_feec_t) :: fs
        type(meridian_point_basis_t) :: pb
        type(ampere_n_t) :: amp
        type(ring_angles_t) :: ra(nsh), ra_d
        complex(dp) :: m(nsh, nsh), phid(nsh), lam(nsh), res(nsh)
        complex(dp), allocatable :: b2(:, :), bd(:), bt(:)
        real(dp) :: asym, scale
        integer :: r, s, ipiv(nsh), info

        r0 = 3.0_dp
        call setup(fm, fs, pb, amp, n)
        allocate (b2(fs%sp(2)%ndof, nsh))
        do s = 1, nsh
            call cylinder_angles(real(mm(s), dp)/n + 0.1_dp*s, ra(s))
            b2(:, s) = sheet_field(fs, pb, fm, amp, ring(s), ra(s), n, mm(s))
        end do
        do s = 1, nsh
            do r = 1, nsh
                m(r, s) = ring_flux_harmonic(fs, pb, fm, ring(r), ra(r), n, &
                    b2(:, s), mm(r))
            end do
        end do
        scale = maxval(abs(m))
        asym = maxval(abs(m - conjg(transpose(m))))/scale
        print '(a,es10.2)', "  reciprocity |M - M^H|/|M|", asym
        call record(asym < 1.0e-10_dp, "sheet response matrix Hermitian")
        call record(all([(real(m(r, r), dp) < 0.0_dp, r=1, nsh)]), &
            "sheet self-inductance negative (shielding sign)")
        ! drive: a sheet on ring 30 with m = 5
        call cylinder_angles(2.3_dp, ra_d)
        bd = sheet_field(fs, pb, fm, amp, 30, ra_d, n, 5)
        do r = 1, nsh
            phid(r) = ring_flux_harmonic(fs, pb, fm, ring(r), ra(r), n, bd, mm(r))
        end do
        lam = -phid
        call zgesv(nsh, 1, m, nsh, ipiv, lam, nsh, info)
        bt = bd + matmul(b2, lam)
        do r = 1, nsh
            res(r) = ring_flux_harmonic(fs, pb, fm, ring(r), ra(r), n, bt, mm(r))
        end do
        print '(a,es10.2)', "  shielded resonant flux / drive", &
            maxval(abs(res))/maxval(abs(phid))
        call record(info == 0 .and. maxval(abs(res)) < &
            1.0e-10_dp*maxval(abs(phid)), "total resonant flux shielded")
    end subroutine reciprocity_and_shielding


    end subroutine run_test_mhd_shielding

    subroutine run_test_mhd_shielding_circle_map(ctx, u, theta, x)
        type(run_test_mhd_shielding_context_t), intent(in) :: ctx
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        x = [ctx%r0 + u*cos(theta), u*sin(theta)]
    end subroutine run_test_mhd_shielding_circle_map

    subroutine run_test_mhd_shielding_circle_map_binding(self, u, theta, x)
        class(run_test_mhd_shielding_circle_map_t), intent(in) :: self
        real(dp), intent(in) :: u, theta
        real(dp), intent(out) :: x(2)
        call run_test_mhd_shielding_circle_map(self%data, u, theta, x)
    end subroutine
end module test_mhd_shielding_callbacks

program test_mhd_shielding
    use test_mhd_shielding_callbacks, only: run_test_mhd_shielding
    implicit none
    call run_test_mhd_shielding()
end program test_mhd_shielding
