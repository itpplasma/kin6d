program test_island_provider
    ! Closure A+D provider (kin6d_island_provider) on the helical-coil island
    ! case (test/support/island_a_fixture.f90).
    !   tangent k nr: implicit JVP against central differences of converged
    !     solves (independent of the linearisation), adjoint pairing of the
    !     VJP, at epsilon = 0 (island appearance) and epsilon = 0.5.
    !   nested k nr: epsilon = 0 equals the axisymmetric N = 0 NN2 solve and
    !     has no n >= 1 response; mode-0 deviation O(epsilon^2) and the
    !     n >= 1 response O(epsilon); parallel defect falls with the ratio.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use island_a_fixture, only: island_a_setup, shell_source
    use nn_vacuum_fixture, only: field, r0
    use kin6d_nonnested_coupled, only: nc_problem_t, nc_layout_t, nc_layout
    use kin6d_island_provider, only: island_provider_t, island_provider_init, &
        island_provider_solve, island_provider_sample, island_provider_jvp, &
        island_provider_vjp, island_provider_provenance, IP_NTHETA, IP_OK, &
        IP_EPSILON
    implicit none
    character(16) :: mode, arg
    integer :: k, nr
    call get_command_argument(1, mode)
    call get_command_argument(2, arg)
    read (arg, *) k
    call get_command_argument(3, arg)
    read (arg, *) nr
    select case (trim(mode))
    case ('tangent')
        call tangent_case()
    case ('nested')
        call nested_case()
    case default
        error stop 'mode'
    end select
contains
    function sample_points() result(xyz)
        real(dp), allocatable :: xyz(:, :)
        real(dp), parameter :: us(4) = [0.45_dp, 0.52_dp, 0.62_dp, 0.9_dp]
        real(dp), parameter :: phis(3) = [0.0_dp, 0.7_dp, 2.1_dp]
        real(dp), parameter :: ths(2) = [0.4_dp, 3.0_dp]
        real(dp) :: rr, zz
        integer :: i, j, l, m
        allocate (xyz(3, size(us)*size(phis)*size(ths)))
        m = 0
        do i = 1, size(us)
            do j = 1, size(ths)
                do l = 1, size(phis)
                    m = m + 1
                    rr = r0 + us(i)*cos(ths(j)); zz = us(i)*sin(ths(j))
                    xyz(:, m) = [rr*cos(phis(l)), rr*sin(phis(l)), zz]
                end do
            end do
        end do
    end function

    subroutine make(prov, eps, ratio, modes)
        type(island_provider_t), intent(inout) :: prov
        real(dp), intent(in) :: eps, ratio
        integer, intent(in) :: modes
        type(nc_problem_t) :: problem
        real(dp), allocatable :: base(:), pert(:)
        integer :: st
        call island_a_setup(nr, k, modes, base, pert, problem)
        call island_provider_init(prov, field, problem, base, pert, shell_source, &
            [1.0_dp, eps, 0.0_dp, ratio], 0.9_dp, 1.0e5_dp, 'helical-coil-annulus', st)
        if (st /= IP_OK) error stop 'init'
        call island_provider_solve(prov, prov%theta, st)
        if (st /= IP_OK) error stop 'solve'
    end subroutine

    subroutine tangent_case()
        type(island_provider_t) :: prov
        real(dp), allocatable :: xyz(:, :), jv(:, :), fd(:, :), vp(:, :), vm(:, :)
        real(dp), allocatable :: w(:, :)
        real(dp) :: theta0(IP_NTHETA), dth(IP_NTHETA), tbar(IP_NTHETA), h
        real(dp) :: err, scale, lhs, rhs, worst, floor
        real(dp), parameter :: eps0(2) = [0.0_dp, 0.5_dp], ratio = 1.0e2_dp
        real(dp), parameter :: typical(IP_NTHETA) = [1.0_dp, 1.0_dp, 1.0e6_dp, 1.0_dp]
        integer :: e, kk, st, c
        xyz = sample_points()
        allocate (jv(4, size(xyz, 2)), fd(4, size(xyz, 2)), vp(4, size(xyz, 2)), &
            vm(4, size(xyz, 2)), w(4, size(xyz, 2)))
        call random_seed(put=[(17 + c, c=1, 64)])
        call random_number(w)
        w = w - 0.5_dp
        w(1:3, :) = w(1:3, :)*1.0e3_dp ! comparable B and p contributions
        do e = 1, size(eps0)
            call make(prov, eps0(e), ratio, 1)
            print '(a)', island_provider_provenance(prov)
            print '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)', 'eps ', eps0(e), &
                ' pressure rows ', prov%report%pressure, ' eps_parallel ', &
                prov%report%epsilon_parallel, ' D3 force ', prov%report%force%normalized
            if (prov%report_status /= 0) error stop 'residual report'
            theta0 = prov%theta
            worst = 0
            do kk = 1, IP_NTHETA
                dth = 0; dth(kk) = 1
                call island_provider_jvp(prov, xyz, dth, jv, st)
                if (st /= IP_OK) error stop 'jvp'
                ! Steps on each datum's natural scale (I_phi in A).
                h = 1.0e-4_dp*max(typical(kk), abs(theta0(kk)))
                call island_provider_solve(prov, theta0 + h*dth, st, tol=1.0e-13_dp)
                if (st /= IP_OK) error stop 'fd plus'
                call island_provider_sample(prov, xyz, vp, st)
                call island_provider_solve(prov, theta0 - h*dth, st, tol=1.0e-13_dp)
                if (st /= IP_OK) error stop 'fd minus'
                call island_provider_sample(prov, xyz, vm, st)
                fd = (vp - vm)/(2*h)
                do c = 1, 2 ! B rows, then p
                    ! The difference quotient carries the solve's relative
                    ! rounding of the sampled value divided by h (floor 1e-11).
                    if (c == 1) then
                        scale = maxval(abs(fd(1:3, :)))
                        err = maxval(abs(jv(1:3, :) - fd(1:3, :)))
                        floor = 1.0e-11_dp*maxval(abs(vp(1:3, :)))/h
                    else
                        scale = maxval(abs(fd(4, :)))
                        err = maxval(abs(jv(4, :) - fd(4, :)))
                        floor = 1.0e-11_dp*maxval(abs(vp(4, :)))/h
                    end if
                    print '(a,es9.1,a,i0,a,a,es10.2,a,es10.2,a,es10.2)', 'eps ', &
                        eps0(e), ' theta_', kk, merge(' B', ' p', c == 1), &
                        ' |JVP-FD|/|FD| ', err/max(scale, tiny(1.0_dp)), ' |FD| ', &
                        scale, ' noise floor ', floor
                    worst = max(worst, err/max(scale, floor/1.0e-5_dp))
                end do
                call island_provider_solve(prov, theta0, st, tol=1.0e-13_dp)
                if (st /= IP_OK) error stop 'restore'
            end do
            if (worst > 1.0e-5_dp) error stop 'JVP disagrees with central differences'
            ! Adjoint pairing <w, J dtheta> = <J^T w, dtheta>.
            dth = [0.3_dp, -0.7_dp, 1.1_dp, 13.0_dp]
            call island_provider_jvp(prov, xyz, dth, jv, st)
            if (st /= IP_OK) error stop 'jvp pairing'
            call island_provider_vjp(prov, xyz, w, tbar, st)
            if (st /= IP_OK) error stop 'vjp'
            lhs = sum(w*jv); rhs = dot_product(tbar, dth)
            print '(a,es9.1,a,2es24.15,a,es10.2)', 'eps ', eps0(e), &
                ' pairing ', lhs, rhs, ' rel ', abs(lhs - rhs)/max(abs(lhs), abs(rhs))
            if (abs(lhs - rhs) > 1.0e-9_dp*sum(abs(tbar*dth))) &
                error stop 'adjoint pairing'
        end do
        print '(a)', 'PASS: closure A+D provider implicit JVP/VJP'
    end subroutine

    subroutine nested_case()
        type(island_provider_t) :: prov, axi
        type(nc_layout_t) :: lay, lay0
        real(dp), allocatable :: x0(:), dev0(:), devn(:)
        real(dp), parameter :: epss(3) = [0.1_dp, 0.05_dp, 0.025_dp]
        real(dp), parameter :: ratios(3) = [1.0e2_dp, 1.0e4_dp, 1.0e6_dp]
        real(dp) :: rel, scale, par(3), slope0, slopen
        integer :: i, st
        ! Existing nested path: the same NN2 closure with N = 0.
        call make(axi, 0.0_dp, ratios(1), 0)
        lay0 = nc_layout(axi%field)
        call make(prov, 0.0_dp, ratios(1), 1)
        lay = nc_layout(prov%field)
        x0 = prov%x
        scale = maxval(abs(axi%x))
        rel = max(maxval(abs(blocks(prov%x, lay, 0) - axi%x)), &
            maxval(abs(blocks(prov%x, lay, 1))))/scale
        print '(a,es10.2)', 'eps 0 vs N=0 solve and n>=1 response, rel ', rel
        if (rel > 1.0e-9_dp) error stop 'epsilon = 0 is not the nested N = 0 state'
        allocate (dev0(size(epss)), devn(size(epss)))
        do i = 1, size(epss)
            prov%theta(IP_EPSILON) = epss(i)
            call island_provider_solve(prov, prov%theta, st)
            if (st /= IP_OK) error stop 'eps solve'
            dev0(i) = maxval(abs(blocks(prov%x - x0, lay, 0)))/scale
            devn(i) = maxval(abs(blocks(prov%x, lay, 1)))/scale
            print '(a,es10.2,a,2es11.3)', 'eps ', epss(i), &
                ' |x0-x0(0)|, |x_n| ', dev0(i), devn(i)
        end do
        slope0 = log(dev0(2)/dev0(3))/log(2.0_dp)
        slopen = log(devn(2)/devn(3))/log(2.0_dp)
        print '(a,2f8.3)', 'orders mode 0, n >= 1 ', slope0, slopen
        if (abs(slope0 - 2) > 0.05_dp .or. abs(slopen - 1) > 0.05_dp) &
            error stop 'nested-limit orders'
        do i = 1, size(ratios)
            axi%theta(4) = ratios(i)
            call island_provider_solve(axi, axi%theta, st)
            if (st /= IP_OK) error stop 'ratio solve'
            par(i) = axi%report%epsilon_parallel
            print '(a,es9.1,a,es10.2,a,es10.2,a,es10.2)', 'nested ratio ', ratios(i), &
                ' eps_parallel ', par(i), ' pressure rows ', axi%report%pressure, &
                ' D3 ', axi%report%force%normalized
        end do
        if (par(2) >= par(1) .or. par(3) >= par(2)) &
            error stop 'parallel defect does not fall with the ratio'
        print '(a)', 'PASS: closure A+D nested limit'
    end subroutine

    function blocks(x, lay, which) result(y)
        ! which = 0: the mode-0 block of every field (N = 0 layout);
        ! which = 1: all n >= 1 blocks, concatenated.
        real(dp), intent(in) :: x(:)
        type(nc_layout_t), intent(in) :: lay
        integer, intent(in) :: which
        real(dp), allocatable :: y(:)
        integer :: f, off(4), sz(4)
        off = [lay%ia, lay%il, lay%ip, lay%iphi]
        sz = [lay%na1, lay%n0, lay%n0, lay%n0]
        allocate (y(0))
        do f = 1, 4
            if (which == 0) then
                y = [y, x(off(f) + 1:off(f) + sz(f))]
            else
                y = [y, x(off(f) + sz(f) + 1:off(f) + lay%nb*sz(f))]
            end if
        end do
        if (which == 0) y = [y, x(lay%iv)]
    end function
end program
