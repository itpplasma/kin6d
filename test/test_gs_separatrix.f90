program test_gs_separatrix
    !! Exact double-null Solov'ev member of the E2.1 generated CF family:
    !! psi=(R^2-1)^2+(5/4-R^2) Z^2, Delta*psi=6 R^2+5/2.
    !! Independent factorization: psi-1/16=(R^2-5/4)(R^2-3/4-Z^2).
    !! Axis (1,0), X-points (sqrt(5)/2,+/-1/sqrt(2)). No PDE solve
    !! or external fixture enters the interpolation/tracing oracle.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_set_num_threads, omp_get_wtime
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs, sample_gs
    use kin6d_gs_free_boundary, only: gs_fb_config_t, gs_fb_lcfs_levels, &
                                      gs_fb_trace_lcfs, GS_FB_XPOINT, GS_FB_OK
    use kin6d_cerfon_freidberg_leaves, only: cerfon_freidberg_jet
    use kin6d_gs_separatrix_leaves, only: gs_lcfs_ray_jet, &
                                          gs_lcfs_ray_jet_jvp, gs_lcfs_ray_jet_vjp
    implicit none
    integer, parameter :: npoints = 256
    real(dp), parameter :: exact_x(2, 2) = reshape([sqrt(1.25_dp), &
                                    -sqrt(0.5_dp), sqrt(1.25_dp), sqrt(0.5_dp)], [2, 2])
    type(gs_state_t) :: state, half
    type(gs_fb_config_t) :: cfg
    real(dp), allocatable :: chain(:, :), serial(:, :), xp(:, :), ga(:), gb(:)
    real(dp), allocatable :: contour(:, :)
    real(dp) :: pa, pb, bound, axis(2), xpoint(2), jet(6), value, br, bz, bt
    real(dp) :: flux_error(3), x_error(3), level_error, fe_gap, order, h, t1, t4
    integer :: degree, rung, nside, n, i, j, status, cell

    call ray_leaf_oracle()
    call omp_set_num_threads(4)
    do degree = 2, 3
        do rung = 1, 3
            nside = 24*2**(rung - 1)
            n = 4*nside
            allocate (contour(2, n))
            do i = 1, n
                call outer_curve(real(i - 1, dp)/n, contour(:, i))
                contour(2, i) = (contour(2, i) + 1)/2
            end do
            h = 2.0_dp/nside
            call initialize_gs(half, contour, 0.5_dp*h**2, status, &
                               element_order=degree)
            call check(status == 0, 'mesh')
            call mirror_mesh(half, state)
            do i = 1, size(state%psi)
                call exact_jet(state%points(:, i), jet)
                state%psi(i) = jet(1)
                call check(abs(jet(1) - exact_flux(state%points(:, i))) < 2e-14_dp, &
                           'independent exact Solovev polynomial')
            end do
            cfg%mode = GS_FB_XPOINT
            cfg%axis_sign = -1
            cfg%axis_guess = [1.0_dp, 0.0_dp]
            cfg%search_radius = 0.15_dp
            bound = huge(1.0_dp)
            do j = 1, 2
                cfg%xpoint_guess = exact_x(:, j)
                call gs_fb_lcfs_levels(state, state%psi, cfg, cfg%axis_guess, &
                           cfg%xpoint_guess, pa, pb, axis, xpoint, cell, ga, gb, status)
                call check(status == GS_FB_OK, 'Solovev critical points')
                bound = min(bound, pb)
            end do
            call gs_fb_trace_lcfs(state, state%psi, axis, bound, -1.0_dp, &
                                  npoints, chain, status, xpoints=xp)
            call check(status == GS_FB_OK, 'double-null LCFS trace')
            call check(size(xp, 2) == 2, 'both double-null X-points returned')
            x_error(rung) = 0
            do j = 1, size(xp, 2)
                x_error(rung) = max(x_error(rung), minval(norm2( &
                                              exact_x - spread(xp(:, j), 2, 2), dim=1)))
                call check(minval(norm2(chain - spread(xp(:, j), 2, npoints), &
                                dim=1)) < 2e-15_dp, 'X-point retained exactly in chain')
            end do
            flux_error(rung) = 0
            fe_gap = 0
            do i = 1, npoints
                flux_error(rung) = max(flux_error(rung), &
                                       abs(exact_flux(chain(:, i)) - 0.0625_dp))
                call sample_gs(state, chain(1, i), chain(2, i), 1.0_dp, &
                               value, br, bz, bt, status)
                call check(status == 0, 'chain FE sample')
                fe_gap = max(fe_gap, abs(value - bound))
                call check(abs(chain(2, i)) <= sqrt(0.5_dp) + 4*h**degree, &
                           'closed component excludes private flux')
            end do
            level_error = abs(bound - 0.0625_dp)
            print '(a,3i7,4es15.6)', 'Solovev p,rung,nx, flux,x,level,FE ', &
                degree, rung, size(xp, 2), flux_error(rung), x_error(rung), &
                level_error, fe_gap
            call check(fe_gap < 2e-13_dp, 'FE flux root near roundoff')
            if (rung == 3) then
                call omp_set_num_threads(1)
                t1 = omp_get_wtime()
                call gs_fb_trace_lcfs(state, state%psi, axis, bound, -1.0_dp, &
                                      npoints, serial, status)
                t1 = omp_get_wtime() - t1
                call check(status == GS_FB_OK, 'serial separatrix')
                call omp_set_num_threads(4)
                t4 = omp_get_wtime()
                call gs_fb_trace_lcfs(state, state%psi, axis, bound, -1.0_dp, &
                                      npoints, chain, status)
                t4 = omp_get_wtime() - t4
                call check(status == GS_FB_OK, 'parallel separatrix')
                call check(maxval(abs(serial - chain)) < 2e-13_dp, 'thread parity')
                print '(a,i2,3es15.6)', 'trace p, 1-thread s, 4-thread s, speedup ', &
                    degree, t1, t4, t1/t4
                ! A regular limited contour has no bounding saddle. Its old
                ! fixed-step/normal-projection path remains exact for P2/P3.
                state%psi = (state%points(1, :) - 1)**2 + state%points(2, :)**2
                call gs_fb_trace_lcfs(state, state%psi, [1.0_dp, 0.0_dp], &
                                   0.04_dp, -1.0_dp, npoints, chain, status, xpoints=xp)
                call check(status == GS_FB_OK, 'limited circle trace')
                call check(size(xp, 2) == 0, 'limited circle has no X-point')
                call check(maxval(abs(norm2(chain - spread([1.0_dp, 0.0_dp], &
                                            2, npoints), dim=1) - 0.2_dp)) < 2e-14_dp, &
                           'limited circle unchanged')
            end if
            deallocate (contour)
        end do
        do rung = 2, 3
            order = log(flux_error(rung - 1)/flux_error(rung))/log(2.0_dp)
            print '(a,2i3,f10.4)', 'flux observed order p,rung ', degree, rung, order
        end do
        call check(flux_error(3) < flux_error(1)/2.0_dp**(2*degree), &
                   'grid-order flux convergence')
        call check(x_error(3) < x_error(1)/2.0_dp**(2*(degree - 1)), &
                   'critical position convergence')
    end do
    print '(a)', 'verified candidate contours; continuous bounds unavailable'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (ok) return
        print '(2a)', 'FAIL: ', message
        error stop 1
    end subroutine check

    pure function exact_flux(point) result(flux)
        real(dp), intent(in) :: point(2)
        real(dp) :: flux
        flux = (point(1)**2 - 1)**2 + (1.25_dp - point(1)**2)*point(2)**2
    end function exact_flux

    subroutine exact_jet(point, jet)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: jet(6)
        call cerfon_freidberg_jet(point(1), point(2), 5.0_dp/17, &
                          2.0_dp/17, -4.0_dp/17, 5.0_dp/34, 1.0_dp/34, 0.0_dp, 0.0_dp, &
                                 0.0_dp, jet(1), jet(2), jet(3), jet(4), jet(5), jet(6))
        jet = 8.5_dp*jet
    end subroutine exact_jet

    subroutine outer_curve(t, point)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: point(2)
        real(dp) :: f
        f = 4*t
        if (f < 1) then
            point = [0.65_dp + 0.75_dp*f, -1.0_dp]
        else if (f < 2) then
            point = [1.4_dp, -1.0_dp + 2*(f - 1)]
        else if (f < 3) then
            point = [1.4_dp - 0.75_dp*(f - 2), 1.0_dp]
        else
            point = [0.65_dp, 1.0_dp - 2*(f - 3)]
        end if
    end subroutine outer_curve

    subroutine mirror_mesh(half, whole)
        ! Reflect the native upper-half mesh to keep the two saddle levels
        ! identical. Unrelated upper/lower triangulation errors would split
        ! this double-null field into two nearby single-null levels.
        type(gs_state_t), intent(in) :: half
        type(gs_state_t), intent(out) :: whole
        integer :: n, ne, ndof
        n = size(half%psi)
        ne = size(half%area)
        ndof = size(half%triangles, 1)
        whole%element_order = half%element_order
        allocate (whole%points(2, 2*n), whole%psi(2*n), whole%boundary(2*n))
        allocate (whole%triangles(ndof, 2*ne), whole%gradient(2, 3, 2*ne))
        allocate (whole%area(2*ne), whole%curved_cell(2*ne))
        whole%points(:, :n) = half%points
        whole%points(:, n + 1:) = half%points
        whole%points(2, n + 1:) = -half%points(2, :)
        whole%triangles(:, :ne) = half%triangles
        whole%triangles(:, ne + 1:) = half%triangles + n
        whole%gradient(:, :, :ne) = half%gradient
        whole%gradient(:, :, ne + 1:) = half%gradient
        whole%gradient(2, :, ne + 1:) = -half%gradient(2, :, :)
        whole%area(:ne) = half%area
        whole%area(ne + 1:) = half%area
        whole%boundary(:n) = half%boundary
        whole%boundary(n + 1:) = half%boundary
        whole%curved_cell = .false.
        whole%ready = .true.
    end subroutine mirror_mesh

    subroutine ray_leaf_oracle()
        real(dp) :: gap, slope, dg, ds, adj(7)
        call gs_lcfs_ray_jet(0.7_dp, 0.2_dp, -1.0_dp, 2.0_dp, &
                             -3.0_dp, 0.6_dp, 0.8_dp, gap, slope)
        call check(abs(gap + 0.5_dp) < 1e-15_dp, 'ray gap independent oracle')
        call check(abs(slope - 1.2_dp) < 1e-15_dp, 'ray slope independent oracle')
        call gs_lcfs_ray_jet_jvp(0.7_dp, 0.2_dp, -1.0_dp, 2.0_dp, &
                              -3.0_dp, 0.6_dp, 0.8_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                                 5.0_dp, 6.0_dp, 7.0_dp, dg, ds)
        call check(abs(dg - 2.5_dp) < 1e-14_dp, 'ray JVP independent gap')
        call check(abs(ds - (-1.0_dp)) < 1e-14_dp, 'ray JVP independent slope')
        call gs_lcfs_ray_jet_vjp(0.7_dp, 0.2_dp, -1.0_dp, 2.0_dp, &
                              -3.0_dp, 0.6_dp, 0.8_dp, 0.3_dp, 0.4_dp, adj(1), adj(2), &
                                 adj(3), adj(4), adj(5), adj(6), adj(7))
        call check(abs(sum(adj*[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                 6.0_dp, 7.0_dp]) - 0.3_dp*dg - 0.4_dp*ds) < 1e-14_dp, 'ray AD pairing')
    end subroutine ray_leaf_oracle
end program test_gs_separatrix
