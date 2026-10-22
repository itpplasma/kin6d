program test_nonnested_leaves
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_core_leaves, only: maxwell_stress_delta, &
        maxwell_stress_delta_jvp, maxwell_stress_delta_vjp, stress_residual_density
    use kin6d_nonnested_leaves, only: segment_potential, segment_potential_jvp, &
        segment_potential_vjp, fieldline_velocity, fieldline_velocity_jvp, &
        fieldline_velocity_vjp, coupled_current, coupled_current_jvp, &
        coupled_current_vjp, anisotropic_flux, island_defect_density, &
        island_defect_density_jvp, island_defect_density_vjp, island_defect_norm, &
        island_defect_norm_jvp, island_defect_norm_vjp
    use kin6d_coils_3d, only: segment_b, coil_mu0
    implicit none
    real(dp) :: x(6), dx(6), g(6), w(9), y(9), yp(9), ym(9), dy(9), c(5)
    real(dp) :: err, dot_error, h, a(3), l(3), p(3), pp(3), ap(3), am(3)
    real(dp) :: grad(3, 3), b(3), bref(3), factor, d, dfd, da, db, dl
    real(dp) :: vx(5), dvx(5), gv(5), v(2), vp(2), vm(2), dv(2), bv(2)
    real(dp) :: cx(12), cdx(12), cg(12), cj(3), cjp(3), cjm(3), cdj(3), cw(3)
    real(dp) :: q(3), bb(3), gpv(3), jxb(3)
    real(dp) :: ix(10), idx(10), ig(10), iy(2), iyp(2), iym(2), idy(2)
    real(dp) :: nx(5), ndx(5), ng(5), cross(3), reference(2)
    integer :: i, j
    h = 1.0e-5_dp
    do i = 1, 4
        x = [0.2_dp*i, -0.07_dp*i, 0.3_dp, 8.0_dp, 0.17_dp, 1.3_dp]
        dx = [0.5_dp, 0.2_dp, -0.4_dp, 0.3_dp, 0.13_dp, -0.1_dp]
        w = [(sin(real(j, dp)), j=1, 9)]
        call stress(x, y)
        call stress(x + h*dx, yp)
        call stress(x - h*dx, ym)
        call maxwell_stress_delta_jvp(x(1), x(2), x(3), x(4), x(5), x(6), &
            dx(1), dx(2), dx(3), dx(4), dx(5), dx(6), dy(1), dy(2), dy(3), &
            dy(4), dy(5), dy(6), dy(7), dy(8), dy(9))
        call maxwell_stress_delta_vjp(x(1), x(2), x(3), x(4), x(5), x(6), &
            w(1), w(2), w(3), w(4), w(5), w(6), w(7), w(8), w(9), &
            g(1), g(2), g(3), g(4), g(5), g(6))
        err = maxval(abs((yp - ym)/(2*h) - dy))
        dot_error = abs(dot_product(w, dy) - dot_product(g, dx))
        if (err > 1.0e-8_dp .or. dot_error > 1.0e-12_dp) error stop 'D3 derivatives'
        call stress_residual_density(6.2_dp, x(1), x(3), x(4) + x(2), x(4), &
            x(5), x(6), c(1), c(2), c(3), c(4), c(5))
        err = maxval(abs(c - [y(5), 6.2_dp*y(1), 6.2_dp*y(3), &
            6.2_dp*y(7), 6.2_dp*y(9)]))
        if (err > 1.0e-12_dp) error stop 'one D3 residual'
    end do
    print '(a)', 'PASS: 3D D3 JVP/VJP and axisymmetric contraction'
    ! I1 independent dot/cross-product oracle plus generated JVP/VJP checks.
    ix = [0.3_dp, 2.1_dp, -0.4_dp, 0.7_dp, 0.2_dp, -0.6_dp, &
        0.1_dp, -0.3_dp, 0.4_dp, 0.5_dp]
    idx = [(sin(real(i, dp)), i=1, 10)]
    call defects(ix, iy)
    cross = [ix(2)*ix(9) - ix(3)*ix(8), ix(3)*ix(7) - ix(1)*ix(9), &
        ix(1)*ix(8) - ix(2)*ix(7)]
    reference = [dot_product(ix(1:3), ix(4:6))**2/sum(ix(1:3)**2), &
        ix(10)**2*sum(cross**2)]
    if (maxval(abs(iy - reference)) > 1.0e-13_dp) error stop 'I1 density oracle'
    call defects(ix + h*idx, iyp)
    call defects(ix - h*idx, iym)
    call island_defect_density_jvp(ix(1), ix(2), ix(3), ix(4), ix(5), ix(6), &
        ix(7), ix(8), ix(9), ix(10), idx(1), idx(2), idx(3), idx(4), idx(5), &
        idx(6), idx(7), idx(8), idx(9), idx(10), idy(1), idy(2))
    call island_defect_density_vjp(ix(1), ix(2), ix(3), ix(4), ix(5), ix(6), &
        ix(7), ix(8), ix(9), ix(10), 0.7_dp, -0.3_dp, ig(1), ig(2), ig(3), &
        ig(4), ig(5), ig(6), ig(7), ig(8), ig(9), ig(10))
    if (maxval(abs(idy - (iyp - iym)/(2*h))) > 1.0e-8_dp) error stop 'I1 density JVP'
    if (abs(0.7_dp*idy(1) - 0.3_dp*idy(2) - dot_product(idx, ig)) &
        > 1.0e-12_dp) error stop 'I1 density VJP'
    nx = [0.7_dp, 0.3_dp, 2.1_dp, 1.2_dp, 0.4_dp]
    ndx = [0.2_dp, -0.1_dp, 0.3_dp, 0.5_dp, -0.3_dp]
    call norms(nx, iy)
    reference = [nx(4)*sqrt(nx(1)/nx(3))/nx(5), sqrt(nx(2)/nx(3))]
    if (maxval(abs(iy - reference)) > 1.0e-13_dp) error stop 'I1 norm oracle'
    call norms(nx + h*ndx, iyp)
    call norms(nx - h*ndx, iym)
    call island_defect_norm_jvp(nx(1), nx(2), nx(3), nx(4), nx(5), ndx(1), &
        ndx(2), ndx(3), ndx(4), ndx(5), idy(1), idy(2))
    call island_defect_norm_vjp(nx(1), nx(2), nx(3), nx(4), nx(5), &
        0.7_dp, -0.3_dp, ng(1), ng(2), ng(3), ng(4), ng(5))
    if (maxval(abs(idy - (iyp - iym)/(2*h))) > 1.0e-8_dp) error stop 'I1 norm JVP'
    if (abs(0.7_dp*idy(1) - 0.3_dp*idy(2) - dot_product(ndx, ng)) &
        > 1.0e-12_dp) error stop 'I1 norm VJP'
    print '(a)', 'PASS: I1 physical defects, JVP/VJP and normalization'
    a = [1.0_dp, -0.4_dp, 0.3_dp]; l = [0.4_dp, 0.7_dp, -0.2_dp]
    p = [0.3_dp, 0.1_dp, 1.2_dp]
    do j = 1, 3
        pp = p; pp(j) = pp(j) + h
        call potential(pp, ap)
        pp(j) = pp(j) - 2*h
        call potential(pp, am)
        grad(:, j) = (ap - am)/(2*h)
    end do
    b = [grad(3, 2) - grad(2, 3), grad(1, 3) - grad(3, 1), &
        grad(2, 1) - grad(1, 2)]
    bref = segment_b(a, l, p)
    err = norm2(b - bref)/norm2(bref)
    print '(a,es12.4)', 'segment-potential curl versus Biot-Savart ', err
    if (err > 1.0e-8_dp) error stop 'potential curl'
    call segment_potential(2.0_dp, 2.1_dp, 0.8_dp, factor)
    call segment_potential_jvp(2.0_dp, 2.1_dp, 0.8_dp, 0.2_dp, -0.3_dp, &
        0.1_dp, d)
    call segment_potential(2 + h*0.2_dp, 2.1_dp - h*0.3_dp, 0.8_dp + h*0.1_dp, da)
    call segment_potential(2 - h*0.2_dp, 2.1_dp + h*0.3_dp, 0.8_dp - h*0.1_dp, db)
    dfd = (da - db)/(2*h)
    call segment_potential_vjp(2.0_dp, 2.1_dp, 0.8_dp, 0.7_dp, da, db, dl)
    if (abs(d - dfd) > 1.0e-9_dp) error stop 'potential JVP'
    if (abs(0.7_dp*d - (0.2_dp*da - 0.3_dp*db + 0.1_dp*dl)) > 1.0e-13_dp) &
        error stop 'potential VJP'
    vx = [6.2_dp, 0.4_dp, 0.3_dp, -0.2_dp, 5.3_dp]
    dvx = [0.2_dp, -0.4_dp, 0.1_dp, 0.3_dp, -0.1_dp]; bv = [0.7_dp, -0.3_dp]
    call velocity(vx, v)
    call velocity(vx + h*dvx, vp)
    call velocity(vx - h*dvx, vm)
    call fieldline_velocity_jvp(vx(1), vx(2), vx(3), vx(4), vx(5), &
        dvx(1), dvx(2), dvx(3), dvx(4), dvx(5), dv(1), dv(2))
    call fieldline_velocity_vjp(vx(1), vx(2), vx(3), vx(4), vx(5), bv(1), bv(2), &
        gv(1), gv(2), gv(3), gv(4), gv(5))
    if (maxval(abs((vp - vm)/(2*h) - dv)) > 1.0e-9_dp) error stop 'velocity JVP'
    if (abs(dot_product(bv, dv) - dot_product(gv, dvx)) > 1.0e-13_dp) &
        error stop 'velocity VJP'
    print '(a)', 'PASS: potential and field-line derivatives'
    ! NN2 coupled current: JVP, VJP pairing, J x B = grad_perp p with the
    ! resistive terms off, and the parallel flux chi_par b.grad p.
    cx = [0.3_dp, 5.1_dp, -0.4_dp, 0.7_dp, 0.05_dp, -0.2_dp, 0.11_dp, &
        -0.03_dp, 0.08_dp, -1.1_dp, 2.5_dp, 0.01_dp]
    cdx = [(cos(1.7_dp*j), j=1, 12)]; cw = [0.4_dp, -0.9_dp, 0.25_dp]
    call current(cx, cj)
    call current(cx + h*cdx, cjp)
    call current(cx - h*cdx, cjm)
    call coupled_current_jvp(cx(1), cx(2), cx(3), cx(4), cx(5), cx(6), cx(7), &
        cx(8), cx(9), cx(10), cx(11), cx(12), cdx(1), cdx(2), cdx(3), cdx(4), &
        cdx(5), cdx(6), cdx(7), cdx(8), cdx(9), cdx(10), cdx(11), cdx(12), &
        cdj(1), cdj(2), cdj(3))
    call coupled_current_vjp(cx(1), cx(2), cx(3), cx(4), cx(5), cx(6), cx(7), &
        cx(8), cx(9), cx(10), cx(11), cx(12), cw(1), cw(2), cw(3), cg(1), cg(2), &
        cg(3), cg(4), cg(5), cg(6), cg(7), cg(8), cg(9), cg(10), cg(11), cg(12))
    if (maxval(abs((cjp - cjm)/(2*h) - cdj)) > 1.0e-8_dp) error stop 'current JVP'
    if (abs(dot_product(cw, cdj) - dot_product(cg, cdx)) > 1.0e-12_dp) &
        error stop 'current VJP'
    bb = cx(1:3); gpv = cx(4:6)
    call coupled_current(bb(1), bb(2), bb(3), gpv(1), gpv(2), gpv(3), 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, cx(11), cx(12), cj(1), cj(2), cj(3))
    jxb = [cj(2)*bb(3) - cj(3)*bb(2), cj(3)*bb(1) - cj(1)*bb(3), &
        cj(1)*bb(2) - cj(2)*bb(1)]
    err = maxval(abs(jxb - (gpv - bb*dot_product(bb, gpv)/dot_product(bb, bb))))
    if (err > 1.0e-14_dp) error stop 'J x B'
    call anisotropic_flux(bb(1), bb(2), bb(3), gpv(1), gpv(2), gpv(3), 7.0_dp, &
        0.5_dp, q(1), q(2), q(3))
    if (abs(dot_product(q, bb) - 7*dot_product(gpv, bb)) > 1.0e-13_dp) &
        error stop 'parallel flux'
    print '(a)', 'PASS: NN2 coupled current and flux leaves'
contains
    subroutine defects(v, y)
        real(dp), intent(in) :: v(10)
        real(dp), intent(out) :: y(2)
        call island_defect_density(v(1), v(2), v(3), v(4), v(5), v(6), v(7), &
            v(8), v(9), v(10), y(1), y(2))
    end subroutine
    subroutine norms(v, y)
        real(dp), intent(in) :: v(5)
        real(dp), intent(out) :: y(2)
        call island_defect_norm(v(1), v(2), v(3), v(4), v(5), y(1), y(2))
    end subroutine
    subroutine current(x, y)
        real(dp), intent(in) :: x(12)
        real(dp), intent(out) :: y(3)
        call coupled_current(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            x(9), x(10), x(11), x(12), y(1), y(2), y(3))
    end subroutine
    subroutine stress(x, y)
        real(dp), intent(in) :: x(6)
        real(dp), intent(out) :: y(9)
        call maxwell_stress_delta(x(1), x(2), x(3), x(4), x(5), x(6), &
            y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), y(9))
    end subroutine
    subroutine potential(p, value)
        real(dp), intent(in) :: p(3)
        real(dp), intent(out) :: value(3)
        real(dp) :: f
        call segment_potential(norm2(p - a), norm2(p - a - l), norm2(l), f)
        value = coil_mu0/(4*acos(-1.0_dp))*f*l/norm2(l)
    end subroutine
    subroutine velocity(x, y)
        real(dp), intent(in) :: x(5)
        real(dp), intent(out) :: y(2)
        call fieldline_velocity(x(1), x(2), x(3), x(4), x(5), y(1), y(2))
    end subroutine
end program
