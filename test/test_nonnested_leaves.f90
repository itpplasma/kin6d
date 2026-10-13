program test_nonnested_leaves
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kin6d_equilibrium_core_leaves, only: maxwell_stress_delta, &
        maxwell_stress_delta_jvp, maxwell_stress_delta_vjp, stress_residual_density
    use kin6d_nonnested_leaves, only: segment_potential, segment_potential_jvp, &
        segment_potential_vjp, fieldline_velocity, fieldline_velocity_jvp, &
        fieldline_velocity_vjp, coupled_current, coupled_current_jvp, &
        coupled_current_vjp, anisotropic_flux
    use kin6d_coils_3d, only: segment_b, coil_mu0
    implicit none
    real(dp) :: x(6), dx(6), g(6), w(9), y(9), yp(9), ym(9), dy(9), c(5)
    real(dp) :: err, dot_error, h, a(3), l(3), p(3), pp(3), ap(3), am(3)
    real(dp) :: grad(3, 3), b(3), bref(3), factor, d, dfd, da, db, dl
    real(dp) :: vx(5), dvx(5), gv(5), v(2), vp(2), vm(2), dv(2), bv(2)
    real(dp) :: cx(12), cdx(12), cg(12), cj(3), cjp(3), cjm(3), cdj(3), cw(3)
    real(dp) :: q(3), bb(3), gpv(3), jxb(3)
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
