module kin6d_vrmhd_fourier
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use milca_gs_symmetry, only: gs_symmetry_t, GS_ROTATION, GS_SCREW
    use kin6d_vrmhd_parameters, only: vrmhd_parameters_t
    use kin6d_vrmhd_point, only: vrmhd_point, vrmhd_point_matrix
    use kin6d_vrmhd_section, only: vrmhd_section_t, vrmhd_section_build, &
        vrmhd_section_basis, vrmhd_section_residual, vrmhd_section_test
    use kin6d_vrmhd_section_leaves, only: vrmhd_orbit_grad_rotation, &
        vrmhd_orbit_grad_screw, vrmhd_orbit_velocity_rotation, &
        vrmhd_orbit_velocity_screw, vrmhd_orbit_curl_rotation, &
        vrmhd_orbit_curl_screw
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t, status_ok
    implicit none
    private
    ! Harmonic block b holds the section layout once for n = 0 and twice,
    ! cos(kappa t) then sin(kappa t), otherwise; kappa = n/orbit_scale.
    type, public :: vrmhd_fourier_t
        type(vrmhd_section_t) :: section
        real(dp) :: orbit_scale = 1
        integer :: ndof = 0
        integer, allocatable :: harmonic(:), offset(:), copies(:)
        logical, allocatable :: fixed(:)
    end type
    public :: vrmhd_fourier_build, vrmhd_fourier_data, vrmhd_fourier_residual
    public :: vrmhd_fourier_jacobian, vrmhd_fourier_wavenumber
contains
    subroutine vrmhd_fourier_build(mm, symmetry, degree, harmonics, orbit_scale, &
            op, status, scalar_degree, velocity_degree, coordinate_coframe)
        type(meridian_mesh_t), intent(in) :: mm
        type(gs_symmetry_t), intent(in) :: symmetry
        integer, intent(in) :: degree, harmonics(:)
        real(dp), intent(in) :: orbit_scale
        type(vrmhd_fourier_t), intent(out) :: op
        integer, intent(out) :: status
        integer, intent(in), optional :: scalar_degree, velocity_degree
        logical, intent(in), optional :: coordinate_coframe
        logical :: closed
        integer :: b, c, ns, nh
        status = 2
        nh = size(harmonics)
        if (nh < 1 .or. any(harmonics < 0) .or. orbit_scale <= 0) return
        do b = 2, nh
            if (any(harmonics(:b - 1) == harmonics(b))) return
        end do
        ! The rotation angle is 2 pi periodic; other scales change the measure.
        if (symmetry%kind == GS_ROTATION .and. abs(orbit_scale - 1) > 0) return
        closed = .false.
        if (present(coordinate_coframe)) closed = coordinate_coframe
        ! Nonzero twist requires closed coordinate forms for Fourier gradients.
        if (symmetry%kind == GS_SCREW .and. abs(symmetry%k) > 0 &
            .and. any(harmonics /= 0) .and. .not. closed) return
        call vrmhd_section_build(mm, symmetry, degree, 1, op%section, status, &
            scalar_degree=scalar_degree, velocity_degree=velocity_degree, &
            coordinate_coframe=closed)
        if (status /= 0) return
        ns = op%section%ndof
        op%orbit_scale = orbit_scale
        allocate (op%harmonic(nh), op%offset(nh), op%copies(nh))
        op%harmonic = harmonics
        op%copies = merge(1, 2, harmonics == 0)
        op%ndof = 0
        do b = 1, nh
            op%offset(b) = op%ndof
            op%ndof = op%ndof + op%copies(b)*ns
        end do
        allocate (op%fixed(op%ndof))
        do b = 1, nh
            do c = 1, op%copies(b)
                op%fixed(op%offset(b) + (c - 1)*ns + 1: &
                    op%offset(b) + c*ns) = op%section%fixed
            end do
        end do
    end subroutine

    pure real(dp) function vrmhd_fourier_wavenumber(op, b) result(kappa)
        type(vrmhd_fourier_t), intent(in) :: op
        integer, intent(in) :: b
        kappa = op%harmonic(b)/op%orbit_scale
    end function

    ! Section data (basis) and the orbit-derivative data per unit kappa.
    subroutine vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, status)
        type(vrmhd_fourier_t), intent(in) :: op
        integer, intent(in) :: t, q
        real(dp), intent(out) :: x(2), weight
        real(dp), intent(out) :: basis(36, op%section%nlocal)
        real(dp), intent(out) :: orbit(36, op%section%nlocal)
        integer, intent(out) :: status
        real(dp) :: k
        integer :: i
        call vrmhd_section_basis(op%section, t, q, x, weight, basis, status)
        if (status /= 0) return
        k = op%section%symmetry%k
        orbit = 0
        do i = 1, op%section%nlocal
            if (op%section%symmetry%kind == GS_ROTATION) then
                call vrmhd_orbit_grad_rotation(x(1), x(2), k, basis(1, i), &
                    orbit(13, i), orbit(14, i), orbit(15, i))
                call vrmhd_orbit_grad_rotation(x(1), x(2), k, basis(11, i), &
                    orbit(28, i), orbit(29, i), orbit(30, i))
                call vrmhd_orbit_grad_rotation(x(1), x(2), k, basis(12, i), &
                    orbit(31, i), orbit(32, i), orbit(33, i))
                call vrmhd_orbit_velocity_rotation(x(1), x(2), k, basis(2, i), &
                    basis(3, i), basis(4, i), orbit(16, i), orbit(17, i), &
                    orbit(18, i), orbit(19, i), orbit(20, i), orbit(21, i), &
                    orbit(22, i), orbit(23, i), orbit(24, i))
                call vrmhd_orbit_curl_rotation(x(1), x(2), k, basis(5, i), &
                    basis(6, i), basis(7, i), orbit(25, i), orbit(26, i), orbit(27, i))
                call vrmhd_orbit_curl_rotation(x(1), x(2), k, basis(8, i), &
                    basis(9, i), basis(10, i), orbit(34, i), orbit(35, i), &
                    orbit(36, i))
            else
                call vrmhd_orbit_grad_screw(x(1), x(2), k, basis(1, i), &
                    orbit(13, i), orbit(14, i), orbit(15, i))
                call vrmhd_orbit_grad_screw(x(1), x(2), k, basis(11, i), &
                    orbit(28, i), orbit(29, i), orbit(30, i))
                call vrmhd_orbit_grad_screw(x(1), x(2), k, basis(12, i), &
                    orbit(31, i), orbit(32, i), orbit(33, i))
                call vrmhd_orbit_velocity_screw(x(1), x(2), k, basis(2, i), &
                    basis(3, i), basis(4, i), orbit(16, i), orbit(17, i), &
                    orbit(18, i), orbit(19, i), orbit(20, i), orbit(21, i), &
                    orbit(22, i), orbit(23, i), orbit(24, i))
                call vrmhd_orbit_curl_screw(x(1), x(2), k, basis(5, i), &
                    basis(6, i), basis(7, i), orbit(25, i), orbit(26, i), orbit(27, i))
                call vrmhd_orbit_curl_screw(x(1), x(2), k, basis(8, i), &
                    basis(9, i), basis(10, i), orbit(34, i), orbit(35, i), &
                    orbit(36, i))
            end if
        end do
    end subroutine

    ! Linearised residual about the invariant state u0: block n = 0 carries
    ! the section residual R(u0) and its JVP; n > 0 blocks carry only the JVP.
    subroutine vrmhd_fourier_residual(op, p, u0, ud0, du, dud, source, r, jv, status)
        type(vrmhd_fourier_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u0(op%section%ndof), ud0(op%section%ndof)
        real(dp), intent(in) :: du(op%ndof), dud(op%ndof)
        real(dp), intent(in) :: source(6, op%section%fs%nq, &
            op%section%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: r(op%ndof), jv(op%ndof)
        integer, intent(out) :: status
        real(dp), allocatable :: r0(:), j0(:)
        integer :: b, ns, lo
        ns = op%section%ndof
        r = 0; jv = 0; status = 0
        allocate (r0(ns), j0(ns))
        do b = 1, size(op%harmonic)
            lo = op%offset(b)
            if (op%copies(b) == 1) then
                call vrmhd_section_residual(op%section, p, u0, ud0, &
                    du(lo + 1:lo + ns), dud(lo + 1:lo + ns), source, r0, j0, status)
                if (status /= 0) return
                r(lo + 1:lo + ns) = r0; jv(lo + 1:lo + ns) = j0
            else
                call harmonic_jvp(op, p, b, u0, ud0, du(lo + 1:lo + 2*ns), &
                    dud(lo + 1:lo + 2*ns), source, jv(lo + 1:lo + 2*ns), status)
                if (status /= 0) return
            end if
        end do
    end subroutine

    subroutine harmonic_jvp(op, p, b, u0, ud0, du, dud, source, jv, status)
        type(vrmhd_fourier_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        integer, intent(in) :: b
        real(dp), intent(in) :: u0(op%section%ndof), ud0(op%section%ndof)
        real(dp), intent(in) :: du(2*op%section%ndof), dud(2*op%section%ndof)
        real(dp), intent(in) :: source(6, op%section%fs%nq, &
            op%section%fs%mm%topology%n_triangles)
        real(dp), intent(out) :: jv(2*op%section%ndof)
        integer, intent(out) :: status
        real(dp) :: basis(36, op%section%nlocal), orbit(36, op%section%nlocal)
        real(dp) :: x(2), weight, v(36), vd(12), dc(36), ds(36), ddc(12), dds(12)
        real(dp) :: tp(19), tq(19), lr(2*op%section%nlocal), kappa, h1, h2, value
        real(dp) :: uc(op%section%nlocal), us(op%section%nlocal)
        integer :: t, q, i, st, fail, ns, nl, row
        ns = op%section%ndof; nl = op%section%nlocal
        kappa = vrmhd_fourier_wavenumber(op, b)
        jv = 0; fail = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(t,q,i,st,basis,orbit,x,weight,v,vd,dc,ds,ddc,dds,tp,tq, &
        !$omp lr,h1,h2,value,uc,us,row) reduction(max:fail)
        do t = 1, op%section%fs%mm%topology%n_triangles
            lr = 0
            do q = 1, op%section%fs%nq
                call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                if (st /= 0) then
                    fail = 1; cycle
                end if
                v = matmul(basis, u0(op%section%map(:, t)))
                vd = matmul(basis(:12, :), ud0(op%section%map(:, t)))
                uc = du(op%section%map(:, t)); us = du(ns + op%section%map(:, t))
                dc = matmul(basis, uc) + kappa*matmul(orbit, us)
                ds = matmul(basis, us) - kappa*matmul(orbit, uc)
                ddc = matmul(basis(:12, :), dud(op%section%map(:, t)))
                dds = matmul(basis(:12, :), dud(ns + op%section%map(:, t)))
                do i = 1, nl
                    row = op%section%row(i)
                    call vrmhd_section_test(row, basis(:, i), tp)
                    call vrmhd_section_test(row, orbit(:, i), tq)
                    call vrmhd_point(p, row, v(:33), dc(:33), vd, ddc, tp, &
                        source(:, q, t), value, h1)
                    call vrmhd_point(p, row, v(:33), ds(:33), vd, dds, &
                        -kappa*tq, source(:, q, t), value, h2)
                    lr(i) = lr(i) + 0.5_dp*weight*(h1 + h2)
                    call vrmhd_point(p, row, v(:33), dc(:33), vd, ddc, &
                        kappa*tq, source(:, q, t), value, h1)
                    call vrmhd_point(p, row, v(:33), ds(:33), vd, dds, tp, &
                        source(:, q, t), value, h2)
                    lr(nl + i) = lr(nl + i) + 0.5_dp*weight*(h1 + h2)
                end do
            end do
            do i = 1, nl
                !$omp atomic update
                jv(op%section%map(i, t)) = jv(op%section%map(i, t)) + lr(i)
                !$omp atomic update
                jv(ns + op%section%map(i, t)) = jv(ns + op%section%map(i, t)) &
                    + lr(nl + i)
            end do
        end do
        !$omp end parallel do
        status = fail
    end subroutine

    ! Sparse R_U + (M U_dot)_U + shift M, block diagonal in the harmonics.
    ! Orbit averages of cos/sin products give the factor 1/2 for n > 0.
    subroutine vrmhd_fourier_jacobian(op, p, u0, ud0, source, shift, jac, status)
        type(vrmhd_fourier_t), intent(in) :: op
        type(vrmhd_parameters_t), intent(in) :: p
        real(dp), intent(in) :: u0(op%section%ndof), ud0(op%section%ndof), shift
        real(dp), intent(in) :: source(6, op%section%fs%nq, &
            op%section%fs%mm%topology%n_triangles)
        type(csc_t), intent(out) :: jac
        integer, intent(out) :: status
        real(dp) :: basis(36, op%section%nlocal), orbit(36, op%section%nlocal)
        real(dp) :: x(2), weight, v(36), vd(12), kmat(19, 45, 6)
        real(dp) :: dp_(45, op%section%nlocal), dq(45, op%section%nlocal)
        real(dp) :: kp(19, op%section%nlocal, 6), kq(19, op%section%nlocal, 6)
        real(dp) :: tp(19), tq(19), kappa, value
        ! Fine Bloch blocks cancel element contributions; retain wide quadrature sums.
        real(wp), dimension(op%section%nlocal, op%section%nlocal) :: pp, qq, pq, qp
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:)
        type(fortsparse_status_t) :: sparse_status
        integer :: t, q, i, j, b, r, st, fail, nl, ns, nt, per, pos, ci, cj
        nl = op%section%nlocal; ns = op%section%ndof
        nt = op%section%fs%mm%topology%n_triangles
        per = sum((op%copies*nl)**2)
        allocate (rows(per*nt), cols(per*nt), vals(per*nt))
        fail = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(t,q,i,j,b,r,st,basis,orbit,x,weight,v,vd,kmat,dp_,dq,kp,kq, &
        !$omp tp,tq,kappa,value,pp,qq,pq,qp,pos,ci,cj) reduction(max:fail)
        do t = 1, nt
            pp = 0; qq = 0; pq = 0; qp = 0
            do q = 1, op%section%fs%nq
                call vrmhd_fourier_data(op, t, q, x, weight, basis, orbit, st)
                if (st /= 0) then
                    fail = 1; cycle
                end if
                v = matmul(basis, u0(op%section%map(:, t)))
                vd = matmul(basis(:12, :), ud0(op%section%map(:, t)))
                call vrmhd_point_matrix(p, v(:33), vd, source(:, q, t), shift, kmat)
                dp_(:33, :) = basis(:33, :); dp_(34:, :) = basis(:12, :)
                dq(:33, :) = orbit(:33, :); dq(34:, :) = 0
                do r = 1, 6
                    kp(:, :, r) = matmul(kmat(:, :, r), dp_)
                    kq(:, :, r) = matmul(kmat(:, :, r), dq)
                end do
                do i = 1, nl
                    r = op%section%row(i)
                    call vrmhd_section_test(r, basis(:, i), tp)
                    call vrmhd_section_test(r, orbit(:, i), tq)
                    pp(i, :) = pp(i, :) + real(weight, wp)* &
                        real(matmul(tp, kp(:, :, r)), wp)
                    qq(i, :) = qq(i, :) + real(weight, wp)* &
                        real(matmul(tq, kq(:, :, r)), wp)
                    pq(i, :) = pq(i, :) + real(weight, wp)* &
                        real(matmul(tp, kq(:, :, r)), wp)
                    qp(i, :) = qp(i, :) + real(weight, wp)* &
                        real(matmul(tq, kp(:, :, r)), wp)
                end do
            end do
            pos = (t - 1)*per
            do b = 1, size(op%harmonic)
                kappa = vrmhd_fourier_wavenumber(op, b)
                do cj = 1, op%copies(b)
                    do j = 1, nl
                        do ci = 1, op%copies(b)
                            do i = 1, nl
                                pos = pos + 1
                                rows(pos) = op%offset(b) + (ci - 1)*ns &
                                    + op%section%map(i, t)
                                cols(pos) = op%offset(b) + (cj - 1)*ns &
                                    + op%section%map(j, t)
                                if (op%copies(b) == 1) then
                                    value = real(pp(i, j), dp)
                                else if (ci == cj) then
                                    value = real(0.5_dp*(pp(i, j) + &
                                        kappa**2*qq(i, j)), dp)
                                else
                                    value = real(0.5_dp*kappa*(pq(i, j) - qp(i, j)), dp)
                                    if (ci == 2) value = -value
                                end if
                                vals(pos) = value
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        status = fail
        if (status /= 0) return
        call csc_from_triplet(op%ndof, op%ndof, rows, cols, vals, jac, sparse_status)
        if (.not. status_ok(sparse_status)) status = 1
    end subroutine
end module kin6d_vrmhd_fourier
