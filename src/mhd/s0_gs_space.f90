module kin6d_s0_gs_space
    ! A conforming radial P1 x poloidal Fourier subspace on a circular torus.
    ! The fields are sampled exclusively through the solved core-state evaluator.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kin6d_equilibrium_state, only: equilibrium_state_t, EQ_VALUES, &
        EQ_BR, EQ_BZ, EQ_BPHI, EQ_P
    use kin6d_equilibrium_chart, only: cylindrical_chart_t
    use kin6d_s0_gs_space_leaves, only: s0_polar_geometry, s0_polar_jet
    use kin6d_ideal_xi, only: xi_gradient_n, XI_OK
    use kin6d_s0_volume, only: s0_assemble_volume
    use fortfem_generated_line_basis, only: generated_line_jet
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: csc_z_t, fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX
    implicit none
    private
    public :: s0_gs_space_t, s0_gs_tables, s0_gs_assemble
    type :: s0_gs_space_t
        integer :: ndof = 0, n = 0
        real(dp) :: center(2) = 0, radius = 0
        integer, allocatable :: dof(:,:)
        real(dp), allocatable :: weight(:,:), b(:,:,:), p(:,:), rho(:,:)
        complex(dp), allocatable :: value(:,:,:,:), gradient(:,:,:,:,:)
    end type s0_gs_space_t
contains
    subroutine s0_gs_tables(eq, center, radius, radial_cells, angular_points, &
            n, density, space, status)
        type(equilibrium_state_t), intent(in) :: eq
        real(dp), intent(in) :: center(2), radius, density
        integer, intent(in) :: radial_cells, angular_points, n
        type(s0_gs_space_t), intent(out) :: space
        integer, intent(out) :: status
        real(dp), allocatable :: nodes(:), weights(:), points(:,:), fields(:,:)
        integer, allocatable :: statuses(:)
        real(dp) :: h, r, theta, jac, phi(2,4), deriv(2,4), u(3), du(3)
        real(dp) :: basis_gradient(2), basis_hessian(2,2)
        real(dp) :: raw(18), pi, angle_weight
        complex(dp) :: xi(3), dr(3), dz(3)
        integer :: e, q, t, g, a, component, local, i, ierr, failures, nq
        status = 1
        if (.not. allocated(eq%field)) return
        if (.not. allocated(eq%chart)) return
        select type (chart => eq%chart)
        type is (cylindrical_chart_t)
        class default
            return
        end select
        if (.not. all(ieee_is_finite([center,radius,density,eq%mu0]))) return
        if (radius <= 0 .or. density <= 0 .or. eq%mu0 <= 0) return
        if (center(1) <= radius) return
        if (radial_cells < 2 .or. angular_points < 8 .or. n == 0) return
        space%ndof = 1+3*(radial_cells-1)
        space%n = n; space%center = center; space%radius = radius
        h = radius/radial_cells; pi = acos(-1.0_dp)
        angle_weight = 2*pi/angular_points
        allocate(nodes(4),weights(4))
        call gauss_legendre_ab(4,0.0_dp,1.0_dp,nodes,weights)
        do g = 1, size(nodes)
            do a = 1, 2
                call generated_line_jet(a,nodes(g),0.0_dp,phi(a,g), &
                    basis_gradient,basis_hessian)
                deriv(a,g) = basis_gradient(1)/h
            end do
        end do
        nq = size(nodes)*angular_points
        allocate(space%dof(6,radial_cells),space%weight(nq,radial_cells), &
            space%b(3,nq,radial_cells),space%p(nq,radial_cells), &
            space%rho(nq,radial_cells),space%value(3,6,nq,radial_cells), &
            space%gradient(3,3,6,nq,radial_cells), &
            points(2,nq*radial_cells),fields(EQ_VALUES,nq*radial_cells), &
            statuses(nq*radial_cells))
        space%dof = 0; space%value = 0; space%gradient = 0
        space%rho = density
        do e = 1, radial_cells
            do a = 1, 2
                i = e+a-2
                if (i == radial_cells) cycle ! Pointwise fixed conducting wall.
                if (i == 0) then
                    ! m=1 axis: u(0)=v(0), w(0)=0. One combined regular jet.
                    space%dof(1,e) = 1
                else
                    space%dof(3*(a-1)+1:3*a,e) = [3*i-1,3*i,3*i+1]
                end if
            end do
        end do
        failures = 0
        !$omp parallel do default(shared) schedule(static) &
        !$omp private(e,t,g,q,r,theta,jac,a,component,local,u,du,raw,xi,dr,dz) &
        !$omp private(i,ierr) &
        !$omp reduction(+:failures)
        do e = 1, radial_cells
            do t = 1, angular_points
                theta = angle_weight*(t-1)
                do g = 1, size(nodes)
                    q = (t-1)*size(nodes)+g; i = (e-1)*nq+q
                    r = h*(e-1+nodes(g))
                    call s0_polar_geometry(r,theta,center(1),center(2), &
                        points(1,i),points(2,i),jac)
                    space%weight(q,e) = 2*pi*angle_weight*h*weights(g)*jac
                    do a = 1, 2
                        do component = 1, 3
                            local = 3*(a-1)+component
                            if (space%dof(local,e) == 0) cycle
                            u = 0; du = 0
                            u(component) = phi(a,g); du(component) = deriv(a,g)
                            if (e == 1 .and. a == 1) then
                                u(2) = phi(a,g); du(2) = deriv(a,g)
                            end if
                            call s0_polar_jet(r,theta,1.0_dp,u(1),u(2),u(3), &
                                du(1),du(2),du(3),raw(1),raw(2),raw(3), &
                                raw(4),raw(5),raw(6),raw(7),raw(8),raw(9), &
                                raw(10),raw(11),raw(12),raw(13),raw(14),raw(15), &
                                raw(16),raw(17),raw(18))
                            xi = cmplx(raw(1:3),raw(4:6),dp)
                            dr = cmplx(raw(7:9),raw(10:12),dp)
                            dz = cmplx(raw(13:15),raw(16:18),dp)
                            space%value(:,local,q,e) = xi
                            call xi_gradient_n(n,points(1,i),xi,dr,dz, &
                                space%gradient(:,:,local,q,e),ierr)
                            if (ierr /= XI_OK) failures = failures+1
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if (failures /= 0) return
        ! One batched D7/core evaluation; no analytical replacement of B or p.
        call eq%field%sample(points,fields,statuses)
        if (any(statuses /= 0)) then
            status = 2
            return
        end if
        space%b(1,:,:) = reshape(fields(EQ_BR,:),[nq,radial_cells])
        space%b(2,:,:) = reshape(fields(EQ_BPHI,:),[nq,radial_cells])
        space%b(3,:,:) = reshape(fields(EQ_BZ,:),[nq,radial_cells])
        space%p = reshape(fields(EQ_P,:),[nq,radial_cells])
        status = 0
    end subroutine s0_gs_tables

    subroutine s0_gs_assemble(eq, space, stiffness, mass, status)
        type(equilibrium_state_t), intent(in) :: eq
        type(s0_gs_space_t), intent(in) :: space
        type(csc_z_t), intent(out) :: stiffness, mass
        type(fortsparse_status_t), intent(out) :: status
        call status_set(status,FORTSPARSE_INVALID_MATRIX,'unprepared S0 GS space')
        if (.not. allocated(space%dof)) return
        call s0_assemble_volume(space%ndof,space%dof,space%weight,space%b, &
            space%p,space%rho,eq%mu0,space%value,space%gradient, &
            stiffness,mass,status)
    end subroutine s0_gs_assemble
end module kin6d_s0_gs_space
