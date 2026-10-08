module kin6d_grad_shafranov
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use triangle_compat, only: tc_result_t, triangulate_compat
   use basis_p1_2d_module, only: basis_p1_2d_t
   use basis_p2_2d_module, only: basis_p2_2d_t
   use fortfem_triangle_level_quadrature, only: triangle_level_workspace_t, &
       initialize_triangle_level_workspace, integrate_triangle_levels, &
       triangle_level_coefficients
   use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
   use fortnum_krylov, only: real_conjugate_gradient_operator, real_gmres_operator, KRYLOV_OK
   use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights, &
       lagrange_combined_jvp, &
       lagrange_combined_vjp, lagrange_nodes_jvp, lagrange_nodes_vjp
   use kin6d_grad_shafranov_leaves, only: gs_source, gs_source_jvp, gs_source_vjp, &
       gs_field, gs_segment_primitive, &
       gs_current_density, gs_current_sheet, gs_smooth_axis_q
   use kin6d_gs_inverse_leaves, only: gs_normalized_profile_response
   implicit none
   private
   public :: gs_state_t, initialize_gs, solve_gs, gs_linear_response, sample_gs
   public :: write_gs_mesh
   public :: gs_profile_response, gs_profile_implicit_response, gs_profile_tangent_cache
   public :: gs_profile_response_cache_t, prepare_gs_profile_response_cache
   public :: gs_profile_ff_load_basis
   public :: gs_profile_load, gs_profile_load_jvp, gs_profile_load_vjp
   public :: solve_gs_profiles, gs_profile_value, gs_profile_integral
   public :: gs_element_jet, gs_axis_candidate, gs_current_jump
   public :: gs_axis_safety_factor, gs_classical_current
   ! Fixed-boundary P1/P2 weak GS in physical (R,Z), R>0. The domain is a
   ! polygonal approximation; nested flux surfaces are not solver inputs.
   ! A result is a candidate, NOT a continuous certificate. Geometry,
   ! discretization and source-model errors remain separate open budgets.
   type :: gs_state_t
      real(dp), allocatable :: points(:, :), psi(:)
      integer, allocatable :: triangles(:, :)
      integer, allocatable :: neighbors(:, :)
      logical, allocatable :: boundary(:)
      real(dp), allocatable :: stiffness(:, :, :), gradient(:, :, :)
      real(dp), allocatable :: area(:), diagonal(:)
      integer :: element_order = 1
      real(dp) :: residual = huge(1.0_dp)
      integer :: iterations = 0
      character(32) :: claim = 'candidate'
   end type gs_state_t
   ! One owned immutable response context for one state/profile, with no globals.
   type :: gs_profile_response_cache_t
      private
      type(gs_state_t) :: state
      real(dp), allocatable :: flux(:), ff(:), matrix(:, :, :), axis_response(:)
      real(dp), allocatable :: axis_shape(:)
      real(dp), allocatable :: load_basis(:, :), load_error(:, :)
      integer :: axis_element = 0
      logical :: ready = .false.
   end type gs_profile_response_cache_t
   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*), info
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      end subroutine dgesv
   end interface
contains
   subroutine initialize_gs(state, contour, max_area, status, element_order)
      type(gs_state_t), intent(out) :: state
      real(dp), intent(in) :: contour(:, :), max_area
      integer, intent(out) :: status
      integer, intent(in), optional :: element_order
      type(tc_result_t) :: mesh
      type(basis_p1_2d_t) :: basis
      type(gauss_quadrature_triangle_t) :: quadrature
      real(dp) :: jac(2, 2), matrix(2, 2), det, grads(2, 3)
      real(dp) :: point(2), shape(6), grad_q(2, 6), weight
      integer, allocatable :: segments(:, :)
      integer :: j, e, k, l, q, pivot(2), stat, n, ndof
      status = 1
      if (present(element_order)) state%element_order = element_order
      if (state%element_order < 1 .or. state%element_order > 2) return
      if (size(contour, 1) /= 2) return
      n = size(contour, 2)
      if (n < 3 .or. max_area <= 0.0_dp) return
      if (.not. ieee_is_finite(max_area)) return
      if (any(.not. ieee_is_finite(contour))) return
      if (minval(contour(1, :)) <= 0.0_dp) return
      allocate (segments(2, n))
      do j = 1, n
         segments(:, j) = [j, mod(j, n) + 1]
      end do
      call triangulate_compat(contour, segments, reshape([real(dp) ::], [2, 0]), &
                           mesh, stat, min_angle=20.0_dp, max_area=max_area, nobisect=0)
      if (stat /= 0 .or. mesh%ntriangles < 1) return
      state%neighbors = mesh%neighbors
      if (state%element_order == 1) then
         state%points = mesh%points
         state%triangles = mesh%triangles
         state%boundary = mesh%pointmarkers /= 0
      else
         call quadratic_mesh(mesh, state, stat)
         if (stat /= 0) return
      end if
      ndof = size(state%triangles, 1)
      n = size(state%points, 2)
      allocate (state%psi(n), state%diagonal(n), &
                state%stiffness(ndof, ndof, mesh%ntriangles), &
                state%gradient(2, 3, mesh%ntriangles), state%area(mesh%ntriangles))
      state%psi = 0.0_dp; state%diagonal = 0.0_dp
      call gs_quadrature(state%element_order, quadrature)
      do e = 1, mesh%ntriangles
         call basis%compute_jacobian(state%points(:, state%triangles(1:3, e)), &
                                     jac, det)
         if (abs(det) <= tiny(det)) return
         matrix = transpose(jac)
         do j = 1, 3
            grads(:, j) = basis%grad(j, 0.0_dp, 0.0_dp)
         end do
         call dgesv(2, 3, matrix, 2, pivot, grads, 2, stat)
         if (stat /= 0) return
         ! Affine geometry gradients, not the represented P2 flux gradients.
         state%gradient(:, :, e) = grads
         state%area(e) = 0.5_dp*abs(det)
         state%stiffness(:, :, e) = 0.0_dp
         do q = 1, quadrature%n_points
            call element_shapes(state, e, quadrature%xi(q), &
                                quadrature%eta(q), shape(:ndof), grad_q(:, :ndof))
            point = matmul(state%points(:, state%triangles(1:3, e)), &
                 [1-quadrature%xi(q)-quadrature%eta(q), &
                  quadrature%xi(q), quadrature%eta(q)])
            weight = 2*state%area(e)*quadrature%weights(q)/point(1)
            do l = 1, ndof
               do k = 1, ndof
                  state%stiffness(k, l, e) = state%stiffness(k, l, e) + &
                       weight*dot_product(grad_q(:, k), grad_q(:, l))
               end do
            end do
         end do
         do l = 1, ndof
            j = state%triangles(l, e)
            state%diagonal(j) = state%diagonal(j) + state%stiffness(l, l, e)
         end do
      end do
      where (state%boundary) state%diagonal = 1.0_dp
      if (minval(state%diagonal) <= 0) return
      status = 0
   end subroutine initialize_gs

   subroutine gs_quadrature(order, quadrature)
      integer, intent(in) :: order
      type(gauss_quadrature_triangle_t), intent(out) :: quadrature
      if (order == 1) then
         call quadrature%init(2)
      else
         ! Positive seven-point rule; degree-five exact for polynomials.
         ! The physical 1/R weight is rational: this is not exact integration.
         call quadrature%init(5)
      end if
   end subroutine gs_quadrature

   subroutine quadratic_mesh(mesh, state, status)
      type(tc_result_t), intent(in) :: mesh
      type(gs_state_t), intent(inout) :: state
      integer, intent(out) :: status
      real(dp), allocatable :: points(:, :)
      logical, allocatable :: boundary(:)
      integer :: e, k, other, other_edge, first, second, count, j, next
      status = 1
      allocate (state%triangles(6, mesh%ntriangles))
      state%triangles = 0
      state%triangles(1:3, :) = mesh%triangles
      allocate (points(2, mesh%npoints+3*mesh%ntriangles), &
                boundary(mesh%npoints+3*mesh%ntriangles))
      points(:, :mesh%npoints) = mesh%points
      boundary(:mesh%npoints) = mesh%pointmarkers /= 0
      count = mesh%npoints
      ! Triangle's neighbour k is across local edge (k,mod(k,3)+1).
      ! Each edge creates one shared midpoint, in O(number of triangles).
      do e = 1, mesh%ntriangles
         do k = 1, 3
            if (state%triangles(k+3, e) /= 0) cycle
            first = mesh%triangles(k, e)
            second = mesh%triangles(mod(k, 3)+1, e)
            other = mesh%neighbors(k, e)
            other_edge = 0
            if (other > 0) then
               if (other > mesh%ntriangles .or. other == e) return
               do j = 1, 3
                  next = mod(j, 3)+1
                  if (mesh%triangles(j, other) /= second) cycle
                  if (mesh%triangles(next, other) /= first) cycle
                  if (mesh%neighbors(j, other) /= e) return
                  other_edge = j
                  exit
               end do
               if (other_edge == 0) return
            end if
            count = count+1
            points(:, count) = 0.5_dp*(mesh%points(:, first)+mesh%points(:, second))
            boundary(count) = other < 1
            state%triangles(k+3, e) = count
            if (other_edge > 0) state%triangles(other_edge+3, other) = count
         end do
      end do
      state%points = points(:, :count)
      state%boundary = boundary(:count)
      status = 0
   end subroutine quadratic_mesh

   subroutine element_shapes(state, e, xi, eta, shape, gradients, hessians)
      type(gs_state_t), intent(in) :: state
      integer, intent(in) :: e
      real(dp), intent(in) :: xi, eta
      real(dp), intent(out) :: shape(:), gradients(:, :)
      real(dp), intent(out), optional :: hessians(:, :, :)
      type(basis_p2_2d_t) :: basis
      type(basis_p1_2d_t) :: affine_basis
      integer :: j
      if (state%element_order == 1) then
         do j = 1, 3
            shape(j) = affine_basis%eval(j, xi, eta)
         end do
         gradients = state%gradient(:, :, e)
         if (present(hessians)) hessians = 0.0_dp
      else
         do j = 1, 6
            shape(j) = basis%eval(j, xi, eta)
            gradients(:, j) = matmul(state%gradient(:, 2:3, e), &
                                    basis%grad(j, xi, eta))
            if (present(hessians)) then
               hessians(:, :, j) = matmul(matmul(state%gradient(:, 2:3, e), &
                    basis%hessian(j, xi, eta)), transpose(state%gradient(:, 2:3, e)))
            end if
         end do
      end if
   end subroutine element_shapes

   subroutine apply_gs(state, input, output)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: input(:)
      real(dp), intent(out) :: output(:)
      integer :: e, k, l, i, j
      output = 0.0_dp
      do e = 1, size(state%triangles, 2)
         do l = 1, size(state%triangles, 1)
            j = state%triangles(l, e)
            if (state%boundary(j)) cycle
            do k = 1, size(state%triangles, 1)
               i = state%triangles(k, e)
               if (state%boundary(i)) cycle
               output(i) = output(i) + state%stiffness(k, l, e)*input(j)
            end do
         end do
      end do
      where (state%boundary) output = input
   end subroutine apply_gs

   subroutine gs_linear_response(state, source, response, status)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: source(:)
      real(dp), allocatable, intent(out) :: response(:)
      integer, intent(out) :: status
      real(dp) :: residual, tolerance
      integer :: iterations
      status = 1
      if (size(source) /= size(state%psi)) return
      allocate (response(size(source)))
      response = 0.0_dp
      tolerance = 1.0e-11_dp
      if (norm2(source) > 0.0_dp) tolerance = tolerance*min(1.0_dp, norm2(source))
      call real_conjugate_gradient_operator(matvec, source, response, &
                   tolerance, max(500, 10*size(source)), status, iterations, residual, &
                                            preconditioner)
   contains
      subroutine matvec(input, output)
         real(dp), intent(in) :: input(:)
         real(dp), intent(out) :: output(:)
         call apply_gs(state, input, output)
      end subroutine matvec
      subroutine preconditioner(input, output)
         real(dp), intent(in) :: input(:)
         real(dp), intent(out) :: output(:)
         output = input/state%diagonal
      end subroutine preconditioner
   end subroutine gs_linear_response

   subroutine solve_gs(state, mu0, pprime, ffprime, status)
      type(gs_state_t), intent(inout) :: state
      real(dp), intent(in) :: mu0, pprime, ffprime
      integer, intent(out) :: status
      real(dp), allocatable :: rhs(:), solution(:), action(:)
      status = 1
      if (mu0 <= 0.0_dp) return
      if (any(.not. ieee_is_finite([mu0, pprime, ffprime]))) return
      allocate (rhs(size(state%psi)), action(size(state%psi)))
      call assemble_gs_rhs(state, mu0, pprime, ffprime, rhs)
      call gs_linear_response(state, rhs, solution, status)
      if (status /= KRYLOV_OK) return
      state%psi = solution
      call apply_gs(state, solution, action)
      state%residual = norm2(action - rhs)/max(norm2(rhs), tiny(1.0_dp))
   end subroutine solve_gs

   subroutine assemble_gs_rhs(state, mu0, pprime, ffprime, rhs)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: mu0, pprime, ffprime
      real(dp), intent(out) :: rhs(:)
      type(gauss_quadrature_triangle_t) :: quadrature
      real(dp) :: point(2), shape(6), gradients(2, 6), source, weight
      integer :: e, q, k, l, i, j, ndof
      rhs = 0.0_dp
      ndof = size(state%triangles, 1)
      call gs_quadrature(state%element_order, quadrature)
      do e = 1, size(state%triangles, 2)
         do q = 1, quadrature%n_points
            call element_shapes(state, e, quadrature%xi(q), &
                 quadrature%eta(q), shape(:ndof), gradients(:, :ndof))
            point = matmul(state%points(:, state%triangles(1:3, e)), &
                 [1-quadrature%xi(q)-quadrature%eta(q), &
                  quadrature%xi(q), quadrature%eta(q)])
            call gs_source(point(1), mu0, pprime, ffprime, source)
            weight = 2*state%area(e)*quadrature%weights(q)
            do k = 1, ndof
               i = state%triangles(k, e)
               rhs(i) = rhs(i) + weight*source*shape(k)
            end do
         end do
         do l = 1, ndof
            j = state%triangles(l, e)
            if (.not. state%boundary(j)) cycle
            do k = 1, ndof
               i = state%triangles(k, e)
               rhs(i) = rhs(i) - state%stiffness(k, l, e)*state%psi(j)
            end do
         end do
      end do
      where (state%boundary) rhs = state%psi
   end subroutine assemble_gs_rhs

    subroutine gs_profile_load(state, mu0, flux, pp, ff, rhs, status, &
                               relative_tolerance, max_evaluations, error_estimate)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:)
        real(dp), intent(out) :: rhs(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(in), optional :: max_evaluations
        real(dp), intent(out), optional :: error_estimate
        call profile_load_products(state, mu0, flux, pp, ff, 1, rhs, status, &
                                   relative_tolerance, max_evaluations, error_estimate)
    end subroutine gs_profile_load

    subroutine gs_profile_load_jvp(state, mu0, flux, pp, ff, dpsi, dmu0, &
                                   dflux, dpp, dff, tangent, status, &
                                   relative_tolerance, max_evaluations)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:), dpsi(:), dmu0
        real(dp), intent(in) :: dflux(:), dpp(:), dff(:)
        real(dp), intent(out) :: tangent(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(in), optional :: max_evaluations
        real(dp) :: rhs(size(state%psi))
        tangent = 0.0_dp
        status = 1
        if (size(tangent) /= size(state%psi)) return
        if (size(dpsi) /= size(state%psi)) return
        if (size(dflux) /= size(flux)) return
        if (size(dpp) /= size(pp) .or. size(dff) /= size(ff)) return
        if (any(.not. ieee_is_finite(dpsi))) return
        if (.not. ieee_is_finite(dmu0)) return
        if (any(.not. ieee_is_finite(dflux))) return
        if (any(.not. ieee_is_finite(dpp))) return
        if (any(.not. ieee_is_finite(dff))) return
        call profile_load_products(state, mu0, flux, pp, ff, 2, rhs, status, &
                            relative_tolerance, max_evaluations, dpsi=dpsi, dmu0=dmu0, &
                                   dflux=dflux, dpp=dpp, dff=dff, tangent=tangent)
    end subroutine gs_profile_load_jvp

    subroutine gs_profile_load_vjp(state, mu0, flux, pp, ff, cotangent, &
                                   gpsi, gmu0, gflux, gpp, gff, status, &
                                   relative_tolerance, max_evaluations)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:), cotangent(:)
        real(dp), intent(out) :: gpsi(:), gmu0, gflux(:), gpp(:), gff(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(in), optional :: max_evaluations
        real(dp) :: rhs(size(state%psi))
        gpsi = 0; gmu0 = 0; gflux = 0; gpp = 0; gff = 0
        status = 1
        if (size(cotangent) /= size(state%psi)) return
        if (size(gpsi) /= size(state%psi)) return
        if (size(gflux) /= size(flux)) return
        if (size(gpp) /= size(pp) .or. size(gff) /= size(ff)) return
        if (any(.not. ieee_is_finite(cotangent))) return
        call profile_load_products(state, mu0, flux, pp, ff, 3, rhs, status, &
                             relative_tolerance, max_evaluations, cotangent=cotangent, &
                                   gpsi=gpsi, gmu0=gmu0, gflux=gflux, gpp=gpp, gff=gff)
    end subroutine gs_profile_load_vjp

    subroutine profile_load_products(state, mu0, flux, pp, ff, mode, rhs, status, &
                                     relative_tolerance, max_evaluations, &
                                     error_estimate, &
                                     dpsi, dmu0, dflux, dpp, dff, tangent, &
                                     cotangent, gpsi, gmu0, gflux, gpp, gff)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:)
        integer, intent(in) :: mode
        real(dp), intent(out) :: rhs(:)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(in), optional :: max_evaluations
        real(dp), intent(out), optional :: error_estimate
        real(dp), intent(in), optional :: dpsi(:), dmu0, dflux(:), dpp(:), dff(:)
        real(dp), intent(out), optional :: tangent(:)
        real(dp), intent(in), optional :: cotangent(:)
        real(dp), intent(out), optional :: gpsi(:), gmu0, gflux(:), gpp(:), gff(:)
        type(triangle_level_workspace_t) :: workspace
        real(dp), allocatable :: assembled_error(:), tangent_error(:)
        real(dp), allocatable :: values(:), errors(:), epsabs(:)
        real(dp), allocatable :: nodes_pp(:), nodes_ff(:), values_pp(:), values_ff(:)
        real(dp), allocatable :: shifted_flux(:), shifted_dflux(:)
        real(dp) :: coefficients(6), tolerance, local_psi(6), local_dpsi(6)
        real(dp) :: vertices(2, 3)
        real(dp) :: source_scale, error_norm, rhs_norm, pp_slope, ff_slope, rmin, rmax
        real(dp) :: activity_scale, dual_scale
        integer :: e, k, l, node, boundary_node, ndof, work_limit, nc, nf
        integer :: ip, iff, ix, im
        character(160) :: message
        status = 1
        rhs = 0.0_dp
        if (present(error_estimate)) error_estimate = huge(1.0_dp)
        if (mu0 <= 0.0_dp) return
        if (.not. ieee_is_finite(mu0)) return
        if (size(rhs) /= size(state%psi)) return
        if (size(flux) < 2) return
        if (size(pp) /= size(flux) .or. size(ff) /= size(flux)) return
        if (any(.not. ieee_is_finite(flux))) return
        if (any(.not. ieee_is_finite(pp))) return
        if (any(.not. ieee_is_finite(ff))) return
        if (any(.not. ieee_is_finite(state%psi))) return
        if (any(flux(2:) <= flux(:size(flux) - 1))) return
        tolerance = 1.0e-11_dp
        if (present(relative_tolerance)) tolerance = relative_tolerance
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) return
        work_limit = 200000
        if (present(max_evaluations)) work_limit = max_evaluations
        if (work_limit < 1) return
        ndof = size(state%triangles, 1)
        nf = size(flux)
        nc = ndof
        if (mode == 2) nc = 2*ndof
        if (mode == 3) nc = 2*ndof + 3*nf + 1
        ip = 2*ndof
        iff = ip + nf
        ix = iff + nf
        im = ix + nf + 1
        call initialize_triangle_level_workspace(workspace, nc, nf, &
                                                 status, max_evaluations=work_limit)
        if (status /= 0) return
        allocate (assembled_error(size(rhs)), values(nc), errors(nc), epsabs(nc), &
                  nodes_pp(nf), nodes_ff(nf), values_pp(nf), values_ff(nf), &
                  shifted_flux(nf))
        assembled_error = 0.0_dp
        if (mode == 2) then
            allocate(tangent_error(size(rhs)), shifted_dflux(nf))
            tangent_error = 0.0_dp
        end if
        pp_slope = maxval(abs((pp(2:) - pp(:nf - 1))/(flux(2:) - flux(:nf - 1))))
        ff_slope = maxval(abs((ff(2:) - ff(:nf - 1))/(flux(2:) - flux(:nf - 1))))
        do e = 1, size(state%triangles, 2)
            local_psi(:ndof) = state%psi(state%triangles(:, e))
            shifted_flux = flux - local_psi(1)
            if (mode == 2) then
                ! Differentiate the primal's common local gauge before composition.
                local_dpsi(:ndof) = dpsi(state%triangles(:, e)) - &
                    dpsi(state%triangles(1, e))
                shifted_dflux = dflux - dpsi(state%triangles(1, e))
            end if
            vertices = state%points(:, state%triangles(1:3, e))
            call triangle_level_coefficients(state%element_order, &
                                             local_psi(:ndof), coefficients, status)
            if (status /= 0) exit
            ! A full cell on a genuine profile kink has no linear derivative.
            if (mode /= 1) then
                if (all(coefficients(2:) == 0.0_dp)) then
                    if (profile_has_kink(flux, pp, ff, coefficients(1))) then
                        status = 4
                        exit
                    end if
                end if
            end if
            rmin = minval(vertices(1, :)); rmax = maxval(vertices(1, :))
            if (rmin <= 0.0_dp) then
                status = 1
                exit
            end if
            source_scale = mu0*rmax*maxval(abs(pp)) + maxval(abs(ff))/rmin
            epsabs = max(tiny(1.0_dp), 0.01_dp*tolerance*source_scale)
            if (mode == 2) then
                activity_scale = maxval(abs(local_dpsi(:ndof))) + &
                                 maxval(abs(shifted_dflux))
                source_scale = abs(dmu0)*rmax*maxval(abs(pp)) + &
                               mu0*rmax*(maxval(abs(dpp)) + pp_slope*activity_scale) + &
                               (maxval(abs(dff)) + ff_slope*activity_scale)/rmin
                epsabs(ndof + 1:) = max(tiny(1.0_dp), 0.01_dp*tolerance*source_scale)
            else if (mode == 3) then
                dual_scale = maxval(abs(cotangent(state%triangles(:, e))))
                epsabs(ndof + 1:2*ndof) = max(tiny(1.0_dp), &
                       0.01_dp*tolerance*dual_scale*(mu0*rmax*pp_slope + ff_slope/rmin))
           epsabs(ip + 1:iff) = max(tiny(1.0_dp), 0.01_dp*tolerance*dual_scale*mu0*rmax)
               epsabs(iff + 1:ix) = max(tiny(1.0_dp), 0.01_dp*tolerance*dual_scale/rmin)
                epsabs(ix + 1:im - 1) = epsabs(ndof + 1)
                epsabs(im) = max(tiny(1.0_dp), &
                                 0.01_dp*tolerance*dual_scale*rmax*maxval(abs(pp)))
            end if
            if (any(.not. ieee_is_finite(epsabs))) then
                status = 1
                exit
            end if
            call integrate_triangle_levels(coefficients, flux, source_integrand, &
                  epsabs, 0.01_dp*tolerance, workspace, values, errors, status, message)
            if (status /= 0) exit
            values = 2*state%area(e)*values
            errors = 2*state%area(e)*errors
            do k = 1, ndof
                node = state%triangles(k, e)
                rhs(node) = rhs(node) + values(k)
                assembled_error(node) = assembled_error(node) + errors(k)
                if (mode == 2) then
                    tangent(node) = tangent(node) + values(ndof + k)
                    tangent_error(node) = tangent_error(node) + errors(ndof + k)
                end if
                if (mode == 3) gpsi(node) = gpsi(node) + values(ndof + k)
            end do
            if (mode == 3) then
                gpp = gpp + values(ip + 1:iff)
                gff = gff + values(iff + 1:ix)
                gflux = gflux + values(ix + 1:im - 1)
                gmu0 = gmu0 + values(im)
            end if
            do l = 1, ndof
                boundary_node = state%triangles(l, e)
                if (.not. state%boundary(boundary_node)) cycle
                do k = 1, ndof
                    node = state%triangles(k, e)
                    rhs(node) = rhs(node) - &
                                state%stiffness(k, l, e)*state%psi(boundary_node)
                    if (mode == 2) tangent(node) = tangent(node) - &
                                            state%stiffness(k, l, e)*dpsi(boundary_node)
                    if (mode == 3) then
                        if (state%boundary(node)) cycle
                        gpsi(boundary_node) = gpsi(boundary_node) - &
                                              state%stiffness(k, l, e)*cotangent(node)
                    end if
                end do
            end do
        end do
        if (status /= 0) then
            call clear_outputs()
            return
        end if
        where (state%boundary)
            rhs = state%psi
            assembled_error = 0.0_dp
        end where
        if (mode == 2) then
            where (state%boundary)
                tangent = dpsi
                tangent_error = 0.0_dp
            end where
        end if
        if (mode == 3) then
            where (state%boundary) gpsi = gpsi + cotangent
        end if
        error_norm = norm2(assembled_error)
        rhs_norm = norm2(rhs)
        if (present(error_estimate)) error_estimate = error_norm
        if (error_norm > tolerance*max(rhs_norm, tiny(1.0_dp))) then
            status = 3 ! Global accumulated quadrature estimator exceeds the budget.
            call clear_outputs()
            return
        end if
        if (mode == 2) then
            if (norm2(tangent_error) > tolerance*max(norm2(tangent), tiny(1.0_dp))) then
                status = 3 ! The directional load retains its own relative budget.
                call clear_outputs()
                return
            end if
        end if
        status = 0
    contains
        subroutine source_integrand(point, output, callback_status)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: callback_status
            real(dp) :: shape(6), r, psi_q, pp_q, ff_q, source, dpsi_q, dpp_q, dff_q
            real(dp) :: dsource, dual, gr, gmu, gp, gf, xp, xf
            real(dp) :: local_dual(6)
            call element_values(state%element_order, point(1), point(2), shape(:ndof))
          r = dot_product(vertices(1, :), [1 - point(1) - point(2), point(1), point(2)])
            ! Shift flux and knots together before interpolating the profile law.
            psi_q = dot_product(shape(:ndof), local_psi(:ndof) - local_psi(1))
            pp_q = gs_profile_value(shifted_flux, pp, psi_q)
            ff_q = gs_profile_value(shifted_flux, ff, psi_q)
            call gs_source(r, mu0, pp_q, ff_q, source)
            output(:ndof) = source*shape(:ndof)
            if (mode == 2) then
                dpsi_q = dot_product(shape(:ndof), local_dpsi(:ndof))
                call profile_tangent(shifted_flux, pp, psi_q, shifted_dflux, &
                    dpp, dpsi_q, dpp_q)
                call profile_tangent(shifted_flux, ff, psi_q, shifted_dflux, &
                    dff, dpsi_q, dff_q)
             call gs_source_jvp(r, mu0, pp_q, ff_q, 0.0_dp, dmu0, dpp_q, dff_q, dsource)
                output(ndof + 1:) = dsource*shape(:ndof)
            else if (mode == 3) then
                local_dual(:ndof) = cotangent(state%triangles(:, e))
                where (state%boundary(state%triangles(:, e))) local_dual(:ndof) = 0
                dual = dot_product(shape(:ndof), local_dual(:ndof))
                call gs_source_vjp(r, mu0, pp_q, ff_q, dual, gr, gmu, gp, gf)
              call profile_adjoint(shifted_flux, pp, psi_q, gp, xp, nodes_pp, values_pp)
              call profile_adjoint(shifted_flux, ff, psi_q, gf, xf, nodes_ff, values_ff)
                output(ndof + 1:2*ndof) = (xp + xf)*shape(:ndof)
                output(ip + 1:iff) = values_pp
                output(iff + 1:ix) = values_ff
                output(ix + 1:im - 1) = nodes_pp + nodes_ff
                output(im) = gmu
            end if
            callback_status = 0
            if (any(.not. ieee_is_finite(output))) callback_status = 1
        end subroutine source_integrand

        subroutine clear_outputs()
            rhs = 0
            if (present(tangent)) tangent = 0
            if (present(gpsi)) gpsi = 0
            if (present(gmu0)) gmu0 = 0
            if (present(gflux)) gflux = 0
            if (present(gpp)) gpp = 0
            if (present(gff)) gff = 0
        end subroutine clear_outputs
    end subroutine profile_load_products

    logical function profile_has_kink(flux, pp, ff, psi) result(kink)
        real(dp), intent(in) :: flux(:), pp(:), ff(:), psi
        real(dp) :: weights(2), left_pp, right_pp, left_ff, right_ff
        integer :: j
        kink = .false.
        do j = 1, size(flux)
            if (psi /= flux(j)) cycle
            left_pp = 0; right_pp = 0; left_ff = 0; right_ff = 0
            if (j > 1) then
                call lagrange_deriv_weights(2,psi,flux(j-1:j),weights)
                left_pp = dot_product(weights,pp(j-1:j))
                left_ff = dot_product(weights,ff(j-1:j))
            end if
            if (j < size(flux)) then
                call lagrange_deriv_weights(2,psi,flux(j:j+1),weights)
                right_pp = dot_product(weights,pp(j:j+1))
                right_ff = dot_product(weights,ff(j:j+1))
            end if
            kink = left_pp /= right_pp .or. left_ff /= right_ff
            return
        end do
    end function profile_has_kink

    subroutine profile_tangent(flux, values, psi, dflux, dvalues, dpsi, tangent)
        real(dp), intent(in) :: flux(:), values(:), psi, dflux(:), dvalues(:), dpsi
        real(dp), intent(out) :: tangent
        real(dp) :: motion
        integer :: j
        if (psi <= flux(1)) then
            tangent = dvalues(1)
            return
        end if
        if (psi >= flux(size(flux))) then
            tangent = dvalues(size(flux))
            return
        end if
        do j = 1, size(flux) - 1
            if (psi <= flux(j + 1)) exit
        end do
        call lagrange_combined_jvp(2, psi, flux(j:j + 1), values(j:j + 1), &
                                   dpsi, dvalues(j:j + 1), tangent)
 call lagrange_nodes_jvp(2, psi, flux(j:j + 1), values(j:j + 1), dflux(j:j + 1), motion)
        tangent = tangent + motion
    end subroutine profile_tangent

    subroutine profile_adjoint(flux, values, psi, dual, psi_bar, flux_bar, values_bar)
        real(dp), intent(in) :: flux(:), values(:), psi, dual
        real(dp), intent(out) :: psi_bar, flux_bar(:), values_bar(:)
        integer :: j
        psi_bar = 0; flux_bar = 0; values_bar = 0
        if (psi <= flux(1)) then
            values_bar(1) = dual
            return
        end if
        if (psi >= flux(size(flux))) then
            values_bar(size(flux)) = dual
            return
        end if
        do j = 1, size(flux) - 1
            if (psi <= flux(j + 1)) exit
        end do
        call lagrange_combined_vjp(2, psi, flux(j:j + 1), values(j:j + 1), &
                                   dual, psi_bar, values_bar(j:j + 1))
call lagrange_nodes_vjp(2, psi, flux(j:j + 1), values(j:j + 1), dual, flux_bar(j:j + 1))
    end subroutine profile_adjoint

    subroutine element_values(order, xi, eta, shape)
        integer, intent(in) :: order
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: shape(:)
        type(basis_p1_2d_t) :: affine
        type(basis_p2_2d_t) :: quadratic
        integer :: j
        do j = 1, size(shape)
            if (order == 1) then
                shape(j) = affine%eval(j, xi, eta)
            else
                shape(j) = quadratic%eval(j, xi, eta)
            end if
        end do
    end subroutine element_values

    subroutine gs_profile_tangent_cache(state, flux, ff, axis_flux, mass, &
            axis_response, status, max_evaluations)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:), ff(:), axis_flux
        real(dp), allocatable, intent(out) :: mass(:, :, :), axis_response(:)
        integer, intent(out) :: status
        integer, intent(in), optional :: max_evaluations
        type(triangle_level_workspace_t) :: workspace
        real(dp), allocatable :: values(:), errors(:), epsabs(:), shifted_flux(:)
        real(dp) :: coefficients(6), local_psi(6), vertices(2, 3), scale, axis_scale
        real(dp) :: zero_profile(size(flux)), local_axis(6)
        integer :: ndof, nc, e, k, l, index, work_limit
        character(160) :: message
        status = 1
        ndof = size(state%triangles, 1)
        allocate(mass(ndof, ndof, size(state%triangles, 2)), &
            axis_response(size(state%psi)))
        mass = 0; axis_response = 0
        if (size(flux) < 2 .or. size(ff) /= size(flux)) return
        if (any(.not. ieee_is_finite(flux))) return
        if (any(.not. ieee_is_finite(ff))) return
        if (any(.not. ieee_is_finite(state%psi))) return
        if (.not. ieee_is_finite(axis_flux)) return
        if (axis_flux == 0.0_dp) return
        if (any(flux(2:) <= flux(:size(flux)-1))) return
        work_limit = 200000
        if (present(max_evaluations)) work_limit = max_evaluations
        if (work_limit < 1) return
        nc = ndof*ndof + ndof
        call initialize_triangle_level_workspace(workspace, nc, size(flux), &
            status, max_evaluations=work_limit)
        if (status /= 0) return
        allocate(values(nc), errors(nc), epsabs(nc), shifted_flux(size(flux)))
        zero_profile = 0
        do e = 1, size(state%triangles, 2)
            local_psi(:ndof) = state%psi(state%triangles(:, e))
            shifted_flux = flux - local_psi(1)
            vertices = state%points(:, state%triangles(1:3, e))
            if (minval(vertices(1, :)) <= 0) then
                status = 1
                exit
            end if
            call triangle_level_coefficients(state%element_order, &
                local_psi(:ndof), coefficients, status)
            if (status /= 0) exit
            if (all(coefficients(2:) == 0)) then
                if (profile_has_kink(flux, zero_profile, ff, coefficients(1))) then
                    status = 4
                    exit
                end if
            end if
            scale = maxval(abs((ff(2:)-ff(:size(ff)-1))/ &
                (flux(2:)-flux(:size(flux)-1))))/minval(vertices(1, :))
            axis_scale = scale*maxval(abs(local_psi(:ndof)))/abs(axis_flux)
            epsabs(:ndof*ndof) = max(tiny(1.0_dp), 1.0e-13_dp*scale)
            epsabs(ndof*ndof+1:) = max(tiny(1.0_dp), 1.0e-13_dp*axis_scale)
            if (any(.not. ieee_is_finite(epsabs))) then
                status = 1
                exit
            end if
            call integrate_triangle_levels(coefficients, flux, cache_integrand, &
                epsabs, 1.0e-13_dp, workspace, values, errors, status, message)
            if (status /= 0) exit
            index = 0
            do l = 1, ndof
                do k = 1, ndof
                    index = index + 1
                    mass(k,l,e) = 2*state%area(e)*values(index)
                end do
            end do
            local_axis(:ndof) = 2*state%area(e)*values(ndof*ndof+1:)
            do k = 1, ndof
                if (state%boundary(state%triangles(k,e))) cycle
                axis_response(state%triangles(k,e)) = &
                    axis_response(state%triangles(k,e)) + local_axis(k)
            end do
        end do
        if (status /= 0) then
            mass = 0; axis_response = 0
        end if
    contains
        subroutine cache_integrand(point, output, callback_status)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: callback_status
            real(dp) :: shape(6), psi_q, r, slope, coefficient, source
            integer :: i, j, offset
            call element_values(state%element_order,point(1),point(2),shape(:ndof))
            psi_q = dot_product(shape(:ndof),local_psi(:ndof)-local_psi(1))
            r = dot_product(vertices(1,:),[1-point(1)-point(2),point(1),point(2)])
            call profile_tangent(shifted_flux,ff,psi_q,zero_profile, &
                zero_profile,1.0_dp,slope)
            call gs_source(r,1.0_dp,0.0_dp,slope,source)
            offset = 0
            do j = 1, ndof
                do i = 1, ndof
                    offset = offset + 1
                    output(offset) = source*shape(i)*shape(j)
                end do
            end do
            call gs_normalized_profile_response(slope,0.0_dp, &
                (psi_q+local_psi(1))/axis_flux,1.0_dp,coefficient)
            call gs_source(r,1.0_dp,0.0_dp,coefficient,source)
            output(ndof*ndof+1:) = -source*shape(:ndof)
            callback_status = 0
            if (any(.not. ieee_is_finite(output))) callback_status = 1
        end subroutine cache_integrand
    end subroutine gs_profile_tangent_cache

    subroutine gs_profile_response(state, flux, ff, solution, residual, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:), ff(:)
        real(dp), allocatable, intent(out) :: solution(:)
        real(dp), intent(out) :: residual
        integer, intent(out) :: status
        real(dp), allocatable :: rhs(:), action(:), pp(:)
        status = 1; residual = huge(1.0_dp)
        if (size(flux) < 2 .or. size(ff) /= size(flux)) return
        if (.not. all(ieee_is_finite(flux))) return
        if (.not. all(ieee_is_finite(ff))) return
        if (any(flux(2:) <= flux(:size(flux)-1))) return
        allocate(rhs(size(state%psi)),action(size(state%psi)),pp(size(flux)))
        pp=0
        call gs_profile_load(state,1.0_dp,flux,pp,ff,rhs,status)
        if(status/=0)return
        call apply_gs(state,state%psi,action)
        residual=norm2(action-rhs)/max(norm2(rhs),tiny(1.0_dp))
        call gs_linear_response(state,rhs,solution,status)
    end subroutine gs_profile_response

    subroutine gs_profile_ff_load_basis(state, flux, loads, estimates, status, &
            relative_tolerance, max_evaluations)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:)
        real(dp), allocatable, intent(out) :: loads(:, :), estimates(:, :)
        integer, intent(out) :: status
        real(dp), intent(in), optional :: relative_tolerance
        integer, intent(in), optional :: max_evaluations
        type(triangle_level_workspace_t) :: workspace
        real(dp), allocatable :: values(:), errors(:), epsabs(:), shifted_flux(:)
        real(dp) :: coefficients(6), local_psi(6), vertices(2, 3), tolerance, rmin
        integer :: ndof, nf, e, k, j, node, index, work_limit
        character(160) :: message
        status = 1
        if (.not. allocated(state%psi)) return
        nf = size(flux)
        allocate (loads(size(state%psi), nf), estimates(size(state%psi), nf))
        loads = 0; estimates = huge(1.0_dp)
        if (nf < 2) return
        if (any(.not. ieee_is_finite(flux))) return
        if (any(flux(2:) <= flux(:nf - 1))) return
        if (.not. allocated(state%triangles)) return
        if (.not. allocated(state%points)) return
        if (.not. allocated(state%boundary)) return
        if (.not. allocated(state%area)) return
        ndof = size(state%triangles, 1)
        if (state%element_order == 1) then
            if (ndof /= 3) return
        else if (state%element_order == 2) then
            if (ndof /= 6) return
        else
            return
        end if
        if (size(state%triangles, 2) < 1) return
        if (any(shape(state%points) /= [2, size(state%psi)])) return
        if (size(state%boundary) /= size(state%psi)) return
        if (size(state%area) /= size(state%triangles, 2)) return
        if (any(state%triangles < 1)) return
        if (any(state%triangles > size(state%psi))) return
        if (any(.not. ieee_is_finite(state%psi))) return
        if (any(.not. ieee_is_finite(state%points))) return
        if (any(.not. ieee_is_finite(state%area))) return
        if (any(state%area <= 0)) return
        tolerance = 1.0e-11_dp
        if (present(relative_tolerance)) tolerance = relative_tolerance
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0 .or. tolerance >= 1) return
        work_limit = 200000
        if (present(max_evaluations)) work_limit = max_evaluations
        if (work_limit < 1) return
        if (nf > huge(1)/ndof) return
        call initialize_triangle_level_workspace(workspace, ndof*nf, nf, &
            status, max_evaluations=work_limit)
        if (status /= 0) return
        allocate (values(ndof*nf), errors(ndof*nf), epsabs(ndof*nf), shifted_flux(nf))
        estimates = 0
        do e = 1, size(state%triangles, 2)
            local_psi(:ndof) = state%psi(state%triangles(:, e))
            shifted_flux = flux - local_psi(1)
            vertices = state%points(:, state%triangles(1:3, e))
            rmin = minval(vertices(1, :))
            if (rmin <= 0) then
                status = 1
                exit
            end if
            call triangle_level_coefficients(state%element_order, &
                local_psi(:ndof), coefficients, status)
            if (status /= 0) exit
            ! Identical per-component scale to a native unit FFprime coefficient.
            epsabs = max(tiny(1.0_dp), 0.01_dp*tolerance/rmin)
            call integrate_triangle_levels(coefficients, flux, integrand, epsabs, &
                0.01_dp*tolerance, workspace, &
                values, errors, status, message)
            if (status /= 0) exit
            do j = 1, nf
                do k = 1, ndof
                    node = state%triangles(k, e)
                    index = (j - 1)*ndof + k
                    loads(node, j) = loads(node, j) + 2*state%area(e)*values(index)
                    estimates(node, j) = estimates(node, j) + &
                        2*state%area(e)*errors(index)
                end do
            end do
        end do
        if (status == 0) then
            do j = 1, nf
                where (state%boundary)
                    loads(:, j) = 0
                    estimates(:, j) = 0
                end where
                if (norm2(estimates(:, j)) > &
                    tolerance*max(norm2(loads(:, j)), tiny(1.0_dp))) then
                    status = 3
                    exit
                end if
            end do
        end if
        if (status /= 0) then
            loads = 0
            estimates = huge(1.0_dp)
        end if
    contains
        subroutine integrand(point, output, callback_status)
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: output(:)
            integer, intent(out) :: callback_status
            real(dp) :: shape(6), r, psi_q, source, profile_weights(nf)
            real(dp) :: zero_profile(nf), nodes_bar(nf), psi_bar
            integer :: column
            call element_values(state%element_order, point(1), point(2), shape(:ndof))
            r = dot_product(vertices(1, :), &
                [1 - point(1) - point(2), point(1), point(2)])
            psi_q = dot_product(shape(:ndof), local_psi(:ndof) - local_psi(1))
            zero_profile = 0
            ! The existing interpolation adjoint supplies the authoritative hat basis.
            call profile_adjoint(shifted_flux, zero_profile, psi_q, 1.0_dp, &
                psi_bar, nodes_bar, profile_weights)
            call gs_source(r, 1.0_dp, 0.0_dp, 1.0_dp, source)
            output = 0
            do column = 1, nf
                if (profile_weights(column) == 0) cycle
                output((column - 1)*ndof + 1:column*ndof) = &
                    source*profile_weights(column)*shape(:ndof)
            end do
            callback_status = 0
            if (any(.not. ieee_is_finite(output))) callback_status = 1
        end subroutine integrand
    end subroutine gs_profile_ff_load_basis



    subroutine prepare_gs_profile_response_cache(state, flux, ff, cache, status)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:), ff(:)
        type(gs_profile_response_cache_t), intent(out) :: cache
        integer, intent(out) :: status
        real(dp), allocatable :: mass(:, :, :)
        real(dp) :: axis(2), axis_flux, hessian(2, 2), local(2), gradients(2, 6)
        integer :: ndof
        status = 1
        if (size(flux) < 2 .or. size(ff) /= size(flux)) return
        if (.not. all(ieee_is_finite(flux))) return
        if (.not. all(ieee_is_finite(ff))) return
        if (any(flux(2:) <= flux(:size(flux) - 1))) return
        if (.not. allocated(state%psi)) return
        if (.not. allocated(state%points)) return
        if (.not. allocated(state%triangles)) return
        if (.not. allocated(state%gradient)) return
        if (.not. allocated(state%stiffness)) return
        if (.not. allocated(state%area)) return
        if (.not. allocated(state%diagonal)) return
        if (.not. allocated(state%boundary)) return
        ndof = size(state%triangles, 1)
        if (state%element_order == 1) then
            if (ndof /= 3) return
        else if (state%element_order == 2) then
            if (ndof /= 6) return
        else
            return
        end if
        if (size(state%points, 1) /= 2) return
        if (size(state%points, 2) /= size(state%psi)) return
        if (size(state%diagonal) /= size(state%psi)) return
        if (size(state%boundary) /= size(state%psi)) return
        if (size(state%area) /= size(state%triangles, 2)) return
        if (any(shape(state%gradient) /= [2, 3, size(state%area)])) return
        if (any(shape(state%stiffness) /= [ndof, ndof, size(state%area)])) return
        if (any(state%triangles < 1)) return
        if (any(state%triangles > size(state%psi))) return
        if (any(.not. ieee_is_finite(state%psi))) return
        if (any(.not. ieee_is_finite(state%points))) return
        if (any(.not. ieee_is_finite(state%gradient))) return
        if (any(.not. ieee_is_finite(state%stiffness))) return
        if (any(.not. ieee_is_finite(state%area))) return
        if (any(.not. ieee_is_finite(state%diagonal))) return
        if (any(state%area <= 0) .or. any(state%diagonal <= 0)) return
        call gs_axis_candidate(state, axis, axis_flux, hessian, &
            cache%axis_element, status)
        if (status /= 0) return
        status = 1
        if (axis_flux >= 0) return
        if (abs(flux(1) - axis_flux) > 1024*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(state%psi)))) return
        ndof = size(state%triangles, 1)
        allocate (cache%axis_shape(ndof))
        local = matmul(transpose(state%gradient(:, 2:3, cache%axis_element)), &
            axis - state%points(:, state%triangles(1, cache%axis_element)))
        call element_shapes(state, cache%axis_element, local(1), local(2), &
            cache%axis_shape, gradients(:, :ndof))
        call gs_profile_tangent_cache(state, flux, ff, axis_flux, mass, &
            cache%axis_response, status)
        if (status /= 0) return
        cache%state = state
        cache%flux = flux
        cache%ff = ff
        cache%matrix = state%stiffness - mass
        call gs_profile_ff_load_basis(state, flux, cache%load_basis, &
            cache%load_error, status)
        if (status /= 0) return
        cache%ready = .true.
    end subroutine prepare_gs_profile_response_cache

    logical function profile_cache_matches(state, flux, ff, cache) result(matches)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:), ff(:)
        type(gs_profile_response_cache_t), intent(in) :: cache
        matches = .false.
        if (.not. cache%ready) return
        if (state%element_order /= cache%state%element_order) return
        if (size(flux) /= size(cache%flux) .or. size(ff) /= size(cache%ff)) return
        if (any(flux /= cache%flux) .or. any(ff /= cache%ff)) return
        if (.not. allocated(state%psi)) return
        if (size(state%psi) /= size(cache%state%psi)) return
        if (any(state%psi /= cache%state%psi)) return
        if (.not. allocated(state%points)) return
        if (any(shape(state%points) /= shape(cache%state%points))) return
        if (any(state%points /= cache%state%points)) return
        if (.not. allocated(state%triangles)) return
        if (any(shape(state%triangles) /= shape(cache%state%triangles))) return
        if (any(state%triangles /= cache%state%triangles)) return
        if (.not. allocated(state%gradient)) return
        if (any(shape(state%gradient) /= shape(cache%state%gradient))) return
        if (any(state%gradient /= cache%state%gradient)) return
        if (.not. allocated(state%stiffness)) return
        if (any(shape(state%stiffness) /= shape(cache%state%stiffness))) return
        if (any(state%stiffness /= cache%state%stiffness)) return
        if (.not. allocated(state%area)) return
        if (size(state%area) /= size(cache%state%area)) return
        if (any(state%area /= cache%state%area)) return
        if (.not. allocated(state%diagonal)) return
        if (size(state%diagonal) /= size(cache%state%diagonal)) return
        if (any(state%diagonal /= cache%state%diagonal)) return
        if (.not. allocated(state%boundary)) return
        if (size(state%boundary) /= size(cache%state%boundary)) return
        if (any(state%boundary .neqv. cache%state%boundary)) return
        matches = .true.
    end function profile_cache_matches

    subroutine gs_profile_implicit_response(state, flux, ff, dff, &
            response, status, cache)
        type(gs_state_t), intent(in) :: state
        real(dp), intent(in) :: flux(:), ff(:), dff(:)
        real(dp), allocatable, intent(out) :: response(:)
        integer, intent(out) :: status
        type(gs_profile_response_cache_t), intent(in), optional :: cache
        type(gs_profile_response_cache_t) :: prepared
        status = 1
        if (size(dff) /= size(flux)) return
        if (.not. all(ieee_is_finite(dff))) return
        if (present(cache)) then
            if (.not. profile_cache_matches(state, flux, ff, cache)) return
            call cached_profile_response(cache, dff, response, status)
        else
            call prepare_gs_profile_response_cache(state, flux, ff, prepared, status)
            if (status /= 0) return
            call cached_profile_response(prepared, dff, response, status)
        end if
    end subroutine gs_profile_implicit_response

    subroutine cached_profile_response(cache, dff, response, status)
        type(gs_profile_response_cache_t), intent(in) :: cache
        real(dp), intent(in) :: dff(:)
        real(dp), allocatable, intent(out) :: response(:)
        integer, intent(out) :: status
        real(dp), allocatable :: rhs(:), pp(:), zero_state(:)
        real(dp) :: residual
        integer :: iterations
        allocate (rhs(size(cache%state%psi)), pp(size(cache%flux)), &
            response(size(cache%state%psi)), &
            zero_state(size(cache%state%psi)))
        pp = 0; response = 0; zero_state = 0
        rhs = matmul(cache%load_basis, dff)
        if (any(.not. ieee_is_finite(rhs))) then
            status = 1
            return
        end if
        ! Conservative column estimates can lose a relative budget on cancellation.
        ! Reintegrate that direction with the original API and unchanged tolerance.
        if (norm2(matmul(cache%load_error, abs(dff))) > &
            1.0e-11_dp*max(norm2(rhs), tiny(1.0_dp))) then
            call gs_profile_load_jvp(cache%state, 1.0_dp, cache%flux, pp, cache%ff, &
                zero_state, 0.0_dp, pp, pp, dff, rhs, status)
            if (status /= 0) return
        end if
        where (cache%state%boundary) rhs = 0
        rhs = rhs/cache%state%diagonal
        call real_gmres_operator(matvec, rhs, response, 1.0e-10_dp, &
            max(500, 5*size(rhs)), min(80, size(rhs)), &
            status, iterations, residual)
    contains
        subroutine matvec(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            real(dp) :: axis_tangent
            integer :: cell, first, second, ni, nj
            output = 0
            axis_tangent = dot_product(cache%axis_shape, &
                input(cache%state%triangles(:, cache%axis_element)))
            do cell = 1, size(cache%state%triangles, 2)
                do second = 1, size(cache%state%triangles, 1)
                    nj = cache%state%triangles(second, cell)
                    if (cache%state%boundary(nj)) cycle
                    do first = 1, size(cache%state%triangles, 1)
                        ni = cache%state%triangles(first, cell)
                        if (cache%state%boundary(ni)) cycle
                        output(ni) = output(ni) + &
                            cache%matrix(first, second, cell)*input(nj)
                    end do
                end do
            end do
            output = output + cache%axis_response*axis_tangent
            where (cache%state%boundary) output = input
            output = output/cache%state%diagonal
        end subroutine matvec
    end subroutine cached_profile_response



   real(dp) function gs_profile_value(flux, values, psi) result(value)
      real(dp), intent(in) :: flux(:), values(:), psi
      real(dp) :: weights(2)
      integer :: j
      ! The declared input law is piecewise linear and constant beyond the
      ! first/last supplied flux. Extension is explicit in the result claim.
      if (psi <= flux(1)) then
         value = values(1)
         return
      end if
      do j = 1, size(flux) - 1
         if (psi <= flux(j + 1)) then
            call lagrange_weights(2,psi,flux(j:j+1),weights)
            value = dot_product(weights,values(j:j+1))
            return
         end if
      end do
      value = values(size(flux))
   end function gs_profile_value

   real(dp) function gs_profile_integral(flux, values, psi) result(value)
      real(dp), intent(in) :: flux(:), values(:), psi
      real(dp) :: upper, width, segment
      integer :: j
      ! Integral from the supplied edge to psi of the executed derivative law.
      ! Each linear segment is integrated exactly; endpoints extend constantly.
      value = 0.0_dp
      if (psi < flux(1)) value = (psi - flux(1))*values(1)
      do j = 1, size(flux) - 1
         upper = max(flux(j), min(psi, flux(j + 1)))
         width = upper - flux(j)
         call gs_segment_primitive(width, values(j), &
            (values(j + 1) - values(j))/(flux(j + 1) - flux(j)), segment)
         value = value + segment
         value = value - 0.5_dp*(flux(j + 1) - flux(j))* &
                 (values(j) + values(j + 1))
      end do
      if (psi > flux(size(flux))) value = value + &
                                          (psi - flux(size(flux)))*values(size(flux))
   end function gs_profile_integral

   subroutine solve_gs_profiles(state, mu0, flux, pp, ff, status, relative_tolerance)
      type(gs_state_t), intent(inout) :: state
      real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: relative_tolerance
      real(dp), allocatable :: rhs(:), solution(:), action(:)
      real(dp) :: tolerance
      integer :: iteration
      status = 1
      if (mu0 <= 0.0_dp) return
      if (.not. ieee_is_finite(mu0)) return
      tolerance = 1.0e-8_dp
      if (present(relative_tolerance)) tolerance = relative_tolerance
      if (.not. ieee_is_finite(tolerance)) return
      if (tolerance <= 0 .or. tolerance >= 1) return
      if (size(flux) < 2) return
      if (size(pp) /= size(flux) .or. size(ff) /= size(flux)) return
      if (any(flux(2:) <= flux(:size(flux) - 1))) return
      if (any(.not. ieee_is_finite(flux))) return
      if (any(.not. ieee_is_finite(pp))) return
      if (any(.not. ieee_is_finite(ff))) return
      allocate (rhs(size(state%psi)), action(size(state%psi)))
      do iteration = 1, 1000
         call gs_profile_load(state, mu0, flux, pp, ff, rhs, status, &
              relative_tolerance=min(1.0e-11_dp,0.01_dp*tolerance))
         if (status /= 0) return
         call gs_linear_response(state, rhs, solution, status)
         if (status /= KRYLOV_OK) return
         where (.not. state%boundary)
            state%psi = 0.3_dp*solution + 0.7_dp*state%psi
         end where
         call gs_profile_load(state, mu0, flux, pp, ff, rhs, status, &
              relative_tolerance=min(1.0e-11_dp,0.01_dp*tolerance))
         if (status /= 0) return
         call apply_gs(state, state%psi, action)
         state%residual = norm2(action - rhs)/max(norm2(rhs), tiny(1.0_dp))
         state%iterations = iteration
         if (state%residual <= tolerance) then
            status = 0
            return
         end if
      end do
      status = 2 ! Candidate did not solve the declared nonlinear model.
   end subroutine solve_gs_profiles

   subroutine gs_element_jet(state, element, r, z, psi, grad, hessian, status)
      type(gs_state_t), intent(in) :: state
      integer, intent(in) :: element
      real(dp), intent(in) :: r, z
      real(dp), intent(out) :: psi, grad(2), hessian(2, 2)
      integer, intent(out) :: status
      real(dp) :: shape(6), gradients(2, 6), hessians(2, 2, 6), local(2), differences(6)
      integer :: ndof, j
      psi = 0.0_dp; grad = 0.0_dp; hessian = 0.0_dp
      status = 1
      if (.not. ieee_is_finite(r) .or. .not. ieee_is_finite(z)) return
      if (r <= 0) return
      if (element < 1 .or. element > size(state%triangles, 2)) return
      local = matmul(transpose(state%gradient(:, 2:3, element)), &
                     [r, z]-state%points(:, state%triangles(1, element)))
      if (minval([1-local(1)-local(2), local]) < -1.0e-12_dp) return
      ndof = size(state%triangles, 1)
      call element_shapes(state, element, local(1), local(2), shape(:ndof), &
                          gradients(:, :ndof), hessians(:, :, :ndof))
      psi = dot_product(shape(:ndof), state%psi(state%triangles(:, element)))
      ! Physical derivatives annihilate a constant gauge. Subtract it before
      ! composition, avoiding cancellation amplified by inverse cell sizes.
      differences(:ndof) = state%psi(state%triangles(:, element)) &
                           - state%psi(state%triangles(1, element))
      grad = matmul(gradients(:, :ndof), differences(:ndof))
      do j = 1, ndof
         hessian = hessian + hessians(:, :, j)*differences(j)
      end do
      status = 0
   end subroutine gs_element_jet

   subroutine sample_gs(state, r, z, f, psi, br, bz, bt, status)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: r, z, f
      real(dp), intent(out) :: psi, br, bz, bt
      integer, intent(out) :: status
      real(dp) :: grad(2), hessian(2, 2)
      integer :: e
      psi = 0.0_dp; br = 0.0_dp; bz = 0.0_dp; bt = 0.0_dp
      status = 1
      if (.not. ieee_is_finite(f)) return
      do e = 1, size(state%triangles, 2)
         call gs_element_jet(state, e, r, z, psi, grad, hessian, status)
         if (status /= 0) cycle
         call gs_field(r, grad(1), grad(2), f, br, bz, bt)
         return
      end do
   end subroutine sample_gs

   subroutine gs_classical_current(state, element, r, z, mu0, current, status)
      type(gs_state_t), intent(in) :: state
      integer, intent(in) :: element
      real(dp), intent(in) :: r, z, mu0
      real(dp), intent(out) :: current
      integer, intent(out) :: status
      real(dp) :: psi, grad(2), hessian(2, 2)
      status = 1; current = 0.0_dp
      if (.not. ieee_is_finite(mu0)) return
      if (mu0 <= 0) return
      call gs_element_jet(state, element, r, z, psi, grad, hessian, status)
      if (status /= 0) return
      call gs_current_density(r, mu0, grad(1), hessian(1, 1), hessian(2, 2), current)
   end subroutine gs_classical_current

   subroutine gs_axis_safety_factor(r, f_axis, hessian, q, status)
      real(dp), intent(in) :: r, f_axis, hessian(2, 2)
      real(dp), intent(out) :: q
      integer, intent(out) :: status
      real(dp) :: determinant
      status = 1; q = 0.0_dp
      if (.not. all(ieee_is_finite([r, f_axis]))) return
      if (.not. all(ieee_is_finite(hessian))) return
      if (r <= 0 .or. abs(f_axis) <= tiny(f_axis)) return
      if (abs(hessian(1, 2)-hessian(2, 1)) > &
          64*epsilon(r)*maxval(abs(hessian))) return
      determinant = hessian(1, 1)*hessian(2, 2)-hessian(1, 2)**2
      if (hessian(1, 1) <= 0 .or. determinant <= 0) return
      call gs_smooth_axis_q(r, f_axis, hessian(1, 1), &
                            hessian(1, 2), hessian(2, 2), q)
      if (.not. ieee_is_finite(q)) return
      status = 0
   end subroutine gs_axis_safety_factor

   subroutine gs_axis_candidate(state, axis, psi_axis, hessian_axis, element, status)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(out) :: axis(2), psi_axis, hessian_axis(2, 2)
      integer, intent(out) :: element, status
      real(dp) :: origin(2), point(2), grad(2), hessian(2, 2), matrix(2, 2)
      real(dp) :: displacement(2, 1), shape(3), value, determinant
      real(dp) :: edge_min, tangent(2), edge_curve, edge_slope, t, tolerance, second_min
      integer :: e, j, k, pivot(2), stat
      axis = 0.0_dp; psi_axis = huge(1.0_dp); hessian_axis = 0.0_dp
      element = 0; status = 3
      if (state%element_order /= 2) return
      status = 4 ! No regular strict interior minimum; edge/vertex cusp possible.
      edge_min = minval(state%psi)
      second_min = huge(1.0_dp)
      do e = 1, size(state%triangles, 2)
         origin = state%points(:, state%triangles(1, e))
         call gs_element_jet(state, e, origin(1), origin(2), value, grad, hessian, stat)
         if (stat /= 0) return
         ! Each edge is a quadratic restriction. Include its interior extrema.
         do j = 1, 3
            k = mod(j, 3)+1
            point = state%points(:, state%triangles(j, e))
            tangent = state%points(:, state%triangles(k, e))-point
            call gs_element_jet(state, e, point(1), point(2), &
                                value, grad, hessian, stat)
            if (stat /= 0) return
            edge_curve = dot_product(tangent, matmul(hessian, tangent))
            edge_slope = dot_product(tangent, grad)
            if (edge_curve <= 0) cycle
            t = -edge_slope/edge_curve
            if (t <= 0 .or. t >= 1) cycle
            point = point+t*tangent
            call gs_element_jet(state, e, point(1), point(2), &
                                value, grad, hessian, stat)
            if (stat /= 0) return
            edge_min = min(edge_min, value)
         end do
         call gs_element_jet(state, e, origin(1), origin(2), value, grad, hessian, stat)
         determinant = hessian(1, 1)*hessian(2, 2)-hessian(1, 2)**2
         if (hessian(1, 1) <= 0 .or. determinant <= 0) cycle
         matrix = hessian; displacement(:, 1) = -grad
         call dgesv(2, 1, matrix, 2, pivot, displacement, 2, stat)
         if (stat /= 0) cycle
         point = origin+displacement(:, 1)
         shape(2:3) = matmul(transpose(state%gradient(:, 2:3, e)), point-origin)
         shape(1) = 1-sum(shape(2:3))
         if (minval(shape) <= 1.0e-10_dp) cycle
         call gs_element_jet(state, e, point(1), point(2), value, grad, hessian, stat)
         if (stat /= 0) cycle
         if (value >= psi_axis) then
            second_min = min(second_min, value)
            cycle
         end if
         second_min = psi_axis
         axis = point; psi_axis = value; hessian_axis = hessian; element = e
      end do
      if (element == 0) return
      tolerance = 256*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(state%psi)))
      if (edge_min <= psi_axis+tolerance) then
         element = 0
         return
      end if
      if (second_min <= psi_axis+tolerance) then
         element = 0; status = 5 ! Multiple candidate axes; do not choose silently.
         return
      end if
      status = 0
   end subroutine gs_axis_candidate

   subroutine gs_current_jump(state, left, right, r, z, normal, mu0, current, status)
      type(gs_state_t), intent(in) :: state
      integer, intent(in) :: left, right
      real(dp), intent(in) :: r, z, normal(2), mu0
      real(dp), intent(out) :: current
      integer, intent(out) :: status
      real(dp) :: psi, gl(2), gr(2), hessian(2, 2), tangent(2), toward(2)
      integer :: j, k, first, second, edge
      status = 1; current = 0.0_dp
      if (.not. all(ieee_is_finite([normal, mu0]))) return
      if (mu0 <= 0) return
      if (abs(norm2(normal)-1) > 1.0e-10_dp) return
      if (left < 1 .or. left > size(state%triangles, 2)) return
      if (right < 1 .or. right > size(state%triangles, 2)) return
      edge = 0
      do j = 1, 3
         if (state%neighbors(j, left) == right) edge = j
      end do
      if (edge == 0) return
      k = mod(edge, 3)+1
      first = state%triangles(edge, left); second = state%triangles(k, left)
      tangent = state%points(:, second)-state%points(:, first)
      toward = sum(state%points(:, state%triangles(1:3, right)), dim=2) &
             - sum(state%points(:, state%triangles(1:3, left)), dim=2)
      if (abs(dot_product(normal, tangent)) > 1.0e-10_dp*norm2(tangent)) return
      if (dot_product(normal, toward) <= 0) return
      if (abs(dot_product(normal, [r, z]-state%points(:, first))) > &
          1.0e-10_dp*norm2(tangent)) return
      call gs_element_jet(state, left, r, z, psi, gl, hessian, status)
      if (status /= 0) return
      call gs_element_jet(state, right, r, z, psi, gr, hessian, status)
      if (status /= 0) return
      call gs_current_sheet(r, mu0, gr(1)-gl(1), gr(2)-gl(2), &
                            normal(1), normal(2), current)
   end subroutine gs_current_jump

   subroutine write_gs_mesh(state, path, status)
      type(gs_state_t), intent(in) :: state
      character(*), intent(in) :: path
      integer, intent(out) :: status
      integer :: unit, j
      open (newunit=unit, file=path, status='replace', action='write', &
            iostat=status)
      if (status /= 0) return
      if (state%element_order == 1) then
         write (unit, '(a)') '# kin6d-gs-mesh-v1: SI R[m] Z[m] psi[Wb/rad] boundary'
      else
         write (unit, '(a)') '# kin6d-gs-mesh-v2: P2 SI R[m] Z[m] psi[Wb/rad] boundary'
      end if
      write (unit, *) size(state%psi), size(state%triangles, 2)
      do j = 1, size(state%psi)
         write (unit, '(3es25.16e3,1x,i1)') state%points(:, j), state%psi(j), &
            merge(1, 0, state%boundary(j))
      end do
      write (unit, '(a)') '# triangle node indices (one-based)'
      do j = 1, size(state%triangles, 2)
         write (unit, '(*(i10))') state%triangles(:, j)
      end do
      close (unit)
   end subroutine write_gs_mesh
end module kin6d_grad_shafranov
