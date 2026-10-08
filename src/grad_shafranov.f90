module kin6d_grad_shafranov
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use triangle_compat, only: tc_result_t, triangulate_compat
   use basis_p1_2d_module, only: basis_p1_2d_t
   use basis_p2_2d_module, only: basis_p2_2d_t
   use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
   use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
   use kin6d_grad_shafranov_leaves, only: gs_source, gs_field, gs_segment_primitive, &
       gs_current_density, gs_current_sheet, gs_smooth_axis_q
   implicit none
   private
   public :: gs_state_t, initialize_gs, solve_gs, gs_linear_response, sample_gs
   public :: write_gs_mesh
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

   subroutine assemble_gs_rhs(state, mu0, pprime, ffprime, rhs, &
                              flux, pp, ff)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: mu0, pprime, ffprime
      real(dp), intent(out) :: rhs(:)
      real(dp), intent(in), optional :: flux(:), pp(:), ff(:)
      type(gauss_quadrature_triangle_t) :: quadrature
      real(dp) :: point(2), shape(6), gradients(2, 6), source, weight
      real(dp) :: psi_q, pp_q, ff_q
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
            pp_q = pprime; ff_q = ffprime
            if (present(flux)) then
               psi_q = dot_product(shape(:ndof), state%psi(state%triangles(:, e)))
               pp_q = gs_profile_value(flux, pp, psi_q)
               ff_q = gs_profile_value(flux, ff, psi_q)
            end if
            call gs_source(point(1), mu0, pp_q, ff_q, source)
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

   real(dp) function gs_profile_value(flux, values, psi) result(value)
      real(dp), intent(in) :: flux(:), values(:), psi
      real(dp) :: weight
      integer :: j
      ! The declared input law is piecewise linear and constant beyond the
      ! first/last supplied flux. Extension is explicit in the result claim.
      if (psi <= flux(1)) then
         value = values(1)
         return
      end if
      do j = 1, size(flux) - 1
         if (psi <= flux(j + 1)) then
            weight = (psi - flux(j))/(flux(j + 1) - flux(j))
            value = (1 - weight)*values(j) + weight*values(j + 1)
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
         call assemble_gs_rhs(state, mu0, 0.0_dp, 0.0_dp, rhs, flux, pp, ff)
         call gs_linear_response(state, rhs, solution, status)
         if (status /= KRYLOV_OK) return
         where (.not. state%boundary)
            state%psi = 0.3_dp*solution + 0.7_dp*state%psi
         end where
         call assemble_gs_rhs(state, mu0, 0.0_dp, 0.0_dp, rhs, flux, pp, ff)
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
