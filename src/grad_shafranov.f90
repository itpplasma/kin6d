module kin6d_grad_shafranov
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use triangle_compat, only: tc_result_t, triangulate_compat
   use basis_p1_2d_module, only: basis_p1_2d_t
   use fortnum_krylov, only: real_conjugate_gradient_operator, KRYLOV_OK
   use kin6d_grad_shafranov_leaves, only: gs_source, gs_field, gs_segment_primitive
   implicit none
   private
   public :: gs_state_t, initialize_gs, solve_gs, gs_linear_response, sample_gs
   public :: write_gs_mesh
   public :: solve_gs_profiles, gs_profile_value, gs_profile_integral
   ! Fixed-boundary P1 weak GS in physical (R,Z), R>0. The domain is a
   ! polygonal approximation; nested flux surfaces are not solver inputs.
   ! A result is a candidate, NOT a continuous certificate. Geometry,
   ! discretization and source-model errors remain separate open budgets.
   type :: gs_state_t
      real(dp), allocatable :: points(:, :), psi(:)
      integer, allocatable :: triangles(:, :)
      logical, allocatable :: boundary(:)
      real(dp), allocatable :: stiffness(:, :, :), gradient(:, :, :)
      real(dp), allocatable :: area(:), diagonal(:)
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
   subroutine initialize_gs(state, contour, max_area, status)
      type(gs_state_t), intent(out) :: state
      real(dp), intent(in) :: contour(:, :), max_area
      integer, intent(out) :: status
      type(tc_result_t) :: mesh
      type(basis_p1_2d_t) :: basis
      real(dp) :: jac(2, 2), matrix(2, 2), det, grads(2, 3)
      real(dp) :: point(2), basis_q(3), weight
      real(dp) :: quadrature(3, 3)
      integer, allocatable :: segments(:, :)
      integer :: j, e, k, l, q, pivot(2), stat, n
      status = 1
      if (size(contour, 1) /= 2) return
      n = size(contour, 2)
      if (n < 3 .or. max_area <= 0.0_dp) return
      if (any(.not. ieee_is_finite(contour))) return
      if (minval(contour(1, :)) <= 0.0_dp) return
      allocate (segments(2, n))
      do j = 1, n
         segments(:, j) = [j, mod(j, n) + 1]
      end do
      call triangulate_compat(contour, segments, reshape([real(dp) ::], [2, 0]), &
                           mesh, stat, min_angle=20.0_dp, max_area=max_area, nobisect=0)
      if (stat /= 0 .or. mesh%ntriangles < 1) return
      state%points = mesh%points
      state%triangles = mesh%triangles
      state%boundary = mesh%pointmarkers /= 0
      allocate (state%psi(mesh%npoints), state%diagonal(mesh%npoints), &
                state%stiffness(3, 3, mesh%ntriangles), &
                state%gradient(2, 3, mesh%ntriangles), state%area(mesh%ntriangles))
      state%psi = 0.0_dp; state%diagonal = 0.0_dp
      quadrature = 1.0_dp/6.0_dp
      do j = 1, 3
         quadrature(j, j) = 2.0_dp/3.0_dp
      end do
      do e = 1, mesh%ntriangles
         call basis%compute_jacobian(state%points(:, state%triangles(:, e)), &
                                     jac, det)
         if (abs(det) <= tiny(det)) return
         matrix = transpose(jac)
         do j = 1, 3
            grads(:, j) = basis%grad(j, 0.0_dp, 0.0_dp)
         end do
         call dgesv(2, 3, matrix, 2, pivot, grads, 2, stat)
         if (stat /= 0) return
         state%gradient(:, :, e) = grads
         state%area(e) = 0.5_dp*abs(det)
         weight = 0.0_dp
         do q = 1, 3
            basis_q = quadrature(:, q)
            point = matmul(state%points(:, state%triangles(:, e)), basis_q)
            weight = weight + state%area(e)/(3*point(1))
         end do
         do l = 1, 3
            do k = 1, 3
               state%stiffness(k, l, e) = weight*dot_product(grads(:, k), &
                                                             grads(:, l))
            end do
            j = state%triangles(l, e)
            state%diagonal(j) = state%diagonal(j) + &
                                state%stiffness(l, l, e)
         end do
      end do
      where (state%boundary) state%diagonal = 1.0_dp
      status = 0
   end subroutine initialize_gs

   subroutine apply_gs(state, input, output)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: input(:)
      real(dp), intent(out) :: output(:)
      integer :: e, k, l, i, j
      output = 0.0_dp
      do e = 1, size(state%triangles, 2)
         do l = 1, 3
            j = state%triangles(l, e)
            if (state%boundary(j)) cycle
            do k = 1, 3
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
      real(dp) :: point(2), shape(3), source, quadrature(3, 3)
      real(dp) :: psi_q, pp_q, ff_q
      integer :: e, q, k, l, i, j
      rhs = 0.0_dp
      quadrature = 1.0_dp/6.0_dp
      do q = 1, 3
         quadrature(q, q) = 2.0_dp/3.0_dp
      end do
      do e = 1, size(state%triangles, 2)
         do q = 1, 3
            shape = quadrature(:, q)
            point = matmul(state%points(:, state%triangles(:, e)), shape)
            pp_q = pprime; ff_q = ffprime
            if (present(flux)) then
               psi_q = dot_product(shape, state%psi(state%triangles(:, e)))
               pp_q = gs_profile_value(flux, pp, psi_q)
               ff_q = gs_profile_value(flux, ff, psi_q)
            end if
            call gs_source(point(1), mu0, pp_q, ff_q, source)
            do k = 1, 3
               i = state%triangles(k, e)
               rhs(i) = rhs(i) + state%area(e)*source*shape(k)/3
            end do
         end do
         do l = 1, 3
            j = state%triangles(l, e)
            if (.not. state%boundary(j)) cycle
            do k = 1, 3
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

   subroutine solve_gs_profiles(state, mu0, flux, pp, ff, status)
      type(gs_state_t), intent(inout) :: state
      real(dp), intent(in) :: mu0, flux(:), pp(:), ff(:)
      integer, intent(out) :: status
      real(dp), allocatable :: rhs(:), solution(:), action(:)
      integer :: iteration
      status = 1
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
         if (state%residual <= 1.0e-8_dp) then
            status = 0
            return
         end if
      end do
      status = 2 ! Candidate did not solve the declared nonlinear model.
   end subroutine solve_gs_profiles

   subroutine sample_gs(state, r, z, f, psi, br, bz, bt, status)
      type(gs_state_t), intent(in) :: state
      real(dp), intent(in) :: r, z, f
      real(dp), intent(out) :: psi, br, bz, bt
      integer, intent(out) :: status
      real(dp) :: point(2), shape(3), grad(2), a(2), b(2), c(2), det
      integer :: e
      psi = 0.0_dp; br = 0.0_dp; bz = 0.0_dp; bt = 0.0_dp
      status = 1
      point = [r, z]
      do e = 1, size(state%triangles, 2)
         a = state%points(:, state%triangles(1, e))
         b = state%points(:, state%triangles(2, e))
         c = state%points(:, state%triangles(3, e))
         det = (b(1) - a(1))*(c(2) - a(2)) - &
               (c(1) - a(1))*(b(2) - a(2))
         ! Basis gradients come from FortFEM. The barycentric reconstruction
         ! uses that same affine shape, without a private geometry solver.
         shape(2) = dot_product(state%gradient(:, 2, e), point - a)
         shape(3) = dot_product(state%gradient(:, 3, e), point - a)
         shape(1) = 1.0_dp - shape(2) - shape(3)
         if (minval(shape) < -1.0e-12_dp) cycle
         if (abs(det) <= tiny(det)) cycle
         psi = dot_product(shape, state%psi(state%triangles(:, e)))
         grad = matmul(state%gradient(:, :, e), &
                       state%psi(state%triangles(:, e)))
         call gs_field(r, grad(1), grad(2), f, br, bz, bt)
         status = 0
         return
      end do
   end subroutine sample_gs

   subroutine write_gs_mesh(state, path, status)
      type(gs_state_t), intent(in) :: state
      character(*), intent(in) :: path
      integer, intent(out) :: status
      integer :: unit, j
      open (newunit=unit, file=path, status='replace', action='write', &
            iostat=status)
      if (status /= 0) return
      write (unit, '(a)') '# kin6d-gs-mesh-v1: SI R[m] Z[m] psi[Wb/rad] boundary'
      write (unit, *) size(state%psi), size(state%triangles, 2)
      do j = 1, size(state%psi)
         write (unit, '(3es25.16e3,1x,i1)') state%points(:, j), state%psi(j), &
            merge(1, 0, state%boundary(j))
      end do
      write (unit, '(a)') '# triangle node indices (one-based)'
      do j = 1, size(state%triangles, 2)
         write (unit, '(3i10)') state%triangles(:, j)
      end do
      close (unit)
   end subroutine write_gs_mesh
end module kin6d_grad_shafranov
