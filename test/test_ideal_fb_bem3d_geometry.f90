module test_bem3d_geometry_charts
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use kin6d_equilibrium_chart, only: ideal_chart_t, cylindrical_chart_t
   implicit none
   type, extends(ideal_chart_t) :: manufactured_chart_t
   contains
      procedure :: jet => manufactured_jet
   end type
   type, extends(cylindrical_chart_t) :: shifted_chart_t
   contains
      procedure :: jet => shifted_jet
   end type
contains
   pure subroutine manufactured_jet(chart, q, x, d1, d2, d3)
      class(manufactured_chart_t), intent(in) :: chart
      real(dp), intent(in) :: q(3)
      real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
      select type (chart)
      class is (manufactured_chart_t)
         x = [q(1) + 0.2_dp*q(2), q(2), q(3) + 0.1_dp*q(1)*q(2)]
         d1 = 0; d2 = 0; d3 = 0
         d1(1, 1) = 1; d1(1, 2) = 0.2_dp; d1(2, 2) = 1
         d1(3, 1) = 0.1_dp*q(2); d1(3, 2) = 0.1_dp*q(1); d1(3, 3) = 1
         d2(3, 1, 2) = 0.1_dp; d2(3, 2, 1) = 0.1_dp
      end select
   end subroutine
   pure subroutine shifted_jet(chart, q, x, d1, d2, d3)
      class(shifted_chart_t), intent(in) :: chart
      real(dp), intent(in) :: q(3)
      real(dp), intent(out) :: x(3), d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
      type(cylindrical_chart_t) :: base
      select type (chart)
      class is (shifted_chart_t)
         call base%jet(q, x, d1, d2, d3)
         x(1) = x(1) + 0.2_dp*q(2); d1(1, 2) = d1(1, 2) + 0.2_dp
      end select
   end subroutine
end module

program test_ideal_fb_bem3d_geometry
   use test_bem3d_geometry_charts
   use kin6d_ideal_free_boundary_bem3d, only: ideal_fb_bem3d_surface_jet
   use kin6d_triangle_basis, only: triangle_basis_jet
   implicit none
   type(cylindrical_chart_t) :: cylindrical
   type(manufactured_chart_t) :: manufactured
   type(shifted_chart_t) :: shifted
   call check(cylindrical)
   call check(manufactured)
   call check(shifted)
   print '(a)', 'P2/P3 swept geometry: full-jet points/tangents/normals and finite differences PASS'
contains
   subroutine check(chart)
      class(ideal_chart_t), intent(in) :: chart
      real(dp) :: edge(2, 0:3), parameters(2, 3), uv(2), q(3), x(3), a(3), b(3)
      real(dp) :: xr(3), ar(3), br(3), d1(3, 3), d1r(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3)
      real(dp) :: plus(3), minus(3), unused_a(3), unused_b(3), shifted_uv(2), t
      real(dp), parameter :: delta = 2e-6_dp, tol = 128*epsilon(1.0_dp)
      integer :: order, panel, node, sample, k, status
      q = [1.7_dp, -0.3_dp, 0.4_dp]
      call chart%jet(q, xr, d1r, d2, d3)
      call chart%first_jet(q, x, d1)
      call require(norm2(x - xr) <= tol*(1 + norm2(xr)), 'first-jet position')
      call require(maxval(abs(d1 - d1r)) <= tol*(1 + maxval(abs(d1r))), 'first-jet Jacobian')
      do order = 2, 3
         do node = 0, order
            t = real(node, dp)/order
            edge(:, node) = [1.3_dp + 0.2_dp*t + 0.1_dp*t*t, &
                             -0.4_dp + 0.3_dp*t - 0.2_dp*t*t + 0.05_dp*t*t*t]
         end do
         do panel = 1, 3
            parameters = reshape([0.0_dp, 0.2_dp, 0.0_dp, 0.9_dp, 1.0_dp, 0.2_dp], [2, 3])
            if (panel == 2) parameters = parameters(:, [3, 2, 1])
            if (panel == 3) parameters = reshape([0.1_dp, -0.3_dp, 0.8_dp, -0.1_dp, &
                                                  0.2_dp, 0.9_dp], [2, 3])
            do sample = 1, 4
               uv = [0.1_dp*sample, 0.15_dp*sample]
               call ideal_fb_bem3d_surface_jet(chart, order, edge(:, :order), parameters, uv, x, a, b, status)
               call require(status == 0, 'geometry status')
               call reference(chart, order, edge(:, :order), parameters, uv, xr, ar, br)
               call require(norm2(x - xr) <= tol*(1 + norm2(xr)), 'surface position')
               call require(norm2(a - ar) <= tol*(1 + norm2(ar)), 'first tangent')
               call require(norm2(b - br) <= tol*(1 + norm2(br)), 'second tangent')
               call require(norm2(cross(a, b) - cross(ar, br)) <= &
                            tol*norm2(ar)*norm2(br), 'oriented normal')
               do k = 1, 2
                  shifted_uv = uv; shifted_uv(k) = uv(k) + delta
                  call ideal_fb_bem3d_surface_jet(chart, order, edge(:, :order), &
                                                  parameters, shifted_uv, plus, unused_a, unused_b, status)
                  shifted_uv(k) = uv(k) - delta
                  call ideal_fb_bem3d_surface_jet(chart, order, edge(:, :order), &
                                                  parameters, shifted_uv, minus, unused_a, unused_b, status)
                  if (k == 1) then
                     call require(norm2((plus - minus)/(2*delta) - a) < 1e-7_dp*(1 + norm2(a)), &
                                  'first tangent finite difference')
                  else
                     call require(norm2((plus - minus)/(2*delta) - b) < 1e-7_dp*(1 + norm2(b)), &
                                  'second tangent finite difference')
                  end if
               end do
            end do
         end do
      end do
   end subroutine
   subroutine reference(chart, order, edge, parameters, uv, x, a, b)
      class(ideal_chart_t), intent(in) :: chart
      integer, intent(in) :: order
      real(dp), intent(in) :: edge(2, 0:order), parameters(2, 3), uv(2)
      real(dp), intent(out) :: x(3), a(3), b(3)
      real(dp) :: parameter(2), deriv(2, 2), q(3), qd(3, 2), w(0:3), dw(0:3)
      real(dp) :: d1(3, 3), d2(3, 3, 3), d3(3, 3, 3, 3), values(10), gradient(2, 10)
      integer :: k, nloc, index(4)
      parameter = parameters(:, 1) + uv(1)*(parameters(:, 2) - parameters(:, 1)) + &
                  uv(2)*(parameters(:, 3) - parameters(:, 1))
      deriv(:, 1) = parameters(:, 2) - parameters(:, 1); deriv(:, 2) = parameters(:, 3) - parameters(:, 1)
      nloc = 6; index = [1, 4, 2, 0]
      if (order == 3) then
         nloc = 10; index = [1, 4, 5, 2]
      end if
      call triangle_basis_jet([parameter(1), 0.0_dp], values(:nloc), gradient(:, :nloc))
      w(:order) = values(index(:order + 1)); dw(:order) = gradient(1, index(:order + 1))
      q(1:2) = matmul(edge, w(:order)); q(3) = parameter(2)
      do k = 1, 2
         qd(1:2, k) = matmul(edge, dw(:order))*deriv(1, k)
         qd(3, k) = deriv(2, k)
      end do
      call chart%jet(q, x, d1, d2, d3)
      a = matmul(d1, qd(:, 1)); b = matmul(d1, qd(:, 2))
   end subroutine
   pure function cross(a, b) result(c)
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: c(3)
      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   end function
   subroutine require(ok, message)
      logical, intent(in) :: ok
      character(*), intent(in) :: message
      if (.not. ok) then
         print '(a)', message
         error stop 'swept BEM geometry equivalence'
      end if
   end subroutine
end program
