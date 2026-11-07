! Element-wise Gauss quadrature of B-spline products (Galerkin integrals).
# include <petsc/finclude/petscsys.h>

module quad
   use types
   use bspline
   implicit none

   private

   public :: quad_x_B
   public :: quad_y_B
   public :: quad_z_B
   public :: quad_x_BB
   public :: quad_y_BB
   public :: quad_z_BB
   public :: quad_x_BBB
   public :: quad_y_BBB
   public :: quad_z_BBB
   public :: quad_x_BBBB
   public :: quad_y_BBBB
   public :: quad_z_BBBB
   public :: quad_xyz_B
   public :: quad_reciprocal_x_B
   public :: quad_reciprocal_xyz_B
contains

   function quad_reciprocal_x_B(grid, B, a0, b0, Bi, derivi, vvx) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, derivi
      PetscReal :: a0, b0, gauss_pt
      PetscScalar, dimension(:), pointer :: vvx

      PetscInt :: i, j, sx
      PetscReal :: func, num, gauss_point

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      if (.not. allocated(B%intx)) allocate(B%intx(deg+1))

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)

         sx = Int(gauss_pt / grid%dx)

         do j=sx, sx+deg
            B%intx(j-sx+1) = xBJK(grid, B, j+1, deg, gauss_pt)
         end do

         func = 1.0 / (DOT_PRODUCT(B%intx, vvx(sx-deg:sx)))

         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               xdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * func
      end do

      end associate
   end function quad_reciprocal_x_B

   function quad_reciprocal_xyz_B(ctx, grid, B, x0, x1, y0, y1, z0, z1, Bi, Bj, Bk, &
         derivi, derivj, derivk, vv, ierr) result (num)
      use linalg
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bk, derivi, derivj, derivk
      PetscReal :: x0, x1, y0, y1, z0, z1
      PetscScalar, dimension(:, :, :), pointer :: vv
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscReal :: num, func
      PetscReal, dimension(3) :: x

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do k = 1, nodes
         do j = 1, nodes
            do i = 1, nodes
               x(1) = 0.5 * (x1 - x0) * xi(i) + 0.5 * (x1 + x0)
               x(2) = 0.5 * (y1 - y0) * xi(j) + 0.5 * (y1 + y0)
               x(3) = 0.5 * (z1 - z0) * xi(k) + 0.5 * (z1 + z0)

               call interpolate(ctx, grid, B, func, x, vv, ierr)

               func = 1.0 / func

               num = num + &
                  0.5 * (x1 - x0) * wts(i) * &
                  0.5 * (y1 - y0) * wts(j) * &
                  0.5 * (z1 - z0) * wts(k) * &
                  xdBIJK(grid, B, derivi, Bi, deg, x(1)) * &
                  ydBIJK(grid, B, derivj, Bj, deg, x(2)) * &
                  zdBIJK(grid, B, derivk, Bk, deg, x(3)) * func
            end do
         end do
      end do

      end associate
   end function quad_reciprocal_xyz_B

   function quad_xyz_B(grid, B, x0, x1, y0, y1, z0, z1, Bi, Bj, Bk, &
         derivi, derivj, derivk, ptr, t) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bk, derivi, derivj, derivk
      PetscReal :: x0, x1, y0, y1, z0, z1, t
      procedure (func3D), pointer :: ptr

      PetscInt :: i, j, k
      PetscReal :: num
      PetscReal, dimension(3) :: x

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do k = 1, nodes
         do j = 1, nodes
            do i = 1, nodes
               x(1) = 0.5 * (x1 - x0) * xi(i) + 0.5 * (x1 + x0)
               x(2) = 0.5 * (y1 - y0) * xi(j) + 0.5 * (y1 + y0)
               x(3) = 0.5 * (z1 - z0) * xi(k) + 0.5 * (z1 + z0)

               num = num + &
                  0.5 * (x1 - x0) * wts(i) * &
                  0.5 * (y1 - y0) * wts(j) * &
                  0.5 * (z1 - z0) * wts(k) * &
                  xdBIJK(grid, B, derivi, Bi, deg, x(1)) * &
                  ydBIJK(grid, B, derivj, Bj, deg, x(2)) * &
                  zdBIJK(grid, B, derivk, Bk, deg, x(3)) * ptr(x(1), x(2), x(3), t)
            end do
         end do
      end do

      end associate
   end function quad_xyz_B

   function quad_x_B(grid, B, a0, b0, Bi, derivi, ptr, t) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, derivi
      PetscReal :: a0, b0, t
      procedure (func1D), pointer :: ptr

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               xdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * ptr(gauss_pt, t)
      end do

      end associate
   end function quad_x_B

   function quad_y_B(grid, B, a0, b0, Bi, derivi, ptr, t) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, derivi
      PetscReal :: a0, b0, t
      procedure (func1D), pointer :: ptr

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               ydBIJK(grid, B, derivi, Bi, deg, gauss_pt) * ptr(gauss_pt, t)
      end do

      end associate
   end function quad_y_B

   function quad_z_B(grid, B, a0, b0, Bi, derivi, ptr, t) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, derivi
      PetscReal :: a0, b0, t
      procedure (func1D), pointer :: ptr

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               zdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * ptr(gauss_pt, t)
      end do

      end associate
   end function quad_z_B

   function quad_x_BB(grid, B, a0, b0, Bi, Bj, derivi, derivj) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, derivi, derivj
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               xdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               xdBIJK(grid, B, derivj, Bj, deg, gauss_pt)
      end do

      end associate
   end function quad_x_BB

   function quad_y_BB(grid, B, a0, b0, Bi, Bj, derivi, derivj) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, derivi, derivj
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               ydBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               ydBIJK(grid, B, derivj, Bj, deg, gauss_pt)
      end do

      end associate
   end function quad_y_BB

   function quad_z_BB(grid, B, a0, b0, Bi, Bj, derivi, derivj) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, derivi, derivj
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal ::num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               zdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               zdBIJK(grid, B, derivj, Bj, deg, gauss_pt)
      end do

      end associate
   end function quad_z_BB

   function quad_x_BBB(grid, B, a0, b0, Bi, Bj, Bl, derivi, derivj, derivl) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, derivi, derivj, derivl
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               xdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               xdBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               xdBIJK(grid, B, derivl, Bl, deg, gauss_pt)
      end do

      end associate
   end function quad_x_BBB

   function quad_y_BBB(grid, B, a0, b0, Bi, Bj, Bl, derivi, derivj, derivl) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, derivi, derivj, derivl
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               ydBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               ydBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               ydBIJK(grid, B, derivl, Bl, deg, gauss_pt)
      end do

      end associate
   end function quad_y_BBB

   function quad_z_BBB(grid, B, a0, b0, Bi, Bj, Bl, derivi, derivj, derivl) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, derivi, derivj, derivl
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               zdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               zdBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               zdBIJK(grid, B, derivl, Bl, deg, gauss_pt)
      end do

      end associate
   end function quad_z_BBB

   function quad_x_BBBB(grid, B, a0, b0, Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               xdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               xdBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               xdBIJK(grid, B, derivl, Bl, deg, gauss_pt) * &
               xdBIJK(grid, B, derivm, Bm, deg, gauss_pt)
      end do

      end associate
   end function quad_x_BBBB

   function quad_y_BBBB(grid, B, a0, b0, Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               ydBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               ydBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               ydBIJK(grid, B, derivl, Bl, deg, gauss_pt) * &
               ydBIJK(grid, B, derivm, Bm, deg, gauss_pt)
      end do

      end associate
   end function quad_y_BBBB

   function quad_z_BBBB(grid, B, a0, b0, Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm) result (num)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: Bi, Bj, Bl, Bm, derivi, derivj, derivl, derivm
      PetscReal :: a0, b0

      PetscInt :: i
      PetscReal :: num, gauss_pt

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg)

      num = 0.0
      do i = 1, nodes
         gauss_pt = 0.5 * (b0 - a0) * xi(i) + 0.5 * (b0 + a0)
         num = num + 0.5 * (b0 - a0) * wts(i) *      &
               zdBIJK(grid, B, derivi, Bi, deg, gauss_pt) * &
               zdBIJK(grid, B, derivj, Bj, deg, gauss_pt) * &
               zdBIJK(grid, B, derivl, Bl, deg, gauss_pt) * &
               zdBIJK(grid, B, derivm, Bm, deg, gauss_pt)
      end do

      end associate
   end function quad_z_BBBB

end module quad

! vim: set ts=3 sw=3 sts=3:
