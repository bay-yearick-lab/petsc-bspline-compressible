! Assembly of 1D Galerkin mass, derivative and triple/quadruple product operators per direction.
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>

module operator
   use types
   use bspline
   use quad
   use gausstable
   implicit none

   private

   public :: mass_x
   public :: mass_y
   public :: mass_z
   public :: mass_init
   public :: linear_x
   public :: linear_y
   public :: linear_z
   public :: linear_init
   public :: nonlinear_x
   public :: nonlinear_y
   public :: nonlinear_z
   public :: nonlinear_init
contains

   subroutine nonlinear_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      call nonlinear_x(ctx, grid, B, ierr)
      call nonlinear_y(ctx, grid, B, ierr)
      call nonlinear_z(ctx, grid, B, ierr)

      ierr = 0
      return
   end subroutine nonlinear_init

   subroutine check_dirichlet_index(type, s, e, m, d, N, minsize, maxsize)
      character(len=*) :: type
      PetscInt :: s, e, m, d, N, minsize, maxsize

      PetscErrorCode :: ierr

      if (type == 'knot') then
         if (s.eq.1) then
            minsize = 0
         else
            minsize = s - d - 1
         end if

         if (e.eq.m) then
            maxsize = N - d + 2
         else
            maxsize = e + d + 2
         end if
      else if (type == 'column') then
         if (s.eq.1) then
            minsize = 0
         else
            minsize = s - d - 1
         end if

         if (e.eq.m) then
            maxsize = N - 1
         else
            maxsize = e + d - 1
         end if
      else
         call PetscPrintf(PETSC_COMM_WORLD, "operator: check_index: No such option.", ierr)
         stop
      end if
   end subroutine check_dirichlet_index

   subroutine check_periodic_index(type, s, e, m, d, N, minsize, maxsize)
      character(len=*) :: type
      PetscInt :: s, e, m, d, N, minsize, maxsize

      PetscErrorCode :: ierr

      if (type == 'knot') then
         if (s.eq.1) then
            minsize = 0
         else
            minsize = s - d - 1
         end if

         if (e.eq.m) then
            maxsize = N + d + 2
         else
            maxsize = e + d + 2
         end if
      else if (type == 'column') then
         if (s.eq.1) then
            minsize = 0
         else
            minsize = s - d - 1
         end if

         if (e.eq.m) then
            maxsize = N - 1
         else
            maxsize = e + d - 1
         end if
      else
         call PetscPrintf(PETSC_COMM_WORLD, "operator: check_index: No such option.", ierr)
         stop
      end if
   end subroutine check_periodic_index

   subroutine mass_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      call mass_x(ctx, grid, B, ierr)
      call mass_y(ctx, grid, B, ierr)
      call mass_z(ctx, grid, B, ierr)

      ierr = 0
      return
   end subroutine mass_init

   subroutine linear_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      call mass_init(ctx, grid, B, ierr)
      call linear_x(ctx, grid, B, ierr)
      call linear_y(ctx, grid, B, ierr)
      call linear_z(ctx, grid, B, ierr)

      ierr = 0
      return
   end subroutine linear_init

   subroutine mass_x(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), allocatable :: val
      PetscScalar, dimension(:), pointer :: FFTx

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         Mx => B%Mx, dx => grid%dx, xs => ctx%xs, xe => ctx%xe, Xk => B%Xk, &
         Nx => grid%Nx, xNc => grid%xNc)

      call VecDuplicate(B%rhsx, B%FFTx, ierr)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(val(2*deg+1))

      ! get min/max indices
      call check_periodic_index('knot', xs, xe, ctx%mx, deg, xNc, kmin, kmax)

      i = xs + deg
      val = 0.0
      do j=i-deg, i+deg
         do k = kmin, kmax
            val(j-xs+1) = val(j-xs+1) + quad_x_BB(grid, B, Xk(k), Xk(k+1), i, j, 0, 0)
         end do
      end do

      ! get array
      call DMDAVecGetArrayF90(ctx%DAx, B%FFTx, FFTx, ierr)
      CHKERRQ(ierr)

      if (xs-1 .eq. 0) then
         FFTx(0:B%deg) = val(B%deg+1:2*deg+1)
      end if

      if (xe-1 .eq. ctx%mx-1) then
         FFTx(xe-B%deg:xe-1) = val(1:B%deg)
      end if

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAx, B%FFTx, FFTx, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      deallocate(val)
      end associate
      ierr = 0
      return
   end subroutine mass_x

   subroutine mass_y(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B

      PetscInt :: i, j, k
      PetscInt :: jmin, jmax, kmin, kmax
      PetscScalar :: val
      PetscErrorCode :: ierr

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         My => B%My, dy => grid%dy, ys => ctx%ys, ye => ctx%ye, Yk => B%Yk, &
         Ny => grid%Ny)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAy, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, My, ierr)
      CHKERRQ(ierr)

      ! get min/max indices
      call check_dirichlet_index('knot', ys, ye, ctx%my, deg, Ny, kmin, kmax)
      call check_dirichlet_index('column', ys, ye, ctx%my, deg, Ny, jmin, jmax)

      ! initialize 1-D matrix
      do j=jmin, jmax
         do i=ys-1, ye-1
            if (i .eq. 0 .or. i .eq. ctx%my-1) then
               if (j .eq. i) then
               call MatSetValue(My, i, i, complex(1.0, 0.0), ADD_VALUES, ierr)
               CHKERRQ(ierr)
               end if
            else if (ABS(i-j) .le. deg) then
               do k=kmin, kmax
                  val = quad_y_BB(grid, B, Yk(k), Yk(k+1), i+1, j+1, 0, 0)
                  call MatSetValue(My, i, j, val, ADD_VALUES, ierr)
                  CHKERRQ(ierr)
               end do
            end if
         end do
      end do
      call MatAssemblyBegin(My, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(My, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine mass_y

   subroutine mass_z(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B

      PetscInt :: i, j, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), allocatable :: val
      PetscScalar, dimension(:), pointer :: FFTz
      PetscErrorCode :: ierr

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         Mz => B%Mz, dz => grid%dz, zs => ctx%zs, ze => ctx%ze, Zk => B%Zk, &
         Nz => grid%Nz, zNc => grid%zNc)

      call VecDuplicate(B%rhsz, B%FFTz, ierr)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(val(2*deg+1))

      ! get min/max indices
      call check_periodic_index('knot', zs, ze, ctx%mz, deg, zNc, kmin, kmax)

      i = zs + deg
      val = 0.0
      do j=i-deg, i+deg
         do k = kmin, kmax
            val(j-zs+1) = val(j-zs+1) + quad_z_BB(grid, B, Zk(k), Zk(k+1), i, j, 0, 0)
         end do
      end do

      ! get array
      call DMDAVecGetArrayF90(ctx%DAz, B%FFTz, FFTz, ierr)
      CHKERRQ(ierr)

      if (zs-1 .eq. 0) then
         FFTz(0:B%deg) = val(B%deg+1:2*deg+1)
      end if

      if (ze-1 .eq. ctx%mz-1) then
         FFTz(ze-B%deg:ze-1) = val(1:B%deg)
      end if

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAz, B%FFTz, FFTz, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      deallocate(val)
      end associate
      ierr = 0
      return
   end subroutine mass_z

   subroutine nonlinear_x(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m
      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         xs => ctx%xs, xe => ctx%xe, Xk => B%Xk, xNc => grid%xNc, Nx => grid%Nx)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 1
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBx(xs-1:xe-1, -deg:deg, -deg:deg))
      allocate(B%BBDx(xs-1:xe-1, -deg:deg, -deg:deg))

      B%BBBx = 0.0
      B%BBDx = 0.0

      do i=xs+deg, xe+deg
         do j=-deg, deg
            if (i+j<1 .or. i+j>Nx+deg) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Nx+deg) cycle
               do k=max(i-deg-1, 0), min(i+deg+1, Nx+deg-1)

                  B%BBBx(i-deg-1, j, l) = B%BBBx(i-deg-1, j, l) + &
                     quad_x_BBB(grid, B, Xk(k), Xk(k+1), i, i+j, i+l, 0, 0, 0)

                  B%BBDx(i-deg-1, j, l) = B%BBDx(i-deg-1, j, l) + &
                     quad_x_BBB(grid, B, Xk(k), Xk(k+1), i, i+j, i+l, 0, 0, 1)
               end do
            end do
         end do
      end do

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 3
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBBx(xs-1:xe-1, -deg:deg, -deg:deg, -deg:deg))
      allocate(B%BBBDx(xs-1:xe-1, -deg:deg, -deg:deg, -deg:deg))

      B%BBBBx = 0.0
      B%BBBDx = 0.0

      do i=xs+deg, xe+deg
         do j=-deg, deg
            if (i+j<1 .or. i+j>Nx+deg) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Nx+deg) cycle
               do m=-deg, deg
                  if (i+m<1 .or. i+m>Nx+deg) cycle
                  do k=max(i-deg-1, 0), min(i+deg+1, Nx+deg-1)

                     B%BBBBx(i-deg-1, j, l, m) = B%BBBBx(i-deg-1, j, l, m) + &
                        quad_x_BBBB(grid, B, Xk(k), Xk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 0)

                     B%BBBDx(i-deg-1, j, l, m) = B%BBBDx(i-deg-1, j, l, m) + &
                        quad_x_BBBB(grid, B, Xk(k), Xk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 1)
                  end do
               end do
            end do
         end do
      end do

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_x

   subroutine nonlinear_y(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m
      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         ys => ctx%ys, ye => ctx%ye, Yk => B%Yk, yNc => grid%yNc, Ny => grid%Ny)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 1
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBy(ys-1:ye-1, -deg:deg, -deg:deg))
      allocate(B%FBBBy(ys-1:ye-1, -deg:deg, -deg:deg))
      allocate(B%BBDy(ys-1:ye-1, -deg:deg, -deg:deg))

      B%BBBy = 0.0
      B%BBDy = 0.0

      do i=ys, ye
         do j=-deg, deg
            if (i+j<1 .or. i+j>Ny) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Ny) cycle
               do k=max(i-deg-1, 0), min(i+deg+1, yNc-1)

                  B%BBBy(i-1, j, l) = B%BBBy(i-1, j, l) + &
                     quad_y_BBB(grid, B, Yk(k), Yk(k+1), i, i+j, i+l, 0, 0, 0)

                  B%BBDy(i-1, j, l) = B%BBDy(i-1, j, l) + &
                     quad_y_BBB(grid, B, Yk(k), Yk(k+1), i, i+j, i+l, 0, 0, 1)
               end do
            end do
         end do
      end do

      B%FBBBy = B%BBBy

      if (ys-1 .eq. 0) then
         B%BBBy(ys-1, :, :) = 0.0
         B%BBBy(ys-1, 0, 0) = 1.0
      end if

      if (ye-1 .eq. ctx%my-1) then
         B%BBBy(ye-1, :, :) = 0.0
         B%BBBy(ye-1, 0, 0) = 1.0
      end if

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 3
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBBy(ys-1:ye-1, -deg:deg, -deg:deg, -deg:deg))
      allocate(B%BBBDy(ys-1:ye-1, -deg:deg, -deg:deg, -deg:deg))

      B%BBBBy = 0.0
      B%BBBDy = 0.0

      do i=ys, ye
         do j=-deg, deg
            if (i+j<1 .or. i+j>Ny) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Ny) cycle
               do m=-deg, deg
                  if (i+m<1 .or. i+m>Ny) cycle
                  do k=max(i-deg-1, 0), min(i+deg+1, yNc-1)

                     B%BBBBy(i-1, j, l, m) = B%BBBBy(i-1, j, l, m) + &
                        quad_y_BBBB(grid, B, Yk(k), Yk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 0)

                     B%BBBDy(i-1, j, l, m) = B%BBBDy(i-1, j, l, m) + &
                        quad_y_BBBB(grid, B, Yk(k), Yk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 1)
                  end do
               end do
            end do
         end do
      end do

      if (ys-1 .eq. 0) then
         B%BBBBy(ys-1, :, :, :) = 0.0
         B%BBBBy(ys-1, 0, 0, 0) = 1.0
      end if

      if (ye-1 .eq. ctx%my-1) then
         B%BBBBy(ye-1, :, :, :) = 0.0
         B%BBBBy(ye-1, 0, 0, 0) = 1.0
      end if

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_y

   subroutine nonlinear_z(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m
      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,  &
         zs => ctx%zs, ze => ctx%ze, Zk => B%Zk, zNc => grid%zNc, Nz => grid%Nz)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 1
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBz(zs-1:ze-1, -deg:deg, -deg:deg))
      allocate(B%BBDz(zs-1:ze-1, -deg:deg, -deg:deg))

      B%BBBz = 0.0
      B%BBDz = 0.0

      do i=zs+deg, ze+deg
         do j=-deg, deg
            if (i+j<1 .or. i+j>Nz+deg) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Nz+deg) cycle
               do k=max(i-deg-1, 0), min(i+deg+1, Nz+deg-1)

                  B%BBBz(i-deg-1, j, l) = B%BBBz(i-deg-1, j, l) + &
                     quad_z_BBB(grid, B, Zk(k), Zk(k+1), i, i+j, i+l, 0, 0, 0)

                  B%BBDz(i-deg-1, j, l) = B%BBDz(i-deg-1, j, l) + &
                     quad_z_BBB(grid, B, Zk(k), Zk(k+1), i, i+j, i+l, 0, 0, 1)
               end do
            end do
         end do
      end do

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 2*deg + 3
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      allocate(B%BBBBz(zs-1:ze-1, -deg:deg, -deg:deg, -deg:deg))
      allocate(B%BBBDz(zs-1:ze-1, -deg:deg, -deg:deg, -deg:deg))

      B%BBBBz = 0.0
      B%BBBDz = 0.0

      do i=zs+deg, ze+deg
         do j=-deg, deg
            if (i+j<1 .or. i+j>Nz+deg) cycle
            do l=-deg, deg
               if (i+l<1 .or. i+l>Nz+deg) cycle
               do m=-deg, deg
                  if (i+m<1 .or. i+m>Nz+deg) cycle
                  do k=max(i-deg-1, 0), min(i+deg+1, Nz+deg-1)

                     B%BBBBz(i-deg-1, j, l, m) = B%BBBBz(i-deg-1, j, l, m) + &
                        quad_z_BBBB(grid, B, Zk(k), Zk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 0)

                     B%BBBDz(i-deg-1, j, l, m) = B%BBBDz(i-deg-1, j, l, m) + &
                        quad_z_BBBB(grid, B, Zk(k), Zk(k+1), i, i+j, i+l, i+m, 0, 0, 0, 1)
                  end do
               end do
            end do
         end do
      end do

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_z

   subroutine linear_x(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B

      PetscInt :: i, j, k, l
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), allocatable :: val, val2, val3
      PetscErrorCode :: ierr

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,      &
         BBx => B%BBx, BDx => B%BDx, BLx => B%BLx, dx => grid%dx, xs => ctx%xs, &
         xe => ctx%xe, Xk => B%Xk, Nx => grid%Nx, xNc => grid%xNc)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAx, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAx, BBx, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAx, BDx, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAx, BLx, ierr)
      CHKERRQ(ierr)

      allocate(val(2*deg+1), val2(2*deg+1), val3(2*deg+1))

      ! get min/max indices
      call check_periodic_index('knot', xs, xe, ctx%mx, deg, xNc, kmin, kmax)

      i = xs + deg
      val = 0.0; val2 = 0.0; val3 = 0.0
      do j=i-deg, i+deg
         do k = kmin, kmax
            val(j-xs+1) = val(j-xs+1) + quad_x_BB(grid, B, Xk(k), Xk(k+1), i, j, 0, 0)
            val2(j-xs+1) = val2(j-xs+1) + quad_x_BB(grid, B, Xk(k), Xk(k+1), i, j, 0, 1)
            val3(j-xs+1) = val3(j-xs+1) + quad_x_BB(grid, B, Xk(k), Xk(k+1), i, j, 0, 2)
         end do
      end do

      do i=xs-1, xe-1
         l = 1
         do k=i-deg, i+deg
            if (k .lt. 0) then
               j = k + ctx%mx
            else if (k .gt. (ctx%mx-1)) then
               j = k - ctx%mx
            else
               j = k
            end if
            call MatSetValue(BBx, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            call MatSetValue(BDx, i, j, val2(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            call MatSetValue(BLx, i, j, val3(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(BBx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BBx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BDx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BDx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BLx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BLx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      deallocate(val, val2, val3)
      end associate
      ierr = 0
      return
   end subroutine linear_x

   subroutine linear_y(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B

      PetscInt :: i, j, k
      PetscInt :: jmin, jmax, kmin, kmax
      PetscScalar :: val
      PetscErrorCode :: ierr

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,      &
         BBy => B%BBy, BDy => B%BDy, BLy => B%BLy, dy => grid%dy, ys => ctx%ys, &
         ye => ctx%ye, Yk => B%Yk, Ny => grid%Ny)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAy, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, BBy, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, BDy, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, BLy, ierr)
      CHKERRQ(ierr)

      ! get min/max indices
      call check_dirichlet_index('knot', ys, ye, ctx%my, deg, Ny, kmin, kmax)
      call check_dirichlet_index('column', ys, ye, ctx%my, deg, Ny, jmin, jmax)

      ! initialize 1-D matrix
      do j=jmin, jmax
         do i=ys-1, ye-1
            if (ABS(i-j) .le. deg) then
               do k=kmin, kmax
                  val = quad_y_BB(grid, B, Yk(k), Yk(k+1), i+1, j+1, 0, 0)
                  call MatSetValue(BBy, i, j, val, ADD_VALUES, ierr)
                  CHKERRQ(ierr)

                  val = quad_y_BB(grid, B, Yk(k), Yk(k+1), i+1, j+1, 0, 1)
                  call MatSetValue(BDy, i, j, val, ADD_VALUES, ierr)
                  CHKERRQ(ierr)

                  val = quad_y_BB(grid, B, Yk(k), Yk(k+1), i+1, j+1, 0, 2)
                  call MatSetValue(BLy, i, j, val, ADD_VALUES, ierr)
                  CHKERRQ(ierr)
               end do
            end if
         end do
      end do
      call MatAssemblyBegin(BBy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BBy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BDy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BDy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BLy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BLy, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine linear_y

  subroutine linear_z(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B

      PetscInt :: i, j, k, l
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), allocatable :: val, val2, val3
      PetscErrorCode :: ierr

      associate (xi => B%xi, wts => B%wts, nodes => B%nodes, deg => B%deg,      &
         BBz => B%BBz, BDz => B%BDz, BLz => B%BLz, dz => grid%dz, zs => ctx%zs, &
         ze => ctx%ze, Zk => B%Zk, Nz => grid%Nz, zNc => grid%zNc)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAz, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAz, BBz, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAz, BDz, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAz, BLz, ierr)
      CHKERRQ(ierr)

      allocate(val(2*deg+1), val2(2*deg+1), val3(2*deg+1))

      ! get min/max indices
      call check_periodic_index('knot', zs, ze, ctx%mz, deg, zNc, kmin, kmax)

      i = zs + deg
      val = 0.0; val2 = 0.0; val3 = 0.0
      do j=i-deg, i+deg
         do k = kmin, kmax
            val(j-zs+1) = val(j-zs+1) + quad_z_BB(grid, B, Zk(k), Zk(k+1), i, j, 0, 0)
            val2(j-zs+1) = val2(j-zs+1) + quad_z_BB(grid, B, Zk(k), Zk(k+1), i, j, 0, 1)
            val3(j-zs+1) = val3(j-zs+1) + quad_z_BB(grid, B, Zk(k), Zk(k+1), i, j, 0, 2)
         end do
      end do

      do i=zs-1, ze-1
         l = 1
         do k=i-deg, i+deg
            if (k .lt. 0) then
               j = k + ctx%mz
            else if (k .gt. (ctx%mz-1)) then
               j = k - ctx%mz
            else
               j = k
            end if
            call MatSetValue(BBz, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            call MatSetValue(BDz, i, j, val2(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            call MatSetValue(BLz, i, j, val3(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(BBz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BBz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BDz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BDz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyBegin(BLz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(BLz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      deallocate(val, val2, val3)
      end associate
      ierr = 0
      return
   end subroutine linear_z

end module operator

! vim: set ts=3 sw=3 sts=3:
