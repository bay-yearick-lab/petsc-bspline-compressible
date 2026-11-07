! Max-norm error against manufactured solutions and exact operator evaluations.
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>

module error
   use types
   use bspline
   use linalg
   implicit none

   private

   public Linf_eval_x
   public Linf_eval_y
   public Linf_eval_z
   public Linf_eval_xyz
   public compute_sol_convection
   public compute_sol_viscous
   public compute_sol_pressure
   public compute_sol_momentum
   public compute_sol_energy_viscous
   public compute_sol_energy_convection
   public compute_sol_energy
contains

   subroutine check_index(type, s, e, m, d, N, minsize, maxsize)
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
            maxsize = N - d - 1
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
   end subroutine check_index

   subroutine Linf_eval_xyz(ctx, grid, B, f, alpha, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure(func3D), pointer :: f
      Vec :: alpha
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, imin, imax
      PetscReal :: dy, Linf
      PetscScalar, dimension(:), allocatable :: val
      PetscScalar, dimension(:, :, :), pointer :: sol

      associate (deg => B%deg, Bx => B%Bx, xs => ctx%xs, xe => ctx%xe, xNc => grid%xNc,      &
         Xk => B%Xk, Nx => grid%Nx, By => B%By, ys => ctx%ys, ye => ctx%ye, Yk => B%Yk,      &
         Ny => grid%Ny, Bz => B%Bz, zs => ctx%zs, ze => ctx%ze, zNc => grid%zNc, Zk => B%Zk, &
         Nz => grid%Nz)

      ! setup x matrix
      call DMSetMatType(ctx%DAx, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAx, Bx, ierr)
      CHKERRQ(ierr)

      allocate(val(deg+1))

      i = xs + deg
      val = 0.0
      do j=i-deg, i
         val(j-xs+1) = xBJK(grid, B, i+1, deg, 0.5 * SUM(Xk(j:j+1)))
      end do

      ! initialize x matrix
      do i=xs-1, xe-1
         l = 1
         do k=i-deg, i
            if (k .lt. 0) then
               j = k + ctx%mx
            else if (k .gt. (ctx%mx-1)) then
               j = k - ctx%mx
            else
               j = k
            end if
            call MatSetValue(Bx, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(Bx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(Bx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! setup y matrix
      call DMSetMatType(ctx%DAy, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, By, ierr)
      CHKERRQ(ierr)

      dy = (grid%Ly - grid%dy) / Ny

      call check_index('column', ys, ye, ctx%my, deg, Ny, imin, imax)

      ! initialize y matrix
      do j=ys-1, ye-1
         do i=imin, imax
            if (ABS(i-j) .le. deg) then
               val(1) = yBJK(grid, B, i+1, deg, (j+1)*dy)
               call MatSetValue(By, j, i, val(1), INSERT_VALUES, ierr)
               CHKERRQ(ierr)
            end if
         end do
      end do
      call MatAssemblyBegin(By, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(By, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! setup z matrix
      call DMSetMatType(ctx%DAz, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAz, Bz, ierr)
      CHKERRQ(ierr)

      i = zs + deg
      val = 0.0
      do j=i-deg, i
         val(j-zs+1) = zBJK(grid, B, i+1, deg, 0.5 * SUM(Zk(j:j+1)))
      end do

      ! initialize 1-D matrix
      do i=zs-1, ze-1
         l = 1
         do k=i-deg, i
            if (k .lt. 0) then
               j = k + ctx%mz
            else if (k .gt. (ctx%mz-1)) then
               j = k - ctx%mz
            else
               j = k
            end if
            call MatSetValue(Bz, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(Bz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(Bz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      ! numerical solution
      call matvec(ctx, B, Bx, By, Bz, alpha, B%temp, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, B%sol, sol, ierr)
      CHKERRQ(ierr)

      ! loop over B_splines
      do k=zs-1, ze-1
         do j=ys-1, ye-1
            do i=xs-1, xe-1
               sol(i, j, k) = f(0.5 * SUM(Xk(i:i+1)), &
                  (j+1) * dy, 0.5 * SUM(Zk(k:k+1)), 0.0)
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, B%sol, sol, ierr)
      CHKERRQ(ierr)

      ! error
      call VecAXPY(B%temp, -1.0 + 0*PETSC_i, B%sol, ierr)
      CHKERRQ(ierr)
      call VecNorm(B%temp, NORM_INFINITY, Linf, ierr)
      CHKERRQ(ierr)

      print*, Linf

      ! Destroy
      call MatDestroy(Bx, ierr)
      CHKERRQ(ierr)
      call MatDestroy(By, ierr)
      CHKERRQ(ierr)
      call MatDestroy(Bz, ierr)
      CHKERRQ(ierr)

      deallocate(val)
      end associate
      ierr = 0
      return
   end subroutine Linf_eval_xyz

   subroutine Linf_eval_x(ctx, grid, B, ptr, alpha, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: alpha
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l
      PetscReal :: Linf
      PetscScalar, dimension(:), allocatable :: val
      PetscScalar, dimension(:), pointer :: vvec

      associate (deg => B%deg, Bx => B%Bx, xs => ctx%xs, xe => ctx%xe, xNc => grid%xNc, &
         Xk => B%Xk, Nx => grid%Nx, tempx => B%tempx, solx => B%solx)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAx, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAx, Bx, ierr)
      CHKERRQ(ierr)

      allocate(val(deg+1))

      i = xs + deg
      val = 0.0
      do j=i-deg, i
         val(j-xs+1) = xBJK(grid, B, i+1, deg, 0.5 * SUM(Xk(j:j+1)))
      end do

      ! initialize 1-D matrix
      do i=xs-1, xe-1
         l = 1
         do k=i-deg, i
            if (k .lt. 0) then
               j = k + ctx%mx
            else if (k .gt. (ctx%mx-1)) then
               j = k - ctx%mx
            else
               j = k
            end if
            call MatSetValue(Bx, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(Bx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(Bx, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      call MatMult(Bx, alpha, tempx, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAx, solx, vvec, ierr)
      CHKERRQ(ierr)

      ! loop over B_splines
      do i=xs-1, xe-1
         vvec(i) = ptr(0.5 * SUM(Xk(i:i+1)), 0.0)
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAx, solx, vvec, ierr)
      CHKERRQ(ierr)

      ! error
      call VecAXPY(tempx, -1.0 + 0*PETSC_i, solx, ierr)
      CHKERRQ(ierr)
      call VecNorm(tempx, NORM_INFINITY, Linf, ierr)
      CHKERRQ(ierr)

      print*, Linf

      deallocate(val)
      end associate
      ierr = 0
      return
   end subroutine Linf_eval_x

   subroutine Linf_eval_y(ctx, grid, B, ptr, alpha, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: alpha
      PetscErrorCode :: ierr

      PetscInt :: i, j, imin, imax
      PetscReal :: dy, Linf
      PetscScalar :: val
      PetscScalar, dimension(:), pointer :: vvec

      associate (deg => B%deg, By => B%By, ys => ctx%ys, ye => ctx%ye, &
         Yk => B%Yk, Ny => grid%Ny, tempy => B%tempy, soly => B%soly)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAy, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAy, By, ierr)
      CHKERRQ(ierr)

      call check_index('column', ys, ye, ctx%my, deg, Ny, imin, imax)

      dy = (grid%Ly - grid%dy) / Ny

      ! initialize 1-D matrix
      do j=ys-1, ye-1
         do i=imin, imax
            if (ABS(i-j) .le. deg) then
               val = yBJK(grid, B, i+1, deg, (j+1) * dy)
               call MatSetValue(By, j, i, val, INSERT_VALUES, ierr)
               CHKERRQ(ierr)
            end if
         end do
      end do
      call MatAssemblyBegin(By, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(By, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      call MatMult(By, alpha, tempy, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAy, soly, vvec, ierr)
      CHKERRQ(ierr)

      ! loop over B_splines
      do i=ys-1, ye-1
         vvec(i) = ptr((i+1) * dy, 0.0)
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAy, soly, vvec, ierr)
      CHKERRQ(ierr)

      ! error
      call VecAXPY(tempy, -1.0 + 0*PETSC_i, soly, ierr)
      CHKERRQ(ierr)
      call VecNorm(tempy, NORM_INFINITY, Linf, ierr)
      CHKERRQ(ierr)

      print*, Linf

      end associate
      ierr = 0
      return
   end subroutine Linf_eval_y

   subroutine Linf_eval_z(ctx, grid, B, ptr, alpha, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: alpha
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l
      PetscReal :: Linf
      PetscScalar, dimension(:), allocatable :: val
      PetscScalar, dimension(:), pointer :: vvec

      associate (deg => B%deg, Bz => B%Bz, zs => ctx%zs, ze => ctx%ze, zNc => grid%zNc, &
         Zk => B%Zk, Nz => grid%Nz, tempz => B%tempz, solz => B%solz)

      ! setup 1-D matrix
      call DMSetMatType(ctx%DAz, MATAIJ, ierr)
      CHKERRQ(ierr)
      call DMCreateMatrix(ctx%DAz, Bz, ierr)
      CHKERRQ(ierr)

      allocate(val(deg+1))

      i = zs + deg
      val = 0.0
      do j=i-deg, i
         val(j-zs+1) = zBJK(grid, B, i+1, deg, 0.5 * SUM(Zk(j:j+1)))
      end do

      ! initialize 1-D matrix
      do i=zs-1, ze-1
         l = 1
         do k=i-deg, i
            if (k .lt. 0) then
               j = k + ctx%mz
            else if (k .gt. (ctx%mz-1)) then
               j = k - ctx%mz
            else
               j = k
            end if
            call MatSetValue(Bz, i, j, val(l), INSERT_VALUES, ierr)
            CHKERRQ(ierr)
            l = l + 1
         end do
      end do
      call MatAssemblyBegin(Bz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(Bz, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)

      call MatMult(Bz, alpha, tempz, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAz, solz, vvec, ierr)
      CHKERRQ(ierr)

      ! loop over B_splines
      do i=zs-1, ze-1
         vvec(i) = ptr(0.5 * SUM(Zk(i:i+1)), 0.0)
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAz, solz, vvec, ierr)
      CHKERRQ(ierr)

      ! error
      call VecAXPY(tempz, -1.0 + 0*PETSC_i, solz, ierr)
      CHKERRQ(ierr)
      call VecNorm(tempz, NORM_INFINITY, Linf, ierr)
      CHKERRQ(ierr)

      print*, Linf

      deallocate(val)
      end associate
      ierr = 0
      return
   end subroutine Linf_eval_z

   subroutine compute_sol_convection(ctx, grid, B, i, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1, s2, s3, s4, s5, s6

      if (i .eq. 1) then
         ! compute solution
         s1(1)%ptr => fcos4sin4; s1(2)%ptr => fsin42; s1(3)%ptr => fcos42
         s2(1)%ptr => fcos2cos4sin4; s2(2)%ptr => fcos4sin42; s2(3)%ptr => fcos2cos42
         s3(1)%ptr => fcos2cos4; s3(2)%ptr => fsin2sin4; s3(3)%ptr => fcos2sin4
         s4(1)%ptr => fcos22cos4; s4(2)%ptr => fsin2sin4cos4; s4(3)%ptr => fcos22sin4
         s5(1)%ptr => fcos4sin4; s5(2)%ptr => fcos42; s5(3)%ptr => fsin8
         s6(1)%ptr => fcos2cos4sin4; s6(2)%ptr => fcos43; s6(3)%ptr => fcos2sin8

         call VecZeroEntries(vector, ierr)

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-2.0, 0.0), temp, ierr)

      else if (i .eq. 2) then
         ! compute solution
         s1(1)%ptr => fcos2sin4; s1(2)%ptr => fsin2cos4; s1(3)%ptr => fcos2cos4
         s2(1)%ptr => fcos22sin4; s2(2)%ptr => fsin2cos42; s2(3)%ptr => fcos22cos4
         s3(1)%ptr => fsin42; s3(2)%ptr => fcos4sin4; s3(3)%ptr => fsin42
         s4(1)%ptr => fcos2sin42; s4(2)%ptr => fsin4cos42; s4(3)%ptr => fcos2sin42
         s5(1)%ptr => fcos42; s5(2)%ptr => fcos4sin4; s5(3)%ptr => fsin8
         s6(1)%ptr => fcos2cos42; s6(2)%ptr => fsin4cos42; s6(3)%ptr => fcos2sin8

         call VecZeroEntries(vector, ierr)

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-2.0, 0.0), temp, ierr)

      else if (i .eq. 3) then
         ! compute solution
         s1(1)%ptr => fsin2cos4; s1(2)%ptr => fsin2sin4; s1(3)%ptr => fcos2cos4
         s2(1)%ptr => fcos2sin2cos4; s2(2)%ptr => fsin2sin4cos4; s2(3)%ptr => fcos22cos4
         s3(1)%ptr => fcos22; s3(2)%ptr => fsin22; s3(3)%ptr => fcos2sin2
         s4(1)%ptr => fcos23; s4(2)%ptr => fsin22cos4; s4(3)%ptr => fcos22sin2
         s5(1)%ptr => fcos2sin4; s5(2)%ptr => fcos2cos4; s5(3)%ptr => fcos2sin4
         s6(1)%ptr => fcos22sin4; s6(2)%ptr => fcos2cos42; s6(3)%ptr => fcos22sin4

         call VecZeroEntries(vector, ierr)

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(2.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(2.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-2.0, 0.0), temp, ierr)
      end if

      ierr = 0
      return
   end subroutine compute_sol_convection

   subroutine compute_sol_viscous(ctx, grid, B, i, mu, Re, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, Re
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1, s2, s3

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      if (i .eq. 1) then
         ! compute solution
         s1(1)%ptr => fcos4; s1(2)%ptr => fsin4; s1(3)%ptr => fcos4
         s2(1)%ptr => fcos4; s2(2)%ptr => fsin4; s2(3)%ptr => fsin4
         s3(1)%ptr => fsin2; s3(2)%ptr => fsin2; s3(3)%ptr => fsin2

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-160.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(4.0/3.0, 0.0), temp, ierr)

      else if (i .eq. 2) then
         ! compute solution
         s1(1)%ptr => fsin4; s1(2)%ptr => fcos4; s1(3)%ptr => fsin4
         s2(1)%ptr => fsin4; s2(2)%ptr => fcos4; s2(3)%ptr => fcos4
         s3(1)%ptr => fcos2; s3(2)%ptr => fcos2; s3(3)%ptr => fsin2

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-160.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0/3.0, 0.0), temp, ierr)

      else if (i .eq. 3) then
         ! compute solution
         s1(1)%ptr => fsin4; s1(2)%ptr => fsin4; s1(3)%ptr => fsin4
         s2(1)%ptr => fsin4; s2(2)%ptr => fsin4; s2(3)%ptr => fcos4
         s3(1)%ptr => fcos2; s3(2)%ptr => fsin2; s3(3)%ptr => fcos2

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(16.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-40.0/3.0, 0.0), temp, ierr)
      end if

      call VecScale(vector, complex(mu/Re, 0.0), ierr)

      ierr = 0
      return
   end subroutine compute_sol_viscous

   subroutine compute_sol_pressure(ctx, grid, B, i, gamma, M, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: gamma, M
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1, s2, s3, s4

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      if (i .eq. 1) then
         ! compute solution
         s1(1)%ptr => fsin2; s1(2)%ptr => fcos4; s1(3)%ptr => fcos2
         s2(1)%ptr => fsin2; s2(2)%ptr => fcos2cos4; s2(3)%ptr => fcos2cos4
         s3(1)%ptr => fsin4; s3(2)%ptr => fcos2; s3(3)%ptr => fcos4
         s4(1)%ptr => fsin6; s4(2)%ptr => fcos2cos4; s4(3)%ptr => fcos2cos4

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-1.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-3.0, 0.0), temp, ierr)

      else if (i .eq. 2) then
         ! compute solution
         s1(1)%ptr => fcos4; s1(2)%ptr => fsin2; s1(3)%ptr => fcos4
         s2(1)%ptr => fcos2cos4; s2(2)%ptr => fsin2; s2(3)%ptr => fcos2cos4
         s3(1)%ptr => fcos2; s3(2)%ptr => fsin4; s3(3)%ptr => fcos2
         s4(1)%ptr => fcos2cos4; s4(2)%ptr => fsin6; s4(3)%ptr => fcos2cos4

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-1.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-3.0, 0.0), temp, ierr)

      else if (i .eq. 3) then
         ! compute solution
         s1(1)%ptr => fcos2; s1(2)%ptr => fcos4; s1(3)%ptr => fsin2
         s2(1)%ptr => fcos2cos4; s2(2)%ptr => fcos2cos4; s2(3)%ptr => fsin2cos4
         s3(1)%ptr => fcos4; s3(2)%ptr => fcos2; s3(3)%ptr => fsin4
         s4(1)%ptr => fcos2cos4; s4(2)%ptr => fcos2cos4; s4(3)%ptr => fcos2sin4

         call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-2.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

         call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
         CHKERRQ(ierr)
         call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)
      end if

      call VecScale(vector, complex(-1.0/(gamma * M**2), 0.0), ierr)

      ierr = 0
      return
   end subroutine compute_sol_pressure

   subroutine compute_sol_momentum(ctx, grid, B, i, &
         mu, gamma, M, Re, temp1, temp2, vector, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, gamma, M, Re
      Vec :: temp1, temp2, vector
      PetscErrorCode, intent(out) :: ierr

      call VecZeroEntries(temp1, ierr)
      call VecZeroEntries(temp2, ierr)
      call VecZeroEntries(vector, ierr)

      call compute_sol_convection(ctx, grid, B, i, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      call compute_sol_pressure(ctx, grid, B, i, gamma, M, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      call compute_sol_viscous(ctx, grid, B, i, mu, Re, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      ierr = 0
      return
   end subroutine compute_sol_momentum

   subroutine compute_sol_energy_viscous(ctx, grid, B, mu, gamma, M, Re, &
         temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: gamma, mu, M, Re
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      type (pointer3D), dimension(3) :: s1, s2, s3, s4
      type (pointer3D), dimension(3) :: s5, s6, s7, s8, s9

      call VecZeroEntries(vector, ierr)

      do i=1, 3
         if (i .eq. 1) then
            ! compute solution
            s1(1)%ptr => fcos42; s1(2)%ptr => fone; s1(3)%ptr => fone
            s2(1)%ptr => fcos42; s2(2)%ptr => fcos8; s2(3)%ptr => fcos8
            s3(1)%ptr => fsin42; s3(2)%ptr => fsin42; s3(3)%ptr => fcos42
            s4(1)%ptr => fcos2sin4; s4(2)%ptr => fsin2sin4; s4(3)%ptr => fsin2cos4
            s5(1)%ptr => fsin2cos4; s5(2)%ptr => fsin2sin4; s5(3)%ptr => fcos2sin4
            s6(1)%ptr => fsin42; s6(2)%ptr => fsin42; s6(3)%ptr => fcos4sin4
            s7(1)%ptr => fcos42; s7(2)%ptr => fcos42; s7(3)%ptr => fsin8

            call VecZeroEntries(temp, ierr)

            call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(64.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-32.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s7, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

         else if (i .eq. 2) then
            ! compute solution
            s1(1)%ptr => fcos42; s1(2)%ptr => fcos42; s1(3)%ptr => fone
            s2(1)%ptr => fcos42; s2(2)%ptr => fcos42; s2(3)%ptr => fcos8
            s3(1)%ptr => fcos2sin4; s3(2)%ptr => fcos2cos4; s3(3)%ptr => fcos2cos4
            s4(1)%ptr => fsin42; s4(2)%ptr => fcos42; s4(3)%ptr => fcos42
            s5(1)%ptr => fsin42; s5(2)%ptr => fsin42; s5(3)%ptr => fone
            s6(1)%ptr => fsin42; s6(2)%ptr => fsin42; s6(3)%ptr => fcos8
            s7(1)%ptr => fcos2sin4; s7(2)%ptr => fsin2sin4; s7(3)%ptr => fsin2sin4
            s8(1)%ptr => fcos42; s8(2)%ptr => fcos42; s8(3)%ptr => fsin8
            s9(1)%ptr => fsin42; s9(2)%ptr => fsin42; s9(3)%ptr => fsin8

            call VecZeroEntries(temp, ierr)

            call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(16.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(32.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-32.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s7, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s8, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s9, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

         else if (i .eq. 3) then
            ! compute solution
            s1(1)%ptr => fcos22; s1(2)%ptr => fcos22; s1(3)%ptr => fcos22
            s2(1)%ptr => fcos2sin4; s2(2)%ptr => fcos2cos4; s2(3)%ptr => fcos2cos4
            s3(1)%ptr => fsin22; s3(2)%ptr => fsin22; s3(3)%ptr => fcos22
            s4(1)%ptr => fcos2sin4; s4(2)%ptr => fcos2sin22; s4(3)%ptr => fcos2
            s5(1)%ptr => fcos2sin4; s5(2)%ptr => fcos2sin22; s5(3)%ptr => fcos6
            s6(1)%ptr => fcos2sin4; s6(2)%ptr => fcos2sin22; s6(3)%ptr => fsin2
            s7(1)%ptr => fcos22; s7(2)%ptr => fsin22; s7(3)%ptr => fsin22
            s8(1)%ptr => fsin2cos4; s8(2)%ptr => fsin2sin4; s8(3)%ptr => fcos2sin4
            s9(1)%ptr => fcos2sin4; s9(2)%ptr => fcos2sin22; s9(3)%ptr => fsin6

            call VecZeroEntries(temp, ierr)

            call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s7, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(16.0/3.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s8, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

            call rhs_init_separable(ctx, grid, B, s9, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
            CHKERRQ(ierr)
            call VecAXPY(vector, complex(-16.0/3.0, 0.0), temp, ierr)
         end if
      end do

      call VecScale(vector, complex((gamma* (gamma - 1.0) * M**2 * mu)/Re, 0.0), ierr)

      ierr = 0
      return
   end subroutine compute_sol_energy_viscous

   subroutine compute_sol_energy_laplacian(ctx, grid, B, &
         gamma, mu, Re, Pr, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: gamma, mu, Re, Pr
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1

      ! compute solution
      s1(1)%ptr => fcos4; s1(2)%ptr => fcos2; s1(3)%ptr => fcos4

      call VecZeroEntries(vector, ierr)

      call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-36.0, 0.0), temp, ierr)

      call VecScale(vector, complex((gamma * mu) / (Re * Pr), 0.0), ierr)

      ierr = 0
      return
   end subroutine compute_sol_energy_laplacian

   subroutine compute_sol_energy_divergence(ctx, grid, B, gamma, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: gamma
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1, s2, s3, s4, s5, s6

      ! compute solution
      s1(1)%ptr => fsin4; s1(2)%ptr => fsin4; s1(3)%ptr => fcos4
      s2(1)%ptr => fcos4sin4; s2(2)%ptr => fcos2sin4; s2(3)%ptr => fcos42
      s3(1)%ptr => fcos2; s3(2)%ptr => fsin2; s3(3)%ptr => fsin2
      s4(1)%ptr => fcos2cos4; s4(2)%ptr => fcos2sin2; s4(3)%ptr => fsin2cos4
      s5(1)%ptr => fsin4; s5(2)%ptr => fsin4; s5(3)%ptr => fsin4
      s6(1)%ptr => fcos4sin4; s6(2)%ptr => fcos2sin4; s6(3)%ptr => fcos4sin4

      call VecZeroEntries(vector, ierr)

      call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-2.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-8.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(-4.0, 0.0), temp, ierr)

      call VecScale(vector, complex(-(gamma - 1.0), 0.0), ierr)

      ierr = 0
      return
   end subroutine compute_sol_energy_divergence

   subroutine compute_sol_energy_convection(ctx, grid, B, temp, vector, ierr)
      use mms
      use rhs
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: temp, vector
      PetscErrorCode, intent(out) :: ierr

      type (pointer3D), dimension(3) :: s1, s2, s3, s4, s5, s6

      ! compute solution
      s1(1)%ptr => fcos4sin4; s1(2)%ptr => fcos2sin4; s1(3)%ptr => fcos42
      s2(1)%ptr => fcos2cos4sin4; s2(2)%ptr => fcos2cos4sin4; s2(3)%ptr => fcos2cos42
      s3(1)%ptr => fcos2cos4; s3(2)%ptr => fsin4; s3(3)%ptr => fcos22sin2
      s4(1)%ptr => fcos22cos4; s4(2)%ptr => fcos4sin4; s4(3)%ptr => fsin2cos23
      s5(1)%ptr => fcos4sin4; s5(2)%ptr => fsin2cos4; s5(3)%ptr => fcos4sin4
      s6(1)%ptr => fcos2cos4sin4; s6(2)%ptr => fsin2cos42; s6(3)%ptr => fcos2cos4sin4

      call VecZeroEntries(vector, ierr)

      call rhs_init_separable(ctx, grid, B, s1, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s2, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s3, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(8.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s4, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s5, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(4.0, 0.0), temp, ierr)

      call rhs_init_separable(ctx, grid, B, s6, temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(vector, complex(2.0, 0.0), temp, ierr)

      ierr = 0
      return
   end subroutine compute_sol_energy_convection

   subroutine compute_sol_energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         temp1, temp2, vector, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: mu, gamma, M, Re, Pr
      Vec :: temp1, temp2, vector
      PetscErrorCode, intent(out) :: ierr

      call VecZeroEntries(temp1, ierr)
      call VecZeroEntries(temp2, ierr)
      call VecZeroEntries(vector, ierr)

      call compute_sol_energy_convection(ctx, grid, B, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      call compute_sol_energy_divergence(ctx, grid, B, gamma, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      call compute_sol_energy_viscous(ctx, grid, B, mu, gamma, M, Re, &
         temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      call compute_sol_energy_laplacian(ctx, grid, B, &
         gamma, mu, Re, Pr, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      ierr = 0
      return
   end subroutine compute_sol_energy

end module error

! vim: set ts=3 sw=3 sts=3:
