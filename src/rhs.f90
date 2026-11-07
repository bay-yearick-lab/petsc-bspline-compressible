! Galerkin right-hand sides (load vectors) and Navier-Stokes residual terms on the DMDA.
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

module rhs
   use types
   use bspline
   use quad
   use mms
   use gausstable
   use linalg
   implicit none

   private

   public :: rhs_x
   public :: rhs_y
   public :: rhs_random_y
   public :: rhs_z
   public :: rhs_xyz
   public :: rhs_init_separable
   public :: rhs_random_init_separable
   public :: rhs_reciprocal_x
   public :: rhs_reciprocal_xyz
   public :: rhs_mass
   public :: rhs_momentum
   public :: rhs_energy
contains

   subroutine rhs_energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         rho, u, v, w, T, invrho, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp1, temp2, temp3, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: mu, gamma, M, Re, Pr
      Vec :: rho, u, v, w, T, invrho
      Vec :: rho_u, rho_v, rho_w
      Vec :: rho_uT, rho_vT, rho_wT
      Vec :: temp1, temp2, temp3, vector
      PetscErrorCode :: ierr

      PetscInt :: i

      call VecZeroEntries(temp1, ierr)
      call VecZeroEntries(temp2, ierr)
      call VecZeroEntries(temp3, ierr)
      call VecZeroEntries(vector, ierr)

      call convection(ctx, grid, B, u, v, w, T, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      do i=1, 3
         call energy_viscous(ctx, grid, B, i, mu, gamma, M, Re, &
            u, v, w, temp1, temp2, temp3, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp3, ierr)
      end do

      call energy_laplacian(ctx, grid, B, gamma, &
         mu, Re, Pr, T, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      ! divide by rho
      call ksp_solve(ctx, grid, B, vector, temp1, 'full', ierr)
      CHKERRQ(ierr)

      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         invrho, temp1, vector, B%local1, B%local2, ierr)

      call energy_divergence(ctx, grid, B, gamma, &
         u, v, w, T, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      ierr = 0
      return
   end subroutine rhs_energy

   subroutine rhs_momentum(ctx, grid, B, i, mu, gamma, M, Re, &
         rho, u, v, w, T, invrho, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp1, temp2, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, gamma, M, Re
      Vec :: rho, u, v, w, T, invrho
      Vec :: rho_u, rho_v, rho_w
      Vec :: rho_uT, rho_vT, rho_wT
      Vec :: temp1, temp2, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(temp1, ierr)
      call VecZeroEntries(temp2, ierr)
      call VecZeroEntries(vector, ierr)

      if (i .eq. 1) then
         call convection(ctx, grid, B, u, v, w, u, rho_u, rho_v, rho_w, &
            rho_uT, rho_vT, rho_wT, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call pressure(ctx, grid, B, i, gamma, M, rho, T, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call viscous(ctx, grid, B, i, mu, Re, u, v, w, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      else if (i .eq. 2) then
         call convection(ctx, grid, B, u, v, w, v, rho_u, rho_v, rho_w, &
            rho_uT, rho_vT, rho_wT, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call pressure(ctx, grid, B, i, gamma, M, rho, T, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call viscous(ctx, grid, B, i, mu, Re, u, v, w, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      else if (i .eq. 3) then
         call convection(ctx, grid, B, u, v, w, w, rho_u, rho_v, rho_w, &
            rho_uT, rho_vT, rho_wT, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call pressure(ctx, grid, B, i, gamma, M, rho, T, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         call viscous(ctx, grid, B, i, mu, Re, u, v, w, temp1, temp2, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)
      end if

      ! divide by rho
      call ksp_solve(ctx, grid, B, vector, temp1, 'full', ierr)
      CHKERRQ(ierr)

      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         invrho, temp1, vector, B%local1, B%local2, ierr)

      ierr = 0
      return
   end subroutine rhs_momentum

   ! All RHS routines takes alpha as input and outputs non-alpha
   !
   subroutine rhs_mass(ctx, grid, B, rho_u, rho_v, rho_w, vector, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: rho_u, rho_v, rho_w
      Vec :: vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)

      ! div(rho * velocity)
      call divergence(ctx, grid, B, rho_u, rho_v, rho_w, vector, ierr)
      CHKERRQ(ierr)

      ! scale by -1.0
      call VecScale(vector, complex(-1.0, 0.0), ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine rhs_mass

   subroutine check_index(s, e, m, d, N, minsize, maxsize)
      PetscInt :: s, e, m, d, N, minsize, maxsize

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
   end subroutine check_index

   subroutine rhs_reciprocal_xyz(ctx, grid, B, vv, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscScalar, dimension(:, :, :), pointer :: vv
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscInt :: kx, ky, kz
      PetscScalar, dimension(:, :, :), pointer :: vvec

      ! initialize gauss points and weights (nodes give exact quadrature)
      B%nodes = B%deg
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, vec, vvec, ierr)
      CHKERRQ(ierr)

      vvec = 0.0
      do k=ctx%zs-1, ctx%ze-1
         do kz=max(k-B%deg, 0), min(k+B%deg, grid%Nz+B%deg-1)
            do j=ctx%ys-1, ctx%ye-1
               do ky=max(j-B%deg, 0), min(j+B%deg, grid%yNc-1)
                  do i=ctx%xs-1, ctx%xe-1
                     do kx=max(i-B%deg, 0), min(i+B%deg, grid%Nx+B%deg-1)

                        vvec(i, j, k) = vvec(i, j, k) + &
                           quad_reciprocal_xyz_B(ctx, grid, B, &
                                                B%Xk(kx), B%Xk(kx+1), &
                                                B%Yk(ky), B%Yk(ky+1), &
                                                B%Zk(kz), B%Zk(kz+1), &
                              i+1+B%deg, j+1, k+1+B%deg, 0, 0, 0, vv, ierr)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine rhs_reciprocal_xyz

   subroutine rhs_xyz(ctx, grid, B, f, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure(func3D), pointer :: f
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscInt :: kx, ky, kz
      PetscScalar, dimension(:, :, :), pointer :: vv

      ! initialize gauss points and weights (nodes give exact quadrature)
      B%nodes = B%deg + 2
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, vec, vv, ierr)
      CHKERRQ(ierr)

      vv = 0.0
      do k=ctx%zs-1, ctx%ze-1
         do kz=max(k-B%deg, 0), min(k+B%deg, grid%Nz+B%deg-1)
            do j=ctx%ys-1, ctx%ye-1
               do ky=max(j-B%deg, 0), min(j+B%deg, grid%yNc-1)
                  do i=ctx%xs-1, ctx%xe-1
                     do kx=max(i-B%deg, 0), min(i+B%deg, grid%Nx+B%deg-1)

                        vv(i, j, k) = vv(i, j, k) + &
                           quad_xyz_B(grid, B, B%Xk(kx), B%Xk(kx+1), &
                                               B%Yk(ky), B%Yk(ky+1), &
                                               B%Zk(kz), B%Zk(kz+1), &
                              i+1+B%deg, j+1, k+1+B%deg, 0, 0, 0, f, 0.0)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, vec, vv, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine rhs_xyz

   subroutine rhs_random_init_separable(ctx, grid, B, f, vec3D, vecx, vecy, vecz, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      type(pointer3D), dimension(3) :: f
      Vec :: vec3D, vecx, vecy, vecz
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:), pointer :: rhsx, rhsy, rhsz
      PetscScalar, dimension(:, :, :), pointer :: rhs

      call VecZeroEntries(vecx, ierr)
      call VecZeroEntries(vecy, ierr)
      call VecZeroEntries(vecz, ierr)
      call VecZeroEntries(vec3D, ierr)

      call rhs_x(ctx, grid, B, f(1)%ptr, vecx, ierr)
      call rhs_random_y(ctx, grid, B, f(2)%ptr, vecy, ierr)
      call rhs_z(ctx, grid, B, f(3)%ptr, vecz, ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, vec3D, rhs, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAx, vecx, rhsx, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAy, vecy, rhsy, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAz, vecz, rhsz, ierr)
      CHKERRQ(ierr)

      do k=ctx%zs-1, ctx%ze-1
         do j=ctx%ys-1, ctx%ye-1
            do i=ctx%xs-1, ctx%xe-1
               rhs(i, j, k) = rhsx(i) * rhsy(j) * rhsz(k)
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, vec3D, rhs, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAx, vecx, rhsx, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAy, vecy, rhsy, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAz, vecz, rhsz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine rhs_random_init_separable

   subroutine rhs_init_separable(ctx, grid, B, f, vec3D, vecx, vecy, vecz, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      type(pointer3D), dimension(3) :: f
      Vec :: vec3D, vecx, vecy, vecz
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:), pointer :: rhsx, rhsy, rhsz
      PetscScalar, dimension(:, :, :), pointer :: rhs

      call VecZeroEntries(vecx, ierr)
      call VecZeroEntries(vecy, ierr)
      call VecZeroEntries(vecz, ierr)
      call VecZeroEntries(vec3D, ierr)

      call rhs_x(ctx, grid, B, f(1)%ptr, vecx, ierr)
      call rhs_y(ctx, grid, B, f(2)%ptr, vecy, ierr)
      call rhs_z(ctx, grid, B, f(3)%ptr, vecz, ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, vec3D, rhs, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAx, vecx, rhsx, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAy, vecy, rhsy, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAz, vecz, rhsz, ierr)
      CHKERRQ(ierr)

      do k=ctx%zs-1, ctx%ze-1
         do j=ctx%ys-1, ctx%ye-1
            do i=ctx%xs-1, ctx%xe-1
               rhs(i, j, k) = rhsx(i) * rhsy(j) * rhsz(k)
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, vec3D, rhs, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAx, vecx, rhsx, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAy, vecy, rhsy, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAz, vecz, rhsz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine rhs_init_separable

   subroutine rhs_reciprocal_x(ctx, grid, B, vvx, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscScalar, dimension(:), pointer :: vvx
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), pointer :: vvec

      associate (nodes => B%nodes, deg => B%deg, xs => ctx%xs, xe => ctx%xe, &
         Xk => B%Xk, Nx => grid%Nx, xNc => grid%xNc)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 15
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAx, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! get looping index
      call check_index(xs, xe, ctx%mx, deg, Nx+2*deg, kmin, kmax)

      vvec = 0.0
      ! calculate \int B(x) rhs(x) dx
      ! loop over B_splines
      do i=xs-1, xe-1
         ! loop over all knot intervals
         do k=kmin, kmax
            vvec(i) = vvec(i) + quad_reciprocal_x_B(grid, B, Xk(k), Xk(k+1), i+1+deg, 0, vvx)
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAx, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine rhs_reciprocal_x

   subroutine rhs_x(ctx, grid, B, ptr, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), pointer :: vvec

      associate (nodes => B%nodes, deg => B%deg, xs => ctx%xs, xe => ctx%xe, &
         Xk => B%Xk, Nx => grid%Nx, xNc => grid%xNc)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 15
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAx, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! get looping index
      call check_index(xs, xe, ctx%mx, deg, Nx+2*deg, kmin, kmax)

      vvec = 0.0
      ! calculate \int B(x) rhs(x) dx
      ! loop over B_splines
      do i=xs-1, xe-1
         ! loop over all knot intervals
         do k=kmin, kmax
            vvec(i) = vvec(i) + quad_x_B(grid, B, Xk(k), Xk(k+1), i+1+deg, 0, ptr, 0.0)
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAx, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine rhs_x

   subroutine rhs_random_y(ctx, grid, B, ptr, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, k
      PetscInt :: kmin, kmax
      PetscReal :: val
      PetscScalar, dimension(:), pointer :: vvec

      associate (nodes => B%nodes, deg => B%deg, ys => ctx%ys, &
         ye => ctx%ye, Yk => B%Yk, Ny => grid%Ny)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 15
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAy, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! get looping index
      call check_index(ys, ye, ctx%my, deg, Ny, kmin, kmax)

      vvec = 0.0
      ! calculate \int B(y) rhs(y) dy
      ! loop over B_splines
      do i=ys-1, ye-1
            ! loop over all knot intervals
            do k=kmin, kmax

               call PetscRandomGetValueReal(B%random_ctx, val, ierr)
               val = 0.1 * val

               vvec(i) = vvec(i) + quad_y_B(grid, B, Yk(k), Yk(k+1), i+1, 0, ptr, val)
            end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAy, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine rhs_random_y

   subroutine rhs_y(ctx, grid, B, ptr, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), pointer :: vvec

      associate (nodes => B%nodes, deg => B%deg, ys => ctx%ys, &
         ye => ctx%ye, Yk => B%Yk, Ny => grid%Ny)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 15
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAy, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! get looping index
      call check_index(ys, ye, ctx%my, deg, Ny, kmin, kmax)

      vvec = 0.0
      ! calculate \int B(y) rhs(y) dy
      ! loop over B_splines
      do i=ys-1, ye-1
         ! loop over all knot intervals
         do k=kmin, kmax
            vvec(i) = vvec(i) + quad_y_B(grid, B, Yk(k), Yk(k+1), i+1, 0, ptr, 0.0)
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAy, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine rhs_y

   subroutine rhs_z(ctx, grid, B, ptr, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      procedure (func1D), pointer :: ptr
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, k
      PetscInt :: kmin, kmax
      PetscScalar, dimension(:), pointer :: vvec

      associate (nodes => B%nodes, deg => B%deg, zs => ctx%zs, ze => ctx%ze, &
         Zk => B%Zk, Nz => grid%Nz, zNc => grid%zNc)

      ! initialize gauss points and weights (nodes give exact quadrature)
      nodes = 15
      call gauss_init(B, ierr)
      CHKERRQ(ierr)

      ! get array
      call DMDAVecGetArrayF90(ctx%DAz, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! get looping index
      call check_index(zs, ze, ctx%mz, deg, Nz+2*deg, kmin, kmax)

      vvec = 0.0
      ! calculate \int B(x) rhs(x) dx
      ! loop over B_splines
      do i=zs-1, ze-1
         ! loop over all knot intervals
         do k=kmin, kmax
            vvec(i) = vvec(i) + quad_z_B(grid, B, Zk(k), Zk(k+1), i+1+deg, 0, ptr, 0.0)
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DAz, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! free gauss points and weights
      call gauss_free(B, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine rhs_z
end module rhs

! vim: set ts=3 sw=3 sts=3:
