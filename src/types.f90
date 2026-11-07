! Derived types: parallel context, grid, B-spline operators and Runge-Kutta coefficients.
# include <petsc/finclude/petscksp.h>

module types
   use petscsys
   use petscdmda
   use petscksp

   type ctx_t
      PetscInt :: xs, ys, zs
      PetscInt :: xe, ye, ze
      PetscInt :: xm, ym, zm
      PetscInt :: gxs, gys, gzs
      PetscInt :: gxe, gye, gze
      PetscInt :: gxm, gym, gzm
      PetscInt :: mx, my, mz
      PetscInt :: procx, procy, procz
      PetscMPIInt :: rank, rankx, ranky, rankz
      MPI_Comm :: commx, commy, commz
      PetscInt, dimension(:), allocatable :: rangex, rangey, rangez
      DM :: DA
      DM :: DAx, DAy, DAz
   end type ctx_t

   type grid_t
      PetscInt :: Nx, Ny, Nz
      PetscInt :: xNc, yNc, zNc
      PetscReal :: Lx, Ly, Lz
      PetscReal :: dx, dy, dz
      PetscReal, dimension(:), allocatable :: x, y, z
   end type grid_t

   type bspline_t
      PetscRandom :: random_ctx
      Mat :: Mx, My, Mz
      Mat :: BBx, BDx, BLx
      Mat :: BBy, BDy, BLy
      Mat :: BBz, BDz, BLz
      Mat :: Bx, By, Bz
      Vec :: const
      Vec :: rhsx, rhsy, rhsz
      Vec :: alphax, alphay, alphaz
      Vec :: solx, soly, solz
      Vec :: tempx, tempy, tempz
      Vec :: rhs, sol, alpha
      Vec :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
      Vec :: local1, local2, local3, local4, local5, local6, local7, local8, local9
      Vec :: local10, local11, local12, local13, local14, local15, local16, local17, local18
      Vec :: natural, vecseq
      VecScatter :: scatter
      KSP :: ksp, fullksp
      PC :: pc
      PetscInt :: deg
      PetscInt :: nodes
      PetscReal, dimension(2) :: old_ybc, new_ybc
      PetscReal, dimension(:), allocatable :: Xk, Yk, Zk
      PetscReal, dimension(:), allocatable :: wts, xi
      PetscReal, dimension(:), allocatable :: intx, inty, intz
      PetscReal, dimension(:, :, :), allocatable :: BBBx, BBBy, FBBBy, BBBz
      PetscReal, dimension(:, :, :), allocatable :: BBDx, BBDy, BBDz
      PetscReal, dimension(:, :, :, :), allocatable :: BBBBx, BBBBy, BBBBz
      PetscReal, dimension(:, :, :, :), allocatable :: BBBDx, BBBDy, BBBDz
      Vec :: FFTx, FFTz
      Mat :: FFTmatx, FFTmatx1, FFTmatz, FFTmatz1
      Vec :: inputx, outputx, backwardx, inputz, outputz, backwardz
      Vec :: inputx1, outputx1, backwardx1, inputz1, outputz1, backwardz1

      ! Navier-stokes variables
      PetscReal :: mu, M, gamma, Re, Pr
      Vec :: tempRho, tempU, tempV, tempW, tempT
      Vec :: Rho1, Rho2, U1, U2, V1, V2, W1, W2, T1, T2
      Vec :: invrho, rho_u, rho_v, rho_w
      Vec :: rho_uu, rho_uv, rho_uw, rho_vv, rho_vw, rho_ww, rho_uT, rho_vT, rho_wT

   end type bspline_t

   type rk_t
      Vec :: temp1, temp2
      PetscInt :: s
      PetscReal, dimension(:), allocatable :: b, c
      PetscReal, dimension(:, :), allocatable :: a
   end type rk_t

   abstract interface
      function func1D(v, t)
         PetscReal :: func1D, v, t
      end function func1D

      function func3D(x, y, z, t)
         PetscReal :: func3D, x, y, z, t
      end function func3D
   end interface

   type pointer3D
      procedure (func1D), pointer, nopass :: ptr
   end type pointer3D

   procedure (func1D), pointer :: f1d_ptr => null()
end module

! vim: set ts=3 sw=3 sts=3:
