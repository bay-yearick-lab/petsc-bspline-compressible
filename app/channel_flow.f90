! Driver: compressible turbulent channel flow (M = 1.5, Re = 3000) with RK3 time stepping.
#include <petsc/finclude/petscsys.h>

module misc
   use types
   use bspline
   use operator
   use rhs
   use mms
   use error
   use linalg
   use rk
   use output
   implicit none

   private

   public :: create_vector
   public :: destroy_vector
   public :: init_vector
   public :: nonlinear_variables
   public :: ns_io
contains

   subroutine create_vector(ctx, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local4, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local5, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local6, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local7, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local8, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local9, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local10, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local11, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local12, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local13, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local14, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local15, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local16, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local17, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local18, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp4, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp5, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp6, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp7, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp8, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp9, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%const, ierr)
      CHKERRQ(ierr)

      call VecDuplicate(B%rhs, B%tempRho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%tempU, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%tempV, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%tempW, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%tempT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%invrho, ierr)
      CHKERRQ(ierr)

      call VecDuplicate(B%rhs, B%Rho1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%U1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%V1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%W1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%T1, ierr)
      CHKERRQ(ierr)

      call VecDuplicate(B%rhs, B%Rho2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%U2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%V2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%W2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%T2, ierr)
      CHKERRQ(ierr)

      call VecDuplicate(B%rhs, B%rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_w, ierr)
      CHKERRQ(ierr)

      call VecDuplicate(B%rhs, B%rho_uu, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_uv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_uw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_vv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_vw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_ww, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_uT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_vT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%rho_wT, ierr)
      CHKERRQ(ierr)

      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%alphaz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine create_vector

   subroutine destroy_vector(ctx, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      call PetscRandomDestroy(B%random_ctx, ierr)
      CHKERRQ(ierr)

      ! create vectors (solution and rhs)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%const, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphay, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphaz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine destroy_vector

   subroutine init_vector(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      type(pointer3D), dimension(3) :: u, rho, temp, one

      ! init function pointers
      u(1)%ptr => fone; u(2)%ptr => fu; u(3)%ptr => fone
      rho(1)%ptr => fone; rho(2)%ptr => fone; rho(3)%ptr => fone
      temp(1)%ptr => fone; temp(2)%ptr => ftemp; temp(3)%ptr => fone
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone

      ! create random context
      call PetscRandomCreate(PETSC_COMM_WORLD,  B%random_ctx, ierr)
      call PetscRandomSetInterval(B%random_ctx, complex(-1.0, 0.0), complex(1.0, 0.0), ierr)
      CHKERRQ(ierr)

      ! initialize constant vector
      call rhs_init_separable(ctx, grid, B, one, B%const, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! u velocity
      call rhs_random_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%rhs, B%U2, 'full', ierr)
      CHKERRQ(ierr)

      ! v velocity
      call VecZeroEntries(B%V2, ierr)
      CHKERRQ(ierr)

      ! w velocity
      call VecZeroEntries(B%W2, ierr)
      CHKERRQ(ierr)

      ! density
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%rhs, B%Rho2, 'full', ierr)
      CHKERRQ(ierr)

      ! temperature
      call rhs_init_separable(ctx, grid, B, temp, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecScale(B%rhs, complex(((B%gamma - 1.0) * B%Pr * B%M**2)/3.0, 0.0), ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(1.0, 0.0), B%const, ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%rhs, B%T2, 'full', ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine init_vector

   ! Finds 1/rho, rho_u, rho_v, rho_w,
   ! rho_uu, rho_uv, rho_uw, rho_vv, rho_vw, rho_ww,
   ! rho_uT, rho_vT, rho_wT
   subroutine nonlinear_variables(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscScalar, dimension(:, :, :), pointer :: vv

      ! global to local
      call DMGlobalToLocalBegin(ctx%DA, B%tempRho, INSERT_VALUES, B%local1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, B%tempRho, INSERT_VALUES, B%local1, ierr)

      call DMDAVecGetArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      call rhs_reciprocal_xyz(ctx, grid, B, vv, B%rhs, ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      ! find invrho
      call ksp_solve(ctx, grid, B, B%rhs, B%invrho, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho u/v/w
      call nonlinear_three(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         B%tempRho, B%tempU, B%temp1, B%local1, B%local2, &
         B%tempRho, B%tempV, B%temp2, B%local3, B%local4, &
         B%tempRho, B%tempW, B%temp3, B%local5, B%local6, ierr)

      ! linear solve
      call ksp_solve(ctx, grid, B, B%temp1, B%rho_u, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp2, B%rho_v, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp3, B%rho_w, 'full', ierr)
      CHKERRQ(ierr)

      ! find other nonlinear terms
      call nonlinear_nine(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         B%rho_u, B%tempU, B%temp1, B%local1, B%local2, &
         B%rho_u, B%tempV, B%temp2, B%local3, B%local4, &
         B%rho_u, B%tempW, B%temp3, B%local5, B%local6, &
         B%rho_v, B%tempV, B%temp4, B%local7, B%local8, &
         B%rho_v, B%tempW, B%temp5, B%local9, B%local10, &
         B%rho_w, B%tempW, B%temp6, B%local11, B%local12, &
         B%rho_u, B%tempT, B%temp7, B%local13, B%local14, &
         B%rho_v, B%tempT, B%temp8, B%local15, B%local16, &
         B%rho_w, B%tempT, B%temp9, B%local17, B%local18, ierr)

      ! linear solve
      call ksp_solve(ctx, grid, B, B%temp1, B%rho_uu, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp2, B%rho_uv, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp3, B%rho_uw, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp4, B%rho_vv, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp5, B%rho_vw, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp6, B%rho_ww, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp7, B%rho_uT, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp8, B%rho_vT, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp9, B%rho_wT, 'full', ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine nonlinear_variables

   subroutine ns_io(ctx, grid, B, n, nout, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: n, nout
      PetscErrorCode :: ierr

      if (mod(n, nout) .eq. 0) then
         ! io
         call io_fortran_scatter(ctx, B, B%Rho2, ierr)
         CHKERRQ(ierr)
         call io_fortran_write(ctx, grid, B, n, "rho", ierr)
         CHKERRQ(ierr)

         call io_fortran_scatter(ctx, B, B%U2, ierr)
         CHKERRQ(ierr)
         call io_fortran_write(ctx, grid, B, n, "u", ierr)
         CHKERRQ(ierr)

         call io_fortran_scatter(ctx, B, B%V2, ierr)
         CHKERRQ(ierr)
         call io_fortran_write(ctx, grid, B, n, "v", ierr)
         CHKERRQ(ierr)

         call io_fortran_scatter(ctx, B, B%W2, ierr)
         CHKERRQ(ierr)
         call io_fortran_write(ctx, grid, B, n, "w", ierr)
         CHKERRQ(ierr)

         call io_fortran_scatter(ctx, B, B%T2, ierr)
         CHKERRQ(ierr)
         call io_fortran_write(ctx, grid, B, n, "temp", ierr)
         CHKERRQ(ierr)
      end if

      ierr = 0
      return
   end subroutine ns_io

end module misc

program channel_flow
   use types
   use lattice
   use environment
   use bspline
   use misc
   use operator
   use rhs
   use mms
   use error
   use linalg
   use rk
   use output
   implicit none

   type(ctx_t) :: ctx
   type(grid_t) :: grid
   type(bspline_t) :: B
   type(rk_t) :: rk
   PetscErrorCode :: ierr

   character(len=30) :: otime
   PetscInt :: k, n, nmax, nout
   PetscBool :: found
   PetscReal :: Tbc
   PetscReal :: finalt, t, dt, scalar, rkt

   call environment_init(ctx, ierr)
   call lattice_init(ctx, grid, B, ierr)
   call grid_init(ctx, grid, ierr)

   ! init bspline
   call bspline_init(ctx, grid, B, ierr)

   ! create vec
   call create_vector(ctx, B, ierr)

   ! initialize operators
   call linear_init(ctx, grid, B, ierr)
   CHKERRQ(ierr)
   call nonlinear_init(ctx, grid, B, ierr)
   CHKERRQ(ierr)

   ! initialize ksp
   call ksp_init(ctx, grid, B, ierr)
   CHKERRQ(ierr)

   ! init rk coefficients
   call rk_init(rk, ierr)

   ! initialize time and time-step
   finalt = 30.0
   t = 0.0
   dt = 1.0e-3
   nmax = INT(finalt / dt)
   nout = 10

   ! runtime overrides: -nsteps <int> -output_every <int>
   call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                           '-nsteps', nmax, found, ierr)
   call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                           '-output_every', nout, found, ierr)
   rkt = 0.0

   ! flow variables
   B%mu = (0.5*(1.25 + 1.0))**(3.0/2.0) * &
      ((1.0 + 110.4/293.15)/((0.5*(1.25 + 1.0) + 110.4/293.15)))
   B%M = 1.5
   B%gamma = 1.4
   B%Re = 3000
   B%Pr = 0.72

   ! Tbc
   Tbc = 1.0

   ! initialize vectors
   call init_vector(ctx, grid, B, ierr)
   CHKERRQ(ierr)

   ! SET BC
   call channel_flow_bc(ctx, grid, B, B%U2, B%V2, B%W2, B%T2, Tbc, ierr)

   ! initialize io
   call io_init(ctx, grid, B, ierr)
   CHKERRQ(ierr)

   ! output files
   call ns_io(ctx, grid, B, 0, nout, ierr)

   ! print time
   write(otime,"(F10.6)") 0*dt
   call PetscPrintf(PETSC_COMM_WORLD, "time = "//otime//"\n", ierr)

   ! number of time-steps
   do n=1, nmax

      ! rk stage
      do k=1, rk%s

         rkt = t + rk%c(k) * dt

         if (k .eq. 1) then
            call VecCopy(B%Rho2, B%Rho1, ierr)
            CHKERRQ(ierr)
            call VecCopy(B%U2, B%U1, ierr)
            CHKERRQ(ierr)
            call VecCopy(B%V2, B%V1, ierr)
            CHKERRQ(ierr)
            call VecCopy(B%W2, B%W1, ierr)
            CHKERRQ(ierr)
            call VecCopy(B%T2, B%T1, ierr)
            CHKERRQ(ierr)
         else
            scalar = (rk%a(k, k-1) - rk%b(k-1)) * dt

            call VecAYPX(B%Rho1, complex(scalar, 0.0), B%Rho2, ierr)
            CHKERRQ(ierr)
            call VecAYPX(B%U1, complex(scalar, 0.0), B%U2, ierr)
            CHKERRQ(ierr)
            call VecAYPX(B%V1, complex(scalar, 0.0), B%V2, ierr)
            CHKERRQ(ierr)
            call VecAYPX(B%W1, complex(scalar, 0.0), B%W2, ierr)
            CHKERRQ(ierr)
            call VecAYPX(B%T1, complex(scalar, 0.0), B%T2, ierr)
            CHKERRQ(ierr)
         end if

         ! SET BC
         call channel_flow_bc(ctx, grid, B, B%U1, B%V1, B%W1, B%T1, Tbc, ierr)

         call VecCopy(B%Rho1, B%tempRho, ierr)
         CHKERRQ(ierr)
         call VecCopy(B%U1, B%tempU, ierr)
         CHKERRQ(ierr)
         call VecCopy(B%V1, B%tempV, ierr)
         CHKERRQ(ierr)
         call VecCopy(B%W1, B%tempW, ierr)
         CHKERRQ(ierr)
         call VecCopy(B%T1, B%tempT, ierr)
         CHKERRQ(ierr)

         ! Find all non-linear alpha's
         call nonlinear_variables(ctx, grid, B, ierr)

         ! Mass equation
         call rhs_mass(ctx, grid, B, B%rho_u, B%rho_v, B%rho_w, B%temp1, ierr)
         CHKERRQ(ierr)
         call ksp_solve(ctx, grid, B, B%temp1, B%tempRho, 'full', ierr)
         CHKERRQ(ierr)

         ! Momentum equation
         ! i = 1
         call rhs_momentum(ctx, grid, B, 1, B%mu, B%gamma, B%M, B%Re, &
               B%tempRho, B%tempU, B%tempV, B%tempW, B%tempT, &
               B%invrho, B%rho_u, B%rho_v, B%rho_w, &
               B%rho_uu, B%rho_uv, B%rho_uw, B%temp1, B%temp2, B%temp3, ierr)

         call bodyforce(ctx, grid, B, B%mu, B%Re, B%tempU, B%tempV, &
            B%temp1, B%temp2, ierr)
         call VecAXPY(B%temp3, complex(1.0, 0.0), B%temp2, ierr)

         call ksp_solve(ctx, grid, B, B%temp3, B%U1, 'full', ierr)
         CHKERRQ(ierr)

         ! i = 2
         call rhs_momentum(ctx, grid, B, 2, B%mu, B%gamma, B%M, B%Re, &
               B%tempRho, B%tempU, B%tempV, B%tempW, B%tempT, &
               B%invrho, B%rho_u, B%rho_v, B%rho_w, &
               B%rho_uv, B%rho_vv, B%rho_vw, B%temp1, B%temp2, B%temp3, ierr)

         call ksp_solve(ctx, grid, B, B%temp3, B%V1, 'full', ierr)
         CHKERRQ(ierr)

         ! i = 3
         call rhs_momentum(ctx, grid, B, 3, B%mu, B%gamma, B%M, B%Re, &
               B%tempRho, B%tempU, B%tempV, B%tempW, B%tempT, &
               B%invrho, B%rho_u, B%rho_v, B%rho_w, &
               B%rho_uw, B%rho_vw, B%rho_ww, B%temp1, B%temp2, B%temp3, ierr)

         call ksp_solve(ctx, grid, B, B%temp3, B%W1, 'full', ierr)
         CHKERRQ(ierr)

         ! Energy equation
         call rhs_energy(ctx, grid, B, B%mu, B%gamma, B%M, B%Re, B%Pr, &
               B%tempRho, B%tempU, B%tempV, B%tempW, B%tempT, &
               B%invrho, B%rho_u, B%rho_v, B%rho_w, &
               B%rho_uT, B%rho_vT, B%rho_wT, B%temp1, B%temp2, B%temp3, B%temp4, ierr)

         call ksp_solve(ctx, grid, B, B%temp4, B%T1, 'full', ierr)
         CHKERRQ(ierr)

         scalar = rk%b(k) * dt
         call VecAXPY(B%Rho2, complex(scalar, 0.0), B%Rho1, ierr)
         CHKERRQ(ierr)
         call VecAXPY(B%U2, complex(scalar, 0.0), B%U1, ierr)
         CHKERRQ(ierr)
         call VecAXPY(B%V2, complex(scalar, 0.0), B%V1, ierr)
         CHKERRQ(ierr)
         call VecAXPY(B%W2, complex(scalar, 0.0), B%W1, ierr)
         CHKERRQ(ierr)
         call VecAXPY(B%T2, complex(scalar, 0.0), B%T1, ierr)
         CHKERRQ(ierr)

         call channel_flow_bc(ctx, grid, B, B%U2, B%V2, B%W2, B%T2, Tbc, ierr)
      end do

      t = n * dt

      ! print time
      write(otime,"(F10.6)") n*dt
      call PetscPrintf(PETSC_COMM_WORLD, "time = "//otime//"\n", ierr)

      ! output files
      call ns_io(ctx, grid, B, n, nout, ierr)
   end do

   ! destroy
   call io_destroy(B, ierr)
   CHKERRQ(ierr)
   call destroy_vector(ctx, B, ierr)
   CHKERRQ(ierr)
   call lattice_free(ctx, ierr)
   CHKERRQ(ierr)
   call environment_free(ierr)
   CHKERRQ(ierr)
end program

! vim: set ts=3 sw=3 sts=3:
