! Manufactured-solution checks for every discrete operator; select one with -test <name>.
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>

module test
   use types
   use bspline
   use operator
   use rhs
   use mms
   use error
   use linalg
   use rk
   implicit none

   private

   public :: test_linearsolve_x
   public :: test_linearsolve_y
   public :: test_linearsolve_z
   public :: test_linearsolve_xyz
   public :: test_operator
   public :: test_dmdagetray
   public :: test_dmda
   public :: test_bspline_1d
   public :: test_nonlinear_sum_BBB
   public :: test_nonlinear_sum_BBBB
   public :: test_divergence
   public :: test_nonlinear_divergence
   public :: test_nonlinear_x
   public :: test_nonlinear_y
   public :: test_nonlinear_z
   public :: test_rk
   public :: test_reciprocal
   public :: test_serial_interpolate_1d
   public :: test_serial_interpolate_3d
   public :: test_momentum_energy
   public :: test_reciprocal_1d
   public :: test_reciprocal_3d
   public :: test_rhs
   public :: test_nonlinear_combined_rhs
   public :: test_petsc_write
   public :: test_fortran_write
   public :: run_test
contains

   ! Run one manufactured-solution check selected by name (PETSc option -test).
   subroutine run_test(name, ctx, grid, B, rk, ierr)
      implicit none

      character(len=*), intent(in) :: name
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      type(rk_t) :: rk
      PetscErrorCode, intent(out) :: ierr

      call PetscPrintf(PETSC_COMM_WORLD, "test: "//trim(name)//"\n", ierr)

      select case (trim(name))
      case ('dmda')
         call test_dmda(ctx, ierr)
      case ('operator')
         call test_operator(ctx, grid, B, ierr)
      case ('linearsolve_x')
         call test_linearsolve_x(ctx, grid, B, ierr)
      case ('linearsolve_y')
         call test_linearsolve_y(ctx, grid, B, ierr)
      case ('linearsolve_z')
         call test_linearsolve_z(ctx, grid, B, ierr)
      case ('linearsolve_xyz')
         call test_linearsolve_xyz(ctx, grid, B, ierr)
      case ('nonlinear_x')
         call test_nonlinear_x(ctx, grid, B, ierr)
      case ('nonlinear_y')
         call test_nonlinear_y(ctx, grid, B, ierr)
      case ('nonlinear_z')
         call test_nonlinear_z(ctx, grid, B, ierr)
      case ('nonlinear_sum_BBB')
         call test_nonlinear_sum_BBB(ctx, grid, B, ierr)
      case ('nonlinear_sum_BBBB')
         call test_nonlinear_sum_BBBB(ctx, grid, B, ierr)
      case ('divergence')
         call test_divergence(ctx, grid, B, ierr)
      case ('nonlinear_divergence')
         call test_nonlinear_divergence(ctx, grid, B, ierr)
      case ('rk')
         call test_rk(ctx, grid, B, rk, ierr)
      case ('reciprocal')
         call test_reciprocal(ctx, grid, B, ierr)
      case ('reciprocal_1d');         call test_reciprocal_1d(ctx, grid, B, ierr)
      case ('reciprocal_3d');         call test_reciprocal_3d(ctx, grid, B, ierr)
      case ('convection')
         call test_convection(ctx, grid, B, ierr)
      case ('viscous')
         call test_viscous(ctx, grid, B, ierr)
      case ('energy_viscous')
         call test_energy_viscous(ctx, grid, B, ierr)
      case ('momentum_energy')
         call test_momentum_energy(ctx, grid, B, ierr)
      case ('interpolate_1d');        call test_serial_interpolate_1d(ctx, grid, B, ierr)
      case ('interpolate_3d');        call test_serial_interpolate_3d(ctx, grid, B, ierr)
      case ('rhs')
         call test_rhs(ctx, grid, B, ierr)
      case ('nonlinear_combined_rhs')
         call test_nonlinear_combined_rhs(ctx, grid, B, ierr)
      case default
         call PetscPrintf(PETSC_COMM_WORLD, "unknown test name\n", ierr)
         ierr = 1
      end select

      return
   end subroutine run_test

   subroutine test_dmda(ctx, ierr)
      implicit none

      type(ctx_t) :: ctx
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: xpoints, xproc
      PetscInt :: ypoints, yproc
      PetscInt :: zpoints, zproc
      PetscInt :: Nxs, Nxm, Nys, Nym, Nzs, Nzm
      PetscInt :: xs, xm, ys, ym, zs, zm
      PetscInt :: rankx, sizex, mygroupx, worldgroup
      PetscInt :: ranky, sizey, mygroupy
      PetscInt :: rankz, sizez, mygroupz
      PetscInt, dimension(:), allocatable :: allrankx, globalrankx
      PetscInt, dimension(:), allocatable :: allranky, globalranky
      PetscInt, dimension(:), allocatable :: allrankz, globalrankz
      MPI_Comm :: globalcomm, commx, commy, commz

      call DMDAGetInfo(ctx%DAx, PETSC_NULL_INTEGER,                                &
                       xpoints, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,            &
                       xproc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,              &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      call DMDAGetInfo(ctx%DAy, PETSC_NULL_INTEGER,                                &
                       ypoints, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,            &
                       yproc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,              &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      call DMDAGetInfo(ctx%DAz, PETSC_NULL_INTEGER,                                &
                       zpoints, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,            &
                       zproc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,              &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      call DMDAGetCorners(ctx%DA, Nxs, Nys, Nzs, Nxm, Nym, Nzm, ierr)

      call DMDAGetCorners(ctx%DAx, xs, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   xm, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      call DMDAGetCorners(ctx%DAy, ys, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   ym, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      call DMDAGetCorners(ctx%DAz, zs, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   zm, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      ! should print to be zeros

      call PetscObjectGetComm(ctx%DA, globalcomm, ierr)
      call PetscObjectGetComm(ctx%DAx, commx, ierr)
      call PetscObjectGetComm(ctx%DAy, commy, ierr)
      call PetscObjectGetComm(ctx%DAz, commz, ierr)

      ! check which ranks are in each of the communicators
      call MPI_Comm_rank(commx, rankx, ierr)
      call MPI_Comm_size(commx, sizex, ierr)
      call MPI_Comm_rank(commy, ranky, ierr)
      call MPI_Comm_size(commy, sizey, ierr)
      call MPI_Comm_rank(commz, rankz, ierr)
      call MPI_Comm_size(commz, sizez, ierr)

      ! get processors in those communicators
      allocate(allrankx(sizex))
      allocate(globalrankx(sizex))
      allocate(allranky(sizey))
      allocate(globalranky(sizey))
      allocate(allrankz(sizez))
      allocate(globalrankz(sizez))

      call MPI_Allgather(rankx, 1, MPI_INT, allrankx, 1, MPI_INT, commx, ierr)
      call MPI_Allgather(ranky, 1, MPI_INT, allranky, 1, MPI_INT, commy, ierr)
      call MPI_Allgather(rankz, 1, MPI_INT, allrankz, 1, MPI_INT, commz, ierr)

      ! group the processors in each communicators
      call MPI_Comm_group(PETSC_COMM_WORLD, worldgroup, ierr)
      call MPI_Comm_group(commx, mygroupx, ierr)
      call MPI_Comm_group(commy, mygroupy, ierr)
      call MPI_Comm_group(commz, mygroupz, ierr)

      ! translate from local comm ranks to global ranks
      call MPI_Group_translate_ranks(mygroupx, sizex, allrankx, worldgroup, globalrankx, ierr)
      call MPI_Group_translate_ranks(mygroupy, sizey, allranky, worldgroup, globalranky, ierr)
      call MPI_Group_translate_ranks(mygroupz, sizez, allrankz, worldgroup, globalrankz, ierr)

      ierr = 0
      return
   end subroutine test_dmda

   subroutine test_operator(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      MPI_Comm :: comm
      PetscBool :: flag

      call test_bspline_1d(ctx, grid, B, ierr)
      call mass_y(ctx, grid, B, ierr)

      ! check symmetry
      call MatIsSymmetric(B%My, 1.0e-14, flag, ierr)
      CHKERRQ(ierr)

      if (flag) then
         call PetscPrintf(comm, "MATRIX IS SYMMETRIC. \n", ierr)
      else
         call PetscPrintf(comm, "MATRIX IS NON-SYMMETRIC. \n", ierr)
      end if
      ierr = 0
      return
   end subroutine test_operator

   subroutine test_bspline_1d(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal, dimension (ctx%xs:ctx%xe) :: xspline_value
      PetscReal, dimension (ctx%ys:ctx%ye) :: yspline_value
      PetscReal, dimension (ctx%zs:ctx%ze) :: zspline_value

      ! x b-splines
      call initxknot(ctx, grid, B, ierr)
      do i=ctx%xs, ctx%xe
         xspline_value(i) = xdBIJK(grid, B, 1, i, B%deg, 0.5 * (B%Xk(i) + B%Xk(i-1)))
      end do

      ! y b-splines
      call inityknot(ctx, grid, B, ierr)
      do i=ctx%ys, ctx%ye
         yspline_value(i) = ydBIJK(grid, B, 2, i, B%deg, 0.5)
      end do

      ! z b-splines
      call initzknot(ctx, grid, B, ierr)
      do i=ctx%zs, ctx%ze
         zspline_value(i) = zdBIJK(grid, B, 1, i, B%deg, 0.5)
      end do

      ierr = 0
      return
   end subroutine test_bspline_1d

   subroutine test_linearsolve_x(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      procedure (func1D), pointer :: ptr

      ptr => f1d

      ! B-spline init
      call test_bspline_1d(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      call VecDuplicate(B%rhsx, B%solx, ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_x(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_x(ctx, grid, B, ptr, B%rhsx, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_x(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, ptr, B%alphax, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_x(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solx, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_linearsolve_x

   subroutine test_linearsolve_y(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      KSP :: ksp
      Mat :: Factor
      PC :: pc
      PetscInt :: i
      MPI_Comm :: commy
      procedure (func1D), pointer :: ptr

      ptr => fcos4

      ! B-spline init
      call test_bspline_1d(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      call VecDuplicate(B%rhsy, B%soly, ierr)
      call VecDuplicate(B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_y(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_y(ctx, grid, B, ptr, B%rhsy, ierr)

      ! get comm from dmda object
      call PetscObjectGetComm(ctx%DAy, commy, ierr)
      CHKERRQ(ierr)

      ! create ksp
      call KSPCreate(ctx%commy, ksp, ierr)
      CHKERRQ(ierr)
      call KSPSetOperators(ksp, B%My, B%My, ierr)
      CHKERRQ(ierr)

      call KSPSetType(ksp, KSPGMRES, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0e-14, 1.0e-14, &
      PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

      call KSPSetFromOptions(ksp, ierr)
      CHKERRQ(ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Solve linear system \n", ierr)
      call KSPSolve(ksp, B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_y(ctx, grid, B, ptr, B%alphay, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call MatDestroy(B%My, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphay, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%soly, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempy, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_linearsolve_y

   subroutine test_linearsolve_z(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Mat :: FFTmat, FFTmat1
      Vec :: input, output, backward
      Vec :: input1, output1, backward1
      PetscInt, dimension(1) :: FFTvec
      procedure (func1D), pointer :: ptr

      ptr => f1d

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      call VecDuplicate(B%rhsz, B%solz, ierr)
      call VecDuplicate(B%rhsz, B%alphaz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_z(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_z(ctx, grid, B, ptr, B%rhsz, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_z(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_z(ctx, grid, B, ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_z(ctx, grid, B, ptr, B%alphaz, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_z(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphaz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_linearsolve_z

   subroutine test_linearsolve_xyz(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      type(pointer3D), dimension(3) :: f
      procedure(func3D), pointer :: solptr

      f(1)%ptr => f1d; f(2)%ptr => f1d; f(3)%ptr => f1d
      solptr => gsin4sin4sin4

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
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

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, f, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      B%new_ybc(1) = f(2)%ptr(0.0, 0.0)
      B%new_ybc(2) = f(2)%ptr(grid%Ly, 0.0)

      call set_ybc(ctx, grid, B, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, B%alpha, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, solptr, B%alpha, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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
   end subroutine test_linearsolve_xyz

   subroutine test_dmdagetray(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:), pointer :: rhsy, tempy
      PetscScalar, dimension(:, :, :), pointer :: rhs

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      CHKERRQ(ierr)

      call VecSet(B%rhsy, 1.0 + 0*PETSC_i, ierr)

      call DMDAVecGetArrayF90(ctx%DA, B%rhs, rhs, ierr)
      call DMDAVecGetArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)

      do k=ctx%zs-1, ctx%ze-1
         do j=ctx%ys-1, ctx%ye-1
            do i=ctx%xs-1, ctx%xe-1
               if (i .eq. 5 .and. k.eq. 2) then
                  rhs(i, j, k) = rhsy(j)
               end if
            end do
         end do
      end do

      call DMDAVecRestoreArrayF90(ctx%DA, B%rhs, rhs, ierr)
      call DMDAVecRestoreArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)

      call DMDAVecGetArrayF90(ctx%DA, B%rhs, rhs, ierr)

      do k=ctx%zs-1, ctx%ze-1
         do i=ctx%xs-1, ctx%xe-1
            if (i .eq. 5 .and. k.eq. 2) then
               call DMDAVecGetArrayF90(ctx%DAy, B%tempy, tempy, ierr)
               tempy = rhs(i, :, k)
               call DMDAVecRestoreArrayF90(ctx%DAy, B%tempy, tempy, ierr)
               call VecView(B%tempy, PETSC_VIEWER_STDOUT_(ctx%commy), ierr)
            end if
         end do
      end do

      print*, ctx%procy

      call DMDAVecRestoreArrayF90(ctx%DA, B%rhs, rhs, ierr)

      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempy, ierr)
      CHKERRQ(ierr)
   end subroutine test_dmdagetray

   subroutine test_nonlinear_z(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, alpha_rho
      type (pointer3D), dimension(3) :: u, rho, rho_u

      u(1)%ptr => fcos4; rho(1)%ptr => fcos2; rho_u(1)%ptr => fcos2cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%alphaz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%solz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_z(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_z(ctx, grid, B, u(1)%ptr, B%rhsz, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_z(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_z(ctx, grid, B, ierr)

      call VecCopy(B%alphaz, alpha_u, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_z(ctx, grid, B, u(1)%ptr, alpha_u, ierr)
      CHKERRQ(ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_z(ctx, grid, B, rho(1)%ptr, B%rhsz, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_z(ctx, grid, B, ierr)

      call VecCopy(B%alphaz, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_z(ctx, grid, B, rho(1)%ptr, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize x nonlinear term
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_z(ctx, grid, B, ierr)

      ! nonlinear sums
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_z(ctx, grid, B, alpha_rho, alpha_u, B%rhsz, ierr)

      ! solve
      call fft_solve_z(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_z(ctx, grid, B, rho_u(1)%ptr, B%alphaz, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_z(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphaz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_rho, ierr)
      CHKERRQ(ierr)

      if (allocated(B%BBBz)) deallocate(B%BBBz)
   end subroutine test_nonlinear_z

   subroutine test_nonlinear_y(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      KSP :: ksp
      Vec :: alpha_u, alpha_rho
      type (pointer3D), dimension(3) :: u, rho, rho_u

      u(1)%ptr => fcos4; rho(1)%ptr => fcos2; rho_u(1)%ptr => fcos2cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%soly, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_y(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_y(ctx, grid, B, u(1)%ptr, B%rhsy, ierr)

      ! create ksp
      call KSPCreate(ctx%commy, ksp, ierr)
      CHKERRQ(ierr)
      call KSPSetOperators(ksp, B%My, B%My, ierr)
      CHKERRQ(ierr)

      call KSPSetType(ksp, KSPGMRES, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0e-14, 1.0e-14, &
      PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

      call KSPSetFromOptions(ksp, ierr)
      CHKERRQ(ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Solve linear system \n", ierr)
      call KSPSolve(ksp, B%rhsy, alpha_u, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_y(ctx, grid, B, u(1)%ptr, alpha_u, ierr)
      CHKERRQ(ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_y(ctx, grid, B, rho(1)%ptr, B%rhsy, ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Solve linear system \n", ierr)
      call KSPSolve(ksp, B%rhsy, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_y(ctx, grid, B, rho(1)%ptr, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize y nonlinear term
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_y(ctx, grid, B, ierr)

      ! nonlinear sums
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_y(ctx, grid, B, alpha_rho, alpha_u, B%rhsy, ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Solve linear system \n", ierr)
      call KSPSolve(ksp, B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_y(ctx, grid, B, rho_u(1)%ptr, B%alphay, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call KSPDestroy(ksp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphay, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%soly, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_rho, ierr)
      CHKERRQ(ierr)

      if (allocated(B%BBBy)) deallocate(B%BBBy)
   end subroutine test_nonlinear_y

   subroutine test_nonlinear_x(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, alpha_rho
      type (pointer3D), dimension(3) :: u, rho, rho_u

      u(1)%ptr => fcos4; rho(1)%ptr => fcos2; rho_u(1)%ptr => fcos2cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%solx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_x(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_x(ctx, grid, B, u(1)%ptr, B%rhsx, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_x(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      call VecCopy(B%alphax, alpha_u, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, u(1)%ptr, alpha_u, ierr)
      CHKERRQ(ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_x(ctx, grid, B, rho(1)%ptr, B%rhsx, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      call VecCopy(B%alphax, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Evaluate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, rho(1)%ptr, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! initialize x nonlinear term
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_x(ctx, grid, B, ierr)

      ! nonlinear sums
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_x(ctx, grid, B, alpha_rho, alpha_u, B%rhsx, ierr)

      ! solve
      call fft_solve_x(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! calculate error
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, rho_u(1)%ptr, B%alphax, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_x(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_rho, ierr)
      CHKERRQ(ierr)

      if (allocated(B%BBBx)) deallocate(B%BBBx)
   end subroutine test_nonlinear_x

   subroutine test_nonlinear_sum_BBB(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, alpha_rho, alpha_rho_u
      type (pointer3D), dimension(3) :: u, rho
      procedure(func3D), pointer :: sol_u, sol_rho, sol_rho_u

      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      rho(1)%ptr => fcos2; rho(2)%ptr => fsin2; rho(3)%ptr => fcos2

      sol_u => gcos4sin4cos4
      sol_rho => gcos2sin2cos2
      sol_rho_u => gcos2cos4sin2sin4cos2cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
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

      ! test u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! test rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_u, B%rhs, B%local1, B%local2, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
   end subroutine test_nonlinear_sum_BBB

   subroutine test_nonlinear_sum_BBBB(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, alpha_v, alpha_rho, alpha_rho_u_v
      type (pointer3D), dimension(3) :: u, v, rho
      procedure(func3D), pointer :: sol_u, sol_v, sol_rho, sol_rho_u_v

      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      rho(1)%ptr => fcos2; rho(2)%ptr => fsin2; rho(3)%ptr => fcos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_rho => gcos2sin2cos2
      sol_rho_u_v => gcos2cos4sin4sin2sin4cos4cos2cos4sin4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u_v, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! test u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! test v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! test rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_x(ctx, grid, B, ierr)
      call nonlinear_y(ctx, grid, B, ierr)
      call nonlinear_z(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_BBBB(ctx, grid, B, B%BBBBx, B%BBBBy, B%BBBBz, &
         alpha_rho, alpha_u, alpha_v, B%rhs, B%local1, B%local2, B%local3, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u_v, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u_v, alpha_rho_u_v, ierr)
      CHKERRQ(ierr)

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local3, ierr)
      CHKERRQ(ierr)

      call VecDestroy(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_rho_u_v, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBBBx)) deallocate(B%BBBBx)
      if (allocated(B%BBBBy)) deallocate(B%BBBBy)
      if (allocated(B%BBBBz)) deallocate(B%BBBBz)
   end subroutine test_nonlinear_sum_BBBB

   subroutine test_nonlinear_divergence(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, alpha_rho, alpha_rho_u, div_alpha_rho_u
      type(pointer3D), dimension(3) :: u, rho
      procedure(func3D), pointer :: sol_u, sol_rho, sol_rho_u
      procedure(func3D), pointer :: sol_div_rho_u

      u(1)%ptr => fcos4; u(2)%ptr => fcos4; u(3)%ptr => fcos4
      rho(1)%ptr => fcos2; rho(2)%ptr => fcos4; rho(3)%ptr => fcos2

      sol_u => gcos4cos4cos4
      sol_rho => gcos2cos4cos2
      sol_rho_u => gcos2cos4cos4cos4cos2cos4
      sol_div_rho_u => div_cos2cos4cos4cos4cos2cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, div_alpha_rho_u, ierr)
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

      ! test u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! test rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_x(ctx, grid, B, ierr)
      call nonlinear_y(ctx, grid, B, ierr)
      call nonlinear_z(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums \n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBBy, B%BBBz, &
         alpha_rho, alpha_u, B%rhs, B%local1, B%local2, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute divergence \n", ierr)
      call divergence(ctx, grid, B, alpha_rho_u, alpha_rho_u, alpha_rho_u, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, div_alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_div_rho_u, div_alpha_rho_u, ierr)
      CHKERRQ(ierr)

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_nonlinear_divergence

   subroutine test_divergence(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: alpha_u, div_alpha
      type (pointer3D), dimension(3) :: u, rho
      procedure(func3D), pointer :: sol_u, sol_div

      u(1)%ptr => fcos4; u(2)%ptr => fcos4; u(3)%ptr => fcos4

      sol_u => gcos4cos4cos4
      sol_div => div_cos4cos4cos4

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, div_alpha, ierr)
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

      ! test u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute divergence \n", ierr)
      call divergence(ctx, grid, B, alpha_u, alpha_u, alpha_u, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, div_alpha, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_div, div_alpha, ierr)
      CHKERRQ(ierr)

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_divergence

   subroutine test_rk(ctx, grid, B, rk, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      type(rk_t) :: rk
      PetscErrorCode, intent(out) :: ierr

      Vec :: x, y
      PetscInt :: k, n, nmax
      PetscReal :: finalt, t, dt
      procedure(func1D), pointer :: ptr, sol

      ptr => expt
      sol => INTexpt

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, x, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, y, ierr)
      CHKERRQ(ierr)

      ! init rk coefficients
      call rk_init(rk, ierr)

      ! initialize time and time-step
      finalt = 1.0
      t = 0.0
      dt = 0.5e-2
      nmax = INT(finalt / dt)

      ! initialize vector
      call VecSet(x, complex(ptr(0.0, t), 0.0), ierr)
      CHKERRQ(ierr)

      do n=1, nmax
         do k=1, rk%s

            if (k .eq. 1) then
               call VecCopy(x, y, ierr)
               CHKERRQ(ierr)
            else
               call VecSet(B%rhs, &
                  complex((rk%a(k, k-1) - rk%b(k-1)) * dt * ptr(0.0, t + rk%c(k-1) * dt), 0.0), ierr)
               CHKERRQ(ierr)

               call VecWAXPY(y, complex(1.0, 0.0), x, B%rhs, ierr)
               CHKERRQ(ierr)
            end if

            call VecSet(B%rhs, &
               complex(rk%b(k) * dt * ptr(0.0, t + rk%c(k) * dt), 0.0), ierr)
            CHKERRQ(ierr)

            call VecAXPY(x, complex(1.0, 0.0), B%rhs, ierr)
            CHKERRQ(ierr)
         end do
         t = n * dt
      end do

      ! set solution
      call VecSet(B%rhs, complex(sol(0.0, finalt), 0.0), ierr)
      CHKERRQ(ierr)

      ! check error
      call VecAXPY(x, complex(-1.0, 0.0), B%rhs, ierr)
      CHKERRQ(ierr)

      call VecNorm(x, NORM_INFINITY, t, ierr)
      CHKERRQ(ierr)

      print*, t
   end subroutine test_rk

   subroutine test_reciprocal(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i, j, k
      PetscReal :: vol
      Vec :: alpha_u, alpha_ru
      type (pointer3D), dimension(3) :: u, ru, one
      procedure(func3D), pointer :: sol_u, sol_ru
      PetscScalar, dimension(:, :, :), pointer :: vvec

      u(1)%ptr => fsin4; u(2)%ptr => fcos4; u(3)%ptr => fsin4
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone

      sol_u => gsin4cos4sin4_2
      sol_ru => rsin4cos4sin4_2

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_ru, ierr)
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

      ! test u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%temp, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'bc', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      vol = grid%dx * grid%dy * grid%dz

      call VecReciprocal(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecScale(alpha_u, complex(vol, 0.0), ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, alpha_u, alpha_ru, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_ru, alpha_ru, ierr)
      CHKERRQ(ierr)

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)

      call VecDestroy(alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDestroy(alpha_ru, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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
   end subroutine test_reciprocal

   subroutine test_convection(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val
      Vec :: conv_alpha, source
      Vec :: alpha_dudx, alpha_dudy, alpha_dudz
      Vec :: alpha_u, alpha_v, alpha_w, alpha_rho
      Vec :: alpha_rho_u, alpha_rho_v, alpha_rho_w
      Vec :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
      type (pointer3D), dimension(3) :: u, v, w, rho, one
      type (pointer3D), dimension(3) :: rho_u, rho_v, rho_w
      procedure(func3D), pointer :: sol_u, sol_v, sol_w, sol_rho
      procedure(func3D), pointer :: sol_rho_u, sol_rho_v, sol_rho_w

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2
      rho(1)%ptr => fcos2; rho(2)%ptr => fcos4; rho(3)%ptr => fcos2
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone
      rho_u(1)%ptr => fcos2cos4; rho_u(2)%ptr => fcos4sin4; rho_u(3)%ptr => fcos2cos4
      rho_v(1)%ptr => fcos2sin4; rho_v(2)%ptr => fcos4cos4; rho_v(3)%ptr => fcos2sin4
      rho_w(1)%ptr => fcos2cos2; rho_w(2)%ptr => fcos4sin2; rho_w(3)%ptr => fcos2cos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2
      sol_rho => gcos2cos4cos2_2
      sol_rho_u => grho_u
      sol_rho_v => grho_v
      sol_rho_w => grho_w

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uu, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, conv_alpha, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, source, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      ! init rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%temp1, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      ! find rho u
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho u\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_u, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      ! find rho v
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho v\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_v, alpha_rho_v, ierr)
      CHKERRQ(ierr)

      ! find rho w
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho w\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_w, alpha_rho_w, ierr)
      CHKERRQ(ierr)

      do i=1, 3
         if (i .eq. 1) then
            ! find convect u
            call PetscPrintf(PETSC_COMM_WORLD, "Convection term u \n", ierr)

            call convection(ctx, grid, B, &
               alpha_u, alpha_v, alpha_w, alpha_u, &
               alpha_rho_u, alpha_rho_v, alpha_rho_w, &
               alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, conv_alpha, B%rhs, ierr)

         else if (i .eq. 2) then
            ! find convect v
            call PetscPrintf(PETSC_COMM_WORLD, "Convection term v \n", ierr)

            call convection(ctx, grid, B, &
               alpha_u, alpha_v, alpha_w, alpha_v, &
               alpha_rho_u, alpha_rho_v, alpha_rho_w, &
               alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, conv_alpha, B%rhs, ierr)

         else if (i .eq. 3) then
            ! find convect w
            call PetscPrintf(PETSC_COMM_WORLD, "Convection term w \n", ierr)

            call convection(ctx, grid, B, &
               alpha_u, alpha_v, alpha_w, alpha_w, &
               alpha_rho_u, alpha_rho_v, alpha_rho_w, &
               alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, conv_alpha, B%rhs, ierr)
         end if

         ! compute solution
         call compute_sol_convection(ctx, grid, B, i, source, conv_alpha, ierr)

         ! compute error
         call VecAXPY(conv_alpha, complex(-1.0, 0.0), B%rhs, ierr)
         call VecNorm(conv_alpha, NORM_INFINITY, val, ierr)

         print*, val
      end do

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_convection

   subroutine test_viscous(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val, mu, Re
      Vec :: visc_alpha, source
      Vec :: alpha_u, alpha_v, alpha_w
      type (pointer3D), dimension(3) :: u, v, w
      type (pointer3D), dimension(3) :: s1, s2, s3, s4, s5
      procedure(func3D), pointer :: sol_u, sol_v, sol_w

      mu = 1.0; Re = 1.0

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, source, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, visc_alpha, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing u \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call VecCopy(B%rhs, visc_alpha, ierr)
      call VecScale(visc_alpha, complex(-16.0, 0.0), ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing v \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing w \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      do i=1, 3
         if (i .eq. 1) then

            ! find viscous u
            call PetscPrintf(PETSC_COMM_WORLD, "Viscous term i=1 \n", ierr)
            call viscous(ctx, grid, B, i, mu, Re, alpha_u, alpha_v, alpha_w, B%temp, B%rhs, ierr)

         else if (i .eq. 2) then

            ! find viscous v
            call PetscPrintf(PETSC_COMM_WORLD, "Viscous term i=2 \n", ierr)
            call viscous(ctx, grid, B, i, mu, Re, alpha_u, alpha_v, alpha_w, B%temp, B%rhs, ierr)

         else if (i .eq. 3) then

            ! find viscous w
            call PetscPrintf(PETSC_COMM_WORLD, "Viscous term i=3 \n", ierr)
            call viscous(ctx, grid, B, i, mu, Re, alpha_u, alpha_v, alpha_w, B%temp, B%rhs, ierr)
         end if

            ! compute solution
            call compute_sol_viscous(ctx, grid, B, i, mu, Re, source, visc_alpha, ierr)

            ! compute error
            call VecAXPY(visc_alpha, complex(-1.0, 0.0), B%rhs, ierr)
            call VecNorm(visc_alpha, NORM_INFINITY, val, ierr)

            print*, val
      end do

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

   end subroutine test_viscous

   subroutine test_energy_viscous(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val, mu, gamma, M, Re
      Vec :: energy_tau, source
      Vec :: alpha_u, alpha_v, alpha_w
      type (pointer3D), dimension(3) :: u, v, w
      procedure(func3D), pointer :: sol_u, sol_v, sol_w

      mu = 1.0; gamma = 1.4; M = 1.0; Re = gamma * (gamma - 1.0)

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, energy_tau, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, source, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      ! init nonlinear operators
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Energy viscous \n", ierr)

      call VecZeroEntries(B%rhs, ierr)
      do i=1, 3
         call energy_viscous(ctx, grid, B, i, mu, gamma, M, Re, &
            alpha_u, alpha_v, alpha_w, B%temp1, B%temp2, B%temp3, ierr)
         call VecAXPY(B%rhs, complex(1.0, 0.0), B%temp3, ierr)
      end do

      ! compute solution
      call compute_sol_energy_viscous(ctx, grid, B, &
         mu, gamma, M, Re, B%temp1, energy_tau, ierr)

      ! compute error
      call VecAXPY(energy_tau, complex(-1.0, 0.0), B%rhs, ierr)
      call VecNorm(energy_tau, NORM_INFINITY, val, ierr)

      print*, val

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_energy_viscous

   subroutine test_momentum_energy(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val, mu, gamma, M, Re, Pr
      Vec :: numerical, solution
      Vec :: alpha_dudx, alpha_dudy, alpha_dudz
      Vec :: alpha_u, alpha_v, alpha_w, alpha_rho, alpha_T
      Vec :: alpha_rho_u, alpha_rho_v, alpha_rho_w
      Vec :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
      Vec :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
      Vec :: alpha_rho_uT, alpha_rho_vT, alpha_rho_wT
      type (pointer3D), dimension(3) :: u, v, w, rho, tem, one
      type (pointer3D), dimension(3) :: rho_u, rho_v, rho_w
      procedure(func3D), pointer :: sol_u, sol_v, sol_w, sol_rho, sol_tem
      procedure(func3D), pointer :: sol_rho_u, sol_rho_v, sol_rho_w

      mu = 1.0; gamma = 1.4; M = 1.0; Re = 1.0; Pr = 1.0

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2
      rho(1)%ptr => fcos2; rho(2)%ptr => fcos4; rho(3)%ptr => fcos2
      tem(1)%ptr => fcos4; tem(2)%ptr => fcos2; tem(3)%ptr => fcos4
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone
      rho_u(1)%ptr => fcos2cos4; rho_u(2)%ptr => fcos4sin4; rho_u(3)%ptr => fcos2cos4
      rho_v(1)%ptr => fcos2sin4; rho_v(2)%ptr => fcos4cos4; rho_v(3)%ptr => fcos2sin4
      rho_w(1)%ptr => fcos2cos2; rho_w(2)%ptr => fcos4sin2; rho_w(3)%ptr => fcos2cos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2
      sol_rho => gcos2cos4cos2_2
      sol_tem => gcos4cos2cos4_2
      sol_rho_u => grho_u
      sol_rho_v => grho_v
      sol_rho_w => grho_w

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_T, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uu, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_ww, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_wT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, numerical, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, solution, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing u \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing v \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing w \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      ! init rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rho \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%temp1, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! init T
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing T \n", ierr)
      call rhs_init_separable(ctx, grid, B, tem, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_T, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_tem, alpha_T, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      ! find rho u
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho u\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_u, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      ! find rho v
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho v\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_v, alpha_rho_v, ierr)
      CHKERRQ(ierr)

      ! find rho w
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho w\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_w, alpha_rho_w, ierr)
      CHKERRQ(ierr)

      ! find rho uu
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uu\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_u, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uu, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uv
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uv\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uv, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uw
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uw\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uw, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vv
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vv\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vv, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vw
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vw\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vw, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho ww
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho ww\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_w, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_ww, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uT, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vT, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho wT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho wT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_w, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_wT, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 1\n", ierr)

      call momentum(ctx, grid, B, 1, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, B%temp1, B%temp2, numerical, ierr)

      ! compute solution
      call compute_sol_momentum(ctx, grid, B, 1, &
         mu, gamma, M, Re, B%temp1, B%temp2, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 2\n", ierr)

      call momentum(ctx, grid, B, 2, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uv, alpha_rho_vv, alpha_rho_vw, B%temp1, B%temp2, numerical, ierr)

      ! compute solution
      call compute_sol_momentum(ctx, grid, B, 2, &
         mu, gamma, M, Re, B%temp1, B%temp2, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 3\n", ierr)

      call momentum(ctx, grid, B, 3, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uw, alpha_rho_vw, alpha_rho_ww, B%temp1, B%temp2, numerical, ierr)

      ! compute solution
      call compute_sol_momentum(ctx, grid, B, 3, &
         mu, gamma, M, Re, B%temp1, B%temp2, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Energy\n", ierr)

      call energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uT, alpha_rho_vT, alpha_rho_wT, B%temp1, B%temp2, B%temp3, numerical, ierr)

      ! compute solution
      call compute_sol_energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         B%temp1, B%temp2, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_momentum_energy

   subroutine test_serial_interpolate_1d(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      KSP :: ksp
      Vec :: localx, localy, localz
      PetscInt :: i, j, sx, sy, sz
      PetscReal :: x, numx
      PetscReal :: y, numy
      PetscReal :: z, numz
      PetscScalar, dimension(:), allocatable :: valx, valy, valz
      PetscScalar, dimension(:), pointer :: vvx, vvy, vvz
      procedure (func1D), pointer :: ptr

      ptr => fcos4

      ! B-spline init
      call test_bspline_1d(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      call VecDuplicate(B%rhsx, B%solx, ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAx, localx, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_x(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_x(ctx, grid, B, ptr, B%rhsx, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_x(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      ! interpolates in x
      x = 1.0 + grid%dx * 0.5

      allocate(valx(B%deg+1))

      sx = Int(x / grid%dx)

      do i=sx, sx+B%deg
         valx(i-sx+1) = xBJK(grid, B, i+1, B%deg, x)
      end do

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAx, B%alphax, INSERT_VALUES ,localx, ierr)
      call DMGlobalToLocalEnd(ctx%DAx, B%alphax, INSERT_VALUES ,localx, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAx, localx, vvx, ierr)
      CHKERRQ(ierr)

      numx = DOT_PRODUCT(valx, vvx(sx-B%deg:sx))

      print*, numx - COS(4.0 * x)

      call DMDAVecRestoreArrayF90(ctx%DAx, localx, vvx, ierr)
      CHKERRQ(ierr)

      !!! y dimension
      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      call VecDuplicate(B%rhsy, B%soly, ierr)
      call VecDuplicate(B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAy, localy, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call linear_y(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_y(ctx, grid, B, ptr, B%rhsy, ierr)

      ! create ksp
      call KSPCreate(ctx%commy, ksp, ierr)
      CHKERRQ(ierr)
      call KSPSetOperators(ksp, B%BBy, B%BBy, ierr)
      CHKERRQ(ierr)

      call KSPSetType(ksp, KSPGMRES, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0e-14, 1.0e-14, &
      PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)

      call KSPSetFromOptions(ksp, ierr)
      CHKERRQ(ierr)

      ! solve linear system
      call PetscPrintf(PETSC_COMM_WORLD, "Solve linear system \n", ierr)
      call KSPSolve(ksp, B%rhsy, B%alphay, ierr)
      CHKERRQ(ierr)

      ! interpolates in y
      y = 1.0 + grid%dy * 0.5

      allocate(valy(B%deg+1))

      sy = Int(y / grid%dy)

      do i=sy, sy+B%deg
         valy(i-sy+1) = yBJK(grid, B, i+1, B%deg, y)
      end do

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAy, B%alphay, INSERT_VALUES ,localy, ierr)
      call DMGlobalToLocalEnd(ctx%DAy, B%alphay, INSERT_VALUES ,localy, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAy, localy, vvy, ierr)
      CHKERRQ(ierr)

      numy = DOT_PRODUCT(valy, vvy(sy:sy+B%deg))

      print*, numy - COS(4.0 * y)

      call DMDAVecRestoreArrayF90(ctx%DAy, localy, vvy, ierr)
      CHKERRQ(ierr)

      !!! z dimension
      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      call VecDuplicate(B%rhsz, B%solz, ierr)
      call VecDuplicate(B%rhsz, B%alphaz, ierr)
      CHKERRQ(ierr)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAz, localz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_z(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_z(ctx, grid, B, ptr, B%rhsz, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_z(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_z(ctx, grid, B, ierr)

      ! interpolates in z
      z = 1.0 + grid%dz * 0.5

      allocate(valz(B%deg+1))

      sz = Int(z / grid%dz)

      do i=sz, sz+B%deg
         valz(i-sz+1) = zBJK(grid, B, i+1, B%deg, z)
      end do

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAz, B%alphaz, INSERT_VALUES ,localz, ierr)
      call DMGlobalToLocalEnd(ctx%DAz, B%alphaz, INSERT_VALUES ,localz, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAz, localz, vvz, ierr)
      CHKERRQ(ierr)

      numz = DOT_PRODUCT(valz, vvz(sz-B%deg:sz))

      print*, numz - COS(4.0 * z)

      call DMDAVecRestoreArrayF90(ctx%DAz, localz, vvz, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_x(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solx, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_serial_interpolate_1d

   subroutine test_serial_interpolate_3d(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscReal :: num
      PetscReal, dimension(3) :: x
      type(pointer3D), dimension(3) :: f
      PetscScalar, dimension(:, :, :), pointer :: vv
      procedure(func3D), pointer :: solptr

      f(1)%ptr => fcos4; f(2)%ptr => fsin4; f(3)%ptr => fcos4
      solptr => gcos4sin4cos4

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
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

      ! create local vectors
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, f, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, B%alpha, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, solptr, B%alpha, ierr)
      CHKERRQ(ierr)

      ! interpolate
      x(1) = 1.0 + grid%dx * 0.5
      x(2) = 1.0 + grid%dy * 0.5
      x(3) = 1.0 + grid%dz * 0.5

      ! global to local
      call DMGlobalToLocalBegin(ctx%DA, B%alpha, INSERT_VALUES, B%local1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, B%alpha, INSERT_VALUES, B%local1, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      call interpolate(ctx, grid, B, num, x, vv, ierr)

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      print*, num - solptr(x(1), x(2), x(3), 0.0)

      ! Destroy mat/vec
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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
   end subroutine test_serial_interpolate_3d

   subroutine test_reciprocal_1d(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      Vec :: localx
      PetscInt :: i, j, sx
      PetscReal :: x, numx
      PetscScalar, dimension(:), allocatable :: valx
      PetscScalar, dimension(:), pointer :: vvx
      procedure (func1D), pointer :: ptr, one, solptr, rsolptr

      ptr => fcos4; one => fone; solptr => fcos4_2; rsolptr => rcos4_2

      ! B-spline init
      call test_bspline_1d(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      call VecDuplicate(B%rhsx, B%solx, ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAx, localx, ierr)
      CHKERRQ(ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_x(ctx, grid, B, ierr)

      ! initialize rhs
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_x(ctx, grid, B, ptr, B%rhsx, ierr)
      call rhs_x(ctx, grid, B, one, B%solx, ierr)

      call VecAXPY(B%rhsx, complex(2.0, 0.0), B%solx, ierr)

      ! initialize FFT
      call PetscPrintf(PETSC_COMM_WORLD, "FFT init \n", ierr)
      call fft_init_x(ctx, grid, B, ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, solptr, B%alphax, ierr)
      CHKERRQ(ierr)

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAx, B%alphax, INSERT_VALUES ,localx, ierr)
      call DMGlobalToLocalEnd(ctx%DAx, B%alphax, INSERT_VALUES ,localx, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAx, localx, vvx, ierr)
      CHKERRQ(ierr)

      call rhs_reciprocal_x(ctx, grid, B, vvx, B%rhsx, ierr)

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DAx, localx, vvx, ierr)
      CHKERRQ(ierr)

      ! solve
      call PetscPrintf(PETSC_COMM_WORLD, "FFT solve \n", ierr)
      call fft_solve_x(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_x(ctx, grid, B, rsolptr, B%alphax, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call fft_destroy_x(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%solx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(localx, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_reciprocal_1d

   subroutine test_reciprocal_3d(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscReal :: num
      PetscReal, dimension(3) :: x
      type(pointer3D), dimension(3) :: f, one
      procedure(func3D), pointer :: solptr, rsolptr
      PetscScalar, dimension(:, :, :), pointer :: vv
      VecScatter :: scatter
      Vec :: seqvec

      f(1)%ptr => fsin4; f(2)%ptr => fcos4; f(3)%ptr => fsin4
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone
      solptr => gsin4cos4sin4_2
      rsolptr => rsin4cos4sin4_2

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
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

      ! create local vectors
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, f, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%sol, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call VecAXPY(B%rhs, complex(2.0, 0.0), B%sol, ierr)

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing mass matrix \n", ierr)
      call mass_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, B%alpha, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, solptr, B%alpha, ierr)
      CHKERRQ(ierr)

      ! global to local
      call DMGlobalToLocalBegin(ctx%DA, B%alpha, INSERT_VALUES, B%local1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, B%alpha, INSERT_VALUES, B%local1, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      call rhs_reciprocal_xyz(ctx, grid, B, vv, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, B%alpha, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, rsolptr, B%alpha, ierr)
      CHKERRQ(ierr)

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alpha, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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
   end subroutine test_reciprocal_3d

   subroutine test_nonlinear_combined_rhs(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val, mu, gamma, M, Re, Pr
      Vec :: numerical, solution
      Vec :: alpha_dudx, alpha_dudy, alpha_dudz
      Vec :: alpha_u, alpha_v, alpha_w, alpha_rho, alpha_T, alpha_invrho
      Vec :: alpha_rho_u, alpha_rho_v, alpha_rho_w
      Vec :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
      Vec :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
      Vec :: alpha_rho_uT, alpha_rho_vT, alpha_rho_wT
      type (pointer3D), dimension(3) :: u, v, w, rho, tem, one
      type (pointer3D), dimension(3) :: rho_u, rho_v, rho_w
      procedure(func3D), pointer :: sol_u, sol_v, sol_w, sol_rho, sol_tem
      procedure(func3D), pointer :: sol_rho_u, sol_rho_v, sol_rho_w, sol_invrho
      procedure(func3D), pointer :: sol_u_momentum, sol_v_momentum
      procedure(func3D), pointer :: sol_w_momentum, sol_energy
      PetscScalar, dimension(:, :, :), pointer :: vv

      mu = 1.0; gamma = 1.4; M = 1.0; Re = 1.0; Pr = 1.0

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2
      rho(1)%ptr => fcos2; rho(2)%ptr => fcos4; rho(3)%ptr => fcos2
      tem(1)%ptr => fcos4; tem(2)%ptr => fcos2; tem(3)%ptr => fcos4
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone
      rho_u(1)%ptr => fcos2cos4; rho_u(2)%ptr => fcos4sin4; rho_u(3)%ptr => fcos2cos4
      rho_v(1)%ptr => fcos2sin4; rho_v(2)%ptr => fcos4cos4; rho_v(3)%ptr => fcos2sin4
      rho_w(1)%ptr => fcos2cos2; rho_w(2)%ptr => fcos4sin2; rho_w(3)%ptr => fcos2cos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2
      sol_rho => gcos2cos4cos2_2
      sol_invrho => rcos2cos4cos2_2
      sol_tem => gcos4cos2cos4_2
      sol_rho_u => grho_u
      sol_rho_v => grho_v
      sol_rho_w => grho_w
      sol_u_momentum => u_momentum
      sol_v_momentum => v_momentum
      sol_w_momentum => w_momentum
      sol_energy => g_energy

      call bspline_init(ctx, grid, B, ierr)

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
      call VecDuplicate(B%rhs, B%sol, ierr)
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
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_T, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_invrho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uu, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_ww, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_wT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, numerical, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, solution, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing u \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing v \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing w \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      ! init rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rho \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%temp1, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! init T
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing T \n", ierr)
      call rhs_init_separable(ctx, grid, B, tem, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_T, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_tem, alpha_T, ierr)
      CHKERRQ(ierr)

      ! global to local
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing 1/rho\n", ierr)
      call DMGlobalToLocalBegin(ctx%DA, alpha_rho, INSERT_VALUES, B%local1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, alpha_rho, INSERT_VALUES, B%local1, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      call rhs_reciprocal_xyz(ctx, grid, B, vv, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_invrho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_invrho, alpha_invrho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho u/v/w\n", ierr)
      call nonlinear_three(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_u, B%temp1, B%local1, B%local2, &
         alpha_rho, alpha_v, B%temp2, B%local3, B%local4, &
         alpha_rho, alpha_w, B%temp3, B%local5, B%local6, ierr)

      call ksp_solve(ctx, grid, B, B%temp1, alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%temp2, alpha_rho_v, 'full', ierr)
      CHKERRQ(ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_v, alpha_rho_v, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%temp3, alpha_rho_w, 'full', ierr)
      CHKERRQ(ierr)
      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_w, alpha_rho_w, ierr)
      CHKERRQ(ierr)

      ! find rho uu
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums OTHERS", ierr)
      call nonlinear_nine(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_u, B%temp1, B%local1, B%local2, &
         alpha_rho_u, alpha_v, B%temp2, B%local3, B%local4, &
         alpha_rho_u, alpha_w, B%temp3, B%local5, B%local6, &
         alpha_rho_v, alpha_v, B%temp4, B%local7, B%local8, &
         alpha_rho_v, alpha_w, B%temp5, B%local9, B%local10, &
         alpha_rho_w, alpha_w, B%temp6, B%local11, B%local12, &
         alpha_rho_u, alpha_T, B%temp7, B%local13, B%local14, &
         alpha_rho_v, alpha_T, B%temp8, B%local15, B%local16, &
         alpha_rho_w, alpha_T, B%temp9, B%local17, B%local18, ierr)

      call ksp_solve(ctx, grid, B, B%temp1, alpha_rho_uu, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp2, alpha_rho_uv, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp3, alpha_rho_uw, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp4, alpha_rho_vv, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp5, alpha_rho_vw, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp6, alpha_rho_ww, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp7, alpha_rho_uT, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp8, alpha_rho_vT, 'full', ierr)
      CHKERRQ(ierr)
      call ksp_solve(ctx, grid, B, B%temp9, alpha_rho_wT, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 1\n", ierr)

      call rhs_momentum(ctx, grid, B, 1, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_u_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 2\n", ierr)

      call rhs_momentum(ctx, grid, B, 2, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uv, alpha_rho_vv, alpha_rho_vw, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_v_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 3\n", ierr)

      call rhs_momentum(ctx, grid, B, 3, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uw, alpha_rho_vw, alpha_rho_ww, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_w_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Energy\n", ierr)

      call rhs_energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uT, alpha_rho_vT, alpha_rho_wT, B%temp1, B%temp2, B%temp3, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_energy, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_nonlinear_combined_rhs

   subroutine test_rhs(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i
      PetscReal :: val, mu, gamma, M, Re, Pr
      Vec :: numerical, solution
      Vec :: alpha_dudx, alpha_dudy, alpha_dudz
      Vec :: alpha_u, alpha_v, alpha_w, alpha_rho, alpha_T, alpha_invrho
      Vec :: alpha_rho_u, alpha_rho_v, alpha_rho_w
      Vec :: alpha_rho_uu, alpha_rho_uv, alpha_rho_uw
      Vec :: alpha_rho_vv, alpha_rho_vw, alpha_rho_ww
      Vec :: alpha_rho_uT, alpha_rho_vT, alpha_rho_wT
      type (pointer3D), dimension(3) :: u, v, w, rho, tem, one
      type (pointer3D), dimension(3) :: rho_u, rho_v, rho_w
      procedure(func3D), pointer :: sol_u, sol_v, sol_w, sol_rho, sol_tem
      procedure(func3D), pointer :: sol_rho_u, sol_rho_v, sol_rho_w, sol_invrho
      procedure(func3D), pointer :: sol_u_momentum, sol_v_momentum
      procedure(func3D), pointer :: sol_w_momentum, sol_energy
      PetscScalar, dimension(:, :, :), pointer :: vv

      mu = 1.0; gamma = 1.4; M = 1.0; Re = 1.0; Pr = 1.0

      ! init function pointers
      u(1)%ptr => fcos4; u(2)%ptr => fsin4; u(3)%ptr => fcos4
      v(1)%ptr => fsin4; v(2)%ptr => fcos4; v(3)%ptr => fsin4
      w(1)%ptr => fcos2; w(2)%ptr => fsin2; w(3)%ptr => fcos2
      rho(1)%ptr => fcos2; rho(2)%ptr => fcos4; rho(3)%ptr => fcos2
      tem(1)%ptr => fcos4; tem(2)%ptr => fcos2; tem(3)%ptr => fcos4
      one(1)%ptr => fone; one(2)%ptr => fone; one(3)%ptr => fone
      rho_u(1)%ptr => fcos2cos4; rho_u(2)%ptr => fcos4sin4; rho_u(3)%ptr => fcos2cos4
      rho_v(1)%ptr => fcos2sin4; rho_v(2)%ptr => fcos4cos4; rho_v(3)%ptr => fcos2sin4
      rho_w(1)%ptr => fcos2cos2; rho_w(2)%ptr => fcos4sin2; rho_w(3)%ptr => fcos2cos2

      sol_u => gcos4sin4cos4
      sol_v => gsin4cos4sin4
      sol_w => gcos2sin2cos2
      sol_rho => gcos2cos4cos2_2
      sol_invrho => rcos2cos4cos2_2
      sol_tem => gcos4cos2cos4_2
      sol_rho_u => grho_u
      sol_rho_v => grho_v
      sol_rho_w => grho_w
      sol_u_momentum => u_momentum
      sol_v_momentum => v_momentum
      sol_w_momentum => w_momentum
      sol_energy => g_energy

      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DA, B%local1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%local1, B%local2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%sol, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp1, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp2, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp3, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_T, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_invrho, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_u, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_v, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_w, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uu, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vv, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vw, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_ww, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_uT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_vT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, alpha_rho_wT, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, numerical, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, solution, ierr)
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

      ! initialize mass matrix
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing linear operators \n", ierr)
      call linear_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! initialize ksp
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing KSP \n", ierr)
      call ksp_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! init u
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing u \n", ierr)
      call rhs_init_separable(ctx, grid, B, u, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_u, alpha_u, ierr)
      CHKERRQ(ierr)

      ! init v
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing v \n", ierr)
      call rhs_init_separable(ctx, grid, B, v, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_v, alpha_v, ierr)
      CHKERRQ(ierr)

      ! init w
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing w \n", ierr)
      call rhs_init_separable(ctx, grid, B, w, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_w, alpha_w, ierr)
      CHKERRQ(ierr)

      ! init rho
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rho \n", ierr)
      call rhs_init_separable(ctx, grid, B, rho, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call rhs_init_separable(ctx, grid, B, one, B%temp1, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho, alpha_rho, ierr)
      CHKERRQ(ierr)

      ! init T
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing T \n", ierr)
      call rhs_init_separable(ctx, grid, B, tem, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecAXPY(B%rhs, complex(2.0, 0.0), B%temp1, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_T, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_tem, alpha_T, ierr)
      CHKERRQ(ierr)

      ! global to local
      call PetscPrintf(PETSC_COMM_WORLD, "Initializing 1/rho\n", ierr)
      call DMGlobalToLocalBegin(ctx%DA, alpha_rho, INSERT_VALUES, B%local1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, alpha_rho, INSERT_VALUES, B%local1, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, B%local1, vv, ierr)
      CHKERRQ(ierr)

      call rhs_reciprocal_xyz(ctx, grid, B, vv, B%rhs, ierr)

      call ksp_solve(ctx, grid, B, B%rhs, alpha_invrho, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_invrho, alpha_invrho, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing nonlinear operators \n", ierr)
      call nonlinear_init(ctx, grid, B, ierr)

      ! find rho u
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho u\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_u, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_u, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_u, alpha_rho_u, ierr)
      CHKERRQ(ierr)

      ! find rho v
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho v\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_v, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_v, alpha_rho_v, ierr)
      CHKERRQ(ierr)

      ! find rho w
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho w\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_w, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Calculate error \n", ierr)
      call Linf_eval_xyz(ctx, grid, B, sol_rho_w, alpha_rho_w, ierr)
      CHKERRQ(ierr)

      ! find rho uu
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uu\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_u, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uu, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uv
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uv\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uv, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uw
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uw\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uw, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vv
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vv\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_v, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vv, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vw
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vw\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vw, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho ww
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho ww\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_w, alpha_w, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_ww, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho uT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho uT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_u, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_uT, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho vT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho vT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_v, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_vT, 'full', ierr)
      CHKERRQ(ierr)

      ! find rho wT
      call PetscPrintf(PETSC_COMM_WORLD, "Nonlinear sums rho wT\n", ierr)
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
         alpha_rho_w, alpha_T, B%rhs, B%local1, B%local2, ierr)
      call ksp_solve(ctx, grid, B, B%rhs, alpha_rho_wT, 'full', ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 1\n", ierr)

      call rhs_momentum(ctx, grid, B, 1, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uu, alpha_rho_uv, alpha_rho_uw, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_u_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 2\n", ierr)

      call rhs_momentum(ctx, grid, B, 2, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uv, alpha_rho_vv, alpha_rho_vw, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_v_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Momentum, i = 3\n", ierr)

      call rhs_momentum(ctx, grid, B, 3, mu, gamma, M, Re, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uw, alpha_rho_vw, alpha_rho_ww, B%temp1, B%temp2, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_w_momentum, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      call PetscPrintf(PETSC_COMM_WORLD, "Energy\n", ierr)

      call rhs_energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         alpha_rho, alpha_u, alpha_v, alpha_w, alpha_T, alpha_invrho, &
         alpha_rho_u, alpha_rho_v, alpha_rho_w, &
         alpha_rho_uT, alpha_rho_vT, alpha_rho_wT, B%temp1, B%temp2, B%temp3, numerical, ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Compute solution\n", ierr)
      call rhs_xyz(ctx, grid, B, sol_energy, solution, ierr)

      ! compute error
      call VecAXPY(solution, complex(-1.0, 0.0), numerical, ierr)
      call VecNorm(solution, NORM_INFINITY, val, ierr)

      print*, val

      ! destroy
      call ksp_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%sol, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%local2, ierr)
      CHKERRQ(ierr)

      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%tempx, ierr)
      CHKERRQ(ierr)

      call MatDestroy(B%My, ierr)
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

      if (allocated(B%BBBx)) deallocate(B%BBBx)
      if (allocated(B%BBDx)) deallocate(B%BBDx)
      if (allocated(B%BBBy)) deallocate(B%BBBy)
      if (allocated(B%BBDy)) deallocate(B%BBDy)
      if (allocated(B%BBBz)) deallocate(B%BBBz)
      if (allocated(B%BBDz)) deallocate(B%BBDz)
   end subroutine test_rhs

   subroutine test_petsc_write(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      MPI_Comm :: comm
      PetscViewer :: view
      type(pointer3D), dimension(3) :: f

      f(1)%ptr => f1d; f(2)%ptr => f1d; f(3)%ptr => f1d

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      CHKERRQ(ierr)

      call PetscPrintf(PETSC_COMM_WORLD, "Initializing rhs \n", ierr)
      call rhs_init_separable(ctx, grid, B, f, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      call PetscViewerCreate(PETSC_COMM_WORLD, view, ierr)
      CHKERRA(ierr)
      call PetscViewerSetType(view, PETSCVIEWERVTK, ierr)
      CHKERRA(ierr)
      call PetscViewerFileSetName(view, "sol.vts", ierr)
      CHKERRA(ierr)
      call VecView(B%rhs, view, ierr)
      CHKERRQ(ierr)
      call PetscViewerDestroy(view, ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_petsc_write

   subroutine test_fortran_write(ctx, grid, B, ierr)
      use output
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: time
      type(pointer3D), dimension(3) :: f

      f(1)%ptr => f1d; f(2)%ptr => f1d; f(3)%ptr => f1d

      ! B-spline init
      call bspline_init(ctx, grid, B, ierr)

      ! create vectors (solution and rhs)
      call DMCreateGlobalVector(ctx%DA, B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhs, B%temp, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAx, B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%alphax, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsx, B%tempx, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAy, B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsy, B%tempy, ierr)
      CHKERRQ(ierr)
      call DMCreateGlobalVector(ctx%DAz, B%rhsz, ierr)
      CHKERRQ(ierr)
      call VecDuplicate(B%rhsz, B%tempz, ierr)
      CHKERRQ(ierr)

      call rhs_init_separable(ctx, grid, B, f, B%rhs, B%rhsx, B%rhsy, B%rhsz, ierr)
      CHKERRQ(ierr)

      ! do this only once
      call io_init(ctx, grid, B, ierr)
      CHKERRQ(ierr)

      ! do this everytime you scatter
      call io_fortran_scatter(ctx, B, B%rhs, ierr)
      CHKERRQ(ierr)

      if (ctx%procx .eq. 1 .and. ctx%procy .eq. 1 .and. ctx%procz .eq. 1) then
         time = 0
      else
         time = 5
      end if

      call io_fortran_write(ctx, grid, B, time, "rho", ierr)
      CHKERRQ(ierr)

      ! Destroy mat/vec
      call io_destroy(B, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhs, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%temp, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsy, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%rhsz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine test_fortran_write

end module test

program test_operators
   use types
   use lattice
   use bspline
   use test
   use environment
   implicit none

   type(ctx_t) :: ctx
   type(grid_t) :: grid
   type(bspline_t) :: B
   type(rk_t) :: rk
   PetscErrorCode :: ierr
   character(len=64) :: name
   PetscBool :: found

   call environment_init(ctx, ierr)
   call lattice_init(ctx, grid, B, ierr)
   call grid_init(ctx, grid, ierr)

   ! select the check at run time: -test <name> (default: nonlinear_combined_rhs)
   name = 'nonlinear_combined_rhs'
   call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                              '-test', name, found, ierr)

   call run_test(name, ctx, grid, B, rk, ierr)

   call lattice_free(ctx, ierr)
   call environment_free(ierr)
end program

! vim: set ts=3 sw=3 sts=3:
