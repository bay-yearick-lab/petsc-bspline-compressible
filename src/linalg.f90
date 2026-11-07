! Tensor-product mass solves (FFTW in x and z, PETSc KSP in y) and Navier-Stokes operator kernels.
# include <petsc/finclude/petscsys.h>
# include <petsc/finclude/petscksp.h>

module linalg
   use types
   implicit none

   private

   public :: matvec
   public :: set_ybc
   public :: channel_flow_bc
   public :: ksp_init
   public :: ksp_solve
   public :: ksp_destroy
   public :: fft_init_x
   public :: fft_solve_x
   public :: fft_destroy_x
   public :: fft_init_z
   public :: fft_solve_z
   public :: fft_destroy_z
   public :: nonlinear_sum_BBB
   public :: nonlinear_three
   public :: nonlinear_nine
   public :: nonlinear_sum_BBBB
   public :: nonlinear_sum_x
   public :: nonlinear_sum_y
   public :: nonlinear_sum_z
   public :: divergence
   public :: convection
   public :: pressure
   public :: viscous
   public :: energy_viscous
   public :: momentum
   public :: energy_laplacian
   public :: energy_divergence
   public :: energy
   public :: d1dx
   public :: d1dy
   public :: d1dz
   public :: d2dx
   public :: d2dy
   public :: d2dz
   public :: dxdy
   public :: dydz
   public :: dxdz
   public :: interpolate
   public :: bodyforce
contains

   subroutine bodyforce(ctx, grid, B, mu, Re, u, v, temp, vector, ierr)
      use bspline
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: mu, Re
      Vec :: u, v, temp, vector
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscReal :: aveup, avedown, mean
      PetscReal :: localup, localdown, globalup, globaldown
      PetscScalar, dimension(:, :, :), pointer :: vv

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      call d1dx(ctx, grid, B, v, vector, ierr)
      call d1dy(ctx, grid, B, u, temp, ierr)

      call VecAXPY(temp, complex(1.0, 0.0), vector, ierr)
      CHKERRQ(ierr)
      call VecScale(temp, complex(mu, 0.0), ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, temp, vv, ierr)
      CHKERRQ(ierr)

      localup = 0.0
      localdown = 0.0
      do k=ctx%zs-1, ctx%ze-1
         do i=ctx%xs-1, ctx%xe-1
            if (ctx%ys-1 .eq. 0) then
               j = 0

               localup = localup + vv(i, j, k)
            else if (ctx%ye-1 .eq. ctx%my-1) then
               j = ctx%my-1

               localdown = localdown + vv(i, j, k)
            end if
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, temp, vv, ierr)
      CHKERRQ(ierr)

      call MPI_Allreduce(localup, globalup, 1, MPI_FLOAT, MPI_SUM, &
                    PETSC_COMM_WORLD, ierr)
      call MPI_Allreduce(localdown, globaldown, 1, MPI_FLOAT, MPI_SUM, &
                    PETSC_COMM_WORLD, ierr)

      aveup = globalup / (grid%xNc * grid%zNc)
      avedown = globaldown / (grid%xNc * grid%zNc)

      mean = - (aveup - avedown) / (2.0 * Re)

      call VecSet(vector, complex(mean, 0.0), ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine bodyforce

   subroutine interpolate(ctx, grid, B, num, x, vv, ierr)
      use bspline
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: num
      PetscReal, dimension(3) :: x
      PetscScalar, dimension(:, :, :), pointer :: vv
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscInt :: sx, sy, sz

      if (.not. allocated(B%intx)) allocate(B%intx(B%deg+1))
      if (.not. allocated(B%inty)) allocate(B%inty(B%deg+1))
      if (.not. allocated(B%intz)) allocate(B%intz(B%deg+1))

      sx = Int(x(1) / grid%dx)
      sy = Int(x(2) / grid%dy)
      sz = Int(x(3) / grid%dz)

      do i=sx, sx+B%deg
         B%intx(i-sx+1) = xBJK(grid, B, i+1, B%deg, x(1))
      end do

      do i=sy, sy+B%deg
         B%inty(i-sy+1) = yBJK(grid, B, i+1, B%deg, x(2))
      end do

      do i=sz, sz+B%deg
         B%intz(i-sz+1) = zBJK(grid, B, i+1, B%deg, x(3))
      end do

      num = 0.0
      do k=sz-B%deg, sz
         do j=sy,sy+B%deg
            do i=sx-B%deg, sx
               num = num + B%intx(i-(sx-B%deg)+1) * &
                           B%inty(j-sy+1) * &
                           B%intz(k-(sz-B%deg)+1) * &
                           vv(i, j, k)
            end do
         end do
      end do

      ierr = 0
      return
   end subroutine interpolate

   subroutine channel_flow_bc(ctx, grid, B, u, v, w, T, Tbc, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: u, v, w, T
      PetscReal :: Tbc
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:, :, :), pointer :: vu, vv, vw, vT

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, u, vu, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, v, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, w, vw, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, T, vT, ierr)
      CHKERRQ(ierr)

      ! reset the ybc
      do k=ctx%zs-1, ctx%ze-1
         do j=ctx%ys-1, ctx%ye-1
            do i=ctx%xs-1, ctx%xe-1
               if (j .eq. 0) then
                  vu(i, j, k) = 0.0
                  vv(i, j, k) = 0.0
                  vw(i, j, k) = 0.0
                  vT(i, j, k) = Tbc
               else if (j .eq. ctx%my-1) then
                  vu(i, j, k) = 0.0
                  vv(i, j, k) = 0.0
                  vw(i, j, k) = 0.0
                  vT(i, j, k) = Tbc
               end if
            end do
         end do
      end do

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, u, vu, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, v, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, w, vw, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, T, vT, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine channel_flow_bc

   subroutine set_ybc(ctx, grid, B, vec, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: vec
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:, :, :), pointer :: vvec

      ! get array
      call DMDAVecGetArrayF90(ctx%DA, vec, vvec, ierr)
      CHKERRQ(ierr)

      ! reset the ybc
      do k=ctx%zs-1, ctx%ze-1
         do i=ctx%xs-1, ctx%xe-1
            if (ctx%ys-1 .eq. 0) then
               j = 0
               vvec(i, j, k) = vvec(i, j, k) / B%old_ybc(1) * B%new_ybc(1)
            else if (ctx%ye-1 .eq. ctx%my-1) then
               j = ctx%my-1
               vvec(i, j, k) = vvec(i, j, k) / B%old_ybc(2) * B%new_ybc(2)
            end if
         end do
      end do

      B%old_ybc(:) = B%new_ybc(:)

      ! restore array
      call DMDAVecRestoreArrayF90(ctx%DA, vec, vvec, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine set_ybc

   subroutine fft_solve_z(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt, dimension(1) :: FFTvec

      associate (FFTmat => B%FFTmatz, input => B%inputz, output => B%outputz, backward => B%backwardz, &
         FFTmat1 => B%FFTmatz1, input1 => B%inputz1, output1 => B%outputz1, backward1 => B%backwardz1)

      ! FFT rhs
      call VecScatterPetscToFFTW(FFTmat1, B%rhsz, input1, ierr)
      CHKERRQ(ierr)
      call MatMult(FFTmat1, input1, output1, ierr)
      CHKERRQ(ierr)

      ! FFT column of periodic matrix
      call VecScatterPetscToFFTW(FFTmat, B%FFTz, input, ierr)
      CHKERRQ(ierr)
      call MatMult(FFTmat, input, output, ierr)
      CHKERRQ(ierr)

      ! divide in FFT space
      call VecPointwiseDivide(output, output1, output, ierr)

      ! inverse FFT
      call MatMultTranspose(FFTmat, output, backward, ierr)
      CHKERRQ(ierr)
      call VecScatterFFTWToPetsc(FFTmat, backward, B%tempz, ierr)
      CHKERRQ(ierr)

      ! scale solution
      call VecScale(B%tempz, 1.0/grid%zNc + 0*PETSC_i, ierr)

      ! copy into answer
      call VecCopy(B%tempz, B%alphaz, ierr)

      end associate
      ierr = 0
      return
   end subroutine fft_solve_z

   subroutine fft_init_z(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt, dimension(1) :: FFTvec

      associate (FFTmat => B%FFTmatz, input => B%inputz, output => B%outputz, backward => B%backwardz, &
         FFTmat1 => B%FFTmatz1, input1 => B%inputz1, output1 => B%outputz1, backward1 => B%backwardz1)

      FFTvec(1) = grid%zNc

      call MatCreateFFT(ctx%commz, 1, FFTvec, MATFFTW, FFTmat, ierr)
      CHKERRQ(ierr)
      call MatCreateFFT(ctx%commz, 1, FFTvec, MATFFTW, FFTmat1, ierr)
      CHKERRQ(ierr)
      call MatCreateVecsFFTW(FFTmat, input, output, backward, ierr)
      CHKERRQ(ierr)
      call MatCreateVecsFFTW(FFTmat1, input1, output1, backward1, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine fft_init_z

   subroutine fft_destroy_z(B, ierr)
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      ! Destroy mat/vec
      call MatDestroy(B%FFTmatz, ierr)
      CHKERRQ(ierr)
      call MatDestroy(B%FFTmatz1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%outputz1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%inputz1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%backwardz1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%inputz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%outputz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%backwardz, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%FFTz, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine fft_destroy_z

   subroutine fft_solve_x(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt, dimension(1) :: FFTvec

      associate (FFTmat => B%FFTmatx, input => B%inputx, output => B%outputx, backward => B%backwardx, &
         FFTmat1 => B%FFTmatx1, input1 => B%inputx1, output1 => B%outputx1, backward1 => B%backwardx1)

      ! FFT rhs
      call VecScatterPetscToFFTW(FFTmat1, B%rhsx, input1, ierr)
      CHKERRQ(ierr)
      call MatMult(FFTmat1, input1, output1, ierr)
      CHKERRQ(ierr)

      ! FFT column of periodic matrix
      call VecScatterPetscToFFTW(FFTmat, B%FFTx, input, ierr)
      CHKERRQ(ierr)
      call MatMult(FFTmat, input, output, ierr)
      CHKERRQ(ierr)

      ! divide in FFT space
      call VecPointwiseDivide(output, output1, output, ierr)

      ! inverse FFT
      call MatMultTranspose(FFTmat, output, backward, ierr)
      CHKERRQ(ierr)
      call VecScatterFFTWToPetsc(FFTmat, backward, B%tempx, ierr)
      CHKERRQ(ierr)

      ! scale solution
      call VecScale(B%tempx, 1.0/grid%xNc + 0*PETSC_i, ierr)

      ! copy into answer
      call VecCopy(B%tempx, B%alphax, ierr)

      end associate
      ierr = 0
      return
   end subroutine fft_solve_x

   subroutine fft_init_x(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      PetscInt, dimension(1) :: FFTvec

      associate (FFTmat => B%FFTmatx, input => B%inputx, output => B%outputx, backward => B%backwardx, &
         FFTmat1 => B%FFTmatx1, input1 => B%inputx1, output1 => B%outputx1, backward1 => B%backwardx1)

      FFTvec(1) = grid%xNc

      call MatCreateFFT(ctx%commx, 1, FFTvec, MATFFTW, FFTmat, ierr)
      CHKERRQ(ierr)
      call MatCreateFFT(ctx%commx, 1, FFTvec, MATFFTW, FFTmat1, ierr)
      CHKERRQ(ierr)
      call MatCreateVecsFFTW(FFTmat, input, output, backward, ierr)
      CHKERRQ(ierr)
      call MatCreateVecsFFTW(FFTmat1, input1, output1, backward1, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine fft_init_x

   subroutine fft_destroy_x(B, ierr)
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      ! Destroy mat/vec
      call MatDestroy(B%FFTmatx, ierr)
      CHKERRQ(ierr)
      call MatDestroy(B%FFTmatx1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%outputx1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%inputx1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%backwardx1, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%inputx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%outputx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%backwardx, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%FFTx, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine fft_destroy_x

   subroutine ksp_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      associate (ksp => B%ksp, fullksp => B%fullksp, pc => B%pc)

      ! y dimension create ksp
      call KSPCreate(ctx%commy, ksp, ierr)
      CHKERRQ(ierr)
      call KSPSetOperators(ksp, B%My, B%My, ierr)
      CHKERRQ(ierr)

      ! Krylov method and preconditioner are runtime options (KSPSetFromOptions),
      ! e.g. -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist

      call KSPSetType(ksp, KSPGMRES, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(ksp, 1.0e-14, 1.0e-14, &
      PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call KSPSetFromOptions(ksp, ierr)
      CHKERRQ(ierr)

      ! y dimension create full ksp
      call KSPCreate(ctx%commy, fullksp, ierr)
      CHKERRQ(ierr)
      call KSPSetOperators(fullksp, B%BBy, B%BBy, ierr)
      CHKERRQ(ierr)

      call KSPSetType(fullksp, KSPCG, ierr)
      CHKERRQ(ierr)
      call KSPSetTolerances(fullksp, 1.0e-14, 1.0e-14, &
      PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
      call KSPSetFromOptions(fullksp, ierr)
      CHKERRQ(ierr)

      ! fft in x and z dimensions
      call fft_init_x(ctx, grid, B, ierr)
      call fft_init_z(ctx, grid, B, ierr)

      end associate
      ierr = 0
      return
   end subroutine ksp_init

   subroutine ksp_solve(ctx, grid, B, rhs3d, x, solver_type, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: rhs3d, x
      character(len=*) :: solver_type
      PetscErrorCode :: ierr

      KSP :: ksp
      PetscInt :: i, j, k
      PetscScalar, dimension(:), pointer :: rhsx, rhsy, rhsz
      PetscScalar, dimension(:), pointer :: alphax, alphay, alphaz
      PetscScalar, dimension(:, :, :), pointer :: rhs

      if (solver_type == 'bc') then
         ksp = B%ksp
      else if (solver_type == 'full') then
         ksp = B%fullksp
      else
         call PetscPrintf(PETSC_COMM_WORLD, 'Erroneous solver_type in ksp_solve routine!' , ierr)
      end if

      ! copy rhs to solution to solve in place and not ruin rhs
      call VecCopy(rhs3d, x, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, x, rhs, ierr)
      CHKERRQ(ierr)

      ! solve x linear system
      do k=ctx%zs-1, ctx%ze-1
         do j=ctx%ys-1, ctx%ye-1

            call DMDAVecGetArrayF90(ctx%DAx, B%rhsx, rhsx, ierr)
            CHKERRQ(ierr)

            rhsx = rhs(:, j, k)
            call DMDAVecRestoreArrayF90(ctx%DAx, B%rhsx, rhsx, ierr)
            CHKERRQ(ierr)

            call fft_solve_x(ctx, grid, B, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAx, B%alphax, alphax, ierr)
            CHKERRQ(ierr)

            rhs(:, j, k) = alphax

            call DMDAVecRestoreArrayF90(ctx%DAx, B%alphax, alphax, ierr)
            CHKERRQ(ierr)
         end do
      end do

      ! solve y linear system
      do k=ctx%zs-1, ctx%ze-1
         do i=ctx%xs-1, ctx%xe-1

            call DMDAVecGetArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)
            CHKERRQ(ierr)

            rhsy = rhs(i, :, k)
            call DMDAVecRestoreArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)
            CHKERRQ(ierr)

            call KSPSolve(ksp, B%rhsy, B%alphay, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAy, B%alphay, alphay, ierr)
            CHKERRQ(ierr)

            rhs(i, :, k) = alphay

            call DMDAVecRestoreArrayF90(ctx%DAy, B%alphay, alphay, ierr)
            CHKERRQ(ierr)
         end do
      end do

      ! solve z linear system
      do j=ctx%ys-1, ctx%ye-1
         do i=ctx%xs-1, ctx%xe-1

            call DMDAVecGetArrayF90(ctx%DAz, B%rhsz, rhsz, ierr)
            CHKERRQ(ierr)

            rhsz = rhs(i, j, :)
            call DMDAVecRestoreArrayF90(ctx%DAz, B%rhsz, rhsz, ierr)
            CHKERRQ(ierr)

            call fft_solve_z(ctx, grid, B, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAz, B%alphaz, alphaz, ierr)
            CHKERRQ(ierr)

            rhs(i, j, :) = alphaz

            call DMDAVecRestoreArrayF90(ctx%DAz, B%alphaz, alphaz, ierr)
            CHKERRQ(ierr)
         end do
      end do

      call DMDAVecRestoreArrayF90(ctx%DA, x, rhs, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine ksp_solve

   subroutine ksp_destroy(B, ierr)
      type(bspline_t) :: B
      PetscErrorCode :: ierr

      ! release the wall-normal KSP and the FFTW plans
      call KSPDestroy(B%ksp, ierr)
      call fft_destroy_x(B, ierr)
      call fft_destroy_z(B, ierr)

      ierr = 0
      return
   end subroutine ksp_destroy

   subroutine matvec(ctx, B, X, Y, Z, v, v2, ierr)
      type(ctx_t) :: ctx
      type(bspline_t) :: B
      Mat :: X, Y, Z
      Vec :: v, v2
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      PetscScalar, dimension(:), pointer :: rhsx, rhsy, rhsz
      PetscScalar, dimension(:), pointer :: tempx, tempy, tempz
      PetscScalar, dimension(:, :, :), pointer :: vv

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye, &
         zs => ctx%zs, ze => ctx%ze)

      call VecCopy(v, v2, ierr)
      CHKERRQ(ierr)

      ! x - dimension
      do k=zs-1, ze-1
         do j=ys-1, ye-1

            call DMDAVecGetArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
            call DMDAVecGetArrayF90(ctx%DAx, B%rhsx, rhsx, ierr)
            CHKERRQ(ierr)

            rhsx = vv(:, j, k)
            call DMDAVecRestoreArrayF90(ctx%DAx, B%rhsx, rhsx, ierr)
            CHKERRQ(ierr)

            call MatMult(X, B%rhsx, B%tempx, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAx, B%tempx, tempx, ierr)
            CHKERRQ(ierr)

            vv(:, j, k) = tempx

            call DMDAVecRestoreArrayF90(ctx%DAx, B%tempx, tempx, ierr)
            CHKERRQ(ierr)
            call DMDAVecRestoreArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
         end do
      end do

      ! y - dimension
      do k=zs-1, ze-1
         do i=xs-1, xe-1

            call DMDAVecGetArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
            call DMDAVecGetArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)
            CHKERRQ(ierr)

            rhsy = vv(i, :, k)
            call DMDAVecRestoreArrayF90(ctx%DAy, B%rhsy, rhsy, ierr)
            CHKERRQ(ierr)

            call MatMult(Y, B%rhsy, B%tempy, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAy, B%tempy, tempy, ierr)
            CHKERRQ(ierr)

            vv(i, :, k) = tempy

            call DMDAVecRestoreArrayF90(ctx%DAy, B%tempy, tempy, ierr)
            CHKERRQ(ierr)
            call DMDAVecRestoreArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
         end do
      end do

      ! z - dimension
      do j=ys-1, ye-1
         do i=xs-1, xe-1

            call DMDAVecGetArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
            call DMDAVecGetArrayF90(ctx%DAz, B%rhsz, rhsz, ierr)
            CHKERRQ(ierr)

            rhsz = vv(i, j, :)
            call DMDAVecRestoreArrayF90(ctx%DAz, B%rhsz, rhsz, ierr)
            CHKERRQ(ierr)

            call MatMult(Z, B%rhsz, B%tempz, ierr)
            CHKERRQ(ierr)

            call DMDAVecGetArrayF90(ctx%DAz, B%tempz, tempz, ierr)
            CHKERRQ(ierr)

            vv(i, j, :) = tempz

            call DMDAVecRestoreArrayF90(ctx%DAz, B%tempz, tempz, ierr)
            CHKERRQ(ierr)
            call DMDAVecRestoreArrayF90(ctx%DA, v2, vv, ierr)
            CHKERRQ(ierr)
         end do
      end do

      end associate
      ierr = 0
      return
   end subroutine matvec

   subroutine nonlinear_sum_z(ctx, grid, B, v, v2, v3, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: v, v2, v3
      PetscErrorCode :: ierr

      PetscInt :: r, s, t
      Vec :: local, local2
      PetscScalar, dimension(:), pointer :: vv, vv2, vv3

      associate (zs => ctx%zs, ze => ctx%ze, d => B%deg)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAz, local, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DAz, local2, ierr)
      CHKERRQ(ierr)

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAz, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalEnd(ctx%DAz, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalBegin(ctx%DAz, v2, INSERT_VALUES ,local2, ierr)
      call DMGlobalToLocalEnd(ctx%DAz, v2, INSERT_VALUES ,local2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAz, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAz, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAz, v3, vv3, ierr)
      CHKERRQ(ierr)

      vv3 = 0.0

      do t=-d, d
         do s=-d, d
            do r=zs-1, ze-1
               vv3(r) = vv3(r) + B%BBBz(r, s, t) * vv(r+s) * vv2(r+t)
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DAz, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAz, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAz, v3, vv3, ierr)
      CHKERRQ(ierr)

      call VecDestroy(local, ierr)
      call VecDestroy(local2, ierr)
      end associate
      ierr = 0
      return
   end subroutine nonlinear_sum_z

   subroutine nonlinear_sum_y(ctx, grid, B, v, v2, v3, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: v, v2, v3
      PetscErrorCode :: ierr

      PetscInt :: l, m, n
      Vec :: local, local2
      PetscScalar, dimension(:), pointer :: vv, vv2, vv3

      associate (ys => ctx%ys, ye => ctx%ye, d => B%deg)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAy, local, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DAy, local2, ierr)
      CHKERRQ(ierr)

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAy, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalEnd(ctx%DAy, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalBegin(ctx%DAy, v2, INSERT_VALUES ,local2, ierr)
      call DMGlobalToLocalEnd(ctx%DAy, v2, INSERT_VALUES ,local2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAy, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAy, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAy, v3, vv3, ierr)
      CHKERRQ(ierr)

      vv3 = 0.0

      do n=-d, d
         do m=-d, d
            do l=ys-1, ye-1
               vv3(l) = vv3(l) + B%BBBy(l, m, n) * vv(l+n) * vv2(l+m)
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DAy, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAy, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAy, v3, vv3, ierr)
      CHKERRQ(ierr)

      call VecDestroy(local, ierr)
      call VecDestroy(local2, ierr)
      end associate
      ierr = 0
      return
   end subroutine nonlinear_sum_y

   subroutine nonlinear_sum_x(ctx, grid, B, v, v2, v3, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: v, v2, v3
      PetscErrorCode :: ierr

      PetscInt :: i, j, k
      Vec :: local, local2
      PetscScalar, dimension(:), pointer :: vv, vv2, vv3

      associate (xs => ctx%xs, xe => ctx%xe, d => B%deg)

      ! create local vectors
      call DMCreateLocalVector(ctx%DAx, local, ierr)
      CHKERRQ(ierr)
      call DMCreateLocalVector(ctx%DAx, local2, ierr)
      CHKERRQ(ierr)

      ! global to local
      call DMGlobalToLocalBegin(ctx%DAx, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalEnd(ctx%DAx, v, INSERT_VALUES ,local, ierr)
      call DMGlobalToLocalBegin(ctx%DAx, v2, INSERT_VALUES ,local2, ierr)
      call DMGlobalToLocalEnd(ctx%DAx, v2, INSERT_VALUES ,local2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DAx, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAx, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DAx, v3, vv3, ierr)
      CHKERRQ(ierr)

      vv3 = 0.0

      do k=-d, d
         do j=-d, d
            do i=xs-1, xe-1
               vv3(i) = vv3(i) + B%BBBx(i, j, k) * vv(i+k) * vv2(i+j)
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DAx, local, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAx, local2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DAx, v3, vv3, ierr)
      CHKERRQ(ierr)

      call VecDestroy(local, ierr)
      call VecDestroy(local2, ierr)
      end associate
      ierr = 0
      return
   end subroutine nonlinear_sum_x

   subroutine nonlinear_nine(ctx, grid, B, X, Y, Z, &
        a1, a2, a3, al1, al2, e1, e2, e3, el1, el2, &
        b1, b2, b3, bl1, bl2, f1, f2, f3, fl1, fl2, &
        c1, c2, c3, cl1, cl2, g1, g2, g3, gl1, gl2, &
        d1, d2, d3, dl1, dl2, h1, h2, h3, hl1, hl2, &
        i1, i2, i3, il1, il2, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal, dimension(ctx%xs-1:ctx%xe-1, -B%deg:B%deg, -B%deg:B%deg) :: X
      PetscReal, dimension(ctx%ys-1:ctx%ye-1, -B%deg:B%deg, -B%deg:B%deg) :: Y
      PetscReal, dimension(ctx%zs-1:ctx%ze-1, -B%deg:B%deg, -B%deg:B%deg) :: Z
      Vec :: a1, a2, a3, al1, al2
      Vec :: b1, b2, b3, bl1, bl2
      Vec :: c1, c2, c3, cl1, cl2
      Vec :: d1, d2, d3, dl1, dl2
      Vec :: e1, e2, e3, el1, el2
      Vec :: f1, f2, f3, fl1, fl2
      Vec :: g1, g2, g3, gl1, gl2
      Vec :: h1, h2, h3, hl1, hl2
      Vec :: i1, i2, i3, il1, il2
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m, n, r, s, t
      PetscScalar, dimension(:, :, :), pointer :: av, av2, av3
      PetscScalar, dimension(:, :, :), pointer :: bv, bv2, bv3
      PetscScalar, dimension(:, :, :), pointer :: cv, cv2, cv3
      PetscScalar, dimension(:, :, :), pointer :: dv, dv2, dv3
      PetscScalar, dimension(:, :, :), pointer :: ev, ev2, ev3
      PetscScalar, dimension(:, :, :), pointer :: fv, fv2, fv3
      PetscScalar, dimension(:, :, :), pointer :: gv, gv2, gv3
      PetscScalar, dimension(:, :, :), pointer :: hv, hv2, hv3
      PetscScalar, dimension(:, :, :), pointer :: iv, iv2, iv3

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye, &
         zs => ctx%zs, ze => ctx%ze, d=>B%deg)

      ! global to local with proper ghost cells
      call DMGlobalToLocalBegin(ctx%DA, a1, INSERT_VALUES, al1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, a1, INSERT_VALUES, al1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, a2, INSERT_VALUES, al2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, a2, INSERT_VALUES, al2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, b1, INSERT_VALUES, bl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, b1, INSERT_VALUES, bl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, b2, INSERT_VALUES, bl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, b2, INSERT_VALUES, bl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, c1, INSERT_VALUES, cl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, c1, INSERT_VALUES, cl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, c2, INSERT_VALUES, cl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, c2, INSERT_VALUES, cl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, d1, INSERT_VALUES, dl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, d1, INSERT_VALUES, dl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, d2, INSERT_VALUES, dl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, d2, INSERT_VALUES, dl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, e1, INSERT_VALUES, el1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, e1, INSERT_VALUES, el1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, e2, INSERT_VALUES, el2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, e2, INSERT_VALUES, el2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, f1, INSERT_VALUES, fl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, f1, INSERT_VALUES, fl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, f2, INSERT_VALUES, fl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, f2, INSERT_VALUES, fl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, g1, INSERT_VALUES, gl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, g1, INSERT_VALUES, gl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, g2, INSERT_VALUES, gl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, g2, INSERT_VALUES, gl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, h1, INSERT_VALUES, hl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, h1, INSERT_VALUES, hl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, h2, INSERT_VALUES, hl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, h2, INSERT_VALUES, hl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, i1, INSERT_VALUES, il1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, i1, INSERT_VALUES, il1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, i2, INSERT_VALUES, il2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, i2, INSERT_VALUES, il2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, al1, av, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, al2, av2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, a3, av3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, bl1, bv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, bl2, bv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, b3, bv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, cl1, cv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, cl2, cv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, c3, cv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, dl1, dv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, dl2, dv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, d3, dv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, el1, ev, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, el2, ev2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, e3, ev3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, fl1, fv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, fl2, fv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, f3, fv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, gl1, gv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, gl2, gv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, g3, gv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, hl1, hv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, hl2, hv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, h3, hv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, il1, iv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, il2, iv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, i3, iv3, ierr)
      CHKERRQ(ierr)

      av3 = 0.0
      bv3 = 0.0
      cv3 = 0.0
      dv3 = 0.0
      ev3 = 0.0
      fv3 = 0.0
      gv3 = 0.0
      hv3 = 0.0
      iv3 = 0.0

      do t=-d, d
         do s=-d, d
            do r=zs-1, ze-1
               do n=-d, d
                  do m=-d, d
                     do l=ys-1, ye-1
                        do k=-d, d
                           do j=-d, d
                              do i=xs-1, xe-1

                                 av3(i, l, r) = av3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    av(i+k, l+n, r+t) * av2(i+j, l+m, r+s)

                                 bv3(i, l, r) = bv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    bv(i+k, l+n, r+t) * bv2(i+j, l+m, r+s)

                                 cv3(i, l, r) = cv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    cv(i+k, l+n, r+t) * cv2(i+j, l+m, r+s)

                                 dv3(i, l, r) = dv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    dv(i+k, l+n, r+t) * dv2(i+j, l+m, r+s)

                                 ev3(i, l, r) = ev3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    ev(i+k, l+n, r+t) * ev2(i+j, l+m, r+s)

                                 fv3(i, l, r) = fv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    fv(i+k, l+n, r+t) * fv2(i+j, l+m, r+s)

                                 gv3(i, l, r) = gv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    gv(i+k, l+n, r+t) * gv2(i+j, l+m, r+s)

                                 hv3(i, l, r) = hv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    hv(i+k, l+n, r+t) * hv2(i+j, l+m, r+s)

                                 iv3(i, l, r) = iv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    iv(i+k, l+n, r+t) * iv2(i+j, l+m, r+s)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, al1, av, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, al2, av2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, a3, av3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, bl1, bv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, bl2, bv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, b3, bv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, cl1, cv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, cl2, cv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, c3, cv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, dl1, dv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, dl2, dv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, d3, dv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, el1, ev, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, el2, ev2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, e3, ev3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, fl1, fv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, fl2, fv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, f3, fv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, gl1, gv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, gl2, gv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, g3, gv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, hl1, hv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, hl2, hv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, h3, hv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, il1, iv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, il2, iv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, i3, iv3, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_nine

   subroutine nonlinear_three(ctx, grid, B, X, Y, Z, &
        u1, u2, u3, ul1, ul2, v1, v2, v3, vl1, vl2, w1, w2, w3, wl1, wl2, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal, dimension(ctx%xs-1:ctx%xe-1, -B%deg:B%deg, -B%deg:B%deg) :: X
      PetscReal, dimension(ctx%ys-1:ctx%ye-1, -B%deg:B%deg, -B%deg:B%deg) :: Y
      PetscReal, dimension(ctx%zs-1:ctx%ze-1, -B%deg:B%deg, -B%deg:B%deg) :: Z
      Vec :: u1, u2, u3, ul1, ul2
      Vec :: v1, v2, v3, vl1, vl2
      Vec :: w1, w2, w3, wl1, wl2
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m, n, r, s, t
      PetscScalar, dimension(:, :, :), pointer :: uv, uv2, uv3
      PetscScalar, dimension(:, :, :), pointer :: vv, vv2, vv3
      PetscScalar, dimension(:, :, :), pointer :: wv, wv2, wv3

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye, &
         zs => ctx%zs, ze => ctx%ze, d=>B%deg)

      ! global to local with proper ghost cells
      call DMGlobalToLocalBegin(ctx%DA, u1, INSERT_VALUES, ul1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, u1, INSERT_VALUES, ul1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, u2, INSERT_VALUES, ul2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, u2, INSERT_VALUES, ul2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, v1, INSERT_VALUES, vl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v1, INSERT_VALUES, vl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, v2, INSERT_VALUES, vl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v2, INSERT_VALUES, vl2, ierr)

      call DMGlobalToLocalBegin(ctx%DA, w1, INSERT_VALUES, wl1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, w1, INSERT_VALUES, wl1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, w2, INSERT_VALUES, wl2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, w2, INSERT_VALUES, wl2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, ul1, uv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, ul2, uv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, u3, uv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, vl1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, vl2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, v3, vv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecGetArrayF90(ctx%DA, wl1, wv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, wl2, wv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, w3, wv3, ierr)
      CHKERRQ(ierr)

      uv3 = 0.0
      vv3 = 0.0
      wv3 = 0.0

      do t=-d, d
         do s=-d, d
            do r=zs-1, ze-1
               do n=-d, d
                  do m=-d, d
                     do l=ys-1, ye-1
                        do k=-d, d
                           do j=-d, d
                              do i=xs-1, xe-1

                                 uv3(i, l, r) = uv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    uv(i+k, l+n, r+t) * uv2(i+j, l+m, r+s)

                                 vv3(i, l, r) = vv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    vv(i+k, l+n, r+t) * vv2(i+j, l+m, r+s)

                                 wv3(i, l, r) = wv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    wv(i+k, l+n, r+t) * wv2(i+j, l+m, r+s)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, ul1, uv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, ul2, uv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, u3, uv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, vl1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, vl2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, v3, vv3, ierr)
      CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(ctx%DA, wl1, wv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, wl2, wv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, w3, wv3, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_three

   subroutine nonlinear_sum_BBB(ctx, grid, B, X, Y, Z, v1, v2, v3, l1, l2, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal, dimension(ctx%xs-1:ctx%xe-1, -B%deg:B%deg, -B%deg:B%deg) :: X
      PetscReal, dimension(ctx%ys-1:ctx%ye-1, -B%deg:B%deg, -B%deg:B%deg) :: Y
      PetscReal, dimension(ctx%zs-1:ctx%ze-1, -B%deg:B%deg, -B%deg:B%deg) :: Z
      Vec :: v1, v2, v3, l1, l2
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, l, m, n, r, s, t
      PetscScalar, dimension(:, :, :), pointer :: vv, vv2, vv3

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye, &
         zs => ctx%zs, ze => ctx%ze, d=>B%deg)

      ! global to local with proper ghost cells
      call DMGlobalToLocalBegin(ctx%DA, v1, INSERT_VALUES, l1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v1, INSERT_VALUES, l1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, v2, INSERT_VALUES, l2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v2, INSERT_VALUES, l2, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, l1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, l2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, v3, vv3, ierr)
      CHKERRQ(ierr)

      vv3 = 0.0

      do t=-d, d
         do s=-d, d
            do r=zs-1, ze-1
               do n=-d, d
                  do m=-d, d
                     do l=ys-1, ye-1
                        do k=-d, d
                           do j=-d, d
                              do i=xs-1, xe-1

                                 vv3(i, l, r) = vv3(i, l, r) +             &
                                    X(i, j, k) * Y(l, m, n) * Z(r, s, t) * &
                                    vv(i+k, l+n, r+t) * vv2(i+j, l+m, r+s)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, l1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, l2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, v3, vv3, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_sum_BBB

   subroutine nonlinear_sum_BBBB(ctx, grid, B, X, Y, Z, v1, v2, v3, v4, l1, l2, l3, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal, dimension(ctx%xs-1:ctx%xe-1, -B%deg:B%deg, -B%deg:B%deg, -B%deg:B%deg) :: X
      PetscReal, dimension(ctx%ys-1:ctx%ye-1, -B%deg:B%deg, -B%deg:B%deg, -B%deg:B%deg) :: Y
      PetscReal, dimension(ctx%zs-1:ctx%ze-1, -B%deg:B%deg, -B%deg:B%deg, -B%deg:B%deg) :: Z
      Vec :: v1, v2, v3, v4, l1, l2, l3
      PetscErrorCode :: ierr

      PetscInt :: i, j, k, p, l, m, n, q, r, s, t, w
      PetscScalar, dimension(:, :, :), pointer :: vv, vv2, vv3, vv4

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye, &
         zs => ctx%zs, ze => ctx%ze, d=>B%deg)

      ! global to local with proper ghost cells
      call DMGlobalToLocalBegin(ctx%DA, v1, INSERT_VALUES, l1, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v1, INSERT_VALUES, l1, ierr)
      call DMGlobalToLocalBegin(ctx%DA, v2, INSERT_VALUES, l2, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v2, INSERT_VALUES, l2, ierr)
      call DMGlobalToLocalBegin(ctx%DA, v3, INSERT_VALUES, l3, ierr)
      call DMGlobalToLocalEnd(ctx%DA, v3, INSERT_VALUES, l3, ierr)

      ! get arrays
      call DMDAVecGetArrayF90(ctx%DA, l1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, l2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, l3, vv3, ierr)
      CHKERRQ(ierr)
      call DMDAVecGetArrayF90(ctx%DA, v4, vv4, ierr)
      CHKERRQ(ierr)

      vv4 = 0.0

      do w=-d, d
         do t=-d, d
            do s=-d, d
               do r=zs-1, ze-1
                  do q=-d, d
                     do n=-d, d
                        do m=-d, d
                           do l=ys-1, ye-1
                              do p=-d, d
                                 do k=-d, d
                                    do j=-d, d
                                       do i=xs-1, xe-1

                                          vv4(i, l, r) = vv4(i, l, r) + &
                                             X(i, j, k, p) * Y(l, m, n, q) * Z(r, s, t, w) * &
                                             vv(i+k, l+n, r+t) * vv2(i+j, l+m, r+s) * vv3(i+p, l+q, r+w)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! restore arrays
      call DMDAVecRestoreArrayF90(ctx%DA, l1, vv, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, l2, vv2, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, l3, vv3, ierr)
      CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(ctx%DA, v4, vv4, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine nonlinear_sum_BBBB

   subroutine divergence(ctx, grid, B, v1, v2, v3, v4, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: v1, v2, v3, v4
      PetscErrorCode :: ierr

      associate (xs => ctx%xs, xe => ctx%xe, ys => ctx%ys, ye => ctx%ye,   &
         zs => ctx%zs, ze => ctx%ze, d=>B%deg, BBx => B%BBx, BBy => B%BBy, &
         BBz => B%BBz, BDx => B%BDx, BDy => B%BDy, BDz => B%BDz, temp => B%temp)

      call VecSet(v4, complex(0.0, 0.0), ierr)
      CHKERRQ(ierr)

      ! BD(x) BB(y) BB(z) alpha
      call matvec(ctx, B, BDx, BBy, BBz, v1, temp, ierr)
      CHKERRQ(ierr)
      call VecAXPY(v4, complex(1.0, 0.0), temp, ierr)
      CHKERRQ(ierr)

      ! BB(x) BD(y) BB(z) alpha
      call matvec(ctx, B, BBx, BDy, BBz, v2, temp, ierr)
      CHKERRQ(ierr)
      call VecAXPY(v4, complex(1.0, 0.0), temp, ierr)
      CHKERRQ(ierr)

      ! BB(x) BB(y) BD(z) alpha
      call matvec(ctx, B, BBx, BBy, BDz, v3, temp, ierr)
      CHKERRQ(ierr)
      call VecAXPY(v4, complex(1.0, 0.0), temp, ierr)
      CHKERRQ(ierr)

      end associate
      ierr = 0
      return
   end subroutine divergence

   ! Convection routine calculates the skew-symmetric convection term.
   !
   ! Input:
   !    u/v/w: alpha for u, v and w
   !    T: alpha for advected scalar
   !    rho_u/v/w: alpha for rho * u/ rho * v/ rho * w (pre-calculated)
   !    rho_uT/vT/wT: alpha for rho * u * T/ rho * v * T/ rho * w * T (pre-calculated)
   !    temp: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine convection(ctx, grid, B, u, v, w, T, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: u, v, w, T
      Vec :: rho_u, rho_v, rho_w
      Vec :: rho_uT, rho_vT, rho_wT
      Vec :: temp, vector
      PetscErrorCode :: ierr

      ! find d(rho ui uj)/dxj
      call divergence(ctx, grid, B, rho_uT, rho_vT, rho_wT, temp, ierr)

      call VecCopy(temp, vector, ierr)

      ! find rho u dT/dx
      call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
         T, rho_u, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! find rho v dT/dy
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
         T, rho_v, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! find rho w dT/dy
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
         T, rho_w, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! -T d(rho u)/dx
      call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
         rho_u, T, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(-1.0, 0.0), temp, ierr)

      ! -T d(rho v)/dy
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
         rho_v, T, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(-1.0, 0.0), temp, ierr)

      ! -T d(rho w)/dz
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
         rho_w, T, temp, B%local1, B%local2, ierr)

      call VecAXPY(vector, complex(-1.0, 0.0), temp, ierr)

      call VecScale(vector, complex(-0.5, 0.0), ierr)

      ierr = 0
      return
   end subroutine convection

   subroutine d1dx(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BD(x) BB(y) BB(z) alpha
      call matvec(ctx, B, B%BDx, B%BBy, B%BBz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d1dx

   subroutine d1dy(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BB(x) BD(y) BB(z) alpha
      call matvec(ctx, B, B%BBx, B%BDy, B%BBz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d1dy

   subroutine d1dz(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BB(x) BB(y) BD(z) alpha
      call matvec(ctx, B, B%BBx, B%BBy, B%BDz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d1dz

   subroutine d2dx(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BL(x) BB(y) BB(z) alpha
      call matvec(ctx, B, B%BLx, B%BBy, B%BBz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d2dx

   subroutine d2dy(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BB(x) BL(y) BB(z) alpha
      call matvec(ctx, B, B%BBx, B%BLy, B%BBz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d2dy

   subroutine d2dz(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(vector, ierr)
      CHKERRQ(ierr)

      ! BB(x) BB(y) BL(z) alpha
      call matvec(ctx, B, B%BBx, B%BBy, B%BLz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine d2dz

   subroutine dxdy(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecSet(vector, complex(0.0, 0.0), ierr)
      CHKERRQ(ierr)

      ! BD(x) BD(y) BB(z) alpha
      call matvec(ctx, B, B%BDx, B%BDy, B%BBz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine dxdy

   subroutine dydz(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecSet(vector, complex(0.0, 0.0), ierr)
      CHKERRQ(ierr)

      ! BB(x) BD(y) BD(z) alpha
      call matvec(ctx, B, B%BBx, B%BDy, B%BDz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine dydz

   subroutine dxdz(ctx, grid, B, alpha, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      Vec :: alpha, vector
      PetscErrorCode :: ierr

      call VecSet(vector, complex(0.0, 0.0), ierr)
      CHKERRQ(ierr)

      ! BD(x) BB(y) BD(z) alpha
      call matvec(ctx, B, B%BDx, B%BBy, B%BDz, alpha, vector, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine dxdz

   ! Viscous routine calculates the viscous term
   !
   ! Input:
   !    i: is the ith index
   !    mu: kinematic viscosity
   !    u/v/w: alpha for u, v and w
   !    temp: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine viscous(ctx, grid, B, i, mu, Re, u, v, w, temp, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, Re
      Vec :: u, v, w
      Vec :: temp, vector
      PetscErrorCode :: ierr

      call VecSet(temp, complex(0.0, 0.0), ierr)
      call VecSet(vector, complex(0.0, 0.0), ierr)

      if (i .eq. 1) then

         ! Laplacian
         call d2dx(ctx, grid, B, u, temp, ierr)
         call VecAXPY(vector, complex(4.0/3.0, 0.0), temp, ierr)
         call d2dy(ctx, grid, B, u, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)
         call d2dz(ctx, grid, B, u, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

         ! + 1/3 (d2v/dxdy + d2w/dxdz)
         call dxdy(ctx, grid, B, v, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)
         call dxdz(ctx, grid, B, w, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)

         ! scale by mu/Re
         call VecScale(vector, complex(mu/Re, 0.0), ierr)

      else if (i .eq. 2) then

         ! Laplacian
         call d2dx(ctx, grid, B, v, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)
         call d2dy(ctx, grid, B, v, temp, ierr)
         call VecAXPY(vector, complex(4.0/3.0, 0.0), temp, ierr)
         call d2dz(ctx, grid, B, v, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

         ! + 1/3 (d2u/dxdy + d2w/dydz)
         call dxdy(ctx, grid, B, u, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)
         call dydz(ctx, grid, B, w, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)

         ! scale by mu/Re
         call VecScale(vector, complex(mu/Re, 0.0), ierr)

      else if (i .eq. 3) then

         ! Laplacian
         call d2dx(ctx, grid, B, w, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)
         call d2dy(ctx, grid, B, w, temp, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)
         call d2dz(ctx, grid, B, w, temp, ierr)
         call VecAXPY(vector, complex(4.0/3.0, 0.0), temp, ierr)

         ! + 1/3 (d2u/dxdz + d2v/dydz)
         call dxdz(ctx, grid, B, u, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)
         call dydz(ctx, grid, B, v, temp, ierr)
         call VecAXPY(vector, complex(1.0/3.0, 0.0), temp, ierr)

         ! scale by mu/Re
         call VecScale(vector, complex(mu/Re, 0.0), ierr)
      end if

      ierr = 0
      return
   end subroutine viscous

   ! Computes sum_j tau_ij dui/dxj
   !
   ! Input:
   !    i: is the ith index
   !    mu: kinematic viscosity
   !    u/v/w: alpha for u, v and w
   !    temp1: temporary vector that cannot be B%temp
   !    temp2: temporary vector that cannot be B%temp
   !    temp3: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine energy_viscous(ctx, grid, B, i, mu, gamma, M, Re, u, v, w, temp1, temp2, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, gamma, M, Re
      Vec :: u, v, w
      Vec :: temp1, temp2, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(temp1, ierr)
      call VecZeroEntries(temp2, ierr)
      call VecZeroEntries(vector, ierr)

      if (i .eq. 1) then

         ! tau_11 du/dx
         call d1dx(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(4.0/3.0, 0.0), temp1, ierr)
         call d1dy(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)
         call d1dz(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
            u, temp1, vector, B%local1, B%local2, ierr)

         ! + tau_12 du/dy
         call VecZeroEntries(temp2, ierr)

         call d1dy(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dx(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
            u, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         ! + tau_13 du/dz
         call VecZeroEntries(temp2, ierr)

         call d1dz(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dx(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
            u, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      else if (i .eq. 2) then

         ! tau_21 dv/dx
         call d1dx(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dy(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
            v, temp1, vector, B%local1, B%local2, ierr)

         ! + tau_22 dv/dy
         call VecZeroEntries(temp2, ierr)

         call d1dy(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(4.0/3.0, 0.0), temp1, ierr)
         call d1dx(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)
         call d1dz(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
            v, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         ! + tau_23 dv/dz
         call VecZeroEntries(temp2, ierr)

         call d1dz(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dy(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
            v, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      else if (i .eq. 3) then

         ! tau_31 dw/dx
         call d1dx(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dz(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
            w, temp1, vector, B%local1, B%local2, ierr)

         ! + tau_32 dw/dy
         call VecZeroEntries(temp2, ierr)

         call d1dy(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)
         call d1dz(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(1.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
            w, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

         ! + tau_33 dw/dz
         call VecZeroEntries(temp2, ierr)

         call d1dz(ctx, grid, B, w, temp1, ierr)
         call VecAXPY(temp2, complex(4.0/3.0, 0.0), temp1, ierr)
         call d1dx(ctx, grid, B, u, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)
         call d1dy(ctx, grid, B, v, temp1, ierr)
         call VecAXPY(temp2, complex(-2.0/3.0, 0.0), temp1, ierr)

         call ksp_solve(ctx, grid, B, temp2, temp1, 'full', ierr)
         CHKERRQ(ierr)

         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
            w, temp1, temp2, B%local1, B%local2, ierr)

         call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)
      end if

      call VecScale(vector, complex((gamma* (gamma - 1.0) * M**2 * mu)/Re, 0.0), ierr)

      ierr = 0
      return
   end subroutine energy_viscous

   ! Computes -1/(gamma * M^2) * d(rho * T)/dx_i
   !
   ! Input:
   !    i: is the ith index
   !    gamma: ratio of specific heats
   !    M: Mach number
   !    rho: density
   !    T: temperature
   !    temp: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine pressure(ctx, grid, B, i, gamma, M, rho, T, temp, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: gamma, M
      Vec :: rho, T
      Vec :: temp, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      if (i .eq. 1) then
         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
            rho, T, vector, B%local1, B%local2, ierr)

         call ksp_solve(ctx, grid, B, vector, temp, 'full', ierr)
         CHKERRQ(ierr)

         call d1dx(ctx, grid, B, temp, vector, ierr)

      else if (i .eq. 2) then
         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
            rho, T, vector, B%local1, B%local2, ierr)

         call ksp_solve(ctx, grid, B, vector, temp, 'full', ierr)
         CHKERRQ(ierr)

         call d1dy(ctx, grid, B, temp, vector, ierr)

      else if (i .eq. 3) then
         call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBBz, &
            rho, T, vector, B%local1, B%local2, ierr)

         call ksp_solve(ctx, grid, B, vector, temp, 'full', ierr)
         CHKERRQ(ierr)

         call d1dz(ctx, grid, B, temp, vector, ierr)
      end if

      ! scale by mu
      call VecScale(vector, complex(-1.0/(gamma * M**2), 0.0), ierr)

      ierr = 0
      return
   end subroutine pressure

   subroutine momentum(ctx, grid, B, i, mu, gamma, M, Re, &
         rho, u, v, w, T, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp1, temp2, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i
      PetscReal :: mu, gamma, M, Re
      Vec :: rho, u, v, w, T
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

      ierr = 0
      return
   end subroutine momentum

   ! Computes (gamma * mu) / (Re * Pr) * d2T/dxjdxj
   !
   ! Input:
   !    gamma: ratio of specific heats
   !    mu: kinematic viscosity
   !    Re: Reynolds number
   !    Pr: Prandtl number
   !    T: temperature
   !    temp: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine energy_laplacian(ctx, grid, B, gamma, mu, Re, Pr, T, temp, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: gamma, mu, Re, Pr
      Vec :: T
      Vec :: temp, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      call d2dx(ctx, grid, B, T, temp, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      call d2dy(ctx, grid, B, T, temp, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      call d2dz(ctx, grid, B, T, temp, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! scale by (gamma * mu) / (Re * Pr)
      call VecScale(vector, complex((gamma * mu)/(Re * Pr), 0.0), ierr)

      ierr = 0
      return
   end subroutine energy_laplacian

   ! Computes -(gamma - 1) * T * duj/dxj
   !
   ! Input:
   !    gamma: ratio of specific heats
   !    u/v/w: velocity
   !    T: temperature
   !    temp: temporary vector that cannot be B%temp
   !
   ! Output:
   !    vector: resultant vector which is not inverted
   !
   subroutine energy_divergence(ctx, grid, B, gamma, u, v, w, T, temp, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: gamma
      Vec :: u, v, w, T
      Vec :: temp, vector
      PetscErrorCode :: ierr

      call VecZeroEntries(temp, ierr)
      call VecZeroEntries(vector, ierr)

      ! T * du/dx
      call nonlinear_sum_BBB(ctx, grid, B, B%BBDx, B%FBBBy, B%BBBz, &
         u, T, temp, B%local1, B%local2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! + T * dv/dy
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%BBDy, B%BBBz, &
         v, T, temp, B%local1, B%local2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! + T * dw/dz
      call nonlinear_sum_BBB(ctx, grid, B, B%BBBx, B%FBBBy, B%BBDz, &
         w, T, temp, B%local1, B%local2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp, ierr)

      ! scale by -(gamma - 1)
      call VecScale(vector, complex(-(gamma - 1.0), 0.0), ierr)

      ierr = 0
      return
   end subroutine energy_divergence

   subroutine energy(ctx, grid, B, mu, gamma, M, Re, Pr, &
         rho, u, v, w, T, rho_u, rho_v, rho_w, &
         rho_uT, rho_vT, rho_wT, temp1, temp2, temp3, vector, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscReal :: mu, gamma, M, Re, Pr
      Vec :: rho, u, v, w, T
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

      call energy_divergence(ctx, grid, B, gamma, &
         u, v, w, T, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      do i=1, 3
         call energy_viscous(ctx, grid, B, i, mu, gamma, M, Re, &
            u, v, w, temp1, temp2, temp3, ierr)
         call VecAXPY(vector, complex(1.0, 0.0), temp3, ierr)
      end do

      call energy_laplacian(ctx, grid, B, gamma, &
         mu, Re, Pr, T, temp1, temp2, ierr)
      call VecAXPY(vector, complex(1.0, 0.0), temp2, ierr)

      ierr = 0
      return
   end subroutine energy

end module linalg

! vim: set ts=3 sw=3 sts=3:
