! Gather of distributed fields to rank 0 in natural ordering and ASCII output.
# include <petsc/finclude/petscsys.h>
# include <petsc/finclude/petscvec.h>

module output
   use types
   implicit none

   private

   public :: io_init
   public :: io_destroy
   public :: io_fortran_scatter
   public :: io_fortran_write
   public :: index_3d_to_1d
contains

   subroutine io_destroy(B, ierr)
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      call VecScatterDestroy(B%scatter, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%natural, ierr)
      CHKERRQ(ierr)
      call VecDestroy(B%vecseq, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine io_destroy

   subroutine io_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      ! do this only once
      call DMDACreateNaturalVector(ctx%DA, B%natural, ierr)
      call VecScatterCreateToZero(B%natural, B%scatter, B%vecseq, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine io_init

   subroutine io_fortran_scatter(ctx, B, vector, ierr)
      type(ctx_t) :: ctx
      type(bspline_t) :: B
      Vec :: vector
      PetscErrorCode, intent(out) :: ierr

      ! scatter as much as you need using scatter context
      call DMDAGlobalToNaturalBegin(ctx%DA, vector, INSERT_VALUES, B%natural, ierr)
      CHKERRQ(ierr)
      call DMDAGlobalToNaturalEnd(ctx%DA, vector, INSERT_VALUES, B%natural, ierr)
      CHKERRQ(ierr)
      call VecScatterBegin(B%scatter, B%natural, B%vecseq, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CHKERRQ(ierr)
      call VecScatterEnd(B%scatter, B%natural, B%vecseq, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CHKERRQ(ierr)

      ierr = 0
      return
   end subroutine io_fortran_scatter

   subroutine io_fortran_write(ctx, grid, B, itime, var, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: itime
      character(len=*) :: var
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i, j, k, v, N
      PetscReal :: dx, dy, dz
      character(10) :: otime
      PetscScalar, dimension(:), pointer :: vv

      if (ctx%rank .eq. 0) then
         write(otime,"(I10.10)")itime

         dx = grid%Lx / (grid%xNc-1)
         dy = grid%Ly / (grid%Ny-1)
         dz = grid%Lz / (grid%zNc-1)

         call VecGetArrayReadF90(B%vecseq, vv, ierr)
         CHKERRQ(ierr)

         open(1,file= var//"_"//otime//".out",form='FORMATTED')
         do k=0, grid%zNc-1
            do j=0, grid%Ny-1
               do i=0, grid%xNc-1

                  v = index_3d_to_1d(grid, i, j, k)
                  write(1,"(4F20.14)") i*dx, j*dy, k*dz, Real(vv(v))
               end do
            end do
         end do
         close(1)

         call VecRestoreArrayReadF90(B%vecseq, vv, ierr)
         CHKERRQ(ierr)
      end if

      ierr = 0
      return
   end subroutine io_fortran_write

   function index_3d_to_1d(grid, i, j, k) result (v)
      type(grid_t) :: grid
      PetscInt :: i, j, k

      PetscInt :: v

      v = i + grid%xNc * (j + grid%Ny * k) + 1
   end function

end module

! vim: set ts=3 sw=3 sts=3:
