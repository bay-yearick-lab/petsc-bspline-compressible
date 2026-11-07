! MPI and PETSc initialization and finalization.
# include <petsc/finclude/petscsys.h>

module environment
   use types
   implicit none

   private

   public :: environment_init
   public :: environment_free

contains

   subroutine environment_init(self, ierr)
      implicit none

      type(ctx_t) :: self
      PetscErrorCode, intent(out) :: ierr

      ! initialize MPI
      call MPI_Init(ierr)

      ! initialize Petsc
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

      if (ierr .ne. 0) then
         print*, 'Unable to initialize PETSc'
      end if

      ierr = 0
      return
   end subroutine environment_init

   subroutine environment_free(ierr)
      implicit none

      PetscErrorCode, intent(out) :: ierr

      ! finalize PETSc
      call PetscFinalize(ierr)

      ! finalize MPI
      call MPI_Finalize(ierr)

      ierr = 0
      return
   end subroutine environment_free

end module environment

! vim: set ts=3 sw=3 sts=3:
