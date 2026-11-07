! Low-storage three-stage third-order Runge-Kutta coefficients.
# include <petsc/finclude/petscsys.h>

module rk
   use types
   implicit none

   private

   public :: rk_init
contains

   subroutine rk_init(rk, ierr)
      type(rk_t) :: rk
      PetscErrorCode, intent(out) :: ierr

      rk%s = 3

      allocate(rk%a(rk%s, rk%s))
      allocate(rk%b(rk%s), rk%c(rk%s))

      rk%a = 0.0; rk%b = 0.0; rk%c = 0.0

      rk%a(2, 1) =  8.0/15.0
      rk%a(3, 1) =  1.0/4.0
      rk%a(3, 2) =  5.0/12.0

      rk%b(1) = 1.0/4.0
      rk%b(3) = 3.0/4.0

      rk%c(2) = 8.0/15.0
      rk%c(3) = 2.0/3.0

      ierr = 0
      return
   end subroutine rk_init

end module

! vim: set ts=3 sw=3 sts=3:
