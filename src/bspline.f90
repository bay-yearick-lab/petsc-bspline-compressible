! Uniform knot vectors and recursive B-spline basis evaluation with derivatives.
# include <petsc/finclude/petscsys.h>

module bspline
   use types
   implicit none

   private

   public :: bspline_init
   public :: initxknot
   public :: inityknot
   public :: initzknot
   public :: xdBIJK
   public :: ydBIJK
   public :: zdBIJK
   public :: xBJK
   public :: yBJK
   public :: zBJK
contains

   subroutine check_index(s, e, m, d, N, minsize, maxsize)
      PetscInt :: s, e, m, d, N, minsize, maxsize

      if (s.eq.1) then
         minsize = s - (d + 2)
      else
         minsize = s - (2 * d + 2)
      end if

      if (e.eq.m) then
         maxsize = N + (2*d + 1)
      else
         maxsize = e + (2*d + 2)
      end if
   end subroutine

   subroutine bspline_init(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      call initxknot(ctx, grid, B, ierr)
      call inityknot(ctx, grid, B, ierr)
      call initzknot(ctx, grid, B, ierr)

      ierr = 0
      return
   end subroutine bspline_init

   subroutine initxknot(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: j, minsize, maxsize

      call check_index(ctx%xs, ctx%xe, ctx%mx, B%deg, grid%Nx, minsize, maxsize)

      allocate(B%Xk(minsize:maxsize))
      B%Xk(ctx%xs-1:ctx%xe-1) = grid%x

      do j = ctx%xs-2, minsize, -1
         B%Xk(j) = B%Xk(j+1) - grid%dx
      end do

      do j = ctx%xe, maxsize
         B%Xk(j) = B%Xk(j-1) + grid%dx
      end do
      return
   end subroutine initxknot

   subroutine inityknot(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: j, minsize, maxsize

      call check_index(ctx%ys, ctx%ye, ctx%my, B%deg, grid%Ny, minsize, maxsize)

      if (ctx%procy .gt. 1) then
         if (ctx%ye .eq. ctx%my) then
            maxsize = maxsize - B%deg

            allocate(B%Yk(minsize:maxsize))

            B%Yk(ctx%ys-1:ctx%ye-1) = grid%y

            B%Yk(ctx%ye:maxsize) = grid%Ly

            do j = ctx%ys-2, minsize, -1
               B%Yk(j) = B%Yk(j+1) - grid%dy
            end do

         else if (ctx%ys .eq. 1) then
            allocate(B%Yk(minsize:maxsize))

            B%Yk(ctx%ys-1:ctx%ye-1) = grid%y

            B%Yk(minsize:ctx%ys-2) = 0.0

            do j = ctx%ye, maxsize
               B%Yk(j) = B%Yk(j-1) + grid%dy
            end do
         end if
      else
         allocate(B%Yk(minsize:maxsize))

         B%Yk(ctx%ys-1:ctx%ye-1) = grid%y

         B%Yk(minsize:ctx%ys-2) = 0.0

         B%Yk(ctx%ye:maxsize) = grid%Ly
      end if
      return
   end subroutine inityknot

   subroutine initzknot(ctx, grid, B, ierr)
      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: j, minsize, maxsize

      call check_index(ctx%zs, ctx%ze, ctx%mz, B%deg, grid%Nz, minsize, maxsize)

      allocate(B%Zk(minsize:maxsize))
      B%Zk(ctx%zs-1:ctx%ze-1) = grid%z

      do j = ctx%zs-2, minsize, -1
         B%Zk(j) = B%Zk(j+1) - grid%dz
      end do

      do j = ctx%ze, maxsize
         B%Zk(j) = B%Zk(j-1) + grid%dz
      end do
      return
   end subroutine initzknot

   function xdBIJK(grid, B, i, j, k, xx) result (dB)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i, j, k, l
      PetscReal :: dB, xx

      dB = 0.0
      do l = j-i, j
         dB = dB + xalpha(B, j, l, i, k) * xBJK(grid, B, l, k-i, xx)
      end do
   end function xdBIJK

   recursive function xalpha(B, j, l, i, k) result (a)
      type(bspline_t) :: B
      PetscInt :: j, l, i, k
      PetscReal :: a

      if (j .eq. l .and. i .eq. 0) then
         a = 1.0
      else if (i .eq. 0) then
         a = 0.0
      else if ((B%Xk(l) - B%Xk(l-k+i-1)) .ne. 0.0) then
         a = REAL(k - i + 1) * (xalpha(B, j, l+1, i-1, k) -  &
                                xalpha(B, j, l  , i-1, k)) / &
             (B%Xk(l) - B%Xk(l-k+i-1))
      else
         a = 0.0
      end if
   end function xalpha

   recursive function xBJK(grid, B, j, k, xx) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: xx, Bl

      if (j .eq. 1 .and. xx .eq. 0.0 .or. &
          j .eq. grid%Nx .and. xx .eq. grid%Lx) then
         Bl = 1.0
      else
         Bl = xBJKr(grid, B, j, k, xx)
      end if
   end function xBJK

   recursive function xBJKr(grid, B, j, k, xx) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: xx, Bl

      PetscReal :: b1, b2

      if (k .eq. 0) then
         if (xx .gt. B%Xk(j-1) .and. xx .lt. B%Xk(j)) then
            Bl = 1.0
         else if (xx .eq. B%Xk(j-1) .or. xx .eq. B%Xk(j)) then
            Bl = 0.5
         else
            Bl = 0.0
         end if
      else
         b1 = xBJK(grid, B, j-1, k-1, xx)
         b2 = xBJK(grid, B, j, k-1, xx)
         Bl = 0.0
         if (b1 .ne. 0.0 .and. B%Xk(j-1) .ne. B%Xk(j-k-1)) then
            Bl = (xx - B%Xk(j-k-1)) / (B%Xk(j-1) - B%Xk(j-k-1)) * b1
         end if
         if (b2 .ne. 0.0 .and. B%Xk(j) .ne. B%Xk(j-k)) then
            Bl = Bl + (B%Xk(j) - xx) / (B%Xk(j) - B%Xk(j-k)) * b2
         end if
      end if
   end function xBJKr

   function ydBIJK(grid, B, i, j, k, yy) result (dB)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i, j, k, l
      PetscReal :: dB, yy

      dB = 0.0
      do l = j-i, j
         dB = dB + yalpha(B, j, l, i, k) * yBJK(grid, B, l, k-i, yy)
      end do
   end function ydBIJK

   recursive function yalpha(B, j, l, i, k) result (a)
      type(bspline_t) :: B
      PetscInt :: j, l, i, k
      PetscReal :: a

      if (j .eq. l .and. i .eq. 0) then
         a = 1.0
      else if (i .eq. 0) then
         a = 0.0
      else if ((B%Yk(l) - B%Yk(l-k+i-1)) .ne. 0.0) then
         a = REAL(k - i + 1) * (yalpha(B, j, l+1, i-1, k) -  &
                                yalpha(B, j, l  , i-1, k)) / &
             (B%Yk(l) - B%Yk(l-k+i-1))
      else
         a = 0.0
      end if
   end function yalpha

   recursive function yBJK(grid, B, j, k, yy) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: yy, Bl

      if (j .eq. 1 .and. yy .eq. 0.0 .or. &
          j .eq. grid%Ny .and. yy .eq. grid%Ly) then
         Bl = 1.0
      else
         Bl = yBJKr(grid, B, j, k, yy)
      end if
   end function yBJK

   recursive function yBJKr(grid, B, j, k, yy) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: yy, Bl

      PetscReal :: b1, b2

      if (k .eq. 0) then
         if (yy .gt. B%Yk(j-1) .and. yy .lt. B%Yk(j)) then
            Bl = 1.0
         else if (yy .eq. B%Yk(j-1) .or. yy .eq. B%Yk(j)) then
            Bl = 0.5
         else
            Bl = 0.0
         end if
      else
         b1 = yBJK(grid, B, j-1, k-1, yy)
         b2 = yBJK(grid, B, j, k-1, yy)
         Bl = 0.0
         if (b1 .ne. 0.0 .and. B%Yk(j-1) .ne. B%Yk(j-k-1)) then
            Bl = (yy - B%Yk(j-k-1)) / (B%Yk(j-1) - B%Yk(j-k-1)) * b1
         end if
         if (b2 .ne. 0.0 .and. B%Yk(j) .ne. B%Yk(j-k)) then
            Bl = Bl + (B%Yk(j) - yy) / (B%Yk(j) - B%Yk(j-k)) * b2
         end if
      end if
   end function yBJKr

   function zdBIJK(grid, B, i, j, k, zz) result (dB)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: i, j, k, l
      PetscReal :: dB, zz

      dB = 0.0
      do l = j-i, j
         dB = dB + zalpha(B, j, l, i, k) * zBJK(grid, B, l, k-i, zz)
      end do
   end function zdBIJK

   recursive function zalpha(B, j, l, i, k) result (a)
      type(bspline_t) :: B
      PetscInt :: j, l, i, k
      PetscReal :: a

      if (j .eq. l .and. i .eq. 0) then
         a = 1.0
      else if (i .eq. 0) then
         a = 0.0
      else if ((B%Zk(l) - B%Zk(l-k+i-1)) .ne. 0.0) then
         a = REAL(k - i + 1) * (zalpha(B, j, l+1, i-1, k) -  &
                                zalpha(B, j, l  , i-1, k)) / &
             (B%Zk(l) - B%Zk(l-k+i-1))
      else
         a = 0.0
      end if
   end function zalpha

   recursive function zBJK(grid, B, j, k, zz) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: zz, Bl

      if (j .eq. 1 .and. zz .eq. 0.0 .or. &
          j .eq. grid%Nz .and. zz .eq. grid%Lz) then
         Bl = 1.0
      else
         Bl = zBJKr(grid, B, j, k, zz)
      end if
   end function zBJK

   recursive function zBJKr(grid, B, j, k, zz) result (Bl)
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscInt :: j, k
      PetscReal :: zz, Bl

      PetscReal :: b1, b2

      if (k .eq. 0) then
         if (zz .gt. B%Zk(j-1) .and. zz .lt. B%Zk(j)) then
            Bl = 1.0
         else if (zz .eq. B%Zk(j-1) .or. zz .eq. B%Zk(j)) then
            Bl = 0.5
         else
            Bl = 0.0
         end if
      else
         b1 = zBJK(grid, B, j-1, k-1, zz)
         b2 = zBJK(grid, B, j, k-1, zz)
         Bl = 0.0
         if (b1 .ne. 0.0 .and. B%Zk(j-1) .ne. B%Zk(j-k-1)) then
            Bl = (zz - B%Zk(j-k-1)) / (B%Zk(j-1) - B%Zk(j-k-1)) * b1
         end if
         if (b2 .ne. 0.0 .and. B%Zk(j) .ne. B%Zk(j-k)) then
            Bl = Bl + (B%Zk(j) - zz) / (B%Zk(j) - B%Zk(j-k)) * b2
         end if
      end if
   end function zBJKr

end module bspline

! vim: set ts=3 sw=3 sts=3:
