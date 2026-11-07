! 3D DMDA over the channel, 1D DMDA line communicators (MPI_Comm_split) and grid coordinates.
# include <petsc/finclude/petscsys.h>

module lattice
   use types
   implicit none

   private

   public :: lattice_init
   public :: grid_init
   public :: lattice_free

contains

   subroutine lattice_init(ctx, grid, B, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      type(bspline_t) :: B
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: color
      MPI_Comm :: comm, size

      B%deg = 3

      grid%Nx = 80
      grid%Ny = 120
      grid%Nz = 80

      grid%xNc = grid%Nx - B%deg
      grid%yNc = grid%Ny - B%deg
      grid%zNc = grid%Nz - B%deg

      grid%Lx = 4.0 * PETSC_PI
      grid%Ly = 2.0
      grid%Lz = 4.0 * PETSC_PI / 3.0

      grid%dx = grid%Lx / grid%xNc
      grid%dy = grid%Ly / grid%yNc
      grid%dz = grid%Lz / grid%zNc

      call DMDACreate3D(PETSC_COMM_WORLD,                                                &
                        DM_BOUNDARY_PERIODIC, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX,                                                &
                        grid%xNc, grid%Ny, grid%zNc,                                     &
                        PETSC_DECIDE, 1, PETSC_DECIDE,                                   &
                        1, B%deg*2,                                                      &
                        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,      &
                        ctx%DA, ierr)
      call DMSetFromOptions(ctx%DA, ierr)
      call DMSetUp(ctx%DA, ierr)

      call DMDAGetInfo(ctx%DA, PETSC_NULL_INTEGER, ctx%mx, ctx%my, ctx%mz,         &
                       ctx%procx, ctx%procy, ctx%procz,                            &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      ! get global rank of current process
      call MPI_Comm_rank(PETSC_COMM_WORLD, ctx%rank, ierr)
      CHKERRQ(ierr)

      ! compute (i, j, k) rank coordinates in the 3D grid mapping of the DM
      ctx%rankx = mod(mod(ctx%rank, ctx%procx * ctx%procy), ctx%procx)
      ctx%ranky = mod(ctx%rank, ctx%procx * ctx%procy) / (ctx%procx)
      ctx%rankz = ctx%rank / (ctx%procx * ctx%procy)

      ! allocate and get process ownership ranges
      if (.not. allocated(ctx%rangex)) allocate(ctx%rangex(0:ctx%procx - 1))
      if (.not. allocated(ctx%rangey)) allocate(ctx%rangey(0:ctx%procy - 1))
      if (.not. allocated(ctx%rangez)) allocate(ctx%rangez(0:ctx%procz - 1))

      call DMDAGetOwnershipRanges(ctx%DA, &
         ctx%rangex, ctx%rangey, ctx%rangez, ierr)
      CHKERRQ(ierr)

      ! comm_split in x-dimension
      color = ctx%procy * ctx%rankz + ctx%ranky
      call MPI_Comm_split(PETSC_COMM_WORLD, color, ctx%rankx, comm, ierr)
      CHKERRQ(ierr)

      call DMDACreate1d(comm, DM_BOUNDARY_PERIODIC, grid%xNc, 1, B%deg, ctx%rangex, ctx%DAx, ierr)
      CHKERRQ(ierr)

      call DMSetFromOptions(ctx%DAx, ierr)
      CHKERRQ(ierr)
      call DMSetUp(ctx%DAx, ierr)
      CHKERRQ(ierr)
      call MPI_Comm_free(comm, ierr)
      CHKERRQ(ierr)

      ! comm_split in y-dimension
      color = ctx%procz * ctx%rankx + ctx%rankz
      call MPI_Comm_split(PETSC_COMM_WORLD, color, ctx%ranky, comm, ierr)
      CHKERRQ(ierr)

      call DMDACreate1d(comm, DM_BOUNDARY_GHOSTED, grid%Ny, 1, B%deg, ctx%rangey, ctx%DAy, ierr)
      CHKERRQ(ierr)

      call DMSetFromOptions(ctx%DAy, ierr)
      CHKERRQ(ierr)
      call DMSetUp(ctx%DAy, ierr)
      CHKERRQ(ierr)
      call MPI_Comm_free(comm, ierr)
      CHKERRQ(ierr)

      ! comm_split in z-dimension
      color = ctx%procx * ctx%ranky + ctx%rankx
      call MPI_Comm_split(PETSC_COMM_WORLD, color, ctx%rankz, comm, ierr)
      CHKERRQ(ierr)

      call DMDACreate1d(comm, DM_BOUNDARY_PERIODIC, grid%zNc, 1, B%deg, ctx%rangez, ctx%DAz, ierr)
      CHKERRQ(ierr)

      call DMSetFromOptions(ctx%DAz, ierr)
      CHKERRQ(ierr)
      call DMSetUp(ctx%DAz, ierr)
      CHKERRQ(ierr)
      call MPI_Comm_free(comm, ierr)
      CHKERRQ(ierr)

      ! get comm from dmda object
      call PetscObjectGetComm(ctx%DAx, ctx%commx, ierr)
      CHKERRQ(ierr)
      call PetscObjectGetComm(ctx%DAy, ctx%commy, ierr)
      CHKERRQ(ierr)
      call PetscObjectGetComm(ctx%DAz, ctx%commz, ierr)
      CHKERRQ(ierr)
      return
   end subroutine lattice_init

   subroutine grid_init(ctx, grid, ierr)
      implicit none

      type(ctx_t) :: ctx
      type(grid_t) :: grid
      PetscErrorCode, intent(out) :: ierr

      PetscInt :: i, j, k
      PetscInt :: xs, xm, ys, ym, zs, zm

      ! assign grid to each process
      call DMDAGetCorners(ctx%DA, ctx%xs, ctx%ys, ctx%zs, ctx%xm, ctx%ym, ctx%zm, ierr)
      CHKERRQ(ierr)

      ! check grid allocation in each 1d dm
      call DMDAGetCorners(ctx%DAx, xs, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   xm, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      call DMDAGetCorners(ctx%DAy, ys, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   ym, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      call DMDAGetCorners(ctx%DAz, zs, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,     &
                                   zm, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)

      ! should be zeros
      if (ctx%xs - xs .ne. 0.0 .or. ctx%xm - xm .ne. 0.0 .or. &
          ctx%ys - ys .ne. 0.0 .or. ctx%ym - ym .ne. 0.0 .or. &
          ctx%zs - zs .ne. 0.0 .or. ctx%zm - zm .ne. 0.0) then
         print*, "grid_init: Error in grid allocation"
      end if

      call DMDAGetGhostCorners(ctx%DA, ctx%gxs, ctx%gys, ctx%gzs,  &
                                       ctx%gxm, ctx%gym, ctx%gzm, ierr)
      CHKERRQ(ierr)

      ctx%xs = ctx%xs + 1
      ctx%ys = ctx%ys + 1
      ctx%zs = ctx%zs + 1
      ctx%gxs = ctx%gxs + 1
      ctx%gys = ctx%gys + 1
      ctx%gzs = ctx%gzs + 1

      ctx%xe = ctx%xs + ctx%xm - 1
      ctx%ye = ctx%ys + ctx%ym - 1
      ctx%ze = ctx%zs + ctx%zm - 1
      ctx%gxe = ctx%gxs + ctx%gxm - 1
      ctx%gye = ctx%gys + ctx%gym - 1
      ctx%gze = ctx%gzs + ctx%gzm - 1

      allocate(grid%x(ctx%xs:ctx%xe))
      allocate(grid%y(ctx%ys:ctx%ye))
      allocate(grid%z(ctx%zs:ctx%ze))

      do i=ctx%xs, ctx%xe
         if (i .eq. 1) then
            grid%x(i) = 0.0
         else if (i .gt. grid%xNc) then
            grid%x(i) = grid%x(i-1) + grid%dx
         else
            grid%x(i) = (i-1) * grid%dx
         end if
      end do

      do j=ctx%ys, ctx%ye
         if (j .eq. 1) then
            grid%y(j) = 0.0
         else if (j .gt. grid%yNc) then
            grid%y(j) = grid%Ly
         else
            grid%y(j) = (j-1) * grid%dy
         end if
      end do

      do k=ctx%zs, ctx%ze
         if (k.eq.1) then
            grid%z(k) = 0.0
         else if (k .gt. grid%zNc) then
            grid%z(k) = grid%z(k-1) + grid%dz
         else
            grid%z(k) = (k-1) * grid%dz
         end if
      end do

      ierr = 0
      return
   end subroutine grid_init

   subroutine lattice_free(ctx, ierr)
      implicit none

      type(ctx_t) :: ctx
      PetscErrorCode, intent(out) :: ierr

      call DMDestroy(ctx%DA, ierr)
      CHKERRQ(ierr)
      call DMDestroy(ctx%DAx, ierr)
      CHKERRQ(ierr)
      call DMDestroy(ctx%DAy, ierr)
      CHKERRQ(ierr)
      call DMDestroy(ctx%DAz, ierr)
      CHKERRQ(ierr)
      return
   end subroutine lattice_free

end module lattice

! vim: set ts=3 sw=3 sts=3:
