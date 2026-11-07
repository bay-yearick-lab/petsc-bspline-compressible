# PETSc solver options

The solver builds every distributed object through PETSc and calls the
`*SetFromOptions` routines, so the Krylov method, preconditioner, direct
solver and domain decomposition can all be changed on the command line or
through an options file (`-options_file examples/options/<file>.opts`).

## Objects created in the code

| Object | Where | Communicator | Compiled-in settings |
|---|---|---|---|
| `DM` 3D (`DMDACreate3d`) | `src/lattice.f90` | `PETSC_COMM_WORLD` | 77 x 120 x 77 coefficients, `DM_BOUNDARY_PERIODIC` / `GHOSTED` / `PERIODIC`, `DMDA_STENCIL_BOX`, stencil width 2p = 6, ranks `PETSC_DECIDE x 1 x PETSC_DECIDE` |
| `DM` 1D (`DMDACreate1d`) x 3 | `src/lattice.f90` | `commx`, `commy`, `commz` from `MPI_Comm_split` | ownership ranges copied from the 3D DMDA, stencil width p = 3 |
| `Mat` (`DMCreateMatrix`, `MATAIJ`) | `src/operator.f90` | line communicator | 1D Galerkin mass `M_y`, `BB`, first-derivative `BD` and stiffness `BL` matrices per direction |
| `Mat` (`MatCreateFFT`, `MATFFTW`) | `src/linalg.f90` | `commx`, `commz` | parallel FFTW transforms that diagonalize the circulant periodic mass matrices |
| `KSP` `ksp` | `src/linalg.f90` (`ksp_init`) | `commy` | `KSPGMRES` on `M_y` with Dirichlet rows, `rtol = atol = 1e-14` |
| `KSP` `fullksp` | `src/linalg.f90` (`ksp_init`) | `commy` | `KSPCG` on the symmetric positive definite mass matrix `BB_y`, `rtol = atol = 1e-14` |

Both `KSP` objects call `KSPSetFromOptions`, so the options below apply to
the wall-normal line solves. Each rank owns complete wall-normal lines
(one rank in y), so the y systems are solved rank-locally with no
communication, while the x and z FFTs run across `commx` and `commz`.

## Krylov method (`-ksp_*`)

| Option | Effect |
|---|---|
| `-ksp_type gmres` | compiled-in default for the boundary-condition system `M_y` |
| `-ksp_type cg` | compiled-in default for the SPD mass matrix `BB_y` |
| `-ksp_type preonly` | apply the preconditioner once (use with a direct factorization) |
| `-ksp_gmres_restart <m>` | GMRES restart length (PETSc default 30) |
| `-ksp_rtol`, `-ksp_atol` | relative and absolute residual tolerances (code sets 1e-14) |
| `-ksp_max_it <n>` | iteration cap |
| `-ksp_monitor`, `-ksp_converged_reason` | per-iteration residuals and termination reason |
| `-ksp_view` | print the full solver configuration |

## Preconditioner (`-pc_*`)

| Option | Effect |
|---|---|
| `-pc_type ilu` | incomplete LU, PETSc default for the rank-local AIJ line matrices |
| `-pc_type icc` | incomplete Cholesky for the SPD system solved with CG |
| `-pc_type jacobi` | diagonal scaling |
| `-pc_type lu`, `-pc_type cholesky` | exact factorization, reused across all lines and time steps |
| `-pc_factor_mat_solver_type superlu_dist` | sparse direct LU from SuperLU_DIST |
| `-mat_superlu_dist_colperm METIS_AT_PLUS_A` | METIS fill-reducing column ordering for SuperLU_DIST |

SuperLU_DIST direct solve:

```
-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist
```

## Domain decomposition

`DMDACreate3d` is called with `PETSC_DECIDE` ranks in x and z and one rank in
y, so `mpirun -n P` produces a P = px x pz pencil decomposition chosen by
PETSc from the coefficient counts. The 1D line DMDAs reuse the 3D ownership
ranges, so a rank's slab of every x-, y- and z-line matches its 3D block
exactly and no redistribution is needed between the 3D field and the line
solves. `-dm_view` prints the layout.

## Profiling

| Option | Effect |
|---|---|
| `-log_view` | time, flop rate and message counts per event (`KSPSolve`, `MatMult`, `VecScatterBegin/End`, `DMGlobalToLocal`) |
| `-memory_view` | memory usage summary |

## Driver options

| Option | Default | Effect |
|---|---|---|
| `-nsteps <n>` | 30000 (t = 30 at dt = 1e-3) | number of RK3 time steps in `channel_flow` |
| `-output_every <n>` | 10 | write fields every n steps |
| `-test <name>` | `nonlinear_combined_rhs` | select a manufactured-solution check in `test_operators` |
