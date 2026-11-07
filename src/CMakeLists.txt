add_library(petb
  environment.f90
  types.f90
  lattice.f90
  bspline.f90
  gausstable.f90
  quad.f90
  operator.f90
  mms.f90
  rhs.f90
  error.f90
  linalg.f90
  rk.f90
  output.f90)
target_compile_options(petb PUBLIC ${PETB_Fortran_FLAGS})
target_include_directories(petb PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(petb PUBLIC PkgConfig::PETSC MPI::MPI_Fortran)
