cmake_minimum_required(VERSION 3.13)
project(petsc_bspline_compressible LANGUAGES C Fortran)

# Requirements: an MPI Fortran compiler (gfortran + MPICH or Open MPI) and
# PETSc 3.10.x configured with complex scalars, FFTW, SuperLU_DIST and ParMETIS
# (see scripts/install_petsc.sh). PETSc is located from PETSC_DIR/PETSC_ARCH.

if(NOT DEFINED PETSC_DIR)
  set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "PETSc source or install directory")
endif()
if(NOT DEFINED PETSC_ARCH)
  set(PETSC_ARCH "$ENV{PETSC_ARCH}" CACHE STRING "PETSc build architecture")
endif()
if(NOT PETSC_DIR)
  message(FATAL_ERROR "Set PETSC_DIR (and PETSC_ARCH for in-tree builds).")
endif()

# PETSc ships a pkg-config file in <PETSC_DIR>/<PETSC_ARCH>/lib/pkgconfig
set(ENV{PKG_CONFIG_PATH}
    "${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:${PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc>=3.10)
message(STATUS "PETSc ${PETSC_VERSION} from ${PETSC_DIR}/${PETSC_ARCH}")

find_package(MPI REQUIRED COMPONENTS Fortran)

# gfortran flags used by the solver: preprocessed sources, free-form long lines,
# implicit none everywhere, double-precision default reals
set(PETB_Fortran_FLAGS -cpp -ffree-line-length-none -fimplicit-none -fdefault-real-8)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  list(APPEND PETB_Fortran_FLAGS -fallow-argument-mismatch)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(BUILD_SHARED_LIBS ON)

enable_testing()
add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(tests)
