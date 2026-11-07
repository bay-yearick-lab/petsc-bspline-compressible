#!/usr/bin/env bash
# Build PETSc 3.10.2 with the configuration this solver needs:
#   complex scalars, MPI, FFTW (parallel FFTs in the periodic directions),
#   SuperLU_DIST (sparse direct LU) and METIS/ParMETIS (fill-reducing ordering).
#
# Usage: MPICC=mpicc MPIF90=mpif90 scripts/install_petsc.sh [PETSC_PREFIX] [debug|opt]
# Optional: SUPERLU_DIST_ROOT, PARMETIS_ROOT point at existing installs
#           (otherwise PETSc downloads and builds them).
set -euf -o pipefail

pkgver=3.10.2
url=https://gitlab.com/petsc/petsc/-/archive/v${pkgver}/petsc-v${pkgver}.tar.gz

prefix_root="${1:-${HOME}/opt/petsc}"
buildtype="${2:-opt}"
petsc_arch="linux-gnu-complex-${buildtype}"
prefix="${prefix_root}/${petsc_arch}"
builddir="build-petsc-${pkgver}-${buildtype}"
MPICC="${MPICC:-mpicc}"
MPIF90="${MPIF90:-mpif90}"
PYTHON="${PYTHON:-python3}"

mkdir -p "${builddir}"
cd "${builddir}"

echo "Downloading PETSc ${pkgver}..."
wget -q -N "${url}"
tar xf "petsc-v${pkgver}.tar.gz"
cd "petsc-v${pkgver}"

export PETSC_DIR="$(pwd)"
export PETSC_ARCH="${petsc_arch}"

if [ "${buildtype}" = "debug" ]; then
   conf_flags=(--with-debugging=1)
else
   conf_flags=(--with-debugging=0
               COPTFLAGS="-O3 -march=native"
               FOPTFLAGS="-O3 -march=native")
fi

if [ -n "${PARMETIS_ROOT:-}" ]; then
   conf_flags+=(--with-metis-dir="${PARMETIS_ROOT}" --with-parmetis-dir="${PARMETIS_ROOT}")
else
   conf_flags+=(--download-metis --download-parmetis)
fi
if [ -n "${SUPERLU_DIST_ROOT:-}" ]; then
   conf_flags+=(--with-superlu_dist-dir="${SUPERLU_DIST_ROOT}")
else
   conf_flags+=(--download-superlu_dist)
fi

"${PYTHON}" ./configure \
   --with-cc="${MPICC}" --with-fc="${MPIF90}" --with-cxx=0 \
   --PETSC_ARCH="${petsc_arch}" \
   --prefix="${prefix}" \
   --with-shared-libraries=1 \
   --with-scalar-type=complex \
   --with-fortran-kernels=1 \
   --with-mpi=1 \
   --download-fftw \
   "${conf_flags[@]}"

make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" all
make PETSC_DIR="${PETSC_DIR}" PETSC_ARCH="${PETSC_ARCH}" install

echo "export PETSC_DIR=${prefix}"
echo "export PETSC_ARCH="
