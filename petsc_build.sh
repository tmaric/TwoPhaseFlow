#!/usr/bin/env bash
set -euo pipefail

# PETSc + petsc4Foam installer for TwoPhaseFlow.
#
# Default behavior:
# 1) Clone/update PETSc in $HOME/petsc
# 2) Configure PETSc with MPI + Helmholtz-relevant packages
# 3) Build PETSc
# 4) Build external/petsc4Foam (after sourcing OpenFOAM environment)
#
# Usage examples:
#   ./petsc_build.sh
#   ./petsc_build.sh --install-deps
#   ./petsc_build.sh --petsc-dir "$HOME/petsc" --petsc-arch arch-linux-c-opt -j 16
#   ./petsc_build.sh --skip-petsc4foam

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PETSC_DIR_DEFAULT="$HOME/petsc"
PETSC_ARCH_DEFAULT="arch-linux-c-opt"
PETSC_BRANCH_DEFAULT="release"
OPENFOAM_BASHRC_DEFAULT="$HOME/openfoam/etc/bashrc"

PETSC_DIR="${PETSC_DIR_DEFAULT}"
PETSC_ARCH="${PETSC_ARCH_DEFAULT}"
PETSC_BRANCH="${PETSC_BRANCH_DEFAULT}"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC_DEFAULT}"
INSTALL_DEPS=0
SKIP_PETSC=0
SKIP_PETSC4FOAM=0
JOBS="$(nproc)"

usage()
{
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --install-deps               Install Ubuntu packages via apt (uses sudo)
  --petsc-dir DIR              PETSc source/install directory (default: ${PETSC_DIR_DEFAULT})
  --petsc-arch ARCH            PETSc arch name (default: ${PETSC_ARCH_DEFAULT})
  --petsc-branch BRANCH        PETSc git branch/tag (default: ${PETSC_BRANCH_DEFAULT})
  --openfoam-bashrc FILE       OpenFOAM bashrc for petsc4Foam build (default: ${OPENFOAM_BASHRC_DEFAULT})
  -j, --jobs N                 Parallel build jobs (default: nproc)
  --skip-petsc                 Skip PETSc clone/configure/build
  --skip-petsc4foam            Skip petsc4Foam build
  -h, --help                   Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-deps) INSTALL_DEPS=1; shift ;;
        --petsc-dir) PETSC_DIR="$2"; shift 2 ;;
        --petsc-arch) PETSC_ARCH="$2"; shift 2 ;;
        --petsc-branch) PETSC_BRANCH="$2"; shift 2 ;;
        --openfoam-bashrc) OPENFOAM_BASHRC="$2"; shift 2 ;;
        -j|--jobs) JOBS="$2"; shift 2 ;;
        --skip-petsc) SKIP_PETSC=1; shift ;;
        --skip-petsc4foam) SKIP_PETSC4FOAM=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 2
            ;;
    esac
done

need_cmd()
{
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "Error: missing command '$1'" >&2
        exit 1
    fi
}

if [[ "${INSTALL_DEPS}" -eq 1 ]]; then
    need_cmd sudo
    need_cmd apt-get
    sudo apt-get update
    sudo apt-get install -y \
        build-essential cmake gfortran git pkg-config python3 \
        openmpi-bin libopenmpi-dev m4 bison flex
fi

need_cmd git
need_cmd python3
need_cmd make

if [[ "${SKIP_PETSC}" -eq 0 ]]; then
    if [[ ! -d "${PETSC_DIR}/.git" ]]; then
        git clone -b "${PETSC_BRANCH}" https://gitlab.com/petsc/petsc.git "${PETSC_DIR}"
    else
        git -C "${PETSC_DIR}" fetch --tags origin
        git -C "${PETSC_DIR}" checkout "${PETSC_BRANCH}"
        git -C "${PETSC_DIR}" pull --ff-only
    fi

    cd "${PETSC_DIR}"

    ./configure \
        --PETSC_ARCH="${PETSC_ARCH}" \
        --with-debugging=0 \
        --with-precision=double \
        --with-64-bit-indices=0 \
        --COPTFLAGS=-O3 \
        --CXXOPTFLAGS=-O3 \
        --FOPTFLAGS=-O3 \
        --with-mpi=1 \
        --download-f2cblaslapack \
        --download-hypre \
        --download-mumps \
        --download-superlu_dist \
        --download-scalapack \
        --download-metis \
        --download-parmetis

    make -j"${JOBS}" all
fi

if [[ "${SKIP_PETSC4FOAM}" -eq 0 ]]; then
    if [[ ! -f "${OPENFOAM_BASHRC}" ]]; then
        echo "Error: OpenFOAM bashrc not found: ${OPENFOAM_BASHRC}" >&2
        echo "Use --openfoam-bashrc to set the correct path." >&2
        exit 1
    fi

    if [[ ! -x "${SCRIPT_DIR}/external/petsc4Foam/Allwmake" ]]; then
        echo "Error: petsc4Foam Allwmake not found at ${SCRIPT_DIR}/external/petsc4Foam/Allwmake" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "${OPENFOAM_BASHRC}"
    export PETSC_DIR
    export PETSC_ARCH
    export PETSC_ARCH_PATH="${PETSC_DIR}/${PETSC_ARCH}"

    cd "${SCRIPT_DIR}/external/petsc4Foam"
    ./Allwmake
fi

cat <<EOF
Done.
PETSC_DIR=${PETSC_DIR}
PETSC_ARCH=${PETSC_ARCH}
PETSC_ARCH_PATH=${PETSC_DIR}/${PETSC_ARCH}
EOF
