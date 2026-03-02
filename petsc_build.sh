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
PETSC_DIR_DEFAULT="${SCRIPT_DIR}/petsc"
PETSC_ARCH_DEFAULT="arch-linux-c-opt"
PETSC_BRANCH_DEFAULT="release"
OPENFOAM_BASHRC_DEFAULT="$HOME/openfoam/etc/bashrc"
TWOPHASEFLOW_BASHRC_DEFAULT="${SCRIPT_DIR}/scripts/bashrc"
PETSC4FOAM_REPO_DEFAULT="https://gitlab.com/petsc/petsc4foam.git"

PETSC_DIR="${PETSC_DIR_DEFAULT}"
PETSC_ARCH="${PETSC_ARCH_DEFAULT}"
PETSC_BRANCH="${PETSC_BRANCH_DEFAULT}"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC_DEFAULT}"
TWOPHASEFLOW_BASHRC="${TWOPHASEFLOW_BASHRC_DEFAULT}"
PETSC4FOAM_REPO="${PETSC4FOAM_REPO_DEFAULT}"
INSTALL_DEPS=0
SKIP_PETSC=0
SKIP_PETSC4FOAM=0
JOBS="$(nproc)"
DOWNLOAD_CMAKE=0

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
  --twophaseflow-bashrc FILE   TwoPhaseFlow bashrc to source (default: ${TWOPHASEFLOW_BASHRC_DEFAULT})
  --petsc4foam-repo URL        petsc4Foam git repo (default: ${PETSC4FOAM_REPO_DEFAULT})
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
        --twophaseflow-bashrc) TWOPHASEFLOW_BASHRC="$2"; shift 2 ;;
        --petsc4foam-repo) PETSC4FOAM_REPO="$2"; shift 2 ;;
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

cmake_version_ok()
{
    local required="3.26.0"
    local current
    if ! command -v cmake >/dev/null 2>&1; then
        return 1
    fi
    current="$(cmake --version | head -n 1 | awk '{print $3}')"
    # Use sort -V for semver-like comparison
    if [[ "$(printf '%s\n' "${required}" "${current}" | sort -V | head -n 1)" == "${required}" ]]; then
        return 0
    fi
    return 1
}

if [[ "${SKIP_PETSC}" -eq 0 ]]; then
    if ! cmake_version_ok; then
        DOWNLOAD_CMAKE=1
        echo "CMake >= 3.26.0 not found; PETSc will download its own CMake." >&2
    fi

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
        --download-parmetis \
        $( [[ "${DOWNLOAD_CMAKE}" -eq 1 ]] && echo --download-cmake )

    make -j"${JOBS}" all
fi

if [[ "${SKIP_PETSC4FOAM}" -eq 0 ]]; then
    PETSC_DIR_SCRIPT="${PETSC_DIR}"
    PETSC_ARCH_SCRIPT="${PETSC_ARCH}"

    if [[ ! -f "${TWOPHASEFLOW_BASHRC}" ]]; then
        echo "Error: TwoPhaseFlow bashrc not found: ${TWOPHASEFLOW_BASHRC}" >&2
        echo "Use --twophaseflow-bashrc to set the correct path." >&2
        exit 1
    fi

    if [[ ! -f "${OPENFOAM_BASHRC}" ]]; then
        echo "Error: OpenFOAM bashrc not found: ${OPENFOAM_BASHRC}" >&2
        echo "Use --openfoam-bashrc to set the correct path." >&2
        exit 1
    fi

    if [[ ! -d "${SCRIPT_DIR}/external/petsc4Foam/.git" ]]; then
        mkdir -p "${SCRIPT_DIR}/external"
        git clone "${PETSC4FOAM_REPO}" "${SCRIPT_DIR}/external/petsc4Foam"
    fi
    if [[ ! -x "${SCRIPT_DIR}/external/petsc4Foam/Allwmake" ]]; then
        echo "Error: petsc4Foam Allwmake not found at ${SCRIPT_DIR}/external/petsc4Foam/Allwmake" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    set +eu
    source "${TWOPHASEFLOW_BASHRC}"
    source "${OPENFOAM_BASHRC}"
    set -eu
    PETSC_DIR="${PETSC_DIR_SCRIPT}"
    PETSC_ARCH="${PETSC_ARCH_SCRIPT}"
    export PETSC_DIR
    export PETSC_ARCH
    export PETSC_ARCH_PATH="${PETSC_DIR_SCRIPT}"
    if [[ ! -f "${PETSC_DIR}/include/petsc.h" ]]; then
        echo "Error: PETSc header not found at ${PETSC_DIR}/include/petsc.h" >&2
        echo "Verify PETSC_DIR is correct and PETSc is built." >&2
        exit 1
    fi
    if [[ ! -f "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" ]]; then
        echo "Error: PETSc arch header not found at ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" >&2
        echo "Verify PETSC_ARCH is correct and PETSc is built." >&2
        exit 1
    fi

    cd "${SCRIPT_DIR}/external/petsc4Foam"
    ./Allwclean
    ./Allwmake
fi

cat <<EOF
Done.
PETSC_DIR=${PETSC_DIR}
PETSC_ARCH=${PETSC_ARCH}
PETSC_ARCH_PATH=${PETSC_ARCH_PATH}
EOF
