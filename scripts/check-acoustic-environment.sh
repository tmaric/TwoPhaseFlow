#!/usr/bin/env bash
set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname "${script_dir}")"
expected_openfoam="v2606"
failures=0

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--expected-openfoam VERSION]

Checks the sourced OpenFOAM/TwoPhaseFlow environment and acoustic executables.
The check is read-only and does not build software or run simulations.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --expected-openfoam)
            if [[ $# -lt 2 ]]; then
                echo "Error: --expected-openfoam requires a version argument" >&2
                exit 2
            fi
            expected_openfoam="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option '$1'" >&2
            usage >&2
            exit 2
            ;;
    esac
done

pass()
{
    printf 'PASS  %s\n' "$1"
}

fail()
{
    printf 'FAIL  %s\n' "$1" >&2
    failures=$((failures + 1))
}

check_command()
{
    local command_name="$1"
    if command -v "${command_name}" >/dev/null 2>&1; then
        pass "command ${command_name}: $(command -v "${command_name}")"
    else
        fail "command ${command_name} is unavailable"
    fi
}

report_git_revision()
{
    local name="$1"
    local directory="$2"

    if ! command -v git >/dev/null 2>&1; then
        return
    fi

    if [[ -d "${directory}/.git" ]]; then
        pass "${name} source revision $(git -C "${directory}" rev-parse --short=12 HEAD)"
    else
        fail "${name} Git checkout missing: ${directory}"
    fi
}

if [[ "${WM_PROJECT_VERSION:-}" == "${expected_openfoam}" ]]; then
    pass "OpenFOAM version ${WM_PROJECT_VERSION}"
else
    fail "OpenFOAM version is '${WM_PROJECT_VERSION:-not sourced}', expected '${expected_openfoam}'"
fi

if [[ "${TPF_PROJECT_DIR:-}" == "${project_dir}" ]]; then
    pass "TwoPhaseFlow environment points to ${project_dir}"
else
    fail "source ${project_dir}/scripts/bashrc after sourcing OpenFOAM"
fi

for command_name in git ldd wmake mpirun python3 cartesian2DMesh; do
    check_command "${command_name}"
done

petsc_dir="${PETSC_DIR:-${project_dir}/petsc}"
petsc_arch="${PETSC_ARCH:-arch-linux-c-opt}"

report_git_revision "PETSc" "${petsc_dir}"
report_git_revision "petsc4Foam" "${project_dir}/external/petsc4Foam"

if [[ -f "${petsc_dir}/include/petsc.h" ]]; then
    pass "PETSc headers in ${petsc_dir}"
else
    fail "PETSc header missing: ${petsc_dir}/include/petsc.h"
fi

if [[ -f "${petsc_dir}/${petsc_arch}/include/petscconf.h" ]]; then
    pass "PETSc architecture ${petsc_arch}"
else
    fail "PETSc architecture header missing for ${petsc_arch}"
fi

if compgen -G "${petsc_dir}/${petsc_arch}/lib/libpetsc.*" >/dev/null; then
    pass "PETSc library in ${petsc_dir}/${petsc_arch}/lib"
else
    fail "PETSc library missing in ${petsc_dir}/${petsc_arch}/lib"
fi

if [[ -z "${FOAM_USER_APPBIN:-}" ]]; then
    fail "FOAM_USER_APPBIN is not defined"
else
    for executable in \
        setAlphaField \
        setPMLFields \
        acousticHelmholtzFoam \
        acousticWaveFoam \
        interFALFlow
    do
        executable_path="${FOAM_USER_APPBIN}/${executable}"
        if [[ ! -x "${executable_path}" ]]; then
            fail "executable missing: ${executable_path}"
            continue
        fi

        if command -v ldd >/dev/null 2>&1; then
            missing_libraries="$(ldd "${executable_path}" 2>/dev/null | awk '/not found/{print $1}')"
            if [[ -n "${missing_libraries}" ]]; then
                fail "${executable} has unresolved libraries: ${missing_libraries//$'\n'/, }"
            else
                pass "executable ${executable}"
            fi
        fi
    done
fi

if (( failures > 0 )); then
    printf '\nEnvironment check failed with %d issue(s).\n' "${failures}" >&2
    exit 1
fi

printf '\nAcoustic environment check passed.\n'
