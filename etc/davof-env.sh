# etc/davof-env.sh -- source AFTER OpenFOAM's etc/bashrc, in every shell that builds or
# runs a TwoPhaseFlow binary of THIS clone:
#
#     source $HOME/OpenFOAM/OpenFOAM-v2606/etc/bashrc && . ./etc/davof-env.sh
#
# Every Make/files in this repository installs through FOAM_USER_APPBIN / FOAM_USER_LIBBIN,
# i.e. by default into the account-wide $HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION, where
# every other TwoPhaseFlow checkout on the same account installs too. This file points them
# at <clone>/platforms/<WM_OPTIONS>/{bin,lib} (git-ignored) instead, so two clones never
# share binaries, and REMOVES the account default from PATH and LD_LIBRARY_PATH: it may hold
# binaries of other projects (2026-09-28: a libVoF.so, libalphaFieldFunctions.so and
# setAlphaField of an unrelated build sat in $HOME/OpenFOAM/tmaric-v2606) that would shadow
# this clone's or OpenFOAM's own. Nothing in that directory is touched. OpenFOAM's
# installation ($WM_PROJECT_DIR) and ThirdParty stay on the paths. Idempotent.
# Sourced by Allwmake and Allwclean; leia's etc/leia-env.sh appends this clone's
# directories itself (LEIA_TPF_DIR), so leia's studies need no reference to this file.
if [ -z "$WM_OPTIONS" ] || [ -z "$WM_PROJECT_VERSION" ]; then
    echo "davof-env.sh: WM_OPTIONS is empty; source OpenFOAM's etc/bashrc first" >&2
    return 1 2>/dev/null || exit 1
fi
TPF_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export TPF_ROOT

# The account default as etc/bashrc set it, captured BEFORE it is overridden. When an
# earlier source (this file, or a leia clone's etc/leia-env.sh) already replaced it by a
# clone root, fall back to the bashrc's own formula
# (etc/bashrc: WM_PROJECT_USER_DIR="$HOME/$WM_PROJECT/${USER:-user}-$WM_PROJECT_VERSION").
_tpf_acct="${WM_PROJECT_USER_DIR:-}"
if [ -z "$_tpf_acct" ] || [ -f "$_tpf_acct/etc/davof-env.sh" ] || [ -f "$_tpf_acct/etc/leia-env.sh" ]; then
    _tpf_acct="$HOME/${WM_PROJECT:-OpenFOAM}/${USER:-user}-$WM_PROJECT_VERSION"
fi

export WM_PROJECT_USER_DIR="$TPF_ROOT"
export FOAM_USER_APPBIN="$TPF_ROOT/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$TPF_ROOT/platforms/$WM_OPTIONS/lib"

# v2512 and v2606 share WM_OPTIONS=linux64GccDPInt32Opt, so a platforms/ built against one
# is loaded under the other without any error from the loader. Allwmake writes the stamp,
# Allwclean removes the directory; this file refuses to run across versions.
_tpf_stamp="$TPF_ROOT/platforms/$WM_OPTIONS/.openfoam-version"
if [ -f "$_tpf_stamp" ] && [ "$(cat "$_tpf_stamp")" != "$WM_PROJECT_VERSION" ]; then
    echo "davof-env.sh: $TPF_ROOT/platforms/$WM_OPTIONS was built with OpenFOAM-$(cat "$_tpf_stamp"), but OpenFOAM-$WM_PROJECT_VERSION is sourced; run ./Allwclean && ./Allwmake under the version you want" >&2
    unset _tpf_stamp _tpf_acct
    return 1 2>/dev/null || exit 1
fi
unset _tpf_stamp

# Drop the account default's platforms/ and any earlier copy of this clone's, keep
# $WM_PROJECT_DIR/ and $WM_THIRD_PARTY_DIR/ entries, then put this clone first.
_tpf_strip() {
    echo "$1" | tr ':' '\n' | awk -v of="$WM_PROJECT_DIR/" -v tp="${WM_THIRD_PARTY_DIR:-/nonexistent}/" \
        -v mine="$TPF_ROOT/platforms/" -v acct="$_tpf_acct/platforms/" '
        length($0) == 0          { next }
        index($0, of) == 1       { print; next }
        index($0, tp) == 1       { print; next }
        index($0, mine) == 1     { next }
        index($0, acct) == 1     { next }
        { print }' | paste -sd:
}
export PATH="$FOAM_USER_APPBIN:$(_tpf_strip "$PATH")"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$(_tpf_strip "$LD_LIBRARY_PATH")"
unset -f _tpf_strip
unset _tpf_acct

# Eigen 3 for the vendored overlap library (src/VoF/davof/thirdParty/overlap, the exact
# sphere-tetrahedron overlap of Strobl, Formella and Poeschel 2016, doi:10.1016/j.jcp.2016.02.003):
# src/VoF/Make/options adds -I$(EIGEN_INC). Keep an EIGEN_INC set by a module (Lichtenberg:
# eigen/5.0.1 puts its include on CPATH as well); on the laptop the headers live in
# $HOME/OpenFOAM/thirdparty/eigen/eigen-3.4.0 (downloaded 2026-10-24, MPL2).
if [ -z "${EIGEN_INC:-}" ]; then
    if [ -d "$HOME/OpenFOAM/thirdparty/eigen/eigen-3.4.0" ]; then
        export EIGEN_INC="$HOME/OpenFOAM/thirdparty/eigen/eigen-3.4.0"
    elif [ -d /usr/include/eigen3 ]; then
        export EIGEN_INC=/usr/include/eigen3
    else
        export EIGEN_INC="${CPATH%%:*}"
    fi
fi
