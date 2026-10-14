"""The rows of the advection benchmarks (shared by the genCases.py of
planeTranslation3D and sphereTranslation3D).

A row is (RECONSCHEME, SETALPHAMETHOD, DAVOFMETHOD, NORMALMODEL, APPLYTO,
SLIVERTOL, MAXANGLE, NONORMALFALLBACK), the keys of system/simulationParameter
that select the reconstruction scheme, the initial state and the davof
scheme's normal model (davofNormalModel.H: normalModel, applyTo all | check,
sliverTol, maxAngle, noNormalFallback keep | ringArea):

    davof-<method>              isoAdvection + davof on the initDavofFields
                                state <method>; every interface cell takes
                                the DAVOF normal m_c (normalModel cellLocal)
    davof-<method>-Youngs<c>    the mixed normal: sliver cells
                                (|m_c| < c h^2) take the Youngs normal
                                fvc::grad(alpha), the others m_c
                                (youngs, check, sliverTol c, maxAngle 180:
                                the slivers-only rule)
    davof-<method>-sliverPlicRDF<c>
                                sliver cells take the plicRDF normal
                                (plicRDF, check, sliverTol c, maxAngle 180)
    davof-<method>-exact        control: every interface cell takes the exact
                                normal of the translated function, alphaf is
                                still transported by the moved plane
                                (exact, all)
    davof-<method>-plicRDF      control: every interface cell takes the
                                plicRDF normal of the same alpha (plicRDF, all);
                                the cells the helper plicRDF does not mark
                                keep m_c (noNormalFallback keep)
    davof-plicRDFnormal-ring    the same control with the ring normal in the
                                cells without a plicRDF normal (plicRDF, all,
                                0, 15, noNormalFallback ringArea), on the
                                first state construction given; the name of
                                leia's tpf_rows.py
    ringArea                    every interface cell takes the area-weighted
                                ring normal sum_k m_k/|sum_k m_k| (ringArea,
                                all), on the first state construction given
    ringArea-check              the ring normal in the slivers (|m_c| < 0.1 h^2)
                                and the cells deviating from it by more than
                                15 degrees (ringArea, check, 0.1, 15)
    ringFit, ringFit-check      the same with the ring quadric fit (ringFit)
    plicRDF, isoAlpha           the comparators on the setAlphaField state

The default set is the cluster list: the six davof variants (cellLocal, Youngs
0.01, 0.1, 1, exact, plicRDF), the four ring rows, plus the two comparators;
davof-plicRDFnormal-ring is generated on request (--rows).
"""

COMPARATORS = {
    "plicRDF": {"RECONSCHEME": "plicRDF", "SETALPHAMETHOD": "setAlphaField"},
    "isoAlpha": {"RECONSCHEME": "isoAlpha", "SETALPHAMETHOD": "setAlphaField"},
}

DEFAULT_YOUNGS_TOLS = [0.01, 0.1, 1.0]

# the ring rows: (normalModel, applyTo, sliverTol, maxAngle, noNormalFallback)
RING_ROWS = {
    "ringArea": ("ringArea", "all", "0", "15", "keep"),
    "ringArea-check": ("ringArea", "check", "0.1", "15", "keep"),
    "ringFit": ("ringFit", "all", "0", "15", "keep"),
    "ringFit-check": ("ringFit", "check", "0.1", "15", "keep"),
}

# the named controls generated on request (genCases.py --rows), the same tuple
CONTROL_ROWS = {
    "davof-plicRDFnormal-ring": ("plicRDF", "all", "0", "15", "ringArea"),
}

NAMED_ROWS = {**RING_ROWS, **CONTROL_ROWS}


def fmt(c):
    """0.01 -> '0.01', 1.0 -> '1'."""
    s = f"{c:g}"
    return s


def davof_base(method):
    return {"RECONSCHEME": "davof", "SETALPHAMETHOD": "davof", "DAVOFMETHOD": method,
            "NORMALMODEL": "cellLocal", "APPLYTO": "all", "SLIVERTOL": "0", "MAXANGLE": "15",
            "NONORMALFALLBACK": "keep"}


def davof_rows(method, youngs_tols=DEFAULT_YOUNGS_TOLS, sliver_plicrdf_tols=(),
               controls=("exact", "plicRDF")):
    """The davof variants on the state construction `method`, in the cluster order."""
    base = davof_base(method)
    rows = {f"davof-{method}": dict(base)}
    for c in youngs_tols:
        rows[f"davof-{method}-Youngs{fmt(c)}"] = dict(
            base, NORMALMODEL="youngs", APPLYTO="check", SLIVERTOL=fmt(c), MAXANGLE="180")
    for c in sliver_plicrdf_tols:
        rows[f"davof-{method}-sliverPlicRDF{fmt(c)}"] = dict(
            base, NORMALMODEL="plicRDF", APPLYTO="check", SLIVERTOL=fmt(c), MAXANGLE="180")
    for ctrl in controls:
        rows[f"davof-{method}-{ctrl}"] = dict(base, NORMALMODEL=ctrl, APPLYTO="all")
    return rows


def ring_rows(method, names=tuple(RING_ROWS)):
    """The named rows (NAMED_ROWS: the ring rows and the controls) on the
    state construction `method`."""
    base = davof_base(method)
    rows = {}
    for name in names:
        model, apply_to, tol, angle, fallback = NAMED_ROWS[name]
        rows[name] = dict(base, NORMALMODEL=model, APPLYTO=apply_to, SLIVERTOL=tol, MAXANGLE=angle,
                          NONORMALFALLBACK=fallback)
    return rows


def all_rows(methods, youngs_tols=DEFAULT_YOUNGS_TOLS, sliver_plicrdf_tols=(),
             controls=("exact", "plicRDF"), comparators=True, rings=True, extras=()):
    """Every row of the given state constructions, the ring rows and the
    named `extras` (CONTROL_ROWS) on the first one, then the comparators."""
    rows = {}
    for m in methods:
        rows.update(davof_rows(m, youngs_tols, sliver_plicrdf_tols, controls))
    if methods:
        names = (tuple(RING_ROWS) if rings else ()) + tuple(extras)
        rows.update(ring_rows(methods[0], names))
    if comparators:
        rows.update(COMPARATORS)
    return rows


def complete(values):
    """Fill the davof keys the comparator rows leave out (the
    simulationParameter needs every key; the other schemes ignore them)."""
    out = {"DAVOFMETHOD": "detrixheAslam", "NORMALMODEL": "cellLocal", "APPLYTO": "all",
           "SLIVERTOL": "0", "MAXANGLE": "15", "NONORMALFALLBACK": "keep"}
    out.update(values)
    return out
