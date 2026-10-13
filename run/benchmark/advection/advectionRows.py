"""The rows of the advection benchmarks (shared by the genCases.py of
planeTranslation3D and sphereTranslation3D).

A row is (RECONSCHEME, SETALPHAMETHOD, DAVOFMETHOD, SLIVERNORMAL, SLIVERTOL,
NORMALSOURCE), the keys of system/simulationParameter that select the
reconstruction scheme, the initial state and the davof scheme's normal
controls (davof.H):

    davof-<method>              isoAdvection + davof on the initDavofFields
                                state <method>; every interface cell takes
                                the DAVOF normal m_c (sliverNormal none)
    davof-<method>-Youngs<c>    the mixed normal: sliver cells
                                (|m_c| < c h^2) take the Youngs normal
                                fvc::grad(alpha), the others m_c
    davof-<method>-sliverPlicRDF<c>
                                sliver cells take the plicRDF normal
    davof-<method>-exact        control: every interface cell takes the exact
                                normal of the translated function, alphaf is
                                still transported by the moved plane
    davof-<method>-plicRDF      control: every interface cell takes the
                                plicRDF normal of the same alpha
    plicRDF, isoAlpha           the comparators on the setAlphaField state

The default set is the cluster list: the six davof variants (none, Youngs
0.01, 0.1, 1, exact, plicRDF) plus the two comparators.
"""

COMPARATORS = {
    "plicRDF": {"RECONSCHEME": "plicRDF", "SETALPHAMETHOD": "setAlphaField"},
    "isoAlpha": {"RECONSCHEME": "isoAlpha", "SETALPHAMETHOD": "setAlphaField"},
}

DEFAULT_YOUNGS_TOLS = [0.01, 0.1, 1.0]


def fmt(c):
    """0.01 -> '0.01', 1.0 -> '1'."""
    s = f"{c:g}"
    return s


def davof_rows(method, youngs_tols=DEFAULT_YOUNGS_TOLS, sliver_plicrdf_tols=(),
               controls=("exact", "plicRDF")):
    """The davof variants on the state construction `method`, in the cluster order."""
    base = {"RECONSCHEME": "davof", "SETALPHAMETHOD": "davof", "DAVOFMETHOD": method,
            "SLIVERNORMAL": "none", "SLIVERTOL": "0", "NORMALSOURCE": "davof"}
    rows = {f"davof-{method}": dict(base)}
    for c in youngs_tols:
        rows[f"davof-{method}-Youngs{fmt(c)}"] = dict(base, SLIVERNORMAL="Youngs", SLIVERTOL=fmt(c))
    for c in sliver_plicrdf_tols:
        rows[f"davof-{method}-sliverPlicRDF{fmt(c)}"] = dict(base, SLIVERNORMAL="plicRDF", SLIVERTOL=fmt(c))
    for ctrl in controls:
        rows[f"davof-{method}-{ctrl}"] = dict(base, NORMALSOURCE=ctrl)
    return rows


def all_rows(methods, youngs_tols=DEFAULT_YOUNGS_TOLS, sliver_plicrdf_tols=(),
             controls=("exact", "plicRDF"), comparators=True):
    """Every row of the given state constructions, then the comparators."""
    rows = {}
    for m in methods:
        rows.update(davof_rows(m, youngs_tols, sliver_plicrdf_tols, controls))
    if comparators:
        rows.update(COMPARATORS)
    return rows


def complete(values):
    """Fill the davof-control keys the comparator rows leave out (the
    simulationParameter needs every key; the other schemes ignore them)."""
    out = {"DAVOFMETHOD": "detrixheAslam", "SLIVERNORMAL": "none", "SLIVERTOL": "0",
           "NORMALSOURCE": "davof"}
    out.update(values)
    return out
