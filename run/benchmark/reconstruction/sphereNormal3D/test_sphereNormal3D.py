"""Sphere-normal reconstruction ladder as a pytest (oftest) test.

    pytest run/benchmark/reconstruction/sphereNormal3D            # N = 20, 40, 80
    pytest run/benchmark/reconstruction/sphereNormal3D -m slow    # + N = 160

Every parametrised test rewrites system/simulationParameter of THIS directory
(oftest.Case_modifiers), runs Allrun in place, and reads the appended columns of
postProcessing/reconstructionError/0/error.dat (copied to error.dat by Allrun):
LNormalDiff1/2/Inf = mean / rms / max of |n - n_exact| over the interface cells,
areaSum = sum of the PLIC areas, areaExact = 4 pi R^2, nInterface. The last test
fits the least-squares order of LNormalDiff2 against h per scheme over the rungs
that ran and writes results_pytest_sphereNormal3D.csv.
"""
import os

import numpy as np
import pandas as pd
import pytest
import oftest
from oftest import run_reset_case  # noqa: F401  (the fixture)

HERE = os.path.dirname(os.path.abspath(__file__))
L, R = 0.004, 0.001


def simMod(n, scheme, setAlpha):
    filemod = {
        "system/simulationParameter": [
            ("N", n),
            ("RECONSCHEME", scheme),
            ("SETALPHAMETHOD", setAlpha),
        ]
    }
    meta = {"scheme": scheme, "setAlpha": setAlpha, "N": n}
    return oftest.Case_modifiers(filemod, HERE, meta)


def _param(n, scheme, setAlpha, slow=False):
    p = simMod(n, scheme, setAlpha)
    return pytest.param(p, marks=pytest.mark.slow) if slow else p


parameters = []
ids = []
for scheme in ("plicRDF", "gradAlpha"):
    for n in (20, 40, 80, 160):
        parameters.append(_param(n, scheme, "cutCellIso", slow=(n == 160)))
        ids.append(f"{scheme}-cutCellIso-N{n}")

results = {"scheme": [], "setAlpha": [], "N": [], "h": [], "LNormalDiff1": [],
           "LNormalDiff2": [], "LNormalDiffInf": [], "areaSum": [], "areaExact": [],
           "nInterface": [], "runtime": []}


def read_error_dat(path):
    names, last = None, None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            if s[1:].split() and s[1:].split()[0] == "Time":
                names = s[1:].split()
            continue
        last = s.split()
    return dict(zip(names, last))


@pytest.mark.parametrize("run_reset_case", parameters, indirect=["run_reset_case"], ids=ids)
def test_reconstruction(run_reset_case):
    log = oftest.path_log()
    assert oftest.case_status(log) == "completed"
    rec = read_error_dat(os.path.join(HERE, "error.dat"))
    meta = run_reset_case.meta_data
    n = meta["N"]
    results["scheme"].append(meta["scheme"])
    results["setAlpha"].append(meta["setAlpha"])
    results["N"].append(n)
    results["h"].append(L/n)
    for k in ("LNormalDiff1", "LNormalDiff2", "LNormalDiffInf", "areaSum", "areaExact"):
        results[k].append(float(rec[k]))
    results["nInterface"].append(int(float(rec["nInterface"])))
    results["runtime"].append(oftest.run_time(log)[0])
    # sanity: a sphere of R/h >= 5 has an interface and a finite normal error
    assert results["nInterface"][-1] > 0
    assert 0 < float(rec["LNormalDiff2"]) < 1


# Least-squares order of LNormalDiff2 each scheme must reach. MEASURED 2026-09-28
# (N = 20/40/80, cutCellIso): plicRDF 1.07; gradAlpha 0.11 -- the centred-gradient
# (Youngs-type) normal of a sharp indicator function does not converge under
# refinement (Pilliod & Puckett, JCP 199, 2004), its rms error stays at 0.066, so
# the ladder records it and asserts completion only (None).
MIN_ORDER = {"plicRDF": 0.9, "gradAlpha": None}


def test_orders():
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(HERE, "results_pytest_sphereNormal3D.csv"), index=False)
    for scheme, g in df.groupby("scheme"):
        g = g.sort_values("h", ascending=False)
        if len(g) < 3:
            continue
        p = np.polyfit(np.log(g["h"]), np.log(g["LNormalDiff2"]), 1)[0]
        print(f"{scheme}: least-squares order of LNormalDiff2 over N={list(g['N'])}: {p:.2f}")
        pmin = MIN_ORDER.get(scheme)
        if pmin is not None:
            assert p >= pmin, f"{scheme}: order {p:.2f} < {pmin}"
