#!/usr/bin/env python3
"""Exact time-harmonic pressure for an axial point source and a fluid sphere.

Convention exp(-i omega t); outgoing waves use h_l = j_l + i y_l, as in
sphere_pressure of advanced.py. Axisymmetry axis = y; sphere of radius a
centred at the origin; point source at (0, d, 0) with d > a, normalised to
1 Pa (incident field) at the sphere centre. Pressure and rho^-1 dp/dr are
continuous at r = a.
"""
import numpy as np
from scipy.special import spherical_jn as jn, spherical_yn as yn, eval_legendre

def hn(l, x, derivative=False):
    return jn(l, x, derivative) + 1j*yn(l, x, derivative)

def coefficients(a, d, kg, kl, rhog, rhol, lmax, rigid=False):
    """Scattered (A_l) and interior (B_l) coefficients.

    Incident expansion for r < d:  S i kg sum (2l+1) j_l(kg r) h_l(kg d) P_l,
    with S = d exp(-i kg d) so that p_inc(0) = 1.
    """
    S = d*np.exp(-1j*kg*d)
    A = np.zeros(lmax+1, complex); B = np.zeros(lmax+1, complex)
    for l in range(lmax+1):
        c = S*1j*kg*(2*l+1)*hn(l, kg*d)
        if rigid:
            A[l] = -c*jn(l, kg*a, True)/hn(l, kg*a, True)
            continue
        M = np.array([[hn(l, kg*a), -jn(l, kl*a)],
                      [kg/rhog*hn(l, kg*a, True), -kl/rhol*jn(l, kl*a, True)]], complex)
        rhs = np.array([-c*jn(l, kg*a), -c*kg/rhog*jn(l, kg*a, True)], complex)
        A[l], B[l] = np.linalg.solve(M, rhs)
    return S, A, B

def pressure(points, a, d, kg, kl, rhog, rhol, lmax=60, rigid=False):
    """Total pressure; host (g) outside and drop (l) inside the sphere."""
    pts = np.asarray(points, float)
    r = np.linalg.norm(pts, axis=1)
    mu = np.divide(pts[:, 1], r, out=np.ones_like(r), where=r > 0)
    R = np.linalg.norm(pts - np.array([0.0, d, 0.0]), axis=1)
    S, A, B = coefficients(a, d, kg, kl, rhog, rhol, lmax, rigid)
    inside = r < a
    p = np.zeros(len(r), complex)
    ext = ~inside
    p[ext] = S*np.exp(1j*kg*R[ext])/R[ext]
    for l in range(lmax+1):
        P = eval_legendre(l, mu)
        p[ext] += A[l]*hn(l, kg*r[ext])*P[ext]
        if not rigid:
            p[inside] += B[l]*jn(l, kl*r[inside])*P[inside]
    if not np.isfinite(p).all():
        raise ValueError("non-finite reference pressure")
    return p

def self_check(a=1e-3, d=2.5e-3, kg=1000.0, cg=343.0, cl=1500.0, rhog=1.2, rhol=1000.0, lmax=60):
    """Interface conditions, incident expansion, homogeneous limit, truncation."""
    kl = kg*cg/cl
    mu = np.cos(np.linspace(0, np.pi, 73))
    S, A, B = coefficients(a, d, kg, kl, rhog, rhol, lmax)
    # Both sides at r = a with analytic radial derivatives; the exterior uses
    # the incident expansion, valid for r < d.
    pe = pi = ve = vi = 0
    for l in range(lmax+1):
        c = S*1j*kg*(2*l+1)*hn(l, kg*d)
        P = eval_legendre(l, mu)
        pe = pe + (c*jn(l, kg*a) + A[l]*hn(l, kg*a))*P
        ve = ve + kg/rhog*(c*jn(l, kg*a, True) + A[l]*hn(l, kg*a, True))*P
        pi = pi + B[l]*jn(l, kl*a)*P
        vi = vi + kl/rhol*B[l]*jn(l, kl*a, True)*P
    ring = np.c_[a*np.sqrt(1-mu**2), a*mu, 0*mu]
    R = np.linalg.norm(ring-np.array([0, d, 0]), axis=1)
    inc = sum(S*1j*kg*(2*l+1)*hn(l, kg*d)*jn(l, kg*a)*eval_legendre(l, mu) for l in range(lmax+1))
    out = {
        "pressureJump": float(np.max(abs(pe-pi))/np.max(abs(pe))),
        "velocityJump": float(np.max(abs(ve-vi))/np.max(abs(ve))),
        "incidentExpansion": float(np.max(abs(inc-S*np.exp(1j*kg*R)/R))/np.max(abs(S*np.exp(1j*kg*R)/R))),
    }
    # Homogeneous limit: no scatterer -> incident field only.
    pts = np.c_[np.linspace(0, 3*a, 20), np.linspace(-3*a, 1.5*a, 20), np.zeros(20)]
    ph = pressure(pts, a, d, kg, kg, rhog, rhog)
    R = np.linalg.norm(pts-np.array([0, d, 0]), axis=1)
    pinc = d*np.exp(-1j*kg*d)*np.exp(1j*kg*R)/R
    out["homogeneousLimit"] = float(np.max(abs(ph-pinc))/np.max(abs(pinc)))
    # Series truncation: lmax 60 vs 80 on a sample of exterior/interior points.
    q = np.c_[np.linspace(0, 4*a, 50), np.linspace(-4*a, 2.2*a, 50), np.zeros(50)]
    out["truncation60vs80"] = float(np.max(abs(pressure(q, a, d, kg, kl, rhog, rhol, 60) -
                                               pressure(q, a, d, kg, kl, rhog, rhol, 80))))
    return out

if __name__ == "__main__":
    import json
    print(json.dumps(self_check(), indent=2))
