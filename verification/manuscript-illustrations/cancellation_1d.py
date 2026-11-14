"""Exact check of the one-dimensional cancellation example (Section 2 of the manuscript).

    dx/dt = w(t),  w = v + dz/dt,  v(t) = U cos(pi t/T),  z(t) = beta t^2 (T-t)/T^2,

integrated with the midpoint rule x_{n+1} = x_n + h w(t_n + h/2), h = T/N, N even.
The script confirms with sympy the closed forms stated in the paper:

    retracing (beta = 0):  x_{N/2} = U h / (2 sin(pi h/(2T))),   x_N = 0,
    with the frame:        x_N = beta h^2 / (4T),

and the closing remark of Section 2: for z(t) = eps sin(2 pi t/T) the midpoint
samples of dz/dt sum to zero on every uniform grid with N >= 2 steps, so this
frame leaves the terminal error at zero although x(T-t) != x(t).

    python3 cancellation_1d.py
"""
import sympy as sp

t, h, T, U, beta = sp.symbols("t h T U beta", positive=True)
n, M = sp.symbols("n M", integer=True, nonnegative=True)

v = U * sp.cos(sp.pi * t / T)
z = beta * t**2 * (T - t) / T**2
assert sp.simplify(z.subs(t, 0)) == 0 and sp.simplify(z.subs(t, T)) == 0


def midpoint_sum(f, steps):
    """h * sum_{n=0}^{steps-1} f(t_n + h/2), summed in closed form in `steps`."""
    return sp.summation(h * f.subs(t, (n + sp.Rational(1, 2)) * h), (n, 0, steps - 1))


# cos part: with theta = pi h/T the midpoint terms telescope,
#   2 sin(theta/2) cos((n+1/2) theta) = sin((n+1) theta) - sin(n theta),
# so h sum_{n<M} U cos((n+1/2) theta) = U h sin(M theta) / (2 sin(theta/2)).
th = sp.symbols("theta", positive=True)
tele = 2 * sp.sin(th / 2) * sp.cos((n + sp.Rational(1, 2)) * th) - (sp.sin((n + 1) * th) - sp.sin(n * th))
print("telescoping identity residual  =", sp.simplify(sp.expand_trig(sp.expand(tele))))
x_M = U * h * sp.sin(M * th) / (2 * sp.sin(th / 2))          # x after M steps
x_M = x_M.subs(th, sp.pi * h / T)

# retracing motion, N = 2M steps, h = T/N: M = N/2 steps reach T/2, 2M reach T
half = sp.simplify(x_M.subs(M, T / (2 * h)))
full = sp.simplify(x_M.subs(M, T / h))
claim_half = U * h / (2 * sp.sin(sp.pi * h / (2 * T)))
print("x_{N/2} - U h/(2 sin(pi h/2T)) =", sp.simplify(half - claim_half))
print("x_N                           =", full)
err_half = sp.series(claim_half - U * T / sp.pi, h, 0, 4).removeO()
print("x_{N/2} - U T/pi              =", sp.simplify(err_half), "+ O(h^4)")

# frame part: the polynomial sum is exact in N; the cos part is the zero above
zdot = sp.diff(z, t)
frame = sp.simplify(midpoint_sum(zdot, 2 * M).subs(h, T / (2 * M)) - (z.subs(t, T) - z.subs(t, 0)))
print("x_N with frame - beta h^2/(4T)  =",
      sp.simplify(frame - (beta * h**2 / (4 * T)).subs(h, T / (2 * M))))

# spot check in exact arithmetic for a few N, independent of the closed forms
for N in (4, 8, 16, 64):
    hh = T / N
    s = sum(hh * (v + zdot).subs(t, (k + sp.Rational(1, 2)) * hh) for k in range(N))
    print(f"N = {N:3d}: x_N - beta h^2/(4T) =",
          sp.nsimplify(sp.simplify(s - beta * hh**2 / (4 * T))))

# sine frame z = eps sin(2 pi t/T): returning, but the midpoint samples of its
# velocity cancel over the period for every N >= 2 (odd N included)
eps = sp.symbols("epsilon", positive=True)
zs = eps * sp.sin(2 * sp.pi * t / T)
zsdot = sp.diff(zs, t)
assert sp.simplify(zs.subs(t, 0)) == 0 and sp.simplify(zs.subs(t, T)) == 0
print("sine frame: z(T/4) - z(3T/4)   =", sp.simplify(zs.subs(t, T / 4) - zs.subs(t, 3 * T / 4)))
for N in range(2, 13):
    hh = T / N
    s = sum(hh * zsdot.subs(t, (k + sp.Rational(1, 2)) * hh) for k in range(N))
    # the sum is eps * (exact trig sum); evaluate the trig sum to 60 digits
    s = sp.N(s.subs({T: 1, eps: 1}), 60)
    assert abs(s) < 1e-50, (N, s)
print("sine frame: midpoint sum of dz/dt = 0 (to 60 digits) for N = 2, ..., 12")
