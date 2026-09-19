"""Chapter 24, Section 5: regression calibration and SIMEX in the blood-pressure study.

Regressors (1, mean of two readings, age); the error variance of the mean reading is estimated
from the replicate pairs.

(a) Regression calibration: replace the mean reading by its moment-based best linear predictor
    of x given (reading, age). In the linear model this reproduces the method-of-moments estimate
    exactly.
(b) SIMEX: add extra error of variance zeta * sigma_ubar^2 for zeta in {0.5, 1, 1.5, 2}, average
    B = 200 naive fits at each zeta, extrapolate to zeta = -1 with a quadratic and with the
    rational function a + b / (c + zeta). The rational extrapolant is exact in the linear model;
    the quadratic one is not, and its large-sample limit is computed exactly.
(c) A pairs bootstrap standard error for the SIMEX estimate.
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch24", "simex", prefix="sx")

# <<setup>>
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x
# <</setup>>

# <<calibration>>
wbar = w.mean(axis=1)
s2_ubar = np.mean((w[:, 0] - w[:, 1]) ** 2) / 4               # error variance of the mean reading
W = np.column_stack([np.ones(n), wbar, age])
S = W.T @ W / n                                               # raw second moments of (1, wbar, age)
Suu = np.diag([0.0, s2_ubar, 0.0])

# calibration: best linear predictor of (1, x, age) from (1, wbar, age), using E(w x^T) = S - Suu
Gamma = np.linalg.solve(S, S - Suu)
x_hat = (W @ Gamma)[:, 1]                                     # predicted long-run pressure
b_rc = np.linalg.lstsq(np.column_stack([np.ones(n), x_hat, age]), y, rcond=None)[0]
b_mom = np.linalg.solve(S - Suu, W.T @ y / n)                 # method of moments, Section 24.3
print("regression calibration:", b_rc.round(4))
print("method of moments:     ", b_mom.round(4))
# <</calibration>>

assert np.allclose(b_rc, b_mom)
b_naive = np.linalg.lstsq(W, y, rcond=None)[0]
gen.num("s2_ubar", s2_ubar, 2)
gen.num("rc_x", b_rc[1], 4)
gen.num("rc_age", b_rc[2], 4)
gen.num("naive_x", b_naive[1], 4)
gen.num("naive_age", b_naive[2], 4)
gen.num("xhat_sd", x_hat.std(), 2)
gen.num("wbar_sd", wbar.std(), 2)


# <<simex>>
def naive_fits(wcol, age, y):
    """Least squares coefficients of y on (1, wcol, age) for a stack of remeasured columns."""
    m, k = wcol.shape                                         # m remeasured data sets of size k
    Xs = np.stack([np.ones((m, k)), wcol, np.broadcast_to(age, (m, k))], axis=2)
    XtX = np.einsum("mki,mkj->mij", Xs, Xs)
    Xty = np.einsum("mki,k->mi", Xs, y)
    return np.linalg.solve(XtX, Xty[:, :, None])[:, :, 0]


def simex(wbar, age, y, s2_ubar, zetas, B, rng):
    """Average naive coefficients after adding error of variance zeta * s2_ubar."""
    out = [naive_fits(wbar[None, :], age, y)[0]]
    for z in zetas[1:]:
        extra = rng.normal(0, np.sqrt(z * s2_ubar), (B, len(y)))
        out.append(naive_fits(wbar + extra, age, y).mean(axis=0))
    return np.array(out)                                      # one row per zeta


def rational(z, a, b, c):
    return a + b / (c + z)


zetas = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
G = simex(wbar, age, y, s2_ubar, zetas, 200, np.random.default_rng(2450))
est_quad, est_rat = [], []
for j in [1, 2]:                                              # coefficient of x, of age
    est_quad.append(np.polyval(np.polyfit(zetas, G[:, j], 2), -1.0))
    p, _ = curve_fit(rational, zetas, G[:, j], p0=[0.0, G[0, j], 1.5], maxfev=20000)
    est_rat.append(rational(-1.0, *p))
print("naive (zeta = 0):   ", G[0, 1:].round(4))
print("SIMEX, quadratic:   ", np.round(est_quad, 4))
print("SIMEX, rational:    ", np.round(est_rat, 4))
# <</simex>>

assert np.isclose(G[0, 1], b_naive[1]) and np.isclose(G[0, 2], b_naive[2])
assert G[0, 1] < est_quad[0] < est_rat[0]                     # quadratic undercorrects
assert abs(est_rat[0] - b_mom[1]) < 0.003
assert abs(est_rat[1] - b_mom[2]) < 0.003
for j, name in [(0, "x"), (1, "age")]:
    gen.num(f"quad_{name}", est_quad[j], 4)
    gen.num(f"rat_{name}", est_rat[j], 4)
gen.num("mom_x", b_mom[1], 4)
gen.num("mom_age", b_mom[2], 4)

# ---- exact large-sample analysis of the extrapolants -----------------------------------------
beta, sxz2, su2bar, pi = 0.08, 108.0, 81.0 / 2, 0.6
def g_x(z):
    return beta * sxz2 / (sxz2 + (1 + z) * su2bar)
def g_age(z):
    return (1 - sxz2 / (sxz2 + (1 + z) * su2bar)) * beta * pi
assert np.isclose(g_x(-1), beta) and np.isclose(g_age(-1), 0)
quad_lim_x = np.polyval(np.polyfit(zetas, g_x(zetas), 2), -1.0)
quad_lim_age = np.polyval(np.polyfit(zetas, g_age(zetas), 2), -1.0)
# reciprocal of the x coefficient's limit is linear in zeta
assert np.allclose(np.diff(1 / g_x(zetas), 2), 0)
gen.num("g0_x", g_x(0), 4)
gen.num("quad_lim_x", quad_lim_x, 4)
gen.num("quad_lim_age", quad_lim_age, 4)
gen.num("quad_lim_frac", (quad_lim_x - g_x(0)) / (beta - g_x(0)), 2)

# ---- pairs bootstrap standard error for SIMEX ------------------------------------------------
# <<bootstrap>>
def simex_rational_x(idx, rng):
    """Rational-extrapolant SIMEX estimate of the x coefficient on the cases idx."""
    wb, a, yy = w[idx].mean(axis=1), age[idx], y[idx]
    s2 = np.mean((w[idx, 0] - w[idx, 1]) ** 2) / 4            # re-estimate the error variance too
    Gb = simex(wb, a, yy, s2, zetas, 50, rng)
    p, _ = curve_fit(rational, zetas, Gb[:, 1], p0=[0.0, Gb[0, 1], 1.5], maxfev=20000)
    return rational(-1.0, *p)


boot_rng = np.random.default_rng(2451)
boot = [simex_rational_x(boot_rng.integers(0, n, n), boot_rng) for _ in range(200)]
print(f"bootstrap standard error of the SIMEX slope: {np.std(boot, ddof=1):.4f}")
# <</bootstrap>>
se_boot = np.std(boot, ddof=1)
assert 0.005 < se_boot < 0.011
gen.num("se_boot", se_boot, 4)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
zz = np.linspace(-1, 2, 200)
for ax, j, truth, title, g in [(axes[0], 1, 0.08, "(a) coefficient of blood pressure", g_x),
                               (axes[1], 2, 0.0, "(b) coefficient of age", g_age)]:
    cq = np.polyfit(zetas, G[:, j], 2)
    pr, _ = curve_fit(rational, zetas, G[:, j], p0=[0.0, G[0, j], 1.5], maxfev=20000)
    ax.plot(zz, np.polyval(cq, zz), color=COLORS["second"], linewidth=1.0, label="quadratic")
    ax.plot(zz, rational(zz, *pr), color=COLORS["accent"], linewidth=1.0, label=r"$a+b/(c+\zeta)$")
    ax.plot(zz, g(zz), color=COLORS["muted"], linestyle=":", linewidth=1.0, label="large-sample limit")
    ax.scatter(zetas, G[:, j], s=16, color=COLORS["ink"], zorder=3, label="simulated")
    ax.scatter([-1, -1], [np.polyval(cq, -1), rational(-1, *pr)], s=22, marker="D",
               color=[COLORS["second"], COLORS["accent"]], zorder=3)
    ax.axhline(truth, color=COLORS["ink"], linestyle="--", linewidth=0.7)
    ax.axvline(0, color=COLORS["grid"], linewidth=0.7, zorder=0)
    ax.set_xlabel(r"added error, $\zeta$")
    ax.set_title(title)
axes[0].legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch24", "simex"))
