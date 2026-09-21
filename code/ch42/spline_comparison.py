"""Chapter 42, Section 4: five trend bases for the same data, fitted by ordinary least squares.

Monthly mean CO2 at Mauna Loa, 1958-2001 (public domain, statsmodels.datasets.co2).
The four seasonal columns of Example 5.3 are in every model; only the trend basis
changes. The last block refits on 1958-1990 and predicts 1991-2001.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<bases>>
import numpy as np
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]
monthly = co2.resample("MS").mean().dropna()
y = monthly.to_numpy()
t = (monthly.index.year + (monthly.index.month - 0.5) / 12).to_numpy()
n = len(y)
lo, hi = t.min(), t.max()
u = lambda s: 2 * (s - lo) / (hi - lo) - 1

K = 6
knots = np.quantile(t, np.arange(1, K + 1) / (K + 1))      # quantile knots
print("knots at", np.round(knots, 2))


def polynomial(s, d):
    return np.vander(u(s), d + 1, increasing=True)


def step(s):                                               # piecewise constant, d = 0
    piece = np.digitize(s, knots)
    return np.column_stack([(piece == k).astype(float) for k in range(K + 1)])


def broken_line(s):                                        # continuous piecewise linear
    return np.column_stack([np.ones_like(s), u(s)]
                           + [np.clip(u(s) - u(k), 0.0, None) for k in knots])


def cubic_spline(s):                                       # truncated power basis, d = 3
    return np.column_stack([u(s) ** j for j in range(4)]
                           + [np.clip(u(s) - u(k), 0.0, None) ** 3 for k in knots])


def with_season(f):
    return lambda s: np.column_stack([f(s), np.cos(2 * np.pi * s), np.sin(2 * np.pi * s),
                                      np.cos(4 * np.pi * s), np.sin(4 * np.pi * s)])


bases = {"quadratic": lambda s: polynomial(s, 2), "degree 9": lambda s: polynomial(s, 9),
         "degree 16": lambda s: polynomial(s, 16), "piecewise constant": step,
         "piecewise linear": broken_line, "cubic spline": cubic_spline}
# <</bases>>

# <<compare>>
for name, f in bases.items():
    X = with_season(f)(t)
    p = np.linalg.matrix_rank(X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = resid @ resid
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)
    print(f"{name:20s} p = {p:2d}   RSS = {rss:8.2f}   sigma = {np.sqrt(rss / (n - p)):.4f}"
          f"   GCV = {n * rss / (n - p) ** 2:.4f}   max leverage {lev.max():.4f}"
          f"   cond = {np.linalg.cond(X):.1e}")
# <</compare>>

# <<beyond>>
train = t < 1991.0
lo, hi = t[train].min(), t[train].max()                    # rescale and reknot on the training range
knots = np.quantile(t[train], np.arange(1, K + 1) / (K + 1))

for name in ("quadratic", "degree 9", "piecewise linear", "cubic spline"):
    X = with_season(bases[name])
    beta, *_ = np.linalg.lstsq(X(t[train]), y[train], rcond=None)
    out = X(t[~train]) @ beta - y[~train]
    print(f"{name:20s} prediction error 1991-2001: root mean square {np.sqrt(np.mean(out ** 2)):9.3f}"
          f"   December 2001 {out[-1]:10.2f}")
# <</beyond>>

# ---- assertions -------------------------------------------------------------
lo, hi = t.min(), t.max()
knots = np.quantile(t, np.arange(1, K + 1) / (K + 1))
table = {}
for name, f in bases.items():
    X = with_season(f)(t)
    p = np.linalg.matrix_rank(X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    lev = np.einsum("ij,ij->i", X @ np.linalg.pinv(X.T @ X), X)
    table[name] = (p, resid @ resid, np.sqrt(resid @ resid / (n - p)),
                   n * (resid @ resid) / (n - p) ** 2, lev.max(), np.linalg.cond(X))
# at the same number of parameters the spline fits better than the polynomial,
# with a smaller maximum leverage
assert table["cubic spline"][0] == table["degree 9"][0]
assert table["cubic spline"][1] < table["degree 9"][1]
assert table["cubic spline"][4] < table["degree 9"][4]
# the piecewise constant fit is much the worst: a step function cannot follow a trend
assert table["piecewise constant"][1] > 20 * table["cubic spline"][1]
# in-sample GCV prefers the degree-16 polynomial to every other basis here
assert min(table, key=lambda k: table[k][3]) == "degree 16"

# the truncated power basis is far worse conditioned than a B-spline basis for the
# same space, and the gap grows with the number of knots
from scipy.interpolate import BSpline

cond_gap = {}
for k_count in (4, 8, 16, 32):
    inner = np.quantile(u(t), np.arange(1, k_count + 1) / (k_count + 1))
    T = np.column_stack([u(t) ** j for j in range(4)]
                        + [np.clip(u(t) - k, 0.0, None) ** 3 for k in inner])
    B = BSpline.design_matrix(np.clip(u(t), -1, 1), np.r_[[-1] * 4, inner, [1] * 4], 3).toarray()
    assert np.linalg.matrix_rank(np.column_stack([T, B])) == T.shape[1]    # same space
    cond_gap[k_count] = (np.linalg.cond(T), np.linalg.cond(B))
assert cond_gap[32][0] > 100 * cond_gap[4][0]
assert cond_gap[32][1] < 2 * cond_gap[4][1]

# a knot outside the range of the data produces a column of zeros
dead = np.clip(u(t) - 1.5, 0.0, None) ** 3
assert np.all(dead == 0.0)

# how far one observation reaches: the change in the fitted curve elsewhere when the
# last response moves by 10 ppm
reach = {}
for name in ("degree 12", "cubic spline", "piecewise constant"):
    f = (lambda s: polynomial(s, 12)) if name == "degree 12" else bases[name]
    X = with_season(f)(t)
    cols = f(t).shape[1]
    coef = (np.linalg.pinv(X.T @ X) @ X[-1])[:cols]        # the trend part of the fit
    w = 10 * (f(t) @ coef)
    reach[name] = (w[-1], w[0])
assert abs(reach["piecewise constant"][1]) < 0.01           # essentially local
assert abs(reach["piecewise constant"][1]) < abs(reach["cubic spline"][1]) < abs(reach["degree 12"][1])

train = t < 1991.0
lo, hi = t[train].min(), t[train].max()
knots = np.quantile(t[train], np.arange(1, K + 1) / (K + 1))
ahead = {}
for name in ("quadratic", "degree 9", "piecewise linear", "cubic spline"):
    X = with_season(bases[name])
    beta, *_ = np.linalg.lstsq(X(t[train]), y[train], rcond=None)
    out = X(t[~train]) @ beta - y[~train]
    ahead[name] = (np.sqrt(np.mean(out ** 2)), out[-1])
# the basis that fits best inside the range is not the one that predicts best outside it
assert ahead["piecewise linear"][0] < ahead["quadratic"][0] < ahead["cubic spline"][0] < ahead["degree 9"][0]

gen = Generated("ch42", "spline_comparison", prefix="cmp")
gen.int("n", n)
gen.int("K", K)
short = {"quadratic": "quad", "degree 9": "d9", "degree 16": "d16",
         "piecewise constant": "pc", "piecewise linear": "pl", "cubic spline": "cs"}
for name, key in short.items():
    p, rss, sigma, gcv, hmax, cond = table[name]
    gen.int(f"p{key}", p)
    gen.num(f"rss{key}", rss, 2)
    gen.num(f"sigma{key}", sigma, 4)
    gen.num(f"gcv{key}", gcv, 4)
    gen.num(f"h{key}", hmax, 4)
    gen.num(f"k{key}", cond, 2, sci=True)
for k_count, (ct, cb) in cond_gap.items():
    gen.num(f"ct{k_count}", ct, 2, sci=True)
    gen.num(f"cb{k_count}", cb, 2)
for name in ("quadratic", "degree 9", "piecewise linear", "cubic spline"):
    gen.num(f"rmse{short[name]}", ahead[name][0], 2)
    gen.num(f"dec{short[name]}", ahead[name][1], 1)
gen.num("reachpoly", reach["degree 12"][1], 3)
gen.num("reachspline", reach["cubic spline"][1], 3)
gen.num("reachstep", reach["piecewise constant"][1], 3)
gen.write()

# ---- figure -----------------------------------------------------------------
lo, hi = t.min(), t.max()
knots = np.quantile(t, np.arange(1, K + 1) / (K + 1))
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
grid = np.linspace(lo, hi, 800)
quad = with_season(lambda s: polynomial(s, 2))
beta_quad, *_ = np.linalg.lstsq(quad(t), y, rcond=None)
base = polynomial(grid, 2) @ beta_quad[:3]
ax = axes[0]
ax.scatter(t, y - quad(t) @ beta_quad, s=2, color=COLORS["muted"], alpha=0.35, linewidths=0)
for name, colour in (("degree 9", COLORS["second"]), ("cubic spline", COLORS["accent"])):
    X = with_season(bases[name])
    beta, *_ = np.linalg.lstsq(X(t), y, rcond=None)
    trend = bases[name](grid) @ beta[:bases[name](grid).shape[1]]
    ax.plot(grid, trend - base, color=colour, label=name)
for k in knots:
    ax.axvline(k, color=COLORS["grid"], linewidth=0.5, zorder=0)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("year")
ax.set_ylabel("ppm above the quadratic")
ax.set_title("(a) ten trend parameters, two ways")
ax.legend(frameon=False, loc="lower center", fontsize=6.5)
ax = axes[1]
window = (t >= 1986.0)
gw = np.linspace(1986.0, hi, 500)
gw = np.sort(np.concatenate([gw, knots[knots > 1986.0] - 1e-7, knots[knots > 1986.0] + 1e-7]))
ax.scatter(t[window], (y - quad(t) @ beta_quad + polynomial(t, 2) @ beta_quad[:3])[window],
           s=3, color=COLORS["muted"], alpha=0.4, linewidths=0)
for name, colour in (("piecewise constant", COLORS["third"]), ("piecewise linear", COLORS["second"]),
                     ("cubic spline", COLORS["accent"])):
    X = with_season(bases[name])
    beta, *_ = np.linalg.lstsq(X(t), y, rcond=None)
    cols = bases[name](gw).shape[1]
    curve = bases[name](gw) @ beta[:cols]
    if name == "piecewise constant":                      # break the line at each jump
        curve = np.where(np.isin(gw, knots + 1e-7), np.nan, curve)
    ax.plot(gw, curve, color=colour, label=name)
for k in knots[knots > 1986.0]:
    ax.axvline(k, color=COLORS["grid"], linewidth=0.5, zorder=0)
ax.set_xlabel("year")
ax.set_ylabel("seasonally adjusted (ppm)")
ax.set_title("(b) the same knots, three degrees")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch42", "spline_comparison"))
