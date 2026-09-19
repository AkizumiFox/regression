"""Chapter 23, Section 6: when resampling fails.

Fixed-seed simulations of six failures:
(a) the bootstrap of a sample maximum puts an atom of about 1 - 1/e on the observed maximum;
(b) with infinite-variance errors, neither the t interval nor the bootstrap interval for a mean
    reaches its nominal coverage, and more data do not help;
(c) a single high-leverage case makes the case-bootstrap distribution of a slope bimodal;
(d) with p/n = 1/2, the residual bootstrap underestimates and the case bootstrap overestimates the
    variance of a coefficient;
(e) with serially correlated errors, resampling residuals independently underestimates the standard
    error of a slope, and resampling blocks of residuals recovers much of it;
(f) with n = 8, the residual-bootstrap percentile interval undercovers.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch23", "failures", prefix="fail")

# ---- (a) the maximum ------------------------------------------------------------------------
# <<maximum>>
rng = np.random.default_rng(2310)
n = 50
u = rng.uniform(size=n)
B = 20_000
boot_max = u[rng.integers(0, n, size=(B, n))].max(axis=1)
atom = np.mean(boot_max == u.max())
print(f"bootstrap P*(max* = max) = {atom:.3f}; 1 - (1 - 1/n)^n = {1 - (1 - 1 / n) ** n:.3f}")
# <</maximum>>
atom_exact = 1 - (1 - 1 / n) ** n
assert abs(atom - atom_exact) < 4 * np.sqrt(atom_exact * (1 - atom_exact) / B)
gen.int("max:n", n)
gen.num("max:atom", atom, 3)
gen.num("max:exact", atom_exact, 3)
gen.num("max:limit", 1 - np.exp(-1), 3)

# ---- (b) heavy tails ------------------------------------------------------------------------
# <<heavy>>
def heavy_coverage(n, reps, B, rng, alpha_tail=1.5):
    """Coverage of 95% t and bootstrap-percentile intervals for the mean of a Pareto law."""
    mean = alpha_tail / (alpha_tail - 1)          # Pareto(alpha) on [1, inf): infinite variance
    hit_t = hit_b = 0
    for _ in range(reps):
        y = (1 - rng.uniform(size=n)) ** (-1 / alpha_tail)
        m, s = y.mean(), y.std(ddof=1)
        hit_t += abs(m - mean) <= stats.t.ppf(0.975, n - 1) * s / np.sqrt(n)
        boot = y[rng.integers(0, n, size=(B, n))].mean(axis=1)
        lo, hi = np.quantile(boot, [0.025, 0.975])
        hit_b += lo <= mean <= hi
    return hit_t / reps, hit_b / reps


rng = np.random.default_rng(2311)
for n_h in (25, 400):
    print(n_h, heavy_coverage(n_h, reps=200, B=499, rng=rng))
# <</heavy>>
rng = np.random.default_rng(2312)
heavy_n = [25, 100, 400, 1600]
heavy = {m: heavy_coverage(m, reps=2000, B=999, rng=rng) for m in heavy_n}
for m, (ct, cb) in heavy.items():
    print(f"n = {m}: t {ct:.3f}, bootstrap {cb:.3f}")
for m in heavy_n:
    assert heavy[m][0] < 0.85 and heavy[m][1] < 0.85  # far below 0.95 at every n
for k in (0, 1):
    assert abs(heavy[1600][k] - heavy[400][k]) < 0.03  # and the coverage has stopped improving
for m, (ct, cb) in heavy.items():
    gen.num(f"heavy:{m}:t", ct, 3)
    gen.num(f"heavy:{m}:boot", cb, 3)

# ---- (c) one high-leverage case -------------------------------------------------------------
# <<leverage>>
rng = np.random.default_rng(2313)
n_l = 15
x_l = np.append(np.linspace(0, 1, n_l - 1), 4.0)    # one case far from the rest
X_l = np.column_stack([np.ones(n_l), x_l])
h_l = np.sum((X_l @ np.linalg.inv(X_l.T @ X_l)) * X_l, axis=1)
y_l = 1 + 0.5 * x_l + 0.3 * rng.normal(size=n_l)
y_l[-1] -= 1.0                                     # and not quite on the line of the others
idx = rng.integers(0, n_l, size=(20_000, n_l))
xs, ys = x_l[idx], y_l[idx]
xc = xs - xs.mean(axis=1, keepdims=True)
b_case = np.sum(xc * ys, axis=1) / np.sum(xc ** 2, axis=1)
absent = ~np.any(idx == n_l - 1, axis=1)
print(f"leverage of the far case {h_l[-1]:.3f}; resamples without it: {absent.mean():.3f}")
# <</leverage>>
absent_exact = (1 - 1 / n_l) ** n_l
assert abs(absent.mean() - absent_exact) < 0.01
b_full = np.polyfit(x_l, y_l, 1)[0]
b_without = np.polyfit(x_l[:-1], y_l[:-1], 1)[0]
# the two modes: resamples with the far case sit near the full-data slope, those without it spread
# around the slope of the other fourteen
assert abs(np.median(b_case[~absent]) - b_full) < 0.05
assert abs(np.median(b_case[absent]) - b_without) < 0.1
assert b_case[absent].std() > 3 * b_case[~absent].std()
gen.int("lev:n", n_l)
gen.num("lev:h", h_l[-1], 3)
gen.num("lev:absent", absent.mean(), 3)
gen.num("lev:absent_exact", absent_exact, 3)
gen.num("lev:b_full", b_full, 3)
gen.num("lev:b_without", b_without, 3)
gen.num("lev:sd_with", b_case[~absent].std(), 3)
gen.num("lev:sd_without", b_case[absent].std(), 3)

# ---- (d) many parameters --------------------------------------------------------------------
# <<dimension>>
def variance_ratios(n, p, reps, B, rng):
    """Bootstrap variance of the first coefficient divided by its true variance (errors N(0, 1))."""
    X = rng.normal(size=(n, p))                   # fixed design
    XtX_inv = np.linalg.inv(X.T @ X)
    true_var = XtX_inv[0, 0]
    H = X @ XtX_inv @ X.T
    out = {"residual": [], "rescaled": [], "case": []}
    for _ in range(reps):
        y = rng.normal(size=n)                    # beta = 0
        e = y - H @ y
        fit = H @ y
        for name, r in [("residual", e), ("rescaled", e * np.sqrt(n / (n - p)))]:
            Ys = fit + r[rng.integers(0, n, size=(B, n))]
            out[name].append(np.var(Ys @ (X @ XtX_inv)[:, 0]) / true_var)
        bs = []
        for _ in range(B):
            i = rng.integers(0, n, size=n)
            bs.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0][0])
        out["case"].append(np.var(bs) / true_var)
    return {k: float(np.median(v)) for k, v in out.items()}


rng = np.random.default_rng(2314)
print(variance_ratios(60, 30, reps=20, B=100, rng=rng))
# <</dimension>>
rng = np.random.default_rng(2315)
dim = variance_ratios(60, 30, reps=200, B=400, rng=rng)
print(dim)
assert abs(dim["residual"] - 0.5) < 0.1
assert abs(dim["rescaled"] - 1.0) < 0.15
assert dim["case"] > 1.5
for k, v in dim.items():
    gen.num(f"dim:{k}", v, 2)

# ---- (e) serial correlation -----------------------------------------------------------------
# <<dependence>>
def dependence_study(n, phi, reps, B, block, rng):
    """True sd of the slope, and the average iid and moving-block residual-bootstrap standard errors."""
    x = np.linspace(0, 1, n)
    X = np.column_stack([np.ones(n), x])
    A = np.linalg.solve(X.T @ X, X.T)
    slopes, se_iid, se_block = [], [], []
    starts_max = n - block + 1
    for _ in range(reps):
        eps = np.empty(n)
        eps[0] = rng.normal() / np.sqrt(1 - phi ** 2)
        for i in range(1, n):
            eps[i] = phi * eps[i - 1] + rng.normal()   # AR(1) errors
        y = 1 + 2 * x + eps
        b = A @ y
        e = y - X @ b
        slopes.append(b[1])
        Ys = X @ b + e[rng.integers(0, n, size=(B, n))]
        se_iid.append(np.std(Ys @ A[1]))
        starts = rng.integers(0, starts_max, size=(B, n // block))
        idx = (starts[:, :, None] + np.arange(block)).reshape(B, -1)
        Ys = X @ b + e[idx]
        se_block.append(np.std(Ys @ A[1]))
    return np.std(slopes), np.mean(se_iid), np.mean(se_block)


rng = np.random.default_rng(2316)
print(dependence_study(100, 0.7, reps=100, B=200, block=10, rng=rng))
# <</dependence>>
rng = np.random.default_rng(2317)
sd_true, se_iid, se_block = dependence_study(100, 0.7, reps=2000, B=400, block=10, rng=rng)
print(f"true sd {sd_true:.3f}, iid bootstrap {se_iid:.3f}, block bootstrap {se_block:.3f}")
assert se_iid < 0.6 * sd_true
assert se_iid < se_block < sd_true
gen.num("dep:sd", sd_true, 3)
gen.num("dep:iid", se_iid, 3)
gen.num("dep:block", se_block, 3)
gen.num("dep:ratio_iid", se_iid / sd_true, 2)
gen.num("dep:ratio_block", se_block / sd_true, 2)

# ---- (f) small n -----------------------------------------------------------------------------
# <<small>>
def small_n_coverage(n, reps, B, rng):
    """Coverage of 95% intervals for a slope with normal errors: t, residual percentile, bootstrap-t."""
    x = np.linspace(0, 1, n)
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    A = XtX_inv @ X.T
    c = np.sqrt(XtX_inv[1, 1])
    hits = {"t": 0, "percentile": 0, "studentized": 0}
    for _ in range(reps):
        y = rng.normal(size=n)                    # true slope 0
        b = A[1] @ y
        e = y - X @ (A @ y)
        s = np.sqrt(e @ e / (n - 2))
        hits["t"] += abs(b) <= stats.t.ppf(0.975, n - 2) * s * c
        Ys = X @ (A @ y) + e[rng.integers(0, n, size=(B, n))]
        bs = Ys @ A[1]
        ss = np.sqrt(np.sum((Ys - (Ys @ A.T) @ X.T) ** 2, axis=1) / (n - 2))
        lo, hi = np.quantile(bs, [0.025, 0.975])
        hits["percentile"] += lo <= 0 <= hi
        tl, th = np.quantile((bs - b) / (ss * c), [0.025, 0.975])
        hits["studentized"] += b - th * s * c <= 0 <= b - tl * s * c
    return {k: v / reps for k, v in hits.items()}


rng = np.random.default_rng(2318)
print(small_n_coverage(8, reps=300, B=499, rng=rng))
# <</small>>
rng = np.random.default_rng(2319)
small = small_n_coverage(8, reps=4000, B=999, rng=rng)
print(small)
# the percentile interval is close to the normal interval with variance (n-p)/n s^2 c^2:
# its coverage is approximately Pr(|t_6| <= 1.96 sqrt(6/8))
approx = 2 * stats.t.cdf(stats.norm.ppf(0.975) * np.sqrt(6 / 8), 6) - 1
mc = np.sqrt(0.05 * 0.95 / 4000)
assert abs(small["t"] - 0.95) < 4 * mc
assert small["percentile"] < 0.9 and abs(small["percentile"] - approx) < 0.03
assert small["studentized"] > small["percentile"] + 0.03
for k, v in small.items():
    gen.num(f"small:{k}", v, 3)
gen.num("small:approx", approx, 3)
gen.num("small:q", stats.norm.ppf(0.975) * np.sqrt(6 / 8), 3)
gen.write()

# ---- figure: the leverage case and heavy tails ------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
bins = np.linspace(np.quantile(b_case, 0.005), np.quantile(b_case, 0.995), 70)
ax.hist([b_case[~absent], b_case[absent]], bins=bins, stacked=True, density=True,
        color=[COLORS["accent"], COLORS["second"]], label=["far case drawn", "far case absent"])
ax.axvline(b_full, color=COLORS["ink"], linewidth=0.7)
ax.set_xlabel("case-bootstrap slope")
ax.set_yticks([])
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax.set_title("(a) one high-leverage case")
ax = axes[1]
ax.plot(heavy_n, [heavy[m][0] for m in heavy_n], "o-", color=COLORS["accent"], markersize=4, label="t interval")
ax.plot(heavy_n, [heavy[m][1] for m in heavy_n], "s-", color=COLORS["second"], markersize=4, label="bootstrap")
ax.axhline(0.95, color=COLORS["ink"], linewidth=0.6, linestyle=":")
ax.set_xscale("log")
ax.set_xticks(heavy_n)
ax.set_xticklabels([str(m) for m in heavy_n])
ax.set_ylim(0.6, 1.0)
ax.set_xlabel("n")
ax.set_ylabel("coverage")
ax.legend(frameon=False, fontsize=7, loc="lower right")
ax.set_title("(b) infinite variance")
fig.tight_layout()
fig.savefig(figure_path("ch23", "failures"))
