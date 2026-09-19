"""Chapter 23, Section 3: coverage of bootstrap confidence intervals for a regression slope.

A fixed skewed design with n = 20 or 40, and two error laws with mean zero and variance one on average:
centred exponential errors with constant variance, and normal errors whose standard deviation is
proportional to the regressor. For each simulated data set, 95% intervals for the slope from the
normal-theory t interval, three residual-bootstrap intervals (percentile, basic, studentized),
the residual-bootstrap BCa interval, the case-bootstrap percentile interval, the HC2 t interval and
the wild bootstrap-t interval. The script uses 2000 data sets and 999 resamples; the runnable cell
uses fewer.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
level = 0.95
lo_q, hi_q = (1 - level) / 2, (1 + level) / 2
beta = np.array([1.0, 1.0])
methods = ["t", "percentile", "basic", "studentized", "BCa", "case", "HC2 t", "wild-t"]


def make_design(n):
    """A skewed fixed design and the quantities every data set reuses."""
    x = np.exp(0.5 * stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n))
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    A = XtX_inv @ X.T                             # beta_hat = A y
    h = np.sum((X @ XtX_inv) * X, axis=1)         # leverages
    sd = x / np.sqrt(np.mean(x ** 2))             # heteroscedastic sd, mean square one
    return dict(n=n, x=x, X=X, A=A, h=h, a=A[1], c=np.sqrt(XtX_inv[1, 1]), sd=sd)


def errors(rng, law, d):
    if law == "skewed":                           # centred exponential, constant variance
        return rng.exponential(size=d["n"]) - 1.0
    return d["sd"] * rng.normal(size=d["n"])      # normal, sd proportional to x


def one_data_set(rng, law, d, B):
    """Return, for each method, whether its 95% interval covers the true slope."""
    n, x, X, A, h, a, c = (d[k] for k in ("n", "x", "X", "A", "h", "a", "c"))
    y = X @ beta + errors(rng, law, d)
    b = a @ y
    e = y - X @ (A @ y)
    s = np.sqrt(e @ e / (n - 2))
    cover = {"t": abs(b - beta[1]) <= stats.t.ppf(hi_q, n - 2) * s * c}
    # residual bootstrap from leverage-adjusted, centred residuals
    r = e / np.sqrt(1 - h)
    r -= r.mean()
    Ystar = X @ (A @ y) + r[rng.integers(0, n, size=(B, n))]
    bstar = Ystar @ a
    Estar = Ystar - (Ystar @ A.T) @ X.T
    tstar = (bstar - b) / (np.sqrt(np.sum(Estar ** 2, axis=1) / (n - 2)) * c)
    q_lo, q_hi = np.quantile(bstar, [lo_q, hi_q])
    t_lo, t_hi = np.quantile(tstar, [lo_q, hi_q])
    cover["percentile"] = q_lo <= beta[1] <= q_hi
    cover["basic"] = 2 * b - q_hi <= beta[1] <= 2 * b - q_lo
    cover["studentized"] = b - t_hi * s * c <= beta[1] <= b - t_lo * s * c
    # BCa: residual-bootstrap bias correction, case-jackknife acceleration (a heuristic pairing)
    z0 = stats.norm.ppf(np.clip(np.mean(bstar < b), 1 / B, 1 - 1 / B))
    loo = b - a * e / (1 - h)                     # slopes with one case deleted
    dev = loo.mean() - loo
    acc = np.sum(dev ** 3) / (6 * np.sum(dev ** 2) ** 1.5)
    zs = z0 + stats.norm.ppf([lo_q, hi_q])
    bca_lo, bca_hi = np.quantile(bstar, stats.norm.cdf(z0 + zs / (1 - acc * zs)))
    cover["BCa"] = bca_lo <= beta[1] <= bca_hi
    # case bootstrap percentile
    idx = rng.integers(0, n, size=(B, n))
    xs, ys = x[idx], y[idx]
    xc = xs - xs.mean(axis=1, keepdims=True)
    sxx = np.sum(xc ** 2, axis=1)
    ok = sxx > 1e-12                              # drop the (rare) resamples with one distinct x
    bcase = np.sum(xc[ok] * ys[ok], axis=1) / sxx[ok]
    cl, ch = np.quantile(bcase, [lo_q, hi_q])
    cover["case"] = cl <= beta[1] <= ch
    # HC2 t interval, and the wild bootstrap-t studentized by HC2
    se2 = np.sqrt(np.sum(a ** 2 * e ** 2 / (1 - h)))
    cover["HC2 t"] = abs(b - beta[1]) <= stats.t.ppf(hi_q, n - 2) * se2
    v = rng.choice([-1.0, 1.0], size=(B, n))
    Ystar = X @ (A @ y) + (e / np.sqrt(1 - h)) * v
    Estar = Ystar - (Ystar @ A.T) @ X.T
    wstar = (Ystar @ a - b) / np.sqrt((Estar ** 2 / (1 - h)) @ a ** 2)
    w_lo, w_hi = np.quantile(wstar, [lo_q, hi_q])
    cover["wild-t"] = b - w_hi * se2 <= beta[1] <= b - w_lo * se2
    return cover


def coverage(law, n, reps, B, seed):
    rng = np.random.default_rng(seed)
    d = make_design(n)
    hits = {m: 0 for m in methods}
    for _ in range(reps):
        for m, ok in one_data_set(rng, law, d, B).items():
            hits[m] += ok
    return {m: hits[m] / reps for m in methods}
# <</design>>


# <<small>>
for law in ["skewed", "heteroscedastic"]:
    cov = coverage(law, n=20, reps=200, B=499, seed=2331)
    print(f"{law:16s}" + "  ".join(f"{m} {v:.2f}" for m, v in cov.items()))
# <</small>>

reps, B = 2000, 999
laws = ["skewed", "heteroscedastic"]
sizes = [20, 40]
table = {(law, n): coverage(law, n, reps, B, seed=2330 + 10 * k + n)
         for k, law in enumerate(laws) for n in sizes}
for (law, n), cov in table.items():
    print(f"{law:16s} n = {n}: " + "  ".join(f"{m} {v:.3f}" for m, v in cov.items()))

# ---- checks -------------------------------------------------------------------------------
mc_se = np.sqrt(level * (1 - level) / reps)
for n in sizes:
    sk, he = table[("skewed", n)], table[("heteroscedastic", n)]
    # constant variance: the t interval and the studentized bootstrap are close to nominal
    assert abs(sk["t"] - level) < 4 * mc_se
    assert abs(sk["studentized"] - level) < 4 * mc_se
    # the percentile interval and BCa undercover a little
    assert sk["percentile"] < sk["studentized"] and sk["BCa"] < sk["studentized"]
    # heteroscedasticity: every method built on the residual bootstrap fails, and does not improve
    for m in ["t", "percentile", "basic", "studentized", "BCa"]:
        assert he[m] < 0.80
    # the heteroscedasticity-robust methods are much better, though still below nominal
    for m in ["case", "HC2 t", "wild-t"]:
        assert he[m] > he["t"] + 0.05 and he[m] < level
for m in ["case", "wild-t"]:
    assert table[("heteroscedastic", 40)][m] > table[("heteroscedastic", 20)][m] - 2 * mc_se
hmax = {n: make_design(n)["h"].max() for n in sizes}

gen = Generated("ch23", "coverage", prefix="cov")
gen.int("reps", reps)
gen.int("B", B)
gen.num("mcse", mc_se, 3)
for n in sizes:
    gen.num(f"hmax:{n}", hmax[n], 3)
for (law, n), cov in table.items():
    for m, v in cov.items():
        gen.num(f"{law}:{n}:{m.replace(' ', '')}", v, 3)
gen.write()

# ---- figure -------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.7), sharey=True)
ypos = np.arange(len(methods))[::-1]
for ax, n in zip(axes, sizes):
    for law, col, mk, off in [("skewed", COLORS["accent"], "o", 0.13),
                              ("heteroscedastic", COLORS["second"], "s", -0.13)]:
        cov = table[(law, n)]
        ax.scatter([cov[m] for m in methods], ypos + off, color=col, marker=mk, s=16, label=law, zorder=3)
    ax.axvline(level, color=COLORS["ink"], linewidth=0.7)
    ax.axvspan(level - 2 * mc_se, level + 2 * mc_se, color=COLORS["grid"], alpha=0.7, linewidth=0)
    ax.set_xlim(0.66, 1.0)
    ax.set_title(f"n = {n}")
    ax.set_xlabel("coverage")
axes[0].set_yticks(ypos)
axes[0].set_yticklabels(methods)
axes[0].legend(frameon=False, loc="lower left", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch23", "coverage"))
