"""Chapter 14, Section 7: lack-of-fit tests without exact replicates.

(a) Stack loss with all three regressors (statsmodels.datasets.stackloss, public domain):
    only one pair of identical rows, so the pure-error test has 1 df. Near-replicate clusters,
    a central-subset (rainbow) test, added squares and RESET instead.
(b) A simulation with a fixed design: the size of every test is exact, and their powers
    against a smooth curvature and against a local bump differ.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from statsmodels.stats.diagnostic import linear_reset

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "lof_alternatives", prefix="alt")


def proj(A):
    """Orthogonal projection onto C(A), via an orthonormal basis of the column space."""
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    U = U[:, s > s[0] * 1e-10]
    return U @ U.T


def f_test(y, X, Zbig):
    """F test of C(X) against a larger C(Zbig) ⊇ C(X); returns F, df1, df2, p."""
    n = len(y)
    M, MZ = proj(X), proj(Zbig)
    r, s = round(np.trace(M)), round(np.trace(MZ))
    sse, sse_z = y @ (y - M @ y), y @ (y - MZ @ y)
    F = ((sse - sse_z) / (s - r)) / (sse_z / (n - s))
    return F, s - r, n - s, stats.f.sf(F, s - r, n - s)


# ---- (a) stack loss, three regressors ---------------------------------------------------------
# <<stackloss>>
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
R = data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), R])
fit = sm.OLS(y, X).fit()

# exact replicates: rows of X that occur more than once
_, rows = np.unique(R, axis=0, return_inverse=True)
print("distinct rows:", rows.max() + 1, "of", n)

# near replicates: cluster the standardized regressor rows (Ward linkage), c clusters
Rs = (R - R.mean(axis=0)) / R.std(axis=0)
clusters = fcluster(linkage(Rs, method="ward"), t=8, criterion="maxclust")
D = np.eye(clusters.max())[clusters - 1]              # cluster indicators
tests = {"near replicates": f_test(y, X, np.column_stack([X, D])),
         "added squares": f_test(y, X, np.column_stack([X, R ** 2]))}

# rainbow: refit on the half of the cases with the smallest leverages
h = np.diag(proj(X))
central = np.argsort(h)[: n // 2 + 1]
E = np.eye(n)[:, np.setdiff1d(np.arange(n), central)]  # one free mean per outer case
tests["rainbow"] = f_test(y, X, np.column_stack([X, E]))

# RESET: add squares and cubes of the fitted values
reset = linear_reset(fit, power=3, test_type="fitted", use_f=True)
tests["RESET"] = (reset.fvalue, reset.df_num, reset.df_denom, reset.pvalue)
for name, (F, d1, d2, p) in tests.items():
    print(f"{name:16s} F = {F:6.2f} on ({d1:.0f}, {d2:.0f}) df,  p = {p:.3f}")
# <</stackloss>>
assert rows.max() + 1 == n - 1                        # one pair of identical rows
# RESET by hand: it is the F test of X against [X, fitted^2, fitted^3]
yh = fit.fittedvalues
F_hand = f_test(y, X, np.column_stack([X, yh ** 2, yh ** 3]))
assert np.isclose(F_hand[0], reset.fvalue)
gen.int("n", n)
gen.int("c", clusters.max())
gen.int("m", len(central))
for key, name in [("nr", "near replicates"), ("sq", "added squares"),
                  ("rb", "rainbow"), ("rs", "RESET")]:
    F, d1, d2, p = tests[name]
    gen.num(f"F_{key}", F, 2)
    gen.int(f"d1_{key}", round(d1))
    gen.int(f"d2_{key}", round(d2))
    gen.num(f"p_{key}", p, 3)

# ---- (b) power with a fixed design ---------------------------------------------------------------
# <<power>>
rng = np.random.default_rng(1471)
n2, reps = 50, 4000
x = (np.arange(n2) + 0.5) / n2                         # fixed design on (0, 1)
X2 = np.column_stack([np.ones(n2), x])
M2 = proj(X2)
knots = (0.25, 0.5, 0.75)
spline = np.column_stack([x ** 2, x ** 3] + [np.clip(x - k, 0, None) ** 3 for k in knots])
cells = np.eye(10)[np.arange(n2) // 5]                 # 10 clusters of 5 neighbouring x
bigger = {"quadratic": np.column_stack([X2, x ** 2]),
          "spline": np.column_stack([X2, spline]),
          "near replicates": np.column_stack([X2, cells])}
proj_big = {k: proj(v) for k, v in bigger.items()}

# kernel smoother of the residuals; T = |S e|^2 / |e|^2, null law free of beta and sigma
bw = 0.06
K = np.exp(-0.5 * ((x[:, None] - x[None, :]) / bw) ** 2)
S = K / K.sum(axis=1, keepdims=True)
null_e = (np.eye(n2) - M2) @ rng.normal(size=(n2, 50_000))
T_null = np.sum((S @ null_e) ** 2, axis=0) / np.sum(null_e ** 2, axis=0)
T_crit = np.quantile(T_null, 0.95)

def powers(mean):
    """Rejection rates at level 0.05 of each test when E(Y) = mean."""
    Y = mean[:, None] + rng.normal(size=(n2, reps))
    e = Y - M2 @ Y
    sse = np.sum(e ** 2, axis=0)
    out = {}
    for name, P in proj_big.items():
        s = round(np.trace(P))
        sse_z = np.sum((Y - P @ Y) ** 2, axis=0)
        F = ((sse - sse_z) / (s - 2)) / (sse_z / (n2 - s))
        out[name] = np.mean(F > stats.f.ppf(0.95, s - 2, n2 - s))
    Fr = np.empty(reps)
    for j in range(reps):                              # RESET: regressors depend on y
        yh = M2 @ Y[:, j]
        W = np.column_stack([yh ** 2, yh ** 3])
        Wt = W - M2 @ W
        g, *_ = np.linalg.lstsq(Wt, e[:, j], rcond=None)
        ss_add = e[:, j] @ (Wt @ g)
        Fr[j] = (ss_add / 2) / ((sse[j] - ss_add) / (n2 - 4))
    out["RESET"] = np.mean(Fr > stats.f.ppf(0.95, 2, n2 - 4))
    T = np.sum((S @ e) ** 2, axis=0) / sse
    out["smoother"] = np.mean(T > T_crit)
    return out

curve = (x - 0.5) ** 2                                 # smooth curvature
bump = np.exp(-0.5 * ((x - 0.8) / 0.05) ** 2)          # a local bump
sizes = powers(1 + 2 * x)
print("sizes:", {k: round(float(v), 3) for k, v in sizes.items()})
reported = {}
for label, shape, a in [("curvature", curve, 6.0), ("bump", bump, 2.5)]:
    reported[label] = powers(1 + 2 * x + a * shape)
    print(f"{label:9s}:", {k: round(float(v), 3) for k, v in reported[label].items()})
# <</power>>
se = np.sqrt(0.05 * 0.95 / reps)
for k, v in sizes.items():
    assert abs(v - 0.05) < 4.5 * se, (k, v)
    gen.num(f"size_{k.split()[0]}", v, 3)
amps_curve = np.linspace(0, 10, 11)
amps_bump = np.linspace(0, 4, 11)
power_curve = {a: powers(1 + 2 * x + a * curve) for a in amps_curve}
power_bump = {a: powers(1 + 2 * x + a * bump) for a in amps_bump}
pc, pb = reported["curvature"], reported["bump"]
assert pc["quadratic"] > pc["near replicates"] and pc["RESET"] > pc["smoother"] - 0.05
assert pb["near replicates"] > pb["quadratic"] and pb["smoother"] > pb["RESET"]
for k in pc:
    gen.num(f"pc_{k.split()[0]}", pc[k], 3)
    gen.num(f"pb_{k.split()[0]}", pb[k], 3)
gen.int("n2", n2)
gen.int("reps", reps)
gen.num("bw", bw, 2)
gen.int("nullreps", null_e.shape[1])
gen.num("acurve", 6.0, 1)
gen.num("abump", 2.5, 1)
gen.write()

# ---- figure ----------------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.6, 2.3), gridspec_kw={"width_ratios": [0.8, 1, 1]})
ax = axes[0]
ax.plot(x, 1 + 2 * x + 6.0 * curve, color=COLORS["accent"], label="curvature")
ax.plot(x, 1 + 2 * x + 2.5 * bump, color=COLORS["second"], label="bump")
ax.plot(x, 1 + 2 * x, color=COLORS["muted"], lw=0.8, ls="--")
ax.set_xlabel(r"$x$")
ax.set_ylabel("mean")
ax.set_title("(a) two departures")
ax.legend(frameon=False, fontsize=7)
styles = {"quadratic": ("-o", COLORS["accent"]), "RESET": ("--s", COLORS["ink"]),
          "spline": ("-^", COLORS["third"]), "near replicates": ("-.D", COLORS["second"]),
          "smoother": (":v", COLORS["thread"])}
for ax, table, amps, title in [(axes[1], power_curve, amps_curve, "(b) curvature"),
                               (axes[2], power_bump, amps_bump, "(c) local bump")]:
    for name, (ls, col) in styles.items():
        ax.plot(amps, [table[a][name] for a in amps], ls, color=col, label=name, lw=1.1,
                ms=2.5)
    ax.axhline(0.05, color=COLORS["grid"], lw=0.6)
    ax.set_ylim(0, 1)
    ax.set_xlabel("size of departure")
    ax.set_title(title)
axes[1].set_ylabel("power")
axes[2].legend(frameon=False, fontsize=6.5, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch14", "lof_power"))
