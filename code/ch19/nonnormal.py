"""Chapter 19, Section 3: least squares inference with non-normal errors.

(a) Coverage of the nominal 95% t interval for a slope, for a spread design and a design with
    one far point at x = 4n, whose leverage falls only slowly (0.95, 0.80, 0.49 at n = 10, 40,
    160), under normal, centred exponential and scaled t(3) errors.
(b) Coverage of the chi-squared interval for sigma^2 against its limit
    2 Phi(z sqrt(2/(2+kappa))) - 1, where kappa is the excess kurtosis.
(c) The linearized F statistic: Var(L) = 2/f1 + 2/f2 + kappa ||d||^2. A balanced one-way F test
    is protected, the two-sample variance-ratio test is not.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch19", "nonnormal", prefix="nn")

# <<laws>>
import numpy as np
from scipy import stats

def errors(law, size, rng):
    """Errors with mean 0 and variance 1."""
    if law == "normal":
        return rng.normal(size=size)
    if law == "exponential":                       # skewed: skewness 2, excess kurtosis 6
        return rng.exponential(size=size) - 1.0
    if law == "t3":                                # heavy tails: variance 3 before scaling
        return rng.standard_t(3, size=size) / np.sqrt(3.0)
    if law == "uniform":                           # short tails: excess kurtosis -1.2
        return rng.uniform(-np.sqrt(3), np.sqrt(3), size=size)

def design(kind, n):
    x = np.arange(1.0, n + 1)
    if kind == "far point":
        x[-1] = 4.0 * n                            # one point far to the right
    return np.column_stack([np.ones(n), x])
# <</laws>>

# <<slope>>
def slope_coverage(kind, law, n, reps, rng):
    """Proportion of nominal 95% t intervals for the slope that cover it."""
    X = design(kind, n)
    G = np.linalg.inv(X.T @ X)
    A = G @ X.T
    E = errors(law, (reps, n), rng)
    B = E @ A.T                                    # beta_hat - beta
    s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n - 2)
    T = B[:, 1] / np.sqrt(s2 * G[1, 1])
    tq = stats.t.ppf(0.975, n - 2)
    return np.mean(np.abs(T) <= tq), np.mean(T > tq), np.mean(T < -tq)
# <</slope>>

# <<small>>
rng = np.random.default_rng(1903)
for kind in ["spread", "far point"]:
    for law in ["normal", "exponential", "t3"]:
        cov, up, lo = slope_coverage(kind, law, 10, 20_000, rng)
        print(f"{kind:9s} {law:11s} n=10: coverage {cov:.3f} (misses above {up:.3f}, below {lo:.3f})")
# <</small>>

rng = np.random.default_rng(1903)
reps = 100_000
ns = [10, 40, 160]
laws = ["normal", "exponential", "t3"]
table = {}
for kind in ["spread", "far point"]:
    X = design(kind, max(ns))
    for n in ns:
        h = np.diag(design(kind, n) @ np.linalg.solve(design(kind, n).T @ design(kind, n), design(kind, n).T))
        table[(kind, "hmax", n)] = h.max()
        for law in laws:
            table[(kind, law, n)] = slope_coverage(kind, law, n, reps, rng)
for n in ns:
    gen.num(f"hmax_spread_{n}", table[("spread", "hmax", n)], 3)
    gen.num(f"hmax_far_{n}", table[("far point", "hmax", n)], 3)
    for kind, tag in [("spread", "sp"), ("far point", "far")]:
        for law in laws:
            c, up, lo = table[(kind, law, n)]
            gen.num(f"cov_{tag}_{law}_{n}", c, 3)
            gen.num(f"up_{tag}_{law}_{n}", up, 3)
            gen.num(f"lo_{tag}_{law}_{n}", lo, 3)
gen.int("reps", reps)
# normal errors: exact, up to simulation error
for kind in ["spread", "far point"]:
    for n in ns:
        assert abs(table[(kind, "normal", n)][0] - 0.95) < 0.004
# spread design: close to nominal by n = 160 for every law
for law in laws:
    assert abs(table[("spread", law, 160)][0] - 0.95) < 0.006
# far point with skewed errors: the tails are unbalanced, less so as the leverage falls
assert table[("far point", "exponential", 10)][2] < 0.1 * table[("far point", "exponential", 10)][1]
assert table[("far point", "exponential", 160)][2] < 0.5 * table[("far point", "exponential", 160)][1]

# ---- (b) the interval for sigma^2 ------------------------------------------------------
# <<sigma>>
def sigma2_coverage(law, n, reps, rng):
    X = design("spread", n)
    E = errors(law, (reps, n), rng)
    B = np.linalg.solve(X.T @ X, X.T @ E.T).T
    sse = np.sum((E - B @ X.T) ** 2, axis=1)     # sigma^2 = 1
    lo, hi = stats.chi2.ppf([0.025, 0.975], n - 2)
    return np.mean((sse >= lo) & (sse <= hi))

def limit(kappa, z=stats.norm.ppf(0.975)):
    """Limiting coverage of the nominal 95% chi-squared interval."""
    return 2 * stats.norm.cdf(z * np.sqrt(2 / (2 + kappa))) - 1

print("limits:", {k: round(limit(k), 3) for k in [0, 6, -1.2]})
# <</sigma>>
kappas = {"normal": 0.0, "exponential": 6.0, "uniform": -1.2}
ns_sig = [10, 40, 160, 640, 2560]
sig_cov = {}
rng = np.random.default_rng(1904)
for law in kappas:
    for n in ns_sig:
        sig_cov[(law, n)] = sigma2_coverage(law, n, 40_000 if n <= 640 else 10_000, rng)
for law, k in kappas.items():
    gen.num(f"sig_limit_{law}", limit(k), 3)
    for n in ns_sig:
        gen.num(f"sig_{law}_{n}", sig_cov[(law, n)], 3)
assert abs(limit(6.0) - 0.673) < 0.001
assert abs(sig_cov[("exponential", 2560)] - limit(6.0)) < 0.03
assert abs(sig_cov[("uniform", 640)] - limit(-1.2)) < 0.01
assert all(abs(sig_cov[("normal", n)] - 0.95) < 0.01 for n in ns_sig)

# ---- (c) quadratic balance --------------------------------------------------------------
def linearized_var(P1, P2, kappa):
    f1, f2 = np.trace(P1), np.trace(P2)
    d = np.diag(P1) / f1 - np.diag(P2) / f2
    return 2 / f1 + 2 / f2 + kappa * d @ d

def oneway(sizes):
    """Projections for the one-way F test of equal means."""
    g = np.repeat(np.arange(len(sizes)), sizes)
    D = (g[:, None] == np.arange(len(sizes))[None, :]).astype(float)
    N = len(g)
    M = D @ np.linalg.solve(D.T @ D, D.T)
    M0 = np.full((N, N), 1 / N)
    return M - M0, np.eye(N) - M

rng = np.random.default_rng(1905)
reps_q = 100_000
res_q = {}
for label, sizes in [("balanced", [8, 8, 8]), ("unbalanced", [2, 4, 18])]:
    P1, P2 = oneway(sizes)
    f1, f2 = round(np.trace(P1)), round(np.trace(P2))
    E = errors("exponential", (reps_q, sum(sizes)), rng)
    q1 = np.einsum("ri,ij,rj->r", E, P1, E) / f1
    q2 = np.einsum("ri,ij,rj->r", E, P2, E) / f2
    L = q1 - q2
    Fstat = q1 / q2
    size = np.mean(Fstat > stats.f.ppf(0.95, f1, f2))
    v_formula = linearized_var(P1, P2, 6.0)
    assert abs(L.var() / v_formula - 1) < 0.03
    res_q[label] = (size, v_formula, 2 / f1 + 2 / f2)
    gen.num(f"size_{label}", size, 3)
    gen.num(f"vL_{label}", v_formula, 3)
    gen.num(f"vL_normal_{label}", 2 / f1 + 2 / f2, 3)
assert np.isclose(res_q["balanced"][1], res_q["balanced"][2])        # balance: kappa drops out
assert res_q["unbalanced"][1] > 1.2 * res_q["unbalanced"][2]

# two independent samples of 20: the variance-ratio test
n1 = n2 = 20
E1 = errors("exponential", (reps_q, n1), rng)
E2 = errors("exponential", (reps_q, n2), rng)
ratio = E1.var(axis=1, ddof=1) / E2.var(axis=1, ddof=1)
lo, hi = stats.f.ppf([0.025, 0.975], n1 - 1, n2 - 1)
size_var = np.mean((ratio < lo) | (ratio > hi))
z = stats.norm.ppf(0.975)
gen.num("size_varratio", size_var, 3)
gen.num("size_varratio_limit", 2 * (1 - stats.norm.cdf(z * np.sqrt(2 / 8))), 3)
assert size_var > 0.2
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
rng = np.random.default_rng(1906)
n_h = 40
X = design("far point", n_h)
G = np.linalg.inv(X.T @ X)
E = errors("exponential", (100_000, n_h), rng)
B = E @ (G @ X.T).T
s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n_h - 2)
T = B[:, 1] / np.sqrt(s2 * G[1, 1])
bins = np.linspace(-6, 6, 97)
ax.hist(T, bins=bins, density=True, color=COLORS["accent"], alpha=0.55, lw=0, label="simulated")
g = np.linspace(-6, 6, 300)
ax.plot(g, stats.t.pdf(g, n_h - 2), color=COLORS["second"], label=f"$t({n_h - 2})$")
ax.set_xlim(-6, 6)
ax.set_xlabel("$t$ statistic for the slope")
ax.set_title("(a) far point, exponential errors")
ax.legend(frameon=False, loc="upper left")
ax = axes[1]
cols = {"normal": COLORS["ink"], "exponential": COLORS["second"], "uniform": COLORS["third"]}
for law in kappas:
    ax.plot(ns_sig, [sig_cov[(law, n)] for n in ns_sig], "o-", ms=3, color=cols[law], label=law)
    ax.axhline(limit(kappas[law]), color=cols[law], lw=0.7, ls="--")
ax.set_xscale("log")
ax.set_ylim(0.6, 1.01)
ax.set_xlabel("$n$")
ax.set_ylabel("coverage")
ax.set_title(r"(b) 95% interval for $\sigma^2$")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch19", "nonnormal"))
