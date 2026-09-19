"""Chapter 20, Section 2: how large leverages are in designed and in observational data.

(1) Under a Gaussian random design, h_ii - 1/n is (1 - 1/n) times a Beta(k/2, (n-1-k)/2)
    variable: a check by simulation, and the implied rates for the 2p/n and 3p/n rules.
(2) The largest leverage as n grows, for regressors with light and heavy tails.
(3) A two-level factorial design, where every leverage is p/n.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style


def leverages(X):
    Q, _ = np.linalg.qr(X)
    return np.sum(Q ** 2, axis=1)


# <<beta>>
rng = np.random.default_rng(2020)
n, k = 30, 3                                       # intercept plus k Gaussian regressors
p = k + 1
reps = 4000
h1 = np.empty(reps)
for r in range(reps):
    Z = rng.normal(size=(n, k)) @ np.array([[1.0, 0.0, 0.0], [0.8, 0.6, 0.0], [0.0, 2.0, 1.0]])
    X = np.column_stack([np.ones(n), Z + 5.0])     # correlated, shifted regressors
    h1[r] = leverages(X)[0]                        # leverage of case 1
B = (h1 - 1 / n) / (1 - 1 / n)                     # should be Beta(k/2, (n-1-k)/2)
law = stats.beta(k / 2, (n - 1 - k) / 2)
print("mean of h_11:", h1.mean(), " theory p/n =", p / n)
print("Kolmogorov-Smirnov p-value against the Beta law:", stats.kstest(B, law.cdf).pvalue)
# <</beta>>

assert stats.kstest(B, law.cdf).pvalue > 0.01
assert abs(h1.mean() - p / n) < 0.003


def tail(n, p, c):
    """P(h_ii > c p / n) under a Gaussian design with an intercept and p - 1 regressors."""
    k = p - 1
    return stats.beta(k / 2, (n - 1 - k) / 2).sf((c * p / n - 1 / n) / (1 - 1 / n))


rates = {(nn, pp, c): tail(nn, pp, c) for nn in (50, 500) for pp in (4, 11) for c in (2, 3)}
for key, v in rates.items():
    print(key, round(v, 4))
# large-n limit: n (h_ii - 1/n) is approximately chi-squared with k degrees of freedom
limit = {c: stats.chi2(3).sf(c * 4 - 1) for c in (2, 3)}
assert abs(rates[(500, 4, 2)] - limit[2]) < 0.003 and abs(rates[(500, 4, 3)] - limit[3]) < 0.002

# ---- largest leverage as n grows ---------------------------------------------------
ns = np.array([25, 50, 100, 200, 400, 800, 1600, 3200, 6400])
laws = {
    "normal": lambda size: rng.normal(size=size),
    "t, 3 df": lambda size: rng.standard_t(3, size=size),
    "Cauchy": lambda size: rng.standard_cauchy(size=size),
}
nrep = 300
med = {name: np.empty(len(ns)) for name in laws}
for name, draw in laws.items():
    for a, m in enumerate(ns):
        mx = np.empty(nrep)
        for r in range(nrep):
            X = np.column_stack([np.ones(m), draw((m, 2))])
            mx[r] = leverages(X).max()
        med[name][a] = np.median(mx)
    print(name, np.round(med[name], 3))
assert np.all(np.diff(med["normal"]) < 0) and med["normal"][-1] < 0.01
assert med["t, 3 df"][-1] < med["t, 3 df"][0] / 5
assert med["Cauchy"].min() > 0.5

# ---- a 2^3 factorial with all interactions: every leverage equals p/n ---------------
pts = np.array(list(itertools.product([-1.0, 1.0], repeat=3)))
cols = [np.ones(8)] + [pts[:, j] for j in range(3)] + [pts[:, 0] * pts[:, 1]]
Xf = np.column_stack(cols)                         # main effects and one interaction
Xf = np.vstack([Xf, Xf])                           # two replicates
hf = leverages(Xf)
assert np.allclose(hf, Xf.shape[1] / Xf.shape[0])

gen = Generated("ch20", "leverage_design", prefix="ld")
gen.int("n", n)
gen.int("k", k)
gen.int("reps", reps)
gen.num("ks_p", stats.kstest(B, law.cdf).pvalue, 2)
for key in [(50, 4, 2), (500, 11, 2)]:
    gen.num("rate_%d_%d_%d" % key, rates[key], 3)
gen.num("limit_2", limit[2], 3)
gen.num("limit_3", limit[3], 3)
gen.num("normal_max_small", med["normal"][0], 3)
gen.num("normal_max_big", med["normal"][-1], 4)
gen.num("t3_max_big", med["t, 3 df"][-1], 3)
gen.num("cauchy_max_small", med["Cauchy"][0], 3)
gen.num("cauchy_max_big", med["Cauchy"][-1], 3)
gen.int("n_small", ns[0])
gen.int("n_big", ns[-1])
gen.num("h_fact", hf[0], 4)
gen.write()

# ---- figure -------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.hist(h1, bins=45, density=True, color=COLORS["accent"], alpha=0.35, linewidth=0)
grid = np.linspace(1 / n + 1e-4, 0.6, 400)
ax.plot(grid, law.pdf((grid - 1 / n) / (1 - 1 / n)) / (1 - 1 / n), color=COLORS["accent"])
for c, ls in [(2, "--"), (3, ":")]:
    ax.axvline(c * p / n, color=COLORS["second"], linestyle=ls, linewidth=0.9)
ax.set_xlim(0, 0.6)
ax.set_xlabel(r"$h_{11}$")
ax.set_ylabel("density")
ax.set_title(r"(a) Gaussian design, $n=%d$, $p=%d$" % (n, p))
ax = axes[1]
styles = {"normal": (COLORS["accent"], "o"), "t, 3 df": (COLORS["third"], "s"), "Cauchy": (COLORS["second"], "^")}
for name, vals in med.items():
    col, mk = styles[name]
    ax.plot(ns, vals, color=col, marker=mk, markersize=3, label=name)
ax.plot(ns, 3 / ns, color=COLORS["muted"], linestyle=":", linewidth=0.8, label=r"$p/n$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$n$")
ax.set_ylabel(r"median of $\max_i h_{ii}$")
ax.set_title("(b) largest leverage, two regressors")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch20", "leverage_design"))
