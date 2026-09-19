"""Chapter 29, Section 2: AIC against BIC.

(a) Consistency. Six independent normal regressors, three of them relevant; all 64 subsets
    (with intercept) are candidates. BIC finds the true subset with probability tending to one;
    AIC's probability tends to Pr(chi2_1 <= 2)^3 = 0.598, since it keeps each irrelevant
    regressor whose squared t statistic exceeds about 2.
(b) Efficiency. The true mean has (here 1000) nonzero coefficients 2/j that decay; the
    candidates are the nested models with the first k regressors, k <= n/4. The loss of the
    AIC choice approaches the best candidate's loss; BIC's does not.
(c) The exact expected Kullback-Leibler discrepancy behind AICc, checked by simulation.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "aic_bic", prefix="ab")
rng = np.random.default_rng(2902)

# ---- (a) probability of choosing the true subset -----------------------------------------
# <<consistency>>
import itertools
import numpy as np
rng = np.random.default_rng(2903)
beta = np.array([1.0, 0.5, 0.25, 0.0, 0.0, 0.0])        # three relevant, three irrelevant
true_set = (0, 1, 2)
subsets = [s for k in range(7) for s in itertools.combinations(range(6), k)]


def choose(X, y, penalty):
    """All-subsets choice by n log(SSE/n) + penalty * (number of columns)."""
    n = len(y)
    best, arg = np.inf, None
    for s in subsets:
        Xs = np.column_stack([np.ones(n), X[:, list(s)]])
        b, *_ = np.linalg.lstsq(Xs, y, rcond=None)
        crit = n * np.log(np.sum((y - Xs @ b) ** 2) / n) + penalty * Xs.shape[1]
        if crit < best:
            best, arg = crit, s
    return arg


def hit_rates(n, reps):
    X = rng.normal(size=(n, 6))                           # the design, fixed for this n
    hits = np.zeros(2)
    for _ in range(reps):
        y = 1 + X @ beta + rng.normal(size=n)
        hits += [choose(X, y, 2.0) == true_set, choose(X, y, np.log(n)) == true_set]
    return hits / reps


for n in (50, 200, 800):
    aic, bic = hit_rates(n, 100)
    print(f"n = {n:4d}: AIC finds the true subset {aic:.2f} of the time, BIC {bic:.2f}")
# <</consistency>>

limit = chi2.cdf(2, 1) ** 3
ns = np.array([25, 50, 100, 200, 400, 800, 1600, 3200])
rates = np.array([hit_rates(n, 1000) for n in ns])
assert abs(rates[-1, 0] - limit) < 0.05                   # AIC settles near 0.598
assert rates[-1, 1] > 0.97 and rates[-1, 1] > rates[2, 1]  # BIC -> 1
assert np.all(rates[-3:, 1] > rates[-3:, 0])
gen.num("limit", limit, 3)
gen.num("pchi1", chi2.sf(2, 1), 3)
gen.num("pchi2", chi2.sf(4, 2), 3)
gen.num("pchi5", chi2.sf(10, 5), 3)
for n, (a, b) in zip(ns, rates):
    gen.num(f"aic_{n}", a, 2)
    gen.num(f"bic_{n}", b, 2)

# ---- (b) efficiency: an infinite-dimensional truth --------------------------------------
J = 1000
coef = 2.0 / np.arange(1, J + 1)                           # every coefficient is nonzero


def loss_ratios(n, reps):
    K = n // 4                                             # candidates: the first k regressors, k <= n/4
    Xall = rng.normal(size=(n, J))
    mu = 1 + Xall @ coef
    Q, _ = np.linalg.qr(np.column_stack([np.ones(n), Xall[:, :K]]))   # nested bases
    coef_mu = Q.T @ mu
    bias2 = mu @ mu - np.cumsum(coef_mu**2)               # ||(I - M_k) mu||^2
    cols = np.arange(1, K + 2)
    E = rng.normal(size=(reps, n))
    QE = E @ Q
    out = np.zeros((reps, 2))
    for r in range(reps):
        y = mu + E[r]
        sse = y @ y - np.cumsum((QE[r] + coef_mu) ** 2)
        loss = np.cumsum(QE[r] ** 2) + bias2               # ||M_k y - mu||^2
        aic = n * np.log(sse / n) + 2 * cols
        bic = n * np.log(sse / n) + np.log(n) * cols
        out[r] = loss[np.argmin(aic)] / loss.min(), loss[np.argmin(bic)] / loss.min()
    return out.mean(axis=0)


ns_b = np.array([100, 200, 400, 800, 1600, 3200])
ratios = np.array([loss_ratios(n, 400) for n in ns_b])
assert np.all(ratios[:, 0] < ratios[:, 1])
assert ratios[-1, 0] < 1.1 and ratios[-1, 1] > 1.4
for n, (a, b) in zip(ns_b, ratios):
    gen.num(f"eff_aic_{n}", a, 2)
    gen.num(f"eff_bic_{n}", b, 2)

# ---- (c) AICc: the exact expected discrepancy ------------------------------------------------
# <<aicc>>
n_c, p_c, reps_c = 20, 5, 100_000
Xc = np.column_stack([np.ones(n_c), rng.normal(size=(n_c, p_c - 1))])
mu_c = Xc @ np.linspace(1, 2, p_c)                         # the model is correct
H = Xc @ np.linalg.solve(Xc.T @ Xc, Xc.T)
Yc = mu_c + rng.normal(size=(reps_c, n_c))
fit_c = Yc @ H
s2ml = np.mean((Yc - fit_c) ** 2, axis=1)                 # ML variance estimate
# discrepancy: expected -2 log-likelihood of an independent copy Y*, at the fitted values
disc = n_c * np.log(2 * np.pi * s2ml) + (n_c + np.sum((fit_c - mu_c) ** 2, axis=1)) / s2ml
target = np.mean(n_c * np.log(2 * np.pi * s2ml)) + n_c * (n_c + p_c) / (n_c - p_c - 2)
print(f"E(discrepancy) {disc.mean():.2f},  formula {target:.2f}")
print(f"penalty in AIC {2 * (p_c + 1)},  in AICc {2 * (p_c + 1) * n_c / (n_c - p_c - 2):.2f}")
# <</aicc>>
assert abs(disc.mean() - target) < 0.15
naive = np.mean(n_c * np.log(2 * np.pi * s2ml) + n_c + 2 * (p_c + 1))
assert disc.mean() - naive > 5                             # AIC's penalty is far too small here
gen.int("n_c", n_c)
gen.int("p_c", p_c)
gen.num("disc", disc.mean(), 2)
gen.num("target", target, 2)
gen.num("pen_aicc", 2 * (p_c + 1) * n_c / (n_c - p_c - 2), 2)
gen.num("gap", disc.mean() - naive, 2)
gen.write()

# ---- figure ------------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
ax.plot(ns, rates[:, 0], "o-", color=COLORS["accent"], markersize=3, label="AIC")
ax.plot(ns, rates[:, 1], "s-", color=COLORS["second"], markersize=3, label="BIC")
ax.axhline(limit, color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax.set_xscale("log")
ax.set_ylim(0, 1.02)
ax.set_xlabel("n")
ax.set_ylabel("Pr(true subset chosen)")
ax.set_title("(a) a finite true model")
ax.legend(frameon=False, loc="lower right")
ax = axes[1]
ax.plot(ns_b, ratios[:, 0], "o-", color=COLORS["accent"], markersize=3, label="AIC")
ax.plot(ns_b, ratios[:, 1], "s-", color=COLORS["second"], markersize=3, label="BIC")
ax.axhline(1, color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax.set_xscale("log")
ax.set_xlabel("n")
ax.set_ylabel("loss / best candidate's loss")
ax.set_title("(b) no finite true model")
ax.set_ylim(0.9, 2.0)
ax.legend(frameon=False, loc="upper left", ncol=2)
fig.tight_layout()
fig.savefig(figure_path("ch29", "aic_bic"))
