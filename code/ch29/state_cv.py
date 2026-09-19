"""Chapter 29, Section 3: cross-validation for the state regressions.

Log violent crime rate on subsets of five regressors (50 states, statecrime data, public
domain). For each size, the best subset by residual sum of squares; ten-fold CV with its
standard error, and the one-standard-error rule. Leave-one-out (PRESS / n) for all 32 subsets.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "state_cv", prefix="scv")

# <<cv>>
import itertools
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y = np.log(data["violent"].to_numpy())
n = len(y)


def design(cols):
    return np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])


def sse(cols):
    X = design(cols)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.sum((y - X @ b) ** 2)


# best subset of each size, by residual sum of squares
best = [min(itertools.combinations(names, k), key=sse) for k in range(len(names) + 1)]

rng = np.random.default_rng(2906)
K = 10
folds = rng.permutation(np.arange(n) % K)                   # a random split into 10 folds of 5
cv, se = [], []
for cols in best:
    X = design(cols)
    fold_err = []
    for f in range(K):
        tr, te = folds != f, folds == f
        b, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        fold_err.append(np.mean((y[te] - X[te] @ b) ** 2))
    cv.append(np.mean(fold_err))
    se.append(np.std(fold_err, ddof=1) / np.sqrt(K))        # the usual (optimistic) standard error
cv, se = np.array(cv), np.array(se)
k_min = int(np.argmin(cv))
k_1se = int(np.min(np.nonzero(cv <= cv[k_min] + se[k_min])[0]))
for k, cols in enumerate(best):
    print(f"{k} regressors {cols}: CV {cv[k]:.4f} (se {se[k]:.4f})")
print("minimum CV at size", k_min, "; one-standard-error rule chooses size", k_1se)
# <</cv>>

assert best[1] == ("single",) and best[2] == ("single", "urban")
assert k_1se <= k_min
gen.int("k_min", k_min)
gen.int("k_1se", k_1se)
for k in range(len(best)):
    gen.num(f"cv{k}", cv[k], 4)
    gen.num(f"se{k}", se[k], 4)

# leave-one-out for all 32 subsets, by the leverage shortcut
press = {}
for k in range(len(names) + 1):
    for cols in itertools.combinations(names, k):
        X = design(cols)
        H = X @ np.linalg.solve(X.T @ X, X.T)
        e = y - H @ y
        press[cols] = np.mean((e / (1 - np.diag(H))) ** 2)
order = sorted(press, key=press.get)
print("smallest leave-one-out errors:", [(c, round(press[c], 4)) for c in order[:3]])
assert order[0] == ("single", "urban")
gen.num("loo_su", press[("single", "urban")], 4)
gen.num("loo_s", press[("single",)], 4)
gen.num("loo_psu", press[("poverty", "single", "urban")], 4)
gen.num("loo_full", press[tuple(names)], 4)
gen.num("loo_null", press[()], 4)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.8))
ks = np.arange(len(best))
ax.errorbar(ks, cv, yerr=se, fmt="o-", color=COLORS["accent"], capsize=3, markersize=4)
ax.axhline(cv[k_min] + se[k_min], color=COLORS["second"], linestyle="--", linewidth=0.8)
ax.text(1.15, cv[k_min] + se[k_min] + 0.003, "minimum + 1 se", fontsize=7, color=COLORS["second"])
ax.set_xlabel("number of regressors (best subset of each size)")
ax.set_ylabel("10-fold CV error")
ax.set_ylim(0.05, 0.25)
fig.tight_layout()
fig.savefig(figure_path("ch29", "cv_state"))
