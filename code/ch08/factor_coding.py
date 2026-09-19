"""Chapter 8, Section 5: one factor, five codings, one fit.

1996 American National Election Study (public domain, statsmodels.datasets.anes96):
respondents' self-placement on a seven-point liberal-conservative scale (selfLR), by their
seven-level party identification (PID). Every coding spans the same column space; the
coefficients are different linear functions of the seven cell means.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].to_numpy()        # 1 = very liberal ... 7 = very conservative
level = anes["PID"].to_numpy().astype(int)    # 0 = strong Dem. ... 6 = strong Rep.
g = 7
Z = np.eye(g)[level]                          # n x 7 indicator (cell-means) matrix
n_k = Z.sum(axis=0)
cell_means = (Z.T @ y) / n_k
print("group sizes:", n_k.astype(int))
print("cell means: ", np.round(cell_means, 3))
# <</data>>

# <<codings>>
def coding_matrices(n_k):
    """Each coding is X = Z K for a nonsingular g x g K; the coefficients are K^{-1} mu."""
    g = len(n_k)
    one = np.ones((g, 1))
    reference = np.vstack([np.zeros(g - 1), np.eye(g - 1)])            # level 0 is baseline
    deviation = np.vstack([np.eye(g - 1), -np.ones(g - 1)])            # sum-to-zero
    weighted = np.vstack([np.eye(g - 1), -n_k[:-1] / n_k[-1]])         # sum n_k * effect = 0
    successive = np.array([[(j + 1 - g if i <= j else j + 1) / g for j in range(g - 1)]
                           for i in range(g)])                         # backward differences
    return {"cell means": np.eye(g),
            "reference": np.hstack([one, reference]),
            "sum to zero": np.hstack([one, deviation]),
            "weighted effect": np.hstack([one, weighted]),
            "successive differences": np.hstack([one, successive])}

for name, K in coding_matrices(n_k).items():
    Xc = Z @ K
    gamma = np.linalg.lstsq(Xc, y, rcond=None)[0]
    print(f"{name:23s}", np.round(gamma, 3))
# <</codings>>

codings = coding_matrices(n_k)
M_ref = Z @ np.diag(1 / n_k) @ Z.T
fits = {}
coefs = {}
for name, K in codings.items():
    assert np.linalg.matrix_rank(K) == g
    Xc = Z @ K
    gamma = np.linalg.lstsq(Xc, y, rcond=None)[0]
    coefs[name] = gamma
    fits[name] = Xc @ gamma
    assert np.allclose(gamma, np.linalg.solve(K, cell_means))       # gamma = K^{-1} mu-hat
    assert np.allclose(fits[name], cell_means[level])
ybar, mbar = y.mean(), cell_means.mean()
cm = cell_means
assert np.allclose(coefs["reference"], np.r_[cm[0], cm[1:] - cm[0]])
assert np.allclose(coefs["sum to zero"], np.r_[mbar, cm[:-1] - mbar])
assert np.allclose(coefs["weighted effect"], np.r_[ybar, cm[:-1] - ybar])
assert np.allclose(coefs["successive differences"], np.r_[mbar, np.diff(cm)])
# the overparameterized model [1, Z] has rank 7, not 8
assert np.linalg.matrix_rank(np.column_stack([np.ones(len(y)), Z])) == g

gen = Generated("ch08", "factor_coding", prefix="fc")
gen.int("n", len(y))
for k in range(g):
    gen.int(f"n{k}", n_k[k])
    gen.num(f"m{k}", cm[k], 3)
for key, name in [("ref", "reference"), ("sum", "sum to zero"), ("wt", "weighted effect"),
                  ("sd", "successive differences")]:
    for i, v in enumerate(coefs[name]):
        gen.num(f"{key}:{i}", v, 3)
gen.num("ybar", ybar, 3)
gen.num("mbar", mbar, 3)
gen.num("sse", np.sum((y - fits["reference"]) ** 2), 1)
# precision: reference = pure independents (n = 37) versus strong Democrats (n = 200)
others = [k for k in range(g) if k not in (0, 3)]
ratios = [np.sqrt((1 / n_k[k] + 1 / n_k[3]) / (1 / n_k[k] + 1 / n_k[0])) for k in others]
assert min(ratios) > 1.5
gen.num("se_ratio_min", min(ratios), 2)
gen.num("se_ratio_max", max(ratios), 2)
gen.write()

# ---- figure: what the intercepts and effects measure ---------------------------------
use_book_style()
labels = ["strong\nDem", "weak\nDem", "ind.\nDem", "ind.", "ind.\nRep", "weak\nRep", "strong\nRep"]
fig, ax = plt.subplots(figsize=(4.6, 2.7))
xs = np.arange(g)
ax.scatter(xs, cm, s=n_k / 3, color=COLORS["accent"], zorder=3)
ax.hlines(cm[0], -0.5, 6.35, color=COLORS["second"], lw=0.9, ls="--")
ax.hlines(mbar, -0.5, 6.35, color=COLORS["third"], lw=0.9)
ax.hlines(ybar, -0.5, 6.35, color=COLORS["thread"], lw=0.9, ls=":")
ax.text(6.45, cm[0], "reference level", va="center", fontsize=7, color=COLORS["second"])
ax.text(6.45, mbar + 0.07, "mean of the means", va="center", fontsize=7, color=COLORS["third"])
ax.text(6.45, ybar - 0.09, "overall mean", va="center", fontsize=7, color=COLORS["thread"])
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=7)
ax.set_xlim(-0.5, 8.3)
ax.set_ylabel("mean self-placement (1-7)")
ax.set_xlabel("party identification")
fig.savefig(figure_path("ch08", "coding_means"))
