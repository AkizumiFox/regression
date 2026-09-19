"""Chapter 17, Section 1: the cell-means model for an unbalanced two-way layout.

Synthetic data: a multi-site trial of three wheat varieties (V1-V3) at three sites (S1-S3),
yield in tonnes per hectare. Each site grew mostly its own favourite variety, so the
counts are far from proportional. The true cell means and the seed are fixed below.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
counts = np.array([[7, 4, 2],        # rows: varieties V1-V3; columns: sites S1-S3
                   [3, 5, 4],
                   [2, 3, 7]])
true_means = np.array([[5.8, 7.2, 7.6],
                       [5.0, 6.5, 7.5],
                       [4.4, 5.9, 7.2]])
rng = np.random.default_rng(20172)
rows = [(f"V{i + 1}", f"S{j + 1}", round(true_means[i, j] + rng.normal(0, 0.5), 1))
        for i in range(3) for j in range(3) for _ in range(counts[i, j])]
trial = pd.DataFrame(rows, columns=["variety", "site", "y"])
# <</data>>

# <<cellmeans>>
cells = trial.groupby(["variety", "site"])["y"]
ybar = cells.mean().unstack().to_numpy()             # 3 x 3 table of cell means
n_ij = cells.size().unstack().to_numpy()
n, m = n_ij.sum(), (n_ij > 0).sum()
sse = ((trial["y"] - cells.transform("mean")) ** 2).sum()
s2 = sse / (n - m)                                   # pooled within-cell variance
se_cell = np.sqrt(s2 / n_ij)
weighted_rows = (n_ij * ybar).sum(axis=1) / n_ij.sum(axis=1)   # = raw variety averages
unweighted_rows = ybar.mean(axis=1)
se_weighted = np.sqrt(s2 / n_ij.sum(axis=1))
se_unweighted = np.sqrt(s2 * (1 / n_ij).sum(axis=1)) / 3
print(np.round(ybar, 3))
print("weighted  :", weighted_rows.round(3), " se", se_weighted.round(3))
print("unweighted:", unweighted_rows.round(3), " se", se_unweighted.round(3))
# <</cellmeans>>

# <<hypotheses>>
def hypothesis_test(L, ybar, n_ij, s2, dfe):
    """F test of L mu = 0 in the cell-means model, all cells occupied (cells in row-major order)."""
    u = L @ ybar.ravel()
    ss = u @ np.linalg.solve(L @ np.diag(1 / n_ij.ravel()) @ L.T, u)
    q = np.linalg.matrix_rank(L)
    F = ss / q / s2
    return ss, q, F, stats.f.sf(F, q, dfe)

C3 = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]])   # any 3 x 2 matrix of contrasts will do
avg3 = np.ones((1, 3)) / 3
L_A = np.kron(C3.T, avg3)                # unweighted variety means equal
L_B = np.kron(avg3, C3.T)                # unweighted site means equal
L_AB = np.kron(C3.T, C3.T)               # all interaction contrasts zero
for name, L in [("variety", L_A), ("site", L_B), ("interaction", L_AB)]:
    ss, q, F, p = hypothesis_test(L, ybar, n_ij, s2, n - m)
    print(f"{name:12s} SS = {ss:6.3f}  df = {q}  F = {F:6.2f}  p = {p:.4f}")
# <</hypotheses>>

res = {name: hypothesis_test(L, ybar, n_ij, s2, n - m)
       for name, L in [("A", L_A), ("B", L_B), ("AB", L_AB)]}

# ---- checks -------------------------------------------------------------------
# the cell-means model is a one-way layout on the nine cells: fitted values are cell averages
W = pd.get_dummies(trial["variety"] + trial["site"]).to_numpy(float)
y = trial["y"].to_numpy()
coef, *_ = np.linalg.lstsq(W, y, rcond=None)
assert np.allclose(np.sort(coef), np.sort(ybar.ravel()))
assert np.isclose(np.sum((y - W @ coef) ** 2), sse)
# weighted row means are the raw variety averages
assert np.allclose(weighted_rows, trial.groupby("variety")["y"].mean().to_numpy())
# the unweighted mean is never more precise than the weighted one (harmonic <= arithmetic)
assert np.all(se_unweighted >= se_weighted - 1e-12)
# the hypothesis does not depend on the choice of contrast matrix
H = np.array([[1.0, 1.0], [-1.0, 1.0], [0.0, -2.0]])    # Helmert-type contrasts
assert np.isclose(hypothesis_test(np.kron(H.T, avg3), ybar, n_ij, s2, n - m)[0], res["A"][0])
# the direct SS equals SSE(reduced) - SSE(full) computed by least squares under the constraint
from scipy.linalg import null_space
K = null_space(L_A)                       # cell means allowed by H_A
Wc = pd.get_dummies(pd.Categorical(trial["variety"] + trial["site"],
                    categories=[f"V{i + 1}S{j + 1}" for i in range(3) for j in range(3)])).to_numpy(float)
th, *_ = np.linalg.lstsq(Wc @ K, y, rcond=None)
assert np.isclose(np.sum((y - Wc @ K @ th) ** 2) - sse, res["A"][0])
# contrast: V1 minus V2, unweighted over sites
d12 = unweighted_rows[0] - unweighted_rows[1]
se12 = np.sqrt(s2 * ((1 / n_ij[0]).sum() + (1 / n_ij[1]).sum())) / 3
t975 = stats.t.ppf(0.975, n - m)
raw12 = weighted_rows[0] - weighted_rows[1]

gen = Generated("ch17", "cell_means", prefix="cm")
gen.int("n", n)
gen.int("m", m)
gen.int("dfe", n - m)
for i in range(3):
    for j in range(3):
        gen.int(f"n:{i + 1}{j + 1}", n_ij[i, j])
        gen.num(f"ybar:{i + 1}{j + 1}", ybar[i, j], 2)
        gen.num(f"se:{i + 1}{j + 1}", se_cell[i, j], 3)
    gen.int(f"nrow:{i + 1}", n_ij[i].sum())
    gen.num(f"w:{i + 1}", weighted_rows[i], 3)
    gen.num(f"u:{i + 1}", unweighted_rows[i], 3)
    gen.num(f"sew:{i + 1}", se_weighted[i], 3)
    gen.num(f"seu:{i + 1}", se_unweighted[i], 3)
for j in range(3):
    gen.num(f"wcol:{j + 1}", (n_ij[:, j] * ybar[:, j]).sum() / n_ij[:, j].sum(), 3)
    gen.num(f"ucol:{j + 1}", ybar[:, j].mean(), 3)
gen.num("sse", sse, 3)
gen.num("s2", s2, 4)
gen.num("s", np.sqrt(s2), 3)
for k, (ss, q, F, p) in res.items():
    gen.num(f"ss:{k}", ss, 3)
    gen.int(f"df:{k}", q)
    gen.num(f"F:{k}", F, 2)
    gen.num(f"p:{k}", p, 4)
gen.num("d12", d12, 3)
gen.num("se12", se12, 3)
gen.num("lo12", d12 - t975 * se12, 2)
gen.num("hi12", d12 + t975 * se12, 2)
gen.num("raw12", raw12, 3)
gen.num("t975", t975, 3)
gen.write()

# ---- figure: cell means with counts --------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.9))
colors = [COLORS["accent"], COLORS["second"], COLORS["third"]]
x = np.arange(3)
for i in range(3):
    ax.plot(x + (i - 1) * 0.10, ybar[i], color=colors[i], linewidth=1.0, zorder=1)
    ax.scatter(x + (i - 1) * 0.10, ybar[i], s=18 * n_ij[i], color=colors[i], alpha=0.85,
               edgecolors="white", linewidths=0.6, zorder=2, label=f"V{i + 1}")
    for j in range(3):
        ax.annotate(str(n_ij[i, j]), (x[j] + (i - 1) * 0.10, ybar[i, j]), xytext=(9, -3),
                    textcoords="offset points", fontsize=7, color=colors[i])
ax.set_xticks(x, ["S1", "S2", "S3"])
ax.set_xlim(-0.4, 2.6)
ax.set_xlabel("site")
ax.set_ylabel("mean yield (t/ha)")
ax.legend(frameon=False, loc="upper left", markerscale=0.5)
fig.tight_layout()
fig.savefig(figure_path("ch17", "cell_means"))
