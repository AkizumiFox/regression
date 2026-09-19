"""Chapter 17, Section 2: Types I, II and III sums of squares for the variety trial,
computed from projections, checked against statsmodels, and matched to the cell-mean
hypotheses of the theorem on what each type tests.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from regbook import COLORS, Generated, figure_path, use_book_style

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

# <<projections>>
def proj(*blocks):
    """Orthogonal projection onto the span of the given column blocks (any rank)."""
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T

y = trial["y"].to_numpy()
n = len(y)
one = np.ones((n, 1))
ZA = pd.get_dummies(trial["variety"]).to_numpy(float)     # variety indicators
ZB = pd.get_dummies(trial["site"]).to_numpy(float)        # site indicators
W = np.column_stack([ZA[:, [i]] * ZB for i in range(3)])  # the nine cell indicators

def ss(P):
    return y @ P @ y

P0, PA, PB, PAB, M = proj(one), proj(one, ZA), proj(one, ZB), proj(one, ZA, ZB), proj(W)
typeI_AB = {"A": ss(PA - P0), "B": ss(PAB - PA), "AB": ss(M - PAB)}   # variety first
typeI_BA = {"B": ss(PB - P0), "A": ss(PAB - PB), "AB": ss(M - PAB)}   # site first
typeII = {"A": ss(PAB - PB), "B": ss(PAB - PA), "AB": ss(M - PAB)}
# <</projections>>

# <<type3>>
def coding(labels, kind):
    """Columns coding a factor: sum-to-zero, or reference level (treatment)."""
    D = pd.get_dummies(labels).to_numpy(float)
    g = D.shape[1]
    C = np.vstack([np.eye(g - 1), -np.ones(g - 1)]) if kind == "sum" else np.eye(g)[:, 1:]
    return D @ C

def last_in(CA, CB):
    """Partial (last-in) sums of squares of A and B with interaction columns = products."""
    CAB = np.column_stack([CA[:, [i]] * CB for i in range(CA.shape[1])])
    full = proj(one, CA, CB, CAB)
    assert np.allclose(full, M)                  # every coding spans the cell-means space
    return {"A": ss(full - proj(one, CB, CAB)), "B": ss(full - proj(one, CA, CAB))}

typeIII = last_in(coding(trial["variety"], "sum"), coding(trial["site"], "sum"))
reference = last_in(coding(trial["variety"], "treatment"), coding(trial["site"], "treatment"))
sse = ss(np.eye(n) - M)
for name, t in [("Type I, variety first", typeI_AB), ("Type I, site first", typeI_BA),
                ("Type II", typeII), ("Type III (sum)", typeIII), ("reference coding", reference)]:
    print(f"{name:22s}", {k: round(float(v), 3) for k, v in t.items()})
# <</type3>>

# <<hypotheses>>
cells = trial.groupby(["variety", "site"])["y"]
ybar = cells.mean().to_numpy()                        # nine cell means, row-major
N = cells.size().to_numpy().reshape(3, 3)             # counts n_ij
npi, npj = N.sum(axis=1), N.sum(axis=0)

def hypothesis_ss(L):
    """Sum of squares of L mu = 0 in the cell-means model (redundant rows allowed)."""
    u = L @ ybar
    return u @ np.linalg.pinv(L @ np.diag(1 / N.ravel()) @ L.T) @ u

# Type II for variety: row i has coefficient n_ij (1[i = k] - n_kj / n_.j) on cell (k, j)
L2A = np.array([[N[i, j] * ((i == k) - N[k, j] / npj[j]) for k in range(3) for j in range(3)]
                for i in range(3)])
# Type I, variety first: the count-weighted variety means are equal
wrow = np.array([[(i == k) * N[k, j] / npi[k] for k in range(3) for j in range(3)]
                 for i in range(3)])
L1 = wrow[[0]] - wrow[1:]
print(f"Type II hypothesis: SS = {hypothesis_ss(L2A):.3f};  Type I hypothesis: SS = {hypothesis_ss(L1):.3f}")
# <</hypotheses>>


sumC = np.vstack([np.eye(2), -np.ones(2)])
helC = np.array([[-1.0, -1.0], [1.0, -1.0], [0.0, 2.0]])
helmert = last_in(ZA @ helC, ZB @ helC)
mixed = last_in(ZA @ np.eye(3)[:, 1:], ZB @ sumC)   # A by reference levels, B by contrasts

# ---- checks against statsmodels ------------------------------------------------
a1 = sm.stats.anova_lm(smf.ols("y ~ C(variety) * C(site)", trial).fit(), typ=1)
a1b = sm.stats.anova_lm(smf.ols("y ~ C(site) * C(variety)", trial).fit(), typ=1)
a2 = sm.stats.anova_lm(smf.ols("y ~ C(variety) * C(site)", trial).fit(), typ=2)
a3 = sm.stats.anova_lm(smf.ols("y ~ C(variety, Sum) * C(site, Sum)", trial).fit(), typ=3)
a3h = sm.stats.anova_lm(smf.ols("y ~ C(variety, Helmert) * C(site, Helmert)", trial).fit(), typ=3)
a3t = sm.stats.anova_lm(smf.ols("y ~ C(variety) * C(site)", trial).fit(), typ=3)
assert np.isclose(a1.loc["C(variety)", "sum_sq"], typeI_AB["A"])
assert np.isclose(a1.loc["C(site)", "sum_sq"], typeI_AB["B"])
assert np.isclose(a1.loc["C(variety):C(site)", "sum_sq"], typeI_AB["AB"])
assert np.isclose(a1b.loc["C(site)", "sum_sq"], typeI_BA["B"])
assert np.isclose(a1b.loc["C(variety)", "sum_sq"], typeI_BA["A"])
assert np.isclose(a2.loc["C(variety)", "sum_sq"], typeII["A"])
assert np.isclose(a2.loc["C(site)", "sum_sq"], typeII["B"])
assert np.isclose(a3.loc["C(variety, Sum)", "sum_sq"], typeIII["A"])
assert np.isclose(a3.loc["C(site, Sum)", "sum_sq"], typeIII["B"])
assert np.isclose(a3.loc["C(variety, Sum):C(site, Sum)", "sum_sq"], typeII["AB"])
assert np.isclose(a3h.loc["C(variety, Helmert)", "sum_sq"], helmert["A"])
assert np.isclose(a3t.loc["C(variety)", "sum_sq"], reference["A"])
assert np.isclose(a3t.loc["C(site)", "sum_sq"], reference["B"])
assert np.isclose(a1.loc["Residual", "sum_sq"], sse)
# Helmert coding is a contrast coding, so it reproduces Type III
assert np.isclose(helmert["A"], typeIII["A"]) and np.isclose(helmert["B"], typeIII["B"])
# the hypothesis of A's last-in SS depends on B's coding only
assert np.isclose(mixed["A"], typeIII["A"])
assert not np.isclose(reference["A"], typeIII["A"])

# ---- the cell-mean hypotheses of each type (cells in row-major order) -----------
assert np.allclose(N, counts)
Dinv = np.diag(1 / N.ravel())
hyp_ss = hypothesis_ss

def cellrow(i, coef):
    """Coefficient vector on the 9 cells with coef (length 3) in row i."""
    v = np.zeros(9)
    v[3 * i:3 * i + 3] = coef
    return v

wrow = list(wrow)
# (b) Type I, B after A: sum_i n_ij (mu_ij - weighted row mean_i) = 0 for every column j
L2B = np.array([sum(N[i, j] * (cellrow(i, np.eye(3)[j]) - wrow[i]) for i in range(3))
                for j in range(3)])
# (c) Type II for A, built a second way from the weighted column means
wcol = [sum(cellrow(k, np.eye(3)[j]) * N[k, j] / npj[j] for k in range(3)) for j in range(3)]
L2A_b = np.array([sum(N[i, j] * (cellrow(i, np.eye(3)[j]) - wcol[j]) for j in range(3))
                  for i in range(3)])
assert np.allclose(L2A, L2A_b)
# (d) Type III: unweighted row means equal
urow = [cellrow(i, np.ones(3) / 3) for i in range(3)]
L3 = np.array([urow[0] - urow[1], urow[0] - urow[2]])
# (e) interaction: the four tetrads of adjacent cells
L4 = np.array([cellrow(i, np.eye(3)[j] - np.eye(3)[j + 1]) - cellrow(i + 1, np.eye(3)[j] - np.eye(3)[j + 1])
               for i in range(2) for j in range(2)])
assert np.linalg.matrix_rank(L2A) == 2 and np.allclose(L2A.sum(axis=0), 0)
assert np.linalg.matrix_rank(L2B) == 2 and np.allclose(L2B.sum(axis=0), 0)
assert np.isclose(hyp_ss(L1), typeI_AB["A"])
assert np.isclose(hyp_ss(L2B), typeI_AB["B"])
assert np.isclose(hyp_ss(L2A), typeII["A"])
assert np.isclose(hyp_ss(L3), typeIII["A"])
assert np.isclose(hyp_ss(L4), typeII["AB"])
# reference coding of B: equality of the variety means at site S1 only
Lref = np.array([cellrow(0, [1, 0, 0]) - cellrow(1, [1, 0, 0]), cellrow(0, [1, 0, 0]) - cellrow(2, [1, 0, 0])])
assert np.isclose(hyp_ss(Lref), reference["A"])
# Helmert coding of B: w orthogonal to its columns is the vector of ones
# general coding theorem: a random coding of B with w = its orthogonal complement
Cb = rng.normal(size=(3, 2))
w = np.linalg.svd(Cb.T)[2][-1]                 # w spans C(Cb)^perp
CA = ZA @ sumC
CB = ZB @ Cb
CAB = np.column_stack([CA[:, [i]] * CB for i in range(2)])
ss_rand = ss(M - proj(one, CB, CAB))
Lw = np.array([cellrow(0, w) - cellrow(1, w), cellrow(0, w) - cellrow(2, w)])
assert np.isclose(ss_rand, hyp_ss(Lw))
# the Type I sequences add up to the model sum of squares, the Type II/III columns do not
for t in (typeI_AB, typeI_BA):
    assert np.isclose(sum(t.values()), ss(M - P0))
assert not np.isclose(sum(typeII.values()), ss(M - P0))
mse = sse / (n - 9)
Fs = {k: v / 2 / mse for k, v in [("I", typeI_AB["A"]), ("Ib", typeI_BA["A"]),
                                  ("II", typeII["A"]), ("III", typeIII["A"]), ("ref", reference["A"])]}
from scipy import stats
ps = {k: stats.f.sf(F, 2, n - 9) for k, F in Fs.items()}
print({k: round(p, 4) for k, p in ps.items()})

gen = Generated("ch17", "sums_of_squares", prefix="ty")
for tag, t in [("1ab", typeI_AB), ("1ba", typeI_BA), ("2", typeII), ("3", typeIII),
               ("ref", reference)]:
    for k, v in t.items():
        gen.num(f"{tag}:{k}", v, 3)
gen.num("model", ss(M - P0), 3)
gen.num("sse", sse, 3)
gen.num("sum2", sum(typeII.values()), 3)
gen.num("sum3", typeIII["A"] + typeIII["B"] + typeII["AB"], 3)
for k in Fs:
    gen.num(f"F:{k}", Fs[k], 2)
    gen.num(f"p:{k}", ps[k], 3)
gen.write()

# ---- figure: the coefficient tables of the first equation of each hypothesis ----
def table(v):
    return v.reshape(3, 3)

panels = [("Type I, variety first", table(wrow[0] - sum(npi[k] * wrow[k] for k in range(3)) / n)),
          ("Type II", table(L2A[0] / npi[0])),
          ("Type III", table(urow[0] - sum(urow) / 3))]
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.2, 2.25))
vmax = max(np.abs(p[1]).max() for p in panels)
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("div", [COLORS["second"], "#FFFFFF", COLORS["accent"]])
for ax, (title, T) in zip(axes, panels):
    ax.imshow(T, cmap=cmap, vmin=-vmax, vmax=vmax)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{T[i, j]:+.2f}".replace("+0.00", "0").replace("-0.00", "0"),
                    ha="center", va="center", fontsize=7.5, color=COLORS["ink"])
    ax.set_xticks(range(3), ["S1", "S2", "S3"])
    ax.set_yticks(range(3), ["V1", "V2", "V3"])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title)
fig.tight_layout()
fig.savefig(figure_path("ch17", "hypothesis_weights"))
