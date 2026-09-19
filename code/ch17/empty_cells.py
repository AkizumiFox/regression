"""Chapter 17, Section 4: empty cells. The variety trial with two combinations never grown,
the interaction test on its estimable part, what Type III becomes, how statsmodels behaves,
and the dimension of the estimable interaction contrasts on random designs.
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

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

# <<empty>>
def proj(*blocks):
    """Orthogonal projection onto the span of the given column blocks (any rank)."""
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T

def rank(*blocks):
    return np.linalg.matrix_rank(np.column_stack(blocks))

missing = ((trial["variety"] == "V1") & (trial["site"] == "S3")) | \
          ((trial["variety"] == "V3") & (trial["site"] == "S1"))
sub = trial[~missing].reset_index(drop=True)          # V1 never grown at S3, V3 never at S1
y = sub["y"].to_numpy()
n = len(y)
one = np.ones((n, 1))
ZA = pd.get_dummies(sub["variety"]).to_numpy(float)
ZB = pd.get_dummies(sub["site"]).to_numpy(float)
W = np.column_stack([ZA[:, [i]] * ZB for i in range(3)])
m = rank(W)                                           # occupied cells
M, PAB, PA, PB = proj(W), proj(one, ZA, ZB), proj(one, ZA), proj(one, ZB)
sse = y @ (np.eye(n) - M) @ y
df_inter = m - rank(one, ZA, ZB)                      # m - a - b + 1 for a connected design
ss_inter = y @ (M - PAB) @ y
F_inter = ss_inter / df_inter / (sse / (n - m))
print(f"occupied cells {m}, interaction df {df_inter}, SS {ss_inter:.3f}, F {F_inter:.2f}")
ss_typeII_A = y @ (PAB - PB) @ y                      # well defined, a - 1 = 2 df
# <</empty>>

# <<type3empty>>
def sum_coding(labels):
    D = pd.get_dummies(labels).to_numpy(float)
    return D[:, :-1] - D[:, [-1]]

CA, CB = sum_coding(sub["variety"]), sum_coding(sub["site"])
CAB = np.column_stack([CA[:, [i]] * CB for i in range(2)])
df_III_A = m - rank(one, CB, CAB)
ss_III_A = y @ (M - proj(one, CB, CAB)) @ y
print(f"last-in SS for variety with sum coding: {abs(ss_III_A):.3f} on {df_III_A} df")
# <</type3empty>>

p_inter = stats.f.sf(F_inter, df_inter, n - m)

# ---- checks -------------------------------------------------------------------
assert m == 7 and df_inter == 2 and n - m == 26
assert abs(ss_III_A) < 1e-9 and df_III_A == 0         # the complete rows: only V2, so nothing to test
# the two tetrads that avoid the empty cells span the estimable interaction contrasts
cells = sub.groupby(["variety", "site"])["y"]
yb = cells.mean()
nn = cells.size()
def cell_vec(coef):
    return np.array([coef.get(k, 0.0) for k in yb.index])
T1 = cell_vec({("V1", "S1"): 1, ("V1", "S2"): -1, ("V2", "S1"): -1, ("V2", "S2"): 1})
T2 = cell_vec({("V2", "S2"): 1, ("V2", "S3"): -1, ("V3", "S2"): -1, ("V3", "S3"): 1})
L = np.array([T1, T2])
u = L @ yb.to_numpy()
assert np.isclose(u @ np.linalg.solve(L @ np.diag(1 / nn.to_numpy()) @ L.T, u), ss_inter)
# estimability in the cell-means model: V1's unweighted mean is not estimable,
# V1 - V2 over the shared sites S1, S2 is
X_cells = W[:, np.abs(W).sum(axis=0) > 0]
full_cells = [(f"V{i + 1}", f"S{j + 1}") for i in range(3) for j in range(3)]
def estimable(coef):
    lam = np.array([coef.get(c, 0.0) for c in full_cells])
    return np.linalg.matrix_rank(np.vstack([W, lam])) == np.linalg.matrix_rank(W)
assert not estimable({("V1", f"S{j}"): 1 / 3 for j in (1, 2, 3)})
assert estimable({("V1", "S1"): 0.5, ("V1", "S2"): 0.5, ("V2", "S1"): -0.5, ("V2", "S2"): -0.5})
shared = 0.5 * (yb[("V1", "S1")] + yb[("V1", "S2")] - yb[("V2", "S1")] - yb[("V2", "S2")])
se_shared = 0.5 * np.sqrt(sse / (n - m) * sum(1 / nn[c] for c in
                          [("V1", "S1"), ("V1", "S2"), ("V2", "S1"), ("V2", "S2")]))
# the additive model is connected and estimates every cell, including the empty ones
add = smf.ols("y ~ C(variety) + C(site)", sub).fit()
grid = pd.DataFrame(full_cells, columns=["variety", "site"])
fill = add.predict(grid).to_numpy().reshape(3, 3)
a2add = sm.stats.anova_lm(add, typ=2)
assert np.isclose(a2add.loc["C(variety)", "sum_sq"], ss_typeII_A)

# one empty cell only: Type III compares the complete rows V2 and V3
sub1 = trial[~((trial["variety"] == "V1") & (trial["site"] == "S3"))].reset_index(drop=True)
y1 = sub1["y"].to_numpy()
o1 = np.ones((len(y1), 1))
A1, B1 = sum_coding(sub1["variety"]), sum_coding(sub1["site"])
AB1 = np.column_stack([A1[:, [i]] * B1 for i in range(2)])
W1 = pd.get_dummies(sub1["variety"] + sub1["site"]).to_numpy(float)
ss1 = y1 @ (proj(W1) - proj(o1, B1, AB1)) @ y1
df1 = np.linalg.matrix_rank(W1) - rank(o1, B1, AB1)
c1 = sub1.groupby(["variety", "site"])["y"]
yb1, nn1 = c1.mean(), c1.size()
l23 = np.array([(1 if v == "V2" else -1 if v == "V3" else 0) / 3 for v, s in yb1.index])
ss23 = (l23 @ yb1.to_numpy()) ** 2 / np.sum(l23 ** 2 / nn1.to_numpy())
assert df1 == 1 and np.isclose(ss1, ss23)
# ... and the same value for every labelling of the levels (sum coding puts -1 on the last level)
for order_v in (["V3", "V1", "V2"], ["V2", "V3", "V1"]):
    for order_s in (["S3", "S1", "S2"], ["S1", "S3", "S2"]):
        Dv = np.column_stack([(sub1["variety"] == v).to_numpy(float) for v in order_v])
        Ds = np.column_stack([(sub1["site"] == s).to_numpy(float) for s in order_s])
        Av, Bs = Dv[:, :-1] - Dv[:, [-1]], Ds[:, :-1] - Ds[:, [-1]]
        ABv = np.column_stack([Av[:, [i]] * Bs for i in range(2)])
        assert np.isclose(y1 @ (proj(W1) - proj(o1, Bs, ABv)) @ y1, ss1)

# ---- what statsmodels 0.15 prints for the two-empty-cell data --------------------
def relabel(d, mapping):
    d = d.copy()
    for col in ["variety", "site"]:
        d[col] = d[col].replace(mapping)
    return d

printed = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for mapping in [{}, {"S1": "S3", "S3": "S1"}, {"V1": "V3", "V3": "V1"}]:
        d = relabel(sub, mapping)
        t1 = sm.stats.anova_lm(smf.ols("y ~ C(variety) * C(site)", d).fit(), typ=1)
        t2 = sm.stats.anova_lm(smf.ols("y ~ C(variety) * C(site)", d).fit(), typ=2)
        t3 = sm.stats.anova_lm(smf.ols("y ~ C(variety, Sum) * C(site, Sum)", d).fit(), typ=3)
        Wd = pd.get_dummies(d["variety"] + d["site"]).to_numpy(float)
        od = np.ones((len(d), 1))
        Pd = proj(od, pd.get_dummies(d["variety"]).to_numpy(float),
                  pd.get_dummies(d["site"]).to_numpy(float))
        yd = d["y"].to_numpy()
        assert np.isclose(yd @ (proj(Wd) - Pd) @ yd, ss_inter)   # the right answer is label-free
        printed.append({"t1_inter_ss": t1.iloc[2]["sum_sq"], "t1_inter_df": t1.iloc[2]["df"],
                        "t2_inter_ss": t2.iloc[2]["sum_sq"], "t3_A_ss": t3.iloc[1]["sum_sq"],
                        "t3_A_df": t3.iloc[1]["df"]})
print(pd.DataFrame(printed).round(3))
for row in printed:
    assert row["t1_inter_df"] == 4                     # 4 df reported where there are 2
    assert row["t3_A_df"] == 2 and row["t3_A_ss"] > 1  # a 'Type III' test where there is none
t1_values = [r["t1_inter_ss"] for r in printed]
assert all(abs(v - ss_inter) > 0.5 for v in t1_values)
t2_values = [r["t2_inter_ss"] for r in printed]
t3_values = [r["t3_A_ss"] for r in printed]
assert max(t2_values) - min(t2_values) > 10           # depends on how the levels are named
assert max(t3_values) - min(t3_values) > 10

# ---- dimension of the estimable interaction contrasts on random designs ----------
def interaction_rank(occ):
    """dim of tables d supported on the occupied cells with zero row and column sums."""
    a, b = occ.shape
    cells_ = np.argwhere(occ)
    T = np.zeros((a + b, len(cells_)))
    for k, (i, j) in enumerate(cells_):
        T[i, k] = 1
        T[a + j, k] = 1
    return len(cells_) - np.linalg.matrix_rank(T)

rng_d = np.random.default_rng(17)
for _ in range(300):
    a, b = rng_d.integers(2, 7, size=2)
    occ = rng_d.random((a, b)) < rng_d.uniform(0.25, 0.9)
    occ[rng_d.integers(a), :] |= ~occ.any(axis=0)       # no empty column
    for i in np.flatnonzero(~occ.any(axis=1)):
        occ[i, rng_d.integers(b)] = True                # no empty row
    adj = np.zeros((a + b, a + b))
    for i, j in np.argwhere(occ):
        adj[i, a + j] = adj[a + j, i] = 1
    c = connected_components(csr_matrix(adj), directed=False)[0]
    assert interaction_rank(occ) == occ.sum() - a - b + c

gen = Generated("ch17", "empty_cells", prefix="em")
gen.int("n", n)
gen.int("m", m)
gen.int("dfe", n - m)
gen.int("dfint", df_inter)
gen.num("ssint", ss_inter, 3)
gen.num("Fint", F_inter, 2)
gen.num("pint", p_inter, 3)
gen.num("mse", sse / (n - m), 4)
gen.num("ss2A", ss_typeII_A, 3)
gen.num("shared", shared, 3)
gen.num("seshared", se_shared, 3)
gen.num("fill13", fill[0, 2], 2)
gen.num("fill31", fill[2, 0], 2)
gen.num("ss1", ss1, 3)
for k, r in enumerate(printed):
    gen.num(f"t1int:{k}", r["t1_inter_ss"], 2)
    gen.num(f"t2int:{k}", r["t2_inter_ss"], 2)
    gen.num(f"t3A:{k}", r["t3_A_ss"], 2)
gen.write()

# ---- figure: the design graph and its two independent cycles --------------------
use_book_style()
fig, ax = plt.subplots(figsize=(3.6, 2.5))
occ = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]], bool)
left = {i: (0.0, 2 - i) for i in range(3)}
right = {j: (2.2, 2 - j) for j in range(3)}
cycle1 = {(0, 0), (0, 1), (1, 0), (1, 1)}
cycle2 = {(1, 1), (1, 2), (2, 1), (2, 2)}
for i, j in np.argwhere(occ):
    (x0, y0), (x1, y1_) = left[i], right[j]
    col = COLORS["accent"] if (i, j) in cycle1 - cycle2 else \
        COLORS["second"] if (i, j) in cycle2 - cycle1 else COLORS["thread"]
    ax.plot([x0, x1], [y0, y1_], color=col, linewidth=1.6, zorder=1)
for i, j in np.argwhere(~occ):
    (x0, y0), (x1, y1_) = left[i], right[j]
    ax.plot([x0, x1], [y0, y1_], color=COLORS["muted"], linewidth=0.8, linestyle=(0, (2, 2)),
            zorder=0)
for k, (x, yy) in left.items():
    ax.scatter([x], [yy], s=260, color="white", edgecolors=COLORS["ink"], zorder=2)
    ax.text(x, yy, f"V{k + 1}", ha="center", va="center", fontsize=8, zorder=3)
for k, (x, yy) in right.items():
    ax.scatter([x], [yy], s=260, color="white", edgecolors=COLORS["ink"], zorder=2)
    ax.text(x, yy, f"S{k + 1}", ha="center", va="center", fontsize=8, zorder=3)
ax.set_xlim(-0.5, 2.7)
ax.set_ylim(-0.5, 2.5)
ax.axis("off")
fig.tight_layout()
fig.savefig(figure_path("ch17", "design_graph"))
