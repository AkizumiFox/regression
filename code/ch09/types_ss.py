"""Chapter 9, Section 3: Type I, II and III sums of squares in an unbalanced two-way layout.

Self-placement on a 1-7 liberal-conservative scale in the 1996 American National
Election Study, classified by party identification (3 levels) and education
(3 levels). Public-domain data shipped with statsmodels (statsmodels.datasets.anes96).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from regbook import Generated

# <<data>>
d = sm.datasets.anes96.load_pandas().data
party = pd.cut(d["PID"], [-1, 1.5, 4.5, 6], labels=["Dem", "Ind", "Rep"]).rename("party")
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"]).rename("educ")
y = d["selfLR"].to_numpy()
n = len(y)
print(pd.crosstab(party, educ))
print(d.groupby([party, educ], observed=True)["selfLR"].mean().unstack().round(3))
# <</data>>


# <<projections>>
def proj(*blocks):
    """Orthogonal projection onto the span of the given column blocks."""
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]                     # orthonormal basis of C(Z), any rank
    return U @ U.T


one = np.ones((n, 1))
A = pd.get_dummies(party).to_numpy(float)           # indicators of party
B = pd.get_dummies(educ).to_numpy(float)            # indicators of education
AB = np.column_stack([A[:, [i]] * B for i in range(3)])   # the nine cells
P0, PA, PB = proj(one), proj(one, A), proj(one, B)
PAB, M = proj(one, A, B), proj(AB)                  # additive model, cell-means model
I = np.eye(n)


def ss(P):
    return y @ P @ y


type1_AB = {"party": ss(PA - P0), "educ": ss(PAB - PA), "inter": ss(M - PAB)}
type1_BA = {"educ": ss(PB - P0), "party": ss(PAB - PB), "inter": ss(M - PAB)}
type2 = {"party": ss(PAB - PB), "educ": ss(PAB - PA), "inter": ss(M - PAB)}
sse = ss(I - M)
# <</projections>>


# <<type3>>
def coded(labels, kind):
    """Columns for a factor: sum-to-zero ("sum") or reference-level ("treatment") coding."""
    D = pd.get_dummies(labels).to_numpy(float)
    return D[:, 1:] - D[:, [0]] if kind == "sum" else D[:, 1:]


type3 = {}
for kind in ["sum", "treatment"]:
    Ca, Cb = coded(party, kind), coded(educ, kind)
    Cab = np.column_stack([Ca[:, [i]] * Cb for i in range(2)])
    full = proj(one, Ca, Cb, Cab)                    # always the cell-means space
    type3[kind] = {"party": ss(full - proj(one, Cb, Cab)),
                   "educ": ss(full - proj(one, Ca, Cab))}
    assert np.allclose(full, M)
for name, table in [("Type I, party first", type1_AB), ("Type I, education first", type1_BA),
                    ("Type II", type2), ("Type III, sum coding", type3["sum"]),
                    ("'Type III', treatment coding", type3["treatment"])]:
    print(f"{name:30s}", {k: round(float(v), 2) for k, v in table.items()})
# <</type3>>

# ---- checks against statsmodels ------------------------------------------------
df = pd.DataFrame({"y": y, "party": party, "educ": educ})
a1 = sm.stats.anova_lm(smf.ols("y ~ C(party) * C(educ)", df).fit(), typ=1)
a1b = sm.stats.anova_lm(smf.ols("y ~ C(educ) * C(party)", df).fit(), typ=1)
a2 = sm.stats.anova_lm(smf.ols("y ~ C(party) * C(educ)", df).fit(), typ=2)
a3 = sm.stats.anova_lm(smf.ols("y ~ C(party, Sum) * C(educ, Sum)", df).fit(), typ=3)
a3t = sm.stats.anova_lm(smf.ols("y ~ C(party) * C(educ)", df).fit(), typ=3)
assert np.isclose(a1.loc["C(party)", "sum_sq"], type1_AB["party"])
assert np.isclose(a1.loc["C(educ)", "sum_sq"], type1_AB["educ"])
assert np.isclose(a1.loc["C(party):C(educ)", "sum_sq"], type1_AB["inter"])
assert np.isclose(a1b.loc["C(educ)", "sum_sq"], type1_BA["educ"])
assert np.isclose(a1b.loc["C(party)", "sum_sq"], type1_BA["party"])
assert np.isclose(a1b.loc["C(educ):C(party)", "sum_sq"], type1_BA["inter"])
assert np.isclose(a2.loc["C(party)", "sum_sq"], type2["party"])
assert np.isclose(a2.loc["C(educ)", "sum_sq"], type2["educ"])
assert np.isclose(a2.loc["C(party):C(educ)", "sum_sq"], type2["inter"])
assert np.isclose(a3.loc["C(party, Sum)", "sum_sq"], type3["sum"]["party"])
assert np.isclose(a3.loc["C(educ, Sum)", "sum_sq"], type3["sum"]["educ"])
assert np.isclose(a3t.loc["C(party)", "sum_sq"], type3["treatment"]["party"])
assert np.isclose(a3t.loc["C(educ)", "sum_sq"], type3["treatment"]["educ"])
assert np.isclose(a1.loc["Residual", "sum_sq"], sse)
# each Type I sequence adds up to the corrected model sum of squares
for t in (type1_AB, type1_BA):
    assert np.isclose(sum(t.values()), ss(M - P0))
assert not np.isclose(type2["party"] + type2["educ"] + type2["inter"], ss(M - P0))

# ---- which hypothesis about the cell means each sum of squares tests ----------
# cell-means parameterization: mu (9 cells, party-major); L' mu = 0 via the hypothesis SS formula
counts = pd.crosstab(party, educ).to_numpy()
cellmeans = d.groupby([party, educ], observed=True)["selfLR"].mean().unstack().to_numpy()


def hyp_ss(Lam):
    """SS of L' mu = 0 in the cell-means model: (L'ybar)' [L' N^{-1} L]^{-1} (L'ybar)."""
    u = Lam.T @ cellmeans.ravel()
    return u @ np.linalg.solve(Lam.T @ np.diag(1 / counts.ravel()) @ Lam, u)


def row_contrasts(weights):
    """Columns expressing 'weighted row means equal': rows 1-2 and 1-3 compared."""
    L = np.zeros((9, 2))
    for k, i in enumerate([1, 2]):
        L[0:3, k] = weights[0]
        L[3 * i:3 * i + 3, k] = -weights[i]
    return L


unweighted = row_contrasts([np.ones(3) / 3] * 3)
weighted = row_contrasts([counts[i] / counts[i].sum() for i in range(3)])
reference = row_contrasts([np.array([1.0, 0, 0])] * 3)
assert np.isclose(hyp_ss(unweighted), type3["sum"]["party"])       # Type III: unweighted means
assert np.isclose(hyp_ss(weighted), type1_AB["party"])             # Type I first: weighted means
assert np.isclose(hyp_ss(reference), type3["treatment"]["party"])  # treatment coding: HS column only

row_unw = cellmeans.mean(axis=1)
row_w = (counts * cellmeans).sum(axis=1) / counts.sum(axis=1)

gen = Generated("ch09", "types_ss", prefix="typ")
gen.int("n", n)
for i, pl in enumerate(["dem", "ind", "rep"]):
    for j, el in enumerate(["hs", "col", "grad"]):
        gen.int(f"n:{pl}:{el}", counts[i, j])
        gen.num(f"m:{pl}:{el}", cellmeans[i, j], 2)
    gen.num(f"unw:{pl}", row_unw[i], 3)
    gen.num(f"w:{pl}", row_w[i], 3)
for key, table in [("1ab", type1_AB), ("1ba", type1_BA), ("2", type2)]:
    for k, v in table.items():
        gen.num(f"{key}:{k}", v, 2)
for kind in ["sum", "treatment"]:
    for k, v in type3[kind].items():
        gen.num(f"3{kind}:{k}", v, 2)
gen.num("sse", sse, 2)
gen.int("dfe", n - 9)
gen.num("mse", sse / (n - 9), 3)
gen.num("model", ss(M - P0), 2)
for i, pl in enumerate(["dem", "ind", "rep"]):
    gen.num(f"hsshare:{pl}", 100 * counts[i, 0] / counts[i].sum(), 0)
gen.write()
