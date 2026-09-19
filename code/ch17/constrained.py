"""Chapter 17, Section 3: the additive model as a constrained cell-means model, and
hierarchies of constrained models in an unbalanced three-way layout.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.linalg import null_space

from regbook import Generated

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

# <<constrained>>
cells = trial.groupby(["variety", "site"])["y"]
ybar = cells.mean().to_numpy()                        # 9 cell means, row-major
n_c = cells.size().to_numpy()
Dinv = np.diag(1 / n_c)
C3 = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]])
G = np.kron(C3.T, C3.T)                               # 4 interaction contrasts: G mu = 0
K = null_space(G)                                     # 9 x 5 basis of additive tables

V_G = K @ np.linalg.inv(K.T @ np.diag(n_c) @ K) @ K.T           # Cov(mu_hat_G) / sigma^2
mu_G = V_G @ np.diag(n_c) @ ybar                                # constrained estimate
mu_G_direct = ybar - Dinv @ G.T @ np.linalg.solve(G @ Dinv @ G.T, G @ ybar)
print("additive fitted cell means:\n", mu_G.reshape(3, 3).round(3))

n, m, g = n_c.sum(), 9, G.shape[0]
sse = ((trial["y"] - cells.transform("mean")) ** 2).sum()
Gy = G @ ybar
ss_constraint = Gy @ np.linalg.solve(G @ Dinv @ G.T, Gy)       # test of additivity
sse_G = sse + ss_constraint
s2_G = sse_G / (n - m + g)                            # residual mean square of the additive model

def constrained_test(L):
    """F test of L mu = 0 within the constrained model G mu = 0."""
    u = L @ mu_G
    V = L @ V_G @ L.T
    ss = u @ np.linalg.pinv(V) @ u
    q = np.linalg.matrix_rank(L @ K)
    F = ss / q / s2_G
    return ss, q, F, stats.f.sf(F, q, n - m + g)

avg3 = np.ones((1, 3)) / 3
print("variety, unweighted:", np.round(constrained_test(np.kron(C3.T, avg3)), 4))
print("site,    unweighted:", np.round(constrained_test(np.kron(avg3, C3.T)), 4))
# <</constrained>>

# <<lsmeans>>
lsmeans = (mu_G.reshape(3, 3)).mean(axis=1)           # estimated unweighted variety means
L_rows = np.kron(np.eye(3), avg3)
se_ls = np.sqrt(np.diag(L_rows @ V_G @ L_rows.T) * s2_G)
d = np.array([1.0, -1.0, 0.0]) @ L_rows
diff12, se_diff12 = d @ mu_G, np.sqrt(d @ V_G @ d * s2_G)
print("least squares means:", lsmeans.round(3), " se", se_ls.round(3))
print(f"V1 - V2: {diff12:.3f} (se {se_diff12:.3f})")
# <</lsmeans>>

ss_A = constrained_test(np.kron(C3.T, avg3))
ss_B = constrained_test(np.kron(avg3, C3.T))

# ---- checks -------------------------------------------------------------------
assert np.allclose(mu_G, mu_G_direct)
assert np.allclose(G @ mu_G, 0)
add = smf.ols("y ~ C(variety) + C(site)", trial).fit()
fitted_cells = add.fittedvalues.groupby([trial["variety"], trial["site"]]).mean().to_numpy()
assert np.allclose(fitted_cells, mu_G)
assert np.isclose(add.ssr, sse_G)
full = smf.ols("y ~ C(variety) * C(site)", trial).fit()
a1 = sm.stats.anova_lm(full, typ=1)
assert np.isclose(a1.loc["C(variety):C(site)", "sum_sq"], ss_constraint)
# identity K (K'DK)^{-1} K' = D^{-1} - D^{-1} G'(G D^{-1} G')^{-1} G D^{-1}
assert np.allclose(V_G, Dinv - Dinv @ G.T @ np.linalg.solve(G @ Dinv @ G.T, G @ Dinv))
# within the additive model any weighting of the row means gives the same hypothesis
a2add = sm.stats.anova_lm(add, typ=2)
a3add = sm.stats.anova_lm(smf.ols("y ~ C(variety, Sum) + C(site, Sum)", trial).fit(), typ=3)
assert np.isclose(ss_A[0], a2add.loc["C(variety)", "sum_sq"])
assert np.isclose(ss_A[0], a3add.loc["C(variety, Sum)", "sum_sq"])
assert np.isclose(ss_B[0], a2add.loc["C(site)", "sum_sq"])
def row_hypothesis(wts):
    """Rows 1-2 and 1-3 compared, row i averaged over sites with weights wts[i]."""
    return np.array([np.r_[wts[0], np.zeros(6)] - np.r_[np.zeros(3), wts[1], np.zeros(3)],
                     np.r_[wts[0], np.zeros(6)] - np.r_[np.zeros(6), wts[2]]])

w_common = np.array([0.5, 0.3, 0.2])                  # the same site weights in every row ...
assert np.isclose(constrained_test(row_hypothesis([w_common] * 3))[0], ss_A[0])  # ... same test
w_own = counts / counts.sum(axis=1, keepdims=True)    # each row weighted by its own counts
ss_own = constrained_test(row_hypothesis(w_own))[0]
assert not np.isclose(ss_own, ss_A[0])                # a different hypothesis, even under additivity
a1add = sm.stats.anova_lm(add, typ=1)                 # ... namely the Type I hypothesis, variety first
assert np.isclose(ss_own, a1add.loc["C(variety)", "sum_sq"])
Lref = np.array([[1, 0, 0, -1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, -1, 0, 0]], float)
assert np.isclose(constrained_test(Lref)[0], ss_A[0])  # so does the S1 column alone
# LS means differences equal the additive model's variety coefficients
assert np.isclose(diff12, -add.params["C(variety)[T.V2]"])
F_add = ss_constraint / g / (sse / (n - m))
p_add = stats.f.sf(F_add, g, n - m)

gen = Generated("ch17", "constrained", prefix="con")
for i in range(3):
    for j in range(3):
        gen.num(f"mu:{i + 1}{j + 1}", mu_G[3 * i + j], 2)
    gen.num(f"ls:{i + 1}", lsmeans[i], 3)
    gen.num(f"sels:{i + 1}", se_ls[i], 3)
gen.num("ssG", ss_constraint, 3)
gen.num("Fadd", F_add, 2)
gen.num("padd", p_add, 3)
gen.num("sseG", sse_G, 3)
gen.num("s2G", s2_G, 4)
gen.int("dfG", n - m + g)
for k, t in [("A", ss_A), ("B", ss_B)]:
    gen.num(f"ss:{k}", t[0], 3)
    gen.num(f"F:{k}", t[2], 2)
    gen.num(f"p:{k}", t[3], 5)
gen.num("ssown", ss_own, 3)
gen.num("d12", diff12, 3)
gen.num("sed12", se_diff12, 3)
gen.write()

# ---- a three-way layout: a backward path through hierarchical models -------------
# <<threeway>>
import statsmodels.formula.api as smf
rng3 = np.random.default_rng(31)
levels = [(i, j, k) for i in range(3) for j in range(2) for k in range(3)]
n3 = rng3.integers(1, 4, size=len(levels))                     # 1 to 3 runs per cell

def true_hardness(i, j, k):
    """Temperature i, hardener j, supplier k: hardener interacts with both, nothing else does."""
    return (40 + [0, 1.0, 2.0][i] + 1.5 * j + [0, 0.5, -0.5][k]
            + j * [0, -1.0, -2.0][i] + j * [0, 1.2, -0.8][k])

cure = pd.DataFrame([(f"t{i + 1}", f"h{j + 1}", f"s{k + 1}", true_hardness(i, j, k) + rng3.normal(0, 0.7))
                     for (i, j, k), r in zip(levels, n3) for _ in range(r)],
                    columns=["temp", "hardener", "supplier", "y"])

def fit(formula):
    f = smf.ols(formula, cure).fit()
    return f.ssr, int(f.df_resid)

sse_full, df_full = fit("y ~ temp * hardener * supplier")
mse_full = sse_full / df_full
path = [("[THS]", "[TH][TS][HS]", "y ~ temp * hardener * supplier", "y ~ (temp + hardener + supplier) ** 2"),
        ("[TH][TS][HS]", "[TH][HS]", "y ~ (temp + hardener + supplier) ** 2", "y ~ temp * hardener + hardener * supplier"),
        ("[TH][HS]", "[TH][S]", "y ~ temp * hardener + hardener * supplier", "y ~ temp * hardener + supplier"),
        ("[TH][HS]", "[T][HS]", "y ~ temp * hardener + hardener * supplier", "y ~ temp + hardener * supplier"),
        ("[TH][HS]", "[HS]", "y ~ temp * hardener + hardener * supplier", "y ~ hardener * supplier")]
steps = []
for big, small, f_big, f_small in path:
    (s1, d1), (s0, d0) = fit(f_big), fit(f_small)
    F = (s0 - s1) / (d0 - d1) / mse_full
    steps.append((big, small, d0 - d1, s0 - s1, F, stats.f.sf(F, d0 - d1, df_full)))
    print(f"{big:13s} -> {small:13s} df = {d0 - d1}  SS = {s0 - s1:6.3f}  F = {F:5.2f}  p = {steps[-1][5]:.4f}")
# <</threeway>>

gen3 = Generated("ch17", "constrained_threeway", prefix="con3")
gen3.int("n", len(cure))
gen3.int("dfe", df_full)
gen3.num("mse", mse_full, 3)
for k, (big, small, df, ss, F, p) in enumerate(steps):
    gen3.int(f"df:{k + 1}", df)
    gen3.num(f"ss:{k + 1}", ss, 3)
    gen3.num(f"F:{k + 1}", F, 2)
    gen3.num(f"p:{k + 1}", p, 4)
# the first two steps keep; the three ways of shrinking [TH][HS] all fail at the 5% level
assert steps[0][5] > 0.05 and steps[1][5] > 0.05
assert all(st[5] < 0.05 for st in steps[2:])
assert steps[0][2] == 4 and steps[1][2] == 4 and steps[2][2] == 2 and steps[3][2] == 2 and steps[4][2] == 4
# the no-three-factor-interaction model is the constrained cell-means model G3 mu = 0,
# with the (3-1)(2-1)(3-1) = 4 three-factor contrasts; its SS comes from the cell means
C2 = np.array([[1.0], [-1.0]])
G3 = np.kron(np.kron(C3.T, C2.T), C3.T)
cells3 = cure.groupby(["temp", "hardener", "supplier"])["y"]
yb3, nn3 = cells3.mean().to_numpy(), cells3.size().to_numpy()
assert len(yb3) == 18 and nn3.min() >= 1
u3 = G3 @ yb3
ss3 = u3 @ np.linalg.solve(G3 @ np.diag(1 / nn3) @ G3.T, u3)
assert np.isclose(ss3, steps[0][3])
gen3.write()
