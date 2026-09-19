"""Chapter 11, Section 6: t tests for single coefficients, the overall F test, tests of groups.

(1) 2009 state data: murder rate on poverty, high-school graduation, single-parent households,
    percentage white and urbanization (public domain, statsmodels.datasets.statecrime).
(2) Longley's employment series (public domain, statsmodels.datasets.longley) regressed on the
    five economic series without the time trend: two coefficients that are individually
    insignificant but jointly significant. The figure shows the two acceptance regions in the
    plane of the two t statistics.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<coefficients>>
data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
names = ["poverty", "hs_grad", "single", "white", "urban"]
y = data["murder"].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), data[names]])
p = X.shape[1]

C = np.linalg.inv(X.T @ X)
b = C @ X.T @ y
sse = np.sum((y - X @ b) ** 2)
s2 = sse / (n - p)
se = np.sqrt(s2 * np.diag(C))
t = b / se
pv = 2 * stats.t.sf(np.abs(t), n - p)
for name, bj, sj, tj, pj in zip(["intercept"] + names, b, se, t, pv):
    print(f"{name:10s} {bj:9.4f} {sj:8.4f} {tj:7.3f} {pj:7.4f}")
# <</coefficients>>

fit = sm.OLS(y, X).fit()
assert np.allclose(fit.params, b) and np.allclose(fit.bse, se) and np.allclose(fit.pvalues, pv)
tcrit = stats.t.ppf(0.975, n - p)

# <<overall>>
sst = np.sum((y - y.mean()) ** 2)
R2 = 1 - sse / sst
F_all = (R2 / (p - 1)) / ((1 - R2) / (n - p))

def extra_F(keep):
    """F test that the coefficients not in `keep` (columns of X) are zero."""
    Xr = X[:, keep]
    br, *_ = np.linalg.lstsq(Xr, y, rcond=None)
    sse_r = np.sum((y - Xr @ br) ** 2)
    q = p - len(keep)
    F = ((sse_r - sse) / q) / s2
    return F, stats.f.sf(F, q, n - p)

F_group, p_group = extra_F([0, 1, 3])          # drop hs_grad, white, urban
print(f"overall: R2 = {R2:.3f}, F = {F_all:.2f} on ({p - 1}, {n - p})")
print(f"hs_grad, white, urban: F = {F_group:.3f}, p = {p_group:.3f}")
# <</overall>>

assert np.isclose(F_all, fit.fvalue)
assert np.isclose(F_group, float(fit.f_test("x2 = 0, x4 = 0, x5 = 0").fvalue))
# t^2 = F for each single coefficient
for j in range(1, p):
    Fj, pj = extra_F([k for k in range(p) if k != j])
    assert np.isclose(Fj, t[j] ** 2) and np.isclose(pj, pv[j])
# the group statistic dominates each of its t^2 / q: q F = max over directions of t(a)^2
q = 3
idx = [2, 4, 5]
Wg = C[np.ix_(idx, idx)]
qF = b[idx] @ np.linalg.solve(Wg, b[idx]) / s2
assert np.isclose(qF, q * F_group)
a_star = np.linalg.solve(Wg, b[idx])
t_star2 = (a_star @ b[idx]) ** 2 / (s2 * a_star @ Wg @ a_star)
assert np.isclose(t_star2, qF) and all(t[j] ** 2 <= qF for j in idx)
assert np.sum(pv[1:] < 0.05) == 1 and pv[3] < 0.01 and fit.f_pvalue < 1e-6 and p_group > 0.5

# ---- (2) Longley without the time trend ----------------------------------------------------
# <<longley>>
lg = sm.datasets.longley.load_pandas()
Z = sm.add_constant(lg.exog.drop(columns="YEAR"))
res = sm.OLS(lg.endog, Z).fit()
joint = res.f_test("UNEMP = 0, POP = 0")
print(res.tvalues[["UNEMP", "POP"]].round(3).to_dict(), res.pvalues[["UNEMP", "POP"]].round(3).to_dict())
print(f"joint F = {float(joint.fvalue):.3f}, p = {float(joint.pvalue):.4f}")
# <</longley>>

cov = res.cov_params().loc[["UNEMP", "POP"], ["UNEMP", "POP"]].to_numpy()
rho = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
t2 = res.tvalues[["UNEMP", "POP"]].to_numpy()
nu = int(res.df_resid)
Fcrit2 = stats.f.ppf(0.95, 2, nu)
tc = stats.t.ppf(0.975, nu)
Rm = np.array([[1, rho], [rho, 1]])
assert np.isclose(t2 @ np.linalg.solve(Rm, t2) / 2, float(joint.fvalue))
assert np.all(np.abs(t2) < tc) and float(joint.pvalue) < 0.05
assert rho < -0.5 and nu == 10

gen = Generated("ch11", "coefficient_tests")
gen.int("n", n)
gen.int("p", p)
gen.int("df", n - p)
for name, bj, sj, tj, pj in zip(["int"] + names, b, se, t, pv):
    gen.num(f"b_{name}", bj, 4)
    gen.num(f"se_{name}", sj, 4)
    gen.num(f"t_{name}", tj, 3)
    gen.num(f"p_{name}", pj, 3)
gen.num("tcrit", tcrit, 3)
gen.num("s2", s2, 3)
gen.num("R2", R2, 3)
gen.num("F_all", F_all, 2)
gen.num("p_all", fit.f_pvalue, 1, sci=True)
gen.num("F_group", F_group, 3)
gen.num("p_group", p_group, 3)
gen.num("lg_t_unemp", t2[0], 3)
gen.num("lg_t_pop", t2[1], 3)
gen.num("lg_p_unemp", res.pvalues["UNEMP"], 3)
gen.num("lg_p_pop", res.pvalues["POP"], 3)
gen.num("lg_F", float(joint.fvalue), 3)
gen.num("lg_p", float(joint.pvalue), 4)
gen.num("lg_rho", rho, 3)
# the same two t statistics with the sign of the correlation reversed (exercise exr-glh-longley-F)
Rm_flip = np.array([[1, -rho], [-rho, 1]])
F_flip = t2 @ np.linalg.solve(Rm_flip, t2) / 2
assert F_flip < 1 < float(joint.fvalue)
gen.num("lg_F_flip", F_flip, 3)
gen.num("lg_tcrit", tc, 3)
gen.num("lg_Fcrit", Fcrit2, 3)
gen.int("lg_df", nu)
gen.write()

# ---- figure: acceptance regions in the plane of the two t statistics -----------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 3.6))
# F <= Fcrit  <=>  t' R^{-1} t <= 2 Fcrit: an ellipse with axes along the eigenvectors of R
vals, vecs = np.linalg.eigh(Rm)
ang = np.linspace(0, 2 * np.pi, 400)
circle = np.stack([np.cos(ang), np.sin(ang)])
ellipse = vecs @ (np.sqrt(2 * Fcrit2 * vals)[:, None] * circle)
ax.fill(ellipse[0], ellipse[1], color=COLORS["accent"], alpha=0.12, linewidth=0)
ax.plot(ellipse[0], ellipse[1], color=COLORS["accent"], label="joint F test accepts")
ax.plot([-tc, tc, tc, -tc, -tc], [-tc, -tc, tc, tc, -tc], color=COLORS["second"],
        label="both t tests accept")
ax.plot(t2[0], t2[1], "o", color=COLORS["ink"], markersize=4)
ax.annotate("Longley\ndata", (t2[0], t2[1]), xytext=(-1.2, -1.75), fontsize=8, ha="right", va="center",
            arrowprops=dict(arrowstyle="-", color=COLORS["muted"], linewidth=0.6))
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.axvline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
lim = 4.2
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.set_xlabel(r"$t$ statistic, unemployment")
ax.set_ylabel(r"$t$ statistic, population")
ax.legend(loc="upper left", frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch11", "joint_vs_individual"))
