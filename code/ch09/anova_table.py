"""Chapter 9, Section 1: the analysis of variance identity for a regression.

Murder rate of the 50 US states (2009) on poverty, single-parent and urban
percentages, as in Chapter 6. Public-domain data shipped with statsmodels
(statsmodels.datasets.statecrime).
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")

# <<table>>
y = data["murder"].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), data["poverty"], data["single"], data["urban"]])
r = np.linalg.matrix_rank(X)

Q, _ = np.linalg.qr(X)                 # orthonormal basis of C(X); first column spans 1
z = Q.T @ y                            # coordinates of My in that basis
ss_mean = z[0] ** 2                    # n * ybar^2
ss_reg = np.sum(z[1:] ** 2)            # ||(M - P1) y||^2, r - 1 squared coordinates
ss_res = y @ y - np.sum(z ** 2)        # ||(I - M) y||^2

rows = [("mean", 1, ss_mean), ("regression", r - 1, ss_reg),
        ("residual", n - r, ss_res), ("total", n, y @ y)]
for name, df, ss in rows:
    print(f"{name:11s} df = {df:3d}   SS = {ss:9.2f}   MS = {ss / df:8.3f}")
F = (ss_reg / (r - 1)) / (ss_res / (n - r))
print(f"corrected total = {ss_reg + ss_res:.2f},  F = {F:.2f}")
# <</table>>

ybar = y.mean()
assert np.isclose(ss_mean, n * ybar ** 2)
assert np.isclose(ss_reg + ss_res, np.sum((y - ybar) ** 2))
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
yhat = X @ beta
assert np.isclose(ss_reg, np.sum((yhat - ybar) ** 2))
# computational formulas
assert np.isclose(ss_res, y @ y - beta @ X.T @ y)
assert np.isclose(ss_reg, beta @ X.T @ y - n * ybar ** 2)
fit = sm.OLS(y, X).fit()
assert np.isclose(fit.ess, ss_reg) and np.isclose(fit.ssr, ss_res) and np.isclose(fit.fvalue, F)
p_value = stats.f.sf(F, r - 1, n - r)

# ---- no intercept: the corrected identity fails by -2 ybar 1'e ---------------
# <<nointercept>>
X0 = X[:, 1:]                                        # drop the intercept column
b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
e0 = y - X0 @ b0
sst = np.sum((y - ybar) ** 2)
ssr_star = np.sum((X0 @ b0 - ybar) ** 2)
sse0 = e0 @ e0
print(f"SST = {sst:.2f}  but  SSR* + SSE = {ssr_star + sse0:.2f}")
print(f"sum of residuals = {e0.sum():.3f},  -2 ybar * sum = {-2 * ybar * e0.sum():.2f}")
# <</nointercept>>
assert np.isclose(sst, ssr_star + sse0 - 2 * ybar * e0.sum())
assert np.isclose(y @ y, np.sum((X0 @ b0) ** 2) + sse0)
r2_uncentred = np.sum((X0 @ b0) ** 2) / (y @ y)

gen = Generated("ch09", "anova_table", prefix="tab")
gen.int("n", n)
gen.int("r", r)
gen.num("ybar", ybar, 3)
gen.num("mean", ss_mean, 2)
gen.num("reg", ss_reg, 2)
gen.num("res", ss_res, 2)
gen.num("sst", ss_reg + ss_res, 2)
gen.num("total", y @ y, 2)
gen.num("msreg", ss_reg / (r - 1), 2)
gen.num("msres", ss_res / (n - r), 3)
gen.num("msmean", ss_mean, 2)
gen.num("F", F, 2)
gen.num("p", p_value, 1, sci=True)
gen.num("R2", ss_reg / (ss_reg + ss_res), 4)
gen.num("ni:ssrstar", ssr_star, 2)
gen.num("ni:sse", sse0, 2)
gen.num("ni:sum", ssr_star + sse0, 2)
gen.num("ni:esum", e0.sum(), 3)
gen.num("ni:corr", -2 * ybar * e0.sum(), 2)
gen.num("ni:R2u", r2_uncentred, 4)
gen.write()
