"""Chapter 11, Section 1: testing a reduced model, the region effect in the 2009 state data.

Response: murder rate per 100,000. Full model: intercept, indicators for three of the four
Census regions (Midwest as reference), poverty, single-parent households and urbanization.
Reduced model: the same without the region indicators. Public-domain data shipped with
statsmodels (statsmodels.datasets.statecrime); the Census region of each state is public
information from the US Census Bureau.
"""
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

from regbook import Generated

# <<setup>>
data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
# Census region of each state, in the data's alphabetical order (N, M, S, W)
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])                    # reduced model
D = np.column_stack([(region == g).astype(float) for g in "NSW"])  # Midwest is the reference
X = np.column_stack([X0, D])                                      # full model
# <</setup>>

assert n == 50 and len(codes) == n
assert {g: int(np.sum(region == g)) for g in "NMSW"} == {"N": 9, "M": 12, "S": 16, "W": 13}
assert region[list(data.index).index("Massachusetts")] == "N"
assert region[list(data.index).index("Texas")] == "S"
assert region[list(data.index).index("Nebraska")] == "M"
assert region[list(data.index).index("Hawaii")] == "W"

# <<ftest>>
def sse(Z):
    """Residual sum of squares of the least squares fit of y on the columns of Z."""
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return np.sum((y - Z @ b) ** 2)

r, r0 = np.linalg.matrix_rank(X), np.linalg.matrix_rank(X0)
sse_full, sse_reduced = sse(X), sse(X0)
F = ((sse_reduced - sse_full) / (r - r0)) / (sse_full / (n - r))
p_value = stats.f.sf(F, r - r0, n - r)
crit = stats.f.ppf(0.95, r - r0, n - r)
print(f"SSE0 = {sse_reduced:.2f}, SSE = {sse_full:.2f}, df = ({r - r0}, {n - r})")
print(f"F = {F:.3f}, p = {p_value:.4f}, 5% critical value = {crit:.3f}")
# <</ftest>>

# the same test by statsmodels' model comparison
df = data.assign(region=region)
fit0 = smf.ols("murder ~ poverty + single + urban", df).fit()
fit1 = smf.ols("murder ~ C(region, Treatment('M')) + poverty + single + urban", df).fit()
table = sm.stats.anova_lm(fit0, fit1)
assert np.isclose(table["F"].iloc[1], F) and np.isclose(table["Pr(>F)"].iloc[1], p_value)
assert (r, r0) == (7, 4) and F < crit and 0.05 < p_value < 0.10
# the Northeast coefficient on its own has a small p-value
assert fit1.pvalues["C(region, Treatment('M'))[T.N]"] < 0.02

# the extra sum of squares is the squared length of (M - M0) y
def proj(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

M, M0 = proj(X), proj(X0)
assert np.isclose(y @ (M - M0) @ y, sse_reduced - sse_full)
s2 = sse_full / (n - r)
partial_r2 = (sse_reduced - sse_full) / sse_reduced

gen = Generated("ch11", "region_test")
gen.int("n", n)
gen.int("r", r)
gen.int("r0", r0)
gen.int("df1", r - r0)
gen.int("df2", n - r)
gen.num("sse0", sse_reduced, 2)
gen.num("sse", sse_full, 2)
gen.num("extra", sse_reduced - sse_full, 2)
gen.num("ms_extra", (sse_reduced - sse_full) / (r - r0), 3)
gen.num("s2", s2, 3)
gen.num("F", F, 3)
gen.num("p", p_value, 3)
gen.num("crit", crit, 3)
gen.num("partial_r2", partial_r2, 3)
gen.num("b_north", fit1.params["C(region, Treatment('M'))[T.N]"], 3)
gen.num("t_north", fit1.tvalues["C(region, Treatment('M'))[T.N]"], 3)
gen.num("p_north", fit1.pvalues["C(region, Treatment('M'))[T.N]"], 3)
gen.write()
