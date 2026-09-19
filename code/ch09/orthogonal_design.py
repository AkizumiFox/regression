"""Chapter 9, Section 4: orthogonal designs.

Grunfeld's investment panel (11 firms x 20 years, one observation per cell) is a
balanced two-way layout, so the firm and year sums of squares do not depend on the
order of fitting. Adding a covariate destroys the orthogonality. Public-domain data
shipped with statsmodels (statsmodels.datasets.grunfeld). The second part checks the
proportional-frequencies criterion on small layouts.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import Generated


def proj(*blocks):
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T


# <<grunfeld>>
g = sm.datasets.grunfeld.load_pandas().data
y = np.log(g["invest"].to_numpy())
n = len(y)
one = np.ones((n, 1))
F = pd.get_dummies(g["firm"]).to_numpy(float)          # 11 firms
T = pd.get_dummies(g["year"]).to_numpy(float)          # 20 years
P0, PF, PT, PFT = proj(one), proj(one, F), proj(one, T), proj(one, F, T)

print("centred spaces orthogonal:", np.abs((PF - P0) @ (PT - P0)).max() < 1e-10)
print(f"firm:  first {y @ (PF - P0) @ y:8.3f}   after year {y @ (PFT - PT) @ y:8.3f}")
print(f"year:  first {y @ (PT - P0) @ y:8.3f}   after firm {y @ (PFT - PF) @ y:8.3f}")
# <</grunfeld>>
ss_firm, ss_year = y @ (PF - P0) @ y, y @ (PFT - PF) @ y
assert np.isclose(ss_firm, y @ (PFT - PT) @ y)
assert np.isclose(ss_year, y @ (PT - P0) @ y)
assert np.allclose(PFT, P0 + (PF - P0) + (PT - P0))
sse_add = y @ (np.eye(n) - PFT) @ y
sst = y @ (np.eye(n) - P0) @ y
assert np.isclose(ss_firm + ss_year + sse_add, sst)

# <<covariate>>
v = np.log(g["value"].to_numpy())[:, None]             # log market value
PFv, PTv, PFTv = proj(one, F, v), proj(one, T, v), proj(one, F, T, v)
print(f"with log value: firm after year, value {y @ (PFTv - PTv) @ y:8.3f}"
      f"   firm after value {y @ (PFv - proj(one, v)) @ y:8.3f}")
# <</covariate>>
ss_firm_cov = y @ (PFTv - PTv) @ y
ss_firm_after_v = y @ (PFv - proj(one, v)) @ y
assert not np.isclose(ss_firm_cov, ss_firm)

# ---- proportional frequencies ---------------------------------------------------


def layout(counts):
    rows, cols = [], []
    for i in range(counts.shape[0]):
        for j in range(counts.shape[1]):
            rows += [i] * counts[i, j]
            cols += [j] * counts[i, j]
    m = len(rows)
    Ai, Bj = np.eye(counts.shape[0])[rows], np.eye(counts.shape[1])[cols]
    o = np.ones((m, 1))
    Q0 = proj(o)
    return np.abs((proj(o, Ai) - Q0) @ (proj(o, Bj) - Q0)).max()


prop = np.array([[2, 4, 6], [1, 2, 3]])                  # n_ij = n_i. n_.j / n
nonprop = np.array([[2, 4, 6], [1, 3, 2]])
assert layout(prop) < 1e-12
assert layout(nonprop) > 1e-3

gen = Generated("ch09", "orthogonal_design", prefix="orth")
gen.int("n", n)
gen.int("firms", F.shape[1])
gen.int("years", T.shape[1])
gen.num("firm", ss_firm, 3)
gen.num("year", ss_year, 3)
gen.num("sse", sse_add, 3)
gen.num("sst", sst, 3)
gen.num("firm_cov", ss_firm_cov, 3)
gen.num("firm_after_v", ss_firm_after_v, 3)
gen.num("nonprop", layout(nonprop), 3)
gen.write()
