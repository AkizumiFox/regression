"""Chapter 9, Section 5: R^2, adjusted R^2 and partial R^2 on the 2009 US state data.

Public-domain data shipped with statsmodels (statsmodels.datasets.statecrime).
"""
from itertools import combinations

import numpy as np
import statsmodels.api as sm

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")

# <<partial>>
y = data["murder"].to_numpy()
names = ["poverty", "single", "urban"]
n = len(y)


def fit(cols):
    X = np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])
    return sm.OLS(y, X).fit()


full = fit(names)
df_res = full.df_resid
for j, name in enumerate(names, start=1):
    reduced = fit([c for c in names if c != name])
    partial_r2 = (reduced.ssr - full.ssr) / reduced.ssr
    t = full.tvalues[j]
    print(f"{name:8s} partial R^2 = {partial_r2:.4f},  t = {t:6.3f},"
          f"  t^2/(t^2 + df) = {t ** 2 / (t ** 2 + df_res):.4f}")
print(f"R^2 = {full.rsquared:.4f}, adjusted R^2 = {full.rsquared_adj:.4f}")
# <</partial>>

gen = Generated("ch09", "partial_r2", prefix="pr2")
gen.num("R2", full.rsquared, 4)
gen.num("R2adj", full.rsquared_adj, 4)
gen.int("dfres", df_res)
for j, name in enumerate(names, start=1):
    reduced = fit([c for c in names if c != name])
    partial_r2 = (reduced.ssr - full.ssr) / reduced.ssr
    t = full.tvalues[j]
    assert np.isclose(partial_r2, t ** 2 / (t ** 2 + df_res))
    # 1 - R^2 = (1 - R^2_reduced)(1 - partial R^2)
    assert np.isclose(1 - full.rsquared, (1 - reduced.rsquared) * (1 - partial_r2))
    # partial R^2 is the squared correlation of the two residual vectors
    others = np.column_stack([np.ones(n)] + [data[c] for c in names if c != name])
    ry = y - others @ np.linalg.lstsq(others, y, rcond=None)[0]
    xj = data[name].to_numpy()
    rx = xj - others @ np.linalg.lstsq(others, xj, rcond=None)[0]
    assert np.isclose(np.corrcoef(ry, rx)[0, 1] ** 2, partial_r2)
    gen.num(f"{name}:partial", partial_r2, 4)
    gen.num(f"{name}:t", t, 3)
    gen.num(f"{name}:R2red", reduced.rsquared, 4)
    gen.num(f"{name}:adjred", reduced.rsquared_adj, 4)

# adjusted R^2 over all seven models: it rises when a regressor with |t| > 1 enters
table = {}
for size in range(1, 4):
    for S in combinations(names, size):
        m = fit(list(S))
        table[S] = (m.rsquared, m.rsquared_adj)
        gen.num("r2:" + "-".join(s[:3] for s in S), m.rsquared, 4)
        gen.num("adj:" + "-".join(s[:3] for s in S), m.rsquared_adj, 4)
urban_t = full.tvalues[3]
assert abs(urban_t) < 1 and full.rsquared_adj < table[("poverty", "single")][1]
assert full.rsquared > table[("poverty", "single")][0]
best = max(table, key=lambda S: table[S][1])
gen.text("best", ", ".join(best))
# the F statistic for adding urban equals t^2, and adj R^2 rises iff F > 1
F_urban = urban_t ** 2
gen.num("F_urban", F_urban, 3)
# sequential product: 1 - R^2 = prod (1 - partial r^2 of each entering variable)
prod = 1.0
cols = []
for name in names:
    before = fit(cols) if cols else None
    sse_before = before.ssr if before else np.sum((y - y.mean()) ** 2)
    cols.append(name)
    prod *= fit(cols).ssr / sse_before
assert np.isclose(prod, 1 - full.rsquared)
# a block: poverty and urban given single
red = fit(["single"])
q = 2
F_block = ((red.ssr - full.ssr) / q) / full.scale
pr2_block = (red.ssr - full.ssr) / red.ssr
assert np.isclose(pr2_block, q * F_block / (q * F_block + df_res))
gen.num("block:F", F_block, 3)
gen.num("block:pr2", pr2_block, 4)
gen.write()
