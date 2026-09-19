"""Chapter 6, Section 6: the fixed-effects (within) estimator is Frisch-Waugh-Lovell.

Grunfeld's investment data: 11 US firms observed annually 1935-1954 (public domain,
statsmodels.datasets.grunfeld).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import Generated

df = sm.datasets.grunfeld.load_pandas().data
y = df["invest"].to_numpy()
X1 = df[["value", "capital"]].to_numpy()

# <<within>>
# (i) least squares with one indicator column per firm ("dummy variables")
D = pd.get_dummies(df["firm"]).to_numpy(dtype=float)
beta_dummies, *_ = np.linalg.lstsq(np.column_stack([X1, D]), y, rcond=None)

# (ii) FWL: projecting onto C(D) replaces each value by its firm mean,
#      so (I - M_D) subtracts firm means -- the "within" transformation
def demean_by(v, groups):
    v = pd.DataFrame(v)
    return (v - v.groupby(groups.to_numpy()).transform("mean")).to_numpy()


X1_within = demean_by(X1, df["firm"])
y_within = demean_by(y, df["firm"]).ravel()
beta_within, *_ = np.linalg.lstsq(X1_within, y_within, rcond=None)

print("with firm indicators:", beta_dummies[:2])
print("within transformation:", beta_within)
# <</within>>

assert np.allclose(beta_dummies[:2], beta_within, rtol=1e-9)
# the within transformation really is I - M_D
M_D = D @ np.linalg.solve(D.T @ D, D.T)
assert np.allclose((np.eye(len(y)) - M_D) @ y, y_within)
beta_pooled, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(y)), X1]), y, rcond=None)

gen = Generated("ch06", "within_estimator", prefix="fe")
gen.int("n", len(y))
gen.int("firms", D.shape[1])
gen.int("years", df["year"].nunique())
gen.num("value", beta_within[0], 4)
gen.num("capital", beta_within[1], 4)
gen.num("pooled_value", beta_pooled[1], 4)
gen.num("pooled_capital", beta_pooled[2], 4)
gen.write()
