"""Chapter 9, Section 4: splitting a sum of squares into single-degree-of-freedom pieces.

One-way layout: liberal-conservative self-placement by education (3 levels) in the
1996 American National Election Study (statsmodels.datasets.anes96, public domain).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import Generated

d = sm.datasets.anes96.load_pandas().data
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"])
y = d["selfLR"].to_numpy()

# <<contrasts>>
counts = educ.value_counts(sort=False).to_numpy().astype(float)
means = d.groupby(educ, observed=True)["selfLR"].mean().to_numpy()
grand = y.mean()
ss_treat = np.sum(counts * (means - grand) ** 2)          # SS(educ | 1), 2 df


def contrast_ss(c):
    """(c' ybar)^2 / sum(c_i^2 / n_i): the SS for the hypothesis c' mu = 0."""
    return (c @ means) ** 2 / np.sum(c ** 2 / counts)


def orthogonal(c, h):
    return np.isclose(np.sum(c * h / counts), 0.0)


n2, n3 = counts[1], counts[2]
c1 = np.array([1.0, -n2 / (n2 + n3), -n3 / (n2 + n3)])   # HS against the rest, weighted
c2 = np.array([0.0, 1.0, -1.0])                          # College against Graduate
c1_plain = np.array([1.0, -0.5, -0.5])                   # HS against the plain average
print("orthogonal:", orthogonal(c1, c2), orthogonal(c1_plain, c2))
print(f"treatment SS {ss_treat:.3f} = {contrast_ss(c1):.3f} + {contrast_ss(c2):.3f}")
print(f"with the plain contrast: {contrast_ss(c1_plain):.3f} + {contrast_ss(c2):.3f}")
# <</contrasts>>
assert orthogonal(c1, c2) and not orthogonal(c1_plain, c2)
assert np.isclose(contrast_ss(c1) + contrast_ss(c2), ss_treat)
assert not np.isclose(contrast_ss(c1_plain) + contrast_ss(c2), ss_treat)
# projection check: the contrast SS is ||P_v y||^2 with v = M rho, v_i = c_k / n_k in group k
codes = educ.cat.codes.to_numpy()
for c in (c1, c2, c1_plain):
    v = (c / counts)[codes]
    assert np.isclose((v @ y) ** 2 / (v @ v), contrast_ss(c))

gen = Generated("ch09", "single_df", prefix="sdf")
for k, lab in enumerate(["hs", "col", "grad"]):
    gen.int(f"n:{lab}", counts[k])
    gen.num(f"m:{lab}", means[k], 3)
gen.num("treat", ss_treat, 3)
gen.num("c1", contrast_ss(c1), 3)
gen.num("c2", contrast_ss(c2), 3)
gen.num("c1plain", contrast_ss(c1_plain), 3)
gen.num("sumplain", contrast_ss(c1_plain) + contrast_ss(c2), 3)
gen.num("w2", n2 / (n2 + n3), 3)
gen.num("w3", n3 / (n2 + n3), 3)
gen.write()
