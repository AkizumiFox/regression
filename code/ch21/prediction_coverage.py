"""Chapter 21, Section 1: prediction intervals when the errors are not normal.

The usual interval yhat_0 +/- t s sqrt(1 + h_00) covers a new observation with probability
tending to Pr(|eps_0| <= z sigma) as n grows (prp-ci-pi-limit, Section 12.3). This
script computes that limit, and its split between the two tails, for several error laws with
variance one; Section 12.3 checks the limit by simulation.
"""
import numpy as np
from scipy import stats

from regbook import Generated

# <<limits>>
laws = {
    "normal": stats.norm(),
    "t(5)": stats.t(5, scale=np.sqrt(3 / 5)),
    "t(3)": stats.t(3, scale=np.sqrt(1 / 3)),
    "uniform": stats.uniform(loc=-np.sqrt(3), scale=2 * np.sqrt(3)),
    "centred exponential": stats.expon(loc=-1.0),
}
for name, law in laws.items():
    assert abs(law.var() - 1) < 1e-9 and abs(law.mean()) < 1e-9
    row = []
    for level in (0.80, 0.95, 0.99):
        z = stats.norm.ppf(0.5 + level / 2)
        below, above = law.cdf(-z), law.sf(z)            # miss low, miss high
        row.append(f"{1 - below - above:.3f} ({below:.3f}, {above:.3f})")
    print(f"{name:20s}", "  ".join(row))
# <</limits>>

gen = Generated("ch21", "prediction_coverage", prefix="pred")
keys = {"normal": "norm", "t(5)": "t5", "t(3)": "t3", "uniform": "unif", "centred exponential": "exp"}
limit = {}
for name, law in laws.items():
    for level in (0.80, 0.95, 0.99):
        z = stats.norm.ppf(0.5 + level / 2)
        below, above = law.cdf(-z), law.sf(z)
        tag = f"{keys[name]}_{int(round(100 * level))}"
        limit[tag] = 1 - below - above
        gen.num(f"cov_{tag}", 1 - below - above, 3)
        gen.num(f"low_{tag}", below, 3)
        gen.num(f"high_{tag}", above, 3)
assert limit["exp_95"] > 0.94 and limit["exp_99"] < 0.985 and limit["t3_99"] < 0.985
assert limit["t3_80"] > 0.88 and limit["unif_80"] < 0.75
# t(3) misses a 99% interval twice as often
assert 2 < (1 - limit["t3_99"]) / 0.01 < 2.2
# all the exponential misses are above
assert 5 < (1 - limit["exp_99"]) / 0.005 < 6

gen.write()
