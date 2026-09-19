"""Chapter 18, Section 5: a balanced incomplete block design.

Synthetic tasting experiment: seven recipes of a sports drink, fourteen tasters, and each
taster can judge only three recipes. The blocks are the lines of two different Fano planes
(the cyclic difference sets {0, 1, 3} and {0, 1, 5} modulo 7), so every pair of recipes is
tasted together by exactly two tasters. Response: a rating on a 0-10 scale.
"""
import itertools

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from regbook import Generated

# <<design>>
t, k = 7, 3
blocks = [sorted((d + s) % t for d in base) for base in ([0, 1, 3], [0, 1, 5]) for s in range(t)]
b = len(blocks)                                       # 14 tasters
N = np.zeros((t, b), dtype=int)                       # incidence: recipe i tasted by taster j
for j, blk in enumerate(blocks):
    N[blk, j] = 1
r = N.sum(axis=1)[0]
lam = (N @ N.T)[0, 1]
print("r =", r, " lambda =", lam, " NN^T =\n", N @ N.T)
# <</design>>

assert np.all(N.sum(axis=1) == r) and np.all(N.sum(axis=0) == k)
assert np.all((N @ N.T)[~np.eye(t, dtype=bool)] == lam)
assert b * k == t * r and lam * (t - 1) == r * (k - 1)

# <<data>>
rng = np.random.default_rng(77)
tau = np.array([0.0, 0.4, -0.3, 0.8, 0.1, -0.6, 0.5])   # recipe effects
taster = rng.normal(0, 1.2, b)                          # tasters use the scale differently
rows = [(i, j) for j, blk in enumerate(blocks) for i in blk]
recipe = np.array([i for i, j in rows])
block = np.array([j for i, j in rows])
y = np.round(6 + tau[recipe] + taster[block] + rng.normal(0, 0.5, len(rows)), 1)
n = len(y)
# <</data>>

# <<intrablock>>
T_tot = np.bincount(recipe, weights=y, minlength=t)          # recipe totals
B_tot = np.bincount(block, weights=y, minlength=b)           # taster totals
Q = T_tot - N @ B_tot / k                                    # adjusted recipe totals
tau_hat = k * Q / (lam * t)                                  # intra-block estimates (sum to zero)
ss_trt_adj = k / (lam * t) * np.sum(Q ** 2)
ss_blk = np.sum(B_tot ** 2) / k - y.sum() ** 2 / n           # blocks, ignoring recipes
ss_tot = np.sum((y - y.mean()) ** 2)
ss_err = ss_tot - ss_blk - ss_trt_adj
df_err = n - b - t + 1
ms_err = ss_err / df_err
F = ss_trt_adj / (t - 1) / ms_err
se_diff = np.sqrt(2 * k * ms_err / (lam * t))
eff = lam * t / (r * k)                                      # efficiency factor
print("raw recipe means     ", (T_tot / r).round(2))
print("adjusted recipe means", (y.mean() + tau_hat).round(2))
print(f"adjusted recipe SS {ss_trt_adj:.3f}, error SS {ss_err:.3f} on {df_err} df,"
      f" F = {F:.2f}, p = {stats.f.sf(F, t - 1, df_err):.2g}")
print(f"se of a difference {se_diff:.3f}; efficiency factor {eff:.4f}")
# <</intrablock>>

# checks against least squares on the full model
df = pd.DataFrame({"y": y, "rec": recipe, "tas": block})
full = smf.ols("y ~ C(tas) + C(rec)", df).fit()
red = smf.ols("y ~ C(tas)", df).fit()
assert np.isclose(full.ssr, ss_err)
assert np.isclose(red.ssr - full.ssr, ss_trt_adj)
coef = np.r_[0.0, [full.params[f"C(rec)[T.{i}]"] for i in range(1, t)]]
assert np.allclose(coef - coef.mean(), tau_hat)              # same contrasts
Lc = np.zeros(len(full.params))
Lc[list(full.params.index).index("C(rec)[T.1]")] = 1
Lc[list(full.params.index).index("C(rec)[T.2]")] = -1
assert np.isclose(np.sqrt(Lc @ full.cov_params().to_numpy() @ Lc), se_diff)
# information matrix C = r I - N N^T / k = (lambda t / k)(I - J/t)
Cmat = r * np.eye(t) - N @ N.T / k
assert np.allclose(Cmat, lam * t / k * (np.eye(t) - np.ones((t, t)) / t))
p_F = stats.f.sf(F, t - 1, df_err)
assert p_F < 0.01
assert np.argmax(T_tot) == 5 and tau_hat[5] < np.median(tau_hat)   # recipe 6: best raw, not adjusted
# the raw means are distorted by which tasters each recipe met
raw_dev = T_tot / r - y.mean()
gap = np.max(np.abs(raw_dev - tau_hat))
se_rcbd = np.sqrt(2 * ms_err / r)
assert np.isclose(se_diff, se_rcbd / np.sqrt(eff))

gen = Generated("ch18", "bibd")
gen.int("b", b)
gen.int("r", r)
gen.int("lam", lam)
gen.int("n", n)
gen.int("dferr", df_err)
for i in range(t):
    gen.num(f"raw{i + 1}", T_tot[i] / r, 2)
    gen.num(f"adj{i + 1}", y.mean() + tau_hat[i], 2)
    gen.num(f"Q{i + 1}", Q[i], 2)
gen.num("grand", y.mean(), 3)
gen.num("sstrt", ss_trt_adj, 2)
gen.num("ssblk", ss_blk, 2)
gen.num("sserr", ss_err, 2)
gen.num("sstot", ss_tot, 2)
gen.num("mserr", ms_err, 4)
gen.num("F", F, 2)
gen.num("p", p_F, 4)
gen.num("sediff", se_diff, 3)
gen.num("sercbd", se_rcbd, 3)
gen.num("eff", eff, 4)
gen.num("gap", gap, 2)
gen.write()
