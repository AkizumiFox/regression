"""Chapter 18, Section 4: a randomized complete block design.

Synthetic experiment: five wood adhesives are compared on six boards of timber (the blocks).
Each board is cut into five strips and the adhesives are assigned to the strips of every board
at random. Response: shear strength of the glued joint in MPa. Used again in Section 6
(missing observations).
"""
import itertools
import math

import numpy as np
from scipy import stats

from regbook import Generated

# <<data>>
rng = np.random.default_rng(606)
t, b = 5, 6                                           # adhesives, boards
plan = np.array([rng.permutation(t) for _ in range(b)])   # which strip of each board gets which adhesive
tau = np.array([0.0, 0.40, 0.65, 0.15, 0.80])         # adhesive effects (unknown in practice)
board = rng.normal(0, 0.75, b)                        # board effects
y = np.round(9 + tau[:, None] + board[None, :] + rng.normal(0, 0.3, (t, b)), 2)
print("adhesive given to strips 1-5 (rows) of each board (columns):")
print(plan.T + 1)
print("shear strength, adhesives in rows, boards in columns:")
print(y)
# <</data>>

# <<anova>>
grand = y.mean()
trt_mean, blk_mean = y.mean(axis=1), y.mean(axis=0)
resid = y - trt_mean[:, None] - blk_mean[None, :] + grand
ss_trt = b * np.sum((trt_mean - grand) ** 2)
ss_blk = t * np.sum((blk_mean - grand) ** 2)
ss_err = np.sum(resid ** 2)
df_err = (t - 1) * (b - 1)
ms_trt, ms_blk, ms_err = ss_trt / (t - 1), ss_blk / (b - 1), ss_err / df_err
F_trt = ms_trt / ms_err
print(f"treatments SS {ss_trt:.2f}, blocks SS {ss_blk:.2f}, error SS {ss_err:.2f} on {df_err} df")
print(f"F for adhesives = {F_trt:.2f}, p = {stats.f.sf(F_trt, t - 1, df_err):.2g}")
print("adhesive means", trt_mean.round(2), " se of a difference", np.sqrt(2 * ms_err / b).round(3))
# <</anova>>

# checks against the projection decomposition
n = t * b
T_ind = np.kron(np.eye(t), np.ones((b, 1)))            # y.ravel() is adhesive-major
B_ind = np.kron(np.ones((t, 1)), np.eye(b))
yv = y.ravel()
def proj(*blocks):
    U, s, _ = np.linalg.svd(np.column_stack(blocks), full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T
P0, PT, PB = proj(np.ones(n)), proj(T_ind), proj(B_ind)
assert np.allclose(PT @ PB, P0)                        # balance: centred spaces orthogonal
assert np.isclose(yv @ (PT - P0) @ yv, ss_trt)
assert np.isclose(yv @ (np.eye(n) - PT - PB + P0) @ yv, ss_err)
p_trt = stats.f.sf(F_trt, t - 1, df_err)
assert p_trt < 1e-4

# <<efficiency>>
s2_crd = ((b - 1) * ms_blk + b * (t - 1) * ms_err) / (b * t - 1)   # estimated CRD error variance
re = s2_crd / ms_err
nu_b, nu_c = df_err, b * t - t                        # error df: blocked, completely randomized
re_fisher = re * (nu_b + 1) * (nu_c + 3) / ((nu_b + 3) * (nu_c + 1))
F_crd = ms_trt / s2_crd                               # what the treatment F would roughly have been
print(f"estimated CRD error variance {s2_crd:.2f} against {ms_err:.2f}:"
      f" relative efficiency {re:.2f} ({re_fisher:.2f} with Fisher's df factor)")
# <</efficiency>>
assert re > 1

# <<randomization>>
def f_stat(yy):
    tm, bm, gm = yy.mean(axis=1), yy.mean(axis=0), yy.mean()
    sst = b * np.sum((tm - gm) ** 2)
    sse = np.sum((yy - tm[:, None] - bm[None, :] + gm) ** 2)
    return (sst / (t - 1)) / (sse / df_err)

rng_perm = np.random.default_rng(7)
draws = 20_000
F_perm = np.empty(draws)
for r in range(draws):
    # re-randomize: permute the adhesive labels independently within every board
    perm = np.array([rng_perm.permutation(t) for _ in range(b)]).T
    F_perm[r] = f_stat(np.take_along_axis(y, perm, axis=0))
p_rand = (1 + np.sum(F_perm >= F_trt)) / (draws + 1)
print(f"randomization p-value {p_rand:.4f}  (F-test p-value {stats.f.sf(F_trt, t - 1, df_err):.2g})")
# <</randomization>>
assert p_rand < 0.001

# Tukey HSD on the adhesive means (balanced one-way comparisons with block-free error)
q_crit = stats.studentized_range.ppf(0.95, t, df_err)
hsd = q_crit * np.sqrt(ms_err / b)
sig = [(i + 1, k + 1) for i, k in itertools.combinations(range(t), 2)
       if abs(trt_mean[i] - trt_mean[k]) > hsd]
print("Tukey HSD", round(hsd, 3), " significant pairs", sig)

# Tukey's one-degree-of-freedom test for nonadditivity (Chapter 16)
a_hat, b_hat = trt_mean - grand, blk_mean - grand
ss_nonadd = (a_hat @ y @ b_hat) ** 2 / (np.sum(a_hat ** 2) * np.sum(b_hat ** 2))
F_nonadd = ss_nonadd / ((ss_err - ss_nonadd) / (df_err - 1))
p_nonadd = stats.f.sf(F_nonadd, 1, df_err - 1)
assert p_nonadd > 0.05

gen = Generated("ch18", "rcbd")
gen.int("t", t)
gen.int("b", b)
gen.int("dferr", df_err)
for i in range(t):
    gen.num(f"mean{i + 1}", trt_mean[i], 3)
gen.num("grand", grand, 3)
gen.num("sstrt", ss_trt, 3)
gen.num("ssblk", ss_blk, 3)
gen.num("sserr", ss_err, 3)
gen.num("sstot", np.sum((y - grand) ** 2), 3)
gen.num("mstrt", ms_trt, 4)
gen.num("msblk", ms_blk, 4)
gen.num("mserr", ms_err, 4)
gen.num("Ftrt", F_trt, 2)
gen.num("Fblk", ms_blk / ms_err, 2)
gen.num("ptrt", p_trt, 1, sci=True)
gen.num("sed", np.sqrt(2 * ms_err / b), 3)
gen.num("s2crd", s2_crd, 4)
gen.num("re", re, 2)
gen.num("refisher", re_fisher, 2)
gen.num("Fcrd", F_crd, 2)
gen.num("prand", p_rand, 5)
gen.int("draws", draws)
gen.int("nexceed", int(np.sum(F_perm >= F_trt)))
gen.text("nplans", "2.99\\times 10^{12}")
assert f"{float(math.factorial(t)) ** b:.2e}" == "2.99e+12"
gen.num("hsd", hsd, 3)
gen.text("sigpairs", ", ".join(f"{i}–{k}" for i, k in sig))
gen.num("Fnonadd", F_nonadd, 2)
gen.num("pnonadd", p_nonadd, 2)
for i in range(t):
    gen.text(f"row{i + 1}", " | ".join(f"{v:.2f}" for v in y[i]))
    gen.text(f"plan{i + 1}", " | ".join(str(v + 1) for v in plan.T[i]))
gen.write()
