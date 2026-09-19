"""Chapter 18, Section 5: a 5 x 5 Latin square.

Synthetic experiment: five keyboard layouts are compared by five typists, each typing in five
sessions. Rows are typists, columns are sessions (practice improves speed over the sessions),
and the Latin square gives each layout once to each typist and once in each session.
Response: words per minute on a standard passage.
"""
import numpy as np
from scipy import stats

from regbook import Generated

# <<design>>
rng = np.random.default_rng(62)
a = 5
cyclic = (np.arange(a)[:, None] + np.arange(a)[None, :]) % a    # a standard square
square = cyclic[rng.permutation(a)][:, rng.permutation(a)]      # permute rows and columns
square = rng.permutation(a)[square]                             # and relabel the layouts
print("layout used by typist (row) in session (column):")
print(square + 1)
# <</design>>

assert all(sorted(square[i]) == list(range(a)) for i in range(a))
assert all(sorted(square[:, j]) == list(range(a)) for j in range(a))

# <<data>>
typist = rng.normal(0, 8, a)                          # typists differ a lot
session = np.array([-3.0, -1.0, 0.0, 1.5, 2.5])       # practice effect over sessions
layout = np.array([0.0, 2.0, -1.0, 3.0, 1.0])
y = np.round(55 + typist[:, None] + session[None, :] + layout[square]
             + rng.normal(0, 1.5, (a, a)), 1)
print(y)
# <</data>>

# <<anova>>
grand = y.mean()
row_m, col_m = y.mean(axis=1), y.mean(axis=0)
trt_m = np.array([y[square == k].mean() for k in range(a)])
ss_row = a * np.sum((row_m - grand) ** 2)
ss_col = a * np.sum((col_m - grand) ** 2)
ss_trt = a * np.sum((trt_m - grand) ** 2)
ss_tot = np.sum((y - grand) ** 2)
ss_err = ss_tot - ss_row - ss_col - ss_trt
df_err = (a - 1) * (a - 2)
ms_err = ss_err / df_err
F_trt = ss_trt / (a - 1) / ms_err
print(f"SS rows {ss_row:.1f}, columns {ss_col:.1f}, layouts {ss_trt:.1f}, error {ss_err:.2f} on {df_err} df")
print(f"F for layouts = {F_trt:.2f}, p = {stats.f.sf(F_trt, a - 1, df_err):.2g}")
print("layout means", trt_m.round(2), " se of a difference", np.sqrt(2 * ms_err / a).round(3))
# <</anova>>

# checks: the residual is y - row - col - trt + 2 grand, and the three centred spaces are orthogonal
resid = y - row_m[:, None] - col_m[None, :] - trt_m[square] + 2 * grand
assert np.isclose(np.sum(resid ** 2), ss_err)
n = a * a
R = np.kron(np.eye(a), np.ones((a, 1)))
C = np.kron(np.ones((a, 1)), np.eye(a))
T = np.eye(a)[square.ravel()]
center = lambda A: A - A.mean(axis=0)
for A1, A2 in [(R, C), (R, T), (C, T)]:
    assert np.allclose(center(A1).T @ center(A2), 0)
X = np.column_stack([np.ones(n), R, C, T])
coef = np.linalg.lstsq(X, y.ravel(), rcond=None)[0]
assert np.isclose(np.sum((y.ravel() - X @ coef) ** 2), ss_err)
assert np.linalg.matrix_rank(X) == 3 * a - 2
p_trt = stats.f.sf(F_trt, a - 1, df_err)
assert p_trt < 0.01

# <<efficiency>>
ms_row, ms_col = ss_row / (a - 1), ss_col / (a - 1)
re_rcbd_rows = (ms_col + (a - 1) * ms_err) / (a * ms_err)       # typists as the only blocks
re_rcbd_cols = (ms_row + (a - 1) * ms_err) / (a * ms_err)       # sessions as the only blocks
re_crd = (ms_row + ms_col + (a - 1) * ms_err) / ((a + 1) * ms_err)
print(f"relative efficiency: against blocking on typists only {re_rcbd_rows:.2f},"
      f" on sessions only {re_rcbd_cols:.2f}, against no blocking {re_crd:.2f}")
# <</efficiency>>
assert re_rcbd_cols > re_rcbd_rows > 1

gen = Generated("ch18", "latin_square")
gen.int("dferr", df_err)
for k in range(a):
    gen.num(f"mean{k + 1}", trt_m[k], 2)
gen.num("ssrow", ss_row, 1)
gen.num("sscol", ss_col, 1)
gen.num("sstrt", ss_trt, 1)
gen.num("sserr", ss_err, 2)
gen.num("sstot", ss_tot, 1)
gen.num("msrow", ms_row, 1)
gen.num("mscol", ms_col, 2)
gen.num("mstrt", ss_trt / (a - 1), 2)
gen.num("mserr", ms_err, 3)
gen.num("Frow", ms_row / ms_err, 1)
gen.num("Fcol", ms_col / ms_err, 1)
gen.num("Ftrt", F_trt, 2)
gen.num("ptrt", p_trt, 4)
gen.num("sed", np.sqrt(2 * ms_err / a), 3)
gen.num("rerows", re_rcbd_rows, 2)
gen.num("recols", re_rcbd_cols, 2)
gen.num("recrd", re_crd, 2)
gen.text("square", r" \\ ".join(" & ".join(str(v + 1) for v in row) for row in square))
gen.write()
