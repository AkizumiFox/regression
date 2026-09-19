"""Chapter 14, Section 6: the pure-error lack-of-fit test.

Stack loss data (Brownlee; statsmodels.datasets.stackloss, public domain): 21 days of a
plant oxidizing ammonia. Stack loss against air flow alone has 7 distinct air-flow values,
so a straight line (and a parabola) can be tested against the cell-means model.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "lack_of_fit", prefix="lof")

# <<lof>>
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
x = data["AIRFLOW"].to_numpy()
n = len(y)
levels, groups = np.unique(x, return_inverse=True)   # c distinct rows of X
c = len(levels)
Z = np.eye(c)[groups]                                 # cell-means model matrix

def sse(A):
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return np.sum((y - A @ coef) ** 2)

ss_pe = sse(Z)                                        # pure error, n - c df
for name, X in [("line", np.column_stack([np.ones(n), x])),
                ("parabola", np.column_stack([np.ones(n), x, x ** 2]))]:
    r = np.linalg.matrix_rank(X)
    ss_lf = sse(X) - ss_pe                            # lack of fit, c - r df
    F = (ss_lf / (c - r)) / (ss_pe / (n - c))
    print(f"{name:8s}: SSE {sse(X):7.2f} = SSLF {ss_lf:7.2f} ({c - r} df)"
          f" + SSPE {ss_pe:6.2f} ({n - c} df);  F = {F:.2f}, p = {stats.f.sf(F, c - r, n - c):.4f}")
# <</lof>>

X1 = np.column_stack([np.ones(n), x])
X2 = np.column_stack([np.ones(n), x, x ** 2])
coef1, *_ = np.linalg.lstsq(X1, y, rcond=None)
coef2, *_ = np.linalg.lstsq(X2, y, rcond=None)
means = np.array([y[groups == g].mean() for g in range(c)])
counts = np.bincount(groups)
# SSLF as a weighted sum over the c groups
fit_levels = coef1[0] + coef1[1] * levels
sslf_weighted = np.sum(counts * (means - fit_levels) ** 2)
assert np.isclose(sslf_weighted, sse(X1) - ss_pe)
# pure error is the within-group sum of squares
assert np.isclose(ss_pe, sum(np.sum((y[groups == g] - means[g]) ** 2) for g in range(c)))
F1 = ((sse(X1) - ss_pe) / (c - 2)) / (ss_pe / (n - c))
F2 = ((sse(X2) - ss_pe) / (c - 3)) / (ss_pe / (n - c))
# the same F from statsmodels' nested-model comparison
cmp = sm.OLS(y, Z).fit().compare_f_test(sm.OLS(y, X1).fit())
assert np.isclose(cmp[0], F1)
s2_line = sse(X1) / (n - 2)
s2_pe = ss_pe / (n - c)
gen.int("n", n)
gen.int("c", c)
gen.int("dfpe", n - c)
gen.num("sspe", ss_pe, 2)
gen.num("mspe", s2_pe, 2)
gen.num("sse1", sse(X1), 2)
gen.num("sslf1", sse(X1) - ss_pe, 2)
gen.num("mslf1", (sse(X1) - ss_pe) / (c - 2), 2)
gen.num("s2line", s2_line, 2)
gen.num("F1", F1, 2)
gen.num("p1", stats.f.sf(F1, c - 2, n - c), 4)
gen.num("sse2", sse(X2), 2)
gen.num("sslf2", sse(X2) - ss_pe, 2)
gen.num("F2", F2, 2)
gen.num("p2", stats.f.sf(F2, c - 3, n - c), 4)
gen.num("fcrit1", stats.f.ppf(0.95, c - 2, n - c), 2)
gen.num("b0", coef1[0], 2)
gen.num("b1", coef1[1], 3)
gen.num("R2line", 1 - sse(X1) / np.sum((y - y.mean()) ** 2), 3)
gen.num("mean70", means[list(levels).index(70.0)], 1)
gen.num("mean62", means[list(levels).index(62.0)], 1)
gen.write()

# ---- figure -------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.6))
jitter = np.zeros(n)
for g in range(c):                                   # spread replicates sideways a little
    idx = np.where(groups == g)[0]
    jitter[idx] = np.linspace(-0.35, 0.35, len(idx)) if len(idx) > 1 else 0
ax.scatter(x + jitter, y, s=12, color=COLORS["accent"], alpha=0.8, linewidths=0, label="days")
ax.scatter(levels, means, s=38, marker="_", color=COLORS["ink"], linewidths=1.5,
           label="cell means")
xs = np.linspace(48, 82, 200)
ax.plot(xs, coef1[0] + coef1[1] * xs, color=COLORS["second"], label="least squares line")
ax.plot(xs, coef2[0] + coef2[1] * xs + coef2[2] * xs ** 2, color=COLORS["third"], ls="--",
        label="parabola")
ax.set_xlabel("air flow")
ax.set_ylabel("stack loss")
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch14", "lack_of_fit"))
