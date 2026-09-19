"""Chapter 8, Section 1: age, period and cohort. A linear dependency built into the design.

A repeated survey is fielded in five periods; respondents fall into eleven age bands, and
the birth cohort is period minus age. The linear model with a term for each of age, period
and cohort has a model matrix of rank 3, not 4: one direction of the coefficient vector is
invisible in the mean. Simulated data with a fixed seed.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<design>>
import numpy as np

ages = np.arange(20, 75, 5)                  # 11 age bands: 20, 25, ..., 70
periods = np.arange(2000, 2025, 5)           # 5 survey waves: 2000, ..., 2020
A, P = np.meshgrid(ages, periods, indexing="ij")
age, period = A.ravel() - 45.0, P.ravel() - 2010.0   # centred at 45 and 2010
cohort = period - age                        # birth year minus 1965, exactly
X = np.column_stack([np.ones(age.size), age, period, cohort])
print("n =", X.shape[0], " p =", X.shape[1], " rank =", np.linalg.matrix_rank(X))
print("X @ (0, 1, -1, 1) =", np.abs(X @ [0, 1, -1, 1]).max())
# <</design>>

rng = np.random.default_rng(8)
truth = np.array([50.0, 0.30, -0.10, 0.20])
y = X @ truth + rng.normal(scale=1.5, size=age.size)


# <<stories>>
def fit_without(j):
    """Least squares after deleting column j (one of age, period, cohort)."""
    keep = [k for k in range(4) if k != j]
    b = np.zeros(4)
    b[keep] = np.linalg.lstsq(X[:, keep], y, rcond=None)[0]
    return b

stories = {"no cohort term": fit_without(3),
           "no period term": fit_without(2),
           "no age term": fit_without(1)}
for name, b in stories.items():
    sse = np.sum((y - X @ b) ** 2)
    print(f"{name:15s} b = {np.round(b, 3)}  SSE = {sse:.3f}",
          f" age+period = {b[1] + b[2]:.3f}  period+cohort = {b[2] + b[3]:.3f}")
# <</stories>>

null = np.array([0.0, 1.0, -1.0, 1.0])
assert np.linalg.matrix_rank(X) == 3
assert np.allclose(X @ null, 0)
fits = [X @ b for b in stories.values()]
assert all(np.allclose(f, fits[0]) for f in fits)
bs = list(stories.values())
# any two stories differ by a multiple of the null vector
for b in bs[1:]:
    d = b - bs[0]
    t = d @ null / (null @ null)
    assert np.allclose(d, t * null)
# the identifiable combinations agree; the individual slopes do not
for b in bs:
    assert np.isclose(b[1] + b[2], bs[0][1] + bs[0][2])
    assert np.isclose(b[2] + b[3], bs[0][2] + bs[0][3])
    assert np.isclose(b[1] - b[3], bs[0][1] - bs[0][3])
assert not np.isclose(bs[0][1], bs[1][1])
# the truth is not on the line of least squares solutions, but its identifiable parts are close
d = truth - bs[0]
assert not np.allclose(d, (d @ null / (null @ null)) * null)
assert abs(truth[1] + truth[2] - (bs[0][1] + bs[0][2])) < 0.02
assert abs(truth[2] + truth[3] - (bs[0][2] + bs[0][3])) < 0.02

# factor version (Exercise C1): one indicator per age band, period and cohort
a_idx = (age + 45 - 20) // 5
p_idx = (period + 2010 - 2000) // 5
c_idx = (cohort + 1965 - 1930) // 5
Za, Zp, Zc = (np.eye(int(k.max()) + 1)[k.astype(int)] for k in (a_idx, p_idx, c_idx))
XF = np.column_stack([np.ones(age.size), Za, Zp, Zc])
nA, nP, nC = Za.shape[1], Zp.shape[1], Zc.shape[1]
assert (nA, nP, nC) == (11, 5, 15)
assert np.linalg.matrix_rank(XF) == nA + nP + nC - 3          # null space: 3 + 1 dimensions
ext = np.r_[0, np.arange(nA) - (nA - 1) / 2, -(np.arange(nP) - (nP - 1) / 2), np.arange(nC) - (nC - 1) / 2]
assert np.allclose(XF @ ext, 0)                               # the linear-trend null vector
second = np.r_[0, [1, -2, 1], np.zeros(nA - 3 + nP + nC)]
assert np.linalg.matrix_rank(np.vstack([XF, second])) == np.linalg.matrix_rank(XF)

gen = Generated("ch08", "apc", prefix="apc")
gen.int("n", X.shape[0])
for key, b in zip(["nocoh", "noper", "noage"], bs):
    for i in range(1, 4):
        gen.num(f"{key}:b{i}", b[i], 3)
gen.num("sse", np.sum((y - fits[0]) ** 2), 2)
gen.num("ap", bs[0][1] + bs[0][2], 3)
gen.num("pc", bs[0][2] + bs[0][3], 3)
gen.num("amc", bs[0][1] - bs[0][3], 3)
gen.write()

# ---- figure: the three stories tell different age profiles ---------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.0), sharey=True)
colors = [COLORS["accent"], COLORS["second"], COLORS["third"]]
labels = list(stories)
grids = [(ages - 45.0, "age"), (periods - 2010.0, "period"), (np.arange(1930, 2005, 5) - 1965.0, "birth cohort")]
offsets = [45, 2010, 1965]
for ax, (grid, title), j, off in zip(axes, grids, [1, 2, 3], offsets):
    for b, c, lab in zip(bs, colors, labels):
        ax.plot(grid + off, b[j] * grid, color=c, label=lab)
    ax.axhline(0, color=COLORS["grid"], lw=0.6, zorder=0)
    ax.set_title(f"{title} term")
    ax.set_xlabel(title)
axes[0].set_ylabel("contribution to the mean")
axes[0].legend(frameon=False, fontsize=6.5, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch08", "apc_stories"))
