"""Chapter 20, Section 5: masking and swamping by a cluster of identical bad cases.

Thirty clean cases follow a line; four identical cases, recorded far to the right and far
below the line, hide one another from every single-case diagnostic. The script checks the
group-deletion formulas and the formulas of the masking proposition, then shows what the
single-case and group diagnostics report.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style


def single_case(X, y):
    """Leverages, residuals, s^2, externally studentized residuals and Cook's distances."""
    n, p = X.shape
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)
    e = y - Q @ (Q.T @ y)
    s2 = e @ e / (n - p)
    s2_del = (e @ e - e ** 2 / (1 - h)) / (n - p - 1)
    t = e / np.sqrt(s2_del * (1 - h))
    cook = e ** 2 * h / (p * s2 * (1 - h) ** 2)
    return h, e, s2, t, cook


# <<data>>
rng = np.random.default_rng(20)
n0, m = 30, 4                                      # clean cases, size of the cluster
x_clean = np.sort(rng.uniform(0, 10, n0))
y_clean = 2 + 0.5 * x_clean + rng.normal(0, 0.5, n0)
x0, y0 = 20.0, 2.0                                 # four identical records
x = np.r_[x_clean, np.full(m, x0)]
y = np.r_[y_clean, np.full(m, y0)]
X = np.column_stack([np.ones(len(x)), x])
n, p = X.shape
cluster = np.arange(n0, n)
# <</data>>

# <<single>>
h, e, s2, t, cook = single_case(X, y)
crit = stats.t.isf(0.05 / (2 * n), n - p - 1)      # Bonferroni outlier test at 5%
print("slope with the cluster:", np.linalg.lstsq(X, y, rcond=None)[0][1])
print("cluster: h =", h[-1].round(3), " t =", t[-1].round(2), " Cook's D =", cook[-1].round(3))
i = int(np.argmax(np.abs(t)))
print(f"largest |t|: case {i + 1} at x = {x[i]:.2f}, t = {t[i]:.2f}; Bonferroni cut-off {crit:.2f}")
# <</single>>

# <<group>>
D = cluster                                        # delete the whole cluster at once
K = np.linalg.inv(X.T @ X)
H_D = X[D] @ K @ X[D].T
pred_D = np.linalg.solve(np.eye(m) - H_D, e[D])    # y_D minus its prediction from the rest
sse = e @ e
sse_D = sse - e[D] @ pred_D
F = (e[D] @ pred_D / m) / (sse_D / (n - p - m))
print("group prediction residuals:", pred_D.round(3))
print(f"F = {F:.1f} on {m} and {n - p - m} degrees of freedom,"
      f" p-value {stats.f.sf(F, m, n - p - m):.1e}")
# <</group>>

# ---- checks of the group-deletion theorem on random subsets -----------------------------
beta = K @ X.T @ y
for trial in range(50):
    d = rng.integers(1, 6)
    Dset = rng.choice(n, size=d, replace=False)
    keep = np.setdiff1d(np.arange(n), Dset)
    b_D, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
    HD = X[Dset] @ K @ X[Dset].T
    if np.linalg.matrix_rank(X[keep]) < p:
        continue
    corr = np.linalg.solve(np.eye(d) - HD, e[Dset])
    assert np.allclose(beta - b_D, K @ X[Dset].T @ corr)
    assert np.allclose(y[Dset] - X[Dset] @ b_D, corr)
    assert np.isclose(np.sum((y[keep] - X[keep] @ b_D) ** 2), sse - e[Dset] @ corr)
    # mean-shift model with the indicators of Dset gives the same F statistic
    Zs = np.column_stack([X, np.eye(n)[:, Dset]])
    bz, *_ = np.linalg.lstsq(Zs, y, rcond=None)
    sse_z = np.sum((y - Zs @ bz) ** 2)
    assert np.isclose(sse_z, sse - e[Dset] @ corr)

# ---- checks of the masking proposition ------------------------------------------------
A = X[:n0].T @ X[:n0]
xv = np.array([1.0, x0])
h1 = xv @ np.linalg.solve(A, xv)                   # leverage of x0 relative to the clean data
b_clean = np.linalg.solve(A, X[:n0].T @ y_clean)
dev = y0 - xv @ b_clean                            # how far the cluster lies from the clean fit
assert np.allclose(h[cluster], h1 / (1 + m * h1)) and h[-1] < 1 / m
assert np.allclose(e[cluster], dev / (1 + m * h1))
assert np.allclose(e[cluster] / (1 - h[cluster]), dev / (1 + (m - 1) * h1))
assert np.allclose(pred_D, dev)
d_single = K @ xv * e[-1] / (1 - h[-1])
d_group = beta - b_clean
assert np.allclose(d_group, m * (1 + (m - 1) * h1) * d_single)
cross = X[:n0] @ np.linalg.solve(A, xv)            # cross-leverages in the clean fit
assert np.allclose(e[:n0], (y_clean - X[:n0] @ b_clean) - m * dev * cross / (1 + m * h1))
assert np.isclose(sse, np.sum((y_clean - X[:n0] @ b_clean) ** 2) + m * dev ** 2 / (1 + m * h1))

# swamping: clean cases flagged, cluster not
assert np.all(np.abs(t[cluster]) < 2) and i < n0 and np.abs(t[i]) > 2
assert np.all(np.abs(t) < crit)
assert F > 50
# a single copy at x0 would have leverage h1 / (1 + h1)
h_single = h1 / (1 + h1)
_, e_one, _, t_one, cook_one = single_case(X[: n0 + 1], y[: n0 + 1])
assert np.abs(t_one[-1]) > crit

# clean-data residual scale and the group-deleted fit
s_clean = np.sqrt(np.sum((y_clean - X[:n0] @ b_clean) ** 2) / (n0 - p))

gen = Generated("ch20", "masking", prefix="mk")
gen.int("n0", n0)
gen.int("m", m)
gen.int("n", n)
gen.num("x0", x0, 0)
gen.num("h1", h1, 3)
gen.num("h_cl", h[-1], 3)
gen.num("h_single", h_single, 3)
gen.num("dev", dev, 3)
gen.num("e_cl", e[-1], 3)
gen.num("pred_single", e[-1] / (1 - h[-1]), 3)
gen.num("t_cl", t[-1], 2)
gen.num("D_cl", cook[-1], 3)
gen.num("cut_D", 4 / n, 3)
gen.num("t_one", t_one[-1], 2)
gen.num("D_one", cook_one[-1], 2)
gen.num("ratio", m * (1 + (m - 1) * h1), 2)
gen.num("slope_all", beta[1], 3)
gen.num("slope_clean", b_clean[1], 3)
gen.num("s_all", np.sqrt(s2), 3)
gen.num("s_clean", s_clean, 3)
gen.int("i_max", i + 1)
gen.num("x_imax", x[i], 2)
gen.num("t_imax", t[i], 2)
gen.int("n_t2", int(np.sum(np.abs(t) > 2)))
gen.int("n_cook_clean", int(np.sum(cook[:n0] > 4 / n)))
gen.num("cook_clean_max", cook[:n0].max(), 3)
gen.num("crit", crit, 2)
gen.num("F", F, 1)
gen.num("F_p", stats.f.sf(F, m, n - p - m), 1, sci=True)
gen.write()

# ---- figure --------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
ax.scatter(x_clean, y_clean, s=9, color=COLORS["accent"], linewidths=0, label="clean cases")
ax.scatter([x0], [y0], s=40, marker="D", color=COLORS["second"], linewidths=0,
           label=f"{m} identical cases")
g = np.linspace(0, 21, 2)
ax.plot(g, beta[0] + beta[1] * g, color=COLORS["second"], linestyle="--", label="fit, all cases")
ax.plot(g, b_clean[0] + b_clean[1] * g, color=COLORS["accent"], label="fit, cluster deleted")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("(a) data and two fits")
ax = axes[1]
idx = np.arange(1, n + 1)
ax.axhline(crit, color=COLORS["muted"], linestyle=":", linewidth=0.8)
ax.axhline(-crit, color=COLORS["muted"], linestyle=":", linewidth=0.8)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6)
ax.vlines(idx[:n0], 0, t[:n0], color=COLORS["accent"], linewidth=1.2)
ax.vlines(idx[n0:], 0, t[n0:], color=COLORS["second"], linewidth=2.0)
ax.set_xlabel("case (clean cases in order of $x$)")
ax.set_ylabel(r"$t_i$")
ax.set_ylim(-4, 4)
ax.set_title("(b) externally studentized residuals")
fig.tight_layout()
fig.savefig(figure_path("ch20", "masking"))
