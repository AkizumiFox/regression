"""Chapter 30, Section 5: forward stagewise regression and the lasso path.

(1) Orthonormal design: every coefficient of incremental forward stagewise (step eps) stays within
    eps of the lasso solution at lambda = the current largest absolute residual correlation.
(2) A correlated design (n = 40, six regressors sharing a common factor) in which one lasso
    coefficient shrinks back towards zero along the path. Forward stagewise with a small step
    tracks the lasso path until that happens, then holds the coefficient (the monotone lasso).
    At every L1 norm the lasso has the smaller residual sum of squares, as its definition demands.
(3) Test error along the lasso, stagewise and componentwise boosting paths.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)


def lasso_cd(X, y, lam, b=None, tol=1e-12):
    G, c = X.T @ X, X.T @ y
    b = np.zeros(X.shape[1]) if b is None else b.copy()
    while True:
        largest = 0.0
        for j in range(len(b)):
            old = b[j]
            b[j] = soft(c[j] - G[j] @ b + G[j, j] * old, lam) / G[j, j]
            largest = max(largest, abs(b[j] - old))
        if largest < tol:
            return b


# <<stagewise>>
def forward_stagewise(X, y, eps, steps):
    """Incremental forward stagewise: move the most correlated coefficient by eps."""
    b = np.zeros(X.shape[1])
    r = y.astype(float).copy()
    path = [b.copy()]
    for _ in range(steps):
        c = X.T @ r
        j = np.argmax(np.abs(c))
        delta = eps * np.sign(c[j])
        b[j] += delta
        r -= delta * X[:, j]
        path.append(b.copy())
    return np.array(path)
# <</stagewise>>


# ---- (1) orthonormal design ----------------------------------------------------------------------
rng = np.random.default_rng(3005)
Q, _ = np.linalg.qr(rng.normal(size=(30, 5)))
yq = Q @ np.array([3.0, -2.0, 1.2, 0.4, 0.0]) + 0.3 * rng.normal(size=30)
zq = Q.T @ yq
eps_q = 0.01
Pq = forward_stagewise(Q, yq, eps_q, 800)
worst_q = 0.0
for b in Pq:
    tau = np.max(np.abs(zq - b))
    if tau < eps_q:
        break
    gap = np.max(np.abs(b - soft(zq, tau)))
    assert gap <= eps_q + 1e-12
    worst_q = max(worst_q, gap)

# ---- (2) a correlated design ---------------------------------------------------------------------
# <<data>>
rng = np.random.default_rng(4)
n, p = 40, 6
z = rng.normal(size=(n, 1))
X_raw = 0.8 * z + 0.6 * rng.normal(size=(n, p))        # a shared factor: correlations near 0.64
center = X_raw.mean(axis=0)
scale = np.linalg.norm(X_raw - center, axis=0)
X = (X_raw - center) / scale                           # centred, unit-length columns
beta = np.array([3.0, -2.0, 2.0, 0.0, 0.0, -1.0])
y = X @ beta + 0.5 * rng.normal(size=n)
y_bar = y.mean()
y = y - y_bar
# <</data>>

lam_max = np.max(np.abs(X.T @ y))
grid = lam_max * np.logspace(0, -3, 400)
lasso_path = []
b = np.zeros(p)
for lam in grid:
    b = lasso_cd(X, y, lam, b)
    lasso_path.append(b.copy())
lasso_path = np.array(lasso_path)
l1_lasso = np.abs(lasso_path).sum(axis=1)
shrinks = np.diff(np.abs(lasso_path), axis=0) < -1e-8
assert shrinks.any()                                    # the lasso path is not monotone here
j_back = int(np.argmax(shrinks.any(axis=0)))            # the coefficient that shrinks back
first_back = int(np.argmax(shrinks.any(axis=1)))
t_back = l1_lasso[first_back]
peak_back = lasso_path[first_back, j_back]

# <<compare>>
eps = 0.002
fs_path = forward_stagewise(X, y, eps, 6000)
l1_fs = np.abs(fs_path).sum(axis=1)
arc_fs = eps * np.arange(len(fs_path))                  # total variation of the stagewise path
print("lasso at the end:     ", np.round(lasso_path[-1], 3))
print("stagewise at the end: ", np.round(fs_path[-1], 3))
# <</compare>>

# stagewise agrees with the lasso until the lasso coefficient turns back
gaps = []
for k in range(0, len(fs_path), 25):
    t = l1_fs[k]
    if t >= t_back - 0.05:
        break
    i = np.searchsorted(l1_lasso, t)
    gaps.append(np.max(np.abs(fs_path[k] - lasso_path[i])))
agree_gap = max(gaps)
assert agree_gap < 0.05
# afterwards stagewise holds the coefficient while the lasso shrinks it
k_mid = np.searchsorted(l1_fs, 0.5 * (t_back + l1_lasso[-1]))
i_mid = np.searchsorted(l1_lasso, l1_fs[k_mid])
hold_gap = abs(fs_path[k_mid, j_back] - lasso_path[i_mid, j_back])
assert hold_gap > 0.1
# until stagewise reaches the least squares region its path is monotone: arc length = L1 norm
k_mono = np.searchsorted(l1_fs, 0.9 * l1_lasso[-1])
assert np.isclose(arc_fs[k_mono], l1_fs[k_mono], atol=1e-9)


def lasso_rss_at(t):
    """Residual sum of squares of the lasso with ||b||_1 = t, by bisection on lambda."""
    lo, hi = 0.0, lam_max
    bb = np.zeros(p)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        bb = lasso_cd(X, y, mid, bb)
        lo, hi = (mid, hi) if np.abs(bb).sum() > t else (lo, mid)
    return np.sum((y - X @ bb) ** 2)


rss_excess = []
for k in np.linspace(100, np.searchsorted(l1_fs, l1_lasso[-1] * 0.98), 12).astype(int):
    rss_fs = np.sum((y - X @ fs_path[k]) ** 2)
    rss_l = lasso_rss_at(l1_fs[k])
    assert rss_fs >= rss_l - 1e-8
    rss_excess.append(rss_fs - rss_l)

# ---- (3) test error along the three paths -------------------------------------------------------
def boost(X, y, nu, steps):
    norms2 = np.sum(X ** 2, axis=0)
    b, r, path = np.zeros(X.shape[1]), y.copy(), []
    for _ in range(steps):
        c = X.T @ r
        j = np.argmax(c ** 2 / norms2)
        s = nu * c[j] / norms2[j]
        b[j] += s
        r -= s * X[:, j]
        path.append(b.copy())
    return np.array(path)


# the design of Section 30.4: n = 60, p = 30, correlation 0.5^|j-k|, five nonzero coefficients, sigma 1.5
rng3 = np.random.default_rng(3006)
n3, p3, sigma3, reps = 60, 30, 1.5, 50
C3 = 0.5 ** np.abs(np.subtract.outer(np.arange(p3), np.arange(p3)))
L3 = np.linalg.cholesky(C3)
beta3 = np.zeros(p3)
beta3[[0, 4, 9, 14, 19]] = [2.0, -1.5, 1.0, 1.0, -0.5]
table = []
for rep in range(reps):
    Xr = rng3.normal(size=(n3, p3)) @ L3.T
    yr = Xr @ beta3 + sigma3 * rng3.normal(size=n3)
    c3, s3, m3 = Xr.mean(axis=0), np.linalg.norm(Xr - Xr.mean(axis=0), axis=0), yr.mean()
    Xs3, ys3 = (Xr - c3) / s3, yr - m3
    def risk(path):                                    # exact risk on a new case from the design law
        d = beta3[None, :] - path / s3
        return np.sum((d @ C3) * d, axis=1) + (path @ (c3 / s3) - m3) ** 2 + sigma3 ** 2
    lm3 = np.max(np.abs(Xs3.T @ ys3))
    lp, bb = [], np.zeros(p3)
    for lam in lm3 * np.logspace(0, -3, 60):
        bb = lasso_cd(Xs3, ys3, lam, bb, tol=1e-9)
        lp.append(bb.copy())
    b_ls3 = np.linalg.lstsq(Xs3, ys3, rcond=None)[0]
    table.append([risk(np.array(lp)).min(), risk(forward_stagewise(Xs3, ys3, 0.05, 1500)).min(),
                  risk(boost(Xs3, ys3, 0.05, 1500)).min(), risk(b_ls3[None, :])[0]])
table = np.array(table)
mean_risk = table.mean(axis=0)
print("mean best risk: lasso, stagewise, boosting; least squares:", np.round(mean_risk, 3))
assert np.all(mean_risk[:3] < mean_risk[3])
assert np.ptp(mean_risk[:3]) < 0.05 * mean_risk[:3].min()
# data set by data set: the largest relative difference from the lasso's best risk
rel_diff = np.abs(table[:, 1:3] - table[:, [0]]) / table[:, [0]]
max_rel_diff = rel_diff.max()
print("largest relative difference from the lasso, over data sets:", round(max_rel_diff, 4))
assert max_rel_diff < 0.1

gen = Generated("ch30", "stagewise")
gen.num("eps_q", eps_q, 2)
gen.num("worst_q", worst_q, 4)
gen.int("n", n)
gen.int("p", p)
gen.num("eps", eps, 3)
gen.int("j_back", j_back + 1)
gen.num("t_back", t_back, 2)
gen.num("peak_back", peak_back, 3)
gen.num("agree_gap", agree_gap, 3)
gen.num("hold_gap", hold_gap, 3)
gen.num("rss_excess_max", max(rss_excess), 3)
gen.num("l1_ls", l1_lasso[-1], 2)
gen.int("reps", reps)
gen.num("risk_lasso", mean_risk[0], 3)
gen.num("risk_fs", mean_risk[1], 3)
gen.num("risk_boost", mean_risk[2], 3)
gen.num("risk_ls", mean_risk[3], 3)
gen.num("max_rel_diff_pct", 100 * max_rel_diff, 1)
gen.num("sd_diff_fs", np.std(table[:, 1] - table[:, 0]) / np.sqrt(reps), 3)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.2, 2.2))
cols = [COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["muted"], COLORS["ink"], COLORS["thread"]]
for ax, (l1, path, title) in zip(axes[:2], [(l1_lasso, lasso_path, "(a) lasso"),
                                            (l1_fs, fs_path, r"(b) forward stagewise")]):
    for j in range(p):
        ax.plot(l1, path[:, j], color=cols[j], linewidth=1.4 if j == j_back else 0.9)
    ax.axvline(t_back, color=COLORS["grid"], linewidth=0.7, linestyle="--", zorder=0)
    ax.set_xlabel(r"$\|\mathbf{b}\|_1$")
    ax.set_title(title)
axes[0].set_ylabel("coefficient")
axes[1].sharey(axes[0])
ax = axes[2]
ax.scatter(table[:, 0], table[:, 1], s=8, color=COLORS["second"], label="stagewise", linewidths=0)
ax.scatter(table[:, 0], table[:, 2], s=8, color=COLORS["third"], marker="^", label=r"boosting", linewidths=0)
lims = [table[:, :3].min() * 0.97, table[:, :3].max() * 1.03]
ax.plot(lims, lims, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel("best risk, lasso")
ax.set_title("(c) best risk, 50 data sets")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch30", "stagewise_lasso"))
