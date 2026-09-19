"""Chapter 28, Section 1: minimum-norm least squares when p > n.

Part 1: a small example (n = 4 observations, p = 6 regressors) comparing the minimum-norm
interpolator with another interpolator, the ridgeless limit and gradient descent from zero.
Part 2: the risk of the minimum-norm fit for isotropic Gaussian regressors when only the first
p of D features are used (the "double descent" curve): the exact formula of the text against
simulation, and the risk of ridge regression with the best penalty.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch28", "min_norm")

# <<small>>
import numpy as np

rng = np.random.default_rng(2801)
n, p = 4, 6
X = np.round(rng.normal(size=(n, p)), 1)            # a 4 x 6 design, rank 4
y = np.array([1.0, 2.0, 0.0, -1.0])

b_plus = X.T @ np.linalg.solve(X @ X.T, y)          # X^T (X X^T)^{-1} y
print("min-norm solution ", np.round(b_plus, 3))
print("residual norm     ", np.linalg.norm(y - X @ b_plus))

# another interpolator: add a null-space vector of X
_, _, Vt = np.linalg.svd(X)
b_other = b_plus + 2.0 * Vt[-1]                     # Vt[-1] spans part of N(X)
print("other solution    ", np.round(b_other, 3))
print("residual norm     ", np.linalg.norm(y - X @ b_other))
print("norms", np.linalg.norm(b_plus), np.linalg.norm(b_other))

x0 = np.ones(p)                                      # a new point
print("predictions at x0:", x0 @ b_plus, x0 @ b_other)
# <</small>>

assert np.linalg.matrix_rank(X) == n
assert np.allclose(b_plus, np.linalg.pinv(X) @ y)
assert np.allclose(X @ b_plus, y) and np.allclose(X @ b_other, y)
assert np.linalg.norm(b_plus) < np.linalg.norm(b_other)
assert abs(Vt[-1] @ b_plus) < 1e-10                   # b_plus is orthogonal to N(X)
gen.num("small_norm_plus", np.linalg.norm(b_plus), 3)
gen.num("small_norm_other", np.linalg.norm(b_other), 3)
gen.num("small_pred_plus", x0 @ b_plus, 3)
gen.num("small_pred_other", x0 @ b_other, 3)
gen.num("small_b1_plus", b_plus[0], 3)
gen.num("small_b1_other", b_other[0], 3)

# <<ridgeless>>
for lam in [1.0, 0.1, 0.01, 0.001]:
    b_ridge = np.linalg.solve(X.T @ X + lam * np.eye(p), X.T @ y)
    print(f"lambda = {lam:6.3f}   distance to min-norm {np.linalg.norm(b_ridge - b_plus):.2e}")

# gradient descent on ||y - X b||^2 / 2 started at zero
step = 1.0 / np.linalg.norm(X, 2) ** 2
b = np.zeros(p)
for t in range(5000):
    b = b + step * X.T @ (y - X @ b)
print("gradient descent: distance to min-norm", np.linalg.norm(b - b_plus))
# <</ridgeless>>

dists = []
for lam in [1.0, 0.1, 0.01, 0.001]:
    b_ridge = np.linalg.solve(X.T @ X + lam * np.eye(p), X.T @ y)
    dists.append(np.linalg.norm(b_ridge - b_plus))
assert all(d2 < d1 for d1, d2 in zip(dists, dists[1:]))
assert np.linalg.norm(b - b_plus) < 1e-8
# started elsewhere, gradient descent finds a different interpolator: b_start's null-space part survives
b2 = Vt[-1].copy()
for t in range(5000):
    b2 = b2 + step * X.T @ (y - X @ b2)
assert np.allclose(b2, b_plus + Vt[-1], atol=1e-8)
gen.num("ridge_dist_001", dists[2], 4)


# ---- Part 2: double descent -------------------------------------------------------------
# <<risk>>
def exact_risk(p, n, beta, sigma2):
    """Excess prediction risk of min-norm least squares on the first p of D isotropic features."""
    kept, left = np.sum(beta[:p] ** 2), np.sum(beta[p:] ** 2)
    noise = sigma2 + left                           # omitted features act as extra noise
    if p == 0:
        return left
    if p <= n - 2:
        return left + noise * p / (n - p - 1)
    if p >= n + 2:
        return left + kept * (1 - n / p) + noise * n / (p - n - 1)
    return np.inf                                   # |p - n| <= 1: not covered by the formula


def simulated_risk(p, n, beta, sigma2, reps, rng, lams=()):
    """Average over designs of the conditional risk given X; also ridge for each penalty."""
    bS, left = beta[:p], np.sum(beta[p:] ** 2)
    noise = sigma2 + left
    out, ridge = [], np.zeros(len(lams))
    for _ in range(reps):
        X = rng.normal(size=(n, p))
        _, d, Vt = np.linalg.svd(X, full_matrices=False)
        c = Vt @ bS                                 # coordinates of beta_S in the row space
        bias2 = bS @ bS - c @ c                     # ||(I - P) beta_S||^2
        out.append(left + bias2 + noise * np.sum(1 / d ** 2))
        if len(lams):
            shrink = d ** 2 / (d ** 2 + np.asarray(lams)[:, None])     # one row per penalty
            ridge += left + bias2 + np.sum(((1 - shrink) * c) ** 2, axis=1) \
                + noise * np.sum(shrink ** 2 / d ** 2, axis=1)
    return np.mean(out), np.std(out) / np.sqrt(reps), ridge / reps
# <</risk>>


# <<risk-small>>
n, D, sigma2 = 20, 100, 0.25
beta = np.full(D, 1 / np.sqrt(D))                  # signal spread evenly, ||beta||^2 = 1
rng = np.random.default_rng(2802)
for p in [5, 10, 15, 25, 40, 100]:
    sim, se, _ = simulated_risk(p, n, beta, sigma2, reps=200, rng=rng)
    print(f"p = {p:3d}: formula {exact_risk(p, n, beta, sigma2):.3f}   simulation {sim:.3f}")
# <</risk-small>>

n, D, sigma2 = 50, 200, 0.25
reps = 400
lams = np.geomspace(1e-2, 1e4, 61)
grid = [p for p in list(range(1, 46)) + list(range(55, 101)) + list(range(110, D + 1, 10))]
settings = {
    "spread": np.full(D, 1 / np.sqrt(D)),
    "decaying": 0.8 ** np.arange(D),
}
settings["decaying"] /= np.linalg.norm(settings["decaying"])
results = {}
rng = np.random.default_rng(2803)
for name, beta in settings.items():
    ps = np.arange(0, D + 1)
    formula = np.array([exact_risk(p, n, beta, sigma2) for p in ps])
    sim_p = [p for p in grid if p % 5 == 0 and abs(p - n) >= 5]
    sims, ses, ridge_best = [], [], []
    for p in sim_p:
        m, se, rr = simulated_risk(p, n, beta, sigma2, reps, rng, lams)
        sims.append(m)
        ses.append(se)
        ridge_best.append(rr.min())
        # the simulation agrees with the exact formula
        assert abs(m - exact_risk(p, n, beta, sigma2)) < 5 * se + 1e-3, (name, p, m)
        # the best ridge is never worse than min-norm least squares
        assert rr.min() <= m + 1e-9
    results[name] = dict(ps=ps, formula=formula, sim_p=np.array(sim_p), sims=np.array(sims),
                         ridge=np.array(ridge_best))

# ---- numbers quoted in the text ---------------------------------------------------------------
sp, dc = results["spread"], results["decaying"]
f_sp, f_dc = sp["formula"], dc["formula"]
under = np.arange(0, n - 1)                        # p = 0, ..., n - 2
over = np.arange(n + 2, D + 1)
best_under_sp = under[np.argmin(f_sp[under])]
best_over_sp = over[np.argmin(f_sp[over])]
best_under_dc = under[np.argmin(f_dc[under])]
gen.int("n", n)
gen.int("D", D)
gen.num("sigma2", sigma2, 2)
gen.num("sp_null", f_sp[0], 3)
gen.int("sp_best_under_p", best_under_sp)
gen.num("sp_best_under", f_sp[best_under_sp], 3)
gen.int("sp_best_over_p", best_over_sp)
gen.num("sp_best_over", f_sp[best_over_sp], 3)
gen.num("sp_p45", f_sp[45], 2)
gen.num("sp_p55", f_sp[55], 2)
gen.int("dc_best_under_p", best_under_dc)
gen.num("dc_best_under", f_dc[best_under_dc], 3)
gen.num("dc_over_D", f_dc[D], 3)
gen.num("sp_ridge_min", sp["ridge"].min(), 3)
gen.num("dc_ridge_min", dc["ridge"].min(), 3)
# the qualitative claims of the text
assert best_under_sp == 0 and best_over_sp == D      # spread signal: interpolating with everything wins
assert f_sp[D] < f_sp[0]
assert 0 < best_under_dc < n - 1 and f_dc[best_under_dc] < f_dc[over].min()   # decaying: a small model wins
assert np.all(np.diff(f_sp[over]) < 0)                # descent after the peak
for r in results.values():                          # tuned ridge has no peak near p = n
    near = (np.abs(r["sim_p"] - n) <= 10)
    assert r["ridge"][near].max() < 1.2 * max(r["ridge"].min(), r["formula"][0])
gen.write()

# ---- figure -----------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), sharey=True)
for ax, (name, r), title in zip(axes, results.items(),
                                ["(a) signal spread over all 200 features",
                                 "(b) signal in the first few features"]):
    ps, f = r["ps"], r["formula"].copy()
    f[~np.isfinite(f)] = np.nan
    ax.plot(ps, f, color=COLORS["accent"], label="min-norm, exact")
    ax.plot(r["sim_p"], r["sims"], "o", ms=2.2, color=COLORS["accent"], alpha=0.8, label="min-norm, simulated")
    ax.plot(r["sim_p"], r["ridge"], color=COLORS["second"], lw=1.1, label="ridge, best penalty")
    ax.axvline(n, color=COLORS["grid"], lw=0.8, zorder=0)
    ax.axhline(r["formula"][0], color=COLORS["muted"], lw=0.7, ls=":", zorder=0)
    ax.set_xlim(0, D)
    ax.set_ylim(0, 3)
    ax.set_xlabel("number of features used, $p$")
    ax.set_title(title)
axes[0].set_ylabel("excess prediction risk")
axes[1].legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch28", "double_descent"))
