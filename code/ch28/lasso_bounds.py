"""Chapter 28, Section 3: prediction and estimation bounds for the lasso.

A fixed Gaussian design with n = 200, p = 1000, independent N(0, 1) entries and columns
scaled to ||x_j||^2 = n (the design of the compatibility example in Section 2); s = 5 nonzero coefficients; normal errors with sigma = 1.
(1) The probability of the event max_j |x_j^T eps| / n <= lambda_0 with
    lambda_0 = sigma sqrt(2 log(2p/delta)/n), delta = 0.05.
(2) The prediction and estimation errors along the lasso path, against the bounds of the text,
    with the compatibility constant of the design computed (certified by a duality gap).
(3) How the error at lambda = 2 lambda_0 grows with p.
The lasso is computed by cyclic coordinate descent, vectorized over replications.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch28", "lasso_bounds")


# <<lasso>>
def soft(z, t):
    """Soft thresholding S(z, t) = sign(z) max(|z| - t, 0)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, Y, lam, B=None, tol=1e-9, max_sweeps=2000):
    """Minimize ||y - X b||^2 / (2n) + lam ||b||_1 for every column y of Y, by coordinate descent."""
    n, p = X.shape
    Y = Y.reshape(n, -1)
    B = np.zeros((p, Y.shape[1])) if B is None else B.copy()
    R = Y - X @ B                                    # residuals, one column per response
    scale = np.sum(X ** 2, axis=0) / n
    active = np.arange(p)
    for sweep in range(max_sweeps):
        change = 0.0
        for j in active:
            z = X[:, j] @ R / n + scale[j] * B[j]
            new = soft(z, lam) / scale[j]
            d = new - B[j]
            if np.any(d != 0):
                R -= np.outer(X[:, j], d)
                B[j] = new
                change = max(change, np.max(np.abs(d)))
        if change < tol:
            if len(active) == p:                     # converged over all coordinates
                break
            active = np.arange(p)                    # check every coordinate once more
        else:
            active = np.flatnonzero(np.any(B != 0, axis=1)) if sweep % 5 else np.arange(p)
    return B


def kkt_violation(X, Y, B, lam):
    """Largest violation of the lasso optimality conditions (zero at a solution)."""
    n = X.shape[0]
    G = X.T @ (Y.reshape(n, -1) - X @ B) / n         # gradient part, should lie in lam * subdifferential
    off = np.where(B == 0, np.maximum(np.abs(G) - lam, 0), np.abs(G - lam * np.sign(B)))
    return off.max()
# <</lasso>>


def compatibility(X, S, radius=3.0, iters=4000):
    """phi^2(S) = min s ||Xv||^2/n / ||v_S||_1^2 over the cone ||v_{S^c}||_1 <= 3 ||v_S||_1.

    For each sign pattern of v_S (up to a global sign) this is a convex problem over
    {v_S in the signed simplex, ||v_{S^c}||_1 <= 3}; it is solved by accelerated projected
    gradient, and the Frank-Wolfe duality gap gives a certified lower bound.
    """
    n, p = X.shape
    s = len(S)
    Sc = np.setdiff1d(np.arange(p), S)
    A, C = X[:, S], X[:, Sc]
    L = 2 * s / n * np.linalg.norm(X, 2) ** 2

    def proj_simplex(u):
        v = np.sort(u)[::-1]
        css = np.cumsum(v) - 1
        k = np.flatnonzero(v - css / np.arange(1, len(u) + 1) > 0)[-1]
        return np.maximum(u - css[k] / (k + 1), 0)

    def proj_l1(w):
        if np.abs(w).sum() <= radius:
            return w
        return np.sign(w) * proj_simplex(np.abs(w) / radius) * radius

    best_val, best_lower = np.inf, np.inf
    for code in range(2 ** (s - 1)):
        sg = np.array([1.0] + [1.0 if (code >> k) & 1 == 0 else -1.0 for k in range(s - 1)])
        u, w = np.full(s, 1 / s), np.zeros(len(Sc))
        yu, yw, t = u.copy(), w.copy(), 1.0
        for _ in range(iters):
            r = A @ (sg * yu) + C @ yw
            gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
            u_new, w_new = proj_simplex(yu - gu / L), proj_l1(yw - gw / L)
            t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
            yu = u_new + (t - 1) / t_new * (u_new - u)
            yw = w_new + (t - 1) / t_new * (w_new - w)
            u, w, t = u_new, w_new, t_new
        r = A @ (sg * u) + C @ w
        val = s / n * r @ r
        gu, gw = 2 * s / n * sg * (A.T @ r), 2 * s / n * (C.T @ r)
        gap = gu @ u - gu.min() + gw @ w + radius * np.abs(gw).max()
        best_val, best_lower = min(best_val, val), min(best_lower, val - gap)
    return best_val, best_lower


# <<small>>
rng = np.random.default_rng(2832)
n, p, s, sigma, delta = 100, 300, 5, 1.0, 0.05
X = rng.normal(size=(n, p))
X /= np.sqrt(np.sum(X ** 2, axis=0) / n)           # columns scaled to ||x_j||^2 = n
beta = np.zeros(p)
beta[:s] = [1.0, -1.0, 1.0, -1.0, 1.0]
lam0 = sigma * np.sqrt(2 * np.log(2 * p / delta) / n)
E = sigma * rng.normal(size=(n, 50))
Y = (X @ beta)[:, None] + E
print("lambda_0 =", round(lam0, 3))
print("frequency of the event:", np.mean(np.max(np.abs(X.T @ E) / n, axis=0) <= lam0))
for ratio in [2.0, 1.0, 0.5]:
    B = lasso_cd(X, Y, ratio * lam0)
    err = np.sum((X @ (B - beta[:, None])) ** 2, axis=0) / n
    print(f"lambda = {ratio} lambda_0: mean prediction error {err.mean():.3f},"
          f" mean number selected {np.mean(np.sum(B != 0, axis=0)):.1f}")
print("oracle sigma^2 s / n =", sigma ** 2 * s / n, "  slow-rate bound at 2 lambda_0:", 6 * lam0 * s)
# <</small>>

# <<setup>>
rng = np.random.default_rng(2830)
n, p, s, sigma, delta = 200, 1000, 5, 1.0, 0.05
X = rng.normal(size=(n, p))
X /= np.sqrt(np.sum(X ** 2, axis=0) / n)           # columns scaled to ||x_j||^2 = n
beta = np.zeros(p)
S = np.arange(s)
beta[S] = [1.0, -1.0, 1.0, -1.0, 1.0]
lam0 = sigma * np.sqrt(2 * np.log(2 * p / delta) / n)
# <</setup>>

reps = 200
E = sigma * rng.normal(size=(n, reps))
Y = (X @ beta)[:, None] + E
event = np.max(np.abs(X.T @ E) / n, axis=0) <= lam0
gen.num("lam0", lam0, 4)
gen.num("event_freq", event.mean(), 3)
gen.int("reps", reps)
assert event.mean() >= 1 - delta - 3 * np.sqrt(delta * (1 - delta) / reps)

# compatibility constant of the design for S
phi2, phi2_lower = compatibility(X, S)
print(f"compatibility constant phi^2(S) = {phi2:.4f} (certified lower bound {phi2_lower:.4f})")
assert phi2 - phi2_lower < 0.02 * phi2
gen.num("phi2", phi2_lower, 5)

# the lasso path
ratios = np.geomspace(0.1, 6, 25)[::-1]
B = None
pred, l1err, l2err, kkt, supp = [], [], [], [], []
for r in ratios:
    B = lasso_cd(X, Y, r * lam0, B)
    kkt.append(kkt_violation(X, Y, B, r * lam0))
    D = B - beta[:, None]
    pred.append(np.mean(np.sum((X @ D) ** 2, axis=0) / n))
    l1err.append(np.mean(np.sum(np.abs(D), axis=0)))
    l2err.append(np.mean(np.sqrt(np.sum(D ** 2, axis=0))))
    supp.append(np.mean(np.sum(B != 0, axis=0)))
pred, l1err, l2err, supp = map(np.array, (pred, l1err, l2err, supp))
assert max(kkt) < 1e-6
B2 = lasso_cd(X, Y, 2 * lam0)
D2 = B2 - beta[:, None]
pred2 = np.sum((X @ D2) ** 2, axis=0) / n
l1_2 = np.sum(np.abs(D2), axis=0)
l2_2 = np.sqrt(np.sum(D2 ** 2, axis=0))
cone = np.sum(np.abs(np.delete(D2, S, axis=0)), axis=0) <= 3 * np.sum(np.abs(D2[S]), axis=0) + 1e-12
lam = 2 * lam0
slow = 3 * lam * np.abs(beta).sum()
fast = 9 * lam ** 2 * s / phi2_lower
l1_bound = 4 * lam * s / phi2_lower
# the theorems hold on the event (they are deterministic given it)
assert np.all(cone[event])
assert np.all(pred2[event] <= min(slow, fast))
assert np.all(l1_2[event] <= l1_bound)
oracle = sigma ** 2 * s / n
k_best = np.argmin(pred)
gen.num("pred_2lam0", pred2.mean(), 4)
gen.num("pred_best", pred[k_best], 4)
gen.num("ratio_best", ratios[k_best], 2)
gen.num("oracle", oracle, 3)
gen.num("slow_bound", slow, 2)
gen.num("fast_bound", fast, 1)
gen.num("l1_2lam0", l1_2.mean(), 3)
gen.num("l1_bound", l1_bound, 1)
gen.num("l2_2lam0", l2_2.mean(), 3)
gen.num("supp_2lam0", np.mean(np.sum(B2 != 0, axis=0)), 1)
gen.num("cone_freq", cone.mean(), 3)
gen.num("pred_small_lam", pred[-1], 3)
gen.num("lam", lam, 3)
gen.num("s_lam2", s * lam ** 2, 2)
gen.num("ratio_small", ratios[-1], 1)
print(f"lambda_0 = {lam0:.4f}, P(event) ~ {event.mean():.3f}")
print(f"prediction error at 2 lambda_0: {pred2.mean():.4f}; oracle {oracle:.3f};"
      f" slow bound {slow:.2f}; fast bound {fast:.1f}")

# growth with p, independent design, lambda = 2 lambda_0(p)
ps = [100, 200, 500, 1000, 2000, 5000]
reps_p = 200
growth = []
rng_p = np.random.default_rng(2831)
X_all = rng_p.normal(size=(n, max(ps)))
X_all /= np.sqrt(np.sum(X_all ** 2, axis=0) / n)
E_all = sigma * rng_p.normal(size=(n, reps_p))
for q in ps:                                         # nested designs: the first q columns
    Xq = X_all[:, :q]
    bq = np.zeros(q)
    bq[:s] = [1.0, -1.0, 1.0, -1.0, 1.0]
    Yq = (Xq @ bq)[:, None] + E_all
    lq = sigma * np.sqrt(2 * np.log(2 * q / delta) / n)
    Bq = lasso_cd(Xq, Yq, 2 * lq)
    err = np.sum((Xq @ (Bq - bq[:, None])) ** 2, axis=0) / n
    growth.append((err.mean(), err.std() / np.sqrt(reps_p), lq))
growth = np.array(growth)
# error grows roughly in proportion to lambda_0^2, i.e. to log(2p/delta)
ratio = growth[:, 0] / growth[:, 2] ** 2
assert ratio.max() / ratio.min() < 1.5
print(growth)
assert np.all(np.diff(growth[:, 0]) > 0)
gen.num("growth_p100", growth[0, 0], 3)
gen.num("growth_p5000", growth[-1, 0], 3)
log_ratio = np.log(2 * 5000 / delta) / np.log(2 * 100 / delta)
err_ratio = growth[-1, 0] / growth[0, 0]
assert abs(err_ratio / log_ratio - 1) < 0.02            # close, not equal
gen.num("log_ratio", log_ratio, 3)
gen.num("err_ratio", err_ratio, 3)
gen.write()

# ---- figure --------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.plot(ratios, pred, color=COLORS["accent"], label="prediction error")
ax.plot(ratios, supp / 100, color=COLORS["third"], lw=0.9, ls="--", label="selected / 100")
ax.axhline(oracle, color=COLORS["muted"], lw=0.8, ls=":", label=r"oracle $\sigma^2 s/n$")
ax.axvline(2, color=COLORS["second"], lw=0.8)
ax.text(1.9, 0.9 * ax.get_ylim()[1], r"$\lambda=2\lambda_0$", color=COLORS["second"], fontsize=8, ha="right")
ax.set_xscale("log")
ax.set_xlabel(r"$\lambda/\lambda_0$")
ax.set_title("(a) error along the lasso path")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.85))
ax = axes[1]
ax.errorbar(ps, growth[:, 0], yerr=2 * growth[:, 1], fmt="o-", ms=3, color=COLORS["accent"],
            label=r"lasso, $\lambda=2\lambda_0$")
ax.set_xscale("log")
ax.set_ylim(0, None)
ax.set_xlabel("number of regressors $p$ ($n=200$, $s=5$)")
ax.set_title("(b) growth with $p$")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch28", "lasso_bounds"))
