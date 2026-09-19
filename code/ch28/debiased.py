"""Chapter 28, Section 5: intervals after selection, and the debiased lasso.

A fixed Gaussian design with n = 200, p = 400 and Toeplitz correlation 0.5^|j-k|, columns
scaled to ||x_j||^2 = n. Four nonzero coefficients: 1 at x_1, 0.15 at x_6, -1 at x_11 and 0.15 at
x_16; normal errors with sigma = 1. Targets of inference: beta_1 (strong), beta_6 (weak) and
beta_2 (zero, but correlated 0.5 with x_1). For 1000 data sets we compute
  (i)   the naive post-selection interval: lasso at lambda = sigma sqrt(2 log p / n), then least
        squares on the selected columns and the usual t interval (reported only when the
        target is selected);
  (ii)  the debiased lasso interval with node-wise lasso rows of Theta and
        sigma_hat^2 = RSS / (n - |selected|);
  (iii) the oracle t interval from least squares on the true support plus the target.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch28", "debiased")

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


# <<debias>>
def nodewise_row(X, j, lam_j):
    """Row theta_j of Theta from the lasso of x_j on the other columns, and tau_j^2."""
    n, p = X.shape
    others = np.delete(np.arange(p), j)
    gamma = lasso_cd(X[:, others], X[:, j], lam_j)[:, 0]
    z = X[:, j] - X[:, others] @ gamma               # lasso residual of x_j on the rest
    tau2 = X[:, j] @ z / n
    theta = np.zeros(p)
    theta[j] = 1.0
    theta[others] = -gamma
    return theta / tau2, tau2


def debiased_interval(X, Y, B, j, theta, level=0.95):
    """Debiased estimate of beta_j and its normal-theory interval, one per column of Y."""
    n = X.shape[0]
    R = Y - X @ B
    b = B[j] + theta @ X.T @ R / n                    # one-step correction
    sigma_hat = np.sqrt(np.sum(R ** 2, axis=0) / (n - np.sum(B != 0, axis=0)))
    G = X.T @ X / n
    se = sigma_hat * np.sqrt(theta @ G @ theta / n)
    z = stats.norm.ppf((1 + level) / 2)
    return b, b - z * se, b + z * se
# <</debias>>


def toeplitz_design(rng, n, p, rho):
    idx = np.arange(p)
    L = np.linalg.cholesky(rho ** np.abs(idx[:, None] - idx[None, :]))
    X = rng.normal(size=(n, p)) @ L.T
    return X / np.sqrt(np.sum(X ** 2, axis=0) / n)


def naive_interval(X, y, selected, j, level=0.95):
    """Least squares on the selected columns and the usual t interval for beta_j (j selected)."""
    A = X[:, selected]
    n, k = A.shape
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    s2 = np.sum((y - A @ coef) ** 2) / (n - k)
    pos = list(selected).index(j)
    se = np.sqrt(s2 * np.linalg.inv(A.T @ A)[pos, pos])
    t = stats.t.ppf((1 + level) / 2, n - k)
    return coef[pos] - t * se, coef[pos] + t * se


# <<setup>>
rng = np.random.default_rng(2850)
n, p, sigma = 200, 400, 1.0
X = toeplitz_design(rng, n, p, rho=0.5)
beta = np.zeros(p)
beta[[0, 5, 10, 15]] = [1.0, 0.15, -1.0, 0.15]
lam = sigma * np.sqrt(2 * np.log(p) / n)
targets = {"strong": 0, "weak": 5, "null": 1}
# <</setup>>

# <<small>>
Ys = (X @ beta)[:, None] + sigma * rng.normal(size=(n, 200))
Bs = lasso_cd(X, Ys, lam)
for name, j in targets.items():
    theta, tau2 = nodewise_row(X, j, 0.1 * np.sqrt(2 * np.log(p) / n))
    b, lo, hi = debiased_interval(X, Ys, Bs, j, theta)
    print(f"{name:6s}: lasso bias {Bs[j].mean() - beta[j]:+.3f}, debiased bias {b.mean() - beta[j]:+.3f},"
          f" coverage {np.mean((lo <= beta[j]) & (beta[j] <= hi)):.3f}")
# <</small>>

reps = 1000
Y = (X @ beta)[:, None] + sigma * rng.normal(size=(n, reps))
B = lasso_cd(X, Y, lam)
assert kkt_violation(X, Y, B, lam) < 1e-6

lam_node = 0.1 * np.sqrt(2 * np.log(p) / n)
cover = {}
for name, j in targets.items():
    theta, tau2 = nodewise_row(X, j, lam_node)
    G = X.T @ X / n
    e_j = np.eye(p)[j]
    mu_j = np.max(np.abs(G @ theta - e_j))
    assert mu_j <= lam_node / tau2 + 1e-8              # the bound from the node-wise KKT conditions
    assert theta @ G @ theta >= (1 - mu_j) ** 2 - 1e-10   # the variance lower bound of the text
    b, lo, hi = debiased_interval(X, Y, B, j, theta)
    deb = np.mean((lo <= beta[j]) & (beta[j] <= hi))
    # naive post-selection
    sel_hits, sel_count = 0, 0
    for r in range(reps):
        selected = np.flatnonzero(B[:, r])
        if j in selected:
            sel_count += 1
            l, h = naive_interval(X, Y[:, r], selected, j)
            sel_hits += l <= beta[j] <= h
    # oracle: least squares on the true support plus the target
    support = sorted(set(np.flatnonzero(beta)) | {j})
    orc_int = np.array([naive_interval(X, Y[:, r], support, j) for r in range(reps)])
    orc = np.mean((orc_int[:, 0] <= beta[j]) & (beta[j] <= orc_int[:, 1]))
    t_stat = (b - beta[j]) / ((hi - lo) / (2 * stats.norm.ppf(0.975)))
    cover[name] = dict(debiased=deb, naive=sel_hits / max(sel_count, 1), selected=sel_count / reps,
                       oracle=orc, t=t_stat, mu=mu_j, tau2=tau2,
                       lasso_bias=np.mean(B[j]) - beta[j], deb_bias=np.mean(b) - beta[j],
                       width_deb=np.mean(hi - lo), width_orc=np.mean(orc_int[:, 1] - orc_int[:, 0]))
    print(name, {k: np.round(v, 3) for k, v in cover[name].items() if k != "t"})

mc = 3 * np.sqrt(0.95 * 0.05 / reps)
for name in targets:
    assert abs(cover[name]["oracle"] - 0.95) < mc
    assert abs(cover[name]["debiased"] - 0.95) < 0.03
assert cover["weak"]["naive"] < 0.95 - mc
assert cover["null"]["naive"] < 0.95 - mc
assert cover["strong"]["lasso_bias"] < -0.05              # the lasso shrinks the strong coefficient
assert abs(cover["strong"]["deb_bias"]) < 0.02
t_strong = cover["strong"]["t"]
assert abs(t_strong.mean()) < 0.25 and abs(t_strong.std() - 1) < 0.1
gen.num("t_mean", t_strong.mean(), 2)
gen.num("t_sd", t_strong.std(), 2)
gen.int("n", n)
gen.int("p", p)
gen.int("reps", reps)
gen.num("lam", lam, 3)
gen.num("lam_node", lam_node, 3)
# the node-wise penalty is a tenth of sqrt(2 log p / n); tau_j^2 comes out below its population
# value 1 / (Sigma^{-1})_{jj} (0.75 for j = 1, 0.6 inside the AR(1) band), which inflates omega_j
tau2_pop = {name: (0.75 if j == 0 else 0.6) for name, j in targets.items()}
for name in targets:
    assert cover[name]["tau2"] < tau2_pop[name]
    gen.num(f"{name}_tau2", cover[name]["tau2"], 2)
l1_err = np.mean(np.sum(np.abs(B - beta[:, None]), axis=0))      # average ||lasso - beta||_1
gen.num("l1_err", l1_err, 2)
for name in targets:
    # the debiased bias is the remainder term, bounded by mu_j ||Delta||_1
    assert abs(cover[name]["deb_bias"]) <= cover[name]["mu"] * l1_err
for name in targets:
    c = cover[name]
    gen.num(f"{name}_deb", c["debiased"], 3)
    gen.num(f"{name}_naive", c["naive"], 3)
    gen.num(f"{name}_sel", c["selected"], 3)
    gen.num(f"{name}_oracle", c["oracle"], 3)
    gen.num(f"{name}_mu", c["mu"], 3)
    gen.num(f"{name}_lasso_bias", c["lasso_bias"], 3)
    gen.num(f"{name}_deb_bias", c["deb_bias"], 3)
    gen.num(f"{name}_width_deb", c["width_deb"], 3)
    gen.num(f"{name}_width_orc", c["width_orc"], 3)
gen.write()

# ---- figure ---------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
t = cover["strong"]["t"]
ax.hist(t, bins=40, density=True, color=COLORS["accent"], alpha=0.55, edgecolor="white", linewidth=0.3)
grid = np.linspace(-4, 4, 200)
ax.plot(grid, stats.norm.pdf(grid), color=COLORS["ink"], lw=1)
ax.set_xlabel("standardized debiased estimate of $\\beta_1$")
ax.set_title("(a) against the N(0, 1) density")
ax = axes[1]
names = ["strong", "weak", "null"]
labels = ["$\\beta_1=1$", "$\\beta_6=0.15$", "$\\beta_2=0$"]
w = 0.26
xs = np.arange(3)
for k, (m, col, lab) in enumerate([("naive", COLORS["second"], "naive, given selected"),
                                   ("debiased", COLORS["accent"], "debiased lasso"),
                                   ("oracle", COLORS["muted"], "oracle")]):
    ax.bar(xs + (k - 1) * w, [cover[nm][m] for nm in names], width=w, color=col, label=lab)
ax.axhline(0.95, color=COLORS["ink"], lw=0.8, ls="--")
ax.set_xticks(xs)
ax.set_xticklabels(labels)
ax.set_ylim(0, 1.42)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_ylabel("coverage of 95% interval")
ax.set_title("(b) coverage")
ax.legend(frameon=False, loc="upper center", fontsize=7, ncol=2, columnspacing=0.8, handlelength=1.2)
fig.tight_layout()
fig.savefig(figure_path("ch28", "debiased"))
