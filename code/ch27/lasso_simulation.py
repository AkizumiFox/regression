"""Chapter 27, Section 4: lasso, ridge and best subset compared by simulation.

Fixed design: n = 50 rows drawn once from N(0, Sigma) with Sigma_ij = 0.5^|i-j|, p = 10, columns
centred and scaled to ||x_j||^2 = n, no intercept. Errors N(0, sigma^2) with sigma = 2. Three
coefficient vectors: sparse and strong, dense and weak, and in between. Each method's tuning
parameter is chosen on an independent copy of the responses at the same design (an idealized
validation set). Reported: E||X (b - beta)||^2 / (p sigma^2), which is 1 for least squares.
The script uses 300 data sets per scenario; the runnable cell uses 20.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
def soft(z, t):
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)


def lasso_cd(X, y, lam, b=None, tol=1e-8, max_sweeps=10_000):
    """Cyclic coordinate descent for (1/2)||y - X b||^2 + lam ||b||_1."""
    p = X.shape[1]
    b = np.zeros(p) if b is None else b.copy()
    r = y - X @ b
    sq = np.sum(X ** 2, axis=0)
    for _ in range(max_sweeps):
        biggest = 0.0
        for j in range(p):
            new = soft(X[:, j] @ r + sq[j] * b[j], lam) / sq[j]
            if new != b[j]:
                r -= X[:, j] * (new - b[j])
                biggest = max(biggest, abs(new - b[j]))
                b[j] = new
        if biggest < tol:
            break
    return b


rng = np.random.default_rng(2027)
n, p, sigma = 50, 10, 2.0
Sigma = 0.5 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
X = rng.multivariate_normal(np.zeros(p), Sigma, size=n)
X = (X - X.mean(axis=0)) / np.sqrt(np.mean((X - X.mean(axis=0)) ** 2, axis=0))   # ||x_j||^2 = n
U, d, Vt = np.linalg.svd(X, full_matrices=False)
subsets = [list(S) for k in range(p + 1) for S in itertools.combinations(range(p), k)]
ridge_grid = np.logspace(-2, 4, 40)
scenarios = {
    "sparse": np.r_[3.0, -2.0, 1.5, np.zeros(p - 3)],
    "dense": 0.5 * (-1.0) ** np.arange(p),
    "intermediate": np.r_[1.0, -1.0, 1.0, -1.0, 1.0, np.zeros(p - 5)],
}


def fits(y):
    """Candidate fitted vectors of each method along its tuning path."""
    c = U.T @ y
    out = {"least squares": [U @ c]}
    out["ridge"] = [U @ (d ** 2 / (d ** 2 + lam) * c) for lam in ridge_grid]
    lam_max = np.max(np.abs(X.T @ y))
    b, lasso = np.zeros(p), []
    for lam in lam_max * np.logspace(0, -3, 30):
        b = lasso_cd(X, y, lam, b)
        lasso.append(X @ b)
    out["lasso"] = lasso
    best = {}                                    # best subset of each size, by residual sum of squares
    for S in subsets:
        fit = X[:, S] @ np.linalg.lstsq(X[:, S], y, rcond=None)[0] if S else np.zeros(n)
        rss = np.sum((y - fit) ** 2)
        if len(S) not in best or rss < best[len(S)][0]:
            best[len(S)] = (rss, fit)
    out["best subset"] = [best[k][1] for k in range(p + 1)]
    return out


def run(beta, reps, seed=1):
    """Mean of ||X b - X beta||^2 / (p sigma^2) for each method, tuned on a validation copy."""
    gen = np.random.default_rng(seed)
    mu = X @ beta
    loss = {m: [] for m in ["least squares", "ridge", "lasso", "best subset"]}
    for _ in range(reps):
        y = mu + sigma * gen.normal(size=n)
        y_val = mu + sigma * gen.normal(size=n)
        for m, cands in fits(y).items():
            chosen = min(cands, key=lambda f: np.sum((y_val - f) ** 2))
            loss[m].append(np.sum((chosen - mu) ** 2) / (p * sigma ** 2))
    return {m: np.array(v) for m, v in loss.items()}
# <</setup>>


# <<quick>>
quick = run(scenarios["sparse"], reps=10)
print({m: round(float(v.mean()), 2) for m, v in quick.items()})
# <</quick>>

reps = 300
losses = {name: run(beta, reps, seed=k + 10) for k, (name, beta) in enumerate(scenarios.items())}
results = {name: {m: v.mean() for m, v in res.items()} for name, res in losses.items()}
diff = losses["intermediate"]["lasso"] - losses["intermediate"]["ridge"]
se_diff = diff.std(ddof=1) / np.sqrt(reps)                 # Monte Carlo standard error of a paired difference
se_ls = losses["sparse"]["least squares"].std(ddof=1) / np.sqrt(reps)
print(f"intermediate: lasso - ridge = {diff.mean():.3f} (standard error {se_diff:.3f})")
for name, res in results.items():
    print(f"{name:13s}", "  ".join(f"{m}: {v:.3f}" for m, v in res.items()))

methods = ["least squares", "ridge", "lasso", "best subset"]
assert abs(results["dense"]["least squares"] - 1) < 0.1 and abs(results["sparse"]["least squares"] - 1) < 0.1
assert results["sparse"]["lasso"] < results["sparse"]["ridge"]
assert results["sparse"]["best subset"] < results["sparse"]["ridge"]
assert results["dense"]["ridge"] < results["dense"]["lasso"] < results["dense"]["best subset"]

use_book_style()
fig, ax = plt.subplots(figsize=(5.0, 2.4))
width = 0.2
cols = [COLORS["muted"], COLORS["third"], COLORS["accent"], COLORS["second"]]
for k, m in enumerate(methods):
    vals = [results[s][m] for s in scenarios]
    ax.bar(np.arange(3) + (k - 1.5) * width, vals, width=width, color=cols[k], label=m)
ax.axhline(1, color=COLORS["ink"], linewidth=0.6)
ax.set_xticks(np.arange(3))
ax.set_xticklabels(list(scenarios))
ax.set_ylabel(r"prediction risk / $p\sigma^2$")
ax.set_ylim(0, 1.35)
ax.legend(frameon=False, fontsize=6.5, ncol=4, loc="upper center")
fig.tight_layout()
fig.savefig(figure_path("ch27", "lasso_simulation"))

gen = Generated("ch27", "lasso_simulation")
gen.int("reps", reps)
gen.num("diff", diff.mean(), 2)
gen.num("se_diff", se_diff, 2)
gen.num("se_ls", se_ls, 2)
for s in scenarios:
    for m in methods:
        gen.num(f"{s}_{m.replace(' ', '_')}", results[s][m], 2)
gen.write()
