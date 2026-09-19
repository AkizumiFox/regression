"""Chapter 27, Section 3: principal component regression on Longley's data.

Same standardization as ridge_longley.py (regressors centred with unit standard deviation,
response centred, intercept unpenalized). For k = 0, ..., 6 components the script computes the
PCR coefficients and leave-one-out prediction error (PCR is least squares on the component
scores, so the deletion formula applies), the t statistics of the components, and compares
with ridge and with the components chosen by their t statistics.
"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<pcr>>
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape
X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
yc = y - y.mean()
U, d, Vt = np.linalg.svd(X, full_matrices=False)
c = U.T @ yc                                   # c_j = u_j^T y, the fitted coordinate of component j


def fit_components(keep):
    """Least squares on the chosen component scores: coefficients and leave-one-out error."""
    keep = list(keep)
    beta = Vt[keep].T @ (c[keep] / d[keep])
    resid = yc - X @ beta
    h = 1 / n + np.sum(U[:, keep] ** 2, axis=1)            # leverages, intercept included
    return beta, np.mean((resid / (1 - h)) ** 2)


resid_ls = yc - U @ c
s = np.sqrt(resid_ls @ resid_ls / (n - p - 1))
t = c / s                                        # t statistic of each component
share = d ** 2 / np.sum(d ** 2)                  # share of the regressors' variance
for k in range(p + 1):
    beta_k, loo_k = fit_components(range(k))
    print(f"k = {k}  LOO = {loo_k:9.0f}   coefficients", np.round(beta_k, 0))
print("variance shares:", np.round(share, 5))
print("component t statistics:", np.round(t, 2))
# <</pcr>>

# the deletion shortcut against brute force, for k = 3
W3 = np.column_stack([np.ones(n), X @ Vt[:3].T])
errs = []
for i in range(n):
    keep = np.arange(n) != i
    bi = np.linalg.lstsq(W3[keep], y[keep], rcond=None)[0]
    errs.append(y[i] - W3[i] @ bi)
assert np.isclose(np.mean(np.square(errs)), fit_components(range(3))[1])
# PCR with k = p is least squares; PCR solves least squares under v_j^T b = 0 for j > k
beta_ls = np.linalg.lstsq(X, yc, rcond=None)[0]
assert np.allclose(fit_components(range(p))[0], beta_ls)
k = 3
Vdrop = Vt[k:].T
K = np.block([[X.T @ X, Vdrop], [Vdrop.T, np.zeros((p - k, p - k))]])
sol = np.linalg.solve(K, np.r_[X.T @ yc, np.zeros(p - k)])
assert np.allclose(sol[:p], fit_components(range(k))[0])

loo = np.array([fit_components(range(k))[1] for k in range(p + 1)])
k_best = int(np.argmin(loo))
chosen = [j for j in range(p) if t[j] ** 2 > 2]        # keep components with t_j^2 > 2
beta_t, loo_t = fit_components(chosen)
# best subset of components by LOO (all 64 subsets)
subsets = [(fit_components(S)[1], S) for r in range(p + 1) for S in itertools.combinations(range(p), r)]
loo_sub, best_sub = min(subsets)
print("LOO-best number of leading components:", k_best, " components with t^2 > 2:", chosen,
      f"LOO {loo_t:.0f};  best subset of components {best_sub} LOO {loo_sub:.0f}")
assert abs(t[4]) > 3 and share[4] < 0.001            # a minor component that matters
assert loo_t < loo[k_best] and chosen == [0, 1, 2, 4, 5]


# honest assessment of the t-rule: nested leave-one-out, with standardization, components, s and
# the selection t_j^2 > 2 all recomputed from the fifteen training years
def pcr_t_predict(Zt, yt, z0):
    nt = len(yt)
    mt, st = Zt.mean(axis=0), Zt.std(axis=0, ddof=1)
    Ut, dt, Vtt = np.linalg.svd((Zt - mt) / st, full_matrices=False)
    ct = Ut.T @ (yt - yt.mean())
    rt = yt - yt.mean() - Ut @ ct
    keep = (ct / np.sqrt(rt @ rt / (nt - p - 1))) ** 2 > 2
    bt = Vtt[keep].T @ (ct[keep] / dt[keep])
    return yt.mean() + ((z0 - mt) / st) @ bt


nested_t = [y[i] - pcr_t_predict(Z[np.arange(n) != i], y[np.arange(n) != i], Z[i]) for i in range(n)]
rmse_t_nested = np.sqrt(np.mean(np.square(nested_t)))
print(f"t-rule: LOO with the selection fixed {np.sqrt(loo_t):.0f}, nested {rmse_t_nested:.0f}")
assert rmse_t_nested > np.sqrt(loo_t)


def ridge_loo(lam):
    f = d ** 2 / (d ** 2 + lam)
    resid = yc - U @ (f * c)
    h = 1 / n + np.sum(U ** 2 * f, axis=1)
    return np.mean((resid / (1 - h)) ** 2), np.sum(f)


lams = np.logspace(-5, 3, 400)
rl = np.array([ridge_loo(l) for l in lams])
ridge_best = rl[:, 0].min()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.6))
ax = axes[0]
idx = np.arange(1, p + 1)
ax.bar(idx - 0.2, share, width=0.4, color=COLORS["accent"], label="share of variance")
ax.set_yscale("log")
ax.set_xlabel("component $j$")
ax.set_ylabel("share of variance", color=COLORS["accent"])
ax2 = ax.twinx()
ax2.bar(idx + 0.2, np.abs(t), width=0.4, color=COLORS["second"], label="$|t_j|$")
ax2.axhline(np.sqrt(2), color=COLORS["muted"], linestyle="--", linewidth=0.8)
ax2.set_ylabel("$|t_j|$", color=COLORS["second"])
ax2.spines["right"].set_visible(True)
ax.set_title("(a) variance and signal by component")
ax = axes[1]
ax.plot(rl[:, 1], np.sqrt(rl[:, 0]), color=COLORS["accent"], label="ridge")
ax.plot(np.arange(p + 1), np.sqrt(loo), "o", color=COLORS["second"], markersize=4, label="PCR, leading $k$")
ax.plot([len(chosen)], [np.sqrt(loo_t)], "s", color=COLORS["third"], markersize=5, label="components with $t_j^2>2$")
ax.set_yscale("log")
ax.set_xlabel("(effective) degrees of freedom")
ax.set_ylabel("LOO root mean squared error")
ax.set_title("(b) leave-one-out error")
ax.legend(frameon=False, fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch27", "pcr_longley"))

gen = Generated("ch27", "pcr")
for j in range(p):
    gen.num(f"share_{j + 1}", 100 * share[j], 3)
    gen.num(f"t_{j + 1}", t[j], 2)
for k in range(p + 1):
    gen.num(f"rmse_{k}", np.sqrt(loo[k]), 0)
gen.int("k_best", k_best)
gen.num("rmse_t", np.sqrt(loo_t), 0)
gen.num("rmse_ridge", np.sqrt(ridge_best), 0)
gen.text("best_sub", ", ".join(str(j + 1) for j in best_sub))
gen.num("rmse_sub", np.sqrt(loo_sub), 0)
gen.num("rmse_t_nested", rmse_t_nested, 0)
for j, name in enumerate(names):
    gen.num(f"k3_{name}", fit_components(range(3))[0][j], 0)
gen.write()
