"""Chapter 34, Section 5: deviance, dispersion and residuals for the gamma/log fit to
Engel's household budget data.

Builds the analysis-of-deviance table by hand from the log-likelihood, compares the two
estimates of the dispersion, and draws normal quantile plots of the Pearson, deviance and
Anscombe residuals.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<anodev>>
import numpy as np
import statsmodels.api as sm
from scipy import stats

engel = sm.datasets.engel.load_pandas().data
y = engel["foodexp"].to_numpy()
u = np.log(engel["income"].to_numpy())
u = u - u.mean()                                   # centred, so the terms are less collinear
n = len(y)

def gamma_deviance(y, mu):
    """D = 2 sum { -log(y/mu) + (y - mu)/mu } for the gamma family."""
    return 2 * np.sum(-np.log(y / mu) + (y - mu) / mu)

models = {"1": np.ones((n, 1)),
          "1 + u": np.column_stack([np.ones(n), u]),
          "1 + u + u^2": np.column_stack([np.ones(n), u, u ** 2])}
fits, dev = {}, {}
for name, Xm in models.items():
    f = sm.GLM(y, Xm, family=sm.families.Gamma(sm.families.links.Log())).fit()
    fits[name], dev[name] = f, gamma_deviance(y, f.fittedvalues)

full = fits["1 + u + u^2"]
mu_full = full.fittedvalues
phi_pearson = np.sum((y - mu_full) ** 2 / mu_full ** 2) / (n - 3)
phi_deviance = dev["1 + u + u^2"] / (n - 3)

names = list(models)
print(f"{'term added':14s} {'df':>4s} {'deviance':>10s} {'change':>9s} {'F':>8s} {'p':>8s}")
prev = None
for name in names:
    df_res = n - models[name].shape[1]
    row = f"{name:14s} {df_res:4d} {dev[name]:10.4f}"
    if prev is not None:
        change = dev[prev] - dev[name]
        F = change / phi_pearson
        row += f" {change:9.4f} {F:8.2f} {stats.f.sf(F, 1, n - 3):8.4f}"
    print(row)
    prev = name
print(f"dispersion: Pearson {phi_pearson:.5f}   deviance {phi_deviance:.5f}")
# <</anodev>>

assert dev["1"] > dev["1 + u"] > dev["1 + u + u^2"] > 0
# The deviance change is exactly twice the change in the gamma log-likelihood with phi = 1.
loglik = lambda mu: np.sum(-y / mu - np.log(mu))
assert abs(2 * (loglik(mu_full) - loglik(fits["1 + u"].fittedvalues))
           - (dev["1 + u"] - dev["1 + u + u^2"])) < 1e-9

drop_u2 = dev["1 + u"] - dev["1 + u + u^2"]
F_u2 = drop_u2 / phi_pearson
p_u2 = stats.f.sf(F_u2, 1, n - 3)
drop_u = dev["1"] - dev["1 + u"]
F_u = drop_u / phi_pearson

# ---- residuals, on data simulated from a Poisson log-linear model ------------
# <<residuals>>
rng = np.random.default_rng(3405)
N = 2000
Xp = np.column_stack([np.ones(N), rng.uniform(-1, 1, N)])
mu_true = np.exp(Xp @ np.array([1.0, 0.8]))                # means between about 1 and 6
yp = rng.poisson(mu_true).astype(float)
pois = sm.GLM(yp, Xp, family=sm.families.Poisson()).fit()
mu_p = pois.fittedvalues

pearson = (yp - mu_p) / np.sqrt(mu_p)                                   # r = (y - mu)/sqrt(V)
d_i = 2 * (np.where(yp > 0, yp * np.log(np.maximum(yp, 1e-300) / mu_p), 0.0) - (yp - mu_p))
deviance_r = np.sign(yp - mu_p) * np.sqrt(d_i)                          # signed root of d_i
anscombe = 1.5 * (yp ** (2 / 3) - mu_p ** (2 / 3)) / mu_p ** (1 / 6)    # A(mu) = 1.5 mu^{2/3}

W = mu_p                                                   # working weights for Poisson/log
root = np.sqrt(W)[:, None] * Xp
H = root @ np.linalg.inv(root.T @ root) @ root.T           # the weighted hat matrix
h = np.diag(H)

for name, r in [("Pearson", pearson), ("deviance", deviance_r), ("Anscombe", anscombe)]:
    print(f"{name:9s} sd {r.std(ddof=2):6.3f}  skewness {stats.skew(r):7.3f}")
print(f"leverages: sum {h.sum():.4f}, largest {h.max():.4f}")
# <</residuals>>

assert np.all(d_i >= -1e-12)                       # every deviance contribution is nonnegative
assert abs(np.sum(d_i) - pois.deviance) < 1e-8
assert stats.skew(pearson) > stats.skew(deviance_r) > stats.skew(anscombe) > -0.1
assert abs(h.sum() - 2) < 1e-8 and h.max() < 1

gen = Generated("ch34", "deviance")
gen.int("n", n)
gen.num("dev0", dev["1"], 2)
gen.num("dev1", dev["1 + u"], 2)
gen.num("dev2", dev["1 + u + u^2"], 2)
gen.num("drop_u", drop_u, 2)
gen.num("drop_u2", drop_u2, 4)
gen.num("F_u", F_u, 1)
gen.num("F_u2", F_u2, 2)
gen.num("p_u2", p_u2, 3)
gen.num("phi_pearson", phi_pearson, 5)
gen.num("phi_deviance", phi_deviance, 5)
gen.num("skew_pearson", float(stats.skew(pearson)), 3)
gen.num("skew_deviance", float(stats.skew(deviance_r)), 3)
gen.num("skew_anscombe", float(stats.skew(anscombe)), 3)
gen.num("hmax", float(h.max()), 4)
gen.int("n_sim", N)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.2), sharey=True)
for ax, (name, r, col) in zip(axes, [("Pearson", pearson, COLORS["accent"]),
                                     ("deviance", deviance_r, COLORS["second"]),
                                     ("Anscombe", anscombe, COLORS["third"])]):
    q = stats.norm.ppf((np.arange(1, N + 1) - 0.5) / N)
    ax.scatter(q, np.sort(r), s=6, color=col, alpha=0.75, linewidths=0)
    lim = [-3.2, 3.2]
    ax.plot(lim, lim, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_xlim(lim)
    ax.set_ylim(-4.2, 4.2)
    ax.set_xlabel("normal quantile")
    ax.set_title(f"{name}, skewness {stats.skew(r):.2f}", fontsize=8)
axes[0].set_ylabel("ordered residual")
fig.tight_layout()
fig.savefig(figure_path("ch34", "poisson_residuals"))
