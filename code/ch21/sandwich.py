"""Chapter 21, Section 4: heteroscedasticity-consistent (sandwich) standard errors.

(a) Engel's food expenditure data (public domain, statsmodels.datasets.engel): classical and
HC0-HC3 standard errors of the slope, and the identity between HC3 and the (uncentred)
jackknife variance.
(b) Coverage of nominal 95% intervals for a slope, by simulation. The regressor values are
the n quantiles of a standard lognormal law (a skewed design whose largest leverage falls
slowly with n); errors are normal with standard deviation 1 or x_i.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<engel>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])
Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ food)
e = food - X @ beta
h = np.sum(Q**2, axis=1)                              # leverages
Rinv = np.linalg.inv(R)                               # (X'X)^{-1} = Rinv Rinv'


def sandwich(omega):
    """(X'X)^{-1} X' diag(omega) X (X'X)^{-1}, computed through Q."""
    meat = (Q * omega[:, None]).T @ Q
    return Rinv @ meat @ Rinv.T


p = X.shape[1]
covs = {
    "classical": (e @ e / (n - p)) * (Rinv @ Rinv.T),
    "HC0": sandwich(e**2),
    "HC1": sandwich(e**2) * n / (n - p),
    "HC2": sandwich(e**2 / (1 - h)),
    "HC3": sandwich(e**2 / (1 - h) ** 2),
}
for name, V in covs.items():
    print(f"{name:9s} slope se {np.sqrt(V[1, 1]):.4f}")
# <</engel>>

fit = sm.OLS(food, X).fit()
for name in ("HC0", "HC1", "HC2", "HC3"):
    assert np.allclose(covs[name], fit.get_robustcov_results(name).cov_params())

# <<jackknife>>
loo = np.array([np.linalg.lstsq(np.delete(X, i, 0), np.delete(food, i),
                                rcond=None)[0]
                for i in range(n)])
# sum over i of (b_(i) - b)(b_(i) - b)'
jack = (loo - beta).T @ (loo - beta)
print("HC3 equals the sum of squared deletion changes:",
      np.allclose(jack, covs["HC3"]))
# <</jackknife>>
assert np.allclose(jack, covs["HC3"])

i_max = int(np.argmax(h))
gen = Generated("ch21", "sandwich", prefix="sw")
for name, V in covs.items():
    gen.num(f"se_{name}", np.sqrt(V[1, 1]), 4)
gen.num("ratio_hc3", np.sqrt(covs["HC3"][1, 1] / covs["classical"][1, 1]), 1)
gen.num("h_max", h[i_max], 3)
gen.num("share_max", (e[i_max] ** 2 * (Q[i_max] @ Rinv.T[:, 1]) ** 2 / (1 - h[i_max]) ** 2)
        / covs["HC3"][1, 1], 2)

# ---- (b) coverage by simulation ---------------------------------------------------------
# <<coverage>>
rng = np.random.default_rng(2107)


def coverage(n, hetero, reps):
    """Coverage of nominal 95% slope intervals, classical and HC0-HC3."""
    # lognormal quantiles
    x = np.exp(stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n))
    Xs = np.column_stack([np.ones(n), x])
    Qs, Rs = np.linalg.qr(Xs)
    hs = np.sum(Qs**2, axis=1)
    # slope estimate - slope = c @ errors
    c = Qs @ np.linalg.inv(Rs).T[:, 1]
    E = rng.normal(size=(reps, n)) * (x if hetero else 1.0)
    err = E @ c
    res = E - (E @ Qs) @ Qs.T
    var = {"classical": (res**2).sum(axis=1) / (n - 2) * (c @ c),
           "HC0": (res**2) @ c**2,
           "HC1": (res**2) @ c**2 * n / (n - 2),
           "HC2": (res**2 / (1 - hs)) @ c**2,
           "HC3": (res**2 / (1 - hs) ** 2) @ c**2}
    q = stats.t.ppf(0.975, n - 2)
    cover = {k: np.mean(np.abs(err) <= q * np.sqrt(v)) for k, v in var.items()}
    return cover, hs.max()


for n_sim in (25, 100, 400):
    cov, hmax = coverage(n_sim, hetero=True, reps=2000)
    print(f"n = {n_sim:3d} (max leverage {hmax:.2f}):",
          {k: round(v, 3) for k, v in cov.items()})
# <</coverage>>

ns = [15, 25, 50, 100, 200, 400, 1000]
res = {hetero: {m: coverage(m, hetero, 20_000) for m in ns} for hetero in (False, True)}
for hetero, tag in ((False, "hom"), (True, "het")):
    for m in (25, 1000):
        cov, hmax = res[hetero][m]
        for k, v in cov.items():
            gen.num(f"{tag}_{k}_{m}", v, 3)
        gen.num(f"hmax_{m}", hmax, 2)
het = res[True]
# classical gets worse with n
assert het[1000][0]["classical"] < het[25][0]["classical"]
assert all(het[m][0]["HC3"] > het[m][0]["HC2"] > het[m][0]["HC0"] for m in ns)
assert all(abs(res[False][m][0]["classical"] - 0.95) < 0.01 for m in ns)
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5), sharey=True)
styles = {"classical": (COLORS["muted"], "-"), "HC0": (COLORS["second"], "-"),
          "HC2": (COLORS["third"], "--"), "HC3": (COLORS["accent"], "-")}
for ax, hetero, title in zip(axes, (False, True), ("(a) constant variance", "(b) sd proportional to $x$")):
    for k, (color, ls) in styles.items():
        ax.plot(ns, [res[hetero][m][0][k] for m in ns], color=color, ls=ls, marker="o",
                markersize=2.5, label=k)
    ax.axhline(0.95, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_xticks([15, 50, 200, 1000], ["15", "50", "200", "1000"])
    ax.minorticks_off()
    ax.set_xlabel("sample size $n$")
    ax.set_title(title)
axes[0].set_ylabel("coverage of nominal 95%")
axes[0].set_ylim(0.2, 1.0)
axes[1].legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch21", "sandwich_coverage"))
