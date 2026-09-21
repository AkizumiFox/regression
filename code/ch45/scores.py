"""Chapter 45, Section 6: comparing three distributional fits to Engel's data honestly.

Everything is out of sample: the same ten folds carry a homoscedastic linear model, a
grid of quantile regressions (rearranged so that the fitted quantiles are ordered)
and the location-scale GAMLSS of Section 45.5. Each fit is judged by the continuous
ranked probability score, computed for all three by the same quantile integral, by
the pinball loss at three levels, and by a probability integral transform histogram.
Differences in score carry a paired standard error.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.interpolate import BSpline
from scipy.optimize import linprog

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data
DEGREE, N_INNER = 3, 8
LEVELS = np.arange(1, 200) / 200.0          # the tau grid the CRPS integral uses
PINBALL = (0.10, 0.50, 0.90)

income = data["income"].to_numpy() / 1000.0
food = data["foodexp"].to_numpy() / 1000.0
n = len(food)
lo, hi = income.min(), income.max()
KNOTS = np.r_[[lo] * (DEGREE + 1),
              np.quantile(income, np.linspace(0, 1, N_INNER + 2)[1:-1]),
              [hi] * (DEGREE + 1)]


def basis(z):
    return np.asarray(BSpline.design_matrix(np.clip(np.atleast_1d(z), lo, hi),
                                            KNOTS, DEGREE).todense())


def difference_penalty(q, order=2):
    D = np.diff(np.eye(q), order, axis=0)
    return D.T @ D


def qreg(X, y, tau):
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]


def gamlss_normal(Bm, Bs, y, lam_mu, lam_sig, iters=100, tol=1e-9):
    Km, Ks = difference_penalty(Bm.shape[1]), difference_penalty(Bs.shape[1])
    b_mu = np.linalg.solve(Bm.T @ Bm + lam_mu * Km, Bm.T @ y)
    b_sig = np.linalg.solve(Bs.T @ Bs, Bs.T @ np.full(len(y), np.log(y.std())))
    for _ in range(iters):
        sigma = np.exp(np.clip(Bs @ b_sig, -8, 4))
        w = 1.0 / sigma ** 2
        new_mu = np.linalg.solve(Bm.T @ (w[:, None] * Bm) + lam_mu * Km, Bm.T @ (w * y))
        u = (y - Bm @ new_mu) ** 2 / sigma ** 2 - 1.0
        new_sig = b_sig + np.linalg.solve(2 * Bs.T @ Bs + lam_sig * Ks,
                                          Bs.T @ u - lam_sig * Ks @ b_sig)
        step = max(np.max(np.abs(new_mu - b_mu)), np.max(np.abs(new_sig - b_sig)))
        b_mu, b_sig = new_mu, new_sig
        if step < tol:
            break
    return b_mu, b_sig


# <<scores>>
def check_loss(u, tau):
    """The check function of Section 45.2."""
    return u * (tau - (u < 0))


def crps_from_quantiles(q, y):
    """CRPS = 2 * integral of the check loss over tau, on the grid LEVELS.

    q has one row per observation and one column per level of LEVELS.
    """
    return 2 * np.mean(check_loss(y[:, None] - q, LEVELS[None, :]), axis=1)


def pit_from_quantiles(q, y):
    """The probability integral transform read off a set of fitted quantiles."""
    return np.mean(q <= y[:, None], axis=1)


# <</scores>>

# <<demo>>
# one household: a normal predictive with mean 600 and scale 120 francs, observed at 700
q_one = 600.0 + 120.0 * stats.norm.ppf(LEVELS)[None, :]
y_one = np.array([700.0])
print(f"CRPS {crps_from_quantiles(q_one, y_one)[0]:.2f} francs,"
      f" PIT {pit_from_quantiles(q_one, y_one)[0]:.3f}")
# <</demo>>

# the closed form of exercise B1 in Section 45.6, as a check on the 199-level grid
z_one = (700.0 - 600.0) / 120.0
exact_one = 120.0 * (z_one * (2 * stats.norm.cdf(z_one) - 1) + 2 * stats.norm.pdf(z_one)
                     - 1 / np.sqrt(np.pi))
assert abs(crps_from_quantiles(q_one, y_one)[0] - exact_one) < 0.01 * exact_one

folds = np.random.default_rng(7).permutation(n) % 10
X = np.column_stack([np.ones(n), income])
out = {name: {"crps": np.empty(n), "pit": np.empty(n), "q": np.empty((n, len(LEVELS)))}
       for name in ("mean model", "quantile model", "GAMLSS")}

for k in range(10):
    train, test = folds != k, folds == k
    ytr, yte = food[train], food[test]

    # (i) the homoscedastic linear model of Part II, with the exact t predictive of
    # Section 12.3: the leverage of the held-out row enters through sqrt(1 + h)
    b = np.linalg.lstsq(X[train], ytr, rcond=None)[0]
    df = train.sum() - 2
    s = np.sqrt(np.sum((ytr - X[train] @ b) ** 2) / df)
    G = np.linalg.inv(X[train].T @ X[train])
    lev = np.einsum("ij,jk,ik->i", X[test], G, X[test])
    q_mean = ((X[test] @ b)[:, None]
              + (s * np.sqrt(1 + lev))[:, None] * stats.t.ppf(LEVELS, df)[None, :])

    # (ii) a quantile regression at every level, sorted so that the quantiles are ordered
    fits = np.array([qreg(X[train], ytr, t) for t in LEVELS])
    q_qr = np.sort(X[test] @ fits.T, axis=1)

    # (iii) the location-scale GAMLSS of Section 45.5
    b_mu, b_sig = gamlss_normal(basis(income[train]), basis(income[train]), ytr, 31.6, 10.0)
    Bt = basis(income[test])
    q_gam = (Bt @ b_mu)[:, None] + np.exp(Bt @ b_sig)[:, None] * stats.norm.ppf(LEVELS)[None, :]

    for name, q in (("mean model", q_mean), ("quantile model", q_qr), ("GAMLSS", q_gam)):
        out[name]["crps"][test] = crps_from_quantiles(q, yte)
        out[name]["pit"][test] = pit_from_quantiles(q, yte)
        out[name]["q"][test] = q

for name in out:
    out[name]["pinball"] = {
        tau: np.mean(check_loss(food - out[name]["q"][:, int(round(tau * 200)) - 1], tau))
        for tau in PINBALL}
    out[name]["cover"] = {
        tau: np.mean(food <= out[name]["q"][:, int(round(tau * 200)) - 1]) for tau in PINBALL}

for name in ("mean model", "quantile model", "GAMLSS"):
    print(f"{name:15s} CRPS {1000 * out[name]['crps'].mean():6.2f} francs   "
          + "  ".join(f"pinball({t}) {1000 * out[name]['pinball'][t]:6.2f}" for t in PINBALL))

# paired standard errors for the differences in CRPS
base = out["mean model"]["crps"]
diffs = {name: out[name]["crps"] - base for name in ("quantile model", "GAMLSS")}
diffs["quantile model - GAMLSS"] = out["quantile model"]["crps"] - out["GAMLSS"]["crps"]
for name, d in diffs.items():
    print(f"{name:24s} CRPS difference {1000 * d.mean():+6.2f}"
          f" (paired se {1000 * d.std(ddof=1) / np.sqrt(n):.2f})")

assert out["GAMLSS"]["crps"].mean() < out["mean model"]["crps"].mean()
assert out["quantile model"]["crps"].mean() < out["mean model"]["crps"].mean()
assert abs(out["mean model"]["cover"][0.10] - 0.10) > abs(out["GAMLSS"]["cover"][0.10] - 0.10)

# a Kolmogorov-Smirnov distance of the PIT values from uniformity
ks = {name: stats.kstest(out[name]["pit"], "uniform").statistic for name in out}

gen = Generated("ch45", "scores")
gen.int("n", n)
gen.int("levels", len(LEVELS))
short = {"mean model": "mean", "quantile model": "qr", "GAMLSS": "gam",
         "quantile model - GAMLSS": "qrgam"}
for name in out:
    key = short[name]
    gen.num(f"crps_{key}", 1000 * out[name]["crps"].mean(), 2)
    gen.num(f"ks_{key}", ks[name], 3)
    for tau in PINBALL:
        t = f"{round(100 * tau):02d}"
        gen.num(f"pin_{key}_{t}", 1000 * out[name]["pinball"][tau], 2)
        gen.num(f"cov_{key}_{t}", out[name]["cover"][tau], 3)
for name, d in diffs.items():
    key = short[name]
    gen.num(f"gap_{key}", 1000 * d.mean(), 2)
    gen.num(f"gapse_{key}", 1000 * d.std(ddof=1) / np.sqrt(n), 2)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.6, 2.0), sharey=True)
edges = np.linspace(0, 1, 11)
for ax, name in zip(axes, ("mean model", "quantile model", "GAMLSS")):
    ax.hist(out[name]["pit"], bins=edges, color=COLORS["accent"], alpha=0.55,
            edgecolor=COLORS["ink"], linewidth=0.5)
    ax.axhline(n / 10, color=COLORS["second"], linestyle="--", linewidth=1.0)
    ax.set_xlabel("PIT value")
    ax.set_title(name)
axes[0].set_ylabel("households")
fig.tight_layout()
fig.savefig(figure_path("ch45", "pit_histograms"))
