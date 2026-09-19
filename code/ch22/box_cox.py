"""Chapter 22, Section 2: the Box-Cox profile likelihood on Engel's data.

The profile log-likelihood of lambda is -(n/2) log SSE_z(lambda) + const, where SSE_z is the
residual sum of squares of the geometric-mean-normalized response. Two mean models: food
expenditure linear in income, and linear in log income. Also: invariance to the units of y,
the constructed variable, and a small simulation of the likelihood-ratio interval.
Public-domain data shipped with statsmodels (statsmodels.datasets.engel).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<profile>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)


def box_cox_z(y, lam):
    """Box-Cox transform normalized by the geometric mean, so every lambda is in the units of y."""
    gm = np.exp(np.mean(np.log(y)))
    if abs(lam) < 1e-12:
        return gm * np.log(y)
    return np.expm1(lam * np.log(y)) / (lam * gm ** (lam - 1))   # (y^lam - 1)/(lam gm^(lam-1))


def sse(v, X):
    coef, *_ = np.linalg.lstsq(X, v, rcond=None)
    r = v - X @ coef
    return r @ r


def profile(lam, y, X):
    """Profile log-likelihood of lambda, up to an additive constant."""
    return -0.5 * len(y) * np.log(sse(box_cox_z(y, lam), X))


def box_cox_fit(y, X, level=0.95):
    """MLE of lambda and the likelihood-ratio interval {lambda: 2(l_max - l_p) <= chi2_1}."""
    lam_hat = optimize.minimize_scalar(lambda l: -profile(l, y, X), bounds=(-3, 3), method="bounded",
                                       options={"xatol": 1e-10}).x
    l_max = profile(lam_hat, y, X)
    cut = stats.chi2.ppf(level, 1) / 2
    g = lambda l: l_max - profile(l, y, X) - cut
    lo = optimize.brentq(g, lam_hat - 3, lam_hat)
    hi = optimize.brentq(g, lam_hat, lam_hat + 3)
    return lam_hat, (lo, hi), l_max


X_lin = np.column_stack([np.ones(n), income])            # food linear in income
X_log = np.column_stack([np.ones(n), np.log(income)])    # food linear in log income
for name, X in [("income", X_lin), ("log income", X_log)]:
    lam_hat, (lo, hi), l_max = box_cox_fit(food, X)
    lr = {lam: 2 * (l_max - profile(lam, food, X)) for lam in (0, 0.5, 1)}
    print(f"{name:10s} lambda_hat {lam_hat:.3f}  95% interval [{lo:.3f}, {hi:.3f}]  LR at 0, 1/2, 1:",
          " ".join(f"{v:.1f}" for v in lr.values()))
# <</profile>>

fits = {name: box_cox_fit(food, X) for name, X in [("lin", X_lin), ("log", X_log)]}
lam_lin, ci_lin, lmax_lin = fits["lin"]
lam_log, ci_log, lmax_log = fits["log"]
lr_lin = {lam: 2 * (lmax_lin - profile(lam, food, X_lin)) for lam in (0, 0.5, 1)}
lr_log = {lam: 2 * (lmax_log - profile(lam, food, X_log)) for lam in (0, 0.5, 1)}
crit = stats.chi2.ppf(0.95, 1)
assert ci_lin[0] > 0 and ci_lin[1] < 1                  # neither log nor raw with income linear
assert ci_log[0] < 0 < ci_log[1]                         # log accepted with log income
assert lr_lin[0] > crit and lr_lin[1] > crit and lr_log[1] > crit and lr_log[0.5] > crit

# units of y do not matter (the model has an intercept)
lam_units, ci_units, _ = box_cox_fit(food / 100, X_log)
assert abs(lam_units - lam_log) < 1e-6 and np.allclose(ci_units, ci_log, atol=1e-6)
prof_diff = [profile(l, food / 100, X_log) - profile(l, food, X_log) for l in (-1, 0, 0.7, 2)]
assert np.ptp(prof_diff) < 1e-8                          # profiles differ by a constant

# residual SD on the normalized scale is in the units of food expenditure, comparable across lambda
s_z = {lam: np.sqrt(sse(box_cox_z(food, lam), X_log) / (n - 2)) for lam in (0, 1)}


# ---- the constructed variable ----------------------------------------------------------------
# <<constructed>>
def constructed(y, lam0, h=1e-5):
    """w = d z(lambda)/d lambda at lam0 (numerical derivative of the normalized transform)."""
    return (box_cox_z(y, lam0 + h) - box_cox_z(y, lam0 - h)) / (2 * h)


def constructed_test(y, X, lam0):
    """Regress z(lam0) on [X, w]: the coefficient of w and its t statistic."""
    w = constructed(y, lam0)
    fit = sm.OLS(box_cox_z(y, lam0), np.column_stack([X, w])).fit()
    return fit.params[-1], fit.tvalues[-1]


for lam0 in (1, 0):
    g, t = constructed_test(food, X_log, lam0)
    print(f"lambda0 = {lam0}: coefficient of w {g:.3f}, t = {t:.2f}, one-step estimate {lam0 - g:.3f}")
# <</constructed>>

g1, t1 = constructed_test(food, X_log, 1)
g0, t0 = constructed_test(food, X_log, 0)
g_hat, t_hat = constructed_test(food, X_log, lam_log)
assert abs(g_hat) < 1e-4                                  # stationary point: w has zero coefficient
assert abs(t1) > 3 and abs(t0) < 1
# closed forms of w, up to columns of X: y(log(y/gm) - 1) at 1, (gm/2) log(y/gm)^2 at 0
gm = np.exp(np.mean(np.log(food)))
for lam0, w_closed in [(1, food * (np.log(food / gm) - 1)), (0, gm / 2 * np.log(food / gm) ** 2)]:
    d = constructed(food, lam0) - w_closed
    assert np.ptp(d) < 1e-4 * np.ptp(w_closed)           # they differ by a constant

# ---- influence of single households on lambda_hat (log-income model) --------------------------
drops = [box_cox_fit(np.delete(food, i), np.delete(X_log, i, axis=0)) for i in range(n)]
lam_drop = np.array([d[0] for d in drops])
i_max = int(np.argmax(np.abs(lam_drop - lam_log)))
assert all(lo <= 0 <= hi for _, (lo, hi), _ in drops)   # the log stays in the interval after any deletion
assert food[i_max] == food.max()                         # the household with the largest food expenditure


# ---- simulation: coverage of the LR interval and size of the constructed-variable test ---------
def simulate(reps, seed):
    rng = np.random.default_rng(seed)
    fit = sm.OLS(np.log(food), X_log).fit()           # the log-log fit: data generated with lambda = 0
    mu, sigma = fit.fittedvalues, np.sqrt(fit.scale)
    cover = reject = 0
    for _ in range(reps):
        y_sim = np.exp(mu + sigma * rng.standard_normal(n))
        lam_s, (lo_s, hi_s), _ = box_cox_fit(y_sim, X_log)
        cover += lo_s <= 0 <= hi_s
        reject += abs(constructed_test(y_sim, X_log, 0)[1]) > stats.t.ppf(0.975, n - 3)
    return cover / reps, reject / reps


# <<simulation>>
coverage, size = simulate(reps=200, seed=2201)
print(f"coverage of the 95% LR interval {coverage:.3f}, size of the 5% constructed-variable test {size:.3f}")
# <</simulation>>
coverage, size = simulate(reps=4000, seed=2202)
assert 0.93 < coverage < 0.96          # the LR interval: close to 95%
assert 0.06 < size < 0.09              # the constructed-variable t test: liberal, since w depends on y

gen = Generated("ch22", "box_cox", prefix="bc")
gen.num("lam_lin", lam_lin, 2)
gen.num("lo_lin", ci_lin[0], 2)
gen.num("hi_lin", ci_lin[1], 2)
gen.num("lam_log", lam_log, 3)
gen.num("lo_log", ci_log[0], 3)
gen.num("hi_log", ci_log[1], 3)
gen.num("lr_lin_0", lr_lin[0], 1)
gen.num("lr_lin_1", lr_lin[1], 1)
gen.num("lr_log_0", lr_log[0], 2)
gen.num("lr_log_half", lr_log[0.5], 1)
gen.num("lr_log_1", lr_log[1], 1)
gen.num("crit", crit, 2)
gen.num("s_z0", s_z[0], 1)
gen.num("s_z1", s_z[1], 1)
gen.num("gm", gm, 1)
gen.num("g1", g1, 3)
gen.num("t1", t1, 2)
gen.num("onestep1", 1 - g1, 3)
gen.num("g0", g0, 3)
gen.num("t0", t0, 2)
gen.num("drop_lam", lam_drop[i_max], 3)
gen.num("drop_income", income[i_max], 0)
gen.num("drop_food", food[i_max], 0)
gen.num("coverage", coverage, 3)
gen.num("size", size, 3)
gen.write()

# ---- figure: the two profile log-likelihoods ---------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
for ax, (X, lam_hat, ci, l_max, title) in zip(axes, [
        (X_lin, lam_lin, ci_lin, lmax_lin, "(a) food linear in income"),
        (X_log, lam_log, ci_log, lmax_log, "(b) food linear in log income")]):
    lams = np.linspace(-0.6, 1.2, 181)
    ax.plot(lams, [profile(l, food, X) - l_max for l in lams], color=COLORS["accent"])
    ax.axhline(-crit / 2, color=COLORS["second"], linewidth=0.8, linestyle="--")
    for v in ci:
        ax.axvline(v, color=COLORS["muted"], linewidth=0.6)
    ax.set_ylim(-12, 1)
    ax.set_xlabel(r"$\lambda$")
    ax.set_title(title)
axes[0].set_ylabel(r"$\ell_p(\lambda)-\ell_p(\hat\lambda)$")
fig.tight_layout()
fig.savefig(figure_path("ch22", "box_cox_profile"))
