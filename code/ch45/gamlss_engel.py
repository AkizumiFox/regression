"""Chapter 45, Section 5: a normal location-scale GAMLSS for Engel's budget survey.

The mean and the log standard deviation each carry their own predictor, and the two
blocks are updated in turn by penalized Fisher scoring. The information matrix of the
normal distribution in this parameterization is block diagonal, so neither block
needs the other's derivatives. The two smoothing parameters are chosen by ten-fold
cross-validation of the predictive log score; choosing them by the in-sample
likelihood does not work, because that likelihood is unbounded, and the script checks
what happens when the scale is left nearly unpenalized.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.interpolate import BSpline

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data
DEGREE, N_INNER = 3, 8
GRID = 10.0 ** np.arange(-2, 5.1, 0.5)
CHOSEN = (10.0 ** 1.5, 10.0 ** 1.0)                # the pair the grid search picks out
WILD = 0.01                                        # a scale penalty small enough to run away

# <<fit>>
income = data["income"].to_numpy() / 1000.0        # thousands of francs
food = data["foodexp"].to_numpy() / 1000.0
n = len(food)

lo, hi = income.min(), income.max()
inner = np.quantile(income, np.linspace(0, 1, N_INNER + 2)[1:-1])
KNOTS = np.r_[[lo] * (DEGREE + 1), inner, [hi] * (DEGREE + 1)]


def basis(z):
    """The cubic B-spline basis of Section 43.3, with knots at income quantiles."""
    z = np.clip(np.atleast_1d(z), lo, hi)
    return np.asarray(BSpline.design_matrix(z, KNOTS, DEGREE).todense())


def difference_penalty(q, order=2):
    """K = D'D for the order-th difference matrix D: the penalty of Section 43.4."""
    D = np.diff(np.eye(q), order, axis=0)
    return D.T @ D


def gamlss_normal(Bm, Bs, y, lam_mu, lam_sig, iters=100, tol=1e-9):
    """Penalized Fisher scoring for y ~ N(mu, sigma^2), mu = Bm b_mu, log sigma = Bs b_sig.

    One block at a time: the mean step is a penalized weighted least squares fit with
    weights 1 / sigma^2, and the scale step is a penalized Fisher scoring step whose
    expected information is the constant 2.
    """
    Km, Ks = difference_penalty(Bm.shape[1]), difference_penalty(Bs.shape[1])
    b_mu = np.linalg.solve(Bm.T @ Bm + lam_mu * Km, Bm.T @ y)
    b_sig = np.linalg.solve(Bs.T @ Bs, Bs.T @ np.full(len(y), np.log(y.std())))
    for _ in range(iters):
        sigma = np.exp(np.clip(Bs @ b_sig, -8, 4))
        w = 1.0 / sigma ** 2                                  # working weights for the mean
        new_mu = np.linalg.solve(Bm.T @ (w[:, None] * Bm) + lam_mu * Km, Bm.T @ (w * y))
        u = (y - Bm @ new_mu) ** 2 / sigma ** 2 - 1.0         # score in the log-sigma direction
        new_sig = b_sig + np.linalg.solve(2 * Bs.T @ Bs + lam_sig * Ks,
                                          Bs.T @ u - lam_sig * Ks @ b_sig)
        step = max(np.max(np.abs(new_mu - b_mu)), np.max(np.abs(new_sig - b_sig)))
        b_mu, b_sig = new_mu, new_sig
        if step < tol:
            break
    return b_mu, b_sig


def neg_log_score(mu, sigma, y):
    """Minus the mean log density: the logarithmic score of Section 45.6."""
    return -np.mean(stats.norm.logpdf(y, mu, sigma))


folds = np.random.default_rng(7).permutation(n) % 10


def cross_validate(design, lam_mu, lam_sig):
    """Ten-fold cross-validated log score of one specification."""
    total = 0.0
    for k in range(10):
        train = folds != k
        Bm, Bs = design(income[train])
        b_mu, b_sig = gamlss_normal(Bm, Bs, food[train], lam_mu, lam_sig)
        Bm, Bs = design(income[~train])
        total += neg_log_score(Bm @ b_mu, np.exp(np.clip(Bs @ b_sig, -8, 4)), food[~train])
    return total / 10


spline_scale = lambda z: (basis(z), basis(z))
Bm, Bs = spline_scale(income)
b_mu, b_sig = gamlss_normal(Bm, Bs, food, *CHOSEN)
mu, sigma = Bm @ b_mu, np.exp(Bs @ b_sig)
cv_chosen = cross_validate(spline_scale, *CHOSEN)
print(f"smoothing parameters: {CHOSEN[0]:.3g} for the mean, {CHOSEN[1]:.3g} for the scale")
print(f"cross-validated log score {cv_chosen:.4f}")
# <</fit>>

# the search the cell above takes on trust: a full grid in both penalties
cv = {(a, b): cross_validate(spline_scale, a, b) for a in GRID for b in GRID}
best = min(cv, key=cv.get)
assert best == CHOSEN                            # the grid confirms the pair quoted in the text
assert abs(cv[best] - cv_chosen) < 1e-12

K = difference_penalty(Bm.shape[1])
w = 1.0 / sigma ** 2
edf_mu = np.trace(np.linalg.solve(Bm.T @ (w[:, None] * Bm) + CHOSEN[0] * K,
                                  Bm.T @ (w[:, None] * Bm)))
edf_sig = np.trace(np.linalg.solve(2 * Bs.T @ Bs + CHOSEN[1] * K, 2 * Bs.T @ Bs))
aic = -2 * np.sum(stats.norm.logpdf(food, mu, sigma)) + 2 * (edf_mu + edf_sig)

# the homoscedastic comparison: the same mean spline, one constant sigma
constant_scale = lambda z: (basis(z), np.ones((len(np.atleast_1d(z)), 1)))
cv_homo = {a: cross_validate(constant_scale, a, 0.0) for a in GRID}
best_homo = min(cv_homo, key=cv_homo.get)
Bm_h, Bs_h = constant_scale(income)
b_mu_h, b_sig_h = gamlss_normal(Bm_h, Bs_h, food, best_homo, 0.0)
sigma_h = float(np.exp(b_sig_h[0]))
w_h = np.full(n, 1 / sigma_h ** 2)
edf_mu_h = np.trace(np.linalg.solve(Bm_h.T @ (w_h[:, None] * Bm_h) + best_homo * K,
                                    Bm_h.T @ (w_h[:, None] * Bm_h)))
aic_homo = -2 * np.sum(stats.norm.logpdf(food, Bm_h @ b_mu_h, sigma_h)) + 2 * (edf_mu_h + 1)

assert cv_chosen < cv_homo[best_homo] - 0.03      # the scale model predicts better
assert sigma.max() > 5 * sigma.min()             # and the fitted scale is far from constant
assert np.all(np.diff(np.exp(basis(np.linspace(lo, np.quantile(income, 0.99), 50)) @ b_sig)) > 0)
assert abs(sigma_h - np.std(food - Bm_h @ b_mu_h)) < 0.03 * sigma_h

# <<unbounded>>
# barely penalize the scale: the in-sample likelihood climbs, the honest score collapses
b_mu_w, b_sig_w = gamlss_normal(Bm, Bs, food, CHOSEN[0], WILD)
wild_loglik = np.sum(stats.norm.logpdf(food, Bm @ b_mu_w, np.exp(np.clip(Bs @ b_sig_w, -8, 4))))
loglik = np.sum(stats.norm.logpdf(food, mu, sigma))
wild_cv = cross_validate(spline_scale, CHOSEN[0], WILD)
print(f"log-likelihood: {loglik:.1f} at the chosen penalty, {wild_loglik:.1f} at"
      f" lambda_sigma = {WILD:g}")
print(f"cross-validated log score: {cv_chosen:.4f} against {wild_cv:.4f}")
# <</unbounded>>

assert wild_loglik > loglik                      # the likelihood rewards the wilder scale
assert wild_cv > cv_chosen + 10                  # and honest prediction is destroyed

xs = np.linspace(lo, np.quantile(income, 0.99), 200)
Bx = basis(xs)
mu_curve, sig_curve = Bx @ b_mu, np.exp(Bx @ b_sig)
TAUS = (0.10, 0.50, 0.90)
implied = {tau: mu_curve + sig_curve * stats.norm.ppf(tau) for tau in TAUS}
assert np.all(np.diff(np.array([implied[t] for t in TAUS]), axis=0) > 0)   # never cross

at = lambda x: (float((basis(x) @ b_mu)[0]), float(np.exp(basis(x) @ b_sig)[0]))

gen = Generated("ch45", "gamlss_engel")
gen.int("n", n)
gen.int("q", Bm.shape[1])
gen.num("lam_mu", CHOSEN[0], 1)
gen.num("lam_sig", CHOSEN[1], 1)
gen.num("edf_mu", edf_mu, 2)
gen.num("edf_sig", edf_sig, 2)
gen.num("edf", edf_mu + edf_sig, 2)
gen.num("cv", cv_chosen, 4)
gen.num("cv_homo", cv_homo[best_homo], 4)
gen.num("aic", aic, 1)
gen.num("aic_homo", aic_homo, 1)
gen.num("loglik", loglik, 1)
gen.num("wild_loglik", wild_loglik, 1)
gen.num("wild_cv", wild_cv, 1)
gen.num("sigma_500", at(0.5)[1] * 1000, 1)
gen.num("sigma_2000", at(2.0)[1] * 1000, 1)
gen.num("sigma_ratio", at(2.0)[1] / at(0.5)[1], 2)
gen.num("sigma_homo", sigma_h * 1000, 1)
gen.num("mu_1000", at(1.0)[0] * 1000, 1)
gen.num("q90_1000", (at(1.0)[0] + at(1.0)[1] * stats.norm.ppf(0.9)) * 1000, 1)
gen.write()

# ---- figure ----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

ax = axes[0]
ax.scatter(income * 1000, food * 1000, s=8, color=COLORS["muted"], alpha=0.5, linewidths=0)
shades = plt.cm.viridis(np.linspace(0.15, 0.85, len(TAUS)))
for tau, colour in zip(TAUS, shades):
    ax.plot(xs * 1000, implied[tau] * 1000, color=colour, linewidth=1.1)
    ax.annotate(f"{tau:g}", (xs[-1] * 1000, implied[tau][-1] * 1000), color=colour,
                fontsize=6.5, xytext=(2, -2), textcoords="offset points")
ax.set_xlim(0, 3100)
ax.set_ylim(0, 2000)
ax.set_xlabel("household income")
ax.set_ylabel("food expenditure")
ax.set_title("(a) quantiles implied by the fit")

ax = axes[1]
ax.plot(xs * 1000, sig_curve * 1000, color=COLORS["accent"])
ax.axhline(sigma_h * 1000, color=COLORS["second"], linestyle="--", linewidth=1.0)
ax.annotate("one constant", (300, sigma_h * 1000), color=COLORS["second"],
            fontsize=7, xytext=(0, 5), textcoords="offset points")
ax.set_xlim(0, 3100)
ax.set_ylim(0, None)
ax.set_xlabel("household income")
ax.set_ylabel(r"fitted $\sigma$ (francs)")
ax.set_title("(b) the fitted scale")
fig.tight_layout()
fig.savefig(figure_path("ch45", "gamlss_engel"))
