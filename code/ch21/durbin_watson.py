"""Chapter 21, Section 5: the Durbin-Watson test on a quarterly Okun's-law regression.

US macroeconomic data 1959Q1-2009Q3 (public domain, statsmodels.datasets.macrodata). The
response is the quarterly change in the unemployment rate (percentage points), the
regressor is real GDP growth (percent, annualized). The script computes the Durbin-Watson
statistic, its exact null distribution given X as a ratio of quadratic forms (by Imhof's
inversion formula, checked by simulation), its exact mean and variance, the Durbin-Watson
bounds, and the Breusch-Godfrey test.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.integrate import quad
from scipy.optimize import brentq
from statsmodels.stats.diagnostic import acorr_breusch_godfrey
from statsmodels.stats.stattools import durbin_watson

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
macro = sm.datasets.macrodata.load_pandas().data
# change in unemployment rate
du = np.diff(macro["unemp"].to_numpy())
# GDP growth, % a year
growth = 400 * np.diff(np.log(macro["realgdp"].to_numpy()))
year = (macro["year"] + (macro["quarter"] - 1) / 4).to_numpy()[1:]
n = len(du)
X = np.column_stack([np.ones(n), growth])
beta, *_ = np.linalg.lstsq(X, du, rcond=None)
e = du - X @ beta
# Durbin-Watson statistic
d = np.sum(np.diff(e) ** 2) / np.sum(e**2)
# lag-one residual autocorrelation
r1 = np.sum(e[1:] * e[:-1]) / np.sum(e**2)
print(f"n = {n}, slope {beta[1]:.4f}, d = {d:.3f}, r1 = {r1:.3f}, "
      f"2(1 - r1) = {2 * (1 - r1):.3f}")
# <</data>>
assert np.isclose(d, durbin_watson(e))


# <<exact>>
def dw_matrix(n):
    """A = D'D, where D takes first differences, so that d = e'Ae / e'e."""
    D = np.diff(np.eye(n), axis=0)
    return D.T @ D


def null_eigenvalues(X):
    """The n - p eigenvalues nu_k of the DW numerator on C(X)-perp."""
    n, p = X.shape
    Qfull, _ = np.linalg.qr(X, mode="complete")
    # orthonormal basis of C(X)-perp
    Qp = Qfull[:, p:]
    return np.linalg.eigvalsh(Qp.T @ dw_matrix(n) @ Qp)


def prob_negative(c):
    """Pr(sum_k c_k xi_k^2 < 0), xi_k iid standard normal (Imhof's formula)."""
    def integrand(u):
        theta = 0.5 * np.sum(np.arctan(c * u))
        rho = np.exp(0.25 * np.sum(np.log1p((c * u) ** 2)))
        return np.sin(theta) / (u * rho)
    return 0.5 - quad(integrand, 0, np.inf, limit=400)[0] / np.pi


nu = null_eigenvalues(X)
# Pr(d <= observed) under H0
p_value = prob_negative(nu - d)
m = len(nu)
mean_d = nu.mean()
var_d = 2 * (m * np.sum(nu**2) - np.sum(nu) ** 2) / (m**2 * (m + 2))
print(f"exact one-sided p-value {p_value:.2e}; "
      f"null mean {mean_d:.4f}, sd {np.sqrt(var_d):.4f}")
# <</exact>>

# simulation checks of the exact distribution and of its moments
rng = np.random.default_rng(2108)
xi2 = rng.chisquare(1, size=(100_000, m))
dsim = (xi2 @ nu) / xi2.sum(axis=1)
assert abs(dsim.mean() - mean_d) < 3e-3 and abs(dsim.var() / var_d - 1) < 0.03
q05 = brentq(lambda t: prob_negative(nu - t) - 0.05, 0.5, 2.0)
assert abs(np.mean(dsim <= q05) - 0.05) < 0.004
assert p_value < 1e-3

# eigenvalues of A and the bounds (Lemma 21.5.1, Theorem 21.5.3)
# lambda_1 = 0 < lambda_2 < ...
lam = 4 * np.sin(np.pi * np.arange(n) / (2 * n)) ** 2
assert np.allclose(np.sort(np.linalg.eigvalsh(dw_matrix(n))), lam)
p = X.shape[1]
assert np.all(nu >= lam[1:n - p + 1] - 1e-10) and np.all(nu <= lam[p:] + 1e-10)
crit_L = brentq(lambda t: prob_negative(lam[1:n - p + 1] - t) - 0.05, 0.5, 2.0)
crit_U = brentq(lambda t: prob_negative(lam[p:] - t) - 0.05, 0.5, 2.0)
assert crit_L < q05 < crit_U

# Breusch-Godfrey with four lags
# <<bg>>
lags = np.column_stack([np.concatenate([np.zeros(k), e[:-k]])
                        for k in range(1, 5)])
Z = np.column_stack([X, lags])
coef, *_ = np.linalg.lstsq(Z, e, rcond=None)
# n R^2
bg = n * (1 - np.sum((e - Z @ coef) ** 2) / np.sum((e - e.mean()) ** 2))
print(f"Breusch-Godfrey (4 lags) = {bg:.2f}, p = {stats.chi2.sf(bg, 4):.1e}")
# <</bg>>
bg_sm = acorr_breusch_godfrey(sm.OLS(du, X).fit(), nlags=4, result_object=False)
assert np.isclose(bg, bg_sm[0])

# a forty-quarter window, 1999Q2-2009Q1, where the bounds test is inconclusive
# <<window>>
w = slice(160, 200)
Xw, yw = X[w], du[w]
ew = yw - Xw @ np.linalg.lstsq(Xw, yw, rcond=None)[0]
d_w = np.sum(np.diff(ew) ** 2) / np.sum(ew**2)
nw_, pw_ = Xw.shape
lam_w = 4 * np.sin(np.pi * np.arange(nw_) / (2 * nw_)) ** 2
dL_w = brentq(lambda t: prob_negative(lam_w[1:nw_ - pw_ + 1] - t) - 0.05,
              0.5, 2.0)
dU_w = brentq(lambda t: prob_negative(lam_w[pw_:] - t) - 0.05, 0.5, 2.0)
p_w = prob_negative(null_eigenvalues(Xw) - d_w)
print(f"1999Q2-2009Q1: d = {d_w:.3f}, bounds ({dL_w:.3f}, {dU_w:.3f}), "
      f"exact p = {p_w:.3f}")
# <</window>>
assert dL_w < d_w < dU_w and p_w < 0.05

gen = Generated("ch21", "durbin_watson", prefix="dw")
gen.int("n", n)
gen.num("slope", beta[1], 4)
gen.num("intercept", beta[0], 3)
gen.num("d", d, 3)
gen.num("r1", r1, 3)
gen.num("approx", 2 * (1 - r1), 3)
gen.num("p", p_value, 1, sci=True)
gen.num("mean_d", mean_d, 4)
gen.num("sd_d", np.sqrt(var_d), 4)
gen.num("q05", q05, 3)
gen.num("crit_L", crit_L, 3)
gen.num("crit_U", crit_U, 3)
gen.num("bg", bg, 2)
gen.num("bg_p", stats.chi2.sf(bg, 4), 1, sci=True)
gen.num("d_w", d_w, 3)
gen.num("dL_w", dL_w, 3)
gen.num("dU_w", dU_w, 3)
gen.num("p_w", p_w, 3)
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.25), gridspec_kw={"width_ratios": [1.5, 1, 1.2]})
ax = axes[0]
ax.plot(year, e, color=COLORS["accent"], linewidth=0.7)
ax.axhline(0, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("year")
ax.set_ylabel("residual")
ax.set_title("(a) residuals in time order")
ax = axes[1]
ax.scatter(e[:-1], e[1:], s=5, color=COLORS["accent"], alpha=0.8, linewidths=0)
ax.set_xlabel(r"$\hat\varepsilon_{t-1}$")
ax.set_ylabel(r"$\hat\varepsilon_t$")
ax.set_title(f"(b) lag plot, $r_1={r1:.2f}$")
ax = axes[2]
grid = np.linspace(0.9, 4.4, 300)
for eig, label, color, ls in [(lam_w[1:nw_ - pw_ + 1], "$d_L$", COLORS["second"], "--"),
                              (null_eigenvalues(Xw), "$d$ given $X$", COLORS["accent"], "-"),
                              (lam_w[pw_:], "$d_U$", COLORS["third"], "--")]:
    cdf = np.array([prob_negative(eig - t) for t in grid])
    ax.plot(grid[1:], np.diff(cdf) / np.diff(grid), color=color, ls=ls, label=label)
ax.axvline(d_w, color=COLORS["ink"], linewidth=0.7)
ax.set_xlabel("value of the statistic")
ax.set_title("(c) null densities, $n=40$")
ax.legend(frameon=False, fontsize=6, loc="upper right", handlelength=1.2)
ax.set_xlim(0.9, 4.4)
fig.tight_layout()
fig.savefig(figure_path("ch21", "durbin_watson"))
