"""Chapter 12, Section 1: t intervals for estimable functions, and an interval for sigma^2.

Murder rate of the 50 US states (2009) on poverty, single-parent households and urbanization,
as in Chapters 6 and 9 (statsmodels.datasets.statecrime, public domain). The simulation checks
the coverage of the t interval and of the chi-squared interval for sigma^2, with normal and with
heavy-tailed errors of the same variance.
"""
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

# <<fit>>
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
r = np.linalg.matrix_rank(X)
G = np.linalg.pinv(X.T @ X)                     # any generalized inverse will do
beta_hat = G @ X.T @ y
s2 = np.sum((y - X @ beta_hat) ** 2) / (n - r)
q = stats.t.ppf(0.975, n - r)                   # upper 2.5% point of t(n - r)


def t_interval(lam):
    """95% interval for the estimable function lam^T beta."""
    est = lam @ beta_hat
    se = np.sqrt(s2 * lam @ G @ lam)
    return est, se, est - q * se, est + q * se


for j, name in enumerate(["intercept", "poverty", "single", "urban"]):
    est, se, lo, hi = t_interval(np.eye(p)[j])
    print(f"{name:9s} {est:8.4f}  se {se:.4f}  95% interval ({lo:.3f}, {hi:.3f})")
# <</fit>>

# <<mean>>
x0 = np.array([1.0, 18.0, 28.0, 70.0])          # a state with high poverty and single parenthood
est0, se0, lo0, hi0 = t_interval(x0)
print(f"mean murder rate at x0: {est0:.3f}, se {se0:.3f}, interval ({lo0:.3f}, {hi0:.3f})")

sse = s2 * (n - r)
lo_s2 = sse / stats.chi2.ppf(0.975, n - r)      # interval for sigma^2 from SSE/sigma^2 ~ chi^2(n-r)
hi_s2 = sse / stats.chi2.ppf(0.025, n - r)
print(f"s^2 = {s2:.3f}; 95% interval for sigma^2 ({lo_s2:.3f}, {hi_s2:.3f})")
# <</mean>>

# ---- checks ---------------------------------------------------------------------------------
assert r == p == 4 and n == 50
fit = sm.OLS(y, X).fit()
ci = fit.conf_int(0.05)
for j in range(p):
    _, se, lo, hi = t_interval(np.eye(p)[j])
    assert np.isclose(se, fit.bse[j]) and np.isclose(lo, ci[j, 0]) and np.isclose(hi, ci[j, 1])
# duality with the t test: d is inside the interval iff |T_d| <= q
est1, se1, lo1, hi1 = t_interval(np.eye(p)[1])
for d in np.linspace(lo1 - 0.2, hi1 + 0.2, 41):
    assert (abs(est1 - d) / se1 <= q) == (lo1 <= d <= hi1)
# invariance: a redundant column and a different generalized inverse give the same interval
Xr = np.column_stack([X, X[:, 1] + X[:, 2]])    # rank still 4, p = 5
A = Xr.T @ Xr
G2 = np.zeros((5, 5))
G2[:4, :4] = np.linalg.inv(A[:4, :4])           # a non-Moore-Penrose generalized inverse
assert np.allclose(A @ G2 @ A, A) and not np.allclose(G2, np.linalg.pinv(A))
lam_r = np.array([1.0, 18.0, 28.0, 70.0, 46.0])  # the same mean response, written in the new coordinates
assert np.isclose(lam_r @ G2 @ Xr.T @ y, est0) and np.isclose(lam_r @ G2 @ lam_r, x0 @ G @ x0)
# FWL form of the standard error of a coefficient (Chapter 6)
others = X[:, [0, 2, 3]]
xt = X[:, 1] - others @ np.linalg.lstsq(others, X[:, 1], rcond=None)[0]
assert np.isclose(np.sqrt(s2 / (xt @ xt)), se1)

# ---- coverage simulation: t interval for the poverty coefficient, and the sigma^2 interval ------
rng = np.random.default_rng(1201)
reps = 40_000
sigma = np.sqrt(s2)                             # use the fitted model as the truth
H = G @ X.T                                     # beta_hat = H y
M = X @ H
lam = np.eye(p)[1]


def coverage(E):
    Y = X @ beta_hat + E
    B = Y @ H.T
    S2 = np.sum((Y - Y @ M) ** 2, axis=1) / (n - r)
    se = np.sqrt(S2 * (lam @ G @ lam))
    cov_t = np.mean(np.abs(B @ lam - lam @ beta_hat) <= q * se)
    sse_sim = S2 * (n - r)
    cov_s2 = np.mean((sse_sim / stats.chi2.ppf(0.975, n - r) <= sigma**2)
                     & (sigma**2 <= sse_sim / stats.chi2.ppf(0.025, n - r)))
    return cov_t, cov_s2


cov_t_norm, cov_s2_norm = coverage(sigma * rng.normal(size=(reps, n)))
nu = 5                                          # t(5) errors rescaled to variance sigma^2
cov_t_heavy, cov_s2_heavy = coverage(sigma * rng.standard_t(nu, size=(reps, n)) / np.sqrt(nu / (nu - 2)))
print(f"coverage, normal errors: t interval {cov_t_norm:.4f}, sigma^2 interval {cov_s2_norm:.4f}")
print(f"coverage, t(5) errors:   t interval {cov_t_heavy:.4f}, sigma^2 interval {cov_s2_heavy:.4f}")
se_mc = np.sqrt(0.95 * 0.05 / reps)
assert abs(cov_t_norm - 0.95) < 4 * se_mc and abs(cov_s2_norm - 0.95) < 4 * se_mc
assert abs(cov_t_heavy - 0.95) < 0.01          # the t interval is robust
assert cov_s2_heavy < 0.85                      # the sigma^2 interval is not

# the factor t * kappa that replaces z when sigma is estimated (expected interval length)
def kappa(nu):
    from scipy.special import gammaln
    return np.sqrt(2 / nu) * np.exp(gammaln((nu + 1) / 2) - gammaln(nu / 2))


rs = np.random.default_rng(12011)
V = rs.chisquare(5, size=400_000)
assert abs(np.mean(np.sqrt(V / 5)) - kappa(5)) < 0.003
assert np.isclose(kappa(1), np.sqrt(2 / np.pi)) and np.isclose(kappa(2), np.sqrt(np.pi) / 2)
tk = {nu: stats.t.ppf(0.975, nu) * kappa(nu) for nu in (5, 20)}
assert tk[5] > tk[20] > stats.norm.ppf(0.975)

gen = Generated("ch12", "estimable_intervals", prefix="est")
gen.num("tk5", tk[5], 2)
gen.num("tk20", tk[20], 2)
gen.num("z", stats.norm.ppf(0.975), 2)
gen.num("kappa1", kappa(1), 3)
gen.num("kappa2", kappa(2), 3)
gen.int("n", n)
gen.int("df", n - r)
gen.num("q", q, 3)
gen.num("s2", s2, 3)
gen.num("s", np.sqrt(s2), 3)
for j, name in enumerate(["int", "pov", "sing", "urb"]):
    est, se, lo, hi = t_interval(np.eye(p)[j])
    gen.num(f"{name}:est", est, 4)
    gen.num(f"{name}:se", se, 4)
    gen.num(f"{name}:lo", lo, 3)
    gen.num(f"{name}:hi", hi, 3)
gen.num("m:est", est0, 3)
gen.num("m:se", se0, 3)
gen.num("m:lo", lo0, 3)
gen.num("m:hi", hi0, 3)
gen.num("m:h0", x0 @ G @ x0, 4)
gen.num("s2:lo", lo_s2, 3)
gen.num("s2:hi", hi_s2, 3)
gen.num("chi:lo", stats.chi2.ppf(0.025, n - r), 3)
gen.num("chi:hi", stats.chi2.ppf(0.975, n - r), 3)
gen.num("cov:tnorm", cov_t_norm, 4)
gen.num("cov:s2norm", cov_s2_norm, 4)
gen.num("cov:theavy", cov_t_heavy, 4)
gen.num("cov:s2heavy", cov_s2_heavy, 4)
gen.int("reps", reps)
gen.write()
