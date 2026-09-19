"""Chapter 14, Section 5: sample partial correlations.

(a) With q conditioning variables, the sample partial correlation from n cases has the law of
    an ordinary sample correlation from n - q cases with the same (partial) correlation.
    Checked by simulation, with skewed conditioning variables to show that only the
    conditional normality of the pair matters.
(b) The 2009 US state data (statsmodels.datasets.statecrime, public domain): the partial
    correlation of violent crime and urbanization given single parenthood and poverty, the
    identity t^2 = (n-p) r^2/(1-r^2), and intervals.
"""
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

from regbook import Generated

gen = Generated("ch14", "partial_correlation", prefix="pc")


def r_cdf(c, rho, n, m=4000):
    """Exact P(r <= c) for an ordinary sample correlation of n bivariate normal pairs."""
    if c < 0:
        return 1 - r_cdf(-c, -rho, n, m)
    t = c * np.sqrt(n - 2) / np.sqrt(1 - c ** 2)
    if rho == 0:
        return stats.t.cdf(t, n - 2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)
    return np.mean(stats.nct.cdf(t, n - 2, rho / np.sqrt(1 - rho ** 2) * np.sqrt(V)))


# ---- (a) simulation ------------------------------------------------------------------------
# <<simulate>>
rng = np.random.default_rng(1451)
n, q, rho, reps = 12, 2, 0.6, 100_000
B = np.array([[1.0, -2.0], [0.5, 3.0]])              # dependence of the pair on Z
C = np.array([[1.0, rho], [rho, 1.0]])               # conditional covariance of the pair
Lc = np.linalg.cholesky(C)

partial = np.empty(reps)
for s in range(reps):
    Zc = rng.exponential(size=(n, q)) ** 2           # skewed conditioning variables
    pair = Zc @ B.T + rng.normal(size=(n, 2)) @ Lc.T
    D = np.column_stack([np.ones(n), Zc])
    coef, *_ = np.linalg.lstsq(D, pair, rcond=None)
    e = pair - D @ coef                              # residuals of both variables
    partial[s] = e[:, 0] @ e[:, 1] / np.sqrt((e[:, 0] @ e[:, 0]) * (e[:, 1] @ e[:, 1]))

# ordinary correlations from n - q pairs with correlation rho
x = rng.normal(size=(reps, n - q))
y = rho * x + np.sqrt(1 - rho ** 2) * rng.normal(size=(reps, n - q))
ordinary = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(x, y)])
print("quantiles, partial :", np.quantile(partial, [0.1, 0.5, 0.9]).round(3))
print("quantiles, ordinary:", np.quantile(ordinary, [0.1, 0.5, 0.9]).round(3))
print("two-sample KS p-value:", round(stats.ks_2samp(partial, ordinary).pvalue, 3))
# <</simulate>>
ks = stats.ks_2samp(partial, ordinary).pvalue
assert ks > 1e-3
for qq in (0.1, 0.5, 0.9):
    assert abs(r_cdf(np.quantile(partial, qq), rho, n - q) - qq) < 0.005
qp = np.quantile(partial, [0.1, 0.5, 0.9])
qo = np.quantile(ordinary, [0.1, 0.5, 0.9])
gen.int("n", n)
gen.int("q", q)
gen.num("rho", rho, 1)
gen.int("reps", reps)
for name, v in zip(("10", "50", "90"), qp):
    gen.num(f"qp{name}", v, 3)
for name, v in zip(("10", "50", "90"), qo):
    gen.num(f"qo{name}", v, 3)
gen.num("ks", ks, 2)

# ---- (b) the state data -----------------------------------------------------------------------
# <<states>>
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["violent"].to_numpy()
x = data["urban"].to_numpy()
Z = np.column_stack([np.ones(len(y)), data["single"], data["poverty"]])
n_s, q_s = len(y), Z.shape[1] - 1                    # q = 2 conditioning variables

def resid(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef

ey, ex = resid(y, Z), resid(x, Z)
r_p = ey @ ex / np.sqrt((ey @ ey) * (ex @ ex))       # sample partial correlation
fit = sm.OLS(y, np.column_stack([Z, x])).fit()       # p = 4 columns
t = fit.tvalues[-1]
df = n_s - 4
print(f"partial r = {r_p:.4f},  t = {t:.3f},  sqrt(df) r/sqrt(1-r^2) = "
      f"{np.sqrt(df) * r_p / np.sqrt(1 - r_p**2):.3f}")
zq = stats.norm.ppf(0.975)
lo_z, hi_z = np.tanh(np.arctanh(r_p) + np.array([-1, 1]) * zq / np.sqrt(n_s - q_s - 3))
lo_x = optimize.brentq(lambda rh: r_cdf(r_p, rh, n_s - q_s) - 0.975, -0.99, 0.99)
hi_x = optimize.brentq(lambda rh: r_cdf(r_p, rh, n_s - q_s) - 0.025, -0.99, 0.99)
print(f"Fisher z interval ({lo_z:.3f}, {hi_z:.3f}); exact ({lo_x:.3f}, {hi_x:.3f})")
# <</states>>
assert np.isclose(t, np.sqrt(df) * r_p / np.sqrt(1 - r_p ** 2))
assert lo_x < r_p < hi_x
r_marg = np.corrcoef(y, x)[0, 1]
gen.int("ns", n_s)
gen.int("df", df)
gen.num("r_p", r_p, 4)
gen.num("r_marg", r_marg, 3)
gen.num("t", t, 3)
gen.num("pval", fit.pvalues[-1], 4)
gen.num("coef", fit.params[-1], 3)
gen.num("lo_z", lo_z, 3)
gen.num("hi_z", hi_z, 3)
gen.num("lo_x", lo_x, 3)
gen.num("hi_x", hi_x, 3)
gen.write()
