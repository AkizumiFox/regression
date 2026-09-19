"""Chapter 14, Section 3: the sampling distribution of R^2 under multivariate normality.

Given the regressors, F = [(n-k-1)/k] R^2/(1-R^2) is F(k, n-k-1, lambda * V) with
lambda = rho^2/(1-rho^2) and V ~ chi^2(n-1) (the centred regressor sum of squares in the
direction of beta). Averaging over V gives the unconditional law of R^2. Checked against
simulation and against the Olkin-Pratt formula for E(R^2). Applied to the 2009 US state
data (statsmodels.datasets.statecrime, public domain): the overall F test and an interval
for the population squared multiple correlation.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "multiple_correlation", prefix="mc")


# <<mixture>>
def r2_cdf(r, rho2, n, k, m=4000):
    """P(R^2 <= r) for iid normal rows: average the conditional noncentral F over V."""
    df2 = n - k - 1
    f = df2 * r / (k * (1 - r))
    if rho2 == 0:
        return stats.f.cdf(f, k, df2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)      # equal-probability nodes
    return np.mean(stats.ncf.cdf(f, k, df2, rho2 / (1 - rho2) * V))


def r2_pdf(r, rho2, n, k, m=4000):
    """Density of R^2, by the same average of noncentral F densities."""
    df2 = n - k - 1
    f = df2 * r / (k * (1 - r))
    jac = df2 / (k * (1 - r) ** 2)                            # df/dr
    if rho2 == 0:
        return stats.f.pdf(f, k, df2) * jac
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)
    return np.mean(stats.ncf.pdf(f[:, None], k, df2, rho2 / (1 - rho2) * V[None, :]), axis=1) * jac
# <</mixture>>


def olkin_pratt_mean(rho2, n, k):
    """E(R^2) = 1 - (n-k-1)/(n-1) (1-rho^2) 2F1(1, 1; (n+1)/2; rho^2)."""
    return 1 - (n - k - 1) / (n - 1) * (1 - rho2) * special.hyp2f1(1, 1, (n + 1) / 2, rho2)


# ---- simulation check -------------------------------------------------------------
rng = np.random.default_rng(1431)
n, k, reps = 25, 4, 100_000
gen.int("n", n)
gen.int("k", k)
gen.int("reps", reps)
sims = {}
for rho2 in (0.0, 0.3, 0.6):
    lam = rho2 / (1 - rho2)
    beta = np.full(k, np.sqrt(lam / k))                      # beta^T beta = lambda, sigma = 1
    X = rng.normal(size=(reps, n, k))
    y = X @ beta + rng.normal(size=(reps, n))
    Xc = X - X.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    Sxy = np.einsum("rij,ri->rj", Xc, yc)
    Sxx = np.einsum("rij,rik->rjk", Xc, Xc)
    b = np.linalg.solve(Sxx, Sxy[..., None])[..., 0]
    R2 = np.einsum("rj,rj->r", b, Sxy) / np.einsum("ri,ri->r", yc, yc)
    sims[rho2] = R2
    mean_mix = olkin_pratt_mean(rho2, n, k)
    se = R2.std() / np.sqrt(reps)
    assert abs(R2.mean() - mean_mix) < 4 * se, (rho2, R2.mean(), mean_mix)
    # mixture cdf against the simulation at three points
    for q in (0.25, 0.5, 0.75):
        c = np.quantile(R2, q)
        assert abs(r2_cdf(c, rho2, n, k) - q) < 0.006, (rho2, q, r2_cdf(c, rho2, n, k))
    tag = f"{int(round(rho2 * 10))}"
    gen.num(f"mean_sim{tag}", R2.mean(), 4)
    gen.num(f"mean_exact{tag}", mean_mix, 4)
    adj = 1 - (1 - R2) * (n - 1) / (n - k - 1)
    gen.num(f"adj_sim{tag}", adj.mean(), 4)
    gen.num(f"sd_sim{tag}", R2.std(), 3)

# ---- power of the overall F test as a function of rho^2 ----------------------------
fcrit = stats.f.ppf(0.95, k, n - k - 1)
r2crit = k * fcrit / (k * fcrit + n - k - 1)
powers = {rho2: 1 - r2_cdf(r2crit, rho2, n, k) for rho2 in (0.1, 0.2, 0.3, 0.5)}
p3 = np.mean(sims[0.3] > r2crit)
assert abs(p3 - powers[0.3]) < 0.01
gen.num("r2crit", r2crit, 3)
for rho2, pw in powers.items():
    gen.num(f"power{int(round(rho2 * 10))}", pw, 3)

# ---- the state data ---------------------------------------------------------------------
# <<states>>
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
regs = ["hs_grad", "poverty", "single", "white", "urban"]
y = data["violent"].to_numpy()
X = np.column_stack([np.ones(len(y))] + [data[c].to_numpy() for c in regs])
fit = sm.OLS(y, X).fit()
n_s, k_s = len(y), len(regs)
R2_obs = fit.rsquared
print(f"R^2 = {R2_obs:.4f}, adjusted {fit.rsquared_adj:.4f}, F = {fit.fvalue:.2f},"
      f" p = {fit.f_pvalue:.2g}")

# 95% interval for rho^2: the values at which the observed R^2 is not extreme
def limit(prob):
    g = lambda rho2: r2_cdf(R2_obs, rho2, n_s, k_s) - prob
    return 0.0 if g(1e-12) < 0 else optimize.brentq(g, 1e-12, 0.999)

lower, upper = limit(0.975), limit(0.025)
print(f"95% interval for rho^2: ({lower:.3f}, {upper:.3f})")
# <</states>>
assert np.isclose(fit.fvalue, (n_s - k_s - 1) / k_s * R2_obs / (1 - R2_obs))
assert np.isclose(r2_cdf(R2_obs, lower, n_s, k_s), 0.975, atol=1e-6)
assert np.isclose(r2_cdf(R2_obs, upper, n_s, k_s), 0.025, atol=1e-6)
assert lower < fit.rsquared_adj < upper
gen.int("ns", n_s)
gen.int("ks", k_s)
gen.num("R2", R2_obs, 4)
gen.num("R2adj", fit.rsquared_adj, 4)
gen.num("F", fit.fvalue, 2)
gen.num("lower", lower, 3)
gen.num("upper", upper, 3)
gen.int("df2s", n_s - k_s - 1)
gen.num("pval", fit.f_pvalue, 1, sci=True)
gen.write()

# ---- figure ----------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.5))
grid = np.linspace(0.002, 0.95, 300)
for rho2, c in zip((0.0, 0.3, 0.6), (COLORS["accent"], COLORS["second"], COLORS["third"])):
    ax.hist(sims[rho2], bins=np.linspace(0, 1, 81), density=True, color=c,
            histtype="step", lw=0.7, alpha=0.7)
    ax.plot(grid, r2_pdf(grid, rho2, n, k, m=800), color=c, label=rf"$\rho^2={rho2:.1f}$")
    ax.axvline(rho2, color=c, lw=0.8, ls="--")
ax.set_xlim(0, 1)
ax.set_xlabel(r"$R^2$  ($n=25$, $k=4$)")
ax.set_ylabel("density")
ax.legend(frameon=False, loc="upper right")
fig.savefig(figure_path("ch14", "r2_distribution"))
