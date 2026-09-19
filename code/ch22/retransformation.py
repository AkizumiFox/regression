"""Chapter 22, Section 5: back to the original scale after a log-log fit to Engel's data.

Medians, means (lognormal correction and Duan's smearing estimate), prediction intervals and a
confidence interval for the mean, at two incomes. A simulation with skewed errors shows the
normal-theory correction failing where smearing does not.
Public-domain data shipped with statsmodels (statsmodels.datasets.engel).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<back>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), np.log(income)])
fit = sm.OLS(np.log(food), X).fit()
b, s2, resid = fit.params, fit.scale, fit.resid
XtX_inv = np.linalg.inv(X.T @ X)
t = stats.t.ppf(0.975, n - 2)

smear = np.mean(np.exp(resid))                    # Duan's smearing factor
lognormal = np.exp(s2 / 2)                        # the normal-theory factor
print(f"smearing factor {smear:.4f}, lognormal factor {lognormal:.4f}")

for x0 in (500, 2000):
    v0 = np.array([1, np.log(x0)])
    eta, h0 = v0 @ b, v0 @ XtX_inv @ v0
    median = np.exp(eta)
    ci_med = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * h0))          # for the median, exact
    pi = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * (1 + h0)))        # for a new household, exact
    print(f"income {x0}: median {median:.1f}, CI [{ci_med[0]:.1f}, {ci_med[1]:.1f}];"
          f" mean {median * smear:.1f} (smearing), {median * lognormal:.1f} (lognormal);"
          f" 95% prediction interval [{pi[0]:.0f}, {pi[1]:.0f}]")
# <</back>>

gen = Generated("ch22", "retransformation", prefix="rt")
gen.num("smear", smear, 4)
gen.num("lognormal", lognormal, 4)
gen.num("s2", s2, 5)
for x0 in (500, 2000):
    v0 = np.array([1, np.log(x0)])
    eta, h0 = v0 @ b, v0 @ XtX_inv @ v0
    median = np.exp(eta)
    ci_med = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * h0))
    pi = np.exp(eta + np.array([-1, 1]) * t * np.sqrt(s2 * (1 + h0)))
    # approximate interval for the lognormal mean exp(eta + sigma^2/2): theta_hat +- z se
    theta = eta + s2 / 2
    se_theta = np.sqrt(s2 * h0 + s2**2 / (2 * (n - 2)))
    ci_mean = np.exp(theta + np.array([-1, 1]) * stats.norm.ppf(0.975) * se_theta)
    assert pi[0] < ci_med[0] < median < ci_med[1] < pi[1]
    assert (pi[1] - median) > (median - pi[0])            # the back-transformed interval is skewed
    k = f"{x0}"
    gen.num(f"eta{k}", eta, 4)
    gen.num(f"h{k}", h0, 5)
    gen.num(f"median{k}", median, 1)
    gen.num(f"cimed_lo{k}", ci_med[0], 1)
    gen.num(f"cimed_hi{k}", ci_med[1], 1)
    gen.num(f"mean_smear{k}", median * smear, 1)
    gen.num(f"mean_ln{k}", median * lognormal, 1)
    gen.num(f"pi_lo{k}", pi[0], 0)
    gen.num(f"pi_hi{k}", pi[1], 0)
    gen.num(f"cimean_lo{k}", ci_mean[0], 1)
    gen.num(f"cimean_hi{k}", ci_mean[1], 1)
    gen.num(f"se_theta{k}", se_theta, 4)
assert 1 < smear < 1.02 and 1 < lognormal < 1.02
gen.num("t", t, 3)
# the share of households below the fitted median curve is about one half
below = np.mean(np.log(food) < fit.fittedvalues)
gen.num("below", below, 3)
assert 0.4 < below < 0.6
# elasticities at the two incomes: the log-log slope, and the linear fit's elasticity at the fitted mean
lin = sm.OLS(food, np.column_stack([np.ones(n), income])).fit()
for x0 in (500, 2000):
    gen.num(f"elast_lin{x0}", lin.params[1] * x0 / (lin.params[0] + lin.params[1] * x0), 2)


# ---- simulation: skewed errors -----------------------------------------------------------------
def simulate(reps, n_sim, seed):
    """Log-linear model with right-skewed errors e = 0.6 (E - 1), E ~ Exp(1): estimate / true mean."""
    rng = np.random.default_rng(seed)
    xs = np.linspace(0, 2, n_sim)
    Xs = np.column_stack([np.ones(n_sim), xs])
    beta, sig = np.array([1.0, 0.5]), 0.6
    x0 = np.array([1.0, 1.0])
    true_mean = np.exp(x0 @ beta) * np.exp(-sig) / (1 - sig)       # E exp(e) = e^(-sig)/(1 - sig)
    out = np.empty((reps, 3))
    for r in range(reps):
        e = sig * (rng.exponential(size=n_sim) - 1)
        z = Xs @ beta + e
        coef, *_ = np.linalg.lstsq(Xs, z, rcond=None)
        res = z - Xs @ coef
        s2_sim = res @ res / (n_sim - 2)
        naive = np.exp(x0 @ coef)
        out[r] = [naive, naive * np.exp(s2_sim / 2), naive * np.mean(np.exp(res))]
    return out.mean(axis=0) / true_mean


# <<skewed>>
ratios = simulate(reps=300, n_sim=100, seed=2205)
print("mean of estimate / true mean:  naive %.3f   lognormal %.3f   smearing %.3f" % tuple(ratios))
# <</skewed>>
ratios = simulate(reps=4000, n_sim=100, seed=2206)
assert ratios[0] < 0.8 and ratios[1] < 0.93 and abs(ratios[2] - 1) < 0.02
gen.num("r_naive", ratios[0], 3)
gen.num("r_ln", ratios[1], 3)
gen.num("r_smear", ratios[2], 3)
gen.num("true_factor", np.exp(-0.6) / 0.4, 3)
gen.num("ln_factor", np.exp(0.18), 3)
gen.write()

# ---- figure: the fit on the original scale ----------------------------------------------------
use_book_style()
grid = np.geomspace(income.min(), income.max(), 200)
V = np.column_stack([np.ones_like(grid), np.log(grid)])
eta_g = V @ b
h_g = np.einsum("ij,jk,ik->i", V, XtX_inv, V)
fig, ax = plt.subplots(figsize=(4.4, 2.9))
ax.fill_between(grid, np.exp(eta_g - t * np.sqrt(s2 * (1 + h_g))), np.exp(eta_g + t * np.sqrt(s2 * (1 + h_g))),
                color=COLORS["grid"], alpha=0.6, linewidth=0, label="95% prediction interval")
ax.scatter(income, food, s=6, color=COLORS["accent"], alpha=0.6, linewidths=0)
ax.plot(grid, np.exp(eta_g), color=COLORS["second"], label="median")
ax.plot(grid, np.exp(eta_g) * smear, color=COLORS["third"], linestyle="--", label="mean (smearing)")
ax.set_xlabel("income")
ax.set_ylabel("food expenditure")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch22", "retransformation"))
