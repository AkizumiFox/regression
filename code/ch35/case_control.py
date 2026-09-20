"""Chapter 35, Section 5: retrospective sampling, and what survives it.

A simulated population of 200,000 with a rare binary outcome. Cases are all taken and controls
sampled at random; the logistic slopes are recovered and only the intercept moves, while the
probit slopes are not recovered. Matched pairs illustrate the incidental-parameter problem and
conditional logistic regression. Fixed seeds throughout.
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import minimize

from regbook import COLORS, Generated, figure_path, use_book_style

warnings.filterwarnings("ignore")

# <<population>>
import numpy as np
import statsmodels.api as sm

rng = np.random.default_rng(20)
N = 200_000
alpha, beta1, beta2 = -5.0, 1.2, 0.8
x1 = rng.normal(size=N)                                # a continuous exposure
x2 = (rng.random(N) < 0.3).astype(float)               # a binary covariate
Xpop = np.column_stack([np.ones(N), x1, x2])
ypop = (rng.random(N) < 1 / (1 + np.exp(-(Xpop @ [alpha, beta1, beta2])))).astype(float)

prospective = sm.GLM(ypop, Xpop, family=sm.families.Binomial()).fit()
cases = np.flatnonzero(ypop == 1)
controls = np.flatnonzero(ypop == 0)
n1 = len(cases)

print("prevalence", ypop.mean().round(4), " cases", n1)
print("whole population:", prospective.params.round(3))
# <</population>>

assert 0.005 < ypop.mean() < 0.03                      # a rare outcome
assert np.allclose(prospective.params, [alpha, beta1, beta2], atol=0.12)


# <<retrospective>>
def case_control(rng, ratio=4, link=None):
    """Take every case and `ratio` controls per case; fit the same model prospectively."""
    keep = np.r_[cases, rng.choice(controls, ratio * n1, replace=False)]
    family = sm.families.Binomial(link=link) if link else sm.families.Binomial()
    return sm.GLM(ypop[keep], Xpop[keep], family=family).fit().params


rng2 = np.random.default_rng(7)
draws = np.array([case_control(rng2) for _ in range(400)])
shift = np.log(1.0 / (4 * n1 / len(controls)))          # log(rho_1 / rho_0)

print("mean case-control estimate:", draws.mean(0).round(3))
print("predicted intercept:", round(prospective.params[0] + shift, 3))
print("slopes in the population:", prospective.params[1:].round(3))
# <</retrospective>>

assert abs(draws[:, 0].mean() - (prospective.params[0] + shift)) < 0.05
assert abs(draws[:, 1].mean() - beta1) < 0.05 and abs(draws[:, 2].mean() - beta2) < 0.06

# ---- the same experiment with a probit link ---------------------------------------------------
rng3 = np.random.default_rng(21)
yprob = (rng3.random(N) < stats.norm.cdf(Xpop @ [-2.4, 0.55, 0.35])).astype(float)
cases_p, controls_p = np.flatnonzero(yprob == 1), np.flatnonzero(yprob == 0)
probit_pop = sm.GLM(yprob, Xpop,
                    family=sm.families.Binomial(sm.families.links.Probit())).fit().params
assert np.allclose(probit_pop, [-2.4, 0.55, 0.35], atol=0.05)

rng4 = np.random.default_rng(8)
probit_draws = []
for _ in range(400):
    keep = np.r_[cases_p, rng4.choice(controls_p, 4 * len(cases_p), replace=False)]
    probit_draws.append(sm.GLM(yprob[keep], Xpop[keep],
                               family=sm.families.Binomial(sm.families.links.Probit())).fit().params)
probit_draws = np.array(probit_draws)
probit_ratio = probit_draws[:, 1].mean() / probit_pop[1]
assert probit_ratio > 1.3                               # the probit slope is not preserved

# the logistic slopes of the same retrospective samples are still right
logit_on_probit = []
for _ in range(100):
    keep = np.r_[cases_p, rng4.choice(controls_p, 4 * len(cases_p), replace=False)]
    logit_on_probit.append(sm.GLM(yprob[keep], Xpop[keep], family=sm.families.Binomial()).fit().params)
logit_pop = sm.GLM(yprob, Xpop, family=sm.families.Binomial()).fit().params
assert abs(np.mean(logit_on_probit, 0)[1] / logit_pop[1] - 1) < 0.06

# ---- noncollapsibility of the odds ratio -------------------------------------------------------
rng5 = np.random.default_rng(31)
ycom = (rng5.random(N) < 1 / (1 + np.exp(-(0.0 + 1.0 * x1 + 4.0 * x2)))).astype(float)
cond_or = sm.GLM(ycom, Xpop, family=sm.families.Binomial()).fit().params
marg_or = sm.GLM(ycom, Xpop[:, :2], family=sm.families.Binomial()).fit().params
assert marg_or[1] < 0.8 * cond_or[1]                    # the marginal slope is nearer zero
# x1 and x2 are independent by construction, so this is not confounding
assert abs(np.corrcoef(x1, x2)[0, 1]) < 0.01

# ---- matched pairs: conditional versus unconditional -------------------------------------------
# <<matched>>
from scipy.optimize import minimize


def matched_pairs(K=3000, beta=1.0, seed=11):
    """K strata, each with one case and one control, and a stratum-specific intercept."""
    rng = np.random.default_rng(seed)
    a = rng.normal(-1.0, 2.0, K)                        # the nuisance intercepts
    xc, xk = np.empty(K), np.empty(K)
    for k in range(K):
        while True:                                     # sample the stratum until it is discordant
            x = rng.normal(size=2)
            pr = 1 / (1 + np.exp(-(a[k] + beta * x)))
            yy = rng.random(2) < pr
            if yy[0] != yy[1]:
                xc[k], xk[k] = (x[0], x[1]) if yy[0] else (x[1], x[0])
                break
    return xc, xk                                       # exposure of the case, of the control


x_case, x_control = matched_pairs()
# the conditional likelihood for 1:1 matching depends on beta only through the difference
d = x_case - x_control
cond = minimize(lambda b: np.sum(np.log1p(np.exp(-b[0] * d))), [0.0]).x[0]
print("conditional estimate:", round(cond, 3))
# <</matched>>

# unconditional maximum likelihood with one intercept per stratum, computed without using the
# analytic profile of prp-bin-conditional(c): for each beta the K stratum intercepts are found
# by Newton's method on the stratum score, and the resulting profile is maximized over beta.
def stratum_alphas(b, x1, x0, steps=80):
    """Solve 1 - logistic(a + b x1) - logistic(a + b x0) = 0 for a, by bisection.

    The left-hand side is strictly decreasing in a, so bisection on a wide bracket is safe.
    """
    def score(a):
        return 1 - 1 / (1 + np.exp(-(a + b * x1))) - 1 / (1 + np.exp(-(a + b * x0)))

    lo = np.full(len(x1), -200.0)
    hi = np.full(len(x1), 200.0)
    for _ in range(steps):
        mid = (lo + hi) / 2
        pos = score(mid) > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return (lo + hi) / 2


def unconditional_profile(b, x1, x0):
    a = stratum_alphas(b, x1, x0)
    return np.sum((a + b * x1) - np.logaddexp(0, a + b * x1) - np.logaddexp(0, a + b * x0))


a_check = stratum_alphas(1.5, x_case, x_control)
p1c = 1 / (1 + np.exp(-(a_check + 1.5 * x_case)))
p0c = 1 / (1 + np.exp(-(a_check + 1.5 * x_control)))
assert np.max(np.abs(1 - p1c - p0c)) < 1e-10           # the intercepts really are maximizers
uncond = minimize(lambda b: -unconditional_profile(b[0], x_case, x_control), [0.5]).x[0]
assert abs(uncond - 2 * cond) < 1e-3                   # exactly twice, as the proposition says
# the observed information of the conditional log-likelihood at its maximum
q_c = 1 / (1 + np.exp(-cond * d))
cond_se = 1 / np.sqrt(np.sum(d ** 2 * q_c * (1 - q_c)))
assert abs(cond - 1.0) < 2.5 * cond_se

gen = Generated("ch35", "cc")
gen.int("N", N)
gen.int("cases", n1)
gen.num("prevalence", ypop.mean(), 4)
gen.num("pop_alpha", prospective.params[0], 3)
gen.num("pop_b1", prospective.params[1], 3)
gen.num("pop_b2", prospective.params[2], 3)
gen.num("cc_alpha", draws[:, 0].mean(), 3)
gen.num("cc_b1", draws[:, 1].mean(), 3)
gen.num("cc_b2", draws[:, 2].mean(), 3)
gen.num("shift", shift, 3)
gen.num("predicted_alpha", prospective.params[0] + shift, 3)
gen.num("cc_sd1", draws[:, 1].std(), 3)
gen.num("probit_pop_b1", probit_pop[1], 3)
gen.num("probit_cc_b1", probit_draws[:, 1].mean(), 3)
gen.num("probit_ratio", probit_ratio, 2)
gen.num("cond_or", cond_or[1], 3)
gen.num("marg_or", marg_or[1], 3)
gen.num("common_prev", ycom.mean(), 3)
gen.int("K", len(d))
gen.num("cond", cond, 3)
gen.num("uncond", uncond, 3)
gen.num("cond_se", cond_se, 3)
gen.write()

# ---- two sampling distributions ------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))

ax = axes[0]
ax.hist(draws[:, 1], bins=30, color=COLORS["accent"], alpha=0.8)
ax.axvline(prospective.params[1], color=COLORS["second"], linestyle="--",
           label="whole population")
ax.set_xlabel(r"retrospective $\hat\beta_1$, logit link")
ax.set_ylabel("replications")
ax.legend(frameon=False, fontsize=7)
ax.set_title("(a) odds ratios survive")

ax = axes[1]
ax.hist(probit_draws[:, 1], bins=30, color=COLORS["accent"], alpha=0.8)
ax.axvline(probit_pop[1], color=COLORS["second"], linestyle="--", label="whole population")
ax.set_xlabel(r"retrospective $\hat\beta_1$, probit link")
ax.set_ylabel("replications")
ax.legend(frameon=False, fontsize=7)
ax.set_title("(b) probit slopes do not")
fig.tight_layout()
fig.savefig(figure_path("ch35", "case_control"))
