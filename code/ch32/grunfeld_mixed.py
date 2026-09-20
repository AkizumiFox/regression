"""Chapter 32, Sections 2 and 5: mixed models for Grunfeld's investment panel.

Public-domain data shipped with statsmodels (statsmodels.datasets.grunfeld): gross
investment, market value and capital stock for 11 US firms over 20 years. Section 6.6
fitted this panel with firm indicators (the within estimator) and with a single
intercept. Here the firm is a random effect, and the script compares the four analyses,
fits the model by a REML routine written from the definition, checks it against
statsmodels' MixedLM, and shows how flat the likelihood becomes when firm-specific
slopes are added.
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import optimize

from regbook import COLORS, Generated, figure_path, use_book_style

warnings.simplefilter("ignore")

# <<data>>
import numpy as np
import pandas as pd
import statsmodels.api as sm

df = sm.datasets.grunfeld.load_pandas().data.sort_values(["firm", "year"])
firms = df["firm"].unique()
y = df["invest"].to_numpy()
X = np.column_stack([np.ones(len(y)), df["value"], df["capital"]])
blocks = [np.where(df["firm"].to_numpy() == f)[0] for f in firms]
n, p, g = len(y), X.shape[1], len(firms)
print(f"{n} observations, {g} firms, {n // g} years each")
# <</data>>

# <<reml>>
def restricted_loglik(theta, Zcols):
    """Restricted log-likelihood; theta holds log sigma^2 and the log variances of Zcols."""
    s2, dvar = np.exp(theta[0]), np.exp(theta[1:])
    logdet, W, XtVy, yVy = 0.0, np.zeros((p, p)), np.zeros(p), 0.0
    for ix in blocks:
        Zi = Zcols[ix]
        Vi = (Zi * dvar) @ Zi.T + s2 * np.eye(len(ix))
        Vinv = np.linalg.inv(Vi)
        logdet += np.linalg.slogdet(Vi)[1]
        W += X[ix].T @ Vinv @ X[ix]
        XtVy += X[ix].T @ Vinv @ y[ix]
        yVy += y[ix] @ Vinv @ y[ix]
    beta = np.linalg.solve(W, XtVy)
    return -0.5 * (logdet + yVy - XtVy @ beta + np.linalg.slogdet(W)[1]
                   - np.linalg.slogdet(X.T @ X)[1] + (n - p) * np.log(2 * np.pi)), beta, W


Z_int = np.ones((n, 1))                                   # random intercept only
fit = optimize.minimize(lambda t: -restricted_loglik(t, Z_int)[0],
                        np.log([2000.0, 5000.0]), method="Nelder-Mead",
                        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 5000})
ll_int, beta_int, W_int = restricted_loglik(fit.x, Z_int)
s2_hat, s2_firm = np.exp(fit.x)
se_int = np.sqrt(np.diag(np.linalg.inv(W_int)))

print(f"sigma^2 {s2_hat:.1f}   firm variance {s2_firm:.1f}")
print(f"value   {beta_int[1]:.4f} ({se_int[1]:.4f})")
print(f"capital {beta_int[2]:.4f} ({se_int[2]:.4f})")
# <</reml>>

# statsmodels finds the same fit
sm_fit = sm.MixedLM.from_formula("invest ~ value + capital", groups="firm", data=df).fit()
assert abs(sm_fit.scale / s2_hat - 1) < 1e-3, (sm_fit.scale, s2_hat)
assert abs(sm_fit.cov_re.iloc[0, 0] / s2_firm - 1) < 1e-3
assert np.allclose(sm_fit.params.to_numpy()[:p], beta_int, rtol=1e-3)
assert np.allclose(sm_fit.bse.to_numpy()[:p], se_int, rtol=1e-2)

# ---- the four analyses of the same panel ------------------------------------------
pooled, *_ = np.linalg.lstsq(X, y, rcond=None)
D = np.zeros((n, g))
for j, ix in enumerate(blocks):
    D[ix, j] = 1.0
within, *_ = np.linalg.lstsq(np.column_stack([X[:, 1:], D]), y, rcond=None)
icc = s2_firm / (s2_firm + s2_hat)
shrink = (n // g) * s2_firm / ((n // g) * s2_firm + s2_hat)
# with 20 years per firm the mixed estimate is far closer to the within estimate
assert abs(beta_int[1] - within[0]) < 0.2 * abs(beta_int[1] - pooled[1])
assert abs(beta_int[2] - within[1]) < 0.1 * abs(beta_int[2] - pooled[2])

# ---- adding random slopes: a much flatter likelihood -------------------------------
Z_slope = np.column_stack([np.ones(n), df["value"] / 1000])
best, all_ll = None, []
for seed in range(30):
    start = np.log([2000.0, 5000.0, 20000.0]) + np.random.default_rng(seed).normal(0, 0.8, 3)
    r = optimize.minimize(lambda t: -restricted_loglik(t, Z_slope)[0], start,
                          method="Nelder-Mead",
                          options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 20000})
    all_ll.append(-r.fun)
    if best is None or r.fun < best.fun:
        best = r
# every one of the thirty random starts reaches the same maximum
assert np.ptp(all_ll) < 1e-3, np.ptp(all_ll)
ll_slope, beta_slope, W_slope = restricted_loglik(best.x, Z_slope)
s2_slope, dvar = np.exp(best.x[0]), np.exp(best.x[1:])
se_slope = np.sqrt(np.diag(np.linalg.inv(W_slope)))
assert ll_slope > ll_int
lrt = 2 * (ll_slope - ll_int)

# On the log scale a start with a negligible slope variance stalls on the boundary: the
# search stops at the random-intercept fit rather than at the maximum.
stall = optimize.minimize(lambda t: -restricted_loglik(t, Z_slope)[0],
                          np.log([2000.0, 5000.0, 1e-6]), method="Nelder-Mead",
                          options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 20000})
ll_stall = -stall.fun
assert abs(ll_stall - ll_int) < 0.01, (ll_stall, ll_int)   # it is the random-intercept fit
assert ll_slope - ll_stall > 10.0

# What the two models estimate. Fitting firm-specific intercepts and firm-specific value
# slopes with a common capital coefficient shows how widely the firms' slopes are spread.
Xsep = np.column_stack([D, D * df["value"].to_numpy()[:, None], X[:, 2]])
sep, *_ = np.linalg.lstsq(Xsep, y, rcond=None)
firm_slopes = sep[g:2 * g]
assert np.std(firm_slopes, ddof=1) > 10 * abs(np.mean(firm_slopes))

# REML AIC and BIC: legitimate here because the fixed effects are the same
aic_int, aic_slope = -2 * ll_int + 2 * 2, -2 * ll_slope + 2 * 3
bic_int = -2 * ll_int + 2 * np.log(n - p)
bic_slope = -2 * ll_slope + 3 * np.log(n - p)
assert aic_slope < aic_int and bic_slope < bic_int

# the firm intercepts: dummy-variable estimates against the mixed model's predictions
alpha_ols = within[2:] - within[2:].mean()
Vinv_blocks = []
u_int = np.empty(g)
for j, ix in enumerate(blocks):
    Vi = s2_firm * np.ones((len(ix), len(ix))) + s2_hat * np.eye(len(ix))
    u_int[j] = s2_firm * np.ones(len(ix)) @ np.linalg.solve(Vi, y[ix] - X[ix] @ beta_int)
assert np.corrcoef(alpha_ols, u_int)[0, 1] > 0.99      # with 20 years per firm, no shrinkage left

gen = Generated("ch32", "grunfeld")
gen.int("n", n)
gen.int("g", g)
gen.int("years", n // g)
gen.num("s2", s2_hat, 1)
gen.num("s2firm", s2_firm, 1)
gen.num("sdfirm", np.sqrt(s2_firm), 1)
gen.num("icc", icc, 3)
gen.num("shrink", shrink, 4)
gen.num("valmix", beta_int[1], 4)
gen.num("capmix", beta_int[2], 4)
gen.num("valse", se_int[1], 4)
gen.num("capse", se_int[2], 4)
gen.num("valpool", pooled[1], 4)
gen.num("cappool", pooled[2], 4)
gen.num("valwithin", within[0], 4)
gen.num("capwithin", within[1], 4)
gen.num("valslope", beta_slope[1], 4)
gen.num("valslopese", se_slope[1], 4)
gen.num("capslope", beta_slope[2], 4)
gen.num("sdslope", np.sqrt(dvar[1]) / 1000, 4)
gen.num("llint", ll_int, 2)
gen.num("llslope", ll_slope, 2)
gen.num("lrt", lrt, 2)
gen.num("aicint", aic_int, 1)
gen.num("aicslope", aic_slope, 1)
gen.num("bicint", bic_int, 1)
gen.num("bicslope", bic_slope, 1)
gen.num("llstall", ll_stall, 2)
gen.num("stallgap", ll_slope - ll_stall, 2)
gen.num("slopemean", np.mean(firm_slopes), 4)
gen.num("slopesd", np.std(firm_slopes, ddof=1), 3)
gen.num("sdols", np.std(alpha_ols, ddof=1), 1)
gen.num("sdblup", np.std(u_int, ddof=1), 1)
gen.write()

# ---- figure ------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
labels = ["one intercept", "firm indicators", "random intercept", "random slopes"]
est = [pooled[1], within[0], beta_int[1], beta_slope[1]]
se = [np.sqrt(np.sum((y - X @ pooled) ** 2) / (n - p) * np.linalg.inv(X.T @ X)[1, 1]),
      np.nan, se_int[1], se_slope[1]]
Xw = np.column_stack([X[:, 1:], D])
sw = np.sum((y - Xw @ within) ** 2) / (n - Xw.shape[1])
se[1] = np.sqrt(sw * np.linalg.inv(Xw.T @ Xw)[0, 0])
for j in range(4):
    ax.plot([est[j] - 1.96 * se[j], est[j] + 1.96 * se[j]], [j, j],
            color=COLORS["accent"], lw=1.0)
    ax.plot([est[j]], [j], "o", ms=4, color=COLORS["accent"], mew=0)
ax.set_yticks(range(4))
ax.set_yticklabels(labels, fontsize=7)
ax.axvline(0, color=COLORS["ink"], lw=0.6, ls="--")
ax.set_xlabel("coefficient of market value")
ax.set_title("(a) four analyses of one panel")

ax = axes[1]
lim = [alpha_ols.min() - 15, alpha_ols.max() + 15]
ax.plot(lim, lim, color=COLORS["grid"], lw=0.8)
ax.plot(alpha_ols, u_int, "o", ms=4, color=COLORS["accent"], mew=0)
ax.set_xlim(lim)
ax.set_xlabel("firm indicator estimate")
ax.set_ylabel("predicted firm effect")
ax.set_title("(b) 11 firm effects, $m=20$")
fig.tight_layout()
fig.savefig(figure_path("ch32", "grunfeld_mixed"))
