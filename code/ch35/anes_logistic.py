"""Chapter 35, Section 1: logistic regression on the 1996 American National Election Studies survey.

Response: the respondent's expected vote (1 = Dole, 0 = Clinton). Regressors: party
identification on a seven-point scale (0 = strong Democrat, 6 = strong Republican), age in
decades, education on a seven-point scale, household income on a 24-point scale. Public-domain
data shipped with statsmodels (statsmodels.datasets.anes96, COPYRIGHT: "This is public domain").
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"]                                   # 1 = Dole, 0 = Clinton

fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
pi_hat = fit.fittedvalues

print(fit.params.round(4).to_string())
print("odds ratio for one step on the party scale:", round(np.exp(fit.params["PID"]), 3))
print("divide by four:", round(fit.params["PID"] / 4, 3))
print("average marginal effect:", round((fit.params["PID"] * pi_hat * (1 - pi_hat)).mean(), 3))
# <</fit>>

n = len(y)
b = fit.params
se = fit.bse
ame = (b["PID"] * pi_hat * (1 - pi_hat)).mean()

# Fisher scoring is IRLS: one hand-coded run must reproduce the software's answer.
Xm = X.to_numpy()
beta = np.zeros(Xm.shape[1])
for _ in range(30):
    eta = Xm @ beta
    mu = 1 / (1 + np.exp(-eta))
    w = mu * (1 - mu)
    z = eta + (y.to_numpy() - mu) / w                       # working response
    beta = np.linalg.solve(Xm.T @ (w[:, None] * Xm), Xm.T @ (w * z))
assert np.allclose(beta, b.to_numpy(), atol=1e-8)
# the score equation X^T (y - mu) = 0 holds at the canonical-link fit
assert np.allclose(Xm.T @ (y.to_numpy() - pi_hat), 0, atol=1e-7)
# with an intercept the fitted probabilities reproduce the observed number of successes
assert np.isclose(pi_hat.sum(), y.sum())

# Wald and likelihood ratio tests for education
X0 = X.drop(columns="educ")
fit0 = sm.GLM(y, X0, family=sm.families.Binomial()).fit()
wald_educ = (b["educ"] / se["educ"]) ** 2
lr_educ = fit0.deviance - fit.deviance
assert abs(wald_educ - lr_educ) < 0.01                      # the two agree to two decimals here

# a profile likelihood interval for the party coefficient
def profile_limit(name, side, level=0.95):
    cut = stats.chi2.ppf(level, 1)
    keep = [c for c in X.columns if c != name]
    lo, hi = b[name], b[name] + side * 5.0
    for _ in range(60):
        mid = (lo + hi) / 2
        off = mid * X[name]
        m = sm.GLM(y, X[keep], family=sm.families.Binomial(), offset=off).fit()
        if 2 * (fit.llf - m.llf) < cut:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

prof = (profile_limit("PID", -1), profile_limit("PID", +1))
wald_ci = (b["PID"] - 1.96 * se["PID"], b["PID"] + 1.96 * se["PID"])
assert max(abs(prof[0] - wald_ci[0]), abs(prof[1] - wald_ci[1])) < 0.01  # nearly quadratic

# <<grouped>>
# the same model in grouped form: one binomial count per distinct covariate pattern
cols = ["PID", "educ"]
cells = anes.groupby(cols)["vote"].agg(["sum", "count"]).reset_index()
Xg = sm.add_constant(cells[cols])
counts = np.column_stack([cells["sum"], cells["count"] - cells["sum"]])
grouped = sm.GLM(counts, Xg, family=sm.families.Binomial()).fit()

Xu = sm.add_constant(anes[cols])
ungrouped = sm.GLM(anes["vote"], Xu, family=sm.families.Binomial()).fit()

print("grouped  ", grouped.params.round(4).to_numpy(), "deviance", round(grouped.deviance, 2))
print("ungrouped", ungrouped.params.round(4).to_numpy(), "deviance", round(ungrouped.deviance, 2))
# <</grouped>>

assert np.allclose(grouped.params.to_numpy(), ungrouped.params.to_numpy(), atol=1e-8)
assert np.allclose(grouped.bse.to_numpy(), ungrouped.bse.to_numpy(), atol=1e-8)
# the two deviances differ by twice the log-likelihood of the grouped saturated model
p_sat = np.clip(cells["sum"] / cells["count"], 1e-12, 1 - 1e-12)
sat_ll = np.sum(cells["sum"] * np.log(p_sat) + (cells["count"] - cells["sum"]) * np.log(1 - p_sat))
assert np.isclose(ungrouped.deviance - grouped.deviance, -2 * sat_ll)

gen = Generated("ch35", "anes")
gen.int("n", n)
gen.int("dole", int(y.sum()))
gen.int("cells", len(cells))
gen.num("b_pid", b["PID"], 4)
gen.num("se_pid", se["PID"], 4)
gen.num("or_pid", np.exp(b["PID"]), 3)
gen.num("or_pid3", np.exp(3 * b["PID"]), 2)
gen.num("quarter_pid", b["PID"] / 4, 3)
gen.num("ame_pid", ame, 3)
gen.num("b_const", b["const"], 3)
gen.num("se_const", se["const"], 3)
gen.num("b_age", b["age"], 4)
gen.num("se_age", se["age"], 4)
gen.num("b_educ", b["educ"], 4)
gen.num("se_educ", se["educ"], 4)
gen.num("b_income", b["income"], 4)
gen.num("se_income", se["income"], 4)
gen.num("z_age", b["age"] / se["age"], 2)
gen.num("p_age", 2 * stats.norm.sf(abs(b["age"] / se["age"])), 3)
gen.num("wald_educ", wald_educ, 3)
gen.num("lr_educ", lr_educ, 3)
gen.num("prof_lo", prof[0], 3)
gen.num("prof_hi", prof[1], 3)
gen.num("wald_lo", wald_ci[0], 3)
gen.num("wald_hi", wald_ci[1], 3)
gen.num("dev", fit.deviance, 2)
gen.int("df_resid", int(fit.df_resid))
gen.num("dev_grouped", grouped.deviance, 2)
gen.int("df_grouped", int(grouped.df_resid))
gen.num("dev_ungrouped", ungrouped.deviance, 2)
gen.int("df_ungrouped", int(ungrouped.df_resid))
gen.int("iters", int(fit.fit_history["iteration"]))
gen.write()

# ---- fitted curve against party identification --------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))

grid = np.linspace(-0.5, 6.5, 200)
means = {c: X[c].mean() for c in ["age", "educ", "income"]}
eta = b["const"] + b["PID"] * grid + sum(b[c] * means[c] for c in means)
curve = 1 / (1 + np.exp(-eta))
obs = anes.groupby("PID")["vote"].agg(["mean", "count"])

ax = axes[0]
ax.plot(grid, curve, color=COLORS["accent"])
ax.scatter(obs.index, obs["mean"], s=obs["count"] * 0.6, color=COLORS["second"],
           zorder=3, linewidths=0)
ax.set_xlabel("party identification (0 = strong Democrat)")
ax.set_ylabel(r"$\Pr(\mathrm{Dole})$")
ax.set_ylim(-0.03, 1.03)
ax.set_title("(a) fit at average covariates")

# the tangent at the steepest point has slope beta/4
x_half = (-b["const"] - sum(b[c] * means[c] for c in means)) / b["PID"]
ax = axes[1]
ax.plot(grid, curve, color=COLORS["accent"])
tan = 0.5 + (b["PID"] / 4) * (grid - x_half)
ok = (tan > -0.03) & (tan < 1.03)
ax.plot(grid[ok], tan[ok], color=COLORS["third"], linestyle="--")
ax.axhline(0.5, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("party identification")
ax.set_ylabel(r"$\Pr(\mathrm{Dole})$")
ax.set_ylim(-0.03, 1.03)
ax.set_title(r"(b) tangent of slope $\hat\beta/4$")
fig.tight_layout()
fig.savefig(figure_path("ch35", "anes_logistic"))
