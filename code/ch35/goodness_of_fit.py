"""Chapter 35, Section 4: goodness of fit for binary data.

Grouped Pearson and deviance statistics on the seven party-identification groups of the 1996
ANES survey; the exact identity that makes the ungrouped deviance useless; a simulation of its
null distribution; the Hosmer-Lemeshow statistic at three group counts; cross-validated
calibration; and the area under the ROC curve. Data: statsmodels.datasets.anes96 (public domain).
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

warnings.filterwarnings("ignore")

# <<grouped>>
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

anes = sm.datasets.anes96.load_pandas().data
cells = anes.groupby("PID")["vote"].agg(["sum", "count"]).reset_index()
counts = np.column_stack([cells["sum"], cells["count"] - cells["sum"]])
Xg = sm.add_constant(cells[["PID"]])

g = sm.GLM(counts, Xg, family=sm.families.Binomial()).fit()
print(cells.to_string(index=False))
print(f"deviance {g.deviance:.3f} on {g.df_resid} df, p = {stats.chi2.sf(g.deviance, g.df_resid):.4f}")
print(f"Pearson  {g.pearson_chi2:.3f} on {g.df_resid} df, "
      f"p = {stats.chi2.sf(g.pearson_chi2, g.df_resid):.4f}")
# <</grouped>>

m = cells["count"].to_numpy()
s_obs = cells["sum"].to_numpy()
pi_g = g.fittedvalues.to_numpy()
pearson = np.sum((s_obs - m * pi_g) ** 2 / (m * pi_g * (1 - pi_g)))
assert np.isclose(pearson, g.pearson_chi2)
assert m.min() >= 30                                    # every group is large
res_g = (s_obs - m * pi_g) / np.sqrt(m * pi_g * (1 - pi_g))
assert np.argmin(res_g) == 2 and np.argmax(res_g) == 4  # the misfit straddles the midpoint

# a quadratic term does not repair the fit
cells["PID2"] = cells["PID"] ** 2
gq = sm.GLM(counts, sm.add_constant(cells[["PID", "PID2"]]), family=sm.families.Binomial()).fit()
assert gq.deviance > g.deviance - 0.1                   # essentially no improvement

# <<entropy>>
# the ungrouped deviance is a function of the fitted probabilities alone
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"].to_numpy()
fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
p = fit.fittedvalues.to_numpy()

entropy = -2 * np.sum(p * np.log(p) + (1 - p) * np.log(1 - p))
print("deviance reported  ", round(fit.deviance, 6))
print("entropy of the fit ", round(entropy, 6))
# <</entropy>>

assert np.isclose(entropy, fit.deviance)                # the identity of the text, exactly

# the same holds for a wrong model fitted to the same responses
wrong = sm.GLM(y, X[["const", "educ"]], family=sm.families.Binomial()).fit()
pw = wrong.fittedvalues.to_numpy()
assert np.isclose(-2 * np.sum(pw * np.log(pw) + (1 - pw) * np.log(1 - pw)), wrong.deviance)


# <<devsim>>
def deviance_null(n, B, seed=0):
    """Ungrouped deviances under three data-generating mechanisms: a correct model with a
    strong slope, a correct model with a weak slope, and a model that omits a real quadratic
    effect. All three fit an intercept and one slope to n observations."""
    rng = np.random.default_rng(seed)
    out = np.empty((B, 3))
    for b in range(B):
        x = rng.normal(size=n)
        Z = np.column_stack([np.ones(n), x])
        for j, eta in enumerate([0.4 + 1.5 * x, 0.4 + 0.3 * x, 0.4 + 1.5 * x - 1.8 * x ** 2]):
            yb = (rng.random(n) < 1 / (1 + np.exp(-eta))).astype(float)
            out[b, j] = sm.GLM(yb, Z, family=sm.families.Binomial()).fit().deviance
    return out


demo = deviance_null(200, 100, seed=1)
print("n = 200 and p = 2, so a chi-squared reference would have mean 198 and "
      f"standard deviation {np.sqrt(2 * 198):.1f}")
for j, label in enumerate(["correct, strong slope", "correct, weak slope  ",
                           "quadratic omitted    "]):
    print(f"{label}: mean {demo[:, j].mean():6.1f}, "
          f"standard deviation {demo[:, j].std():4.1f}")
# <</devsim>>

sim = deviance_null(200, 2000, seed=0)
# two correct models with the same n and p give deviances hundreds apart
assert sim[:, 1].mean() - sim[:, 0].mean() > 50
# and the wrong model's deviance is smaller than the correct weak-slope model's
assert sim[:, 2].mean() < sim[:, 1].mean()
assert sim[:, 0].std() < np.sqrt(2 * 198) and sim[:, 2].std() < np.sqrt(2 * 198)
overlap = (sim[:, 2] > np.quantile(sim[:, 1], 0.05)).mean()

# <<hl>>
def hosmer_lemeshow(y, p, groups):
    """The Hosmer-Lemeshow statistic and its nominal p-value with the given number of groups."""
    order = np.argsort(p)
    chunks = np.array_split(order, groups)
    stat = 0.0
    for c in chunks:
        e = p[c].sum()
        stat += (y[c].sum() - e) ** 2 / (e * (1 - e / len(c)))
    return stat, stats.chi2.sf(stat, groups - 2)


for gr in (8, 10, 12):
    stat, pv = hosmer_lemeshow(y, p, gr)
    print(f"{gr:2d} groups: statistic {stat:6.3f}, p = {pv:.3f}")
# <</hl>>

hl = {gr: hosmer_lemeshow(y, p, gr) for gr in (8, 10, 12)}
assert max(v[1] for v in hl.values()) / min(v[1] for v in hl.values()) > 2   # p-values disagree


# ---- cross-validated calibration and discrimination ------------------------------------------
def cv_predict(X, y, folds=10, seed=3):
    rng = np.random.default_rng(seed)
    who = rng.permutation(len(y)) % folds
    out = np.empty(len(y))
    for k in range(folds):
        tr, te = who != k, who == k
        f = sm.GLM(y[tr], X[tr], family=sm.families.Binomial()).fit()
        out[te] = f.predict(X[te])
    return out


Xa = X.to_numpy()
cv = cv_predict(Xa, y)
eta_cv = np.log(cv / (1 - cv))
cal = sm.GLM(y, sm.add_constant(eta_cv), family=sm.families.Binomial()).fit()
# in sample the calibration slope is exactly one; out of sample it need not be
eta_in = np.log(p / (1 - p))
cal_in = sm.GLM(y, sm.add_constant(eta_in), family=sm.families.Binomial()).fit()
assert np.isclose(cal_in.params[1], 1.0, atol=1e-6) and abs(cal_in.params[0]) < 1e-6


def auc(y, score):
    """Area under the ROC curve = Mann-Whitney probability that a case outranks a control."""
    r = stats.rankdata(score)
    n1 = y.sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1))


def roc(y, score):
    thresh = np.sort(np.unique(score))[::-1]
    tpr = [(score[y == 1] >= t).mean() for t in thresh]
    fpr = [(score[y == 0] >= t).mean() for t in thresh]
    return np.r_[0, fpr, 1], np.r_[0, tpr, 1]


pid_only = sm.GLM(y, X[["const", "PID"]], family=sm.families.Binomial()).fit()
auc_full, auc_pid = auc(y, p), auc(y, pid_only.fittedvalues.to_numpy())
assert auc_full > auc_pid and auc_full - auc_pid < 0.02      # three extra regressors buy little
assert np.isclose(auc_full, roc_area := np.trapezoid(*roc(y, p)[::-1]), atol=1e-3)

gen = Generated("ch35", "fit")
gen.num("dev_grouped", g.deviance, 3)
gen.int("df_grouped", int(g.df_resid))
gen.num("p_dev", stats.chi2.sf(g.deviance, g.df_resid), 4)
gen.num("pearson", g.pearson_chi2, 3)
gen.num("p_pearson", stats.chi2.sf(g.pearson_chi2, g.df_resid), 4)
gen.num("dev_quad", gq.deviance, 3)
gen.num("res2", res_g[2], 3)
gen.num("res4", res_g[4], 3)
gen.num("obs2", s_obs[2] / m[2], 3)
gen.num("fit2", pi_g[2], 3)
gen.num("obs4", s_obs[4] / m[4], 3)
gen.num("fit4", pi_g[4], 3)
gen.int("m_min", int(m.min()))
gen.int("m_max", int(m.max()))
gen.num("dev_ungrouped", fit.deviance, 2)
gen.int("df_ungrouped", int(fit.df_resid))
gen.num("dev_wrong", wrong.deviance, 2)
gen.int("df_wrong", int(wrong.df_resid))
gen.num("sim_mean", sim[:, 0].mean(), 1)
gen.num("sim_sd", sim[:, 0].std(), 1)
gen.num("sim_mean_weak", sim[:, 1].mean(), 1)
gen.num("sim_sd_weak", sim[:, 1].std(), 1)
gen.num("sim_mean_bad", sim[:, 2].mean(), 1)
gen.num("sim_sd_bad", sim[:, 2].std(), 1)
gen.num("chisq_sd", np.sqrt(2 * 198), 1)
gen.num("overlap", overlap, 3)
for gr in (8, 10, 12):
    gen.num(f"hl{gr}", hl[gr][0], 3)
    gen.num(f"hlp{gr}", hl[gr][1], 3)
gen.num("cal_slope", cal.params[1], 3)
gen.num("cal_int", cal.params[0], 3)
gen.num("auc_full", auc_full, 3)
gen.num("auc_pid", auc_pid, 3)
gen.num("auc_gap", auc_full - auc_pid, 3)
gen.write()

# ---- three views of fit ----------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.3))

ax = axes[0]
ax.hist(sim[:, 0], bins=40, density=True, color=COLORS["accent"], alpha=0.75,
        label="correct, strong")
ax.hist(sim[:, 1], bins=40, density=True, color=COLORS["third"], alpha=0.55,
        label="correct, weak")
ax.hist(sim[:, 2], bins=40, density=True, color=COLORS["second"], alpha=0.55,
        label="quadratic omitted")
grid = np.linspace(150, 300, 300)
ax.plot(grid, stats.chi2.pdf(grid, 198), color=COLORS["ink"], linewidth=0.9,
        label=r"$\chi^2(198)$")
ax.set_xlabel("ungrouped deviance")
ax.set_ylabel("density")
ax.set_yticks([])
ax.legend(frameon=False, fontsize=6)
ax.set_title("(a) null distribution")

ax = axes[1]
order = np.argsort(cv)
chunks = np.array_split(order, 10)
ax.plot([0, 1], [0, 1], color=COLORS["grid"], linewidth=0.8)
xs = [cv[c].mean() for c in chunks]
ys = [y[c].mean() for c in chunks]
lo = [stats.beta.ppf(0.025, y[c].sum() + 0.5, len(c) - y[c].sum() + 0.5) for c in chunks]
hi = [stats.beta.ppf(0.975, y[c].sum() + 0.5, len(c) - y[c].sum() + 0.5) for c in chunks]
ax.errorbar(xs, ys, yerr=[np.array(ys) - lo, hi - np.array(ys)], fmt="o", ms=3,
            color=COLORS["accent"], linewidth=0.8)
ax.set_xlabel("cross-validated probability")
ax.set_ylabel("observed proportion")
ax.set_title("(b) calibration")

ax = axes[2]
f1, t1 = roc(y, p)
f2, t2 = roc(y, pid_only.fittedvalues.to_numpy())
ax.plot(f1, t1, color=COLORS["accent"], label=f"four regressors ({auc_full:.3f})")
ax.plot(f2, t2, color=COLORS["second"], linestyle="--", label=f"party only ({auc_pid:.3f})")
ax.plot([0, 1], [0, 1], color=COLORS["grid"], linewidth=0.8)
ax.set_xlabel("false positive rate")
ax.set_ylabel("true positive rate")
ax.legend(frameon=False, fontsize=6, loc="lower right")
ax.set_title("(c) ROC curves")
fig.tight_layout()
fig.savefig(figure_path("ch35", "fit_diagnostics"))
