"""Chapter 37, Section 2: offsets, exposure, grouping, and what a free exposure coefficient costs.

Part 1 (real data): the RAND Health Insurance Experiment person-years
(statsmodels.datasets.randhie, public domain) are collapsed to the distinct
patterns of five categorical covariates. The grouped Poisson fit with offset
log(person-years) reproduces the person-level fit exactly.

Part 2 (simulation): counts with a known exposure. Three analyses -- offset fixed
at one, a free coefficient on log t, and no exposure term at all -- are compared
when the true mean is proportional to exposure and when it is not.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# ---- Part 1: grouping a Poisson fit ----------------------------------------
# <<grouping>>
import numpy as np
import pandas as pd
import statsmodels.api as sm


def scoring(X, y, offset=0.0, tol=1e-11, maxit=60):
    """Fisher scoring for a Poisson log-linear model with a fixed offset."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 1e-6))
    for _ in range(maxit):
        mu = np.exp(X @ beta + offset)
        z = X @ beta + (y - mu) / mu                 # working response
        XW = X * mu[:, None]                         # weights w_i = mu_i
        step = np.linalg.solve(X.T @ XW, XW.T @ z)
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


data = sm.datasets.randhie.load_pandas().data.copy()
data["band"] = pd.cut(data["disea"], [-0.1, 6, 10, 15, 60], labels=[0, 1, 2, 3])
cells = ["idp", "hlthg", "hlthf", "hlthp", "band"]   # all categorical
y = data["mdvis"].to_numpy(float)


def design(frame):
    """Intercept, the plan and health indicators, and three disease-band indicators."""
    columns = [np.ones(len(frame))] + [frame[v].to_numpy(float)
                                       for v in ["idp", "hlthg", "hlthf", "hlthp"]]
    band = frame["band"].to_numpy(int)
    return np.column_stack(columns + [(band == k).astype(float) for k in (1, 2, 3)])


grouped = data.groupby(cells, as_index=False, observed=True).agg(
    visits=("mdvis", "sum"), person_years=("mdvis", "size"))
s = grouped["visits"].to_numpy(float)                # total visits in the group
t = grouped["person_years"].to_numpy(float)          # exposure: person-years
X, Xg = design(data), design(grouped)

beta_person = scoring(X, y)                          # one row per person-year
beta_group = scoring(Xg, s, offset=np.log(t))        # one row per pattern, offset log t

print("groups:", len(s), " largest coefficient difference:",
      np.abs(beta_person - beta_group).max())
print("fitted rate, free-plan person-years in excellent health:",
      np.exp(beta_group[0]))
# <</grouping>>

assert np.allclose(beta_person, beta_group, atol=1e-8)
mu_g = t * np.exp(Xg @ beta_group)
assert np.allclose(Xg.T @ (s - mu_g), 0.0, atol=1e-6)

# the same fit read as rates, and the deviance of the grouped data
dev_group = 2 * np.sum(np.where(s > 0, s * np.log(np.where(s > 0, s, 1) / mu_g), 0.0)
                       - (s - mu_g))
df_group = len(s) - Xg.shape[1]

# is the lack of fit in the mean model or in the variance function? the grouped design is
# main effects only, so refit with every two-factor interaction among the three factors
# (plan, self-rated health, disease band) -- the homogeneous-association model
health = np.column_stack([grouped[v].to_numpy(float) for v in ["hlthg", "hlthf", "hlthp"]])
bands = np.column_stack([(grouped["band"].to_numpy(int) == k).astype(float) for k in (1, 2, 3)])
plan = grouped["idp"].to_numpy(float)
Xint = np.column_stack([Xg, plan[:, None] * health, plan[:, None] * bands]
                       + [health[:, [i]] * bands[:, [j]] for i in range(3) for j in range(3)])
assert np.linalg.matrix_rank(Xint) == Xint.shape[1]
beta_int = scoring(Xint, s, offset=np.log(t))
mu_int = t * np.exp(Xint @ beta_int)
dev_int = 2 * np.sum(np.where(s > 0, s * np.log(np.where(s > 0, s, 1) / mu_int), 0.0)
                     - (s - mu_int))
df_int = len(s) - Xint.shape[1]
dev_drop = dev_group - dev_int
assert df_int == 9 and dev_drop > 0
assert dev_drop > 0.8 * dev_group          # most of the lack of fit is the mean model

# estimating the exposure coefficient instead of fixing it
Xpsi = np.column_stack([Xg, np.log(t)])
beta_psi = scoring(Xpsi, s)
mu_psi = np.exp(Xpsi @ beta_psi)
se_psi = np.sqrt(np.diag(np.linalg.inv(Xpsi.T @ (Xpsi * mu_psi[:, None]))))[-1]
dev_psi = 2 * np.sum(np.where(s > 0, s * np.log(np.where(s > 0, s, 1) / mu_psi), 0.0)
                     - (s - mu_psi))
lr_psi = dev_group - dev_psi
assert lr_psi > 0

gen = Generated("ch37", "offsets")
gen.int("groups", len(s))
gen.int("n", len(y))
gen.num("rate_base", np.exp(beta_group[0]), 4)
gen.num("rr_idp", np.exp(beta_group[1]), 4)
gen.num("rr_hlthp", np.exp(beta_group[4]), 4)
gen.num("rr_band3", np.exp(beta_group[7]), 4)
gen.num("dev_group", dev_group, 2)
gen.int("df_group", df_group)
gen.num("dev_int", dev_int, 2)
gen.int("df_int", df_int)
gen.num("dev_drop", dev_drop, 1)
gen.int("df_drop", df_group - df_int)
gen.num("psi", beta_psi[-1], 4)
gen.num("psi_se", se_psi, 4)
gen.num("psi_lo", beta_psi[-1] - 1.96 * se_psi, 4)
gen.num("psi_hi", beta_psi[-1] + 1.96 * se_psi, 4)
gen.num("lr_psi", lr_psi, 3)
gen.num("smallest_fitted", mu_g.min(), 1)
gen.num("smallest_group", t.min(), 0)
gen.num("largest_group", t.max(), 0)
gen.write()

# ---- Part 2: simulation ----------------------------------------------------
rng = np.random.default_rng(20370202)
reps, m = 2000, 300


def simulate(psi_true, rng):
    """One data set: log mu = 0.4 + 0.5 x + psi_true * log t, with t correlated with x."""
    x = rng.normal(size=m)
    t = np.exp(1.5 + 0.6 * x + 0.4 * rng.normal(size=m))     # exposure, correlated with x
    mu = np.exp(0.4 + 0.5 * x + psi_true * np.log(t))
    return x, t, rng.poisson(mu).astype(float)


def three_fits(x, t, y):
    """Slope of x under: offset log t; free coefficient on log t; no exposure term."""
    X0 = np.column_stack([np.ones(m), x])
    b_off = scoring(X0, y, offset=np.log(t))[1]
    b_free = scoring(np.column_stack([X0, np.log(t)]), y)[1]
    b_none = scoring(X0, y)[1]
    return b_off, b_free, b_none


results = {}
for psi_true in (1.0, 0.7):
    rows = [three_fits(*simulate(psi_true, rng)) for _ in range(reps)]
    results[psi_true] = np.array(rows)

prop = results[1.0]
assert abs(prop[:, 0].mean() - 0.5) < 0.01          # offset model unbiased
assert abs(prop[:, 1].mean() - 0.5) < 0.01          # free coefficient unbiased
assert prop[:, 1].std() > prop[:, 0].std()          # but less precise
assert prop[:, 2].mean() > 0.8                      # dropping exposure inflates the slope

nonprop = results[0.7]
assert abs(nonprop[:, 1].mean() - 0.5) < 0.02       # free coefficient still unbiased
assert nonprop[:, 0].mean() < 0.42                  # fixing psi = 1 biases the slope

gen2 = Generated("ch37", "offset_simulation")
gen2.int("reps", reps)
gen2.int("m", m)
for tag, res in (("prop", prop), ("nonprop", nonprop)):
    for k, label in enumerate(["off", "free", "none"]):
        gen2.num(f"{tag}_{label}_mean", res[:, k].mean(), 3)
        gen2.num(f"{tag}_{label}_sd", res[:, k].std(ddof=1), 3)
gen2.num("efficiency", prop[:, 1].std(ddof=1) / prop[:, 0].std(ddof=1), 3)
gen2.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4), sharey=True)
titles = [r"(a) true exposure exponent $\psi=1$", r"(b) true exposure exponent $\psi=0.7$"]
for ax, (psi_true, res), title in zip(axes, results.items(), titles):
    bins = np.linspace(0.25, 1.30, 70)
    for k, (label, colour) in enumerate([("offset fixed at 1", "accent"),
                                         ("free coefficient", "third"),
                                         ("no exposure term", "second")]):
        ax.hist(res[:, k], bins=bins, histtype="step", color=COLORS[colour],
                linewidth=1.0, label=label)
    ax.axvline(0.5, color=COLORS["ink"], linewidth=0.7, linestyle=(0, (3, 2)))
    ax.set_title(title)
    ax.set_xlabel(r"estimated slope of $x$")
axes[0].set_ylabel("frequency")
axes[0].legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch37", "offset_simulation"))
