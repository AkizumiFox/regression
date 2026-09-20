"""Chapter 37, Section 3: log-linear models for a contingency table and the multinomial fit.

Data: the 1996 American National Election Studies sample shipped with statsmodels
(statsmodels.datasets.anes96, public domain), 944 respondents, cross-classified by
education (four levels) and party identification (three levels).

Fits three nested pairs of models -- independence, a linear education score, and the
saturated model -- once as Poisson log-linear models for the twelve cell counts and
once as baseline-category multinomial logit models for party given education, and
checks that each pair gives the same fitted counts and the same association
parameters.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import Generated

# <<table>>
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
party = pd.cut(anes["PID"], [-0.5, 1.5, 4.5, 6.5],
               labels=["Democrat", "Independent", "Republican"])
education = pd.cut(anes["educ"], [0.5, 2.5, 3.5, 4.5, 7.5],
                   labels=["no diploma", "high school", "some college", "degree"])
table = pd.crosstab(education, party)
print(table)
# <</table>>

counts = table.to_numpy(float)
r, c = counts.shape
n = counts.sum()
score = np.arange(1.0, r + 1)                 # education score u_i = 1, 2, 3, 4

# <<loglinear>>
def poisson_mle(X, y, tol=1e-12, maxit=80):
    """Fisher scoring for log mu = X beta."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        XW = X * mu[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


def dummies(levels, k):
    """Indicator columns for levels 1, ..., k-1 of a factor taking the values 0..k-1."""
    return np.column_stack([(levels == j).astype(float) for j in range(1, k)])


rows, cols = np.divmod(np.arange(counts.size), c)
y = counts.ravel()
main = np.column_stack([np.ones(counts.size), dummies(rows, r), dummies(cols, c)])
linear = np.column_stack([main, score[rows][:, None] * dummies(cols, c)])
saturated = np.column_stack([main] + [dummies(rows, r)[:, [i]] * dummies(cols, c)[:, [j]]
                                      for i in range(r - 1) for j in range(c - 1)])

fits = {}
for name, X in [("independence", main), ("linear", linear), ("saturated", saturated)]:
    beta = poisson_mle(X, y)
    fits[name] = (beta, np.exp(X @ beta).reshape(counts.shape))

print("independence fit:\n", np.round(fits["independence"][1], 2))
print("row totals match:", np.allclose(fits["independence"][1].sum(1), counts.sum(1)))
print("column totals match:", np.allclose(fits["independence"][1].sum(0), counts.sum(0)))
# <</loglinear>>

# the independence fit is the outer product of the observed margins
assert np.allclose(fits["independence"][1], np.outer(counts.sum(1), counts.sum(0)) / n)
assert np.allclose(fits["saturated"][1], counts, atol=1e-6)

# <<tests>>
def deviance(observed, fitted):
    return 2 * np.sum(observed * np.log(observed / fitted))


for name, dfree in [("independence", 6), ("linear", 4)]:
    dev = deviance(counts, fits[name][1])
    pear = np.sum((counts - fits[name][1]) ** 2 / fits[name][1])
    print(f"{name:13s} G2 = {dev:6.3f}  X2 = {pear:6.3f}  df = {dfree}  "
          f"p = {stats.chi2.sf(dev, dfree):.4f}")
# <</tests>>

dev_ind = deviance(counts, fits["independence"][1])
dev_lin = deviance(counts, fits["linear"][1])
pear_ind = np.sum((counts - fits["independence"][1]) ** 2 / fits["independence"][1])
pear_lin = np.sum((counts - fits["linear"][1]) ** 2 / fits["linear"][1])
assert dev_lin < dev_ind
assert stats.chi2.sf(dev_ind, 6) < 0.05 < stats.chi2.sf(dev_lin, 4)

# ---- the same three models as multinomial logits ---------------------------
# <<multinomial>>
def multinomial_logit(counts, design, tol=1e-12, maxit=100):
    """Baseline-category logits: log(pi_ik / pi_i0) = design_i . gamma_k, k = 1..c-1."""
    r, c = counts.shape
    q = design.shape[1]
    gamma = np.zeros((c - 1) * q)
    total = counts.sum(1)
    for _ in range(maxit):
        eta = np.column_stack([np.zeros(r)] +
                              [design @ gamma[k * q:(k + 1) * q] for k in range(c - 1)])
        pi = np.exp(eta) / np.exp(eta).sum(1, keepdims=True)
        mu = total[:, None] * pi
        u = np.concatenate([design.T @ (counts[:, k + 1] - mu[:, k + 1])
                            for k in range(c - 1)])
        info = np.zeros((len(gamma), len(gamma)))
        for k in range(c - 1):
            for l in range(c - 1):
                w = total * pi[:, k + 1] * ((k == l) - pi[:, l + 1])
                info[k * q:(k + 1) * q, l * q:(l + 1) * q] = design.T @ (design * w[:, None])
        step = np.linalg.solve(info, u)
        gamma = gamma + step
        if np.max(np.abs(step)) < tol:
            break
    return gamma, mu


designs = {"independence": np.ones((r, 1)),
           "linear": np.column_stack([np.ones(r), score]),
           "saturated": np.column_stack([np.ones(r), dummies(np.arange(r), r)])}
for name, design in designs.items():
    gamma, mu = multinomial_logit(counts, design)
    print(f"{name:13s} largest difference in fitted counts: "
          f"{np.abs(mu - fits[name][1]).max():.2e}")
# <</multinomial>>

for name, design in designs.items():
    gamma, mu = multinomial_logit(counts, design)
    assert np.allclose(mu, fits[name][1], atol=1e-6)

# the association parameters agree: the education slopes of the multinomial fit are the
# log-linear score-by-party coefficients
gamma_lin, _ = multinomial_logit(counts, designs["linear"])
beta_lin = fits["linear"][0]
assert np.allclose(gamma_lin[[1, 3]], beta_lin[-2:], atol=1e-7)

# a saturated log-linear interaction parameter is a log odds ratio: the last one compares
# the "degree" row with "no diploma" and the Republican column with the Democrat column
beta_sat = fits["saturated"][0]
lam = beta_sat[main.shape[1] + (r - 2) * (c - 1) + (c - 2)]
odds_ratio = (counts[r - 1, c - 1] * counts[0, 0]) / (counts[r - 1, 0] * counts[0, c - 1])
assert np.isclose(np.exp(lam), odds_ratio)

# where the independence model fails: Pearson residuals, cell by cell
resid = (counts - fits["independence"][1]) / np.sqrt(fits["independence"][1])
worst = np.unravel_index(np.abs(resid).argmax(), resid.shape)
assert worst == (0, c - 1)                    # no diploma, Republican
assert np.argmax(resid) == 0                  # the largest positive one is in the same row
assert fits["independence"][1].min() > 10     # every expected count is comfortably large

gen = Generated("ch37", "contingency")
gen.int("n", int(n))
gen.num("dev_ind", dev_ind, 3)
gen.num("pear_ind", pear_ind, 3)
gen.num("p_ind", stats.chi2.sf(dev_ind, 6), 4)
gen.num("dev_lin", dev_lin, 3)
gen.num("pear_lin", pear_lin, 3)
gen.num("p_lin", stats.chi2.sf(dev_lin, 4), 4)
gen.num("dev_drop", dev_ind - dev_lin, 3)
gen.num("p_drop", stats.chi2.sf(dev_ind - dev_lin, 2), 5)
gen.num("slope_ind", gamma_lin[1], 4)
gen.num("slope_rep", gamma_lin[3], 4)
gen.num("rr_rep", np.exp(gamma_lin[3]), 4)
gen.num("odds_ratio", odds_ratio, 3)
gen.num("lambda11", lam, 4)
gen.num("resid_worst", resid[worst], 2)
gen.num("obs_worst", counts[worst], 0)
gen.num("fit_worst", fits["independence"][1][worst], 2)
gen.num("resid_00", resid[0, 0], 2)
gen.num("min_expected", fits["independence"][1].min(), 2)
gen.num("fit_00", fits["independence"][1][0, 0], 2)
gen.num("fit_lin_00", fits["linear"][1][0, 0], 2)
gen.int("obs_00", int(counts[0, 0]))
for k, label in enumerate(["dem", "ind", "rep"]):
    gen.int("total_" + label, int(counts.sum(0)[k]))
for k, label in enumerate(["nodip", "hs", "some", "deg"]):
    gen.int("row_" + label, int(counts.sum(1)[k]))
gen.write()
