"""Chapter 34, Section 3: the likelihood equations and the information.

Checks the score equations at the fit, the equality of observed and expected information
under the canonical link (and their inequality otherwise), and the fitted-marginal
property: under the canonical link the fitted means reproduce the weighted column totals
of X exactly.

Data: statsmodels.datasets.anes96 (public domain; a subset of the 1996 American National
Election Studies). Response: the respondent's reported vote (1 = Dole). Regressors: an
indicator for each of the seven party-identification levels, and age.
"""
import numpy as np
import statsmodels.api as sm

from regbook import Generated

# <<fit>>
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["vote"].to_numpy(float)                       # 1 if the respondent voted for Dole
pid = anes["PID"].to_numpy(int)                        # party identification, 0 (Dem) to 6 (Rep)
D = np.eye(7)[pid]                                     # one indicator column per level
age = (anes["age"].to_numpy(float) - 45.0) / 10.0
X = np.column_stack([D, age])                          # no separate intercept: D already sums to 1

logit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
mu = logit.fittedvalues

print("level   n   observed Dole votes   fitted")
for k in range(7):
    inside = pid == k
    print(f"  {k}   {inside.sum():4d}      {y[inside].sum():8.0f}   {mu[inside].sum():12.6f}")
print(f"age column: X^T y = {age @ y:.6f},  X^T mu = {age @ mu:.6f}")
# <</fit>>

# The score equations are exactly X^T (y - mu) = 0 for the canonical link.
assert np.allclose(X.T @ (y - mu), 0.0, atol=1e-8)
for k in range(7):
    assert abs(y[pid == k].sum() - mu[pid == k].sum()) < 1e-8

# Observed and expected information agree under the canonical link.
W = mu * (1 - mu)
expected = (X.T * W) @ X
hess_obs = -logit.model.hessian(logit.params)
assert np.allclose(expected, hess_obs, rtol=1e-8)
assert np.allclose(np.sqrt(np.diag(np.linalg.inv(expected))), logit.bse, atol=1e-8)

# With a probit link neither statement holds.
probit = sm.GLM(y, X, family=sm.families.Binomial(sm.families.links.Probit())).fit()
mu_p = probit.fittedvalues
resid_p = X.T @ (y - mu_p)
gap_ind = np.abs(resid_p[:7]).max()                    # the seven indicator columns
gap_age = abs(resid_p[7])
assert gap_ind > 1e-3 and gap_age > 0.1                # the column totals are not reproduced
from scipy import stats
eta_p = X @ probit.params
d = stats.norm.pdf(eta_p)
W_p = d ** 2 / (mu_p * (1 - mu_p))
expected_p = (X.T * W_p) @ X
hess_obs_p = -probit.model.hessian(probit.params)
rel = np.abs(expected_p - hess_obs_p).max() / np.abs(expected_p).max()
assert rel > 1e-3                                      # observed and expected differ

gen = Generated("ch34", "score")
gen.int("n", len(y))
gen.int("dole0", int(y[pid == 0].sum()))
gen.int("dole6", int(y[pid == 6].sum()))
gen.int("n0", int((pid == 0).sum()))
gen.int("n6", int((pid == 6).sum()))
gen.num("age_total", float(age @ y), 4)
gen.num("probit_gap_ind", float(gap_ind), 4)
gen.num("probit_gap_age", float(gap_age), 4)
gen.num("info_rel", float(rel), 4)
gen.num("beta_age", float(logit.params[-1]), 4)
gen.num("se_age", float(logit.bse[-1]), 4)
gen.write()
