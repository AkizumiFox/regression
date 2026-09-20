"""Chapter 35, Section 3: separation in logistic regression, and three remedies.

The 13 respondents of the 1996 ANES survey whose schooling stopped at or before the eighth grade
(statsmodels.datasets.anes96, public domain) are completely separated by party identification:
everyone below the midpoint of the scale expected to vote Clinton, everyone above it Dole.
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import linprog

from regbook import COLORS, Generated, figure_path, use_book_style

warnings.filterwarnings("ignore")

# <<data>>
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
sub = anes[anes["educ"] == 1]                      # schooling stopped at or before grade eight
X = sm.add_constant(sub[["PID"]]).to_numpy()
y = sub["vote"].to_numpy()

print(pd.crosstab(sub["PID"], sub["vote"]))
for maxiter in (5, 10, 25, 50):
    f = sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=maxiter, tol=1e-16)
    print(f"{maxiter:3d} steps: slope {f.params[1]:10.3f}   standard error {f.bse[1]:14.3f}"
          f"   deviance {f.deviance:.2e}")
# <</data>>

n = len(y)
s = 2 * y - 1


def separating_direction(X, y):
    """Maximize the smallest margin min_i s_i x_i'b over the cube |b_j| <= 1.

    A strictly positive value certifies complete separation; a value of zero together with a
    nonzero direction along which no margin is negative certifies quasi-complete separation.
    """
    S = (2 * y - 1)[:, None] * X
    p = X.shape[1]
    # variables (b, t): maximize t subject to s_i x_i'b - t >= 0
    A = np.column_stack([-S, np.ones(len(y))])
    r = linprog(c=np.r_[np.zeros(p), -1.0], A_ub=A, b_ub=np.zeros(len(y)),
                bounds=[(-1, 1)] * p + [(0, None)], method="highs")
    return r.x[p], r.x[:p]


def quasi_direction(X, y):
    """Maximize the sum of margins subject to none being negative; zero means overlap."""
    S = (2 * y - 1)[:, None] * X
    r = linprog(c=-S.sum(0), A_ub=-S, b_ub=np.zeros(len(y)),
                bounds=[(-1, 1)] * X.shape[1], method="highs")
    return -r.fun, r.x


margin, direction = separating_direction(X, y)
assert margin > 1e-8                                # completely separated
assert np.all(s * (X @ direction) > 0)              # the certificate really separates
full_margin, _ = separating_direction(
    sm.add_constant(anes[["PID"]]).to_numpy(), anes["vote"].to_numpy())
full_quasi, _ = quasi_direction(
    sm.add_constant(anes[["PID"]]).to_numpy(), anes["vote"].to_numpy())
assert full_margin < 1e-8 and full_quasi < 1e-8     # the whole survey overlaps


def loglik(beta, X, y):
    eta = X @ beta
    return -np.sum(np.log1p(np.exp(-(2 * y - 1) * eta)))


# the log-likelihood along the separating ray climbs to zero and never attains it
ray = np.array([loglik(t * direction, X, y) for t in [1, 4, 16, 64, 256]])
assert np.all(np.diff(ray) > 0) and ray[-1] < 0


# <<firth>>
def firth(X, y, steps=200):
    """Maximize the log-likelihood penalized by (1/2) log det(X'WX): Firth's estimate."""
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        mu = 1 / (1 + np.exp(-X @ beta))
        w = mu * (1 - mu)
        F = X.T @ (w[:, None] * X)                          # the expected information
        h = w * np.einsum("ij,jk,ik->i", X, np.linalg.inv(F), X)
        beta = beta + np.linalg.solve(F, X.T @ (y - mu + h * (0.5 - mu)))
    return beta


beta_firth = firth(X, y)
print("Firth estimate:", beta_firth.round(3))
# <</firth>>

mu_f = 1 / (1 + np.exp(-X @ beta_firth))
w_f = mu_f * (1 - mu_f)
F_f = X.T @ (w_f[:, None] * X)
h_f = w_f * np.einsum("ij,jk,ik->i", X, np.linalg.inv(F_f), X)
assert np.allclose(X.T @ (y - mu_f + h_f * (0.5 - mu_f)), 0, atol=1e-8)
se_firth = np.sqrt(np.diag(np.linalg.inv(F_f)))


def penalized_loglik(beta, X, y):
    mu = 1 / (1 + np.exp(-X @ beta))
    w = np.clip(mu * (1 - mu), 1e-300, None)
    sign, logdet = np.linalg.slogdet(X.T @ (w[:, None] * X))
    return loglik(beta, X, y) + 0.5 * logdet


assert penalized_loglik(beta_firth, X, y) > penalized_loglik(beta_firth * 1.2, X, y)
assert penalized_loglik(beta_firth, X, y) > penalized_loglik(beta_firth * 0.8, X, y)


def profile(b1, X, y, penalized=False):
    """Maximize over the intercept with the slope held at b1."""
    f = penalized_loglik if penalized else loglik
    lo, hi = -400.0, 400.0
    for _ in range(200):                                     # golden-free ternary search
        a, c = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if f(np.array([a, b1]), X, y) < f(np.array([c, b1]), X, y):
            lo = a
        else:
            hi = c
    return f(np.array([(lo + hi) / 2, b1]), X, y)


cut = stats.chi2.ppf(0.95, 1)
# the ordinary profile likelihood has supremum 0 under complete separation
assert profile(60.0, X, y) > -1e-6
lo, hi = 0.0, 40.0
for _ in range(80):                                          # lower profile limit
    mid = (lo + hi) / 2
    (lo, hi) = (mid, hi) if 2 * (0.0 - profile(mid, X, y)) > cut else (lo, mid)
prof_lower = (lo + hi) / 2
assert profile(prof_lower, X, y) < 0

# the penalized profile likelihood does have a maximum, and a two-sided interval
top = penalized_loglik(beta_firth, X, y)
limits = []
for side in (-1, 1):
    lo, hi = beta_firth[1], beta_firth[1] + side * 40.0
    for _ in range(80):
        mid = (lo + hi) / 2
        (lo, hi) = (mid, hi) if 2 * (top - profile(mid, X, y, True)) < cut else (lo, mid)
    limits.append((lo + hi) / 2)
assert limits[0] < beta_firth[1] < limits[1]


# <<ridge>>
def ridge_logistic(X, y, lam, steps=200):
    """Maximize the log-likelihood minus (lam/2) times the squared norm of the slopes."""
    d = np.ones(X.shape[1])
    d[0] = 0.0                                               # leave the intercept unpenalized
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        mu = 1 / (1 + np.exp(-X @ beta))
        w = mu * (1 - mu)
        H = X.T @ (w[:, None] * X) + lam * np.diag(d)
        beta = beta + np.linalg.solve(H, X.T @ (y - mu) - lam * d * beta)
    return beta


for lam in (0.1, 1.0, 10.0):
    print(f"ridge, lambda = {lam:4.1f}: {ridge_logistic(X, y, lam).round(3)}")
# <</ridge>>

ridge = {lam: ridge_logistic(X, y, lam) for lam in (0.1, 1.0, 10.0)}
assert ridge[0.1][1] > ridge[1.0][1] > ridge[10.0][1] > 0     # heavier penalty, smaller slope

# ---- quasi-complete separation: a factor with an empty cell ---------------------------------
phd = anes[anes["educ"] == 7]
dummies = pd.get_dummies(phd["PID"].astype(int), prefix="pid").astype(float)
Xq = sm.add_constant(dummies.drop(columns="pid_4"))           # PID = 4, the baseline category
fq = sm.GLM(phd["vote"].to_numpy(), Xq.to_numpy(),
            family=sm.families.Binomial()).fit(maxiter=100, tol=1e-14)
q_margin, _ = separating_direction(Xq.to_numpy(), phd["vote"].to_numpy())
q_sum, q_dir = quasi_direction(Xq.to_numpy(), phd["vote"].to_numpy())
assert q_margin < 1e-8 and q_sum > 1e-8                       # quasi-complete, not complete
params = pd.Series(fq.params, index=Xq.columns)
finite = np.abs(params) < 10
assert set(params.index[~finite]) == {"pid_0", "pid_2", "pid_6"}
assert finite.sum() == 4                                      # the other four are ordinary
q_cells = phd.groupby("PID")["vote"].agg(["sum", "count"])

gen = Generated("ch35", "separation")
gen.int("n", n)
gen.int("dole", int(y.sum()))
gen.num("slope10", sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=10, tol=1e-16).params[1], 2)
gen.num("se10", sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=10, tol=1e-16).bse[1], 1)
gen.num("slope25", sm.GLM(y, X, family=sm.families.Binomial()).fit(maxiter=25, tol=1e-16).params[1], 2)
gen.num("firth_const", beta_firth[0], 3)
gen.num("firth_slope", beta_firth[1], 3)
gen.num("firth_se", se_firth[1], 3)
gen.num("firth_wald", beta_firth[1] / se_firth[1], 2)
gen.num("prof_lower", prof_lower, 3)
gen.num("firth_lo", limits[0], 3)
gen.num("firth_hi", limits[1], 3)
gen.num("ridge_01", ridge[0.1][1], 3)
gen.num("ridge_1", ridge[1.0][1], 3)
gen.num("ridge_10", ridge[10.0][1], 3)
gen.int("phd_n", len(phd))
gen.int("phd_infinite", int((~finite).sum()))
gen.num("phd_pid1", params["pid_1"], 3)
gen.num("phd_se1", fq.bse[list(Xq.columns).index("pid_1")], 3)
gen.num("phd_dev", fq.deviance, 3)
gen.write()

# ---- profile log-likelihoods ----------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.6))
grid = np.linspace(0.05, 12, 120)
plain = np.array([profile(b, X, y) for b in grid])
pen = np.array([profile(b, X, y, True) for b in grid])
ax.plot(grid, plain, color=COLORS["accent"], label="log-likelihood")
ax.plot(grid, pen - top, color=COLORS["second"], linestyle="--", label="Firth penalized")
ax.axhline(-cut / 2, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.axvline(beta_firth[1], color=COLORS["second"], linewidth=0.6, linestyle=":", zorder=0)
ax.set_xlabel(r"slope $\beta_1$ on the party scale")
ax.set_ylabel("profile, relative to its supremum")
ax.set_ylim(-6, 0.4)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch35", "separation_profile"))
