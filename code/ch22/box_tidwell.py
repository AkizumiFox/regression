"""Chapter 22, Section 3: transforming a regressor on Engel's data (Box-Tidwell), and the joint
profile likelihood of a response power and a regressor power.

Response log food expenditure, regressor income^(alpha) (Box-Cox transform of income).
Public-domain data shipped with statsmodels (statsmodels.datasets.engel).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<tidwell>>
df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
y = np.log(food)


def bc(x, a):
    """Box-Cox transform (x^a - 1)/a, with log x at a = 0."""
    return np.log(x) if abs(a) < 1e-12 else np.expm1(a * np.log(x)) / a


def fit(y, x, a):
    return sm.OLS(y, np.column_stack([np.ones(len(y)), bc(x, a)])).fit()


def box_tidwell_step(y, x, a0):
    """One Box-Tidwell (Gauss-Newton) step: a1 = a0 + gamma/beta1."""
    beta1 = fit(y, x, a0).params[1]                    # slope of x^(a0) in the null fit
    w = x**a0 * np.log(x) / a0 if a0 != 0 else np.log(x) ** 2 / 2   # d bc(x, a)/da, up to C(X0)
    aug = sm.OLS(y, np.column_stack([np.ones(len(y)), bc(x, a0), w])).fit()
    return a0 + aug.params[2] / beta1, aug.tvalues[2]


a, path = 1.0, [1.0]
for _ in range(8):
    a, t = box_tidwell_step(y, income, a)
    path.append(a)
print("Box-Tidwell iterates from alpha = 1:", " ".join(f"{v:.4f}" for v in path))
# <</tidwell>>

sse = lambda a: fit(y, income, a).ssr
res = optimize.minimize_scalar(sse, bounds=(-2, 2), method="bounded", options={"xatol": 1e-10})
a_hat = res.x
assert abs(path[-1] - a_hat) < 1e-5                       # the iteration converges to the minimizer of SSE
cut = stats.chi2.ppf(0.95, 1)
g = lambda a: n * np.log(sse(a) / sse(a_hat)) - cut
a_lo, a_hi = optimize.brentq(g, a_hat - 2, a_hat), optimize.brentq(g, a_hat, a_hat + 2)
_, t_at1 = box_tidwell_step(y, income, 1.0)
_, t_at0 = box_tidwell_step(y, income, 0.0)
assert a_lo < a_hi < 0 and abs(t_at1) > 10 and 2 < abs(t_at0) < 3     # log income rejected, narrowly

# the evidence against alpha = 0 rests largely on the richest household
rich = int(np.argmax(income))
keep = np.arange(n) != rich
sse_k = lambda a: fit(y[keep], income[keep], a).ssr
a_hat_k = optimize.minimize_scalar(sse_k, bounds=(-2, 2), method="bounded", options={"xatol": 1e-10}).x
lr0_k = (n - 1) * np.log(sse_k(0.0) / sse_k(a_hat_k))
lr0 = n * np.log(sse(0.0) / sse(a_hat))
assert lr0 > cut > lr0_k

# exact size of the constructed-variable t test at alpha0 (w is not random): simulate under the log fit
rng = np.random.default_rng(2203)
fit0 = fit(y, income, 0.0)
tt = [box_tidwell_step(fit0.fittedvalues + np.sqrt(fit0.scale) * rng.standard_normal(n), income, 0.0)[1]
      for _ in range(4000)]
size_bt = np.mean(np.abs(tt) > stats.t.ppf(0.975, n - 3))
assert 0.04 < size_bt < 0.06                             # exact t(n-3) under the null


# ---- joint profile of (lambda, alpha): response Box-Cox, regressor Box-Tidwell ----------------
def z(yv, lam):
    gm = np.exp(np.mean(np.log(yv)))
    return gm * np.log(yv) if abs(lam) < 1e-12 else np.expm1(lam * np.log(yv)) / (lam * gm ** (lam - 1))


def joint_profile(lam, a, rows=slice(None)):
    xv, fv = income[rows], food[rows]
    Xa = np.column_stack([np.ones(len(xv)), bc(xv, a)])
    coef, *_ = np.linalg.lstsq(Xa, z(fv, lam), rcond=None)
    r = z(fv, lam) - Xa @ coef
    return -0.5 * len(xv) * np.log(r @ r)


def joint_mle(rows=slice(None)):
    opt = optimize.minimize(lambda v: -joint_profile(*v, rows), x0=[0.3, 0.3], method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-10})
    return opt.x, -opt.fun


(lam_j, a_j), l_j = joint_mle()
lr_00 = 2 * (l_j - joint_profile(0, 0))
lr_11 = 2 * (l_j - joint_profile(1, 1))
lam_best = optimize.minimize_scalar(lambda l: -joint_profile(l, 1), bounds=(-2, 2), method="bounded",
                                    options={"xatol": 1e-8}).x   # the best lambda with income untransformed
assert abs(lam_best - 0.631) < 1e-3                             # agrees with bc:lam_lin of box_cox.py
lr_lin = 2 * (l_j - joint_profile(lam_best, 1))
_, l_jk = joint_mle(keep)
lr_00_k = 2 * (l_jk - joint_profile(0, 0, keep))
chi2_2 = stats.chi2.ppf(0.95, 2)
assert chi2_2 < lr_00 < 3 * chi2_2 < lr_lin < lr_11 and lr_00_k < lr_00

gen = Generated("ch22", "box_tidwell", prefix="bt")
gen.num("a_hat", a_hat, 3)
gen.num("a_lo", a_lo, 3)
gen.num("a_hi", a_hi, 3)
gen.num("step1", path[1], 3)
gen.num("step2", path[2], 3)
gen.num("t_at1", t_at1, 2)
gen.num("t_at0", t_at0, 2)
gen.num("size_bt", size_bt, 3)
gen.num("t_crit", stats.t.ppf(0.975, n - 3), 2)
gen.num("lam_j", lam_j, 2)
gen.num("a_j", a_j, 2)
gen.num("lr_00", lr_00, 1)
gen.num("lr_11", lr_11, 1)
gen.num("lr_lin", lr_lin, 1)
gen.num("lr_00_k", lr_00_k, 2)
gen.num("chi2_2", chi2_2, 2)
gen.num("a_hat_k", a_hat_k, 3)
gen.num("lr0", lr0, 2)
gen.num("lr0_k", lr0_k, 2)
gen.num("rich_income", income[rich], 0)
gen.num("rich_food", food[rich], 0)
gen.write()

# ---- figure: contours of the joint profile log-likelihood ------------------------------------
use_book_style()
lams = np.linspace(-0.8, 1.1, 131)
alphas = np.linspace(-1.4, 1.4, 131)
L = np.array([[joint_profile(l, a) - l_j for l in lams] for a in alphas])
fig, ax = plt.subplots(figsize=(3.6, 3.0))
cs = ax.contour(lams, alphas, L, levels=[-40, -20, -10, -chi2_2 / 2],
                colors=[COLORS["muted"]] * 3 + [COLORS["second"]], linewidths=[0.6, 0.6, 0.6, 1.2])
ax.plot([lam_j], [a_j], "o", color=COLORS["accent"], markersize=4)
ax.plot([0, 1], [0, 1], "s", color=COLORS["ink"], markersize=3.5)
ax.annotate("log-log", (0, 0), xytext=(6, -10), textcoords="offset points", fontsize=8)
ax.annotate("raw", (1, 1), xytext=(-18, 5), textcoords="offset points", fontsize=8)
ax.set_xlabel(r"response power $\lambda$")
ax.set_ylabel(r"income power $\alpha$")
fig.tight_layout()
fig.savefig(figure_path("ch22", "joint_profile"))
