"""Chapter 27, Section 5: the risk of the James-Stein estimator, and shrinkage in a regression.

(1) Z ~ N_p(theta, I) with p = 10. Exact risk of Z (p), of James-Stein (Poisson mixture formula),
    of a fixed proportional shrinker Z/(1 + lam) and of the oracle linear shrinker; positive-part
    James-Stein by simulation. Simulation also checks the exact formula.
(2) Longley's regression: the James-Stein shrinkage factor with estimated sigma^2 and its relation
    to the F statistic, and a simulation at Longley's design confirming domination in the
    prediction loss (b - beta)^T X^T X (b - beta).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import integrate, optimize
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<risk>>
p = 10


def js_risk(gamma, p=p):
    """Exact risk of (1 - (p-2)/||Z||^2) Z when Z ~ N_p(theta, I) and gamma = ||theta||^2."""
    k = np.arange(0, 400)
    w = stats.poisson.pmf(k, gamma / 2)            # ||Z||^2 ~ chi^2(p, gamma): Poisson(gamma/2) mixture
    return p - (p - 2) ** 2 * np.sum(w / (p - 2 + 2 * k))


def js(z, positive=False):
    factor = 1 - (len(z) - 2) / np.sum(z ** 2)
    return (max(factor, 0.0) if positive else factor) * z


rng = np.random.default_rng(27)
for gamma in [0.0, 5.0, 20.0, 80.0]:
    theta = np.sqrt(gamma / p) * np.ones(p)
    Zs = theta + rng.normal(size=(2_000, p))
    loss_js = np.mean([np.sum((js(z) - theta) ** 2) for z in Zs])
    loss_pp = np.mean([np.sum((js(z, True) - theta) ** 2) for z in Zs])
    print(f"||theta||^2 = {gamma:5.1f}: exact JS risk {js_risk(gamma):.3f}, "
          f"simulated {loss_js:.3f}, positive part {loss_pp:.3f}  (MLE: {p})")
# <</risk>>

assert np.isclose(js_risk(0.0), 2.0)                   # risk 2 at theta = 0
gammas = np.linspace(0, 100, 101)
risk_js = np.array([js_risk(g) for g in gammas])
assert np.all(risk_js < p) and np.all(np.diff(risk_js) > 0)
check_rng = np.random.default_rng(5)
for gamma in [3.0, 30.0]:
    theta = np.sqrt(gamma) * np.eye(p)[0]               # the risk depends on theta only via its length
    Zs = theta + check_rng.normal(size=(100_000, p))
    losses = np.array([np.sum((js(z) - theta) ** 2) for z in Zs])
    assert abs(losses.mean() - js_risk(gamma)) < 4 * losses.std() / np.sqrt(len(losses))

# positive part, simulated on the same grid with common random numbers
base = np.random.default_rng(11).normal(size=(4000, p))
risk_pp, risk_js_sim = [], []
for g in gammas:
    theta = np.sqrt(g) * np.eye(p)[0]
    Zs = theta + base
    risk_pp.append(np.mean([np.sum((js(z, True) - theta) ** 2) for z in Zs]))
    risk_js_sim.append(np.mean([np.sum((js(z) - theta) ** 2) for z in Zs]))
risk_pp = np.array(risk_pp)
assert np.all(risk_pp <= np.array(risk_js_sim) + 1e-12)   # pathwise-averaged: positive part is better

lam_fixed = 0.5
risk_ridge = (p + lam_fixed ** 2 * gammas) / (1 + lam_fixed) ** 2
risk_oracle = p * gammas / (p + gammas)
assert np.all(risk_oracle <= risk_js + 1e-12) and np.all(risk_js - risk_oracle <= 2 + 1e-12)
crossing = gammas[np.argmax(risk_ridge > p)]
assert np.isclose(crossing, 5 * p + 1)                  # fixed shrinkage exceeds p once gamma > 5p


def fixed_risk(g):
    return (p + lam_fixed ** 2 * g) / (1 + lam_fixed) ** 2


# fixed shrinkage is worse than James-Stein at the origin, better only on a middle interval
assert fixed_risk(0.0) > js_risk(0.0)
fixed_lo = optimize.brentq(lambda g: fixed_risk(g) - js_risk(g), 0.1, 20.0)
fixed_hi = optimize.brentq(lambda g: fixed_risk(g) - js_risk(g), 20.0, 60.0)
mid = np.linspace(fixed_lo + 1e-6, fixed_hi - 1e-6, 50)
assert all(fixed_risk(g) < js_risk(g) for g in mid)
assert all(fixed_risk(g) > js_risk(g) for g in np.r_[np.linspace(0, fixed_lo - 1e-6, 20),
                                                      np.linspace(fixed_hi + 1e-6, 100, 20)])

# exact risk of the positive part at theta = 0: E[(1 - (p-2)/U)_+^2 U] with U ~ chi^2(p)
pp_0_exact = integrate.quad(lambda u: (u - (p - 2)) ** 2 / u * stats.chi2.pdf(u, p), p - 2, np.inf)[0]
assert abs(pp_0_exact - risk_pp[0]) < 0.05                # agrees with the simulation

# ---- a regression: Longley ------------------------------------------------------------------
# <<longley>>
data = sm.datasets.longley.load_pandas().data
Zr = data.drop(columns="TOTEMP").to_numpy()
y = data["TOTEMP"].to_numpy()
n, k = Zr.shape
Xc = Zr - Zr.mean(axis=0)                              # centred slopes; the intercept is ybar
yc = y - y.mean()
b = np.linalg.lstsq(Xc, yc, rcond=None)[0]
m = n - k - 1                                          # residual degrees of freedom
sse = np.sum((yc - Xc @ b) ** 2)
fit_ss = b @ Xc.T @ Xc @ b                             # b^T X^T X b, the regression sum of squares
F = (fit_ss / k) / (sse / m)
factor = 1 - (k - 2) / (m + 2) * sse / fit_ss
print(f"F = {F:.1f}; James-Stein factor {factor:.5f} = 1 - (k-2) m / ((m+2) k F)")
# <</longley>>
assert np.isclose(factor, 1 - (k - 2) * m / ((m + 2) * k * F))

# domination at Longley's design, with a weak true signal so that shrinkage matters
S = Xc.T @ Xc
sim = np.random.default_rng(3)
beta_true = b / 60                                      # a much weaker signal in the same direction
L = np.linalg.cholesky(S)
loss_ls, loss_js = [], []
for _ in range(20_000):
    yy = Xc @ beta_true + np.sqrt(sse / m) * sim.normal(size=n)
    yy = yy - yy.mean()
    bb = np.linalg.solve(S, Xc.T @ yy)
    ss = np.sum((yy - Xc @ bb) ** 2)
    bj = (1 - (k - 2) / (m + 2) * ss / (bb @ S @ bb)) * bb
    loss_ls.append(np.sum((L.T @ (bb - beta_true)) ** 2))
    loss_js.append(np.sum((L.T @ (bj - beta_true)) ** 2))
ratio = np.mean(loss_js) / np.mean(loss_ls)
print(f"Longley design, weak signal: risk ratio JS / LS = {ratio:.3f}")
assert ratio < 1 and abs(np.mean(loss_ls) / (sse / m) - k) < 0.1
gamma_weak = beta_true @ S @ beta_true / (sse / m)

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.6))
ax.axhline(p, color=COLORS["muted"], linestyle="--", linewidth=0.9, label="$\\mathbf{Z}$ (least squares)")
ax.plot(gammas, risk_js, color=COLORS["accent"], label="James–Stein (exact)")
ax.plot(gammas, risk_pp, color=COLORS["third"], label="positive part (simulated)")
ax.plot(gammas, risk_ridge, color=COLORS["second"], label=f"fixed shrinkage $\\mathbf{{Z}}/(1+{lam_fixed})$")
ax.plot(gammas, risk_oracle, color=COLORS["ink"], linestyle=":", label="oracle linear shrinkage")
ax.set_ylim(0, 14)
ax.set_xlabel(r"$\|\mathbf{\theta}\|^2/\sigma^2$")
ax.set_ylabel(r"risk / $\sigma^2$")
ax.legend(frameon=False, fontsize=6.5, loc="lower right")
fig.tight_layout()
fig.savefig(figure_path("ch27", "james_stein_risk"))

gen = Generated("ch27", "james_stein")
gen.num("risk_5", js_risk(5.0), 3)
gen.num("risk_20", js_risk(20.0), 3)
gen.num("risk_80", js_risk(80.0), 3)
gen.num("pp_0", pp_0_exact, 3)
gen.num("fixed_0", fixed_risk(0.0), 2)
gen.num("fixed_lo", fixed_lo, 1)
gen.num("fixed_hi", fixed_hi, 1)
gen.num("crossing", crossing, 0)
gen.num("F", F, 1)
gen.num("factor", factor, 5)
gen.int("m", m)
gen.num("ratio", ratio, 3)
gen.num("gamma_weak", gamma_weak, 1)
gen.write()
