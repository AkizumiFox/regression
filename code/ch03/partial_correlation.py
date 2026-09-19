"""Chapter 3, Section 6: a sample partial correlation computed four ways.

Longley's macroeconomic series (US, 1947-1962; public domain, shipped with
statsmodels): total employment and unemployment, adjusted for the year."""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.longley.load_pandas().data

# <<pcor>>
emp = data["TOTEMP"].to_numpy() / 1000       # thousands -> millions of persons
unemp = data["UNEMP"].to_numpy() / 1000
year = data["YEAR"].to_numpy()
n = len(emp)
Z = np.column_stack([np.ones(n), year])      # the variables we adjust for


def resid(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef


# 1. correlation of residuals (detrended series)
r_resid = np.corrcoef(resid(emp, Z), resid(unemp, Z))[0, 1]

# 2. the one-variable recursion from ordinary correlations
R = np.corrcoef(np.vstack([emp, unemp, year]))
r12, r13, r23 = R[0, 1], R[0, 2], R[1, 2]
r_formula = (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2))

# 3. from the inverse of the correlation (or covariance) matrix
W = np.linalg.inv(R)
r_precision = -W[0, 1] / np.sqrt(W[0, 0] * W[1, 1])

# 4. from the t statistic of unemployment in the regression of emp on (1, unemp, year)
fit = sm.OLS(emp, np.column_stack([Z, unemp])).fit()
t = fit.tvalues[2]
r_t = np.sign(t) * np.sqrt(t**2 / (t**2 + n - 3))

print(f"marginal r = {r12:.4f}")
print(f"partial r:  {r_resid:.10f} {r_formula:.10f} {r_precision:.10f} {r_t:.10f}")
# <</pcor>>

assert np.allclose([r_formula, r_precision, r_t], r_resid, atol=1e-10)
assert r12 > 0 > r_resid
# the covariance-matrix version of method 3 gives the same number
Wc = np.linalg.inv(np.cov(np.vstack([emp, unemp, year])))
assert np.isclose(-Wc[0, 1] / np.sqrt(Wc[0, 0] * Wc[1, 1]), r_resid)
# FWL: the slope of the residual regression is the multiple-regression coefficient
eu, ee = resid(unemp, Z), resid(emp, Z)
assert np.isclose(eu @ ee / (eu @ eu), fit.params[2])

# Proposition (recursion), general S: partial correlations given S + {k} from those given S
rng = np.random.default_rng(306)
A = rng.standard_normal((5, 5))
Sig = A @ A.T


def pcov(Sig, idx, S):
    """partial covariance matrix of the variables idx given the variables S"""
    if not S:
        return Sig[np.ix_(idx, idx)]
    return (Sig[np.ix_(idx, idx)]
            - Sig[np.ix_(idx, S)] @ np.linalg.pinv(Sig[np.ix_(S, S)]) @ Sig[np.ix_(S, idx)])


def pcorr(C):
    d = 1 / np.sqrt(np.diag(C))
    return C * np.outer(d, d)


for S in ([], [3], [3, 4]):
    k = 2 if S else 4
    P = pcorr(pcov(Sig, [0, 1, k], S))
    rec = (P[0, 1] - P[0, 2] * P[1, 2]) / np.sqrt((1 - P[0, 2]**2) * (1 - P[1, 2]**2))
    assert np.isclose(pcorr(pcov(Sig, [0, 1], S + [k]))[0, 1], rec)
# residual form: covariance of the adjusted residuals equals the recursion's numerator
C = pcov(Sig, [0, 1, 2], [3])
assert np.isclose(pcov(Sig, [0, 1], [3, 2])[0, 1], C[0, 1] - C[0, 2] * C[1, 2] / C[2, 2])
# precision-matrix row of the response: (1, -beta) / sigma^2
Om = np.linalg.inv(Sig)
beta = np.linalg.solve(Sig[1:, 1:], Sig[1:, 0])
s2 = Sig[0, 0] - Sig[0, 1:] @ beta
assert np.allclose(Om[0], np.concatenate([[1.0], -beta]) / s2)

gen = Generated("ch03", "partial_correlation", prefix="pcor")
gen.int("n", n)
gen.int("first", year.min()); gen.int("last", year.max())
gen.num("r12", r12, 3)
gen.num("r13", r13, 3)
gen.num("r23", r23, 3)
gen.num("partial", r_resid, 3)
gen.num("t", t, 2)
gen.num("slope", fit.params[2], 3)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
sc = ax.scatter(unemp, emp, c=year, cmap="viridis", s=14, linewidths=0)
ax.set_xlabel("unemployed (millions)")
ax.set_ylabel("employed (millions)")
ax.set_title(f"(a) raw series: $r={r12:.2f}$")
cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.06)
cb.ax.tick_params(labelsize=7)
ax = axes[1]
ax.scatter(eu, ee, c=year, cmap="viridis", s=14, linewidths=0)
xs = np.linspace(eu.min(), eu.max(), 2)
ax.plot(xs, fit.params[2] * xs, color=COLORS["second"])
ax.axhline(0, color=COLORS["grid"], lw=0.6, zorder=0)
ax.axvline(0, color=COLORS["grid"], lw=0.6, zorder=0)
ax.set_xlabel("unemployed, detrended")
ax.set_ylabel("employed, detrended")
ax.set_title(f"(b) adjusted for year: $r={r_resid:.2f}$")
fig.tight_layout()
fig.savefig(figure_path("ch03", "partial_correlation"))
