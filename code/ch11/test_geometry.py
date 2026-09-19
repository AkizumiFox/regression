"""Chapter 11, Section 2: the F statistic as an angle, and its rejection region as a cone.

Uses the region test of Section 11.1 (murder rate, 2009 state data). The data enter the test
only through the two lengths ||(M - M0) y|| and ||(I - M) y||; F is a function of the angle
between (I - M0) y and the test space. The figure plots these two lengths for the observed data
and for simulated data sets under the reduced model and under an alternative.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(1102)

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])
D = np.column_stack([(region == g).astype(float) for g in "NSW"])
X = np.column_stack([X0, D])

# <<angle>>
def proj(Z):
    """Orthogonal projection onto C(Z), from an orthonormal basis."""
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

M, M0 = proj(X), proj(X0)
q, nu = 3, n - 7                               # r - r0 and n - r
u = (M - M0) @ y                               # component in the test space
e = y - M @ y                                  # residual vector
w = y - M0 @ y                                 # what the reduced model leaves unexplained
theta = np.arccos(np.linalg.norm(u) / np.linalg.norm(w))
F = (nu / q) / np.tan(theta) ** 2
theta_crit = np.arctan(np.sqrt(nu / (q * stats.f.ppf(0.95, q, nu))))
print(f"angle = {np.degrees(theta):.1f} degrees, F = {F:.3f}")
print(f"reject when the angle is below {np.degrees(theta_crit):.1f} degrees")
# <</angle>>

assert np.allclose(u + e, w) and abs(u @ e) < 1e-9          # orthogonal split of (I - M0) y
F_direct = (u @ u / q) / (e @ e / nu)
assert np.isclose(F, F_direct) and np.isclose(F, 2.440, atol=5e-4)
assert np.isclose(np.cos(theta) ** 2, q * F / (q * F + nu))   # partial R^2 = cos^2
assert theta > theta_crit

# canonical form: orthonormal bases of C(X0), the test space and the residual space
Qfull, _ = np.linalg.qr(np.column_stack([X0, D, rng.standard_normal((n, n - 7))]))
Q0, Q1, Q2 = Qfull[:, :4], Qfull[:, 4:7], Qfull[:, 7:]
z1, z2 = Q1.T @ y, Q2.T @ y
assert np.isclose(z1 @ z1, u @ u) and np.isclose(z2 @ z2, e @ e)
assert np.allclose(Q1 @ Q1.T, M - M0)

# ---- simulated clouds: null and an alternative, in (||u||, ||e||) coordinates ----------
b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
sigma = np.sqrt(e @ e / nu)
mu0 = X0 @ b0                                  # a mean in the reduced model
mu1 = mu0 - 2.0 * D[:, 0]                      # Northeast lowered by 2 per 100,000
reps = 4000
Z = rng.standard_normal((reps, n)) * sigma
clouds = {}
for name, mu in (("null", mu0), ("alt", mu1)):
    Y = mu + Z
    U = Y @ (M - M0)
    E = Y - Y @ M
    clouds[name] = (np.linalg.norm(U, axis=1), np.linalg.norm(E, axis=1))
rates = {}
for name, (a, b) in clouds.items():
    rates[name] = np.mean(np.arctan2(b, a) < theta_crit)
assert abs(rates["null"] - 0.05) < 0.012 and rates["alt"] > 0.3

gen = Generated("ch11", "test_geometry")
gen.num("theta", np.degrees(theta), 1)
gen.num("theta_crit", np.degrees(theta_crit), 1)
gen.num("len_u", np.linalg.norm(u), 3)
gen.num("len_e", np.linalg.norm(e), 3)
gen.num("cos2", np.cos(theta) ** 2, 3)
gen.num("rate_null", rates["null"], 3)
gen.num("rate_alt", rates["alt"], 3)
gen.int("reps", reps)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 3.2))
xmax, ymax = 8.0, 14.0
xs = np.array([0, xmax])
ax.fill_between(xs, 0, xs * np.tan(theta_crit), color=COLORS["second"], alpha=0.12, linewidth=0)
ax.plot(xs, xs * np.tan(theta_crit), color=COLORS["second"], linewidth=0.9)
a, b = clouds["null"]
ax.scatter(a[:600], b[:600], s=3, color=COLORS["muted"], alpha=0.5, linewidths=0, label="simulated, reduced model true")
a, b = clouds["alt"]
ax.scatter(a[:600], b[:600], s=3, color=COLORS["third"], alpha=0.6, linewidths=0, label="simulated, Northeast lower by 2")
ax.plot([0, np.linalg.norm(u)], [0, np.linalg.norm(e)], color=COLORS["ink"], linewidth=0.8)
ax.plot(np.linalg.norm(u), np.linalg.norm(e), "o", color=COLORS["ink"], markersize=4, label="observed data")
ax.text(xmax * 0.97, 0.8, "rejection region, 5% level", color=COLORS["second"], ha="right", fontsize=8)
ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
ax.set_xlabel(r"$\|(\mathbf{M}-\mathbf{M}_0)\mathbf{y}\|$, length in the test space")
ax.set_ylabel(r"$\|(\mathbf{I}-\mathbf{M})\mathbf{y}\|$, residual length")
leg = ax.legend(loc="upper left", frameon=False)
for handle in leg.legend_handles[:2]:
    handle.set_sizes([12])
fig.tight_layout()
fig.savefig(figure_path("ch11", "test_cone"))
