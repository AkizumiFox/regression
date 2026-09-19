"""Chapter 12, Section 2: a confidence ellipse for two coefficients, its shadows, and the F test.

Murder rate on poverty, single-parent households and urbanization for the 50 US states
(statsmodels.datasets.statecrime, public domain). The ellipse is the 95% region for the
poverty and single-parent coefficients.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

# <<ellipse>>
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
C = np.linalg.inv(X.T @ X)
beta_hat = C @ X.T @ y
s2 = np.sum((y - X @ beta_hat) ** 2) / (n - p)

L = np.eye(p)[:, [1, 2]]                        # Lambda: the poverty and single coefficients
q = L.shape[1]
phi_hat = L.T @ beta_hat
W = L.T @ C @ L                                  # Cov(phi_hat) / sigma^2
Fq = stats.f.ppf(0.95, q, n - p)


def F_stat(d):
    """F statistic for the hypothesis Lambda^T beta = d."""
    u = phi_hat - d
    return u @ np.linalg.solve(W, u) / (q * s2)


def in_ellipse(d):
    return F_stat(d) <= Fq


corr = W[0, 1] / np.sqrt(W[0, 0] * W[1, 1])
print(f"estimates {phi_hat}, correlation of the estimates {corr:.3f}")
print(f"F({q}, {n - p}) 95% point {Fq:.3f}; F statistic for (0, 0): {F_stat(np.zeros(2)):.2f}")
# <</ellipse>>

# <<shadows>>
tq = stats.t.ppf(0.975, n - p)
se = np.sqrt(s2 * np.diag(W))
t_box = np.column_stack([phi_hat - tq * se, phi_hat + tq * se])            # one-at-a-time t intervals
k = np.sqrt(q * Fq)
shadow = np.column_stack([phi_hat - k * se, phi_hat + k * se])            # projections of the ellipse
print("t intervals:\n", t_box.round(3))
print("shadows of the ellipse:\n", shadow.round(3), f"\nwider by the factor {k / tq:.3f}")
# <</shadows>>

# two instructive points: a corner of the t box, and a point on the long axis outside the box
A_pt = t_box[:, 1].copy()                        # upper corner (both coefficients at their upper limits)
evals, evecs = np.linalg.eigh(W)
long_axis = evecs[:, np.argmax(evals)]
long_axis *= np.sign(long_axis[0])
B_pt = phi_hat + 0.95 * np.sqrt(q * Fq * s2 * evals.max()) * long_axis
tA = (phi_hat - A_pt) / se
tB = (phi_hat - B_pt) / se

# ---- checks ---------------------------------------------------------------------------------
# (1) the F statistic equals the reduced-model form (SSE0 - SSE)/q / s^2
def reduced_F(d):
    y_star = y - X[:, 1] * d[0] - X[:, 2] * d[1]          # offset the constrained coefficients
    X0 = X[:, [0, 3]]
    sse0 = np.sum((y_star - X0 @ np.linalg.lstsq(X0, y_star, rcond=None)[0]) ** 2)
    return (sse0 - s2 * (n - p)) / q / s2


for d in (np.zeros(2), A_pt, B_pt, phi_hat + np.array([0.1, -0.2])):
    assert np.isclose(F_stat(d), reduced_F(d))
# (2) the corner of the box is outside the ellipse, the long-axis point inside it but outside the box
assert not in_ellipse(A_pt) and np.all(np.abs(tA) <= tq + 1e-12)
assert in_ellipse(B_pt) and np.any(np.abs(tB) > tq)
assert not in_ellipse(np.zeros(2))
# (3) shadow = support function: max of a^T phi over the ellipse is a^T phi_hat + sqrt(qF s^2 a^T W a)
Wc = np.linalg.cholesky(q * Fq * s2 * W)
theta = np.linspace(0, 2 * np.pi, 20001)
boundary = phi_hat[:, None] + Wc @ np.vstack([np.cos(theta), np.sin(theta)])
assert np.allclose(boundary.max(axis=1), shadow[:, 1], atol=1e-6)
assert np.allclose(boundary.min(axis=1), shadow[:, 0], atol=1e-6)
for a in (np.array([1.0, -1.0]), np.array([2.0, 1.0])):
    assert np.isclose((a @ boundary).max(), a @ phi_hat + np.sqrt(q * Fq * s2 * a @ W @ a), atol=1e-6)
# (4) semi-axes of the ellipse
semi = np.sqrt(q * Fq * s2 * evals)
area = np.pi * semi.prod()
assert np.isclose(area, np.pi * q * Fq * s2 * np.sqrt(np.linalg.det(W)))
box_area = np.prod(2 * tq * se)
# (4b) area ratio for q = 2 depends on the estimates only through their correlation
ratio0 = np.pi * Fq / (2 * tq**2)                  # ratio when the estimates are uncorrelated
area_ratio = ratio0 * np.sqrt(1 - corr**2)
assert np.isclose(area / box_area, area_ratio)
rho_equal = np.sqrt(1 - 1 / ratio0**2)             # |rho| at which the two areas agree
assert ratio0 > 1 and area_ratio < 1 and abs(corr) > rho_equal
# (5) coverage of the ellipse and of the t box by simulation (the fitted model as truth)
rng = np.random.default_rng(1202)
reps = 40_000
Y = X @ beta_hat + np.sqrt(s2) * rng.normal(size=(reps, n))
B = Y @ (C @ X.T).T
S2 = np.sum((Y - B @ X.T) ** 2, axis=1) / (n - p)
D = B[:, [1, 2]] - phi_hat
Fsim = np.einsum("ij,jk,ik->i", D, np.linalg.inv(W), D) / (q * S2)
cov_ellipse = np.mean(Fsim <= Fq)
Tsim = np.abs(D) / np.sqrt(S2[:, None] * np.diag(W))
cov_box = np.mean(np.all(Tsim <= tq, axis=1))
print(f"simulated coverage: ellipse {cov_ellipse:.4f}, box of t intervals {cov_box:.4f}")
assert abs(cov_ellipse - 0.95) < 4 * np.sqrt(0.95 * 0.05 / reps)
assert cov_box < 0.93

gen = Generated("ch12", "ellipse", prefix="ell")
gen.num("pov", phi_hat[0], 4)
gen.num("sing", phi_hat[1], 4)
gen.num("corr", corr, 3)
gen.num("F", Fq, 3)
gen.num("tq", tq, 3)
gen.num("k", k, 3)
gen.num("ratio", k / tq, 3)
gen.num("F00", F_stat(np.zeros(2)), 2)
for j, name in enumerate(["pov", "sing"]):
    gen.num(f"box:{name}:lo", t_box[j, 0], 3)
    gen.num(f"box:{name}:hi", t_box[j, 1], 3)
    gen.num(f"sh:{name}:lo", shadow[j, 0], 3)
    gen.num(f"sh:{name}:hi", shadow[j, 1], 3)
gen.num("A:pov", A_pt[0], 3)
gen.num("A:sing", A_pt[1], 3)
gen.num("A:F", F_stat(A_pt), 3)
gen.num("A:p", stats.f.sf(F_stat(A_pt), q, n - p), 4)
gen.num("B:pov", B_pt[0], 3)
gen.num("B:sing", B_pt[1], 3)
gen.num("B:F", F_stat(B_pt), 3)
gen.num("B:tpov", abs(tB[0]), 3)
gen.num("B:tsing", abs(tB[1]), 3)
gen.num("area", area, 4)
gen.num("boxarea", box_area, 4)
gen.num("arearatio", area_ratio, 3)
gen.num("ratio0", ratio0, 3)
gen.num("rhoequal", rho_equal, 3)
gen.num("cov:ell", cov_ellipse, 4)
gen.num("cov:box", cov_box, 4)
gen.int("reps", reps)
gen.write()

# ---- figure ----------------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 3.4))
ax.plot(*boundary, color=COLORS["accent"], label="95% confidence ellipse")
ax.add_patch(plt.Rectangle(t_box[:, 0], *(t_box[:, 1] - t_box[:, 0]), fill=False,
                           edgecolor=COLORS["second"], linewidth=1.0, label="95% $t$ intervals"))
ax.add_patch(plt.Rectangle(shadow[:, 0], *(shadow[:, 1] - shadow[:, 0]), fill=False,
                           edgecolor=COLORS["muted"], linewidth=0.8, linestyle="--",
                           label="shadows of the ellipse"))
ax.plot(*phi_hat, "o", color=COLORS["ink"], markersize=3)
ax.plot(*A_pt, "s", color=COLORS["second"], markersize=4)
ax.annotate("A", A_pt, textcoords="offset points", xytext=(4, 3), fontsize=8)
ax.plot(*B_pt, "^", color=COLORS["third"], markersize=4)
ax.annotate("B", B_pt, textcoords="offset points", xytext=(4, -9), fontsize=8)
ax.set_xlabel("coefficient of poverty")
ax.set_ylabel("coefficient of single-parent households")
ax.legend(loc="lower left", frameon=False)
pad = 0.04
ax.set_xlim(shadow[0, 0] - pad, shadow[0, 1] + pad)
ax.set_ylim(shadow[1, 0] - 0.12, shadow[1, 1] + pad)
fig.tight_layout()
fig.savefig(figure_path("ch12", "confidence_ellipse"))
