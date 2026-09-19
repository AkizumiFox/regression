"""Chapter 19, Section 6: collinearity inflates the variance of coefficients, not of fits.

(a) Variance inflation factors for Longley's (1967) regression (public domain,
    statsmodels.datasets.longley), computed as 1/(1 - R_j^2) and as diagonal entries of the
    inverse correlation matrix.
(b) Two regressors with correlation 0.98: the standard deviation of the fitted mean at a new
    point is small along the data's direction and large across it (figure).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch19", "collinearity", prefix="col")

# <<vif>>
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
Z = data.drop(columns="TOTEMP").to_numpy()           # six regressors
n, k = Z.shape

def r_squared(v, W):
    W1 = np.column_stack([np.ones(len(v)), W])
    fitted = W1 @ np.linalg.lstsq(W1, v, rcond=None)[0]
    return 1 - np.sum((v - fitted) ** 2) / np.sum((v - v.mean()) ** 2)

vif = np.array([1 / (1 - r_squared(Z[:, j], np.delete(Z, j, axis=1))) for j in range(k)])
R = np.corrcoef(Z, rowvar=False)
print(dict(zip(data.columns[1:], np.round(vif, 1))))
print("diagonal of R^{-1}:", np.round(np.diag(np.linalg.inv(R)), 1))
print("eigenvalues of R:  ", np.round(np.linalg.eigvalsh(R), 5))
# <</vif>>
assert np.allclose(vif, np.diag(np.linalg.inv(R)))
lam = np.linalg.eigvalsh(R)
assert np.isclose(lam.sum(), k)
names = list(data.columns[1:])
jmax = int(np.argmax(vif))
gen.text("vif_max_name", names[jmax])
gen.num("vif_max", vif[jmax], 1)
gen.num("vif_min", vif.min(), 2)
gen.text("vif_min_name", names[int(np.argmin(vif))])
gen.num("lam_min", lam[0], 5)
gen.num("lam_max", lam[-1], 3)
gen.num("sum_inv_lam", np.sum(1 / lam), 1)
gen.num("sum_vif", vif.sum(), 1)
gen.int("n", n)

# <<prediction>>
rng = np.random.default_rng(1906)
m, rho = 40, 0.98
x1 = rng.normal(size=m)
x2 = rho * x1 + np.sqrt(1 - rho**2) * rng.normal(size=m)
X = np.column_stack([np.ones(m), x1, x2])
G = np.linalg.inv(X.T @ X)                           # Cov(beta_hat) / sigma^2
along = np.array([1.0, 1.0, 1.0])                    # a new point that follows the pattern
across = np.array([1.0, 1.0, -1.0])                  # same distance, against the pattern
print("sd of beta_1, beta_2 (units of sigma):", np.round(np.sqrt(np.diag(G)[1:]), 3))
print("sd of the fitted mean, along :", round(np.sqrt(along @ G @ along), 3))
print("sd of the fitted mean, across:", round(np.sqrt(across @ G @ across), 3))
print("average over the design points:", round(np.trace(X @ G @ X.T) / m, 3), "= p/n")
# <</prediction>>
r12 = np.corrcoef(x1, x2)[0, 1]
assert np.isclose(np.trace(X @ G @ X.T), 3)
sd_along = np.sqrt(along @ G @ along)
sd_across = np.sqrt(across @ G @ across)
assert sd_across > 5 * sd_along
# eigen-decomposition formula for the variance of a linear function
lamG, V = np.linalg.eigh(X.T @ X)
assert np.isclose(across @ G @ across, np.sum((V.T @ across) ** 2 / lamG))
gen.num("r12", r12, 4)
gen.int("m", m)
gen.num("sd_b1", np.sqrt(G[1, 1]), 3)
gen.num("sd_b2", np.sqrt(G[2, 2]), 3)
gen.num("sd_along", sd_along, 3)
gen.num("sd_across", sd_across, 3)
gen.num("sd_design", np.sqrt(3 / m), 3)
gen.num("vif_sim", 1 / (1 - r12**2), 1)
gen.write()

# ---- figure -----------------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(3.4, 3.0))
g = np.linspace(-3.6, 3.6, 289)
A, B = np.meshgrid(g, g)
P = np.stack([np.ones_like(A), A, B], axis=-1)
sd = np.sqrt(np.einsum("...i,ij,...j->...", P, G, P))
cs = ax.contour(A, B, sd, levels=[0.2, 0.4, 0.8, 1.6, 3.2], colors=COLORS["accent"], linewidths=0.8)
ax.clabel(cs, levels=cs.levels[2:], fmt="%.1f", fontsize=7)
ax.scatter(x1, x2, s=8, color=COLORS["ink"], zorder=3)
ax.plot([1], [1], "o", color=COLORS["third"], ms=5)
ax.plot([1], [-1], "s", color=COLORS["second"], ms=5)
ax.annotate("along", (1, 1), xytext=(1.3, 0.5), fontsize=8, color=COLORS["third"])
ax.annotate("across", (1, -1), xytext=(1.3, -1.5), fontsize=8, color=COLORS["second"])
ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
lim = (min(x1.min(), x2.min()) - 0.4, max(x1.max(), x2.max()) + 0.4)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(figure_path("ch19", "collinearity"))
