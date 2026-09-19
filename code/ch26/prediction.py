"""Chapter 26, Section 4: collinearity harms coefficients, not predictions that follow the data.

The consumption regression on the US macroeconomic data (public domain,
statsmodels.datasets.macrodata). Exact standard deviations of the fitted mean at an
observed quarter (1985Q1), at a convex combination of observed quarters, and at a point that breaks the
income-population pattern; their decomposition over the canonical directions of the scaled
design; and a parametric simulation of coefficients and predictions (figure).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch26", "prediction")

# <<setup>>
import numpy as np
import statsmodels.api as sm
macro = sm.datasets.macrodata.load_pandas().data
names = ["realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
X = sm.add_constant(macro[names].to_numpy())
y = macro["realcons"].to_numpy()
n, p = X.shape
fit = sm.OLS(y, X).fit()
s = np.sqrt(fit.scale)                                # estimate of sigma
G = np.linalg.inv(X.T @ X)

x_on = X[104]                                         # 1985Q1, a design point
x_mix = 0.5 * (X[40] + X[160])                        # midway between 1969Q1 and 1999Q1
x_off = x_on.copy()
x_off[2] *= 0.90                                      # 1985Q1 with population 10% lower
for label, x0 in [("1985Q1", x_on), ("midway", x_mix), ("off pattern", x_off)]:
    print(f"{label:13s} fitted {x0 @ fit.params:8.1f}   sd {s * np.sqrt(x0 @ G @ x0):6.1f}")
# <</setup>>

sd_on = s * np.sqrt(x_on @ G @ x_on)
sd_mix = s * np.sqrt(x_mix @ G @ x_mix)
sd_off = s * np.sqrt(x_off @ G @ x_off)
Q = np.linalg.qr(X)[0]
h = np.sum(Q**2, axis=1)
# convex combination of design rows: Var <= sigma^2 max h_ii; any x0 = X'w: Var <= sigma^2 ||w||^2
assert x_mix @ G @ x_mix <= h.max()
w = np.zeros(n); w[40] = w[160] = 0.5
assert x_mix @ G @ x_mix <= w @ w
# average over the design rows is p/n
assert np.isclose(np.mean(np.einsum("ij,jk,ik->i", X, G, X)), p / n)

# <<decompose>>
d = np.linalg.norm(X, axis=0)
U, mu, Vt = np.linalg.svd(X / d, full_matrices=False) # scaled design, as in Section 26.3

def shares(x0):
    """Contributions of the canonical directions to Var(x0' beta_hat), as fractions."""
    c = (Vt @ (x0 / d)) ** 2 / mu**2                  # (x0_scaled' v_k)^2 / mu_k^2
    return c / c.sum()

np.set_printoptions(precision=3, suppress=True)
print("condition indices:", mu[0] / mu)
print("shares, 1985Q1     :", shares(x_on))
print("shares, off pattern :", shares(x_off))
# <</decompose>>
sh_on, sh_off = shares(x_on), shares(x_off)
# the variance formula of the eigen-decomposition
for x0 in (x_on, x_off):
    c = (Vt @ (x0 / d)) ** 2 / mu**2
    assert np.isclose(c.sum(), x0 @ G @ x0)
assert sh_off[-1] > 0.8 and sh_off[-1] > 3 * sh_on[-1]

# <<simulate>>
rng = np.random.default_rng(2604)
reps = 2000
Ystar = X @ fit.params + s * rng.normal(size=(reps, n))      # parametric replicates, same design
Bstar = np.linalg.solve(X.T @ X, X.T @ Ystar.T).T
print("sd of income, population coefficients:", Bstar[:, 1:3].std(axis=0))
print("correlation of the two:", np.corrcoef(Bstar[:, 1], Bstar[:, 2])[0, 1])
print("sd of predictions, 1985Q1 / off:", (Bstar @ x_on).std(), (Bstar @ x_off).std())
# <</simulate>>
corr_b = np.corrcoef(Bstar[:, 1], Bstar[:, 2])[0, 1]
corr_exact = G[1, 2] / np.sqrt(G[1, 1] * G[2, 2])
assert abs(corr_b - corr_exact) < 0.02
assert abs((Bstar @ x_off).std() / sd_off - 1) < 0.05
assert abs((Bstar @ x_on).std() / sd_on - 1) < 0.05

gen.num("s", s, 1)
gen.num("sd_on", sd_on, 1)
gen.num("sd_mix", sd_mix, 1)
gen.num("sd_off", sd_off, 1)
gen.num("ratio_off", sd_off / sd_on, 1)
gen.num("h_on", h[104], 3)
gen.num("h_max", h.max(), 3)
gen.num("sd_bound", s * np.sqrt(h.max()), 1)
gen.num("pop_on", X[104, 2], 1)
gen.num("pop_off", x_off[2], 1)
gen.num("share_off", 100 * sh_off[-1], 1)
gen.num("share_on", 100 * sh_on[-1], 1)
gen.num("corr_exact", corr_exact, 3)
gen.num("corr_b", corr_b, 3)
gen.num("se_dpi", fit.bse[1], 4)
gen.num("se_pop", fit.bse[2], 3)
gen.num("sd_avg", s * np.sqrt(p / n), 1)
gen.int("reps", reps)
gen.write()

# ---- figure -------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
ax.scatter(Bstar[:, 1], Bstar[:, 2], s=2, color=COLORS["accent"], alpha=0.35, linewidths=0)
ax.plot([fit.params[1]], [fit.params[2]], "o", color=COLORS["second"], ms=4)
ax.set_xlabel("income coefficient")
ax.set_ylabel("population coefficient")
ax.set_title("(a) coefficients")
ax = axes[1]
bins = np.linspace(-130, 130, 66)
ax.hist(Bstar @ x_on - x_on @ fit.params, bins=bins, color=COLORS["third"], alpha=0.8,
        label="1985Q1")
ax.hist(Bstar @ x_off - x_off @ fit.params, bins=bins, color=COLORS["second"], alpha=0.6,
        label="population 10% lower")
ax.set_xlabel("fitted mean minus its centre")
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.legend(frameon=False, loc="upper left", fontsize=7)
ax.set_title("(b) predictions")
fig.tight_layout()
fig.savefig(figure_path("ch26", "prediction_along_across"))
