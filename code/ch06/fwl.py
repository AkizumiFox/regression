"""Chapter 6, Section 6: the Frisch-Waugh-Lovell theorem on the 2009 US state data.

Response: murder rate per 100,000. Regressors: poverty rate, percentage of
single-parent households, percentage urban. Public-domain data shipped with
statsmodels (statsmodels.datasets.statecrime).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")        # an extreme point; see the exercises of Section 6.6

# <<fwl>>
y = data["murder"].to_numpy()
x1 = data["poverty"].to_numpy()                       # the coefficient we care about
X2 = np.column_stack([np.ones(len(y)), data["single"], data["urban"]])
X = np.column_stack([x1, X2])


def resid(v, Z):
    """Residual of v after projecting onto C(Z)."""
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef


beta_full, *_ = np.linalg.lstsq(X, y, rcond=None)     # regress y on everything

y_tilde = resid(y, X2)                                # (I - M2) y
x_tilde = resid(x1, X2)                               # (I - M2) x1
beta_fwl = (x_tilde @ y_tilde) / (x_tilde @ x_tilde)  # simple regression, no intercept

print(f"full regression:  {beta_full[0]:.6f}")
print(f"FWL regression:   {beta_fwl:.6f}")
# <</fwl>>

assert np.isclose(beta_full[0], beta_fwl, rtol=1e-10)
# residuals of the two regressions coincide
assert np.allclose(y - X @ beta_full, y_tilde - beta_fwl * x_tilde, atol=1e-10)

n, p = X.shape
rss = np.sum((y - X @ beta_full) ** 2)
se_full = np.sqrt(rss / (n - p) * np.linalg.inv(X.T @ X)[0, 0])
se_naive = np.sqrt(rss / (n - 1) / (x_tilde @ x_tilde))
assert np.isclose(se_naive * np.sqrt((n - 1) / (n - p)), se_full)
fit_simple = np.polyfit(x1, y, 1)

gen = Generated("ch06", "fwl")
gen.int("n", n)
gen.int("p", p)
gen.num("beta", beta_full[0], 4)
gen.num("beta_single", beta_full[2], 4)
gen.num("beta_urban", beta_full[3], 4)
gen.num("beta_simple", fit_simple[0], 4)
gen.num("se_full", se_full, 4)
gen.num("se_naive", se_naive, 4)
gen.num("se_ratio", se_full / se_naive, 4)
gen.write()

# ---- added-variable plot vs marginal plot -----------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
ax.scatter(x1, y, s=10, color=COLORS["accent"], alpha=0.8, linewidths=0)
xs = np.linspace(x1.min(), x1.max(), 2)
ax.plot(xs, np.polyval(fit_simple, xs), color=COLORS["second"])
ax.set_xlabel("poverty (%)")
ax.set_ylabel("murder rate")
ax.set_title("(a) marginal: slope %.3f" % fit_simple[0])
ax = axes[1]
ax.scatter(x_tilde, y_tilde, s=10, color=COLORS["accent"], alpha=0.8, linewidths=0)
xs = np.linspace(x_tilde.min(), x_tilde.max(), 2)
ax.plot(xs, beta_fwl * xs, color=COLORS["second"])
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.axvline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"poverty, residualized on single, urban")
ax.set_ylabel(r"murder, residualized")
ax.set_title("(b) added-variable: slope %.3f" % beta_fwl)
fig.tight_layout()
fig.savefig(figure_path("ch06", "fwl_added_variable"))
