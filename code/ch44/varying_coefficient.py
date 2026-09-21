"""Chapter 44, Section 3: a varying-coefficient model on the 2009 US state data.

Chapter 6 fitted the murder rate on the poverty rate, the single-parent percentage and
the urban percentage, and reported one number for the coefficient of poverty. Here that
coefficient is allowed to vary smoothly with how urban the state is:

    murder = beta_0 + beta_1 single + {beta_2 + f(urban)} (poverty - mean) + f_0(urban) + e,

a P-spline term multiplied by the poverty column. The smoothing parameters come from the
restricted-likelihood updates of Section 44.5. Public-domain data shipped with statsmodels
(statsmodels.datasets.statecrime).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from psplines import Fit, reml_lambdas, spline_term
from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
import numpy as np
import statsmodels.api as sm

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
urban = data["urban"].to_numpy()                      # per cent living in urban areas
poverty = data["poverty"].to_numpy()
single = data["single"].to_numpy()
pov_c = poverty - poverty.mean()                      # centred, so beta_2 is the average slope

Z0, K0, design0 = spline_term(urban, n_basis=8)                       # smooth main effect
Z1, K1, design1 = spline_term(urban, n_basis=8, weight=pov_c)         # varying coefficient
X = np.column_stack([np.ones(len(y)), single, pov_c])
lambdas, fit, sigma2 = reml_lambdas(y, [(Z0, K0), (Z1, K1)], X)

print("smoothing parameters", lambdas.round(2))
print(f"edf: total {fit.edf_total:.2f}, "
      f"smooth main effect {fit.edf[0]:.2f}, varying coefficient {fit.edf[1]:.2f}")
print("average poverty slope", round(fit.beta[2], 4))
# <</fit>>

# <<slope>>
grid = np.linspace(urban.min(), urban.max(), 200)
row = np.zeros((len(grid), len(fit.gamma)))
row[:, 2] = 1.0                                       # the constant part beta_2 ...
row[:, X.shape[1] + Z0.shape[1]:] = design1(grid, np.ones(len(grid)))   # ... plus f(urban)
slope = row @ fit.gamma
cov = sigma2 * np.linalg.inv(fit.A)                   # the Bayesian covariance of Section 44.5
se = np.sqrt(np.einsum("ij,jk,ik->i", row, cov, row))
print(f"slope at 20 per cent urban {np.interp(20, grid, slope):.3f}, "
      f"at 60 per cent {np.interp(60, grid, slope):.3f}, "
      f"at 90 per cent {np.interp(90, grid, slope):.3f}")
# <</slope>>

if __name__ == "__main__":
    n = len(y)
    # the like-for-like comparison: the same smooth main effect of urbanization, but a
    # constant poverty slope. (Chapter 6's fit, with urbanization entering linearly, is
    # the third number below.)
    _, fit_const, _ = reml_lambdas(y, [(Z0, K0)], X)
    b_const = fit_const.beta[2]
    X_ch6 = np.column_stack([np.ones(n), single, urban, pov_c])
    b_ch6 = np.linalg.lstsq(X_ch6, y, rcond=None)[0][3]
    rss_const = np.sum((y - fit_const.fitted) ** 2)
    rss_vary = np.sum((y - fit.fitted) ** 2)
    assert rss_vary < rss_const
    assert fit.edf[1] > 2.0                            # the coefficient function is not a line
    lo, hi = slope - 2 * se, slope + 2 * se
    covers_zero = np.mean(lo <= 0)                     # fraction of the range whose band holds 0
    assert np.all(se > 0)
    assert np.isclose(np.sum(fit.parts[1] / pov_c), 0.0, atol=1e-7)    # sum_i f(urban_i) = 0

    # the difference between the coefficient at 60 and at 20 per cent urban
    contrast = np.zeros(len(fit.gamma))
    contrast[X.shape[1] + Z0.shape[1]:] = (design1([60.0], [1.0]) - design1([20.0], [1.0]))[0]
    diff = contrast @ fit.gamma
    se_diff = np.sqrt(contrast @ cov @ contrast)
    print(f"slope(60) - slope(20) = {diff:.3f} with standard error {se_diff:.3f}")
    assert abs(diff) < 2.5 * se_diff                   # not resolved by 50 states

    gen = Generated("ch44", "varying")
    gen.int("n", n)
    gen.num("lam0", lambdas[0], 1)
    gen.num("lam1", lambdas[1], 1)
    gen.num("edf_total", fit.edf_total, 2)
    gen.num("edf0", fit.edf[0], 2)
    gen.num("edf1", fit.edf[1], 2)
    gen.num("beta_const", b_const, 3)
    gen.num("beta_ch6", b_ch6, 4)
    gen.num("beta_avg", fit.beta[2], 3)
    gen.num("single", fit.beta[1], 3)
    gen.num("slope20", float(np.interp(20, grid, slope)), 3)
    gen.num("slope60", float(np.interp(60, grid, slope)), 3)
    gen.num("slope90", float(np.interp(90, grid, slope)), 3)
    gen.num("se60", float(np.interp(60, grid, se)), 3)
    gen.num("covers", 100 * covers_zero, 0)
    gen.num("diff", diff, 3)
    gen.num("se_diff", se_diff, 3)
    gen.num("ratio", diff / se_diff, 2)
    gen.num("sigma", np.sqrt(sigma2), 3)
    gen.write()

    # ---- figure -------------------------------------------------------------
    use_book_style()
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
    ax = axes[0]
    ax.fill_between(grid, lo, hi, color=COLORS["accent"], alpha=0.18, linewidth=0)
    ax.plot(grid, slope, color=COLORS["accent"])
    ax.axhline(b_const, color=COLORS["second"], linestyle="--",
               label="constant slope, smooth urban")
    ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
    ax.set_xlabel("urban (%)")
    ax.set_ylabel("coefficient of poverty")
    ax.set_title("(a) the varying coefficient")
    ax.legend(frameon=False, loc="lower right")
    ax = axes[1]
    ax.scatter(poverty, y, s=10, color=COLORS["muted"], alpha=0.6, linewidths=0)
    povs = np.linspace(poverty.min(), poverty.max(), 50)
    for level, colour in zip((25.0, 55.0, 90.0), ("third", "accent", "thread")):
        rows = np.zeros((len(povs), len(fit.gamma)))
        rows[:, 0], rows[:, 1] = 1.0, single.mean()
        rows[:, 2] = povs - poverty.mean()
        rows[:, X.shape[1]:X.shape[1] + Z0.shape[1]] = design0(np.full(len(povs), level))
        rows[:, X.shape[1] + Z0.shape[1]:] = design1(np.full(len(povs), level),
                                                     povs - poverty.mean())
        ax.plot(povs, rows @ fit.gamma, color=COLORS[colour], label=f"urban {level:.0f}%")
    ax.set_xlabel("poverty (%)")
    ax.set_ylabel("murder rate")
    ax.set_title("(b) what the fit says")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(figure_path("ch44", "varying_coefficient"))
