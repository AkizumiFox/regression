"""Chapter 44, Section 6: a generalized additive model for physician visits.

Chapter 37 fitted a Poisson log-linear model to the RAND Health Insurance Experiment
(statsmodels.datasets.randhie, public domain, 20190 person-years) with every regressor
entering linearly. Here the disease index and the log participation incentive payment get
smooth terms instead, fitted by penalized IRLS, with the smoothing parameters chosen by
the restricted-likelihood updates of Section 44.5 applied to the working problem and with
the estimated dispersion of Chapter 38 in place of the nominal one.

The script reports the effective degrees of freedom of the two smooths, the dispersion,
the quasi-F comparison with the log-linear fit, and draws the two estimated curves with
bands and two checks of the fitted model.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from psplines import irls_reml, penalized_irls, poisson_deviance, spline_term
from regbook import COLORS, Generated, figure_path, use_book_style

# <<fit>>
import numpy as np
import statsmodels.api as sm

data = sm.datasets.randhie.load_pandas().data          # 20190 person-years, public domain
y = data["mdvis"].to_numpy(float)                      # outpatient physician visits in the year
disea = data["disea"].to_numpy(float)                  # a general disease index
lpi = data["lpi"].to_numpy(float)                      # log participation incentive payment
linear = ["lncoins", "idp", "physlm", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in linear])

Z1, K1, design1 = spline_term(disea, n_basis=12)       # the two smooth terms
Z2, K2, design2 = spline_term(lpi, n_basis=12)
lambdas, fit = irls_reml(y, [(Z1, K1), (Z2, K2)], X, family="poisson", dispersion=True)

print("smoothing parameters", lambdas.round(2))
print(f"edf: disease index {fit.edf[0]:.2f}, incentive payment {fit.edf[1]:.2f}, "
      f"model {fit.edf_total:.2f}")
print(f"dispersion {fit.phi:.2f}")
for name, b in zip(linear, fit.beta[1:]):
    print(f"  {name:8s} {b:8.4f}   rate ratio {np.exp(b):.4f}")
# <</fit>>

# <<compare>>
X_linear = np.column_stack([X, disea, lpi])            # the same regressors, both smooths linear
fit_linear = penalized_irls(y, [], X_linear, [], family="poisson")
D_gam = poisson_deviance(y, fit.mu)
D_linear = poisson_deviance(y, fit_linear.mu)
ddf = fit.edf_total - fit_linear.edf_total
F = (D_linear - D_gam) / ddf / fit.phi                 # a quasi-F, not an exact null statistic
print(f"deviance {D_linear:.0f} (linear) against {D_gam:.0f} (smooth)")
print(f"quasi-F {F:.2f} on {ddf:.2f} and {len(y) - fit.edf_total:.0f} degrees of freedom")
# <</compare>>

if __name__ == "__main__":
    n = len(y)
    assert D_gam < D_linear and ddf > 0
    assert fit.phi > 5.0                               # heavily overdispersed, as in Chapter 37
    assert np.abs(X.T @ (y - fit.mu)).max() < 1e-4     # the unpenalized score equations hold
    p_value = stats.f.sf(F, ddf, n - fit.edf_total)
    assert p_value < 1e-10

    cov = fit.phi * np.linalg.inv(fit.A)               # the Bayesian covariance, scaled by phi
    grids, curves, ses = {}, {}, {}
    for name, x, design, sl in [("disea", disea, design1, slice(X.shape[1],
                                                               X.shape[1] + Z1.shape[1])),
                                ("lpi", lpi, design2, slice(X.shape[1] + Z1.shape[1],
                                                            X.shape[1] + Z1.shape[1] +
                                                            Z2.shape[1]))]:
        g = np.linspace(x.min(), x.max(), 200)
        rows = np.zeros((len(g), len(fit.gamma)))
        rows[:, sl] = design(g)
        grids[name], curves[name] = g, rows @ fit.gamma
        ses[name] = np.sqrt(np.einsum("ij,jk,ik->i", rows, cov, rows))
    assert np.all(ses["disea"] > 0)

    # the curve is least certain where the data are thinnest
    q = np.quantile(disea, [0.05, 0.95])
    se_mid = float(np.interp(np.median(disea), grids["disea"], ses["disea"]))
    se_top = float(np.interp(disea.max(), grids["disea"], ses["disea"]))
    assert se_top > se_mid

    gen = Generated("ch44", "gam")
    gen.int("n", n)
    gen.int("nbasis", 12)
    gen.num("lam1", lambdas[0], 2)
    gen.num("lam2", lambdas[1], 2)
    gen.num("edf1", fit.edf[0], 2)
    gen.num("edf2", fit.edf[1], 2)
    gen.num("edf_total", fit.edf_total, 2)
    gen.num("phi", fit.phi, 2)
    gen.num("dev_gam", D_gam, 0)
    gen.num("dev_linear", D_linear, 0)
    gen.num("ddf", ddf, 2)
    gen.num("F", F, 2)
    gen.num("p", p_value, 1, sci=True)
    gen.num("beta_coins", fit.beta[1], 4)
    gen.num("rr_coins", np.exp(fit.beta[1]), 4)
    gen.num("beta_physlm", fit.beta[3], 4)
    gen.num("beta_disea_linear", fit_linear.beta[-2], 4)
    gen.num("beta_coins_linear", fit_linear.beta[1], 4)
    # Chapter 37's own model: the same list without the participation incentive payment
    fit_ch37 = penalized_irls(y, [], np.column_stack([X, disea]), [], family="poisson")
    gen.num("beta_coins_ch37", fit_ch37.beta[1], 4)
    gen.num("se_mid", se_mid, 3)
    gen.num("se_top", se_top, 3)
    gen.num("range_disea", float(curves["disea"].max() - curves["disea"].min()), 2)
    gen.int("n_top", int(np.sum(disea > 40)))
    gen.num("pct_top", 100 * float(np.mean(disea > 40)), 1)
    gen.num("disea_max", float(disea.max()), 1)
    gen.write()

    # ---- figure: the two curves and two checks ------------------------------
    use_book_style()
    fig, axes = plt.subplots(1, 4, figsize=(6.8, 1.95))
    for ax, name, label, title in [(axes[0], "disea", "disease index", "(a) smooth of disease"),
                                   (axes[1], "lpi", "log incentive payment",
                                    "(b) smooth of incentive")]:
        g, c, s = grids[name], curves[name], ses[name]
        ax.fill_between(g, c - 2 * s, c + 2 * s, color=COLORS["accent"], alpha=0.18, linewidth=0)
        ax.plot(g, c, color=COLORS["accent"])
        ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
        ax.plot(np.quantile(disea if name == "disea" else lpi, np.linspace(0, 1, 30)),
                np.full(30, -1.3), "|", color=COLORS["muted"], markersize=3)
        ax.set_xlabel(label)
        ax.set_ylabel("effect on log mean")
        ax.set_ylim(-1.5, 1.6)
        ax.set_title(title, fontsize=8.5)

    edges = np.quantile(fit.mu, np.linspace(0, 1, 21))
    bin_id = np.clip(np.searchsorted(edges, fit.mu) - 1, 0, 19)
    means = np.array([fit.mu[bin_id == b].mean() for b in range(20)])
    obs = np.array([y[bin_id == b].mean() for b in range(20)])
    var = np.array([y[bin_id == b].var() for b in range(20)])
    ax = axes[2]
    ax.scatter(means, obs, s=10, color=COLORS["third"], linewidths=0)
    lims = [means.min(), means.max()]
    ax.plot(lims, lims, color=COLORS["muted"], linewidth=0.7)
    ax.set_xlabel("fitted mean")
    ax.set_ylabel("observed mean")
    ax.set_title("(c) calibration", fontsize=8.5)
    ax = axes[3]
    ax.scatter(means, var, s=10, color=COLORS["second"], linewidths=0)
    ax.plot(lims, lims, color=COLORS["muted"], linewidth=0.7, label="Poisson")
    ax.plot(lims, fit.phi * np.array(lims), color=COLORS["accent"], linewidth=0.9,
            label="$\\hat\\phi\\,\\mu$")
    ax.set_xlabel("fitted mean")
    ax.set_ylabel("observed variance")
    ax.set_title("(d) mean and variance", fontsize=8.5)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(figure_path("ch44", "gam_visits"))
