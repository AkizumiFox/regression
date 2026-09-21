"""Chapter 44, Section 5: one structured additive fit with four kinds of term.

Simulated data (fixed seed) with a linear effect, a smooth effect of a continuous
covariate, a coefficient that varies smoothly with a second covariate, an areal spatial
effect on an 8 x 8 map, and a random intercept for 40 clusters. Every term is a pair
(design, penalty); the smoothing parameters come from the restricted-likelihood updates,
and the variance components they imply are compared with the values used to generate the
data. The script checks that

  * each term is recovered (the correlation with the truth and the error are reported),
  * the variance component of the random-intercept term is close to its true value,
  * the effective degrees of freedom add up over the terms,

and draws the four fitted terms against the truth.
"""
import matplotlib.pyplot as plt
import numpy as np

from psplines import markov_block, random_block, reml_lambdas, spline_term
from regbook import COLORS, Generated, figure_path, use_book_style
from spatial import lattice


# <<simulate>>
def simulate(n=1200, side=8, n_cluster=40, seed=44041):
    """A response built from a linear term, a smooth term, a varying coefficient,
    an areal spatial effect and a random intercept."""
    rng = np.random.default_rng(seed)
    coords, neighbours = lattice(side)
    d = side * side
    z1, z2 = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    u = rng.normal(size=n)                                   # the modified covariate
    x = rng.binomial(1, 0.4, n).astype(float)
    region = rng.integers(0, d, n)
    cluster = rng.integers(0, n_cluster, n)

    f1 = 1.5 * np.sin(2 * np.pi * z1) - 0.4
    g = 1.2 * (z2 - 0.5) ** 2 * 4 - 0.4                      # the varying coefficient
    a, b = coords[:, 0] / side, coords[:, 1] / side
    spatial = 1.2 * np.cos(3.0 * a + 1.4 * b) - 0.3 * b
    spatial -= spatial.mean()
    effect = rng.normal(scale=0.5, size=n_cluster)           # the random intercept
    effect -= effect.mean()

    truth = (2.0 + 0.8 * x + (f1 - f1.mean()) + (g - g.mean()) * u
             + spatial[region] + effect[cluster])
    y = truth + rng.normal(scale=0.7, size=n)
    return dict(y=y, z1=z1, z2=z2, u=u, x=x, region=region, cluster=cluster,
                neighbours=neighbours, f1=f1 - f1.mean(), g=g - g.mean(),
                spatial=spatial, effect=effect, d=d, n_cluster=n_cluster)
# <</simulate>>


if __name__ == "__main__":
    D = simulate()
    y, n = D["y"], len(D["y"])

    # <<fit>>
    Z1, K1, _ = spline_term(D["z1"], n_basis=15)                       # smooth effect
    Z2, K2, design2 = spline_term(D["z2"], n_basis=15, weight=D["u"])  # varying coefficient
    Z3, K3 = markov_block(D["region"], D["neighbours"], D["d"])        # areal spatial effect
    Z4, K4 = random_block(D["cluster"], D["n_cluster"])                # random intercept
    blocks = [(Z1, K1), (Z2, K2), (Z3, K3), (Z4, K4)]
    X = np.column_stack([np.ones(n), D["x"]])

    lambdas, fit, sigma2 = reml_lambdas(y, blocks, X)
    names = ["smooth", "varying", "spatial", "random"]
    for name, lam, e in zip(names, lambdas, fit.edf):
        print(f"{name:8s} lambda {lam:9.2f}   variance {sigma2 / lam:7.4f}   edf {e:6.2f}")
    print(f"residual sd {np.sqrt(sigma2):.3f}, total edf {fit.edf_total:.2f}, "
          f"linear effect {fit.beta[1]:.3f}")
    # <</fit>>

    assert np.isclose(fit.edf_total, X.shape[1] + sum(fit.edf))
    assert abs(np.sqrt(sigma2) - 0.7) < 0.05
    assert abs(sigma2 / lambdas[3] - 0.25) < 0.12            # the random-effect variance

    truth = [D["f1"], D["g"] * D["u"], D["spatial"][D["region"]], D["effect"][D["cluster"]]]
    corr = [float(np.corrcoef(t, p)[0, 1]) for t, p in zip(truth, fit.parts)]
    assert min(corr) > 0.8
    rmse = [float(np.sqrt(np.mean((t - p) ** 2))) for t, p in zip(truth, fit.parts)]

    gen = Generated("ch44", "star")
    gen.int("n", n)
    gen.int("d", D["d"])
    gen.int("nclust", D["n_cluster"])
    gen.int("ncoef", len(fit.gamma))
    gen.num("sigma", np.sqrt(sigma2), 3)
    gen.num("beta_x", fit.beta[1], 3)
    gen.num("edf_total", fit.edf_total, 2)
    for name, lam, e, c, r in zip(names, lambdas, fit.edf, corr, rmse):
        gen.num(f"lam_{name}", lam, 2)
        gen.num(f"var_{name}", sigma2 / lam, 4)
        gen.num(f"edf_{name}", e, 2)
        gen.num(f"corr_{name}", c, 3)
        gen.num(f"rmse_{name}", r, 3)
    gen.write()

    # ---- figure: the four fitted terms --------------------------------------
    use_book_style()
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 1.95))
    o = np.argsort(D["z1"])
    ax = axes[0]
    ax.plot(D["z1"][o], D["f1"][o], color=COLORS["second"], linestyle="--", label="truth")
    ax.plot(D["z1"][o], fit.parts[0][o], color=COLORS["accent"], label="fitted")
    ax.set_xlabel("$z_1$")
    ax.set_title("(a) smooth term", fontsize=8.5)
    ax.legend(frameon=False, fontsize=7)
    ax = axes[1]
    grid = np.linspace(0, 1, 100)
    coefs = design2(grid, np.ones(len(grid))) @ fit.gamma[X.shape[1] + Z1.shape[1]:
                                                          X.shape[1] + Z1.shape[1] + Z2.shape[1]]
    ax.plot(grid, np.interp(grid, D["z2"][np.argsort(D["z2"])],
                            D["g"][np.argsort(D["z2"])]),
            color=COLORS["second"], linestyle="--")
    ax.plot(grid, coefs, color=COLORS["accent"])
    ax.set_xlabel("$z_2$")
    ax.set_title("(b) varying coefficient", fontsize=8.5)
    ax = axes[2]
    hat_spatial = np.array([Z3[np.argmax(D["region"] == s)] @
                            fit.gamma[X.shape[1] + Z1.shape[1] + Z2.shape[1]:
                                      X.shape[1] + Z1.shape[1] + Z2.shape[1] + Z3.shape[1]]
                            for s in range(D["d"])])
    ax.scatter(D["spatial"], hat_spatial, s=8, color=COLORS["third"], linewidths=0)
    lims = [-1.4, 1.4]
    ax.plot(lims, lims, color=COLORS["muted"], linewidth=0.7)
    ax.set_xlabel("true district effect")
    ax.set_title("(c) spatial term", fontsize=8.5)
    ax = axes[3]
    hat_effect = np.array([Z4[np.argmax(D["cluster"] == k)] @ fit.gamma[-Z4.shape[1]:]
                           for k in range(D["n_cluster"])])
    ax.scatter(D["effect"], hat_effect, s=8, color=COLORS["thread"], linewidths=0)
    lims = [-1.1, 1.1]
    ax.plot(lims, lims, color=COLORS["muted"], linewidth=0.7)
    ax.set_xlabel("true cluster effect")
    ax.set_title("(d) random intercept", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(figure_path("ch44", "star_terms"))
