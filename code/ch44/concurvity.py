"""Chapter 44, Section 2: concurvity, the additive analogue of collinearity.

The second covariate is a noisy function of the first, and the noise level is a knob: as
it shrinks, the two model spaces come close to sharing a direction. For a fixed design the
fitted components are linear in the response, so the exact covariance of each component is
available and no simulation is needed for the variances. The script

  * computes the worst-case concurvity index (the largest squared canonical correlation
    between the two centred model spaces),
  * computes the exact pointwise standard deviation of the first fitted component,
  * checks that it grows as the dependence tightens while the standard deviation of the
    fitted mean barely moves,

and draws replicate fits under weak and strong concurvity.
"""
import matplotlib.pyplot as plt
import numpy as np

from psplines import Fit, spline_block
from regbook import COLORS, Generated, figure_path, use_book_style


# <<setup>>
def design(n, noise, seed=44012):
    """z2 is a noisy monotone function of z1: small noise means strong concurvity."""
    rng = np.random.default_rng(seed)
    z1 = rng.uniform(0, 1, n)
    z2 = np.clip(z1**1.5 + noise * rng.normal(size=n), 0.001, 0.999)
    return z1, z2


def concurvity(Z1, Z2):
    """The largest squared canonical correlation between the two column spaces."""
    Q1, _ = np.linalg.qr(Z1)
    Q2, _ = np.linalg.qr(Z2)
    return np.linalg.svd(Q1.T @ Q2, compute_uv=False)[0] ** 2
# <</setup>>


if __name__ == "__main__":
    n, lam, sigma = 300, 5.0, 0.4
    noises = (0.30, 0.10, 0.03)
    kappa, sd_comp, sd_fit, edf = {}, {}, {}, {}
    target = 0.3                                      # the point at which the sd is reported

    # <<variance>>
    for noise in noises:
        z1, z2 = design(n, noise)
        blocks = [spline_block(z1), spline_block(z2)]
        fit = Fit(np.zeros(n), blocks, np.ones((n, 1)), [lam, lam])
        i = int(np.argmin(np.abs(z1 - target)))       # an observation near z1 = 0.3
        kappa[noise] = concurvity(blocks[0][0], blocks[1][0])
        sd_comp[noise] = sigma * np.linalg.norm(fit.comp[0][i])   # sd of the component
        sd_fit[noise] = sigma * np.linalg.norm(fit.S[i])          # sd of the fitted mean
        edf[noise] = fit.edf_total
        print(f"noise {noise}: concurvity {kappa[noise]:.3f}, "
              f"sd of component {sd_comp[noise]:.3f}, sd of fit {sd_fit[noise]:.3f}")
    # <</variance>>

    assert kappa[0.30] < kappa[0.10] < kappa[0.03] < 1.0
    assert sd_comp[0.30] < sd_comp[0.10] < sd_comp[0.03]
    assert sd_fit[0.03] < 1.6 * sd_fit[0.30]          # the fitted mean is far less affected
    inflation = sd_comp[0.03] / sd_comp[0.30]
    fit_inflation = sd_fit[0.03] / sd_fit[0.30]

    gen = Generated("ch44", "concurvity")
    gen.int("n", n)
    gen.num("lam", lam, 1)
    gen.num("sigma", sigma, 1)
    for noise in noises:
        tag = str(noise).replace(".", "")
        gen.num(f"kappa{tag}", kappa[noise], 3)
        gen.num(f"sd{tag}", sd_comp[noise], 3)
        gen.num(f"sdfit{tag}", sd_fit[noise], 3)
        gen.num(f"edf{tag}", edf[noise], 2)
    gen.num("inflation", inflation, 1)
    gen.num("fitinflation", fit_inflation, 2)
    gen.write()

    # ---- figure: replicate component fits under weak and strong concurvity ----
    use_book_style()
    rng = np.random.default_rng(7788)
    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.2))
    for ax, noise, title in [(axes[0], 0.30, "(a) weak concurvity"),
                             (axes[1], 0.03, "(b) strong concurvity")]:
        z1, z2 = design(n, noise)
        blocks = [spline_block(z1), spline_block(z2)]
        fit = Fit(np.zeros(n), blocks, np.ones((n, 1)), [lam, lam])
        f1 = 1.2 * np.sin(3.0 * z1) - 0.9 * z1
        f2 = 4.8 * (z2 - 0.5) ** 2
        f1, f2 = f1 - f1.mean(), f2 - f2.mean()
        o = np.argsort(z1)
        for _ in range(25):
            y = 1.0 + f1 + f2 + sigma * rng.normal(size=n)
            ax.plot(z1[o], (fit.comp[0] @ y)[o], color=COLORS["accent"],
                    alpha=0.25, linewidth=0.7)
        ax.plot(z1[o], f1[o], color=COLORS["second"], linestyle="--")
        ax.set_xlabel("$z_1$")
        ax.set_ylabel("first component")
        ax.set_ylim(-2.6, 2.6)
        ax.set_title(title)
    ax = axes[2]
    ks, ss = [], []
    for noise in np.geomspace(0.01, 0.6, 24):
        z1, z2 = design(n, noise)
        blocks = [spline_block(z1), spline_block(z2)]
        fit = Fit(np.zeros(n), blocks, np.ones((n, 1)), [lam, lam])
        ks.append(concurvity(blocks[0][0], blocks[1][0]))
        j = int(np.argmin(np.abs(z1 - target)))
        ss.append(sigma * np.linalg.norm(fit.comp[0][j]))
    ax.plot(ks, ss, color=COLORS["third"], marker="o", markersize=2.5)
    ax.set_xlabel("concurvity index")
    ax.set_ylabel("sd of the component")
    ax.set_title("(c) the price of concurvity")
    fig.tight_layout()
    fig.savefig(figure_path("ch44", "concurvity"))
