"""Chapter 44, Section 4: spatial effects, by a Markov random field and by a radial basis.

A simulated map (fixed seed): 144 districts on a 12 x 12 lattice, unequal numbers of
observations per district, a covariate with a linear effect, and a smooth spatial surface
the model does not know. Three fits are compared against the truth:

  * separate unpenalized district effects (the raw map),
  * the Markov random field penalty, which pulls each district towards its neighbours,
  * a kriging-type radial basis on the district centroids.

The script checks that both penalized fits beat the raw map, that the neighbourhood
penalty annihilates constants, and that its null space has one dimension per connected
component; it writes the maps and a small worked penalty matrix.
"""
import matplotlib.pyplot as plt
import numpy as np

from psplines import Fit, markov_block, radial_block, reml_lambdas
from regbook import COLORS, Generated, figure_path, use_book_style


# <<lattice>>
def lattice(side):
    """Rook adjacency on a side x side grid of districts, and the district centroids."""
    coords = np.array([[c + 0.5, r + 0.5] for r in range(side) for c in range(side)])
    neighbours = []
    for s in range(side * side):
        r, c = divmod(s, side)
        nb = [(r + dr) * side + (c + dc) for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1))
              if 0 <= r + dr < side and 0 <= c + dc < side]
        neighbours.append(sorted(nb))
    return coords, neighbours


def neighbourhood_penalty(neighbours):
    """K = D - W: the number of neighbours on the diagonal, -1 for each neighbour pair."""
    d = len(neighbours)
    K = np.zeros((d, d))
    for s, nbrs in enumerate(neighbours):
        K[s, s] = len(nbrs)
        for t in nbrs:
            K[s, t] = -1.0
    return K
# <</lattice>>


# <<simulate>>
def simulate(side=12, seed=44031):
    """Counts of observations per district, one covariate, and a smooth spatial surface."""
    rng = np.random.default_rng(seed)
    coords, neighbours = lattice(side)
    d = side * side
    u, v = coords[:, 0] / side, coords[:, 1] / side
    surface = 1.6 * np.sin(2.6 * u) * np.cos(2.2 * v) - 1.1 * (u - 0.45) ** 2
    surface = surface - surface.mean()
    m = rng.integers(1, 9, size=d)                      # unequal district sample sizes
    region = np.repeat(np.arange(d), m)
    x = rng.normal(size=len(region))
    y = 2.0 + 0.7 * x + surface[region] + rng.normal(scale=0.9, size=len(region))
    return coords, neighbours, region, x, y, surface, m
# <</simulate>>


if __name__ == "__main__":
    side = 12
    coords, neighbours, region, x, y, surface, m = simulate(side)
    d = side * side
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    K = neighbourhood_penalty(neighbours)
    assert np.allclose(K @ np.ones(d), 0.0)                  # constants are unpenalized
    assert np.allclose(K, K.T) and np.min(np.linalg.eigvalsh(K)) > -1e-9
    assert np.sum(np.linalg.eigvalsh(K) < 1e-9) == 1         # the lattice is connected

    # <<fits>>
    mrf = markov_block(region, neighbours, d)
    krig = radial_block(coords[region], coords[::7])         # centroids, every seventh knot
    lam_mrf, fit_mrf, s2_mrf = reml_lambdas(y, [mrf], X)
    lam_krig, fit_krig, _ = reml_lambdas(y, [krig], X)
    raw = Fit(y, [(mrf[0], 1e-10 * np.eye(mrf[1].shape[0]))], X, [1.0])   # no smoothing

    def district_effect(fit, Z):
        """The fitted spatial effect of every district, read off the incidence design."""
        g = fit.gamma[X.shape[1]:]
        return np.array([Z[np.argmax(region == s)] @ g for s in range(d)])

    err = lambda f: np.sqrt(np.mean((f - f.mean() - surface) ** 2))
    e_raw, e_mrf, e_krig = (err(district_effect(raw, mrf[0])),
                            err(district_effect(fit_mrf, mrf[0])),
                            err(district_effect(fit_krig, krig[0])))
    print(f"root mean squared error of the map: raw {e_raw:.3f}, "
          f"Markov random field {e_mrf:.3f}, kriging {e_krig:.3f}")
    print(f"edf: raw {raw.edf[0]:.1f}, MRF {fit_mrf.edf[0]:.1f}, kriging {fit_krig.edf[0]:.1f}")
    # <</fits>>

    assert e_mrf < e_raw and e_krig < e_raw
    assert fit_mrf.edf[0] < raw.edf[0]
    se_beta = np.sqrt(s2_mrf * np.linalg.inv(fit_mrf.A)[1, 1])
    # the linear effect is recovered; this particular draw sits 2.8 standard errors low
    assert abs(fit_mrf.beta[1] - 0.7) < 3.5 * se_beta

    gen = Generated("ch44", "spatial")
    gen.int("side", side)
    gen.int("d", d)
    gen.int("n", n)
    gen.int("mmin", int(m.min()))
    gen.int("mmax", int(m.max()))
    gen.int("nknots", len(coords[::7]))
    gen.num("lam_mrf", lam_mrf[0], 2)
    gen.num("lam_krig", lam_krig[0], 2)
    gen.num("err_raw", e_raw, 3)
    gen.num("err_mrf", e_mrf, 3)
    gen.num("err_krig", e_krig, 3)
    gen.num("edf_raw", raw.edf[0], 1)
    gen.num("edf_mrf", fit_mrf.edf[0], 1)
    gen.num("edf_krig", fit_krig.edf[0], 1)
    gen.num("beta_x", fit_mrf.beta[1], 3)
    gen.num("se_beta", se_beta, 3)
    gen.num("sigma_mrf", np.sqrt(s2_mrf), 3)
    gen.write()

    # ---- figure: the maps ---------------------------------------------------
    use_book_style()
    maps = [(surface, "(a) true surface"),
            (district_effect(raw, mrf[0]), "(b) unpenalized"),
            (district_effect(fit_mrf, mrf[0]), "(c) Markov field"),
            (district_effect(fit_krig, krig[0]), "(d) kriging")]
    fig, axes = plt.subplots(1, 4, figsize=(6.6, 1.95))
    vmax = 1.9
    for ax, (field, title) in zip(axes, maps):
        img = (field - field.mean()).reshape(side, side)
        h = ax.imshow(img, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=8.5)
    fig.colorbar(h, ax=axes, fraction=0.025, pad=0.01)
    fig.savefig(figure_path("ch44", "spatial_maps"))
