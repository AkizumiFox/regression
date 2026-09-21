"""Chapter 44, Section 2: the additive model, backfitting and its Gauss-Seidel reading.

Simulated data (fixed seed): two continuous covariates whose dependence is a knob, two
smooth components, normal errors. The P-spline terms of Chapter 43 are the building blocks
(code/ch44/psplines.py). The script checks that

  * backfitting and one direct solve of the penalized normal equations agree,
  * the backfitting change decays geometrically, at a rate that worsens as the two
    covariates become more strongly dependent,
  * the effective degrees of freedom of the whole fit split over the components,

and draws the fitted components and the convergence curves.
"""
import matplotlib.pyplot as plt
import numpy as np

from psplines import spline_block
from regbook import COLORS, Generated, figure_path, use_book_style


# <<data>>
def simulate(n, rho, seed=20440):
    """Two dependent covariates on [0, 1] and an additive mean with normal errors."""
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n, 2)) @ np.array([[1.0, rho], [0.0, np.sqrt(1 - rho**2)]])
    z1, z2 = 0.5 * (1 + np.tanh(0.8 * u[:, 0])), 0.5 * (1 + np.tanh(0.8 * u[:, 1]))
    f1 = 16.0 * (z1 - 0.5) ** 3 - 0.6 * np.sin(4.5 * z1)
    f2 = 1.4 * np.exp(-25 * (z2 - 0.35) ** 2) - 0.5 * z2
    f1, f2 = f1 - f1.mean(), f2 - f2.mean()
    y = 3.0 + f1 + f2 + rng.normal(scale=0.4, size=n)
    return z1, z2, f1, f2, y
# <</data>>


# <<backfit>>
def backfit(y, blocks, lambdas, tol=1e-12, maxit=600):
    """Backfitting: cycle over the terms, refitting each to the current partial residual."""
    A = [Z.T @ Z + lam * K for (Z, K), lam in zip(blocks, lambdas)]
    f = [np.zeros_like(y) for _ in blocks]
    b0, history = y.mean(), []
    for _ in range(maxit):
        change = 0.0
        for j, (Z, _) in enumerate(blocks):
            r = y - b0 - sum(f[k] for k in range(len(blocks)) if k != j)
            new = Z @ np.linalg.solve(A[j], Z.T @ r)      # smooth the partial residual
            change = max(change, np.abs(new - f[j]).max())
            f[j] = new
        history.append(change)
        if change < tol:
            break
    return b0, f, np.array(history)


def direct_fit(y, blocks, lambdas):
    """One solve of the penalized normal equations, with an intercept column in front."""
    Z = np.column_stack([np.ones(len(y))] + [Zj for Zj, _ in blocks])
    K = np.zeros((Z.shape[1], Z.shape[1]))
    start = 1
    for (_, Kj), lam in zip(blocks, lambdas):
        d = Kj.shape[0]
        K[start:start + d, start:start + d] = lam * Kj
        start += d
    C = np.linalg.solve(Z.T @ Z + K, Z.T)                # gamma = C y
    edf = [np.trace(Zj @ C[s:s + Zj.shape[1]]) for Zj, s in
           zip([Zj for Zj, _ in blocks], np.cumsum([1] + [Zj.shape[1] for Zj, _ in blocks]))]
    return Z @ (C @ y), np.trace(Z @ C), edf
# <</backfit>>


if __name__ == "__main__":
    n, lam = 400, 8.0
    z1, z2, f1_true, f2_true, y = simulate(n, rho=0.5)
    blocks = [spline_block(z1), spline_block(z2)]
    lambdas = [lam, lam]

    # <<compare>>
    b0, f, history = backfit(y, blocks, lambdas)
    fitted, edf_total, edf = direct_fit(y, blocks, lambdas)
    print("cycles:", len(history))
    print("largest disagreement:", np.abs(fitted - (b0 + f[0] + f[1])).max())
    print(f"edf {edf_total:.3f} = 1 + {edf[0]:.3f} + {edf[1]:.3f}")
    # <</compare>>

    assert np.abs(fitted - (b0 + f[0] + f[1])).max() < 1e-9
    assert np.isclose(edf_total, 1.0 + edf[0] + edf[1])

    # geometric decay, slower when the covariates are more strongly dependent
    cycles, rates = {}, {}
    for rho in (0.0, 0.5, 0.9):
        zz1, zz2, _, _, yy = simulate(n, rho=rho)
        _, _, h = backfit(yy, [spline_block(zz1), spline_block(zz2)], lambdas)
        cycles[rho] = len(h)
        rates[rho] = (h[5] / h[1]) ** (1 / 4)             # the geometric rate of the decay
        assert np.all(np.diff(h[:len(h) - 1]) < 0)        # the change shrinks every cycle
    assert cycles[0.0] < cycles[0.5] < cycles[0.9]
    assert rates[0.0] < rates[0.5] < rates[0.9] < 1.0

    gen = Generated("ch44", "backfitting")
    gen.int("n", n)
    gen.int("nbasis", 15)
    gen.num("lam", lam, 1)
    gen.num("edf_total", edf_total, 2)
    gen.num("edf1", edf[0], 2)
    gen.num("edf2", edf[1], 2)
    gen.num("agree", np.abs(fitted - (b0 + f[0] + f[1])).max(), 1, sci=True)
    for rho in (0.0, 0.5, 0.9):
        tag = str(rho).replace(".", "")
        gen.int(f"cycles{tag}", cycles[rho])
        gen.num(f"rate{tag}", rates[rho], 3)
    gen.write()

    # ---- figure: components and convergence ---------------------------------
    use_book_style()
    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.2))
    for ax, z, ftrue, fhat, name, title in [
            (axes[0], z1, f1_true, f[0], "z_1", "(a) first component"),
            (axes[1], z2, f2_true, f[1], "z_2", "(b) second component")]:
        o = np.argsort(z)
        ax.scatter(z, y - y.mean(), s=4, color=COLORS["muted"], alpha=0.35, linewidths=0)
        ax.plot(z[o], ftrue[o], color=COLORS["second"], linestyle="--", label="truth")
        ax.plot(z[o], fhat[o], color=COLORS["accent"], label="backfitted")
        ax.set_xlabel(f"${name}$")
        ax.set_ylabel("component")
        ax.set_ylim(-2.2, 2.2)
        ax.set_title(title)
    axes[0].legend(frameon=False, loc="upper left")
    ax = axes[2]
    for rho, colour in zip((0.0, 0.5, 0.9), ("third", "accent", "second")):
        zz1, zz2, _, _, yy = simulate(n, rho=rho)
        _, _, h = backfit(yy, [spline_block(zz1), spline_block(zz2)], lambdas)
        ax.semilogy(np.arange(1, len(h) + 1), h, color=COLORS[colour], label=f"$\\rho={rho}$")
    ax.set_xlim(0, 40)
    ax.set_xlabel("backfitting cycle")
    ax.set_ylabel("largest change")
    ax.set_title("(c) convergence")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_path("ch44", "backfitting"))
