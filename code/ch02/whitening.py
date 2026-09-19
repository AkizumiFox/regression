"""Chapter 2, Section 2: standardizing a random vector (whitening).

A correlated bivariate sample is transformed by the symmetric inverse square root
Sigma^{-1/2} and by the inverse Cholesky factor L^{-1}. Both produce identity covariance;
they differ by an orthogonal matrix, and Mahalanobis distances are the same for both.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(2)

# <<whiten>>
mu = np.array([3.0, 1.0])
Sigma = np.array([[4.0, 2.4],
                  [2.4, 2.25]])                   # correlation 2.4 / (2 * 1.5) = 0.8

lam, U = np.linalg.eigh(Sigma)                    # spectral decomposition
W_sym = U @ np.diag(lam ** -0.5) @ U.T            # Sigma^{-1/2}
L = np.linalg.cholesky(Sigma)                     # Sigma = L L', L lower triangular
W_chol = np.linalg.inv(L)                         # L^{-1}

n = 400
Y = mu + rng.standard_normal((n, 2)) @ L.T        # any draws with Cov = Sigma will do
Z_sym = (Y - mu) @ W_sym.T
Z_chol = (Y - mu) @ W_chol.T

Q = W_chol @ np.linalg.inv(W_sym)                 # L^{-1} Sigma^{1/2}
print("Q'Q =", np.round(Q.T @ Q, 12))
# <</whiten>>

I2 = np.eye(2)
assert np.allclose(W_sym @ Sigma @ W_sym.T, I2)
assert np.allclose(W_chol @ Sigma @ W_chol.T, I2)
assert np.allclose(W_sym, W_sym.T)
assert np.allclose(Q.T @ Q, I2) and np.allclose(Z_chol, Z_sym @ Q.T)
# Mahalanobis distance is the squared length after either whitening
d2 = np.einsum("ij,jk,ik->i", Y - mu, np.linalg.inv(Sigma), Y - mu)
assert np.allclose(d2, (Z_sym**2).sum(1)) and np.allclose(d2, (Z_chol**2).sum(1))
# sample covariances of the whitened draws are close to I
assert np.allclose(np.cov(Z_sym, rowvar=False), I2, atol=0.2)
theta = np.degrees(np.arctan2(Q[1, 0], Q[0, 0]))
# Q is a rotation (not a reflection), through the angle reported in the text
assert np.isclose(np.linalg.det(Q), 1.0)
t = np.radians(theta)
assert np.allclose(Q, [[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
assert np.isclose(abs(theta), 22.5, atol=0.05)

gen = Generated("ch02", "whitening", prefix="wh")
gen.num("lam1", lam[1], 3)
gen.num("lam2", lam[0], 3)
for (i, j), v in np.ndenumerate(W_sym):
    gen.num(f"ws{i+1}{j+1}", v, 3)
for (i, j), v in np.ndenumerate(W_chol):
    gen.num(f"wc{i+1}{j+1}", v, 3)
for (i, j), v in np.ndenumerate(L):
    gen.num(f"L{i+1}{j+1}", v, 2)
gen.num("angle", abs(theta), 1)
gen.int("n", n)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.05),
                         gridspec_kw={"width_ratios": [15 / 11, 1, 1]})
mark = np.argsort(d2)[-3:]                        # three far-out points, tracked across panels
mcol = [COLORS["second"], COLORS["third"], COLORS["thread"]]
panels = [(Y, "(a) $\\mathbf{Y}$"), (Z_sym, "(b) $\\mathbf{\\Sigma}^{-1/2}(\\mathbf{Y}-\\mathbf{\\mu})$"),
          (Z_chol, "(c) $\\mathbf{L}^{-1}(\\mathbf{Y}-\\mathbf{\\mu})$")]
for k, (ax, (P, title)) in enumerate(zip(axes, panels)):
    ax.scatter(P[:, 0], P[:, 1], s=4, color=COLORS["accent"], alpha=0.45, linewidths=0)
    for c, i in zip(mcol, mark):
        ax.scatter(P[i, 0], P[i, 1], s=18, color=c, zorder=3, linewidths=0)
    for r in (1, 2, 3):                           # Mahalanobis contours d = r
        if k == 0:
            e = Ellipse(mu, 2 * r * np.sqrt(lam[1]), 2 * r * np.sqrt(lam[0]),
                        angle=np.degrees(np.arctan2(U[1, 1], U[0, 1])))
        else:
            e = Ellipse((0, 0), 2 * r, 2 * r)
        e.set(fill=False, edgecolor=COLORS["ink"], linewidth=0.6, zorder=2)
        ax.add_patch(e)
    ax.set_aspect("equal")
    ax.set_title(title)
    lim = 4.2 if k else None
    if k:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    else:
        ax.set_xlim(mu[0] - 7.5, mu[0] + 7.5)
        ax.set_ylim(mu[1] - 5.5, mu[1] + 5.5)
fig.tight_layout(w_pad=0.6)
fig.savefig(figure_path("ch02", "whitening"))
