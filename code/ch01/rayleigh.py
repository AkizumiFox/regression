"""Chapter 1, Section 7: the Rayleigh quotient of a 2 x 2 symmetric matrix.

The quotient x'Ax / x'x depends only on the direction of x. Along the unit circle it
oscillates between the smallest and largest eigenvalues, reaching them at the eigenvectors.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<rayleigh>>
A = np.array([[2.5, 1.0],
              [1.0, 1.5]])
lam, Q = np.linalg.eigh(A)                     # ascending: lam[0] <= lam[1]
theta = np.linspace(0, np.pi, 721)
U = np.vstack([np.cos(theta), np.sin(theta)])  # unit vectors u(theta), as columns
R = np.einsum("it,ij,jt->t", U, A, U)          # Rayleigh quotient u'Au

print("eigenvalues:", lam)
print("range of the quotient on the grid:", R.min(), R.max())
# <</rayleigh>>

assert lam[0] - 1e-12 <= R.min() and R.max() <= lam[1] + 1e-12
assert np.isclose(R.max(), lam[1], atol=1e-4) and np.isclose(R.min(), lam[0], atol=1e-4)
angle = lambda v: np.mod(np.arctan2(v[1], v[0]), np.pi)
a_max, a_min = angle(Q[:, 1]), angle(Q[:, 0])
assert abs(theta[np.argmax(R)] - a_max) < 0.01
assert abs(theta[np.argmin(R)] - a_min) < 0.01
assert np.isclose(abs(a_max - a_min), np.pi / 2)                  # eigenvectors orthogonal
# the quotient is a weighted average of eigenvalues: R = lam1 c1^2 + lam2 c2^2
C = Q.T @ U
assert np.allclose(R, lam[0] * C[0] ** 2 + lam[1] * C[1] ** 2)
# spectral decomposition and the trace/determinant identities
assert np.allclose(Q @ np.diag(lam) @ Q.T, A)
assert np.isclose(lam.sum(), np.trace(A)) and np.isclose(lam.prod(), np.linalg.det(A))

gen = Generated("ch01", "rayleigh", prefix="ray")
gen.num("lmax", lam[1], 3)
gen.num("lmin", lam[0], 3)
gen.num("amax", np.degrees(a_max), 1)
gen.num("amin", np.degrees(a_min), 1)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5), gridspec_kw={"width_ratios": [1, 1.35]})
ax = axes[0]
t = np.linspace(0, 2 * np.pi, 400)
circle = np.vstack([np.cos(t), np.sin(t)])
level = Q @ np.diag(1 / np.sqrt(lam)) @ Q.T @ circle              # {x : x'Ax = 1}
ax.plot(circle[0], circle[1], color=COLORS["grid"], linewidth=1.0)
ax.plot(level[0], level[1], color=COLORS["accent"])
for j, col, name in [(1, COLORS["second"], r"$\bf q_1$"), (0, COLORS["third"], r"$\bf q_2$")]:
    q = Q[:, j] * np.sign(Q[0, j] if abs(Q[0, j]) > 1e-9 else 1) * (1 if j == 1 else -1)
    ax.annotate("", xy=q, xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=col, linewidth=1.1))
    ax.text(*(1.18 * q), name, color=col, ha="center", va="center")
ax.set_aspect("equal")
ax.set_xlim(-1.35, 1.35)
ax.set_ylim(-1.35, 1.35)
ax.set_xticks([-1, 0, 1])
ax.set_yticks([-1, 0, 1])
ax.set_title(r"(a) $\mathbf{x}^{\top}\mathbf{A}\mathbf{x}=1$ and the unit circle")
ax = axes[1]
deg = np.degrees(theta)
ax.plot(deg, R, color=COLORS["accent"])
ax.axhline(lam[1], color=COLORS["second"], linestyle="--", linewidth=0.8)
ax.axhline(lam[0], color=COLORS["third"], linestyle="--", linewidth=0.8)
ax.axvline(np.degrees(a_max), color=COLORS["second"], linewidth=0.6, alpha=0.6)
ax.axvline(np.degrees(a_min), color=COLORS["third"], linewidth=0.6, alpha=0.6)
ax.text(181, lam[1], r"$\lambda_1$", color=COLORS["second"], va="center")
ax.text(181, lam[0], r"$\lambda_2$", color=COLORS["third"], va="center")
ax.set_xlim(0, 180)
ax.set_xticks([0, 45, 90, 135, 180])
ax.set_xlabel(r"direction $\theta$ of $\mathbf{u}=(\cos\theta,\sin\theta)^{\top}$ (degrees)")
ax.set_ylabel(r"$\mathbf{u}^{\top}\mathbf{A}\mathbf{u}$")
ax.set_title("(b) Rayleigh quotient")
fig.tight_layout()
fig.savefig(figure_path("ch01", "rayleigh"))
