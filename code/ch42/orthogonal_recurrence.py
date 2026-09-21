"""Chapter 42, Section 2: the three-term recurrence for the orthogonal polynomials of a design.

Annual flow of the Nile at Aswan, 1871-1970 (10^8 m^3), the series of Example 9.4
(exm-ss-nile-polynomials). The recurrence reproduces the polynomials that chapter built
by a QR factorization, at a fraction of the work and without ever forming a power of the
design vector beyond the first.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<recurrence>>
import numpy as np
import statsmodels.api as sm

nile = sm.datasets.nile.load_pandas().data
year = nile["year"].to_numpy().astype(float)
y = nile["volume"].to_numpy()
n = len(y)
u = 2 * (year - year.min()) / (year.max() - year.min()) - 1     # design rescaled to [-1, 1]


def orthogonal_polynomials(t, d):
    """Monic polynomials p_0, ..., p_d orthogonal over the design t, by the recurrence."""
    P = np.zeros((len(t), d + 1))
    P[:, 0] = 1.0
    a, b, sq = np.zeros(d), np.zeros(d), [float(len(t))]
    for j in range(d):
        a[j] = (t * P[:, j] @ P[:, j]) / sq[j]
        P[:, j + 1] = (t - a[j]) * P[:, j]
        if j > 0:
            b[j] = sq[j] / sq[j - 1]
            P[:, j + 1] -= b[j] * P[:, j - 1]
        sq.append(P[:, j + 1] @ P[:, j + 1])
    return P, a, b, np.array(sq)


P, a, b, sq = orthogonal_polynomials(u, 8)
Q = P / np.sqrt(sq)                                             # orthonormal columns
print("a_j:", np.round(a[:5], 6))
print("b_j:", np.round(b[1:5], 5))
print("largest inner product between distinct columns:",
      f"{np.max(np.abs(Q.T @ Q - np.eye(9))):.2e}")
# <</recurrence>>

# <<sequential>>
coef = Q.T @ y                                    # gamma_j = <q_j, y>, unchanged by the degree
seq_ss = coef[1:] ** 2                            # sequential sum of squares for each degree
sse = np.sum((y - y.mean()) ** 2) - np.cumsum(seq_ss)
for j in range(1, 7):
    F = seq_ss[j - 1] / (sse[j - 1] / (n - j - 1))
    print(f"degree {j}: SS = {seq_ss[j - 1]:9.0f}   SSE = {sse[j - 1]:9.0f}   F = {F:7.2f}")
# <</sequential>>

# <<conditioning>>
for d in (2, 4, 6, 8):
    powers = {"raw years": np.vander(year, d + 1, increasing=True),
              "rescaled": np.vander(u, d + 1, increasing=True),
              "orthogonal": orthogonal_polynomials(u, d)[0]}
    line = "   ".join(f"{k} {np.linalg.cond(v):.2e}" for k, v in powers.items())
    print(f"degree {d}:  {line}   orthonormal {np.linalg.cond(Q[:, :d + 1]):.2f}")
# <</conditioning>>

# ---- assertions -------------------------------------------------------------
assert np.allclose(Q.T @ Q, np.eye(9), atol=1e-12)
# each p_j is monic of degree j: its coefficient on u^j is 1
for j in range(9):
    c = np.polyfit(u, P[:, j], j)
    assert np.isclose(c[0], 1.0, atol=1e-8)
# the design is symmetric about its centre, so every a_j vanishes
assert np.max(np.abs(a)) < 1e-12
# the recurrence reproduces the Gram-Schmidt basis of Example 9.4 (up to sign)
V = np.vander(u, 9, increasing=True)
Qqr, R = np.linalg.qr(V)
Qqr = Qqr * np.sign(np.diag(R))
qr_gap = np.max(np.abs(np.abs(Q) - np.abs(Qqr)))
assert qr_gap < 1e-12
# and the sequential sums of squares of that example
assert np.allclose(seq_ss[:4], [613893, 309415, 1894, 134142], atol=1.0)
# sequential SS is the drop in residual sum of squares, whatever basis is used
for j in range(1, 7):
    fit = lambda Z: np.sum((y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]) ** 2)
    assert np.isclose(seq_ss[j - 1], fit(V[:, :j]) - fit(V[:, :j + 1]), rtol=1e-6)

# the discrete recurrence coefficients approach those of the Legendre polynomials
legendre_b = np.array([j ** 2 / (4 * j ** 2 - 1) for j in range(1, 5)])
assert np.all(np.abs(b[1:5] - legendre_b) < 0.01)

conds = {}
for d in (2, 4, 6, 8):
    conds[d] = (np.linalg.cond(np.vander(year, d + 1, increasing=True)),
                np.linalg.cond(np.vander(u, d + 1, increasing=True)),
                np.linalg.cond(orthogonal_polynomials(u, d)[0]),
                np.linalg.cond(Q[:, :d + 1]))
    assert conds[d][3] < 1.0 + 1e-8 < conds[d][2] < conds[d][1] < conds[d][0]

F_stats = {}
for j in range(1, 7):
    F_stats[j] = seq_ss[j - 1] / (sse[j - 1] / (n - j - 1))
p6 = 1 - stats.f.cdf(F_stats[6], 1, n - 7)

gen = Generated("ch42", "orthogonal_recurrence", prefix="ort")
gen.int("n", n)
for j in range(1, 5):
    gen.num(f"b{j}", b[j] if j < 5 else 0.0, 5)
    gen.num(f"leg{j}", legendre_b[j - 1], 5)
for j in range(1, 7):
    gen.num(f"ss{j}", seq_ss[j - 1], 0)
    gen.num(f"F{j}", F_stats[j], 2)
gen.num("p6", p6, 3)
gen.num("qrgap", qr_gap, 1, sci=True)
gen.num("orth", np.max(np.abs(Q.T @ Q - np.eye(9))), 1, sci=True)
for d in (2, 4, 6, 8):
    gen.num(f"kraw{d}", conds[d][0], 2, sci=True)
    gen.num(f"ksca{d}", conds[d][1], 2, sci=True)
    gen.num(f"kmon{d}", conds[d][2], 2, sci=True)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
for j, colour in ((1, COLORS["third"]), (2, COLORS["accent"]),
                  (3, COLORS["second"]), (4, COLORS["thread"])):
    ax.plot(year, Q[:, j], color=colour, label=f"$q_{j}$")
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel("year")
ax.set_title("(a) orthonormal polynomials of the design")
ax.legend(frameon=False, loc="upper center", ncol=4, fontsize=6.5)
ax = axes[1]
degrees = np.arange(1, 13)
curves = {"raw years": [], "rescaled to $[-1,1]$": [], "monic orthogonal": []}
for d in degrees:
    curves["raw years"].append(np.linalg.cond(np.vander(year, d + 1, increasing=True)))
    curves["rescaled to $[-1,1]$"].append(np.linalg.cond(np.vander(u, d + 1, increasing=True)))
    curves["monic orthogonal"].append(np.linalg.cond(orthogonal_polynomials(u, d)[0]))
for (name, vals), colour in zip(curves.items(),
                                [COLORS["second"], COLORS["accent"], COLORS["third"]]):
    ax.semilogy(degrees, vals, color=colour, marker="o", markersize=2.5, label=name)
ax.semilogy(degrees, [np.linalg.cond(orthogonal_polynomials(u, d)[0]
                                     / np.sqrt(orthogonal_polynomials(u, d)[3])) for d in degrees],
            color=COLORS["thread"], marker="o", markersize=2.5, label="orthonormal")
ax.set_ylim(0.5, 1e60)
ax.set_yticks([1e0, 1e15, 1e30, 1e45, 1e60])
ax.set_xlabel("degree")
ax.set_ylabel("condition number")
ax.set_title("(b) conditioning of four bases")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
fig.tight_layout()
fig.savefig(figure_path("ch42", "orthogonal_recurrence"))
