"""Chapter 8, Section 7: contrasts among seven ordered groups.

1996 ANES (public domain, statsmodels.datasets.anes96): self-placement (selfLR) by party
identification (PID, seven ordered levels, unequal group sizes). Polynomial contrasts,
their sums of squares, what happens to orthogonality with unequal sizes, and the difference
between a coding matrix and the contrasts its coefficients estimate.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].to_numpy()
level = anes["PID"].to_numpy().astype(int)       # 0 = strong Dem. ... 6 = strong Rep.
g = 7
n_k = np.bincount(level).astype(float)
means = np.bincount(level, weights=y) / n_k
# <</data>>


# <<polynomial>>
def orthogonal_polynomials(scores, weights):
    """Columns 1..g-1: polynomials in the scores, orthogonal in sum_k w_k u_k v_k."""
    V = np.vander(scores, len(scores), increasing=True)      # 1, s, s^2, ...
    W = np.sqrt(weights)[:, None]
    Q, _ = np.linalg.qr(W * V)                               # weighted Gram-Schmidt
    P = Q / W                                                # back to the original scale
    P = P / np.abs(P).max(axis=0)                            # cosmetic rescaling
    return P[:, 1:]                                          # drop the constant column

def contrast_ss(c, means, n_k):
    """Estimate c'mu-hat and its sum of squares (c'mu-hat)^2 / sum(c_k^2 / n_k)."""
    psi = c @ means
    return psi, psi ** 2 / np.sum(c ** 2 / n_k)

scores = np.arange(g, dtype=float)
C_unw = orthogonal_polynomials(scores, np.ones(g))           # the classical table (equal n)
C_wtd = n_k[:, None] * orthogonal_polynomials(scores, n_k)   # c_k = n_k p(k)
ss_between = np.sum(n_k * (means - y.mean()) ** 2)
for name, C in [("unweighted", C_unw), ("weighted", C_wtd)]:
    ss = [contrast_ss(C[:, j], means, n_k)[1] for j in range(g - 1)]
    print(f"{name:10s} SS by degree: {np.round(ss, 2)}  total {sum(ss):.2f}",
          f"(between groups {ss_between:.2f})")
# <</polynomial>>

ss_unw = np.array([contrast_ss(C_unw[:, j], means, n_k)[1] for j in range(g - 1)])
ss_wtd = np.array([contrast_ss(C_wtd[:, j], means, n_k)[1] for j in range(g - 1)])
# contrasts: columns sum to zero (unweighted) or n-weighted sum to zero (weighted)
assert np.allclose(C_unw.sum(axis=0), 0) and np.allclose(C_wtd.sum(axis=0), 0)
# the unweighted linear contrast is proportional to -3, ..., 3
lin = C_unw[:, 0] / C_unw[-1, 0] * 3
assert np.allclose(lin, scores - 3)
quad = C_unw[:, 1] / C_unw[0, 1] * 5
assert np.allclose(quad, [5, 0, -3, -4, -3, 0, 5])
# orthogonality: weighted contrasts give an exact decomposition of the between-group SS
assert np.isclose(ss_wtd.sum(), ss_between)
assert not np.isclose(ss_unw.sum(), ss_between, rtol=1e-3)
G = C_wtd.T @ np.diag(1 / n_k) @ C_wtd                     # sum_k c_k d_k / n_k
assert np.allclose(G - np.diag(np.diag(G)), 0, atol=1e-9)
# the between-group SS equals the SS for the full set of g-1 contrasts
Zm = np.eye(g)[level]
M = Zm @ np.diag(1 / n_k) @ Zm.T
assert np.isclose(y @ M @ y - len(y) * y.mean() ** 2, ss_between)

# <<helmert>>
def helmert(g):
    """Coding matrix whose column k compares level k+1 with the mean of levels 1..k."""
    C = np.zeros((g, g - 1))
    for k in range(1, g):
        C[:k, k - 1] = -1.0
        C[k, k - 1] = k
    return C

H = helmert(g)
K = np.column_stack([np.ones(g), H])
gamma = np.linalg.solve(K, means)                  # Helmert-coded coefficients
psi = H.T @ means                                  # the contrasts c_k' mu-hat
print("coefficients:", np.round(gamma[1:], 4))
print("contrasts:   ", np.round(psi, 4))
print("ratio:       ", np.round(psi / gamma[1:], 1))   # = c_k' c_k = k(k+1)
# <</helmert>>

assert np.allclose(psi / gamma[1:], [k * (k + 1) for k in range(1, g)])
assert np.isclose(gamma[0], means.mean())
# a non-orthogonal coding: reference coding's columns are not the contrasts it estimates
R = np.vstack([np.zeros(g - 1), np.eye(g - 1)])
gamma_ref = np.linalg.solve(np.column_stack([np.ones(g), R]), means)
assert np.allclose(gamma_ref[1:], means[1:] - means[0])

gen = Generated("ch08", "contrasts", prefix="ct")
for j in range(3):
    psi_j, ss_j = contrast_ss(np.array([scores - 3, [5, 0, -3, -4, -3, 0, 5],
                                        [-1, 1, 1, 0, -1, -1, 1]][j], float), means, n_k)
    gen.num(f"psi{j + 1}", psi_j, 3)
    gen.num(f"ss{j + 1}", ss_j, 2)
    assert np.isclose(ss_j, ss_unw[j])
gen.num("ssbetween", ss_between, 2)
gen.num("ssunwsum", ss_unw.sum(), 2)
gen.num("sswtd1", ss_wtd[0], 2)
gen.num("sswtdrest", ss_wtd[1:].sum(), 2)
gen.num("share_lin", ss_unw[0] / ss_between * 100, 1)
for k in range(3):
    gen.num(f"helm:g{k + 1}", gamma[k + 1], 4)
    gen.num(f"helm:psi{k + 1}", psi[k], 4)
gen.num("v_lin", np.sum((scores - 3) ** 2 / n_k), 4)
gen.write()

# ---- figure: cell means with the linear and quadratic polynomial fits ----------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.6))
xs = np.linspace(0, 6, 200)
for deg, col, ls in [(1, COLORS["second"], "--"), (2, COLORS["third"], "-")]:
    coef = np.polyfit(scores, means, deg, w=np.sqrt(n_k))
    ax.plot(xs, np.polyval(coef, xs), color=col, ls=ls, lw=1.0,
            label="linear trend" if deg == 1 else "quadratic trend")
ax.scatter(scores, means, s=n_k / 3, color=COLORS["accent"], zorder=3, label="cell means")
ax.set_xticks(scores)
ax.set_xticklabels(["strong\nDem", "weak\nDem", "ind.\nDem", "ind.", "ind.\nRep", "weak\nRep", "strong\nRep"], fontsize=7)
ax.set_ylabel("mean self-placement")
ax.legend(frameon=False, fontsize=7, loc="upper left")
fig.savefig(figure_path("ch08", "pid_trend"))
