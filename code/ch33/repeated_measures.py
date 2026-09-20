"""Chapter 33, Section 3: repeated measures, sphericity and the corrections.

Twenty-four patients have grip strength measured at six weekly visits. The data are
simulated with first-order autoregressive within-patient errors, so the sphericity
condition fails and the uncorrected univariate F test is not exact.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

N, m = 24, 6
mu_true = np.array([48.0, 48.5, 48.9, 49.2, 49.4, 49.5])   # kg, a slow improvement
sigma, rho_true = 2.0, 0.6

lags = np.abs(np.subtract.outer(np.arange(m), np.arange(m)))
Sigma_true = sigma**2 * rho_true**lags               # AR(1) within a patient

rng = np.random.default_rng(3303)
Yobs = np.round(mu_true + rng.multivariate_normal(np.zeros(m), Sigma_true, N), 2)


# <<contrast>>
def contrast_basis(m):
    """Orthonormal basis U of the space of contrasts among m repeated measures."""
    Q, _ = np.linalg.qr(np.column_stack([np.ones(m), np.eye(m)[:, 1:]]))
    return Q[:, 1:]                                  # columns orthonormal and orthogonal to 1


def epsilon(Psi):
    """Box's measure of departure from sphericity, 1/(m-1) <= eps <= 1."""
    k = Psi.shape[0]
    return np.trace(Psi)**2 / (k * np.sum(Psi * Psi))


U = contrast_basis(m)
Z = Yobs @ U                                         # N x (m-1) contrasts, one row per patient
zbar = Z.mean(axis=0)
Sz = np.cov(Z, rowvar=False, ddof=1)                 # Wishart, (N-1) degrees of freedom

SS_time = N * zbar @ zbar
SS_err = (N - 1) * np.trace(Sz)
F_uncorrected = N * (zbar @ zbar) / np.trace(Sz)
eps_hat = epsilon(Sz)
eps_hf = min(1.0, (N * (m - 1) * eps_hat - 2)
             / ((m - 1) * ((N - 1) - (m - 1) * eps_hat)))
d1, d2 = m - 1, (N - 1) * (m - 1)
p_un = stats.f.sf(F_uncorrected, d1, d2)
p_gg = stats.f.sf(F_uncorrected, d1 * eps_hat, d2 * eps_hat)
p_hf = stats.f.sf(F_uncorrected, d1 * eps_hf, d2 * eps_hf)
T2 = N * zbar @ np.linalg.solve(Sz, zbar)            # Hotelling's statistic
F_mv = (N - m + 1) / ((N - 1) * (m - 1)) * T2
p_mv = stats.f.sf(F_mv, m - 1, N - m + 1)
print(f"F = {F_uncorrected:.3f} on ({d1}, {d2}) df,  p = {p_un:.4f}")
print(f"Greenhouse-Geisser eps = {eps_hat:.3f}, p = {p_gg:.4f}")
print(f"Huynh-Feldt eps = {eps_hf:.3f}, p = {p_hf:.4f}")
print(f"multivariate F({m - 1}, {N - m + 1}) = {F_mv:.3f}, p = {p_mv:.4f}")
# <</contrast>>

assert np.allclose(U.T @ U, np.eye(m - 1)) and np.allclose(U.T @ np.ones(m), 0)
# the two sums of squares are the usual repeated-measures ones
grand = Yobs.mean()
assert np.isclose(SS_time, N * np.sum((Yobs.mean(axis=0) - grand)**2))
resid = Yobs - Yobs.mean(axis=0) - Yobs.mean(axis=1)[:, None] + grand
assert np.isclose(SS_err, np.sum(resid**2))
# the basis does not matter
rng_b = np.random.default_rng(1)
Q2, _ = np.linalg.qr(np.column_stack([np.ones(m), rng_b.standard_normal((m, m - 1))]))
Z2 = Yobs @ Q2[:, 1:]
assert np.isclose(epsilon(np.cov(Z2, rowvar=False, ddof=1)), eps_hat)
assert np.isclose(N * Z2.mean(axis=0) @ Z2.mean(axis=0), SS_time)
# epsilon lies between its lower bound and 1
assert 1 / (m - 1) <= eps_hat <= 1

eps_true = epsilon(U.T @ Sigma_true @ U)
# a covariance of Huynh-Feldt form has epsilon exactly 1, without being equicorrelated
a_vec = np.array([0.4, 0.9, 1.5, 2.2, 1.1, 0.6])
Sigma_hf = np.add.outer(a_vec, a_vec) / 2 + 1.0 * np.eye(m)
assert np.isclose(epsilon(U.T @ Sigma_hf @ U), 1.0)
assert np.min(np.linalg.eigvalsh(Sigma_hf)) > 0
assert not np.allclose(Sigma_hf - np.diag(np.diag(Sigma_hf)),
                       Sigma_hf[0, 1] * (np.ones((m, m)) - np.eye(m)))


# ---- level of the four tests -------------------------------------------------
def levels(Sigma, B=20000, seed=5, N=N):
    """Simulated rejection rates under the hypothesis of no change over time."""
    m = Sigma.shape[0]
    U = contrast_basis(m)
    Psi = U.T @ Sigma @ U
    L = np.linalg.cholesky(Psi)
    rng = np.random.default_rng(seed)
    d1, d2 = m - 1, (N - 1) * (m - 1)
    out = np.zeros(4)
    for _ in range(B):
        Z = rng.standard_normal((N, m - 1)) @ L.T
        zb = Z.mean(axis=0)
        S = np.cov(Z, rowvar=False, ddof=1)
        F = N * (zb @ zb) / np.trace(S)
        e = epsilon(S)
        ehf = min(1.0, (N * d1 * e - 2) / (d1 * ((N - 1) - d1 * e)))
        T2 = N * zb @ np.linalg.solve(S, zb)
        out += [F > stats.f.ppf(0.95, d1, d2),
                F > stats.f.ppf(0.95, d1 * e, d2 * e),
                F > stats.f.ppf(0.95, d1 * ehf, d2 * ehf),
                (N - m + 1) / d2 * T2 > stats.f.ppf(0.95, d1, N - m + 1)]
    return out / B


rhos = [0.0, 0.3, 0.6, 0.8]
lev = np.array([levels(sigma**2 * r**lags) for r in rhos])
eps_rho = [epsilon(U.T @ (sigma**2 * r**lags) @ U) for r in rhos]
for r, e, row in zip(rhos, eps_rho, lev):
    print(f"rho = {r:.1f}  eps = {e:.3f}  uncorrected {row[0]:.3f}  GG {row[1]:.3f}  "
          f"HF {row[2]:.3f}  multivariate {row[3]:.3f}")
lev_hf = levels(Sigma_hf)
print("Huynh-Feldt covariance:", lev_hf.round(3))

assert abs(lev[0, 0] - 0.05) < 0.005                 # exact when the errors are independent
assert lev[-1, 0] > 0.08                             # badly inflated under strong correlation
assert abs(lev_hf[0] - 0.05) < 0.006                 # exact under sphericity without symmetry
assert all(abs(lev[:, 3] - 0.05) < 0.006)            # the multivariate test is always exact
assert lev[-1, 1] < lev[-1, 0] and lev[-1, 2] < lev[-1, 0]

gen = Generated("ch33", "repeated_measures")
gen.int("N", N)
gen.int("m", m)
gen.num("rho_true", rho_true, 1)
gen.num("eps_true", eps_true, 3)
gen.num("eps_hat", eps_hat, 3)
gen.num("eps_hf", eps_hf, 3)
gen.num("eps_lower", 1 / (m - 1), 3)
gen.num("F", F_uncorrected, 3)
gen.num("p_un", p_un, 4)
gen.num("p_gg", p_gg, 4)
gen.num("p_hf", p_hf, 4)
gen.num("F_mv", F_mv, 3)
gen.num("T2", T2, 3)
gen.num("p_mv", p_mv, 4)
gen.num("df_gg1", d1 * eps_hat, 2)
gen.num("df_gg2", d2 * eps_hat, 1)
gen.num("df_hf1", d1 * eps_hf, 2)
gen.num("df_hf2", d2 * eps_hf, 1)
gen.int("d1", d1)
gen.int("d2", d2)
for j in range(m):
    gen.num(f"visit{j + 1}", Yobs.mean(axis=0)[j], 3)
for r, e, row in zip(rhos, eps_rho, lev):
    tag = str(int(10 * r))
    gen.num("eps_rho" + tag, e, 3)
    gen.num("lev_un" + tag, row[0], 3)
    gen.num("lev_gg" + tag, row[1], 3)
    gen.num("lev_hf" + tag, row[2], 3)
    gen.num("lev_mv" + tag, row[3], 3)
gen.num("lev_hfcov", lev_hf[0], 3)
gen.write()

# ---- figure ------------------------------------------------------------------
rng_f = np.random.default_rng(99)
Psi = U.T @ Sigma_true @ U
Lf = np.linalg.cholesky(Psi)
Fs = np.empty(40000)
for i in range(Fs.size):
    Zf = rng_f.standard_normal((N, m - 1)) @ Lf.T
    zb = Zf.mean(axis=0)
    Fs[i] = N * (zb @ zb) / np.trace(np.cov(Zf, rowvar=False, ddof=1))

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
grid = np.linspace(0.01, 6, 400)
ax.hist(Fs, bins=80, range=(0, 6), density=True, color=COLORS["grid"],
        edgecolor="none", label="simulated")
ax.plot(grid, stats.f.pdf(grid, d1, d2), color=COLORS["second"],
        label=f"F({d1}, {d2})")
ax.plot(grid, stats.f.pdf(grid, d1 * eps_true, d2 * eps_true), "--",
        color=COLORS["accent"], label=r"F($\epsilon d_1$, $\epsilon d_2$)")
ax.axvline(stats.f.ppf(0.95, d1, d2), color=COLORS["ink"], linewidth=0.6)
ax.set_xlabel("F statistic")
ax.set_ylabel("density")
ax.set_title(f"(a) null distribution, $\\epsilon$ = {eps_true:.2f}")
ax.legend(fontsize=6.5, frameon=False)
ax = axes[1]
labels = ["uncorrected", "Greenhouse–Geisser", "Huynh–Feldt", "multivariate"]
marks = ["o-", "s-", "d-", "^-"]
cols = [COLORS["second"], COLORS["accent"], COLORS["thread"], COLORS["third"]]
for j in range(4):
    ax.plot(rhos, lev[:, j], marks[j], color=cols[j], markersize=3, label=labels[j])
ax.axhline(0.05, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"autoregressive parameter $\rho$")
ax.set_ylabel("rejection rate")
ax.set_title("(b) level of a 5% test")
ax.legend(fontsize=6.5, frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch33", "sphericity"))
