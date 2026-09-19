"""Chapter 4, Section 3: when is a quadratic form chi-squared?

Y ~ N_n(mu, Sigma) with AR(1) covariance and a trending mean. The GLS-centred form
Y^T A Y, A = W - W 1 (1^T W 1)^{-1} 1^T W with W = Sigma^{-1}, has A Sigma idempotent,
so it is chi^2(n-1, mu^T A mu). The ordinary centred form Y^T C Y is not chi-squared:
its law is a weighted sum of noncentral chi^2(1)'s, computed exactly by numerical
inversion of the characteristic function (Gil-Pelaez) and compared with simulation
and with a two-moment scaled chi-squared approximation.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate, stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(31415)

n, rho = 8, 0.6
idx = np.arange(n)
Sigma = rho ** np.abs(idx[:, None] - idx[None, :])     # AR(1) correlation matrix
mu = 0.4 * (idx - idx.mean())                          # linear trend, mean zero
one = np.ones(n)

# <<gls>>
W = np.linalg.inv(Sigma)
A = W - np.outer(W @ one, W @ one) / (one @ W @ one)   # GLS-centred form
AS = A @ Sigma
print("A Sigma idempotent:", np.allclose(AS @ AS, AS))
df, gamma = np.trace(AS), mu @ A @ mu                  # chi^2(df, gamma)

L = np.linalg.cholesky(Sigma)                          # Y = mu + L Z
Z = rng.standard_normal((200_000, n))
Y = mu + Z @ L.T
q_gls = np.einsum("ij,jk,ik->i", Y, A, Y)
print(df, gamma, stats.kstest(q_gls, stats.ncx2(df, gamma).cdf).pvalue)
# <</gls>>

assert np.allclose(AS @ AS, AS) and np.isclose(df, n - 1)
assert np.linalg.matrix_rank(A) == n - 1
assert stats.kstest(q_gls, stats.ncx2(round(df), gamma).cdf).pvalue > 1e-3


# <<canonical>>
def canonical(A, mu, L):
    """Y^T A Y = sum_j lam_j (W_j + nu_j)^2 + c (all lam_j nonzero here), W ~ N(0, I)."""
    lam, Qm = np.linalg.eigh(L.T @ A @ L)
    b = Qm.T @ L.T @ A @ mu
    keep = np.abs(lam) > 1e-10
    assert np.allclose(b[~keep], 0)                    # no linear terms in this example
    c = mu @ A @ mu - np.sum(b[keep] ** 2 / lam[keep])
    return lam[keep], b[keep] / lam[keep], c


def tail_prob(x, lam, nu, c=0.0):
    """P(sum lam_j (W_j + nu_j)^2 + c > x) by Gil-Pelaez inversion."""
    def integrand(u):
        z = 1 - 2j * lam * u
        phi = np.prod(z ** -0.5 * np.exp(1j * u * lam * nu**2 / z))
        return (np.exp(-1j * u * (x - c)) * phi).imag / u
    edges = np.linspace(0, 400, 401)                   # integrand decays like u^(-1-k/2)
    val = sum(integrate.quad(integrand, a, b)[0] for a, b in zip(edges[:-1], edges[1:]))
    return 0.5 + val / np.pi
# <</canonical>>

C = np.eye(n) - np.outer(one, one) / n                 # ordinary centring
lam, nu, c = canonical(C, mu, L)
assert np.isclose(c, 0, atol=1e-10)
assert np.allclose(np.sort(lam), np.sort(np.linalg.eigvals(C @ Sigma).real)[1:])
q_ols = np.einsum("ij,jk,ik->i", Y, C, Y)

# sanity check of the inversion on a case with a known answer
assert np.isclose(tail_prob(9.0, np.ones(4), np.array([1.0, 1.0, 0.5, 0.0])),
                  stats.ncx2.sf(9.0, 4, 2.25), atol=1e-6)

# exact tail probabilities vs simulation
cuts = np.quantile(q_ols, [0.5, 0.9, 0.99])
for x in cuts:
    p_exact = tail_prob(x, lam, nu)
    p_sim = np.mean(q_ols > x)
    assert abs(p_exact - p_sim) < 5 * np.sqrt(p_sim * (1 - p_sim) / len(q_ols)) + 1e-4

# two-moment (Satterthwaite) approximation a * chi^2(nu_s): match mean and variance
m = np.trace(C @ Sigma) + mu @ C @ mu
v = 2 * np.trace(C @ Sigma @ C @ Sigma) + 4 * mu @ C @ Sigma @ C @ mu
a_s, nu_s = v / (2 * m), 2 * m**2 / v
x95 = np.quantile(q_ols, 0.95)
p_exact95 = tail_prob(x95, lam, nu)
p_satt95 = stats.chi2.sf(x95 / a_s, nu_s)
p_naive95 = stats.ncx2.sf(x95, n - 1, mu @ C @ mu)     # pretending Sigma = I
assert abs(p_exact95 - 0.05) < 0.003

# the naive error is not monotone in rho, but strong correlation makes it severe
rho_high = 0.9
L_high = np.linalg.cholesky(rho_high ** np.abs(idx[:, None] - idx[None, :]))
q_high = np.einsum("ij,jk,ik->i", mu + Z @ L_high.T, C, mu + Z @ L_high.T)
p_naive_high = stats.ncx2.sf(np.quantile(q_high, 0.95), n - 1, mu @ C @ mu)
assert p_naive_high > 0.15 and abs(p_naive_high - 0.05) > 3 * abs(p_naive95 - 0.05)

gen = Generated("ch04", "quadform_chisq", prefix="qfc")
gen.int("n", n)
gen.num("rho", rho, 1)
gen.int("df", round(df))
gen.num("gamma", gamma, 3)
gen.num("gammaols", mu @ C @ mu, 3)
gen.num("lammax", lam.max(), 3)
gen.num("lammin", lam.min(), 3)
gen.num("x95", x95, 2)
gen.num("pexact", p_exact95, 4)
gen.num("psatt", p_satt95, 4)
gen.num("pnaive", p_naive95, 3)
gen.num("rhohigh", rho_high, 1)
gen.num("pnaivehigh", p_naive_high, 3)
gen.num("satta", a_s, 3)
gen.num("sattnu", nu_s, 2)
gen.text("reps", "200{,}000")
gen.write()

# ---- figure ---------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
bins = np.linspace(0, np.quantile(q_gls, 0.995), 60)
ax.hist(q_gls, bins=bins, density=True, color=COLORS["accent"], alpha=0.35, linewidth=0)
xs = np.linspace(0.01, bins[-1], 400)
ax.plot(xs, stats.ncx2.pdf(xs, n - 1, gamma), color=COLORS["accent"])
ax.set_title(rf"(a) GLS-centred: $\chi^2({n-1},\,{gamma:.2f})$")
ax.set_xlabel("quadratic form")
ax.set_ylabel("density")
ax = axes[1]
bins = np.linspace(0, np.quantile(q_ols, 0.995), 60)
ax.hist(q_ols, bins=bins, density=True, color=COLORS["second"], alpha=0.35, linewidth=0)
xs = np.linspace(0.01, bins[-1], 400)
ax.plot(xs, stats.chi2.pdf(xs / a_s, nu_s) / a_s, color=COLORS["second"],
        label="two-moment fit")
ax.plot(xs, stats.ncx2.pdf(xs, n - 1, mu @ C @ mu), color=COLORS["muted"],
        linestyle="--", label=r"$\chi^2$ if $\Sigma=I$")
ax.set_title("(b) ordinary centring")
ax.set_xlabel("quadratic form")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch04", "quadform_chisq"))
