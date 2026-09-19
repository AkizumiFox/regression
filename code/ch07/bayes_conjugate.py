"""Chapter 7, Section 6: the conjugate normal-inverse-gamma analysis of the stack loss data.

Prior: beta | sigma^2 ~ N(m0, sigma^2 V0), sigma^2 ~ inverse gamma(a0, b0).
Posterior: the same family, with V_n = (V0^{-1} + X^T X)^{-1}, m_n = V_n (V0^{-1} m0 + X^T y),
a_n = a0 + n/2, b_n = b0 + (y^T y + m0^T V0^{-1} m0 - m_n^T V_n^{-1} m_n)/2.
The script checks the formulas against the unnormalized log posterior, checks the
shrinkage form of m_n, and checks the flat-prior limit against least squares and t.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape
XtX, Xty = X.T @ X, X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)
sse = np.sum((y - X @ beta_hat) ** 2)
# <</data>>


# <<posterior>>
def nig_posterior(m0, V0, a0, b0):
    """Conjugate update of the normal-inverse-gamma prior NIG(m0, V0, a0, b0)."""
    P0 = np.linalg.inv(V0)                         # prior precision (per unit sigma^2)
    Vn = np.linalg.inv(P0 + XtX)
    mn = Vn @ (P0 @ m0 + Xty)
    an = a0 + n / 2
    bn = b0 + (y @ y + m0 @ P0 @ m0 - mn @ (P0 + XtX) @ mn) / 2
    return mn, Vn, an, bn


# prior: vague intercept; slopes centred on 1, 1, 0 with sd 0.1, 0.2, 0.1 in units of sigma
m0 = np.array([0.0, 1.0, 1.0, 0.0])
V0 = np.diag([100.0**2, 0.1**2, 0.2**2, 0.1**2])
a0, b0 = 3.0, 18.0                                 # prior mean of sigma^2 is b0/(a0-1) = 9
mn, Vn, an, bn = nig_posterior(m0, V0, a0, b0)
scale = np.sqrt(np.diag(Vn) * bn / an)             # marginal posterior: t with 2 a_n d.f.
q = stats.t.ppf(0.975, 2 * an)
for j, name in enumerate(["intercept", "air flow", "water temp", "acid conc"]):
    print(f"{name:10s} LS {beta_hat[j]:8.4f}   posterior mean {mn[j]:8.4f}   "
          f"95% interval ({mn[j] - q * scale[j]:.3f}, {mn[j] + q * scale[j]:.3f})")
print(f"posterior mean of sigma^2 = {bn / (an - 1):.3f}")
# <</posterior>>


def log_nig_density(beta, s2, m, V, a, b):
    """Normalized log density of NIG(m, V, a, b) at (beta, sigma^2)."""
    k = len(m)
    d = beta - m
    _, logdetV = np.linalg.slogdet(V)
    log_normal = -k / 2 * np.log(2 * np.pi * s2) - logdetV / 2 - d @ np.linalg.solve(V, d) / (2 * s2)
    return log_normal + stats.invgamma.logpdf(s2, a, scale=b)


def log_prior_times_lik(beta, s2):
    r = y - X @ beta
    return (log_nig_density(beta, s2, m0, V0, a0, b0)
            - n / 2 * np.log(2 * np.pi * s2) - r @ r / (2 * s2))


# the posterior differs from prior x likelihood by a constant (the log marginal likelihood)
rng = np.random.default_rng(76)
diffs = []
for _ in range(50):
    b = beta_hat + rng.normal(size=p) * np.array([5, 0.3, 0.5, 0.3])
    s2 = rng.uniform(2, 30)
    diffs.append(log_prior_times_lik(b, s2) - log_nig_density(b, s2, mn, Vn, an, bn))
diffs = np.array(diffs)
assert np.ptp(diffs) < 1e-8

# the equivalent form of b_n, and the shrinkage (matrix-weighted average) form of m_n
d = beta_hat - m0
bn_alt = b0 + (sse + d @ np.linalg.solve(V0 + np.linalg.inv(XtX), d)) / 2
assert np.isclose(bn, bn_alt)
W = Vn @ np.linalg.inv(V0)
assert np.allclose(mn, W @ m0 + (np.eye(p) - W) @ beta_hat)
assert np.allclose(np.eye(p) - W, Vn @ XtX)
# the two pieces of b_n: lack of fit at m_n, and the prior-data conflict term
rss_mn = np.sum((y - X @ mn) ** 2)
conflict = (mn - m0) @ np.linalg.solve(V0, mn - m0)
assert np.isclose(bn, b0 + (rss_mn + conflict) / 2)
assert rss_mn > sse and conflict > 0
sigma2_noconf = (b0 + rss_mn / 2) / (an - 1)      # what the posterior mean would be without the conflict term
assert sigma2_noconf < b0 / (a0 - 1) < bn / (an - 1)

# Monte Carlo check of the marginal t: draw sigma^2 then beta
S2 = stats.invgamma.rvs(an, scale=bn, size=200_000, random_state=rng)
B1 = mn[1] + np.sqrt(S2 * Vn[1, 1]) * rng.normal(size=S2.size)
assert abs(B1.mean() - mn[1]) < 0.003
assert abs(np.quantile(B1, 0.975) - (mn[1] + q * scale[1])) < 0.01

# flat-prior limit: V0^{-1} -> 0, a0 = -p/2, b0 = 0 reproduces least squares and the t interval
mf, Vf, af, bf = nig_posterior(np.zeros(p), 1e10 * np.eye(p), -p / 2, 0.0)
s2 = sse / (n - p)
assert np.allclose(mf, beta_hat, atol=1e-5)
assert np.isclose(2 * af, n - p)
assert np.isclose(bf / af, s2, rtol=1e-6)
se = np.sqrt(np.diag(Vf) * bf / af)
assert np.allclose(se, np.sqrt(s2 * np.diag(np.linalg.inv(XtX))), rtol=1e-5)

# shrinkage path: multiply the prior precision of the slopes by tau
taus = np.logspace(-3, 3, 61)
path = []
for tau in taus:
    V0t = V0.copy()
    V0t[1:, 1:] /= tau
    path.append(nig_posterior(m0, V0t, a0, b0)[0])
path = np.array(path)
assert np.allclose(path[0, 1:], beta_hat[1:], atol=5e-3)
assert np.allclose(path[-1, 1:], m0[1:], atol=5e-3)
assert path[:, 2].min() < m0[2] - 0.03            # water temperature overshoots its prior mean

gen = Generated("ch07", "bayes_conjugate", prefix="bay")
for j in range(p):
    gen.num(f"ls{j}", beta_hat[j], 4)
    gen.num(f"mn{j}", mn[j], 4)
    gen.num(f"lo{j}", mn[j] - q * scale[j], 3)
    gen.num(f"hi{j}", mn[j] + q * scale[j], 3)
    gen.num(f"w{j}", W[j, j], 3)
gen.num("wt_min", path[:, 2].min(), 3)
gen.num("an", an, 1)
gen.int("df", round(2 * an))
gen.num("bn", bn, 3)
gen.num("sigma2_post", bn / (an - 1), 3)
gen.num("q", q, 3)
gen.num("rss_mn", rss_mn, 3)
gen.num("conflict", conflict, 3)
gen.num("sigma2_noconf", sigma2_noconf, 3)
gen.num("sse", sse, 3)
tq = stats.t.ppf(0.975, n - p)
se1 = np.sqrt(s2 * np.linalg.inv(XtX)[1, 1])
gen.num("tlo1", beta_hat[1] - tq * se1, 3)
gen.num("thi1", beta_hat[1] + tq * se1, 3)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.9, 2.5))
ax = axes[0]
g = np.linspace(0.2, 1.5, 400)
prior_scale = np.sqrt(V0[1, 1] * b0 / a0)
ax.plot(g, stats.t.pdf(g, 2 * a0, loc=m0[1], scale=prior_scale), color=COLORS["muted"],
        ls="--", label="prior")
ax.plot(g, stats.t.pdf(g, n - p, loc=beta_hat[1], scale=se1), color=COLORS["second"],
        label="flat-prior posterior (LS)")
ax.plot(g, stats.t.pdf(g, 2 * an, loc=mn[1], scale=scale[1]), color=COLORS["accent"],
        label="posterior")
ax.set_xlabel(r"air-flow coefficient $\beta_1$")
ax.set_ylabel("density")
ax.set_ylim(0, 6.2)
ax.set_title("(a) prior and two posteriors")
ax.legend(frameon=False, fontsize=6.5, loc="upper right")
ax = axes[1]
names = ["air flow", "water temp", "acid conc"]
colors = [COLORS["accent"], COLORS["second"], COLORS["third"]]
for j in range(1, p):
    ax.semilogx(taus, path[:, j], color=colors[j - 1], label=names[j - 1])
    ax.plot(taus[0], beta_hat[j], "o", ms=3, color=colors[j - 1])
    ax.plot(taus[-1], m0[j], "s", ms=3, color=colors[j - 1])
ax.axvline(1, color=COLORS["grid"], lw=0.7, zorder=0)
ax.set_xlabel(r"prior precision multiplier $\tau$")
ax.set_ylabel("posterior mean")
ax.set_title("(b) from least squares to the prior")
ax.legend(frameon=False, fontsize=6.5, loc="center right")
fig.tight_layout()
fig.savefig(figure_path("ch07", "bayes_conjugate"))
