"""Chapter 40, Section 3: approximations to the integrated likelihood of a GLMM.

Adaptive Gauss-Hermite quadrature, its non-adaptive form, the Laplace approximation (one
adaptive node) and penalized quasi-likelihood are compared on simulated clustered binary data:
m clusters of n observations, a random intercept, one covariate varying within the cluster.
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize, minimize_scalar

from regbook import COLORS, Generated, figure_path, use_book_style

M, N = 150, 5                       # clusters, observations per cluster
BETA = np.array([-0.5, 1.0])
TAU = 1.0
REPS = 200


# <<quadrature>>
def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


def cluster_mode(y, X, beta, tau):
    """Mode and curvature of log f(y_i | u) + log p(u), cluster by cluster."""
    eta0, u = X @ beta, np.zeros(y.shape[0])
    for _ in range(100):
        mu = expit(eta0 + u[:, None])
        step = ((y - mu).sum(1) - u / tau**2) / (-(mu * (1 - mu)).sum(1) - 1 / tau**2)
        u = u - step
        if np.max(np.abs(step)) < 1e-12:
            break
    mu = expit(eta0 + u[:, None])
    return u, 1.0 / np.sqrt((mu * (1 - mu)).sum(1) + 1 / tau**2)


def loglik(y, X, beta, tau, K, adapt=True):
    """Integrated log-likelihood by K-node Gauss-Hermite quadrature.

    Adaptive: centre and scale the nodes at each cluster's own mode and curvature.
    K = 1 with adapt=True is exactly the Laplace approximation."""
    node, w = hermgauss(K)
    if adapt:
        centre, scale = cluster_mode(y, X, beta, tau)
    else:
        centre, scale = np.zeros(y.shape[0]), np.full(y.shape[0], tau)
    pts = centre[:, None] + np.sqrt(2) * scale[:, None] * node            # (m, K)
    eta = (X @ beta)[:, :, None] + pts[:, None, :]                        # (m, n, K)
    cond = (y[:, :, None] * eta - np.logaddexp(0.0, eta)).sum(1)          # (m, K)
    prior = -0.5 * (pts / tau) ** 2 - np.log(tau) - 0.5 * np.log(2 * np.pi)
    terms = np.log(w) + node**2 + cond + prior
    return float(np.sum(np.log(np.sqrt(2) * scale) + np.logaddexp.reduce(terms, axis=1)))


def fit_ml(y, X, K, start=None):
    p = X.shape[2]
    obj = lambda par: -loglik(y, X, par[:p], np.exp(par[p]), K)
    x0 = np.r_[np.zeros(p), 0.0] if start is None else start
    r = minimize(obj, x0, method="Nelder-Mead",
                 options=dict(maxiter=8000, maxfev=8000, xatol=1e-8, fatol=1e-10))
    return r.x[:p], float(np.exp(r.x[p])), -r.fun
# <</quadrature>>


# <<pql>>
def weighted_lmm(z, X, W):
    """REML fit of a weighted random-intercept linear mixed model: z_ij = x_ij'b + u_i + e_ij
    with Var(e_ij) = 1/W_ij and Var(u_i) = t2. Returns b, t2 and the predicted u."""
    def solve(t2):
        s, WX = W.sum(1), W[..., None] * X
        den = 1 + t2 * s
        XtVX = (np.einsum("mnp,mnq->pq", WX, X)
                - t2 * np.einsum("mp,mq->pq", WX.sum(1), WX.sum(1) / den[:, None]))
        XtVz = (np.einsum("mnp,mn->p", WX, z)
                - t2 * np.einsum("mp,m->p", WX.sum(1), (W * z).sum(1) / den))
        return np.linalg.solve(XtVX, XtVz), XtVX, den

    def neg_reml(log_t2):
        t2 = np.exp(log_t2)
        b, XtVX, den = solve(t2)
        r = z - X @ b
        q = np.sum(W * r * r) - t2 * np.sum((W * r).sum(1) ** 2 / den)
        ld = np.sum(np.log(den)) - np.sum(np.log(W))
        return 0.5 * (ld + np.linalg.slogdet(XtVX)[1] + q)

    opt = minimize_scalar(neg_reml, bounds=(-12.0, 6.0), method="bounded",
                          options=dict(xatol=1e-10))
    t2 = float(np.exp(opt.x))
    b, _, den = solve(t2)
    r = z - X @ b
    return b, t2, t2 * (W * r).sum(1) / den


def fit_pql(y, X, maxit=300):
    """Breslow-Clayton penalized quasi-likelihood: a linear mixed model on the working
    response, refitted until it stops changing."""
    beta, u = np.zeros(X.shape[2]), np.zeros(y.shape[0])
    for _ in range(maxit):
        eta = X @ beta + u[:, None]
        mu = expit(eta)
        v = np.clip(mu * (1 - mu), 1e-8, None)          # working weights
        z = eta + (y - mu) / v                          # working response
        b_new, t2, u_new = weighted_lmm(z, X, v)
        done = max(np.max(np.abs(b_new - beta)), np.max(np.abs(u_new - u))) < 1e-10
        beta, u = b_new, u_new
        if done:
            break
    return beta, float(np.sqrt(t2)), u
# <</pql>>


# <<demo>>
M_demo, N_demo = 60, 5
rng_demo = np.random.default_rng(40_007)
x_demo = rng_demo.normal(size=(M_demo, N_demo))
X_demo = np.stack([np.ones((M_demo, N_demo)), x_demo], axis=2)
u_demo = rng_demo.normal(0, TAU, M_demo)
y_demo = (rng_demo.uniform(size=(M_demo, N_demo))
          < expit(X_demo @ BETA + u_demo[:, None])).astype(float)

for K in [1, 5, 15]:
    b, t_hat, ll = fit_ml(y_demo, X_demo, K)
    print(f"AGHQ K = {K:2d}:  beta = {b.round(4)}   tau = {t_hat:.4f}   loglik = {ll:.4f}")
# <</demo>>


# <<demopql>>
b_pql, t_pql, _ = fit_pql(y_demo, X_demo)
print(f"PQL        :  beta = {b_pql.round(4)}   tau = {t_pql:.4f}")
# <</demopql>>


def simulate(rng):
    x = rng.normal(size=(M, N))
    X = np.stack([np.ones((M, N)), x], axis=2)
    u = rng.normal(0, TAU, M)
    y = (rng.uniform(size=(M, N)) < expit(X @ BETA + u[:, None])).astype(float)
    return y, X


rng = np.random.default_rng(40_003)
rows = {"PQL": [], "Laplace": [], "AGHQ": []}
for _ in range(REPS):
    y, X = simulate(rng)
    b_p, t_p, _ = fit_pql(y, X)
    b_l, t_l, _ = fit_ml(y, X, 1)
    b_a, t_a, _ = fit_ml(y, X, 9, start=np.r_[b_l, np.log(t_l)])
    rows["PQL"].append([b_p[1], t_p])
    rows["Laplace"].append([b_l[1], t_l])
    rows["AGHQ"].append([b_a[1], t_a])
means = {k: np.mean(np.array(v), axis=0) for k, v in rows.items()}
for k, v in means.items():
    print(f"{k:8s} mean beta1 = {v[0]:.4f}   mean tau = {v[1]:.4f}")

# every method is biased downwards, and the order of the bias is PQL < Laplace < AGHQ
assert means["PQL"][0] < means["Laplace"][0] < means["AGHQ"][0] < BETA[1] * 1.05
assert means["PQL"][1] < means["Laplace"][1] < means["AGHQ"][1] < TAU * 1.02
assert means["AGHQ"][1] > 0.97 * TAU and means["PQL"][1] < 0.85 * TAU

# ---- accuracy of the quadrature itself, at one data set ----------------------
y, X = simulate(np.random.default_rng(40_004))
beta0, tau0, _ = fit_ml(y, X, 21)
exact = loglik(y, X, beta0, tau0, 61)
Ks = np.arange(1, 22, 2)
err_ad = np.array([abs(loglik(y, X, beta0, tau0, K) - exact) for K in Ks])
err_gh = np.array([abs(loglik(y, X, beta0, tau0, K, adapt=False) - exact) for K in Ks])
assert err_ad[2] < err_gh[2] / 10          # K = 5: adaptive is an order of magnitude better
assert err_ad[-1] < 1e-8
assert err_ad[6] < 1e-9                    # K = 13: the adaptive rule has converged
assert err_gh[9] > 1e-6                    # K = 19: the plain rule has not

# one adaptive node is the Laplace approximation, to machine precision
u_hat, sigma = cluster_mode(y, X, beta0, tau0)
eta = X @ beta0 + u_hat[:, None]
lap = np.sum((y * eta - np.logaddexp(0.0, eta)).sum(1)
             - 0.5 * (u_hat / tau0) ** 2 - np.log(tau0) + np.log(sigma))
assert abs(lap - loglik(y, X, beta0, tau0, 1)) < 1e-10

gen = Generated("ch40", "approximations")
gen.int("m", M)
gen.int("n", N)
gen.int("reps", REPS)
gen.num("beta_pql", means["PQL"][0], 4)
gen.num("beta_lap", means["Laplace"][0], 4)
gen.num("beta_aghq", means["AGHQ"][0], 4)
gen.num("tau_pql", means["PQL"][1], 4)
gen.num("tau_lap", means["Laplace"][1], 4)
gen.num("tau_aghq", means["AGHQ"][1], 4)
gen.num("bias_pql_pct", 100 * (means["PQL"][0] / BETA[1] - 1), 1)
gen.num("bias_tau_pql_pct", 100 * (means["PQL"][1] / TAU - 1), 1)
gen.num("err_ad5", err_ad[2], 5, sci=True)
gen.num("err_gh5", err_gh[2], 5, sci=True)
gen.num("err_ad13", err_ad[6], 3, sci=True)
gen.num("err_gh19", err_gh[9], 3, sci=True)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.4))
labels = ["PQL", "Laplace", "AGHQ"]
colours = [COLORS["second"], COLORS["thread"], COLORS["accent"]]

for ax, j, truth, name in [(axes[0], 0, BETA[1], r"$\hat\beta_1$"),
                           (axes[1], 1, TAU, r"$\hat\tau$")]:
    data = [np.array(rows[k])[:, j] for k in labels]
    bp = ax.boxplot(data, tick_labels=labels, widths=0.55, showfliers=False,
                    medianprops=dict(color=COLORS["ink"]))
    for patch, colour in zip(bp["boxes"], colours):
        patch.set_color(colour)
    ax.axhline(truth, color=COLORS["grid"], linewidth=0.8, zorder=0)
    ax.set_title(f"(a) {name}" if j == 0 else f"(b) {name}")
    ax.tick_params(axis="x", labelsize=7)

ax = axes[2]
ax.semilogy(Ks, np.maximum(err_ad, 1e-14), "o-", color=COLORS["accent"], markersize=3,
            label="adaptive")
ax.semilogy(Ks, np.maximum(err_gh, 1e-14), "s--", color=COLORS["second"], markersize=3,
            label="non-adaptive")
ax.set_xlabel("nodes $K$")
ax.set_ylabel("error in log-likelihood")
ax.set_title("(c) quadrature error")
ax.legend(frameon=False, fontsize=7)

fig.tight_layout()
fig.savefig(figure_path("ch40", "approximations"))
