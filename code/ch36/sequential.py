"""Chapter 36, Section 4: sequential (continuation-ratio) models.

Three things are checked here.

1. The factorization: the multinomial likelihood of a sequential model is a product of
   c - 1 binomial likelihoods, so with cut-specific slopes the maximum likelihood fit
   is c - 1 separate binary logistic regressions, and with a common slope it is one
   binary logistic regression on the expanded ("one row per stage") data set.
2. The Laara-Matthews identity: with the complementary log-log link the sequential and
   cumulative models are the same family, the cutpoints being related by a cumulative sum
   on the exponential scale.
3. Sequential and cumulative models are genuinely different families under the logit link;
   a simulation from a sequential mechanism is fitted with both.

The self-placement data are the public-domain 1996 American National Election Study
shipped with statsmodels; see code/ch36/ordinal.py.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].astype(int).to_numpy()
X = np.column_stack([anes["age"] / 10.0, anes["educ"], anes["income"]])
n, p = X.shape
c = 7


def seq_probs(theta, beta, X, F=lambda z: 1.0 / (1.0 + np.exp(-z))):
    """Category probabilities of the sequential model with P(Y = r | Y >= r) = F(theta_r - x'beta)."""
    delta = F(theta[None, :] - (X @ beta)[:, None])
    surv = np.column_stack([np.ones(len(X)), np.cumprod(1.0 - delta, axis=1)])
    return np.column_stack([delta * surv[:, :-1], surv[:, -1]])


def seq_loglik(theta, beta, X, y, F=lambda z: 1.0 / (1.0 + np.exp(-z))):
    P = seq_probs(theta, beta, X, F)
    return float(np.sum(np.log(P[np.arange(len(y)), y - 1])))
# <</data>>


# <<separate>>
# unrestricted sequential model: one binary logistic fit per stage, on those still at risk
stage_beta, stage_ll = [], 0.0
for r in range(1, c):
    at_risk = y >= r
    Z = np.column_stack([np.ones(at_risk.sum()), -X[at_risk]])
    fit_r = sm.Logit((y[at_risk] == r).astype(float), Z).fit(disp=0)
    stage_beta.append(np.asarray(fit_r.params))
    stage_ll += fit_r.llf
    print(f"stage {r}: at risk {at_risk.sum():4d}, stop {int((y == r).sum()):4d}, "
          f"theta {fit_r.params[0]:7.4f}, beta {np.round(fit_r.params[1:], 4)}")

print(f"sum of the {c - 1} binary log-likelihoods: {stage_ll:.4f}")
# <</separate>>

stage_beta = np.array(stage_beta)
theta_free, beta_free = stage_beta[:, 0], stage_beta[:, 1:]

# the same number is the multinomial log-likelihood of the sequential fit
delta_free = 1.0 / (1.0 + np.exp(-(theta_free[None, :] - np.einsum("ij,rj->ir", X, beta_free))))
surv_free = np.column_stack([np.ones(n), np.cumprod(1.0 - delta_free, axis=1)])
P_free = np.column_stack([delta_free * surv_free[:, :-1], surv_free[:, -1]])
ll_free = float(np.sum(np.log(P_free[np.arange(n), y - 1])))
assert np.isclose(ll_free, stage_ll, atol=1e-8)
assert np.allclose(P_free.sum(axis=1), 1.0)

# <<expanded>>
# common slope: one binary logistic regression on the expanded data set
rows_X, rows_y = [], []
for i in range(n):
    for r in range(1, min(y[i], c - 1) + 1):        # stages this respondent reached
        stage = np.zeros(c - 1)
        stage[r - 1] = 1.0
        rows_X.append(np.concatenate([stage, -X[i]]))
        rows_y.append(float(y[i] == r))
Xe, ye = np.array(rows_X), np.array(rows_y)

fit_seq = sm.Logit(ye, Xe).fit(disp=0)
theta_seq, beta_seq = fit_seq.params[:c - 1], fit_seq.params[c - 1:]
se_seq = fit_seq.bse[c - 1:]
ll_seq = seq_loglik(theta_seq, beta_seq, X, y)

print(f"expanded data: {len(ye)} rows from {n} respondents")
print("cutpoints", np.round(theta_seq, 4))
print("slopes   ", np.round(beta_seq, 4), " se", np.round(se_seq, 4))
print(f"binary log-likelihood {fit_seq.llf:.4f} = multinomial log-likelihood {ll_seq:.4f}")
# <</expanded>>

assert np.isclose(fit_seq.llf, ll_seq, atol=1e-8)
aic_seq = -2 * ll_seq + 2 * (c - 1 + p)

# <<cloglog>>
# Laara and Matthews: with a complementary log-log link the sequential model IS a cumulative model
cll = lambda z: 1.0 - np.exp(-np.exp(z))
fit_cll = sm.GLM(ye, Xe, family=sm.families.Binomial(sm.families.links.CLogLog())).fit()
theta_t, beta_c = fit_cll.params[:c - 1], fit_cll.params[c - 1:]
theta_cum = np.log(np.cumsum(np.exp(theta_t)))                   # the cumulative model's cutpoints

P_step = seq_probs(theta_t, beta_c, X, F=cll)                    # built stage by stage
gamma_cum = cll(theta_cum[None, :] - (X @ beta_c)[:, None])      # built from cumulative probabilities
P_cum = np.diff(np.column_stack([np.zeros(n), gamma_cum, np.ones(n)]), axis=1)

print("sequential cutpoints ", np.round(theta_t, 4))
print("cumulative cutpoints ", np.round(theta_cum, 4))
print("largest difference between the two sets of fitted probabilities:",
      float(np.max(np.abs(P_step - P_cum))))
# <</cloglog>>

assert np.allclose(P_step, P_cum, atol=1e-12)
assert np.all(np.diff(theta_cum) > 0)

# <<simulate>>
# a mechanism that really is sequential, fitted by both models
def simulate_sequential(n, theta, beta, rng):
    X = rng.normal(size=(n, len(beta)))
    delta = 1.0 / (1.0 + np.exp(-(theta[None, :] - (X @ beta)[:, None])))
    stop = rng.random(delta.shape) < delta
    y = np.where(stop.any(axis=1), stop.argmax(axis=1) + 1, len(theta) + 1)
    return X, y


rng = np.random.default_rng(20360436)
theta0, beta0 = np.array([-0.6, 0.0, 0.5]), np.array([1.0])
Xs, ys = simulate_sequential(4000, theta0, beta0, rng)
print("simulated counts:", np.bincount(ys)[1:])
# <</simulate>>

# sequential fit of the simulated data
rows_X, rows_y = [], []
for i in range(len(ys)):
    for r in range(1, min(ys[i], 3) + 1):
        stage = np.zeros(3)
        stage[r - 1] = 1.0
        rows_X.append(np.concatenate([stage, -Xs[i]]))
        rows_y.append(float(ys[i] == r))
fit_s = sm.Logit(np.array(rows_y), np.array(rows_X)).fit(disp=0)
ll_s = seq_loglik(fit_s.params[:3], fit_s.params[3:], Xs, ys)
aic_s = -2 * ll_s + 2 * 4

# cumulative logit fit of the same data, by Fisher scoring as in Section 36.2
def cum_fit_simple(X, y, c, maxit=100, tol=1e-11):
    n, p = X.shape
    Y = np.zeros((n, c))
    Y[np.arange(n), y - 1] = 1.0
    psi = np.concatenate([np.linspace(-1.0, 1.0, c - 1), np.zeros(p)])
    for _ in range(maxit):
        G = 1.0 / (1.0 + np.exp(-(psi[:c - 1][None, :] - (X @ psi[c - 1:])[:, None])))
        lam = G * (1.0 - G)
        dG = np.zeros((n, c - 1, c - 1 + p))
        for r in range(c - 1):
            dG[:, r, r] = lam[:, r]
            dG[:, r, c - 1:] = -lam[:, [r]] * X
        P = np.diff(np.column_stack([np.zeros(n), G, np.ones(n)]), axis=1)
        D = np.concatenate([dG[:, :1, :], np.diff(dG, axis=1)], axis=1)
        Pm = P[:, :c - 1]
        Sigma = np.einsum("ir,rs->irs", Pm, np.eye(c - 1)) - np.einsum("ir,is->irs", Pm, Pm)
        A = np.linalg.solve(Sigma, D)
        step = np.linalg.solve(np.einsum("irk,irl->kl", D, A),
                               np.einsum("irk,ir->k", A, Y[:, :c - 1] - Pm))
        psi = psi + step
        if np.max(np.abs(step)) < tol:
            break
    G = 1.0 / (1.0 + np.exp(-(psi[:c - 1][None, :] - (X @ psi[c - 1:])[:, None])))
    P = np.diff(np.column_stack([np.zeros(n), G, np.ones(n)]), axis=1)
    return psi, float(np.sum(Y * np.log(P)))


psi_c, ll_c = cum_fit_simple(Xs, ys, 4)
aic_c = -2 * ll_c + 2 * 4
assert abs(fit_s.params[3] - beta0[0]) < 3 * fit_s.bse[3]        # the sequential fit recovers beta
assert aic_s < aic_c - 10                                        # and fits much better

gen = Generated("ch36", "sequential")
gen.int("n", n)
gen.int("rows", len(ye))
for r in range(c - 1):
    gen.num(f"theta{r + 1}", theta_seq[r], 4)
    gen.int(f"risk{r + 1}", int((y >= r + 1).sum()))
for j, key in enumerate(["age", "educ", "income"]):
    gen.num(f"beta_{key}", beta_seq[j], 4)
    gen.num(f"se_{key}", se_seq[j], 4)
gen.num("ll_seq", ll_seq, 4)
gen.num("ll_free", ll_free, 4)
gen.num("aic_seq", aic_seq, 2)
gen.num("aic_free", -2 * ll_free + 2 * (c - 1) * (1 + p), 2)
gen.int("nsim", len(ys))
gen.num("beta_sim_seq", fit_s.params[3], 4)
gen.num("se_sim_seq", fit_s.bse[3], 4)
gen.num("beta_sim_cum", psi_c[3], 4)
gen.num("aic_sim_seq", aic_s, 2)
gen.num("aic_sim_cum", aic_c, 2)
gen.write()

# ---- figure: each model is a curved surface in the other's coordinates -------
use_book_style()
grid = np.linspace(-2.5, 2.5, 200)[:, None]
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))

P_seq = seq_probs(theta0, beta0, grid)
Gam = np.cumsum(P_seq, axis=1)[:, :3]
ax = axes[0]
for r in range(3):
    ax.plot(grid[:, 0], np.log(Gam[:, r] / (1 - Gam[:, r])), color=COLORS["accent"])
    ax.annotate(f"$r={r + 1}$", (grid[-1, 0], np.log(Gam[-1, r] / (1 - Gam[-1, r]))),
                textcoords="offset points", xytext=(3, -3), fontsize=7, color=COLORS["muted"])
ax.set_xlabel("$x$")
ax.set_ylabel(r"logit $P(Y \leq r)$")
ax.set_title("(a) a sequential model")
ax.set_xlim(-2.5, 3.2)

Gc = 1.0 / (1.0 + np.exp(-(np.array([-0.6, 0.3, 1.2])[None, :] - grid @ beta0[:, None])))
Pc = np.diff(np.column_stack([np.zeros(len(grid)), Gc, np.ones(len(grid))]), axis=1)
Del = Pc[:, :3] / (1.0 - np.column_stack([np.zeros(len(grid)), Gc[:, :2]]))
ax = axes[1]
for r in range(3):
    ax.plot(grid[:, 0], np.log(Del[:, r] / (1 - Del[:, r])), color=COLORS["second"])
    ax.annotate(f"$r={r + 1}$", (grid[-1, 0], np.log(Del[-1, r] / (1 - Del[-1, r]))),
                textcoords="offset points", xytext=(3, -3), fontsize=7, color=COLORS["muted"])
ax.set_xlabel("$x$")
ax.set_ylabel(r"logit $P(Y = r \mid Y \geq r)$")
ax.set_title("(b) a proportional odds model")
ax.set_xlim(-2.5, 3.2)
fig.tight_layout()
fig.savefig(figure_path("ch36", "sequential_vs_cumulative"))
