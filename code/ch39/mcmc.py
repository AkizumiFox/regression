"""Chapter 39, Section 4: Metropolis-Hastings and a data-augmentation Gibbs sampler for a GLM.

Data: statsmodels.datasets.anes96 (public domain). The model is the logistic
regression of the intended vote on party identification, age, education and income
that Chapter 35 fits by maximum likelihood, now with a normal prior and a fully
Bayesian analysis; the probit version of the same model is sampled by the Albert-Chib
data-augmentation scheme.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
anes = sm.datasets.anes96.load_pandas().data
y = anes["vote"].to_numpy(float)                          # 1 if the respondent expected to vote Dole
X = np.column_stack([np.ones(len(y)), anes["PID"], anes["age"] / 10.0,
                     anes["educ"], anes["income"]])
n, p = X.shape
TAU = 10.0                                                # prior: beta ~ N(0, TAU^2 I)

mle = sm.GLM(y, X, family=sm.families.Binomial()).fit()
w = mle.fittedvalues * (1 - mle.fittedvalues)             # the working weights at the mode
V_prop = np.linalg.inv(X.T @ (w[:, None] * X) + np.eye(p) / TAU**2)   # Laplace covariance
CHOL = np.linalg.cholesky(V_prop)
print("maximum likelihood:", mle.params.round(4))
# <</data>>


def log_post(beta):
    eta = X @ beta
    return np.sum(y * eta - np.logaddexp(0.0, eta)) - np.sum(beta**2) / (2 * TAU**2)


# <<metropolis>>
def metropolis(n_draws, seed, start, step=2.38 / np.sqrt(5)):
    """Random-walk Metropolis with the Laplace covariance as the proposal shape."""
    rng = np.random.default_rng(seed)
    beta, lp = start.copy(), log_post(start)
    out, accepted = np.empty((n_draws, p)), 0
    for t in range(n_draws):
        prop = beta + step * (CHOL @ rng.standard_normal(p))
        lp_prop = log_post(prop)
        if np.log(rng.random()) < lp_prop - lp:           # the acceptance probability
            beta, lp, accepted = prop, lp_prop, accepted + 1
        out[t] = beta
    return out, accepted / n_draws


M = 8000                                                  # draws per chain
SD = np.sqrt(np.diag(V_prop))                             # posterior scale, coordinatewise
SIGNS = np.array([[1, 1, 1, 1, 1], [-1, -1, -1, -1, -1],
                  [1, -1, 1, -1, 1], [-1, 1, -1, 1, -1]], float)
starts = [mle.params + 8 * SD * g for g in SIGNS]         # over-dispersed in every coordinate
raw, rates = zip(*[metropolis(M, 3920 + c, s) for c, s in enumerate(starts)])
raw = np.array(raw)                                       # the chains as generated
chains = raw[:, M // 2:, :]                               # discard the first half
print(f"acceptance rate {np.mean(rates):.3f}")
# <</metropolis>>

assert 0.15 < np.mean(rates) < 0.45


# <<diagnostics>>
def split_rhat(draws):
    """Gelman-Rubin statistic on the two halves of every chain."""
    c, N = draws.shape
    s = draws.reshape(2 * c, N // 2)
    means, variances = s.mean(axis=1), s.var(axis=1, ddof=1)
    B, W = means.var(ddof=1), variances.mean()
    return np.sqrt(((N // 2 - 1) / (N // 2) * W + B) / W)


def ess(draws):
    """Effective sample size, with Geyer's initial positive sequence."""
    c, N = draws.shape
    centred = draws - draws.mean(axis=1, keepdims=True)
    var = draws.var(ddof=1)
    rho = np.array([np.mean([np.dot(x[:N - t], x[t:]) / N for x in centred]) / var
                    for t in range(1, min(N, 400))])
    pairs = rho[: 2 * (len(rho) // 2)].reshape(-1, 2).sum(axis=1)
    k = np.argmax(pairs < 0) if np.any(pairs < 0) else len(pairs)
    return c * N / (1 + 2 * pairs[:k].sum())


names = ["intercept", "party id", "age", "education", "income"]
for j in range(p):
    d = chains[:, :, j]
    print(f"{names[j]:10s} mean {d.mean():8.4f}  sd {d.std():7.4f}  "
          f"(MLE {mle.params[j]:8.4f}, se {mle.bse[j]:6.4f})  "
          f"R-hat {split_rhat(d):.4f}  ESS {ess(d):7.0f}")
# <</diagnostics>>

# <<short>>
short, short_rates = zip(*[metropolis(1500, 3960 + c, mle.params) for c in range(2)])
short = np.array(short)[:, 750:, :]
print(f"two short chains: acceptance {np.mean(short_rates):.3f}, "
      f"party id mean {short[:, :, 1].mean():.4f}, R-hat {split_rhat(short[:, :, 1]):.4f}, "
      f"ESS {ess(short[:, :, 1]):.0f}")
# <</short>>

rhat = np.array([split_rhat(chains[:, :, j]) for j in range(p)])
ess_mh = np.array([ess(chains[:, :, j]) for j in range(p)])
assert rhat.max() < 1.02
post_mean, post_sd = chains.reshape(-1, p).mean(axis=0), chains.reshape(-1, p).std(axis=0)
assert np.max(np.abs(post_mean - mle.params) / mle.bse) < 0.2


# <<gibbs>>
def albert_chib(n_draws, seed, start):
    """Gibbs for the probit model: latent normals, then a normal draw for beta."""
    rng = np.random.default_rng(seed)
    V = np.linalg.inv(X.T @ X + np.eye(p) / TAU**2)
    L = np.linalg.cholesky(V)
    beta, out = start.copy(), np.empty((n_draws, p))
    for t in range(n_draws):
        eta = X @ beta
        lower = stats.norm.cdf(-eta)                      # P(latent < 0)
        u = np.where(y > 0, lower + (1 - lower) * rng.random(n), lower * rng.random(n))
        z = eta + stats.norm.ppf(np.clip(u, 1e-12, 1 - 1e-12))
        mean = V @ (X.T @ z)
        beta = mean + L @ rng.standard_normal(p)
        out[t] = beta
    return out


G = 4000
probit_mle = sm.GLM(y, X, family=sm.families.Binomial(sm.families.links.Probit())).fit()
gibbs = np.array([albert_chib(G, 3930 + c, np.zeros(p)) for c in range(4)])[:, G // 2:, :]
print(f"probit, party id: posterior mean {gibbs[:, :, 1].mean():.4f} "
      f"(maximum likelihood {probit_mle.params[1]:.4f}), "
      f"ESS {ess(gibbs[:, :, 1]):.0f} of {gibbs[:, :, 1].size}")
# <</gibbs>>

# <<gibbsshort>>
gibbs_short = np.array([albert_chib(600, 3970 + c, np.zeros(p)) for c in range(2)])[:, 300:, :]
print(f"short Gibbs run: probit party id mean {gibbs_short[:, :, 1].mean():.4f}, "
      f"R-hat {split_rhat(gibbs_short[:, :, 1]):.4f}")
# <</gibbsshort>>

ess_gibbs = np.array([ess(gibbs[:, :, j]) for j in range(p)])
rhat_gibbs = np.array([split_rhat(gibbs[:, :, j]) for j in range(p)])
assert rhat_gibbs.max() < 1.02
eff_mh, eff_gibbs = ess_mh[1] / chains[:, :, 1].size, ess_gibbs[1] / gibbs[:, :, 1].size
assert eff_gibbs > 1.3 * eff_mh

# the logit-to-probit scale factor, empirically
ratio = post_mean[1] / gibbs[:, :, 1].mean()

# <<ppc>>
draws = chains.reshape(-1, p)
rng = np.random.default_rng(3940)
idx = rng.choice(len(draws), 4000, replace=False)
group = anes["PID"].to_numpy() == 2                       # the independent Democrats
observed = y[group].sum()
replicated = np.array([rng.binomial(1, special.expit(X @ draws[i]))[group].sum()
                       for i in idx])
p_value = np.mean(replicated <= observed)
print(f"independent Democrats: {int(observed)} of {int(group.sum())} expected to vote Dole, "
      f"against {replicated.mean():.1f} replicated; "
      f"posterior predictive p-value {p_value:.3f}")
# <</ppc>>

assert p_value < 0.05

# <<ppcshort>>
group = anes["PID"].to_numpy() == 2                       # the independent Democrats
observed = y[group].sum()
short4 = np.array([metropolis(2000, 3980 + c, mle.params)[0] for c in range(4)])[:, 1000:, :]
draws_short = short4.reshape(-1, p)
rng2 = np.random.default_rng(3990)
idx2 = rng2.choice(len(draws_short), 2000, replace=False)
rep_short = np.array([rng2.binomial(1, special.expit(X @ draws_short[i]))[group].sum()
                      for i in idx2])
print(f"short run: {int(observed)} of {int(group.sum())} against {rep_short.mean():.1f} "
      f"replicated; posterior predictive p-value {np.mean(rep_short <= observed):.3f}")
# <</ppcshort>>

# a quadratic in party identification does not explain the misfit
Xquad = np.column_stack([X, anes["PID"].to_numpy() ** 2])
lr_quad = mle.deviance - sm.GLM(y, Xquad, family=sm.families.Binomial()).fit().deviance
assert lr_quad < 0.1

# the same signal from the likelihood: party identification as a seven-level factor
D = np.column_stack([anes["PID"].to_numpy() == k for k in range(1, 7)]).astype(float)
Xf = np.column_stack([np.ones(n), D, anes["age"] / 10.0, anes["educ"], anes["income"]])
lr_factor = mle.deviance - sm.GLM(y, Xf, family=sm.families.Binomial()).fit().deviance
print(f"factor versus linear party identification: deviance drop {lr_factor:.2f} on 5 d.f., "
      f"p = {stats.chi2.sf(lr_factor, 5):.3f}")
assert lr_factor > 10

# R-hat as a function of how long the four chains are run: at each length the usual
# recipe is applied from scratch, discarding the first half of what was generated.
# This is the curve panel (b) draws, and it starts at the over-dispersed values.
lengths = np.arange(200, M + 1, 200)
rhat_curve = np.array([split_rhat(raw[:, k // 2:k, 1]) for k in lengths])
stays_below = np.array([bool(np.all(rhat_curve[k:] < 1.01)) for k in range(len(lengths))])
cross = int(lengths[np.argmax(stays_below)])
below105 = int(lengths[np.argmax(rhat_curve < 1.05)])
below102 = int(lengths[np.argmax(rhat_curve < 1.02)])
print(f"R-hat curve: {rhat_curve[0]:.2f} at {lengths[0]} draws per chain, "
      f"below 1.05 from {below105}, below 1.02 from {below102}, "
      f"permanently below 1.01 from {cross}")
assert stays_below.any() and rhat_curve[0] > 2.0 and rhat_curve[-1] < 1.01

gen = Generated("ch39", "mcmc")
gen.int("cross", cross)
gen.int("below105", below105)
gen.int("below102", below102)
gen.num("rhatstart", float(rhat_curve[0]), 2)
gen.int("n", n)
gen.int("p", p)
gen.int("M", M)
gen.int("G", G)
gen.int("kept", chains.shape[0] * chains.shape[1])
gen.num("tau", TAU, 0)
gen.num("rate", float(np.mean(rates)), 3)
gen.num("rhatmax", float(rhat.max()), 4)
gen.num("essmh", float(ess_mh[1]), 0)
gen.num("essgibbs", float(ess_gibbs[1]), 0)
gen.int("keptgibbs", int(gibbs[:, :, 1].size))
gen.num("effmh", float(eff_mh), 3)
gen.num("effgibbs", float(eff_gibbs), 3)
gen.num("meanpid", post_mean[1], 4)
gen.num("sdpid", post_sd[1], 4)
gen.num("mlepid", mle.params[1], 4)
gen.num("sepid", mle.bse[1], 4)
gen.num("probitpid", float(gibbs[:, :, 1].mean()), 4)
gen.num("probitmle", probit_mle.params[1], 4)
gen.num("ratio", ratio, 3)
gen.int("middle", int(group.sum()))
gen.num("lrfactor", lr_factor, 2)
gen.num("lrquad", lr_quad, 3)
gen.num("plrfactor", float(stats.chi2.sf(lr_factor, 5)), 3)
gen.int("observed", int(observed))
gen.num("expected", float(replicated.mean()), 1)
gen.num("ppc", p_value, 3)
gen.write()

# ---- Figure: traces, R-hat against chain length, autocorrelation -----------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.2))
ax = axes[0]
for c in range(4):
    ax.plot(raw[c, :1200, 1], linewidth=0.5, alpha=0.85,
            color=[COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["thread"]][c])
ax.set_xlabel("iteration")
ax.set_ylabel(r"$\beta_{\mathrm{party}}$")
ax.set_title("(a) four chains")

ax = axes[1]
ax.semilogy(lengths, rhat_curve - 1.0, color=COLORS["accent"])
ax.axhline(0.01, color=COLORS["muted"], linestyle=":")
ax.set_xlabel("draws per chain generated")
ax.set_ylabel(r"$\hat R-1$")
ax.set_title("(b) convergence")

ax = axes[2]
for d, lab, col, ls in [(chains[:, :, 1], "Metropolis", COLORS["accent"], "-"),
                        (gibbs[:, :, 1], "Gibbs (probit)", COLORS["second"], "--")]:
    cen = d - d.mean(axis=1, keepdims=True)
    v = d.var(ddof=1)
    acf = [np.mean([np.dot(x[: x.size - t], x[t:]) / x.size for x in cen]) / v
           for t in range(60)]
    ax.plot(acf, color=col, linestyle=ls, label=lab)
ax.axhline(0, color=COLORS["grid"], linewidth=0.6)
ax.set_xlabel("lag")
ax.set_ylabel("autocorrelation")
ax.set_title("(c) mixing")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch39", "mcmc_diagnostics"))

# ---- Figure: the posterior predictive check --------------------------------
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3))
ax = axes[0]
ax.hist(replicated, bins=np.arange(replicated.min() - 0.5, replicated.max() + 1.5),
        color=COLORS["grid"], edgecolor="white", linewidth=0.3)
ax.axvline(observed, color=COLORS["second"])
ax.set_xlabel("Dole votes among the independent Democrats")
ax.set_ylabel("replications")
ax.set_title("(a) posterior predictive check")

ax = axes[1]
pid = np.arange(7)
prop = np.array([y[anes["PID"] == k].mean() for k in pid])
rows = np.column_stack([np.ones(7), pid, np.full(7, anes["age"].mean() / 10),
                        np.full(7, anes["educ"].mean()), np.full(7, anes["income"].mean())])
curves = special.expit(draws[idx] @ rows.T)
band = np.quantile(curves, [0.025, 0.975], axis=0)
ax.fill_between(pid, band[0], band[1], color=COLORS["grid"], alpha=0.85, linewidth=0)
ax.plot(pid, curves.mean(axis=0), color=COLORS["accent"])
ax.scatter(pid, prop, s=12, color=COLORS["second"], zorder=3)
ax.set_xlabel("party identification")
ax.set_ylabel("probability of a Dole vote")
ax.set_title("(b) fit at the average covariates")
fig.tight_layout()
fig.savefig(figure_path("ch39", "posterior_predictive"))
