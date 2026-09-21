"""Chapter 39, Section 5: priors as regularizers, and a prior that cures separation.

Part one revisits the separated subgroup of Chapter 35 (the thirteen respondents of
the 1996 American National Election Study whose schooling stopped at or before the
eighth grade) and computes the posterior of the party-identification slope under a
flat, a normal and a Cauchy prior.  The quadrature runs over the whole plane, by the
substitution b = s tan(u): under separation the likelihood tends to a positive
constant along the separating cone, so a Cauchy posterior has a polynomial tail and
any fixed box would throw away a large part of its mass.  Part two checks the two
correspondences between priors and penalties on the full survey.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<separated>>
anes = sm.datasets.anes96.load_pandas().data
sub = anes[anes["educ"] == 1]                          # schooling stopped at or before grade eight
Xs = np.column_stack([np.ones(len(sub)), sub["PID"]])
ys = sub["vote"].to_numpy(float)


def loglik(b0, b1):
    """Logistic log-likelihood at one intercept and a whole grid of slopes."""
    eta = b0 + np.outer(b1, Xs[:, 1])
    return np.sum(ys * eta - np.logaddexp(0.0, eta), axis=1)


def log_prior(b0, b1, kind, scale=2.5):
    if kind == "flat":
        return 0.0 * b0
    if kind == "normal":
        return -(b0 / 10.0) ** 2 / 2 - (b1 / scale) ** 2 / 2
    return -np.log1p((b0 / 10.0) ** 2) - np.log1p((b1 / scale) ** 2)   # Cauchy


# With b = s tan(u) and u in (-pi/2, pi/2) the plane becomes the bounded square
# (-pi/2, pi/2)^2, the Jacobian s/cos^2(u) carries the tails, and nothing is truncated.
# Gauss-Legendre nodes in u are used because under a Cauchy prior the integrand of the
# posterior mean tends to a nonzero limit at the ends of the interval, which a rule that
# ignores the endpoints would miss.
def posterior_of_slope(kind, k=1200, s0=30.0, s1=6.0):
    """Marginal posterior of the slope, by quadrature over the whole plane."""
    t, wt = np.polynomial.legendre.leggauss(k)
    u, wu = t * np.pi / 2, wt * np.pi / 2
    b0, b1 = s0 * np.tan(u), s1 * np.tan(u)
    j0, j1 = s0 / np.cos(u) ** 2, s1 / np.cos(u) ** 2
    logq = np.array([loglik(a, b1) for a in b0]) + log_prior(b0[:, None], b1[None, :], kind)
    q = np.exp(logq - logq.max())
    marg = ((j0 * wu)[:, None] * q).sum(axis=0)        # integrate out the intercept
    dens = marg / np.sum(marg * j1 * wu)               # density with respect to b1
    return u, b1, dens, j1, wu


def summarize(u, b1, dens, j1, wu, probs=(0.025, 0.5, 0.975)):
    """Posterior mean and quantiles of the slope from the quadrature grid."""
    f = dens * j1                                      # the density in the u coordinate
    cdf = np.concatenate([[0.0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(u))])
    return np.sum(b1 * f * wu), np.interp(probs, cdf / cdf[-1], b1)


post_n, post_c = posterior_of_slope("normal"), posterior_of_slope("cauchy")
mean_n, (lo_n, med_n, hi_n) = summarize(*post_n)
mean_c, (lo_c, med_c, hi_c) = summarize(*post_c)
print(f"posterior mean of the slope: normal {mean_n:.3f}, Cauchy {mean_c:.3f}")
# <</separated>>

assert 0 < lo_n < med_n < hi_n and 0 < lo_c < med_c < hi_c
assert hi_c > 3 * hi_n                                  # the Cauchy tail leaves far more room

# The quadrature is converged: refining the grid and changing the substitution scale
# must not move the mean or the upper endpoint in the third decimal place.
for k, s0, s1 in [(2000, 30.0, 6.0), (1200, 45.0, 10.0)]:
    m, (_, _, h) = summarize(*posterior_of_slope("cauchy", k=k, s0=s0, s1=s1))
    assert abs(m - mean_c) < 1e-3 and abs(h - hi_c) < 5e-3, (k, s0, s1)


# Firth's estimate on the same subgroup, recomputed here so the figure does not
# hard-code the value Chapter 35 reports
def firth(X, y, steps=300):
    beta = np.zeros(X.shape[1])
    for _ in range(steps):
        mu = special.expit(X @ beta)
        w = mu * (1 - mu)
        F = X.T @ (w[:, None] * X)
        h = w * np.einsum("ij,jk,ik->i", X, np.linalg.inv(F), X)
        beta = beta + np.linalg.solve(F, X.T @ (y - mu + h * (0.5 - mu)))
    return beta


firth_slope = firth(Xs, ys)[1]
assert np.isclose(round(firth_slope, 3), 1.362), firth_slope   # Chapter 35, Example 35.3.3

# the flat prior gives an improper posterior: the mass inside a square box centred at
# the origin grows without bound, like the area of the separating cone
boxes = [10.0, 20.0, 40.0, 80.0]
masses = []
for half in boxes:
    b = np.linspace(-half, half, 1200)
    logq = np.array([loglik(a, b) for a in b])
    masses.append(np.trapezoid(np.trapezoid(np.exp(logq), b, axis=1), b))
masses = np.array(masses)
growth = masses[1:] / masses[:-1]
print("flat prior, mass inside a square box of half-width", boxes, ":", np.round(masses, 1))
assert np.all(growth > 1.8)                             # doubling the box more than doubles the mass

# <<penalties>>
Xf = np.column_stack([anes["PID"], anes["age"], anes["educ"], anes["income"]])
Xf = (Xf - Xf.mean(0)) / Xf.std(0)                      # standardize, as any penalty requires
Xf = np.column_stack([np.ones(len(Xf)), Xf])
yf = anes["vote"].to_numpy(float)


def map_estimate(penalty, lam, steps=30):
    """Posterior mode: maximize the log-likelihood minus lam times a penalty on the slopes.

    A normal prior gives the ridge penalty and is solved by Newton's method; a Laplace
    prior gives the lasso penalty and is solved by coordinate descent on the working
    least squares problem of each Newton step.
    """
    beta = np.zeros(Xf.shape[1])
    for _ in range(steps):
        pi = special.expit(Xf @ beta)
        w = np.maximum(pi * (1 - pi), 1e-8)
        z = Xf @ beta + (yf - pi) / w                   # the working response
        if penalty == "ridge":
            P = np.diag([0.0] + [1.0] * (Xf.shape[1] - 1))
            beta = np.linalg.solve(Xf.T @ (w[:, None] * Xf) + lam * P, Xf.T @ (w * z))
            continue
        for _ in range(40):                             # coordinate descent for the lasso
            for j in range(Xf.shape[1]):
                r = z - Xf @ beta + Xf[:, j] * beta[j]
                num, den = np.sum(w * Xf[:, j] * r), np.sum(w * Xf[:, j] ** 2)
                beta[j] = (num / den if j == 0
                           else np.sign(num) * max(abs(num) - lam, 0.0) / den)
    return beta


unpenalized = map_estimate("ridge", 0.0)                # the maximum likelihood fit
lams = np.array([1.0, 4.0, 16.0, 64.0, 256.0])
ridge_path = np.array([map_estimate("ridge", l) for l in lams])
lasso_path = np.array([map_estimate("lasso", l) for l in lams])
print(f"unpenalized slopes: {np.round(unpenalized[1:], 3)}")
print("normal prior with sd tau = 1/sqrt(lam)  ->  ridge;  Laplace prior with b = 1/lam  ->  lasso")
for k, l in enumerate(lams):
    print(f"  lam {l:6.1f}  ridge {np.round(ridge_path[k, 1:], 3)}   "
          f"lasso {np.round(lasso_path[k, 1:], 3)}")
# <</penalties>>

nonzero = np.sum(np.abs(lasso_path[-1, 1:]) > 1e-3)
assert nonzero < 4 and np.all(np.abs(ridge_path[-1, 1:]) > 1e-3)

# the two correspondences, checked against an independent implementation
lam_check = 16.0
tau = 1.0 / np.sqrt(lam_check)


def neg_log_posterior(beta, lam, penalty):
    eta = Xf @ beta
    pen = np.sum(beta[1:] ** 2) / 2 if penalty == "ridge" else np.sum(np.abs(beta[1:]))
    return -np.sum(yf * eta - np.logaddexp(0.0, eta)) + lam * pen


for penalty in ("ridge", "lasso"):
    b = map_estimate(penalty, lam_check)
    best = optimize.minimize(neg_log_posterior, b, args=(lam_check, penalty),
                             method="Nelder-Mead",
                             options={"xatol": 1e-9, "fatol": 1e-9, "maxiter": 50000})
    assert neg_log_posterior(b, lam_check, penalty) <= best.fun + 1e-6, penalty

# the induced prior on the probability scale: logit-normal, bimodal exactly when sigma^2 > 2
pgrid = np.linspace(1e-4, 1 - 1e-4, 2000)


def logit_normal_density(sigma):
    t = special.logit(pgrid)
    return np.exp(-(t / sigma) ** 2 / 2) / (sigma * np.sqrt(2 * np.pi) * pgrid * (1 - pgrid))


def stationary_points(sigma):
    """Roots of t/sigma^2 = 2 expit(t) - 1: the modes of the induced prior."""
    t = np.linspace(-40, 40, 400000)      # a grid that avoids the exact root at zero
    g = t / sigma**2 - (2 * special.expit(t) - 1)
    return int(np.sum(np.diff(np.sign(g)) != 0))


modes = {s: stationary_points(s) for s in (1.0, np.sqrt(2.0), 1.5, 4.0)}
print("stationary points of the induced prior on the probability:", modes)
assert modes[1.0] == 1 and modes[np.sqrt(2.0)] == 1 and modes[1.5] == 3 and modes[4.0] == 3

gen = Generated("ch39", "priors")
gen.int("nsub", len(sub))
gen.num("meannormal", mean_n, 3)
gen.num("meancauchy", mean_c, 3)
gen.num("mednormal", med_n, 3)
gen.num("medcauchy", med_c, 3)
gen.num("lonormal", lo_n, 3)
gen.num("hinormal", hi_n, 3)
gen.num("locauchy", lo_c, 3)
gen.num("hicauchy", hi_c, 3)
gen.num("hiratio", hi_c / hi_n, 1)
# half the mass of a N(0, 10^2) prior on the logit lies outside this probability interval
tail10 = float(special.expit(-stats.norm.ppf(0.75) * 10.0))
assert 0.001 < tail10 < 0.002
gen.num("tail10", tail10, 4)
gen.num("firth", firth_slope, 3)
gen.num("growth", float(growth.min()), 2)
gen.num("mass10", masses[0], 1)
gen.num("mass80", masses[-1], 1)
gen.num("tau", tau, 3)
gen.num("lamcheck", lam_check, 0)
gen.num("lamquote", lams[3], 0)
gen.num("taucheck", float(1 / np.sqrt(lams[3])), 3)
gen.int("nonzero", int(nonzero))
gen.num("unpen_pid", unpenalized[1], 3)
for k, name in enumerate(["pid", "age", "educ", "income"]):
    gen.num("ridge_" + name, ridge_path[3, k + 1], 3)
    gen.num("lasso_" + name, lasso_path[3, k + 1], 3)
gen.write()

# ---- (a) the induced prior on the probability, (b) posteriors under three priors ----
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
for s, col, ls in [(1.0, COLORS["accent"], "-"), (np.sqrt(2.0), COLORS["third"], "--"),
                   (4.0, COLORS["second"], "-.")]:
    ax.plot(pgrid, logit_normal_density(s), color=col, linestyle=ls,
            label=r"$\sigma=%.2f$" % s)
ax.set_ylim(0, 4)
ax.set_xlabel(r"probability $\pi$")
ax.set_ylabel("prior density")
ax.set_title("(a) what a normal prior implies")
ax.legend(frameon=False)

ax = axes[1]
ax.plot(post_n[1], post_n[2], color=COLORS["accent"], label=r"normal, sd $2.5$")
ax.plot(post_c[1], post_c[2], color=COLORS["second"], linestyle="--", label=r"Cauchy, scale $2.5$")
ax.axvline(firth_slope, color=COLORS["muted"], linestyle=":", label="Firth")
ax.set_xlim(-1, 26)                                     # wide enough to show the Cauchy tail
ax.set_xlabel(r"slope $\beta_1$")
ax.set_ylabel("posterior density")
ax.set_title("(b) the separated subgroup")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(figure_path("ch39", "prior_separation"))
