"""Chapter 36, Sections 2 and 3: cumulative logit models, proportional odds and its checks.

Response: self-placement on the seven-point liberal-conservative scale of the 1996
American National Election Study (1 = extremely liberal, 7 = extremely conservative).
Regressors: age in decades, the seven-point education scale, the twenty-four-point
household income scale. Public-domain data shipped with statsmodels
(statsmodels.datasets.anes96).

The cumulative link model is fitted by Fisher scoring on the multinomial likelihood,
written so that any subset of the regressors may be given cut-specific slopes; that
one routine gives the proportional odds fit, the partial proportional odds fit, the
fully non-proportional fit and the score test of proportional odds.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
NAMES = ["age (decades)", "education", "income"]

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].astype(int).to_numpy()            # 1 = extremely liberal ... 7 = extremely conservative
X = np.column_stack([anes["age"] / 10.0, anes["educ"], anes["income"]])
n, p = X.shape
c = 7
Y = np.zeros((n, c))
Y[np.arange(n), y - 1] = 1.0

print("category counts:", Y.sum(axis=0))
# <</data>>


# <<fit>>
def cum_gamma(psi, X, c, free):
    """Cumulative probabilities Lambda(theta_r - x^T beta_r) and their gradients.

    psi = [theta (c-1), beta (p), delta ((c-2) x len(free))]; the columns of X listed
    in `free` get cut-specific slopes beta_r = beta + delta_r, with delta_1 = 0.
    """
    n, p = X.shape
    q, nf, fl = len(psi), len(free), list(free)
    theta, beta = psi[:c - 1], psi[c - 1:c - 1 + p]
    delta = psi[c - 1 + p:].reshape(c - 2, nf)
    G = np.empty((n, c - 1))
    dG = np.zeros((n, c - 1, q))
    for r in range(c - 1):
        beta_r = beta.copy()
        if r >= 1 and nf:
            beta_r[fl] = beta_r[fl] + delta[r - 1]
        g = 1.0 / (1.0 + np.exp(-(theta[r] - X @ beta_r)))
        lam = (g * (1.0 - g))[:, None]              # the logistic density at the cut
        G[:, r] = g
        dG[:, r, r] = lam[:, 0]
        dG[:, r, c - 1:c - 1 + p] = -lam * X
        if r >= 1 and nf:
            off = c - 1 + p + (r - 1) * nf
            dG[:, r, off:off + nf] = -lam * X[:, fl]
    return G, dG


def cum_probs(psi, X, c, free):
    """Category probabilities pi_ir and their derivatives D[i, r, k] = d pi_ir / d psi_k."""
    G, dG = cum_gamma(psi, X, c, free)
    P = np.diff(np.column_stack([np.zeros(len(X)), G, np.ones(len(X))]), axis=1)
    D = np.concatenate([dG[:, :1, :], np.diff(dG, axis=1)], axis=1)
    return P, D


def cum_score(psi, X, Y, c, free):
    """Score and expected information of the multinomial likelihood."""
    P, D = cum_probs(psi, X, c, free)
    Pm = P[:, :c - 1]
    Sigma = np.einsum("ir,rs->irs", Pm, np.eye(c - 1)) - np.einsum("ir,is->irs", Pm, Pm)
    A = np.linalg.solve(Sigma, D)                   # Sigma_i^{-1} D_i, one observation at a time
    U = np.einsum("irk,ir->k", A, Y[:, :c - 1] - Pm)
    J = np.einsum("irk,irl->kl", D, A)
    return U, J


def cum_loglik(psi, X, Y, c, free):
    P, _ = cum_probs(psi, X, c, free)
    return float(np.sum(Y * np.log(np.clip(P, 1e-300, None))))


def cum_fit(X, Y, c, free=(), tol=1e-10, maxit=100):
    """Fisher scoring for the cumulative logit model."""
    n, p = X.shape
    psi = np.zeros((c - 1) + p + (c - 2) * len(free))
    psi[:c - 1] = np.linspace(-1.5, 1.5, c - 1)     # start from ordered cutpoints, no slopes
    for it in range(1, maxit + 1):
        U, J = cum_score(psi, X, Y, c, free)
        step = np.linalg.solve(J, U)
        psi = psi + step
        if np.max(np.abs(step)) < tol:
            break
    return psi, J, it


psi, J, iters = cum_fit(X, Y, c)
theta, beta = psi[:c - 1], psi[c - 1:]
se = np.sqrt(np.diag(np.linalg.inv(J)))
ll_po = cum_loglik(psi, X, Y, c, ())

print("cutpoints", np.round(theta, 4))
for j, name in enumerate(NAMES):
    print(f"{name:>16s}  {beta[j]:8.4f}  ({se[c - 1 + j]:.4f})   z = {beta[j] / se[c - 1 + j]:6.3f}")
print(f"log-likelihood {ll_po:.4f} after {iters} Fisher scoring steps")
# <</fit>>

# ---- check against statsmodels OrderedModel ---------------------------------
from statsmodels.miscmodels.ordinal_model import OrderedModel

sm_fit = OrderedModel(y - 1, X, distr="logit").fit(method="bfgs", disp=0)
assert np.isclose(sm_fit.llf, ll_po, atol=1e-6)
assert np.allclose(np.asarray(sm_fit.params)[:p], beta, atol=1e-5)
sm_cuts = sm_fit.model.transform_threshold_params(np.asarray(sm_fit.params)[p:])[1:-1]
assert np.allclose(sm_cuts, theta, atol=1e-3)   # statsmodels stops its BFGS search earlier
assert np.all(np.diff(theta) > 0)                   # the fitted cutpoints are ordered

# the fitted probabilities are proper probabilities summing to one
P_po, _ = cum_probs(psi, X, c, ())
assert np.all(P_po > 0) and np.allclose(P_po.sum(axis=1), 1.0)

# <<collapse>>
# collapsing adjacent categories leaves the slopes estimating the same parameters
y3 = np.where(y <= 3, 1, np.where(y == 4, 2, 3))    # liberal / moderate / conservative
Y3 = np.zeros((n, 3))
Y3[np.arange(n), y3 - 1] = 1.0
psi3, J3, _ = cum_fit(X, Y3, 3)
se3 = np.sqrt(np.diag(np.linalg.inv(J3)))

y2 = np.where(y <= 4, 1, 2)                         # dichotomized at the middle category
Y2 = np.zeros((n, 2))
Y2[np.arange(n), y2 - 1] = 1.0
psi2, J2, _ = cum_fit(X, Y2, 2)
se2 = np.sqrt(np.diag(np.linalg.inv(J2)))

print("seven categories :", np.round(beta, 4), " se", np.round(se[c - 1:], 4))
print("three categories :", np.round(psi3[2:], 4), " se", np.round(se3[2:], 4))
print("two categories   :", np.round(psi2[1:], 4), " se", np.round(se2[1:], 4))
# <</collapse>>

assert np.all(se2[1:] > se3[2:]) and np.all(se3[2:] > se[c - 1:])   # collapsing costs information

# <<cuts>>
# the same slopes estimated one cut at a time: c - 1 separate binary logistic fits
Xc = sm.add_constant(X)
cut_beta = np.empty((c - 1, p))
cut_se = np.empty((c - 1, p))
for r in range(1, c):
    fit_r = sm.Logit((y <= r).astype(float), Xc).fit(disp=0)
    cut_beta[r - 1] = -np.asarray(fit_r.params)[1:]    # logit P(Y <= r) = theta_r - x^T beta_r
    cut_se[r - 1] = np.asarray(fit_r.bse)[1:]

for r in range(c - 1):
    print(f"cut {r + 1}: " + "  ".join(f"{b:7.4f} ({s:.4f})" for b, s in zip(cut_beta[r], cut_se[r])))
# <</cuts>>

assert np.allclose(cut_beta[3], psi2[1:], atol=1e-6)    # the cut-4 fit is the dichotomized fit

# <<score>>
# score test of proportional odds: all three slopes free to vary across cuts
free = (0, 1, 2)
psi_free = np.zeros((c - 1) + p + (c - 2) * len(free))
psi_free[:c - 1 + p] = psi                          # the proportional odds fit, with delta = 0
U, J_free = cum_score(psi_free, X, Y, c, free)
score_stat = float(U @ np.linalg.solve(J_free, U))
df = (c - 2) * len(free)

psi_np, _, _ = cum_fit(X, Y, c, free=free)
ll_np = cum_loglik(psi_np, X, Y, c, free)
lr_stat = 2 * (ll_np - ll_po)

print(f"score {score_stat:.3f} on {df} df, p = {stats.chi2.sf(score_stat, df):.3g}")
print(f"LR    {lr_stat:.3f} on {df} df, p = {stats.chi2.sf(lr_stat, df):.3g}")
# <</score>>

assert np.max(np.abs(U[:c - 1 + p])) < 1e-6         # the nuisance score vanishes at the restricted fit
assert abs(score_stat - lr_stat) < 0.05 * lr_stat   # the two statistics agree to within 5%

# one covariate at a time
one_score, one_lr = {}, {}
for j, name in enumerate(NAMES):
    psi_j = np.zeros((c - 1) + p + (c - 2))
    psi_j[:c - 1 + p] = psi
    Uj, Jj = cum_score(psi_j, X, Y, c, (j,))
    one_score[name] = float(Uj @ np.linalg.solve(Jj, Uj))
    psi_jf, _, _ = cum_fit(X, Y, c, free=(j,))
    one_lr[name] = 2 * (cum_loglik(psi_jf, X, Y, c, (j,)) - ll_po)
    assert abs(one_score[name] - one_lr[name]) < 0.1 * one_lr[name]

# <<partial>>
# partial proportional odds: education alone gets cut-specific slopes
psi_pp, J_pp, _ = cum_fit(X, Y, c, free=(1,))
ll_pp = cum_loglik(psi_pp, X, Y, c, (1,))
educ_slopes = np.concatenate([[psi_pp[c - 1 + 1]], psi_pp[c - 1 + 1] + psi_pp[c - 1 + p:]])

aic = lambda ll, k: -2 * ll + 2 * k
print("education slope by cut:", np.round(educ_slopes, 4))
print("AIC  proportional odds %.2f | partial %.2f | unrestricted %.2f"
      % (aic(ll_po, c - 1 + p), aic(ll_pp, len(psi_pp)), aic(ll_np, len(psi_np))))
# <</partial>>

P_pp, _ = cum_probs(psi_pp, X, c, (1,))
assert np.all(P_pp > 0)                             # no ordering violation in this sample
assert aic(ll_pp, len(psi_pp)) < aic(ll_po, c - 1 + p)

# the shape of the education slopes across the cuts, as the text describes it
educ_cuts = cut_beta[:, 1]
assert np.all(educ_cuts[:3] < -0.25)            # near -0.3 at the three lowest cuts
assert np.all(np.abs(educ_cuts[3:5]) < 0.10)    # near zero at cuts 4 and 5
assert educ_cuts[5] < -0.25                     # and back near -0.3 at cut 6

# ---- the concordance probability of exr-mlt-po-rank -------------------------


def concordance(theta, beta):
    """P(Y2 > Y1) + P(Y2 = Y1)/2 under proportional odds with a binary regressor."""
    def cat(xb):
        g = 1.0 / (1.0 + np.exp(-(np.asarray(theta, float) - xb)))
        return np.diff(np.concatenate([[0.0], g, [1.0]]))
    p1, p2 = cat(0.0), cat(beta)
    M = np.outer(p2, p1)                        # M[b, a] = P(Y2 = b) P(Y1 = a)
    return float(np.tril(M, -1).sum() + 0.5 * np.trace(M))


def latent_concordance(beta, lim=60.0, m=400001):
    """P(U2 > U1) = E{Lambda(beta + eps)} for standard logistic errors."""
    u = np.linspace(-lim, lim, m)
    dens = np.exp(-u) / (1.0 + np.exp(-u)) ** 2
    return float(np.trapezoid(dens / (1.0 + np.exp(-(beta + u))), u))


rank_latent = latent_concordance(1.0)
rank_near = concordance([-0.5, 0.5], 1.0)
rank_far = concordance([-3.0, 3.0], 1.0)

assert abs(latent_concordance(0.0) - 0.5) < 1e-9        # no effect, no concordance
assert rank_near > rank_far + 0.05                      # not a function of beta alone
assert (concordance([-0.5, 0.5], 2.0) > rank_near
        > concordance([-0.5, 0.5], 0.5) > concordance([-0.5, 0.5], 0.0))
assert abs(concordance([-0.5, 0.5], 0.0) - 0.5) < 1e-12  # symmetric when beta = 0

gen = Generated("ch36", "ordinal")
gen.int("n", n)
gen.int("c", c)
gen.int("iters", iters)
for r in range(c):
    gen.int(f"count{r + 1}", int(Y[:, r].sum()))
for r in range(c - 1):
    gen.num(f"theta{r + 1}", theta[r], 4)
for j, key in enumerate(["age", "educ", "income"]):
    gen.num(f"beta_{key}", beta[j], 4)
    gen.num(f"se_{key}", se[c - 1 + j], 4)
    gen.num(f"z_{key}", beta[j] / se[c - 1 + j], 3)
    gen.num(f"beta3_{key}", psi3[2 + j], 4)
    gen.num(f"se3_{key}", se3[2 + j], 4)
    gen.num(f"beta2_{key}", psi2[1 + j], 4)
    gen.num(f"se2_{key}", se2[1 + j], 4)
    gen.num(f"cutfit_{key}1", cut_beta[0, j], 4)
    gen.num(f"cutfit_{key}3", cut_beta[2, j], 4)
    gen.num(f"cutfit_{key}4", cut_beta[3, j], 4)
    gen.num(f"score_{key}", one_score[NAMES[j]], 3)
    gen.num(f"lr_{key}", one_lr[NAMES[j]], 3)
    gen.num(f"p_{key}", stats.chi2.sf(one_score[NAMES[j]], c - 2), 2, sci=True)
gen.num("or_age", np.exp(beta[0]), 4)
gen.num("or_educ", np.exp(beta[1]), 4)
gen.num("ll", ll_po, 4)
gen.num("aic_po", aic(ll_po, c - 1 + p), 2)
gen.num("aic_pp", aic(ll_pp, len(psi_pp)), 2)
gen.num("aic_np", aic(ll_np, len(psi_np)), 2)
gen.num("score", score_stat, 3)
gen.num("lr", lr_stat, 3)
gen.int("df", df)
gen.num("p_score", stats.chi2.sf(score_stat, df), 2, sci=True)
gen.num("p_lr", stats.chi2.sf(lr_stat, df), 2, sci=True)
gen.num("ll_pp", ll_pp, 4)
gen.num("lr_pp", 2 * (ll_pp - ll_po), 3)
gen.num("p_pp", stats.chi2.sf(2 * (ll_pp - ll_po), c - 2), 2, sci=True)
for r in range(c - 1):
    gen.num(f"educ_cut{r + 1}", educ_slopes[r], 4)
gen.num("rank_latent", rank_latent, 4)
gen.num("rank_near", rank_near, 4)
gen.num("rank_far", rank_far, 4)
gen.write()

# ---- figure: parallel cumulative logits and the category probabilities -------
use_book_style()
ages = np.linspace(2.0, 9.0, 200)
Xg = np.column_stack([ages, np.full_like(ages, 4.0), np.full_like(ages, 17.0)])
Gg, _ = cum_gamma(psi, Xg, c, ())
Pg, _ = cum_probs(psi, Xg, c, ())

fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6))
ax = axes[0]
for r in range(c - 1):
    ax.plot(10 * ages, np.log(Gg[:, r] / (1 - Gg[:, r])), color=COLORS["accent"], linewidth=1.0)
    ax.annotate(f"$r={r + 1}$", (10 * ages[-1], np.log(Gg[-1, r] / (1 - Gg[-1, r]))),
                textcoords="offset points", xytext=(3, -3), fontsize=7, color=COLORS["muted"])
ax.set_xlabel("age (years)")
ax.set_ylabel(r"logit $P(Y \leq r)$")
ax.set_title("(a) parallel cumulative logits")
ax.set_xlim(20, 97)

ax = axes[1]
shades = plt.cm.RdYlBu_r(np.linspace(0.08, 0.92, c))
ax.stackplot(10 * ages, Pg.T, colors=shades, edgecolor="white", linewidth=0.3)
mids = np.cumsum(Pg[-1]) - Pg[-1] / 2
for r in range(1, c):                               # keep thin bands' labels from colliding
    mids[r] = max(mids[r], mids[r - 1] + 0.075)
for r in range(c):
    ax.annotate(str(r + 1), (10 * ages[-1], mids[r]), textcoords="offset points",
                xytext=(3, -3), fontsize=7, color=COLORS["muted"])
ax.set_xlabel("age (years)")
ax.set_ylabel("fitted probability")
ax.set_title("(b) category probabilities")
ax.set_xlim(20, 93)
ax.set_ylim(0, 1)
fig.tight_layout()
fig.savefig(figure_path("ch36", "ordinal_fit"))

# ---- figure: the slopes cut by cut ------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(5.8, 2.2))
cuts = np.arange(1, c)
for j, ax in enumerate(axes):
    ax.axhline(beta[j], color=COLORS["second"], linewidth=1.0)
    ax.fill_between([0.5, c - 0.5],
                    beta[j] - 1.96 * se[c - 1 + j], beta[j] + 1.96 * se[c - 1 + j],
                    color=COLORS["second"], alpha=0.15, linewidth=0)
    ax.errorbar(cuts, cut_beta[:, j], yerr=1.96 * cut_se[:, j], fmt="o", markersize=3,
                color=COLORS["accent"], elinewidth=0.9, capsize=2)
    ax.axhline(0.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
    ax.set_xlabel("cut $r$")
    ax.set_xlim(0.5, c - 0.5)
    ax.set_title(NAMES[j])
axes[0].set_ylabel(r"$\hat\beta_r$")
fig.tight_layout()
fig.savefig(figure_path("ch36", "ordinal_cuts"))
