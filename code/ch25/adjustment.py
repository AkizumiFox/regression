"""Chapter 25, Section 6: regression adjustment in a completely randomized experiment.

A fixed synthetic finite population of 400 units with one covariate. The treatment
effect grows with the covariate (slope 2 under treatment, 0 under control). For several
treated fractions p, the three estimators (difference in means, Fisher's ANCOVA, Lin's
interacted regression) are computed over many randomizations and their variances are
compared with the exact variance of the fixed-slope estimator of the adjustment theorem of Section 25.6.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<population>>
import numpy as np

rng = np.random.default_rng(2506)
n = 400
x = rng.normal(size=n)
x = (x - x.mean()) / x.std(ddof=1)                    # centred, S_xx = 1
y0 = rng.normal(size=n)                               # control: no dependence on x
y1 = 1 + 2 * x + rng.normal(size=n)                   # treated: slope 2
tau = np.mean(y1 - y0)
beta1, beta0 = np.polyfit(x, y1, 1)[0], np.polyfit(x, y0, 1)[0]

def estimators(T):
    """Difference in means, ANCOVA and Lin's estimator for assignment rows T (reps x n)."""
    n1 = T.sum(axis=1, keepdims=True); n0 = n - n1
    yobs = np.where(T, y1, y0)
    xb1, xb0 = (T * x).sum(1, keepdims=True) / n1, (~T * x).sum(1, keepdims=True) / n0
    yb1, yb0 = (T * yobs).sum(1, keepdims=True) / n1, (~T * yobs).sum(1, keepdims=True) / n0
    dx = np.where(T, x - xb1, x - xb0)                # deviations from the arm means
    dy = np.where(T, yobs - yb1, yobs - yb0)
    b1 = (T * dx * dy).sum(1) / (T * dx ** 2).sum(1)  # within-arm slopes
    b0 = (~T * dx * dy).sum(1) / (~T * dx ** 2).sum(1)
    bp = (dx * dy).sum(1) / (dx ** 2).sum(1)          # pooled within-arm slope
    xb1, xb0, yb1, yb0 = xb1[:, 0], xb0[:, 0], yb1[:, 0], yb0[:, 0]
    unadj = yb1 - yb0
    ancova = unadj - bp * (xb1 - xb0)
    lin = (yb1 - b1 * xb1) - (yb0 - b0 * xb0)         # x is centred, so x-bar = 0
    return unadj, ancova, lin
# <</population>>

# <<quick>>
p, reps = 0.8, 2000
T = np.array([rng.permutation(n) < p * n for _ in range(reps)])
for name, est in zip(["difference in means", "ANCOVA", "Lin (interacted)"], estimators(T)):
    print(f"{name:20s} mean {est.mean():6.3f}  (tau = {tau:.3f})   n * variance {n * est.var():6.2f}")
# <</quick>>

# identities with least squares, for one assignment
t = T[0]
yobs = np.where(t, y1, y0)
one = np.ones(n)
u, a, l = (e[0] for e in estimators(T[:1]))
assert np.isclose(np.linalg.lstsq(np.column_stack([one, t, x]), yobs, rcond=None)[0][1], a)
assert np.isclose(np.linalg.lstsq(np.column_stack([one, t, x, t * x]), yobs, rcond=None)[0][1], l)


def exact_var(c, n1):
    """Exact randomization variance of the estimator with fixed common coefficient c."""
    n0 = n - n1
    u = y1 + (n1 / n0) * y0 - (n / n0) * c * x
    return (1 / n1 - 1 / n) * np.var(u, ddof=1)


def theory(p):
    """n times the variance of the fixed-slope estimators that the three estimators approach."""
    n1 = int(round(p * n)); n0 = n - n1
    c_opt = (n0 * beta1 + n1 * beta0) / n
    b_anc = (n1 * beta1 + n0 * beta0) / n
    return n * exact_var(0.0, n1), n * exact_var(b_anc, n1), n * exact_var(c_opt, n1)


# the variance identity V(c) = V(0) + n/(n0 n1) S_xx [(c - c_opt)^2 - c_opt^2]
for n1 in [80, 200, 320]:
    n0 = n - n1
    c_opt = (n0 * beta1 + n1 * beta0) / n
    for c in [-1.0, 0.3, 2.0]:
        rhs = exact_var(0, n1) + n / (n0 * n1) * ((c - c_opt) ** 2 - c_opt ** 2)
        assert np.isclose(exact_var(c, n1), rhs)

ps = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
sim = []
for pp in ps:
    T = np.array([rng.permutation(n) < pp * n for _ in range(20_000)])
    ests = estimators(T)
    sim.append([n * e.var() for e in ests])
    for e in ests:
        assert abs(e.mean() - tau) < 5 * e.std() / np.sqrt(20_000) + 0.01
sim = np.array(sim)
th = np.array([theory(pp) for pp in ps])
assert np.allclose(sim[:, 0], th[:, 0], rtol=0.05)      # exact for the difference in means
assert np.allclose(sim[:, 1:], th[:, 1:], rtol=0.08)     # asymptotic for the adjusted ones
# ANCOVA is worse than no adjustment exactly when b(b - 2 c_opt) > 0; here for p > 2/3
assert th[-1, 1] > th[-1, 0] and th[0, 1] < th[0, 0]
assert np.all(th[:, 2] <= th[:, 0] + 1e-12) and np.all(th[:, 2] <= th[:, 1] + 1e-12)

# ---- standard errors at p = 0.8 -------------------------------------------------
reps_se = 10_000
cover = np.zeros(3)
n1 = int(0.8 * n)
for r in range(reps_se):
    t = rng.permutation(n) < n1
    yobs = np.where(t, y1, y0)
    # difference in means with Neyman's variance estimate
    d = yobs[t].mean() - yobs[~t].mean()
    se_n = np.sqrt(yobs[t].var(ddof=1) / n1 + yobs[~t].var(ddof=1) / (n - n1))
    # ANCOVA with the classical least squares standard error
    Xa = np.column_stack([one, t, x])
    G = np.linalg.inv(Xa.T @ Xa)
    ba = G @ Xa.T @ yobs
    ra = yobs - Xa @ ba
    se_a = np.sqrt(ra @ ra / (n - 3) * G[1, 1])
    # Lin's estimator with the HC2 sandwich standard error
    Xl = np.column_stack([one, t, x, t * x])
    G = np.linalg.inv(Xl.T @ Xl)
    bl = G @ Xl.T @ yobs
    rl = yobs - Xl @ bl
    h = np.einsum("ij,jk,ik->i", Xl, G, Xl)
    meat = (Xl * (rl ** 2 / (1 - h))[:, None]).T @ Xl
    se_l = np.sqrt((G @ meat @ G)[1, 1])
    for j, (est, se) in enumerate([(d, se_n), (ba[1], se_a), (bl[1], se_l)]):
        cover[j] += abs(est - tau) <= 1.959964 * se
cover /= reps_se
assert cover[0] > 0.93 and cover[2] > 0.93
assert cover[1] < 0.95                                # classical ANCOVA SE undercovers

gen = Generated("ch25", "adjustment")
gen.int("n", n)
gen.num("tau", tau, 3)
gen.num("betaone", beta1, 3)
gen.num("betazero", beta0, 3)
for pp, s_, t_ in zip(ps, sim, th):
    tag = f"{int(round(pp * 10))}"
    for name, sv, tv in zip(["unadj", "ancova", "lin"], s_, t_):
        gen.num(f"sim{name}{tag}", sv, 2)
        gen.num(f"th{name}{tag}", tv, 2)
gen.num("coverunadj", cover[0], 3)
gen.num("coverancova", cover[1], 3)
gen.num("coverlin", cover[2], 3)
gen.int("repsse", reps_se)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.8, 2.8))
grid = np.linspace(0.1, 0.9, 81)
curves = np.array([theory(pp) for pp in grid])
for j, (label, col) in enumerate([("difference in means", COLORS["ink"]),
                                  ("ANCOVA (common slope)", COLORS["second"]),
                                  ("Lin (separate slopes)", COLORS["accent"])]):
    ax.plot(grid, curves[:, j], color=col, label=label)
    ax.plot(ps, sim[:, j], "o", color=col, markersize=3.5)
ax.axvline(2 / 3, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("fraction treated $p$")
ax.set_ylabel(r"$n\times$ variance")
ax.set_ylim(bottom=0)
ax.legend(frameon=False, loc="upper center")
fig.tight_layout()
fig.savefig(figure_path("ch25", "adjustment_variance"))
