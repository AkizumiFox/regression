"""Chapter 41, Section 2: when dropping the incomplete records is safe, and what
the quick repairs cost.

Part one simulates four missingness mechanisms for the response of a linear model
and records the sampling behaviour of the complete-case slope. Part two measures
the variance deflation of mean imputation and of regression imputation without
noise. Part three builds a twelve-row data set on which pairwise deletion returns
a correlation matrix that is not nonnegative definite.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<scenarios>>
rng = np.random.default_rng(4101)
n, reps = 300, 2000
beta = np.array([1.0, 0.8, 0.9])                      # intercept, coefficient of x, of w


def one_run(mechanism, rng):
    x, w = rng.normal(size=n), rng.normal(size=n)
    y = beta[0] + beta[1] * x + beta[2] * w + rng.normal(size=n)
    if mechanism == "MCAR":
        eta = np.full(n, -0.6)
    elif mechanism == "MAR on x":
        eta = -0.6 + 1.6 * x                          # depends on a regressor only
    elif mechanism == "MAR on x and w":
        eta = -0.6 + 1.2 * x + 1.5 * w                # depends on observed variables only
    else:                                             # MNAR
        eta = -0.6 + 1.2 * (y - beta[0])              # depends on the value itself
    missing = rng.uniform(size=n) < 1 / (1 + np.exp(-eta))
    obs = ~missing
    short = np.column_stack([np.ones(obs.sum()), x[obs]])            # y on x
    long = np.column_stack([np.ones(obs.sum()), x[obs], w[obs]])     # y on x and w
    b_short, *_ = np.linalg.lstsq(short, y[obs], rcond=None)
    b_long, *_ = np.linalg.lstsq(long, y[obs], rcond=None)
    return b_short[1], b_long[1], obs.sum()


names = ["MCAR", "MAR on x", "MAR on x and w", "MNAR"]
out = {m: np.array([one_run(m, rng) for _ in range(reps)]) for m in names}
for m in names:
    s, l, k = out[m].T
    print(f"{m:16s} complete cases {k.mean():5.0f}   y ~ x: {s.mean():.4f}"
          f"   y ~ x + w: {l.mean():.4f}")
# <</scenarios>>

bias = {m: (out[m][:, 0].mean() - beta[1], out[m][:, 1].mean() - beta[1]) for m in names}
mcse = {m: out[m][:, 0].std(ddof=1) / np.sqrt(reps) for m in names}
for m in ["MCAR", "MAR on x"]:                       # selection on the regressors is harmless
    assert abs(bias[m][0]) < 4 * mcse[m]
assert bias["MAR on x and w"][0] < -8 * mcse["MAR on x and w"]      # the short fit is biased
assert abs(bias["MAR on x and w"][1]) < 4 * out["MAR on x and w"][:, 1].std(ddof=1) / np.sqrt(reps)
assert bias["MNAR"][0] < -8 * mcse["MNAR"]
assert bias["MNAR"][1] < -8 * out["MNAR"][:, 1].std(ddof=1) / np.sqrt(reps)

# <<single>>
# mean imputation: the filled-in column has too little spread
rng2 = np.random.default_rng(4102)
m_frac, trials = 0.3, 4000
var_imputed, var_cc = [], []
for _ in range(trials):
    z = rng2.normal(0.0, 1.0, n)
    keep = rng2.uniform(size=n) > m_frac                       # MCAR
    filled = np.where(keep, z, z[keep].mean())                 # every gap gets the same value
    var_imputed.append(filled.var(ddof=1))
    var_cc.append(z[keep].var(ddof=1))
print(f"mean imputation: E(s^2) = {np.mean(var_imputed):.4f};"
      f" complete cases {np.mean(var_cc):.4f}")
# <</single>>

k_obs = int(round(n * (1 - m_frac)))
predicted = (k_obs - 1) / (n - 1)                     # the deflation factor of the section
assert abs(np.mean(var_imputed) - predicted) < 0.02
assert abs(np.mean(var_cc) - 1.0) < 0.02

# regression imputation without noise: the imputed points sit exactly on the line
rng3 = np.random.default_rng(4103)
rho_hat, sd_hat = [], []
for _ in range(trials):
    u = rng3.normal(size=n)
    v = 0.6 * u + np.sqrt(1 - 0.36) * rng3.normal(size=n)       # true correlation 0.6
    keep = rng3.uniform(size=n) > m_frac
    b0, b1 = np.polyfit(u[keep], v[keep], 1)[::-1]    # intercept, slope
    v_fill = np.where(keep, v, b0 + b1 * u)           # the gaps land exactly on the line
    rho_hat.append(np.corrcoef(u, v_fill)[0, 1])
    sd_hat.append(v_fill.std(ddof=1))
print(f"regression imputation: mean correlation {np.mean(rho_hat):.4f} (true 0.6),"
      f" mean sd {np.mean(sd_hat):.4f} (true 1)")
assert np.mean(rho_hat) > 0.66 and np.mean(sd_hat) < 0.96

# ---- pairwise deletion can break nonnegative definiteness --------------------
# <<pairwise>>
nan = np.nan
u = np.array([-1.5, -0.5, 0.5, 1.5])
d = np.array([0.5, -1.5, 1.5, -0.5])                  # orthogonal to u, same length
pair = np.column_stack([u, 0.9 * u + np.sqrt(1 - 0.81) * d])             # correlation +0.9
flip = np.column_stack([u, -(0.9 * u + np.sqrt(1 - 0.81) * d)])          # correlation -0.9
D = np.vstack([
    np.column_stack([pair, np.full(4, nan)]),         # pattern A: z3 absent
    np.column_stack([pair[:, 0], np.full(4, nan), pair[:, 1]]),          # pattern B: z2 absent
    np.column_stack([np.full(4, nan), flip]),         # pattern C: z1 absent
])
R = np.ones((3, 3))
for i in range(3):
    for j in range(i + 1, 3):
        both = ~np.isnan(D[:, i]) & ~np.isnan(D[:, j])
        R[i, j] = R[j, i] = np.corrcoef(D[both, i], D[both, j])[0, 1]
print("available-case correlation matrix\n", R)
print("eigenvalues", np.round(np.linalg.eigvalsh(R), 3))
# <</pairwise>>

assert np.linalg.eigvalsh(R).min() < -0.5

gen = Generated("ch41", "complete_case")
gen.int("n", n)
gen.int("reps", reps)
for key, m in zip(["mcar", "marx", "marxw", "mnar"], names):
    gen.num(f"short_{key}", out[m][:, 0].mean(), 3)
    gen.num(f"long_{key}", out[m][:, 1].mean(), 3)
    gen.int(f"cc_{key}", int(round(out[m][:, 2].mean())))
gen.num("var_imputed", float(np.mean(var_imputed)), 3)
gen.num("var_cc", float(np.mean(var_cc)), 3)
gen.num("deflation", predicted, 3)
gen.num("rho_fill", float(np.mean(rho_hat)), 3)
gen.num("sd_fill", float(np.mean(sd_hat)), 3)
gen.num("eig_min", float(np.linalg.eigvalsh(R).min()), 3)
gen.write()

# ---- figure: the four mechanisms --------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(5.6, 2.4))
positions, centres, data, colours = [], [], [], []
for k, m in enumerate(names):
    positions += [2 * k + 0.75, 2 * k + 1.35]
    centres.append(2 * k + 1.05)
    data += [out[m][:, 0], out[m][:, 1]]
    colours += [COLORS["accent"], COLORS["third"]]
bp = ax.boxplot(data, positions=positions, widths=0.5, showfliers=False,
                patch_artist=True, medianprops=dict(color=COLORS["ink"], linewidth=0.8))
for patch, c in zip(bp["boxes"], colours):
    patch.set_facecolor(c)
    patch.set_alpha(0.45)
    patch.set_linewidth(0.5)
ax.axhline(beta[1], color=COLORS["second"], linewidth=0.9, zorder=0)
ax.set_xticks(centres)
ax.set_xticklabels(names)
ax.set_ylabel(r"complete-case estimate of $\beta_1$")
handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, alpha=0.45, linewidth=0.5,
                         edgecolor=COLORS["ink"]) for c in [COLORS["accent"], COLORS["third"]]]
ax.legend(handles, [r"fit $y\sim x$", r"fit $y\sim x+w$"], frameon=False,
          fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch41", "complete_case_bias"))
