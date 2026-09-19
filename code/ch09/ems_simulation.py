"""Chapter 9, Section 6: expected mean squares, checked by simulation.

The design is the unbalanced 3 x 3 party-by-education layout of the 1996 American
National Election Study (statsmodels.datasets.anes96, public domain). The true cell
means are set to the observed ones and sigma = 1.1. Each simulated data set is split
into the sequential (Type I) pieces mean, party, education | party, interaction and
residual, plus the Type III party sum of squares, and the average of every mean square
is compared with sigma^2 + ||P mu||^2 / rank(P). Normal and skewed errors are used.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

d = sm.datasets.anes96.load_pandas().data
party = pd.cut(d["PID"], [-1, 1.5, 4.5, 6], labels=["Dem", "Ind", "Rep"])
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"])
cell = party.cat.codes.to_numpy() * 3 + educ.cat.codes.to_numpy()
n = len(cell)


def basis(*blocks):
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    return U[:, s > 1e-10 * s[0]]


# <<pieces>>
one = np.ones((n, 1))
A = np.eye(3)[party.cat.codes.to_numpy()]
B = np.eye(3)[educ.cat.codes.to_numpy()]
Cell = np.eye(9)[cell]
U0, UA, UAB, UC = basis(one), basis(one, A), basis(one, A, B), basis(Cell)


def complement(U_big, U_small):
    """Orthonormal basis of C(U_small)-perp within C(U_big)."""
    R = U_big - U_small @ (U_small.T @ U_big)
    return basis(R)


pieces = {"mean": U0, "party": complement(UA, U0),
          "educ | party": complement(UAB, UA), "interaction": complement(UC, UAB)}
Ub = basis(one, B, np.column_stack([(A[:, [i]] - A[:, [0]]) * (B[:, [j]] - B[:, [0]])
                                     for i in (1, 2) for j in (1, 2)]))
pieces["party (Type III)"] = complement(UC, Ub)       # unweighted-means hypothesis

sigma = 1.1
mu = d.groupby(cell)["selfLR"].mean().to_numpy()[cell]  # true means: observed cell means
for name, U in pieces.items():
    r = U.shape[1]
    gamma = np.sum((U.T @ mu) ** 2) / sigma ** 2        # noncentrality ||P mu||^2 / sigma^2
    print(f"{name:17s} rank {r}  gamma = {gamma:9.2f}"
          f"  E(MS) = {sigma ** 2 * (1 + gamma / r):9.3f}")
# <</pieces>>

y_obs = d["selfLR"].to_numpy()
assert np.isclose(np.sum((pieces["party (Type III)"].T @ y_obs) ** 2),
                  sm.stats.anova_lm(sm.OLS.from_formula(
                      "y ~ C(p, Sum) * C(e, Sum)",
                      pd.DataFrame({"y": y_obs, "p": party, "e": educ})).fit(), typ=3).iloc[1, 0])

# ---- simulation -------------------------------------------------------------
# <<simulate>>
rng = np.random.default_rng(99)
reps, chunk = 40_000, 4_000
dfe = n - 9
sums = {"normal": [], "skewed": []}
for kind in sums:
    out = []
    for _ in range(reps // chunk):
        if kind == "normal":
            E = sigma * rng.normal(size=(n, chunk))
        else:                                             # centred exponential, variance sigma^2
            E = sigma * (rng.exponential(size=(n, chunk)) - 1.0)
        Y = mu[:, None] + E
        ms = {name: np.sum((U.T @ Y) ** 2, axis=0) / U.shape[1] for name, U in pieces.items()}
        ms["residual"] = (np.sum(Y ** 2, axis=0) - np.sum((UC.T @ Y) ** 2, axis=0)) / dfe
        out.append(pd.DataFrame(ms))
    sums[kind] = pd.concat(out, ignore_index=True)
print(pd.DataFrame({k: v.mean() for k, v in sums.items()}).round(3))
# <</simulate>>

gen = Generated("ch09", "ems_simulation", prefix="ems")
gen.num("sigma2", sigma ** 2, 2)
gen.int("reps", reps)
keys = {"mean": "mean", "party": "party", "educ | party": "educ", "interaction": "inter",
        "party (Type III)": "party3", "residual": "res"}
theory = {}
for name, key in keys.items():
    if name == "residual":
        r, gamma = dfe, 0.0
    else:
        U = pieces[name]
        r, gamma = U.shape[1], np.sum((U.T @ mu) ** 2) / sigma ** 2
    ems = sigma ** 2 * (1 + gamma / r)
    var_normal = 2 * sigma ** 4 * (r + 2 * gamma) / r ** 2  # Var(MS) = sigma^4 (2r + 4 gamma)/r^2
    theory[name] = (r, gamma, ems)
    for kind in ("normal", "skewed"):
        sim = sums[kind][name]
        se = sim.std() / np.sqrt(reps)
        assert abs(sim.mean() - ems) < 4.5 * se, (name, kind, sim.mean(), ems, se)
        gen.num(f"{key}:{kind}", sim.mean(), 3)
    assert abs(sums["normal"][name].var() / var_normal - 1) < 0.05
    gen.int(f"{key}:r", r)
    gen.num(f"{key}:gamma", gamma, 2)
    gen.num(f"{key}:ems", ems, 3)
    gen.num(f"{key}:var", var_normal, 4 if var_normal < 1 else 2)
    gen.num(f"{key}:simvar", sums["normal"][name].var(), 4 if var_normal < 1 else 2)
# independence under normality, but not under skewed errors
corr_normal = np.corrcoef(sums["normal"]["interaction"], sums["normal"]["residual"])[0, 1]
corr_skewed = np.corrcoef(sums["skewed"]["interaction"], sums["skewed"]["residual"])[0, 1]
assert abs(corr_normal) < 0.03
gen.num("corr:normal", corr_normal, 3)
gen.num("corr:skewed", corr_skewed, 3)
gen.num("var:skewed:res", sums["skewed"]["residual"].var(), 4)

# skewed-error theory for the quadratic forms MS_A = Y'AY (interaction) and MS_B = Y'BY
# (residual): centred exponential errors have mu3 = 2 sigma^3 and mu4 = 9 sigma^4, and
# Cov(Y'AY, Y'BY) = (mu4 - 3 sigma^4) sum a_ii b_ii + 2 sigma^4 tr(AB)
#                   + 2 mu3 sum [(A mu)_i b_ii + (B mu)_i a_ii] + 4 sigma^2 mu'ABmu,
# where here AB = 0 and B mu = 0.
mu3, mu4 = 2 * sigma ** 3, 9 * sigma ** 4
Uint = pieces["interaction"]
r_int = Uint.shape[1]
a_diag = np.sum(Uint ** 2, axis=1) / r_int
b_diag = (1 - np.sum(UC ** 2, axis=1)) / dfe
A_mu = Uint @ (Uint.T @ mu) / r_int
var_a = ((mu4 - 3 * sigma ** 4) * np.sum(a_diag ** 2) + 2 * sigma ** 4 / r_int
         + 4 * sigma ** 2 * np.sum(A_mu ** 2) + 4 * mu3 * np.sum(A_mu * a_diag))
var_b = (mu4 - 3 * sigma ** 4) * np.sum(b_diag ** 2) + 2 * sigma ** 4 / dfe
cov_ab = (mu4 - 3 * sigma ** 4) * np.sum(a_diag * b_diag) + 2 * mu3 * np.sum(A_mu * b_diag)
corr_theory = cov_ab / np.sqrt(var_a * var_b)
var_res_skewed = sums["skewed"]["residual"].var()
ratio = var_res_skewed / (2 * sigma ** 4 / dfe)
assert 3.5 < ratio < 4.5                                  # about 8 sigma^4/(n-r) against 2 sigma^4/(n-r)
assert abs(var_res_skewed / var_b - 1) < 0.1
assert 0.015 < corr_theory < 0.035 and abs(corr_skewed - corr_theory) < 4 / np.sqrt(reps)
assert abs(corr_skewed) > 4 / np.sqrt(reps)               # the dependence is detectable
gen.num("var:skewed:res:theory", var_b, 4)
gen.num("var:skewed:ratio", ratio, 2)
gen.num("corr:skewed:theory", corr_theory, 3)

# F ratio for the interaction: noncentral F(4, n - 9, gamma)
F = sums["normal"]["interaction"] / sums["normal"]["residual"]
r_int, g_int, _ = theory["interaction"]
crit = stats.f.ppf(0.95, r_int, dfe)
power_theory = stats.ncf.sf(crit, r_int, dfe, g_int)
power_sim = np.mean(F > crit)
assert abs(power_sim - power_theory) < 4 * np.sqrt(power_theory * (1 - power_theory) / reps)
gen.num("crit", crit, 3)
gen.num("power:theory", power_theory, 4)
gen.num("power:sim", power_sim, 4)
gen.num("Fmean", F.mean(), 3)
gen.num("Fmean:theory", dfe * (r_int + g_int) / (r_int * (dfe - 2)), 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.4))
grid = np.linspace(0.01, 25, 500)
ax.hist(F, bins=np.linspace(0, 25, 101), density=True, color=COLORS["accent"], alpha=0.3,
        label="simulated $F$ for interaction")
ax.plot(grid, stats.ncf.pdf(grid, r_int, dfe, g_int), color=COLORS["accent"],
        label=rf"$F(4,{dfe},{g_int:.1f})$")
ax.plot(grid, stats.f.pdf(grid, r_int, dfe), color=COLORS["second"], label=rf"central $F(4,{dfe})$")
ax.axvline(crit, color=COLORS["muted"], lw=0.8, ls="--")
ax.text(crit + 0.3, 0.62, "5% point", color=COLORS["muted"], fontsize=7)
ax.set_xlabel("F")
ax.set_ylabel("density")
ax.set_ylim(0, 0.8)
ax.legend(frameon=False)
fig.savefig(figure_path("ch09", "ems_f_ratio"))
