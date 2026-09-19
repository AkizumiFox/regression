"""Chapter 13, Section 5: screening 100 regression coefficients with the
Benjamini-Hochberg procedure.

A fixed design with n = 250 observations, an intercept and 100 regressors with
independent standard normal entries (fixed seed). m1 of the 100 coefficients are
nonzero, each chosen so that its t statistic has noncentrality about 3. The t tests
share s and the design is not orthogonal, so the p-values are dependent.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

q = 0.05

# <<design>>
rng = np.random.default_rng(20260919)
n, m = 250, 100
X = np.column_stack([np.ones(n), rng.normal(size=(n, m))])
Qx, Rx = np.linalg.qr(X)
nu = n - m - 1                                        # 149 error degrees of freedom
c_diag = np.diag(np.linalg.inv(Rx.T @ Rx))[1:]        # [(X'X)^{-1}]_jj for the slopes

def t_statistics(Y):
    """t statistics of the 100 slopes for each column of Y."""
    coef = np.linalg.solve(Rx, Qx.T @ Y)
    resid = Y - X @ coef
    s2 = np.sum(resid ** 2, axis=0) / nu
    return coef[1:] / np.sqrt(np.outer(c_diag, s2))
# <</design>>

# <<procedures>>
def bh(p, q):
    """Benjamini-Hochberg step-up at level q; returns a boolean rejection vector."""
    m = len(p)
    order = np.argsort(p)
    below = np.nonzero(p[order] <= q * np.arange(1, m + 1) / m)[0]
    reject = np.zeros(m, dtype=bool)
    if below.size:
        reject[order[: below[-1] + 1]] = True
    return reject

def holm(p, alpha):
    m = len(p)
    order = np.argsort(p)
    ok = np.cumprod(p[order] <= alpha / (m - np.arange(m))).astype(bool)
    reject = np.zeros(m, dtype=bool)
    reject[order] = ok
    return reject

procedures = {
    "unadjusted": lambda p: p <= q,
    "Bonferroni": lambda p: p <= q / m,
    "Holm": lambda p: holm(p, q),
    "BH": lambda p: bh(p, q),
}
# <</procedures>>

# <<one>>
m1 = 10                                               # the first 10 slopes are nonzero
beta = np.zeros(m + 1)
beta[1:m1 + 1] = 3.0 * np.sqrt(c_diag[:m1])           # noncentrality 3 when sigma = 1
y = X @ beta + rng.normal(size=n)
p = 2 * stats.t.sf(np.abs(t_statistics(y[:, None])[:, 0]), nu)
for name, rule in procedures.items():
    rej = rule(p)
    print(f"{name:10s} rejects {rej.sum():3d}, of which {rej[m1:].sum()} are false discoveries")
# <</one>>

gen = Generated("ch13", "bh_screening", prefix="bh")
gen.int("n", n)
gen.int("m", m)
gen.int("nu", nu)
gen.int("mone", m1)
one = {name: rule(p) for name, rule in procedures.items()}
for name, key in (("unadjusted", "t"), ("Bonferroni", "bon"), ("Holm", "holm"), ("BH", "bh")):
    gen.int(f"rej{key}", one[name].sum())
    gen.int(f"false{key}", one[name][m1:].sum())
assert one["Bonferroni"].sum() <= one["Holm"].sum() <= one["BH"].sum() <= one["unadjusted"].sum()
order = np.argsort(p)
k_bh = one["BH"].sum()
gen.num("pkbh", p[order][k_bh - 1], 4)
gen.num("cutkbh", q * k_bh / m, 4)
gen.num("bonfcut", q / m, 4)
p_one, order_one = p.copy(), order.copy()

# ---- simulation of error rates and power -----------------------------------------
reps = 2000
m1s = [0, 5, 10, 25, 50]
summary = {name: {"FDR": [], "FWER": [], "power": []} for name in procedures}
for m1_ in m1s:
    b = np.zeros(m + 1)
    b[1:m1_ + 1] = 3.0 * np.sqrt(c_diag[:m1_])
    Y = (X @ b)[:, None] + rng.normal(size=(n, reps))
    P = 2 * stats.t.sf(np.abs(t_statistics(Y)), nu)
    for name, rule in procedures.items():
        fdp, anyfalse, pw = [], [], []
        for r in range(reps):
            rej = rule(P[:, r])
            V, R = rej[m1_:].sum(), rej.sum()
            fdp.append(V / max(R, 1))
            anyfalse.append(V > 0)
            pw.append(rej[:m1_].mean() if m1_ else np.nan)
        summary[name]["FDR"].append(np.mean(fdp))
        summary[name]["FWER"].append(np.mean(anyfalse))
        summary[name]["power"].append(np.nanmean(pw) if m1_ else np.nan)
        se = np.std(fdp) / np.sqrt(reps)
        if name == "BH":
            assert np.mean(fdp) <= q * (m - m1_) / m + 4 * se + 1e-3
        if name in ("Bonferroni", "Holm"):
            assert np.mean(anyfalse) <= q + 4 * np.sqrt(q * (1 - q) / reps)
    if m1_ == 0:                                       # complete null: FDR = FWER
        assert np.isclose(summary["BH"]["FDR"][-1], summary["BH"]["FWER"][-1])

for i, m1_ in enumerate(m1s):
    for name, key in (("unadjusted", "t"), ("Bonferroni", "bon"), ("Holm", "holm"), ("BH", "bh")):
        gen.num(f"fdr{key}{m1_}", summary[name]["FDR"][i], 3)
        gen.num(f"fwer{key}{m1_}", summary[name]["FWER"][i], 3)
        if m1_:
            gen.num(f"pow{key}{m1_}", summary[name]["power"][i], 3)
gen.int("reps", reps)

# ---- exactness under independence: FDR = q m0 / m ---------------------------------
rng2 = np.random.default_rng(7)
reps2, m0 = 20_000, 90
Zs = rng2.normal(size=(reps2, m))
Zs[:, m0:] += 3.0
Ps = 2 * stats.norm.sf(np.abs(Zs))
fdp = np.empty(reps2)
for r in range(reps2):
    rej = bh(Ps[r], q)
    fdp[r] = rej[:m0].sum() / max(rej.sum(), 1)
fdr_ind = fdp.mean()
assert abs(fdr_ind - q * m0 / m) < 4 * fdp.std() / np.sqrt(reps2)
gen.num("fdrind", fdr_ind, 4)
gen.num("fdrindtheory", q * m0 / m, 4)
gen.int("repsind", reps2)
harmonic = np.sum(1 / np.arange(1, m + 1))          # Benjamini-Yekutieli factor
assert 5.18 < harmonic < 5.19
gen.num("harmonic", harmonic, 3)
gen.num("byq", q / harmonic, 4)
gen.write()

# ---- figures -----------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.7))
K = 20
ranks = np.arange(1, K + 1)
ps = p_one[order_one][:K]
false = order_one[:K] >= m1
ax.plot(ranks, q * ranks / m, color=COLORS["accent"], label=r"BH line $kq/m$")
ax.axhline(q / m, color=COLORS["second"], linestyle="--", label=r"Bonferroni $q/m$")
ax.plot(ranks[~false], ps[~false], "o", color=COLORS["ink"], markersize=3.5,
        label="nonzero coefficient")
ax.plot(ranks[false], ps[false], "o", mfc="white", color=COLORS["ink"], markersize=3.5,
        label="zero coefficient")
ax.axvline(k_bh + 0.5, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.set_xlabel("rank $k$ of the $p$-value")
ax.set_yscale("log")
ax.set_ylabel("ordered $p$-value (log scale)")
ax.legend(frameon=False, loc="lower right")
fig.savefig(figure_path("ch13", "bh_stepup"))

fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.7))
styles = {"unadjusted": (COLORS["muted"], ":", "o"), "Bonferroni": (COLORS["second"], "-", "s"),
          "Holm": (COLORS["third"], "--", "^"), "BH": (COLORS["accent"], "-", "D")}
for name, st in styles.items():
    col, ls, mk = st
    axes[0].plot(m1s, summary[name]["FDR"], color=col, linestyle=ls, marker=mk, markersize=3,
                 label=name)
    axes[1].plot(m1s[1:], summary[name]["power"][1:], color=col, linestyle=ls, marker=mk,
                 markersize=3, label=name)
axes[0].axhline(q, color=COLORS["grid"], linewidth=0.8, zorder=0)
axes[0].set_yscale("log")
axes[0].set_xlabel("nonzero coefficients $m_1$")
axes[0].set_ylabel("false discovery rate")
axes[0].set_title("(a) FDR (log scale)")
axes[1].set_xlabel("nonzero coefficients $m_1$")
axes[1].set_ylabel("average power")
axes[1].set_title("(b) power")
fig.tight_layout(rect=(0, 0.08, 1, 1))
handles, labs = axes[1].get_legend_handles_labels()
fig.legend(handles, labs, frameon=False, ncol=4, loc="lower center")
fig.savefig(figure_path("ch13", "fdr_simulation"))
