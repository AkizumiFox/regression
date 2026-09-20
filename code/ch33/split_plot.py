"""Chapter 33, Section 1: a balanced split-plot experiment.

Four blocks (fields), three irrigation regimes assigned to the whole plots of each
block, four barley varieties assigned to the four subplots of each whole plot.
The yields are simulated from the split-plot model with whole-plot standard
deviation 0.35 and subplot standard deviation 0.25 tonnes per hectare.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<simulate>>
r, a, b = 4, 3, 4                                   # blocks, irrigation regimes, varieties
n = r * a * b
block_eff = np.array([-0.30, 0.10, 0.05, 0.15])     # field-to-field fertility
irrig_eff = np.array([0.00, 0.45, 0.70])            # flood, sprinkler, drip
variety_eff = np.array([-0.30, -0.10, 0.15, 0.25])
inter = np.array([[0.10, -0.05, -0.10, 0.05],       # irrigation x variety
                  [-0.05, 0.10, 0.00, -0.05],
                  [-0.05, -0.05, 0.10, 0.00]])
sigma_w, sigma_s = 0.35, 0.25                       # whole-plot and subplot errors

rng = np.random.default_rng(3301)
eta = rng.normal(0, sigma_w, (r, a))                # one error per whole plot
e = rng.normal(0, sigma_s, (r, a, b))               # one error per subplot
mean = (5.0 + block_eff[:, None, None] + irrig_eff[None, :, None]
        + variety_eff[None, None, :] + inter[None, :, :])
Y = np.round(mean + eta[:, :, None] + e, 2)         # yields, t/ha
y = Y.reshape(-1)                                   # stacked block, then regime, then variety
# <</simulate>>


def Jbar(k):
    return np.full((k, k), 1.0 / k)


def Cen(k):
    return np.eye(k) - Jbar(k)


def kron3(E1, E2, E3):
    return np.kron(np.kron(E1, E2), E3)


# <<anova>>
P = {"blocks": kron3(Cen(r), Jbar(a), Jbar(b)),
     "irrigation": kron3(Jbar(r), Cen(a), Jbar(b)),
     "whole-plot error": kron3(Cen(r), Cen(a), Jbar(b)),
     "varieties": kron3(Jbar(r), Jbar(a), Cen(b)),
     "irrigation x variety": kron3(Jbar(r), Cen(a), Cen(b)),
     "subplot error": kron3(Cen(r), np.eye(a), Cen(b))}
df = {k: int(round(np.trace(M))) for k, M in P.items()}
SS = {k: y @ M @ y for k, M in P.items()}
MS = {k: SS[k] / df[k] for k in P}

F_irrigation = MS["irrigation"] / MS["whole-plot error"]
F_variety = MS["varieties"] / MS["subplot error"]
F_inter = MS["irrigation x variety"] / MS["subplot error"]
for name, F, d in [("irrigation", F_irrigation, df["whole-plot error"]),
                   ("varieties", F_variety, df["subplot error"]),
                   ("irrigation x variety", F_inter, df["subplot error"])]:
    p = stats.f.sf(F, df[name], d)
    print(f"{name:22s} F({df[name]},{d}) = {F:7.3f}   p = {p:.4f}")
# <</anova>>

# the six projections plus the grand mean project onto orthogonal subspaces summing to I
P0 = kron3(Jbar(r), Jbar(a), Jbar(b))
allP = [P0] + list(P.values())
assert np.allclose(sum(allP), np.eye(n))
for i, Mi in enumerate(allP):
    assert np.allclose(Mi @ Mi, Mi) and np.allclose(Mi, Mi.T)
    for Mj in allP[i + 1:]:
        assert np.allclose(Mi @ Mj, 0)
assert sum(df.values()) + 1 == n
assert df == {"blocks": r - 1, "irrigation": a - 1, "whole-plot error": (r - 1) * (a - 1),
              "varieties": b - 1, "irrigation x variety": (a - 1) * (b - 1),
              "subplot error": (r - 1) * a * (b - 1)}

# the covariance matrix and the two strata ------------------------------------
P1 = kron3(np.eye(r), np.eye(a), Jbar(b))           # projection onto whole-plot means
V = sigma_s**2 * np.eye(n) + b * sigma_w**2 * P1
Zw = np.kron(np.eye(r * a), np.ones((b, 1)))        # whole-plot indicator matrix
assert np.allclose(V, sigma_s**2 * np.eye(n) + sigma_w**2 * Zw @ Zw.T)
lam_w, lam_s = sigma_s**2 + b * sigma_w**2, sigma_s**2
for name in ["blocks", "irrigation", "whole-plot error"]:
    assert np.allclose(V @ P[name], lam_w * P[name])          # whole-plot stratum
for name in ["varieties", "irrigation x variety", "subplot error"]:
    assert np.allclose(V @ P[name], lam_s * P[name])          # subplot stratum

# ordinary least squares is BLUE here: C(VX) is contained in C(X)  (Kruskal)
M = P0 + P["blocks"] + P["irrigation"] + P["varieties"] + P["irrigation x variety"]
assert np.allclose(V @ M, M @ V)                    # so V M = M V M has columns in C(M)
assert np.allclose(M @ V @ M, V @ M)

# expected mean squares, checked against a simulation --------------------------
mu = mean[:, :, :].reshape(-1)
ems = {k: (np.trace(P[k] @ V) + mu @ P[k] @ mu) / df[k] for k in P}
sim = {k: 0.0 for k in P}
B = 20000
rs = np.random.default_rng(77)
L = np.linalg.cholesky(V)
for _ in range(B):
    ysim = mu + L @ rs.standard_normal(n)
    for k in P:
        sim[k] += (ysim @ P[k] @ ysim) / df[k] / B
for k in P:
    assert abs(sim[k] - ems[k]) < 0.05 * ems[k], (k, sim[k], ems[k])
assert np.isclose(ems["whole-plot error"], lam_w)
assert np.isclose(ems["subplot error"], lam_s)

# standard errors of the two kinds of comparison ------------------------------
irr_means = Y.mean(axis=(0, 2))
var_means = Y.mean(axis=(0, 1))
se_irrig = np.sqrt(2 * MS["whole-plot error"] / (r * b))
se_variety = np.sqrt(2 * MS["subplot error"] / (r * a))
MS_pooled = ((SS["whole-plot error"] + SS["subplot error"])
             / (df["whole-plot error"] + df["subplot error"]))
se_pooled_irrig = np.sqrt(2 * MS_pooled / (r * b))
se_pooled_variety = np.sqrt(2 * MS_pooled / (r * a))
print("irrigation means", irr_means.round(3), "se of a difference", round(se_irrig, 4))
print("variety means   ", var_means.round(3), "se of a difference", round(se_variety, 4))

# a simple effect of A at one level of C belongs to neither stratum -------------
c_simple = np.zeros(a)
c_simple[2], c_simple[0] = 1.0, -1.0                # drip minus flood
k_fixed = 2                                         # at one variety only
rho_simple = np.kron(np.kron(np.ones(r) / r, c_simple), np.eye(b)[k_fixed])
var_simple = rho_simple @ V @ rho_simple
assert np.isclose(var_simple, 2 * (sigma_w**2 + sigma_s**2) / r)
assert np.isclose(var_simple, 2 * (lam_w + (b - 1) * lam_s) / (r * b))
assert np.isclose(rho_simple @ Y.reshape(-1),
                  Y[:, 2, k_fixed].mean() - Y[:, 0, k_fixed].mean())
var_simple_hat = 2 * (MS["whole-plot error"] + (b - 1) * MS["subplot error"]) / (r * b)
se_simple = np.sqrt(var_simple_hat)
c1, c2 = 2 / (r * b), 2 * (b - 1) / (r * b)         # weights of the two mean squares
nu_simple = var_simple_hat**2 / (
    (c1 * MS["whole-plot error"])**2 / df["whole-plot error"]
    + (c2 * MS["subplot error"])**2 / df["subplot error"])
assert df["whole-plot error"] < nu_simple < df["subplot error"]
print("simple effect: se", round(se_simple, 4), "on", round(nu_simple, 2), "df")

gen = Generated("ch33", "split_plot")
gen.int("n", n)
gen.int("r", r)
gen.int("a", a)
gen.int("b", b)
for k in P:
    key = k.replace(" ", "_").replace("x", "by")
    gen.int("df_" + key, df[k])
    gen.num("ss_" + key, SS[k], 3)
    gen.num("ms_" + key, MS[k], 3)
gen.num("f_irrigation", F_irrigation, 3)
gen.num("f_variety", F_variety, 3)
gen.num("f_inter", F_inter, 3)
gen.num("p_irrigation", stats.f.sf(F_irrigation, df["irrigation"], df["whole-plot error"]), 4)
gen.num("p_variety", stats.f.sf(F_variety, df["varieties"], df["subplot error"]), 4)
gen.num("p_inter", stats.f.sf(F_inter, df["irrigation x variety"], df["subplot error"]), 4)
gen.num("f_irrigation_pooled", MS["irrigation"] / MS_pooled, 3)
gen.num("p_irrigation_pooled",
        stats.f.sf(MS["irrigation"] / MS_pooled, df["irrigation"],
                   df["whole-plot error"] + df["subplot error"]), 2, sci=True)
gen.num("ms_pooled", MS_pooled, 3)
for j, name in enumerate(["flood", "sprinkler", "drip"]):
    gen.num("mean_" + name, irr_means[j], 4)         # exact to 4 dp: 16 yields to 2 dp each
gen.num("diff_drip_flood", irr_means[2] - irr_means[0], 4)
assert abs(round(irr_means[2], 4) - round(irr_means[0], 4)
           - round(irr_means[2] - irr_means[0], 4)) < 1e-12
for k, name in enumerate(["v1", "v2", "v3", "v4"]):
    gen.num("mean_" + name, var_means[k], 3)
gen.num("se_irrig", se_irrig, 4)
gen.num("se_variety", se_variety, 4)
gen.num("se_pooled_irrig", se_pooled_irrig, 4)
gen.num("se_pooled_variety", se_pooled_variety, 4)
gen.num("se_ratio_irrig", se_irrig / se_pooled_irrig, 3)
gen.num("se_ratio_variety", se_variety / se_pooled_variety, 3)
gen.num("se_ratio_variety_inv", se_pooled_variety / se_variety, 3)
gen.num("sigma_w_hat", np.sqrt(max(MS["whole-plot error"] - MS["subplot error"], 0) / b), 4)
gen.num("sigma_s_hat", np.sqrt(MS["subplot error"]), 4)
gen.num("sigma_w", sigma_w, 2)
gen.num("sigma_s", sigma_s, 2)
gen.num("t_drip", (irr_means[2] - irr_means[0]) / se_irrig, 3)
gen.num("se_simple", se_simple, 4)
gen.num("nu_simple", nu_simple, 2)
gen.write()

# ---- the layout ------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(5.6, 2.5))
shade = [COLORS["accent"], COLORS["second"], COLORS["third"]]
rng_layout = np.random.default_rng(11)
for i in range(r):
    order = rng_layout.permutation(a)               # irrigation randomized within the block
    for pos, j in enumerate(order):
        x0, y0 = i * 2.6 + pos * 0.8, 0.0
        ax.add_patch(Rectangle((x0, y0), 0.76, 4.0, facecolor=shade[j],
                               alpha=0.18, edgecolor=shade[j], linewidth=1.4))
        for pos2, k in enumerate(rng_layout.permutation(b)):
            ax.text(x0 + 0.38, 3.5 - pos2, "ABCD"[k], ha="center", va="center", fontsize=7.5)
            if pos2:
                ax.plot([x0, x0 + 0.76], [4.0 - pos2, 4.0 - pos2],
                        color=COLORS["grid"], linewidth=0.5)
        ax.text(x0 + 0.38, 4.15, "FSD"[j], ha="center", va="bottom", fontsize=7.5,
                color=shade[j])
    ax.text(i * 2.6 + 1.2, -0.45, f"block {i + 1}", ha="center", va="top", fontsize=8)
ax.set_xlim(-0.35, (r - 1) * 2.6 + a * 0.8 + 0.35)
ax.set_ylim(-1.35, 4.9)
ax.text((r - 1) * 1.3 + a * 0.4, -1.05,
        "whole plots: F flood, S sprinkler, D drip      subplots: varieties A–D",
        ha="center", va="center", fontsize=8, color=COLORS["muted"])
ax.axis("off")
fig.tight_layout()
fig.savefig(figure_path("ch33", "split_plot_layout"))
