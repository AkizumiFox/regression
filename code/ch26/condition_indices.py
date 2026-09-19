"""Chapter 26, Section 3: condition indices and variance-decomposition proportions.

The design of a consumption regression on the quarterly US macroeconomic data (public domain,
statsmodels.datasets.macrodata): intercept, real disposable income, population, CPI, M1, the
Treasury bill rate and unemployment. Scaled to unit column lengths without centring, as
Belsley, Kuh and Welsch recommend; compared with the centred-and-scaled version. Checks the
locating inequality of the variance-decomposition theorem, the coefficient-of-variation bound,
and the deletion bound. Figure: the table of proportions.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch26", "condition_indices", prefix="ci")

# <<bkw>>
import numpy as np
import statsmodels.api as sm
macro = sm.datasets.macrodata.load_pandas().data
names = ["const", "realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
X = sm.add_constant(macro[names[1:]].to_numpy())
X_s = X / np.linalg.norm(X, axis=0)                   # unit-length columns, not centred

U, mu, Vt = np.linalg.svd(X_s, full_matrices=False)   # X_s = U diag(mu) V'
eta = mu[0] / mu                                      # condition indices
phi = Vt.T ** 2 / mu**2                               # phi[j, k] = v_jk^2 / mu_k^2
prop = (phi / phi.sum(axis=1, keepdims=True)).T       # row k: proportions for index eta_k

print("index   " + " ".join(f"{nm[:7]:>7s}" for nm in names))
for k in range(len(mu)):
    print(f"{eta[k]:7.1f} " + " ".join(f"{p:7.3f}" for p in prop[k]))
# <</bkw>>

p = len(names)
assert np.allclose(prop.sum(axis=0), 1)
assert np.isclose(np.sum(mu**2), p) and 1 <= mu[0] ** 2 <= p
# the diagonal of (X_s'X_s)^{-1} is sum_k v_jk^2/mu_k^2 and is at most kappa^2
Ginv_s = np.linalg.inv(X_s.T @ X_s)
assert np.allclose(np.diag(Ginv_s), phi.sum(axis=1))
assert np.all(np.diag(Ginv_s) <= 1 / mu[-1] ** 2) and 1 / mu[-1] ** 2 <= eta[-1] ** 2
# locating inequality, for S = the last one, two and three indices
for s in (1, 2, 3):
    S = list(range(p - s, p))
    mu_S, mu_L = mu[S].max(), mu[: p - s].min()
    w = (Vt.T[:, S] ** 2).sum(axis=1)
    lower = w / (w + (1 - w) * mu_S**2 / mu_L**2)
    assert np.all(prop[S].sum(axis=0) >= lower - 1e-12)
# each near-dependency: || X_s v_k || = mu_k
assert np.allclose(np.linalg.norm(X_s @ Vt.T, axis=0), mu)

# <<centred>>
Z = X[:, 1:] - X[:, 1:].mean(axis=0)
Z_s = Z / np.linalg.norm(Z, axis=0)                   # centred, then unit length
mu_c = np.linalg.svd(Z_s, compute_uv=False)
print("scaled condition number, not centred:", round(eta[-1], 1))
print("scaled condition number, centred    :", round(mu_c[0] / mu_c[-1], 1))
cv = X[:, 1:].std(axis=0) / X[:, 1:].mean(axis=0)
print("coefficients of variation:", dict(zip(names[1:], np.round(cv, 3))))
# <</centred>>

# every near-dependency visible after centring is visible before: lambda_min(R) >= mu_p^2
R = Z_s.T @ Z_s
assert np.linalg.eigvalsh(R)[0] >= mu[-1] ** 2
# the coefficient-of-variation bound mu_p^2 <= CV^2/(CV^2 + 2) for every column
assert np.all(mu[-1] ** 2 <= cv**2 / (cv**2 + 2))
j_pop = names.index("pop") - 1
eta_cv = np.sqrt((cv[j_pop] ** 2 + 2) / cv[j_pop] ** 2)

# <<auxiliary>>
def r_squared(target, others):
    fit = sm.OLS(macro[target], sm.add_constant(macro[others])).fit()
    return fit.rsquared

print("pop on realdpi     : R^2 =", round(r_squared("pop", ["realdpi"]), 4))
print("cpi on realdpi, m1 : R^2 =", round(r_squared("cpi", ["realdpi", "m1"]), 4))
# <</auxiliary>>
r2_pop = r_squared("pop", ["realdpi"])
r2_cpi = r_squared("cpi", ["realdpi", "m1"])

# ---- deletion: one case changes lambda_min by at most the factor (1 - h_ii) --------------
XtX = X_s.T @ X_s
lam_min = np.linalg.eigvalsh(XtX)[0]
Q = np.linalg.qr(X_s)[0]
h = np.sum(Q**2, axis=1)
ratios = []
for i in range(len(X_s)):
    Xi = np.delete(X_s, i, axis=0)
    lam_i = np.linalg.eigvalsh(Xi.T @ Xi)[0]
    assert lam_i >= (1 - h[i]) * lam_min - 1e-15
    ratios.append(lam_min / lam_i)
    # Sherman-Morrison form of the variance increase, for the worst coefficient
    Gi = np.linalg.inv(Xi.T @ Xi)
    assert np.all(np.diag(Gi) <= np.diag(Ginv_s) / (1 - h[i]) + 1e-9)
ratios = np.array(ratios)
i_max = int(np.argmax(h))
quarter = f"{int(macro.year[i_max])}Q{int(macro.quarter[i_max])}"

gen.int("n", len(X))
for k in range(p):
    gen.num(f"eta{k + 1}", eta[k], 1)
for k in range(p):
    for j, nm in enumerate(names):
        gen.num(f"prop_{k + 1}_{nm}", prop[k, j], 3)
gen.num("eta_centred", mu_c[0] / mu_c[-1], 1)
gen.num("cv_pop", cv[j_pop], 3)
gen.num("eta_cv", eta_cv, 1)
gen.num("r2_pop", r2_pop, 4)
gen.num("r2_cpi", r2_cpi, 4)
gen.num("mu_max2", mu[0] ** 2, 3)
v_last = Vt[-1] * np.sign(Vt[-1, 2])
gen.num("v7_const", v_last[0], 3)
gen.num("v7_dpi", v_last[1], 3)
gen.num("v7_pop", v_last[2], 3)
gen.num("mu7", mu[-1], 4)
gen.num("h_max", h[i_max], 3)
gen.text("h_max_quarter", quarter)
gen.num("h_factor", 1 / (1 - h[i_max]), 3)
gen.num("ratio_max", ratios.max(), 3)
gen.write()

# ---- figure: proportions as a dot table ------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(5.2, 2.9))
labels = ["intercept", "income", "population", "CPI", "M1", "T-bill", "unemployment"]
for k in range(p):
    for j in range(p):
        if prop[k, j] >= 0.005:
            col = COLORS["second"] if (eta[k] >= 30 and prop[k, j] >= 0.5) else COLORS["accent"]
            ax.scatter(j, k, s=260 * prop[k, j], color=col, linewidths=0)
        if prop[k, j] >= 0.1:
            ax.annotate(f"{prop[k, j]:.3f}", (j, k), xytext=(10, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=6.5, color=COLORS["ink"])
ax.set_xticks(range(p))
ax.set_xticklabels(labels, rotation=30, ha="right")
ax.set_yticks(range(p))
ax.set_yticklabels([f"{e:.1f}" for e in eta])
ax.set_ylabel("condition index")
ax.set_ylim(p - 0.4, -0.6)
ax.set_xlim(-0.6, p - 0.2)
for spine in ("left", "bottom"):
    ax.spines[spine].set_visible(False)
ax.tick_params(length=0)
ax.axhline(3.5, color=COLORS["grid"], linewidth=0.8, zorder=0)
ax.annotate("index 30", (-0.55, 3.5), xytext=(0, 2), textcoords="offset points",
            ha="left", va="bottom", fontsize=6.5, color=COLORS["muted"])
fig.tight_layout()
fig.savefig(figure_path("ch26", "variance_proportions"))
