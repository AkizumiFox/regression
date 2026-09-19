"""Chapter 1, Section 8: singular values and the condition number of a real design matrix.

US quarterly macroeconomic series 1959-2009 (statsmodels.datasets.macrodata, public
domain): real GDP, real investment, real government spending, real disposable income,
population and the unemployment rate, with an intercept.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.macrodata.load_pandas().data
names = ["realgdp", "realinv", "realgovt", "realdpi", "pop", "unemp"]
labels = ["GDP", "investment", "government", "disp. income", "population", "unemployment"]
Z = data[names].to_numpy()
n = len(Z)

# <<svd>>
X = np.column_stack([np.ones(n), Z])                          # as recorded
X_unit = X / np.linalg.norm(X, axis=0)                        # unit-length columns
Xc = np.column_stack([np.ones(n), Z - Z.mean(axis=0)])      # regressors centred
X_std = Xc / np.linalg.norm(Xc, axis=0)                       # ... then length 1

for name, M in [("raw", X), ("unit length", X_unit), ("centred", X_std)]:
    U, s, Vt = np.linalg.svd(M, full_matrices=False)          # M = U diag(s) V'
    print(f"{name:12s} kappa = {s[0] / s[-1]:10.1f}   singular values / s_1:",
          np.array2string(s / s[0], precision=4))
# <</svd>>

sing = {}
for key, M in [("raw", X), ("unit", X_unit), ("std", X_std)]:
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    assert np.allclose(U @ np.diag(s) @ Vt, M)
    assert np.allclose(U.T @ U, np.eye(len(s))) and np.allclose(Vt @ Vt.T, np.eye(len(s)))
    assert np.isclose(s[0] / s[-1], np.linalg.cond(M))
    # squared singular values are the eigenvalues of M'M, so kappa(M'M) = kappa(M)^2
    ev = np.linalg.eigvalsh(M.T @ M)[::-1]
    assert np.allclose(ev, s ** 2, rtol=1e-6)
    assert np.isclose(ev[0] / ev[-1], (s[0] / s[-1]) ** 2, rtol=1e-6)
    # spectral and Frobenius norms
    assert np.isclose(np.linalg.norm(M, 2), s[0])
    assert np.isclose(np.linalg.norm(M, "fro"), np.sqrt(np.sum(s ** 2)))
    sing[key] = (s, Vt)

# Moore-Penrose inverse from the SVD satisfies the four Penrose conditions
U, s, Vt = np.linalg.svd(X_std, full_matrices=False)
Xp = Vt.T @ np.diag(1 / s) @ U.T
assert np.allclose(X_std @ Xp @ X_std, X_std) and np.allclose(Xp @ X_std @ Xp, Xp)
assert np.allclose(X_std @ Xp, (X_std @ Xp).T) and np.allclose(Xp @ X_std, (Xp @ X_std).T)
assert np.allclose(Xp, np.linalg.pinv(X_std))

# Eckart-Young (spectral norm): truncation is the best rank-k approximation
rng = np.random.default_rng(7)
for k in range(1, len(s)):
    Xk = U[:, :k] @ np.diag(s[:k]) @ Vt[:k]
    assert np.isclose(np.linalg.norm(X_std - Xk, 2), s[k])
    for _ in range(20):
        B = rng.normal(size=(n, k)) @ rng.normal(size=(k, X_std.shape[1]))
        assert np.linalg.norm(X_std - B, 2) >= s[k] - 1e-12

# the near-dependence: smallest right singular vector of the centred design
s_std, Vt_std = sing["std"]
v_min = Vt_std[-1]
v_min = v_min * np.sign(v_min[1])                   # sign convention: GDP loading positive
assert abs(v_min[0]) < 1e-8                         # intercept is orthogonal after centring
big = np.argsort(-np.abs(v_min))[:2]
assert set(big) == {1, 4}                           # GDP and disposable income
assert np.isclose(np.linalg.norm(X_std @ v_min), s_std[-1])
r_gdp_dpi = np.corrcoef(Z[:, 0], Z[:, 3])[0, 1]

gen = Generated("ch01", "svd_condition", prefix="svd")
gen.int("n", n)
gen.int("p", X.shape[1])
gen.num("kraw", sing["raw"][0][0] / sing["raw"][0][-1], 2, sci=True)
gen.num("kunit", sing["unit"][0][0] / sing["unit"][0][-1], 0)
gen.num("kstd", s_std[0] / s_std[-1], 0)
gen.num("kstdsq", (s_std[0] / s_std[-1]) ** 2, 0)
gen.num("smin", s_std[-1], 4)
gen.num("vgdp", v_min[1], 2)
gen.num("vdpi", v_min[4], 2)
gen.num("r", r_gdp_dpi, 4)
gen.write()

# ---- figure ------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4), gridspec_kw={"width_ratios": [1.05, 1]})
ax = axes[0]
idx = np.arange(1, X.shape[1] + 1)
for key, lab, col, mk in [("raw", "as recorded", COLORS["second"], "o"),
                          ("unit", "unit-length columns", COLORS["accent"], "s"),
                          ("std", "centred, unit length", COLORS["third"], "^")]:
    s_k = sing[key][0]
    ax.semilogy(idx, s_k / s_k[0], marker=mk, markersize=3.5, color=col, label=lab)
ax.set_xlabel("index $i$")
ax.set_ylabel(r"$\sigma_i/\sigma_1$")
ax.set_title("(a) singular values")
ax.legend(frameon=False, loc="lower left")
ax.set_xticks(idx)
ax = axes[1]
ypos = np.arange(len(labels))[::-1]
ax.barh(ypos, v_min[1:], color=COLORS["accent"], height=0.6)
ax.axvline(0, color=COLORS["ink"], linewidth=0.6)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlabel("loading")
ax.set_title("(b) smallest right singular vector")
ax.set_xlim(-1, 1)
fig.tight_layout()
fig.savefig(figure_path("ch01", "svd_condition"))
