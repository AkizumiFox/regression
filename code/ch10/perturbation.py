"""Chapter 10, Section 5: sensitivity of least squares.

(1) The perturbation bound of the section, checked on random problems.
(2) A problem on which the kappa^2 ||e|| term is attained.
(3) Accuracy of three solvers against kappa, with zero and with large residuals.
(4) Longley's data (statsmodels.datasets.longley, public domain): exact solution in
    rational arithmetic, digits achieved by each method, column scaling, and the effect of
    rounding the data to the digits in which they were published.
"""
from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch10", "perturbation", prefix="pert")
u = np.finfo(float).eps / 2


def ls(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


# ---- (1) the bound of thm-cmp-perturbation on random problems ------------------------------
# <<bound>>
def perturbation_bound(X, y, eps):
    """Right-hand side of the bound for ||db||/||b|| when ||dX|| <= eps||X||, ||dy|| <= eps||y||."""
    s = np.linalg.svd(X, compute_uv=False)
    kappa = s[0] / s[-1]
    b = ls(X, y)
    eta = np.linalg.norm(y - X @ b) / (s[0] * np.linalg.norm(b))
    q = kappa * eps
    return q / (1 - q) * (2 + eta) + kappa * q / (1 - q) ** 2 * eta
# <</bound>>


rng = np.random.default_rng(105)
worst_ratio = 0.0
for _ in range(300):
    n, p = 40, 5
    U0, _ = np.linalg.qr(rng.normal(size=(n, p)))
    V0, _ = np.linalg.qr(rng.normal(size=(p, p)))
    X = U0 @ np.diag(np.logspace(0, -rng.uniform(0, 4), p)) @ V0.T
    y = X @ rng.normal(size=p) + rng.uniform(0, 2) * rng.normal(size=n) * np.linalg.norm(X, 2)
    eps = 10.0 ** rng.uniform(-10, -6)
    E = rng.normal(size=(n, p))
    E *= eps * np.linalg.norm(X, 2) / np.linalg.norm(E, 2)
    f = rng.normal(size=n)
    f *= eps * np.linalg.norm(y) / np.linalg.norm(f)
    b, bt = ls(X, y), ls(X + E, y + f)
    actual = np.linalg.norm(bt - b) / np.linalg.norm(b)
    bound = perturbation_bound(X, y, eps)
    assert actual <= bound * (1 + 1e-6)
    worst_ratio = max(worst_ratio, actual / bound)
gen.num("bound_worst", worst_ratio, 2)

# residual bound: ||de||/||y|| <= eps (1 + 2 kappa / (1 - kappa eps))
for _ in range(100):
    X = rng.normal(size=(30, 4)) @ np.diag([1, 1, 1e-2, 1e-3])
    y = rng.normal(size=30)
    eps = 1e-8
    E = rng.normal(size=X.shape)
    E *= eps * np.linalg.norm(X, 2) / np.linalg.norm(E, 2)
    kap = np.linalg.cond(X)
    de = (y - (X + E) @ ls(X + E, y)) - (y - X @ ls(X, y))
    assert np.linalg.norm(de) / np.linalg.norm(y) <= eps * (1 + 2 * kap / (1 - kap * eps))

# ---- (2) the kappa^2 term is attained ------------------------------------------------------
# <<attained>>
delta, eps = 1e-3, 1e-10
X = np.array([[1.0, 0.0], [0.0, delta], [0.0, 0.0]])
y = np.array([0.0, delta, 1.0])            # b = (0, 1), residual e = (0, 0, 1)
E = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, eps]])   # ||E|| = eps ||X||
b = ls(X, y)
b_tilde = ls(X + E, y)
print("kappa =", np.linalg.cond(X))
print("relative change in b:", np.linalg.norm(b_tilde - b) / np.linalg.norm(b))
print("kappa^2 * eps       :", np.linalg.cond(X) ** 2 * eps)
# <</attained>>

rel = np.linalg.norm(b_tilde - b) / np.linalg.norm(b)
k2e = np.linalg.cond(X) ** 2 * eps
assert 0.9 < rel / k2e < 1.1
assert rel <= perturbation_bound(X, y, eps)
gen.num("att_rel", rel, 2, sci=True)
gen.num("att_k2e", k2e, 2, sci=True)
gen.num("att_kappa", np.linalg.cond(X), 0)


# ---- (3) four solvers against kappa, zero and large residual --------------------------------
# <<solvers>>
def by_cholesky(X, y):
    L = np.linalg.cholesky(X.T @ X)
    return np.linalg.solve(L.T, np.linalg.solve(L, X.T @ y))


def by_qr(X, y):
    Q, R = np.linalg.qr(X)
    return np.linalg.solve(R, Q.T @ y)


def by_svd(X, y):
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return Vt.T @ ((U.T @ y) / s)
# <</solvers>>


def test_problem(kappa, resid_size, n=100, p=10, seed=0):
    r = np.random.default_rng(seed)
    U0, _ = np.linalg.qr(r.normal(size=(n, n)))
    V0, _ = np.linalg.qr(r.normal(size=(p, p)))
    X = U0[:, :p] @ np.diag(np.logspace(0, -np.log10(kappa), p)) @ V0.T
    beta = r.normal(size=p)
    beta /= np.linalg.norm(beta)
    e = U0[:, p:] @ r.normal(size=n - p)                  # orthogonal to C(X)
    e *= resid_size / np.linalg.norm(e)                   # ||X||_2 = ||beta|| = 1, so eta = resid_size
    return X, X @ beta + e, beta


kappas = 10.0 ** np.arange(1, 13)
results = {}
for resid in (0.0, 1.0):
    rows = []
    for kap in kappas:
        X, y, beta = test_problem(kap, resid)
        eta = np.linalg.norm(y - X @ beta) / (np.linalg.norm(X, 2) * np.linalg.norm(beta))
        assert abs(eta - resid) < 1e-8                       # the stated eta of the panel
        errs = []
        for f in (by_cholesky, by_qr, by_svd):
            try:
                errs.append(np.linalg.norm(f(X, y) - beta) / np.linalg.norm(beta))
            except np.linalg.LinAlgError:
                errs.append(np.nan)                          # Cholesky breaks down
        rows.append(errs)
    results[resid] = np.array(rows)

z, big = results[0.0], results[1.0]
# zero residual: QR and SVD track kappa u; Cholesky tracks kappa^2 u
mask = kappas ** 2 * u < 1e-2
assert np.all(z[:, 1] < 50 * kappas * u) and np.all(z[:, 2] < 50 * kappas * u)
assert np.all(z[mask & (kappas >= 1e3), 0] > 10 * z[mask & (kappas >= 1e3), 1])
# large residual: even QR and SVD show the kappa^2 term
assert np.all(big[kappas >= 1e4, 1] > 10 * kappas[kappas >= 1e4] * u)
assert np.all(big[:, 1] < 50 * (kappas * u + kappas ** 2 * u))
i6 = np.where(kappas == 1e6)[0][0]
gen.num("z_chol6", z[i6, 0], 1, sci=True)
gen.num("z_qr6", z[i6, 1], 1, sci=True)
gen.num("b_qr6", big[i6, 1], 1, sci=True)
gen.num("b_svd6", big[i6, 2], 1, sci=True)

# ---- (4) Longley --------------------------------------------------------------------------
# <<longley>>
lon = sm.datasets.longley.load_pandas()
XL = sm.add_constant(lon.exog.to_numpy())
yL = lon.endog.to_numpy()
names = ["const"] + list(lon.exog.columns)
print("kappa(X) as recorded:", f"{np.linalg.cond(XL):.3e}")
# <</longley>>


def exact_ls(X, y):
    """Solve the normal equations exactly in rational arithmetic (decimal data)."""
    Xf = [[Fraction(repr(float(v))) for v in row] for row in X]
    yf = [Fraction(repr(float(v))) for v in y]
    p = len(Xf[0])
    A = [[sum(r[i] * r[j] for r in Xf) for j in range(p)] + [sum(r[i] * t for r, t in zip(Xf, yf))] for i in range(p)]
    for k in range(p):                                   # Gauss-Jordan elimination
        piv = A[k][k]
        A[k] = [a / piv for a in A[k]]
        for i in range(p):
            if i != k and A[i][k] != 0:
                f = A[i][k]
                A[i] = [a - f * b for a, b in zip(A[i], A[k])]
    return np.array([float(A[i][p]) for i in range(p)])


b_exact = exact_ls(XL, yL)


def digits(b):
    return -np.log10(np.max(np.abs(b - b_exact) / np.abs(b_exact)))


methods = {
    "normal equations (Cholesky)": by_cholesky(XL, yL),
    "Householder QR": by_qr(XL, yL),
    "SVD": by_svd(XL, yL),
    "statsmodels OLS": sm.OLS(yL, XL).fit().params,
}
dig = {k: digits(v) for k, v in methods.items()}
for k, v in dig.items():
    print(f"{k:28s} correct digits (worst coefficient): {v:.1f}")
assert dig["Householder QR"] > 7 and dig["SVD"] > 7
assert dig["normal equations (Cholesky)"] < dig["Householder QR"] - 3

# column scaling
D0 = 1 / np.linalg.norm(XL, axis=0)
Xs = XL * D0
Zc = lon.exog.to_numpy() - lon.exog.to_numpy().mean(axis=0)
Xc = np.column_stack([np.ones(len(yL)), Zc])
Xcs = Xc / np.linalg.norm(Xc, axis=0)
k_raw, k_unit, k_cs = np.linalg.cond(XL), np.linalg.cond(Xs), np.linalg.cond(Xcs)
b_scaled_ne = by_cholesky(Xs, yL) * D0                    # solve in scaled variables, map back
dig_scaled_ne = digits(b_scaled_ne)
print(f"kappa: raw {k_raw:.2e}, unit columns {k_unit:.2e}, centred+unit {k_cs:.2e}")
print(f"normal equations after scaling columns: {dig_scaled_ne:.1f} digits")
# scaling the columns changes kappa enormously but not the accuracy: both methods behave as
# if the columns were already equilibrated (the scaled condition number is what counts)
assert abs(dig_scaled_ne - dig["normal equations (Cholesky)"]) < 1.5
assert dig["normal equations (Cholesky)"] > 16 - 2 * np.log10(k_unit) - 2
assert dig["Householder QR"] > 16 - np.log10(k_unit) - 2
assert dig["normal equations (Cholesky)"] > 16 - 2 * np.log10(k_raw) + 10     # raw kappa is far too pessimistic
for k, v in [("raw", k_raw), ("unit", k_unit), ("cs", k_cs)]:
    gen.num(f"k_{k}", v, 2, sci=True)
gen.num("dig_ne", dig["normal equations (Cholesky)"], 1)
gen.num("dig_qr", dig["Householder QR"], 1)
gen.num("dig_svd", dig["SVD"], 1)
gen.num("dig_sm", dig["statsmodels OLS"], 1)
gen.num("dig_ne_scaled", dig_scaled_ne, 1)
gen.num("b_gnp_exact", b_exact[2], 6)

# van der Sluis: unit columns are within sqrt(p) of the best diagonal scaling
for _ in range(500):
    Dr = np.exp(rng.normal(scale=3, size=XL.shape[1]))
    assert k_unit <= np.sqrt(XL.shape[1]) * np.linalg.cond(XL * Dr) * (1 + 1e-10)

# ---- rounding the published data: statistical versus numerical error ---------------------
# <<rounding>>
rng = np.random.default_rng(106)
half_unit = np.array([0.0, 0.05, 0.5, 0.5, 0.5, 0.5, 0.0])   # half a unit in the last digit
fitL = sm.OLS(yL, XL).fit()                                  # year and constant are exact
reps = 1000
shifts = []
for _ in range(reps):
    Xp = XL + rng.uniform(-1, 1, size=XL.shape) * half_unit
    shifts.append(ls(Xp, yL) - fitL.params)
spread = np.std(shifts, axis=0)
ratio = spread / fitL.bse
for name, r in zip(names, ratio):
    print(f"{name:8s} sd of coefficient under rounding / standard error = {r:.2f}")
# <</rounding>>

gen.num("round_min", ratio.min(), 3)
gen.num("round_max", ratio.max(), 3)
gen.int("round_reps", reps)
ne_in_se = np.max(np.abs(methods["normal equations (Cholesky)"] - b_exact) / fitL.bse)
# numerical error of even the normal equations << rounding of the data << sampling error
assert ne_in_se < 1e-4 < ratio.min() and ratio.max() < 0.1
gen.num("ne_in_se", ne_in_se, 0, sci=True)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6), sharey=True)
labels = ["normal equations", "Householder QR", "SVD"]
cols = [COLORS["second"], COLORS["accent"], COLORS["third"]]
marks = ["o-", "s-", "^-"]
for ax, (resid, R) in zip(axes, results.items()):
    for j in range(3):
        ax.loglog(kappas, R[:, j], marks[j], ms=3, color=cols[j], label=labels[j])
    ax.loglog(kappas, kappas * u, ":", color=COLORS["muted"], lw=0.8, label=r"$\kappa u$")
    ax.loglog(kappas, np.minimum(kappas ** 2 * u, 1e2), "--", color=COLORS["muted"], lw=0.8, label=r"$\kappa^2 u$")
    ax.set_xlabel(r"condition number $\kappa(\mathbf{X})$")
    ax.set_title("(a) zero residual" if resid == 0 else r"(b) $\eta = 1$")
axes[0].set_ylabel(r"relative error in $\hat{\boldsymbol{\beta}}$")
axes[0].set_ylim(1e-17, 1e2)
axes[0].legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(figure_path("ch10", "residual_size"))
