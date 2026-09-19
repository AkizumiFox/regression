"""Chapter 6, Section 10: computing the projection. Normal equations lose twice the digits QR does."""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<solvers>>
def fit_normal_equations(X, y):
    return np.linalg.solve(X.T @ X, X.T @ y)


def fit_qr(X, y):
    Q, R = np.linalg.qr(X)                  # X = QR, Q has orthonormal columns
    return np.linalg.solve(R, Q.T @ y)      # triangular system R b = Q'y
# <</solvers>>


rng = np.random.default_rng(10)
n = 200
t = np.linspace(0, 1, n)
rows = []
# <<experiment>>
for degree in range(1, 14):
    X = np.vander(t, degree + 1, increasing=True)          # 1, t, t^2, ..., t^degree
    beta = rng.normal(size=degree + 1)
    y = X @ beta                                           # exact fit: the true answer is beta
    kappa = np.linalg.cond(X)
    err_ne = np.linalg.norm(fit_normal_equations(X, y) - beta) / np.linalg.norm(beta)
    err_qr = np.linalg.norm(fit_qr(X, y) - beta) / np.linalg.norm(beta)
    rows.append((degree, kappa, err_ne, err_qr))
    print(f"degree {degree:2d}  cond(X) = {kappa:9.2e}   normal eq. {err_ne:8.1e}   QR {err_qr:8.1e}")
# <</experiment>>

rows = np.array(rows)
eps = np.finfo(float).eps
# QR error tracks kappa * eps; normal equations track kappa^2 * eps (until they fail entirely)
mid = rows[(rows[:, 1] > 1e3) & (rows[:, 1] ** 2 * eps < 1e-2)]
assert len(mid) >= 3
assert np.all(mid[:, 2] > 10 * mid[:, 3])
assert np.all(rows[:, 3] < 1e3 * rows[:, 1] * eps)

# Longley's macroeconomic data: the classic ill-conditioned regression
lon = sm.datasets.longley.load_pandas()
XL = sm.add_constant(lon.exog.to_numpy())
kL = np.linalg.cond(XL)
gen = Generated("ch06", "qr_vs_normal", prefix="qr")
gen.num("kappa_longley", kL, 2, sci=True)
gen.num("kappa_longley_sq", kL ** 2, 2, sci=True)
gen.int("digits_longley_qr", max(0, 16 - int(np.ceil(np.log10(kL)))))
gen.int("digits_longley_ne", max(0, 16 - int(np.ceil(2 * np.log10(kL)))))
last = rows[-1]
gen.int("deg_last", last[0])
gen.num("kappa_last", last[1], 1, sci=True)
gen.num("err_ne_last", last[2], 1, sci=True)
gen.num("err_qr_last", last[3], 1, sci=True)
d8 = rows[rows[:, 0] == 8][0]
gen.num("kappa8", d8[1], 1, sci=True)
gen.num("err_ne8", d8[2], 1, sci=True)
gen.num("err_qr8", d8[3], 1, sci=True)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(3.8, 2.7))
k = rows[:, 1]
ax.loglog(k, rows[:, 2], "o-", ms=3, color=COLORS["second"], label="normal equations")
ax.loglog(k, rows[:, 3], "s-", ms=3, color=COLORS["accent"], label="QR")
kk = np.logspace(0, np.log10(k.max()), 50)
ax.loglog(kk, kk * eps, ":", color=COLORS["accent"], lw=0.8, label=r"$\kappa\,\epsilon$")
ax.loglog(kk, np.minimum(kk ** 2 * eps, 1e3), ":", color=COLORS["second"], lw=0.8, label=r"$\kappa^2\epsilon$")
ax.set_xlabel(r"condition number $\kappa(\mathbf{X})$")
ax.set_ylabel(r"relative error in $\hat{\boldsymbol{\beta}}$")
ax.set_ylim(1e-17, 1e3)
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch06", "qr_vs_normal"))
