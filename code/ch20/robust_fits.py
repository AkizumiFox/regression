"""Chapter 20, Section 6: robust fits as a diagnostic.

Huber and bisquare M-estimates by iteratively reweighted least squares (IRLS), least trimmed
squares (LTS) by concentration steps, and what each reports on (1) the masking data of
Section 5 and (2) Engel's 1857 household budgets (public domain, statsmodels.datasets.engel).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import integrate, stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<irls>>
def huber_weight(u, k=1.345):
    return np.minimum(1.0, k / np.maximum(np.abs(u), 1e-12))


def bisquare_weight(u, c=4.685):
    return np.where(np.abs(u) < c, (1 - (u / c) ** 2) ** 2, 0.0)


def mad_scale(r):
    return np.median(np.abs(r)) / 0.6745           # consistent for sigma at the normal


def irls(X, y, weight, beta=None, scale=None, tol=1e-10, max_iter=500):
    """M-estimate by iteratively reweighted least squares.

    Starts from least squares unless `beta` is given. If `scale` is None it is re-estimated
    from the current residuals at every step; otherwise it is held fixed."""
    if beta is None:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for _ in range(max_iter):
        r = y - X @ beta
        s = mad_scale(r) if scale is None else scale
        w = weight(r / s)
        sw = np.sqrt(w)
        new, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        if np.max(np.abs(new - beta)) < tol * (1 + np.max(np.abs(beta))):
            beta = new
            break
        beta = new
    r = y - X @ beta
    s = mad_scale(r) if scale is None else scale
    return beta, weight(r / s), s
# <</irls>>


def huber_rho(u, k=1.345):
    a = np.abs(u)
    return np.where(a <= k, u ** 2 / 2, k * a - k ** 2 / 2)


def bisquare_rho(u, c=4.685):
    return np.where(np.abs(u) <= c, c ** 2 / 6 * (1 - (1 - (u / c) ** 2) ** 3), c ** 2 / 6)


# <<lts>>
def lts(X, y, q, n_starts=500, seed=0):
    """Least trimmed squares: minimize the sum of the q smallest squared residuals.

    Random elemental starts (p cases), each followed by concentration steps: fit least
    squares to the current q cases, then keep the q cases with the smallest residuals."""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    best, best_obj = None, np.inf
    for _ in range(n_starts):
        idx = rng.choice(n, size=p, replace=False)
        if np.linalg.matrix_rank(X[idx]) < p:
            continue
        beta = np.linalg.solve(X[idx], y[idx])
        for _ in range(100):                       # concentration steps
            keep = np.argsort((y - X @ beta) ** 2)[:q]
            new, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            if np.allclose(new, beta):
                break
            beta = new
        obj = np.sort((y - X @ beta) ** 2)[:q].sum()
        if obj < best_obj:
            best, best_obj = beta, obj
    return best, best_obj
# <</lts>>

# ---- efficiency at the normal of the two tuning constants --------------------------------
def efficiency(psi, dpsi, kink):
    """Asymptotic efficiency at the normal of the location M-estimate: (E psi')^2 / E psi^2."""
    pts = [-kink, kink]
    num = integrate.quad(lambda u: dpsi(u) * stats.norm.pdf(u), -12, 12, points=pts)[0] ** 2
    den = integrate.quad(lambda u: psi(u) ** 2 * stats.norm.pdf(u), -12, 12, points=pts)[0]
    return num / den


k, c = 1.345, 4.685
eff_huber = efficiency(lambda u: np.clip(u, -k, k), lambda u: float(abs(u) <= k), k)
eff_bisq = efficiency(lambda u: u * (1 - (u / c) ** 2) ** 2 * (abs(u) <= c),
                      lambda u: ((1 - (u / c) ** 2) * (1 - 5 * (u / c) ** 2)) * (abs(u) <= c), c)
assert abs(eff_huber - 0.95) < 0.002 and abs(eff_bisq - 0.95) < 0.002

# ---- (1) the masking data of Section 5 -------------------------------------------------
# <<masking>>
rng = np.random.default_rng(20)
n0, m = 30, 4
x_clean = np.sort(rng.uniform(0, 10, n0))
y_clean = 2 + 0.5 * x_clean + rng.normal(0, 0.5, n0)
x = np.r_[x_clean, np.full(m, 20.0)]
y = np.r_[y_clean, np.full(m, 2.0)]
X = np.column_stack([np.ones(len(x)), x])
n, p = X.shape

b_ls, *_ = np.linalg.lstsq(X, y, rcond=None)
b_hub, w_hub, _ = irls(X, y, huber_weight)
b_lts, _ = lts(X, y, q=(n + p + 1) // 2)
s_lts = mad_scale(y - X @ b_lts)
b_bis, w_bis, _ = irls(X, y, bisquare_weight, beta=b_lts, scale=s_lts)
for name, b in [("least squares", b_ls), ("Huber", b_hub), ("LTS", b_lts), ("bisquare from LTS", b_bis)]:
    print(f"{name:18s} intercept {b[0]:6.3f}  slope {b[1]:6.3f}")
print("scaled LTS residuals of the cluster:", ((y - X @ b_lts) / s_lts)[-m:].round(1))
# <</masking>>

b_clean = np.linalg.lstsq(X[:n0], y_clean, rcond=None)[0]
assert b_hub[1] < 0.1 and abs(b_lts[1] - 0.5) < 0.1 and abs(b_bis[1] - b_clean[1]) < 0.02
assert np.all(w_bis[-m:] == 0)
z_lts = (y - X @ b_lts) / s_lts
assert np.all(np.abs(z_lts[-m:]) > 10) and np.all(np.abs(z_lts[:n0]) < 3)
b_bis_ls, _, _ = irls(X, y, bisquare_weight)       # redescending, started badly
assert b_bis_ls[1] < 0.1

# IRLS never increases the objective when the scale is held fixed
for weight, rho in [(huber_weight, huber_rho), (bisquare_weight, bisquare_rho)]:
    s_fix = mad_scale(y - X @ b_ls)
    beta = b_ls.copy()
    objs = [rho((y - X @ beta) / s_fix).sum()]
    for _ in range(30):
        w = weight((y - X @ beta) / s_fix)
        sw = np.sqrt(w)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        objs.append(rho((y - X @ beta) / s_fix).sum())
    assert np.all(np.diff(objs) <= 1e-9)

# concentration steps never increase the trimmed sum of squares
hh = (n + p + 1) // 2
beta = np.linalg.solve(X[[0, n - 1]], y[[0, n - 1]])
vals = []
for _ in range(10):
    keep = np.argsort((y - X @ beta) ** 2)[:hh]
    vals.append(np.sort((y - X @ beta) ** 2)[:hh].sum())
    beta, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
assert np.all(np.diff(vals) <= 1e-9)

# one bad leverage point breaks least squares and Huber alike
x_one = np.r_[x_clean, 50.0]
X_one = np.column_stack([np.ones(n0 + 1), x_one])
y_one = np.r_[y_clean, 0.0]
b_one_ls = np.linalg.lstsq(X_one, y_one, rcond=None)[0]
b_one_hub = irls(X_one, y_one, huber_weight)[0]
assert b_one_ls[1] < 0 and b_one_hub[1] < 0

# ---- (2) Engel's household budgets ----------------------------------------------------
# <<engel>>
engel = sm.datasets.engel.load_pandas().data
inc = engel["income"].to_numpy()
food = engel["foodexp"].to_numpy()
XE = np.column_stack([np.ones(len(inc)), inc])
fits = {
    "least squares": sm.OLS(food, XE).fit(),
    "Huber": sm.RLM(food, XE, M=sm.robust.norms.HuberT()).fit(),
    "bisquare": sm.RLM(food, XE, M=sm.robust.norms.TukeyBiweight()).fit(),
}
for name, f in fits.items():
    print(f"{name:14s} slope {f.params[1]:.4f}")
w = fits["bisquare"].weights
print("households with bisquare weight 0 (row numbers from 1):", np.where(w == 0)[0] + 1)
# <</engel>>

# statsmodels RLM rescales by the MAD about zero at each step: our IRLS agrees
b_hub_E, w_hub_E, _ = irls(XE, food, huber_weight)
assert np.allclose(b_hub_E, fits["Huber"].params, rtol=1e-5)
zero = np.where(w == 0)[0]
QE, _ = np.linalg.qr(XE)
hE = np.sum(QE ** 2, axis=1)
top = int(np.argmax(hE))
assert top in zero

lfit = {
    "ls": sm.OLS(np.log(food), np.column_stack([np.ones(len(inc)), np.log(inc)])).fit(),
    "hub": sm.RLM(np.log(food), np.column_stack([np.ones(len(inc)), np.log(inc)]), M=sm.robust.norms.HuberT()).fit(),
    "bis": sm.RLM(np.log(food), np.column_stack([np.ones(len(inc)), np.log(inc)]), M=sm.robust.norms.TukeyBiweight()).fit(),
}
assert lfit["bis"].weights.min() > 0.1

gen = Generated("ch20", "robust_fits", prefix="rb")
gen.num("eff_huber", eff_huber, 3)
gen.num("eff_bisq", eff_bisq, 3)
gen.num("slope_ls", b_ls[1], 3)
gen.num("slope_hub", b_hub[1], 3)
gen.num("slope_lts", b_lts[1], 3)
gen.num("slope_bis", b_bis[1], 3)
gen.num("slope_bis_ls", b_bis_ls[1], 3)
gen.num("slope_clean", b_clean[1], 3)
gen.num("int_bis", b_bis[0], 3)
gen.num("int_clean", b_clean[0], 3)
gen.num("s_lts", s_lts, 3)
gen.num("z_cluster", z_lts[-1], 1)
gen.num("z_clean_max", np.abs(z_lts[:n0]).max(), 2)
gen.num("w_hub_cluster", w_hub[-1], 3)
gen.num("one_ls", b_one_ls[1], 3)
gen.num("one_hub", b_one_hub[1], 3)
gen.int("h_lts", hh)
gen.num("E_ls", fits["least squares"].params[1], 3)
gen.num("E_hub", fits["Huber"].params[1], 3)
gen.num("E_bis", fits["bisquare"].params[1], 3)
gen.int("E_nzero", len(zero))
gen.text("E_zero", ", ".join(str(i + 1) for i in zero))
gen.int("E_top", top + 1)
gen.num("E_top_inc", inc[top], 0)
gen.num("E_top_food", food[top], 0)
gen.num("E_top_whub", fits["Huber"].weights[top], 3)
gen.num("L_ls", lfit["ls"].params[1], 3)
gen.num("L_hub", lfit["hub"].params[1], 3)
gen.num("L_bis", lfit["bis"].params[1], 3)
gen.num("L_wmin", lfit["bis"].weights.min(), 2)
gen.write()

# ---- figure ----------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.9, 2.6))
ax = axes[0]
ax.scatter(x_clean, y_clean, s=9, color=COLORS["accent"], linewidths=0)
ax.scatter([20.0], [2.0], s=40, marker="D", color=COLORS["second"], linewidths=0)
g = np.linspace(0, 21, 2)
for b, lab, col, ls in [(b_ls, "least squares", COLORS["second"], "--"),
                        (b_hub, "Huber", COLORS["thread"], ":"),
                        (b_lts, "LTS", COLORS["third"], "-."),
                        (b_bis, "bisquare from LTS", COLORS["accent"], "-")]:
    ax.plot(g, b[0] + b[1] * g, color=col, linestyle=ls, label=lab)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("(a) masking data")
ax = axes[1]
ax.scatter(inc, food, s=6, color=COLORS["accent"], alpha=0.6, linewidths=0)
ax.scatter(inc[zero], food[zero], s=26, facecolors="none", edgecolors=COLORS["second"], linewidths=1.0,
           label="bisquare weight 0")
g = np.linspace(0, inc.max() * 1.02, 2)
for key, col, ls in [("least squares", COLORS["second"], "--"), ("Huber", COLORS["thread"], ":"),
                     ("bisquare", COLORS["accent"], "-")]:
    b = fits[key].params
    ax.plot(g, b[0] + b[1] * g, color=col, linestyle=ls, label=key)
ax.set_xlabel("household income")
ax.set_ylabel("food expenditure")
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("(b) Engel's budgets")
fig.tight_layout()
fig.savefig(figure_path("ch20", "robust_fits"))
