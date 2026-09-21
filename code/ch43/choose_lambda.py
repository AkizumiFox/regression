"""Chapter 43, Section 5: cross-validation, GCV and corrected AIC for a penalized spline.

One data set for the criterion curves, then a simulation that compares the three rules,
and finally the same simulation with AR(1) errors, where all three break down.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


# <<criteria>>
def smoother(B, P, lam):
    """The smoother matrix of the penalized fit with basis B and penalty matrix P."""
    return B @ np.linalg.solve(B.T @ B + lam * P, B.T)


def criteria(B, P, y, lam):
    """Leave-one-out CV by the leverage shortcut, GCV, and the corrected AIC."""
    n = len(y)
    S = smoother(B, P, lam)
    resid = y - S @ y
    df = np.trace(S)
    cv = np.mean((resid / (1 - np.diag(S))) ** 2)
    gcv = np.mean(resid**2) / (1 - df / n) ** 2
    aicc = np.log(np.mean(resid**2)) + 1 + 2 * (df + 1) / (n - df - 2)
    return df, cv, gcv, aicc
# <</criteria>>


def select(B, P, y, lams, which):
    values = np.array([criteria(B, P, y, lam)[which] for lam in lams])
    return lams[int(np.argmin(values))]


n, sigma, K, d, k = 120, 0.25, 25, 3, 2
x = (np.arange(1, n + 1) - 0.5) / n
kn = uniform_knots(0.0, 1.0, K, d)
B = bspline_basis(x, kn, d)
Dk = np.diff(np.eye(B.shape[1]), n=k, axis=0)
P = Dk.T @ Dk
rng = np.random.default_rng(4343)
y = f_true(x) + sigma * rng.normal(size=n)

# --- the shortcut is exact: compare with n refits -----------------------------
lam0 = 1e-4
S0 = smoother(B, P, lam0)
loo_short = (y - S0 @ y) / (1 - np.diag(S0))
loo_brute = np.empty(n)
for i in range(n):
    keep = np.arange(n) != i
    Bi, yi = B[keep], y[keep]
    gi = np.linalg.solve(Bi.T @ Bi + lam0 * P, Bi.T @ yi)
    loo_brute[i] = y[i] - B[i] @ gi
shortcut_gap = np.max(np.abs(loo_short - loo_brute))
print("largest gap between the shortcut and n refits:", shortcut_gap)
assert shortcut_gap < 1e-8

# --- the criterion curves on this data set ------------------------------------
lams = np.geomspace(1e-9, 1e3, 121)
table = np.array([criteria(B, P, y, lam) for lam in lams])
dfs, cvs, gcvs, aiccs = table.T
pick = {"CV": lams[np.argmin(cvs)], "GCV": lams[np.argmin(gcvs)], "AICc": lams[np.argmin(aiccs)]}
df_pick = {name: criteria(B, P, y, lam)[0] for name, lam in pick.items()}
print("degrees of freedom chosen:", {k_: round(v, 2) for k_, v in df_pick.items()})

# --- simulation: independent errors, then AR(1) errors ------------------------
def one_run(err):
    yy = f_true(x) + err
    out = {}
    for name, which in (("CV", 1), ("GCV", 2), ("AICc", 3)):
        lam = select(B, P, yy, lams, which)
        fit = smoother(B, P, lam) @ yy
        out[name] = (np.trace(smoother(B, P, lam)), np.mean((fit - f_true(x)) ** 2))
    return out


reps = 200
rho = 0.7
res = {name: {"df": [], "mse": []} for name in ("CV", "GCV", "AICc")}
res_ar = {name: {"df": [], "mse": []} for name in ("CV", "GCV", "AICc")}
for _ in range(reps):
    e = sigma * rng.normal(size=n)
    for name, (df_, mse_) in one_run(e).items():
        res[name]["df"].append(df_)
        res[name]["mse"].append(mse_)
    w = rng.normal(size=n) * sigma * np.sqrt(1 - rho**2)
    a = np.empty(n)
    a[0] = sigma * rng.normal()
    for i in range(1, n):
        a[i] = rho * a[i - 1] + w[i]
    for name, (df_, mse_) in one_run(a).items():
        res_ar[name]["df"].append(df_)
        res_ar[name]["mse"].append(mse_)

med = {name: float(np.median(res[name]["df"])) for name in res}
mse = {name: float(np.mean(res[name]["mse"])) for name in res}
med_ar = {name: float(np.median(res_ar[name]["df"])) for name in res_ar}
mse_ar = {name: float(np.mean(res_ar[name]["mse"])) for name in res_ar}
wild = {name: float(np.mean(np.array(res[name]["df"]) > 15)) for name in res}
q90 = {name: float(np.quantile(res[name]["df"], 0.9)) for name in res}
assert med_ar["GCV"] > 2.5 * med["GCV"]
assert mse_ar["GCV"] > 3 * mse["GCV"]
assert med["AICc"] < med["GCV"] + 1e-9

# an AR(1) data set to draw
a = np.empty(n)
a[0] = sigma * rng.normal()
for i in range(1, n):
    a[i] = rho * a[i - 1] + sigma * np.sqrt(1 - rho**2) * rng.normal()
y_ar = f_true(x) + a
lam_ar = select(B, P, y_ar, lams, 2)
fit_ar = smoother(B, P, lam_ar) @ y_ar
df_ar = np.trace(smoother(B, P, lam_ar))

gen = Generated("ch43", "choose_lambda")
gen.int("n", n)
gen.int("K", K)
gen.int("reps", reps)
gen.num("rho", rho, 1)
gen.num("sigma", sigma, 2)
gen.num("shortcut_gap", shortcut_gap, 1, sci=True)
for name in ("CV", "GCV", "AICc"):
    key = name.lower()
    gen.num(f"df_{key}", df_pick[name], 2)
    gen.num(f"med_{key}", med[name], 2)
    gen.num(f"mse_{key}", mse[name], 5)
    gen.num(f"medar_{key}", med_ar[name], 2)
    gen.num(f"msear_{key}", mse_ar[name], 5)
    gen.num(f"wild_{key}", 100 * wild[name], 1)
gen.num("df_ar", df_ar, 2)
gen.write()

# ---- figure -----------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4))
ax = axes[0]
ax.plot(dfs, cvs / cvs.min(), color=COLORS["accent"], label="CV")
ax.plot(dfs, gcvs / gcvs.min(), color=COLORS["second"], linestyle="--", label="GCV")
ax.plot(dfs, np.exp(aiccs) / np.exp(aiccs).min(), color=COLORS["third"], linestyle="-.", label="AICc")
for name, col in (("CV", COLORS["accent"]), ("GCV", COLORS["second"]), ("AICc", COLORS["third"])):
    ax.axvline(df_pick[name], color=col, linewidth=0.7, linestyle=":")
ax.set_xlabel("effective degrees of freedom")
ax.set_ylabel("criterion / its minimum")
ax.set_xlim(2, 30)
ax.set_ylim(0.995, 1.10)
ax.set_title("(a) three criteria, one data set")
ax.legend(frameon=False, loc="upper center", handlelength=1.6)
ax = axes[1]
ax.scatter(x, y_ar, s=6, color=COLORS["muted"], alpha=0.8, linewidths=0)
ax.plot(x, f_true(x), color=COLORS["ink"], linewidth=1.0, label="true $f$")
ax.plot(x, fit_ar, color=COLORS["second"], linewidth=1.1, label=f"GCV fit, df $={df_ar:.0f}$")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_ylim(-0.8, 2.9)
ax.set_title("(b) AR(1) errors, $\\rho=0.7$")
ax.legend(frameon=False, loc="upper left", handlelength=1.4)
fig.tight_layout()
fig.savefig(figure_path("ch43", "smoothing_criteria"))
