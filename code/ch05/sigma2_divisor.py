"""Chapter 5, Section 6: why the residual sum of squares is divided by n - p.

(i) With n = 30 and p = 1, ..., 25 columns, E(SSE / n) = sigma^2 (n - p) / n, while
    SSE / (n - p) is unbiased whatever p is (Monte Carlo, normal errors).
(ii) The variance of s^2 depends on the fourth moment of the errors (Chapter 2); here
    for centred exponential errors, which are skewed:
    Var(SSE) = (mu4 - 3 sigma^4) sum_i (1 - h_ii)^2 + 2 sigma^4 (n - p).
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch05", "sigma2_divisor", prefix="s2")

# <<divisor>>
rng = np.random.default_rng(506)
n, sigma2, reps = 30, 1.0, 20_000
Xfull = np.column_stack([np.ones(n), rng.normal(size=(n, 24))])

ps = np.arange(1, 26)
mean_n, mean_np = [], []
for p in ps:
    X = Xfull[:, :p]
    R = np.eye(n) - X @ np.linalg.solve(X.T @ X, X.T)     # I - H
    E = rng.standard_normal((reps, n))                      # errors, sigma^2 = 1
    SSE = np.einsum("ri,ij,rj->r", E, R, E)                 # e^T (I - H) e
    mean_n.append(SSE.mean() / n)
    mean_np.append(SSE.mean() / (n - p))
print("p     E(SSE/n)  E(SSE/(n-p))")
for p in (1, 5, 15, 25):
    print(f"{p:2d}    {mean_n[p - 1]:.3f}     {mean_np[p - 1]:.3f}")
# <</divisor>>

mean_n, mean_np = np.array(mean_n), np.array(mean_np)
assert np.all(np.abs(mean_np - 1) < 0.03)
assert np.allclose(mean_n, (n - ps) / n, atol=0.02)

# ---- (ii) variance of SSE with skewed errors ----------------------------------
# <<variance>>
p = 6
X = Xfull[:, :p]
H = X @ np.linalg.solve(X.T @ X, X.T)
R = np.eye(n) - H
mu4 = 9.0                                    # fourth moment of a centred Exp(1) variable
var_normal = 2 * (n - p)
var_exp = (mu4 - 3) * np.sum(np.diag(R) ** 2) + 2 * (n - p)
reps2 = 200_000
E = rng.exponential(size=(reps2, n)) - 1.0
SSE = np.einsum("ri,ij,rj->r", E, R, E)
print(f"Var(SSE): normal errors {var_normal:.2f}, exponential errors {var_exp:.2f},"
      f" simulated {SSE.var():.2f}")
# <</variance>>
assert abs(SSE.mean() / (n - p) - 1) < 0.01
assert abs(SSE.var() / var_exp - 1) < 0.03

gen.int("n", n)
gen.text("reps", f"{reps:,}")
gen.text("reps2", f"{reps2:,}")
for p_ in (5, 15, 25):
    gen.num(f"n{p_}", mean_n[p_ - 1], 3)
    gen.num(f"np{p_}", mean_np[p_ - 1], 3)
gen.int("p", p)
gen.num("varnormal", var_normal, 2)
gen.num("varexp", var_exp, 2)
gen.num("varsim", SSE.var(), 2)
gen.num("sumdiag", np.sum(np.diag(R) ** 2), 3)
gen.num("sds2normal", np.sqrt(var_normal) / (n - p), 3)
gen.num("sds2exp", np.sqrt(var_exp) / (n - p), 3)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 2.5))
ax.plot(ps, (n - ps) / n, color=COLORS["muted"], linewidth=0.8, linestyle="--",
        label=r"$(n-p)/n$")
ax.plot(ps, mean_n, "o", markersize=3, color=COLORS["second"], label=r"SSE$/n$ (simulated mean)")
ax.plot(ps, mean_np, "s", markersize=3, color=COLORS["accent"],
        label=r"SSE$/(n-p)$ (simulated mean)")
ax.axhline(1.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.set_xlabel(r"number of columns $p$ (with $n=30$)")
ax.set_ylabel(r"mean of estimate / $\sigma^2$")
ax.set_ylim(0, 1.15)
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(figure_path("ch05", "sigma2_divisor"))
