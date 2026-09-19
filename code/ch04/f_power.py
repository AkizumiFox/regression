"""Chapter 4, Section 4: the F statistic for nested models and its power.

(1) Power P(F(r, s, gamma) > F_{0.95}(r, s)) as a function of gamma for several r,
    checking monotonicity in gamma (increasing) and in r (decreasing).
(2) Testing curvature: full model 1, x, x^2 against the straight line 1, x on an
    equally spaced design. The noncentrality is ||(M - M0) X beta||^2 / sigma^2 and the
    power from the noncentral F is checked by simulation.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(1104)
alpha, s = 0.05, 12


def power(gamma, q, s, alpha=0.05):
    crit = stats.f.ppf(1 - alpha, q, s)
    return stats.ncf.sf(crit, q, s, np.maximum(gamma, 1e-12))


gammas = np.linspace(0, 30, 301)
gamma_ref, rs = 10.0, (1, 3, 6)                        # r = numerator degrees of freedom
curves = {q: power(gammas, q, s) for q in rs}
for q, pw in curves.items():
    assert np.isclose(pw[0], alpha, atol=1e-6) and np.all(np.diff(pw) > 0)
assert np.all(curves[1][1:] > curves[3][1:]) and np.all(curves[3][1:] > curves[6][1:])

# <<curvature>>
n, sigma, beta0, beta1, beta2, reps = 15, 1.0, 1.0, 0.5, 8.0, 100_000
x = np.linspace(0, 1, n)
X0 = np.column_stack([np.ones(n), x])                  # straight line
X = np.column_stack([X0, x**2])                        # adds curvature


def proj(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T


M, M0 = proj(X), proj(X0)
mean = X @ np.array([beta0, beta1, beta2])
gamma = mean @ (M - M0) @ mean / sigma**2              # noncentrality
q, df_resid = 1, n - 3
crit = stats.f.ppf(0.95, q, df_resid)
exact_power = stats.ncf.sf(crit, q, df_resid, gamma)

Y = mean + sigma * rng.standard_normal((reps, n))
num = np.einsum("ij,jk,ik->i", Y, M - M0, Y) / q
den = np.einsum("ij,jk,ik->i", Y, np.eye(n) - M, Y) / df_resid
print(gamma, exact_power, np.mean(num / den > crit))
# <</curvature>>

sim_power = np.mean(num / den > crit)
assert df_resid == s and np.isclose(power(gamma, 1, s), exact_power)  # point lies on the q=1 curve
assert abs(sim_power - exact_power) < 5 * np.sqrt(exact_power * (1 - exact_power) / len(Y))
assert np.isclose(np.trace(M - M0), 1) and np.allclose((M - M0) @ (np.eye(n) - M), 0)
# n needed for 80% power at the same curvature (equally spaced on [0, 1])
n80 = None
for m in range(4, 200):
    xm = np.linspace(0, 1, m)
    Xm0 = np.column_stack([np.ones(m), xm])
    Xm = np.column_stack([Xm0, xm**2])
    mm = Xm @ np.array([beta0, beta1, beta2])
    g = mm @ (proj(Xm) - proj(Xm0)) @ mm / sigma**2
    if stats.ncf.sf(stats.f.ppf(0.95, 1, m - 3), 1, m - 3, g) >= 0.8:
        n80, g80 = m, g
        break
assert n80 is not None and n80 > n

gen = Generated("ch04", "f_power", prefix="fp")
gen.int("s", s)
gen.num("gammaref", gamma_ref, 0)
gen.int("rone", rs[0])
gen.int("rthree", rs[1])
gen.int("rsix", rs[2])
gen.num("pone", power(gamma_ref, rs[0], s), 3)
gen.num("pthree", power(gamma_ref, rs[1], s), 3)
gen.num("psix", power(gamma_ref, rs[2], s), 3)
gen.text("reps", f"{reps:,}".replace(",", "{,}"))
gen.num("betazero", beta0, 0)
gen.num("betaone", beta1, 1)
gen.int("n", n)
gen.num("beta", beta2, 1)
gen.num("gamma", gamma, 3)
gen.num("crit", crit, 3)
gen.num("power", exact_power, 3)
gen.num("simpower", sim_power, 3)
gen.int("neighty", n80)
gen.num("geighty", g80, 2)
gen.write()

# ---- figure ---------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.4))
for (q_, pw), col in zip(curves.items(), [COLORS["accent"], COLORS["third"], COLORS["second"]]):
    ax.plot(gammas, pw, color=col, label=f"$r={q_}$")
ax.axhline(alpha, color=COLORS["muted"], linewidth=0.6, linestyle=":")
ax.plot([gamma], [exact_power], "o", color=COLORS["ink"], markersize=3.5)
ax.annotate("curvature test", (gamma, exact_power), textcoords="offset points",
            xytext=(6, -12), fontsize=8, color=COLORS["ink"])
ax.set_xlabel(r"noncentrality $\gamma$")
ax.set_ylabel("power")
ax.set_xlim(0, 30)
ax.set_ylim(0, 1)
ax.legend(frameon=False, loc="lower right", title=rf"$F(r,{s},\gamma)$", title_fontsize=8)
fig.tight_layout()
fig.savefig(figure_path("ch04", "f_power"))
