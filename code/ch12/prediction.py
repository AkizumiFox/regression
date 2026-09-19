"""Chapter 12, Section 3: prediction intervals, leverage of a new point, and the role of normality.

Part 1 uses the state murder-rate regression of Sections 12.1-12.2 (statsmodels.datasets.statecrime,
public domain). Part 2 is a fixed-seed simulation of the coverage of intervals for the mean response
and of prediction intervals when the errors are not normal.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

# <<predict>>
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
C = np.linalg.inv(X.T @ X)
beta_hat = C @ X.T @ y
s = np.sqrt(np.sum((y - X @ beta_hat) ** 2) / (n - p))
q = stats.t.ppf(0.975, n - p)
h = np.einsum("ij,jk,ik->i", X, C, X)            # leverages of the 50 states


def intervals(x0):
    """Leverage h0, 95% interval for the mean response, 95% prediction interval."""
    h0 = x0 @ C @ x0
    fit = x0 @ beta_hat
    ci = (fit - q * s * np.sqrt(h0), fit + q * s * np.sqrt(h0))
    pi = (fit - q * s * np.sqrt(1 + h0), fit + q * s * np.sqrt(1 + h0))
    return h0, fit, ci, pi


x_centre = X.mean(axis=0)                        # the average state
x_hidden = np.array([1.0, 20.0, 19.0, 60.0])     # each coordinate inside its observed range
for name, x0 in [("centre", x_centre), ("hidden", x_hidden)]:
    h0, fit, ci, pi = intervals(x0)
    print(f"{name}: h0 = {h0:.4f}, fit {fit:.3f}, mean ({ci[0]:.3f}, {ci[1]:.3f}), "
          f"new state ({pi[0]:.3f}, {pi[1]:.3f})")
print(f"largest leverage in the data: {h.max():.4f}")
# <</predict>>

h_c, fit_c, ci_c, pi_c = intervals(x_centre)
h_h, fit_h, ci_h, pi_h = intervals(x_hidden)
assert np.isclose(h_c, 1 / n)                    # at the centroid the leverage is 1/n
assert np.isclose((pi_c[1] - pi_c[0]) / (ci_c[1] - ci_c[0]), np.sqrt((1 + h_c) / h_c))
lo_rng, hi_rng = X.min(axis=0), X.max(axis=0)
assert np.all((x_hidden[1:] >= lo_rng[1:]) & (x_hidden[1:] <= hi_rng[1:]))
assert h_h > h.max()                             # hidden extrapolation
# the leverage of a new point via the centred form (Chapter 6)
Xc = X[:, 1:] - X[:, 1:].mean(axis=0)
d0 = x_hidden[1:] - X[:, 1:].mean(axis=0)
assert np.isclose(h_h, 1 / n + d0 @ np.linalg.solve(Xc.T @ Xc, d0))
# enlarging the model never decreases the variance of a prediction
X_small = X[:, :3]
C_small = np.linalg.inv(X_small.T @ X_small)
rng = np.random.default_rng(1203)
for _ in range(200):
    x0 = np.concatenate([[1.0], rng.uniform(lo_rng[1:], hi_rng[1:])])
    assert x0 @ C @ x0 >= x0[:3] @ C_small @ x0[:3] - 1e-12
h_small = x_hidden[:3] @ C_small @ x_hidden[:3]

# ---- simulation: coverage when the errors are not normal ---------------------------------------
# <<coverage>>
laws = {
    "normal": lambda size: rng.normal(size=size),
    "exponential": lambda size: rng.exponential(size=size) - 1.0,                  # skewed
    "uniform": lambda size: rng.uniform(-np.sqrt(3), np.sqrt(3), size=size),       # short tails
}                                                                                  # all: mean 0, variance 1
x0_sim = 0.8


def simulate(law, n_sim, reps=20_000):
    x = np.linspace(0.0, 1.0, n_sim)
    Xs = np.column_stack([np.ones(n_sim), x])
    Hs = np.linalg.solve(Xs.T @ Xs, Xs.T)                        # beta_hat = Hs y
    z0 = np.array([1.0, x0_sim])
    h0 = z0 @ np.linalg.solve(Xs.T @ Xs, z0)
    E = laws[law]((reps, n_sim))                                  # true beta = 0, sigma = 1
    fit = E @ (Hs.T @ z0)
    s_sim = np.sqrt(np.sum((E - (E @ Hs.T) @ Xs.T) ** 2, axis=1) / (n_sim - 2))
    new = laws[law](reps)                                         # the new observation's error
    out = {}
    for level in (0.80, 0.95):
        t = stats.t.ppf(0.5 + level / 2, n_sim - 2)
        out[("mean", level)] = np.mean(np.abs(fit) <= t * s_sim * np.sqrt(h0))
        out[("new", level)] = np.mean(np.abs(new - fit) <= t * s_sim * np.sqrt(1 + h0))
        out[("above", level)] = np.mean(new - fit > t * s_sim * np.sqrt(1 + h0))
        out[("below", level)] = np.mean(new - fit < -t * s_sim * np.sqrt(1 + h0))
    return out


table = {(law, m): simulate(law, m) for law in laws for m in (20, 400)}
for (law, m), out in table.items():
    print(f"{law:12s} n = {m:3d}: " + ", ".join(f"{k[0]} {k[1]:.2f}: {v:.3f}" for k, v in out.items()))
# <</coverage>>

# the limits of the prediction-interval coverage as n grows: Pr(|eps| <= z sigma)
def limit(law, level):
    z = stats.norm.ppf(0.5 + level / 2)
    if law == "normal":
        return level
    if law == "exponential":                     # eps = E - 1, E ~ Exp(1)
        return stats.expon.cdf(1 + z) - stats.expon.cdf(1 - z)
    return min(1.0, z / np.sqrt(3))


lims = {(law, lev): limit(law, lev) for law in laws for lev in (0.80, 0.95)}
reps_sim = 20_000
tol = 4 * np.sqrt(0.25 / reps_sim)
for (law, m), out in table.items():
    for level in (0.80, 0.95):
        if law == "normal":                      # exact for every n
            assert abs(out[("mean", level)] - level) < tol and abs(out[("new", level)] - level) < tol
        if m == 400:                             # the mean interval recovers; the prediction interval does not
            assert abs(out[("mean", level)] - level) < 0.012
            assert abs(out[("new", level)] - lims[(law, level)]) < 0.012
assert lims[("exponential", 0.80)] > 0.88 and lims[("uniform", 0.80)] < 0.75
# with exponential errors the 95% prediction interval misses on one side only, in the limit
z95 = stats.norm.ppf(0.975)
assert 1 - z95 < 0 and np.isclose(lims[("exponential", 0.95)], 1 - np.exp(-(1 + z95)))

gen = Generated("ch12", "prediction", prefix="pred")
gen.int("n", n)
gen.num("hc", h_c, 4)
gen.num("fitc", fit_c, 3)
gen.num("cic:lo", ci_c[0], 3)
gen.num("cic:hi", ci_c[1], 3)
gen.num("pic:lo", pi_c[0], 3)
gen.num("pic:hi", pi_c[1], 3)
gen.num("hh", h_h, 3)
gen.num("fith", fit_h, 3)
gen.num("cih:lo", ci_h[0], 3)
gen.num("cih:hi", ci_h[1], 3)
gen.num("pih:lo", pi_h[0], 3)
gen.num("pih:hi", pi_h[1], 3)
gen.num("hmax", h.max(), 3)
gen.num("lenratio", (pi_c[1] - pi_c[0]) / (ci_c[1] - ci_c[0]), 2)
gen.num("pilen", pi_c[1] - pi_c[0], 3)
gen.num("pov:min", lo_rng[1], 1)
gen.num("pov:max", hi_rng[1], 1)
gen.num("sing:min", lo_rng[2], 1)
gen.num("sing:max", hi_rng[2], 1)
gen.text("hmaxstate", data.index[h.argmax()])
gen.num("hsmall", h_small, 3)
gen.num("widthratio", np.sqrt(1 + h_h) / np.sqrt(1 + h_c), 3)
gen.num("ciratio", np.sqrt(h_h / h_c), 2)
for law in laws:
    for m in (20, 400):
        for kind in ("mean", "new"):
            for level in (0.80, 0.95):
                gen.num(f"cov:{law}:{m}:{kind}:{int(level * 100)}", table[(law, m)][(kind, level)], 3)
    for level in (0.80, 0.95):
        gen.num(f"lim:{law}:{int(level * 100)}", lims[(law, level)], 3)
gen.num("exp:above", table[("exponential", 400)][("above", 0.95)], 3)
gen.num("exp:below", table[("exponential", 400)][("below", 0.95)], 4)
gen.int("reps", reps_sim)
gen.write()

# ---- figure: contours of the leverage of a new state ---------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.2, 3.3))
g1 = np.linspace(7, 23, 200)
g2 = np.linspace(16, 34, 200)
P, S = np.meshgrid(g1, g2)
Z = np.stack([np.ones_like(P), P, S, np.full_like(P, X[:, 3].mean())], axis=-1)
H0 = np.einsum("...j,jk,...k->...", Z, C, Z)
cs = ax.contour(P, S, H0, levels=[0.05, 0.1, 0.2, 0.3], colors=COLORS["accent"], linewidths=0.7)
ax.clabel(cs, fontsize=7, fmt="%.2f")
ax.contour(P, S, H0, levels=[h.max()], colors=COLORS["second"], linewidths=1.0, linestyles="--")
ax.scatter(X[:, 1], X[:, 2], s=9, color=COLORS["ink"], linewidths=0, alpha=0.8)
ax.plot(*x_hidden[1:3], "^", color=COLORS["third"], markersize=6)
ax.annotate("new state", x_hidden[1:3], textcoords="offset points", xytext=(-8, -24), fontsize=8,
            arrowprops=dict(arrowstyle="-", color=COLORS["muted"], linewidth=0.6))
ax.add_patch(plt.Rectangle((lo_rng[1], lo_rng[2]), hi_rng[1] - lo_rng[1], hi_rng[2] - lo_rng[2],
                           fill=False, edgecolor=COLORS["muted"], linewidth=0.6, linestyle=":"))
ax.set_xlabel("poverty (%)")
ax.set_ylabel("single-parent households (%)")
fig.tight_layout()
fig.savefig(figure_path("ch12", "prediction_leverage"))
