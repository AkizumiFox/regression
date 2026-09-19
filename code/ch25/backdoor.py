"""Chapter 25, Section 4: adjustment sets in a linear structural causal model.

An eight-variable linear SCM with unit-variance disturbances. For several candidate
adjustment sets Z, the population coefficient of X in the best linear predictor of Y from
(X, Z) is computed exactly from the covariance matrix, with its large-sample standard
deviation sigma_{Y.XZ} / sigma_{X.Z} / sqrt(n). A simulation checks both.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<model>>
import numpy as np

names = ["W1", "W3", "I", "P", "W2", "X", "M", "Y"]      # a topological order
k = {v: j for j, v in enumerate(names)}
B = np.zeros((8, 8))                                      # B[j, i]: coefficient of i in j's equation
for child, parent, coef in [("W2", "W1", 0.8), ("X", "W1", 0.7), ("X", "W3", 0.6), ("X", "I", 0.9),
                            ("M", "X", 0.8), ("Y", "X", 0.5), ("Y", "M", 0.6), ("Y", "W2", 0.7),
                            ("Y", "W3", 0.8), ("Y", "P", 1.0)]:
    B[k[child], k[parent]] = coef
A = np.linalg.inv(np.eye(8) - B)                          # V = A U
Sigma = A @ A.T                                           # unit disturbance variances
tau = A[k["Y"], k["X"]]                                   # total effect: sum over directed paths

def population_fit(Z):
    """Coefficient of X in L(Y | X, Z) and the large-sample sd factor sigma_{Y.XZ}/sigma_{X.Z}."""
    R = [k["X"]] + [k[z] for z in Z]
    coef = np.linalg.solve(Sigma[np.ix_(R, R)], Sigma[R, k["Y"]])
    res_y = Sigma[k["Y"], k["Y"]] - Sigma[k["Y"], R] @ coef
    Zi = [k[z] for z in Z]
    res_x = Sigma[k["X"], k["X"]] - (Sigma[k["X"], Zi] @ np.linalg.solve(Sigma[np.ix_(Zi, Zi)], Sigma[Zi, k["X"]])
                                     if Z else 0.0)
    return coef[0], np.sqrt(res_y / res_x)

sets = [[], ["W3"], ["W1", "W3"], ["W2", "W3"], ["W2", "W3", "P"], ["W1", "W3", "I"],
        ["W1", "W3", "M"], ["W1", "W2", "W3", "I", "P"]]
print(f"total effect of X on Y: {tau:.3f}")
for Z in sets:
    c, f = population_fit(Z)
    print(f"Z = {{{', '.join(Z)}}}".ljust(26) + f"coefficient {c:6.3f}   sd factor {f:5.3f}")
# <</model>>

valid = [False, False, True, True, True, True, False, True]
fits = [population_fit(Z) for Z in sets]
for (c, f), ok in zip(fits, valid):
    assert np.isclose(c, tau) == ok
assert np.isclose(tau, 0.5 + 0.8 * 0.6)
# precision: an outcome predictor helps, a pure cause of X hurts
f_base = fits[3][1]
assert fits[4][1] < f_base and fits[5][1] > fits[2][1]
# coefficient of W3 given (X, W1, W3) is its direct effect, not its total effect
R = [k["X"], k["W1"], k["W3"]]
coef_w3 = np.linalg.solve(Sigma[np.ix_(R, R)], Sigma[R, k["Y"]])[2]
total_w3 = A[k["Y"], k["W3"]]
assert np.isclose(coef_w3, 0.8) and np.isclose(total_w3, 0.8 + 0.6 * tau)
# the partial covariances behind the precision proposition: sigma_{Y I . X W1 W3} = 0, sigma_{X P . W2 W3} = 0
def partial_cov(a, b, Z):
    Zi = [k[z] for z in Z]
    return Sigma[k[a], k[b]] - Sigma[k[a], Zi] @ np.linalg.solve(Sigma[np.ix_(Zi, Zi)], Sigma[Zi, k[b]])
assert abs(partial_cov("Y", "I", ["X", "W1", "W3"])) < 1e-12
assert abs(partial_cov("X", "P", ["W2", "W3"])) < 1e-12
# the three sets of Exercise A1 of Section 25.4 (none is valid)
extra = {tag: population_fit(Z)[0] for tag, Z in
         [("wonewtwo", ["W1", "W2"]), ("wthreep", ["W3", "P"]), ("wtwowthreem", ["W2", "W3", "M"])]}
assert np.isclose(extra["wtwowthreem"], 0.5)
assert extra["wonewtwo"] > tau and extra["wthreep"] > tau

# ---- simulation ---------------------------------------------------------------
rng = np.random.default_rng(2504)
n, reps = 500, 4000
est = np.empty((reps, len(sets)))
for r in range(reps):
    V = rng.normal(size=(n, 8)) @ A.T
    for s, Z in enumerate(sets):
        Xm = np.column_stack([np.ones(n), V[:, k["X"]]] + [V[:, k[z]] for z in Z])
        est[r, s] = np.linalg.lstsq(Xm, V[:, k["Y"]], rcond=None)[0][1]
for s, (c, f) in enumerate(fits):
    assert abs(est[:, s].mean() - c) < 4 * f / np.sqrt(n) / np.sqrt(reps) + 0.002
    assert abs(est[:, s].std() / (f / np.sqrt(n)) - 1) < 0.06

gen = Generated("ch25", "backdoor")
gen.num("tau", tau, 2)
tags = ["empty", "wthree", "wonewthree", "wtwowthree", "withp", "withi", "withm", "all"]
for tag, (c, f), s in zip(tags, fits, range(len(sets))):
    gen.num(f"c{tag}", c, 3)
    gen.num(f"f{tag}", f, 2)
    gen.num(f"sim{tag}", est[:, s].mean(), 3)
    gen.num(f"sd{tag}", est[:, s].std(), 3)
gen.num("coefwthree", coef_w3, 2)
for tag, c in extra.items():
    gen.num(f"x{tag}", c, 3)
gen.num("totalwthree", total_w3, 3)
gen.int("n", n)
gen.int("reps", reps)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(5.0, 2.9))
labels = ["{" + ", ".join(Z) + "}" if Z else "{ } (none)" for Z in sets]
ypos = np.arange(len(sets))[::-1]
for yv, (c, f), ok, s in zip(ypos, fits, valid, range(len(sets))):
    col = COLORS["accent"] if ok else COLORS["second"]
    half = 1.96 * f / np.sqrt(n)
    ax.plot([c - half, c + half], [yv, yv], color=col, linewidth=2.2, alpha=0.55, solid_capstyle="butt")
    ax.plot(c, yv, "o", color=col, markersize=4)
    ax.plot(est[:, s].mean(), yv, "x", color=COLORS["ink"], markersize=4)
ax.axvline(tau, color=COLORS["third"], linestyle="--", linewidth=0.9)
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_xlabel(r"coefficient of $X$ given $Z$  (population value, $\pm1.96$ sd at $n=500$)")
fig.tight_layout()
fig.savefig(figure_path("ch25", "adjustment_sets"))
