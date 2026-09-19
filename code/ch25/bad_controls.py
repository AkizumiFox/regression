"""Chapter 25, Section 5: bad controls in linear structural causal models.

Three models in which adjusting for a third variable C biases the least squares estimate
of the total effect of X on Y: C a mediator, C a common effect (collider) of X and Y, and
C a pre-treatment collider (M-bias). The bias formulas of the section are compared with
simulation, and the sampling distributions are plotted.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<models>>
import numpy as np

rng = np.random.default_rng(2505)

def mediator(n):          # X -> C -> Y and X -> Y; total effect 0.5 + 1.2 * 1.0 = 1.7
    x = rng.normal(size=n); c = 1.2 * x + rng.normal(size=n)
    return x, c, 0.5 * x + 1.0 * c + rng.normal(size=n)

def collider(n):          # X -> Y, and X -> C <- Y; total effect 1
    x = rng.normal(size=n); y = 1.0 * x + rng.normal(size=n)
    return x, 0.5 * x + 1.0 * y + rng.normal(size=n), y

def m_bias(n):            # X <- H1 -> C <- H2 -> Y, and X -> Y; total effect 1
    h1, h2 = rng.normal(size=(2, n))
    x = 1.5 * h1 + rng.normal(size=n)
    c = 1.5 * h1 + 1.5 * h2 + rng.normal(size=n)
    return x, c, 1.0 * x + 1.5 * h2 + rng.normal(size=n)

def slopes(x, c, y):
    """Coefficient of x without and with adjustment for c."""
    one = np.ones(len(x))
    b_short = np.linalg.lstsq(np.column_stack([one, x]), y, rcond=None)[0][1]
    b_long = np.linalg.lstsq(np.column_stack([one, x, c]), y, rcond=None)[0][1]
    return b_short, b_long

# <</models>>

# <<quick>>
n, reps = 400, 500
for name, model in [("mediator", mediator), ("collider", collider), ("M-bias", m_bias)]:
    est = np.array([slopes(*model(n)) for _ in range(reps)])
    print(f"{name:9s} mean without C {est[:, 0].mean():.3f}   mean with C {est[:, 1].mean():.3f}")
# <</quick>>

reps = 4000
results = {}
for name, model in [("mediator", mediator), ("collider", collider), ("M-bias", m_bias)]:
    results[name] = np.array([slopes(*model(n)) for _ in range(reps)])

# <<selection>>
x, c, y = collider(200_000)
keep = c > 1                                          # only units with a large C are observed
slope_all = np.polyfit(x, y, 1)[0]
slope_sel = np.polyfit(x[keep], y[keep], 1)[0]
print(f"slope of Y on X: all units {slope_all:.3f}, units with C > 1 {slope_sel:.3f}")
# <</selection>>

# formulas of the theorem
theta, delta, gamma = 0.5, 1.2, 1.0
tau_med = theta + gamma * delta
tau_col, a, b = 1.0, 0.5, 1.0
bias_col = -(a + b * tau_col) * b / (b ** 2 + 1)
aa = bb = c1 = c2 = 1.5
D = (aa ** 2 + 1) * (c1 ** 2 + c2 ** 2 + 1) - (aa * c1) ** 2
bias_m = -aa * bb * c1 * c2 / D
tol = 4 * 0.1 / np.sqrt(reps)
assert abs(results["mediator"][:, 0].mean() - tau_med) < tol
assert abs(results["mediator"][:, 1].mean() - theta) < tol
assert abs(results["collider"][:, 0].mean() - tau_col) < tol
assert abs(results["collider"][:, 1].mean() - (tau_col + bias_col)) < tol
assert abs(results["M-bias"][:, 0].mean() - 1.0) < tol
assert abs(results["M-bias"][:, 1].mean() - (1.0 + bias_m)) < tol
assert abs(slope_all - 1) < 0.02 and slope_sel < 0.8

gen = Generated("ch25", "bad_controls")
gen.num("taumed", tau_med, 1)
gen.num("biascol", bias_col, 2)
gen.num("colcoef", tau_col + bias_col, 2)
gen.num("D", D, 4)
gen.num("biasm", bias_m, 3)
gen.num("mcoef", 1 + bias_m, 3)
for name, tag in [("mediator", "med"), ("collider", "col"), ("M-bias", "m")]:
    gen.num(f"{tag}short", results[name][:, 0].mean(), 3)
    gen.num(f"{tag}long", results[name][:, 1].mean(), 3)
    gen.num(f"{tag}sdshort", results[name][:, 0].std(), 3)
    gen.num(f"{tag}sdlong", results[name][:, 1].std(), 3)
gen.num("slopeall", slope_all, 3)
gen.num("slopesel", slope_sel, 3)
gen.num("fracsel", keep.mean(), 2)
gen.int("n", n)
gen.int("reps", reps)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.2), sharey=False)
truth = {"mediator": tau_med, "collider": tau_col, "M-bias": 1.0}
for ax, (name, r) in zip(axes, results.items()):
    lo, hi = r.min(), r.max()
    bins = np.linspace(lo, hi, 60)
    ax.hist(r[:, 0], bins=bins, color=COLORS["accent"], alpha=0.65, label="without C")
    ax.hist(r[:, 1], bins=bins, color=COLORS["second"], alpha=0.65, label="with C")
    ax.axvline(truth[name], color=COLORS["ink"], linestyle="--", linewidth=0.9)
    panel = {"mediator": "mediator", "collider": "common effect", "M-bias": "M-bias"}[name]
    ax.set_title(f"({'abc'[list(results).index(name)]}) {panel}")
    ax.set_yticks([])
    ax.set_xlabel("coefficient of X")
axes[0].legend(frameon=False, loc="upper center")
fig.tight_layout()
fig.savefig(figure_path("ch25", "bad_controls"))
