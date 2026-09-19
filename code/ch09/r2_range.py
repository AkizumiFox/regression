"""Chapter 9, Section 5: R^2 depends on the spread of the regressor.

The same straight line and the same errors give very different R^2 when the
regressor is observed over a narrow or a wide range."""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<range>>
rng = np.random.default_rng(9)
n, beta0, beta1, sigma = 50, 2.0, 0.5, 1.0
eps = sigma * rng.normal(size=n)                        # the same errors in both designs


def fit(x):
    y = beta0 + beta1 * x + eps
    X = np.column_stack([np.ones(n), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    sse = np.sum((y - X @ b) ** 2)
    r2 = 1 - sse / np.sum((y - y.mean()) ** 2)
    return y, b, sse / (n - 2), r2


for label, x in [("narrow", np.linspace(4, 6, n)), ("wide", np.linspace(0, 10, n))]:
    y, b, s2, r2 = fit(x)
    Sxx = np.sum((x - x.mean()) ** 2)
    expected = beta1 ** 2 * Sxx / (beta1 ** 2 * Sxx + (n - 1) * sigma ** 2)
    print(f"{label:6s}: slope {b[1]:.3f}, s^2 {s2:.3f}, R^2 {r2:.3f}, design value {expected:.3f}")
# <</range>>

gen = Generated("ch09", "r2_range", prefix="rng")
out = {}
for label, x in [("narrow", np.linspace(4, 6, n)), ("wide", np.linspace(0, 10, n))]:
    y, b, s2, r2 = fit(x)
    Sxx = np.sum((x - x.mean()) ** 2)
    expected = beta1 ** 2 * Sxx / (beta1 ** 2 * Sxx + (n - 1) * sigma ** 2)
    out[label] = (x, y, b, r2, s2)
    gen.num(f"{label}:slope", b[1], 3)
    gen.num(f"{label}:s2", s2, 3)
    gen.num(f"{label}:r2", r2, 3)
    gen.num(f"{label}:design", expected, 3)
    gen.num(f"{label}:Sxx", Sxx, 1)
# both designs span the same C(X) = span(1, x), so the residuals and the residual
# mean squares are identical; R^2 is not
assert np.isclose(out["narrow"][4], out["wide"][4])
assert abs(out["narrow"][3] - out["wide"][3]) > 0.4
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.3), sharey=True)
for ax, label in zip(axes, ["narrow", "wide"]):
    x, y, b, r2, _ = out[label]
    ax.scatter(x, y, s=8, color=COLORS["accent"], alpha=0.8, linewidths=0)
    xs = np.array([0, 10])
    ax.plot(xs, beta0 + beta1 * xs, color=COLORS["muted"], lw=0.8, ls="--")
    ax.plot(xs, b[0] + b[1] * xs, color=COLORS["second"])
    ax.set_xlim(-0.3, 10.3)
    ax.set_xlabel("x")
    ax.set_title(f"({'a' if label == 'narrow' else 'b'}) {label} range: $R^2={r2:.2f}$")
axes[0].set_ylabel("y")
fig.tight_layout()
fig.savefig(figure_path("ch09", "r2_range"))
