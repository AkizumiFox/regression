"""Chapter 29, Section 2: selection in an orthonormal design is hard thresholding.

Risk of the estimator Z 1{|Z| > c} of theta, Z ~ N(theta, 1), for the thresholds chosen by
least squares (c = 0), C_p and AIC (c = sqrt 2), a 5% pretest (c = 1.96) and BIC with
n = 1000 (c = sqrt(log 1000)). The closed form is checked by simulation.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch29", "thresholding_risk", prefix="thr")

# <<risk>>
import numpy as np
from scipy.stats import norm


def hard_risk(theta, c):
    """E(Z 1{|Z| > c} - theta)^2 for Z ~ N(theta, 1)."""
    D = norm.cdf(c - theta) - norm.cdf(-c - theta)        # Pr(|Z| <= c): coordinate dropped
    return theta**2 * D + 1 - D + (c - theta) * norm.pdf(c - theta) + (c + theta) * norm.pdf(c + theta)


thetas = np.linspace(0, 6, 601)
for label, c in [("C_p / AIC", np.sqrt(2)), ("5% pretest", 1.96), ("BIC, n = 1000", np.sqrt(np.log(1000)))]:
    r = hard_risk(thetas, c)
    print(f"{label:14s} c = {c:.3f}: risk at 0 = {r[0]:.3f}, "
          f"maximum {r.max():.3f} at theta = {thetas[r.argmax()]:.2f}")
# <</risk>>

rng = np.random.default_rng(2904)
Z = rng.normal(size=1_000_000)
for c in (np.sqrt(2), 1.96, 2.6):
    for th in (0.0, 1.5, 3.0):
        est = (th + Z) * (np.abs(th + Z) > c)
        assert abs(np.mean((est - th) ** 2) - hard_risk(th, c)) < 0.01
assert np.allclose(hard_risk(thetas, 0.0), 1.0)
cs = {"aic": np.sqrt(2), "pre": 1.96, "bic": np.sqrt(np.log(1000))}
for key, c in cs.items():
    r = hard_risk(thetas, c)
    assert r.max() > 1 and r[0] < 1 and abs(r[-1] - 1) < 0.01
    gen.num(f"{key}_c", c, 3)
    gen.num(f"{key}_r0", r[0], 3)
    gen.num(f"{key}_max", r.max(), 3)
    gen.num(f"{key}_argmax", thetas[r.argmax()], 2)
    gen.num(f"{key}_cross", thetas[np.argmax(r > 1)], 2)     # where the risk first exceeds 1
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.axhline(1, color=COLORS["muted"], linewidth=0.8, linestyle="--")
ax.plot(thetas, thetas**2, color=COLORS["grid"], linewidth=1.0)
ax.plot(thetas, hard_risk(thetas, np.sqrt(2)), color=COLORS["accent"], label=r"$C_p$/AIC, $c=\sqrt{2}$")
ax.plot(thetas, hard_risk(thetas, 1.96), color=COLORS["third"], label=r"5% pretest, $c=1.96$")
ax.plot(thetas, hard_risk(thetas, np.sqrt(np.log(1000))), color=COLORS["second"],
        label=r"BIC ($n=1000$), $c=2.63$")
ax.text(4.55, 0.86, "least squares", fontsize=7, color=COLORS["muted"])
ax.text(0.35, 1.95, r"drop ($\theta^2$)", fontsize=7, color=COLORS["muted"])
ax.set_ylim(0, 4.2)
ax.set_xlim(0, 6)
ax.set_xlabel(r"$\theta$ (in standard errors)")
ax.set_ylabel("risk")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(figure_path("ch29", "pretest_risk"))
