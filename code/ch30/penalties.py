"""Chapter 30, Section 1: penalties and the univariate thresholding rules they induce.

For one coordinate, minimize (1/2)(z - theta)^2 + P(theta). The script checks the closed-form
rules (soft thresholding for the lasso, proportional shrinkage for ridge, hard thresholding for the
l0 penalty, firm thresholding for MCP, and the SCAD rule of Fan and Li) against brute-force
minimization on a fine grid, and compares the posterior mode and the posterior mean of theta
given z ~ N(theta, 1) under a Laplace prior.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from regbook import COLORS, Generated, figure_path, use_book_style

lam, gam, a_scad = 1.0, 3.0, 3.7


# <<rules>>
def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)

def mcp_penalty(t, lam, gam):
    t = np.abs(t)
    return np.where(t <= gam * lam, lam * t - t ** 2 / (2 * gam), gam * lam ** 2 / 2)

def mcp_rule(z, lam, gam):                    # firm thresholding
    az = np.abs(z)
    mid = np.sign(z) * (az - lam) / (1 - 1 / gam)
    return np.where(az <= lam, 0.0, np.where(az <= gam * lam, mid, z))

def scad_penalty(t, lam, a):
    t = np.abs(t)
    quad = (2 * a * lam * t - t ** 2 - lam ** 2) / (2 * (a - 1))
    return np.where(t <= lam, lam * t, np.where(t <= a * lam, quad, (a + 1) * lam ** 2 / 2))

def scad_rule(z, lam, a):
    az = np.abs(z)
    mid = ((a - 1) * z - np.sign(z) * a * lam) / (a - 2)
    return np.where(az <= 2 * lam, soft(z, lam), np.where(az <= a * lam, mid, z))

def hard_rule(z, lam):                        # penalty (lam^2 / 2) * 1{theta != 0}
    return np.where(np.abs(z) > lam, z, 0.0)
# <</rules>>


def brute(z, pen):
    grid = np.linspace(-8, 8, 160001)
    return grid[np.argmin(0.5 * (z - grid) ** 2 + pen(grid))]


zs = np.linspace(-5.5, 5.5, 45)
checks = {
    "lasso": (lambda t: lam * np.abs(t), lambda z: soft(z, lam)),
    "ridge": (lambda t: 0.5 * lam * t ** 2, lambda z: z / (1 + lam)),
    "hard": (lambda t: 0.5 * lam ** 2 * (np.abs(t) > 1e-12), lambda z: hard_rule(z, lam)),
    "mcp": (lambda t: mcp_penalty(t, lam, gam), lambda z: mcp_rule(z, lam, gam)),
    "scad": (lambda t: scad_penalty(t, lam, a_scad), lambda z: scad_rule(z, lam, a_scad)),
}
for name, (pen, rule) in checks.items():
    for z in zs:
        if name == "hard" and abs(abs(z) - lam) < 1e-9:
            continue                          # a tie between 0 and z at |z| = lam
        assert abs(brute(z, pen) - rule(z)) < 2e-4, (name, z)

# continuity of the MCP and SCAD rules, and unbiasedness for large |z|
assert np.isclose(mcp_rule(gam * lam, lam, gam), gam * lam)
assert np.isclose(scad_rule(2 * lam, lam, a_scad), soft(2 * lam, lam))
assert np.isclose(scad_rule(a_scad * lam, lam, a_scad), a_scad * lam)

# ---- Laplace prior: posterior mode is soft thresholding, the posterior mean is never zero ------
# <<laplace>>
def laplace_posterior_mean(z, rate):
    """E(theta | z) when z ~ N(theta, 1) and theta has density (rate/2) exp(-rate |theta|)."""
    w = lambda t: np.exp(-0.5 * (z - t) ** 2 - rate * abs(t))
    num = integrate.quad(lambda t: t * w(t), -30, 30, points=[0.0, z])[0]
    den = integrate.quad(w, -30, 30, points=[0.0, z])[0]
    return num / den

rate = 1.5
for z in [0.5, 1.0, 2.0, 4.0]:
    print(f"z = {z}: posterior mode {float(soft(z, rate)):.3f}, "
          f"posterior mean {laplace_posterior_mean(z, rate):.3f}")
# <</laplace>>

pm = {z: laplace_posterior_mean(z, rate) for z in [0.5, 1.0, 2.0, 4.0]}
assert all(v > 0 for v in pm.values())       # never exactly zero
assert soft(1.0, rate) == 0.0
# for large z the posterior mean approaches z - rate, the mode
assert abs(laplace_posterior_mean(8.0, rate) - (8.0 - rate)) < 1e-3

gen = Generated("ch30", "penalties")
gen.num("rate", rate, 1)
gen.num("pm_half", pm[0.5], 3)
gen.num("pm_one", pm[1.0], 3)
gen.num("pm_two", pm[2.0], 3)
gen.num("pm_four", pm[4.0], 3)
gen.num("mode_two", float(soft(2.0, rate)), 3)
gen.num("mode_four", float(soft(4.0, rate)), 3)
gen.write()

# ---- figure -------------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.2))
t = np.linspace(-4.2, 4.2, 841)
ax = axes[0]
ax.plot(t, lam * np.abs(t), color=COLORS["accent"], label="lasso")
ax.plot(t, scad_penalty(t, lam, a_scad), color=COLORS["second"], label="SCAD")
ax.plot(t, mcp_penalty(t, lam, gam), color=COLORS["third"], label="MCP")
ax.plot(t, 0.5 * lam ** 2 * (np.abs(t) > 0), color=COLORS["muted"], linestyle="--", label=r"$\ell_0$")
ax.set_ylim(0, 4.4)
ax.set_xlabel(r"$\theta$")
ax.set_title(r"(a) penalties, $\lambda=1$")
ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=6.5, columnspacing=0.8, handlelength=1.4)
ax = axes[1]
z = np.linspace(-5, 5, 1001)
ax.plot(z, z, color=COLORS["grid"], linewidth=0.8)
ax.plot(z, soft(z, lam), color=COLORS["accent"], label="soft")
ax.plot(z, scad_rule(z, lam, a_scad), color=COLORS["second"], label="SCAD")
ax.plot(z, mcp_rule(z, lam, gam), color=COLORS["third"], label="MCP")
for piece, lab in [(z < -lam, "hard"), (np.abs(z) <= lam, None), (z > lam, None)]:   # drawn with its jumps
    ax.plot(z[piece], hard_rule(z[piece], lam), color=COLORS["muted"], linestyle="--", label=lab)
ax.set_xlabel(r"$z$")
ax.set_title(r"(b) thresholding rules")
ax.legend(frameon=False, loc="upper left", fontsize=7)
ax = axes[2]
z = np.linspace(-5, 5, 201)
ax.plot(z, z, color=COLORS["grid"], linewidth=0.8)
ax.plot(z, soft(z, rate), color=COLORS["accent"], label="mode")
ax.plot(z, [laplace_posterior_mean(v, rate) for v in z], color=COLORS["thread"], label="mean")
ax.set_xlabel(r"$z$")
ax.set_title(r"(c) Laplace prior, rate 1.5")
ax.legend(frameon=False, loc="upper left", fontsize=7)
fig.tight_layout()
fig.savefig(figure_path("ch30", "thresholding"))
