"""Chapter 43, Section 1: local polynomial regression, bandwidth, and boundary bias.

Simulated design on [0, 1] with a known mean function, so that bias and variance can be
computed by Monte Carlo and compared with the asymptotic expansion of Theorem 43.1.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style

# <<truth>>
A1, C1, M1 = 1.0, 25.0, 0.35          # the bump
A2, C2, M2 = -0.75, 50.0, 0.80        # the dip


def f_true(x):
    """The mean function: a straight line plus a bump and a narrower dip."""
    return 2 * x + A1 * np.exp(-C1 * (x - M1) ** 2) + A2 * np.exp(-C2 * (x - M2) ** 2)


def f_second(x):
    """Its second derivative, needed for the asymptotic bandwidth."""
    g1 = (4 * C1**2 * (x - M1) ** 2 - 2 * C1) * np.exp(-C1 * (x - M1) ** 2)
    g2 = (4 * C2**2 * (x - M2) ** 2 - 2 * C2) * np.exp(-C2 * (x - M2) ** 2)
    return A1 * g1 + A2 * g2
# <</truth>>


# <<fit>>
def epanechnikov(u):
    """K(u) = (3/4)(1 - u^2) on [-1, 1], zero outside."""
    return np.where(np.abs(u) < 1, 0.75 * (1 - u**2), 0.0)


def local_poly_weights(x0, x, h, degree):
    """The weights l_i(x0) of the local polynomial fit: f-hat(x0) = sum_i l_i(x0) y_i."""
    w = epanechnikov((x - x0) / h) / h                 # kernel weights K_h(x_i - x0)
    U = np.vander(x - x0, degree + 1, increasing=True)  # columns (x_i - x0)^j
    A = U.T @ (w[:, None] * U)                          # the weighted moment matrix
    e1 = np.zeros(degree + 1)
    e1[0] = 1.0
    return np.linalg.solve(A, e1) @ U.T * w             # row vector of weights


def local_poly(grid, x, y, h, degree):
    """Local polynomial estimate of degree `degree` on a grid of target points."""
    return np.array([local_poly_weights(x0, x, h, degree) @ y for x0 in grid])
# <</fit>>


n, sigma = 200, 0.25
x = (np.arange(1, n + 1) - 0.5) / n          # the regular design of assumption (A)
rng = np.random.default_rng(20250921)
y = f_true(x) + sigma * rng.normal(size=n)

# --- the weights reproduce polynomials of degree <= q -------------------------
for degree in (0, 1, 2):
    w = local_poly_weights(0.4, x, 0.1, degree)
    for j in range(degree + 1):
        assert np.isclose(w @ (x - 0.4) ** j, 1.0 if j == 0 else 0.0, atol=1e-12)

# --- asymptotic optimal bandwidth at an interior point ------------------------
MU2, RK = 0.2, 0.6                            # \int u^2 K and \int K^2 for Epanechnikov
assert np.isclose(np.trapezoid(epanechnikov(np.linspace(-1, 1, 200001)) * np.linspace(-1, 1, 200001) ** 2,
                               np.linspace(-1, 1, 200001)), MU2, atol=1e-6)
assert np.isclose(np.trapezoid(epanechnikov(np.linspace(-1, 1, 200001)) ** 2, np.linspace(-1, 1, 200001)),
                  RK, atol=1e-6)

X0 = 0.35                                     # the top of the bump
h_star = (sigma**2 * RK / (MU2**2 * f_second(X0) ** 2)) ** 0.2 * n ** (-0.2)
print("asymptotically optimal bandwidth at x = 0.35:", round(h_star, 4))

# --- Monte Carlo bias and variance at X0 over a grid of bandwidths ------------
reps = 2000
E = rng.normal(size=(reps, n)) * sigma
hs = np.geomspace(0.02, 0.30, 29)
bias_mc, var_mc = np.empty(len(hs)), np.empty(len(hs))
for k, h in enumerate(hs):
    w = local_poly_weights(X0, x, h, 1)
    fits = (f_true(x) + E) @ w
    bias_mc[k] = fits.mean() - f_true(X0)
    var_mc[k] = fits.var(ddof=1)
mse_mc = bias_mc**2 + var_mc
h_mc = hs[np.argmin(mse_mc)]
assert abs(np.log(h_mc / h_star)) < 0.35      # the asymptotic bandwidth is in the right place

# asymptotic curves
bias_as = 0.5 * MU2 * f_second(X0) * hs**2
var_as = sigma**2 * RK / (n * hs)
k_ref = np.argmin(np.abs(hs - h_star))
ratio_bias = bias_mc[k_ref] / bias_as[k_ref]
ratio_var = var_mc[k_ref] / var_as[k_ref]
assert 0.8 < ratio_bias < 1.2 and 0.8 < ratio_var < 1.2

# --- boundary: local constant against local linear ---------------------------
hb = 0.08
grid_b = np.linspace(0.0, 0.25, 51)
mean_nw = np.array([((f_true(x) + E) @ local_poly_weights(x0, x, hb, 0)).mean() for x0 in grid_b])
mean_ll = np.array([((f_true(x) + E) @ local_poly_weights(x0, x, hb, 1)).mean() for x0 in grid_b])
bias_nw0 = mean_nw[0] - f_true(grid_b[0])
bias_ll0 = mean_ll[0] - f_true(grid_b[0])
assert abs(bias_nw0) > 8 * abs(bias_ll0)

# order h versus order h^2 at the left endpoint: halve h and watch the biases
def endpoint_bias(h, degree):
    w = local_poly_weights(0.0, x, h, degree)
    return w @ f_true(x) - f_true(0.0)


print("bias at the left endpoint, local constant:", round(endpoint_bias(0.08, 0), 4))
print("bias at the left endpoint, local linear:  ", round(endpoint_bias(0.08, 1), 4))
drop_nw = endpoint_bias(0.08, 0) / endpoint_bias(0.04, 0)
drop_ll = endpoint_bias(0.08, 1) / endpoint_bias(0.04, 1)
assert 1.7 < drop_nw < 2.3 and 3.2 < drop_ll < 4.8

gen = Generated("ch43", "local_polynomial")
gen.int("n", n)
gen.num("sigma", sigma, 2)
gen.num("x0", X0, 2)
gen.num("fpp", f_second(X0), 2)
gen.num("h_star", h_star, 4)
gen.num("h_mc", h_mc, 4)
gen.num("ratio_bias", ratio_bias, 2)
gen.num("ratio_var", ratio_var, 2)
gen.num("hb", hb, 2)
gen.num("bias_nw0", bias_nw0, 4)
gen.num("bias_ll0", bias_ll0, 4)
gen.num("ep_nw", endpoint_bias(0.08, 0), 4)
gen.num("ep_ll", endpoint_bias(0.08, 1), 4)
gen.num("drop_nw", drop_nw, 2)
gen.num("drop_ll", drop_ll, 2)
gen.num("h_ratio", (10.0) ** 0.2, 3)
gen.write()

# ---- figure: bandwidths, and the boundary -----------------------------------
use_book_style()
grid = np.linspace(0, 1, 301)
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
ax.scatter(x, y, s=5, color=COLORS["muted"], alpha=0.6, linewidths=0)
ax.plot(grid, f_true(grid), color=COLORS["ink"], linewidth=1.0, label="true $f$")
for h, col, ls in ((0.012, COLORS["second"], "-"), (0.08, COLORS["accent"], "-"), (0.30, COLORS["third"], "--")):
    ax.plot(grid, local_poly(grid, x, y, h, 1), color=col, linestyle=ls, label=f"$h={h:g}$")
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_ylim(-0.9, 3.0)
ax.set_title("(a) local linear fits")
ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.4, columnspacing=1.0)
ax = axes[1]
ax.axhline(0.0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.plot(grid_b, mean_nw - f_true(grid_b), color=COLORS["second"], linestyle="--", label="local constant")
ax.plot(grid_b, mean_ll - f_true(grid_b), color=COLORS["accent"], label="local linear")
ax.axvline(hb, color=COLORS["muted"], linewidth=0.7)
ax.text(hb + 0.004, 0.085, "$x=h$", fontsize=7, color=COLORS["muted"])
ax.set_xlabel("$x$")
ax.set_ylabel("bias of the fit")
ax.set_title("(b) bias near the left boundary")
ax.legend(frameon=False, loc="upper right", handlelength=1.4)
fig.tight_layout()
fig.savefig(figure_path("ch43", "local_linear_boundary"))

# ---- figure: bias, variance and mean squared error against h -----------------
fig, ax = plt.subplots(figsize=(3.4, 2.4))
ax.plot(hs, bias_mc**2, color=COLORS["second"], label="squared bias")
ax.plot(hs, bias_as**2, color=COLORS["second"], linestyle=":", linewidth=0.9)
ax.plot(hs, var_mc, color=COLORS["accent"], label="variance")
ax.plot(hs, var_as, color=COLORS["accent"], linestyle=":", linewidth=0.9)
ax.plot(hs, mse_mc, color=COLORS["ink"], label="mean squared error")
ax.axvline(h_star, color=COLORS["muted"], linewidth=0.7)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.3], ["0.02", "0.05", "0.1", "0.2", "0.3"])
ax.minorticks_off()
ax.set_ylim(1e-5, 3e-1)
ax.set_xlabel("bandwidth $h$")
ax.set_title("at $x=0.35$, $n=200$")
ax.legend(frameon=False, loc="upper left", handlelength=1.4)
fig.tight_layout()
fig.savefig(figure_path("ch43", "bias_variance_bandwidth"))
