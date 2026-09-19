"""Chapter 19, Section 4: what one contaminated observation does to least squares.

(a) Stack loss data (Brownlee 1965; public domain, statsmodels.datasets.stackloss): the effect of
    each observation on the coefficients, from the deletion formula, checked against refits;
    the mean-shift coefficient e_k / (1 - h_kk).
(b) Adding Delta to one response moves beta_hat by Delta (X'X)^{-1} x_k and the fitted vector
    by a vector of length |Delta| sqrt(h_kk).
(c) Figure: a shifted response at a low- and a high-leverage point; a response held fixed while
    its x value is moved away (a bad leverage point drags the slope towards zero).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch19", "outliers", prefix="out")

# <<stackloss>>
import numpy as np
import statsmodels.api as sm
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]].to_numpy()])
n, p = X.shape

G = np.linalg.inv(X.T @ X)
beta_hat = G @ X.T @ y
e_hat = y - X @ beta_hat
h = np.einsum("ij,jk,ik->i", X, G, X)             # leverages

# the change in beta_hat from deleting case k, without refitting
change = (G @ X.T) * (e_hat / (1 - h))            # column k is beta_hat - beta_hat_(k)
shift = h * e_hat**2 / (1 - h) ** 2               # ||X (beta_hat - beta_hat_(k))||^2
k = int(np.argmax(shift))
print("most influential case:", k + 1, " leverage", round(h[k], 3), " residual", round(e_hat[k], 3))
print("beta_hat           ", np.round(beta_hat, 3))
print("beta_hat without it", np.round(beta_hat - change[:, k], 3))
print("mean-shift estimate e_k/(1-h_kk):", round(e_hat[k] / (1 - h[k]), 3))
# <</stackloss>>

keep = np.arange(n) != k
b_del = np.linalg.lstsq(X[keep], y[keep], rcond=None)[0]
assert np.allclose(beta_hat - change[:, k], b_del)
# mean-shift model: add an indicator of case k
Xa = np.column_stack([X, (np.arange(n) == k).astype(float)])
ba = np.linalg.lstsq(Xa, y, rcond=None)[0]
assert np.allclose(ba[:p], b_del) and np.isclose(ba[p], e_hat[k] / (1 - h[k]))
assert np.allclose(np.linalg.norm(X @ change, axis=0) ** 2, shift)
assert k == 20                                    # the last case of the stack loss data
order = np.argsort(shift)[::-1]
gen.int("n", n)
gen.int("k", k + 1)
gen.num("h_k", h[k], 3)
gen.num("e_k", e_hat[k], 3)
gen.num("delta_k", e_hat[k] / (1 - h[k]), 3)
gen.num("pbar", p / n, 3)
gen.num("h_max", h.max(), 3)
gen.int("k_hmax", int(np.argmax(h)) + 1)
gen.num("e_hmax", e_hat[np.argmax(h)], 3)
gen.num("shift_k", np.sqrt(shift[k]), 3)
gen.num("shift_2nd", np.sqrt(shift[order[1]]), 3)
gen.int("k_2nd", int(order[1]) + 1)
gen.num("shift_median", np.sqrt(np.median(shift)), 3)
for j, name in enumerate(["b0", "b1", "b2", "b3"]):
    gen.num(name, beta_hat[j], 3)
    gen.num(name + "_del", b_del[j], 3)
s2 = e_hat @ e_hat / (n - p)
s2_del = np.sum((y[keep] - X[keep] @ b_del) ** 2) / (n - 1 - p)
gen.num("s", np.sqrt(s2), 3)
gen.num("s_del", np.sqrt(s2_del), 3)

# <<contaminate>>
Delta = 10.0
for j in [int(np.argmin(h)), int(np.argmax(h))]:
    y_c = y.copy()
    y_c[j] += Delta                                        # contaminate one response
    b_c = G @ X.T @ y_c
    print(f"case {j + 1:2d}: h = {h[j]:.3f}, ||X(b_c - b)|| = {np.linalg.norm(X @ (b_c - beta_hat)):.3f},"
          f" |Delta| sqrt(h) = {Delta * np.sqrt(h[j]):.3f}")
# <</contaminate>>
for j in [int(np.argmin(h)), int(np.argmax(h))]:
    y_c = y.copy(); y_c[j] += Delta
    b_c = G @ X.T @ y_c
    assert np.allclose(b_c - beta_hat, Delta * G @ X[j])
    assert np.isclose(np.linalg.norm(X @ (b_c - beta_hat)), Delta * np.sqrt(h[j]))
    # the residual at j absorbs only (1 - h) of the shift
    assert np.isclose((y_c - X @ b_c)[j] - e_hat[j], Delta * (1 - h[j]))
gen.num("Delta", Delta, 0)
gen.num("h_min", h.min(), 3)
gen.int("k_hmin", int(np.argmin(h)) + 1)
gen.num("move_min", Delta * np.sqrt(h.min()), 3)
gen.num("move_max", Delta * np.sqrt(h.max()), 3)
gen.num("resid_min", Delta * (1 - h.min()), 3)
gen.num("resid_max", Delta * (1 - h.max()), 3)
gen.write()

# ---- figure: a straight line ------------------------------------------------------------
rng = np.random.default_rng(1904)
xs = np.concatenate([np.linspace(1, 10, 14), [22.0]])
ys = 2 + 0.5 * xs + 0.6 * rng.normal(size=len(xs))
Xs = np.column_stack([np.ones_like(xs), xs])
hs = np.diag(Xs @ np.linalg.solve(Xs.T @ Xs, Xs.T))
fit = lambda xx, yy: np.polyfit(xx, yy, 1)
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
ax = axes[0]
grid = np.linspace(0, 23, 2)
ax.scatter(xs, ys, s=10, color=COLORS["ink"], zorder=3)
ax.plot(grid, np.polyval(fit(xs, ys), grid), color=COLORS["muted"], lw=1, label="clean")
for j, col, lab in [(6, COLORS["accent"], "shift at low leverage"),
                    (len(xs) - 1, COLORS["second"], "shift at high leverage")]:
    yc = ys.copy(); yc[j] += 6.0
    ax.scatter([xs[j]], [yc[j]], s=18, facecolors="none", edgecolors=col, zorder=4)
    ax.plot(grid, np.polyval(fit(xs, yc), grid), color=col, lw=1, label=lab)
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.set_title(r"(a) one response moved up by 6")
ax.legend(frameon=False, loc="upper left", fontsize=7)
ax = axes[1]
j = 6                                                      # a point in the middle of the cloud
far = np.linspace(xs[j], 200, 300)
slopes = []
for xf in far:
    xm = xs[:-1].copy(); xm[j] = xf
    slopes.append(fit(xm, ys[:-1])[0])
ax.plot(far, slopes, color=COLORS["accent"])
ax.axhline(fit(xs[:-1], ys[:-1])[0], color=COLORS["muted"], lw=0.8, ls=":")
ax.axhline(0, color=COLORS["grid"], lw=0.8)
ax.set_xlabel("$x$ value of the moved case")
ax.set_ylabel("fitted slope")
ax.set_title("(b) a bad leverage point")
fig.tight_layout()
fig.savefig(figure_path("ch19", "outliers"))
assert abs(slopes[-1]) < 0.2 * abs(fit(xs[:-1], ys[:-1])[0])
assert hs[-1] > 0.5
