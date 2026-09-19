"""Chapter 23, Sections 1-2: residual, case and wild bootstraps for the Engel food-expenditure line.

Engel's 1857 Belgian household budgets (statsmodels.datasets.engel, public domain): annual food
expenditure against annual income for 235 households. The spread of food expenditure grows with
income, so the three bootstraps estimate different things: the residual bootstrap reproduces the
classical standard error, the case and wild bootstraps reproduce heteroscedasticity-consistent ones.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

data = sm.datasets.engel.load_pandas().data

# <<setup>>
y = data["foodexp"].to_numpy()
x = data["income"].to_numpy()
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape
XtX_inv = np.linalg.inv(X.T @ X)
beta_hat = XtX_inv @ X.T @ y
e_hat = y - X @ beta_hat
h = np.sum((X @ XtX_inv) * X, axis=1)            # leverages
s2 = e_hat @ e_hat / (n - p)
se_classical = np.sqrt(s2 * XtX_inv[1, 1])
print(f"slope {beta_hat[1]:.4f}, classical se {se_classical:.5f}, max leverage {h.max():.3f}")
# <</setup>>

# <<residual>>
rng = np.random.default_rng(2301)
B = 4000
A = XtX_inv @ X.T                                 # beta_hat = A y
fit = X @ beta_hat
e_c = e_hat - e_hat.mean()                        # centred residuals
r = e_hat / np.sqrt(1 - h)                        # leverage-adjusted residuals
r_c = r - r.mean()


def residual_bootstrap(resid, B):
    """Slopes from y* = X beta_hat + e*, with e* drawn with replacement from resid."""
    idx = rng.integers(0, n, size=(B, n))
    y_star = fit + resid[idx]
    return (y_star @ A.T)[:, 1]


slopes_res = residual_bootstrap(e_c, B)
slopes_mod = residual_bootstrap(r_c, B)
print(f"residual bootstrap se: raw {slopes_res.std():.5f}, leverage-adjusted {slopes_mod.std():.5f}")
# <</residual>>

# <<pairs>>
def pairs_bootstrap(B):
    """Slopes refitted to n cases drawn with replacement from the (x_i, y_i)."""
    idx = rng.integers(0, n, size=(B, n))
    xs, ys = x[idx], y[idx]
    xc = xs - xs.mean(axis=1, keepdims=True)
    return np.sum(xc * ys, axis=1) / np.sum(xc ** 2, axis=1)


slopes_pairs = pairs_bootstrap(B)
print(f"case bootstrap se {slopes_pairs.std():.5f}")
# <</pairs>>

# <<wild>>
golden = (1 + np.sqrt(5)) / 2
mammen_values = np.array([1 - golden, golden])    # -0.618 and 1.618
mammen_probs = np.array([golden / np.sqrt(5), 1 - golden / np.sqrt(5)])


def wild_bootstrap(resid, B, weights="rademacher"):
    """Slopes from y*_i = x_i' beta_hat + resid_i v_i, with independent multipliers v_i."""
    if weights == "rademacher":
        v = rng.choice([-1.0, 1.0], size=(B, n))
    else:
        v = rng.choice(mammen_values, p=mammen_probs, size=(B, n))
    y_star = fit + resid * v
    return (y_star @ A.T)[:, 1]


slopes_wild = wild_bootstrap(e_hat, B)
slopes_mammen = wild_bootstrap(e_hat, B, weights="mammen")
slopes_wild3 = wild_bootstrap(e_hat / (1 - h), B)
se_hc0 = np.sqrt(A[1] ** 2 @ e_hat ** 2)
se_hc3 = np.sqrt(A[1] ** 2 @ (e_hat / (1 - h)) ** 2)
print(f"wild se: Rademacher {slopes_wild.std():.5f}, Mammen {slopes_mammen.std():.5f}, "
      f"HC3-type {slopes_wild3.std():.5f}")
print(f"sandwich se: HC0 {se_hc0:.5f}, HC3 {se_hc3:.5f}")
# <</wild>>

# ---- checks -----------------------------------------------------------------------------
# exact bootstrap moments (Proposition bs-residual-moments and bs-wild-moments)
se_res_exact = np.sqrt(np.mean(e_c ** 2) * XtX_inv[1, 1])
se_mod_exact = np.sqrt(np.mean(r_c ** 2) * XtX_inv[1, 1])
assert np.isclose(se_res_exact, se_classical * np.sqrt((n - p) / n))
se_hc2 = np.sqrt(A[1] ** 2 @ (e_hat ** 2 / (1 - h)))


def mc_ok(sample, target, tol_sd=4.0):
    """The Monte Carlo standard deviation is within tol_sd standard errors of the target."""
    se_of_sd = target / np.sqrt(2 * (len(sample) - 1))   # rough; kurtosis near 3 for these sums
    return abs(sample.std() - target) < tol_sd * se_of_sd * 1.5


assert mc_ok(slopes_res, se_res_exact)
assert mc_ok(slopes_mod, se_mod_exact)
assert mc_ok(slopes_wild, se_hc0)
assert mc_ok(slopes_mammen, se_hc0)
assert mc_ok(slopes_wild3, se_hc3)
assert abs(slopes_res.mean() - beta_hat[1]) < 4 * se_res_exact / np.sqrt(B)
# the Mammen weights: mean 0, variance 1, third moment 1
for k, target in [(1, 0.0), (2, 1.0), (3, 1.0)]:
    assert np.isclose(mammen_probs @ mammen_values ** k, target)
# the one-step case bootstrap has covariance exactly X' diag(e^2) X (Proposition bs-case-linear)
W_cov = (np.eye(n) - 1.0 / n)                     # covariance of multinomial(n, 1/n) counts
lin_cov = (X * e_hat[:, None]).T @ W_cov @ (X * e_hat[:, None])
assert np.allclose(lin_cov, (X * e_hat[:, None] ** 2).T @ X)
# heteroscedasticity: the case and wild bootstraps agree with the sandwich, not with the classical se
assert se_hc0 > 1.4 * se_classical
assert abs(slopes_pairs.std() / se_hc0 - 1) < 0.2
assert abs(slopes_res.std() / se_classical - 1) < 0.05
# the spread of the residuals grows with income
lo = x < np.median(x)
spread_ratio = e_hat[~lo].std() / e_hat[lo].std()
assert spread_ratio > 2

# one household dominates the wild-bootstrap sum sum_i a_i e_i v_i (the Lindeberg ratio of Theorem
# bs-wild-consistency), which is why the wild distribution has two modes
terms = (A[1] * e_hat) ** 2
share = terms.max() / terms.sum()
top = np.argmax(terms)
assert top == np.argmax(x) and share > 0.5
shift = abs(A[1, top] * e_hat[top])              # the modes sit near beta_hat -/+ shift
assert abs(np.median(slopes_wild[slopes_wild > beta_hat[1]]) - (beta_hat[1] + shift)) < 0.01
assert abs(np.median(slopes_wild[slopes_wild < beta_hat[1]]) - (beta_hat[1] - shift)) < 0.01

gen = Generated("ch23", "engel_bootstrap", prefix="engel")
gen.num("share", share, 2)
gen.num("income_top", x[top], 0)
gen.num("shift", shift, 3)
gen.num("mode_lo", beta_hat[1] - shift, 3)
gen.num("mode_hi", beta_hat[1] + shift, 3)
gen.int("n", n)
gen.int("B", B)
gen.num("slope", beta_hat[1], 4)
gen.num("intercept", beta_hat[0], 2)
gen.num("hmax", h.max(), 3)
gen.num("se_classical", se_classical, 5)
gen.num("se_res_exact", se_res_exact, 5)
gen.num("se_res", slopes_res.std(), 5)
gen.num("se_mod_exact", se_mod_exact, 5)
gen.num("se_mod", slopes_mod.std(), 5)
gen.num("se_pairs", slopes_pairs.std(), 5)
gen.num("se_hc0", se_hc0, 5)
gen.num("se_hc2", se_hc2, 5)
gen.num("se_hc3", se_hc3, 5)
gen.num("se_wild", slopes_wild.std(), 5)
gen.num("se_mammen", slopes_mammen.std(), 5)
gen.num("se_wild3", slopes_wild3.std(), 5)
gen.num("factor", np.sqrt((n - p) / n), 4)
gen.num("spread_ratio", spread_ratio, 2)
gen.num("hc0_ratio", se_hc0 / se_classical, 2)
gen.write()

# ---- figure -------------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5))
ax = axes[0]
ax.scatter(x, y, s=7, color=COLORS["accent"], alpha=0.7, linewidths=0)
xs = np.linspace(x.min(), x.max(), 2)
ax.plot(xs, beta_hat[0] + beta_hat[1] * xs, color=COLORS["second"])
ax.set_xlabel("income")
ax.set_ylabel("food expenditure")
ax.set_title("(a) Engel's households")
ax = axes[1]
bins = np.linspace(beta_hat[1] - 0.12, beta_hat[1] + 0.12, 61)
for s_, lab, col in [(slopes_res, "residual", COLORS["muted"]),
                     (slopes_pairs, "case", COLORS["accent"]),
                     (slopes_wild, "wild", COLORS["second"])]:
    ax.hist(s_, bins=bins, histtype="step", density=True, color=col, label=lab, linewidth=1.0)
ax.axvline(beta_hat[1], color=COLORS["ink"], linewidth=0.6)
ax.set_xlabel("bootstrap slope")
ax.set_yticks([])
ax.legend(frameon=False, loc="upper left")
ax.set_title("(b) three bootstrap distributions")
fig.tight_layout()
fig.savefig(figure_path("ch23", "engel_bootstrap"))
