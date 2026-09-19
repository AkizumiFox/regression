"""Chapter 21, Section 1: normal probability plots of studentized residuals with simulated envelopes.

A simulated regression with n = 60 observations, an intercept and two regressors. The same
design is used with three error laws: normal, t with 3 degrees of freedom (heavy tails) and
centred exponential (skewed), all with variance one. The envelope is the pointwise 2.5% and
97.5% points of the ordered externally studentized residuals over 1999 data sets simulated
from the fitted model with normal errors; their law does not depend on beta or sigma.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

# <<setup>>
rng = np.random.default_rng(2104)
n, p = 60, 3
X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
Q, _ = np.linalg.qr(X)
h = np.sum(Q**2, axis=1)                          # leverages


def ext_studentized(y):
    """Externally studentized residuals t_i = e_i / (s_(i) sqrt(1 - h_ii))."""
    e = y - Q @ (Q.T @ y)
    s2 = e @ e / (n - p)
    # squared internally studentized residuals
    r2 = e**2 / (s2 * (1 - h))
    return np.sign(e) * np.sqrt(r2 * (n - p - 1) / (n - p - r2))


# plotting positions
blom = stats.norm.ppf((np.arange(1, n + 1) - 0.375) / (n + 0.25))
sims = np.sort([ext_studentized(rng.normal(size=n)) for _ in range(1999)],
               axis=1)
# pointwise envelope
lower, upper = np.percentile(sims, [2.5, 97.5], axis=0)
# <</setup>>

# <<plots>>
laws = {
    "normal": rng.normal(size=n),
    "t(3)": rng.standard_t(3, size=n) / np.sqrt(3.0),
    "centred exponential": rng.exponential(size=n) - 1.0,
}
beta = np.array([1.0, 2.0, -1.0])
for name, eps in laws.items():
    t = np.sort(ext_studentized(X @ beta + eps))
    outside = np.sum((t < lower) | (t > upper))
    print(f"{name:20s}: {outside} of {n} points outside the envelope, "
          f"Shapiro-Wilk p = {stats.shapiro(t).pvalue:.3f}")
# <</plots>>

# the externally studentized residuals are t(n - p - 1) under normal errors
assert stats.kstest(sims[:, n // 2], lambda v: v) is not None
flat = sims.ravel()
assert abs(np.mean(flat**2) - (n - p - 1) / (n - p - 3)) < 0.03

gen = Generated("ch21", "normal_plots", prefix="nplot")
gen.int("n", n)
gen.int("p", p)
results = {}
for key, (name, eps) in zip(["normal", "t3", "exp"], laws.items()):
    t = np.sort(ext_studentized(X @ beta + eps))
    results[key] = t
    gen.int(f"outside_{key}", np.sum((t < lower) | (t > upper)))
    gen.num(f"sw_{key}", stats.shapiro(t).pvalue, 3)
assert np.sum((results["normal"] < lower) | (results["normal"] > upper)) <= 3
assert stats.shapiro(results["exp"]).pvalue < 0.01
gen.write()

# ---- figure -------------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 3, figsize=(6.2, 2.25), sharey=True)
titles = ["(a) normal errors", "(b) $t(3)$ errors", "(c) centred exponential"]
for ax, key, title in zip(axes, ["normal", "t3", "exp"], titles):
    ax.fill_between(blom, lower, upper, color=COLORS["grid"], alpha=0.8, linewidth=0,
                    label="95% envelope")
    ax.plot(blom, np.median(sims, axis=0), color=COLORS["muted"], linewidth=0.7)
    t = results[key]
    out = (t < lower) | (t > upper)
    ax.scatter(blom[~out], t[~out], s=7, color=COLORS["accent"], linewidths=0)
    ax.scatter(blom[out], t[out], s=9, color=COLORS["second"], linewidths=0)
    ax.set_title(title)
    ax.set_xlabel("normal score")
axes[0].set_ylabel("studentized residual")
fig.tight_layout()
fig.savefig(figure_path("ch21", "normal_plots"))
