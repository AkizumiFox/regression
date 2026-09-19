"""Chapter 9, Section 5: the null distribution of R^2 is Beta((r-1)/2, (n-r)/2),
so R^2 is biased upwards, while the adjusted R^2 has mean exactly zero."""
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from regbook import COLORS, Generated, figure_path, use_book_style

rng = np.random.default_rng(20260919)
n, reps = 30, 100_000
gen = Generated("ch09", "r2_null", prefix="r2n")
gen.int("n", n)
gen.int("reps", reps)
results = {}
# <<simulate>>
for k in (1, 5, 10):                                   # number of regressors besides 1
    X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
    r = k + 1
    Q, _ = np.linalg.qr(X)
    Y = 3.0 + rng.normal(size=(n, reps))               # no regressor matters
    Z = Q.T @ Y
    ssr = np.sum(Z[1:] ** 2, axis=0)
    sst = np.sum(Y ** 2, axis=0) - Z[0] ** 2
    R2 = ssr / sst
    R2adj = 1 - (1 - R2) * (n - 1) / (n - r)
    print(f"k = {k:2d}: mean R^2 = {R2.mean():.4f} (theory {k / (n - 1):.4f}),"
          f" mean adjusted R^2 = {R2adj.mean():+.4f}")
    results[k] = R2
# <</simulate>>
    se = R2.std() / np.sqrt(reps)
    assert abs(R2.mean() - k / (n - 1)) < 4 * se
    assert abs(R2adj.mean()) < 4 * R2adj.std() / np.sqrt(reps)
    ks = stats.kstest(R2, stats.beta(k / 2, (n - r) / 2).cdf)
    assert ks.pvalue > 1e-3
    gen.num(f"mean{k}", R2.mean(), 4)
    gen.num(f"theory{k}", k / (n - 1), 4)
    gen.num(f"adj{k}", R2adj.mean(), 4)
    gen.num(f"q95_{k}", stats.beta(k / 2, (n - r) / 2).ppf(0.95), 3)
    gen.num(f"neg{k}", np.mean(R2adj < 0), 3)
gen.write()

use_book_style()
fig, ax = plt.subplots(figsize=(4.4, 2.4))
grid = np.linspace(0.0005, 0.9995, 400)
for k, c in zip((1, 5, 10), (COLORS["accent"], COLORS["second"], COLORS["third"])):
    ax.hist(results[k], bins=np.linspace(0, 1, 81), density=True, color=c, alpha=0.25)
    ax.plot(grid, stats.beta(k / 2, (n - k - 1) / 2).pdf(grid), color=c,
            label=f"{k} regressor{'s' if k > 1 else ''}")
ax.set_ylim(0, 8)
ax.set_xlim(0, 0.8)
ax.set_xlabel(r"$R^2$")
ax.set_ylabel("density")
ax.legend(frameon=False)
fig.savefig(figure_path("ch09", "r2_null"))
