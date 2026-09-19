"""Chapter 14, Section 2: the prediction error of an estimated best linear predictor.

(a) Lemma: if W ~ Wishart_k(m, Sigma), E(W^{-1}) = Sigma^{-1} / (m - k - 1).
(b) Multivariate normal rows: the expected squared error of the fitted predictor at a new
    random point is sigma^2 (1 + 1/n)(n - 2)/(n - k - 2), while the in-sample residual mean
    square SSE/n has expectation sigma^2 (n - k - 1)/n. Figure: the optimism as k grows.
(c) RAND Health Insurance Experiment (public domain, statsmodels.datasets.randhie): the
    same comparison with the 20,190 rows as the population, and non-normal regressors.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

gen = Generated("ch14", "prediction_error", prefix="pe")

# ---- (a) the mean of an inverse Wishart matrix ---------------------------------
rng = np.random.default_rng(1421)
m, k = 20, 4
A = rng.normal(size=(k, k))
Sigma = A @ A.T + k * np.eye(k)
L = np.linalg.cholesky(Sigma)
Z = rng.normal(size=(50_000, m, k)) @ L.T            # rows iid N(0, Sigma)
W = np.einsum("rij,rik->rjk", Z, Z)                  # W = Z^T Z ~ Wishart_k(m, Sigma)
EWinv = np.linalg.inv(W).mean(axis=0)
target = np.linalg.inv(Sigma) / (m - k - 1)
rel = np.max(np.abs(EWinv - target)) / np.max(np.abs(target))
assert rel < 0.02
gen.num("wishart_rel", rel, 4)

# ---- (b) prediction error of the fitted predictor, normal rows -----------------
# <<msep>>
rng = np.random.default_rng(1422)
n, sigma2, reps = 30, 1.0, 20_000

def simulate(k):
    """Average in-sample and out-of-sample squared error over many training samples."""
    ins, out = np.empty(reps), np.empty(reps)
    beta = np.full(k, 0.5)
    for s in range(reps):
        X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
        y = X @ np.r_[1.0, beta] + rng.normal(size=n)
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        ins[s] = np.mean((y - X @ b) ** 2)           # SSE / n
        # new point X0 ~ N(0, I), Y0 = 1 + beta^T X0 + e0: exact conditional error
        out[s] = sigma2 + (b[0] - 1.0) ** 2 + np.sum((b[1:] - beta) ** 2)
    return ins.mean(), out.mean()

results = {}
for k in (1, 5, 10, 20):
    ins, out = simulate(k)
    results[k] = (ins, out)
    theory = sigma2 * (1 + 1 / n) * (n - 2) / (n - k - 2)
    print(f"k = {k:2d}: in-sample {ins:.3f} (theory {sigma2 * (n - k - 1) / n:.3f}),"
          f" new point {out:.3f} (theory {theory:.3f})")
# <</msep>>
for k in (1, 5, 10, 20):
    gen.num(f"in{k}", results[k][0], 3)
    gen.num(f"out{k}", results[k][1], 3)
    gen.num(f"thin{k}", sigma2 * (n - k - 1) / n, 3)
    gen.num(f"thout{k}", sigma2 * (1 + 1 / n) * (n - 2) / (n - k - 2), 3)
    gen.num(f"thfix{k}", sigma2 * (n + k + 1) / n, 3)
ks = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
sim_in, sim_out = {}, {}
for k in ks:
    sim_in[k], sim_out[k] = results[k] if k in results else simulate(k)
    theory_out = sigma2 * (1 + 1 / n) * (n - 2) / (n - k - 2)
    theory_in = sigma2 * (n - k - 1) / n
    assert abs(sim_out[k] / theory_out - 1) < 0.03, (k, sim_out[k], theory_out)
    assert abs(sim_in[k] / theory_in - 1) < 0.02
gen.int("n", n)
gen.int("reps", reps)

# exercise exr-cor-omitted-prediction: n = 30, k = 10, sigma_1^2 = 1.2 sigma^2
def small_err(k1, ratio=1.2):
    return ratio * (1 + 1 / n) * (n - 2) / (n - k1 - 2)
full_err = (1 + 1 / n) * (n - 2) / (n - 10 - 2)
k1_max = max(k1 for k1 in range(1, 10) if small_err(k1) < full_err)
assert small_err(4) < full_err and k1_max == 6 and small_err(7) > full_err
gen.num("exsmall4", small_err(4), 3)
gen.num("exsmall6", small_err(6), 3)
gen.num("exsmall7", small_err(7), 3)
gen.int("exk1max", k1_max)

# ---- (c) RAND data: the population is the empirical distribution of the rows ----
# <<rand>>
data = sm.datasets.randhie.load_pandas().data
cols = ["lncoins", "idp", "lpi", "fmde", "physlm", "disea"]
Xpop = np.column_stack([np.ones(len(data)), data[cols].to_numpy()])
ypop = np.log1p(data["mdvis"].to_numpy())            # log(1 + number of visits)
b_pop, *_ = np.linalg.lstsq(Xpop, ypop, rcond=None)  # the population projection
sigma2_pop = np.mean((ypop - Xpop @ b_pop) ** 2)     # its mean squared error

rng = np.random.default_rng(1423)
n_train, draws = 40, 4000
k_rand = len(cols)
ins, out = np.empty(draws), np.empty(draws)
for s in range(draws):
    idx = rng.integers(0, len(ypop), n_train)        # iid draws from the population
    b, *_ = np.linalg.lstsq(Xpop[idx], ypop[idx], rcond=None)
    ins[s] = np.mean((ypop[idx] - Xpop[idx] @ b) ** 2)
    out[s] = np.mean((ypop - Xpop @ b) ** 2)         # exact error over the population
print(f"population projection error {sigma2_pop:.4f}")
print(f"in-sample {ins.mean():.4f}, new case {out.mean():.4f}")
print(f"normal-theory values: {sigma2_pop * (n_train - k_rand - 1) / n_train:.4f},"
      f" {sigma2_pop * (1 + 1 / n_train) * (n_train - 2) / (n_train - k_rand - 2):.4f}")
# <</rand>>
th_in = sigma2_pop * (n_train - k_rand - 1) / n_train
th_out = sigma2_pop * (1 + 1 / n_train) * (n_train - 2) / (n_train - k_rand - 2)
assert ins.mean() < sigma2_pop < out.mean()
assert abs(out.mean() / th_out - 1) < 0.1
gen.int("N", len(ypop))
gen.int("ntrain", n_train)
gen.int("draws", draws)
gen.int("krand", k_rand)
gen.num("sig2", sigma2_pop, 4)
gen.num("rin", ins.mean(), 4)
gen.num("rout", out.mean(), 4)
gen.num("rthin", th_in, 4)
gen.num("rthout", th_out, 4)
gen.num("r2pop", 1 - sigma2_pop / ypop.var(), 3)
gen.write()

# ---- figure -----------------------------------------------------------------------
use_book_style()
fig, ax = plt.subplots(figsize=(4.6, 2.6))
kk = np.linspace(0, 25, 200)
ax.plot(kk, sigma2 * (1 + 1 / n) * (n - 2) / (n - kk - 2), color=COLORS["second"],
        label="new random point")
ax.plot(kk, sigma2 * (n + kk + 1) / n, color=COLORS["third"], label="new response, same design")
ax.plot(kk, sigma2 * (n - kk - 1) / n, color=COLORS["accent"], label=r"in-sample, SSE$/n$")
ax.plot(ks, [sim_out[k] for k in ks], "o", ms=3, color=COLORS["second"])
ax.plot(ks, [sim_in[k] for k in ks], "o", ms=3, color=COLORS["accent"])
ax.axhline(sigma2, color=COLORS["muted"], lw=0.8, ls="--")
ax.text(23.6, 0.82, r"$\sigma^2$", color=COLORS["muted"], fontsize=8)
ax.set_ylim(0, 4)
ax.set_xlim(0, 25)
ax.set_xlabel(r"number of regressors $k$ ($n=30$)")
ax.set_ylabel(r"expected squared error / $\sigma^2$")
ax.legend(frameon=False, loc="upper left")
fig.savefig(figure_path("ch14", "prediction_error"))
