"""Chapter 36, Section 1: baseline-category logits on the travel mode choice data.

210 travellers between Sydney and Melbourne choose one of four modes (air, train,
bus, car). The individual-level regressors are household income and the size of the
travelling party. Public-domain data shipped with statsmodels
(statsmodels.datasets.modechoice); statsmodels takes them from the online data
archive of Greene's econometrics text.

The script fits the baseline-category logit model by Newton-Raphson in numpy, checks
it against statsmodels MNLogit, and checks the invariance of the fit to the choice of
baseline category.
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

# <<data>>
MODES = ["air", "train", "bus", "car"]          # car is category 4, the baseline

raw = sm.datasets.modechoice.load_pandas().data
chosen = raw[raw["choice"] == 1].sort_values("individual").reset_index(drop=True)
y = chosen["mode"].astype(int).to_numpy()       # 1 air, 2 train, 3 bus, 4 car
n, c = len(y), 4
X = np.column_stack([np.ones(n), chosen["hinc"] / 10.0, chosen["psize"]])
Y = np.zeros((n, c))
Y[np.arange(n), y - 1] = 1.0                    # indicator matrix, baseline last

print("chosen mode counts:", Y.sum(axis=0))
# <</data>>


# <<fit>>
def probs(X, B):
    """Fitted probabilities of the baseline-category logit model; B is p x (c-1)."""
    eta = np.column_stack([X @ B, np.zeros(len(X))])     # baseline has eta = 0
    eta = eta - eta.max(axis=1, keepdims=True)           # stabilize the exponentials
    E = np.exp(eta)
    return E / E.sum(axis=1, keepdims=True)


def loglik(X, Y, B):
    return float(np.sum(Y * np.log(probs(X, B))))


def mnlogit(X, Y, tol=1e-10, maxit=50):
    """Newton-Raphson (= Fisher scoring, the link is canonical) for baseline logits."""
    n, p = X.shape
    c = Y.shape[1]
    B = np.zeros((p, c - 1))
    for it in range(1, maxit + 1):
        P = probs(X, B)
        U = (X.T @ (Y[:, :c - 1] - P[:, :c - 1])).T.ravel()       # score, stacked by category
        J = np.zeros(((c - 1) * p, (c - 1) * p))                  # expected = observed information
        for r in range(c - 1):
            for s in range(c - 1):
                w = P[:, r] * ((r == s) - P[:, s])
                J[r * p:(r + 1) * p, s * p:(s + 1) * p] = X.T @ (w[:, None] * X)
        step = np.linalg.solve(J, U)
        B = B + step.reshape(c - 1, p).T
        if np.max(np.abs(step)) < tol:
            break
    return B, J, it


B, J, iters = mnlogit(X, Y)
se = np.sqrt(np.diag(np.linalg.inv(J))).reshape(c - 1, X.shape[1]).T

for r in range(c - 1):
    print(f"{MODES[r]:>6s} vs car:  " + "  ".join(f"{b:7.4f} ({s:.4f})" for b, s in zip(B[:, r], se[:, r])))
print(f"log-likelihood {loglik(X, Y, B):.4f} after {iters} Newton steps")
# <</fit>>

ll = loglik(X, Y, B)

# ---- check against statsmodels MNLogit (its baseline is the smallest code) ----
ycode = np.where(y == 4, 0, y)                  # 0 car, 1 air, 2 train, 3 bus
sm_fit = sm.MNLogit(ycode, X).fit(disp=0)
assert np.allclose(sm_fit.params, B, atol=1e-6)
assert np.isclose(sm_fit.llf, ll, atol=1e-8)
assert np.allclose(np.asarray(sm_fit.bse), se, rtol=1e-5)

# ---- the null model and the likelihood ratio statistic -----------------------
B0, _, _ = mnlogit(X[:, :1], Y)
ll0 = loglik(X[:, :1], Y, B0)
lr = 2 * (ll - ll0)
assert np.isclose(ll0, sm_fit.llnull, atol=1e-8)
assert np.allclose(probs(X[:, :1], B0)[0], Y.mean(axis=0), atol=1e-8)   # null fit = sample shares

# <<baseline>>
# changing the baseline: subtract the train column (category 2) from every column
Btrain = B - B[:, [1]]
Btrain = np.delete(np.column_stack([Btrain, -B[:, 1]]), 1, axis=1)      # drop train, add car
Ytrain = Y[:, [0, 2, 3, 1]]                                            # put train last
Bcheck, _, _ = mnlogit(X, Ytrain)

print("air vs train, refitted :", np.round(Bcheck[:, 0], 6))
print("air vs train, subtracted:", np.round(Btrain[:, 0], 6))
print("fitted probabilities identical:", np.allclose(probs(X, Bcheck)[:, [0, 3, 1, 2]], probs(X, B)))
# <</baseline>>

assert np.allclose(Bcheck, Btrain, atol=1e-7)
assert np.allclose(probs(X, Bcheck)[:, [0, 3, 1, 2]], probs(X, B), atol=1e-10)

# <<derivative>>
# the derivative of a fitted probability is not the coefficient
Bfull = np.column_stack([B, np.zeros(X.shape[1])])          # append the baseline's column of zeros
x0 = np.array([1.0, 1.0, 1.0])                              # income $10,000, party of one
p0 = probs(x0[None, :], B)[0]
dp0 = p0 * (Bfull[1] - p0 @ Bfull[1])                       # d pi_r / d (income in $10,000s)

print("fitted probabilities", np.round(p0, 4))
print("derivatives         ", np.round(dp0, 4))
print("income coefficients ", np.round(Bfull[1], 4))
# <</derivative>>

assert np.isclose(dp0.sum(), 0.0, atol=1e-12)
assert B[1, 2] < 0 < dp0[2]          # bus: negative coefficient, rising fitted probability
assert dp0[3] > 0                    # car: zero coefficient, rising fitted probability
assert int((dp0 > dp0[3]).sum()) == 1        # only air's derivative exceeds the car's
assert dp0[0] > dp0[3] > dp0[2] > dp0[1]     # air, car, bus, train in that order

gen = Generated("ch36", "mode")
gen.int("n", n)
for r, m in enumerate(MODES):
    gen.int(f"count_{m}", int(Y[:, r].sum()))
for r, m in enumerate(MODES[:3]):
    gen.num(f"const_{m}", B[0, r], 4)
    gen.num(f"hinc_{m}", B[1, r], 4)
    gen.num(f"psize_{m}", B[2, r], 4)
    gen.num(f"se_const_{m}", se[0, r], 4)
    gen.num(f"se_hinc_{m}", se[1, r], 4)
    gen.num(f"se_psize_{m}", se[2, r], 4)
gen.num("ll", ll, 4)
gen.num("ll0", ll0, 4)
gen.num("lr", lr, 3)
gen.int("iters", iters)
gen.num("or_train", np.exp(B[1, 1]), 4)
gen.num("air_vs_train", B[1, 0] - B[1, 1], 4)
for r, m in enumerate(MODES):
    gen.num(f"p0_{m}", p0[r], 4)
    gen.num(f"dp0_{m}", dp0[r], 4)
gen.write()

# ---- figure: fitted probabilities against household income -------------------
use_book_style()
grid = np.linspace(2.0, 72.0, 200) / 10.0
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.5), sharey=True)
palette = [COLORS["accent"], COLORS["second"], COLORS["third"], COLORS["thread"]]
for ax, size in zip(axes, [1, 3]):
    Xg = np.column_stack([np.ones_like(grid), grid, np.full_like(grid, float(size))])
    Pg = probs(Xg, B)
    for r, m in enumerate(MODES):
        ax.plot(10 * grid, Pg[:, r], color=palette[r], label=m)
    ax.set_xlabel("household income ($1000s)")
    ax.set_title(f"party of {size}")
axes[0].set_ylabel("fitted probability")
axes[0].set_ylim(0, 0.88)
axes[1].legend(frameon=False, loc="upper left", ncol=2, columnspacing=1.0, handlelength=1.4)
fig.tight_layout()
fig.savefig(figure_path("ch36", "mode_probabilities"))
