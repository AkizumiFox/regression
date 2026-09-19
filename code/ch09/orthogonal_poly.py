"""Chapter 9, Section 4: orthogonal polynomials split a trend sum of squares into
order-free linear, quadratic, cubic and quartic pieces.

Annual flow of the Nile at Aswan, 1871-1970 (10^8 m^3). Public-domain data shipped
with statsmodels (statsmodels.datasets.nile).
"""
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from regbook import COLORS, Generated, figure_path, use_book_style

nile = sm.datasets.nile.load_pandas().data
year = nile["year"].to_numpy()
y = nile["volume"].to_numpy()
n = len(y)

# <<basis>>
t = year - year.mean()                              # centring only helps rounding
V = np.vander(t / t.max(), 5, increasing=True)      # 1, t, t^2, t^3, t^4 (scaled)
Q, R = np.linalg.qr(V)                              # Gram-Schmidt in matrix form
Q = Q * np.sign(np.diag(R))                         # make each q_k agree in sign with t^k
coef = Q.T @ y                                      # coordinates of y on q_0, ..., q_4
seq_ss = coef[1:] ** 2                              # sequential SS: linear, ..., quartic
for k, s in enumerate(seq_ss, start=1):
    print(f"degree {k}: sequential SS = {s:10.0f}")
print(f"residual after degree 4: {y @ y - np.sum(coef ** 2):10.0f}")
# <</basis>>

# checks: q_k are orthonormal, nested spans equal those of the powers, and the
# sequential SS for degree k is SSE_{k-1} - SSE_k whatever basis is used
assert np.allclose(Q.T @ Q, np.eye(5))


def sse(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return np.sum((y - Z @ b) ** 2)


for k in range(1, 5):
    assert np.isclose(seq_ss[k - 1], sse(V[:, :k]) - sse(V[:, :k + 1]))
# order-free: entering the cubic column before the linear one changes nothing
assert np.isclose(sse(Q[:, [0, 3]]) - sse(Q[:, [0, 3, 1]]), seq_ss[0])
# with raw powers the order matters
raw_lin_first = sse(V[:, :1]) - sse(V[:, :2])
raw_lin_after_cubic = sse(V[:, [0, 2, 3]]) - sse(V[:, :4])
assert not np.isclose(raw_lin_first, raw_lin_after_cubic)

# coefficients: stable in the orthogonal basis, not in the raw one
lin_raw = [np.linalg.lstsq(V[:, :k + 1], y, rcond=None)[0][1] for k in range(1, 5)]
lin_orth = [np.linalg.lstsq(Q[:, :k + 1], y, rcond=None)[0][1] for k in range(1, 5)]
assert np.allclose(lin_orth, coef[1])
sst = np.sum((y - y.mean()) ** 2)
sse1 = sse(V[:, :2])
mse3 = sse(V[:, :4]) / (n - 4)

gen = Generated("ch09", "orthogonal_poly", prefix="poly")
gen.int("n", n)
for k, s in enumerate(seq_ss, start=1):
    gen.num(f"ss{k}", s, 0)
gen.num("sst", sst, 0)
gen.num("res4", y @ y - np.sum(coef ** 2), 0)
gen.num("share1", seq_ss[0] / sst, 3)
gen.num("share3", seq_ss[:3].sum() / sst, 3)
gen.num("lin_after_cubic", raw_lin_after_cubic, 0)
gen.num("F2", seq_ss[1] / mse3, 2)
gen.num("F3", seq_ss[2] / mse3, 2)
gen.num("mse3", mse3, 0)
for k in range(4):
    gen.num(f"linraw{k + 1}", lin_raw[k], 1)
gen.num("linorth", coef[1], 1)
gen.write()

# ---- figure -------------------------------------------------------------------
use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4))
ax = axes[0]
colors = [COLORS["accent"], COLORS["second"], COLORS["third"]]
for k in range(1, 4):
    ax.plot(year, Q[:, k], color=colors[k - 1], label=f"$q_{k}$")
ax.axhline(0, color=COLORS["grid"], lw=0.6, zorder=0)
ax.set_xlabel("year")
ax.set_title("(a) orthonormal polynomials")
ax.set_ylim(-0.3, 0.38)
ax.legend(frameon=False, ncol=3, loc="upper center", fontsize=7, handlelength=1.2)
ax = axes[1]
ax.scatter(year, y, s=6, color=COLORS["muted"], linewidths=0)
for k in range(1, 4):
    ax.plot(year, Q[:, :k + 1] @ coef[:k + 1], color=colors[k - 1], label=f"degree {k}")
ax.set_xlabel("year")
ax.set_ylabel(r"flow ($10^8\,$m$^3$)")
ax.set_title("(b) Nile flow and nested fits")
ax.set_ylim(400, 1600)
ax.legend(frameon=False, ncol=3, loc="upper center", fontsize=7, handlelength=1.2, columnspacing=0.8)
fig.tight_layout()
fig.savefig(figure_path("ch09", "orthogonal_poly"))
