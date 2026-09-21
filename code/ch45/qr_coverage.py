"""Chapter 45, Section 3: coverage of three confidence intervals for a quantile slope.

Data are generated from a location-scale model whose errors are independent but far
from identically distributed, so that the textbook iid formula is wrong by
construction. The table in the text reports how often each nominal 95% interval for
the slope of the tau-quantile covers the truth.
"""
import numpy as np
from scipy import stats
from scipy.optimize import linprog

from regbook import Generated

REPS, BOOT_REPS, B, N = 800, 150, 99, 150
TAUS = (0.50, 0.90)


def qreg(X, y, tau):
    """Quantile regression by linear programming (as in Section 45.2)."""
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]


def hall_sheather(n, tau, alpha=0.05):
    z, za = stats.norm.ppf(tau), stats.norm.ppf(1 - alpha / 2)
    return n ** (-1 / 3) * za ** (2 / 3) * (1.5 * stats.norm.pdf(z) ** 2 / (2 * z ** 2 + 1)) ** (1 / 3)


def sparsity(resid, tau, n):
    h = min(hall_sheather(n, tau), tau, 1 - tau)
    hi, lo = np.quantile(resid, [min(tau + h, 1.0), max(tau - h, 0.0)])
    return (hi - lo) / (2 * h)


def ses(X, y, tau):
    """Standard errors for the slope: the iid formula and Powell's sandwich."""
    beta = qreg(X, y, tau)
    n = len(y)
    r = y - X @ beta
    s = sparsity(r, tau, n)
    iid = np.sqrt(tau * (1 - tau) * s ** 2 * np.diag(np.linalg.inv(X.T @ X)))
    h = s * hall_sheather(n, tau)
    dens = np.exp(-0.5 * (r / h) ** 2) / (h * np.sqrt(2 * np.pi))
    D1 = np.linalg.inv(X.T @ (dens[:, None] * X))
    sand = np.sqrt(tau * (1 - tau) * np.diag(D1 @ (X.T @ X) @ D1))
    return beta, iid, sand


def sample(rng, n=N):
    """y = 1 + x + (0.2 + 1.2 x) e with x uniform on (0, 2) and e standard normal."""
    x = rng.uniform(0, 2, n)
    y = 1 + x + (0.2 + 1.2 * x) * rng.standard_normal(n)
    return np.column_stack([np.ones(n), x]), y


def true_slope(tau):
    """The tau-quantile of the model above is 1 + x + (0.2 + 1.2 x) z_tau."""
    return 1 + 1.2 * stats.norm.ppf(tau)


rng = np.random.default_rng(20451)
hits = {tau: {"iid": 0, "sandwich": 0} for tau in TAUS}
widths = {tau: {"iid": [], "sandwich": []} for tau in TAUS}
for _ in range(REPS):
    X, y = sample(rng)
    for tau in TAUS:
        beta, iid, sand = ses(X, y, tau)
        for name, se in (("iid", iid), ("sandwich", sand)):
            hits[tau][name] += abs(beta[1] - true_slope(tau)) <= 1.96 * se[1]
            widths[tau][name].append(2 * 1.96 * se[1])

boot_hits = {tau: 0 for tau in TAUS}
boot_width = {tau: [] for tau in TAUS}
for _ in range(BOOT_REPS):
    X, y = sample(rng)
    idx = rng.integers(0, N, (B, N))
    for tau in TAUS:
        beta = qreg(X, y, tau)
        draws = np.array([qreg(X[i], y[i], tau) for i in idx])
        se = draws[:, 1].std(ddof=1)
        boot_hits[tau] += abs(beta[1] - true_slope(tau)) <= 1.96 * se
        boot_width[tau].append(2 * 1.96 * se)

cov = {}
for tau in TAUS:
    for name in ("iid", "sandwich"):
        cov[(tau, name)] = hits[tau][name] / REPS
    cov[(tau, "bootstrap")] = boot_hits[tau] / BOOT_REPS
    print(f"tau = {tau}:  iid {cov[(tau, 'iid')]:.3f}   "
          f"sandwich {cov[(tau, 'sandwich')]:.3f}   bootstrap {cov[(tau, 'bootstrap')]:.3f}")

assert cov[(0.50, "iid")] < 0.92                 # the iid formula undercovers
assert cov[(0.50, "sandwich")] > 0.92
assert cov[(0.50, "bootstrap")] > 0.92

gen = Generated("ch45", "qr_coverage")
gen.int("reps", REPS)
gen.int("boot_reps", BOOT_REPS)
gen.int("B", B)
gen.int("n", N)
for tau in TAUS:
    key = f"{round(100 * tau):02d}"
    for name in ("iid", "sandwich", "bootstrap"):
        gen.num(f"cov_{name}_{key}", cov[(tau, name)], 3)
    gen.num(f"width_iid_{key}", np.mean(widths[tau]["iid"]), 3)
    gen.num(f"width_sand_{key}", np.mean(widths[tau]["sandwich"]), 3)
    gen.num(f"width_boot_{key}", np.mean(boot_width[tau]), 3)
    gen.num(f"true_{key}", true_slope(tau), 3)
gen.write()
