"""Chapter 30, Section 2: the elastic net, its grouping effect, and the naive-versus-rescaled question.

Criterion (columns centred and scaled to unit length, response centred):
    (1/2)||y - X b||^2 + lam1 ||b||_1 + (lam2/2) ||b||^2.
Solved by cyclic coordinate descent. Simulated data: three strongly correlated regressors that all
matter, and seven independent noise regressors.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


# <<solver>>
def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)

def cd_quadratic(G, c, lam1, b=None, tol=1e-12, max_sweeps=100000):
    """Minimize (1/2) b'Gb - c'b + lam1 ||b||_1 by cyclic coordinate descent."""
    p = len(c)
    b = np.zeros(p) if b is None else b.copy()
    for _ in range(max_sweeps):
        largest = 0.0
        for j in range(p):
            old = b[j]
            partial = c[j] - G[j] @ b + G[j, j] * old      # c_j - sum_{k != j} G_jk b_k
            b[j] = soft(partial, lam1) / G[j, j]
            largest = max(largest, abs(b[j] - old))
        if largest < tol:
            return b
    raise RuntimeError("coordinate descent did not converge")

def elastic_net(X, y, lam1, lam2, b=None):
    """The (naive) elastic net: (1/2)||y - Xb||^2 + lam1 ||b||_1 + (lam2/2)||b||^2."""
    G = X.T @ X + lam2 * np.eye(X.shape[1])
    return cd_quadratic(G, X.T @ y, lam1, b)
# <</solver>>


# <<data>>
rng = np.random.default_rng(3002)
n, p = 50, 10
z = rng.normal(size=n)
X = rng.normal(size=(n, p))
X[:, :3] = z[:, None] + 0.3 * rng.normal(size=(n, 3))   # three near-copies of one signal
X -= X.mean(axis=0)
X /= np.linalg.norm(X, axis=0)                          # centred, unit-length columns
beta = np.r_[2.0, 2.0, 2.0, np.zeros(p - 3)]
y = X @ beta + 0.5 * rng.normal(size=n)
y -= y.mean()
# <</data>>

R = X.T @ X
rho_12, rho_13, rho_23 = R[0, 1], R[0, 2], R[1, 2]
assert min(rho_12, rho_13, rho_23) > 0.85


def kkt_ok(b, lam1, lam2, tol=1e-8):
    g = X.T @ (y - X @ b) - lam2 * b                    # must lie in lam1 * subdifferential
    act = b != 0
    return (np.allclose(g[act], lam1 * np.sign(b[act]), atol=tol)
            and np.all(np.abs(g[~act]) <= lam1 + tol))


# <<compare>>
lam1 = 0.8
b_lasso = elastic_net(X, y, lam1, 0.0)
b_naive = elastic_net(X, y, lam1, 1.0)
print("lasso         ", np.round(b_lasso[:5], 3))
print("naive EN      ", np.round(b_naive[:5], 3))
print("rescaled EN   ", np.round(2.0 * b_naive[:5], 3))     # (1 + lam2) * naive
# <</compare>>

lam2 = 1.0
assert kkt_ok(b_lasso, lam1, 0.0) and kkt_ok(b_naive, lam1, lam2)
assert np.sum(b_lasso[:3] != 0) < 3                     # the lasso drops a member of the group
assert np.all(b_naive[:3] > 0)                          # the elastic net keeps all three

# grouping bound |b_i - b_j| <= ||y|| sqrt(2 (1 - rho_ij)) / lam2 for same-sign nonzero pairs
worst = 0.0
for i in range(p):
    for j in range(i + 1, p):
        if b_naive[i] * b_naive[j] > 0:
            bound = np.linalg.norm(y) * np.sqrt(2 * (1 - R[i, j])) / lam2
            sharper = np.linalg.norm(y - X @ b_naive) * np.linalg.norm(X[:, i] - X[:, j]) / lam2
            assert abs(b_naive[i] - b_naive[j]) <= sharper + 1e-10 <= bound + 1e-10
            worst = max(worst, abs(b_naive[i] - b_naive[j]) / bound)
spread_naive = np.ptp(b_naive[:3])
bound_12 = np.linalg.norm(y) * np.sqrt(2 * (1 - rho_12)) / lam2
sharp_12 = np.linalg.norm(y - X @ b_naive) * np.sqrt(2 * (1 - rho_12)) / lam2
diff_12 = abs(b_naive[0] - b_naive[1])

# augmented data: the naive elastic net is a lasso on (X*, y*)
Xs = np.vstack([X, np.sqrt(lam2) * np.eye(p)])
ys = np.r_[y, np.zeros(p)]
b_aug = cd_quadratic(Xs.T @ Xs, Xs.T @ ys, lam1)
assert np.allclose(b_aug, b_naive, atol=1e-9)

# rescaled elastic net = lasso with the Gram matrix shrunk towards the identity
G_shrunk = (R + lam2 * np.eye(p)) / (1 + lam2)
b_resc = cd_quadratic(G_shrunk, X.T @ y, lam1)
assert np.allclose(b_resc, (1 + lam2) * b_naive, atol=1e-9)

# identical columns: the elastic net splits equally, the lasso's solutions are not unique
Xd = np.column_stack([X[:, 0], X[:, 0], X[:, 3]])
bd = elastic_net(Xd, y, 0.3, 0.2)
assert np.isclose(bd[0], bd[1])
bl = elastic_net(Xd, y, 0.3, 0.0)
other = bl.copy()
other[0], other[1] = bl[0] + bl[1], 0.0                 # move all the weight to one copy
obj = lambda b: 0.5 * np.sum((y - Xd @ b) ** 2) + 0.3 * np.sum(np.abs(b))
assert np.isclose(obj(bl), obj(other))

# orthonormal design: naive EN = soft / (1 + lam2); rescaled EN = lasso
Q, _ = np.linalg.qr(rng.normal(size=(n, 4)))
zq = Q.T @ y
assert np.allclose(elastic_net(Q, y, 0.3, 0.7), soft(zq, 0.3) / 1.7, atol=1e-10)

# ---- paths ------------------------------------------------------------------------------------
lam_max = np.max(np.abs(X.T @ y))
grid = lam_max * np.logspace(np.log10(2.0), -2.5, 90)   # the elastic net enters at 2 lam_max
alpha = 0.5                                             # lam1 = alpha lam, lam2 = (1 - alpha) lam
paths = {"lasso": [], "enet": []}
b1 = np.zeros(p)
b2 = np.zeros(p)
for lam in grid:
    b1 = elastic_net(X, y, lam, 0.0, b1)
    b2 = elastic_net(X, y, alpha * lam, (1 - alpha) * lam, b2)
    paths["lasso"].append(b1.copy())
    paths["enet"].append(b2.copy())
paths = {k: np.array(v) for k, v in paths.items()}
# at every lambda on the grid, the elastic net keeps the three correlated columns together
first_in = [np.argmax(paths["enet"][:, j] != 0) for j in range(3)]
first_in_lasso = [np.argmax(paths["lasso"][:, j] != 0) for j in range(3)]
assert max(first_in) - min(first_in) < max(first_in_lasso) - min(first_in_lasso)
entry_lasso = grid[first_in_lasso]
entry_enet = grid[first_in]
lasso_entry_ratio = max(entry_lasso) / min(entry_lasso)       # first partner in / last partner in
enet_entry_ratio = max(entry_enet) / min(entry_enet)
print("entry ratios", lasso_entry_ratio, enet_entry_ratio, first_in, first_in_lasso)
noise_max = np.max(np.abs(b_naive[3:]))

gen = Generated("ch30", "elastic_net")
gen.int("n", n)
gen.int("p", p)
gen.num("rho_min", min(rho_12, rho_13, rho_23), 2)
gen.num("rho_max", max(rho_12, rho_13, rho_23), 2)
gen.num("lam1", lam1, 1)
for j in range(3):
    gen.num(f"lasso_{j + 1}", b_lasso[j], 3)
    gen.num(f"naive_{j + 1}", b_naive[j], 3)
    gen.num(f"resc_{j + 1}", 2 * b_naive[j], 3)
gen.int("lasso_active", int(np.sum(b_lasso != 0)))
gen.int("naive_active", int(np.sum(b_naive != 0)))
gen.num("spread_naive", spread_naive, 3)
gen.num("bound_12", bound_12, 3)
gen.num("sharp_12", sharp_12, 3)
gen.num("diff_12", diff_12, 3)
gen.num("rho_12", rho_12, 3)
gen.num("rnorm", np.linalg.norm(y - X @ b_naive), 3)
gen.num("ynorm", np.linalg.norm(y), 3)
gen.num("worst", worst, 3)
gen.num("lasso_entry_ratio", lasso_entry_ratio, 0)
gen.num("enet_entry_ratio", enet_entry_ratio, 2)
gen.num("noise_max", noise_max, 3)
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.4), sharey=True)
x_axis = np.log10(grid)
for ax, key, title in [(axes[0], "lasso", "(a) lasso"),
                       (axes[1], "enet", r"(b) elastic net, $\lambda_1=\lambda_2=\lambda/2$")]:
    for j in range(p):
        color = [COLORS["accent"], COLORS["second"], COLORS["third"]][j] if j < 3 else COLORS["muted"]
        ax.plot(x_axis, paths[key][:, j], color=color, linewidth=1.2 if j < 3 else 0.7)
    ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
    ax.invert_xaxis()
    ax.set_xlabel(r"$\log_{10}\lambda$")
    ax.set_title(title)
axes[0].set_ylabel("coefficient")
fig.tight_layout()
fig.savefig(figure_path("ch30", "en_paths"))
