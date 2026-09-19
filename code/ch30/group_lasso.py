"""Chapter 30, Section 3: the group lasso with factors as groups.

Simulated data: four factors (4, 3, 5 and 3 levels) and three numeric regressors; the response
depends on factor A, a little on factor B, and on the first numeric regressor. Each factor is coded
by indicator columns, centred, and orthonormalized within its group, so that the block coordinate
update is exact block soft thresholding. The group lasso is compared with the lasso applied to
the individual indicator columns.
"""
import matplotlib.pyplot as plt
import numpy as np

from regbook import COLORS, Generated, figure_path, use_book_style


# <<data>>
rng = np.random.default_rng(3003)
n = 120
levels = {"A": 4, "B": 3, "C": 5, "D": 3}
factors = {f: rng.integers(0, k, size=n) for f, k in levels.items()}
numeric = rng.normal(size=(n, 3))
effect_A = np.array([0.0, 1.0, -0.5, 0.8])
effect_B = np.array([0.0, 0.3, -0.3])
y = effect_A[factors["A"]] + effect_B[factors["B"]] + 0.6 * numeric[:, 0] + rng.normal(size=n)
y -= y.mean()
# <</data>>


# <<groups>>
def indicator_columns(codes, k, reference=0):
    """Indicator columns for the levels other than the reference level."""
    return np.column_stack([(codes == l).astype(float) for l in range(k) if l != reference])

def orthonormal_group(Z):
    """Centre the columns, then an orthonormal basis of their span (Z = Q R)."""
    Z = Z - Z.mean(axis=0)
    Q, _ = np.linalg.qr(Z)
    return Q

def build(reference=0):
    blocks = [orthonormal_group(indicator_columns(factors[f], k, reference))
              for f, k in levels.items()]
    blocks += [orthonormal_group(numeric[:, [j]]) for j in range(3)]
    return blocks

def group_lasso(blocks, y, lam, b=None, tol=1e-12, max_sweeps=100000):
    """(1/2)||y - sum_g X_g b_g||^2 + lam sum_g sqrt(p_g) ||b_g||, with X_g'X_g = I."""
    b = [np.zeros(Xg.shape[1]) for Xg in blocks] if b is None else [v.copy() for v in b]
    r = y - sum(Xg @ bg for Xg, bg in zip(blocks, b))
    for _ in range(max_sweeps):
        largest = 0.0
        for g, Xg in enumerate(blocks):
            s = Xg.T @ r + b[g]                        # block least squares on the partial residual
            w = np.sqrt(Xg.shape[1])
            norm_s = np.linalg.norm(s)
            new = max(0.0, 1 - lam * w / norm_s) * s if norm_s > 0 else 0 * s
            r += Xg @ (b[g] - new)
            largest = max(largest, np.max(np.abs(new - b[g])))
            b[g] = new
        if largest < tol:
            return b
    raise RuntimeError("block coordinate descent did not converge")
# <</groups>>

blocks = build()
names = list(levels) + ["x1", "x2", "x3"]
sizes = [Xg.shape[1] for Xg in blocks]
for Xg in blocks:
    assert np.allclose(Xg.T @ Xg, np.eye(Xg.shape[1]))


def kkt_ok(b, lam, tol=1e-8):
    r = y - sum(Xg @ bg for Xg, bg in zip(blocks, b))
    for Xg, bg in zip(blocks, b):
        grad, w = Xg.T @ r, np.sqrt(Xg.shape[1])
        if np.linalg.norm(bg) > 0:
            if not np.allclose(grad, lam * w * bg / np.linalg.norm(bg), atol=tol):
                return False
        elif np.linalg.norm(grad) > lam * w + tol:
            return False
    return True


# the zero solution holds exactly for lam >= lam_max
lam_max = max(np.linalg.norm(Xg.T @ y) / np.sqrt(Xg.shape[1]) for Xg in blocks)
assert all(np.linalg.norm(v) == 0 for v in group_lasso(blocks, y, lam_max))
assert any(np.linalg.norm(v) > 0 for v in group_lasso(blocks, y, 0.99 * lam_max))

# <<fit>>
lam = 2.0
b = group_lasso(blocks, y, lam)
for name, Xg, bg in zip(names, blocks, b):
    print(f"{name}: {Xg.shape[1]} columns, ||X_g b_g|| = {np.linalg.norm(Xg @ bg):.3f}")
# <</fit>>
assert kkt_ok(b, lam)
selected = [nm for nm, bg in zip(names, b) if np.linalg.norm(bg) > 0]
contrib = {nm: np.linalg.norm(Xg @ bg) for nm, Xg, bg in zip(names, blocks, b)}

# invariance: a different reference level gives the same fitted contribution of every group
b_other = group_lasso(build(reference=1), y, lam)
fits = [Xg @ bg for Xg, bg in zip(blocks, b)]
fits_other = [Xg @ bg for Xg, bg in zip(build(reference=1), b_other)]
assert all(np.allclose(f1, f2, atol=1e-8) for f1, f2 in zip(fits, fits_other))

# the lasso on the raw (centred) indicator columns depends on the reference level
def raw_design(reference):
    cols = [indicator_columns(factors[f], k, reference) for f, k in levels.items()] + [numeric]
    Z = np.column_stack(cols)
    return Z - Z.mean(axis=0)

def lasso_cd(Z, y, lam, tol=1e-12):
    G, c = Z.T @ Z, Z.T @ y
    bb = np.zeros(Z.shape[1])
    while True:
        largest = 0.0
        for j in range(len(bb)):
            old = bb[j]
            part = c[j] - G[j] @ bb + G[j, j] * old
            bb[j] = np.sign(part) * max(abs(part) - lam, 0.0) / G[j, j]
            largest = max(largest, abs(bb[j] - old))
        if largest < tol:
            return bb

Z0, Z1 = raw_design(0), raw_design(1)
lam_l = 6.0
fit0 = Z0 @ lasso_cd(Z0, y, lam_l)
fit1 = Z1 @ lasso_cd(Z1, y, lam_l)
coding_gap = np.linalg.norm(fit0 - fit1) / np.linalg.norm(fit0)
assert coding_gap > 0.01

# orthonormal groups (all columns mutually orthonormal): block soft thresholding
Qall, _ = np.linalg.qr(rng.normal(size=(n, 7)))
oblocks = [Qall[:, :3], Qall[:, 3:5], Qall[:, 5:]]
ob = group_lasso(oblocks, y, 0.8)
for Xg, bg in zip(oblocks, ob):
    zg = Xg.T @ y
    expected = max(0.0, 1 - 0.8 * np.sqrt(Xg.shape[1]) / np.linalg.norm(zg)) * zg
    assert np.allclose(bg, expected, atol=1e-10)

# ---- paths ------------------------------------------------------------------------------------
grid = lam_max * np.logspace(0.05, -2, 70)
gl_path = []
bb = None
for lm in grid:
    bb = group_lasso(blocks, y, lm, bb)
    gl_path.append([np.linalg.norm(Xg @ bg) for Xg, bg in zip(blocks, bb)])
gl_path = np.array(gl_path)
entry_order = [names[g] for g in np.argsort([np.argmax(gl_path[:, g] > 0) if np.any(gl_path[:, g] > 0)
                                              else 999 for g in range(len(blocks))])]
lam_max_raw = np.max(np.abs(Z0.T @ y))
grid_raw = lam_max_raw * np.logspace(0.05, -2, 70)
raw_path = []
for lm in grid_raw:
    raw_path.append(lasso_cd(Z0, y, lm))
raw_path = np.array(raw_path)
group_of = sum([[f] * (k - 1) for f, k in levels.items()], []) + ["x1", "x2", "x3"]
cols_A = [j for j, g in enumerate(group_of) if g == "A"]
entry_A = [grid_raw[np.argmax(raw_path[:, j] != 0)] for j in cols_A]
A_entry_ratio = max(entry_A) / min(entry_A)
null_cols = [j for j, g in enumerate(group_of) if g in ("C", "D")]
# the lasso admits a null-factor column (C or D) at a worse fit than the group lasso reaches with
# both null factors still excluded
null_groups = [names.index("C"), names.index("D")]
i_null_gl = min(int(np.argmax(gl_path[:, g] > 0)) if np.any(gl_path[:, g] > 0) else len(grid)
                for g in null_groups)
b_clean = group_lasso(blocks, y, grid[i_null_gl - 1])
rss_group_clean = np.sum((y - sum(Xg @ bg for Xg, bg in zip(blocks, b_clean))) ** 2)
i_null_l = int(np.argmax(np.any(raw_path[:, null_cols] != 0, axis=1)))
rss_lasso_null = np.sum((y - Z0 @ raw_path[i_null_l]) ** 2)
assert rss_lasso_null > rss_group_clean
print("A entry ratio", A_entry_ratio, "RSS: lasso at first null column", rss_lasso_null,
      "group lasso with null factors out", rss_group_clean)

gen = Generated("ch30", "group_lasso")
gen.int("n", n)
gen.num("lam", lam, 1)
gen.num("lam_max", lam_max, 3)
gen.text("selected", ", ".join(selected))
for nm in names:
    gen.num(f"contrib_{nm}", contrib[nm], 3)
gen.num("coding_gap", 100 * coding_gap, 1)
gen.num("lam_l", lam_l, 1)
gen.text("entry_order", ", ".join(entry_order[:4]))
gen.num("A_entry_ratio", A_entry_ratio, 1)
gen.num("rss_lasso_null", rss_lasso_null, 1)
gen.num("rss_group_clean", rss_group_clean, 1)
gen.int("null_cols", len(null_cols))
gen.write()

use_book_style()
fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.5))
palette = {"A": COLORS["accent"], "B": COLORS["second"], "C": COLORS["thread"], "D": COLORS["third"],
           "x1": COLORS["ink"], "x2": COLORS["muted"], "x3": COLORS["muted"]}
ax = axes[0]
for g, nm in enumerate(names):
    ax.plot(np.log10(grid), gl_path[:, g], color=palette[nm], label=nm if nm not in ("x2", "x3") else None,
            linestyle="--" if nm.startswith("x") else "-")
ax.invert_xaxis()
ax.set_xlabel(r"$\log_{10}\lambda$")
ax.set_ylabel(r"$\|\mathbf{X}_g\hat{\mathbf{b}}_g\|$")
ax.set_title("(a) group lasso: whole factors")
ax.legend(frameon=False, fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 0.3), ncol=3, columnspacing=0.8, handlelength=1.5)
ax = axes[1]
for j, nm in enumerate(group_of):
    ax.plot(np.log10(grid_raw), raw_path[:, j], color=palette[nm], linewidth=0.9,
            linestyle="--" if nm.startswith("x") else "-")
ax.axhline(0, color=COLORS["grid"], linewidth=0.6, zorder=0)
ax.invert_xaxis()
ax.set_xlabel(r"$\log_{10}\lambda$")
ax.set_ylabel("coefficient")
ax.set_title("(b) lasso on indicator columns")
fig.tight_layout()
fig.savefig(figure_path("ch30", "group_paths"))
