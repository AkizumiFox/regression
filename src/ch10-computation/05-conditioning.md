# Conditioning and perturbation bounds

[Sections 10.2](02-cholesky.html) to [10.4](04-svd.html) were about algorithms. This section is about the problem itself: if the data
\( (\X,\y) \) change a little, how much can the least squares solution and the residuals change? The
answer does not depend on the algorithm. Combined with backward stability through
@eq-cmp-rule, it predicts the accuracy of any computed solution.

## The perturbation theorem

Throughout this section \( \X \) is \( n\times p \) with full column rank, \( \hbeta=\X^+\y \) is its least squares
solution (assumed nonzero) and \( \he=\y-\X\hbeta \) its residual. Write
\( \kappa=\kappa(\X)=\sigma_1/\sigma_p \). The perturbed data are \( \tilde{\X}=\X+\bm E \) and
\( \tilde{\y}=\y+\bm f \), with
\[
\norm{\bm E}_2\le\epsilon\norm{\X}_2,\qquad\norm{\bm f}\le\epsilon\norm{\y}.
\]{#eq-cmp-normwise}

The solution and residual for the perturbed data are \( \tilde{\bb} \) and \( \tilde{\be} \). A single
number measures how badly the model fits relative to the size of the problem:
\[
\eta=\frac{\norm{\he}}{\norm{\X}_2\norm{\hbeta}} .
\]
Since \( \norm{\X\hbeta}\le\norm{\X}_2\norm{\hbeta} \), \( \eta\le\tan\theta \), where \( \theta \) is the
angle between \( \y \) and \( \C(\X) \). A model that fits well has small \( \eta \).

::: {#thm-cmp-perturbation}
[Perturbation of least squares]

Under @eq-cmp-normwise with \( \kappa\epsilon<1 \), the matrix \( \tilde{\X} \) has full column rank, and
\[
\frac{\norm{\tilde{\bb}-\hbeta}}{\norm{\hbeta}}\le
\frac{\kappa\epsilon}{1-\kappa\epsilon}\,(2+\eta)+\frac{\kappa^2\epsilon}{(1-\kappa\epsilon)^2}\,\eta ,
\]{#eq-cmp-perturbation}

\[
\frac{\norm{\tilde{\be}-\he}}{\norm{\y}}\le\epsilon\Bigl(1+\frac{2\kappa}{1-\kappa\epsilon}\Bigr).
\]{#eq-cmp-residual-perturbation}

To first order in \( \epsilon \), the coefficient bound is \( \epsilon\bigl(2\kappa+\kappa(\kappa+1)\eta\bigr) \).
:::

::: {.proof}
By @prp-cmp-weyl, \( \sigma_p(\tilde{\X})\ge\sigma_p-\norm{\bm E}_2\ge\sigma_p(1-\kappa\epsilon)>0 \), so
\( \tilde{\X} \) has full column rank and \( \tilde{\X}^+\tilde{\X}=\I \). Hence
\[
\tilde{\bb}-\hbeta=\tilde{\X}^+\tilde{\y}-\tilde{\X}^+\tilde{\X}\hbeta=\tilde{\X}^+(\tilde{\y}-\tilde{\X}\hbeta)
=\tilde{\X}^+(\he+\bm f-\bm E\hbeta).
\]
The residual term simplifies because \( \X\T\he=\bzero \):
\( \tilde{\X}^+\he=(\tilde{\X}\T\tilde{\X})^{-1}(\X+\bm E)\T\he=(\tilde{\X}\T\tilde{\X})^{-1}\bm E\T\he \).
So
\[
\tilde{\bb}-\hbeta=\tilde{\X}^+(\bm f-\bm E\hbeta)+(\tilde{\X}\T\tilde{\X})^{-1}\bm E\T\he .
\]{#eq-cmp-exact-perturbation}

Now \( \norm{\tilde{\X}^+}_2=1/\sigma_p(\tilde{\X}) \) and \( \norm{(\tilde{\X}\T\tilde{\X})^{-1}}_2=1/\sigma_p(\tilde{\X})^2 \)
(@thm-mat-svd). Also \( \norm{\y}\le\norm{\X\hbeta}+\norm{\he}\le\sigma_1\norm{\hbeta}+\norm{\he} \), so
\( \norm{\bm f-\bm E\hbeta}\le\epsilon(\norm{\y}+\sigma_1\norm{\hbeta})\le\epsilon\sigma_1\norm{\hbeta}(2+\eta) \).
Therefore
\[
\norm{\tilde{\bb}-\hbeta}\le\frac{\epsilon\sigma_1\norm{\hbeta}(2+\eta)}{\sigma_p(1-\kappa\epsilon)}
+\frac{\epsilon\sigma_1\norm{\he}}{\sigma_p^2(1-\kappa\epsilon)^2} .
\]
Dividing by \( \norm{\hbeta} \) and using \( \norm{\he}=\eta\sigma_1\norm{\hbeta} \) gives @eq-cmp-perturbation.

For the residual, \( \tilde{\be}=(\I-\tilde{\M})\tilde{\y} \), where \( \tilde{\M} \) projects onto \( \C(\tilde{\X}) \).
Write \( \tilde{\y}=\tilde{\X}\hbeta+(\bm f-\bm E\hbeta)+\he \). The first term is annihilated by
\( \I-\tilde{\M} \), and \( \tilde{\M}\he=(\tilde{\X}^+)\T\tilde{\X}\T\he=(\tilde{\X}^+)\T\bm E\T\he \) as before.
So
\[
\tilde{\be}-\he=(\I-\tilde{\M})(\bm f-\bm E\hbeta)-(\tilde{\X}^+)\T\bm E\T\he .
\]
The first term has norm at most \( \epsilon(\norm{\y}+\sigma_1\norm{\hbeta})\le\epsilon(1+\kappa)\norm{\y} \),
because \( \norm{\hbeta}=\norm{\X^+\y}\le\norm{\y}/\sigma_p \). The second has norm at most
\( \epsilon\sigma_1\norm{\he}/(\sigma_p(1-\kappa\epsilon))\le\epsilon\kappa\norm{\y}/(1-\kappa\epsilon) \).
Adding and using \( 1+\kappa+\kappa/(1-\kappa\epsilon)\le1+2\kappa/(1-\kappa\epsilon) \) gives the
bound @eq-cmp-residual-perturbation.
:::

The exact identity @eq-cmp-exact-perturbation shows where each term comes from. Perturbing the data
along \( \C(\X) \) moves the solution through \( \tilde{\X}^+ \), with amplification \( 1/\sigma_p \). This is the
familiar \( \kappa \) of linear systems (@eq-mat-perturbation). The second term is new. A perturbation
\( \bm E \) that tilts a column of \( \X \) towards the residual direction changes the subspace, and the
residual, however large, then "leaks" into the solution through \( (\tilde{\X}\T\tilde{\X})^{-1} \), with
amplification \( 1/\sigma_p^2 \). The residuals themselves, and hence the fitted values, are only
\( \kappa \)-sensitive.

The \( \kappa^2 \) term is not an artefact of the proof.

::: {#exm-cmp-kappa-squared}
[The squared condition number is attained]

Let \( \delta=10^{-3} \) and
\[
\X=\begin{pmatrix}1&0\\0&\delta\\0&0\end{pmatrix},\qquad \y=\begin{pmatrix}0\\\delta\\1\end{pmatrix},
\]
so that \( \hbeta=(0,1)\T \), \( \he=(0,0,1)\T \), \( \kappa=1000 \) and \( \eta=1 \).
Perturb only the \( (3,2) \) entry of \( \X \), by \( \epsilon=10^{-10} \). This tilts the weak second column
towards the residual. Its coefficient becomes \( (\delta^2+\epsilon)/(\delta^2+\epsilon^2) \), and the relative
change in \( \hbeta \) is \( 1.00\times 10^{-4} \), equal to \( \kappa^2\epsilon=1.00\times 10^{-4} \). A perturbation of
relative size \( 10^{-10} \) has changed the answer in the fourth digit.
:::

```{.python .run #cell-perturbation-attained}
import numpy as np
import statsmodels.api as sm

def ls(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]

delta, eps = 1e-3, 1e-10
X = np.array([[1.0, 0.0], [0.0, delta], [0.0, 0.0]])
y = np.array([0.0, delta, 1.0])            # b = (0, 1), residual e = (0, 0, 1)
E = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, eps]])   # ||E|| = eps ||X||
b = ls(X, y)
b_tilde = ls(X + E, y)
print("kappa =", np.linalg.cond(X))
print("relative change in b:", np.linalg.norm(b_tilde - b) / np.linalg.norm(b))
print("kappa^2 * eps       :", np.linalg.cond(X) ** 2 * eps)
```

The script also draws 300 random problems with random perturbations and checks
@eq-cmp-perturbation on each. The worst ratio of the actual change to the bound is
\( 0.31 \).

## What the theorem says about algorithms

A backward stable algorithm (Householder QR, Givens QR, modified Gram–Schmidt on the augmented
matrix, the SVD) returns the exact solution for data perturbed as in @eq-cmp-normwise with
\( \epsilon\approx c\,u \), where \( c \) grows modestly with the dimensions. Its relative error is therefore
about
\[
c\,u\,(\kappa+\kappa^2\eta).
\]
The normal equations incur a relative error of about \( c\,u\,\kappa^2 \) *whatever* \( \eta \) is,
because forming \( \X\T\X \) already commits it (@prp-cmp-gram-rounding). Two regimes follow:

- **Small residual** (\( \eta\ll1/\kappa \)): QR gains a full factor \( \kappa \) over the normal equations.
  This is the regime of the experiment in [Section 6.10](../ch06-projections/10-computation.html), where
  \( \y=\X\bbeta \) exactly.
- **Large residual** (\( \eta \) of order one): every method, however stable, loses about
  \( 2\log_{10}\kappa \) digits, because the problem is that sensitive. The normal equations are then
  worse only by a modest factor.

[Figure 10.5.1](#fig-cmp-residual-size) shows both regimes on test problems with prescribed condition
numbers. With zero residual and \( \kappa=10^6 \), the normal equations give relative error
\( 1.1\times 10^{-5} \) and QR \( 7.8\times 10^{-12} \). With \( \eta=1 \), that is \( \norm{\he}=\norm{\X}_2\norm{\bbeta} \), QR gives
\( 1.9\times 10^{-7} \) and the SVD \( 1.9\times 10^{-7} \). The \( \kappa^2 \) term has caught up
with the stable methods.

::: {when-format="html"}
![**Figure 10.5.1.** Relative error of the coefficients computed by the normal equations (Cholesky),
Householder QR and the SVD, against \( \kappa(\X) \). (a) Consistent data, \( \he=\bzero \): QR and SVD
track \( \kappa u \), the normal equations track \( \kappa^2u \). (b) A large residual, \( \eta=1 \): all
three track \( \kappa^2u \), as @thm-cmp-perturbation predicts.](residual_size.svg){#fig-cmp-residual-size width=100%}
:::

::: {when-format="pdf"}
![Relative error of the coefficients computed by the normal equations (Cholesky),
Householder QR and the SVD, against \( \kappa(\X) \). (a) Consistent data, \( \he=\bzero \): QR and SVD
track \( \kappa u \), the normal equations track \( \kappa^2u \). (b) A large residual, \( \eta=1 \): all
three track \( \kappa^2u \), as @thm-cmp-perturbation predicts.](residual_size.pdf){width=100%}
:::

```{.python .run #cell-perturbation-solvers}
def by_cholesky(X, y):
    L = np.linalg.cholesky(X.T @ X)
    return np.linalg.solve(L.T, np.linalg.solve(L, X.T @ y))

def by_qr(X, y):
    Q, R = np.linalg.qr(X)
    return np.linalg.solve(R, Q.T @ y)

def by_svd(X, y):
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    return Vt.T @ ((U.T @ y) / s)
```

In regression \( \eta \) is rarely tiny, but QR keeps an advantage: its \( \kappa^2\eta \) term is a property of
the problem, not of the method.

## Column scaling

The normwise condition number depends on the units of the columns: measuring one regressor in
millimetres instead of metres can multiply \( \kappa(\X) \) by up to \( 1000 \) (@exr-mat-cond-scaling),
although the fit is unchanged. Equalizing the column lengths is essentially the best diagonal scaling.

::: {#prp-cmp-van-der-sluis}
[Equilibrating the columns]

Let \( \X \) have full column rank and let \( \bD_0=\diag(1/\norm{\x_1},\dots,1/\norm{\x_p}) \). For every
nonsingular diagonal \( \bD \),
\[
\kappa(\X\bD_0)\le\sqrt p\,\kappa(\X\bD).
\]
:::

::: {.proof}
The columns of \( \X\bD_0 \) have unit length, so \( \sigma_1(\X\bD_0)\le\norm{\X\bD_0}_F=\sqrt p \). For the
smallest singular value, let \( \bz \) be any nonzero vector and put \( \bw=\bD^{-1}\bD_0\bz \). Then
\( \norm{\X\bD_0\bz}=\norm{\X\bD\bw}\ge\sigma_p(\X\bD)\norm{\bw} \) and
\( \norm{\bw}\ge\min_j\bigl(1/(\lvert d_j\rvert\norm{\x_j})\bigr)\norm{\bz} \). The largest column length of \( \X\bD \),
\( \max_j\lvert d_j\rvert\norm{\x_j} \), is at most \( \sigma_1(\X\bD) \), because the columns are images of unit
vectors. So \( \sigma_p(\X\bD_0)\ge\sigma_p(\X\bD)/\sigma_1(\X\bD)=1/\kappa(\X\bD) \), and
\( \kappa(\X\bD_0)\le\sqrt p\,\kappa(\X\bD) \).
:::

This is due to van der Sluis (1969). What makes it more than a remark about units is the form of the
backward errors. Householder QR has a *columnwise* backward error (@thm-cmp-householder(c)). Rescaling
the columns rescales the perturbations with them, so the computed solution is as accurate as if the
columns had been equilibrated first. The rounding errors in \( \X\T\X \) are entrywise small relative to
\( \lvert\X\rvert\T\lvert\X\rvert \) (@prp-cmp-gram-rounding), which has the same property. Both methods therefore
behave according to the condition number of the *scaled* matrix, \( \kappa(\X\bD_0) \), whatever the units.
The rule of thumb should be applied to \( \kappa(\X\bD_0) \), not to \( \kappa(\X) \) as recorded.

::: {#exm-cmp-longley}
[Longley's data, exactly]

Longley (1967) built his macroeconomic data to test regression programs. With an intercept and six regressors, \( \kappa(\X)=4.86\times 10^{9} \)
([Section 6.10](../ch06-projections/10-computation.html); the data are described in @exm-rv-longley). The data are published as decimals, so the exact least squares solution can be
computed in rational arithmetic. The script does so and counts the correct significant digits of the
worst coefficient for each method:

| Method | Correct digits |
|---|---|
| Normal equations (Cholesky) | 7.2 |
| Householder QR | 10.9 |
| SVD | 10.8 |
| statsmodels `OLS` | 10.8 |
| Normal equations after scaling the columns to unit length | 7.0 |

The raw \( \kappa \) predicts about \( 6 \) digits from QR and none from the normal equations
([Section 6.10](../ch06-projections/10-computation.html)). Both predictions are far too pessimistic. After
scaling the columns to unit length, \( \kappa(\X\bD_0)=4.33\times 10^{4} \). Then \( \kappa u\approx5\times10^{-12} \)
predicts about eleven digits from QR, and \( \kappa^2u\approx2\times10^{-7} \) about seven from the normal
equations. That is what is observed. Scaling the columns explicitly before forming the normal equations
changes nothing, as the argument above predicts. Centring as well reduces the condition number to
\( 1.11\times 10^{2} \). That is a reparameterization rather than a scaling (@thm-proj-reparam), and what
remains is genuine collinearity among the regressors (Chapter 26).
:::

```{.python .run #cell-perturbation-longley}
lon = sm.datasets.longley.load_pandas()
XL = sm.add_constant(lon.exog.to_numpy())
yL = lon.endog.to_numpy()
names = ["const"] + list(lon.exog.columns)
print("kappa(X) as recorded:", f"{np.linalg.cond(XL):.3e}")
```

## Numerical error and statistical error

The perturbation theorem applies to any perturbation of the data, including the rounding of recorded
values, which is far above \( u \). Longley's GNP deflator is given to one decimal place and the other series to whole units.
Perturbing each regressor entry (the year is exact) uniformly within half a unit of its last recorded digit, \( 1000 \) times, the
coefficients vary with standard deviations between \( 0.005 \) and \( 0.024 \) of their
standard errors. The numerical error of the normal equations, the least accurate method, is at most
\( 4\times 10^{-8} \) standard errors. Sampling error dominates, then comes the rounding of the published
data, and the arithmetic is last by far. Beaton, Rubin and Barone (1976) made this point about
Longley's data. A problem can be numerically delicate and statistically ordinary, or the reverse.

```{.python .run #cell-perturbation-rounding}
rng = np.random.default_rng(106)
half_unit = np.array([0.0, 0.05, 0.5, 0.5, 0.5, 0.5, 0.0])   # half a unit in the last digit
fitL = sm.OLS(yL, XL).fit()                                  # year and constant are exact
reps = 1000
shifts = []
for _ in range(reps):
    Xp = XL + rng.uniform(-1, 1, size=XL.shape) * half_unit
    shifts.append(ls(Xp, yL) - fitL.params)
spread = np.std(shifts, axis=0)
ratio = spread / fitL.bse
for name, r in zip(names, ratio):
    print(f"{name:8s} sd of coefficient under rounding / standard error = {r:.2f}")
```

::: {.idea}
Least squares has a \( \kappa \) sensitivity, and a \( \kappa^2\eta \) one when the residual is large.
Stable algorithms add nothing beyond this; the normal equations add \( \kappa^2 \) regardless. Measure
\( \kappa \) after equilibrating the columns.
:::

## Exercises

### A. Check your understanding

::: {#exr-cmp-kappa-squared-exact}
[A1]

In @exm-cmp-kappa-squared, find \( \tilde{\bb} \) exactly and show that the relative change is
\( \epsilon(1-\epsilon)/(\delta^2+\epsilon^2) \). For which \( \epsilon \) is it close to \( \kappa^2\epsilon \)?
:::

### B. Practice

::: {#exr-cmp-fitted-perturbation}
[B1]

Deduce from @eq-cmp-residual-perturbation that the fitted values satisfy
\( \norm{\tilde{\X}\tilde{\bb}-\X\hbeta}\le\epsilon\bigl(2+2\kappa/(1-\kappa\epsilon)\bigr)\norm{\y} \). Why are fitted
values never \( \kappa^2 \)-sensitive, even though coefficients can be?
:::

::: {.solution}
\( \tilde{\X}\tilde{\bb}=\tilde{\y}-\tilde{\be} \) and \( \X\hbeta=\y-\he \), so the difference is
\( \bm f-(\tilde{\be}-\he) \), with norm at most
\( \epsilon\norm{\y}+\epsilon(1+2\kappa/(1-\kappa\epsilon))\norm{\y} \). The \( \kappa^2 \) term in the coefficients
comes from \( (\tilde{\X}\T\tilde{\X})^{-1}\bm E\T\he \), a change of coordinates *within* the slightly
rotated column space. Multiplying by \( \tilde{\X} \) turns the factor \( 1/\sigma_p^2 \) into \( 1/\sigma_p \). A
direction that the columns barely span is a direction in which large coefficient changes produce only
small changes in the fit.
:::

::: {#exr-cmp-consistent-bound}
[B2]

Show that for consistent data (\( \he=\bzero \)) the bound reduces to
\( 2\kappa\epsilon/(1-\kappa\epsilon) \). Construct a perturbation \( \bm f \) of \( \y \) alone that changes \( \hbeta \) by a
relative amount \( \kappa\epsilon \), so that the bound is attained up to a factor of two.
:::

### C. Going deeper

::: {#exr-cmp-sharpness}
[C1]

Take data with \( \hbeta=\bv_p \) and \( \he\ne\bzero \), and the perturbation
\( \bm E=\epsilon\sigma_1\,\he\bv_p\T/\norm{\he} \), with \( \bm f=\bzero \). Use @eq-cmp-exact-perturbation to
show that, to first order in \( \epsilon \), \( \tilde{\bb}-\hbeta=\epsilon\sigma_1\norm{\he}\sigma_p^{-2}\bv_p \), so that
the relative change is \( \epsilon\kappa^2\eta \). Compare with the first-order bound
\( \epsilon\bigl(2\kappa+\kappa(\kappa+1)\eta\bigr) \).
:::
