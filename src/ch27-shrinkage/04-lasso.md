# The lasso

Ridge regression and principal components shrink, but neither selects (@exr-shr-pcr-not-selection).
Subset selection does, but it is a combinatorial search whose output
jumps when the data change slightly. Tibshirani (1996) replaced the squared penalty of ridge by
the sum of the absolute coefficients. The resulting **lasso** ("least absolute shrinkage and
selection operator") is convex, like ridge, and sets some coefficients exactly to zero, like
subset selection; Chen, Donoho and Saunders (1998) introduced the same criterion as basis pursuit
denoising. The theory for \( p>n \) is in [Chapter 28](../ch28-high-dimensional/index.html), and generalizations in [Chapter 30](../ch30-regularization-boosting/index.html).

## Definition and optimality conditions

The conventions of [Section 27.2](02-ridge.html) remain: \( \X \) is centred with columns on a common
scale, \( \y \) is centred, and the intercept is \( \bar y \) (@prp-shr-centring). Here \( \X \) may have
any rank.

::: {#def-shr-lasso}
[The lasso]

For \( \lambda\ge0 \), a **lasso estimate** \( \hbeta^{\mathrm{L}}_\lambda \) minimizes
\[
L_\lambda(\bb)=\tfrac12\norm{\y-\X\bb}^2+\lambda\norm{\bb}_1,\qquad \norm{\bb}_1=\sum_j\lvert b_j\rvert .
\]{#eq-shr-lasso}

This is the **penalized form**; the **constrained form** minimizes \( \norm{\y-\X\bb}^2 \) subject to
\( \norm{\bb}_1\le t \). A vector \( \mathbf{s} \) is a **subgradient** of \( \norm{\cdot}_1 \) at \( \bb \), written
\( \mathbf{s}\in\partial\norm{\bb}_1 \), if \( s_j=\operatorname{sign}(b_j) \) when \( b_j\ne0 \) and
\( s_j\in[-1,1] \) when \( b_j=0 \). The **Karush–Kuhn–Tucker (KKT) conditions** for @eq-shr-lasso are
\[
\X\T(\y-\X\bb)=\lambda\mathbf{s}\quad\text{for some }\mathbf{s}\in\partial\norm{\bb}_1 ;
\]{#eq-shr-kkt}

that is, \( \lvert\x_j\T(\y-\X\bb)\rvert\le\lambda \) for every \( j \), with
\( \x_j\T(\y-\X\bb)=\lambda\operatorname{sign}(b_j) \) whenever \( b_j\ne0 \).
:::

The factor \( \tfrac12 \) removes a \( 2 \) from the KKT conditions; [Chapter 28](../ch28-high-dimensional/index.html) writes the criterion as
\( (2n)^{-1}\norm{\y-\X\bb}^2+\lambda\norm{\bb}_1 \), whose \( \lambda \) is the present one divided by \( n \).
The KKT conditions say that no regressor's inner product with the residual exceeds \( \lambda \) in
absolute value, and that only regressors at the bound are active.

::: {#thm-shr-lasso-kkt}
[Lasso solutions]

Let \( \lambda>0 \).

::: {.enumerate options="label=(\alph*)"}
1. \( L_\lambda \) has at least one minimizer, and the set of minimizers is convex.

2. \( \bb \) minimizes \( L_\lambda \) iff it satisfies the KKT conditions @eq-shr-kkt.

3. All minimizers have the same fitted vector \( \X\bb \) and the same norm \( \norm{\bb}_1 \). If \( \X \) has
   full column rank, the minimizer is unique.

4. \( \bzero \) is a minimizer iff \( \lambda\ge\lambda_{\max}=\norm{\X\T\y}_\infty=\max_j\lvert\x_j\T\y\rvert \).

5. Every minimizer \( \hat{\bb} \) also minimizes \( \norm{\y-\X\bb}^2 \) subject to \( \norm{\bb}_1\le\norm{\hat{\bb}}_1 \).
:::
:::

::: {.proof}
(a) \( L_\lambda \) is continuous and \( L_\lambda(\bb)\ge\lambda\norm{\bb}_1 \), so the set
\( \{\bb:L_\lambda(\bb)\le L_\lambda(\bzero)\} \) is closed and bounded, and \( L_\lambda \) attains its minimum
there. Both terms of \( L_\lambda \) are convex, so if \( \bb,\bb' \) are minimizers, every point
between them has value at most the minimum, and is a minimizer too.

(b) *Sufficiency.* Let \( \bb \) satisfy @eq-shr-kkt with subgradient \( \mathbf{s} \), and let \( \mathbf{c} \) be
arbitrary. Expanding the square,
\[
\begin{aligned}
\tfrac12\norm{\y-\X\mathbf{c}}^2&=\tfrac12\norm{\y-\X\bb}^2-(\mathbf{c}-\bb)\T\X\T(\y-\X\bb)\\
&\qquad+\tfrac12\norm{\X(\mathbf{c}-\bb)}^2\\
&\ge\tfrac12\norm{\y-\X\bb}^2-\lambda\mathbf{s}\T(\mathbf{c}-\bb).
\end{aligned}
\]
Since \( \lvert s_j\rvert\le1 \) and \( s_jb_j=\lvert b_j\rvert \), we have \( \norm{\mathbf{c}}_1\ge\mathbf{s}\T\mathbf{c}=\norm{\bb}_1+\mathbf{s}\T(\mathbf{c}-\bb) \).
Adding \( \lambda \) times this to the previous inequality gives \( L_\lambda(\mathbf{c})\ge L_\lambda(\bb) \).
*Necessity.* Let \( \bb \) be a minimizer, \( \br=\y-\X\bb \) and \( g_j=\x_j\T\br \). For fixed \( j \) consider
\( \phi(t)=L_\lambda(\bb+t\mathbf{e}_j) \), which is minimized at \( t=0 \). Its one-sided derivatives at zero are
\[
\phi'(0^+)=-g_j+\lambda\,\sigma_+,\qquad \phi'(0^-)=-g_j+\lambda\,\sigma_-,
\]
with \( \sigma_\pm=\operatorname{sign}(b_j) \) if \( b_j\ne0 \), and \( \sigma_+=1 \), \( \sigma_-=-1 \) if \( b_j=0 \). A minimum at
zero requires \( \phi'(0^+)\ge0\ge\phi'(0^-) \). If \( b_j\ne0 \) this gives \( g_j=\lambda\operatorname{sign}(b_j) \); if
\( b_j=0 \) it gives \( -\lambda\le g_j\le\lambda \). These are the KKT conditions with \( s_j=g_j/\lambda \).

(c) Let \( \bb,\bb' \) be minimizers with \( \X\bb\ne\X\bb' \), and let \( \mathbf{m}=(\bb+\bb')/2 \). The function
\( \bu\mapsto\norm{\y-\bu}^2 \) is strictly convex, so
\( \norm{\y-\X\mathbf{m}}^2<\tfrac12\norm{\y-\X\bb}^2+\tfrac12\norm{\y-\X\bb'}^2 \), while
\( \norm{\mathbf{m}}_1\le\tfrac12\norm{\bb}_1+\tfrac12\norm{\bb'}_1 \). Then \( L_\lambda(\mathbf{m})<L_\lambda(\bb) \), a
contradiction. Equal fitted vectors give equal residual sums of squares, and equal values of
\( L_\lambda \) then force equal \( \ell_1 \) norms. If \( \X \) has full column rank, \( \X\bb=\X\bb' \) implies
\( \bb=\bb' \).

(d) At \( \bb=\bzero \) every \( s_j \) may be chosen in \( [-1,1] \), so @eq-shr-kkt holds iff
\( \lvert\x_j\T\y\rvert\le\lambda \) for all \( j \). Apply (b).

(e) This is the remark on penalized and constrained forms in [Section 27.2](02-ridge.html) with
\( P=\norm{\cdot}_1 \).
:::

By (b) and (c), all solutions share the inner products \( \x_j\T(\y-\X\bb) \) and vanish outside the
set of \( j \) at which these equal \( \pm\lambda \), a set determined by the data even when the solution is
not. The converse of (e) is Lagrangian duality for convex problems (Boyd and
Vandenberghe 2004, chapter 5), which we shall not need.

The constrained forms show why the lasso produces exact zeros and ridge does not. The constrained
estimate is where the smallest contour of \( \norm{\y-\X\bb}^2 \), an ellipsoid, touches the constraint
region. The ridge ball has a single normal direction at each boundary point, so the contact lies on
a coordinate subspace only by accident. The lasso's cross-polytope has vertices and edges in
coordinate subspaces, each with a whole cone of normal directions, so for an open set of data the
contact is at such a face and some coordinates are exactly zero. Algebraically, at \( b_j=0 \) the
subgradient may be anything in \( [-1,1] \).

## The orthonormal case: soft thresholding

Under an orthonormal design the lasso, ridge regression and subset selection are all explicit.
Define the **soft** and **hard
thresholding** functions
\[
S_\lambda(z)=\operatorname{sign}(z)\,(\lvert z\rvert-\lambda)_+ ,\qquad H_\lambda(z)=z\,\mathbf{1}\{\lvert z\rvert>\lambda\}.
\]

::: {#thm-shr-lasso-orthonormal}
[Three kinds of shrinkage under an orthonormal design]

Suppose \( \X\T\X=\I \), so that \( \hbeta=\X\T\y \), and let \( \lambda>0 \).

::: {.enumerate options="label=(\alph*)"}
1. *(Lasso.)* The unique minimizer of @eq-shr-lasso is given by \( \hat{\beta}^{\mathrm{L}}_{\lambda,j}=S_\lambda(\hat{\beta}_j) \).

2. *(Ridge.)* The minimizer of @eq-shr-ridge-criterion is \( \hbeta/(1+\lambda) \).

3. *(Subset selection.)* The minimizers of \( \tfrac12\norm{\y-\X\bb}^2+\tfrac12\lambda^2\#\{j:b_j\ne0\} \)
   are given by \( b_j=H_\lambda(\hat{\beta}_j) \) (with either choice when \( \lvert\hat{\beta}_j\rvert=\lambda \)), and for
   each \( k \) the best subset of size \( k \) consists of the \( k \) largest \( \lvert\hat{\beta}_j\rvert \).
:::
:::

::: {.proof}
Since \( \X\T(\y-\X\hbeta)=\bzero \),
\( \y-\X\bb=(\y-\X\hbeta)+\X(\hbeta-\bb) \) is an orthogonal sum, and \( \norm{\X\mathbf{w}}=\norm{\mathbf{w}} \), so
\[
\norm{\y-\X\bb}^2=\norm{\y-\X\hbeta}^2+\sum_j(\hat{\beta}_j-b_j)^2 .
\]
Each criterion therefore separates into one-dimensional problems. (a) Minimize
\( \tfrac12(\hat{\beta}_j-b)^2+\lambda\lvert b\rvert \) over \( b \). For \( b>0 \) the derivative is
\( b-\hat{\beta}_j+\lambda \), which vanishes at \( b=\hat{\beta}_j-\lambda \), admissible iff \( \hat{\beta}_j>\lambda \);
symmetrically for \( b<0 \); otherwise the convex function is minimized at \( b=0 \). This is
\( S_\lambda(\hat{\beta}_j) \), and it is unique because the problem is strictly convex. (b) Minimize
\( (\hat{\beta}_j-b)^2+\lambda b^2 \), a quadratic with minimum at \( \hat{\beta}_j/(1+\lambda) \). (c) The \( j \)th
problem costs \( \tfrac12\hat{\beta}_j^2 \) with \( b_j=0 \), and at best \( \tfrac12\lambda^2 \), at \( b_j=\hat{\beta}_j \),
with \( b_j\ne0 \). So \( b_j=\hat{\beta}_j \) iff \( \hat{\beta}_j^2>\lambda^2 \). For a fixed subset \( A \) of size \( k \) the
residual sum of squares is \( \norm{\y-\X\hbeta}^2+\sum_{j\notin A}\hat{\beta}_j^2 \), smallest when \( A \) holds
the \( k \) largest \( \lvert\hat{\beta}_j\rvert \).
:::

[Figure 27.4.1](#fig-shr-lasso)(a) draws the three rules. Ridge shrinks every estimate by the same
factor; hard thresholding keeps large estimates and kills small ones, with a jump at \( \pm\lambda \);
soft thresholding kills small estimates and moves the others towards zero by \( \lambda \), continuous
like ridge and sparse like subset selection, at the price of a bias \( \lambda \) on large coefficients. Soft thresholding of independent normal estimates was
analysed by Donoho and Johnstone (1994), whose threshold \( \sigma\sqrt{2\log p} \) reappears in [Chapter 28](../ch28-high-dimensional/index.html). The lasso is also a
posterior *mode* under independent Laplace priors (Tibshirani 1996; @exr-shr-laplace).

## The lasso path

As \( \lambda \) decreases from \( \lambda_{\max} \) the lasso moves from \( \bzero \) towards a least squares
fit, and the KKT conditions show that it does so along straight lines.

::: {#prp-shr-lasso-path}
[Piecewise linearity]

Suppose that for all \( \lambda \) in an interval \( [\lambda_1,\lambda_2] \) with \( \lambda_1>0 \) the lasso solution
is unique, has the same active set \( A=\{j:\hat{\beta}^{\mathrm{L}}_{\lambda,j}\ne0\} \) and the same signs
\( \mathbf{s}_A \), and that \( \X_A \), the columns of \( \X \) in \( A \), has full column rank. Then on that interval
\[
\hbeta^{\mathrm{L}}_{\lambda,A}=(\X_A\T\X_A)^{-1}(\X_A\T\y-\lambda\mathbf{s}_A),
\]
an affine function of \( \lambda \), and the other coordinates are zero.
:::

::: {.proof}
The KKT equations for \( j\in A \) read \( \X_A\T(\y-\X_A\hbeta^{\mathrm{L}}_{\lambda,A})=\lambda\mathbf{s}_A \), and \( \X_A\T\X_A \) is
invertible.
:::

So the path is made of line segments joined where a variable enters or leaves the active set, and
the whole path costs about as much as one least squares fit (Osborne, Presnell and Turlach 2000;
Efron, Hastie, Johnstone and Tibshirani 2004; [Section 30.5](../ch30-regularization-boosting/05-boosting-as-regularization.html) returns to this). A variable leaves when
its affine path reaches zero, which happens with correlated regressors but never under an
orthonormal design (@exr-shr-orthonormal-path).

## Computing the lasso by coordinate descent

As a function of one coordinate \( b_j \), with the others fixed, @eq-shr-lasso is a one-dimensional
lasso, solved by soft thresholding. Cycling through the coordinates gives a simple, fast algorithm, the "shooting" method of Fu (1998),
popularized by Friedman, Hastie, Höfling and Tibshirani (2007).

::: {#lem-shr-coordinate}
[Coordinate updates]

Fix \( \bb \) and \( j \) with \( \x_j\ne\bzero \), and let \( \br_{(j)}=\y-\sum_{k\ne j}\x_kb_k \) be the partial residual.

::: {.enumerate options="label=(\alph*)"}
1. The minimizer of \( L_\lambda \) over \( b_j \), with the other coordinates fixed, is
   \( S_\lambda(\x_j\T\br_{(j)})/\norm{\x_j}^2 \).

2. If \( \bb \) is unchanged by the update of every coordinate, then \( \bb \) satisfies the KKT conditions
   @eq-shr-kkt, and is therefore a lasso solution.
:::
:::

::: {.proof}
(a) With \( z=\x_j\T\br_{(j)} \) and \( q=\norm{\x_j}^2 \),
\( \tfrac12\norm{\br_{(j)}-\x_jb_j}^2=\tfrac12q\,(b_j-z/q)^2+\text{const} \). Dividing the one-dimensional
criterion by \( q \) gives \( \tfrac12(b_j-z/q)^2+(\lambda/q)\lvert b_j\rvert \), minimized by
\( S_{\lambda/q}(z/q)=S_\lambda(z)/q \) as in @thm-shr-lasso-orthonormal(a).
(b) A fixed point has \( b_j=S_\lambda(z_j)/q_j \) for every \( j \), which by the one-dimensional case of
@thm-shr-lasso-kkt(b) means \( z_j-q_jb_j\in\lambda\,\partial\lvert b_j\rvert \). But
\( z_j-q_jb_j=\x_j\T(\y-\X\bb) \). So @eq-shr-kkt holds, and @thm-shr-lasso-kkt(b) applies.
:::

::: {.algorithm}
**Cyclic coordinate descent for the lasso.** Start from \( \bb=\bzero \) (or a warm start) and
\( \br=\y-\X\bb \). Repeat until no coordinate changes by more than a tolerance: for
\( j=1,\dots,p \), compute \( z=\x_j\T\br+\norm{\x_j}^2b_j \), set
\( b_j^{\text{new}}=S_\lambda(z)/\norm{\x_j}^2 \), update \( \br\leftarrow\br-\x_j(b_j^{\text{new}}-b_j) \) and
\( b_j\leftarrow b_j^{\text{new}} \).
:::

Each update lowers \( L_\lambda \) or leaves it unchanged. Tseng (2001) proved, for criteria that are
smooth and convex plus a sum of convex functions of single coordinates, that every limit point of
cyclic coordinate descent is a minimizer; we do not reproduce the proof. The lemma gives what a
computation needs: at a fixed point the KKT conditions hold, and their violation can be measured.
In practice the lasso is computed on a decreasing grid of \( \lambda \) from \( \lambda_{\max} \), each
solution the warm start for the next (Friedman, Hastie and Tibshirani 2010). Highly correlated
regressors slow convergence, because single-coordinate moves zigzag along a narrow valley.

::: {#exm-shr-longley-lasso}
[The lasso for Longley's data]

For the standardized Longley data, \( \lambda_{\max}=51813 \), and the variables enter in the order
GNP, UNEMP, ARMED, YEAR, POP; the GNP deflator has not entered by \( \lambda_{\max}/1000 \).
[Figure 27.4.1](#fig-shr-lasso)(b) shows the path against the ratio of \( \norm{\hbeta^{\mathrm{L}}_\lambda}_1 \)
to the \( \ell_1 \) norm of least squares. GNP, the regressor most correlated with employment, grows
to over \( 3000 \); when YEAR enters, the two compete for the same direction, and GNP falls back to
exactly zero and leaves the active set. Coordinate descent needs \( 42 \) sweeps at
\( \lambda_{\max}/2 \), \( 157 \) at \( \lambda_{\max}/50 \) and \( 3479 \) at \( \lambda_{\max}/1000 \) for a tolerance of
\( 10^{-10} \): the zigzag effect of collinearity.

Leave-one-out cross-validation, refitting the path without each year in turn, chooses
\( \lambda=130 \), about \( 0.0025\lambda_{\max} \), with \( 3 \) nonzero coefficients: \( -1303 \) for
unemployment, \( -484 \) for the armed forces and \( 4475 \) for the year. Its minimized root mean
squared error is \( 373 \), but as with ridge the minimum is optimistic; nested leave-one-out gives
\( 428 \), level with least squares (\( 425 \)). What the lasso offers here is a simpler fit of equal
accuracy, not a better one, and its choice of YEAR over GNP is exactly the kind that collinearity
makes fragile.

The script checks the solutions at five values of \( \lambda \) three ways: the KKT conditions hold to
within \( 10^{-6}\lambda \); a brute-force search over all \( 3^6=729 \) active sets and signs, solving
the equations of @prp-shr-lasso-path, finds the same solution; and a bound-constrained optimizer
on the split \( \bb=\bu-\bv \), \( \bu,\bv\ge\bzero \), finds no smaller criterion value.
:::

::: {when-format="html"}
![**Figure 27.4.1.** (a) Lasso (soft thresholding), subset selection (hard thresholding) and
ridge regression (proportional shrinkage) as functions of the least squares coefficient under an
orthonormal design, with \( \lambda=1 \) (@thm-shr-lasso-orthonormal). (b) The lasso path for the
standardized Longley data; the vertical line marks the leave-one-out choice of \( \lambda \).](lasso_path.svg){#fig-shr-lasso width=100%}
:::

::: {when-format="pdf"}
![(a) Lasso (soft thresholding), subset selection (hard thresholding) and
ridge regression (proportional shrinkage) as functions of the least squares coefficient under an
orthonormal design, with \( \lambda=1 \) (@thm-shr-lasso-orthonormal). (b) The lasso path for the
standardized Longley data; the vertical line marks the leave-one-out choice of \( \lambda \).](lasso_path.pdf){width=100%}
:::

```{.python .run #cell-lasso-cd-cd}
import numpy as np
import statsmodels.api as sm
def soft(z, t):
    """Soft thresholding: sign(z) * max(|z| - t, 0)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)

def lasso_cd(X, y, lam, b=None, tol=1e-10, max_sweeps=100_000):
    """Cyclic coordinate descent for (1/2)||y - X b||^2 + lam ||b||_1."""
    p = X.shape[1]
    b = np.zeros(p) if b is None else b.copy()
    r = y - X @ b                                   # current residual
    sq = np.sum(X ** 2, axis=0)                     # ||x_j||^2
    for sweep in range(1, max_sweeps + 1):
        biggest = 0.0
        for j in range(p):
            z = X[:, j] @ r + sq[j] * b[j]          # x_j^T (partial residual without j)
            new = soft(z, lam) / sq[j]
            if new != b[j]:
                r -= X[:, j] * (new - b[j])
                biggest = max(biggest, abs(new - b[j]) * np.sqrt(sq[j]))
                b[j] = new
        if biggest < tol:
            break
    return b, sweep

def kkt_violation(X, y, b, lam):
    """Largest violation of the KKT conditions X^T (y - X b) = lam * s, s in the subdifferential."""
    g = X.T @ (y - X @ b)
    active = b != 0
    return max(np.max(np.abs(g[active] - lam * np.sign(b[active])), initial=0.0),
               np.max(np.abs(g[~active]) - lam, initial=0.0))
```

The second cell computes the path on a grid of \( 61 \) values of \( \lambda \) with warm starts.

```{.python .run #cell-lasso-cd-longley}
data = sm.datasets.longley.load_pandas().data
names = ["GNPDEFL", "GNP", "UNEMP", "ARMED", "POP", "YEAR"]
Z = data[names].to_numpy()
y = data["TOTEMP"].to_numpy()
n, p = Z.shape
X = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
yc = y - y.mean()

lam_max = np.max(np.abs(X.T @ yc))                  # smallest lam with b = 0
lams = lam_max * np.logspace(0, -3, 61)
path, b = [], np.zeros(p)
for lam in lams:                                     # warm starts along a decreasing grid
    b, _ = lasso_cd(X, yc, lam, b, tol=1e-6)
    path.append(b.copy())
path = np.array(path)
first = {names[j]: np.argmax(path[:, j] != 0) for j in range(p) if np.any(path[:, j] != 0)}
entry = sorted(first, key=first.get)                  # variables in order of entry
print(f"lambda_max = {lam_max:.0f}; order of entry:", entry)
print("GNP coefficient along the path:", np.round(path[::6, 1], 0))
```

## Lasso, ridge and subset selection compared

@thm-shr-lasso-orthonormal suggests that hard thresholding suits a few large coefficients among
zeros, proportional shrinkage many moderate ones, and soft thresholding lies between. A simulation
with correlated regressors tests this.

::: {#exm-shr-lasso-simulation}
[Three coefficient patterns]

Fix a design with \( n=50 \) and \( p=10 \), its rows drawn once from a normal distribution with
correlations \( 0.5^{\lvert j-k\rvert} \), and columns scaled to \( \norm{\x_j}^2=n \); the errors are
normal with \( \sigma=2 \). The coefficients are *sparse*, \( (3,-2,1.5,0,\dots,0) \); *dense*, ten
coefficients \( \pm0.5 \) with alternating signs; or *intermediate*, five coefficients \( \pm1 \) and five
zeros. Each method's tuning constant (\( \lambda \), or the subset size, with the best subsets found
among all \( 1024 \)) minimizes the error on an independent copy of the responses, an idealized
validation set. The table gives the average of
\( \norm{\X\tilde{\bbeta}-\X\bbeta}^2/(p\sigma^2) \) over \( 300 \) data sets; least squares has expectation
one (@exr-shr-smoother-mse), and the Monte Carlo standard error of each entry is about \( 0.02 \).

| Pattern | least squares | ridge | lasso | best subset |
|---|---|---|---|---|
| sparse | 1.02 | 0.95 | 0.74 | 0.46 |
| dense | 0.99 | 0.73 | 0.83 | 1.03 |
| intermediate | 0.99 | 0.84 | 0.86 | 0.90 |

Best subset selection is far ahead in the sparse case and worst in the dense case, where it cannot
shrink and finds no zeros; ridge is the reverse; the lasso is never best and never worst. In the
intermediate case the paired difference between the lasso and ridge is \( 0.02 \) with a standard
error of \( 0.01 \). Hastie, Tibshirani and Tibshirani (2020) reached similar conclusions in a much
larger study, adding that the ranking depends strongly on the signal-to-noise ratio
([Figure 27.4.2](#fig-shr-lasso-sim) displays the table).
:::

::: {when-format="html"}
![**Figure 27.4.2.** Prediction risk at the design points, relative to \( p\sigma^2 \), for least squares,
ridge, the lasso and best subset selection, each tuned on a validation copy of the responses,
under three patterns of true coefficients (\( 300 \) data sets each).](lasso_simulation.svg){#fig-shr-lasso-sim width=80%}
:::

::: {when-format="pdf"}
![Prediction risk at the design points, relative to \( p\sigma^2 \), for least squares,
ridge, the lasso and best subset selection, each tuned on a validation copy of the responses,
under three patterns of true coefficients (\( 300 \) data sets each).](lasso_simulation.pdf){width=80%}
:::

The cell below reruns the sparse scenario with \( 10 \) data sets; each one searches all \( 1024 \)
subsets, so in the browser it may take half a minute. The script uses \( 300 \) per scenario.

```{.python .run #cell-lasso-simulation-setup}
import itertools
import numpy as np
def soft(z, t):
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)

def lasso_cd(X, y, lam, b=None, tol=1e-8, max_sweeps=10_000):
    """Cyclic coordinate descent for (1/2)||y - X b||^2 + lam ||b||_1."""
    p = X.shape[1]
    b = np.zeros(p) if b is None else b.copy()
    r = y - X @ b
    sq = np.sum(X ** 2, axis=0)
    for _ in range(max_sweeps):
        biggest = 0.0
        for j in range(p):
            new = soft(X[:, j] @ r + sq[j] * b[j], lam) / sq[j]
            if new != b[j]:
                r -= X[:, j] * (new - b[j])
                biggest = max(biggest, abs(new - b[j]))
                b[j] = new
        if biggest < tol:
            break
    return b

rng = np.random.default_rng(2027)
n, p, sigma = 50, 10, 2.0
Sigma = 0.5 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
X = rng.multivariate_normal(np.zeros(p), Sigma, size=n)
X = (X - X.mean(axis=0)) / np.sqrt(np.mean((X - X.mean(axis=0)) ** 2, axis=0))   # ||x_j||^2 = n
U, d, Vt = np.linalg.svd(X, full_matrices=False)
subsets = [list(S) for k in range(p + 1) for S in itertools.combinations(range(p), k)]
ridge_grid = np.logspace(-2, 4, 40)
scenarios = {
    "sparse": np.r_[3.0, -2.0, 1.5, np.zeros(p - 3)],
    "dense": 0.5 * (-1.0) ** np.arange(p),
    "intermediate": np.r_[1.0, -1.0, 1.0, -1.0, 1.0, np.zeros(p - 5)],
}

def fits(y):
    """Candidate fitted vectors of each method along its tuning path."""
    c = U.T @ y
    out = {"least squares": [U @ c]}
    out["ridge"] = [U @ (d ** 2 / (d ** 2 + lam) * c) for lam in ridge_grid]
    lam_max = np.max(np.abs(X.T @ y))
    b, lasso = np.zeros(p), []
    for lam in lam_max * np.logspace(0, -3, 30):
        b = lasso_cd(X, y, lam, b)
        lasso.append(X @ b)
    out["lasso"] = lasso
    best = {}                                    # best subset of each size, by residual sum of squares
    for S in subsets:
        fit = X[:, S] @ np.linalg.lstsq(X[:, S], y, rcond=None)[0] if S else np.zeros(n)
        rss = np.sum((y - fit) ** 2)
        if len(S) not in best or rss < best[len(S)][0]:
            best[len(S)] = (rss, fit)
    out["best subset"] = [best[k][1] for k in range(p + 1)]
    return out

def run(beta, reps, seed=1):
    """Mean of ||X b - X beta||^2 / (p sigma^2) for each method, tuned on a validation copy."""
    gen = np.random.default_rng(seed)
    mu = X @ beta
    loss = {m: [] for m in ["least squares", "ridge", "lasso", "best subset"]}
    for _ in range(reps):
        y = mu + sigma * gen.normal(size=n)
        y_val = mu + sigma * gen.normal(size=n)
        for m, cands in fits(y).items():
            chosen = min(cands, key=lambda f: np.sum((y_val - f) ** 2))
            loss[m].append(np.sum((chosen - mu) ** 2) / (p * sigma ** 2))
    return {m: np.array(v) for m, v in loss.items()}

quick = run(scenarios["sparse"], reps=10)
print({m: round(float(v.mean()), 2) for m, v in quick.items()})
```

## Exercises

### A. Check your understanding

::: {#exr-shr-lambda-max}
[A1]

Suppose \( \max_j\lvert\x_j\T\y\rvert \) is attained at a single \( j^* \). Show that for \( \lambda \) slightly below
\( \lambda_{\max} \) the lasso solution has one nonzero coordinate, \( j^* \), with the sign of \( \x_{j^*}\T\y \)
and absolute value \( (\lvert\x_{j^*}\T\y\rvert-\lambda)/\norm{\x_{j^*}}^2 \).
:::

::: {.solution}
Try \( \bb=b\mathbf{e}_{j^*} \) with \( b=\operatorname{sign}(g)(\lvert g\rvert-\lambda)/q \), \( g=\x_{j^*}\T\y \),
\( q=\norm{\x_{j^*}}^2 \). Then \( \x_{j^*}\T(\y-\X\bb)=g-qb=\lambda\operatorname{sign}(g) \), as required. For
\( k\ne j^* \), \( \lvert\x_k\T(\y-\X\bb)\rvert\le\lvert\x_k\T\y\rvert+\lvert b\rvert\lvert\x_k\T\x_{j^*}\rvert \). The first term is
strictly less than \( \lambda_{\max} \), and \( b\to0 \) as \( \lambda\to\lambda_{\max} \), so the bound is below
\( \lambda \) for \( \lambda \) close enough to \( \lambda_{\max} \). By @thm-shr-lasso-kkt(b) \( \bb \) is a solution; it is
unique among solutions supported on \( \{j^*\} \), and every solution is supported on the set where the
inner product with the residual equals \( \pm\lambda \), which is \( \{j^*\} \).
:::

::: {#exr-shr-threshold-numbers}
[A2]

With an orthonormal design and \( \hbeta=(3,-0.5,1.2)\T \), compute the lasso, ridge and
hard-thresholding estimates of @thm-shr-lasso-orthonormal with \( \lambda=1 \).
:::

### B. Practice

::: {#exr-shr-duplicate}
[B1]

Suppose \( \x_1=\x_2 \) and let \( \hat{\bb} \) be a lasso solution with \( \hat b_1+\hat b_2=\gamma\ne0 \). Show that
every \( \bb \) with \( b_1+b_2=\gamma \), \( b_1,b_2 \) of the sign of \( \gamma \), and the other coordinates as in
\( \hat{\bb} \), is also a solution, whereas ridge regression gives \( b_1=b_2 \).
:::

::: {.solution}
Such a \( \bb \) has \( \X\bb=\X\hat{\bb} \) and \( \lvert b_1\rvert+\lvert b_2\rvert=\lvert\gamma\rvert\le\lvert\hat b_1\rvert+\lvert\hat b_2\rvert \), so
\( L_\lambda(\bb)\le L_\lambda(\hat{\bb}) \), with equality by optimality. For ridge, replacing \( (b_1,b_2) \) by
their average leaves \( \X\bb \) unchanged and does not increase \( b_1^2+b_2^2 \), strictly decreasing it
unless \( b_1=b_2 \); since the ridge solution is unique, it has \( b_1=b_2 \). The elastic net of [Chapter 30](../ch30-regularization-boosting/index.html)
combines the two penalties to get the ridge behaviour with lasso sparsity.
:::

::: {#exr-shr-laplace}
[B2]

Let \( \Y\mid\bbeta\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \sigma^2 \) known, and let the \( \beta_j \) be independent
with density \( (\lambda/2\sigma^2)\exp(-\lambda\lvert\beta_j\rvert/\sigma^2) \). Show that the posterior mode is the
lasso estimate. Is the posterior mean sparse?
:::

::: {.solution}
The log posterior is \( -\{\tfrac12\norm{\y-\X\bbeta}^2+\lambda\norm{\bbeta}_1\}/\sigma^2+\text{const} \), maximized
at the minimizer of @eq-shr-lasso. The posterior mean is an average over a continuous
distribution with positive density everywhere, so every coordinate is nonzero with probability
one: the sparsity of the lasso is a property of the mode, not of the Bayesian model.
:::

::: {#exr-shr-orthonormal-path}
[B3]

Under an orthonormal design, describe the lasso path: show that it is piecewise linear with
kinks at \( \lambda=\lvert\hat{\beta}_j\rvert \), that variables enter in decreasing order of \( \lvert\hat{\beta}_j\rvert \),
and that none ever leaves.
:::

### C. Going deeper

::: {#exr-shr-lasso-df}
[C1]

In the orthonormal case with normal errors, show that
\( \sigma^{-2}\sum_j\Cov(S_\lambda(\hat{\beta}_j),\hat{\beta}_j)=\E\,\#\{j:\lvert\hat{\beta}_j\rvert>\lambda\} \): the expected number of
nonzero coefficients plays the role of degrees of freedom (compare @exr-shr-df-covariance). Hint:
integrate by parts, one coordinate at a time.
:::

::: {#exr-shr-sparsity-bound}
[C2]

Show that for every \( \lambda>0 \) there is a lasso solution with at most \( \min(n,p) \) nonzero coordinates.
Hint: if the active columns of a solution are linearly dependent, move along a direction in their
null space that keeps the signs and the fitted vector until a coordinate hits zero.
:::
