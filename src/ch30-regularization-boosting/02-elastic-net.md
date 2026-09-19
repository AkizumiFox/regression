# The elastic net

When several regressors are near copies of one another (repeated measurements, genes in one pathway), the lasso
tends to keep one and drop the rest, and which one can change with a small perturbation of the data. With two
identical columns its solution is not even unique (Exercise B1 of [Section 27.4](../ch27-shrinkage/04-lasso.html)). Ridge spreads the coefficient evenly over
the copies but never sets anything to zero. Zou and Hastie (2005) proposed using both penalties.

## Definition and basic properties

::: {#def-reg-elastic-net}
[Elastic net]

For \( \lambda_1,\lambda_2\ge0 \), the **elastic net** criterion is
\[
Q(\bb)=\tfrac12\norm{\y-\X\bb}^2+\lambda_1\norm{\bb}_1+\tfrac12\lambda_2\norm{\bb}^2,
\]{#eq-reg-enet}

and its minimizer \( \hbeta^{\mathrm{en}} \) is the (naive) **elastic net estimate**. An equivalent
parameterization is \( \lambda_1=\alpha\lambda \), \( \lambda_2=(1-\alpha)\lambda \), with a mixing
parameter \( \alpha\in[0,1] \) and an overall level \( \lambda\ge0 \).
:::

With \( \lambda_2=0 \) this is the lasso, and with \( \lambda_1=0 \) ridge. Here the columns of \( \X \) are centred with unit length and \( \y \) is centred, so \( \rho_{ij}=\x_i\T\x_j \) is a
sample correlation.

::: {#prp-reg-elastic-net}
[Elastic net: uniqueness, augmented data, grouping]

Let \( \lambda_2>0 \).

::: {.enumerate options="label=(\alph*)"}
1. The criterion @eq-reg-enet has a unique minimizer, whatever the rank of \( \X \). It is characterized by
   \[
   \x_j\T(\y-\X\hbeta)-\lambda_2\hat\beta_j=\lambda_1\operatorname{sign}(\hat\beta_j)\ \ (\hat\beta_j\neq0),
   \qquad \lvert\x_j\T(\y-\X\hbeta)\rvert\le\lambda_1\ \ (\hat\beta_j=0).
   \]{#eq-reg-enet-kkt}

2. *(Augmented data.)* Let \( \X^*=\begin{pmatrix}\X\\\sqrt{\lambda_2}\,\I_p\end{pmatrix} \) and
   \( \y^*=\begin{pmatrix}\y\\\bzero\end{pmatrix} \). Then \( \hbeta^{\mathrm{en}} \) is the lasso estimate, with parameter
   \( \lambda_1 \), for the data \( (\X^*,\y^*) \). The matrix \( \X^* \) has full column rank \( p \).

3. *(Orthonormal design.)* If \( \X\T\X=\I \), then \( \hat\beta_j^{\mathrm{en}}=\operatorname{sign}(z_j)(\lvert z_j\rvert-\lambda_1)_+/(1+\lambda_2) \),
   with \( \mathbf{z}=\X\T\y \).

4. *(Grouping effect.)* If \( \hat\beta_i\hat\beta_j>0 \), then
   \[
   \lvert\hat\beta_i-\hat\beta_j\rvert\le\frac{\norm{\y-\X\hbeta}\,\norm{\x_i-\x_j}}{\lambda_2}\le\frac{\norm{\y}\sqrt{2(1-\rho_{ij})}}{\lambda_2}.
   \]{#eq-reg-grouping}

5. *(Identical columns.)* If \( \x_i=\x_j \), then \( \hat\beta_i=\hat\beta_j \).
:::

:::

::: {.proof}
(a) The penalty \( \lambda_1\norm{\bb}_1+\tfrac12\lambda_2\norm{\bb}^2 \) is strictly convex, because \( \norm{\bb}^2 \) is
strictly convex and \( \norm{\bb}_1 \) is convex. It is at least \( \tfrac12\lambda_2\norm{\bb}^2 \), so \( Q\to\infty \). Existence and
uniqueness follow as in @thm-reg-existence(a) and (c). By @lem-reg-optimality applied to the whole penalty, with
the gradient \( \lambda_2\hbeta \) of the smooth part moved to the left-hand side, and @lem-reg-norm-subgradient, the
minimizer is characterized by @eq-reg-enet-kkt.

(b) \( \norm{\y^*-\X^*\bb}^2=\norm{\y-\X\bb}^2+\lambda_2\norm{\bb}^2 \). So \( \tfrac12\norm{\y^*-\X^*\bb}^2+\lambda_1\norm{\bb}_1 \) is @eq-reg-enet. The last \( p \) rows of \( \X^* \) are \( \sqrt{\lambda_2}\,\I_p \), which has rank \( p \).

(c) With \( \X\T\X=\I \), @eq-reg-enet equals \( \sum_j\{\tfrac12(1+\lambda_2)b_j^2-z_jb_j+\lambda_1\lvert b_j\rvert\} \) plus a
constant. Dividing the \( j \)th term by \( 1+\lambda_2 \) gives
\( \tfrac12\{b_j-z_j/(1+\lambda_2)\}^2+\{\lambda_1/(1+\lambda_2)\}\lvert b_j\rvert \) plus a constant. By the univariate lasso solution (@thm-shr-lasso-orthonormal), this is minimized by soft thresholding \( z_j/(1+\lambda_2) \) at
\( \lambda_1/(1+\lambda_2) \), which is the stated formula.

(d) Both coefficients are nonzero with a common sign \( s \), so @eq-reg-enet-kkt gives
\( \x_i\T\br-\lambda_2\hat\beta_i=\lambda_1s=\x_j\T\br-\lambda_2\hat\beta_j \), where \( \br=\y-\X\hbeta \). Subtracting,
\( \lambda_2(\hat\beta_i-\hat\beta_j)=(\x_i-\x_j)\T\br \), and Cauchy–Schwarz gives the first inequality. For the second,
\( \norm{\x_i-\x_j}^2=2-2\rho_{ij} \) for unit-length columns, and \( \tfrac12\norm{\br}^2\le Q(\hbeta)\le Q(\bzero)=\tfrac12\norm{\y}^2 \).

(e) Swapping coordinates \( i \) and \( j \) of \( \hbeta \) changes neither \( \X\hbeta \) nor the penalty. So the swapped
vector is also a minimizer, and by uniqueness it equals \( \hbeta \).
:::

By (b), the elastic net is a lasso on \( n+p \) observations with a full-rank design, so it can select all \( p \)
regressors even when \( p>n \); a unique lasso solution has at most \( \min(n,p) \) nonzero coefficients (@exr-reg-lasso-at-most-n). Part (d) is the grouping effect of Zou and Hastie (2005): as \( \rho_{ij}\to1 \), coefficients of
the same sign are forced together at a rate set by \( \lambda_2 \). No such bound holds for the lasso.

::: {#exm-reg-grouping}
[Three near-copies of one signal]

The script simulates \( n=50 \) observations of \( p=10 \) regressors. The first three are noisy copies
of a common variable, with pairwise correlations between \( 0.93 \) and \( 0.95 \), and each
has coefficient 2 (on the unit-length scale). The other seven are independent noise. With \( \lambda_1=0.8 \):

| | \( \hat\beta_1 \) | \( \hat\beta_2 \) | \( \hat\beta_3 \) | nonzero |
|---|---|---|---|---|
| lasso (\( \lambda_2=0 \)) | 2.335 | 0.000 | 2.522 | 2 |
| elastic net (\( \lambda_2=1 \)) | 1.248 | 1.103 | 1.236 | 4 |
| rescaled elastic net | 2.496 | 2.207 | 2.471 | 4 |

The lasso drops the second copy. The elastic net keeps all three with similar coefficients (and one noise regressor, at
\( 0.081 \)). For the pair \( (1,2) \), with correlation \( 0.927 \), the spread \( 0.144 \) is well inside
@eq-reg-grouping, whose middle term is \( 1.577 \) (residual norm \( 4.129 \)) and right-hand side \( 2.532 \)
(\( \norm{\y}=6.630 \)). Along the lasso path ([Figure 30.2.1](#fig-reg-en-paths)) the second copy enters only after
\( \lambda \) has fallen by a factor of about \( 26 \); along the elastic net path all three enter at the same grid value.
:::

::: {when-format="html"}
![**Figure 30.2.1.** Coefficient paths for the data of @exm-reg-grouping, plotted against
\( \log_{10}\lambda \), with \( \lambda \) decreasing to the right. The three correlated regressors are coloured and
the seven noise regressors are grey. (a) The lasso. (b) The elastic net with \( \alpha=\tfrac12 \), that is,
\( \lambda_1=\lambda_2=\lambda/2 \).](en_paths.svg){#fig-reg-en-paths width=100%}
:::

::: {when-format="pdf"}
![Coefficient paths for the data of @exm-reg-grouping, plotted against
\( \log_{10}\lambda \), with \( \lambda \) decreasing to the right. The three correlated regressors are coloured and
the seven noise regressors are grey. (a) The lasso. (b) The elastic net with \( \alpha=\tfrac12 \), that is,
\( \lambda_1=\lambda_2=\lambda/2 \).](en_paths.pdf){width=100%}
:::

## Computing it: coordinate descent

**Cyclic coordinate descent**, used for the lasso in @lem-shr-coordinate, minimizes over one coordinate at a time, the others fixed, and sweeps until nothing
changes.

::: {#prp-reg-coordinate}
[Coordinate update]

Fix all coordinates of \( \bb \) except \( b_j \), and let \( \br_{(j)}=\y-\sum_{k\neq j}\x_kb_k \) be the partial residual.
The value of \( b_j \) that minimizes @eq-reg-enet is
\[
b_j=\frac{\operatorname{sign}(\x_j\T\br_{(j)})\bigl(\lvert\x_j\T\br_{(j)}\rvert-\lambda_1\bigr)_+}{\norm{\x_j}^2+\lambda_2}.
\]{#eq-reg-cd}

:::

::: {.proof}
As a function of \( b_j \) alone, @eq-reg-enet equals
\( \tfrac12(\norm{\x_j}^2+\lambda_2)b_j^2-b_j\x_j\T\br_{(j)}+\lambda_1\lvert b_j\rvert \) plus a constant. Divided by
\( \norm{\x_j}^2+\lambda_2 \), this is a univariate lasso problem of the form in the proof of @prp-reg-elastic-net(c).
Its solution is @eq-reg-cd.
:::

Because the nonsmooth part is separable, every limit point of the iterates is a minimizer (Tseng, 2001; cited, not proved),
so the iterates converge when the minimizer is unique, as for the elastic net with \( \lambda_2>0 \). For a coupling penalty such as \( \sum_j\lvert b_j-b_{j-1}\rvert \), coordinate descent can stall. Whole paths are
computed on a decreasing grid of \( \lambda \), each fit starting from the last (Friedman, Hastie, Höfling and Tibshirani, 2007;
Friedman, Hastie and Tibshirani, 2010). The listing works through the Gram matrix, as the next proposition requires.

```{.python .run #cell-elastic-net-compare}
import numpy as np

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

lam1 = 0.8
b_lasso = elastic_net(X, y, lam1, 0.0)
b_naive = elastic_net(X, y, lam1, 1.0)
print("lasso         ", np.round(b_lasso[:5], 3))
print("naive EN      ", np.round(b_naive[:5], 3))
print("rescaled EN   ", np.round(2.0 * b_naive[:5], 3))     # (1 + lam2) * naive
```

## Double shrinkage and the rescaled elastic net

By @prp-reg-elastic-net(c), under an orthonormal design the elastic net soft-thresholds and then divides by
\( 1+\lambda_2 \): surviving coefficients are shrunk twice. Zou and Hastie (2005) proposed reporting
\( (1+\lambda_2)\hbeta^{\mathrm{en}} \) instead, which in general is a lasso with a modified Gram matrix.

::: {#prp-reg-rescaled}
[The rescaled elastic net]

Let the columns of \( \X \) have unit length and let \( \lambda_2>0 \). Then \( \tilde{\bbeta}=(1+\lambda_2)\hbeta^{\mathrm{en}} \) is the
unique minimizer of
\[
\tfrac12\,\mathbf{c}\T\bSigma_{\lambda_2}\mathbf{c}-\y\T\X\mathbf{c}+\lambda_1\norm{\mathbf{c}}_1,
\qquad \bSigma_{\lambda_2}=\frac{\X\T\X+\lambda_2\I}{1+\lambda_2}.
\]{#eq-reg-rescaled}

:::

::: {.proof}
Up to the constant \( \tfrac12\y\T\y \), @eq-reg-enet is
\( \tfrac12\bb\T(\X\T\X+\lambda_2\I)\bb-\y\T\X\bb+\lambda_1\norm{\bb}_1 \). Substitute \( \bb=\mathbf{c}/(1+\lambda_2) \) and multiply by
\( 1+\lambda_2>0 \). The result is @eq-reg-rescaled. The substitution is a bijection, so minimizers correspond, and the
minimizer of @eq-reg-enet is unique by @prp-reg-elastic-net(a).
:::

The lasso has \( \X\T\X \) in place of \( \bSigma_{\lambda_2} \), whose correlations are divided by \( 1+\lambda_2 \): the rescaled elastic
net is a lasso computed as if the regressors were less correlated. The rescaling is a heuristic; in @exm-reg-grouping, where the
true coefficients are 2, the naive estimates are about half that and the rescaled ones somewhat too large. The glmnet software
(Friedman, Hastie and Tibshirani, 2010) fits the naive criterion.

## Exercises

### A. Check your understanding

::: {#exr-reg-en-orthonormal}
[A1]

Under an orthonormal design with \( \mathbf{z}=\X\T\y=(3,0.5,-1.2)\T \), \( \lambda_1=1 \) and \( \lambda_2=1 \), compute the lasso,
ridge, naive elastic net and rescaled elastic net estimates.
:::

::: {.solution}
Lasso (\( \lambda=1 \)): \( (2,0,-0.2) \). Ridge (\( \lambda=1 \)): \( \mathbf{z}/2=(1.5,0.25,-0.6) \). Naive elastic net: the lasso
divided by \( 2 \), \( (1,0,-0.1) \). Rescaled: multiply by \( 1+\lambda_2=2 \), which gives the lasso, \( (2,0,-0.2) \).
:::

### B. Practice

::: {#exr-reg-grouping-opposite}
[B1]

Show that if \( \hat\beta_i\hat\beta_j<0 \), then \( \lvert\hat\beta_i+\hat\beta_j\rvert\le\norm{\y}\sqrt{2(1+\rho_{ij})}/\lambda_2 \). Interpret this for
two regressors with correlation close to \( -1 \).
:::

::: {.solution}
With signs \( s \) and \( -s \), @eq-reg-enet-kkt gives \( \x_i\T\br-\lambda_2\hat\beta_i=\lambda_1s \) and
\( \x_j\T\br-\lambda_2\hat\beta_j=-\lambda_1s \). Adding, \( \lambda_2(\hat\beta_i+\hat\beta_j)=(\x_i+\x_j)\T\br \), and
\( \norm{\x_i+\x_j}^2=2+2\rho_{ij} \). Bound \( \norm{\br}\le\norm{\y} \) as before. If \( \rho_{ij}\approx-1 \), then \( \x_j\approx-\x_i \), and
the bound forces \( \hat\beta_j\approx-\hat\beta_i \). The two regressors are treated as one regressor with a sign flip.
:::

::: {#exr-reg-en-identical}
[B2]

Let \( \X=[\x,\x] \) with \( \norm{\x}=1 \). Find the elastic net estimate in closed form and compare it with the ridge and
lasso estimates.
:::

### C. Going deeper

::: {#exr-reg-lasso-at-most-n}
[C1]

Suppose \( p>n \) and the lasso estimate \( \hbeta \) is unique with active set \( A \). Use the KKT conditions to show
that the columns \( \X_A \) are linearly independent, and deduce that \( \lvert A\rvert\le n \). (Hint: if
\( \X_A\mathbf{d}=\bzero \) with \( \mathbf{d}\neq\bzero \), move along \( \hbeta+\tau(\mathbf{d},\bzero) \) and use that the fit and the \( \ell_1 \) norm are
constant among minimizers.)
:::
