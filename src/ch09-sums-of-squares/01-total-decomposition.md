# Decomposing the total sum of squares

Every sum of squares in a linear model is the squared length of a projection of the data.
[Chapter 6](../ch06-projections/index.html) showed this for the residual sum of squares,
\( \norm{(\I-\M)\y}^2 \), and for the drop in residual sum of squares between nested models,
\( \norm{(\M-\Mo)\y}^2 \) (@eq-proj-extra-ss). This section turns these observations into a
general framework. An analysis of variance is a choice of mutually orthogonal subspaces that
together fill the model space. Once that choice is made, the whole table follows from
Pythagoras: the sums of squares add up because the pieces are orthogonal, and the degrees of
freedom add up because dimensions of orthogonal subspaces add.

Throughout the chapter the model is \( \E(\Y)=\X\bbeta \) with \( \Cov(\Y)=\sigma^2\I \)
(@def-lm-linear-model). \( \X \) is \( n\times p \) of rank \( r \), which may be less than \( p \),
and \( \M \) is the orthogonal projection onto \( \C(\X) \). The identities of the chapter hold for
every data vector \( \y \). Probability enters only where a covariance, an expectation or a
distribution is stated, and systematically in [Section 9.6](06-expected-mean-squares.html).

## Orthogonal decompositions of the model space

::: {#def-ss-anova-decomposition}
[Analysis of variance decomposition]

An **analysis of variance decomposition** of the model space is a list of subspaces
\( \mathcal V_1,\dots,\mathcal V_k \) of \( \C(\X) \) that are mutually orthogonal and whose sum is
\( \C(\X) \):
\[
\C(\X)=\mathcal V_1\dirsum\mathcal V_2\dirsum\dots\dirsum\mathcal V_k,
\qquad \mathcal V_i\perp\mathcal V_j\ (i\neq j).
\]
The subspace \( \mathcal V_{k+1}=\C(\X)\perpc \) is the **error space**. Writing \( \bP_i \) for the
orthogonal projection onto \( \mathcal V_i \), the **sum of squares** for \( \mathcal V_i \) is
\( \norm{\bP_i\y}^2=\y\T\bP_i\y \), its **degrees of freedom** are \( r_i=\dim\mathcal V_i \), and its
**mean square** is \( \norm{\bP_i\y}^2/r_i \).
:::

The same space can be decomposed in many ways, and the choice is the statistician's: the
constants, a regressor of interest, a block of nuisance regressors, the linear and quadratic parts
of a trend. The theorem collects what is true of all decompositions.

::: {#thm-ss-decomposition}
[The analysis of variance identity]

Let \( \mathcal V_1,\dots,\mathcal V_k \) be an analysis of variance decomposition of \( \C(\X) \),
with projections \( \bP_1,\dots,\bP_k \) and dimensions \( r_1,\dots,r_k \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \bP_i\bP_j=\bzero \) for \( i\neq j \), \( \bP_1+\dots+\bP_k=\M \), and
   \[
\I=\bP_1+\bP_2+\dots+\bP_k+(\I-\M),
\]
   a sum of \( k+1 \) mutually orthogonal projections;

2. for every \( \y\in\Real^n \) the vectors \( \bP_1\y,\dots,\bP_k\y,(\I-\M)\y \) are mutually
   orthogonal and
   \[
\norm{\y}^2=\norm{\bP_1\y}^2+\dots+\norm{\bP_k\y}^2+\norm{(\I-\M)\y}^2 ;
\]{#eq-ss-anova-identity}

3. \( r_i=\rank(\bP_i)=\tr(\bP_i) \), \( r_1+\dots+r_k=r \), and the error space has dimension
   \( \tr(\I-\M)=n-r \);

4. if the columns of \( \Q_i \) (\( n\times r_i \)) form an orthonormal basis of \( \mathcal V_i \), then
   \( \norm{\bP_i\y}^2=\norm{\Q_i\T\y}^2 \) is a sum of \( r_i \) squares, and together with an
   orthonormal basis \( \Q_{k+1} \) of the error space the matrix \( \bU=[\Q_1,\dots,\Q_{k+1}] \) is
   orthogonal;

5. the subspaces \( \mathcal W_j=\mathcal V_1\dirsum\dots\dirsum\mathcal V_j \) form a chain
   \( \mathcal W_1\subseteq\dots\subseteq\mathcal W_k=\C(\X) \) with projections
   \( \bm N_j=\bP_1+\dots+\bP_j \), and \( \bP_j=\bm N_j-\bm N_{j-1} \) (with \( \bm N_0=\bzero \)).
   Conversely, every chain \( \{\bzero\}=\mathcal W_0\subseteq\mathcal W_1\subseteq\dots\subseteq\mathcal W_k=\C(\X) \)
   gives an analysis of variance decomposition with \( \mathcal V_j=\mathcal W_{j-1}\perpc\cap\mathcal W_j \).
:::

:::

::: {.proof}
(a) For any \( \bu \), \( \bP_j\bu\in\mathcal V_j\subseteq\mathcal V_i\perpc \), so
\( \bP_i\bP_j\bu=\bzero \) by @prp-proj-trace-rank(c). Applying @thm-proj-sum repeatedly, the sum
\( \bP_1+\dots+\bP_k \) is the orthogonal projection onto
\( \mathcal V_1\dirsum\dots\dirsum\mathcal V_k=\C(\X) \), and a projection is determined by its range,
so the sum is \( \M \). Each \( \bP_i \) satisfies \( \M\bP_i=\bP_i \), because its columns lie in
\( \C(\X) \), and therefore \( (\I-\M)\bP_i=\bzero \). The projection \( \I-\M \) onto \( \C(\X)\perpc \) (@thm-proj-sym-idem)
is thus orthogonal to every \( \bP_i \).

(b) Multiply \( \I=\sum_i\bP_i+(\I-\M) \) by \( \y \). For two different pieces,
\( (\bP_i\y)\T(\bP_j\y)=\y\T\bP_i\bP_j\y=0 \), using symmetry and (a). Pythagoras gives the identity.

(c) Rank equals trace for projections (@prp-proj-trace-rank(b)), and the trace is additive, so
\( \sum_ir_i=\tr(\M)=r \) by @thm-proj-M-formula, and \( \tr(\I-\M)=n-r \).

(d) By @prp-proj-orthonormal-formula, \( \bP_i=\Q_i\Q_i\T \), so
\( \y\T\bP_i\y=\norm{\Q_i\T\y}^2 \). Columns taken from different \( \Q_i \) are orthogonal because
the subspaces are, so \( \bU \) has orthonormal columns. It has
\( r_1+\dots+r_k+(n-r)=n \) of them, so it is square and hence orthogonal.

(e) The \( \bm N_j \) are projections onto \( \mathcal W_j \) by @thm-proj-sum, and the chain is
nested by construction. Conversely, given a chain, @thm-proj-nested(b) shows that
\( \bm N_j-\bm N_{j-1} \) is the projection onto \( \mathcal V_j=\mathcal W_{j-1}\perpc\cap\mathcal W_j \).
For \( i<j \) we have \( \mathcal V_i\subseteq\mathcal W_i\subseteq\mathcal W_{j-1} \), which is
orthogonal to \( \mathcal V_j \). The projections telescope to \( \bm N_k=\M \), so by
@thm-proj-sum the orthogonal sum of the \( \mathcal V_j \) is \( \C(\M)=\C(\X) \).
:::

Part (e) says that decompositions and chains are the same thing. Listing the rows of an
analysis of variance table from the top is the same as fitting a sequence of ever larger
models and recording what each step adds. [Section 9.3](03-sequential-partial.html) builds on this.

## Degrees of freedom are dimensions

“Degrees of freedom” can seem like a rule of thumb (“subtract one for every parameter
estimated”). Part (d) of @thm-ss-decomposition gives the phrase an exact meaning. Rotate the data into the
orthonormal coordinates \( \bz=\bU\T\y \). The rotation preserves length, and it splits
\( \bz \) into consecutive blocks \( \Q_1\T\y,\dots,\Q_{k+1}\T\y \) of sizes
\( r_1,\dots,r_k,n-r \). The sum of squares for \( \mathcal V_i \) is the sum of the squares of the
\( r_i \) coordinates in its block. *A sum of squares has \( r_i \) degrees of freedom because it is
a sum of \( r_i \) squares of orthonormal linear combinations of the data, and it cannot be
written with fewer.*

The rotation also explains why degrees of freedom govern the size of a sum of squares. If
\( \Cov(\Y)=\sigma^2\I \), the coordinates \( \bU\T\Y \) are uncorrelated with variance \( \sigma^2 \)
(@thm-rv-linear), so a block of \( r_i \) of them has expected sum of squares \( r_i\sigma^2 \) plus
the squared length of its mean ([Section 9.6](06-expected-mean-squares.html)). The divisor
\( n-1 \) of the sample variance is the case \( \Real^n=\spn(\bone)\dirsum\bone\perpc \): the \( n \)
deviations \( y_i-\bar{y} \) lie in an \( (n-1) \)-dimensional subspace, and the Helmert matrix of
@exm-mat-helmert writes their sum of squares as \( n-1 \) squared orthonormal contrasts.

## The regression analysis of variance table

The most common decomposition comes from the chain
\( \spn(\bone)\subseteq\C(\X) \), available whenever \( \bone\in\C(\X) \), which is the case for every
model with an intercept. Write \( \bP_1=n^{-1}\bone\bone\T \) for the projection onto the constant
vectors (@exm-proj-line). The chain gives the three-way decomposition of @eq-proj-three-way,
\[
\I=\bP_1+(\M-\bP_1)+(\I-\M),
\]
with ranks \( 1 \), \( r-1 \) and \( n-r \). The squared lengths are
\[
\norm{\bP_1\y}^2=n\bar{y}^2,\qquad
\text{SSR}=\norm{(\M-\bP_1)\y}^2=\sum_i(\hat{y}_i-\bar{y})^2,\qquad
\text{SSE}=\norm{(\I-\M)\y}^2=\sum_i(y_i-\hat{y}_i)^2 .
\]
They are usually arranged as follows.

| Source | df | Sum of squares | Mean square |
|---|---|---|---|
| Mean | \( 1 \) | \( n\bar{y}^2 \) | \( n\bar{y}^2 \) |
| Regression (after the mean) | \( r-1 \) | \( \text{SSR} \) | \( \text{SSR}/(r-1) \) |
| Residual | \( n-r \) | \( \text{SSE} \) | \( s^2=\text{SSE}/(n-r) \) |
| Total | \( n \) | \( \y\T\y \) | |

This is the **uncorrected** table. Its total is the raw sum of squares \( \sum_iy_i^2 \), with
\( n \) degrees of freedom. The row for the mean is rarely of interest, because it records only
that the average response is not zero. Deleting it gives the **corrected** table, whose total is
\( \text{SST}=\sum_i(y_i-\bar{y})^2=\norm{(\I-\bP_1)\y}^2 \) on \( n-1 \) degrees of freedom. The
corrected identity \( \text{SST}=\text{SSR}+\text{SSE} \) is @eq-proj-corrected-decomposition.
“Correcting for the mean” is thus a projection. We first project the data onto
\( \bone\perpc \) and then decompose what is left.

The ratio of the two mean squares,
\[
F=\frac{\text{SSR}/(r-1)}{\text{SSE}/(n-r)} ,
\]
compares the variation along the regression directions with the variation in the error space.
Under normality and when no regressor matters, it has the \( F(r-1,n-r) \) distribution
(@thm-qf-nested-f). [Section 9.6](06-expected-mean-squares.html) explains why this ratio is the
right comparison, and [Chapter 11](../ch11-general-linear-hypothesis/index.html) develops the test.

A few formulas make the table easy to compute from any least squares solution.

::: {#prp-ss-computing-formulas}
[Computing the sums of squares]

Let \( \hbeta \) be any least squares estimate. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \norm{\M\y}^2=\hbeta\T\X\T\y \) and \( \text{SSE}=\y\T\y-\hbeta\T\X\T\y \);

2. if \( \bone\in\C(\X) \), then \( \text{SSR}=\hbeta\T\X\T\y-n\bar{y}^2 \);

3. if \( \X=[\bone,\X_1] \) has full column rank and \( \tilde{\X}_1=(\I-\bP_1)\X_1 \) is the matrix of
   centred regressors, then \( \text{SSR}=\hbeta_1\T\tilde{\X}_1\T\tilde{\X}_1\hbeta_1=\hbeta_1\T\tilde{\X}_1\T\y \),
   where \( \hbeta_1 \) holds the slopes.
:::

:::

::: {.proof}
(a) \( \X\hbeta=\M\y \) (@thm-proj-ls-projection), so
\( \norm{\M\y}^2=\y\T\M\y=\y\T\X\hbeta=\hbeta\T\X\T\y \). The second formula follows from
\( \text{SSE}=\y\T\y-\y\T\M\y \) (@prp-proj-fit-residual). (b) \( \M-\bP_1 \) is a projection (@thm-proj-nested), so
\( \text{SSR}=\y\T\M\y-\y\T\bP_1\y=\norm{\M\y}^2-n\bar{y}^2 \); now use (a). (c) By @lem-proj-fwl-split with
\( \X_2=\bone \), \( \M-\bP_1 \) is the projection onto \( \C(\tilde{\X}_1) \), and by
@thm-proj-fwl the slopes are the coefficients of \( (\M-\bP_1)\y \) in terms of the columns of
\( \tilde{\X}_1 \). So \( (\M-\bP_1)\y=\tilde{\X}_1\hbeta_1 \) and
\( \text{SSR}=\norm{\tilde{\X}_1\hbeta_1}^2 \). The last equality uses
\( \tilde{\X}_1\T\y=\tilde{\X}_1\T(\M-\bP_1)\y=\tilde{\X}_1\T\tilde{\X}_1\hbeta_1 \).
:::

In (a), \( \hbeta\T\X\T\y \) is the same for every solution of the normal equations, even though
\( \hbeta \) is not unique when \( r<p \). That is no accident. It equals \( \norm{\M\y}^2 \), which
depends only on the model space.

::: {#exm-ss-crime-anova}
[An analysis of variance table for murder rates]

Return to the regression of the 2009 murder rate of the \( 50 \) US states on poverty,
single parenthood and urbanization, with an intercept (@exm-proj-fwl-crime). Here
\( r=4 \) and \( \bar{y}=4.514 \). The listing computes every sum of squares from the
coordinates of \( \y \) in an orthonormal basis of \( \C(\X) \), as in @thm-ss-decomposition(d):
the first coordinate gives the mean, the next three the regression.

| Source | df | Sum of squares | Mean square |
|---|---|---|---|
| Mean | 1 | 1018.81 | 1018.81 |
| Regression | 3 | 182.80 | 60.93 |
| Residual | 46 | 101.96 | 2.217 |
| Total (uncorrected) | 50 | 1303.57 | |
| Total (corrected) | 49 | 284.76 | |

The regression mean square is \( 27.49 \) times the residual mean square. If none of the three
regressors were related to the murder rate and the errors were normal, a ratio this large
would occur with probability about \( 2.5\times 10^{-10} \). The
proportion of the corrected total that the regression accounts for is
\( 182.80/284.76=0.6419 \), the \( R^2 \) of @thm-proj-r2-cosine.
:::

```{.python .run #cell-anova-table-table}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")

y = data["murder"].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), data["poverty"], data["single"], data["urban"]])
r = np.linalg.matrix_rank(X)

Q, _ = np.linalg.qr(X)                 # orthonormal basis of C(X); first column spans 1
z = Q.T @ y                            # coordinates of My in that basis
ss_mean = z[0] ** 2                    # n * ybar^2
ss_reg = np.sum(z[1:] ** 2)            # ||(M - P1) y||^2, r - 1 squared coordinates
ss_res = y @ y - np.sum(z ** 2)        # ||(I - M) y||^2

rows = [("mean", 1, ss_mean), ("regression", r - 1, ss_reg),
        ("residual", n - r, ss_res), ("total", n, y @ y)]
for name, df, ss in rows:
    print(f"{name:11s} df = {df:3d}   SS = {ss:9.2f}   MS = {ss / df:8.3f}")
F = (ss_reg / (r - 1)) / (ss_res / (n - r))
print(f"corrected total = {ss_reg + ss_res:.2f},  F = {F:.2f}")
```

The QR factorization supplies the orthonormal basis of @thm-ss-decomposition(d); its first
column is \( \pm\bone/\sqrt n \) because Gram–Schmidt starts with the intercept
(@prp-proj-gram-schmidt). No \( n\times n \) projection matrix is formed.

## Models without an intercept

The corrected table depends on \( \bone\in\C(\X) \). Without it, the chain
\( \spn(\bone)\subseteq\C(\X) \) is not available, and the corrected “identity” fails.

::: {#prp-ss-no-intercept}
For any model matrix \( \X \), with \( \hY=\M\y \) and \( \he=(\I-\M)\y \),
\[
\sum_i(y_i-\bar{y})^2=\sum_i(\hat{y}_i-\bar{y})^2+\sum_i\hat{\varepsilon}_i^2-2\bar{y}\sum_i\hat{\varepsilon}_i .
\]
The last term vanishes for every \( \y \) iff \( \bone\in\C(\X) \).
:::

::: {.proof}
Write \( \y-\bar{y}\bone=(\hY-\bar{y}\bone)+\he \) and expand the squared length. The cross term is
\( 2(\hY-\bar{y}\bone)\T\he=2\hY\T\he-2\bar{y}\bone\T\he=-2\bar{y}\,\bone\T\he \), since
\( \hY\perp\he \). If \( \bone\in\C(\X) \), then \( \bone\T\he=0 \) (@prp-proj-fit-residual(b)).
If \( \bone\notin\C(\X) \), take \( \y=\bone \). Then \( \bar{y}=1 \) and
\( \bone\T\he=\bone\T(\I-\M)\bone=\norm{(\I-\M)\bone}^2>0 \), so the last term is not zero.
:::

Without an intercept the valid decomposition is the uncorrected one,
\( \norm{\y}^2=\norm{\M\y}^2+\text{SSE} \), whose regression row has \( r \) degrees of freedom. That is
appropriate only when the mean response is known to vanish when every regressor is zero.

::: {#exm-ss-no-intercept}
[Forcing the murder regression through the origin]

Refit the model of @exm-ss-crime-anova without the intercept. The residual sum of squares
rises from \( 101.96 \) to \( 177.73 \), so the fit is much worse. Yet the residuals now sum to
\( -8.332 \), and the would-be corrected identity is off by
\( -2\bar{y}\sum_i\hat{\varepsilon}_i=75.22 \): the corrected total \( 284.76 \) is not
\( 31.81+177.73=209.54 \). The uncentred ratio \( \norm{\hY}^2/\norm{\y}^2 \), which some
software reports as \( R^2 \) for such models, is \( 0.8637 \), much larger than the
\( 0.6419 \) of the better model with an intercept (@exr-proj-uncentred-r2).
:::

```{.python .run #cell-anova-table-nointercept}
X0 = X[:, 1:]                                        # drop the intercept column
b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
e0 = y - X0 @ b0
ybar = y.mean()
sst = np.sum((y - ybar) ** 2)
ssr_star = np.sum((X0 @ b0 - ybar) ** 2)
sse0 = e0 @ e0
print(f"SST = {sst:.2f}  but  SSR* + SSE = {ssr_star + sse0:.2f}")
print(f"sum of residuals = {e0.sum():.3f},  -2 ybar * sum = {-2 * ybar * e0.sum():.2f}")
```

## Coarsening and refining

**Coarsening** merges rows: replacing \( \mathcal V_1 \) and \( \mathcal V_2 \) by
\( \mathcal V_1\dirsum\mathcal V_2 \) gives another decomposition, whose sums of squares and degrees of
freedom are sums of the old ones. The regression row is a coarsening of the table that gives each
regressor its own row ([Section 9.3](03-sequential-partial.html)). **Refining** splits a row: any
subspace of dimension \( r_i \) can be cut into \( r_i \) orthogonal lines, each carrying one degree of
freedom. A refinement is useful only when its lines answer questions of interest. [Section 9.4](04-orthogonal-designs.html) chooses them with
orthogonal contrasts and orthogonal polynomials.

::: {.idea}
An analysis of variance table is a list of mutually orthogonal subspaces that fill
\( \Real^n \). Its sums of squares are the squared lengths of the projections of \( \y \), its
degrees of freedom are the dimensions, and both columns add up because of Pythagoras.
:::

## Exercises

### A. Check your understanding

::: {#exr-ss-simple-anova}
[A1]

For simple regression, \( \X=[\bone,\x] \) with \( \x \) not constant, show that
\( \text{SSR}=\hat{\beta}_1^2S_{xx}=S_{xy}^2/S_{xx} \) and that the \( F \) ratio of the corrected
table equals \( t^2 \), where \( t=\hat{\beta}_1\sqrt{S_{xx}}/s \).
:::

::: {#exr-ss-fill-table}
[A2]

A regression with an intercept has a corrected table in which the regression row has
\( 3 \) degrees of freedom and sum of squares \( 120 \), the residual mean square is \( 2.5 \), and
the corrected total sum of squares is \( 172.5 \). Find \( n \), the residual degrees of freedom,
the \( F \) ratio and \( R^2 \).
:::

### B. Practice

::: {#exr-ss-canonical-form}
[B1]

Let \( \bU=[\Q_1,\dots,\Q_{k+1}] \) be the orthogonal matrix of @thm-ss-decomposition(d) and
\( \bz=\bU\T\Y \). Suppose \( \E(\Y)=\bmu\in\C(\X) \) and \( \Cov(\Y)=\sigma^2\I \). Show that the
coordinates of \( \bz \) are uncorrelated with variance \( \sigma^2 \), that the last \( n-r \) of them
have mean zero, and deduce \( \E(\text{SSE})=(n-r)\sigma^2 \). Show that if in addition
\( \Y \) is normal, the coordinates are independent.
:::

::: {.solution}
By @thm-rv-linear, \( \Cov(\bz)=\bU\T(\sigma^2\I)\bU=\sigma^2\I \) and \( \E(\bz)=\bU\T\bmu \). The
last block is \( \Q_{k+1}\T\bmu \), and \( \bmu\in\C(\X) \) is orthogonal to the error space, so it
is zero. Then \( \text{SSE}=\norm{\Q_{k+1}\T\Y}^2 \) is a sum of \( n-r \) squares of variables
with mean zero and variance \( \sigma^2 \), so its expectation is \( (n-r)\sigma^2 \). Under
normality \( \bz \) is normal (@thm-mvn-linear) with diagonal covariance, and uncorrelated
jointly normal variables are independent (@thm-mvn-independence).
:::

::: {#exr-ss-negative-r2}
[B2]

For a model without an intercept, the quantity \( 1-\text{SSE}/\text{SST} \) is sometimes reported
as \( R^2 \). Show that it is negative iff the model fits worse than the intercept-only model.
Construct a data set with \( n=4 \) and one regressor for which it is below \( -1 \).
:::

### C. Going deeper

::: {#exr-ss-max-single}
[C1]

Let \( \mathcal V \) be a subspace with projection \( \bP \). Show that
\[
\norm{\bP\y}^2=\max_{\bv\in\mathcal V,\ \bv\neq\bzero}\frac{(\bv\T\y)^2}{\bv\T\bv},
\]
with the maximum attained at \( \bv=\bP\y \) when \( \bP\y\neq\bzero \). Interpret: the sum of
squares of a row with \( r_i \) degrees of freedom is the largest one-degree-of-freedom sum of
squares obtainable by refining that row. [Chapter 13](../ch13-multiplicity/02-scheffe.html) uses this fact to build Scheffé's
simultaneous intervals.
:::

::: {.solution}
For \( \bv\in\mathcal V \), \( \bv\T\y=\bv\T\bP\y \), so by Cauchy–Schwarz (@prp-mat-cauchy-schwarz)
\( (\bv\T\y)^2\le\bv\T\bv\,\norm{\bP\y}^2 \), with equality iff
\( \bv \) is a multiple of \( \bP\y \). Taking \( \bv=\bP\y\in\mathcal V \) attains the bound. A line
\( \spn(\bv) \) inside \( \mathcal V \) carries the one-degree-of-freedom sum of squares
\( (\bv\T\y)^2/\bv\T\bv \) (@exm-proj-line), so no refinement of the row gives a line with a larger
sum of squares than the whole row.
:::

