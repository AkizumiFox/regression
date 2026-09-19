# Collinearity, previewed

Collinearity, near-linear dependence among the columns of \( \X \), violates none of the assumptions of Part III.
As long as \( \X \) has full rank, every result of Parts II and III holds. What collinearity does is make some of
them unhelpful: intervals that are correct but very wide. [Chapter 26](../ch26-collinearity/index.html) develops diagnostics, and [Chapter 27](../ch27-shrinkage/index.html) remedies that trade bias for variance.

## Variance inflation

For a single coefficient, the Frisch–Waugh–Lovell theorem gave @eq-proj-vif-preview:
\[
\Var(\hat{\beta}_j)=\frac{\sigma^2}{S_{jj}}\cdot\frac1{1-R_j^2},
\]
where \( S_{jj} \) is the centred sum of squares of \( \x_j \) and \( R_j^2 \) is the coefficient of determination from regressing
\( \x_j \) on the other regressors. The first factor is the variance \( \x_j \) would give if it were orthogonal to the other
regressors. The second, the **variance inflation factor** \( \text{VIF}_j=1/(1-R_j^2) \), is the price of the part of
\( \x_j \) that the others can already explain. A variance inflation factor of \( 100 \) means the standard error is ten times what an
orthogonal design with the same spread would give.

The factors have a compact matrix form. Let \( \mathbf{R} \) be the correlation matrix of the \( k \) regressors other than the
intercept. Then \( \text{VIF}_j=[\mathbf{R}^{-1}]_{jj} \) (@exr-dep-vif-inverse), and the eigenvalues of \( \mathbf{R} \) show where the
problem lies.

::: {#prp-dep-collinear-directions}
[Where collinearity costs variance]

Let \( \X \) have full rank and let \( \X\T\X=\sum_{\ell=1}^p\lambda_\ell\bv_\ell\bv_\ell\T \) with orthonormal eigenvectors \( \bv_\ell \) and
eigenvalues \( \lambda_\ell>0 \).

::: {.enumerate options="label=(\alph*)"}
1. For every \( \mathbf{a}\in\Real^p \), \( \Var(\mathbf{a}\T\hbeta)=\sigma^2\sum_\ell (\mathbf{a}\T\bv_\ell)^2/\lambda_\ell \).

2. The fitted values at the design points have average variance \( \sigma^2p/n \), whatever the eigenvalues.

3. With \( \mathbf{R} \) as above, \( \sum_j\text{VIF}_j=\tr\mathbf{R}^{-1}=\sum_\ell 1/\lambda_\ell(\mathbf{R}) \), and every
   \( \text{VIF}_j\le1/\lambda_{\min}(\mathbf{R}) \).
:::

:::

::: {.proof}
(a) By @thm-mat-spectral, \( (\X\T\X)^{-1}=\sum_\ell \lambda_\ell^{-1}\bv_\ell\bv_\ell\T \), and
\( \Var(\mathbf{a}\T\hbeta)=\sigma^2\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a} \) (@cor-lm-linear-functions). (b) \( \Cov(\hY)=\sigma^2\M \), and
\( \tr\M=p \) (@prp-proj-trace-rank). (c) The first equality is \( \text{VIF}_j=[\mathbf{R}^{-1}]_{jj} \) summed over \( j \). The
trace is the sum of the eigenvalues of \( \mathbf{R}^{-1} \), which are \( 1/\lambda_\ell(\mathbf{R}) \). Each diagonal entry of
\( \mathbf{R}^{-1} \) is at most its largest eigenvalue, by @thm-mat-extremal-rayleigh applied to a standard basis vector.
:::

Part (a) locates the damage. A small eigenvalue inflates the variance of \( \mathbf{a}\T\hbeta \) only through the
component of \( \mathbf{a} \) along its eigenvector. Individual coefficients usually have such components, and so do
predictions at points that break the pattern of the regressors; predictions that follow the pattern do not, and by (b)
fitted values at the design points are as precise on average as in any design. Collinearity is a problem about
*attributing* an effect to one of several regressors that moved together, not about fitting.

::: {#exm-dep-longley}
[Longley's regression]

Longley's (1967) data, which @exm-cmp-longley used to test the accuracy of least squares algorithms, regress total
employment on six economic series over \( n=16 \) years. The variance inflation factors range from
\( 3.59 \) for ARMED (the size of the armed forces) to \( 1788.5 \) for GNP. The eigenvalues of the
correlation matrix run from \( 4.603 \) down to \( 0.00038 \), and the factors sum to
\( 3119.4 \), which is the sum of the reciprocal eigenvalues, as @prp-dep-collinear-directions(c)
requires. Four of the series (the price deflator, GNP, population and the year itself) rose almost in step over the sixteen
years, each with correlation above \( 0.99 \) with the year. The data can say much about their joint movement and
almost nothing about any one of them with the others held fixed.
:::

```{.python .run #cell-collinearity-vif}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
Z = data.drop(columns="TOTEMP").to_numpy()           # six regressors
n, k = Z.shape

def r_squared(v, W):
    W1 = np.column_stack([np.ones(len(v)), W])
    fitted = W1 @ np.linalg.lstsq(W1, v, rcond=None)[0]
    return 1 - np.sum((v - fitted) ** 2) / np.sum((v - v.mean()) ** 2)

vif = np.array([1 / (1 - r_squared(Z[:, j], np.delete(Z, j, axis=1))) for j in range(k)])
R = np.corrcoef(Z, rowvar=False)
print(dict(zip(data.columns[1:], np.round(vif, 1))))
print("diagonal of R^{-1}:", np.round(np.diag(np.linalg.inv(R)), 1))
print("eigenvalues of R:  ", np.round(np.linalg.eigvalsh(R), 5))
```

::: {#exm-dep-along-across}
[Predicting along and across the pattern]

Draw \( m=40 \) pairs of regressors with sample correlation \( 0.9868 \), so that each has a
variance inflation factor of \( 38.1 \), and fit \( \E(Y)=\beta_0+\beta_1x_1+\beta_2x_2 \). In units of \( \sigma \),
the two slopes have standard deviations \( 0.904 \) and \( 0.907 \). At the new
point \( (x_1,x_2)=(1,1) \), which follows the pattern of the data, the fitted mean has standard deviation
\( 0.228 \), close to the root mean square over the design points,
\( \sqrt{3/40}=0.274 \). At \( (1,-1) \), which is just as far from the origin but against the
pattern, it is \( 1.809 \). [Figure 19.6.1](#fig-dep-collinearity) shows the contours of this
standard deviation. They are long ellipses aligned with the data cloud.
:::

```{.python .run #cell-collinearity-prediction}
rng = np.random.default_rng(1906)
m, rho = 40, 0.98
x1 = rng.normal(size=m)
x2 = rho * x1 + np.sqrt(1 - rho**2) * rng.normal(size=m)
X = np.column_stack([np.ones(m), x1, x2])
G = np.linalg.inv(X.T @ X)                           # Cov(beta_hat) / sigma^2
along = np.array([1.0, 1.0, 1.0])                    # a new point that follows the pattern
across = np.array([1.0, 1.0, -1.0])                  # same distance, against the pattern
print("sd of beta_1, beta_2 (units of sigma):", np.round(np.sqrt(np.diag(G)[1:]), 3))
print("sd of the fitted mean, along :", round(np.sqrt(along @ G @ along), 3))
print("sd of the fitted mean, across:", round(np.sqrt(across @ G @ across), 3))
print("average over the design points:", round(np.trace(X @ G @ X.T) / m, 3), "= p/n")
```

::: {when-format="html"}
![**Figure 19.6.1.** Standard deviation of the fitted mean, in units of \( \sigma \), over the regressor plane for
@exm-dep-along-across, with the \( 40 \) design points. The contours are elongated along the data: prediction
at the circle, which follows the pattern, is precise; prediction at the square, which breaks it, is
not.](collinearity.svg){#fig-dep-collinearity width=55%}
:::

::: {when-format="pdf"}
![Standard deviation of the fitted mean, in units of \( \sigma \), over the regressor plane for
@exm-dep-along-across, with the \( 40 \) design points. The contours are elongated along the data: prediction
at the circle, which follows the pattern, is precise; prediction at the square, which breaks it, is
not.](collinearity.pdf){width=55%}
:::

## How collinearity interacts with the other departures

Collinearity causes no bias, but it amplifies the departures of the earlier sections. By
@eq-dep-two-regressor-ratio the stakes of omitting \( x_2 \) are \( r^2\lvert\gamma-1\rvert \), so the choice between short and
long models matters most exactly when the data are least able to make it. Perturbations of the regressors, from
rounding or measurement error, are governed by the condition number of \( \X \) (@thm-cmp-perturbation), which is large
for the same reason. And a single case can have high leverage without being extreme in any one regressor, because
leverage measures distance in the metric of the regressors' covariance (@prp-proj-leverage-mahalanobis).

::: {.idea}
Collinearity is a shortage of information, not a failure of the model, and it magnifies every other departure.
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-vif-inverse}
[A1]

Let \( \mathbf{Z} \) be the \( n\times k \) matrix of centred regressors, each scaled to unit length, so that \( \mathbf{Z}\T\mathbf{Z}=\mathbf{R} \). Use
@eq-proj-partitioned-inverse to show that \( [\mathbf{R}^{-1}]_{jj}=1/(1-R_j^2) \).
:::

::: {.solution}
By @eq-proj-partitioned-inverse with \( \X_1=\mathbf{z}_j \) and \( \X_2 \) the other columns,
\( [\mathbf{R}^{-1}]_{jj}=1/\norm{(\I-\M_{(j)})\mathbf{z}_j}^2 \), where \( \M_{(j)} \) projects onto the other centred columns. Since
\( \mathbf{z}_j \) is centred with unit length, the squared residual length is \( 1-R_j^2 \), as in @exr-proj-vif (the intercept only
removes means, and the columns are already centred).
:::

::: {#exr-dep-vif-one}
[A2]

Show that \( \text{VIF}_j\ge1 \), with equality iff the centred \( \x_j \) is orthogonal to every other centred regressor.
:::

### B. Practice

::: {#exr-dep-sum-difference}
[B1]

For two centred regressors with \( \norm{\x_1}^2=\norm{\x_2}^2=S \) and correlation \( r \), show that
\( \Var(\hat{\beta}_1+\hat{\beta}_2)=2\sigma^2/\bigl(S(1+r)\bigr) \) and \( \Var(\hat{\beta}_1-\hat{\beta}_2)=2\sigma^2/\bigl(S(1-r)\bigr) \). Relate the
result to @prp-dep-collinear-directions(a) and to [Figure 19.6.1](#fig-dep-collinearity).
:::

::: {.solution}
The slope block of \( \X\T\X \) is \( S\begin{pmatrix}1&r\\r&1\end{pmatrix} \), with eigenvalues \( S(1+r) \) and \( S(1-r) \) and eigenvectors
\( (1,1)\T/\sqrt2 \) and \( (1,-1)\T/\sqrt2 \). The intercept is orthogonal to the centred slopes and plays no role. With
\( \mathbf{a}=(1,\pm1)\T \), \( \mathbf{a}\T\bv \) is \( \sqrt2 \) on one eigenvector and \( 0 \) on the other, so part (a) gives \( 2\sigma^2/\bigl(S(1\pm r)\bigr) \).
When \( r \) is near one, the sum, the combined effect of two regressors that move together, is well determined, and the
difference is not. In the figure, the new point \( (1,1) \) asks for the sum and \( (1,-1) \) for the difference.
:::

::: {#exr-dep-vif-bounds}
[B2]

Show that \( \sum_j\text{VIF}_j\ge k \), and that \( \max_j\text{VIF}_j\ge\bigl(\sum_j\text{VIF}_j\bigr)/k \). Deduce that a correlation matrix
with a small eigenvalue must have at least one large variance inflation factor, and conversely.
:::

### C. Going deeper

::: {#exr-dep-convex-hull}
[C1]

Show that if \( \x_0 \) is a convex combination of the rows of \( \X \), then \( \Var(\x_0\T\hbeta)\le\sigma^2\max_ih_{ii} \). Is the
converse true: does a small prediction variance mean \( \x_0 \) lies inside the data? Compare with @exm-ci-hidden-extrapolation.
:::
