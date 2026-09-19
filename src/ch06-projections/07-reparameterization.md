# Reparameterization, centring and R²

## Same subspace, same fit

The key idea of [Section 6.1](01-nearest-point.html) now becomes a theorem.

::: {#thm-proj-reparam}
[Invariance under reparameterization]

Let \( \X \) (\( n\times p \)) and \( \W \) (\( n\times q \)) satisfy \( \C(\X)=\C(\W) \). Then the
two linear models \( \E(\Y)=\X\bbeta \) and \( \E(\Y)=\W\bgamma \) have the same
projection \( \M \), the same fitted values and residuals for every \( \y \), and the
same residual sum of squares. There are a \( p\times q \) matrix \( \bT \) with
\( \W=\X\bT \) and a \( q\times p \) matrix \( \bS \) with \( \X=\W\bS \). If \( \hat{\bgamma} \) is a
least squares estimate for \( \W \), then \( \bT\hat{\bgamma} \) is one for \( \X \), and
\( \bS\hbeta \) is one for \( \W \) whenever \( \hbeta \) is one for \( \X \). If both matrices
have full column rank, then \( p=q \), \( \bT=\bS^{-1} \), and \( \hat{\bgamma}=\bS\hbeta \).
:::

::: {.proof}
A projection is determined by its range (@thm-proj-sym-idem), so both models
have projection \( \M \), and fitted values and residuals follow from @thm-proj-ls-projection.
Each column of \( \W \) lies in \( \C(\X) \), so \( \W=\X\bT \)
for some \( \bT \), and likewise \( \X=\W\bS \). If \( \W\hat{\bgamma}=\M\y \), then
\( \X(\bT\hat{\bgamma})=\M\y \), so \( \bT\hat{\bgamma} \) is a least squares estimate for \( \X \).
The other direction is symmetric. With full column rank, \( \X=\X\bT\bS \) gives
\( \X(\I-\bT\bS)=\bzero \) and so \( \bT\bS=\I \), and similarly \( \bS\bT=\I \). Uniqueness
of least squares estimates then gives \( \hat{\bgamma}=\bS\hbeta \).
:::

Many operations that change the printed coefficient table leave the fit exactly
unchanged:

- **Changing units** of a regressor (metres to centimetres, dollars to
          thousands of dollars) rescales one column: \( \W=\X\bD \) with \( \bD \) diagonal.
          The coefficient is divided by the scale factor, and nothing else changes.

- **Centring** regressors, \( \x_j\mapsto\x_j-\bar x_j\bone \), when the model has
          an intercept. The slopes are unchanged, and the intercept becomes the fitted
          value at the mean of the regressors instead of at zero, which is often
          far outside the data.

- **Recoding a factor**. Reference-level dummies, sum-to-zero effect
          coding, and cell-means coding without an intercept all span the space of
          vectors that are constant within levels.

- **Replacing polynomial terms by orthogonal polynomials**
          (Chapter 42), which span the same space with much
          better numerical behaviour ([Section 6.10](10-computation.html)).

Operations that *do* change the fit are exactly those that change
\( \C(\X) \). Examples are dropping the intercept when no other combination of
columns produces \( \bone \), replacing a regressor by its logarithm, and adding an
interaction.

## Centring as a projection

Suppose \( \bone\in\C(\X) \). By @thm-proj-nested with \( \X_0=\bone \), the fitted vector
splits as
\[
\M\y=\bar y\bone+(\M-n^{-1}\bone\bone\T)\y ,
\]
and the second projection is onto \( \C(\X)\cap\bone\perpc \), the vectors in the model
space with entries summing to zero. By @lem-proj-fwl-split with \( \X_2=\bone \),
that subspace is spanned by the centred columns
\( (\I-n^{-1}\bone\bone\T)\X \). Applying Pythagoras to the orthogonal pieces
\( \hY-\bar y\bone \) and \( \he \) of \( \y-\bar y\bone \) gives the
*corrected* decomposition
\[
\underbrace{\sum_i(y_i-\bar y)^2}_{\text{SST}}
=\underbrace{\sum_i(\hat y_i-\bar y)^2}_{\text{SSR}}
+\underbrace{\sum_i(y_i-\hat y_i)^2}_{\text{SSE}} ,
\]{#eq-proj-corrected-decomposition}

with dimensions \( n-1=(r-1)+(n-r) \). [Chapter 9](../ch09-sums-of-squares/index.html) builds the analysis of
variance on this identity.

## \( R^2 \) is a squared cosine

The **coefficient of determination** of a model with \( \bone\in\C(\X) \) is
\[
R^2=\frac{\text{SSR}}{\text{SST}}=1-\frac{\text{SSE}}{\text{SST}} .
\]
It depends on the model only through \( \C(\X) \), so by @thm-proj-reparam it
is invariant under every reparameterization. Its geometric meaning is simple.

::: {#thm-proj-r2-cosine}
Assume \( \bone\in\C(\X) \) and \( \y\notin\spn(\bone) \). Let \( \theta \) be the angle between the
centred response \( \y-\bar y\bone \) and the centred fitted vector \( \hY-\bar y\bone \).
Then
\[
R^2=\cos^2\theta=\operatorname{corr}(\y,\hY)^2 ,
\]
where \( \operatorname{corr} \) is the sample correlation. For simple regression on one
regressor \( \x \), \( R^2=\operatorname{corr}(\y,\x)^2 \).
:::

::: {.proof}
Write \( \bu=\y-\bar y\bone \) and \( \bv=\hY-\bar y\bone \). Then \( \bu=\bv+\he \) with
\( \bv\perp\he \), so \( \bu\T\bv=\norm{\bv}^2 \) and
\[
\cos^2\theta=\frac{(\bu\T\bv)^2}{\norm{\bu}^2\norm{\bv}^2}
=\frac{\norm{\bv}^2}{\norm{\bu}^2}=\frac{\text{SSR}}{\text{SST}} .
\]
(If \( \bv=\bzero \) the angle is undefined, but \( R^2=0 \) and the correlation is taken to be zero.) Since
\( \hY \) has mean \( \bar y \) (@prp-proj-fit-residual(b)), \( \cos\theta \) is the
sample correlation of \( \y \) and \( \hY \). In simple regression
\( \bv=(S_{xy}/S_{xx})(\x-\bar{x}\bone) \) by @exm-proj-simple-nested, a scalar
multiple of the centred regressor. Rescaling a vector by a nonzero constant does
not change a squared cosine, so \( \cos^2\theta \) is the squared correlation of \( \y \)
and \( \x \). When \( S_{xy}=0 \) both sides are zero.
:::

::: {when-format="html"}
![**Figure 6.7.1.** \( R^2 \) as a squared cosine, for \( n=3 \) and one regressor. Centred vectors
lie in the plane \( \bone\perpc \), so the whole configuration is two-dimensional. The
centred fitted vector is the projection of the centred response onto the line
spanned by the centred regressor.](r2_cosine.svg){#fig-proj-r2 width=52%}
:::

::: {when-format="pdf"}
![\( R^2 \) as a squared cosine, for \( n=3 \) and one regressor. Centred vectors
lie in the plane \( \bone\perpc \), so the whole configuration is two-dimensional. The
centred fitted vector is the projection of the centred response onto the line
spanned by the centred regressor.](r2_cosine.pdf){width=52%}
:::

So \( R^2 \) measures how closely the model space lines up with the centred
data. [Figure 6.7.1](07-reparameterization.html#fig-proj-r2) draws the case \( n=3 \). For the
three-regressor model of @exm-proj-fwl-crime, \( R^2=0.6419 \), meaning the centred
response makes an angle of \( 36.8^\circ \) with the model space. The
simple regression on poverty alone has \( R^2=0.3599 \), the square of
the correlation \( 0.5999 \).

```{.python .run #cell-r2-cosine-cosine}
import numpy as np
import statsmodels.api as sm

def r_squared_as_cosine(y, X):
    """X must contain the intercept column. Returns (R^2 from sums of squares, cos^2 of angle)."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    yc, yhc = y - y.mean(), y_hat - y.mean()          # centre both vectors
    r2 = 1 - np.sum((y - y_hat) ** 2) / np.sum(yc ** 2)
    cos = (yc @ yhc) / (np.linalg.norm(yc) * np.linalg.norm(yhc))
    return r2, cos ** 2
```

Three consequences follow directly from the geometry.

1. *\( R^2 \) never decreases when the model space grows.* If
           \( \C(\X_0)\subseteq\C(\X) \) and both contain \( \bone \), then
           \( \text{SSE}_0-\text{SSE}=\norm{(\M-\Mo)\y}^2\ge0 \) by @eq-proj-extra-ss.
           Adding a regressor of pure noise increases \( R^2 \) almost surely. This is why
           \( R^2 \) cannot be used for model selection (Chapter 29).

2. *\( R^2=1 \) once the model space contains \( \y \)*, which always happens when
           \( \rank(\X)=n \). A model with as many free parameters as observations
           interpolates, and its \( R^2 \) says nothing.

3. *Without an intercept, the definition changes.* If \( \bone\notin\C(\X) \),
           the corrected decomposition @eq-proj-corrected-decomposition fails,
           because \( \hY-\bar y\bone \) and \( \he \) need not be orthogonal. Software then
           usually reports the *uncentred*
           \( R_0^2=\norm{\hY}^2/\norm{\y}^2 \), the squared cosine of the angle
           between \( \y \) and \( \C(\X) \) itself. It is typically much larger and cannot be
           compared with the centred \( R^2 \) of a model with an intercept
           (@exr-proj-uncentred-r2).

## Exercises

### A. Check your understanding

::: {#exr-proj-baseline-change}
[A1]

A factor with three levels is coded with an intercept and indicators for levels 2 and 3.
Write down the matrices \( \bT \) and \( \bS \) of @thm-proj-reparam that convert to the
coding in which level 3, rather than level 1, is the omitted baseline, and to sum-to-zero
effect coding. Express each new coefficient in terms of the old ones.
:::

::: {#exr-proj-r2-affine}
[A2]

Show that, for a model with \( \bone\in\C(\X) \), \( R^2 \) is unchanged when \( \y \) is replaced by
\( a\bone+b\y \) with \( b\ne0 \), and when any column of \( \X \) is replaced by an affine
function of itself. Is it unchanged when a regressor is replaced by its logarithm?
:::

### B. Practice

::: {#exr-proj-uncentred-r2}
[B1]

Show that if \( \bone\in\C(\X) \) then
\( R_0^2=\norm{\hY}^2/\norm{\y}^2=(n\bar y^2+\text{SSR})/(n\bar y^2+\text{SST})\ge R^2 \).
Construct a data set and a model without an intercept for which \( R_0^2>0.95 \) but the
model's residual sum of squares exceeds that of the intercept-only model.
:::

::: {.solution}
With \( \bone\in\C(\X) \), the fitted vector splits
orthogonally as \( \bar y\bone+(\hY-\bar y\bone) \), so \( \norm{\hY}^2=n\bar y^2+\text{SSR} \), and
likewise \( \norm{\y}^2=n\bar y^2+\text{SST} \). For \( a\ge0 \) and \( 0\le\text{SSR}\le\text{SST} \),
\( (a+\text{SSR})/(a+\text{SST})\ge\text{SSR}/\text{SST} \). For the example, take responses
near \( 100 \) with small noise and a regressor unrelated to them, and fit through the origin.
The fitted line captures the large mean, so \( R_0^2 \) is near \( 1 \), but its residuals exceed
those around \( \bar y \).
:::

### C. Going deeper

::: {#exr-proj-r2-null-distribution}
[C1]

Suppose \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \) and \( \bone\in\C(\X) \) with \( \rank(\X)=r \).
Using the results of [Chapter 4](../ch04-quadratic-forms/index.html), show that
\( R^2\sim\mathrm{Beta}\bigl((r-1)/2,(n-r)/2\bigr) \), so that \( \E(R^2)=(r-1)/(n-1) \) even though
no regressor is related to the response. Use this to motivate the adjusted
\( R^2=1-(1-R^2)(n-1)/(n-r) \).
:::
