# Generalized least squares

The model of this chapter keeps the mean structure of every earlier chapter and
loosens the second moments:
\[
\E(\Y)=\X\bbeta,\qquad \Cov(\Y)=\sigma^2\V ,
\]{#eq-ggm-model}

where \( \X \) is \( n\times p \) of rank \( r \), the matrix \( \V \) is known and
nonnegative definite, and \( \sigma^2>0 \) is unknown. This section and the next assume
\( \V \) nonsingular; [Section 31.3](03-singular-covariance.html) removes that
assumption.

## What is assumed and what is not

The split of the covariance into a known matrix \( \V \) and an unknown scalar
\( \sigma^2 \) is a convention, not extra knowledge. Only the product \( \sigma^2\V \)
enters @eq-ggm-model, so \( (\sigma^2,\V) \) and \( (\sigma^2/c,\,c\V) \) describe the
same model for every \( c>0 \); it is usual to normalize \( \V \) so that
\( \tr(\V)/n=1 \), making \( \sigma^2 \) the variance of a single observation on a
reference scale. What the model really assumes is that the *shape* of the
covariance — the ratios of variances and all the correlations — is known.

That is exactly true only when the design forces the shape, as when an
observation is the mean of \( m_i \) replicates and has variance \( \sigma^2/m_i \)
([Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html)),
or when the construction of the data does, as under the linear transformations of
[Section 31.3](03-singular-covariance.html). Otherwise \( \V \) belongs to a family
\( \V(\boldsymbol{\uptheta}) \) estimated first and then treated as known, and this section's
exactness is a fiction whose price is computed in
[Section 31.4](04-feasible-gls.html).

## The estimator and the whitened model

Generalized least squares was defined in @def-proj-gls as projection in the
\( \V^{-1} \) inner product; we repeat it in the notation of this part.

::: {#def-ggm-gls}
[Generalized least squares]

Under @eq-ggm-model with \( \V \) positive definite, a **generalized least squares
estimate** is any \( \hbeta_{\V} \) minimizing
\[
Q_{\V}(\bb)=(\y-\X\bb)\T\V^{-1}(\y-\X\bb),
\]
equivalently any solution of the **generalized normal equations**
\( \X\T\V^{-1}\X\bb=\X\T\V^{-1}\y \). The fitted vector is \( \X\hbeta_{\V}=\A_{\V}\y \),
where
\[
\A_{\V}=\X(\X\T\V^{-1}\X)\ginv\X\T\V^{-1}
\]{#eq-ggm-gls-projection}

for any generalized inverse, and the **generalized residual** is
\( \he_{\V}=\y-\X\hbeta_{\V} \).
:::

Everything algebraic about \( \A_{\V} \) is @thm-proj-A-projection with \( \A=\V^{-1} \):
the choice of generalized inverse is immaterial, it is idempotent with
column space \( \C(\X) \), and \( \V^{-1}\A_{\V} \) is symmetric although \( \A_{\V} \)
is not — in the Euclidean geometry it is an oblique projection (@prp-proj-oblique).
Everything statistical comes from one change of variables. Let \( \V^{-1/2} \) be the symmetric positive definite square root of \( \V^{-1} \) (@thm-mat-square-root)
and set
\[
\Y_*=\V^{-1/2}\Y,\qquad \X_*=\V^{-1/2}\X .
\]{#eq-ggm-whitened}

By @thm-rv-linear, \( \E(\Y_*)=\X_*\bbeta \) and \( \Cov(\Y_*)=\sigma^2\I \): the
**whitened model** satisfies the assumptions of Part II exactly. Since
\( Q_{\V}(\bb)=\norm{\Y_*-\X_*\bb}^2 \), generalized least squares for \( (\Y,\X) \) *is*
ordinary least squares for \( (\Y_*,\X_*) \). Because \( \V^{-1/2} \) is nonsingular,
\( \rank(\X_*)=r \) and \( \C(\X_*\T)=\C(\X\T) \), so by @thm-est-characterization the two
models have the same estimable functions: nothing about identifiability changes.
Any nonsingular \( \bL \) with \( \V=\bL\bL\T \) serves equally well, by
@prp-rv-whitening(b); the Cholesky factor is the usual choice in computation.

Whitening is rarely carried out by forming \( \V^{-1/2} \); structured \( \V \) almost
always has a shortcut. If the observations fall into \( G \) clusters of
size \( m \) that share a random shock, so that \( \V=(1-\rho)\I+\rho\Z\Z\T \) with
\( \Z \) the \( n\times G \) matrix of cluster indicators, then @thm-mat-woodbury gives
\[
\V^{-1}=\frac{1}{1-\rho}\Bigl(\I-\frac{\rho}{1+(m-1)\rho}\Z\Z\T\Bigr),
\]{#eq-ggm-equicorr-inverse}

and the whitening is the **quasi-demeaning**
\[
y_{kj}\ \longmapsto\ \frac{y_{kj}-c\,\bar y_{k}}{\sqrt{1-\rho}},
\qquad c=1-\sqrt{\frac{1-\rho}{1+(m-1)\rho}},
\]{#eq-ggm-quasi}

where \( \bar y_k \) is the mean of cluster \( k \) (@exr-ggm-quasi). At \( \rho=0 \),
\( c=0 \) and nothing happens; at \( \rho=0.5 \) with \( m=12 \), \( c=0.723 \);
as \( \rho\to1 \), \( c\to1 \) and each cluster is fully centred, throwing away
everything the shared shock contaminated. Generalized least squares interpolates,
and [Chapter 32](../ch32-linear-mixed-models/index.html) will recognize
@eq-ggm-quasi as the transformation a random intercept induces.

## Exact inference when the shape is known

::: {#thm-ggm-inference}
[Inference under a known covariance shape]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V) \) with \( \V \) known and positive
definite, \( r=\rank(\X)<n \). Put
\[
\text{SSE}_{\V}=\he_{\V}\T\V^{-1}\he_{\V}=\Y\T\V^{-1}(\I-\A_{\V})\Y,
\qquad s_{\V}^2=\frac{\text{SSE}_{\V}}{n-r}.
\]{#eq-ggm-sse}

::: {.enumerate options="label=(\alph*)"}
1. For \( \blambda\in\C(\X\T) \), the statistic \( \blambda\T\hbeta_{\V} \) is the same for
   every solution of the generalized normal equations, and
   \[
\blambda\T\hbeta_{\V}\sim\Normal\bigl(\blambda\T\bbeta,\ \sigma^2\blambda\T(\X\T\V^{-1}\X)\ginv\blambda\bigr).
\]
   It is the best linear unbiased estimator of \( \blambda\T\bbeta \) (@cor-opt-aitken).

2. \( \text{SSE}_{\V}/\sigma^2\sim\chi^2(n-r) \), independently of \( \X\hbeta_{\V} \).

3. \( \E(s_{\V}^2)=\sigma^2 \) and \( \Var(s_{\V}^2)=2\sigma^4/(n-r) \).

4. For estimable \( \blambda\T\bbeta \) with \( \blambda\ne\bzero \),
   \[
\frac{\blambda\T\hbeta_{\V}-\blambda\T\bbeta}{s_{\V}\sqrt{\blambda\T(\X\T\V^{-1}\X)\ginv\blambda}}\sim t(n-r).
\]

5. Let \( \C(\X_0)\subseteq\C(\X) \) with \( \rank(\X_0)=r_0<r \), and let
   \( \text{SSE}_{\V,0} \) be @eq-ggm-sse computed in the reduced model. Then
   \[
F=\frac{(\text{SSE}_{\V,0}-\text{SSE}_{\V})/(r-r_0)}{s_{\V}^2}\sim F(r-r_0,\,n-r,\,\gamma),
\]{#eq-ggm-f}

   with \( \gamma=\norm{(\M_*-\M_{*0})\X_*\bbeta}^2/\sigma^2 \), where \( \M_* \) and
   \( \M_{*0} \) are the orthogonal projections onto \( \C(\V^{-1/2}\X) \) and
   \( \C(\V^{-1/2}\X_0) \). The noncentrality is zero iff \( \X\bbeta\in\C(\X_0) \).
:::
:::

::: {.proof}
Pass to the whitened model @eq-ggm-whitened, which satisfies the normal linear
model @eq-opt-normal-model with model matrix \( \X_* \) of rank \( r \). As shown above,
\( \hbeta_{\V} \) is exactly a least squares estimate for \( (\Y_*,\X_*) \), and the two
models have the same estimable functions. Moreover
\[
\text{SSE}_{\V}=\norm{\Y_*-\X_*\hbeta_{\V}}^2=\norm{(\I-\M_*)\Y_*}^2
\]
is the residual sum of squares of the whitened problem, because
\( \V^{-1/2}\he_{\V}=(\I-\M_*)\Y_* \) by @thm-proj-ls-projection, and
\( \X\T\V^{-1}\X=\X_*\T\X_* \).

Now apply @thm-opt-sampling to \( (\Y_*,\X_*) \). Part (a) of that theorem gives (a)
here, since \( \blambda\T(\X_*\T\X_*)\ginv\blambda=\blambda\T(\X\T\V^{-1}\X)\ginv\blambda \)
and invariance of the estimate follows from @thm-proj-invariant-functions applied
in the whitened model. Parts (b) and (c) of @thm-opt-sampling give (b)
here, noting that \( \X\hbeta_{\V}=\V^{1/2}\M_*\Y_* \) is a function of \( \M_*\Y_* \),
and part (d) gives (c).
Part (d) follows from (a) and (b) exactly as @thm-ci-estimable-interval does in the
Euclidean case. For (e), \( \C(\X_0)\subseteq\C(\X) \) implies
\( \C(\V^{-1/2}\X_0)\subseteq\C(\V^{-1/2}\X) \) with ranks \( r_0<r \), so
@def-glh-reduced-model holds for the whitened pair and @thm-glh-f-test applies;
\( \text{SSE}_{\V,0}=\norm{(\I-\M_{*0})\Y_*}^2 \) by the same computation. The
noncentrality @eq-glh-noncentrality is
\( \norm{(\M_*-\M_{*0})\X_*\bbeta}^2/\sigma^2 \), and it vanishes iff
\( \X_*\bbeta\in\C(\X_{*0}) \), that is, iff \( \X\bbeta\in\C(\X_0) \) since
\( \V^{-1/2} \) is nonsingular.
:::

The theorem is not a new theory: it is Parts II and III read in another
coordinate system, and the proof is nothing but the dictionary. Under normality
\( \hbeta_{\V} \) also maximizes the likelihood (@exr-opt-gls-mle), and \( s_{\V}^2 \)
has smallest variance among unbiased estimates of \( \sigma^2 \), by
[Section 7.4](../ch07-optimality/04-minimum-variance.html) applied to the whitened
model.

## What changes, and what does not

The dictionary is exact, but several habits from Parts II and III do not survive
translation.

::: {.warning}
Under @eq-ggm-model with \( \V\ne\I \):

- The generalized residual is **not** orthogonal to \( \C(\X) \); what holds is
  \( \X\T\V^{-1}\he_{\V}=\bzero \). The entries of \( \he_{\V} \) need not sum to zero
  even when the model has an intercept.
- \( \text{SSE}_{\V} \) is **not** \( \sum_i(y_i-\hat y_i)^2 \), and the ordinary
  residual sum of squares of the generalized fit is necessarily *larger* than
  that of the ordinary fit, which minimizes it. A model comparison using
  \( \sum_i(y_i-\hat y_i)^2 \) always prefers ordinary least squares and proves
  nothing.
- \( R^2 \) has no agreed meaning: the squared-cosine reading of
  @thm-proj-r2-cosine holds in the \( \V^{-1} \) geometry, not the Euclidean one.
- The diagonal entries of \( \A_{\V} \) still sum to \( r \), but they are not
  leverages in the sense of @prp-proj-leverage and can fall outside \( [0,1] \) (@exr-ggm-oblique-leverage).
:::

What survives is everything inference needs: unbiasedness, optimality, the
\( \chi^2 \), \( t \) and \( F \) distributions, the degrees of freedom \( n-r \), and
\( \gamma \) as a squared distance from the reduced model space, measured now in the
\( \V^{-1} \) metric.

::: {#exm-ggm-elnino}
[A trend in sea surface temperature]

The El Niño data record the monthly averaged sea surface temperature of a
region of the eastern Pacific for the \( 61 \) complete years \( 1950 \) to
\( 2010 \), \( 732 \) observations in \( 61 \) blocks of \( 12 \). Fit a
linear trend in time together with three harmonic pairs at one, two and three
cycles per year, so \( p=8 \). Warm and cold years are warm and cold in every
month, so the natural working model is a random annual shock: equicorrelation
\( \rho \) inside each year, independence across years. Take \( \rho=0.5 \) as known
for the moment.

Ordinary least squares gives a trend of \( 0.1349 \) degrees per decade with the
usual standard error \( 0.0223 \); generalized least squares gives
\( 0.1347 \) with standard error \( 0.0512 \), larger by a factor of
\( 2.30 \). The estimates agree to three decimals, the standard errors not at
all. The annual cycle tells the same story in reverse: the coefficient of
\( \sin(2\pi\,\text{month}/12) \) is \( 2.7629 \) under both fits, with
standard error \( 0.0554 \) from ordinary and \( 0.0354 \) from
generalized least squares — this time the usual standard error is too *large*.
[Section 31.2](02-ols-blue.html) explains both facts.

The \( F \) test of "no trend" is \( 36.7 \) if the annual shock is ignored
and \( 6.92 \) once it is allowed for, on \( 1 \) and \( 724 \) degrees of freedom
(\( p=0.0087 \)): still clear evidence, but with a tenth of the apparent weight.
The test that the second and third harmonics are absent moves the other way, from
\( 10.4 \) to \( 25.6 \), because the annual shock inflates \( s^2 \) and so
deflates every \( F \) whose numerator lives inside years.
[Figure 31.1.1](#fig-ggm-elnino)(a) shows the annual mean residuals, far too
spread out for independent errors, and panel (b) the two confidence bands.
:::

```{.python .run #cell-elnino-gls-data}
import numpy as np
import statsmodels.api as sm
from scipy import stats

frame = sm.datasets.elnino.load_pandas().data
years = frame["YEAR"].to_numpy()
temps = frame.iloc[:, 1:].to_numpy()                  # 61 years by 12 months
G, m = temps.shape
n = G * m
y = temps.ravel()                                     # the month varies fastest
year = np.repeat(years, m)
month = np.tile(np.arange(1, m + 1), G)
time = (year - 1980.0) + (month - 6.5) / 12.0         # decimal years, centred at 1980

cols = [np.ones(n), time / 10.0]                      # the trend is per decade
for k in (1, 2, 3):
    cols += [np.cos(2 * np.pi * k * month / m), np.sin(2 * np.pi * k * month / m)]
X = np.column_stack(cols)
p = X.shape[1]
```

```{.python .run #cell-elnino-gls-gls}
rho = 0.5                                             # taken as known for now


def gls_fit(X, y, rho, m):
    """GLS with equicorrelation rho inside blocks of m consecutive observations.

    Whitening is the quasi-demeaning y -> (y - c ybar_g)/sqrt(1 - rho) with
    c = 1 - sqrt((1 - rho)/(1 + (m - 1) rho)); ordinary least squares on the whitened
    data is generalized least squares on the original data.
    """
    c = 1.0 - np.sqrt((1 - rho) / (1 + (m - 1) * rho))

    def whiten(a):
        b = a.reshape(len(a) // m, m, -1) if a.ndim > 1 else a.reshape(-1, m)
        return ((b - c * b.mean(axis=1, keepdims=True)) / np.sqrt(1 - rho)).reshape(a.shape)

    Xw, yw = whiten(X), whiten(y)
    b = np.linalg.solve(Xw.T @ Xw, Xw.T @ yw)
    rss = float((yw - Xw @ b) @ (yw - Xw @ b))         # the V^{-1} residual sum of squares
    s2 = rss / (len(y) - X.shape[1])
    return b, s2 * np.linalg.inv(Xw.T @ Xw), rss, c


b_ols = np.linalg.solve(X.T @ X, X.T @ y)
rss_ols = float((y - X @ b_ols) @ (y - X @ b_ols))
V_ols = rss_ols / (n - p) * np.linalg.inv(X.T @ X)     # the standard errors usually reported
b_gls, V_gls, rss_gls, cq = gls_fit(X, y, rho, m)

print(f"trend per decade: OLS {b_ols[1]:+.4f}  GLS {b_gls[1]:+.4f}")
print(f"standard error:   OLS {np.sqrt(V_ols[1, 1]):.4f}  GLS {np.sqrt(V_gls[1, 1]):.4f}")
print(f"quasi-demeaning constant c = {cq:.4f},  s_GLS = {np.sqrt(rss_gls / (n - p)):.4f}")
```

::: {when-format="html"}
![**Figure 31.1.1.** El Niño sea surface temperatures, 1950–2010. (a) The annual means
of the ordinary least squares residuals, with the \( \pm2 \) standard deviation lines
that independent errors would give (dashed). Whole years are warm or cold together.
(b) The estimated trend component with nominal 95% bands, from ordinary least squares
with the usual covariance and from generalized least squares with
\( \rho=0.5 \).](elnino_gls.svg){#fig-ggm-elnino width=100%}
:::

::: {when-format="pdf"}
![El Niño sea surface temperatures, 1950–2010. (a) The annual means
of the ordinary least squares residuals, with the \( \pm2 \) standard deviation lines
that independent errors would give (dashed). Whole years are warm or cold together.
(b) The estimated trend component with nominal 95% bands, from ordinary least squares
with the usual covariance and from generalized least squares with
\( \rho=0.5 \).](elnino_gls.pdf){width=100%}
:::

The \( F \) test of @thm-ggm-inference(e) is computed exactly as in Part III: fit
the full and the reduced model in the same geometry and compare their
\( \V^{-1} \) residual sums of squares. The listing does this for both hypotheses of
@exm-ggm-elnino, using \( \rho=0 \) to recover the ordinary analysis.

```{.python .run #cell-elnino-gls-ftest}
def f_test(X, y, drop, rho, m):
    """F test that the coefficients listed in `drop` are zero, in the rho-geometry."""
    keep = [j for j in range(X.shape[1]) if j not in drop]
    full, null = gls_fit(X, y, rho, m), gls_fit(X[:, keep], y, rho, m)
    df1, df2 = len(drop), len(y) - X.shape[1]
    F = (null[2] - full[2]) / df1 / (full[2] / df2)
    return F, stats.f.sf(F, df1, df2)


F_trend_gls, p_trend_gls = f_test(X, y, [1], rho, m)
F_trend_ols, p_trend_ols = f_test(X, y, [1], 0.0, m)   # rho = 0 is ordinary least squares
F_harm_gls, p_harm_gls = f_test(X, y, [4, 5, 6, 7], rho, m)
F_harm_ols, p_harm_ols = f_test(X, y, [4, 5, 6, 7], 0.0, m)
print(f"trend:      F = {F_trend_ols:7.2f} (rho = 0)   F = {F_trend_gls:6.2f} (rho = 0.5)")
print(f"harmonics:  F = {F_harm_ols:7.2f} (rho = 0)   F = {F_harm_gls:6.2f} (rho = 0.5)")
```

::: {.idea}
Generalized least squares is one change of coordinates away from ordinary least
squares, and everything in Parts II and III transfers — provided *all* of it
does: the estimate, the residual sum of squares, the degrees of freedom and the
sums of squares of every hypothesis. Mixing geometries, as in a generalized
estimate with an ordinary standard error, is the one way to get this wrong.
:::

## Exercises

### A. Check your understanding

::: {#exr-ggm-scale}
[A1]

Show that replacing \( (\sigma^2,\V) \) by \( (\sigma^2/c,\,c\V) \) for \( c>0 \) leaves
\( \hbeta_{\V} \), \( \A_{\V} \) and every statistic in @thm-ggm-inference(d), (e)
unchanged, and multiplies \( \text{SSE}_{\V} \) and \( s_{\V}^2 \) by \( 1/c \). Conclude
that no normalization of \( \V \) affects any inference about \( \bbeta \).
:::

::: {.solution}
The criterion \( Q_{c\V}(\bb)=Q_{\V}(\bb)/c \) is minimized at the same \( \bb \), and
\( (c\V)^{-1}=\V^{-1}/c \) cancels from @eq-ggm-gls-projection. Then
\( \text{SSE}_{c\V}=\text{SSE}_{\V}/c \) and
\( \blambda\T(\X\T(c\V)^{-1}\X)\ginv\blambda=c\,\blambda\T(\X\T\V^{-1}\X)\ginv\blambda \), so the
factors \( c \) cancel in the \( t \) statistic; in \( F \) they cancel between numerator
and denominator.
:::

::: {#exr-ggm-residual-sum}
[A2]

Show that \( \bone\T\he_{\V}=0 \) for every \( \y \) iff \( \bone\in\C(\V^{-1}\X) \). Give an
example with an intercept in which this fails.
:::

### B. Practice

::: {#exr-ggm-quasi}
[B1]

Let \( \V_0=(1-\rho)\I_m+\rho\bone\bone\T \) with \( -1/(m-1)<\rho<1 \). Verify
@eq-ggm-equicorr-inverse by @thm-mat-woodbury, and show that
\( \bL_0=(\I_m-c\bone\bone\T/m)/\sqrt{1-\rho} \) satisfies
\( \bL_0\V_0\bL_0\T=\I_m \) for \( c=1-\sqrt{(1-\rho)/(1+(m-1)\rho)} \). Deduce @eq-ggm-quasi.
:::

::: {.solution}
Write \( \bP=\bone\bone\T/m \), which is idempotent, so \( \V_0=(1-\rho)\I+\rho m\bP \).
For any \( c \),
\( (\I-c\bP)\V_0(\I-c\bP)=(1-\rho)\I+\bigl[(1-\rho)(c^2-2c)+\rho m(1-c)^2\bigr]\bP \),
using \( \bP^2=\bP \). Setting \( a=1-\rho \) and \( b=\rho m \), the bracket vanishes iff
\( (a+b)c^2-2(a+b)c+b=0 \), that is \( c^2-2c+b/(a+b)=0 \), whose root in \( [0,1) \) is
\( c=1-\sqrt{a/(a+b)} \). Dividing by \( 1-\rho \) gives \( \bL_0\V_0\bL_0\T=\I \). The
Woodbury identity with \( \bU=\bone \) gives @eq-ggm-equicorr-inverse directly.
Applying \( \bL_0 \) block by block is @eq-ggm-quasi, since
\( \bone\bone\T\y_k/m=\bar y_k\bone \).
:::

::: {#exr-ggm-oblique-leverage}
[B2]

Show that \( \tr\A_{\V}=r \) for every positive definite \( \V \). Then take \( n=2 \),
\( p=1 \),
\[
\X=\begin{pmatrix}1\\2\end{pmatrix},\qquad
\V=\begin{pmatrix}1&0.9\\0.9&1\end{pmatrix},
\]
and show that the diagonal of \( \A_{\V} \) is \( (-4/7,\,11/7) \). Why does this not
contradict @prp-proj-leverage?
:::

::: {#exr-ggm-gls-vs-ols-rss}
[B3]

Show that \( \he_{\V}\T\he_{\V}\ge\he\T\he \), where \( \he=(\I-\M)\y \) is the ordinary
residual, with equality iff \( \X\hbeta_{\V}=\M\y \). Show also that
\( \he\T\V^{-1}\he\ge\he_{\V}\T\V^{-1}\he_{\V} \). Explain why neither inequality says
anything about which fit is better.
:::

::: {#exr-ggm-two-groups}
[B4]

Observations \( Y_{kj} \), \( j=1,\dots,m_k \), \( k=1,\dots,G \), have mean \( \mu \) and
covariance \( \sigma^2\V \) with \( \V=(1-\rho)\I+\rho\Z\Z\T \), \( \Z \) the cluster
indicator matrix. Using @eq-ggm-equicorr-inverse cluster by cluster, show that
\[
\hat\mu_{\V}=\frac{\sum_k w_k\bar Y_k}{\sum_k w_k},\qquad
w_k=\frac{m_k}{1-\rho+m_k\rho},
\]{#eq-ggm-cluster-weights}

with \( \Var(\hat\mu_{\V})=\sigma^2/\sum_kw_k \). Interpret the weights as
\( \rho\to0 \) and \( \rho\to1 \).
:::

::: {.solution}
With \( \V \) block diagonal, \( \V^{-1} \) is block diagonal with blocks
\( \bigl(\I-\frac{\rho}{1-\rho+m_k\rho}\bone\bone\T\bigr)/(1-\rho) \). Hence
\( \bone\T\V^{-1}\bone=\sum_k\bigl(m_k-\rho m_k^2/(1-\rho+m_k\rho)\bigr)/(1-\rho)=\sum_kw_k \),
and \( \bone\T\V^{-1}\y=\sum_kw_k\bar y_k \) by the same computation. The generalized
normal equation is \( (\bone\T\V^{-1}\bone)\mu=\bone\T\V^{-1}\y \). As \( \rho\to0 \),
\( w_k\to m_k \) and \( \hat\mu_{\V}\to\bar y \); as \( \rho\to1 \), \( w_k\to1 \) and every
cluster counts once, however large, because the extra observations only repeat
the same shock.
:::

### C. Going deeper

::: {#exr-ggm-blup}
[C1]

Suppose \( (\Y\T,Y_0)\T \) is jointly distributed with \( \E(\Y)=\X\bbeta \),
\( \E(Y_0)=\x_0\T\bbeta \), \( \Cov(\Y)=\sigma^2\V \), \( \Cov(\Y,Y_0)=\sigma^2\mathbf{c} \) and
\( \Var(Y_0)=\sigma^2v_0 \), with \( \V \) positive definite and \( \rank(\X)=p \). Among
predictors \( \mathbf{a}\T\Y \) with \( \E(\mathbf{a}\T\Y-Y_0)=0 \) for all \( \bbeta \), show that the
one minimizing \( \E(\mathbf{a}\T\Y-Y_0)^2 \) is
\[
\hat Y_0=\x_0\T\hbeta_{\V}+\mathbf{c}\T\V^{-1}(\Y-\X\hbeta_{\V}),
\]
with prediction error variance
\( \sigma^2\bigl[v_0-\mathbf{c}\T\V^{-1}\mathbf{c}+\mathbf{d}\T(\X\T\V^{-1}\X)^{-1}\mathbf{d}\bigr] \),
\( \mathbf{d}=\x_0-\X\T\V^{-1}\mathbf{c} \). This is Goldberger's formula;
[Chapter 32](../ch32-linear-mixed-models/index.html) recognizes it as the best linear
unbiased predictor of a random effect.
:::

::: {.solution}
Unbiasedness requires \( \X\T\mathbf{a}=\x_0 \), and
\( \E(\mathbf{a}\T\Y-Y_0)^2=\sigma^2(\mathbf{a}\T\V\mathbf{a}-2\mathbf{a}\T\mathbf{c}+v_0) \). By @prp-mat-lagrange
with the positive definite \( \V \), the minimizer is
\( \mathbf{a}=\V^{-1}\mathbf{c}+\V^{-1}\X(\X\T\V^{-1}\X)^{-1}\mathbf{d} \), and collecting terms gives
\( \mathbf{a}\T\Y=\x_0\T\hbeta_{\V}+\mathbf{c}\T\V^{-1}(\Y-\X\hbeta_{\V}) \); the minimum value is
the stated variance. Setting \( \mathbf{c}=\bzero \) recovers the ordinary prediction of
[Chapter 12](../ch12-intervals-and-bands/index.html). A correlated future observation is
predicted by correcting the fitted mean with the part of the current residual
that predicts it.
:::

::: {#exr-ggm-mle}
[C2]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V(\boldsymbol{\uptheta})) \) with \( \rank(\X)=p \). Show that for
each fixed \( \boldsymbol{\uptheta} \) the log-likelihood is maximized over \( (\bbeta,\sigma^2) \) at
\( \hbeta_{\V(\boldsymbol{\uptheta})} \) and
\( \hat\sigma^2(\boldsymbol{\uptheta})=\text{SSE}_{\V(\boldsymbol{\uptheta})}/n \), and that the profile
log-likelihood is
\[
\ell_p(\boldsymbol{\uptheta})=-\tfrac{n}{2}\log\hat\sigma^2(\boldsymbol{\uptheta})-\tfrac12\log\det\V(\boldsymbol{\uptheta})+\text{const}.
\]
Why does the determinant term make the normalization of \( \V \) matter here, when
@exr-ggm-scale says it does not matter for inference about \( \bbeta \)?
:::
