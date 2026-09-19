# The coefficient of determination and its limits

The coefficient of determination is the most widely reported summary of a regression, and
probably the most widely misread. It is a ratio of two rows of the corrected table,
\[
R^2=\frac{\text{SSR}}{\text{SST}}=1-\frac{\text{SSE}}{\text{SST}},
\]
defined when \( \bone\in\C(\X) \) and \( \y \) is not constant. [Section 6.7](../ch06-projections/07-reparameterization.html)
showed that it is a squared cosine (@thm-proj-r2-cosine) that cannot decrease as the model space
grows. This section adds \( R^2 \) as a maximal correlation, its null distribution, the adjusted
and partial versions, and what it does not measure.

## A maximal correlation

The square root \( R=\sqrt{R^2} \) is called the **multiple correlation coefficient**. The name
comes from the following property, the sample counterpart of the population result @prp-rv-multiple-correlation.

::: {#thm-ss-r2-max-correlation}
[\( R \) is the largest correlation with a linear combination]

Suppose \( \bone\in\C(\X) \) and \( \y \) is not constant. For every \( \bb \) such that \( \X\bb \) is not
constant,
\[
\operatorname{corr}(\y,\X\bb)^2\le R^2.
\]
If \( R^2>0 \), equality holds iff \( \X\bb=a\hY+c\bone \) for some \( a\neq0 \) and some \( c \), and the
correlation \( \operatorname{corr}(\y,\X\bb) \) is maximized, with value \( R \), by the fitted values
themselves. If \( R^2=0 \), every such correlation is zero.
:::

::: {.proof}
Let \( \bu=\y-\bar{y}\bone \) and \( \mathcal V=\C(\X)\cap\bone\perpc \), whose projection is
\( \M-\bP_1 \) (@thm-proj-nested), and let \( \bv=\X\bb-\overline{\X\bb}\,\bone \). Since
\( \bone\in\C(\X) \), \( \bv\in\mathcal V \), so \( \bu\T\bv=\bu\T(\M-\bP_1)\bv=((\M-\bP_1)\bu)\T\bv \). The
sample correlation is \( \bu\T\bv/(\norm{\bu}\norm{\bv}) \), and by Cauchy–Schwarz (@prp-mat-cauchy-schwarz)
\[
\operatorname{corr}(\y,\X\bb)^2=\frac{\bigl(((\M-\bP_1)\bu)\T\bv\bigr)^2}{\norm{\bu}^2\norm{\bv}^2}
\le\frac{\norm{(\M-\bP_1)\bu}^2}{\norm{\bu}^2}=\frac{\text{SSR}}{\text{SST}},
\]
with equality iff \( \bv \) and \( (\M-\bP_1)\bu=\hY-\bar{y}\bone \) are linearly dependent. If
\( R^2>0 \) the second vector is nonzero, and since \( \bv\neq\bzero \) dependence means
\( \bv=a(\hY-\bar{y}\bone) \) with \( a\neq0 \), which is the stated condition; for \( a>0 \) the
correlation itself is \( +R \). If \( R^2=0 \), then \( (\M-\bP_1)\bu=\bzero \) and \( \bu\T\bv=0 \).
:::

Among all linear combinations of the regressors, the fitted values are the one most correlated
with the response. When the regressors are pairwise uncorrelated, \( R^2 \) is the sum of the squared simple correlations
(@exr-ss-uncorrelated-regressors). In general it can be smaller or larger than that sum
(@exr-ss-suppression-r2).

## \( R^2 \) when nothing is going on

Since \( R^2 \) never decreases as regressors are added, a model with many regressors has a sizeable
\( R^2 \) even if none of them matters. The theory of [Chapter 4](../ch04-quadratic-forms/index.html)
says how sizeable.

::: {#thm-ss-r2-null}
[The distribution of \( R^2 \)]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \) with \( \bmu\in\C(\X) \), \( \bone\in\C(\X) \) and
\( 1<r=\rank(\X)<n \). Let
\( F=[\text{SSR}/(r-1)]/[\text{SSE}/(n-r)] \) and \( \gamma=\norm{(\M-\bP_1)\bmu}^2/\sigma^2 \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( R^2=\dfrac{(r-1)F}{(r-1)F+n-r} \), and \( F\sim F(r-1,n-r,\gamma) \);

2. if \( \bmu\in\spn(\bone) \), so that no regressor matters, then
   \( R^2\sim\mathrm{Beta}\bigl(\tfrac{r-1}{2},\tfrac{n-r}{2}\bigr) \) and
   \( \E(R^2)=\dfrac{r-1}{n-1} \).
:::

:::

::: {.proof}
(a) \( R^2/(1-R^2)=\text{SSR}/\text{SSE}=(r-1)F/(n-r) \), and solving for \( R^2 \) gives the formula.
The projections \( \M-\bP_1 \) and \( \I-\M \) are orthogonal with ranks \( r-1 \) and \( n-r \), and
\( (\I-\M)\bmu=\bzero \), so by @thm-qf-orthogonal-projections \( \text{SSR}/\sigma^2\sim\chi^2(r-1,\gamma) \)
and \( \text{SSE}/\sigma^2\sim\chi^2(n-r) \), independently. The ratio \( F \) is therefore
\( F(r-1,n-r,\gamma) \) (@def-qf-noncentral-f).

(b) Now \( \gamma=0 \), and \( R^2=U/(U+V) \) with \( U=\text{SSR}/\sigma^2\sim\chi^2(r-1) \) and
\( V=\text{SSE}/\sigma^2\sim\chi^2(n-r) \) independent. By @exr-qf-beta this is
\( \mathrm{Beta}\bigl(\tfrac{r-1}{2},\tfrac{n-r}{2}\bigr) \), whose mean is
\( \tfrac{r-1}{2}\big/\tfrac{n-1}{2}=(r-1)/(n-1) \).
:::

Part (b) answers the exercise @exr-proj-r2-null-distribution of [Chapter 6](../ch06-projections/index.html). With \( k=r-1 \) regressors
of pure noise, \( R^2 \) is on average \( k/(n-1) \). A model with \( 10 \) useless regressors fitted to
\( 30 \) observations has expected \( R^2 \) about one third. [Figure 9.5.1](#fig-ss-r2-null) compares the
Beta densities with \( 100{,}000 \) simulated data sets with \( n=30 \) and \( k=1,5,10 \)
regressors drawn at random. The simulated means are \( 0.0345 \), \( 0.1725 \) and
\( 0.3445 \), against the exact \( 0.0345 \), \( 0.1724 \) and
\( 0.3448 \). The upper \( 5\% \) points of the three null distributions are
\( 0.130 \), \( 0.353 \) and \( 0.556 \), so with ten regressors an \( R^2 \) of
one half is unremarkable. Part (a) shows that \( R^2 \) and \( F \) carry the same information:
testing \( R^2 \) against its null distribution is the overall \( F \) test of the regression.

::: {when-format="html"}
![**Figure 9.5.1.** Null distribution of \( R^2 \) for \( n=30 \) observations and \( 1 \), \( 5 \) or \( 10 \)
regressors unrelated to the response. Histograms: \( 100{,}000 \) simulated data sets each.
Curves: the Beta densities of @thm-ss-r2-null.](r2_null.svg){#fig-ss-r2-null width=78%}
:::

::: {when-format="pdf"}
![Null distribution of \( R^2 \) for \( n=30 \) observations and \( 1 \), \( 5 \) or \( 10 \)
regressors unrelated to the response. Histograms: \( 100{,}000 \) simulated data sets each.
Curves: the Beta densities of @thm-ss-r2-null.](r2_null.pdf){width=78%}
:::

```{.python .run #cell-r2-null-simulate}
import numpy as np
rng = np.random.default_rng(20260919)
n, reps = 30, 100_000
results = {}
for k in (1, 5, 10):                                   # number of regressors besides 1
    X = np.column_stack([np.ones(n), rng.normal(size=(n, k))])
    r = k + 1
    Q, _ = np.linalg.qr(X)
    Y = 3.0 + rng.normal(size=(n, reps))               # no regressor matters
    Z = Q.T @ Y
    ssr = np.sum(Z[1:] ** 2, axis=0)
    sst = np.sum(Y ** 2, axis=0) - Z[0] ** 2
    R2 = ssr / sst
    R2adj = 1 - (1 - R2) * (n - 1) / (n - r)
    print(f"k = {k:2d}: mean R^2 = {R2.mean():.4f} (theory {k / (n - 1):.4f}),"
          f" mean adjusted R^2 = {R2adj.mean():+.4f}")
    results[k] = R2
```

## Adjusted \( R^2 \)

The bias of \( R^2 \) under the null suggests a correction. Write \( s^2=\text{SSE}/(n-r) \) and
\( s_y^2=\text{SST}/(n-1) \) for the unbiased estimators of the error variance in the fitted model and
in the intercept-only model. Replacing sums of squares by these mean squares defines the
**adjusted coefficient of determination**
\[
\bar{R}^2=1-\frac{s^2}{s_y^2}=1-(1-R^2)\frac{n-1}{n-r}.
\]{#eq-ss-adjusted-r2}

::: {#prp-ss-adjusted-r2}
[Properties of the adjusted \( R^2 \)]

Let \( \bone\in\C(\X) \), \( 1\le r=\rank(\X)<n \), and let \( F \) be the overall \( F \) ratio of @thm-ss-r2-null.

::: {.enumerate options="label=(\alph*)"}
1. \( \bar{R}^2\le R^2 \), with equality iff \( r=1 \) or \( R^2=1 \).

2. For \( r>1 \), \( \bar{R}^2<0 \) iff \( F<1 \).

3. If a model with rank \( r_0 \) is enlarged to one of rank \( r=r_0+q \) containing it, then
   \( \bar{R}^2 \) increases iff the \( F \) ratio of the added terms,
   \( [\text{SS}(\text{added}\mid\text{old})/q]/s^2 \), exceeds \( 1 \). For a single added column
   this means \( |t|>1 \).

4. Among models for the same response, \( \bar{R}^2 \) is largest for the model with the smallest
   residual mean square \( s^2 \).

5. Under the normal model with \( \bmu\in\spn(\bone) \), \( \E(\bar{R}^2)=0 \).
:::

:::

::: {.proof}
(a) \( 1-\bar{R}^2=(1-R^2)(n-1)/(n-r) \), and \( (n-1)/(n-r)\ge1 \) with equality iff \( r=1 \), so
\( 1-\bar{R}^2\ge1-R^2 \), with equality iff \( r=1 \) or \( 1-R^2=0 \).

(b) \( \bar{R}^2<0 \) iff \( s^2>s_y^2 \), that is,
\( \text{SSE}(n-1)>(\text{SSR}+\text{SSE})(n-r) \), that is, \( \text{SSE}(r-1)>\text{SSR}(n-r) \), which
for \( r>1 \) is \( F<1 \).

(c) \( s_y^2 \) does not depend on the model, so \( \bar{R}^2 \) increases iff \( s^2 \) decreases. Write
\( \text{SSE}_0=\text{SSE}+D \) with \( D \) the extra sum of squares. Then
\( \text{SSE}/(n-r)<(\text{SSE}+D)/(n-r_0) \) iff \( \text{SSE}\,(n-r_0)<(\text{SSE}+D)(n-r) \), iff
\( q\,\text{SSE}<D(n-r) \), iff \( (D/q)/s^2>1 \). For one column the ratio is \( t^2 \)
(@cor-ss-single-column).

(d) is immediate from @eq-ss-adjusted-r2.

(e) By @thm-ss-r2-null(b), \( \E(1-R^2)=1-(r-1)/(n-1)=(n-r)/(n-1) \), so
\( \E(1-\bar{R}^2)=1 \).
:::

Part (e) is the point of the adjustment. In the simulation the averages of \( \bar{R}^2 \) are
\( 0.0000 \), \( 0.0001 \) and \( -0.0005 \), but \( \bar{R}^2 \) is negative in a
proportion \( 0.674 \) of the data sets with one regressor: it is not a proportion of
anything. Part (c) shows what selection by \( \bar{R}^2 \) amounts to: a regressor is admitted when
\( |t|>1 \), far below the threshold near two of a \( 5\% \) test, so the chosen models tend to be too
large (Chapter 29).

::: {#exm-ss-crime-adjusted}
[Adjusted \( R^2 \) for the state data]

For all seven models built from poverty, single parenthood and urbanization (each with an
intercept):

| Regressors | \( R^2 \) | \( \bar{R}^2 \) |
|---|---|---|
| poverty | 0.3599 | 0.3466 |
| single | 0.5719 | 0.5629 |
| urban | 0.0221 | 0.0017 |
| poverty, single | 0.6410 | 0.6257 |
| poverty, urban | 0.4626 | 0.4398 |
| single, urban | 0.5844 | 0.5668 |
| all three | 0.6419 | 0.6186 |

\( R^2 \) is largest for the full model, as it must be. \( \bar{R}^2 \) prefers poverty and single
parenthood without urbanization. Adding urbanization to that model raises \( R^2 \) slightly but
lowers \( \bar{R}^2 \), because its \( t \) statistic in the full model is only \( 0.349 \), so its
\( F \) ratio is \( 0.121<1 \).
:::

## Partial \( R^2 \)

Applied to what a reduced model leaves unexplained, the idea of \( R^2 \) gives a measure for a term.

::: {#def-ss-partial-r2}
[Partial coefficient of determination]

Let \( \C(\X_1)\subseteq\C(\X) \) with \( \bone\in\C(\X_1) \). The **partial coefficient of determination**
of \( \X \) over \( \X_1 \) is
\[
R^2_{\X\mid\X_1}=\frac{\text{SSE}(\X_1)-\text{SSE}(\X)}{\text{SSE}(\X_1)}
=\frac{\text{SS}(\X\mid\X_1)}{\text{SSE}(\X_1)} ,
\]
the proportion of the residual variation of the reduced model that the larger model explains.
:::

::: {#thm-ss-partial-r2}
[Partial \( R^2 \), \( F \) and \( t \)]

With the notation of @def-ss-partial-r2, let \( r_1=\rank(\X_1) \), \( r=\rank(\X) \), \( q=r-r_1>0 \),
\( s^2=\text{SSE}(\X)/(n-r) \), and let \( F=[\text{SS}(\X\mid\X_1)/q]/s^2 \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( 1-R^2=(1-R_1^2)(1-R^2_{\X\mid\X_1}) \), where \( R^2 \) and \( R_1^2 \) belong to the two models;

2. \( R^2_{\X\mid\X_1}=\dfrac{qF}{qF+n-r} \);

3. if \( \X=[\X_1,\x_j] \) adds one column, then \( R^2_{\X\mid\X_1}=\dfrac{t_j^2}{t_j^2+n-r} \), and it
   equals the squared sample correlation between the residuals \( (\I-\M_1)\y \) and \( (\I-\M_1)\x_j \);

4. for a chain \( \C(\bone)=\C(\X_0)\subseteq\C(\X_1)\subseteq\dots\subseteq\C(\X_k) \),
   \( 1-R^2_k=\prod_{j=1}^k\bigl(1-R^2_{\X_j\mid\X_{j-1}}\bigr) \).
:::

:::

::: {.proof}
(a) \( 1-R^2=\text{SSE}(\X)/\text{SST}=\bigl(\text{SSE}(\X_1)/\text{SST}\bigr)\bigl(\text{SSE}(\X)/\text{SSE}(\X_1)\bigr) \),
and the two factors are \( 1-R_1^2 \) and \( 1-R^2_{\X\mid\X_1} \).

(b) Write \( D=\text{SS}(\X\mid\X_1)=qFs^2 \) and \( \text{SSE}(\X)=(n-r)s^2 \). Then
\( R^2_{\X\mid\X_1}=D/(D+\text{SSE}(\X))=qF/(qF+n-r) \).

(c) With \( q=1 \), \( F=t_j^2 \) by @cor-ss-single-column, and (b) gives the first formula. For the
correlation, let \( \tilde{\y}=(\I-\M_1)\y \) and \( \tilde{\x}=(\I-\M_1)\x_j \). Then
\( \text{SSE}(\X_1)=\norm{\tilde{\y}}^2 \), and by @lem-ss-residualized
\( \text{SS}(\X\mid\X_1)=(\tilde{\x}\T\tilde{\y})^2/\norm{\tilde{\x}}^2 \). Their ratio is the squared
cosine of the angle between \( \tilde{\y} \) and \( \tilde{\x} \). Both have mean zero because
\( \bone\in\C(\X_1) \), so the cosine is their sample correlation (@exr-proj-partial-correlation).

(d) Apply (a) repeatedly along the chain, with \( R_0^2=0 \) for the intercept-only model.
:::

Part (c) identifies the partial \( R^2 \) of one regressor with a squared **sample partial
correlation**, the sample counterpart of @def-mvn-partial-correlation; Chapter 14 studies its sampling distribution. Within one
model, ranking regressors by \( |t| \) and by partial \( R^2 \) is the same thing.

::: {#exm-ss-crime-partial}
[Partial \( R^2 \) in the state data]

In the three-regressor model, with \( n-r=46 \):

| Regressor | \( t \) | partial \( R^2 \) | \( R^2 \) without it |
|---|---|---|---|
| poverty | 2.718 | 0.1383 | 0.5844 |
| single | 4.799 | 0.3337 | 0.4626 |
| urban | 0.349 | 0.0026 | 0.6410 |

Each partial \( R^2 \) equals \( t^2/(t^2+46) \), and each row satisfies part (a): for poverty,
\( 1-0.6419=(1-0.5844)(1-0.1383) \). Poverty explains \( 14\% \) of the variation that single parenthood
and urbanization leave unexplained, although on its own it explains \( 36\% \) of the total
(\( R^2=0.3599 \)). For the pair poverty and urbanization over single parenthood, the \( F \) ratio
on \( 2 \) and \( 46 \) degrees of freedom is \( 4.502 \) and the partial \( R^2 \) is
\( 0.1637 \).
:::

```{.python .run #cell-partial-r2-partial}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")

y = data["murder"].to_numpy()
names = ["poverty", "single", "urban"]
n = len(y)

def fit(cols):
    X = np.column_stack([np.ones(n)] + [data[c].to_numpy() for c in cols])
    return sm.OLS(y, X).fit()

full = fit(names)
df_res = full.df_resid
for j, name in enumerate(names, start=1):
    reduced = fit([c for c in names if c != name])
    partial_r2 = (reduced.ssr - full.ssr) / reduced.ssr
    t = full.tvalues[j]
    print(f"{name:8s} partial R^2 = {partial_r2:.4f},  t = {t:6.3f},"
          f"  t^2/(t^2 + df) = {t ** 2 / (t ** 2 + df_res):.4f}")
print(f"R^2 = {full.rsquared:.4f}, adjusted R^2 = {full.rsquared_adj:.4f}")
```

## What \( R^2 \) does not measure

A high \( R^2 \) is often read as “a good model” and a low one as “a useless model”; neither reading is safe.

**It is not a property of the relationship alone.** \( R^2 \) depends on the spread of the
regressors, which is a property of the design. In simple regression with fixed \( x \) values,
the expectations of the two sums of squares are (@exr-ss-design-r2)
\[
\E(\text{SSR})=\sigma^2+\beta_1^2S_{xx},\qquad \E(\text{SST})=(n-1)\sigma^2+\beta_1^2S_{xx},
\]
so \( R^2 \) is roughly \( \beta_1^2S_{xx}/\bigl(\beta_1^2S_{xx}+(n-1)\sigma^2\bigr) \). In [Figure 9.5.2](#fig-ss-r2-range) the line
\( \E(y)=2+0.5x \) and the same \( 50 \) errors with \( \sigma=1 \) are combined with \( x \) values spread
evenly over \( [4,6] \) and over \( [0,10] \). The two designs span the same \( \C(\X) \), so the
residuals are identical and the residual mean square is \( 1.101 \) in both fits, but \( R^2 \) is \( 0.026 \) in one case and
\( 0.632 \) in the other. The approximate values from the formula are
\( 0.081 \) and \( 0.689 \). Restricting the range of a regressor lowers \( R^2 \) without changing the relationship.

::: {when-format="html"}
![**Figure 9.5.2.** One line, two designs. The dashed line is the true mean \( 2+0.5x \), the solid line the
least squares fit. The errors and the residual mean square are the same in both panels.](r2_range.svg){#fig-ss-r2-range width=100%}
:::

::: {when-format="pdf"}
![One line, two designs. The dashed line is the true mean \( 2+0.5x \), the solid line the
least squares fit. The errors and the residual mean square are the same in both panels.](r2_range.pdf){width=100%}
:::

```{.python .run #cell-r2-range-range}
import numpy as np
rng = np.random.default_rng(9)
n, beta0, beta1, sigma = 50, 2.0, 0.5, 1.0
eps = sigma * rng.normal(size=n)                        # the same errors in both designs

def fit(x):
    y = beta0 + beta1 * x + eps
    X = np.column_stack([np.ones(n), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    sse = np.sum((y - X @ b) ** 2)
    r2 = 1 - sse / np.sum((y - y.mean()) ** 2)
    return y, b, sse / (n - 2), r2

for label, x in [("narrow", np.linspace(4, 6, n)), ("wide", np.linspace(0, 10, n))]:
    y, b, s2, r2 = fit(x)
    Sxx = np.sum((x - x.mean()) ** 2)
    expected = beta1 ** 2 * Sxx / (beta1 ** 2 * Sxx + (n - 1) * sigma ** 2)
    print(f"{label:6s}: slope {b[1]:.3f}, s^2 {s2:.3f}, R^2 {r2:.3f}, design value {expected:.3f}")
```

**It is not a test of the model.** A line fitted to points on a smooth curve can have a high
\( R^2 \) while missing the curvature; residual plots and lack-of-fit tests (Chapter 14) reveal it.
Conversely, a correct model has a low \( R^2 \) whenever the error variance is large, and its
coefficients can still be precisely estimated and important.

**It is not comparable across responses.** Models for \( y \) and \( \log y \), or for individuals and
for group averages, have different denominators \( \text{SST} \). Averaging replicates inflates \( R^2 \)
(@exr-ss-averaging).

**It is not comparable without an intercept** (@prp-ss-no-intercept, @exm-ss-no-intercept).

**It is not a measure of prediction or of cause.** Computed on the fitting data, it overstates
predictive accuracy for new data (Chapter 29), and variation “explained” is linear association,
not the effect of changing a regressor (Chapter 25).

## Exercises

### A. Check your understanding

::: {#exr-ss-r2-f-numbers}
[A1]

A regression with an intercept, \( n=20 \) and \( r=6 \) has \( R^2=0.40 \). Compute \( \bar{R}^2 \) and the
overall \( F \) ratio. Is \( \bar{R}^2 \) positive? Is it consistent with part (b) of @prp-ss-adjusted-r2?
:::

### B. Practice

::: {#exr-ss-design-r2}
[B1]

In simple regression with fixed \( x \) values, \( \E(\Y)=\beta_0\bone+\beta_1\x \) and
\( \Cov(\Y)=\sigma^2\I \), use @thm-rv-quadform-mean to show that
\( \E(\text{SSR})=\sigma^2+\beta_1^2S_{xx} \) and \( \E(\text{SST})=(n-1)\sigma^2+\beta_1^2S_{xx} \).
Show that the ratio of these expectations tends to \( 1 \) as \( S_{xx}\to\infty \) with \( n \) fixed,
and to \( 1/(n-1) \) as \( S_{xx}\to0 \), and interpret the second limit.
:::

::: {.solution}
\( \text{SSR}=\y\T(\M-\bP_1)\y \), where \( \M-\bP_1 \) projects onto the line spanned by
\( \x-\bar{x}\bone \) (@exm-proj-simple-nested) and has trace \( 1 \). With
\( \bmu=\beta_0\bone+\beta_1\x \), \( (\M-\bP_1)\bmu=\beta_1(\x-\bar{x}\bone) \), whose squared length is
\( \beta_1^2S_{xx} \). @thm-rv-quadform-mean gives \( \E(\text{SSR})=\sigma^2+\beta_1^2S_{xx} \).
Similarly \( \text{SST}=\y\T(\I-\bP_1)\y \), with trace \( n-1 \) and
\( (\I-\bP_1)\bmu=\beta_1(\x-\bar{x}\bone) \). The ratio
\( (\sigma^2+\beta_1^2S_{xx})/((n-1)\sigma^2+\beta_1^2S_{xx}) \) tends to \( 1 \) as \( S_{xx}\to\infty \)
and to \( 1/(n-1) \) as \( S_{xx}\to0 \). The second limit is the null value of \( \E(R^2) \) for \( r=2 \)
(@thm-ss-r2-null): with almost no spread in \( x \), the slope contributes almost nothing and the
regression sum of squares is almost pure noise.
:::

::: {#exr-ss-averaging}
[B2]

Suppose each of \( g \) distinct values \( x_1,\dots,x_g \) carries \( m \) observations. Show that the
least squares line fitted to all \( gm \) observations is the same as the line fitted to the \( g \)
group means, and that its \( R^2 \) for the individual data is at most its \( R^2 \) for the means, with
equality iff every observation equals its group mean.
:::

::: {.solution}
The model space for the individual data is spanned by \( \bone \) and \( \x \), both constant within
groups, so it lies inside the space of group-constant vectors. By @thm-proj-nested, projecting
\( \y \) onto it equals projecting the vector of group means (repeated \( m \) times) onto it, which is
the fit to the means with each point weighted equally. Write \( \text{SSR}_g \) and \( \text{SST}_g \) for
the sums over the \( g \) means, and \( W \) for the within-group sum of squares. For the individual
data \( \text{SSR}=m\,\text{SSR}_g \) and \( \text{SST}=m\,\text{SST}_g+W \). Hence
\( R^2=m\,\text{SSR}_g/(m\,\text{SST}_g+W)\le\text{SSR}_g/\text{SST}_g \), with equality iff \( W=0 \).
:::

::: {#exr-ss-chain-product}
[B3]

Using part (d) of @thm-ss-partial-r2 and the sequential sums of squares of
@exm-proj-sequential-order (order A), express \( 1-R^2 \) for the three-regressor model as a product
of three factors, and identify each factor as one minus a partial \( R^2 \).
:::

### C. Going deeper

::: {#exr-ss-r2-monotone-power}
[C1]

Under the conditions of @thm-ss-r2-null, show that for every \( c\in(0,1) \) the probability
\( \Pr(R^2>c) \) is a strictly increasing function of the noncentrality \( \gamma \). *Hint:* part (a)
of the theorem and @thm-qf-f-power.
:::

::: {#exr-ss-suppression-r2}
[C2]

For two regressors with sample correlations \( r_{y1}=0.5 \), \( r_{y2}=0 \) and \( r_{12}=0.6 \), compute
\( R^2 \) and compare it with \( r_{y1}^2+r_{y2}^2 \). Use the formulas of @exr-ss-two-regressor-formula.
Explain geometrically why a regressor that is uncorrelated with the response can raise \( R^2 \).
:::

