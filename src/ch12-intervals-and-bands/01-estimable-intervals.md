# Intervals for estimable functions

A point estimate needs a statement of its precision. A confidence interval gives one exactly: a range
computed from the data that contains the unknown quantity with a probability fixed in advance. This
section builds the interval for one estimable function, relates it to the \( t \) test, and records two
facts that are easy to overlook: an interval for a non-estimable function is useless, and the
interval for \( \sigma^2 \), unlike the one for a coefficient, depends delicately on normality.

Throughout the chapter we work in the normal linear model
\[
\Y\sim\Normal_n(\X\bbeta,\sigma^2\I),\qquad \bbeta\in\Real^p,\ \sigma^2>0,
\]
of @eq-opt-normal-model, with \( r=\rank(\X)<n \). We write \( \G \) for any generalized inverse of
\( \X\T\X \), \( \hbeta \) for any least squares estimate, \( \text{SSE}=\norm{(\I-\M)\Y}^2 \) and
\( s^2=\text{SSE}/(n-r) \). For \( 0<\gamma<1 \), \( t_{\nu,\gamma} \), \( F_\gamma(q,\nu) \) and
\( \chi^2_\gamma(\nu) \) denote the upper \( \gamma \) points of \( t(\nu) \), \( F(q,\nu) \) and
\( \chi^2(\nu) \).

## Confidence sets

::: {#def-ci-confidence-set}
[Confidence set]

Let \( \theta=g(\bbeta,\sigma^2) \) be a quantity with values in a set \( \Theta \). A
**confidence set** for \( \theta \) with **confidence level** \( 1-\alpha \) is a rule that assigns to
each data vector \( \y \) a subset \( C(\y)\subseteq\Theta \), such that
\[
\Pr_{\bbeta,\sigma^2}\{\theta\in C(\Y)\}\ \ge\ 1-\alpha
\]
for every \( \bbeta\in\Real^p \) and \( \sigma^2>0 \).
The left side is the **coverage probability** at \( (\bbeta,\sigma^2) \). The set is **exact** if
the coverage equals \( 1-\alpha \) at every parameter value. When \( C(\y) \) is an interval it is a
**confidence interval**.
:::

The probability is taken over \( \Y \) with the parameters fixed: the set is random and the target is not.
A realized interval either contains \( \theta \) or it does not; \( 1-\alpha \) describes the procedure.

The standard way to build an exact set is a **pivot**: a function \( Q(\Y,\theta) \) whose distribution
is the same at every parameter value. If \( \Pr\{Q(\Y,\theta)\in A\}=1-\alpha \) for a fixed set
\( A \), then \( C(\y)=\{\theta_0:Q(\y,\theta_0)\in A\} \) is exact, because the events
\( \{\theta\in C(\Y)\} \) and \( \{Q(\Y,\theta)\in A\} \) are the same event. All the sets of this
chapter are built this way, from pivots supplied by [Chapter 7](../ch07-optimality/index.html).

## The \( t \) interval

::: {#thm-ci-estimable-interval}
[The \( t \) interval for an estimable function]

Let \( \blambda\in\C(\X\T) \), \( \blambda\ne\bzero \), and put
\( \text{se}(\blambda\T\hbeta)=s\sqrt{\blambda\T\G\blambda} \). Let
\[
I(\Y)=\bigl[\,\blambda\T\hbeta-t_{n-r,\alpha/2}\,\text{se}(\blambda\T\hbeta),\ \
\blambda\T\hbeta+t_{n-r,\alpha/2}\,\text{se}(\blambda\T\hbeta)\,\bigr].
\]{#eq-ci-t-interval}

::: {.enumerate options="label=(\alph*)"}
1. \( I(\Y) \) does not depend on the choice of \( \hbeta \) or of \( \G \).

2. \( I(\Y) \) is an exact \( 1-\alpha \) confidence interval for \( \blambda\T\bbeta \).

3. Its length \( L=2t_{n-r,\alpha/2}\,s\sqrt{\blambda\T\G\blambda} \) has expectation
   \( \E L=2t_{n-r,\alpha/2}\,\kappa_{n-r}\,\sigma\sqrt{\blambda\T\G\blambda} \), where
   \[
\kappa_\nu=\E\sqrt{V/\nu}=\sqrt{\frac2\nu}\;\frac{\Gamma\bigl((\nu+1)/2\bigr)}{\Gamma(\nu/2)}<1,\qquad V\sim\chi^2(\nu).
\]

4. The one-sided intervals \( \bigl(-\infty,\ \blambda\T\hbeta+t_{n-r,\alpha}\,\text{se}\bigr] \) and
   \( \bigl[\blambda\T\hbeta-t_{n-r,\alpha}\,\text{se},\ \infty\bigr) \) are exact \( 1-\alpha \) confidence
   intervals for \( \blambda\T\bbeta \).
:::
:::

::: {.proof}
(a) Write \( \blambda=\X\T\boldsymbol{\uprho} \). Then \( \blambda\T\hbeta=\boldsymbol{\uprho}\T\X\hbeta=\boldsymbol{\uprho}\T\M\Y \) and
\( \blambda\T\G\blambda=\boldsymbol{\uprho}\T\X\G\X\T\boldsymbol{\uprho}=\boldsymbol{\uprho}\T\M\boldsymbol{\uprho} \) by @thm-proj-M-formula. Neither
involves the choices, and neither does \( s \).

(b) By @cor-opt-t, \( T=(\blambda\T\hbeta-\blambda\T\bbeta)/\text{se}(\blambda\T\hbeta)\sim t(n-r) \) at every
\( (\bbeta,\sigma^2) \). The event \( \blambda\T\bbeta\in I(\Y) \) is the event \( \lvert T\rvert\le t_{n-r,\alpha/2} \),
whose probability is \( 1-\alpha \). So part (b) restates @eq-opt-t-interval as a confidence statement.

(c) \( L \) is a constant multiple of \( s=\sigma\sqrt{V/(n-r)} \) with \( V=\text{SSE}/\sigma^2\sim\chi^2(n-r) \) (@thm-opt-sampling(b)).
With the \( \chi^2(\nu) \) density of @eq-qf-chisq-density,
\[
\E\sqrt V=\int_0^\infty\frac{v^{1/2}\,v^{\nu/2-1}e^{-v/2}}{2^{\nu/2}\Gamma(\nu/2)}\,dv
=\frac{2^{(\nu+1)/2}\Gamma\bigl((\nu+1)/2\bigr)}{2^{\nu/2}\Gamma(\nu/2)},
\]
because the integrand is a multiple of the \( \chi^2(\nu+1) \) density. Dividing by \( \sqrt\nu \) gives
\( \kappa_\nu \). By Jensen's inequality \( \E\sqrt{V/\nu}<\sqrt{\E(V/\nu)}=1 \), with strict inequality since
\( V \) is not constant.

(d) The events are \( \{T\ge-t_{n-r,\alpha}\} \) and \( \{T\le t_{n-r,\alpha}\} \), each of probability
\( 1-\alpha \).
:::

Part (a) matters when \( \X \) is rank deficient: the interval is a property of the model, not of the
software's choice of generalized inverse. The degrees of freedom are \( n-r \), not \( n-p \).

Part (c) shows what the interval costs. Its expected length is proportional to
\( \sigma\sqrt{\blambda\T\G\blambda} \), the standard deviation of the estimate. Estimating \( \sigma \)
adds the factor \( t_{n-r,\alpha/2}\kappa_{n-r} \) in place of the normal point \( z_{\alpha/2} \) that
would be used if \( \sigma \) were known. For \( \alpha=0.05 \) that factor is \( 2.45 \) at
\( n-r=5 \) and \( 2.06 \) at \( n-r=20 \), against \( 1.96 \). Beyond about twenty residual degrees of freedom,
precision improves only through \( \blambda\T\G\blambda \), that is, through the design.

::: {#exm-ci-state-intervals}
[Murder rates: coefficients and a mean]

Return to the regression of the 2009 murder rate of the 50 states on the poverty rate, the
percentage of single-parent households and the percentage urban, with an intercept (@exm-proj-fwl-crime).
Here \( n=50 \), \( r=4 \), \( s^2=2.217 \), and
\( t_{46,0.025}=2.013 \). The \( 95\% \) intervals for the coefficients are

| Coefficient | Estimate | Standard error | 95% interval |
|---|---|---|---|
| intercept | \( -9.0946 \) | \( 1.5555 \) | \( (-12.226,\ -5.964) \) |
| poverty | \( 0.2538 \) | \( 0.0934 \) | \( (0.066,\ 0.442) \) |
| single parent | \( 0.3980 \) | \( 0.0829 \) | \( (0.231,\ 0.565) \) |
| urban | \( 0.0046 \) | \( 0.0131 \) | \( (-0.022,\ 0.031) \) |

The data are consistent with anything from \( 0.066 \) to \( 0.442 \) murders per 100,000 for each
point of poverty, with the other regressors held fixed in the sense of @thm-proj-fwl. The urban interval straddles
zero but is narrow, which says more than "not significant": any such association is at most \( 0.031 \)
per percentage point.

A mean response is an estimable function too. For a hypothetical state with poverty \( 18\% \),
single-parent households \( 28\% \) and \( 70\% \) urban, \( \x_0=(1,18,28,70)\T \), the estimated mean
murder rate is \( 6.937 \) with standard error \( 0.412 \), and the \( 95\% \) interval is
\( (6.107,\ 7.767) \). Here \( \x_0\T\G\x_0=0.0767 \). [Section 12.3](03-prediction.html) interprets this
number as the leverage of the new point.
:::

```{.python .run #cell-estimable-intervals-fit}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
r = np.linalg.matrix_rank(X)
G = np.linalg.pinv(X.T @ X)                     # any generalized inverse will do
beta_hat = G @ X.T @ y
s2 = np.sum((y - X @ beta_hat) ** 2) / (n - r)
q = stats.t.ppf(0.975, n - r)                   # upper 2.5% point of t(n - r)


def t_interval(lam):
    """95% interval for the estimable function lam^T beta."""
    est = lam @ beta_hat
    se = np.sqrt(s2 * lam @ G @ lam)
    return est, se, est - q * se, est + q * se


for j, name in enumerate(["intercept", "poverty", "single", "urban"]):
    est, se, lo, hi = t_interval(np.eye(p)[j])
    print(f"{name:9s} {est:8.4f}  se {se:.4f}  95% interval ({lo:.3f}, {hi:.3f})")
```

## What the standard error measures

For a coefficient in a full-rank model with an intercept, @eq-proj-vif-preview gives
\[
\text{se}(\hat{\beta}_j)=\frac{s}{\norm{(\I-\M_{(j)})\x_j}}=\frac{s}{\sqrt{S_{jj}(1-R_j^2)}},
\]
where \( \M_{(j)} \) projects onto the span of the other columns. The interval is short when \( \x_j \) has a large component orthogonal to the other
regressors; collinearity lengthens it through \( 1-R_j^2 \).

In general \( \blambda\T\G\blambda=\norm{\M\boldsymbol{\uprho}}^2 \), with \( \blambda=\X\T\boldsymbol{\uprho} \), is the smallest squared length of a
vector \( \mathbf{a} \) with \( \E(\mathbf{a}\T\Y)=\blambda\T\bbeta \) (@thm-opt-gauss-markov).

## Duality with tests

An interval and a family of tests carry the same information, in any model.

::: {#prp-ci-duality}
[Tests and confidence sets]

For each \( \theta_0\in\Theta \), let \( A(\theta_0) \) be the acceptance region of a test of
\( H:\theta=\theta_0 \) with \( \Pr_{\bbeta,\sigma^2}\{\Y\in A(\theta_0)\}\ge1-\alpha \) whenever
\( g(\bbeta,\sigma^2)=\theta_0 \). Then
\[
C(\y)=\{\theta_0\in\Theta:\y\in A(\theta_0)\}
\]
is a \( 1-\alpha \) confidence set for \( \theta \), and it is exact if every test has size exactly
\( \alpha \) at every parameter value in its hypothesis. Conversely, if \( C \) is a \( 1-\alpha \) confidence set,
the rule "reject \( H:\theta=\theta_0 \) iff \( \theta_0\notin C(\y) \)" is a level-\( \alpha \) test for every
\( \theta_0 \).
:::

::: {.proof}
For every parameter value, with \( \theta=g(\bbeta,\sigma^2) \), the events \( \{\theta\in C(\Y)\} \) and
\( \{\Y\in A(\theta)\} \) coincide, so they have the same probability. Both directions follow.
:::

For an estimable function, the two-sided \( t \) test of \( \blambda\T\bbeta=d \) accepts iff
\( \lvert\blambda\T\hbeta-d\rvert\le t_{n-r,\alpha/2}\,\text{se} \), which is the condition
\( d\in I(\y) \). So the interval @eq-ci-t-interval is exactly the set of values \( d \) that the \( t \) test
of @thm-glh-t-test does not reject. The \( p \)-value of the test of \( d \) is the smallest \( \alpha \)
for which the \( 1-\alpha \) interval excludes \( d \).

The interval is the more informative report. A test of \( \beta_j=0 \) reports one bit, confounded with
the sample size; the interval shows location and resolution. @exr-ci-equivalence uses it to answer "is the
effect negligible?", which a test of \( \beta_j=0 \) cannot answer at all.

## Functions that are not estimable

What if \( \blambda\notin\C(\X\T) \)? The recipe above fails because \( \blambda\T\hbeta \) depends on the
choice of \( \hbeta \). The failure is not a defect of the recipe. No procedure can do better.

::: {#prp-ci-nonestimable}
[No useful interval for a non-estimable function]

Let \( \blambda\notin\C(\X\T) \) and let \( C \) be a \( 1-\alpha \) confidence set for \( \blambda\T\bbeta \) with
\( \alpha<1 \). Then at every parameter value \( (\bbeta,\sigma^2) \) and for every real \( v \),
\( \Pr_{\bbeta,\sigma^2}\{v\in C(\Y)\}\ge1-\alpha \). If the sets \( C(\y) \) are such that
\( \{(\y,v):v\in C(\y)\} \) is measurable, the expected length of \( C(\Y) \) is infinite at every
parameter value.
:::

::: {.proof}
By @thm-est-characterization there is \( \bw\in\Null(\X) \) with \( \blambda\T\bw\ne0 \). Fix \( (\bbeta,\sigma^2) \)
and \( v \), and put \( \bbeta'=\bbeta+c\bw \) with \( c=(v-\blambda\T\bbeta)/\blambda\T\bw \). Then
\( \X\bbeta'=\X\bbeta \), so \( \Y \) has the same distribution under \( (\bbeta',\sigma^2) \) as under
\( (\bbeta,\sigma^2) \), while \( \blambda\T\bbeta'=v \). Hence
\( \Pr_{\bbeta,\sigma^2}\{v\in C(\Y)\}=\Pr_{\bbeta',\sigma^2}\{\blambda\T\bbeta'\in C(\Y)\}\ge1-\alpha \). For the
length, Fubini's theorem gives, for every \( a>0 \),
\[
\E\bigl[\text{length}\bigl(C(\Y)\cap[-a,a]\bigr)\bigr]=\int_{-a}^a\Pr\{v\in C(\Y)\}\,dv\ \ge\ 2a(1-\alpha),
\]
and the left side is at most the expected length of \( C(\Y) \). Let \( a\to\infty \).
:::

Thus the algebraic fact of
[Chapter 8](../ch08-estimability/index.html) becomes a statement about inference.

## An interval for \( \sigma^2 \)

::: {#prp-ci-sigma}
[Interval for the error variance]

The interval
\[
\Bigl[\ \frac{\text{SSE}}{\chi^2_{\alpha/2}(n-r)}\ ,\ \frac{\text{SSE}}{\chi^2_{1-\alpha/2}(n-r)}\ \Bigr]
\]
is an exact \( 1-\alpha \) confidence interval for \( \sigma^2 \), and the square roots of its endpoints give one for
\( \sigma \).
:::

::: {.proof}
\( \text{SSE}/\sigma^2\sim\chi^2(n-r) \) at every parameter value (@thm-opt-sampling(b)), so it is a pivot. The
event \( \chi^2_{1-\alpha/2}(n-r)\le\text{SSE}/\sigma^2\le\chi^2_{\alpha/2}(n-r) \) has probability \( 1-\alpha \),
and it is the event that \( \sigma^2 \) lies in the interval, which is the interval of @exr-opt-sigma-interval. Taking square roots is monotone.
:::

For the murder-rate regression, \( \text{SSE}=46s^2 \) with \( s^2=2.217 \), the chi-squared points are
\( 29.160 \) and \( 66.617 \), and the \( 95\% \) interval for \( \sigma^2 \) is
\( (1.531,\ 3.497) \). The interval is not symmetric about \( s^2 \), and it is not the shortest
one with this coverage (@exr-ci-shortest-sigma).

The two pivots of this section react differently to nonnormal errors. The \( t \) statistic needs
normality of a weighted average of errors, which a central limit theorem restores, and \( s\to\sigma \), which
holds for any error law. The chi-squared pivot needs the whole error distribution, because the variance of SSE
depends on the fourth moment (@prp-opt-var-sse). Take the fitted murder-rate model as the truth, with
\( \sigma^2=s^2 \), and simulate \( 40000 \) data sets. With normal errors the two \( 95\% \) intervals cover in
proportions \( 0.9514 \) (poverty coefficient) and \( 0.9491 \) (\( \sigma^2 \)). With errors from a
\( t(5) \) distribution rescaled to the same variance, the coefficient interval still covers
\( 0.9498 \) of the time, but the interval for \( \sigma^2 \) covers only \( 0.8165 \). More data do not
repair this (@exr-opt-sigma-interval).

```{.python .run #cell-estimable-intervals-mean}
x0 = np.array([1.0, 18.0, 28.0, 70.0])          # a state with high poverty and single parenthood
est0, se0, lo0, hi0 = t_interval(x0)
print(f"mean murder rate at x0: {est0:.3f}, se {se0:.3f}, interval ({lo0:.3f}, {hi0:.3f})")

sse = s2 * (n - r)
lo_s2 = sse / stats.chi2.ppf(0.975, n - r)      # interval for sigma^2 from SSE/sigma^2 ~ chi^2(n-r)
hi_s2 = sse / stats.chi2.ppf(0.025, n - r)
print(f"s^2 = {s2:.3f}; 95% interval for sigma^2 ({lo_s2:.3f}, {hi_s2:.3f})")
```

## Exercises

### A. Check your understanding

::: {#exr-ci-replicate}
[A1]

A design is replicated \( m \) times: \( \X \) is replaced by \( m \) stacked copies of itself. For a fixed
estimable \( \blambda\T\bbeta \), how does \( \blambda\T\G\blambda \) change? Roughly how many replicates are needed to
halve the expected length of the \( 95\% \) interval when \( n-r \) is already large, and why is the answer
slightly less than \( 4 \) when \( n-r \) is small?
:::

::: {.solution}
The stacked matrix has \( \X\T\X \) replaced by \( m\X\T\X \), and \( \G/m \) is a generalized inverse of
\( m\X\T\X \), so \( \blambda\T\G\blambda \) is divided by \( m \). By @thm-ci-estimable-interval(c) the expected
length is proportional to \( t_{\nu,\alpha/2}\kappa_\nu/\sqrt m \), with \( \nu=mn-r \). For large \( \nu \) the factor
\( t\kappa \) is nearly constant and \( m=4 \) halves the length; for small \( \nu \), replication also lowers \( t\kappa \), so
slightly less suffices.
:::

### B. Practice

::: {#exr-ci-kappa}
[B1]

Show that \( \kappa_\nu\to1 \) as \( \nu\to\infty \), and that \( s/\kappa_{n-r} \) is an unbiased estimator of
\( \sigma \). Compute \( \kappa_\nu \) for \( \nu=1 \) and \( \nu=2 \).
:::

::: {.solution}
\( V/\nu\to1 \) in probability by the law of large numbers, and \( \sqrt{V/\nu} \) is bounded in mean square, hence uniformly
integrable, so \( \E\sqrt{V/\nu}\to1 \). Unbiasedness:
\( \E s=\sigma\E\sqrt{V/(n-r)}=\sigma\kappa_{n-r} \). For \( \nu=1 \),
\( \kappa_1=\sqrt2\,\Gamma(1)/\Gamma(1/2)=\sqrt{2/\pi}\approx0.798 \). For \( \nu=2 \),
\( \kappa_2=\Gamma(3/2)/\Gamma(1)=\sqrt\pi/2\approx0.886 \).
:::

::: {#exr-ci-equivalence}
[B2]

Let \( \delta>0 \) be a margin below which \( \lvert\blambda\T\bbeta\rvert \) is considered negligible. The *two one-sided
tests* procedure declares equivalence when both \( H_1:\blambda\T\bbeta\le-\delta \) and \( H_2:\blambda\T\bbeta\ge\delta \)
are rejected by one-sided \( t \) tests at level \( \alpha \). Show that this happens iff the \( 1-2\alpha \) interval
@eq-ci-t-interval lies inside \( (-\delta,\delta) \). Show that the procedure has level \( \alpha \): when
\( \lvert\blambda\T\bbeta\rvert\ge\delta \), the probability of declaring equivalence is at most \( \alpha \).
:::

::: {.solution}
\( H_1 \) is rejected iff \( (\blambda\T\hbeta+\delta)/\text{se}>t_{n-r,\alpha} \), that is iff the lower end of the
\( 1-2\alpha \) interval, \( \blambda\T\hbeta-t_{n-r,\alpha}\text{se} \), exceeds \( -\delta \). Symmetrically \( H_2 \) is rejected iff
the upper end is below \( \delta \). For the level: if \( \blambda\T\bbeta\le-\delta \), declaring equivalence requires
rejecting \( H_1 \). Write \( \tau=\sigma\sqrt{\blambda\T\G\blambda} \), \( Z=(\blambda\T\hbeta-\blambda\T\bbeta)/\tau\sim\Normal(0,1) \) and
\( V=\text{SSE}/\sigma^2\sim\chi^2(\nu) \), \( \nu=n-r \), independent of \( Z \) (@thm-opt-sampling). The test statistic is
\( (Z+\eta)/\sqrt{V/\nu} \) with \( \eta=(\blambda\T\bbeta+\delta)/\tau \). For every realization of \( (Z,V) \) it increases in \( \eta \), so
the probability that it exceeds \( t_{\nu,\alpha} \) increases in \( \eta \). At \( \eta=0 \) that probability is \( \alpha \), so it is at most \( \alpha \)
whenever \( \blambda\T\bbeta\le-\delta \). The case \( \blambda\T\bbeta\ge\delta \) is symmetric.
:::

::: {#exr-ci-known-covariance}
[B3]

Suppose \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\V) \) with \( \V \) known and positive definite and \( \X \) of full column rank.
Using the generalized least squares estimator \( \tilde{\bbeta}=(\X\T\V^{-1}\X)^{-1}\X\T\V^{-1}\Y \) and
\( \tilde{s}^2=(\Y-\X\tilde{\bbeta})\T\V^{-1}(\Y-\X\tilde{\bbeta})/(n-p) \), derive an exact \( 1-\alpha \) interval for
\( \blambda\T\bbeta \). Write it out for a straight line \( \E Y_i=\beta_0+\beta_1x_i \) with \( \V=\diag(1/w_i) \).
:::

::: {.solution}
Let \( \V=\bL\bL\T \) (Cholesky) and \( \Z=\bL^{-1}\Y\sim\Normal_n(\bL^{-1}\X\bbeta,\sigma^2\I) \). The ordinary least squares
estimate for \( (\Z,\bL^{-1}\X) \) is \( \tilde{\bbeta} \), and its residual sum of squares is \( (n-p)\tilde{s}^2 \). Applying
@thm-ci-estimable-interval to the transformed model gives
\( \blambda\T\tilde{\bbeta}\pm t_{n-p,\alpha/2}\,\tilde{s}\sqrt{\blambda\T(\X\T\V^{-1}\X)^{-1}\blambda} \). For the weighted straight
line, with \( \bar{x}_w=\sum w_ix_i/\sum w_i \), \( \bar{y}_w \) defined likewise and \( S_w=\sum w_i(x_i-\bar{x}_w)^2 \), the slope
estimate is \( \tilde{\beta}_1=\sum w_i(x_i-\bar{x}_w)(y_i-\bar{y}_w)/S_w \) and its interval is
\( \tilde{\beta}_1\pm t_{n-2,\alpha/2}\,\tilde{s}/\sqrt{S_w} \).
:::

::: {#exr-ci-oneway-contrast}
[B4]

In a one-way layout with \( a \) groups of sizes \( n_1,\dots,n_a \) and the overparameterized model
\( \E Y_{ij}=\mu+\alpha_i \), show that the interval for \( \alpha_1-\alpha_2 \) is
\( \bar{y}_{1\cdot}-\bar{y}_{2\cdot}\pm t_{n-a,\alpha/2}\,s\sqrt{1/n_1+1/n_2} \), with \( n=\sum n_i \). What does
@prp-ci-nonestimable say about an interval for \( \alpha_1 \) alone?
:::

::: {.solution}
The rank is \( r=a \). With \( \blambda=\vect{e}_2-\vect{e}_3 \) (coordinates \( \mu,\alpha_1,\dots,\alpha_a \)), \( \blambda\T\bbeta=\alpha_1-\alpha_2 \)
is estimable, and \( \blambda\T\hbeta=\bar{y}_{1\cdot}-\bar{y}_{2\cdot} \) for every least squares estimate. Its variance is
\( \sigma^2(1/n_1+1/n_2) \) since the group means are independent, so \( \blambda\T\G\blambda=1/n_1+1/n_2 \). Substitute in
@eq-ci-t-interval with \( n-r=n-a \). The function \( \alpha_1 \) is not estimable (adding \( c \) to every \( \alpha_i \) and
subtracting it from \( \mu \) leaves \( \X\bbeta \) unchanged), so every \( 95\% \) confidence set for it covers every real number
with probability at least \( 0.95 \).
:::

### C. Going deeper

::: {#exr-ci-shortest-t}
[C1]

Among intervals \( [\blambda\T\hbeta-a\,\text{se},\ \blambda\T\hbeta+b\,\text{se}] \) with \( a,b\ge0 \) and coverage exactly \( 1-\alpha \),
show that the expected length is smallest when \( a=b=t_{n-r,\alpha/2} \). (Use that the \( t \) density is symmetric and
strictly decreasing in \( \lvert x\rvert \).)
:::

::: {.solution}
Coverage is \( F(b)-F(-a) \), with \( F \) the \( t(n-r) \) distribution function, and the expected length is proportional to
\( a+b \). If \( a<b \), shifting \( [-a,b] \) slightly left gains density \( f(a) \) and loses \( f(b)<f(a) \), raising the coverage, so a
strictly shorter interval has coverage \( 1-\alpha \). Hence \( a=b \).
:::

::: {#exr-ci-shortest-sigma}
[C2]

Intervals for \( \sigma^2 \) of the form \( [\text{SSE}/c_2,\ \text{SSE}/c_1] \) with \( \Pr(c_1\le\chi^2(\nu)\le c_2)=1-\alpha \) all have
exact coverage. Show that the one with the smallest expected length satisfies \( c_1^2f(c_1)=c_2^2f(c_2) \), where \( f \) is the
\( \chi^2(\nu) \) density, and that the equal-tailed choice of @prp-ci-sigma does not satisfy this condition.
:::
