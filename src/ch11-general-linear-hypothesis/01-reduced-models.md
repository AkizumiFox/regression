# Testing reduced models

Part II fitted linear models without asking whether a fitted term was needed. Part III adds
normal errors, and with them exact answers to that question. Does the murder rate of a state depend on its
region once poverty and family structure are accounted for? Is a straight line enough, or is
there curvature? Do two groups share a regression line? Is a coefficient equal to the value a
theory predicts? Each asks whether a *reduced model*, a special case of the model already
fitted, is adequate. This section sets up that comparison and derives its test. The rest of the
chapter looks at the same test from other directions.

## Hypotheses are subspaces

Throughout the chapter the data follow the normal linear model
\[
\Y\sim\Normal_n(\X\bbeta,\sigma^2\I),\qquad \bbeta\in\Real^p,\ \sigma^2>0,
\]
of @eq-opt-normal-model, with \( r=\rank(\X)<n \). We call it the **full model**. Its content is
that the mean vector \( \bmu=\E(\Y) \) lies in the \( r \)-dimensional subspace \( \C(\X) \)
([Section 6.1](../ch06-projections/01-nearest-point.html)). A reduced model says more: the
mean lies in a smaller subspace.

::: {#def-glh-reduced-model}
[Reduced model and linear hypothesis]

Let \( \X_0 \) be an \( n\times p_0 \) matrix with \( \C(\X_0)\subseteq\C(\X) \) and
\( r_0=\rank(\X_0)<r \). The **reduced model** is \( \E(\Y)=\X_0\mathbf{c} \),
\( \mathbf{c}\in\Real^{p_0} \). Inside the full model it corresponds to the **linear hypothesis**
\[
H_0:\ \X\bbeta\in\C(\X_0)
\qquad\text{against}\qquad
H_1:\ \X\bbeta\in\C(\X),\ \X\bbeta\notin\C(\X_0).
\]
The number \( q=r-r_0 \) is the number of **degrees of freedom of the hypothesis**, and we write
\( \M \) and \( \Mo \) for the orthogonal projections onto \( \C(\X) \) and \( \C(\X_0) \).
:::

Three features of the definition deserve comment. First, the hypothesis is a statement about the
mean vector, not about coefficients. Coefficients enter only through the column spaces, so
recoding the regressors of either model changes nothing (@thm-proj-reparam). Statements about
coefficients, such as \( \beta_2=\beta_3 \), are translated into subspaces in
[Section 11.3](03-testable-hypotheses.html). Second, the models must be *nested*. Two models
neither of which contains the other, say one with poverty and one with urbanization, are not
compared by the tests of this chapter. Third, both models share the same \( \sigma^2 \) and the
same error law. The alternative is not "anything other than the reduced model". It is "the
full model, but not the reduced one".

The most common reduced models come from three operations on \( \X \):

- *Deleting columns.* Dropping a block \( \X_2 \) from \( \X=[\X_1,\X_2] \) tests whether those
  regressors are needed given the others. With a single column this is the familiar question
  about one coefficient. With all columns except the intercept it is the "overall" test of
  [Section 11.6](06-coefficients.html).

- *Merging columns.* Replacing \( \x_2 \) and \( \x_3 \) by their sum \( \x_2+\x_3 \) forces equal
  coefficients. Merging the indicator columns of several groups forces those groups to share a
  mean.

- *Fixing a component.* Moving a known multiple of a column into an offset, for example testing
  that a slope equals \( 1 \), gives a reduced model of the form \( \E(\Y)=\X_0\mathbf{c}+\X\bb_0 \)
  with \( \bb_0 \) known. It is handled by subtracting \( \X\bb_0 \) from the data
  ([Section 11.3](03-testable-hypotheses.html)).

## The F statistic

If the reduced model is adequate, forcing the fit into \( \C(\X_0) \) should cost little: the
residual sum of squares \( \text{SSE}_0=\norm{(\I-\Mo)\y}^2 \) should exceed
\( \text{SSE}=\norm{(\I-\M)\y}^2 \) by no more than noise would produce. The difference is the
extra sum of squares of [Section 9.2](../ch09-sums-of-squares/02-extra-sums.html), and by
@eq-proj-extra-ss it is a squared length,
\[
\text{SSE}_0-\text{SSE}=\norm{(\M-\Mo)\y}^2=\y\T(\M-\Mo)\y ,
\]
the squared length of the projection of \( \y \) onto the orthogonal complement of \( \C(\X_0) \) within
\( \C(\X) \). How large is "more than noise"? Noise contributes \( \sigma^2 \) per dimension (@thm-ss-expected-mean-squares), and \( \M-\Mo \) has rank \( q \), so we compare the
**hypothesis mean square** \( (\text{SSE}_0-\text{SSE})/q \) with the unbiased estimate
\( s^2=\text{SSE}/(n-r) \) of \( \sigma^2 \) from the full model. The ratio
\[
F=\frac{(\text{SSE}_0-\text{SSE})/(r-r_0)}{\text{SSE}/(n-r)}
 =\frac{\y\T(\M-\Mo)\y/q}{\y\T(\I-\M)\y/(n-r)}
\]{#eq-glh-F}

is the **\( F \) statistic** for the reduced model. It is scale-free: multiplying \( \y \) by a
constant, or changing the units of the response, does not change it. Its exact distribution is
the content of the first theorem of the chapter.

::: {#thm-glh-f-test}
[The \( F \) test of a reduced model]

Assume the normal linear model @eq-opt-normal-model and the setting of @def-glh-reduced-model,
so that \( 0\le r_0<r<n \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \text{SSE}_0-\text{SSE} \) and \( \text{SSE} \) are independent, with
   \[
\frac{\text{SSE}_0-\text{SSE}}{\sigma^2}\sim\chi^2(r-r_0,\gamma),\qquad
\frac{\text{SSE}}{\sigma^2}\sim\chi^2(n-r);
\]

2. the statistic @eq-glh-F has the noncentral \( F \) distribution
   \[
F\sim F(r-r_0,\ n-r,\ \gamma),\qquad
\gamma=\frac{\norm{(\M-\Mo)\X\bbeta}^2}{\sigma^2};
\]{#eq-glh-noncentrality}

3. \( \gamma=\norm{(\I-\Mo)\X\bbeta}^2/\sigma^2=\min_{\mathbf{c}}\norm{\X\bbeta-\X_0\mathbf{c}}^2/\sigma^2 \),
   the squared distance from the mean vector to the reduced model space in units of
   \( \sigma^2 \). In particular \( \gamma=0 \) iff \( H_0 \) holds;

4. under \( H_0 \), \( F\sim F(r-r_0,n-r) \), the central \( F \) distribution, whatever the
   values of \( \bbeta \) and \( \sigma^2 \).
:::

:::

::: {.proof}
By @thm-proj-nested, \( \M-\Mo \) and \( \I-\M \) are symmetric idempotent with ranks \( r-r_0 \) and
\( n-r \), and \( (\M-\Mo)(\I-\M)=\bzero \). @thm-qf-orthogonal-projections, applied with
\( \bmu=\X\bbeta \), shows that \( \norm{(\M-\Mo)\Y}^2/\sigma^2 \) and \( \norm{(\I-\M)\Y}^2/\sigma^2 \) are
independent with distributions \( \chi^2(r-r_0,\norm{(\M-\Mo)\X\bbeta}^2/\sigma^2) \) and
\( \chi^2(n-r,\norm{(\I-\M)\X\bbeta}^2/\sigma^2) \). The second noncentrality is zero because
\( \M\X\bbeta=\X\bbeta \). With @eq-proj-extra-ss this is (a). Part (b) follows from (a) and
@def-qf-noncentral-f, since \( \sigma^2 \) cancels in the ratio. (It is also
@thm-qf-nested-f with \( \bP=\M \) and \( \bP_0=\Mo \).) For (c), \( \M\X\bbeta=\X\bbeta \) gives
\( (\M-\Mo)\X\bbeta=(\I-\Mo)\X\bbeta \), which is the difference between \( \X\bbeta \) and its nearest
point in \( \C(\X_0) \) (@thm-proj-projection-theorem). Its length is therefore the distance from
\( \X\bbeta \) to \( \C(\X_0) \), which is zero iff \( \X\bbeta\in\C(\X_0) \). Part (d) is (b) with
\( \gamma=0 \).
:::

Part (d) is what makes a test possible. The reduced model is a *composite* hypothesis: it
leaves \( \mathbf{c} \) and \( \sigma^2 \) unspecified. Yet the null distribution of \( F \) is the same at
every point of it. The unknown coefficients cancel because \( (\M-\Mo)\X_0\mathbf{c}=\bzero \) and
\( (\I-\M)\X_0\mathbf{c}=\bzero \), and the unknown scale cancels in the ratio. A statistic with this
property is called *pivotal* under the hypothesis.

Part (c) says what the test detects. The alternative is measured by one number, the squared
distance of the true mean from the reduced model, relative to the noise variance. Two
alternatives with the same \( \gamma \) are equally hard to detect, however different their
coefficients. [Section 11.5](05-power.html) is built on this.

## The test and its p-value

Large values of \( F \) are evidence against \( H_0 \): under the alternative the numerator mean
square estimates \( \sigma^2+\sigma^2\gamma/q \) while the denominator still estimates \( \sigma^2 \) (@eq-ss-ems). The **\( F \) test** at level \( \alpha \) rejects \( H_0 \) when
\[
F>F_\alpha(r-r_0,\,n-r),
\]
where \( F_\alpha(a,b) \) is the upper \( \alpha \) point of the central \( F(a,b) \) distribution.
By @thm-glh-f-test(d) its probability of rejecting is exactly \( \alpha \) at every point of the
reduced model. It rejects with probability greater than \( \alpha \) at every point of the
alternative, by (c) and @thm-qf-f-power. Instead of a fixed level, one usually reports the
**p-value**
\[
p(\y)=\Pr\{F(r-r_0,n-r)\ge F_{\text{obs}}\},
\]
the probability under the reduced model of a statistic at least as extreme as the one observed.

::: {#prp-glh-p-value}
[The p-value is uniform under the hypothesis]

Under \( H_0 \), the p-value \( p(\Y) \) has the uniform distribution on \( (0,1) \). For every
\( \alpha\in(0,1) \), rejecting when \( p(\Y)\le\alpha \) is the \( F \) test at level \( \alpha \).
:::

::: {.proof}
Let \( \Phi_F \) be the distribution function of \( F(q,n-r) \). It is continuous and strictly
increasing on \( (0,\infty) \), since the distribution has a positive density there, and
\( p(\y)=1-\Phi_F(F_{\text{obs}}) \). Under \( H_0 \) the statistic \( F \) has distribution function
\( \Phi_F \), so for \( u\in(0,1) \),
\( \Pr\{p(\Y)\le u\}=\Pr\{\Phi_F(F)\ge1-u\}=\Pr\{F\ge\Phi_F^{-1}(1-u)\}=u \). Finally,
\( p(\y)\le\alpha \) iff \( \Phi_F(F_{\text{obs}})\ge1-\alpha \) iff \( F_{\text{obs}}\ge F_\alpha(q,n-r) \),
and the event \( F=F_\alpha \) has probability zero.
:::

The p-value is a statement about the data under \( H_0 \). It is *not* the probability that
\( H_0 \) is true, and it is not the probability that the result is "due to chance". A small
p-value says that data like these would be unusual if the reduced model held.

::: {#exm-glh-region}
[Do regions matter?]

In the 2009 data on \( 50 \) US states used in [Chapter 6](../ch06-projections/index.html) and
[Chapter 9](../ch09-sums-of-squares/index.html), regress the murder rate on the poverty rate, the
percentage of single-parent households and the percentage living in urbanized areas, with an
intercept. This is the reduced model, with \( r_0=4 \). The full model adds indicators of three of
the four Census regions (Northeast, South, West, with the Midwest as reference), so \( r=7 \). The
question is whether the regions differ in murder rate once the three covariates are accounted for.
The listing gives

| Source | Degrees of freedom | Sum of squares | Mean square |
|:---|:---:|:---:|:---:|
| Regions, given the covariates | \( 3 \) | \( 14.83 \) | \( 4.944 \) |
| Residual, full model | \( 43 \) | \( 87.13 \) | \( 2.026 \) |
| Residual, reduced model | \( 46 \) | \( 101.96 \) | |

so \( F=4.944/2.026=2.440 \) on \( 3 \) and \( 43 \) degrees of freedom. The
\( 5\% \) critical value is \( 2.822 \) and the p-value is \( 0.077 \). At the
conventional level the reduced model is not rejected. The data do not show convincingly that the
regions differ, given the covariates.

Yet in the full model the Northeast coefficient is \( -1.736 \) murders per 100,000,
with a \( t \) statistic of \( -2.519 \) and a p-value of \( 0.016 \). One of the
three contrasts looks clearly nonzero while the joint test does not reject. There is no
contradiction: [Section 11.6](06-coefficients.html) shows how a joint test can dilute one strong
effect among weak ones, and [Chapter 13](../ch13-multiplicity/index.html) explains why a contrast singled out after looking at
the estimates cannot be judged by its own p-value.
:::

```{.python .run #cell-region-test-ftest}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
# Census region of each state, in the data's alphabetical order (N, M, S, W)
codes = "SWWSWWNSSSWWMMMMSSNSNMMSMWMWNNWNSMMSWNNSMSSWNSWSMW"
region = np.array(list(codes))
y = data["murder"].to_numpy()
n = len(y)
covariates = np.column_stack([data["poverty"], data["single"], data["urban"]])
X0 = np.column_stack([np.ones(n), covariates])                    # reduced model
D = np.column_stack([(region == g).astype(float) for g in "NSW"])  # Midwest is the reference
X = np.column_stack([X0, D])                                      # full model

def sse(Z):
    """Residual sum of squares of the least squares fit of y on the columns of Z."""
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return np.sum((y - Z @ b) ** 2)

r, r0 = np.linalg.matrix_rank(X), np.linalg.matrix_rank(X0)
sse_full, sse_reduced = sse(X), sse(X0)
F = ((sse_reduced - sse_full) / (r - r0)) / (sse_full / (n - r))
p_value = stats.f.sf(F, r - r0, n - r)
crit = stats.f.ppf(0.95, r - r0, n - r)
print(f"SSE0 = {sse_reduced:.2f}, SSE = {sse_full:.2f}, df = ({r - r0}, {n - r})")
print(f"F = {F:.3f}, p = {p_value:.4f}, 5% critical value = {crit:.3f}")
```

## What the test assumes

@thm-glh-f-test rests on three assumptions.

*The full model must hold.* The denominator \( s^2 \) estimates \( \sigma^2 \) only if
\( \E(\Y)\in\C(\X) \). If the full model omits something, \( \E(s^2) \) exceeds \( \sigma^2 \) by
\( \norm{(\I-\M)\bmu}^2/(n-r) \) (@thm-ss-expected-mean-squares(b)), the denominator is inflated,
and \( F \) is pushed towards zero. The test then loses power, and a nonsignificant result says even
less than usual. The test of a reduced model is a test *within* the full model. It cannot certify
the full model itself. [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) tests the full model against a nonparametric alternative when
there are replicate observations.

*The errors must be independent with a common variance.* Correlated or heteroscedastic errors
change both numerator and denominator, and the null distribution is no longer \( F \). [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html)
and [Chapter 31](../ch31-general-gauss-markov/index.html) treat these cases.

*Normality* gives the exact distribution, and it is the least critical of the three. The null
distribution is exactly \( F \) for every spherically symmetric error law
([Section 11.2](02-comparing-projections.html)), and approximately so for independent non-normal errors
unless the sample is small or some points have high leverage (@thm-dep-nonnormal).

::: {.warning}
[What a nonsignificant \( F \) does not show]

Failing to reject the reduced model is *not* evidence that the reduced model is true. The test
controls the probability of wrongly rejecting \( H_0 \). It says nothing on its own about the
probability of wrongly accepting it, which depends on the unknown noncentrality \( \gamma \) and may
be large. In @exm-glh-region, \( p=0.077 \) is compatible with no regional differences.
It is equally compatible with sizeable regional differences.
[Section 11.5](05-power.html) computes that the test had only about an even chance of detecting a
Northeast shortfall of \( 1.5 \) murders per 100,000, about the size estimated in
[Section 11.3](03-testable-hypotheses.html) (\( 1.379 \) below the other three regions). To support "the reduced model is
adequate", one must show that every alternative of practical importance would have been rejected.
That calls for a power calculation made before the data were seen, or better, a confidence region
for the omitted effects that excludes all important values ([Chapter 12](../ch12-intervals-and-bands/index.html)), or a test of equivalence in
which the roles of hypothesis and alternative are reversed. Nor does a nonsignificant \( F \) license
dropping the terms and then reporting the reduced fit as if it had been specified in advance.
:::

::: {.warning}
[What a significant \( F \) does not show]

A significant \( F \) says that the mean vector is not in \( \C(\X_0) \), assuming the full model holds.
It does not say which of the \( q \) dimensions of the test space carries the departure, nor that every
omitted coefficient is nonzero. It does not say that the departure is large enough to matter: with
large \( n \), a tiny \( \gamma/n \) is detected with certainty ([Section 11.5](05-power.html)). And it
does not say that the full model is correct. The statistic only compares two models, and a third,
better model may exist.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-invariance-shift}
[A1]

Show that the \( F \) statistic @eq-glh-F does not change when \( \y \) is replaced by
\( a\y+\X_0\mathbf{c} \) for any \( a\ne0 \) and any vector \( \mathbf{c} \). Interpret the result for the
region test of @exm-glh-region when murder rates are expressed per million instead of per 100,000.
:::

::: {#exr-glh-reduced-examples}
[A2]

In the model \( \E(Y_i)=\beta_0+\beta_1x_{i1}+\beta_2x_{i2}+\beta_3x_{i3} \), write down a matrix
\( \X_0 \) (or an offset and an \( \X_0 \)) for each reduced model and give \( q \):
(i) \( \beta_2=\beta_3=0 \); (ii) \( \beta_1=\beta_2=\beta_3 \); (iii) \( \beta_1=2 \);
(iv) \( \beta_2+\beta_3=0 \).
:::

### B. Practice

::: {#exr-glh-known-sigma}
[B1]

Suppose \( \sigma^2=\sigma_0^2 \) is known. Show that \( X^2=(\text{SSE}_0-\text{SSE})/\sigma_0^2 \) has the
\( \chi^2(q,\gamma) \) distribution and that the test rejecting when \( X^2 \) exceeds the upper \( \alpha \)
point of \( \chi^2(q) \) has level \( \alpha \). Show that \( qF/X^2\to1 \) in probability as \( n\to\infty \)
with \( q \) fixed and the full model true, and explain why the \( F \) test and this test then behave
alike.
:::

::: {.solution}
The distribution is @thm-glh-f-test(a) with \( \sigma^2=\sigma_0^2 \); under \( H_0 \) it is \( \chi^2(q) \), so
the level is exact. Now \( qF/X^2=\sigma_0^2/s^2 \), and \( s^2\to\sigma^2=\sigma_0^2 \) in probability, since
\( \E(s^2)=\sigma^2 \) and \( \Var(s^2)=2\sigma^4/(n-r)\to0 \) (@thm-opt-sampling(d)) and Chebyshev's
inequality applies. So \( qF \) and \( X^2 \) differ by a factor tending to one, and
\( qF_\alpha(q,n-r) \) tends to the upper \( \alpha \) point of \( \chi^2(q) \) for the same reason. With
few residual degrees of freedom the two tests differ, and only the \( F \) test has the right level
when \( \sigma^2 \) must be estimated.
:::

::: {#exr-glh-small-F}
[B2]

Suppose the full model is wrong: \( \E(\Y)=\bmu \) with \( (\I-\M)\bmu=\boldsymbol{\updelta}\ne\bzero \), but
\( (\M-\Mo)\bmu=\bzero \). Show that the numerator and denominator of @eq-glh-F are still independent,
that \( \E(\text{numerator})=\sigma^2 \) and \( \E(\text{denominator})=\sigma^2+\norm{\boldsymbol{\updelta}}^2/(n-r) \), and
that \( \Pr(F\le c) \ge \Pr\{F(q,n-r)\le c\} \) for every \( c>0 \). Explain why an \( F \) statistic far
*below* one, with a p-value near one, is a hint that the full model is inadequate.
:::

::: {.solution}
Independence and the distributions follow from @thm-qf-orthogonal-projections as in the proof of
@thm-glh-f-test: the numerator is \( \sigma^2\chi^2(q)/q \), and \( \text{SSE}/\sigma^2\sim\chi^2(n-r,\norm{\boldsymbol{\updelta}}^2/\sigma^2) \),
with mean \( n-r+\norm{\boldsymbol{\updelta}}^2/\sigma^2 \) (@thm-qf-ncchisq(b)). For the probability, condition on
the numerator \( U \): \( F\le c \) iff \( \text{SSE}/\sigma^2\ge (n-r)U/(qc\sigma^2) \), and a noncentral
\( \chi^2 \) exceeds any fixed value with at least the central probability (@prp-qf-ncchisq-monotone). Averaging over \( U \) gives the inequality. So omitted structure in the
residual space makes small values of \( F \) more likely. Under the full model, \( \Pr(F\le c) \) is
small for small \( c \), and an observed \( F \) in the far lower tail is surprising.
:::

::: {#exr-glh-chain}
[B3]

Let \( \C(\X_0)\subseteq\C(\X_1)\subseteq\C(\X) \) with ranks \( r_0<r_1<r<n \), and let
\( s^2=\text{SSE}/(n-r) \) come from the largest model. Show that
\[
F_1=\frac{\y\T(\M_1-\Mo)\y/(r_1-r_0)}{s^2}
\]
has the \( F(r_1-r_0,n-r,\gamma_1) \) distribution with \( \gamma_1=\norm{(\M_1-\Mo)\X\bbeta}^2/\sigma^2 \)
whenever the largest model holds, and that its numerator is independent of that of the test of
\( \X_1 \) within \( \X \). Why is this preferable to using the residual mean square of \( \X_1 \)?
:::

::: {#exr-glh-chow}
[B4]

Two samples of sizes \( n_1 \) and \( n_2 \) follow \( \E(\Y_k)=\X_k\bbeta_k \), \( k=1,2 \), with a common
\( \sigma^2 \) and \( \X_k \) of full column rank \( p \). Write the model for both samples together as a
full model with block-diagonal model matrix, and the hypothesis \( \bbeta_1=\bbeta_2 \) as a reduced
model. Show that
\[
F=\frac{\bigl(\text{SSE}_{\text{pooled}}-\text{SSE}_1-\text{SSE}_2\bigr)/p}{(\text{SSE}_1+\text{SSE}_2)/(n_1+n_2-2p)},
\]
where \( \text{SSE}_k \) comes from fitting sample \( k \) alone and \( \text{SSE}_{\text{pooled}} \) from one fit to
all the data. Apply the test to the question whether two firms of the Grunfeld panel
(`sm.datasets.grunfeld`) share their investment equation, and comment on the assumption of a
common \( \sigma^2 \).
:::

### C. Going deeper

::: {#exr-glh-chow-few}
[C1]

In @exr-glh-chow suppose instead that \( n_2\le p \) and \( \rank(\X_2)=n_2 \), so that \( \bbeta_2 \)
cannot be estimated from the second sample. Show that the full model then fits the second sample
exactly, identify the test space, and show that the \( F \) statistic of @thm-glh-f-test becomes
\[
F=\frac{(\text{SSE}_{\text{pooled}}-\text{SSE}_1)/n_2}{\text{SSE}_1/(n_1-p)},
\]
a test of whether the second sample is predicted by the fit to the first, on \( n_2 \) and
\( n_1-p \) degrees of freedom.
:::

