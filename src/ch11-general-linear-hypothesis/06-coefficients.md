# Testing single coefficients and groups of coefficients

The tests most often printed by regression software are instances of the general theory: a
\( t \) test for each coefficient, an \( F \) test that all slopes vanish, and on request an \( F \) test for a
group of coefficients. This section derives each from the preceding sections, explains what each
one tests, and studies how a group test and the individual tests inside it can disagree, in both
directions.

## The t test

When the hypothesis consists of a single estimable function, \( q=1 \), the \( F \) statistic is the square of
a statistic with a sign, and the sign allows one-sided alternatives.

::: {#thm-glh-t-test}
[The \( t \) test for an estimable function]

Assume @eq-opt-normal-model with \( r<n \). Let \( \blambda\in\C(\X\T) \), \( \blambda\ne\bzero \), let \( d\in\Real \), let
\( \G \) be a generalized inverse of \( \X\T\X \), and put
\[
T=\frac{\blambda\T\hbeta-d}{s\sqrt{\blambda\T\G\blambda}},\qquad s^2=\frac{\text{SSE}}{n-r}.
\]{#eq-glh-t}

::: {.enumerate options="label=(\alph*)"}
1. \( T\sim t(n-r,\delta) \) with \( \delta=(\blambda\T\bbeta-d)/(\sigma\sqrt{\blambda\T\G\blambda}) \).

2. \( T^2 \) equals the statistic \( F_H \) of @eq-glh-general-F for \( H:\blambda\T\bbeta=d \), and
   \( \gamma=\delta^2 \). The two-sided test that rejects when \( \lvert T\rvert>t_{n-r,\alpha/2} \) is the \( F \) test of
   \( H \), with the same p-value.

3. The one-sided test that rejects when \( T>t_{n-r,\alpha} \) has size exactly \( \alpha \) for the hypothesis
   \( \blambda\T\bbeta\le d \): its rejection probability is at most \( \alpha \) whenever \( \blambda\T\bbeta\le d \), with
   equality when \( \blambda\T\bbeta=d \). Its power increases strictly with \( \blambda\T\bbeta \).
:::

:::

::: {.proof}
(a) is @cor-opt-t. (b) With \( \bLambda=\blambda \), \( \W=\blambda\T\G\blambda>0 \) is a scalar, so
\( F_H=(\blambda\T\hbeta-d)^2/(\blambda\T\G\blambda\,s^2)=T^2 \), and @eq-glh-general-noncentrality gives
\( \gamma=\delta^2 \). Under \( H \), \( T\sim t(n-r) \) is symmetric about zero, so
\( \Pr\{\lvert T\rvert>t\}=\Pr\{T^2>t^2\} \) for every \( t>0 \). Hence \( t_{n-r,\alpha/2}^2=F_\alpha(1,n-r) \), the two
rejection regions coincide, and the two-sided \( t \) p-value \( \Pr\{\lvert T\rvert\ge\lvert T_{\text{obs}}\rvert\} \) equals
the \( F \) p-value \( \Pr\{T^2\ge T^2_{\text{obs}}\} \). (c) By @exr-qf-t-square, the probability that a
\( t(n-r,\delta) \) variable exceeds a fixed constant is strictly increasing in \( \delta \), and \( \delta \) is an
increasing function of \( \blambda\T\bbeta \). At \( \delta=0 \) the probability is \( \alpha \).
:::

For a single coefficient \( \beta_j \) in a full-rank model, @eq-glh-t becomes the familiar
\( t_j=\hat{\beta}_j/\{s\sqrt{c_{jj}}\} \), with \( c_{jj} \) the \( j \)th diagonal entry of \( (\X\T\X)^{-1} \). By
@cor-ss-single-column, \( t_j^2 \) is the \( F \) statistic for dropping column \( j \) from the model. The \( t \) test
of \( \beta_j=0 \) therefore compares the model with all regressors to the model with all *except* \( \x_j \). It
asks whether \( \x_j \) adds anything *given the others*. By the Frisch–Waugh–Lovell theorem (@thm-proj-fwl), what it measures is the association of \( \y \) with the part of \( \x_j \) orthogonal to the other
columns. The same coefficient can be "significant" in one model and not in another, and there is no
contradiction: the two tests answer different questions.

The one-sided test is appropriate only when the direction was specified before seeing the data, for
example when theory says a coefficient cannot be negative, or when only an increase would lead to action.
Choosing the direction after seeing the sign of the estimate doubles the actual size.

## Tests and intervals are one computation

The two-sided test of \( \blambda\T\bbeta=d \) at level \( \alpha \) *accepts* exactly when
\( \lvert\blambda\T\hbeta-d\rvert\le t_{n-r,\alpha/2}\,s\sqrt{\blambda\T\G\blambda} \). Collecting the values of \( d \) that are
accepted gives the interval
\[
\blambda\T\hbeta\pm t_{n-r,\alpha/2}\,s\sqrt{\blambda\T\G\blambda},
\]
which is the \( 1-\alpha \) confidence interval @eq-opt-t-interval. Every interval of this kind is the set of
hypothesized values that a test does not reject, and every test can be read off an interval: reject
\( \blambda\T\bbeta=d \) iff \( d \) lies outside. The interval reports more. It shows which values are ruled out,
and so how precisely the data determine \( \blambda\T\bbeta \), which the verdict of a single test does not.
For several functions at once, inverting the \( F \) test of @thm-glh-general-f gives an ellipsoid.
[Chapter 12](../ch12-intervals-and-bands/index.html) develops both and the duality between them.

## The overall F test

The most common reduced model is the one with no regressors at all. Let \( \X \) contain the intercept column
\( \bone \) and have rank \( r \). The reduced model \( \E(\Y)=\beta_0\bone \) has \( \text{SSE}_0=\text{SST}=\sum_i(y_i-\bar{y})^2 \),
and \( \text{SST}-\text{SSE}=\text{SSR} \) is the regression sum of squares. By @thm-glh-f-test and
@thm-ss-r2-null,
\[
F=\frac{\text{SSR}/(r-1)}{\text{SSE}/(n-r)}=\frac{R^2/(r-1)}{(1-R^2)/(n-r)}\ \sim\ F(r-1,\,n-r,\,\gamma),
\qquad \gamma=\frac{\norm{(\M-n^{-1}\bone\bone\T)\X\bbeta}^2}{\sigma^2} .
\]
The **overall \( F \) test** asks whether any linear combination of the regressors helps to predict the
response. The answer is almost always yes in observational data with a sensible choice of regressors, and a
significant overall \( F \) is correspondingly weak information. It does not show that any particular regressor
matters, that the model fits well, or that it predicts usefully: with \( n \) large, a model with
\( R^2=0.02 \) has a highly significant overall \( F \). A nonsignificant overall \( F \) is more informative. It says
that the data give no evidence that the regressors, taken together, relate to the response at all.

## Groups of coefficients

Between one coefficient and all of them lies the test of a group: whether the terms of a factor, the
powers of a polynomial beyond the linear, or a set of related regressors can be dropped together. The
statistic is the extra sum of squares of the group, divided by its degrees of freedom and by \( s^2 \) (@thm-glh-f-test). How does it relate to the \( t \) statistics of its members?

::: {#prp-glh-max-t}
[The group statistic is the largest single statistic]

Let \( H:\bLambda\T\bbeta=\bm d \) be testable with \( \bLambda \) of full column rank \( q \), let \( F_H \) be its statistic
@eq-glh-general-F, and for \( \bm a\in\Real^q \), \( \bm a\ne\bzero \), let \( T(\bm a) \) be the \( t \) statistic @eq-glh-t for
the single function \( \bm a\T\bLambda\T\bbeta=\bm a\T\bm d \). Then
\[
qF_H=\max_{\bm a\ne\bzero}T(\bm a)^2,
\]
and the maximum is attained at \( \bm a\propto\W^{-1}(\bLambda\T\hbeta-\bm d) \), where \( \W=\bLambda\T\G\bLambda \). In
particular \( t_j^2\le qF_H \) for each coefficient \( \beta_j \) in a tested group \( \beta_j=0 \), \( j\in J \), with \( q=\lvert J\rvert \).
:::

::: {.proof}
Write \( \bu=\bLambda\T\hbeta-\bm d \). The function \( \bm a\T\bLambda\T\bbeta \) is estimable with
\( \blambda=\bLambda\bm a\ne\bzero \) (the columns of \( \bLambda \) are independent), and
\( \blambda\T\G\blambda=\bm a\T\W\bm a \), so \( T(\bm a)^2=(\bm a\T\bu)^2/(s^2\,\bm a\T\W\bm a) \). The matrix \( \W \) is positive
definite (@thm-ss-hypothesis(a)) and has a positive definite square root (@thm-mat-square-root). By the
Cauchy–Schwarz inequality (@prp-mat-cauchy-schwarz) applied to \( \W^{1/2}\bm a \) and \( \W^{-1/2}\bu \),
\[
(\bm a\T\bu)^2\le(\bm a\T\W\bm a)(\bu\T\W^{-1}\bu),
\]
with equality iff \( \W^{1/2}\bm a \) is proportional to \( \W^{-1/2}\bu \), that is, \( \bm a\propto\W^{-1}\bu \) (or \( \bu=\bzero \)).
Dividing by \( s^2\,\bm a\T\W\bm a \) gives \( T(\bm a)^2\le\bu\T\W^{-1}\bu/s^2=qF_H \), with equality at the stated
\( \bm a \). For the coefficient \( \beta_j \), take \( \bm a \) to be the corresponding coordinate vector.
:::

The proposition is the geometric heart of Scheffé's method for multiple comparisons, and it explains the two
ways in which a group test and its members can disagree.

*Dilution.* A single large \( t_j \) need not make the group significant. The group test rejects when
\( qF_H>qF_\alpha(q,n-r) \), and \( qF_\alpha(q,n-r) \) grows with \( q \) (@exr-glh-critical-growth). A \( t_j^2 \) just above
\( F_\alpha(1,n-r) \) can therefore leave \( qF_H \) below the group's critical value, especially when the other members
of the group contribute nothing. This is what happened in @exm-glh-region: the Northeast contrast had
\( \lvert t\rvert=2.519 \), but the three-degree-of-freedom test of all regional differences was not significant.

*Collinearity.* Conversely, the group can be highly significant while every member has a small \( t \). The
direction \( \bm a \) that attains the maximum is generally not a coordinate direction. When the estimates are strongly
correlated, the data can say clearly that *some* combination of the coefficients is nonzero without being able to
attribute it to any single one.

::: {#exm-glh-state-coefficients}
[Five regressors for the murder rate]

Regress the murder rate of the \( 50 \) states on poverty, high-school graduation, single-parent households,
percentage white and urbanization, with an intercept (\( p=6 \), \( n-p=44 \)).

| Regressor | Estimate | Standard error | \( t \) | p-value |
|:---|:---:|:---:|:---:|:---:|
| intercept | \( -10.0586 \) | \( 17.6785 \) | \( -0.569 \) | \( 0.572 \) |
| poverty | \( 0.2378 \) | \( 0.1412 \) | \( 1.684 \) | \( 0.099 \) |
| high-school graduates | \( -0.0047 \) | \( 0.1530 \) | \( -0.030 \) | \( 0.976 \) |
| single-parent households | \( 0.4280 \) | \( 0.1372 \) | \( 3.120 \) | \( 0.003 \) |
| white | \( 0.0108 \) | \( 0.0270 \) | \( 0.398 \) | \( 0.692 \) |
| urban | \( 0.0046 \) | \( 0.0145 \) | \( 0.315 \) | \( 0.754 \) |

The overall test gives \( R^2=0.643 \) and \( F=15.87 \) on \( 5 \) and \( 44 \) degrees of freedom,
\( p=6.4\times10^{-9} \). Only single parenthood has a \( t \) statistic beyond the critical value
\( 2.015 \). The group of three demographic and educational regressors, high-school graduation,
percentage white and urbanization, has \( F=0.096 \) on \( 3 \) and \( 44 \) degrees of freedom,
\( p=0.962 \): the data give no evidence that these three add anything once poverty and single
parenthood are in the model. Every \( t_j^2 \) of the group is at most \( 3F=0.29 \), as @prp-glh-max-t requires.

This group was chosen here because its members have small \( t \) statistics. A group selected that way
will tend to have a small \( F \), and its large p-value makes dropping the group look safer than it is. Groups should be
defined by the scientific question, before the individual results are seen.
:::

```{.python .run #cell-coefficient-tests-coefficients}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.statecrime.load_pandas().data
data.index = data.index.str.strip()
data = data.drop(index="District of Columbia")
names = ["poverty", "hs_grad", "single", "white", "urban"]
y = data["murder"].to_numpy()
n = len(y)
X = np.column_stack([np.ones(n), data[names]])
p = X.shape[1]

C = np.linalg.inv(X.T @ X)
b = C @ X.T @ y
sse = np.sum((y - X @ b) ** 2)
s2 = sse / (n - p)
se = np.sqrt(s2 * np.diag(C))
t = b / se
pv = 2 * stats.t.sf(np.abs(t), n - p)
for name, bj, sj, tj, pj in zip(["intercept"] + names, b, se, t, pv):
    print(f"{name:10s} {bj:9.4f} {sj:8.4f} {tj:7.3f} {pj:7.4f}")
```

```{.python .run #cell-coefficient-tests-overall}
sst = np.sum((y - y.mean()) ** 2)
R2 = 1 - sse / sst
F_all = (R2 / (p - 1)) / ((1 - R2) / (n - p))

def extra_F(keep):
    """F test that the coefficients not in `keep` (columns of X) are zero."""
    Xr = X[:, keep]
    br, *_ = np.linalg.lstsq(Xr, y, rcond=None)
    sse_r = np.sum((y - Xr @ br) ** 2)
    q = p - len(keep)
    F = ((sse_r - sse) / q) / s2
    return F, stats.f.sf(F, q, n - p)

F_group, p_group = extra_F([0, 1, 3])          # drop hs_grad, white, urban
print(f"overall: R2 = {R2:.3f}, F = {F_all:.2f} on ({p - 1}, {n - p})")
print(f"hs_grad, white, urban: F = {F_group:.3f}, p = {p_group:.3f}")
```

::: {#exm-glh-longley}
[Jointly significant, individually not]

Longley's data record US employment for the \( 16 \) years 1947 to 1962, with five economic series: the GNP
deflator, GNP, unemployment, the size of the armed forces and population. Regress employment on the five
series and an intercept, leaving out the time trend, so that \( n-p=10 \). For unemployment and
population, the \( t \) statistics are \( -0.921 \) and \( -1.222 \), with p-values
\( 0.379 \) and \( 0.250 \). Neither comes near the critical value \( 2.228 \). The
test that *both* coefficients are zero gives \( F=7.175 \) on \( 2 \) and \( 10 \) degrees of freedom,
\( p=0.0117 \).

[Figure 11.6.1](#fig-glh-joint) shows why. The two estimates have correlation \( -0.840 \).
Under the hypothesis, their \( t \) statistics are then likely to have opposite signs, and the \( F \) test's
acceptance region is an ellipse stretched along the line \( t_1=-t_2 \). The observed pair has both statistics
negative, a combination the ellipse excludes although each coordinate on its own is unremarkable. The
data are clear that unemployment and population *together* matter, and unclear about how to divide the
credit. With \( 16 \) consecutive years of trending series, the errors are also unlikely to be independent, a
caution taken up in Chapter 21.
:::

::: {when-format="html"}
![**Figure 11.6.1.** Acceptance regions for the hypothesis that the unemployment and population
coefficients in Longley's regression are both zero, drawn in the plane of the two \( t \) statistics. The
square accepts when both individual \( t \) tests accept at level \( 0.05 \). The ellipse accepts when the joint
\( F \) test does. The black point is the data.](joint_vs_individual.svg){#fig-glh-joint width=60%}
:::

::: {when-format="pdf"}
![Acceptance regions for the hypothesis that the unemployment and population
coefficients in Longley's regression are both zero, drawn in the plane of the two \( t \) statistics. The
square accepts when both individual \( t \) tests accept at level \( 0.05 \). The ellipse accepts when the joint
\( F \) test does. The black point is the data.](joint_vs_individual.pdf){width=60%}
:::

```{.python .run #cell-coefficient-tests-longley}
lg = sm.datasets.longley.load_pandas()
Z = sm.add_constant(lg.exog.drop(columns="YEAR"))
res = sm.OLS(lg.endog, Z).fit()
joint = res.f_test("UNEMP = 0, POP = 0")
print(res.tvalues[["UNEMP", "POP"]].round(3).to_dict(), res.pvalues[["UNEMP", "POP"]].round(3).to_dict())
print(f"joint F = {float(joint.fvalue):.3f}, p = {float(joint.pvalue):.4f}")
```

The figure also shows the opposite disagreement. Points inside the ellipse but outside the square have
one significant \( t \) and a nonsignificant \( F \). Neither verdict is wrong: the individual tests ask about the
coordinate directions, and the joint test asks about all directions at once, paying for that breadth with a
larger critical value.

::: {.warning}
[Several t tests are not a test of the group]

Rejecting the group hypothesis "all \( \beta_j=0 \), \( j\in J \)" whenever at least one \( \lvert t_j\rvert \) exceeds
\( t_{n-r,\alpha/2} \) is not a level-\( \alpha \) test: when the \( k \) coefficient estimates are independent and \( n-r \) is large,
its size is about \( 1-(1-\alpha)^k \), about \( 0.23 \) for \( k=5 \) at \( \alpha=0.05 \). (The \( t \) statistics
share the estimate \( s \), so they are not exactly independent.) Reporting the smallest of several p-values as
if it were the only one overstates the evidence in the same way. The \( F \) test of the group, or a
multiple-comparison procedure designed for the family ([Chapter 13](../ch13-multiplicity/index.html)),
controls the familywise error rate. Choosing regressors by their \( t \) statistics and then testing the chosen model on
the same data invalidates all the p-values of the final fit. Chapter 29 returns to this problem.
:::

## Exercises

### A. Check your understanding

::: {#exr-glh-poverty-t}
[A1]

Using the table of @exm-glh-state-coefficients, compute the \( F \) statistic for dropping poverty from the
five-regressor model and its p-value, without refitting.
:::

::: {#exr-glh-one-sided-choice}
[A2]

Show that if the direction of a one-sided \( t \) test is chosen after seeing the sign of \( T \), and the test then
rejects when \( \lvert T\rvert>t_{n-r,\alpha} \), its size is \( 2\alpha \).
:::

### B. Practice

::: {#exr-glh-equal-coefficients}
[B1]

In a full-rank model, show that the \( t \) statistic for \( \beta_j=\beta_k \) is
\[
T=\frac{\hat{\beta}_j-\hat{\beta}_k}{s\sqrt{c_{jj}+c_{kk}-2c_{jk}}},
\]
where \( c_{jk} \) are the entries of \( (\X\T\X)^{-1} \). Show that it equals the \( t \) statistic of the coefficient of
\( \x_j \) in the reparameterized model with regressors \( \x_j \) and \( \x_j+\x_k \) in place of \( \x_j \) and \( \x_k \).
:::

::: {#exr-glh-critical-growth}
[B2]

Show that \( qF_\alpha(q,\nu) \) is strictly increasing in \( q \) for fixed \( \nu \) and \( \alpha \). Deduce that a single \( t \)
statistic with \( t^2 \) just above \( F_\alpha(1,\nu) \) does not force the rejection of a group containing it.
:::

::: {.solution}
\( qF(q,\nu) \) has the law of \( U_q/(V/\nu) \) with \( U_q\sim\chi^2(q) \) independent of \( V\sim\chi^2(\nu) \). Taking
\( U_{q+1}=U_q+W \) with \( W\sim\chi^2(1) \) independent of \( U_q \) and \( V \), we have \( U_{q+1}/(V/\nu)>U_q/(V/\nu) \) with
probability one, so \( \Pr\{U_{q+1}/(V/\nu)>c\}>\Pr\{U_q/(V/\nu)>c\} \) for every \( c>0 \) (strictly, since
\( \Pr\{U_q/(V/\nu)\le c<U_{q+1}/(V/\nu)\}>0 \)). Taking \( c=qF_\alpha(q,\nu) \), where the right side equals \( \alpha \),
shows that the upper \( \alpha \) point of \( U_{q+1}/(V/\nu) \) exceeds \( qF_\alpha(q,\nu) \). For the deduction, if
\( t_j^2 \) lies between \( F_\alpha(1,\nu) \) and \( qF_\alpha(q,\nu) \) and the other members contribute nothing, then
\( qF_H \) can equal \( t_j^2<qF_\alpha(q,\nu) \), and the group is not rejected.
:::

::: {#exr-glh-longley-F}
[B3]

For two coefficients with \( t \) statistics \( t_1,t_2 \) whose estimates have correlation \( \rho \), show that the
\( F \) statistic for both being zero is
\[
F=\frac{t_1^2+t_2^2-2\rho t_1t_2}{2(1-\rho^2)} .
\]
Check the value \( 7.175 \) of @exm-glh-longley from \( t_1=-0.921 \),
\( t_2=-1.222 \) and \( \rho=-0.840 \), and find the value of \( F \) had the correlation been \( +0.840 \).
:::

::: {.solution}
Let \( \bm D=\diag(\text{se}_1,\text{se}_2) \) and \( \bm R \) the correlation matrix of the estimates, so that the
estimated covariance is \( \bm D\bm R\bm D \). With \( \bLambda \) selecting the two coefficients,
\( 2F=\hat{\bbeta}_J\T(\bm D\bm R\bm D)^{-1}\hat{\bbeta}_J=\bm t\T\bm R^{-1}\bm t \), where
\( \bm t=\bm D^{-1}\hat{\bbeta}_J \), and \( \bm R^{-1}=(1-\rho^2)^{-1}\bigl(\begin{smallmatrix}1&-\rho\\-\rho&1\end{smallmatrix}\bigr) \).
Numerically, \( t_1^2+t_2^2=2.341 \), \( -2\rho t_1t_2=1.890 \) and \( 2(1-\rho^2)=0.589 \), so
\( F=7.175 \), as in the example, up to the rounding of these three terms. With \( \rho=+0.840 \) the middle term
changes sign and \( F=0.765 \), about \( (2.341-1.890)/0.589 \): the same two
\( t \) statistics would then be entirely unremarkable.
:::

### C. Going deeper

::: {#exr-glh-scheffe-preview}
[C1]

Show that \( \max\,T(\blambda)^2 \) over all nonzero estimable \( \blambda \), for the hypotheses
\( \blambda\T\bbeta=0 \), equals \( rF_* \), where \( F_* \) is the statistic for the reduced model \( \E(\Y)=\bzero \). Deduce that
\( \Pr\{T(\blambda)^2\le rF_\alpha(r,n-r)\text{ for all estimable }\blambda\}=1-\alpha \) when \( \bbeta=\bzero \). What
changes if \( \bbeta\ne\bzero \) and \( T(\blambda) \) is centred at \( \blambda\T\bbeta \)?
([Chapter 13](../ch13-multiplicity/index.html) builds simultaneous inference on this.)
:::

::: {#exr-glh-one-sided-lrt}
[C2]

Find the likelihood ratio test of \( \blambda\T\bbeta\le d \) against \( \blambda\T\bbeta>d \) for estimable \( \blambda \), and show
that it rejects for large values of \( T \) in @eq-glh-t. (Maximize the likelihood under the inequality
constraint: if \( \blambda\T\hbeta\le d \) the unrestricted maximum is feasible; otherwise the maximum is on the
boundary.)
:::

::: {.solution}
If \( \blambda\T\hbeta\le d \), the unrestricted maximizer satisfies the constraint and \( \Lambda=1 \). Otherwise the
log-likelihood, maximized over \( \sigma^2 \), is a decreasing function of \( \norm{\y-\X\bb}^2 \), which is a convex quadratic
in \( \bb \) with minimum at \( \hbeta \). Its minimum over the half-space \( \blambda\T\bb\le d \) therefore lies on the
boundary \( \blambda\T\bb=d \), where it equals \( \text{SSE}+(\blambda\T\hbeta-d)^2/\blambda\T\G\blambda \) (@eq-opt-constrained-sse).
So \( \Lambda=\{1+T^2/(n-r)\}^{-n/2} \) when \( T>0 \) and \( \Lambda=1 \) when \( T\le0 \). This is a nonincreasing function of
\( T \), strictly decreasing for \( T>0 \), so \( \Lambda\le k<1 \) iff \( T\ge c \) for some \( c>0 \), which is the one-sided \( t \)
test.
:::
