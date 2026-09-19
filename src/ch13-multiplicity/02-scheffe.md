# Scheffé's projection method

Scheffé's method answers the hardest version of the problem: the family is every estimable
function in a subspace, and the analyst may choose which ones to look at after seeing the data. The
intervals have simultaneous coverage exactly \( 1-\alpha \), because the largest standardized error in
the family is the length of a projection, and Scheffé's method *is* the \( F \) test, read one direction
at a time.

Throughout this section the model is the normal linear model @eq-opt-normal-model,
\( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \), with \( r=\rank(\X)<n \). Write \( \bmu=\X\bbeta \),
\( \nu=n-r \) for the error degrees of freedom, \( s^2=\text{SSE}/\nu \),
\( \G=(\X\T\X)\ginv \) for any generalized inverse, and \( F_\alpha(q,\nu) \) for the upper
\( \alpha \) point of \( F(q,\nu) \). For an estimable \( \blambda\T\bbeta \) the standard error is
\( \text{se}(\blambda\T\hbeta)=s\sqrt{\blambda\T\G\blambda} \) (@cor-opt-t).

## Linear functions of the mean as vectors

An estimable function is a linear function of the mean vector, and a linear function of
\( \bmu\in\C(\X) \) can be represented by a vector of \( \C(\X) \). If \( \blambda\T\bbeta \) is
estimable, then \( \blambda=\X\T\bm\rho \) for some vector
\( \bm\rho\in\Real^n \) (@thm-est-characterization), and with \( \bv=\M\bm\rho \),
\[
\blambda\T\bbeta=\bm\rho\T\X\bbeta=\bv\T\bmu,\qquad
\blambda\T\hbeta=\bm\rho\T\M\Y=\bv\T\Y,\qquad
\blambda\T\G\blambda=\bm\rho\T\X\G\X\T\bm\rho=\bm\rho\T\M\bm\rho=\norm{\bv}^2 .
\]{#eq-mc-function-vector}

The second identity is @thm-proj-invariant-functions, and the third uses \( \M=\X\G\X\T \) (@thm-proj-M-formula).
The vector \( \bv \) does not depend on the choice of \( \bm\rho \), since two
choices differ by some \( \bm\delta \) with \( \X\T\bm\delta=\bzero \), hence \( \M\bm\delta=\bzero \). So the
estimate is the inner product of \( \bv \) with the data, the estimand its inner product with the
mean, and the standard error is \( s\norm{\bv} \).

A *family* of estimable functions closed under linear combination is then a subspace
\( \mathcal V\subseteq\C(\X) \). It is usually described in one of two ways. If the family is
\( \{\blambda\T\bbeta:\blambda\in\C(\bLambda)\} \) for a \( p\times q \) matrix \( \bLambda=\X\T\bT \), then \( \bv=\M\bT\bm a \)
and \( \mathcal V=\C(\M\bT) \), the test space of \( \bLambda\T\bbeta=\bzero \) (@thm-proj-constraint-space); if
\( \bLambda \) has rank \( q \), so does \( \M\bT \), because \( \X\T\M\bT=\bLambda \). If instead the family is the set of
functions about which a reduced model \( \C(\X_0)\subseteq\C(\X) \) makes a claim, those vanishing on
\( \C(\X_0) \), then \( \mathcal V=\C(\M-\Mo) \), of dimension \( r-r_0 \) (@thm-proj-nested).

## The theorem

::: {#lem-mc-projection-max}
[The largest standardized component]

Let \( \mathcal V \) be a subspace of \( \Real^n \) with orthogonal projection \( \bP_{\mathcal V} \), and let
\( \bu\in\Real^n \). Then
\[
\max_{\bv\in\mathcal V,\ \bv\ne\bzero}\frac{(\bv\T\bu)^2}{\norm{\bv}^2}=\norm{\bP_{\mathcal V}\bu}^2 ,
\]
and if \( \bP_{\mathcal V}\bu\ne\bzero \) the maximum is attained exactly at the nonzero multiples of
\( \bP_{\mathcal V}\bu \).
:::

::: {.proof}
For \( \bv\in\mathcal V \), \( \bv=\bP_{\mathcal V}\bv \), and since \( \bP_{\mathcal V} \) is symmetric,
\( \bv\T\bu=\bv\T\bP_{\mathcal V}\bu \). By Cauchy–Schwarz (@prp-mat-cauchy-schwarz),
\( (\bv\T\bP_{\mathcal V}\bu)^2\le\norm{\bv}^2\norm{\bP_{\mathcal V}\bu}^2 \), with equality iff \( \bv \) is
a multiple of \( \bP_{\mathcal V}\bu \). That multiple lies in \( \mathcal V \), so the bound is attained.
If \( \bP_{\mathcal V}\bu=\bzero \) both sides are zero.
:::

The lemma is @exr-ss-max-single of [Chapter 9](../ch09-sums-of-squares/index.html), where it described the sum of squares of a row of an analysis of
variance table as the largest one-degree-of-freedom sum of squares inside it. Applied to
\( \bu=\Y-\bmu \), it describes the worst standardized error in a family.

::: {#thm-mc-scheffe}
[Scheffé's simultaneous intervals]

Assume @eq-opt-normal-model with \( r<n \), and let \( 0<\alpha<1 \).

::: {.enumerate options="label=(\alph*)"}
1. Let \( \mathcal V\subseteq\C(\X) \) be a subspace of dimension \( q\ge1 \). Then
   \[
   \max_{\bv\in\mathcal V,\ \bv\ne\bzero}\frac{\bigl(\bv\T\Y-\bv\T\bmu\bigr)^2}{q\,s^2\norm{\bv}^2}
   =\frac{\norm{\bP_{\mathcal V}(\Y-\bmu)}^2/q}{s^2}\sim F(q,\nu),
   \]
   and consequently
   \[
   \Pr\Bigl(\lvert\bv\T(\Y-\bmu)\rvert\le\sqrt{qF_\alpha(q,\nu)}\,s\norm{\bv}\ \ \forall\,\bv\in\mathcal V\Bigr)=1-\alpha .
   \]

2. Let \( \bLambda \) be \( p\times q \) of rank \( q \) with \( \C(\bLambda)\subseteq\C(\X\T) \). With probability
   exactly \( 1-\alpha \), for every \( \blambda\in\C(\bLambda) \) simultaneously,
   \[
   \blambda\T\bbeta\in\Bigl[\blambda\T\hbeta-\sqrt{qF_\alpha(q,\nu)}\;\text{se}(\blambda\T\hbeta),\ \
   \blambda\T\hbeta+\sqrt{qF_\alpha(q,\nu)}\;\text{se}(\blambda\T\hbeta)\Bigr].
   \]{#eq-mc-scheffe-interval}

   Moreover, with \( \bm\theta=\bLambda\T\bbeta \), \( \hat{\bm\theta}=\bLambda\T\hbeta \) and
   \( \W=\bLambda\T\G\bLambda \),
   \[
   \max_{\blambda\in\C(\bLambda),\ \blambda\ne\bzero}\frac{(\blambda\T\hbeta-\blambda\T\bbeta)^2}{\text{se}(\blambda\T\hbeta)^2}
   =\frac{(\hat{\bm\theta}-\bm\theta)\T\W^{-1}(\hat{\bm\theta}-\bm\theta)}{s^2},
   \]{#eq-mc-scheffe-max}

   attained at \( \blambda=\bLambda\W^{-1}(\hat{\bm\theta}-\bm\theta) \).
:::

:::

::: {.proof}
*(a)* The first equality is @lem-mc-projection-max with \( \bu=\Y-\bmu \). The map \( \bv\mapsto\X\T\bv \) is
one to one on \( \C(\X) \) (Step 1 of the proof of @thm-ci-working-hotelling), so
\( \mathcal L=\{\X\T\bv:\bv\in\mathcal V\} \) is a subspace of \( \C(\X\T) \) of dimension \( q \). For
\( \x=\X\T\bv \), @eq-mc-function-vector with \( \bm\rho=\bv \) gives
\( \x\T\hbeta-\x\T\bbeta=\bv\T(\Y-\bmu) \) and \( \x\T\G\x=\norm{\bv}^2 \). So the event in the display is
the event @eq-ci-wh-band with \( k=q \), and it has probability exactly \( 1-\alpha \) by
@thm-ci-working-hotelling. That event is also the event that the maximum is at most
\( F_\alpha(q,\nu) \). Since this holds for every \( \alpha\in(0,1) \), the maximum has the
distribution function of \( F(q,\nu) \).

*(b)* Write \( \bLambda=\X\T\bT \) and \( \mathcal V=\C(\M\bT) \), of dimension \( q \). Every
\( \blambda\in\C(\bLambda) \) is \( \bLambda\bm a=\X\T(\bT\bm a) \), and by @eq-mc-function-vector
\( \blambda\T\hbeta-\blambda\T\bbeta=\bv\T\Y-\bv\T\bmu \) and \( \text{se}(\blambda\T\hbeta)=s\norm{\bv} \)
with \( \bv=\M\bT\bm a\in\mathcal V \). As \( \bm a \) ranges over \( \Real^q \), \( \bv \) ranges over all of
\( \mathcal V \). So the event in (b) is the event in (a).
Finally, @eq-mc-scheffe-max is @prp-glh-max-t for the hypothesis \( \bLambda\T\bbeta=\bm\theta \), which is
testable because \( \bm\theta \) is the true value: its statistic @eq-glh-general-F is
\( (\hat{\bm\theta}-\bm\theta)\T\W^{-1}(\hat{\bm\theta}-\bm\theta)/(q\,s^2) \), and each \( T(\bm a)^2 \) is the ratio on the
left of @eq-mc-scheffe-max at \( \blambda=\bLambda\bm a \), so the maximum is attained at
\( \bm a\propto\W^{-1}(\hat{\bm\theta}-\bm\theta) \).
:::

The coverage is exact, not a bound: some function in the family, depending on the data, always lies
on the edge of the event. Only the subspace matters, so dependent columns of \( \bLambda \) may be
dropped first. For \( q=1 \) the multiplier is \( t_{\nu,\alpha/2} \), since \( F_\alpha(1,\nu)=t_{\nu,\alpha/2}^2 \),
and the interval is that of @thm-ci-estimable-interval.

Part (b) is the observation-space form of a fact met twice already. [Section 11.6](../ch11-general-linear-hypothesis/06-coefficients.html)
proved @eq-mc-scheffe-max, with \( \bm d \) in place of \( \bm\theta \), as the statement that a group \( F \)
statistic is the largest single \( t \) statistic in the group, and [Section 12.2](../ch12-intervals-and-bands/02-ellipsoids.html)
showed that the intervals @eq-mc-scheffe-interval are the *shadows* of the confidence ellipsoid for
\( \bm\theta=\bLambda\T\bbeta \) (@thm-ci-ellipsoid). The ellipsoid covers \( \bm\theta \) iff every shadow covers
the corresponding \( \bm a\T\bm\theta \), which is why an infinite family costs no more than one ellipsoid.
The coverage statement itself is the Working–Hotelling band of @thm-ci-working-hotelling for the
subspace \( \mathcal L=\C(\bLambda) \): a band for a regression surface and a family of intervals for
estimable functions are one object, written in the coordinates of \( \C(\X\T) \) or of \( \C(\X) \).
[Figure 13.2.1](#fig-mc-shadow) extends Figure 12.2.1, which drew the same ellipse for two slopes of the
crime regression, by adding the Bonferroni rectangle and the tangent lines that bound the Scheffé
interval for a combination \( a_1\beta_1+a_2\beta_2 \).

::: {when-format="html"}
![**Figure 13.2.1.** The coefficients of poverty and single parenthood in the crime regression. The
Scheffé intervals (dashed box, multiplier \( 2.530 \)) are the shadows of the 95% confidence ellipse;
the Bonferroni rectangle (multiplier \( 2.317 \)) covers only the two coefficients. The tangent lines
bound the Scheffé interval for the sum of the slopes. The ellipse is that of Figure 12.2.1.](scheffe_shadow.svg){#fig-mc-shadow width=100%}
:::

::: {when-format="pdf"}
![The coefficients of poverty and single parenthood in the crime regression. The
Scheffé intervals (dashed box, multiplier \( 2.530 \)) are the shadows of the 95% confidence ellipse;
the Bonferroni rectangle (multiplier \( 2.317 \)) covers only the two coefficients. The tangent lines
bound the Scheffé interval for the sum of the slopes. The ellipse is that of Figure 12.2.1.](scheffe_shadow.pdf){width=100%}
:::

::: {#exm-mc-scheffe-crime}
[Three slopes of a crime regression]

Return to the regression of the murder rate on poverty, single parenthood and urbanization
for the \( 50 \) states of @exm-proj-fwl-crime, with \( \nu=46 \) error degrees of freedom. The
poverty coefficient is \( 0.2538 \) with standard error \( 0.0934 \), so \( t=2.718 \). Suppose the
family of interest is every linear combination of the three slopes, \( q=3 \). The \( 95\% \)
intervals for the poverty coefficient are:

| method | multiplier | interval |
|---|---|---|
| unadjusted \( t \) | \( 2.013 \) | \( [0.066,\ 0.442] \) |
| Bonferroni, three slopes | \( 2.485 \) | \( [0.022,\ 0.486] \) |
| Scheffé, all combinations | \( 2.902 \) | \( [-0.017,\ 0.525] \) |

The poverty coefficient is significant on its own and among three planned comparisons, but not
once every combination of the slopes is admitted. The \( F \) statistic for the hypothesis that all
three slopes are zero is \( 27.49 \), and the largest \( \lvert t\rvert \) over all combinations is
\( 9.081=\sqrt{3\times27.49} \), in agreement with @eq-mc-scheffe-max.
:::

```{.python .run #cell-scheffe-crime-fit}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")

y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
n, p = X.shape
nu = n - p                                            # error degrees of freedom
XtX_inv = np.linalg.inv(X.T @ X)
beta_hat = XtX_inv @ X.T @ y
s2 = np.sum((y - X @ beta_hat) ** 2) / nu
se = np.sqrt(s2 * np.diag(XtX_inv))
```

```{.python .run #cell-scheffe-crime-multipliers}
alpha, slopes = 0.05, [1, 2, 3]                       # poverty, single, urban
k = len(slopes)
mult = {
    "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
    "Bonferroni": stats.t.ppf(1 - alpha / (2 * k), nu),
    "Scheffe": np.sqrt(k * stats.f.ppf(1 - alpha, k, nu)),   # q = 3 slopes
}
for name, c in mult.items():
    lo, hi = beta_hat[slopes] - c * se[slopes], beta_hat[slopes] + c * se[slopes]
    print(f"{name:10s} c = {c:.3f}  poverty [{lo[0]:.3f}, {hi[0]:.3f}]")
```

```{.python .run #cell-scheffe-crime-snoop}
L = np.eye(p)[:, slopes]                              # Lambda: the three slopes
W = L.T @ XtX_inv @ L
d = L.T @ beta_hat
F = d @ np.linalg.solve(W, d) / (k * s2)              # F test that all slopes are 0
a_star = np.linalg.solve(W, d)                        # the most significant combination
t_star = (a_star @ d) / np.sqrt(s2 * a_star @ W @ a_star)
print(f"F = {F:.2f};  largest |t| over all combinations = {t_star:.3f} = sqrt(3F)")
```

## Scheffé's method is the F test

Let \( \bm d\in\Real^q \) and consider the hypothesis \( \bLambda\T\bbeta=\bm d \), with \( \bLambda \) as in @thm-mc-scheffe(b).
Because \( \bLambda\T \) has rank \( q \), the hypothesis is consistent for every
\( \bm d \). By @thm-glh-general-f its \( F \) statistic is
\[
F=\frac{(\hat{\bm\theta}-\bm d)\T\W^{-1}(\hat{\bm\theta}-\bm d)}{q\,s^2} .
\]
The hypothesis implies a single-degree-of-freedom statement for every
\( \blambda=\bLambda\bm a \) in the family, namely \( \blambda\T\bbeta=\bm a\T\bm d \). Call such a statement
**Scheffé-rejected** if its hypothesized value lies outside the Scheffé interval @eq-mc-scheffe-interval.

::: {#cor-mc-scheffe-f}
[The \( F \) test and Scheffé's intervals]

In the setting of @thm-mc-scheffe(b), for every \( \bm d\in\Real^q \):

::: {.enumerate options="label=(\alph*)"}
1. \( \displaystyle\max_{\bm a\ne\bzero}\frac{(\bm a\T\hat{\bm\theta}-\bm a\T\bm d)^2}{s^2\,\bm a\T\W\bm a}=qF \),
   attained at \( \bm a=\W^{-1}(\hat{\bm\theta}-\bm d) \);

2. the level-\( \alpha \) \( F \) test rejects \( \bLambda\T\bbeta=\bm d \) iff at least one statement
   \( \blambda\T\bbeta=\bm a\T\bm d \), \( \blambda=\bLambda\bm a \), is Scheffé-rejected;

3. if \( \bLambda\T\bbeta=\bm d \) is true, the probability that any of these statements is
   Scheffé-rejected is exactly \( \alpha \). For any \( \bbeta \), the probability of Scheffé-rejecting at
   least one *true* statement \( \blambda\T\bbeta=c \), \( \blambda\in\C(\bLambda) \), is at most \( \alpha \).
:::

:::

::: {.proof}
*(a)* For \( \blambda=\bLambda\bm a \), the ratio is \( (\bm a\T\bm w)^2/(s^2\bm a\T\W\bm a) \) with
\( \bm w=\hat{\bm\theta}-\bm d \), and its maximum is \( \bm w\T\W^{-1}\bm w/s^2=qF \), attained at \( \bm a=\W^{-1}\bm w \), by
@cor-mat-generalized-rayleigh(c) with \( \W \) positive definite (@thm-ss-hypothesis(a)); this is the maximum
identity of [Section 11.6](../ch11-general-linear-hypothesis/06-coefficients.html). *(b)* The statement for \( \bm a \) is Scheffé-rejected iff its squared
\( t \) statistic exceeds \( qF_\alpha(q,\nu) \), so some statement is rejected iff the maximum in (a)
exceeds \( qF_\alpha(q,\nu) \), that is, iff \( F>F_\alpha(q,\nu) \). *(c)* The first claim follows from (b)
and the exact level of the \( F \) test (@thm-glh-general-f). The second is @prp-mc-intervals-to-tests(a)
applied to the simultaneous intervals of @thm-mc-scheffe.
:::

When the \( F \) test rejects, *some* function in the family is Scheffé-significant, though perhaps not
one of interest; conversely, no function is Scheffé-significant unless the \( F \) test rejects, so no
separate protecting \( F \) test is needed. For a reduced model the same holds with
\( \mathcal V=\C(\M-\Mo) \): by @lem-mc-projection-max the \( F \) statistic of @thm-glh-f-test is the
maximum over \( \bv\in\mathcal V \) of \( (\bv\T\Y)^2/\bigl((r-r_0)s^2\norm{\bv}^2\bigr) \) (@exr-mc-reduced-family). Taking \( \C(\bLambda)=\C(\X\T) \), so that \( q=r \) and the family contains the mean response \( \x\T\bbeta \)
at every \( \x\in\C(\X\T) \), gives the band \( \x\T\hbeta\pm\sqrt{rF_\alpha(r,\nu)}\,s\sqrt{\x\T\G\x} \), the
Working–Hotelling band of @thm-ci-working-hotelling.

## Smaller families give shorter intervals

The Scheffé multiplier depends on the family only through its dimension \( q \), and
\( qF_\alpha(q,\nu) \) is strictly increasing in \( q \) for fixed \( \nu \) and \( \alpha \), as a \( \chi^2 \)
coupling shows (an exercise of [Section 11.6](../ch11-general-linear-hypothesis/06-coefficients.html)).
So one should use the smallest subspace containing every function of interest, whose dimension can be
smaller than the length of a list spanning it (@exr-mc-span-three). For a few planned functions,
Bonferroni's method ([Section 13.4](04-bonferroni-holm.html)) is usually shorter still, as in
@exm-mc-scheffe-crime; Scheffé's multiplier is earned when functions are picked after the estimates are seen.

## All contrasts in a one-way layout

In the one-way model \( \E(y_{kj})=\mu+\alpha_k \) with \( g \) levels, \( n_k\ge1 \) observations at level
\( k \) and \( n=\sum_kn_k \), the estimable functions of the effects alone are the contrasts
\( \psi=\sum_kc_k\mu_k \) with \( \sum_kc_k=0 \) (@thm-est-oneway-contrasts). They form a space of dimension
\( g-1 \), the test space of the hypothesis of equal means. The estimate is
\( \hat{\psi}=\sum_kc_k\bar{y}_k \) with standard error \( s\bigl(\sum_kc_k^2/n_k\bigr)^{1/2} \), and
\( \nu=n-g \). Scheffé's intervals are
\[
\sum_kc_k\bar{y}_k\;\pm\;\sqrt{(g-1)F_\alpha(g-1,n-g)}\;s\Bigl(\sum_k\frac{c_k^2}{n_k}\Bigr)^{1/2}
\qquad\text{for every contrast }\bm c,
\]{#eq-mc-scheffe-contrasts}

with simultaneous coverage exactly \( 1-\alpha \), balanced or not. The most significant contrast
has \( c_k=n_k(\bar{y}_k-\bar{y}) \) (@exr-est-max-contrast), and its squared \( t \) statistic is
\( (g-1)F \), where \( F \) is the one-way \( F \) statistic.

::: {#exm-mc-scheffe-education}
[Education and the perceived position of a candidate]

The 1996 American National Election Study asked respondents to place Bill Clinton on a
seven-point scale from extremely liberal (1) to extremely conservative (7). Group the
\( n=944 \) respondents by education, coded from 1 (least) to 7 (most). The group sizes are
\( 13,\ 52,\ 248,\ 187,\ 90,\ 227,\ 127 \), and the mean placements are
\[
4.08,\ 3.58,\ 3.21,\ 2.95,\ 2.79,\ 2.66,\ 2.65 .
\]
More educated respondents place Clinton further to the left. The pooled standard deviation is
\( s=1.354 \) on \( \nu=937 \) degrees of freedom, and the one-way \( F \) statistic is \( 7.88 \), with
\( p \)-value \( 2.7\times 10^{-8} \). The between-groups sum of squares is \( 86.72 \).

Having seen the means, one might ask how the two least educated groups compare with the two
most educated, the contrast \( \psi=\tfrac12(\mu_1+\mu_2)-\tfrac12(\mu_6+\mu_7) \). Its estimate is
\( 1.176 \) with standard error \( 0.223 \), so \( t=5.27 \). Because the contrast was suggested by the
data, the unadjusted \( t \) interval is not legitimate. The Scheffé interval @eq-mc-scheffe-contrasts,
with multiplier \( \sqrt{6F_{0.05}(6,937)}=3.557 \), is
\( [0.383,\ 1.969] \), and it is legitimate however the contrast was found. The most significant
contrast of all has \( t=6.88=\sqrt{6\times7.88} \), and its sum of squares is the whole
between-groups sum of squares.
:::

```{.python .run #cell-education-data}
import numpy as np
import statsmodels.api as sm
from scipy import stats
anes = sm.datasets.anes96.load_pandas().data

y = anes["ClinLR"].to_numpy()
level = anes["educ"].to_numpy().astype(int) - 1       # 0, ..., 6
g = 7
n_k = np.bincount(level).astype(float)
means = np.bincount(level, weights=y) / n_k
n = len(y)
nu = n - g                                            # error degrees of freedom
s2 = np.sum((y - means[level]) ** 2) / nu             # pooled variance estimate
s = np.sqrt(s2)
between = np.sum(n_k * (means - y.mean()) ** 2)
F = between / (g - 1) / s2
print("group sizes", n_k.astype(int))
print("means      ", means.round(3))
print(f"s = {s:.4f} on {nu} df;  F = {F:.3f},  p = {stats.f.sf(F, g - 1, nu):.2g}")
```

```{.python .run #cell-education-scheffe}
alpha = 0.05
c_sch = np.sqrt((g - 1) * stats.f.ppf(1 - alpha, g - 1, nu))   # Scheffe multiplier, q = g - 1

def contrast(c):
    """Estimate, standard error and sum of squares of the contrast sum_k c_k mu_k."""
    est = c @ means
    v = np.sum(c ** 2 / n_k)
    return est, s * np.sqrt(v), est ** 2 / v

c_max = n_k * (means - y.mean())                      # the most significant contrast
c_lohi = np.array([0.5, 0.5, 0, 0, 0, -0.5, -0.5])    # levels 1-2 against levels 6-7
for name, c in (("maximal", c_max), ("1-2 vs 6-7", c_lohi)):
    est, se, ss = contrast(c)
    print(f"{name:11s} estimate {est:8.3f}  t = {est / se:6.2f}  SS = {ss:8.2f}"
          f"  Scheffe interval [{est - c_sch * se:.3f}, {est + c_sch * se:.3f}]")
print(f"between-groups SS = {between:.2f};  Scheffe multiplier = {c_sch:.3f}")
```

## Exercises

### A. Check your understanding

::: {#exr-mc-scheffe-one}
[A1]

Show that for \( q=1 \) Scheffé's interval @eq-mc-scheffe-interval is the \( t \) interval of @thm-ci-estimable-interval,
and explain why no multiplicity adjustment appears.
:::

::: {#exr-mc-span-three}
[A2]

In @exm-mc-scheffe-education, an analyst decides before seeing the data to examine the
three contrasts \( \tfrac12(\mu_1+\mu_2)-\tfrac12(\mu_6+\mu_7) \), \( \mu_1-\mu_7 \) and \( \mu_2-\mu_6 \).
Show that they span a space of dimension \( 2 \), and state the Scheffé multiplier for
this family.
:::

### B. Practice

::: {#exr-mc-reduced-family}
[B1]

Let \( \C(\X_0)\subseteq\C(\X) \) with ranks \( r_0<r \), and \( \mathcal V=\C(\M-\Mo) \). Show that a vector
\( \bv\in\C(\X) \) lies in \( \mathcal V \) iff \( \bv\T\bmu=0 \) for every \( \bmu\in\C(\X_0) \). Show that the maximum
over \( \bv\in\mathcal V \) of the squared \( t \) statistic for \( \bv\T\bmu=0 \) equals \( (r-r_0)F \), where \( F \) is
the statistic of @thm-glh-f-test.
:::

### C. Going deeper

::: {#exr-mc-scheffe-shortest}
[C1]

Show that Scheffé's multiplier cannot be reduced: if \( c<\sqrt{qF_\alpha(q,\nu)} \), the intervals
\( \blambda\T\hbeta\pm c\,\text{se}(\blambda\T\hbeta) \), \( \blambda\in\C(\bLambda) \), have simultaneous coverage less
than \( 1-\alpha \). So among intervals of the form estimate \( \pm \) constant \( \times \) standard error,
valid for the whole family, Scheffé's are the shortest.
:::

::: {.solution}
By @eq-mc-scheffe-max, all intervals cover iff the maximal squared \( t \) statistic, which has
\( q \) times an \( F(q,\nu) \) distribution, is at most \( c^2 \). So the coverage is
\( \Pr\{F(q,\nu)\le c^2/q\} \), which is less than \( 1-\alpha \) when \( c^2/q<F_\alpha(q,\nu) \), because the
\( F \) distribution has a positive density.
:::
