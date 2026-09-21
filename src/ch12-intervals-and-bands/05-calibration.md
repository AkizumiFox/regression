# Inverse prediction and calibration

Prediction runs from \( x \) to \( y \). Calibration runs the other way: an instrument is calibrated on specimens
of known value, a new specimen is measured, and its value is wanted. The estimate reads the fitted line
backwards; the interval reads the prediction band backwards, and it is not always an interval.

## Ratios of estimable functions

Inverse prediction belongs to a wider problem: a confidence set for a *ratio* of linear functions of the coefficients,
such as the zero crossing \( -\beta_0/\beta_1 \) of a line or the turning point \( -\beta_1/(2\beta_2) \) of a quadratic. A delta-method
interval ignores that the denominator is uncertain and may be near zero. Fieller's approach tests, for each candidate value, a
hypothesis that is linear in the coefficients.

::: {#prp-ci-fieller}
[Fieller's confidence set for a ratio]

Let \( \blambda_1,\blambda_2\in\C(\X\T) \) be linearly independent, put \( a=\blambda_1\T\bbeta \), \( b=\blambda_2\T\bbeta \), their estimates
\( \hat{a}=\blambda_1\T\hbeta \), \( \hat{b}=\blambda_2\T\hbeta \), and \( v_{ij}=\blambda_i\T\G\blambda_j \), so that \( \mathbf{V}=(v_{ij}) \) is positive definite.
Suppose \( b\ne0 \) and let \( \theta=a/b \). With \( t=t_{n-r,\alpha/2} \), let
\[
\begin{aligned}
\mathcal F(\Y)=\bigl\{\theta_0\in\Real:\ &(\hat{a}-\theta_0\hat{b})^2\\
&\le t^2s^2\,(v_{11}-2\theta_0v_{12}+\theta_0^2v_{22})\bigr\}.
\end{aligned}
\]{#eq-ci-fieller}

::: {.enumerate options="label=(\alph*)"}
1. \( \mathcal F(\Y) \) is an exact \( 1-\alpha \) confidence set for \( \theta \), at every parameter value with \( b\ne0 \).

2. Let \( Q=(\hat{a},\hat{b})\mathbf{V}^{-1}(\hat{a},\hat{b})\T \). Outside an event of probability zero, \( \mathcal F(\Y) \) is

   - a bounded closed interval containing \( \hat{\theta}=\hat{a}/\hat{b} \), if \( \hat{b}^2>t^2s^2v_{22} \);
   - the union of two disjoint closed rays, if \( \hat{b}^2<t^2s^2v_{22} \) and \( Q>t^2s^2 \);
   - the whole real line, if \( \hat{b}^2<t^2s^2v_{22} \) and \( Q<t^2s^2 \).
:::
:::

::: {.proof}
(a) Fix \( \theta_0 \). The vector \( \blambda_1-\theta_0\blambda_2 \) lies in \( \C(\X\T) \) and is not zero, so
\( (\blambda_1-\theta_0\blambda_2)\T\bbeta=a-\theta_0b \) is estimable, with estimate \( \hat{a}-\theta_0\hat{b} \) and
\( (\blambda_1-\theta_0\blambda_2)\T\G(\blambda_1-\theta_0\blambda_2)=v_{11}-2\theta_0v_{12}+\theta_0^2v_{22} \). The inequality in
@eq-ci-fieller says that the \( t \) interval of @thm-ci-estimable-interval for \( a-\theta_0b \) contains \( 0 \). At \( \theta_0=\theta \),
\( a-\theta b=0 \), so the probability of this event is exactly \( 1-\alpha \).

(b) The defining inequality is \( q(\theta_0)=A\theta_0^2-2B\theta_0+K\le0 \), with
\[
A=\hat{b}^2-t^2s^2v_{22},\qquad B=\hat{a}\hat{b}-t^2s^2v_{12},\qquad K=\hat{a}^2-t^2s^2v_{11}.
\]
Expanding and using \( (\hat{a},\hat{b})\mathbf{V}^{-1}(\hat{a},\hat{b})\T=(v_{22}\hat{a}^2-2v_{12}\hat{a}\hat{b}+v_{11}\hat{b}^2)/\det\mathbf{V} \),
\[
\begin{aligned}
B^2-AK&=t^2s^2\bigl(v_{22}\hat{a}^2-2v_{12}\hat{a}\hat{b}+v_{11}\hat{b}^2-t^2s^2\det\mathbf{V}\bigr)\\
&=t^2s^2\det(\mathbf{V})\,(Q-t^2s^2).
\end{aligned}
\]
If \( A>0 \), then \( q(\hat{\theta})=-t^2s^2(v_{11}-2\hat{\theta} v_{12}+\hat{\theta}^2v_{22})<0 \), so the upward parabola \( q \) has two real
roots, and \( \{q\le0\} \) is the closed interval between them, which contains \( \hat{\theta} \). If \( A<0 \), the parabola opens downward. When
\( B^2-AK<0 \) it is negative everywhere and \( \mathcal F(\Y)=\Real \). When \( B^2-AK>0 \) it is nonpositive exactly outside the open interval between its
roots. The events \( A=0 \) and \( Q=t^2s^2 \) have probability zero, because \( \hat{b} \) and \( Q/s^2 \) have continuous distributions.
:::

The first case says that the \( t \) test rejects \( b=0 \): the denominator is away from zero and the ratio is confined. Otherwise
very large \( \lvert\theta\rvert \) cannot be excluded, and \( Q/(2s^2) \), the \( F \) statistic for \( a=b=0 \), decides between two rays and
the whole line.
@exr-ci-fieller-geometry gives a picture: \( \mathcal F(\Y) \) is the set of slopes of lines through the origin that meet a certain
confidence ellipse for \( (a,b) \).

## The calibration problem

A calibration experiment measures \( n \) standards with known values \( x_1,\dots,x_n \), not all equal, and records readings
\[
Y_i=\beta_0+\beta_1x_i+\varepsilon_i,\qquad i=1,\dots,n .
\]
A new specimen with unknown value \( x_0 \) is then measured \( m\ge1 \) times,
\( Y_{0j}=\beta_0+\beta_1x_0+\varepsilon_{0j} \), and all \( n+m \) errors are independent \( \Normal(0,\sigma^2) \). Let \( \bar{Y}_0 \) be the mean of
the new readings, and let \( \hat{\beta}_0 \), \( \hat{\beta}_1 \) and \( s^2 \) (with \( n-2 \) degrees of freedom) come from the standards alone. The
**classical estimator** reads the fitted line backwards:
\[
\hat{x}_0=\frac{\bar{Y}_0-\hat{\beta}_0}{\hat{\beta}_1}=\bar{x}+\frac{\bar{Y}_0-\bar{Y}}{\hat{\beta}_1}.
\]
It is the maximum likelihood estimator (@exr-ci-calibration-mle) and, being a ratio, inherits a pathology of ratios.

::: {#prp-ci-ratio-no-mean}
[The classical estimator has no mean]

If \( \Pr(\bar{Y}_0\ne\bar{Y})>0 \), which always holds here, then \( \E\lvert\hat{x}_0\rvert=\infty \).
:::

::: {.proof}
\( \bar{Y} \) and \( \hat{\beta}_1 \) are jointly normal and uncorrelated (@cor-lm-simple-moments), hence independent (@thm-mvn-independence),
and \( \bar{Y}_0 \) is independent of both. So \( g=\bar{Y}_0-\bar{Y} \) is independent of \( \hat{\beta}_1 \) and
\( \E\lvert\hat{x}_0-\bar{x}\rvert=\E\lvert g\rvert\,\E\lvert1/\hat{\beta}_1\rvert \), with \( \E\lvert g\rvert>0 \). The density \( f \) of
\( \hat{\beta}_1 \) is continuous and positive at \( 0 \), so \( f\ge c>0 \) on some \( [-\epsilon,\epsilon] \) and
\( \E\lvert1/\hat{\beta}_1\rvert\ge c\int_{-\epsilon}^{\epsilon}\lvert u\rvert^{-1}\,du=\infty \). Finally
\( \lvert\hat{x}_0\rvert\ge\lvert\hat{x}_0-\bar{x}\rvert-\lvert\bar{x}\rvert \).
:::

So a delta-method "standard error" of \( \hat{x}_0 \) approximates a quantity that does not exist. That is harmless when
\( \hat{\beta}_1 \) is many standard errors from zero, and harmful exactly when the slope is poorly determined.

::: {#thm-ci-calibration}
[Calibration set]

In the calibration model, let \( t=t_{n-2,\alpha/2} \) and
\[
\begin{aligned}
\mathcal C=\Bigl\{x\in\Real:\ &\bigl(\bar{Y}_0-\hat{\beta}_0-\hat{\beta}_1x\bigr)^2\\
&\le t^2s^2\Bigl(\frac1m+\frac1n+\frac{(x-\bar{x})^2}{S_{xx}}\Bigr)\Bigr\}.
\end{aligned}
\]{#eq-ci-calibration}

::: {.enumerate options="label=(\alph*)"}
1. \( \mathcal C \) is an exact \( 1-\alpha \) confidence set for \( x_0 \), at every \( (\beta_0,\beta_1,\sigma^2,x_0) \), including \( \beta_1=0 \).

2. \( \mathcal C \) is the set of \( x \) at which \( \bar{Y}_0 \) lies inside the \( 1-\alpha \) prediction interval of
   @thm-ci-prediction-interval(c) for the mean of \( m \) new readings at \( x \).

3. Let \( T_1=\hat{\beta}_1\sqrt{S_{xx}}/s \) be the \( t \) statistic for the slope and \( g=\bar{Y}_0-\bar{Y} \). Outside an event of probability
   zero, \( \mathcal C \) is

   - a bounded interval containing \( \hat{x}_0 \), if \( T_1^2>t^2 \);
   - the union of two disjoint closed rays, if \( T_1^2<t^2 \) and \( T_1^2+g^2/\{s^2(1/m+1/n)\}>t^2 \);
   - the whole real line, if \( T_1^2+g^2/\{s^2(1/m+1/n)\}<t^2 \).
:::
:::

::: {.proof}
(a) and (b). For fixed \( x \), \( \bar{Y}_0-\hat{\beta}_0-\hat{\beta}_1x \) is the difference between the mean of the new readings and the
estimated mean response at \( x \). At \( x=x_0 \), @thm-ci-prediction-interval(c), with \( \x_0=(1,x_0)\T \) and
\( h_0=1/n+(x_0-\bar{x})^2/S_{xx} \) (@eq-lm-var-mean-response), shows that
\[
\frac{\bar{Y}_0-\hat{\beta}_0-\hat{\beta}_1x_0}{s\sqrt{1/m+1/n+(x_0-\bar{x})^2/S_{xx}}}\sim t(n-2),
\]
whatever the parameters are. This is a pivot, and \( x_0\in\mathcal C \) is the event that its absolute value is at most \( t \). Nothing
in the argument requires \( \beta_1\ne0 \). Part (b) restates the definition.

(c) Put \( d=x-\bar{x} \) and \( c=1/m+1/n \). Since \( \bar{Y}_0-\hat{\beta}_0-\hat{\beta}_1x=g-\hat{\beta}_1d \), the defining inequality is
\( Ad^2-2\hat{\beta}_1g\,d+(g^2-t^2s^2c)\le0 \) with \( A=\hat{\beta}_1^2-t^2s^2/S_{xx} \). Its discriminant is
\[
\begin{aligned}
\hat{\beta}_1^2g^2-A(g^2-t^2s^2c)&=t^2s^2\Bigl(\frac{g^2}{S_{xx}}+cA\Bigr)\\
&=\frac{t^2s^4c}{S_{xx}}\Bigl(T_1^2+\frac{g^2}{s^2c}-t^2\Bigr).
\end{aligned}
\]
The sign of \( A \) is that of \( T_1^2-t^2 \). If \( A>0 \), the quadratic is negative at \( d=g/\hat{\beta}_1 \), that is at \( x=\hat{x}_0 \), where it
equals \( -t^2s^2(c+d^2/S_{xx}) \); so it has two real roots and \( \mathcal C \) is the interval between them. If \( A<0 \), the argument of
@prp-ci-fieller(b) applies: two rays if the discriminant is positive and the whole line if it is negative. Equalities have probability zero.
:::

Part (b) is the practical recipe ([Figure 12.5.1](#fig-ci-calibration)): draw the prediction band for the mean of \( m \)
readings and the horizontal line at \( \bar{Y}_0 \), and read off where the line is inside the band. A steep calibration line crosses the
band in a short interval. A flat one lets both hyperbolic edges pass the horizontal line, giving two rays or everything. The set is
bounded iff the slope test rejects, as it does for any useful instrument.

::: {#prp-ci-calibration-unbounded}
[Unbounded calibration sets cannot be avoided]

Let \( \mathcal D \) be any confidence set for \( x_0 \) with coverage at least \( 1-\alpha>0 \) at every parameter value, including those with
\( \beta_1=0 \). Then at every parameter value with \( \beta_1=0 \), \( \Pr\{v\in\mathcal D\}\ge1-\alpha \) for every real \( v \), and if
\( \mathcal D \) is measurable in the sense of @prp-ci-nonestimable, its expected length is infinite.
:::

::: {.proof}
When \( \beta_1=0 \), the joint distribution of the standards and the new readings does not involve \( x_0 \). So for every \( v \), the
probability that \( v\in\mathcal D \) is the same as it would be if \( x_0 \) were \( v \), and that probability is at least \( 1-\alpha \). The rest
is the Fubini argument in the proof of @prp-ci-nonestimable.
:::

By continuity the same happens approximately for small \( \beta_1 \). A
procedure that always reports a bounded interval, like the delta method, pays with coverage below \( 1-\alpha \) somewhere. The Fieller
set's unbounded answers are the honest answer to a question the data cannot settle.

::: {#exm-ci-assay}
[Calibrating an assay]

In a simulated assay, nine standards of concentration \( 0,1,\dots,8 \) are each measured three times (\( n=27 \)). The
fitted line is \( 0.0338+0.1161\,x \), with \( s=0.0131 \) and slope \( t \) statistic \( 118.7 \). A new
specimen is read \( m=3 \) times, with mean reading \( 0.6521 \). The classical estimate is \( \hat{x}_0=5.327 \), and the
\( 95\% \) calibration interval is \( (5.184,\ 5.471) \) (\( t_{25,0.025}=2.060 \)). For an assay this precise the
delta-method interval, \( \hat{x}_0\pm t\,(s/\lvert\hat{\beta}_1\rvert)\sqrt{1/m+1/n+(\hat{x}_0-\bar{x})^2/S_{xx}} \), agrees to three
decimals: \( (5.184,\ 5.471) \).

Now consider a weak assay with the same design, in which the fitted slope is \( 0.0021 \) and its \( t \) statistic is only
\( 1.67 \). The mean reading of the standards is \( 0.0475 \). A specimen with mean reading \( 0.0525 \), close to that
mean, gives the whole real line: nothing about its concentration can be concluded. A specimen with mean reading \( 0.0775 \)
gives two rays, \( (-\infty,\ -55.3]\cup[7.6,\ \infty) \): the middle of the scale is excluded, but not extreme
concentrations of either sign.
[Figure 12.5.1](#fig-ci-calibration) shows both assays.

In \( 20000 \) simulated experiments at \( x_0=5.3 \), both intervals cover in \( 0.9517 \) for the precise assay. With the slope
three standard errors from zero (\( \beta_1=0.00335 \)) the Fieller set again covers in \( 0.9517 \), since its event involves only the
errors, and is bounded in \( 0.817 \) of the experiments, close to the power \( 0.822 \) of the slope \( t \) test; the delta-method
interval covers in \( 0.9827 \), wastefully. With the slope two standard errors from zero and \( x_0=12 \), beyond the standards,
the Fieller set covers in \( 0.9528 \) and the delta-method interval in only \( 0.9337 \).
:::

::: {when-format="html"}
![**Figure 12.5.1.** Calibration inverts the \( 95\% \) prediction band for the mean of \( m=3 \) readings (shaded); bars mark the calibration sets for the observed readings (horizontal lines). (a) A precise assay, magnified in the inset. (b) A weak assay: whole line (lower bar) and two rays (upper bar).](calibration.svg){#fig-ci-calibration width=100%}
:::

::: {when-format="pdf"}
![Calibration inverts the \( 95\% \) prediction band for the mean of \( m=3 \) readings (shaded); bars mark the calibration sets for the observed readings (horizontal lines). (a) A precise assay, magnified in the inset. (b) A weak assay: whole line (lower bar) and two rays (upper bar).](calibration.pdf){width=100%}
:::

```{.python .run #cell-calibration-assay}
import numpy as np
from scipy import stats

rng = np.random.default_rng(1206)
x = np.repeat(np.arange(9.0), 3)                   # standards, in triplicate
n = len(x)
y = 0.04 + 0.115 * x + 0.015 * rng.normal(size=n)  # absorbance readings
m = 3
y0 = 0.04 + 0.115 * 5.3 + 0.015 * rng.normal(size=m)   # the unknown sample, true x0 = 5.3


def calibrate(x, y, y0bar, m, level=0.95):
    """Fieller set {x0 : |y0bar - b0 - b1 x0| <= t s sqrt(1/m + 1/n + (x0 - xbar)^2 / Sxx)}."""
    n = len(x)
    xbar, Sxx = x.mean(), np.sum((x - x.mean()) ** 2)
    b1 = np.sum((x - xbar) * (y - y.mean())) / Sxx
    s2 = np.sum((y - y.mean() - b1 * (x - xbar)) ** 2) / (n - 2)
    t2 = stats.t.ppf(0.5 + level / 2, n - 2) ** 2
    g = y0bar - y.mean()
    # quadratic in d = x0 - xbar:  A d^2 - 2 B d + K <= 0
    A = b1**2 - t2 * s2 / Sxx
    B = b1 * g
    K = g**2 - t2 * s2 * (1 / m + 1 / n)
    disc = B**2 - A * K
    if disc < 0:
        return "whole line", (-np.inf, np.inf), xbar + g / b1
    r1, r2 = sorted([(B - np.sqrt(disc)) / A, (B + np.sqrt(disc)) / A])
    if A > 0:
        return "interval", (xbar + r1, xbar + r2), xbar + g / b1
    return "two rays", (xbar + r1, xbar + r2), xbar + g / b1   # (-inf, first] and [second, inf)


shape, (lo, hi), x0_hat = calibrate(x, y, y0.mean(), m)
print(f"estimate {x0_hat:.3f}; 95% calibration set: {shape} ({lo:.3f}, {hi:.3f})")
```

::: {.remark}
[Classical or inverse estimation?]

The *inverse* estimator regresses the \( x_i \) on the readings and predicts from \( \bar{Y}_0 \):
\( \tilde{x}_0=\bar{x}+(S_{xy}/S_{yy})(\bar{Y}_0-\bar{Y}) \), with \( S_{xy}=\sum(x_i-\bar{x})(Y_i-\bar{Y}) \) and
\( S_{yy}=\sum(Y_i-\bar{Y})^2 \). Since \( \bar{Y}_0-\bar{Y}=\hat{\beta}_1(\hat{x}_0-\bar{x}) \) and \( \hat{\beta}_1=S_{xy}/S_{xx} \),
\[
\tilde{x}_0-\bar{x}=\frac{S_{xy}^2}{S_{xx}S_{yy}}(\hat{x}_0-\bar{x})=R^2(\hat{x}_0-\bar{x}),
\]
where \( R \) is the sample correlation of the calibration pairs. So \( \tilde{x}_0 \) shrinks \( \hat{x}_0 \) toward \( \bar{x} \) by the factor
\( R^2\in[0,1] \). It has finite moments and can beat \( \hat{x}_0 \) in mean squared error near the centre of the standards
(Krutchkoff 1967); away from the centre it does worse (@exr-ci-inverse-estimator). Hoadley (1970) showed that it is a Bayes
estimator under a particular prior for \( x_0 \).
:::

## Exercises

### A. Check your understanding

::: {#exr-ci-x-intercept}
[A1]

For a straight line, write down Fieller's \( 1-\alpha \) confidence set for the point \( \xi=-\beta_0/\beta_1 \) where the line crosses zero, and
state the condition under which it is a bounded interval. Show that it is the calibration set of @thm-ci-calibration with \( \bar{Y}_0 \) replaced by
\( 0 \) and \( 1/m \) replaced by \( 0 \).
:::

::: {.solution}
With \( \blambda_1=(1,0)\T \) and \( \blambda_2=(0,1)\T \), \( \xi=-a/b \). Apply @prp-ci-fieller to \( \theta=a/b \) and change sign, or directly: \( \xi_0 \) is in the set
iff the \( t \) interval for \( \beta_0+\beta_1\xi_0 \) contains \( 0 \), that is iff
\( (\hat{\beta}_0+\hat{\beta}_1\xi_0)^2\le t^2s^2\{1/n+(\xi_0-\bar{x})^2/S_{xx}\} \). This is @eq-ci-calibration with \( \bar{Y}_0=0 \) and the \( 1/m \) term
absent, because \( 0 \) is a known number rather than an average of new readings. It is bounded iff \( T_1^2>t^2 \).
:::

::: {#exr-ci-calibration-mle}
[A2]

Show that \( \hat{x}_0 \), together with the least squares estimates from the standards, maximizes the joint likelihood of the \( n+m \) readings
over \( (\beta_0,\beta_1,\sigma^2,x_0) \) when \( \hat{\beta}_1\ne0 \). What is the maximum likelihood estimate of \( \sigma^2 \)?
:::

::: {.solution}
For fixed \( \sigma^2 \) the likelihood is maximized by minimizing
\( \sum_i(Y_i-\beta_0-\beta_1x_i)^2+\sum_j(Y_{0j}-\beta_0-\beta_1x_0)^2 \). For any \( \beta_1\ne0 \) the second sum is minimized over \( x_0 \) by choosing
\( \beta_0+\beta_1x_0=\bar{Y}_0 \), leaving \( \sum_j(Y_{0j}-\bar{Y}_0)^2 \), which is free of the parameters. The first sum is then minimized by least
squares on the standards, and \( x_0=(\bar{Y}_0-\hat{\beta}_0)/\hat{\beta}_1=\hat{x}_0 \). The maximum likelihood estimate of \( \sigma^2 \) is
\( \{\text{SSE}+\sum_j(Y_{0j}-\bar{Y}_0)^2\}/(n+m) \).
:::

### B. Practice

::: {#exr-ci-calibration-pooled}
[B1]

When \( m\ge2 \), the new readings carry their own information about \( \sigma^2 \). Let
\( s_p^2=\{\text{SSE}+\sum_j(Y_{0j}-\bar{Y}_0)^2\}/(n+m-3) \). Show that replacing \( s \) by \( s_p \) and \( t_{n-2,\alpha/2} \) by \( t_{n+m-3,\alpha/2} \) in
@eq-ci-calibration gives another exact \( 1-\alpha \) set. Why is it shorter on average, and why is it more often bounded?
:::

::: {.solution}
The within-specimen sum of squares is \( \sigma^2\chi^2(m-1) \), independent of \( \bar{Y}_0 \) (the sample mean and sample variance of a normal sample
are independent by @cor-qf-sample-variance) and of the standards. So \( \{\text{SSE}+\sum_j(Y_{0j}-\bar{Y}_0)^2\}/\sigma^2\sim\chi^2(n+m-3) \), independent of the numerator of the
pivot, and the pivot with \( s_p \) is \( t(n+m-3) \). More degrees of freedom give a smaller multiplier and a less variable estimate of \( \sigma \), so
the set is shorter on average. The bounded case requires \( \hat{\beta}_1^2S_{xx}>t^2s_p^2 \), and a smaller \( t \) makes this easier.
:::

::: {#exr-ci-stationary}
[B2]

In the quadratic regression \( \E Y_i=\beta_0+\beta_1x_i+\beta_2x_i^2 \), use @prp-ci-fieller to give a confidence set for the location
\( x^*=-\beta_1/(2\beta_2) \) of the turning point. When is it a bounded interval, and what does that condition mean for the curvature?
:::

::: {.solution}
Take \( \blambda_1=(0,1,0)\T \), \( \blambda_2=(0,0,2)\T \) and \( \theta=\beta_1/(2\beta_2) \), so \( x^*=-\theta \). The set is
\( \{x_0:(\hat{\beta}_1+2\hat{\beta}_2x_0)^2\le t^2s^2(c_{11}+4c_{12}x_0+4c_{22}x_0^2)\} \), where \( c_{ij} \) are the entries of \( (\X\T\X)^{-1} \)
belonging to \( \beta_1,\beta_2 \). It is the set of \( x_0 \) at which the estimated derivative \( \beta_1+2\beta_2x_0 \) is not significantly different from
zero. It is bounded iff \( \hat{\beta}_2^2>t^2s^2c_{22} \), that is iff the curvature is significant: without significant curvature, the fitted curve
may be nearly straight and its turning point anywhere.
:::

::: {#exr-ci-inverse-estimator}
[B3]

Study the two estimators of the Remark "Classical or inverse estimation?" in an idealized large calibration experiment:
the standards are so numerous that \( \hat{\beta}_0 \), \( \hat{\beta}_1 \) and \( R^2 \) may be replaced by constants \( \beta_0 \), \( \beta_1\ne0 \) and
\( \rho^2\in(0,1) \), while \( \bar{Y}_0 \) is still the mean of \( m \) new readings. Let \( v=\sigma^2/(m\beta_1^2) \).

1. Show that \( \hat{x}_0 \) is unbiased with mean squared error \( v \), and that \( \tilde{x}_0=\bar{x}+\rho^2(\hat{x}_0-\bar{x}) \) has mean squared error
   \( \rho^4v+(1-\rho^2)^2(x_0-\bar{x})^2 \).
2. Show that \( \tilde{x}_0 \) has the smaller mean squared error iff \( (x_0-\bar{x})^2<v(1+\rho^2)/(1-\rho^2) \).
3. Let \( \tau^2=S_{xx}/n \), so that in this limit \( \rho^2=\beta_1^2\tau^2/(\beta_1^2\tau^2+\sigma^2) \). Show that the condition of part 2 reads
   \( (x_0-\bar{x})^2<v+2\tau^2/m \). With \( m=1 \) and a precise instrument, roughly how far from \( \bar{x} \), in units of the spread \( \tau \) of the
   standards, does the inverse estimator stay ahead?
:::

::: {.solution}
1. Here \( \hat{x}_0=x_0+\bar{\varepsilon}_0/\beta_1 \) with \( \bar{\varepsilon}_0\sim\Normal(0,\sigma^2/m) \), so it is unbiased with variance \( v \). Then
   \( \tilde{x}_0-x_0=\rho^2(\hat{x}_0-x_0)-(1-\rho^2)(x_0-\bar{x}) \): variance \( \rho^4v \) plus squared bias \( (1-\rho^2)^2(x_0-\bar{x})^2 \).
2. The inequality \( \rho^4v+(1-\rho^2)^2(x_0-\bar{x})^2<v \) is \( (1-\rho^2)^2(x_0-\bar{x})^2<(1-\rho^2)(1+\rho^2)v \); divide by \( (1-\rho^2)^2>0 \).
3. \( (1+\rho^2)/(1-\rho^2)=(2\beta_1^2\tau^2+\sigma^2)/\sigma^2 \), so \( v(1+\rho^2)/(1-\rho^2)=v+2\tau^2/m \). When \( \sigma \) is small, \( v\approx0 \)
   and the inverse estimator wins only for \( \lvert x_0-\bar{x}\rvert<\sqrt2\,\tau \), a little more than one spread from the centre; beyond that its
   bias dominates.
:::

### C. Going deeper

::: {#exr-ci-fieller-geometry}
[C1]

In the setting of @prp-ci-fieller, let \( E_t=\{\boldsymbol{\upphi}\in\Real^2:((\hat{a},\hat{b})\T-\boldsymbol{\upphi})\T\mathbf{V}^{-1}((\hat{a},\hat{b})\T-\boldsymbol{\upphi})\le t^2s^2\} \), an
ellipse with the \( t \) constant in place of \( 2F_\alpha(2,n-r) \). Show that \( \theta_0\in\mathcal F(\Y) \) iff the line \( \{(\phi_1,\phi_2):\phi_1=\theta_0\phi_2\} \)
meets \( E_t \). Use the picture to explain the three cases: the set is the whole line iff \( E_t \) contains the origin, and it is bounded iff \( E_t \)
does not meet the axis \( \phi_2=0 \).
:::

::: {.solution}
The line is \( \{\boldsymbol{\upphi}:\mathbf{c}\T\boldsymbol{\upphi}=0\} \) with \( \mathbf{c}=(1,-\theta_0)\T \). By @cor-ci-shadows (with the constant \( t^2s^2 \)), the values of
\( \mathbf{c}\T\boldsymbol{\upphi} \) over \( E_t \) form the interval \( \hat{a}-\theta_0\hat{b}\pm ts\sqrt{\mathbf{c}\T\mathbf{V}\mathbf{c}} \), and
\( \mathbf{c}\T\mathbf{V}\mathbf{c}=v_{11}-2\theta_0v_{12}+\theta_0^2v_{22} \). So the line meets \( E_t \) iff this interval contains \( 0 \), which is @eq-ci-fieller. Every
line through the origin meets \( E_t \) iff the origin is in \( E_t \), which is \( Q\le t^2s^2 \). Lines of large slope approach the axis \( \phi_2=0 \)
(the \( \phi_1 \)-axis in these coordinates); the set is unbounded iff lines of arbitrarily large slope meet \( E_t \), iff \( E_t \) meets that axis, which by
@cor-ci-shadows with \( \mathbf{c}=(0,1)\T \) is \( \hat{b}^2\le t^2s^2v_{22} \).
:::

::: {#exr-ci-intersection}
[C2]

Two straight lines \( \E Y=\alpha_k+\beta_kx \), \( k=1,2 \), are fitted to two independent samples in one model with a common \( \sigma^2 \), as in a
two-phase regression. Use @prp-ci-fieller to construct a confidence set for the \( x \)-coordinate \( (\alpha_1-\alpha_2)/(\beta_2-\beta_1) \) of their
intersection. When is it bounded? Interpret the unbounded cases in terms of parallel lines.
:::
