# The classical measurement error model

This section sets up the two versions of the classical model, shows that its normal version cannot identify
the slope from the observed data alone, and separates classical error from Berkson error, which does no harm
to a linear regression.

## Structural and functional models

In a **structural** model the true values are random draws, as in @thm-eiv-attenuation. In a **functional** model
they are fixed unknown constants, one nuisance parameter per observation, which cannot be conditioned on because
they are not observed. Attenuation is the same in both.

::: {#prp-eiv-functional}
[Attenuation with fixed true values]

Let \( \mathbf{t}_1,\mathbf{t}_2,\dots\in\Real^p \) be fixed vectors with first coordinate \( 1 \), such that
\( n^{-1}\sum_{i=1}^n\mathbf{t}_i\mathbf{t}_i\T\to\A \), a positive definite matrix. Let
\( y_i=\mathbf{t}_i\T\boldsymbol{\uptheta}+\varepsilon_i \) and \( \mathbf{s}_i=\mathbf{t}_i+\tilde{\bu}_i \), where the pairs
\( (\tilde{\bu}_i,\varepsilon_i) \) are independent and identically distributed, \( \tilde{\bu}_i \) and \( \varepsilon_i \) are independent with
mean zero and finite variances, the first coordinate of \( \tilde{\bu}_i \) is zero (so that \( \mathbf{s}_i \) keeps the
constant and \( \tilde{\bSigma} \) has zero first row and column), and \( \Cov(\tilde{\bu}_i)=\tilde{\bSigma} \). Then the least squares
coefficient vector from regressing \( y \) on \( \mathbf{s} \) converges in probability to
\( (\A+\tilde{\bSigma})^{-1}\A\boldsymbol{\uptheta} \).
:::

::: {.proof}
Stack the vectors as rows of \( \mathbf{T} \), \( \mathbf{S} \) and \( \tilde{\bU} \). Then
\[
\tfrac1n\mathbf{S}\T\mathbf{S}=\tfrac1n\mathbf{T}\T\mathbf{T}+\tfrac1n\bigl(\mathbf{T}\T\tilde{\bU}+\tilde{\bU}\T\mathbf{T}\bigr)+\tfrac1n\tilde{\bU}\T\tilde{\bU} .
\]
The last term tends to \( \tilde{\bSigma} \) in probability by the weak law of large numbers. An entry
\( n^{-1}\sum_it_{ij}\tilde{u}_{ik} \) of \( n^{-1}\mathbf{T}\T\tilde{\bU} \) has mean zero and variance
\( \tilde{\Sigma}_{kk}\,n^{-2}\sum_it_{ij}^2\to0 \), so it tends to zero by Chebyshev's inequality. Hence \( n^{-1}\mathbf{S}\T\mathbf{S}\to\A+\tilde{\bSigma} \), which is positive definite. In the same way,
\( n^{-1}\mathbf{S}\T\y=n^{-1}\mathbf{S}\T\mathbf{T}\boldsymbol{\uptheta}+n^{-1}\mathbf{T}\T\be+n^{-1}\tilde{\bU}\T\be \). Here
\( n^{-1}\mathbf{S}\T\mathbf{T}\to\A \) by the argument just given, \( n^{-1}\mathbf{T}\T\be\to\bzero \) by Chebyshev's
inequality, and \( n^{-1}\tilde{\bU}\T\be\to\bzero \) by the weak law, since \( \E(\tilde{\bu}_i\varepsilon_i)=\bzero \) by
independence. Matrix inversion is continuous at a nonsingular matrix, which finishes the proof.
:::

For a straight line the slope tends to \( \beta s_x^2/(s_x^2+\sigma_u^2) \), where \( s_x^2 \) is the limiting variance of
the fixed values: what matters is the *spread* of the true values relative to the error. Because the number of
parameters grows with \( n \), maximum likelihood can fail in the functional model (@exm-opt-neyman-scott, @exr-eiv-functional-ml).

## What the data cannot tell you

In the **normal structural model**, \( X\sim\Normal(\mu_x,\sigma_x^2) \), \( U\sim\Normal(0,\sigma_u^2) \) and
\( \varepsilon\sim\Normal(0,\sigma^2) \) are independent, \( W=X+U \) and \( Y=\alpha+\beta X+\varepsilon \). Then
\( (W,Y) \) is bivariate normal (@thm-mvn-linear) with mean \( (\mu_x,\alpha+\beta\mu_x) \) and covariance matrix
\[
\begin{pmatrix}\sigma_{ww}&\sigma_{wy}\\\sigma_{wy}&\sigma_{yy}\end{pmatrix}
=\begin{pmatrix}\sigma_x^2+\sigma_u^2&\beta\sigma_x^2\\\beta\sigma_x^2&\beta^2\sigma_x^2+\sigma^2\end{pmatrix}.
\]{#eq-eiv-normal-structural}

A bivariate normal distribution is fixed by five numbers, but the model has six parameters
\( (\mu_x,\alpha,\beta,\sigma_x^2,\sigma_u^2,\sigma^2) \), so data on \( (W,Y) \) cannot determine all six.

::: {#prp-eiv-nonidentified}
[The slope is identified only up to an interval]

Let \( (W,Y) \) be bivariate normal with covariance matrix entries \( \sigma_{ww},\sigma_{wy},\sigma_{yy} \), where
\( \sigma_{ww}\sigma_{yy}>\sigma_{wy}^2 \). If \( \sigma_{wy}\ne0 \), a number \( b \) is the slope of some normal
structural model with \( \sigma_x^2>0 \), \( \sigma_u^2\ge0 \) and \( \sigma^2\ge0 \) that produces this distribution iff \( b \)
lies in the closed interval with end points
\[
b_1=\frac{\sigma_{wy}}{\sigma_{ww}},\qquad b_2=\frac{\sigma_{yy}}{\sigma_{wy}} .
\]
If \( \sigma_{wy}=0 \), the only such slope is \( b=0 \).
:::

::: {.proof}
A normal distribution is determined by its mean and covariance matrix, so a model reproduces the distribution iff it
reproduces both. The means can always be matched by \( \mu_x=\E W \) and \( \alpha=\E Y-b\mu_x \). By
@eq-eiv-normal-structural the covariances are matched iff
\[
\sigma_x^2+\sigma_u^2=\sigma_{ww},\qquad b\,\sigma_x^2=\sigma_{wy},\qquad b^2\sigma_x^2+\sigma^2=\sigma_{yy}.
\]
If \( \sigma_{wy}=0 \), the middle equation with \( \sigma_x^2>0 \) forces \( b=0 \), and then \( \sigma_x^2\in(0,\sigma_{ww}] \) and
\( \sigma^2=\sigma_{yy} \) complete a model. If \( \sigma_{wy}\ne0 \), then \( b\ne0 \) and the equations force
\( \sigma_x^2=\sigma_{wy}/b \), \( \sigma_u^2=\sigma_{ww}-\sigma_{wy}/b \) and \( \sigma^2=\sigma_{yy}-b\sigma_{wy} \). These are admissible
iff \( b \) has the sign of \( \sigma_{wy} \) (so \( \sigma_x^2>0 \)), \( \lvert b\rvert\ge\lvert\sigma_{wy}\rvert/\sigma_{ww}=\lvert b_1\rvert \) (so
\( \sigma_u^2\ge0 \)) and \( \lvert b\rvert\le\sigma_{yy}/\lvert\sigma_{wy}\rvert=\lvert b_2\rvert \) (so \( \sigma^2\ge0 \)). Since \( b_1 \) and \( b_2 \) share the sign of
\( \sigma_{wy} \), these conditions say that \( b \) lies between \( b_1 \) and \( b_2 \). The interval is not empty, because
\( \lvert b_1\rvert\le\lvert b_2\rvert \) is equivalent to \( \sigma_{wy}^2\le\sigma_{ww}\sigma_{yy} \).
:::

The end point \( b_1 \) is the least squares slope of \( Y \) on \( W \), correct if there is no measurement error.
The end point \( b_2 \) is the inverted slope of \( W \) on \( Y \), correct if there is no equation error. Their ratio is the
squared correlation, so the interval is short only when the correlation is high.

::: {#exm-eiv-identified-set}
[The interval for the blood-pressure study]

For the first reading in the running example, the sample covariances are \( s_{ww}=216.21 \),
\( s_{wy}=10.786 \) and \( s_{yy}=1.889 \). The compatible slopes form the interval
\( [0.0499,\,0.1752] \), and the squared correlation is \( 0.285 \). It contains the true \( 0.08 \), whose implied variances \( \sigma_x^2=134.8 \),
\( \sigma_u^2=81.4 \) and \( \sigma^2=1.026 \) are close to the generating \( 144 \), \( 81 \) and \( 1 \), but
the data cannot prefer them to those of any other slope in the interval
([Figure 24.2.1](#fig-eiv-identified-set)).
:::

::: {when-format="html"}
![**Figure 24.2.1.** (a) Direct and reverse regression lines and slopes in between (grey), each fitting a normal
structural model exactly; dashed, the true slope. (b) The reliability \( b_1/b \) and squared correlation \( b/b_2 \)
that slope \( b \) would require.](identified_set.svg){#fig-eiv-identified-set width=100%}
:::

::: {when-format="pdf"}
![(a) Direct and reverse regression lines and slopes in between (grey), each fitting a normal structural model
exactly; dashed, the true slope. (b) The reliability \( b_1/b \) and squared correlation \( b/b_2 \) that slope \( b \)
would require.](identified_set.pdf){width=100%}
:::

```{.python .run #cell-identification-interval}
import numpy as np

rng = np.random.default_rng(2401)
n = 300
age = rng.normal(55, 10, n)
x = 130 + 0.6 * (age - 55) + rng.normal(0, np.sqrt(108), n)   # long-run blood pressure, sd 12
w = x[:, None] + rng.normal(0, 9, (n, 2))                      # two clinic readings
y = -4 + 0.08 * x + 0.0 * age + rng.normal(0, 1.0, n)          # age has no effect given x

S = np.cov(w[:, 0], y)                       # sample covariance matrix of (reading, outcome)
s_ww, s_wy, s_yy = S[0, 0], S[0, 1], S[1, 1]
b_direct = s_wy / s_ww                       # assumes no measurement error (sigma_u^2 = 0)
b_reverse = s_yy / s_wy                      # assumes no equation error (sigma^2 = 0)
print(f"slopes compatible with the data: [{b_direct:.4f}, {b_reverse:.4f}]")

for b in [b_direct, 0.08, b_reverse]:        # the variances each slope would imply
    sx2 = s_wy / b
    print(f"b = {b:.4f}: sigma_x^2 = {sx2:7.2f}, sigma_u^2 = {s_ww - sx2:7.2f}, "
          f"sigma^2 = {s_yy - b * s_wy:.4f}")
```

The sign of \( \beta \), whether \( \beta=0 \), and the lower bound \( \lvert b_1\rvert \) on \( \lvert\beta\rvert \) are identified.
Any information that pins down one variance identifies everything: \( \sigma_u^2 \), \( \lambda \), the ratio
\( \sigma^2/\sigma_u^2 \), or a bound \( \lambda\ge\lambda_0 \) that narrows the interval (@exr-eiv-bound). The rest of the
chapter is about where such information comes from.

::: {.remark}
[Nonnormal true values]

Reiersøl (1950) showed that when \( \beta\ne0 \) and the errors are normal and independent of \( X \), the slope is identified
exactly when \( X \) is *not* normal; higher moments then carry the missing information. We do not prove this.
Estimators based on higher moments are highly variable unless \( X \) is far from normal, and rarely replace
external information.
:::

## Berkson error

Sometimes the truth scatters around the recorded value instead: a furnace set to \( w \) degrees, every resident
of a district assigned the district monitor's pollution level, a solution prepared at a nominal concentration.
Then
\[
X=W+U,\qquad U \text{ independent of } W,
\]
which is **Berkson error** (Berkson 1950). The difference from classical error is which variable the error is
independent of. It cannot be settled from data on \( (W,Y) \), only from how the recorded value was produced. The
scalar case with fixed settings appeared in [Section 19.5](../ch19-theory-of-departures/05-random-regressors.html);
the proposition below allows several regressors and random recorded values.

::: {#prp-eiv-berkson}
[Berkson error in linear regression]

Let \( Y=\alpha+\bbeta\T\x+\varepsilon \) and \( \x=\bw+\bu \), where \( \bw \) is fixed or random and, given \( \bw \), the errors
\( \bu \) and \( \varepsilon \) are independent, with means zero, covariance matrix \( \bSigma_{uu} \) and variance \( \sigma^2 \), not depending on \( \bw \).
Then, given \( \bw \),
\[
\E(Y\mid\bw)=\alpha+\bbeta\T\bw,\qquad \Var(Y\mid\bw)=\sigma^2+\bbeta\T\bSigma_{uu}\bbeta ,
\]
and if \( \bu \) and \( \varepsilon \) are normal, so is \( Y \) given \( \bw \). Independent observations of this kind therefore
satisfy the linear model in the recorded \( \bw \), with the same coefficients and error variance
\( \sigma^2+\bbeta\T\bSigma_{uu}\bbeta \).
:::

::: {.proof}
Given \( \bw \), \( Y=\alpha+\bbeta\T\bw+(\bbeta\T\bu+\varepsilon) \), and the composite error is a linear combination of the
independent \( \bu \) and \( \varepsilon \) (@thm-rv-linear).
:::

So least squares on the settings is unbiased, and under normality the \( t \) and \( F \) procedures are exact,
with \( s^2 \) estimating the enlarged error variance. Nothing needs correcting.

::: {#exm-eiv-berkson}
[Concentrations set on a dial]

Five solutions are prepared at each nominal concentration \( 10,20,\dots,60 \) mg/L (\( n=30 \)); the achieved
concentration adds \( \Normal(0,4^2) \) error, and the response is \( 2+0.5x+\varepsilon \), \( \varepsilon\sim\Normal(0,3^2) \).
Regressing on the nominal concentrations in
\( 20000 \) simulated experiments gives a mean slope of \( 0.4994 \) and a mean \( s^2 \) of
\( 13.01 \), matching \( \sigma^2+\beta^2\sigma_u^2=9+0.25\cdot16=13 \). The nominal 95% \( t \) intervals for
the slope covered \( 0.5 \) in a proportion \( 0.952 \) of the experiments.
:::

```{.python .run #cell-identification-berkson}
from scipy import stats

rng_b = np.random.default_rng(2410)
settings = np.repeat([10.0, 20.0, 30.0, 40.0, 50.0, 60.0], 5)      # dial settings, n = 30
Xd = np.column_stack([np.ones(30), settings])
reps = 20_000
achieved = settings + rng_b.normal(0, 4, (reps, 30))               # x = w + u, u independent of w
Y = 2 + 0.5 * achieved + rng_b.normal(0, 3, (reps, 30))
B = np.linalg.solve(Xd.T @ Xd, Xd.T @ Y.T)                          # least squares on the settings
s2 = np.sum((Y.T - Xd @ B) ** 2, axis=0) / 28
se = np.sqrt(s2 * np.linalg.inv(Xd.T @ Xd)[1, 1])
covered = np.abs(B[1] - 0.5) <= stats.t.ppf(0.975, 28) * se
print(f"mean slope {B[1].mean():.4f}, mean s^2 {s2.mean():.2f} (sigma^2 + beta^2 sigma_u^2 = 13), "
      f"coverage {covered.mean():.3f}")
```

::: {.warning}
Berkson error is harmless only because the model is linear in \( \x \). For nonlinear \( g \),
\( \E\{g(\x)\mid\bw\}\ne g(\bw) \) in general (@exr-eiv-berkson-quadratic). Real measurements often mix Berkson and
classical components, which must then be handled separately.
:::

## Exercises

### A. Check your understanding

::: {#exr-eiv-classify}
[A1]

Classify each regressor error as classical, Berkson or neither, and say whether it biases the slope: (a) weight
from a noisy scale; (b) a prescribed dose, when patients take a random fraction of it with mean one; (c) income
rounded down to the nearest 10,000.
:::

::: {#exr-eiv-bound}
[A2]

In @prp-eiv-nonidentified, suppose it is known only that the reliability is at least \( \lambda_0\in(0,1] \). Show
that the compatible slopes are those between \( b_1 \) and \( \min(\lvert b_2\rvert,\lvert b_1\rvert/\lambda_0) \) in magnitude, with the sign of
\( b_1 \). Evaluate the interval for @exm-eiv-identified-set with \( \lambda_0=0.5 \).
:::

::: {.solution}
The implied reliability of slope \( b \) is \( \sigma_x^2/\sigma_{ww}=b_1/b \), so \( \lambda\ge\lambda_0 \) iff
\( \lvert b\rvert\le\lvert b_1\rvert/\lambda_0 \). Intersect with the interval of the proposition. For the example,
\( b_1/\lambda_0=0.0499/0.5=0.0998<0.1752 \), so the interval shrinks to \( [0.0499,\,0.0998] \), which
still contains \( 0.08 \).
:::

### B. Practice

::: {#exr-eiv-berkson-quadratic}
[B1]

Let \( Y=\beta_0+\beta_1X+\beta_2X^2+\beta_3X^3+\varepsilon \) and \( X=W+U \) with Berkson error \( U\sim\Normal(0,\sigma_u^2) \)
independent of \( W \) and \( \varepsilon \). Compute \( \E(Y\mid W) \). Which coefficients of a cubic regression on \( W \) are
biased, and by how much? What happens for a quadratic?
:::

::: {.solution}
With \( \E U=\E U^3=0 \) and \( \E U^2=\sigma_u^2 \),
\( \E(X^2\mid W)=W^2+\sigma_u^2 \) and \( \E(X^3\mid W)=W^3+3W\sigma_u^2 \). So
\[
\E(Y\mid W)=(\beta_0+\beta_2\sigma_u^2)+(\beta_1+3\beta_3\sigma_u^2)W+\beta_2W^2+\beta_3W^3 .
\]
The cubic regression on \( W \) estimates the intercept with bias \( \beta_2\sigma_u^2 \) and the linear coefficient with bias
\( 3\beta_3\sigma_u^2 \). The top two coefficients are unbiased. For a quadratic (\( \beta_3=0 \)) only the intercept is
biased. The conditional variance also changes with \( W \), so the errors are heteroscedastic.
:::

::: {#exr-eiv-reverse-limit}
[B2]

In the normal structural model, show that the reciprocal of the least squares slope of \( W \) on \( Y \) converges to
\( b_2=\beta+\sigma^2/(\beta\sigma_x^2) \), so that it overstates \( \lvert\beta\rvert \) exactly when there is equation
error. Evaluate \( b_2 \) for the population of the running example and compare it with @exm-eiv-identified-set.
:::

### C. Going deeper

::: {#exr-eiv-functional-ml}
[C1]

In the functional model \( w_i=x_i+u_i \), \( y_i=\alpha+\beta x_i+\varepsilon_i \), with \( u_i,\varepsilon_i \) independent
\( \Normal(0,\sigma^2) \) (equal, unknown variances) and \( x_1,\dots,x_n \) unknown constants:

::: {.enumerate options="label=(\alph*)"}
1. show that, for fixed \( (\alpha,\beta) \), maximizing the likelihood over each \( x_i \) leaves
   \( d_i^2=(y_i-\alpha-\beta w_i)^2/(1+\beta^2) \), the squared perpendicular distance from \( (w_i,y_i) \) to the line;

2. deduce that the maximum likelihood line is the orthogonal regression line of [Section 24.3](03-moment-correction.html), and that
   \( \hat{\sigma}^2_{\text{ML}}=(2n)^{-1}\sum_i\hat{d}_i^2 \);

3. show that \( \hat{\sigma}^2_{\text{ML}}\to\sigma^2/2 \) in probability, assuming the line is estimated consistently.
   Compare @exm-opt-neyman-scott.
:::

:::

::: {.solution}
*(a)* The log likelihood is \( -n\log(2\pi\sigma^2)-\sum_i\{(w_i-x_i)^2+(y_i-\alpha-\beta x_i)^2\}/(2\sigma^2) \). For each \( i \) the
bracket is the squared distance from \( (w_i,y_i) \) to the point \( (x_i,\alpha+\beta x_i) \) of the line, and its minimum
over \( x_i \) is the squared perpendicular distance \( d_i^2 \). (b) The profile log likelihood is
\( -n\log(2\pi\sigma^2)-\sum_id_i^2/(2\sigma^2) \). It is maximized over the line by minimizing \( \sum_id_i^2 \), which is orthogonal
regression, and over \( \sigma^2 \) by \( \sum_i\hat{d}_i^2/(2n) \). (c) At the true line,
\( y_i-\alpha-\beta w_i=\varepsilon_i-\beta u_i \) has variance \( \sigma^2(1+\beta^2) \), so \( n^{-1}\sum_id_i^2\to\sigma^2 \) by the law
of large numbers. A consistent estimate of the line gives the same
limit, so \( \hat{\sigma}^2_{\text{ML}}\to\sigma^2/2 \): each nuisance \( x_i \) uses up half of its observation's two degrees of freedom.
:::
