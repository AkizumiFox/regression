# Estimating the error variance

The covariance formula \( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1} \) is of little use until
\( \sigma^2 \) is estimated. Least squares does not estimate it: the criterion
\( \norm{\y-\X\bb}^2 \) involves only \( \bbeta \). But the size of the residuals clearly
carries information about the size of the errors. This section finds the right way to
turn one into an estimate of the other, and shows why the obvious way is wrong.

## The residual sum of squares

Since \( \sigma^2=\Var(\varepsilon_i)=\E(\varepsilon_i^2) \), a natural first idea is to replace
each unobserved error by its residual and average: \( \text{SSE}/n \). This underestimates
\( \sigma^2 \). The fitted values were chosen to make the residuals as small as possible, so
the residuals are systematically smaller than the errors. The next theorem measures
exactly how much smaller.

::: {#thm-lm-sigma2}
[Unbiased estimation of the error variance]

In the linear model @def-lm-linear-model with \( \rank(\X)=p<n \),
\[
\E(\text{SSE})=\sigma^2(n-p).
\]
Hence
\[
s^2=\frac{\text{SSE}}{n-p}=\frac{\norm{\Y-\X\hbeta}^2}{n-p}
\]{#eq-lm-s2}

is an unbiased estimator of \( \sigma^2 \).
:::

::: {.proof}
By @prp-lm-fit-algebra(e) the residual vector is \( \he=(\I-\bH)\be \), whatever \( \bbeta \) is.
Since \( \I-\bH \) is symmetric and idempotent,
\( \text{SSE}=\he\T\he=\be\T(\I-\bH)\be \), a quadratic form in the error vector, which has mean
\( \bzero \) and covariance \( \sigma^2\I_n \). By @cor-rv-quadform-special(b),
\[
\E(\text{SSE})=\sigma^2\tr(\I-\bH)+\bzero\T(\I-\bH)\bzero=\sigma^2(n-p),
\]
using \( \tr(\I-\bH)=n-p \) from @prp-lm-fit-algebra(a).
:::

The computation is the case \( \boldsymbol{\theta}=\X\bbeta \) of @exm-rv-rss-bias in
[Chapter 2](../ch02-random-vectors/index.html), specialized to full rank.

The same result follows from @prp-lm-fit-moments without quadratic forms. Each residual has
mean zero and variance \( \sigma^2(1-h_{ii}) \), so
\( \E(\text{SSE})=\sum_i\E(\hat{\varepsilon}_i^2)=\sigma^2\sum_i(1-h_{ii})=\sigma^2(n-p) \). This version
shows *where* the shortfall comes from: each residual is less variable than its error by the
factor \( 1-h_{ii} \), and the leverages add up to \( p \).

The divisor \( n-p \) is the number of **residual degrees of freedom**. The name reflects the
constraints \( \X\T\he=\bzero \) of @prp-lm-fit-algebra(b): the \( n \) residuals satisfy \( p \)
independent linear equations, so only \( n-p \) of them can vary freely.
[Chapter 6](../ch06-projections/index.html) makes this precise: the residual vector lies in a
subspace of dimension \( n-p \), and the degrees of freedom of a sum of squares is always the
dimension of the space it lives in. For the straight line \( p=2 \), so
\( s^2=\text{SSE}/(n-2) \). For the sample variance of a single sample, the model matrix is
\( \bone \), \( p=1 \), and @eq-lm-s2 is the familiar divisor \( n-1 \).

::: {#cor-lm-cov-estimate}
[Estimated covariance and standard errors]

Under the conditions of @thm-lm-sigma2, \( s^2(\X\T\X)^{-1} \) is an unbiased estimator of
\( \Cov(\hbeta) \), and \( s^2\mathbf{a}\T(\X\T\X)^{-1}\mathbf{a} \) is an unbiased estimator of
\( \Var(\mathbf{a}\T\hbeta) \). The **standard error** of \( \hat{\beta}_j \) is
\[
\operatorname{se}(\hat{\beta}_j)=s\sqrt{\bigl[(\X\T\X)^{-1}\bigr]_{jj}} .
\]
:::

::: {.proof}
\( (\X\T\X)^{-1} \) and \( \mathbf{a} \) are constants, so the expectations are
\( \E(s^2)(\X\T\X)^{-1}=\sigma^2(\X\T\X)^{-1} \) and similarly for the scalar.
:::

The standard error is not unbiased for the standard deviation of \( \hat{\beta}_j \). By Jensen's
inequality \( \E(s)\le\sqrt{\E(s^2)}=\sigma \), with strict inequality unless \( s \) is
constant. The bias is small when \( n-p \) is large, and it matters little, since the
standard error is used mainly through the ratio \( (\hat{\beta}_j-\beta_j)/\operatorname{se}(\hat{\beta}_j) \),
whose exact distribution under normality is the \( t \) distribution with \( n-p \) degrees of
freedom ([Chapter 7](../ch07-optimality/index.html), @thm-opt-sampling).

::: {#exm-lm-money-se}
[Standard errors for money demand]

For the money demand model of @exm-lm-money, \( n=55 \) and \( p=4 \), so
\[
s^2=\frac{0.09253}{51}=0.001814,\qquad s=0.04259 .
\]
Since the response is a logarithm, \( s \) says that money holdings typically deviate from
the fitted relation by about \( 4 \) percent. The standard errors of the four coefficients are
\( 0.5811 \), \( 0.0940 \), \( 0.3282 \) and
\( 0.6911 \). The income coefficient \( 1.2958 \) is almost
fourteen standard errors from zero and the bond rate coefficient \( -2.6163 \)
about eight, but the deposit rate coefficient \( 0.6186 \) is smaller than its
standard error. Whether the deposit rate belongs in the model at all is a testing question
for [Chapter 11](../ch11-general-linear-hypothesis/index.html).

These standard errors rest on (L2) and (L3). For a quarterly economic series (L3) is
doubtful. The correlation between consecutive residuals is \( 0.62 \), which
suggests that the errors are positively correlated in time, and then the formula
\( \sigma^2(\X\T\X)^{-1} \) usually understates the true variability. Chapter 21 shows how to
detect and allow for this.
:::

```{.python .run #cell-money-demand-sigma2}
import numpy as np
import statsmodels.api as sm
dk = sm.datasets.danish_data.load_pandas().data
y = dk["lrm"].to_numpy()
X = np.column_stack([np.ones(len(y)), dk["lry"], dk["ibo"], dk["ide"]])
n, p = X.shape
XtX = X.T @ X
Xty = X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)       # solve X^T X b = X^T y
fitted = X @ beta_hat
resid = y - fitted
SSE = resid @ resid

s2 = SSE / (n - p)                          # unbiased estimate of sigma^2
cov_hat = s2 * np.linalg.inv(XtX)           # estimated Cov(beta_hat)
se = np.sqrt(np.diag(cov_hat))
print(f"s^2 = {s2:.6f}, s = {np.sqrt(s2):.5f}")
print("standard errors:", np.round(se, 4))
print("SSE / n would give", SSE / n)
```

## Dividing by n

The estimator \( \hat{\sigma}^2=\text{SSE}/n \) is not unbiased, but it is not arbitrary either:
under normality it is the maximum likelihood estimator (@thm-opt-mle). Its expectation is
\[
\E(\hat{\sigma}^2)=\frac{n-p}{n}\,\sigma^2=\sigma^2-\frac{p}{n}\sigma^2 .
\]
The bias is negligible when \( p \) is a small fraction of \( n \), and severe when it is not. For
the money demand data the two estimates of \( \sigma \) are \( 0.04259 \) and
\( 0.04102 \), a difference of about 4 percent. But a model with many
columns relative to the number of cases, such as a two-way layout with a few replicates per
cell or a regression with many indicator variables, can make \( \text{SSE}/n \) badly
misleading.

::: {#exm-lm-divisor}
[How bad is dividing by n?]

Take \( n=30 \) cases and model matrices with an intercept and \( p-1 \) further columns of
standard normal numbers, for \( p=1,\dots,25 \). For each \( p \),
20,000 error vectors are simulated with \( \sigma^2=1 \). [Figure 5.6.1](#fig-lm-divisor)
plots the average of \( \text{SSE}/n \) and of \( \text{SSE}/(n-p) \). The first follows the line
\( (n-p)/n \): at \( p=5 \) it averages \( 0.834 \), at \( p=15 \) it averages
\( 0.500 \), half the true value, and at \( p=25 \) only
\( 0.167 \). The unbiased estimator averages \( 1.000 \),
\( 1.001 \) and \( 0.999 \) at the same values of \( p \).
:::

::: {when-format="html"}
![**Figure 5.6.1.** Average of two estimators of \( \sigma^2=1 \) over 20,000 simulated
data sets, with \( n=30 \) and \( p \) columns. Dividing the residual sum of squares by \( n \)
underestimates \( \sigma^2 \) by the factor \( (n-p)/n \); dividing by \( n-p \) is unbiased for
every \( p \).](sigma2_divisor.svg){#fig-lm-divisor width=72%}
:::

::: {when-format="pdf"}
![Average of two estimators of \( \sigma^2=1 \) over 20,000 simulated
data sets, with \( n=30 \) and \( p \) columns. Dividing the residual sum of squares by \( n \)
underestimates \( \sigma^2 \) by the factor \( (n-p)/n \); dividing by \( n-p \) is unbiased for
every \( p \).](sigma2_divisor.pdf){width=72%}
:::

```{.python .run #cell-sigma2-divisor-divisor}
import numpy as np
rng = np.random.default_rng(506)
n, sigma2, reps = 30, 1.0, 20_000
Xfull = np.column_stack([np.ones(n), rng.normal(size=(n, 24))])

ps = np.arange(1, 26)
mean_n, mean_np = [], []
for p in ps:
    X = Xfull[:, :p]
    R = np.eye(n) - X @ np.linalg.solve(X.T @ X, X.T)     # I - H
    E = rng.standard_normal((reps, n))                      # errors, sigma^2 = 1
    SSE = np.einsum("ri,ij,rj->r", E, R, E)                 # e^T (I - H) e
    mean_n.append(SSE.mean() / n)
    mean_np.append(SSE.mean() / (n - p))
print("p     E(SSE/n)  E(SSE/(n-p))")
for p in (1, 5, 15, 25):
    print(f"{p:2d}    {mean_n[p - 1]:.3f}     {mean_np[p - 1]:.3f}")
```

## The variance of s²

Unbiasedness says where \( s^2 \) is centred, not how much it varies. Its variance
depends on the fourth moment of the errors, so the second-moment assumptions are not enough
to find it. With independent errors, [Chapter 2](../ch02-random-vectors/index.html)
supplies the formula.

::: {#prp-lm-var-s2}
[Variance of s²]

In the linear model with \( \rank(\X)=p<n \), suppose the errors \( \varepsilon_1,\dots,\varepsilon_n \) are
independent with mean \( 0 \), variance \( \sigma^2 \), common third moment
\( \mu_3=\E(\varepsilon_i^3) \) and common fourth moment \( \mu_4=\E(\varepsilon_i^4)<\infty \). Then
\[
\Var(\text{SSE})=2\sigma^4(n-p)+(\mu_4-3\sigma^4)\sum_{i=1}^n(1-h_{ii})^2 .
\]
For normal errors, \( \mu_4=3\sigma^4 \) and \( \Var(s^2)=2\sigma^4/(n-p) \).
:::

::: {.proof}
Apply @thm-rv-quadform-variance to the vector \( \be \), whose means \( \theta_i \) are all zero, with
\( \A=\I-\bH \). The terms involving \( \boldsymbol{\theta} \) vanish (in particular \( \mu_3 \) enters only
multiplied by \( \boldsymbol{\theta}=\bzero \), so its value is irrelevant). The diagonal of \( \A \) has entries
\( 1-h_{ii} \), and \( \tr(\A^2)=\tr(\A)=n-p \) because \( \A \) is idempotent. For the normal case,
the fourth moment of a \( \Normal(0,\sigma^2) \) variable is \( 3\sigma^4 \), and
\( \Var(s^2)=\Var(\text{SSE})/(n-p)^2 \).
:::

The kurtosis term can be large. Heavy-tailed errors, with \( \mu_4 \) much larger than
\( 3\sigma^4 \), make \( s^2 \) much less reliable than the normal-theory formula suggests, even
though it remains unbiased. Since \( 0\le1-h_{ii}\le1 \), the sum in the kurtosis term is
at most \( \sum_i(1-h_{ii})=n-p \), which gives the following.

::: {#cor-lm-s2-consistent}
[Consistency of s²]

Under the conditions of @prp-lm-var-s2,
\[
\Var(s^2)\le\frac{2\sigma^4+\lvert\mu_4-3\sigma^4\rvert}{n-p}.
\]
So if \( n-p\to\infty \) along a sequence of models, \( s^2\to\sigma^2 \) in mean square and in
probability.
:::

::: {.proof}
Divide the formula of @prp-lm-var-s2 by \( (n-p)^2 \) and bound the sum by \( n-p \). Since
\( \E(s^2)=\sigma^2 \), \( \E(s^2-\sigma^2)^2=\Var(s^2)\to0 \), and Chebyshev's inequality gives
convergence in probability.
:::

::: {#exm-lm-kurtosis}
[Skewed errors inflate the variance of s²]

With \( n=30 \) and the first \( p=6 \) columns of the model matrix of
@exm-lm-divisor, the leverages give \( \sum_i(1-h_{ii})^2=19.431 \). For normal
errors with \( \sigma^2=1 \), \( \Var(\text{SSE})=2(n-p)=48.00 \). For centred
exponential errors, which have \( \sigma^2=1 \) and \( \mu_4=9 \), @prp-lm-var-s2 gives
\( 48+6\times19.431=164.59 \), and a simulation with 200,000 data sets
gives \( 165.19 \). The standard deviation of \( s^2 \) is
\( 0.289 \) under normal errors and \( 0.535 \) under exponential
errors, nearly twice as large, although \( s^2 \) is unbiased in both cases.
:::

```{.python .run #cell-sigma2-divisor-variance}
p = 6
X = Xfull[:, :p]
H = X @ np.linalg.solve(X.T @ X, X.T)
R = np.eye(n) - H
mu4 = 9.0                                    # fourth moment of a centred Exp(1) variable
var_normal = 2 * (n - p)
var_exp = (mu4 - 3) * np.sum(np.diag(R) ** 2) + 2 * (n - p)
reps2 = 200_000
E = rng.exponential(size=(reps2, n)) - 1.0
SSE = np.einsum("ri,ij,rj->r", E, R, E)
print(f"Var(SSE): normal errors {var_normal:.2f}, exponential errors {var_exp:.2f},"
      f" simulated {SSE.var():.2f}")
```

Under normality much more is true: \( \text{SSE}/\sigma^2 \) has the chi-squared
distribution with \( n-p \) degrees of freedom and is independent of \( \hbeta \)
(@thm-opt-sampling, using @thm-qf-chisq and @thm-qf-indep-linear). That is what turns
\( \hbeta \), \( s^2 \) and the standard errors into exact confidence intervals and tests.
Among unbiased estimators of \( \sigma^2 \), \( s^2 \) then has the smallest variance
(@thm-opt-umvue).

## When the mean is wrong

@thm-lm-sigma2 assumed that \( \E(\Y)=\X\bbeta \). If instead \( \E(\Y)=\boldsymbol{\theta} \) is not of this
form, @exm-rv-rss-bias in [Chapter 2](../ch02-random-vectors/index.html) shows that
\[
\E(\text{SSE})=\sigma^2(n-p)+\norm{(\I-\bH)\boldsymbol{\theta}}^2 .
\]
The extra term is the squared distance from the true mean vector to the nearest mean vector the
model can produce. So an incomplete model *overestimates* \( \sigma^2 \): variation that the
model fails to capture is counted as error. This is one reason to prefer a model that
includes all relevant regressors, and it is the basis of the lack-of-fit tests of
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) (@thm-cor-lack-of-fit), which compare \( s^2 \) with an estimate of \( \sigma^2 \) that does not depend on the
form of the mean.

## Exercises

### A. Check your understanding

::: {#exr-lm-s2-origin}
[A1]

For regression through the origin, \( \E(Y_i)=\beta x_i \), what is the unbiased estimator of
\( \sigma^2 \)? Verify \( \E(\text{SSE})=\sigma^2(n-1) \) directly, by writing
\( \text{SSE}=\sum_i\varepsilon_i^2-(\sum_ix_i\varepsilon_i)^2/\sum_ix_i^2 \).
:::

::: {#exr-lm-se-elnino}
[A2]

For @exm-lm-elnino, compute \( s \) and the standard error of the slope. Is the slope more
than two standard errors away from \( 1 \), the value that would mean December simply
tracks August one for one?
:::

### B. Practice

::: {#exr-lm-sse-formulas}
[B1]

Show that \( \text{SSE}=\y\T\y-\hbeta\T\X\T\X\hbeta \), and that for a model with an intercept
\( \text{SSE}=S_{yy}-\hbeta_1\T\X_c\T\X_c\hbeta_1 \) in the notation of @prp-lm-centred. Explain why
the second form is numerically preferable when the responses are large and nearly equal, as
in @exm-lm-co2.
:::

::: {.solution}
By @prp-lm-fit-algebra(c), \( \text{SSE}=\y\T\y-\hbeta\T\X\T\y \), and \( \X\T\y=\X\T\X\hbeta \) by
the normal equations. For the centred form, the model with columns \( [\bone,\X_c] \) has the
same fitted values (the proof of @prp-lm-centred), with coefficients \( (\bar{y},\hbeta_1) \). So
\( \text{SSE}=\y\T\y-n\bar{y}^2-\hbeta_1\T\X_c\T\X_c\hbeta_1 \), because the cross-product matrix is block diagonal, and
\( \y\T\y-n\bar{y}^2=S_{yy} \). When the responses are about \( 340 \) with variation of a few
units, \( \y\T\y \) and \( \hbeta\T\X\T\y \) are both near \( n\cdot340^2 \) and agree in their leading
digits; their difference loses those digits to cancellation. The centred form subtracts
numbers of the size of the variation itself.
:::

::: {#exr-lm-mle-mse}
[B2]

Under normality \( \text{SSE}/\sigma^2 \) has mean \( n-p \) and variance \( 2(n-p) \). Among
estimators \( c\cdot\text{SSE} \), find the \( c \) that minimizes the mean squared error
\( \E(c\,\text{SSE}-\sigma^2)^2 \). Compare with \( 1/(n-p) \) and \( 1/n \).
:::

::: {.solution}
Write \( \text{SSE}=\sigma^2Q \) with \( \E Q=m=n-p \) and \( \E Q^2=2m+m^2 \). Then
\( \E(cQ-1)^2\sigma^4=\sigma^4\bigl(c^2(m^2+2m)-2cm+1\bigr) \), minimized at
\( c=m/(m^2+2m)=1/(n-p+2) \). The minimum-MSE estimator divides by \( n-p+2 \): it is biased
downward but less variable. The ordering of divisors is \( n-p<n-p+2 \) and, when \( p\le2 \),
\( n\le n-p+2 \), so for small \( p \) the maximum likelihood divisor \( n \) lies between the unbiased and
the minimum-MSE choice. Unbiasedness is a convention, not an optimality property; it is
kept because it makes \( s^2 \) combine cleanly in the \( t \) and \( F \) statistics.
:::

### C. Going deeper

::: {#exr-lm-quadratic-unbiased}
[C1]

Show that a quadratic form \( \Y\T\A\Y \) with \( \A \) symmetric is unbiased for \( \sigma^2 \) for
every \( \bbeta \) and \( \sigma^2 \) iff \( \X\T\A\X=\bzero \) and \( \tr\A=1 \). Verify that
\( \A=(\I-\bH)/(n-p) \) qualifies, and find another \( \A \) that does. (Under normality,
@thm-opt-umvue shows that \( s^2 \) has the smallest variance among all unbiased estimators.)
:::
