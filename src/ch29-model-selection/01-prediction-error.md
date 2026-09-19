# Prediction error and optimism

[Chapter 19](../ch19-theory-of-departures/index.html) compared a short regression with a long one and
found that neither wins everywhere. Leaving out regressors whose coefficients are small buys a drop in
variance that can outweigh the bias it causes, and @thm-dep-mse located the break-even point exactly:
at the design points the short fit is better iff the noncentrality @eq-dep-gamma of the omitted block
is less than the number of columns omitted. The noncentrality is unknown, so the theorem does not
choose a model. It does say what a good choice should aim at, namely the error of the fitted values as
predictions, and this chapter builds that into procedures.

The residual mean square describes how well a fit reproduces the data it was fitted to, and every
fitting procedure bends towards those data. This section makes the size of that bend exact for any
procedure, and turns it into a definition of the degrees of freedom a procedure spends.

## In-sample prediction error

Throughout the chapter the data follow
\[
\Y=\bmu+\be,\qquad \E(\be)=\bzero,\qquad \Cov(\be)=\sigma^2\I_n ,
\]{#eq-sel-model}

where the mean vector \( \bmu\in\Real^n \) is unrestricted. A candidate linear model \( \C(\X) \) may
or may not contain \( \bmu \). That is the point of model selection: we want procedures whose
properties do not depend on a candidate being right. A **fitting procedure** is any function
\( \hat{\bmu}=\hat{\bmu}(\Y) \) from \( \Real^n \) to \( \Real^n \) with \( \E\norm{\hat{\bmu}(\Y)}^2<\infty \):
least squares on a fixed model, least squares on a model chosen from the data, ridge regression,
thresholding.

::: {#def-sel-prediction-error}
[Training error, prediction error, optimism]

Let \( \Y^*=\bmu+\be^* \), where \( \be^* \) is independent of \( \be \) and has the same distribution: new
responses at the same design points. For a fitting procedure \( \hat{\bmu} \):

::: {.enumerate options="label=(\alph*)"}
1. the **training error** is \( \overline{\text{err}}=n^{-1}\norm{\Y-\hat{\bmu}}^2 \);

2. the **in-sample prediction error** is
   \( \text{Err}_{\text{in}}=n^{-1}\E\bigl(\norm{\Y^*-\hat{\bmu}}^2\bigm|\Y\bigr) \);

3. the **optimism** is \( \text{Err}_{\text{in}}-\overline{\text{err}} \), and its mean
   \( \omega=\E(\text{Err}_{\text{in}}-\overline{\text{err}}) \) is the **expected optimism**.
:::

:::

Expanding \( \Y^*-\hat{\bmu}=\be^*+(\bmu-\hat{\bmu}) \) and using the independence of \( \be^* \) and \( \Y \),
\[
\text{Err}_{\text{in}}=\sigma^2+\frac1n\norm{\hat{\bmu}-\bmu}^2 .
\]{#eq-sel-err-in}

So the in-sample prediction error is the irreducible \( \sigma^2 \) plus the average squared error of the
fitted values as estimates of the mean. Ranking procedures by \( \E(\text{Err}_{\text{in}}) \) is the same as
ranking them by the risk \( \E\norm{\hat{\bmu}-\bmu}^2 \), which is the criterion of @thm-dep-mse(c).
*In-sample* refers to the design points: the responses are new, the regressor values old. New
regressor values are taken up at the end of the section.

## The covariance penalty

::: {#thm-sel-optimism}
[Optimism is a covariance]

Under @eq-sel-model, let \( \hat{\bmu} \) be a fitting procedure.

::: {.enumerate options="label=(\alph*)"}
1. The expected optimism is
   \[
   \omega=\frac2n\sum_{i=1}^n\Cov(\hat{\mu}_i,Y_i).
   \]

2. If \( \hat{\bmu}=\bS\Y \) for a fixed \( n\times n \) matrix \( \bS \), then \( \omega=2\sigma^2\tr(\bS)/n \) and
   \( \E\norm{\hat{\bmu}-\bmu}^2=\norm{(\I-\bS)\bmu}^2+\sigma^2\tr(\bS\bS\T) \).

3. For least squares on a model matrix \( \X \) of rank \( p \), whether or not \( \bmu\in\C(\X) \),
   \[
   \begin{aligned}
   \E(\overline{\text{err}})&=\sigma^2\,\frac{n-p}n+\frac{\norm{(\I-\M)\bmu}^2}n,\\
   \E(\text{Err}_{\text{in}})&=\sigma^2\,\frac{n+p}n+\frac{\norm{(\I-\M)\bmu}^2}n,
   \end{aligned}
   \]
   and \( \omega=2\sigma^2p/n \).
:::

:::

::: {.proof}
(a) Write \( \Y-\hat{\bmu}=(\Y-\bmu)-(\hat{\bmu}-\bmu) \). Then
\[
\norm{\Y-\hat{\bmu}}^2=\norm{\be}^2+\norm{\hat{\bmu}-\bmu}^2-2\be\T(\hat{\bmu}-\bmu).
\]
Take expectations and compare with @eq-sel-err-in: \( n\,\E(\text{Err}_{\text{in}})=n\sigma^2+\E\norm{\hat{\bmu}-\bmu}^2 \),
so \( n\omega=2\E\bigl[\be\T(\hat{\bmu}-\bmu)\bigr]=2\sum_i\E\bigl[(Y_i-\mu_i)\hat{\mu}_i\bigr] \), since \( \bmu \)
is constant and \( \E(\be)=\bzero \). Because \( \E(Y_i-\mu_i)=0 \), the last expectation is
\( \Cov(\hat{\mu}_i,Y_i) \). The moments exist by the Cauchy–Schwarz inequality.
(b) \( \Cov(\bS\Y,\Y)=\sigma^2\bS \), whose diagonal sums to \( \sigma^2\tr(\bS) \). For the risk,
\( \bS\Y-\bmu=\bS\be-(\I-\bS)\bmu \) has mean \( -(\I-\bS)\bmu \) and covariance \( \sigma^2\bS\bS\T \), and
@thm-rv-quadform-mean with \( \A=\I \) gives the sum of the squared mean and the trace.
(c) Take \( \bS=\M \), with \( \tr(\M)=p \) and \( \M\M\T=\M \) (@prp-proj-trace-rank). Then
\( \E(\text{Err}_{\text{in}})=\sigma^2+n^{-1}\{\norm{(\I-\M)\bmu}^2+\sigma^2p\} \), and subtracting \( \omega \) gives
\( \E(\overline{\text{err}}) \).
:::

Part (c) with \( \bmu\in\C(\X) \) is @prp-cor-optimism(a). What is new is that the optimism of least squares does not
depend on \( \bmu \): a wrong model has larger training and prediction errors, but both grow by the same squared bias, so
a penalty depending only on \( p \) is a sensible correction even for wrong models. Part (a) asks nothing of the procedure
but finite second moments: a fit is optimistic to the extent that each fitted value follows its own response.

::: {.idea}
The training error of a procedure understates its prediction error by \( 2/n \) times the total
covariance between fitted values and responses. For least squares on \( p \) columns the covariance is
\( \sigma^2p \), right model or wrong.
:::

::: {#exm-sel-wrong-mean}
[Optimism of a wrong model]

Take \( n=30 \) equally spaced points \( t_i \) on \( [-1,1] \), a curved mean \( \mu_i=1+t_i+1.5t_i^2 \),
\( \sigma=1 \), and fit by least squares a straight line in \( t \) together with two columns of pure noise, so
\( p=4 \). The model is wrong: its squared bias per observation is
\( \norm{(\I-\M)\bmu}^2/n=0.2054 \). Over \( 20000 \) simulated data sets the average optimism is
\( 0.2688 \), against \( 2\sigma^2p/n=0.2667 \) from @thm-sel-optimism(c). Ridge regression (@thm-shr-ridge) on the same
columns with penalty \( \lambda=5 \) is a linear smoother with \( \tr(\bS)=3.207 \); its simulated optimism is
\( 0.2163 \), against \( 2\times3.207/30=0.2138 \) from part (b).
:::

```{.python .run #cell-optimism-penalty}
import numpy as np
rng = np.random.default_rng(2901)
n, sigma, reps = 30, 1.0, 20_000

t = np.linspace(-1, 1, n)
X = np.column_stack([np.ones(n), t, rng.normal(size=(n, 2))])   # line + 2 irrelevant columns
mu = 1 + t + 1.5 * t**2                                         # the true mean is curved
p = X.shape[1]
M = X @ np.linalg.solve(X.T @ X, X.T)                           # projection onto C(X)

Y = mu + sigma * rng.normal(size=(reps, n))                     # one data set per row
fit = Y @ M                                                     # M is symmetric
train = np.mean((Y - fit) ** 2, axis=1)                         # training error
err_in = sigma**2 + np.mean((fit - mu) ** 2, axis=1)            # in-sample prediction error
print(f"average optimism {np.mean(err_in - train):.4f},  2 sigma^2 p / n = {2 * sigma**2 * p / n:.4f}")
```

```{.python .run #cell-optimism-ridge}
lam = 5.0
S = X @ np.linalg.solve(X.T @ X + lam * np.eye(p), X.T)         # ridge smoother
fit_r = Y @ S.T
train_r = np.mean((Y - fit_r) ** 2, axis=1)
err_in_r = sigma**2 + np.mean((fit_r - mu) ** 2, axis=1)
omega_r = np.mean(err_in_r - train_r)
print(f"ridge: optimism {omega_r:.4f},  2 sigma^2 tr(S) / n = {2 * np.trace(S) / n:.4f}")
```

## Degrees of freedom of a procedure

@thm-sel-optimism suggests measuring the complexity of any procedure by the covariance it
creates.

::: {#def-sel-df}
[Degrees of freedom]

The **degrees of freedom** of a fitting procedure \( \hat{\bmu} \) under @eq-sel-model is
\[
\text{df}(\hat{\bmu})=\frac1{\sigma^2}\sum_{i=1}^n\Cov(\hat{\mu}_i,Y_i),
\]
so that \( \omega=2\sigma^2\,\text{df}(\hat{\bmu})/n \).
:::

For least squares on a fixed model \( \text{df}=p \), the number of free coefficients, and for a linear
smoother \( \text{df}=\tr(\bS) \). For ridge regression with singular values \( d_j \) of \( \X \) this is
\( \sum_jd_j^2/(d_j^2+\lambda) \), the effective degrees of freedom of @thm-shr-ridge, which falls from
\( \rank(\X) \) to \( 0 \) as \( \lambda \) grows. If \( \text{df} \) is known,
\[
\widehat{\text{Err}}_{\text{in}}=\overline{\text{err}}+\frac{2\sigma^2}n\,\text{df}(\hat{\bmu})
\]{#eq-sel-covariance-penalty}

is an unbiased estimate of \( \E(\text{Err}_{\text{in}}) \) by @thm-sel-optimism(a). [Section 29.2](02-cp-aic-bic.html)
turns this into Mallows' \( C_p \). For a procedure that searches, the covariance is not a count and has
to be computed. Under normal errors it is computed by Stein's lemma,
@lem-shr-stein(b), proved in [Section 27.5](../ch27-shrinkage/05-james-stein.html) for the James–Stein theorem:
if \( \Z\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \) and each \( g_i \) is absolutely continuous in \( z_i \) with the
integrability stated there, then \( \E[(\Z-\boldsymbol{\uptheta})\T\mathbf{g}(\Z)]=\sigma^2\E\sum_i\partial g_i(\Z)/\partial z_i \).

::: {#prp-sel-divergence}
[Degrees of freedom as a divergence]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \). Suppose that for each \( i \), for almost every value of the other
coordinates, \( y_i\mapsto\hat{\mu}_i(\y) \) is absolutely continuous, and that
\( \E\lvert\partial\hat{\mu}_i/\partial y_i\rvert \) and \( \E\hat{\mu}_i^2 \) are finite. Then
\[
\text{df}(\hat{\bmu})=\E\sum_{i=1}^n\frac{\partial\hat{\mu}_i}{\partial y_i}(\Y),
\]
and \( \overline{\text{err}}+2\sigma^2n^{-1}\sum_i\partial\hat{\mu}_i/\partial y_i \) is an unbiased estimate of
\( \E(\text{Err}_{\text{in}}) \) that does not involve \( \bmu \).
:::

::: {.proof}
Since \( \E(Y_i-\mu_i)=0 \), \( \sum_i\Cov(\hat{\mu}_i,Y_i)=\E[(\Y-\bmu)\T\hat{\bmu}(\Y)] \), and each term is integrable by the
Cauchy–Schwarz inequality. @lem-shr-stein(b) with \( \mathbf{g}=\hat{\bmu} \) turns this into
\( \sigma^2\E\sum_i\partial\hat{\mu}_i/\partial y_i \). The unbiasedness statement is then @thm-sel-optimism(a).
:::

The estimate in the proposition is **Stein's unbiased risk estimate** (Stein 1981). For a linear smoother
the divergence is \( \tr(\bS) \), and nothing is new. The interesting cases are adaptive.

**Soft thresholding.** Take an orthonormal design, \( \X\T\X=\I_p \), and put \( \Z=\X\T\Y\sim\Normal_p(\X\T\bmu,\sigma^2\I) \).
For a fit \( \hat{\bmu}=\X\hat{\bbeta}(\Z) \),
\( \sum_i\Cov(\hat{\mu}_i,Y_i)=\tr\Cov(\X\hat{\bbeta},\Y)=\tr\Cov(\hat{\bbeta},\X\T\Y)=\sum_j\Cov(\hat{\beta}_j,Z_j) \), so the
degrees of freedom can be computed coordinate by coordinate in the sequence model
\( Z_j\sim\Normal(\theta_j,\sigma^2) \). The lasso in this design is soft thresholding,
\( \hat{\beta}_j=\operatorname{sign}(Z_j)(\lvert Z_j\rvert-\lambda)_+ \) (@thm-shr-lasso-orthonormal). It is Lipschitz,
with derivative \( 1\{\lvert Z_j\rvert>\lambda\} \) except at \( \pm\lambda \), so by @prp-sel-divergence its
degrees of freedom are the expected number of nonzero coefficients, as an exercise of
[Section 27.4](../ch27-shrinkage/04-lasso.html) found by direct integration. The shrinkage exactly pays for the search.

**Hard thresholding.** Subset selection in the same design keeps or drops each coordinate whole,
\( \hat{\beta}_j=Z_j1\{\lvert Z_j\rvert>\lambda\} \) (@thm-shr-lasso-orthonormal again). This jumps by \( \lambda \) at \( \pm\lambda \),
so Stein's lemma does not apply, and the degrees of freedom exceed the expected count.

::: {#prp-sel-hard-df}
[Degrees of freedom of hard thresholding]

Let \( Z\sim\Normal(\theta,\sigma^2) \) and \( \lambda>0 \). Then
\[
\begin{aligned}
\frac{\Cov\bigl(Z1\{\lvert Z\rvert>\lambda\},Z\bigr)}{\sigma^2}
&=\Pr(\lvert Z\rvert>\lambda)\\
&\quad+\frac\lambda\sigma\Bigl[\phi\Bigl(\frac{\lambda-\theta}\sigma\Bigr)+\phi\Bigl(\frac{\lambda+\theta}\sigma\Bigr)\Bigr],
\end{aligned}
\]
where \( \phi \) is the standard normal density.
:::

::: {.proof}
Write \( Z=\theta+\sigma W \) with \( W\sim\Normal(0,1) \), and put \( a=(\lambda-\theta)/\sigma \), \( b=(-\lambda-\theta)/\sigma \),
so that \( \lvert Z\rvert>\lambda \) iff \( W>a \) or \( W<b \). The covariance is
\( \E[\sigma W(\theta+\sigma W)1\{W>a\text{ or }W<b\}] \). The truncated moments of the standard normal are
\( \E[W1\{W>a\}]=\phi(a) \), \( \E[W1\{W<b\}]=-\phi(b) \), \( \E[W^21\{W>a\}]=1-\Phi(a)+a\phi(a) \) and
\( \E[W^21\{W<b\}]=\Phi(b)-b\phi(b) \), each by one integration by parts. Collecting terms,
\[
\sigma\theta\bigl(\phi(a)-\phi(b)\bigr)+\sigma^2\bigl(1-\Phi(a)+\Phi(b)+a\phi(a)-b\phi(b)\bigr).
\]
Since \( \sigma\theta+\sigma^2a=\sigma\lambda \) and \( -\sigma\theta-\sigma^2b=\sigma\lambda \), this is
\( \sigma^2\Pr(\lvert Z\rvert>\lambda)+\sigma\lambda\bigl(\phi(a)+\phi(b)\bigr) \), and \( \phi(b)=\phi((\lambda+\theta)/\sigma) \).
:::

The first term is the expected number of coefficients kept. The second is the price of the search:
the fit follows \( Z \) most closely just where \( Z \) decides whether to keep the coordinate. At the
threshold \( \lambda=\sigma\sqrt2 \) (which [Section 29.2](02-cp-aic-bic.html) shows is where \( C_p \) and AIC
put it) and \( \theta=0 \), the expected count is \( 0.157 \) but the degrees of freedom are \( 0.572 \): the search
term \( 0.415 \) is more than twice the count. At \( \theta=\sigma \) the degrees of freedom are \( 0.896 \), and at
\( \theta=3\sigma \) they are \( 1.104 \), more than one full parameter for a single coordinate.

::: {#exm-sel-subset-df}
[Best subset of a given size]

Keep the \( k \) largest of \( m=20 \) values \( \lvert Z_j\rvert \) and set the rest to zero, with \( Z_j\sim\Normal(\theta_j,1) \)
independent. For a subset fixed in advance the degrees of freedom would be \( k \). [Figure 29.1.1](#fig-sel-subset-df)
shows them by simulation. When all \( \theta_j=0 \), the best single coordinate already costs \( 4.91 \) degrees of
freedom, the best five cost \( 13.94 \) and the best ten \( 18.27 \). When five means equal \( 6 \) and the rest are
zero, choosing the best five is almost free (\( 5.06 \)), because the five are obvious; choosing a sixth
means choosing among fifteen pure-noise coordinates, and the degrees of freedom jump to
\( 9.47 \). Searching among many candidates of equal merit is expensive, and choosing between clear winners and
clear losers is cheap.
:::

::: {when-format="html"}
![**Figure 29.1.1.** Degrees of freedom of the best subset of size \( k \) among \( 20 \) orthonormal
coordinates, by simulation (\( 20000 \) data sets). The dashed line is \( \text{df}=k \), correct only
when the subset is fixed in advance.](subset_df.svg){#fig-sel-subset-df width=70%}
:::

::: {when-format="pdf"}
![Degrees of freedom of the best subset of size \( k \) among \( 20 \) orthonormal
coordinates, by simulation (\( 20000 \) data sets). The dashed line is \( \text{df}=k \), correct only
when the subset is fixed in advance.](subset_df.pdf){width=70%}
:::

```{.python .run #cell-optimism-subset}
m, reps_s = 20, 20_000
Zs = rng.normal(size=(reps_s, m))


def subset_df(theta, k):
    """Monte Carlo df of 'keep the k largest |Y_j|' when Y ~ N(theta, I)."""
    Yb = theta + Zs
    keep = np.argsort(-np.abs(Yb), axis=1)[:, :k]
    fit = np.zeros_like(Yb)
    np.put_along_axis(fit, keep, np.take_along_axis(Yb, keep, axis=1), axis=1)
    return np.mean(np.sum(fit * Zs, axis=1))                   # sum_j Cov(fit_j, Y_j)


null = np.zeros(m)
sparse = np.r_[np.full(5, 6.0), np.zeros(m - 5)]               # five large means
for k in (1, 5, 10):
    print(f"k = {k:2d}: df null {subset_df(null, k):5.2f}, df sparse {subset_df(sparse, k):5.2f}")
```

So a penalty of \( 2\sigma^2p/n \) is right for a model chosen in advance and too small for a model chosen by
search, a point [Section 29.5](05-selection-bias.html) returns to.

## New regressor values

The in-sample error keeps the design fixed. For a new case with its own regressors the error is larger, and for least
squares with \( k \) normal regressors and a correct model @thm-cor-prediction-error gives it exactly; since
\( \E(s^2)=\sigma^2 \), \( s^2(1+1/n)(n-2)/(n-k-2) \) estimates it without bias. Breiman and Spector (1992) found that
fixed-design estimates understate the error when the design is random, and recommended five- or ten-fold
cross-validation. Cross-validation
([Section 29.3](03-cross-validation.html)) estimates the random-design error directly.

## Exercises

### A. Check your understanding

::: {#exr-sel-interpolate}
[A1]

The procedure \( \hat{\bmu}(\Y)=\Y \) interpolates the data. Compute its training error, its in-sample
prediction error and its degrees of freedom, and check @thm-sel-optimism(a).
:::

### B. Practice

::: {#exr-sel-ridge-df}
[B1]

Let \( \bS_\lambda=\X(\X\T\X+\lambda\I)^{-1}\X\T \) be the ridge smoother, whose trace
\( \sum_jd_j^2/(d_j^2+\lambda) \) was found in @thm-shr-ridge. Show that
\( \tr(\bS_\lambda\bS_\lambda\T)<\tr(\bS_\lambda) \) for every \( \lambda>0 \) with \( \X\ne\mathbf{0} \), and that for a symmetric
smoother \( \bS \) with eigenvalues in \( [0,1] \), \( \tr(\bS\bS\T)=\tr(\bS) \) iff \( \bS \) is an orthogonal projection.
:::

::: {.solution}
A symmetric \( \bS \) with eigenvalues \( s_j\in[0,1] \) has \( \tr(\bS\bS\T)=\tr(\bS^2)=\sum_js_j^2\le\sum_js_j=\tr(\bS) \), with
equality iff every \( s_j\in\{0,1\} \), that is, iff \( \bS \) is idempotent, an orthogonal projection (@thm-proj-sym-idem).
The ridge smoother is symmetric with eigenvalues \( d_j^2/(d_j^2+\lambda) \), which lie in \( (0,1) \) for \( d_j>0 \), and
at least one \( d_j \) is positive. So the variance term of @thm-sel-optimism(b), \( \sigma^2\tr(\bS\bS\T) \), is
strictly smaller than the optimism term \( \sigma^2\tr(\bS) \).
:::

::: {#exr-sel-moving-average}
[B2]

Arrange \( n \) observations on a circle, and let \( \hat{\mu}_i \) be the average of \( Y_{i-q},\dots,Y_{i+q} \)
(indices modulo \( n \), \( 2q+1\le n \)). Find the degrees of freedom and the expected optimism, and explain why
a wider window is "simpler".
:::

### C. Going deeper

::: {#exr-sel-sure-choose}
[C1]

In the sequence model \( Z_j\sim\Normal(\theta_j,\sigma^2) \), \( j=1,\dots,p \), independent, soft thresholding at
\( \lambda \) has Stein's unbiased risk estimate
\[
\text{SURE}(\lambda)=\sum_{j=1}^p\min(Z_j^2,\lambda^2)-p\sigma^2+2\sigma^2\#\{j:\lvert Z_j\rvert>\lambda\}
\]
by @prp-sel-divergence, because \( Z_j-\hat{\theta}_j \) is \( Z_j \) when \( \lvert Z_j\rvert\le\lambda \) and \( \pm\lambda \)
otherwise. Let \( 0=a_0\le a_1\le\dots\le a_p \) be \( 0 \) and the ordered \( \lvert Z_j\rvert \). Show that
\( \text{SURE} \) is nondecreasing on each interval \( [a_k,a_{k+1}) \), so it is minimized over \( \lambda\ge0 \) at one of
\( a_0,\dots,a_p \). Compute \( \text{SURE}(0) \), and explain why the minimized value is biased downwards.
(Minimizing SURE is the SureShrink rule of Donoho and Johnstone 1995.)
:::

::: {.solution}
On \( [a_k,a_{k+1}) \) exactly \( p-k \) of the \( \lvert Z_j\rvert \) exceed \( \lambda \) (for distinct values), so
\( \text{SURE}(\lambda)=\sum_{\lvert Z_j\rvert\le a_k}Z_j^2+(p-k)\lambda^2-p\sigma^2+2\sigma^2(p-k) \), which is nondecreasing in
\( \lambda \) and smallest at \( a_k \). At \( \lambda=0 \) nothing is shrunk and \( \text{SURE}(0)=p\sigma^2 \), the exact risk of
\( \Z \). Each \( \text{SURE}(a_k) \) is unbiased for its own risk, so their minimum has expectation at most the smallest
risk: the optimism of a minimized estimate studied in [Section 29.5](05-selection-bias.html).
:::

::: {#exr-sel-max-df}
[C2]

In @exm-sel-subset-df with all \( \theta_j=0 \) and \( k=1 \), show that the degrees of freedom equal
\( \E\max_jZ_j^2 \) (with \( \sigma=1 \)). Deduce that they are at least \( 1 \), and use
\( \E\max_jZ_j^2\ge(\E\max_j\lvert Z_j\rvert)^2 \) to explain why they grow without bound as \( m\to\infty \).
:::

