# Resampling cases versus resampling residuals

Let \( \theta=\mathbf{c}\T\bbeta \) be estimated by \( \hat{\theta}=\mathbf{c}\T\hbeta \). An interval or a test for \( \theta \)
needs the distribution of a **root** such as \( \hat{\theta}-\theta \), or a studentized version of it. That distribution
depends on the mechanism \( P \) that produced the data. In Part III, \( P \) lay in the normal linear model, where the
studentized root is \( t(n-r) \) whatever the parameters. Without that assumption the distribution is unknown.

The **bootstrap** (Efron 1979) replaces \( P \) by an estimate \( \hat{P} \). Draw a data set \( \Y^* \) from \( \hat{P} \),
compute \( \hat{\theta}^* \) from it, and use the distribution of \( \hat{\theta}^*-\hat{\theta} \), given the data, as the
approximation to that of \( \hat{\theta}-\theta \). Here \( \hat{\theta} \) plays the true value, because it is the value of
\( \theta \) under \( \hat{P} \). Probabilities, expectations and covariances under \( \hat{P} \) with the data held fixed are
written \( {\Pr}_* \), \( \E_* \) and \( \Cov_* \). In practice \( B \) resamples approximate the bootstrap distribution. Their
Monte Carlo error shrinks as \( B \) grows; the statistical error, from the gap between \( \hat{P} \) and \( P \), does not, and
it is the subject of this chapter.

## Two ways to estimate the data-generating mechanism

In the **fixed-design model** the rows \( \x_{(i)}\T \) of \( \X \) are constants and \( \Y=\X\bbeta+\be \) with independent,
identically distributed errors from an unknown law \( F \). Estimate \( \bbeta \) by \( \hbeta \) and \( F \) by the empirical
distribution of the residuals. In the **random-design model** the cases \( (\x_{(i)},Y_i) \) are independent draws from an
unknown joint distribution. Estimate it by the empirical distribution of the \( n \) cases.

::: {#def-bs-bootstrap}
[Residual and case bootstraps]

Let \( \X \) be \( n\times p \) of rank \( p \), with least squares estimate \( \hbeta \) and residuals
\( \hat{\varepsilon}_1,\dots,\hat{\varepsilon}_n \) with mean \( \bar{\hat{\varepsilon}} \).

::: {.enumerate options="label=(\alph*)"}
1. The **residual bootstrap** draws \( \varepsilon^*_1,\dots,\varepsilon^*_n \) independently and uniformly
   from the centred residuals \( \tilde{\varepsilon}_i=\hat{\varepsilon}_i-\bar{\hat{\varepsilon}} \), sets
   \( \Y^*=\X\hbeta+\be^* \), and computes \( \hbeta^*=(\X\T\X)^{-1}\X\T\Y^* \). The design is the same in every resample.

2. The **case bootstrap** (or pairs bootstrap) draws \( n \) indices \( I_1,\dots,I_n \) independently and uniformly
   from \( \{1,\dots,n\} \) and computes the least squares estimate \( \hbeta^* \) from the cases \( (\x_{(I_k)},Y_{I_k}) \),
   \( k=1,\dots,n \). The design changes from resample to resample.
:::

Variants of (a) resample the **leverage-adjusted residuals**
\( r_i=\hat{\varepsilon}_i/\sqrt{1-h_{ii}} \), centred, or the residuals rescaled by \( \sqrt{n/(n-p)} \).
:::

With an intercept the residuals already sum to zero (@prp-proj-fit-residual(b)). Without one, centring prevents a bias
in the bootstrap world (@exr-bs-no-intercept).

## The moments of the residual bootstrap

::: {#prp-bs-residual-moments}
[Moments of the residual bootstrap]

In the residual bootstrap, put \( \hat{\sigma}_*^2=n^{-1}\sum_i\tilde{\varepsilon}_i^2 \), the variance of the resampled law. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \E_*(\hbeta^*)=\hbeta \) and \( \Cov_*(\hbeta^*)=\hat{\sigma}_*^2(\X\T\X)^{-1} \);

2. if \( \bone\in\C(\X) \), then \( \hat{\sigma}_*^2=\text{SSE}/n=(n-p)s^2/n \);

3. if the data follow the linear model with \( \Cov(\be)=\sigma^2\I \), then \( \E(r_i^2)=\sigma^2 \) for each \( i \) with \( h_{ii}<1 \).
:::

:::

::: {.proof}
(a) Given the data, \( \hbeta^*=\hbeta+(\X\T\X)^{-1}\X\T\be^* \), where the \( \varepsilon^*_i \) are independent with
\( \E_*(\varepsilon^*_i)=n^{-1}\sum_j\tilde{\varepsilon}_j=0 \) and \( \Var_*(\varepsilon^*_i)=\hat{\sigma}_*^2 \). Apply @thm-lm-moments in the
bootstrap world. (b) The residuals sum to zero, so \( \tilde{\varepsilon}_i=\hat{\varepsilon}_i \). (c) \( \he=(\I-\M)\be \) by
@prp-proj-fit-residual(e), so \( \E(\hat{\varepsilon}_i)=0 \) and \( \Var(\hat{\varepsilon}_i)=\sigma^2(1-h_{ii}) \).
:::

By (b) the plain residual bootstrap shrinks the usual variance by the factor \( (n-p)/n \): residuals are shorter than errors
because least squares fits part of the noise (@thm-lm-sigma2). The factor matters only when \( p/n \) is not small
([Section 23.6](06-failures.html)). Rescaling by \( \sqrt{n/(n-p)} \) removes it exactly when the model has an
intercept. The leverage-adjusted residuals remove it observation by observation. They are \( s \) times the internally
studentized residuals of [Chapter 20](../ch20-residuals-leverage-influence/index.html) (@def-res-residuals).

## Consistency for a linear contrast

Consider designs \( \X_n \) of rank \( p \), with \( p \) fixed, and errors
\( \varepsilon_{n1},\dots,\varepsilon_{nn} \) drawn independently from a fixed law \( F \) with mean zero and variance
\( \sigma^2\in(0,\infty) \). For \( \mathbf{c}_n\ne\bzero \) put
\[
\begin{aligned}
&\mathbf{a}_n=\X_n(\X_n\T\X_n)^{-1}\mathbf{c}_n,\qquad
\mathbf{c}_n\T(\hbeta-\bbeta)=\mathbf{a}_n\T\be,\\
&\norm{\mathbf{a}_n}^2=\mathbf{c}_n\T(\X_n\T\X_n)^{-1}\mathbf{c}_n ,
\end{aligned}
\]
so the contrast estimate is a weighted sum of the errors with standard deviation \( \sigma\norm{\mathbf{a}_n} \). The condition we
need is that no single weight dominates:
\[
m_n=\max_i\frac{\lvert a_{ni}\rvert}{\norm{\mathbf{a}_n}}\to0 .
\]{#eq-bs-no-dominant}

The Cauchy–Schwarz bound \( a_{ni}^2\le h_{ii}\norm{\mathbf{a}_n}^2 \) in the proof of @thm-dep-nonnormal shows that
@eq-bs-no-dominant holds for every contrast when the largest leverage tends to zero. The limit theorem used there is the
Lindeberg–Feller theorem of [Chapter 19](../ch19-theory-of-departures/index.html); for a triangular array \( U_{ni} \) its condition holds in
particular when \( \sum_i\E U_{ni}^4\to0 \) (Lyapunov), since the Lindeberg sum is at most \( \eta^{-2}\sum_i\E U_{ni}^4 \). We add one more
result from probability, used without proof (van der Vaart 1998, chapter 2).

::: {.remark}
[Pólya's theorem]

If distribution functions \( G_n \) converge pointwise to a continuous distribution function \( G \), then
\( \sup_x\lvert G_n(x)-G(x)\rvert\to0 \).
:::

A bootstrap distribution is a conditional distribution that changes with the data, so we need a conditional form of the
Lindeberg–Feller theorem.

::: {#lem-bs-conditional-clt}
[A conditional central limit theorem]

For each \( n \), given the data \( \mathcal{D}_n \), let \( U^*_{n1},\dots,U^*_{nn} \) be independent with
\( \E_*(U^*_{ni})=0 \) and \( \sum_i\E_*(U^{*2}_{ni})=1 \). If for every \( \eta>0 \)
\[
L^*_n(\eta)=\sum_i\E_*\bigl[U^{*2}_{ni}\,\mathbf{1}\{\lvert U^*_{ni}\rvert>\eta\}\bigr]\to0\quad\text{in probability},
\]
then \( \Delta_n=\sup_x\bigl\lvert{\Pr}_*\bigl(\sum_iU^*_{ni}\le x\bigr)-\Phi(x)\bigr\rvert\to0 \) in probability.
:::

::: {.proof}
Random variables converge to zero in probability if and only if every subsequence has a further subsequence along which they
converge almost surely. Given a subsequence, a diagonal argument gives a further subsequence along which
\( L^*_n(1/k)\to0 \) almost surely for every \( k=1,2,\dots \) at once. Since \( L^*_n(\eta) \) decreases in \( \eta \), on this event
\( L^*_n(\eta)\to0 \) for every \( \eta>0 \). For each outcome in the event the conditional laws satisfy the Lindeberg condition, so
\( \Delta_n\to0 \) by the Lindeberg–Feller and Pólya theorems.
:::

::: {#lem-bs-residual-tails}
[Residuals have uniformly small tails]

In the setting above, \( \hat{\sigma}_*^2\to\sigma^2 \) in probability, and for every sequence of constants \( K_n\to\infty \),
\( n^{-1}\sum_i\tilde{\varepsilon}_i^2\,\mathbf{1}\{\lvert\tilde{\varepsilon}_i\rvert>K_n\}\to0 \) in probability.
:::

::: {.proof}
Write \( \tilde{\varepsilon}_i=\varepsilon_i-\delta_i \) with \( \boldsymbol{\updelta}=\M\be+\bar{\hat{\varepsilon}}\bone \).
Here \( \E\norm{\M\be}^2=\sigma^2p \), and \( \bar{\hat{\varepsilon}}=n^{-1}\bone\T(\I-\M)\be \) has
\( \E(\bar{\hat{\varepsilon}}^2)=\sigma^2n^{-2}\bone\T(\I-\M)\bone\le\sigma^2/n \). So
\( \E\norm{\boldsymbol{\updelta}}^2\le2\sigma^2(p+1) \), and \( n^{-1}\norm{\boldsymbol{\updelta}}^2\to0 \) in probability.

First claim: \( \hat{\sigma}_*^2=n^{-1}\norm{\be}^2-n^{-1}\norm{\M\be}^2-\bar{\hat{\varepsilon}}^2 \). The first term tends to
\( \sigma^2 \) by the weak law of large numbers (the errors of the \( n \)th model are an independent sample from \( F \)), and
the others tend to zero.

Second claim: for \( K>0 \) and reals \( u,d \), with \( v=u-d \),
\[
v^2\,\mathbf{1}\{\lvert v\rvert>K\}\le2u^2\,\mathbf{1}\{\lvert u\rvert>K/2\}+4d^2 .
\]
If \( \lvert u\rvert>K/2 \), use \( v^2\le2u^2+2d^2 \). If \( \lvert u\rvert\le K/2 \) and \( \lvert v\rvert>K \), then
\( \lvert d\rvert\ge\lvert v\rvert-\lvert u\rvert>K/2\ge\lvert u\rvert \), so \( \lvert v\rvert\le2\lvert d\rvert \). Otherwise the
left side is zero. With \( u=\varepsilon_i \), \( d=\delta_i \) and \( K=K_n \), averaging gives
\[
\begin{aligned}
\frac1n\sum_i\tilde{\varepsilon}_i^2\,\mathbf{1}\{\lvert\tilde{\varepsilon}_i\rvert>K_n\}
&\le\frac2n\sum_i\varepsilon_i^2\,\mathbf{1}\{\lvert\varepsilon_i\rvert>K_n/2\}\\
&\quad+\frac4n\norm{\boldsymbol{\updelta}}^2 .
\end{aligned}
\]
The first term has expectation \( 2\E\bigl[\varepsilon^2\mathbf{1}\{\lvert\varepsilon\rvert>K_n/2\}\bigr]\to0 \) by dominated
convergence, and the second tends to zero in probability.
:::

::: {#thm-bs-residual}
[Consistency of the residual bootstrap]

In the fixed-design setting above, assume @eq-bs-no-dominant. Let
\[
\begin{aligned}
G_n(x)&=\Pr\bigl(\mathbf{c}_n\T(\hbeta-\bbeta)\le x\norm{\mathbf{a}_n}\bigr),\\
\hat{G}_n(x)&={\Pr}_*\bigl(\mathbf{c}_n\T(\hbeta^*-\hbeta)\le x\norm{\mathbf{a}_n}\bigr),
\end{aligned}
\]
where \( \hbeta^* \) comes from the residual bootstrap. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \sup_x\lvert G_n(x)-\Phi(x/\sigma)\rvert\to0 \);

2. \( \sup_x\lvert\hat{G}_n(x)-\Phi(x/\sigma)\rvert\to0 \) in probability;

3. \( \sup_x\lvert\hat{G}_n(x)-G_n(x)\rvert\to0 \) in probability.
:::

:::

::: {.proof}
Put \( w_i=a_{ni}/\norm{\mathbf{a}_n} \), so \( \sum_iw_i^2=1 \) and \( \max_i\lvert w_i\rvert=m_n\to0 \).

(a) \( \mathbf{c}_n\T(\hbeta-\bbeta)/\norm{\mathbf{a}_n}=\sum_iw_i\varepsilon_i \). When \( \max_ih_{ii}\to0 \) the claim is
@thm-dep-nonnormal(a) with Pólya's theorem. The same Lindeberg argument needs only @eq-bs-no-dominant: with
\( U_{ni}=w_i\varepsilon_i/\sigma \), the Lindeberg sum is at most
\( \E\bigl[\sigma^{-2}\varepsilon^2\mathbf{1}\{\lvert\varepsilon\rvert>\eta\sigma/m_n\}\bigr]\to0 \), by dominated convergence.

(b) Similarly \( \mathbf{c}_n\T(\hbeta^*-\hbeta)/\norm{\mathbf{a}_n}=\sum_iw_i\varepsilon^*_i \). Take
\( U^*_{ni}=w_i\varepsilon^*_i/\hat{\sigma}_* \). As in (a),
\[
L^*_n(\eta)\le\frac{1}{n\hat{\sigma}_*^2}\sum_i\tilde{\varepsilon}_i^2\,\mathbf{1}\{\lvert\tilde{\varepsilon}_i\rvert>\eta\hat{\sigma}_*/m_n\}.
\]
By @lem-bs-residual-tails, \( \hat{\sigma}_*\to\sigma \) in probability. On the event \( \hat{\sigma}_*>\sigma/2 \), whose
probability tends to one, the right side is at most \( 4\sigma^{-2}n^{-1}\sum_i\tilde{\varepsilon}_i^2\mathbf{1}\{\lvert\tilde{\varepsilon}_i\rvert>K_n\} \)
with \( K_n=\eta\sigma/(2m_n)\to\infty \), which tends to zero in probability by the lemma. So @lem-bs-conditional-clt gives
\( \sup_x\lvert\hat{G}_n(x)-\Phi(x/\hat{\sigma}_*)\rvert\to0 \) in probability. Finally
\( \sup_x\lvert\Phi(x/\hat{\sigma}_*)-\Phi(x/\sigma)\rvert \) is a continuous function of \( \hat{\sigma}_*/\sigma \) that vanishes at
\( 1 \), so it tends to zero in probability.

(c) follows from (a) and (b) by the triangle inequality.
:::

The theorem holds for the rescaled residuals too, since they differ by the factor \( \sqrt{n/(n-p)}\to1 \); for the
leverage-adjusted ones the proof adapts when \( \max_ih_{ii}\to0 \), with small changes we omit. It needs neither normal errors nor
moments beyond the second. Freedman (1981) proved it in a stronger metric, together with a statement for the case bootstrap.
It does not say that the bootstrap beats the normal approximation \( \Phi(x/s) \), which is consistent under the same
conditions. The bootstrap earns its keep through studentization ([Section 23.3](03-confidence-intervals.html)) and through its
flexibility for statistics with no convenient standard error.

::: {.warning}
The residual bootstrap treats the residuals as exchangeable. It gives \( \Cov_*(\hbeta^*)=\hat{\sigma}_*^2(\X\T\X)^{-1} \) whatever
the data look like, which is the classical formula, so it inherits that formula's failure under heteroscedasticity (@thm-dep-covariance, @exr-bs-heteroscedastic-residual). Wu (1986) made this point forcefully.
:::

## The case bootstrap

The case bootstrap assumes only that the cases are independent draws from a population. By @prp-proj-consistency, least squares
then estimates the population projection coefficient whether or not the regression function is linear, and the case bootstrap
estimates the variability of \( \hbeta \) around that target.

Let \( W_i \) be the number of times case \( i \) is drawn, so that \( (W_1,\dots,W_n) \) is multinomial with \( n \) trials and equal
probabilities, and put \( \mathbf{W}=\diag(W_i) \). When \( \X\T\mathbf{W}\X \) is invertible,
\( \hbeta^*=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\y \), and since \( \y=\X\hbeta+\he \),
\[
\hbeta^*-\hbeta=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\he .
\]{#eq-bs-case-exact}

The matrix \( n^{-1}\X\T\mathbf{W}\X \) has expectation \( n^{-1}\X\T\X \) and is close to it for large \( n \). Replacing it by its
expectation gives a linear approximation whose covariance is exact.

::: {#prp-bs-case-linear}
[The case bootstrap, linearized, is the sandwich]

Let \( \mathbf{L}^*=(\X\T\X)^{-1}\X\T\mathbf{W}\he \). Then \( \E_*(\mathbf{L}^*)=\bzero \) and
\[
\Cov_*(\mathbf{L}^*)=(\X\T\X)^{-1}\Bigl(\sum_i\hat{\varepsilon}_i^2\,\x_{(i)}\x_{(i)}\T\Bigr)(\X\T\X)^{-1},
\]
the sandwich covariance estimator HC0 of [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (@thm-het-sandwich).
:::

::: {.proof}
For multinomial counts with \( n \) trials and probabilities \( 1/n \), \( \E(W_i)=1 \), \( \Var(W_i)=1-1/n \) and
\( \Cov(W_i,W_j)=-1/n \) for \( i\ne j \). So \( \E_*(\mathbf{W}\he)=\he \) and
\( \Cov_*(\mathbf{W}\he)=\diag(\he)(\I-n^{-1}\bone\bone\T)\diag(\he) \). Hence \( \E_*(\mathbf{L}^*)=(\X\T\X)^{-1}\X\T\he=\bzero \) by
@prp-proj-fit-residual(a), and
\( \Cov_*(\X\T\mathbf{W}\he)=\sum_i\hat{\varepsilon}_i^2\x_{(i)}\x_{(i)}\T-n^{-1}(\X\T\he)(\X\T\he)\T=\sum_i\hat{\varepsilon}_i^2\x_{(i)}\x_{(i)}\T \).
:::

To first order, then, the case bootstrap estimates the variance of \( \hbeta \) by the sandwich, not the classical formula.
The price is a varying design: a case of high leverage appears in some resamples several times and in others not at all, and
with indicator columns some resamples are singular (@exr-bs-case-singular).

::: {#exm-bs-engel}
[Engel's food expenditure]

Return to Engel's data on the income and food expenditure of \( 235 \) households (@exm-proj-engel). The least squares
slope is \( 0.4852 \), and the scatter fans out ([Figure 23.1.1](#fig-bs-engel)(a)): the residual standard deviation is
\( 2.36 \) times larger for the richer half of the households.

The classical standard error of the slope is \( 0.01437 \). The residual bootstrap with \( B=4000 \) resamples
gives \( 0.01433 \), near its exact value \( 0.01431 \), the classical value times
\( \sqrt{(n-2)/n}=0.9957 \) (@prp-bs-residual-moments); with leverage-adjusted residuals the exact value is
\( 0.01484 \). The case bootstrap gives \( 0.05032 \), more than three times as large and close to the HC0 value
\( 0.05177 \) of @prp-bs-case-linear; [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) computed
this and the other sandwich standard errors for the same slope. The residual bootstrap's answer rests on a constant
variance that the data contradict.
:::

::: {when-format="html"}
![**Figure 23.1.1.** (a) Engel's households, with the least squares line. (b) Bootstrap distributions of the slope: residual (grey), case
(blue) and wild (orange) resampling.](engel_bootstrap.svg){#fig-bs-engel width=100%}
:::

::: {when-format="pdf"}
![(a) Engel's households, with the least squares line. (b) Bootstrap distributions of the slope: residual (grey), case
(blue) and wild (orange) resampling.](engel_bootstrap.pdf){width=100%}
:::

```{.python .run #cell-engel-bootstrap-setup}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.engel.load_pandas().data

y = data["foodexp"].to_numpy()
x = data["income"].to_numpy()
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape
XtX_inv = np.linalg.inv(X.T @ X)
beta_hat = XtX_inv @ X.T @ y
e_hat = y - X @ beta_hat
h = np.sum((X @ XtX_inv) * X, axis=1)            # leverages
s2 = e_hat @ e_hat / (n - p)
se_classical = np.sqrt(s2 * XtX_inv[1, 1])
print(f"slope {beta_hat[1]:.4f}, classical se {se_classical:.5f}, max leverage {h.max():.3f}")

rng = np.random.default_rng(2301)
B = 4000
A = XtX_inv @ X.T                                 # beta_hat = A y
fit = X @ beta_hat
e_c = e_hat - e_hat.mean()                        # centred residuals
r = e_hat / np.sqrt(1 - h)                        # leverage-adjusted residuals
r_c = r - r.mean()

def residual_bootstrap(resid, B):
    """Slopes from y* = X beta_hat + e*, with e* drawn with replacement from resid."""
    idx = rng.integers(0, n, size=(B, n))
    y_star = fit + resid[idx]
    return (y_star @ A.T)[:, 1]

slopes_res = residual_bootstrap(e_c, B)
slopes_mod = residual_bootstrap(r_c, B)
print(f"residual bootstrap se: raw {slopes_res.std():.5f}, leverage-adjusted {slopes_mod.std():.5f}")

def pairs_bootstrap(B):
    """Slopes refitted to n cases drawn with replacement from the (x_i, y_i)."""
    idx = rng.integers(0, n, size=(B, n))
    xs, ys = x[idx], y[idx]
    xc = xs - xs.mean(axis=1, keepdims=True)
    return np.sum(xc * ys, axis=1) / np.sum(xc ** 2, axis=1)

slopes_pairs = pairs_bootstrap(B)
print(f"case bootstrap se {slopes_pairs.std():.5f}")
```

The choice between the two is a choice of model. The residual bootstrap suits a linear mean with exchangeable errors, even
for sampled cases: when the error law does not depend on \( \X \), conditioning on the design is justified as in
[Section 14.1](../ch14-correlation-lack-of-fit-prediction/01-random-regressors.html), where @thm-cor-conditional makes it
exact under normality. It is then the more accurate, since it adds no noise from a varying design. Under heteroscedasticity or
a nonlinear mean the case bootstrap is safer. The wild bootstrap of the next section keeps the design fixed and still respects
heteroscedasticity.

## Exercises

### A. Check your understanding

::: {#exr-bs-no-intercept}
[A1]

In a regression through the origin, \( \E(Y_i)=\beta x_i \), the residuals need not sum to zero. Show that if the
residual bootstrap resampled the uncentred residuals, then \( \E_*(\hat{\beta}^*)-\hat{\beta}=\bar{\hat{\varepsilon}}\sum_ix_i/\sum_ix_i^2 \).
:::

::: {.solution}
Here \( \hat{\beta}^*-\hat{\beta}=\sum_ix_i\varepsilon^*_i/\sum_ix_i^2 \), and each uncentred draw has
\( \E_*(\varepsilon^*_i)=n^{-1}\sum_j\hat{\varepsilon}_j=\bar{\hat{\varepsilon}} \). Taking expectations gives the stated bias, which is
not zero unless \( \bar{\hat{\varepsilon}}=0 \) or \( \sum_ix_i=0 \).
:::

::: {#exr-bs-distinct}
[A2]

A case-bootstrap resample of size \( n \) contains each original case with probability \( 1-(1-1/n)^n \). Show that the
expected number of distinct cases in a resample is \( n\bigl(1-(1-1/n)^n\bigr) \), about \( 0.632n \) for large \( n \).
:::

### B. Practice

::: {#exr-bs-heteroscedastic-residual}
[B1]

Suppose the errors are independent with variances \( \sigma_i^2 \), and the model has an intercept. Show that the residual
bootstrap variance of \( \mathbf{c}\T\hbeta^* \) has expectation \( \norm{\mathbf{a}}^2\,n^{-1}\sum_i(1-h_{ii})\sigma_i^2 \),
whereas \( \Var(\mathbf{c}\T\hbeta)=\sum_ia_i^2\sigma_i^2 \). Construct a two-group example in which the ratio of the two is
as small as you like.
:::

::: {.solution}
By @prp-bs-residual-moments the bootstrap variance is \( \norm{\mathbf{a}}^2\text{SSE}/n \), and
\( \E(\text{SSE})=\tr\bigl((\I-\M)\bSigma\bigr)=\sum_i(1-h_{ii})\sigma_i^2 \) for diagonal \( \bSigma \). The true variance is
\( \mathbf{a}\T\bSigma\mathbf{a} \). Take two groups of sizes \( n_1 \) and \( n_2 \) and the contrast of the group means, so
\( a_i=1/n_1 \) in group 1 and \( -1/n_2 \) in group 2. With \( \sigma_i^2=\sigma_1^2 \) in the small group and \( 0 \) in the
large one, the true variance is \( \sigma_1^2/n_1 \) and the expected bootstrap variance is
\( (1/n_1+1/n_2)(n_1-1)\sigma_1^2/n \). With \( n_1 \) fixed the ratio of the second to the first tends to
\( 0 \) as \( n_2\to\infty \).
:::

::: {#exr-bs-case-singular}
[B2]

In a one-way layout with \( g \) groups of \( m \) observations each, find the probability that a case-bootstrap resample
omits some group entirely, so that the cell-means model cannot be fitted. Evaluate it for \( g=5 \), \( m=4 \), and
describe a resampling scheme that avoids the problem.
:::
