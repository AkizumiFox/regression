# The wild bootstrap

The residual bootstrap keeps the design but moves every residual to a random position, destroying any link between the size
of an error and where it occurred. The case bootstrap keeps that link but changes the design. The **wild bootstrap** keeps both:
each residual stays with its own row of \( \X \), and only a random multiplier changes (Wu 1986; Liu 1988; Mammen 1993).

::: {#def-bs-wild}
[Wild bootstrap]

Let \( \X \) have rank \( p \), with least squares estimate \( \hbeta \), residuals \( \hat{\varepsilon}_i \) and leverages
\( h_{ii}<1 \). Choose transformed residuals \( u_i \), one of
\[
u_i=\hat{\varepsilon}_i,\qquad u_i=\frac{\hat{\varepsilon}_i}{\sqrt{1-h_{ii}}},\qquad u_i=\frac{\hat{\varepsilon}_i}{1-h_{ii}} ,
\]
and a law for the **multipliers**, with mean \( 0 \) and variance \( 1 \). The wild bootstrap draws \( v_1,\dots,v_n \)
independently from that law, independently of the data, sets
\[
Y^*_i=\x_{(i)}\T\hbeta+u_i\,v_i ,\qquad i=1,\dots,n,
\]
and computes \( \hbeta^*=(\X\T\X)^{-1}\X\T\Y^* \). Two standard multiplier laws are the **Rademacher** law,
\( v=\pm1 \) with probability \( 1/2 \) each, and **Mammen's** two-point law
\[
\Pr\Bigl(v=\frac{1-\sqrt5}{2}\Bigr)=\frac{\sqrt5+1}{2\sqrt5},\qquad
\Pr\Bigl(v=\frac{1+\sqrt5}{2}\Bigr)=\frac{\sqrt5-1}{2\sqrt5} .
\]
:::

In the bootstrap world the error of case \( i \) is \( u_iv_i \), with mean zero and its own variance \( u_i^2 \). One residual is a poor
estimate of a variance, but \( \mathbf{c}\T\hbeta \) needs only a weighted average of the variances, as for the sandwich.

::: {#prp-bs-wild-moments}
[Moments of the wild bootstrap]

In the wild bootstrap, with \( \mathbf{U}=\diag(u_1,\dots,u_n) \) and \( \mathbf{a}=\X(\X\T\X)^{-1}\mathbf{c} \):

::: {.enumerate options="label=(\alph*)"}
1. \( \E_*(\hbeta^*)=\hbeta \) and
   \[
\Cov_*(\hbeta^*)=(\X\T\X)^{-1}\X\T\mathbf{U}^2\X(\X\T\X)^{-1},
\]{#eq-bs-wild-cov}

   which is the sandwich estimator HC0, HC2 or HC3 (@thm-het-sandwich) for the three choices of \( u_i \);

2. \( \E_*\bigl[(\mathbf{c}\T\hbeta^*-\mathbf{c}\T\hbeta)^3\bigr]=\E(v^3)\sum_ia_i^3u_i^3 \).
:::

:::

::: {.proof}
\( \hbeta^*-\hbeta=(\X\T\X)^{-1}\X\T\mathbf{U}\mathbf{v} \), where \( \mathbf{v}=(v_1,\dots,v_n)\T \) has mean \( \bzero \)
and covariance \( \I \). So \( \E_*(\hbeta^*)=\hbeta \) and \( \Cov_*(\hbeta^*)=(\X\T\X)^{-1}\X\T\mathbf{U}\mathbf{U}\X(\X\T\X)^{-1} \).
For (b), \( \mathbf{c}\T(\hbeta^*-\hbeta)=\sum_ia_iu_iv_i \) is a sum of independent terms with mean zero, and the third moment of
such a sum is the sum of the third moments, because every cross term contains some \( v_i \) to the first power.
:::

So the wild bootstrap reproduces the sandwich exactly, not only to first order as the case bootstrap does (@prp-bs-case-linear). Part (b) explains the choice of multipliers. The true third moment is
\( \sum_ia_i^3\E(\varepsilon_i^3) \), not zero for skewed errors. Rademacher multipliers have \( \E(v^3)=0 \) and give a symmetric
bootstrap distribution; Mammen's law is the two-point law with third moment \( 1 \) (@exr-bs-mammen) and carries the skewness of the
residuals. Davidson and Flachaire (2008) nonetheless favour Rademacher weights for tests, because matching the third moment costs
accuracy in the fourth.

## Consistency under heteroscedasticity

Take the fixed designs of [Section 23.1](01-cases-and-residuals.html), with contrasts \( \mathbf{c}_n \), weights \( \mathbf{a}_n \),
\( w_i=a_{ni}/\norm{\mathbf{a}_n} \) and @eq-bs-no-dominant. The errors
\( \varepsilon_{n1},\dots,\varepsilon_{nn} \) are independent with mean zero and variances \( \sigma_{ni}^2 \), and for some
constants \( 0<\underline{\sigma}^2 \) and \( K<\infty \),
\[
\sigma_{ni}^2\ge\underline{\sigma}^2,\qquad \E\varepsilon_{ni}^4\le K\qquad\text{for all } n,i .
\]{#eq-bs-hetero-conditions}

The variance of the standardized contrast is \( \tau_n^2=\sum_iw_i^2\sigma_{ni}^2 \), which lies between
\( \underline{\sigma}^2 \) and \( \sqrt K \), since \( \sigma_{ni}^2\le(\E\varepsilon_{ni}^4)^{1/2} \).

::: {#thm-bs-wild-consistency}
[Consistency of the wild bootstrap]

Under @eq-bs-no-dominant and @eq-bs-hetero-conditions, let \( \hbeta^* \) come from the wild bootstrap with
\( u_i=\hat{\varepsilon}_i \) and multipliers from any fixed law with mean \( 0 \) and variance \( 1 \), and let \( G_n \) and
\( \hat{G}_n \) be defined as in @thm-bs-residual, with this \( \hbeta^* \). Then \( \sup_x\lvert\hat{G}_n(x)-G_n(x)\rvert\to0 \) in
probability, and both are uniformly close to \( \Phi(x/\tau_n) \).
:::

::: {.proof}
Write \( \mathbf{d}=\M\be \), so that \( \hat{\varepsilon}_i=\varepsilon_i-d_i \), and put \( m_n=\max_i\lvert w_i\rvert \). Four steps.

*Step 1: the bootstrap variance.* Let \( S_n^2=\sum_iw_i^2\hat{\varepsilon}_i^2 \), the bootstrap variance of
\( \mathbf{c}_n\T(\hbeta^*-\hbeta)/\norm{\mathbf{a}_n} \) by @prp-bs-wild-moments. We show \( S_n^2-\tau_n^2\to0 \) in probability.
Expand \( S_n^2=\sum_iw_i^2\varepsilon_i^2-2\sum_iw_i^2\varepsilon_id_i+\sum_iw_i^2d_i^2 \). The first sum has mean \( \tau_n^2 \)
and variance \( \sum_iw_i^4\Var(\varepsilon_i^2)\le Km_n^2\sum_iw_i^2=Km_n^2\to0 \). The third is at most
\( m_n^2\norm{\mathbf{d}}^2 \), and \( \E\norm{\mathbf{d}}^2=\tr\bigl(\M\Cov(\be)\bigr)\le p\sqrt K \), so it tends to zero in probability. The
middle sum is at most \( 2\bigl(\sum_iw_i^2\varepsilon_i^2\bigr)^{1/2}\bigl(\sum_iw_i^2d_i^2\bigr)^{1/2} \) by Cauchy–Schwarz, a bounded
factor times one that tends to zero.

*Step 2: no dominant term.* \( \max_i\lvert w_i\hat{\varepsilon}_i\rvert\le\max_i\lvert w_i\varepsilon_i\rvert+m_n\norm{\mathbf{d}} \). The
second term tends to zero in probability, and for \( \eta>0 \)
\[
\Pr\bigl(\max_i\lvert w_i\varepsilon_i\rvert>\eta\bigr)\le\sum_i\frac{w_i^4\E\varepsilon_i^4}{\eta^4}\le\frac{Km_n^2}{\eta^4} ,
\]
which tends to zero. So \( \max_i\lvert w_i\hat{\varepsilon}_i\rvert\to0 \) in probability.

*Step 3: the bootstrap limit.* Given the data, \( U^*_{ni}=w_i\hat{\varepsilon}_iv_i/S_n \) are independent with mean zero and
variances summing to one. With \( g(t)=\E\bigl[v^2\mathbf{1}\{\lvert v\rvert>t\}\bigr] \), which decreases to zero as
\( t\to\infty \),
\[
\begin{aligned}
L^*_n(\eta)&=\sum_i\frac{w_i^2\hat{\varepsilon}_i^2}{S_n^2}\,g\Bigl(\frac{\eta S_n}{\lvert w_i\hat{\varepsilon}_i\rvert}\Bigr)\\
&\le g\Bigl(\frac{\eta S_n}{\max_i\lvert w_i\hat{\varepsilon}_i\rvert}\Bigr).
\end{aligned}
\]
By Steps 1 and 2 the argument of \( g \) tends to infinity in probability, because \( S_n^2\ge\underline{\sigma}^2/2 \)
with probability tending to one. So \( L^*_n(\eta)\to0 \) in probability, and @lem-bs-conditional-clt gives
\( \sup_x\lvert{\Pr}_*(\sum_iw_i\hat{\varepsilon}_iv_i\le x)-\Phi(x/S_n)\rvert\to0 \) in probability. Since \( S_n/\tau_n\to1 \) in
probability, the same holds with \( \Phi(x/\tau_n) \), as at the end of the proof of @thm-bs-residual.

*Step 4: the true distribution.* \( \mathbf{c}_n\T(\hbeta-\bbeta)/\norm{\mathbf{a}_n}=\sum_iw_i\varepsilon_i \). With
\( U_{ni}=w_i\varepsilon_i/\tau_n \), Lyapunov's condition holds:
\( \sum_i\E U_{ni}^4\le Km_n^2/\underline{\sigma}^4\to0 \). So \( \sum_iw_i\varepsilon_i/\tau_n\to\Normal(0,1) \), and by Pólya's theorem
\( \sup_x\lvert\Pr(\sum_iw_i\varepsilon_i\le x)-\Phi(x/\tau_n)\rvert\to0 \). The triangle inequality completes the proof.
:::

An argument like Step 1 shows why the residual bootstrap fails here. Its variance for the standardized contrast is
\( \hat{\sigma}_*^2 \), whose expectation \( n^{-1}\sum_i(1-h_{ii})\sigma_{ni}^2 \) is an unweighted average of the error variances (@exr-bs-heteroscedastic-residual), not the weighted average \( \tau_n^2 \). So it settles on a normal law with the
variance the classical formula gives. The wild bootstrap succeeds because each resampled error keeps the
variance of its own case.

::: {#exm-bs-engel-wild}
[The wild bootstrap for Engel's data]

For the line of @exm-bs-engel, \( 4000 \) wild resamples with \( u_i=\hat{\varepsilon}_i \) give standard errors
\( 0.05119 \) (Rademacher) and \( 0.05217 \) (Mammen), Monte Carlo estimates of the HC0 value \( 0.05177 \).
With \( u_i=\hat{\varepsilon}_i/(1-h_{ii}) \) they give \( 0.06653 \), estimating HC3, \( 0.06639 \).

The shape is more interesting. The wild distribution in
[Figure 23.1.1](01-cases-and-residuals.html#fig-bs-engel)(b) has two modes, at \( 0.439 \) and
\( 0.531 \). The cause is the richest household, whose leverage \( 0.255 \) was examined in
[Chapter 20](../ch20-residuals-leverage-influence/index.html) and whose share of the HC3 variance was measured in
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html). Its term \( a_i^2\hat{\varepsilon}_i^2 \) accounts
for a fraction \( 0.78 \) of the HC0 variance \( \sum_ja_j^2\hat{\varepsilon}_j^2 \), so the condition behind Step 2 of the
proof is badly violated. The bootstrap slope is dominated by the contribution \( a_i\hat{\varepsilon}_iv_i \) of this one
household, which with Rademacher multipliers is \( \pm0.046 \), and the modes sit at \( \hat{\beta}_1\mp0.046 \). The
sampling distribution of the slope is dominated by that household's error too, and one observation cannot reveal the law it
came from; the remedy is the one given in those chapters.
:::

```{.python .run #cell-engel-bootstrap-wild}
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
rng = np.random.default_rng(2301)
B = 4000
A = XtX_inv @ X.T                                 # beta_hat = A y
fit = X @ beta_hat

golden = (1 + np.sqrt(5)) / 2
mammen_values = np.array([1 - golden, golden])    # -0.618 and 1.618
mammen_probs = np.array([golden / np.sqrt(5), 1 - golden / np.sqrt(5)])

def wild_bootstrap(resid, B, weights="rademacher"):
    """Slopes from y*_i = x_i' beta_hat + resid_i v_i, with independent multipliers v_i."""
    if weights == "rademacher":
        v = rng.choice([-1.0, 1.0], size=(B, n))
    else:
        v = rng.choice(mammen_values, p=mammen_probs, size=(B, n))
    y_star = fit + resid * v
    return (y_star @ A.T)[:, 1]

slopes_wild = wild_bootstrap(e_hat, B)
slopes_mammen = wild_bootstrap(e_hat, B, weights="mammen")
slopes_wild3 = wild_bootstrap(e_hat / (1 - h), B)
se_hc0 = np.sqrt(A[1] ** 2 @ e_hat ** 2)
se_hc3 = np.sqrt(A[1] ** 2 @ (e_hat / (1 - h)) ** 2)
print(f"wild se: Rademacher {slopes_wild.std():.5f}, Mammen {slopes_mammen.std():.5f}, "
      f"HC3-type {slopes_wild3.std():.5f}")
print(f"sandwich se: HC0 {se_hc0:.5f}, HC3 {se_hc3:.5f}")
```

## Bootstrap or sandwich?

For a standard error the wild bootstrap adds nothing to the sandwich it reproduces. It earns its place through the
*studentized* statistic \( (\mathbf{c}\T\hbeta-\mathbf{c}\T\bbeta)/\widehat{\text{se}}_{\text{HC}} \), whose noisy denominator
makes its small-sample distribution far from \( t \) when some leverages are large; the simulations of
[Section 23.3](03-confidence-intervals.html) show the wild bootstrap-\( t \) beating the HC2 \( t \) interval. For a test of \( \mathbf{c}\T\bbeta=d \), the resamples should come from the restricted fit, so that
the bootstrap approximates the null distribution even when the null hypothesis is false (Davidson and Flachaire 2008).

## Exercises

### A. Check your understanding

::: {#exr-bs-rademacher-exact}
[A1]

With Rademacher multipliers and \( u_i=\hat{\varepsilon}_i \), how many distinct wild resamples are there? Show that the bootstrap
variance of \( \mathbf{c}\T\hbeta^* \) computed by averaging over all of them equals the HC0 variance exactly.
:::

::: {.solution}
There are \( 2^n \) sign vectors, equally likely. Averaging \( (\sum_ia_i\hat{\varepsilon}_iv_i)^2 \) over them, every cross term
\( a_ia_j\hat{\varepsilon}_i\hat{\varepsilon}_jv_iv_j \) with \( i\ne j \) averages to zero, because half of the sign vectors give
\( v_iv_j=1 \) and half give \( -1 \). What remains is \( \sum_ia_i^2\hat{\varepsilon}_i^2 \), the \( (\mathbf{c},\mathbf{c}) \) entry of
HC0. The bootstrap mean is \( \mathbf{c}\T\hbeta \), so this is the variance.
:::

::: {#exr-bs-wild-two-sample}
[A2]

Two groups of sizes \( n_1 \) and \( n_2 \) are compared through the difference of their means, fitted by least squares with
one indicator column per group. Show that the wild bootstrap with \( u_i=\hat{\varepsilon}_i \) gives the variance
\( \hat{\sigma}_1^2/n_1+\hat{\sigma}_2^2/n_2 \), where \( \hat{\sigma}_k^2 \) is the within-group variance with divisor \( n_k \), and that
\( u_i=\hat{\varepsilon}_i/\sqrt{1-h_{ii}} \) gives the Welch variance \( s_1^2/n_1+s_2^2/n_2 \).
:::

### B. Practice

::: {#exr-bs-mammen}
[B1]

Let \( v \) take two values \( a<0<b \). Show that \( \E(v)=0 \) and \( \E(v^2)=1 \) force \( ab=-1 \) and \( \Pr(v=a)=b/(b-a) \), and that then
\( \E(v^3)=a+b \). Deduce that Mammen's law of @def-bs-wild is the only two-point law with mean \( 0 \), variance \( 1 \) and third moment \( 1 \),
and compute its fourth moment.
:::

::: {.solution}
With \( q=\Pr(v=a) \), \( qa+(1-q)b=0 \) gives \( q=b/(b-a) \). Then
\( \E(v^2)=(ba^2-ab^2)/(b-a)=-ab \), so \( ab=-1 \). Similarly \( \E(v^3)=(ba^3-ab^3)/(b-a)=-ab(a+b)=a+b \). Third moment one means
\( a+b=1 \) and \( ab=-1 \), so \( a \) and \( b \) are the roots of \( t^2-t-1=0 \), namely \( (1\mp\sqrt5)/2 \), and
\( q=b/(b-a)=(1+\sqrt5)/(2\sqrt5) \). The fourth moment is \( (ba^4-ab^4)/(b-a)=-ab(a^2+ab+b^2)=(a+b)^2-ab=2 \), compared
with \( 1 \) for Rademacher weights and \( 3 \) for the normal.
:::

### C. Going deeper

::: {#exr-bs-wild-cluster}
[C1]

Suppose the cases fall into clusters and errors are correlated within clusters but independent between them. Show that
the wild bootstrap of @def-bs-wild gives the HC0 covariance, which ignores the correlation. Propose a version in which all
cases in a cluster share one multiplier, and show that its bootstrap covariance is the cluster-robust sandwich
\( (\X\T\X)^{-1}\sum_g\X_g\T\he_g\he_g\T\X_g(\X\T\X)^{-1} \).
:::
