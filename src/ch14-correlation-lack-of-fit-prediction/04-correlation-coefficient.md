# The correlation coefficient and Fisher's z

With a single regressor, the multiple correlation reduces to the ordinary correlation \( \rho \) of
the pair \( (X,Y) \), and \( R^2 \) to the square of the sample correlation
\[
r=\frac{\sum_i(X_i-\bar X)(Y_i-\bar Y)}{\sqrt{\sum_i(X_i-\bar X)^2\sum_i(Y_i-\bar Y)^2}} .
\]
Keeping the sign, this section develops two routes to inference about \( \rho \): the exact
distribution of \( r \), and Fisher's transformation to a nearly
normal variable whose variance does not depend on \( \rho \). Throughout,
\( (X_i,Y_i) \), \( i=1,\dots,n \), are independent draws from a bivariate normal distribution with correlation
\( \lvert\rho\rvert<1 \), and \( n\ge3 \).

## The exact distribution of r

::: {#prp-cor-r-exact}
[Exact law of the sample correlation]

Let \( T=r\sqrt{n-2}/\sqrt{1-r^2} \) and \( \delta=\rho/\sqrt{1-\rho^2} \).

::: {.enumerate options="label=(\alph*)"}
1. Given the \( X \) values, \( T\sim t\bigl(n-2,\delta\sqrt V\bigr) \), where \( V=\sum_i(X_i-\bar X)^2/\sigma_X^2 \).

2. \( V\sim\chi^2(n-1) \), and for \( -1<c<1 \),
   \[
\Pr(r\le c)=\E\Bigl[\Pr\bigl\{t(n-2,\delta\sqrt V)\le c\sqrt{n-2}/\sqrt{1-c^2}\bigr\}\Bigr].
\]{#eq-cor-r-cdf}

3. If \( \rho=0 \), then \( T\sim t(n-2) \) and \( r^2\sim\mathrm{Beta}\bigl(\tfrac12,\tfrac{n-2}2\bigr) \), independently of the \( X \) values.

4. The law of \( r \) is stochastically increasing in \( \rho \): for each \( c\in(-1,1) \), \( \Pr(r\le c) \) is a strictly
   decreasing function of \( \rho \).
:::

:::

::: {.proof}
(a) Let \( \hat{\beta}_1 \) be the least squares slope of \( Y \) on \( X \) and \( t_1 \) its \( t \) statistic for the value zero.
By @thm-ss-partial-r2(c) (see also @exr-proj-partial-correlation), with only the intercept as the other column, \( r^2=t_1^2/(t_1^2+n-2) \),
and \( t_1 \) has the sign of \( \hat{\beta}_1 \), which is the sign of \( r \). Solving, \( t_1=T \). Given the \( X \) values
the pairs follow a normal simple regression with slope \( \beta_1=\rho\sigma_Y/\sigma_X \) and error variance
\( \sigma^2=\sigma_Y^2(1-\rho^2) \) (@eq-mvn-random-x-model), so by @cor-opt-t with \( d=0 \), \( T \) is noncentral \( t \) with
\( n-2 \) degrees of freedom and noncentrality
\( \beta_1\sqrt{S_{xx}}/\sigma=\rho\sqrt{S_{xx}}/(\sigma_X\sqrt{1-\rho^2})=\delta\sqrt V \).

(b) \( V\sim\chi^2(n-1) \) by @cor-qf-sample-variance. Since \( T \) is an increasing function of \( r \), \( r\le c \)
iff \( T\le c\sqrt{n-2}/\sqrt{1-c^2} \), and @eq-cor-r-cdf is the expectation of the conditional
probability.

(c) With \( \delta=0 \) the conditional law is \( t(n-2) \) for every design; apply @lem-cor-conditioning. Then
\( r^2=T^2/(T^2+n-2) \) is the \( k=1 \) case of @thm-cor-r2-distribution(b).

(d) Write a \( t(\nu,\theta) \) variable as \( (Z+\theta)/\sqrt{W/\nu} \) (@def-qf-noncentral-t). For \( \theta_1<\theta_2 \) and any \( u \),
\( \Pr\{t(\nu,\theta_1)\le u\}-\Pr\{t(\nu,\theta_2)\le u\}=\Pr\{u\sqrt{W/\nu}-\theta_2<Z\le u\sqrt{W/\nu}-\theta_1\}>0 \).
Since \( \delta=\rho/\sqrt{1-\rho^2} \) is strictly increasing in \( \rho \), for every \( c \) and every \( V>0 \) the
conditional probability in @eq-cor-r-cdf strictly decreases in \( \rho \), and so does its expectation.
:::

Fisher (1915) found the density of \( r \) as an infinite series; for computation the mixture
@eq-cor-r-cdf, an average of noncentral \( t \) probabilities over quantiles of \( V \), is enough. Part (c) is
the familiar test of \( \rho=0 \), the \( t \) test of a zero slope, symmetric in the two variables. Part (d)
makes exact intervals possible, by the inversion used for \( \rho^2 \) in [Section 14.3](03-multiple-correlation.html).

```{.python .run #cell-fisher-z-exact}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

def r_cdf(c, rho, n, m=4000):
    """Exact P(r <= c) for n iid bivariate normal pairs with correlation rho."""
    if c < 0:                                                # law of -r is that for -rho
        return 1 - r_cdf(-c, -rho, n, m)
    t = c * np.sqrt(n - 2) / np.sqrt(1 - c ** 2)            # increasing in c
    if rho == 0:
        return stats.t.cdf(t, n - 2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)     # equal-probability nodes
    delta = rho / np.sqrt(1 - rho ** 2)
    return np.mean(stats.nct.cdf(t, n - 2, delta * np.sqrt(V)))
```

When \( \rho\ne0 \) the distribution of \( r \) is skewed towards zero, because \( r \) cannot pass \( \pm1 \). For
\( \rho=0.8 \) and \( n=15 \), \( 200000 \) simulated samples give a mean of \( 0.7886 \), a median of
\( 0.811 \), a skewness of \( -1.40 \) and a variance of \( 0.01195 \) (panel (a) of
[Figure 14.4.1](#fig-cor-fisher-z)). The variance depends strongly on \( \rho \): the standard deviation of
\( r \) is \( 0.267 \) at \( \rho=0 \) and \( 0.032 \) at \( \rho=0.95 \) (panel (c)). An interval \( r\pm\text{constant}\times\text{standard error} \) would inherit both problems.

## Fisher's transformation

Fisher (1921) proposed to work instead with
\[
z=\tanh^{-1}r=\frac12\log\frac{1+r}{1-r},\qquad \zeta=\tanh^{-1}\rho .
\]{#eq-cor-fisher-z}

To see why this helps we need the large-sample behaviour of \( r \), and for that three standard tools
from probability theory, used here without proof (see van der Vaart 1998, chapters 2 and 3).

::: {.remark}
[Three limit theorems]

(i) *Central limit theorem.* If \( \mathbf{W}_1,\mathbf{W}_2,\dots \) are independent copies of a random vector with
mean \( \bmu \) and covariance matrix \( \boldsymbol{\Gamma} \), then \( \sqrt n(\bar{\mathbf{W}}_n-\bmu)\to\Normal(\bzero,\boldsymbol{\Gamma}) \) in distribution.
(ii) *Slutsky's lemma.* If \( \mathbf{U}_n\to\mathbf{U} \) in distribution and \( \mathbf{V}_n\to\mathbf{c} \) in probability, then
\( \mathbf{U}_n+\mathbf{V}_n\to\mathbf{U}+\mathbf{c} \), and \( V_nU_n\to cU \) for scalars. (iii) *Delta method.* If
\( \sqrt n(\mathbf{T}_n-\boldsymbol{\theta})\to\Normal(\bzero,\boldsymbol{\Gamma}) \) and \( g \) is differentiable at \( \boldsymbol{\theta} \) with gradient
\( \nabla g \), then \( \sqrt n\bigl(g(\mathbf{T}_n)-g(\boldsymbol{\theta})\bigr)\to\Normal\bigl(0,\nabla g\T\boldsymbol{\Gamma}\nabla g\bigr) \).
:::

The limiting variance of \( r \) needs fourth moments of the bivariate normal distribution.

::: {#lem-cor-bivariate-moments}
[Fourth moments of a standard bivariate normal pair]

Let \( (X,Y) \) be bivariate normal with means \( 0 \), variances \( 1 \) and correlation \( \rho \). Then
\( \E X^4=3 \), \( \E X^3Y=3\rho \) and \( \E X^2Y^2=1+2\rho^2 \). Consequently the covariance matrix of
\( (X^2,Y^2,XY) \) is
\[
\boldsymbol{\Gamma}=\begin{pmatrix}2&2\rho^2&2\rho\\2\rho^2&2&2\rho\\2\rho&2\rho&1+\rho^2\end{pmatrix}.
\]
:::

::: {.proof}
Write \( Y=\rho X+\sqrt{1-\rho^2}\,U \) with \( U\sim\Normal(0,1) \) independent of \( X \); the pair has the right means,
variances and covariance, and it is normal by @thm-mvn-linear. Using \( \E X^2=1 \), \( \E X^4=3 \) (the fourth moment of a
standard normal) and the vanishing of odd moments,
\( \E X^3Y=\rho\E X^4=3\rho \) and \( \E X^2Y^2=\rho^2\E X^4+(1-\rho^2)\E X^2\E U^2=3\rho^2+1-\rho^2=1+2\rho^2 \). The
covariances follow by subtracting products of means, \( \E X^2=\E Y^2=1 \) and \( \E XY=\rho \): for instance
\( \Var(XY)=1+2\rho^2-\rho^2 \) and \( \Cov(X^2,XY)=\E X^3Y-\rho=2\rho \).
:::

::: {#thm-cor-fisher-z}
[Fisher's z transformation]

As \( n\to\infty \),
\[
\sqrt n\,(r-\rho)\to\Normal\bigl(0,(1-\rho^2)^2\bigr)
\qquad\text{and}\qquad
\sqrt n\,(z-\zeta)\to\Normal(0,1)
\]
in distribution. The same limits hold with \( \sqrt{n-3} \) in place of \( \sqrt n \). The limiting variance of \( z \)
does not depend on \( \rho \): \( \tanh^{-1} \) is a variance-stabilizing transformation for \( r \).
:::

::: {.proof}
The correlation is unchanged by the maps \( X\mapsto a+bX \), \( Y\mapsto c+dY \) with \( b,d>0 \), so we may assume
\( \E X=\E Y=0 \) and \( \Var X=\Var Y=1 \). Let \( \bar{\mathbf{W}}_n \) be the average of \( \mathbf{W}_i=(X_i^2,Y_i^2,X_iY_i) \), whose
mean is \( \boldsymbol{\theta}=(1,1,\rho) \) and whose covariance is \( \boldsymbol{\Gamma} \) of @lem-cor-bivariate-moments. With divisor \( n \),
the sample variances and covariance are
\( \bar{\mathbf{W}}_n-(\bar X^2,\bar Y^2,\bar X\bar Y) \). Since \( \sqrt n\,\bar X \) converges in distribution and \( \bar X\to0 \) in
probability, \( \sqrt n\,\bar X^2\to0 \) in probability, and similarly for the other two products, so by
Slutsky's lemma the centring does not affect the limit. The divisor cancels in \( r=g(\text{sample moments}) \),
where \( g(a,b,c)=c/\sqrt{ab} \), and \( \nabla g(\boldsymbol{\theta})=(-\rho/2,-\rho/2,1) \). By the central limit theorem and the
delta method, \( \sqrt n(r-\rho) \) is asymptotically normal with variance
\[
\nabla g\T\boldsymbol{\Gamma}\nabla g
=\frac{\rho^2}{4}(2+2+4\rho^2)-2\cdot\frac{\rho}{2}(2\rho+2\rho)+(1+\rho^2)
=\rho^2+\rho^4-4\rho^2+1+\rho^2=(1-\rho^2)^2 .
\]
The derivative of \( \tanh^{-1} \) at \( \rho \) is \( 1/(1-\rho^2) \), so a second application of the delta method gives
variance \( 1 \) for \( z \). Finally \( \sqrt{n-3}/\sqrt n\to1 \), and Slutsky's lemma allows the substitution.
:::

The transformation solves \( g'(\rho)=1/(1-\rho^2) \), the condition for a limiting standard deviation of one (@exr-cor-variance-stabilizing). Why the approximation is so good for small \( n \) is another matter: Fisher (1921) showed that the mean of \( z \) exceeds \( \zeta \) by
about \( \rho/(2(n-1)) \) and that its variance is close to \( 1/(n-3) \), and the distribution of \( z \) is nearly
symmetric even when that of \( r \) is strongly skewed. In the simulation with \( \rho=0.8 \) and \( n=15 \), \( z \)
has mean \( 1.1295 \) against \( \zeta=1.0986 \) and \( \zeta+\rho/(2(n-1))=1.1272 \), variance
\( 0.0811 \) against \( 1/(n-3)=0.0833 \), and skewness \( 0.01 \). Its standard deviation is
\( 0.288 \) at \( \rho=0 \) and \( 0.283 \) at \( \rho=0.95 \) (panel (c)).

::: {when-format="html"}
![**Figure 14.4.1.** The sample correlation of \( n=15 \) normal pairs. (a) \( r \) for \( \rho=0.8 \), with the exact
density, obtained by differentiating @eq-cor-r-cdf. (b) Fisher's \( z \), with the normal approximation. (c) The standard deviations of
\( r \) and \( z \) across \( \rho \); dashed: \( 1/\sqrt{n-3} \).](fisher_z.svg){#fig-cor-fisher-z width=100%}
:::

::: {when-format="pdf"}
![The sample correlation of \( n=15 \) normal pairs. (a) \( r \) for \( \rho=0.8 \), with the exact
density, obtained by differentiating @eq-cor-r-cdf. (b) Fisher's \( z \), with the normal approximation. (c) The standard deviations of
\( r \) and \( z \) across \( \rho \); dashed: \( 1/\sqrt{n-3} \).](fisher_z.pdf){width=100%}
:::

## Intervals and tests for a correlation

The approximation \( z\approx\Normal(\zeta,1/(n-3)) \) gives the interval
\[
\tanh\Bigl(z-\frac{z_{\alpha/2}}{\sqrt{n-3}}\Bigr)\ \le\ \rho\ \le\ \tanh\Bigl(z+\frac{z_{\alpha/2}}{\sqrt{n-3}}\Bigr),
\]{#eq-cor-z-interval}

with \( z_{\alpha/2} \) the upper \( \alpha/2 \) point of \( \Normal(0,1) \). It is asymmetric about \( r \), as it should be, and
it never leaves \( (-1,1) \). The exact interval inverts @eq-cor-r-cdf, using part (d) of
@prp-cor-r-exact exactly as for \( \rho^2 \): its limits are the values of \( \rho \) at which the observed \( r \)
is the \( 97.5\% \) and the \( 2.5\% \) point.

How good is the approximation? For \( n=10 \) pairs, the coverage of the \( 95\% \) interval
@eq-cor-z-interval over \( 200000 \) samples is \( 0.949 \), \( 0.949 \) and \( 0.952 \) at \( \rho=0 \),
\( 0.5 \) and \( 0.9 \). The interval \( r\pm1.96(1-r^2)/\sqrt n \), which uses the asymptotic standard error of
@thm-cor-fisher-z directly, covers only \( 0.839 \), \( 0.838 \) and \( 0.834 \).

```{.python .run #cell-fisher-z-coverage}
def sample_r(rho, n, reps, rng):
    x = rng.normal(size=(reps, n))
    y = rho * x + np.sqrt(1 - rho ** 2) * rng.normal(size=(reps, n))
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    return np.sum(xc * yc, axis=1) / np.sqrt(np.sum(xc ** 2, axis=1) * np.sum(yc ** 2, axis=1))

rng = np.random.default_rng(1443)
n_c, reps_c, zq = 10, 200_000, stats.norm.ppf(0.975)
for rho_c in (0.0, 0.5, 0.9):
    r_c = sample_r(rho_c, n_c, reps_c, rng)
    z_c = np.arctanh(r_c)
    fisher = np.abs(z_c - np.arctanh(rho_c)) <= zq / np.sqrt(n_c - 3)
    naive = np.abs(r_c - rho_c) <= zq * (1 - r_c ** 2) / np.sqrt(n_c)
    print(f"rho = {rho_c}: Fisher z covers {fisher.mean():.4f}, naive covers {naive.mean():.4f}")
```

Because \( z \) has an approximately known variance, it also compares correlations. With
independent samples of \( n_1 \) and \( n_2 \) cases, the hypothesis \( \rho_1=\rho_2 \) is tested by referring
\[
\frac{z_1-z_2}{\sqrt{1/(n_1-3)+1/(n_2-3)}}
\]
to the standard normal distribution. The \( t \) statistic of @prp-cor-r-exact(c) remains the exact test of
\( \rho=0 \); for any other null value \( \rho_0 \), either \( (z-\tanh^{-1}\rho_0)\sqrt{n-3} \) or the exact law @eq-cor-r-cdf is used.

::: {#exm-cor-states-correlation}
[Poverty and graduation]

Across the \( 50 \) states, the correlation between the poverty rate and the high-school
graduation rate is \( r=-0.7613 \). Its \( t \) statistic for \( \rho=0 \) is \( -8.13 \) on
\( 48 \) degrees of freedom, far out in the tail. The Fisher interval @eq-cor-z-interval is
\( (-0.858, -0.613) \) and the exact interval is \( (-0.855, -0.609) \), almost the same, and both
are noticeably asymmetric about \( r \).
:::

```{.python .run #cell-fisher-z-states}
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
x, y = data["poverty"].to_numpy(), data["hs_grad"].to_numpy()
n_s = len(x)
r_s = np.corrcoef(x, y)[0, 1]
t_s = r_s * np.sqrt(n_s - 2) / np.sqrt(1 - r_s ** 2)
lo_z, hi_z = np.tanh(np.arctanh(r_s) + np.array([-1, 1]) * zq / np.sqrt(n_s - 3))
# exact interval: the rho at which the observed r sits at the 97.5% and 2.5% points
lo_x = optimize.brentq(lambda rh: r_cdf(r_s, rh, n_s) - 0.975, -0.99, 0.99)
hi_x = optimize.brentq(lambda rh: r_cdf(r_s, rh, n_s) - 0.025, -0.99, 0.99)
print(f"r = {r_s:.4f}, t = {t_s:.2f} on {n_s - 2} df")
print(f"Fisher z interval ({lo_z:.3f}, {hi_z:.3f}); exact interval ({lo_x:.3f}, {hi_x:.3f})")
```

::: {.warning}
Everything in this section uses bivariate normality, and not only through the conditional linear
model. The asymptotic variance \( (1-\rho^2)^2/n \) comes from the fourth moments of @lem-cor-bivariate-moments. For other distributions it can be quite different, and then neither \( z \)
nor the exact interval is reliable (@exr-cor-nonnormal-r). The \( t \) test of \( \rho=0 \) is more robust:
by @thm-cor-conditional it needs only that \( Y \) given \( X \) follow a normal linear model with constant
variance, and when \( X \) and \( Y \) are independent this holds as soon as \( Y \) is normal.
:::

## Exercises

### A. Check your understanding

::: {#exr-cor-z-interval-numbers}
[A1]

A sample of \( n=28 \) pairs has \( r=0.62 \). Compute the \( 95\% \) Fisher interval for \( \rho \).
:::

::: {.solution}
\( z=\tanh^{-1}0.62=0.7250 \) and \( 1.96/\sqrt{25}=0.392 \), so \( \zeta \) lies in \( (0.333,1.117) \), and \( \rho \) in
\( (\tanh0.333,\tanh1.117)=(0.321,0.806) \). The interval is longer below \( r \) than above.
:::

::: {#exr-cor-two-correlations}
[A2]

Two independent samples give \( r_1=0.45 \) with \( n_1=60 \) and \( r_2=0.20 \) with \( n_2=80 \). Test \( \rho_1=\rho_2 \) at level
\( 0.05 \).
:::

### B. Practice

::: {#exr-cor-variance-stabilizing}
[B1]

Suppose \( \sqrt n(T_n-\theta)\to\Normal(0,s(\theta)^2) \) for every \( \theta \) in an interval, with \( s \) positive and
continuous. Show that if \( g'(\theta)=1/s(\theta) \), then \( \sqrt n(g(T_n)-g(\theta))\to\Normal(0,1) \). Solve
\( g'(\rho)=1/(1-\rho^2) \) and recover @eq-cor-fisher-z.
:::

::: {.solution}
By the delta method the limiting variance is \( g'(\theta)^2s(\theta)^2=1 \). With partial fractions,
\( 1/(1-\rho^2)=\tfrac12[1/(1+\rho)+1/(1-\rho)] \), whose antiderivative vanishing at zero is
\( \tfrac12[\log(1+\rho)-\log(1-\rho)]=\tanh^{-1}\rho \).
:::

::: {#exr-cor-nonnormal-r}
[B2]

Let \( X \) and \( U \) be independent standard normal variables and \( Y=XU \). Show that \( X \) and \( Y \) are
uncorrelated but not independent, and that \( \sqrt n\,r\to\Normal(0,3) \). What happens to the level of the
test that rejects \( \rho=0 \) when \( \lvert T\rvert>t_{n-2,0.025} \) in large samples? Which assumption of
@thm-cor-conditional fails?
:::

::: {.solution}
\( \E XY=\E X^2\E U=0 \), but \( \E(X^2Y^2)=\E X^4\E U^2=3\ne\E X^2\E Y^2=1 \), so they are dependent. Repeat the
proof of @thm-cor-fisher-z at \( \rho=0 \): the gradient is \( (0,0,1) \), so the limiting variance of \( \sqrt n\,r \) is
\( \Var(XY)/(\Var X\Var Y)=3 \). Since \( T\approx\sqrt n\,r \) for large \( n \), the test rejects with probability tending
to \( 2\{1-\Phi(1.96/\sqrt3)\}=0.258 \), not \( 0.05 \). Given \( X \), \( Y \) is normal with mean zero but variance \( X^2 \): the
errors are not homoscedastic, so the conditional model of @def-cor-random-regressors fails.
:::

### C. Going deeper

::: {#exr-cor-r-general-variance}
[C1]

For a general pair \( (X,Y) \) with finite fourth moments, standardized to means \( 0 \) and variances \( 1 \), show that
\( \sqrt n(r-\rho) \) is asymptotically normal with variance
\( \bigl(1+\tfrac{\rho^2}{2}\bigr)\E X^2Y^2+\tfrac{\rho^2}{4}(\E X^4+\E Y^4)-\rho(\E X^3Y+\E XY^3) \). Check that it reduces to
\( (1-\rho^2)^2 \) under normality, and explain why no transformation of \( r \) stabilizes the variance across
all distributions.
:::
