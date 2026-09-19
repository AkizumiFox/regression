# Stein's paradox and James–Stein estimation

The shrinkage estimators so far beat least squares only for *some* parameter values, which is
unavoidable for linear estimators (@exr-shr-no-linear-dominance). Stein (1956) discovered that a
*nonlinear* shrinkage can do better everywhere: when three or more means are estimated together
under total squared error, the estimator that is at once least squares, maximum likelihood and best
unbiased is inadmissible. James and Stein (1961) gave an explicit estimator that dominates it. This
section proves their theorem and explains it as an empirical Bayes phenomenon.

## The normal means problem

Let \( \Z\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \) with \( \sigma^2 \) known, and measure an estimator
\( \hat{\boldsymbol{\uptheta}}(\Z) \) by its **risk** \( R(\boldsymbol{\uptheta},\hat{\boldsymbol{\uptheta}})=\E\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}}^2 \),
its total mean squared error. The estimator \( \Z \) has risk \( p\sigma^2 \) for every
\( \boldsymbol{\uptheta} \). An estimator **dominates** another if its risk is never larger and is
smaller for some \( \boldsymbol{\uptheta} \); an estimator dominated by none is **admissible**.

This is regression in canonical form. If \( \X\T\X=\I \), then \( \hbeta\sim\Normal_p(\bbeta,\sigma^2\I) \) (@thm-opt-sampling(a)). For a general full-rank design, \( \Z=(\X\T\X)^{1/2}\hbeta \) is
\( \Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \) with \( \boldsymbol{\uptheta}=(\X\T\X)^{1/2}\bbeta \), and for any
estimator \( \tilde{\bbeta} \) with \( \hat{\boldsymbol{\uptheta}}=(\X\T\X)^{1/2}\tilde{\bbeta} \),
\[
\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}}^2=(\tilde{\bbeta}-\bbeta)\T\X\T\X(\tilde{\bbeta}-\bbeta)=\norm{\X\tilde{\bbeta}-\X\bbeta}^2,
\]
the prediction loss at the design points, with the square root of @thm-mat-square-root as the
change of coordinates.

## Stein's lemma

The proof rests on an integration-by-parts identity for the normal distribution.

::: {#lem-shr-stein}
[Stein's lemma]

::: {.enumerate options="label=(\alph*)"}
1. Let \( Z\sim\Normal(\mu,\sigma^2) \), and let \( g:\Real\to\Real \) be absolutely continuous with
   \( \E\lvert g'(Z)\rvert<\infty \). Then \( \E\lvert(Z-\mu)g(Z)\rvert<\infty \) and
   \[
   \E\bigl[(Z-\mu)\,g(Z)\bigr]=\sigma^2\,\E\,g'(Z).
   \]

2. Let \( \Z\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \), and let \( \mathbf{g}=(g_1,\dots,g_p):\Real^p\to\Real^p \) be
   such that, for each \( i \) and almost every value of the other coordinates, \( g_i \) is absolutely
   continuous in \( z_i \), with \( \E\lvert\partial g_i(\Z)/\partial z_i\rvert<\infty \) and
   \( \E\lvert(Z_i-\theta_i)g_i(\Z)\rvert<\infty \). Then
   \[
   \E\bigl[(\Z-\boldsymbol{\uptheta})\T\mathbf{g}(\Z)\bigr]=\sigma^2\,\E\sum_{i=1}^p\frac{\partial g_i}{\partial z_i}(\Z).
   \]{#eq-shr-stein}

:::
:::

::: {.proof}
(a) Write \( Z=\mu+\sigma W \) with \( W\sim\Normal(0,1) \) and \( h(w)=g(\mu+\sigma w) \), so that
\( h'(w)=\sigma g'(\mu+\sigma w) \). The claim becomes \( \E[Wh(W)]=\E\,h'(W) \), so we may take \( \mu=0 \),
\( \sigma=1 \), \( g=h \). Let \( \varphi \) be the standard normal density. Since \( \varphi'(t)=-t\varphi(t) \),
\[
\begin{aligned}
\varphi(w)&=\int_w^\infty t\varphi(t)\,dt\quad(w\ge0),\\
\varphi(w)&=-\int_{-\infty}^wt\varphi(t)\,dt\quad(w\le0).
\end{aligned}
\]
Hence
\[
\begin{aligned}
\E\,h'(W)&=\int_0^\infty h'(w)\int_w^\infty t\varphi(t)\,dt\,dw\\
&\qquad-\int_{-\infty}^0h'(w)\int_{-\infty}^wt\varphi(t)\,dt\,dw .
\end{aligned}
\]
The double integrals converge absolutely: replacing \( h' \) by \( \lvert h'\rvert \) and \( t \) by \( \lvert t\rvert \) gives
\( \int\lvert h'(w)\rvert\varphi(w)\,dw=\E\lvert h'(W)\rvert<\infty \). By Fubini's theorem we may reverse the order of
integration. The first term becomes \( \int_0^\infty t\varphi(t)\int_0^th'(w)\,dw\,dt=\int_0^\infty t\varphi(t)\{h(t)-h(0)\}\,dt \),
and the second \( -\int_{-\infty}^0t\varphi(t)\int_t^0h'(w)\,dw\,dt=\int_{-\infty}^0t\varphi(t)\{h(t)-h(0)\}\,dt \), using
absolute continuity of \( h \). Their sum is \( \E[W\{h(W)-h(0)\}]=\E[Wh(W)] \), since \( \E W=0 \). The same bound
shows \( \E\lvert W\{h(W)-h(0)\}\rvert\le\E\lvert h'(W)\rvert \), so \( \E\lvert Wh(W)\rvert<\infty \).
(b) The coordinates of \( \Z \) are independent. Fix \( i \) and condition on the other coordinates
\( \Z_{-i} \). By Fubini, \( \E\bigl[\lvert\partial g_i/\partial z_i\rvert\mid\Z_{-i}\bigr]<\infty \) for almost every value
of \( \Z_{-i} \), and for those values part (a), applied to \( z_i\mapsto g_i(\Z) \) with \( Z_i\sim\Normal(\theta_i,\sigma^2) \),
gives \( \E[(Z_i-\theta_i)g_i(\Z)\mid\Z_{-i}]=\sigma^2\E[\partial g_i(\Z)/\partial z_i\mid\Z_{-i}] \). Take expectations and sum over \( i \).
:::

The lemma turns an expectation involving the unknown \( \boldsymbol{\uptheta} \) into the expectation
of a function of \( \Z \) alone. Stein (1981) built a theory of risk estimation on it (@exr-shr-sure),
and it computes covariance degrees of freedom (@exr-shr-df-covariance) for
nonlinear fits such as the lasso (@exr-shr-lasso-df).

## The James–Stein theorem

::: {#thm-shr-james-stein}
[James–Stein]

Let \( \Z\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \) with \( p\ge3 \) and \( \sigma^2 \) known, and for a constant \( c \)
let \( \hat{\boldsymbol{\uptheta}}_c=(1-c\sigma^2/\norm{\Z}^2)\Z \). Then \( \E\norm{\Z}^{-2}<\infty \) and:

::: {.enumerate options="label=(\alph*)"}
1. \( \E\norm{\hat{\boldsymbol{\uptheta}}_c-\boldsymbol{\uptheta}}^2=p\sigma^2-c\,\{2(p-2)-c\}\,\sigma^4\,\E\norm{\Z}^{-2} \). This is
   below \( p\sigma^2 \) for every \( \boldsymbol{\uptheta} \) iff \( 0<c<2(p-2) \), and smallest at \( c=p-2 \). The
   **James–Stein estimator** is
   \[
   \hat{\boldsymbol{\uptheta}}^{\mathrm{JS}}=\Bigl(1-\frac{(p-2)\sigma^2}{\norm{\Z}^2}\Bigr)\Z,
   \]
   with risk
   \[
   \E\norm{\hat{\boldsymbol{\uptheta}}^{\mathrm{JS}}-\boldsymbol{\uptheta}}^2=p\sigma^2-(p-2)^2\sigma^4\,\E\norm{\Z}^{-2}<p\sigma^2 .
   \]{#eq-shr-js}

2. With \( \gamma=\norm{\boldsymbol{\uptheta}}^2/\sigma^2 \) and \( K\sim\text{Poisson}(\gamma/2) \),
   \( \sigma^2\E\norm{\Z}^{-2}=\E\{1/(p-2+2K)\} \). The risk of \( \hat{\boldsymbol{\uptheta}}^{\mathrm{JS}} \) depends on
   \( \boldsymbol{\uptheta} \) only through \( \gamma \), equals \( 2\sigma^2 \) at \( \boldsymbol{\uptheta}=\bzero \), increases strictly in
   \( \gamma \), and tends to \( p\sigma^2 \) as \( \gamma\to\infty \).
:::
:::

::: {.proof}
*Finiteness.* The density of \( \Z \) is at most \( (2\pi\sigma^2)^{-p/2} \), so
\[
\begin{aligned}
\E\norm{\Z}^{-2}&\le(2\pi\sigma^2)^{-p/2}\int_{\norm{\bz}\le1}\norm{\bz}^{-2}\,d\bz+1\\
&=(2\pi\sigma^2)^{-p/2}A_p\int_0^1r^{p-3}\,dr+1<\infty,
\end{aligned}
\]
where \( A_p \) is the surface area of the unit sphere in \( \Real^p \) and the integral is finite because
\( p\ge3 \).
(a) Apply @eq-shr-stein to \( \mathbf{g}(\bz)=\bz/\norm{\bz}^2 \). For \( p\ge3 \) the other coordinates
\( \Z_{-i} \) are nonzero with probability one, and then \( g_i \) is a smooth function of \( z_i \) with
\[
\frac{\partial g_i}{\partial z_i}=\frac1{\norm{\bz}^2}-\frac{2z_i^2}{\norm{\bz}^4},\qquad
\Bigl\lvert\frac{\partial g_i}{\partial z_i}\Bigr\rvert\le\frac3{\norm{\bz}^2},\qquad
\sum_i\frac{\partial g_i}{\partial z_i}=\frac{p-2}{\norm{\bz}^2}.
\]
The bound makes \( \E\lvert\partial g_i/\partial z_i\rvert \) finite, and
\( \lvert(Z_i-\theta_i)g_i(\Z)\rvert\le\lvert Z_i-\theta_i\rvert/\norm{\Z} \) has finite expectation by the
Cauchy–Schwarz inequality. So
\( \E[(\Z-\boldsymbol{\uptheta})\T\Z/\norm{\Z}^2]=(p-2)\sigma^2\E\norm{\Z}^{-2} \). Now expand
\[
\norm{\hat{\boldsymbol{\uptheta}}_c-\boldsymbol{\uptheta}}^2=\norm{\Z-\boldsymbol{\uptheta}}^2-2c\sigma^2\frac{(\Z-\boldsymbol{\uptheta})\T\Z}{\norm{\Z}^2}+\frac{c^2\sigma^4}{\norm{\Z}^2}
\]
and take expectations: \( p\sigma^2-2c(p-2)\sigma^4\E\norm{\Z}^{-2}+c^2\sigma^4\E\norm{\Z}^{-2} \). The factor
\( c\{2(p-2)-c\} \) is positive iff \( 0<c<2(p-2) \) and largest at \( c=p-2 \), where it equals \( (p-2)^2 \).
(b) \( \Z/\sigma\sim\Normal_p(\boldsymbol{\uptheta}/\sigma,\I) \), so \( U=\norm{\Z}^2/\sigma^2\sim\chi^2(p,\gamma) \) (@def-qf-noncentral-chisq), and by @thm-qf-ncchisq(d), given \( K=k \), \( U\sim\chi^2(p+2k) \). For \( m\ge3 \),
\( \E\{1/\chi^2(m)\}=1/(m-2) \) (@exr-shr-inverse-moment), which gives the formula. At \( \gamma=0 \),
\( K=0 \) and the risk is \( p\sigma^2-(p-2)\sigma^2=2\sigma^2 \). Since
\( \E(1/U)=\int_0^\infty\Pr(U<1/s)\,ds \) and each \( \Pr(U<1/s) \) decreases strictly in \( \gamma \) (@prp-qf-ncchisq-monotone), \( \E(1/U) \) decreases strictly and the risk increases. As \( \gamma\to\infty \),
\( K\to\infty \) in probability and \( 1/(p-2+2K)\le1/(p-2) \), so \( \E\{1/(p-2+2K)\}\to0 \) by dominated convergence.
:::

The shrinkage factor is close to one when \( \norm{\Z}^2 \) is large compared with its null expectation
\( p\sigma^2 \), and small when the data are consistent with \( \boldsymbol{\uptheta}\approx\bzero \). Nothing
is special about the origin: shrinking towards any point fixed in advance also dominates \( \Z \),
with the largest gains near that point (@exr-shr-js-target).

The paradox is that the coordinates need not be related. If \( \theta_1 \) is a physical constant,
\( \theta_2 \) a crop yield and \( \theta_3 \) a batting average, the James–Stein estimate of the crop yield
depends on the other two measurements. The resolution is that the theorem concerns the *total*
risk; a single coordinate's risk can exceed \( \sigma^2 \) (@exr-shr-coordinate-risk).

::: {#exm-shr-js-risk}
[The risk function for ten means]

For \( p=10 \) and \( \sigma=1 \), the exact risk from @thm-shr-james-stein(b) is \( 2 \) at \( \boldsymbol{\uptheta}=\bzero \),
\( 4.784 \) at \( \norm{\boldsymbol{\uptheta}}^2=5 \), \( 7.587 \) at \( 20 \) and \( 9.257 \) at \( 80 \), always below the
constant risk \( 10 \) of \( \Z \) ([Figure 27.5.1](#fig-shr-js)). The fixed shrinkage
\( \Z/(1+0.5) \), which is ridge regression with \( \lambda=0.5 \) in an orthonormal design, has risk
\( (p+0.25\norm{\boldsymbol{\uptheta}}^2)/2.25 \) (@exr-shr-orthonormal-ridge). It beats James–Stein only on a middle range,
\( 5.9<\norm{\boldsymbol{\uptheta}}^2<36.6 \). At the origin it is worse (\( 10/2.25=4.44 \) against \( 2 \)), it
exceeds \( p \) once \( \norm{\boldsymbol{\uptheta}}^2>5p=50 \), and it grows without bound.
The dotted curve is the risk \( p\norm{\boldsymbol{\uptheta}}^2/(p+\norm{\boldsymbol{\uptheta}}^2) \) of the best linear
shrinkage \( f\Z \) for each \( \boldsymbol{\uptheta} \), the oracle of @prp-shr-canonical(c) with a common
factor; James–Stein stays within \( 2 \) of it without knowing \( \boldsymbol{\uptheta} \).
:::

::: {when-format="html"}
![**Figure 27.5.1.** Risk of estimators of ten normal means with \( \sigma=1 \), against
\( \norm{\boldsymbol{\uptheta}}^2 \): the unbiased estimator \( \Z \), James–Stein (exact, from the Poisson
mixture), its positive part (simulated, \( 4000 \) draws per point), a fixed proportional shrinkage,
and the best linear shrinkage for each \( \boldsymbol{\uptheta} \).](james_stein_risk.svg){#fig-shr-js width=72%}
:::

::: {when-format="pdf"}
![Risk of estimators of ten normal means with \( \sigma=1 \), against
\( \norm{\boldsymbol{\uptheta}}^2 \): the unbiased estimator \( \Z \), James–Stein (exact, from the Poisson
mixture), its positive part (simulated, \( 4000 \) draws per point), a fixed proportional shrinkage,
and the best linear shrinkage for each \( \boldsymbol{\uptheta} \).](james_stein_risk.pdf){width=72%}
:::

```{.python .run #cell-james-stein-risk}
import numpy as np
from scipy import stats
p = 10

def js_risk(gamma, p=p):
    """Exact risk of (1 - (p-2)/||Z||^2) Z when Z ~ N_p(theta, I) and gamma = ||theta||^2."""
    k = np.arange(0, 400)
    w = stats.poisson.pmf(k, gamma / 2)            # ||Z||^2 ~ chi^2(p, gamma): Poisson(gamma/2) mixture
    return p - (p - 2) ** 2 * np.sum(w / (p - 2 + 2 * k))

def js(z, positive=False):
    factor = 1 - (len(z) - 2) / np.sum(z ** 2)
    return (max(factor, 0.0) if positive else factor) * z

rng = np.random.default_rng(27)
for gamma in [0.0, 5.0, 20.0, 80.0]:
    theta = np.sqrt(gamma / p) * np.ones(p)
    Zs = theta + rng.normal(size=(2_000, p))
    loss_js = np.mean([np.sum((js(z) - theta) ** 2) for z in Zs])
    loss_pp = np.mean([np.sum((js(z, True) - theta) ** 2) for z in Zs])
    print(f"||theta||^2 = {gamma:5.1f}: exact JS risk {js_risk(gamma):.3f}, "
          f"simulated {loss_js:.3f}, positive part {loss_pp:.3f}  (MLE: {p})")
```

## The positive part

When \( \norm{\Z}^2<(p-2)\sigma^2 \) the James–Stein factor is negative, and the estimator reverses the
sign of every coordinate, which is absurd. Truncating the factor at zero helps.

::: {#prp-shr-positive-part}
[Positive-part James–Stein]

In the setting of @thm-shr-james-stein, the estimator
\[
\hat{\boldsymbol{\uptheta}}^{\mathrm{JS}+}=\Bigl(1-\frac{(p-2)\sigma^2}{\norm{\Z}^2}\Bigr)_+\Z
\]
has strictly smaller risk than \( \hat{\boldsymbol{\uptheta}}^{\mathrm{JS}} \) for every \( \boldsymbol{\uptheta} \).
:::

::: {.proof}
Write \( \hat{\boldsymbol{\uptheta}}^{\mathrm{JS}}=h\Z \) with \( h=1-(p-2)\sigma^2/\norm{\Z}^2 \), and let
\( \psi=\max(-h,0)\ge0 \), a function of \( \norm{\Z} \). Where \( h\ge0 \) the estimators agree. Where \( h<0 \),
the positive part is \( \bzero \), and
\[
\norm{h\Z-\boldsymbol{\uptheta}}^2-\norm{\boldsymbol{\uptheta}}^2=h^2\norm{\Z}^2-2h\,\boldsymbol{\uptheta}\T\Z=\psi^2\norm{\Z}^2+2\psi\,\boldsymbol{\uptheta}\T\Z .
\]
So the difference of risks is \( \E[\psi^2\norm{\Z}^2]+2\E[\psi\,\boldsymbol{\uptheta}\T\Z] \). The first term is positive, because
\( \Pr(h<0)=\Pr\{\norm{\Z}^2<(p-2)\sigma^2\}>0 \). For the second, choose an orthogonal \( \Q \) with
\( \Q\boldsymbol{\uptheta}=m\mathbf{e}_1 \), \( m=\norm{\boldsymbol{\uptheta}} \). Then \( \mathbf{w}=\Q\Z\sim\Normal_p(m\mathbf{e}_1,\sigma^2\I) \),
\( \norm{\mathbf{w}}=\norm{\Z} \) and \( \boldsymbol{\uptheta}\T\Z=mw_1 \). Given \( R=\sum_{i\ge2}w_i^2 \), which is independent of
\( w_1\sim\Normal(m,\sigma^2) \),
\[
\begin{aligned}
\E[\psi\,w_1\mid R]&=\int_0^\infty\psi\bigl(\sqrt{t^2+R}\bigr)\,t\,\bigl\{\varphi_\sigma(t-m)-\varphi_\sigma(t+m)\bigr\}\,dt\\
&\ge0,
\end{aligned}
\]
where \( \varphi_\sigma \) is the \( \Normal(0,\sigma^2) \) density: the integrand pairs \( t \) with \( -t \), and
\( (t-m)^2\le(t+m)^2 \) for \( t,m\ge0 \). (The expectations are finite since \( \psi\le(p-2)\sigma^2/\norm{\Z}^2 \).) Hence
\( \E[\psi\,\boldsymbol{\uptheta}\T\Z]\ge0 \) and the risk difference is positive.
:::

The positive-part estimator is itself inadmissible, but improvements on it are tiny (Lehmann and
Casella 1998, chapter 5). For \( p=10 \) its risk at the origin is
\( \E[\{1-(p-2)/U\}_+^2U]=1.258 \) with \( U\sim\chi^2(10) \), against \( 2 \) for James–Stein, and the two
curves in [Figure 27.5.1](#fig-shr-js) merge as negative factors become rare.

## The empirical Bayes explanation

Efron and Morris (1973) explained the paradox by drawing the means themselves from a normal
distribution.

::: {#prp-shr-empirical-bayes}
[James–Stein as empirical Bayes]

Let \( \boldsymbol{\uptheta}\sim\Normal_p(\bzero,\omega^2\I) \) and \( \Z\mid\boldsymbol{\uptheta}\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \),
\( p\ge3 \), and put \( B=\sigma^2/(\sigma^2+\omega^2) \).

::: {.enumerate options="label=(\alph*)"}
1. The posterior mean is \( \E(\boldsymbol{\uptheta}\mid\Z)=(1-B)\Z \), and its Bayes risk is
   \( \E\norm{(1-B)\Z-\boldsymbol{\uptheta}}^2=p\sigma^2(1-B) \).

2. Marginally \( \Z\sim\Normal_p(\bzero,(\sigma^2+\omega^2)\I) \), and \( (p-2)\sigma^2/\norm{\Z}^2 \) is an unbiased
   estimator of \( B \).

3. The Bayes risk of \( \hat{\boldsymbol{\uptheta}}^{\mathrm{JS}} \) is \( p\sigma^2(1-B)+2\sigma^2B \).
:::
:::

::: {.proof}
(a) \( (\boldsymbol{\uptheta},\Z) \) is jointly normal with \( \Cov(\boldsymbol{\uptheta},\Z)=\omega^2\I \) and
\( \Cov(\Z)=(\sigma^2+\omega^2)\I \). By @thm-mvn-conditional,
\( \boldsymbol{\uptheta}\mid\Z\sim\Normal_p\bigl(\omega^2(\sigma^2+\omega^2)^{-1}\Z,\ \sigma^2\omega^2(\sigma^2+\omega^2)^{-1}\I\bigr) \); the mean is
\( (1-B)\Z \), and the Bayes risk is the expected trace of the posterior covariance,
\( p\sigma^2\omega^2/(\sigma^2+\omega^2)=p\sigma^2(1-B) \).
(b) \( \Z=\boldsymbol{\uptheta}+(\Z-\boldsymbol{\uptheta}) \) is a sum of independent normal vectors. So
\( \norm{\Z}^2/(\sigma^2+\omega^2)\sim\chi^2(p) \), and \( \E\{(p-2)\sigma^2/\norm{\Z}^2\}=(p-2)\sigma^2/\{(\sigma^2+\omega^2)(p-2)\}=B \)
by @exr-shr-inverse-moment.
(c) Average @eq-shr-js over the prior: the Bayes risk is \( p\sigma^2-(p-2)^2\sigma^4\E\norm{\Z}^{-2} \) with the
marginal expectation \( \E\norm{\Z}^{-2}=1/\{(\sigma^2+\omega^2)(p-2)\} \), which gives
\( p\sigma^2-(p-2)\sigma^2B=p\sigma^2(1-B)+2\sigma^2B \).
:::

James–Stein plugs an unbiased estimate of \( B \), based on the spread of the \( p \) observations,
into the Bayes rule, and part (c) shows that not knowing \( \omega^2 \) costs at most \( 2\sigma^2 \) in total.
This is the **empirical Bayes** idea: the ensemble of parameters supplies its own prior.
@thm-shr-james-stein says more, that the gain is guaranteed even if no such prior exists.

In an orthonormal design ridge regression is the Bayes rule \( \Z/(1+\lambda) \) with
\( \omega^2=\sigma^2/\lambda \) (@thm-shr-ridge(e)). A fixed \( \lambda \) is a guess at the prior variance, and a
bad guess is punished without limit, as the fixed-shrinkage curve in [Figure 27.5.1](#fig-shr-js)
shows; choosing \( \lambda \) from the data is the regression version of James–Stein. With a general design the \( g \)-prior of
@prp-opt-shrinkage(c) shrinks \( \hbeta \) by the common factor \( g/(1+g) \), and estimating that factor
gives the next corollary.

## Shrinking a regression with unknown variance

In regression \( \sigma^2 \) is unknown but estimated independently of \( \hbeta \), which is enough.

::: {#cor-shr-js-regression}
[James–Stein for regression]

Assume the normal linear model with \( \X \) of full column rank \( p\ge3 \), let \( \text{SSE} \) have \( m \) degrees of
freedom (so \( \text{SSE}/\sigma^2\sim\chi^2(m) \), independent of \( \hbeta \)), and for \( c>0 \) let
\[
\tilde{\bbeta}_c=\Bigl(1-\frac{c\,\text{SSE}}{\hbeta\T\X\T\X\hbeta}\Bigr)\hbeta .
\]
Then, with \( \Z=(\X\T\X)^{1/2}\hbeta \),
\[
\E\bigl[(\tilde{\bbeta}_c-\bbeta)\T\X\T\X(\tilde{\bbeta}_c-\bbeta)\bigr]
=p\sigma^2-c\,m\,\{2(p-2)-c(m+2)\}\,\sigma^4\,\E\norm{\Z}^{-2}.
\]
So \( \tilde{\bbeta}_c \) dominates least squares in this loss iff
\( 0<c<2(p-2)/(m+2) \), and the best constant is \( c=(p-2)/(m+2) \). In terms of the \( F \) statistic
\( F=(\hbeta\T\X\T\X\hbeta/p)/(\text{SSE}/m) \) for the hypothesis \( \bbeta=\bzero \), the shrinkage factor is
\( 1-cm/(pF) \).
:::

::: {.proof}
By @thm-opt-sampling, \( \hbeta\sim\Normal_p(\bbeta,\sigma^2(\X\T\X)^{-1}) \) and \( \text{SSE} \) is independent of it.
So \( \Z\sim\Normal_p(\boldsymbol{\uptheta},\sigma^2\I) \) with \( \boldsymbol{\uptheta}=(\X\T\X)^{1/2}\bbeta \), and the loss is
\( \norm{(1-a/\norm{\Z}^2)\Z-\boldsymbol{\uptheta}}^2 \) with \( a=c\,\text{SSE} \). Conditionally on \( \text{SSE} \), \( a \) is a
constant and the computation in the proof of @thm-shr-james-stein(a), with \( c\sigma^2 \) replaced by \( a \),
gives the conditional risk \( p\sigma^2-\{2a(p-2)\sigma^2-a^2\}\E\norm{\Z}^{-2} \). Now \( \E(a)=cm\sigma^2 \) and
\( \E(a^2)=c^2\sigma^4\E\{\chi^2(m)^2\}=c^2m(m+2)\sigma^4 \) (@thm-qf-ncchisq(b)). Taking expectations gives the
formula, and \( c\{2(p-2)-c(m+2)\} \) is positive iff \( 0<c<2(p-2)/(m+2) \) and largest at the midpoint. The
last statement is \( \text{SSE}/(\hbeta\T\X\T\X\hbeta)=m/(pF) \).
:::

The factor \( 1-(p-2)m/\{(m+2)pF\} \) is close to \( 1-1/F \) when \( p \) and \( m \) are large: a regression
that clearly fits is barely shrunk, and one indistinguishable from noise is shrunk hard. The loss
concerns fitted values, not individual coefficients. Shrinking towards the fit of a reduced model,
using its \( F \) statistic, works in the same way when the full model has at least three more
parameters (@exr-shr-js-subspace).

::: {#exm-shr-longley-js}
[James–Stein for Longley's regression]

For Longley's six slopes there are \( m=9 \) residual degrees of freedom and the overall \( F \) statistic is \( 330.3 \). The James–Stein factor with
\( c=(p-2)/(m+2) \) is \( 0.99835 \): the regression fits far too well for shrinkage towards zero to
help. The script therefore also simulates responses at Longley's design with the least squares
estimates divided by \( 60 \) as the truth, a weak signal whose
noncentrality \( \bbeta\T\X\T\X\bbeta/\sigma^2 \) is about \( 0.6 \), with \( \sigma \) equal to Longley's \( s \). Over \( 20000 \)
simulated data sets the ratio of the James–Stein prediction risk to that of least squares is
\( 0.501 \).
:::

```{.python .run #cell-james-stein-longley}
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data
Zr = data.drop(columns="TOTEMP").to_numpy()
y = data["TOTEMP"].to_numpy()
n, k = Zr.shape
Xc = Zr - Zr.mean(axis=0)                              # centred slopes; the intercept is ybar
yc = y - y.mean()
b = np.linalg.lstsq(Xc, yc, rcond=None)[0]
m = n - k - 1                                          # residual degrees of freedom
sse = np.sum((yc - Xc @ b) ** 2)
fit_ss = b @ Xc.T @ Xc @ b                             # b^T X^T X b, the regression sum of squares
F = (fit_ss / k) / (sse / m)
factor = 1 - (k - 2) / (m + 2) * sse / fit_ss
print(f"F = {F:.1f}; James-Stein factor {factor:.5f} = 1 - (k-2) m / ((m+2) k F)")
```

## Exercises

### A. Check your understanding

::: {#exr-shr-inverse-moment}
[A1]

Show that \( \E\{1/\chi^2(m)\}=1/(m-2) \) for \( m\ge3 \) and that the expectation is infinite for
\( m\le2 \). Deduce that the risk of the James–Stein estimator at \( \boldsymbol{\uptheta}=\bzero \) is \( 2\sigma^2 \)
for every \( p\ge3 \), however large \( p \) is.
:::

::: {.solution}
With the \( \chi^2(m) \) density \( x^{m/2-1}e^{-x/2}/\{2^{m/2}\Gamma(m/2)\} \),
\[
\E\,\frac1{\chi^2(m)}=\frac{\int_0^\infty x^{m/2-2}e^{-x/2}\,dx}{2^{m/2}\Gamma(m/2)}
=\frac{2^{m/2-1}\Gamma(m/2-1)}{2^{m/2}\Gamma(m/2)}=\frac1{2(m/2-1)}=\frac1{m-2},
\]
where the integral converges at zero iff \( m/2-2>-1 \), that is, \( m>2 \). At \( \boldsymbol{\uptheta}=\bzero \),
\( \norm{\Z}^2/\sigma^2\sim\chi^2(p) \), so @eq-shr-js gives \( p\sigma^2-(p-2)^2\sigma^2/(p-2)=2\sigma^2 \).
:::

### B. Practice

::: {#exr-shr-js-target}
[B1]

(a) Show that for a fixed \( \boldsymbol{\uptheta}_0 \), the estimator
\( \boldsymbol{\uptheta}_0+(1-(p-2)\sigma^2/\norm{\Z-\boldsymbol{\uptheta}_0}^2)(\Z-\boldsymbol{\uptheta}_0) \) has risk
\( p\sigma^2-(p-2)^2\sigma^4\E\norm{\Z-\boldsymbol{\uptheta}_0}^{-2} \). (b) Lindley's estimator shrinks towards the grand
mean: \( \hat{\theta}_i=\bar Z+(1-(p-3)\sigma^2/S)(Z_i-\bar Z) \) with \( S=\sum_i(Z_i-\bar Z)^2 \). Show that it dominates
\( \Z \) when \( p\ge4 \).
:::

::: {.solution}
(a) \( \Z-\boldsymbol{\uptheta}_0\sim\Normal_p(\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0,\sigma^2\I) \), and the estimation error is that of
the James–Stein estimator of \( \boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0 \) from \( \Z-\boldsymbol{\uptheta}_0 \); apply @eq-shr-js.
(b) Choose an orthogonal \( \Q \) whose first row is \( \bone\T/\sqrt p \). In the coordinates \( \mathbf{w}=\Q\Z \), the first coordinate is \( \sqrt p\bar Z \) and the remaining
\( p-1 \) are a \( \Normal_{p-1} \) vector with squared length \( S \). The estimator leaves the first coordinate alone
and applies James–Stein in dimension \( p-1\ge3 \) to the others, so its risk is \( \sigma^2 \) plus less
than \( (p-1)\sigma^2 \) (@thm-shr-james-stein).
:::

::: {#exr-shr-coordinate-risk}
[B2]

Let \( p \) be large and \( \boldsymbol{\uptheta}=(\sqrt p\,\sigma,0,\dots,0)\T \). Using the approximation
\( \norm{\Z}^2\approx\norm{\boldsymbol{\uptheta}}^2+p\sigma^2 \), show that the James–Stein estimate of \( \theta_1 \) has mean
squared error of roughly \( p\sigma^2/4 \), far above the \( \sigma^2 \) of \( Z_1 \). Why does this not
contradict @thm-shr-james-stein?
:::

::: {.solution}
With \( \norm{\Z}^2\approx2p\sigma^2 \), the factor is about \( 1-(p-2)/(2p)\approx\tfrac12 \), so
\( \hat{\theta}_1\approx Z_1/2 \), with bias \( -\theta_1/2 \) and squared bias \( p\sigma^2/4 \). Halving the other \( p-1 \) coordinates,
whose means are zero, cuts each of their risks from \( \sigma^2 \) to about \( \sigma^2/4 \). The total is roughly \( p\sigma^2/4+(p-1)\sigma^2/4\approx p\sigma^2/2<p\sigma^2 \).
One coordinate pays heavily for the others' gains; the theorem only controls the sum.
:::

::: {#exr-shr-eb-ratio}
[B3]

In @prp-shr-empirical-bayes, show that the ratio of the Bayes risk of James–Stein to that of the Bayes
rule is \( 1+2\sigma^2/(p\omega^2) \). When is the cost of estimating \( \omega^2 \) negligible?
:::

### C. Going deeper

::: {#exr-shr-js-subspace}
[C1]

Let \( \bP \) be the orthogonal projection onto a fixed \( q \)-dimensional subspace of \( \Real^p \), with
\( p-q\ge3 \). Show that \( \bP\Z+\bigl(1-(p-q-2)\sigma^2/\norm{(\I-\bP)\Z}^2\bigr)(\I-\bP)\Z \) dominates \( \Z \), and
translate this into shrinking a regression fit towards a nested reduced model using its \( F \) statistic.
:::

::: {#exr-shr-sure}
[C2]

Let \( \hat{\boldsymbol{\uptheta}}=\Z+\mathbf{g}(\Z) \) with \( \mathbf{g} \) satisfying the conditions of @lem-shr-stein. Show that
\[
\hat R=p\sigma^2+\norm{\mathbf{g}(\Z)}^2+2\sigma^2\sum_i\frac{\partial g_i}{\partial z_i}(\Z)
\]
is unbiased for the risk of \( \hat{\boldsymbol{\uptheta}} \) (Stein's unbiased risk estimate). Compute it for
James–Stein and for soft thresholding \( \hat\theta_i=S_\lambda(Z_i) \); how could it choose \( \lambda \)?
:::
