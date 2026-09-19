# Bayesian thread I: the conjugate linear model

The previous sections judged estimators by their behaviour over repeated samples, for each fixed value of the
parameters. The Bayesian approach treats \( (\bbeta,\sigma^2) \) as uncertain quantities with a probability
distribution. It combines a **prior** distribution, expressing what is known before the data, with the
likelihood of [Section 7.3](03-maximum-likelihood.html), and the result is a **posterior** distribution. For the
normal linear model there is a family of priors for which the posterior is available in closed form, and it
connects directly to the rest of the chapter. The posterior mean is a matrix-weighted combination of the prior mean
and the least squares estimate. It is biased in the frequentist sense, which is how it escapes the
Gauss–Markov bound. And as the prior becomes flat, the posterior reproduces least squares and the
\( t \) intervals of [Section 7.5](05-sampling-distributions.html) exactly. This is the first of several Bayesian
threads in the book. [Section 12.6](../ch12-intervals-and-bands/06-bayes.html) continues it with credible regions and predictive distributions, Chapter 27
with shrinkage priors, and Chapter 39 with generalized linear models.

## The normal-inverse-gamma family

::: {#def-opt-nig}
[Inverse gamma and normal-inverse-gamma distributions]

::: {.enumerate options="label=(\alph*)"}
1. A positive random variable \( V \) has the **inverse gamma** distribution \( \text{IG}(a,b) \), with
   \( a,b>0 \), if \( 1/V \) has the gamma distribution with shape \( a \) and rate \( b \). Its density is
   \[
\pi(v)=\frac{b^a}{\Gamma(a)}\,v^{-a-1}e^{-b/v},\qquad v>0 .
\]

2. The pair \( (\bbeta,\sigma^2) \) has the **normal-inverse-gamma** distribution \( \text{NIG}(\mathbf{m},\V,a,b) \), with
   \( \mathbf{m}\in\Real^p \), \( \V \) positive definite and \( a,b>0 \), if \( \sigma^2\sim\text{IG}(a,b) \) and, given \( \sigma^2 \),
   \( \bbeta\sim\Normal_p(\mathbf{m},\sigma^2\V) \). Its joint density is
   \[
\begin{aligned}
\pi(\bbeta,\sigma^2)&=\frac{b^a}{\Gamma(a)(2\pi)^{p/2}\det(\V)^{1/2}}\,(\sigma^2)^{-(a+p/2)-1}\\
&\qquad\times\exp\Bigl\{-\frac{2b+(\bbeta-\mathbf{m})\T\V^{-1}(\bbeta-\mathbf{m})}{2\sigma^2}\Bigr\}.
\end{aligned}
\]{#eq-opt-nig-density}

:::
:::

The prior on \( \bbeta \) is scaled by \( \sigma^2 \): prior uncertainty about the coefficients is measured in units
of the noise level. This is a modelling choice, and it is what makes the family conjugate. It also makes
\( \V \) dimensionless relative to the likelihood's \( (\X\T\X)^{-1} \), which is the comparison that matters below.
The marginal distribution of \( \bbeta \) is a multivariate \( t \).

::: {#def-opt-mvt}
[Multivariate \( t \)]

Let \( \Z\sim\Normal_p(\bzero,\bS) \) and \( W\sim\chi^2(\nu) \) be independent, with \( \bS \) positive definite. The law of
\( \mathbf{m}+\Z/\sqrt{W/\nu} \) is the **multivariate \( t \) distribution** \( t_\nu(\mathbf{m},\bS) \).
:::

For \( p=1 \) and \( \bS=1 \), this is the \( t(\nu) \) distribution of @def-qf-noncentral-t shifted by \( \mathbf{m} \). For every
\( \mathbf{c}\ne\bzero \), \( \mathbf{c}\T\bu \) with \( \bu\sim t_\nu(\mathbf{m},\bS) \) is \( \mathbf{c}\T\mathbf{m}+\sqrt{\mathbf{c}\T\bS\mathbf{c}}\,T \) with
\( T\sim t(\nu) \), because \( \mathbf{c}\T\Z/\sqrt{\mathbf{c}\T\bS\mathbf{c}}\sim\Normal(0,1) \) is independent of \( W \).

::: {#lem-opt-nig-marginal}
[Marginals of the normal-inverse-gamma]

If \( (\bbeta,\sigma^2)\sim\text{NIG}(\mathbf{m},\V,a,b) \), then \( 2b/\sigma^2\sim\chi^2(2a) \),
\( \E(\sigma^2)=b/(a-1) \) for \( a>1 \), and
\[
\bbeta\sim t_{2a}\bigl(\mathbf{m},\ (b/a)\V\bigr).
\]
In particular, for every \( \blambda\ne\bzero \),
\( (\blambda\T\bbeta-\blambda\T\mathbf{m})/\sqrt{(b/a)\blambda\T\V\blambda}\sim t(2a) \).
:::

::: {.proof}
Put \( W=2b/\sigma^2 \). Changing variables \( v=2b/w \) in the inverse gamma density, with \( \lvert dv/dw\rvert=2b/w^2 \),
gives the density of \( W \):
\[
\frac{b^a}{\Gamma(a)}\Bigl(\frac{2b}{w}\Bigr)^{-a-1}e^{-w/2}\frac{2b}{w^2}
=\frac{w^{a-1}e^{-w/2}}{2^a\Gamma(a)},\qquad w>0,
\]
which is the \( \chi^2(2a) \) density (@eq-qf-chisq-density). Then
\( \E\sigma^2=2b\,\E(1/W)=2b/(2a-2) \), since \( \E(1/W)=1/(\nu-2) \) for \( W\sim\chi^2(\nu) \), \( \nu>2 \). Next let
\( \Z_0=(\bbeta-\mathbf{m})/\sigma \). Given \( \sigma^2 \), \( \Z_0\sim\Normal_p(\bzero,\V) \), a law that does not depend on
\( \sigma^2 \). So \( \Z_0 \) is independent of \( \sigma^2 \), and hence of \( W \), and
\[
\bbeta=\mathbf{m}+\sigma\Z_0=\mathbf{m}+\sqrt{\frac{2b}W}\,\Z_0=\mathbf{m}+\frac{\sqrt{b/a}\,\Z_0}{\sqrt{W/(2a)}} .
\]
Since \( \sqrt{b/a}\,\Z_0\sim\Normal_p(\bzero,(b/a)\V) \), this is \( t_{2a}(\mathbf{m},(b/a)\V) \) by @def-opt-mvt. The last
statement is the remark after @def-opt-mvt.
:::

## The conjugate update

::: {#thm-opt-bayes-conjugate}
[Conjugate posterior]

Let \( \y \) follow the normal linear model @eq-opt-normal-model, with \( \X \) of any rank, and let the prior be
\( (\bbeta,\sigma^2)\sim\text{NIG}(\mathbf{m}_0,\V_0,a_0,b_0) \). Then the posterior is
\( \text{NIG}(\mathbf{m}_n,\V_n,a_n,b_n) \), where
\[
\begin{aligned}
\V_n&=(\V_0^{-1}+\X\T\X)^{-1},\qquad \mathbf{m}_n=\V_n(\V_0^{-1}\mathbf{m}_0+\X\T\y),\\
a_n&=a_0+\frac n2,
\end{aligned}
\]{#eq-opt-posterior-params}

\[
\begin{aligned}
b_n&=b_0+\tfrac12\bigl(\y\T\y+\mathbf{m}_0\T\V_0^{-1}\mathbf{m}_0-\mathbf{m}_n\T\V_n^{-1}\mathbf{m}_n\bigr)\\
&=b_0+\tfrac12\bigl\{\norm{\y-\X\mathbf{m}_n}^2+(\mathbf{m}_n-\mathbf{m}_0)\T\V_0^{-1}(\mathbf{m}_n-\mathbf{m}_0)\bigr\}.
\end{aligned}
\]{#eq-opt-posterior-b}

Consequently \( \bbeta\mid\y\sim t_{2a_n}\bigl(\mathbf{m}_n,(b_n/a_n)\V_n\bigr) \), \( \sigma^2\mid\y\sim\text{IG}(a_n,b_n) \), and
\( \E(\bbeta\mid\y)=\mathbf{m}_n \).
:::

::: {.proof}
Write \( \mathbf{P}_0=\V_0^{-1} \) and \( \mathbf{P}_n=\mathbf{P}_0+\X\T\X \). The matrix \( \mathbf{P}_n \) is positive definite, because
\( \mathbf{P}_0 \) is and \( \X\T\X \) is nonnegative definite. This holds whatever the rank of \( \X \). Expanding both sides and
using \( \mathbf{P}_n\mathbf{m}_n=\mathbf{P}_0\mathbf{m}_0+\X\T\y \) gives the completion of the square
\[
\norm{\y-\X\bbeta}^2+(\bbeta-\mathbf{m}_0)\T\mathbf{P}_0(\bbeta-\mathbf{m}_0)=(\bbeta-\mathbf{m}_n)\T\mathbf{P}_n(\bbeta-\mathbf{m}_n)+C,
\]
with \( C=\y\T\y+\mathbf{m}_0\T\mathbf{P}_0\mathbf{m}_0-\mathbf{m}_n\T\mathbf{P}_n\mathbf{m}_n \). Indeed, the left side is
\( \bbeta\T\mathbf{P}_n\bbeta-2\bbeta\T(\X\T\y+\mathbf{P}_0\mathbf{m}_0)+\y\T\y+\mathbf{m}_0\T\mathbf{P}_0\mathbf{m}_0 \), and
\( (\bbeta-\mathbf{m}_n)\T\mathbf{P}_n(\bbeta-\mathbf{m}_n)=\bbeta\T\mathbf{P}_n\bbeta-2\bbeta\T\mathbf{P}_n\mathbf{m}_n+\mathbf{m}_n\T\mathbf{P}_n\mathbf{m}_n \).
Putting \( \bbeta=\mathbf{m}_n \) in the identity shows that
\( C=\norm{\y-\X\mathbf{m}_n}^2+(\mathbf{m}_n-\mathbf{m}_0)\T\mathbf{P}_0(\mathbf{m}_n-\mathbf{m}_0)\ge0 \), which is the second form of @eq-opt-posterior-b.

By Bayes' theorem the posterior density is proportional, as a function of \( (\bbeta,\sigma^2) \), to prior times
likelihood. Using @eq-opt-nig-density and dropping factors free of the parameters,
\[
\begin{aligned}
\pi(\bbeta,\sigma^2\mid\y)&\propto(\sigma^2)^{-(a_0+p/2)-1}
e^{-\{2b_0+(\bbeta-\mathbf{m}_0)\T\mathbf{P}_0(\bbeta-\mathbf{m}_0)\}/(2\sigma^2)}\\
&\qquad\times(\sigma^2)^{-n/2}e^{-\norm{\y-\X\bbeta}^2/(2\sigma^2)} .
\end{aligned}
\]
By the completion of the square, the exponent is \( -\{2b_0+C+(\bbeta-\mathbf{m}_n)\T\mathbf{P}_n(\bbeta-\mathbf{m}_n)\}/(2\sigma^2) \),
and \( 2b_0+C=2b_n \). The power of \( \sigma^2 \) is \( -(a_0+n/2+p/2)-1=-(a_n+p/2)-1 \). So the posterior is
proportional to the \( \text{NIG}(\mathbf{m}_n,\V_n,a_n,b_n) \) density @eq-opt-nig-density. Two densities that are
proportional are equal, since both integrate to one. The marginals follow from @lem-opt-nig-marginal. The
posterior mean of \( \bbeta \) is \( \mathbf{m}_n \): by the tower property, \( \E(\bbeta\mid\y)=\E\{\E(\bbeta\mid\sigma^2,\y)\mid\y\}=\mathbf{m}_n \).
:::

The update has a simple structure. **Precisions add**: the posterior precision (per unit \( \sigma^2 \)) is the
prior precision \( \V_0^{-1} \) plus the data precision \( \X\T\X \). **Shape parameters count observations**: each
observation adds \( 1/2 \) to \( a \). **Scale parameters accumulate squared discrepancies**: \( b_n \) adds half of the
residual sum of squares at the posterior mean and half of the prior-weighted distance between the posterior and
prior means. Because \( \mathbf{P}_n \) is positive definite even when \( \X \) is rank deficient, the posterior is proper for
every design. The prior supplies information about the directions of \( \Null(\X) \) that the data cannot, and there
the posterior equals the prior (@exr-opt-bayes-nonidentified).

## The posterior mean as shrinkage

::: {#prp-opt-shrinkage}
[Posterior mean as a combination]

In @thm-opt-bayes-conjugate, let \( \X \) have full column rank and put \( \W=\V_n\V_0^{-1} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{m}_n=\W\mathbf{m}_0+(\I-\W)\hbeta \), with \( \I-\W=\V_n\X\T\X \). Equivalently
   \( \V_n^{-1}\mathbf{m}_n=\V_0^{-1}\mathbf{m}_0+(\X\T\X)\hbeta \): the posterior mean is the precision-weighted average of the
   prior mean and the least squares estimate.

2. \( b_n=b_0+\tfrac12\bigl\{\text{SSE}+(\hbeta-\mathbf{m}_0)\T\bigl(\V_0+(\X\T\X)^{-1}\bigr)^{-1}(\hbeta-\mathbf{m}_0)\bigr\} \).

3. *(\( g \)-prior.)* If \( \V_0=g(\X\T\X)^{-1} \) with \( g>0 \), then \( \W=\I/(1+g) \) and
   \( \mathbf{m}_n=(\mathbf{m}_0+g\hbeta)/(1+g) \).

4. *(Ridge.)* If \( \mathbf{m}_0=\bzero \) and \( \V_0=k^{-1}\I \) with \( k>0 \), then \( \mathbf{m}_n=(\X\T\X+k\I)^{-1}\X\T\y \).
:::
:::

::: {.proof}
(a) \( \X\T\y=\X\T\X\hbeta \) by the normal equations (@thm-proj-normal-equations), so
\( \mathbf{m}_n=\V_n\V_0^{-1}\mathbf{m}_0+\V_n\X\T\X\hbeta \), and \( \V_n\V_0^{-1}+\V_n\X\T\X=\V_n\V_n^{-1}=\I \).
(b) Write \( \mathbf{P}_0=\V_0^{-1} \) and \( \mathbf{K}=\X\T\X \). The completion of the square in the proof of
@thm-opt-bayes-conjugate shows that \( 2(b_n-b_0) \) is the minimum over \( \bbeta \) of
\( \norm{\y-\X\bbeta}^2+(\bbeta-\mathbf{m}_0)\T\mathbf{P}_0(\bbeta-\mathbf{m}_0) \), because the remaining term
\( (\bbeta-\mathbf{m}_n)\T\mathbf{P}_n(\bbeta-\mathbf{m}_n) \) is nonnegative and vanishes at \( \mathbf{m}_n \). By @eq-proj-distance-split,
\( \norm{\y-\X\bbeta}^2=\text{SSE}+(\bbeta-\hbeta)\T\mathbf{K}(\bbeta-\hbeta) \). Put \( \mathbf{d}=\hbeta-\mathbf{m}_0 \) and
\( \bbeta=\mathbf{m}_0+\bu \). The quantity to minimize becomes \( \text{SSE}+(\bu-\mathbf{d})\T\mathbf{K}(\bu-\mathbf{d})+\bu\T\mathbf{P}_0\bu \), whose
minimum over \( \bu \), attained at \( \bu=(\mathbf{K}+\mathbf{P}_0)^{-1}\mathbf{K}\mathbf{d} \), is
\( \text{SSE}+\mathbf{d}\T\{\mathbf{K}-\mathbf{K}(\mathbf{K}+\mathbf{P}_0)^{-1}\mathbf{K}\}\mathbf{d} \). By @thm-mat-woodbury,
\( \mathbf{K}-\mathbf{K}(\mathbf{K}+\mathbf{P}_0)^{-1}\mathbf{K}=(\mathbf{K}^{-1}+\mathbf{P}_0^{-1})^{-1}=(\V_0+(\X\T\X)^{-1})^{-1} \).
(c) With \( \V_0^{-1}=\X\T\X/g \), \( \V_n^{-1}=(1+1/g)\X\T\X \), so
\( \W=\V_n\V_0^{-1}=\{g/(1+g)\}(\X\T\X)^{-1}\X\T\X/g=\I/(1+g) \). (d) Substitute into @eq-opt-posterior-params.
:::

Part (a) is the Bayesian counterpart of Gauss–Markov. The posterior mean is linear in \( \y \), and it is biased
unless \( \mathbf{m}_0=\bbeta \): \( \E(\mathbf{m}_n)=\W\mathbf{m}_0+(\I-\W)\bbeta \). So the Gauss–Markov theorem does not apply to it,
and it can have smaller mean squared error than least squares when the prior mean is close to the truth.
Part (d) identifies ridge regression, whose mean squared error was studied in @exr-opt-ridge, as a posterior
mean. Chapter 27 develops the connection. When the prior is weak relative to the data, \( \W \) is close to \( \mathbf{0} \)
and the posterior mean is close to \( \hbeta \). When it is strong, \( \W \) is close to \( \I \) and the data barely move it.

In general \( \W \) is a full matrix, not a multiple of the identity. Shrinkage then happens along the eigenvectors
of the problem, not coordinate by coordinate. A single coefficient's posterior mean need not lie between its
prior mean and its least squares estimate, as the next example shows. Only the \( g \)-prior of part (c) shrinks
every coordinate by the same factor. Its prior covariance copies the shape of the sampling covariance of \( \hbeta \),
and that makes it convenient for model comparison. Zellner (1986) introduced it.

## The flat-prior limit

What happens when the prior carries no information? Letting \( \V_0^{-1}\to\mathbf{0} \) in @eq-opt-posterior-params
gives \( \mathbf{m}_n\to\hbeta \) when \( \X \) has full rank. The shape parameter needs more care. The NIG density
@eq-opt-nig-density carries the factor \( (\sigma^2)^{-(a_0+p/2)-1} \), so letting only \( a_0,b_0\to0 \) leads to the
prior \( (\sigma^2)^{-p/2-1} \), with \( a_n=n/2 \) and \( t \) marginals on \( n \) degrees of freedom rather than
\( n-p \). The improper prior \( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \), uniform in \( \bbeta \) and in \( \log\sigma \),
is the limit \( \V_0^{-1}\to\mathbf{0} \), \( a_0\to-p/2 \), \( b_0\to0 \). A negative shape \( a_0=-p/2 \) lies outside the
proper NIG family, so this prior cannot be reached inside @thm-opt-bayes-conjugate and has to be handled directly. The
result is clean.

::: {#cor-opt-flat-prior}
[The flat prior reproduces least squares]

Let \( \rank(\X)=p<n \), let \( \text{SSE}>0 \), and take the improper prior \( \pi(\bbeta,\sigma^2)\propto1/\sigma^2 \). Then the
posterior is proper and equals \( \text{NIG}\bigl(\hbeta,(\X\T\X)^{-1},(n-p)/2,\text{SSE}/2\bigr) \). In particular,
for every \( \blambda\ne\bzero \),
\[
\frac{\blambda\T\bbeta-\blambda\T\hbeta}{s\sqrt{\blambda\T(\X\T\X)^{-1}\blambda}}\ \Big|\ \y\ \sim\ t(n-p),
\]
so the \( 1-\alpha \) equal-tailed credible interval for \( \blambda\T\bbeta \) is the \( t \) interval @eq-opt-t-interval.
:::

::: {.proof}
By @eq-proj-distance-split, \( \norm{\y-\X\bbeta}^2=\text{SSE}+(\bbeta-\hbeta)\T\X\T\X(\bbeta-\hbeta) \). So prior times
likelihood is proportional to
\[
\begin{aligned}
&(\sigma^2)^{-1-n/2}\exp\Bigl\{-\frac{\text{SSE}+(\bbeta-\hbeta)\T\X\T\X(\bbeta-\hbeta)}{2\sigma^2}\Bigr\}\\
&\qquad=(\sigma^2)^{-(a+p/2)-1}\exp\Bigl\{-\frac{2b+(\bbeta-\hbeta)\T\X\T\X(\bbeta-\hbeta)}{2\sigma^2}\Bigr\},
\end{aligned}
\]
with \( a=(n-p)/2>0 \) and \( b=\text{SSE}/2>0 \). This is proportional to the density @eq-opt-nig-density of
\( \text{NIG}(\hbeta,(\X\T\X)^{-1},a,b) \), which is integrable. So the posterior is proper and equal to it. By
@lem-opt-nig-marginal with \( 2a=n-p \) and \( b/a=\text{SSE}/(n-p)=s^2 \), the stated ratio has the \( t(n-p) \) distribution.
:::

The frequentist and Bayesian statements now read alike, but they say different things. The \( t \) interval of
[Section 7.5](05-sampling-distributions.html) has probability \( 1-\alpha \) of covering the fixed \( \blambda\T\bbeta \)
*before the data are seen*. The credible interval has posterior probability \( 1-\alpha \) of containing the random
\( \blambda\T\bbeta \) *given the observed data*. That the two coincide is a special property of this model and this prior.
It is one reason why \( 1/\sigma^2 \) is regarded as the natural "noninformative" prior for the normal linear model.

## A worked example

::: {#exm-opt-stackloss-bayes}
[A prior for the stack loss plant]

Return to the stack loss regression of @exm-opt-stackloss-mle. Suppose an engineer, before seeing these
\( 21 \) days, believed from the chemistry of the process that stack loss should rise roughly one for one with air
flow and with cooling water temperature, and should not depend on acid concentration. This belief is hypothetical and is
used only for illustration. We encode it as
\[
\mathbf{m}_0=(0,\ 1,\ 1,\ 0)\T,\qquad \V_0=\diag(100^2,\ 0.1^2,\ 0.2^2,\ 0.1^2),\qquad a_0=3,\ b_0=18 .
\]
So the prior standard deviations of the slopes are \( 0.1\sigma \), \( 0.2\sigma \) and \( 0.1\sigma \), the intercept is left
essentially free, and the prior mean of \( \sigma^2 \) is \( b_0/(a_0-1)=9 \).

The posterior has \( a_n=13.5 \), so the marginal posteriors are \( t \) with \( 27 \)
degrees of freedom, and \( b_n=112.523 \). The posterior means of the coefficients are
\[
\mathbf{m}_n=(-41.1490,\ 0.7616,\ 1.1606,\ {-}0.1371)\T,
\]
against the least squares estimates \( (-39.9197,\ 0.7156,\ 1.2953,\ {-}0.1521) \).
Each slope has moved toward its prior mean. The intercept has moved too, though its prior is vague, because it
is strongly correlated with the slopes when the regressors are not centred. The \( 95\% \) credible interval for the
air-flow coefficient, using the \( t(27) \) quantile \( 2.052 \), is
\( (0.550,\ 0.973) \), narrower than the least squares \( t \) interval
\( (0.431,\ 1.000) \), because the prior adds information. For water temperature and acid
concentration the intervals are \( (0.594,\ 1.727) \) and
\( (-0.392,\ 0.118) \). The posterior mean of \( \sigma^2 \) is
\( 9.002 \), almost exactly the prior mean \( 9 \), although \( s^2=10.519 \). The two parts of
@eq-opt-posterior-b explain why. The lack of fit at the posterior mean is
\( \norm{\y-\X\mathbf{m}_n}^2=180.667 \), and on its own, with \( b_0 \), it would give
\( (b_0+180.667/2)/(a_n-1)=8.667 \): the shape \( a_n=a_0+n/2 \) counts all \( n \) observations, while the residuals
have spent some of them on \( \bbeta \). The prior–data conflict term
\( (\mathbf{m}_n-\mathbf{m}_0)\T\V_0^{-1}(\mathbf{m}_n-\mathbf{m}_0)=8.378 \) adds to the scale and makes up the difference, so
the data barely move \( \sigma^2 \) from where the prior put it.

[Figure 7.6.1](#fig-opt-bayes)(a) shows the three marginal densities for the air-flow coefficient. Panel (b) traces the
posterior means of the slopes as the prior precision of the slopes is multiplied by a factor \( \tau \), from
\( 10^{-3} \) (essentially least squares) to \( 10^{3} \) (essentially the prior mean). The water-temperature coefficient
does not move monotonically. It overshoots its prior mean of \( 1 \), down to \( 0.946 \), before
returning. Matrix shrinkage acts on combinations of coefficients, and a single coordinate can pass beyond its prior
mean on the way.
:::

```{.python .run #cell-bayes-conjugate-data}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]]])
n, p = X.shape
XtX, Xty = X.T @ X, X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)
sse = np.sum((y - X @ beta_hat) ** 2)
```

```{.python .run #cell-bayes-conjugate-posterior}
def nig_posterior(m0, V0, a0, b0):
    """Conjugate update of the normal-inverse-gamma prior NIG(m0, V0, a0, b0)."""
    P0 = np.linalg.inv(V0)                         # prior precision (per unit sigma^2)
    Vn = np.linalg.inv(P0 + XtX)
    mn = Vn @ (P0 @ m0 + Xty)
    an = a0 + n / 2
    bn = b0 + (y @ y + m0 @ P0 @ m0 - mn @ (P0 + XtX) @ mn) / 2
    return mn, Vn, an, bn

# prior: vague intercept; slopes centred on 1, 1, 0 with sd 0.1, 0.2, 0.1 in units of sigma
m0 = np.array([0.0, 1.0, 1.0, 0.0])
V0 = np.diag([100.0**2, 0.1**2, 0.2**2, 0.1**2])
a0, b0 = 3.0, 18.0                                 # prior mean of sigma^2 is b0/(a0-1) = 9
mn, Vn, an, bn = nig_posterior(m0, V0, a0, b0)
scale = np.sqrt(np.diag(Vn) * bn / an)             # marginal posterior: t with 2 a_n d.f.
q = stats.t.ppf(0.975, 2 * an)
for j, name in enumerate(["intercept", "air flow", "water temp", "acid conc"]):
    print(f"{name:10s} LS {beta_hat[j]:8.4f}   posterior mean {mn[j]:8.4f}   "
          f"95% interval ({mn[j] - q * scale[j]:.3f}, {mn[j] + q * scale[j]:.3f})")
print(f"posterior mean of sigma^2 = {bn / (an - 1):.3f}")
```

The script behind this example also checks the theorem numerically. The log posterior computed from
@eq-opt-posterior-params differs from log prior plus log-likelihood by the same constant at \( 50 \) random parameter
values. The two forms of \( b_n \) agree. A Monte Carlo sample drawn from \( \sigma^2 \) and then \( \bbeta\mid\sigma^2 \)
reproduces the marginal \( t \) quantiles. And a nearly flat prior reproduces least squares and the \( t \) interval.

::: {when-format="html"}
![**Figure 7.6.1.** Conjugate analysis of the stack loss regression. (a) Marginal prior, flat-prior posterior (which
is the least squares \( t \) distribution) and posterior for the air-flow coefficient. (b) Posterior means of the three
slopes as the prior precision of the slopes is scaled by \( \tau \). Circles mark least squares, squares mark the prior
means, and the grey line marks the prior of the example (\( \tau=1 \)).](bayes_conjugate.svg){#fig-opt-bayes width=100%}
:::

::: {when-format="pdf"}
![Conjugate analysis of the stack loss regression. (a) Marginal prior, flat-prior posterior (which
is the least squares \( t \) distribution) and posterior for the air-flow coefficient. (b) Posterior means of the three
slopes as the prior precision of the slopes is scaled by \( \tau \). Circles mark least squares, squares mark the prior
means, and the grey line marks the prior of the example (\( \tau=1 \)).](bayes_conjugate.pdf){width=100%}
:::

::: {.remark}
[How the Bayesian answer relates to the optimality results]

Seen from the frequentist side, the posterior mean is a linear estimator with bias \( \W(\mathbf{m}_0-\bbeta) \) and covariance
\( \sigma^2(\I-\W)(\X\T\X)^{-1}(\I-\W)\T \), which is smaller than that of \( \hbeta \) in the nonnegative definite order.
Whether its mean squared error is smaller depends on how far \( \mathbf{m}_0 \) is from \( \bbeta \), in the metric that
\( \W \) defines. Averaged over the prior, the posterior mean is the estimator with the smallest expected squared error,
because a conditional mean minimizes mean squared error. Different criteria single out different estimators. The
Gauss–Markov theorem, the Lehmann–Scheffé theorem and the conjugate posterior (@thm-opt-bayes-conjugate) answer three different
questions, and in the flat-prior limit their answers agree.
:::

## Exercises

### A. Check your understanding

::: {#exr-opt-ig-moments}
[A1]

Show that the \( \text{IG}(a,b) \) distribution has mean \( b/(a-1) \) for \( a>1 \) and mode \( b/(a+1) \). What are the mean and
mode of \( \sigma^2\mid\y \) under the flat prior of @cor-opt-flat-prior, and how do they compare with \( s^2 \) and with the
MLE of @thm-opt-mle?
:::

::: {#exr-opt-bayes-known-variance}
[A2]

Suppose \( \sigma^2 \) is known and the prior is \( \bbeta\sim\Normal_p(\mathbf{m}_0,\sigma^2\V_0) \). Show that the posterior is
\( \Normal_p(\mathbf{m}_n,\sigma^2\V_n) \) with \( \mathbf{m}_n \) and \( \V_n \) as in @eq-opt-posterior-params.
:::

### B. Practice

::: {#exr-opt-g-prior}
[B1]

Under the \( g \)-prior of @prp-opt-shrinkage(c) with \( \mathbf{m}_0=\bzero \), show that
\( b_n=b_0+\tfrac12\bigl\{\text{SSE}+\hbeta\T\X\T\X\hbeta/(1+g)\bigr\} \). Show that as \( g\to\infty \) the posterior mean tends to \( \hbeta \),
and that for finite \( g \) the frequentist mean squared error of \( \mathbf{m}_n \) is
\( \{g^2\sigma^2\tr(\X\T\X)^{-1}+\norm{\bbeta}^2\}/(1+g)^2 \).
:::

::: {.solution}
By @prp-opt-shrinkage(b), with \( \V_0+(\X\T\X)^{-1}=(1+g)(\X\T\X)^{-1} \), the quadratic term is
\( \hbeta\T\X\T\X\hbeta/(1+g) \). As \( g\to\infty \), \( \mathbf{m}_n=g\hbeta/(1+g)\to\hbeta \). For the mean squared error,
\( \mathbf{m}_n-\bbeta=\{g(\hbeta-\bbeta)-\bbeta\}/(1+g) \), and \( \hbeta-\bbeta \) has mean zero and covariance
\( \sigma^2(\X\T\X)^{-1} \). So the expected squared norm is \( \{g^2\sigma^2\tr(\X\T\X)^{-1}+\norm{\bbeta}^2\}/(1+g)^2 \). It is
smaller than \( \sigma^2\tr(\X\T\X)^{-1} \), the value for \( \hbeta \), iff
\( \norm{\bbeta}^2<(1+2g)\sigma^2\tr(\X\T\X)^{-1} \).
:::

::: {#exr-opt-bayes-nonidentified}
[B2]

Let \( \X \) be rank deficient, \( \mathbf{m}_0=\bzero \) and \( \V_0=\tau^{-1}\I \). For \( \bv\in\Null(\X) \), show that
\( \bv\T\mathbf{m}_n=0 \) and \( \bv\T\V_n\bv=\norm{\bv}^2/\tau \). So the conditional posterior of \( \bv\T\bbeta \) given \( \sigma^2 \) equals its
prior: the data do not update nonidentified directions. What happens to \( \blambda\T\mathbf{m}_n \) for \( \blambda\in\C(\X\T) \) as
\( \tau\to0 \)?
:::

::: {.solution}
\( \V_n^{-1}\bv=(\tau\I+\X\T\X)\bv=\tau\bv \), so \( \V_n\bv=\bv/\tau \). Then
\( \bv\T\mathbf{m}_n=\bv\T\V_n\X\T\y=(\V_n\bv)\T\X\T\y=(\X\bv)\T\y/\tau=0 \), and \( \bv\T\V_n\bv=\norm{\bv}^2/\tau \). These are the prior
mean and variance (per unit \( \sigma^2 \)) of \( \bv\T\bbeta \). For the second question, \( \mathbf{m}_n=(\X\T\X+\tau\I)^{-1}\X\T\y \) converges as
\( \tau\to0 \) to the minimum-norm least squares estimate \( \X^+\y \) (@prp-proj-min-norm). This follows from the spectral
decomposition of \( \X\T\X \): eigenvectors with eigenvalue \( d>0 \) get the factor \( 1/(d+\tau)\to1/d \), and those with eigenvalue
\( 0 \) get zero. Hence \( \blambda\T\mathbf{m}_n\to\blambda\T\X^+\y=\blambda\T\hbeta \), the least squares estimate of the estimable function.
:::

### C. Going deeper

::: {#exr-opt-predictive}
[C1]

Under the posterior of @thm-opt-bayes-conjugate, let \( Y_0=\x_0\T\bbeta+\varepsilon_0 \) be a future observation, with
\( \varepsilon_0\sim\Normal(0,\sigma^2) \) independent of everything given \( \sigma^2 \). Show that the posterior predictive
distribution of \( Y_0 \) is \( t_{2a_n}\bigl(\x_0\T\mathbf{m}_n,\ (b_n/a_n)(1+\x_0\T\V_n\x_0)\bigr) \). Compare it with the frequentist
prediction interval of @exr-opt-prediction-t under the flat prior.
:::

::: {.solution}
Given \( \sigma^2 \) and \( \y \), \( \x_0\T\bbeta\sim\Normal(\x_0\T\mathbf{m}_n,\sigma^2\x_0\T\V_n\x_0) \), and \( \varepsilon_0 \) is independent of it.
So \( Y_0\mid\sigma^2,\y\sim\Normal(\x_0\T\mathbf{m}_n,\sigma^2(1+\x_0\T\V_n\x_0)) \). Thus \( (Y_0,\sigma^2)\mid\y \) is normal-inverse-gamma
with \( p=1 \), mean \( \x_0\T\mathbf{m}_n \) and \( \V=1+\x_0\T\V_n\x_0 \), and @lem-opt-nig-marginal gives the \( t_{2a_n} \) law. Under the flat
prior, \( \mathbf{m}_n=\hbeta \), \( \V_n=(\X\T\X)^{-1} \), \( 2a_n=n-p \) and \( b_n/a_n=s^2 \). The equal-tailed predictive interval is
then exactly the prediction interval of @exr-opt-prediction-t.
:::

::: {#exr-opt-sequential}
[C2]

Split the data into two batches \( (\y_1,\X_1) \) and \( (\y_2,\X_2) \). Show that updating the prior with the first batch, and
then using that posterior as the prior for the second batch, gives the same posterior as a single update with all the data.
Which quantities in @eq-opt-posterior-params make this obvious?
:::
