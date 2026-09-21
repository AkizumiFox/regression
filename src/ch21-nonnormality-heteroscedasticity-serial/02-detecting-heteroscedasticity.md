# Detecting heteroscedasticity

The errors are **heteroscedastic** when their variances differ, \( \Var(\varepsilon_i)=\sigma_i^2 \). Least squares
remains unbiased, but its covariance becomes
\( (\X\T\X)^{-1}\X\T\bSigma\X(\X\T\X)^{-1} \) with \( \bSigma=\diag(\sigma_1^2,\dots,\sigma_n^2) \), and
\( s^2(\X\T\X)^{-1} \) no longer estimates it (@thm-dep-covariance). Least squares also stops being the best
linear unbiased estimator. How far the usual standard errors go wrong depends on how the variances line up
with the regressors: the damage is greatest when the observations with the most influence on a coefficient
also have the largest variances. [Chapter 15](../ch15-anova-subspaces/index.html) showed the same thing for the
one-way layout, where unequal variances bias the \( F \) test only in unbalanced
designs (@thm-aov-heteroscedastic).

Unequal variances have familiar sources. A response that is an average of \( m_i \) observations has variance
proportional to \( 1/m_i \). Counts, sizes and amounts of money tend to vary more when they are larger. And if a
coefficient itself varies from unit to unit, \( y_i=\x_{(i)}\T(\bbeta+\bu_i) \) with random \( \bu_i \), the variance
of \( y_i \) is a quadratic function of \( \x_{(i)} \) (@exr-het-random-coefficients).

## Residual plots

Under heteroscedasticity the residuals \( \he=(\I-\M)\be \) have covariance \( (\I-\M)\bSigma(\I-\M) \), so
\[
\E\hat{\varepsilon}_i^2=(1-h_{ii})^2\sigma_i^2+\sum_{j\ne i}m_{ij}^2\sigma_j^2 .
\]{#eq-het-residual-second-moment}

When \( h_{ii} \) is small the first term dominates, and the squared residual tracks the variance of its own
observation. At high leverage it does not: a point with \( h_{ii} \) near one has a small residual whatever its
variance. Even with constant variance the residuals have unequal variances \( \sigma^2(1-h_{ii}) \), so plots
for heteroscedasticity should use the studentized residuals \( r_i=\hat{\varepsilon}_i/(s\sqrt{1-h_{ii}}) \),
which have constant variance under the model.

The standard displays are:

- \( r_i \) against the fitted values. A fan that opens to the right says that variance grows with the mean.
- \( r_i \) against each regressor, and against variables outside the model that might affect precision,
  such as time, the instrument or the interviewer.
- \( \lvert r_i\rvert \), or \( \log\lvert r_i\rvert \), against the same variables. Absolute values fold the two
  halves of the fan onto each other and make a trend in spread easier to see.

The logarithmic plot also estimates a power law. If \( \sigma_i=\sigma v_i^{\theta} \) for a positive variable
\( v_i \), then \( \varepsilon_i=\sigma v_i^\theta\zeta_i \) with \( \zeta_i \) of constant law, and
\( \log\lvert\varepsilon_i\rvert=\log\sigma+\theta\log v_i+\log\lvert\zeta_i\rvert \). So the slope of
\( \log\lvert r_i\rvert \) on \( \log v_i \) estimates \( \theta \), with the residuals standing in for the errors.

::: {#exm-het-engel-plots}
[Engel's food expenditure data]

[Section 6.9](../ch06-projections/09-inner-products.html) fitted Engel's data on the annual income and food
expenditure of \( 235 \) Belgian households, and noted that the spread of expenditure grows with income. The
least squares slope is \( 0.4852 \). [Figure 21.2.1](#fig-het-engel) shows the residual plots. Panel (a) is a
fan. Panel (b) plots \( \lvert r_i\rvert \) against income on logarithmic scales. The slope of the fitted line,
\( 0.96 \), suggests a standard deviation roughly proportional to income, the weighting used
in @exm-proj-engel. The most extreme point, at the far right of panel (a), is the household with the largest
income. Its leverage is \( 0.255 \), about thirty times the average \( p/n \), and its residual is large and
negative.
:::

::: {when-format="html"}
![**Figure 21.2.1.** Residuals from the least squares fit to Engel's data. (a) Residuals against fitted
values. (b) Absolute studentized residuals against income, both on logarithmic scales, with the least
squares line.](engel_residuals.svg){#fig-het-engel width=100%}
:::

::: {when-format="pdf"}
![Residuals from the least squares fit to Engel's data. (a) Residuals against fitted
values. (b) Absolute studentized residuals against income, both on logarithmic scales, with the least
squares line.](engel_residuals.pdf){width=100%}
:::

```{.python .run #cell-engel-hetero-fit}
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])
beta, *_ = np.linalg.lstsq(X, food, rcond=None)
fitted = X @ beta
e = food - fitted                                     # OLS residuals
Q, _ = np.linalg.qr(X)
h = np.sum(Q**2, axis=1)
s2 = e @ e / (n - 2)
# internally studentized residuals
r = e / np.sqrt(s2 * (1 - h))
```

Plots show the form of the variance function, which weighting will need. A formal test helps when the
pattern is unclear or many models must be screened.

## Limit theorems we need

Most results in the rest of the chapter are large-sample results. Three standard limit theorems, the central
limit theorem, Slutsky's lemma and the delta method, were stated in
[Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html), and
[Section 19.3](../ch19-theory-of-departures/03-non-normal-errors.html) recalled three more: the
Lindeberg–Feller central limit theorem, the Cramér–Wold device and the continuous mapping theorem. We add
two pieces of notation and one law of large numbers.

::: {.remark}
[Orders in probability and the weak law]

A sequence \( T_n \) is **bounded in probability**, written \( T_n=O_p(1) \), if for every \( \eta>0 \) there is a
\( K \) with \( \Pr(\lvert T_n\rvert>K)<\eta \) for all \( n \). If \( \E\lvert T_n\rvert \) is bounded, then
\( T_n=O_p(1) \) by Markov's inequality. We write \( T_n=o_p(1) \) if \( T_n\to0 \) in probability. Products and
sums obey \( O_p(1)\,o_p(1)=o_p(1) \) and \( O_p(1)+o_p(1)=O_p(1) \). The *weak law of large numbers* says that if
\( W_1,W_2,\dots \) are independent copies of a random variable with finite mean, then \( \bar W_n\to\E W_1 \) in
probability. Proofs are in van der Vaart (1998, chapter 2).
:::

The regression version of the central limit theorem concerns a weighted sum \( \mathbf{c}\T\boldsymbol{\upxi} \) of
independent variables, with weights determined by the design. It holds as long as no single weight
dominates. @lem-dep-weighted-clt proved it for identically distributed variables. The sandwich estimators of
[Section 21.4](04-sandwich-estimators.html) also need it for variables whose variances differ, and for that a
moment slightly higher than the second must be bounded, in the form due to Lyapunov.

::: {#lem-het-clt}
[A central limit theorem for weighted sums]

For each \( n \) let \( \xi_1,\dots,\xi_n \) be independent with mean zero, and let \( \mathbf{c}=\mathbf{c}_n\in\Real^n \) be
nonzero weights with
\[
\max_i\frac{c_i^2}{\norm{\mathbf{c}}^2}\to0 .
\]{#eq-het-no-dominant-weight}

Suppose either (a) the \( \xi_i \) have a common distribution, not depending on \( n \), with variance
\( \tau^2\in(0,\infty) \), or (b) for some
constants \( \delta,K,\tau_{\min}>0 \), \( \Var\xi_i\ge\tau_{\min}^2 \) and \( \E\lvert\xi_i\rvert^{2+\delta}\le K \) for all
\( i \) and \( n \). Then
\[
\frac{\sum_ic_i\xi_i}{\bigl(\sum_ic_i^2\Var\xi_i\bigr)^{1/2}}\to\Normal(0,1)\quad\text{in distribution.}
\]
:::

::: {.proof}
Case (a) is the lemma on weighted sums of errors (@lem-dep-weighted-clt) in
[Section 19.3](../ch19-theory-of-departures/03-non-normal-errors.html), applied with \( \xi_i \) in place of the
errors and the unit weights \( \mathbf{c}/\norm{\mathbf{c}} \). For case (b), let \( v_n=\sum_ic_i^2\Var\xi_i \) and
\( X_{ni}=c_i\xi_i/\sqrt{v_n} \), so that \( \sum_i\Var X_{ni}=1 \), and write
\( \mu_n=\max_i\lvert c_i\rvert/\norm{\mathbf{c}}\to0 \). We check Lindeberg's condition. Here
\( v_n\ge\tau_{\min}^2\norm{\mathbf{c}}^2 \), and on \( \{\lvert X_{ni}\rvert>\eta\} \) we have
\( X_{ni}^2\le\lvert X_{ni}\rvert^{2+\delta}/\eta^\delta \). Hence the Lindeberg sum satisfies
\[
\begin{aligned}
\sum_i\E\bigl[X_{ni}^2\mathbf{1}\{\lvert X_{ni}\rvert>\eta\}\bigr]
&\le\frac{K\sum_i\lvert c_i\rvert^{2+\delta}}{\eta^\delta v_n^{1+\delta/2}}\\
&\le\frac{K\mu_n^\delta\norm{\mathbf{c}}^2\norm{\mathbf{c}}^\delta}{\eta^\delta\tau_{\min}^{2+\delta}\norm{\mathbf{c}}^{2+\delta}}
=\frac{K\mu_n^\delta}{\eta^\delta\tau_{\min}^{2+\delta}}\to0 ,
\end{aligned}
\]
using \( \sum_i\lvert c_i\rvert^{2+\delta}\le\max_i\lvert c_i\rvert^\delta\sum_ic_i^2 \). The Lindeberg–Feller theorem
gives the conclusion.
:::

For the coefficient vector \( \mathbf{c}=\X(\X\T\X)^{-1}\blambda \) of a linear function \( \blambda\T\hbeta \), the
condition @eq-het-no-dominant-weight is implied by \( \max_ih_{ii}\to0 \), since \( c_i^2\le h_{ii}\norm{\mathbf{c}}^2 \), as
the proof of @thm-dep-nonnormal shows. This is the condition behind that theorem.

## The Breusch–Pagan score test

To test for heteroscedasticity we need an alternative. A **variance function** model supposes that
\[
\sigma_i^2=h(\alpha_0+\bz_i\T\boldsymbol{\upalpha}),
\]{#eq-het-variance-function}

where \( \bz_i\in\Real^q \) is a vector of observed variables, possibly some or all of the regressors, and \( h \)
is a positive function. Constant variance is the hypothesis \( \boldsymbol{\upalpha}=\bzero \). Common choices are
\( h(t)=e^t \), which keeps variances positive, and \( h(t)=t \). The score test of
[Section 11.4](../ch11-general-linear-hypothesis/04-likelihood-ratio.html) needs only the fit under the
hypothesis, which is ordinary least squares, and it turns out not to depend on \( h \) at all.

Write \( \tilde{\sigma}^2=\norm{\he}^2/n \) for the maximum likelihood estimate of \( \sigma^2 \) under the hypothesis,
\( u_i=\hat{\varepsilon}_i^2/\tilde{\sigma}^2 \) for the scaled squared residuals, and \( \tilde{\Z} \) for the
\( n\times q \) matrix whose rows are \( (\bz_i-\bar{\bz})\T \). Let \( R^2 \) be the coefficient of determination from
regressing \( \hat{\varepsilon}_i^2 \) on an intercept and \( \bz_i \), and let
\( \hat\kappa=n\sum_i\hat{\varepsilon}_i^4/(\sum_i\hat{\varepsilon}_i^2)^2 \) be the sample kurtosis of the residuals. In
this section \( \kappa \) denotes the raw kurtosis \( \mu_4/\sigma^4 \), which equals \( 3 \) for normal errors;
[Section 19.3](../ch19-theory-of-departures/03-non-normal-errors.html) used the same letter for the excess
kurtosis \( \mu_4/\sigma^4-3 \).

::: {#thm-het-breusch-pagan}
[Breusch–Pagan and Koenker]

Let \( \tilde{\Z} \) have full column rank \( q \).

::: {.enumerate options="label=(\alph*)"}
1. *(Score statistic.)* Suppose \( \varepsilon_i\sim\Normal(0,\sigma_i^2) \) independently, with
   @eq-het-variance-function, \( h \) continuously differentiable, and \( h'\ne0 \) at the null value of
   \( \alpha_0 \). The score statistic for \( \boldsymbol{\upalpha}=\bzero \) is
   \[
\text{BP}=\tfrac12(\mathbf{u}-\bone)\T\tilde{\Z}(\tilde{\Z}\T\tilde{\Z})^{-1}\tilde{\Z}\T(\mathbf{u}-\bone),
\]
   half the explained sum of squares from regressing \( u_i \) on an intercept and \( \bz_i \). It does not depend
   on \( h \).

2. *(Studentized form.)* Let \( \text{K}=nR^2 \). Then \( \text{BP}=\tfrac12(\hat\kappa-1)\,\text{K} \).

3. *(Null distributions.)* Suppose instead only that the errors are independent and identically distributed
   with mean \( 0 \), variance \( \sigma^2 \), finite fourth moment \( \mu_4 \) and kurtosis
   \( \kappa=\mu_4/\sigma^4>1 \), that \( \rank(\X)=p \) and \( q \) are fixed, that the \( \bz_i \) are bounded, and that
   \( n^{-1}\tilde{\Z}\T\tilde{\Z}\to\Q_z \), positive definite. Then, as \( n\to\infty \),
   \[
\text{K}\to\chi^2(q),\qquad \text{BP}\to\tfrac12(\kappa-1)\,\chi^2(q)\qquad\text{in distribution.}
\]
:::

:::

::: {.proof}
(a) The log-likelihood is
\( \ell=-\tfrac12\sum_i\log(2\pi\sigma_i^2)-\tfrac12\sum_i(y_i-\x_{(i)}\T\bbeta)^2/\sigma_i^2 \). Put
\( \boldsymbol{\uptheta}=(\alpha_0,\boldsymbol{\upalpha}\T)\T \), \( \bz_i^*=(1,\bz_i\T)\T \) and
\( \eta_i=\alpha_0+\bz_i\T\boldsymbol{\upalpha} \). Since \( \partial\ell/\partial\sigma_i^2=(\varepsilon_i^2/\sigma_i^2-1)/(2\sigma_i^2) \)
and \( \partial\sigma_i^2/\partial\boldsymbol{\uptheta}=h'(\eta_i)\bz_i^* \),
\[
\begin{gathered}
\frac{\partial\ell}{\partial\boldsymbol{\uptheta}}=\sum_i\frac{h'(\eta_i)}{2\sigma_i^2}\Bigl(\frac{\varepsilon_i^2}{\sigma_i^2}-1\Bigr)\bz_i^*,\\
\frac{\partial\ell}{\partial\bbeta}=\sum_i\frac{\varepsilon_i}{\sigma_i^2}\x_{(i)} .
\end{gathered}
\]
Under normality \( \varepsilon_i^2/\sigma_i^2 \) has variance \( 2 \), and \( \E\varepsilon_i(\varepsilon_i^2-\sigma_i^2)=0 \), so
the information matrix is block diagonal between \( \bbeta \) and \( \boldsymbol{\uptheta} \), and its
\( \boldsymbol{\uptheta} \) block is \( \sum_ih'(\eta_i)^2\bz_i^*\bz_i^{*\top}/(2\sigma_i^4) \). Under the hypothesis
\( \eta_i=\alpha_0 \) for all \( i \), the restricted maximum likelihood estimates are \( \hbeta \) and \( \tilde{\sigma}^2 \), and
\( h(\tilde{\alpha}_0)=\tilde{\sigma}^2 \). Write \( g=h'(\tilde{\alpha}_0) \), let \( \Z \) be the \( n\times q \) matrix with rows \( \bz_i\T \), and put
\( \Z^*=[\bone,\Z] \). At the restricted
estimate the \( \bbeta \) score vanishes, and
\[
\mathbf{S}_\theta=\frac{g}{2\tilde{\sigma}^2}\Z^{*\top}(\mathbf{u}-\bone),\qquad
\boldsymbol{\mathcal I}_{\theta\theta}=\frac{g^2}{2\tilde{\sigma}^4}\Z^{*\top}\Z^* .
\]
By block diagonality the score statistic is
\( \mathbf{S}_\theta\T\boldsymbol{\mathcal I}_{\theta\theta}^{-1}\mathbf{S}_\theta
=\tfrac12(\mathbf{u}-\bone)\T\Z^*(\Z^{*\top}\Z^*)^{-1}\Z^{*\top}(\mathbf{u}-\bone) \), in which \( g \) has cancelled. The
matrix here is the projection onto \( \C(\Z^*)=\C(\bone)\dirsum\C(\tilde{\Z}) \), and \( \bone\T(\mathbf{u}-\bone)=0 \)
because \( \sum_iu_i=n \), so only the projection onto \( \C(\tilde{\Z}) \) contributes (@lem-proj-fwl-split).

(b) Let \( \mathbf{w} \) have entries \( \hat{\varepsilon}_i^2 \). Its explained sum of squares on \( [\bone,\Z] \) is
\( \text{ESS}_w=(\mathbf{w}-\tilde{\sigma}^2\bone)\T\tilde{\M}(\mathbf{w}-\tilde{\sigma}^2\bone) \), where \( \tilde{\M} \) projects onto
\( \C(\tilde{\Z}) \), and its total sum of squares is
\( \text{TSS}_w=\sum_i(\hat{\varepsilon}_i^2-\tilde{\sigma}^2)^2=\sum_i\hat{\varepsilon}_i^4-n\tilde{\sigma}^4 \). Since
\( \mathbf{u}-\bone=(\mathbf{w}-\tilde{\sigma}^2\bone)/\tilde{\sigma}^2 \), (a) gives
\( \text{BP}=\text{ESS}_w/(2\tilde{\sigma}^4) \), while \( \text{K}=n\,\text{ESS}_w/\text{TSS}_w \). Their ratio is
\( \text{TSS}_w/(2n\tilde{\sigma}^4)=(\hat\kappa-1)/2 \).

(c) Put \( \mathbf{a}=\M\be \), so that \( \hat{\varepsilon}_i=\varepsilon_i-a_i \) and
\( \E\norm{\mathbf{a}}^2=\sigma^2p \) (@prp-het-normal-check). Let \( \lvert\tilde z_{ik}\rvert\le C \). By (b),
\[
\begin{gathered}
\text{K}=\frac{\mathbf{T}_n\T(n^{-1}\tilde{\Z}\T\tilde{\Z})^{-1}\mathbf{T}_n}{\text{TSS}_w/n},\\
\mathbf{T}_n=n^{-1/2}\tilde{\Z}\T\mathbf{w}=n^{-1/2}\tilde{\Z}\T(\mathbf{w}-\sigma^2\bone),
\end{gathered}
\]
using \( \tilde{\Z}\T\bone=\bzero \). We show three things.

*The leading term.* Let \( \mathbf{T}_n^0=n^{-1/2}\sum_i\tilde{\bz}_i(\varepsilon_i^2-\sigma^2) \). The
\( \xi_i=\varepsilon_i^2-\sigma^2 \) are independent and identically distributed with variance \( \mu_4-\sigma^4>0 \). For fixed
\( \mathbf{b}\in\Real^q \) the weights \( c_i=\mathbf{b}\T\tilde{\bz}_i \) are bounded and \( \norm{\mathbf{c}}^2/n\to\mathbf{b}\T\Q_z\mathbf{b}>0 \),
so @eq-het-no-dominant-weight holds, and @lem-het-clt(a) shows that
\( \mathbf{b}\T\mathbf{T}_n^0\to\Normal\bigl(0,(\mu_4-\sigma^4)\mathbf{b}\T\Q_z\mathbf{b}\bigr) \). By the Cramér–Wold device,
\( \mathbf{T}_n^0\to\Normal_q\bigl(\bzero,(\mu_4-\sigma^4)\Q_z\bigr) \).

*The residuals do not matter.* \( \hat{\varepsilon}_i^2-\varepsilon_i^2=a_i^2-2a_i\varepsilon_i \), so the \( k \)th coordinate of
\( \mathbf{T}_n-\mathbf{T}_n^0 \) is \( n^{-1/2}\sum_i\tilde z_{ik}a_i^2-2n^{-1/2}\be\T\B\be \), with
\( \B=\tfrac12(\bD\M+\M\bD) \) and \( \bD=\diag(\tilde z_{1k},\dots,\tilde z_{nk}) \). The first term is at most
\( Cn^{-1/2}\norm{\mathbf{a}}^2 \), whose mean \( C\sigma^2p/\sqrt n \) tends to zero. For the second,
\( \lvert\E\be\T\B\be\rvert=\sigma^2\lvert\tr(\bD\M)\rvert\le C\sigma^2p \) (@thm-rv-quadform-mean), and by
@thm-rv-quadform-variance
\[
\begin{aligned}
\Var(\be\T\B\be)&\le\lvert\mu_4-3\sigma^4\rvert\sum_ib_{ii}^2+2\sigma^4\tr(\B^2)\\
&\le(\lvert\mu_4-3\sigma^4\rvert+2\sigma^4)\,C^2p,
\end{aligned}
\]
because \( b_{ii}=\tilde z_{ik}h_{ii} \), \( \sum_ih_{ii}^2\le\sum_ih_{ii}=p \), and
\( \tr(\B^2)\le\norm{\bD\M}_F^2=\tr(\M\bD^2\M)\le C^2\tr\M \). So \( n^{-1/2}\be\T\B\be \) has mean and variance tending to
zero, and \( \mathbf{T}_n-\mathbf{T}_n^0=o_p(1) \).

*The denominator.* \( \tilde{\sigma}^2=(\norm{\be}^2-\norm{\mathbf{a}}^2)/n\to\sigma^2 \) by the weak law and
\( \norm{\mathbf{a}}^2=O_p(1) \). For the fourth powers, Minkowski's inequality in the norm
\( \norm{\mathbf{v}}_4=(\sum_iv_i^4)^{1/4} \) gives
\( \bigl\lvert\norm{\he}_4-\norm{\be}_4\bigr\rvert\le\norm{\mathbf{a}}_4\le\norm{\mathbf{a}}=O_p(1) \), while
\( n^{-1}\norm{\be}_4^4\to\mu_4 \) by the weak law. Hence \( n^{-1}\sum_i\hat{\varepsilon}_i^4\to\mu_4 \),
\( \text{TSS}_w/n\to\mu_4-\sigma^4 \) and \( \hat\kappa\to\kappa \).

Combining these by Slutsky's lemma and continuous mapping, \( \text{K} \) converges in distribution to
\( \mathbf{T}\T\Q_z^{-1}\mathbf{T}/(\mu_4-\sigma^4) \) with \( \mathbf{T}\sim\Normal_q(\bzero,(\mu_4-\sigma^4)\Q_z) \), which is
\( \chi^2(q) \) by @cor-qf-mahalanobis. The statement about BP follows from (b).
:::

Part (a) is the test of Breusch and Pagan (1979). Cook and Weisberg (1983) derived the same statistic, and
used the regression of \( u_i \) on \( \bz_i \) as a diagnostic plot. Part (c) exposes a weakness. BP is
calibrated for normal errors, which have \( \kappa=3 \). With heavier tails its null distribution is inflated
by the factor \( (\kappa-1)/2 \), and the test rejects too often, detecting kurtosis rather than
heteroscedasticity. Koenker (1981) proposed the studentized form \( nR^2 \), whose null limit is the same for
every error law with a fourth moment. It is the version to use.

::: {#exm-het-bp-size}
[Size of the two tests]

Take a straight-line regression with \( n=100 \) equally spaced values of \( x \), constant variance, and
\( \bz_i=x_i \). In \( 40{,}000 \) simulated data sets the 5% tests rejected as follows:

| error law | \( \kappa \) | Breusch–Pagan | its limit | Koenker |
|---|---|---|---|---|
| normal | 3 | 0.048 | 0.050 | 0.050 |
| Laplace | 6 | 0.185 | 0.215 | 0.047 |
| uniform | 1.8 | 0.003 | 0.002 | 0.049 |

The Laplace law makes the Breusch–Pagan test reject nearly four times too often; the uniform law makes it
almost never reject. The studentized test holds its level in all three cases.
:::

```{.python .run #cell-bp-size-size}
import numpy as np
from scipy import stats

rng = np.random.default_rng(2105)
n = 100
x = np.linspace(0, 10, n)
X = np.column_stack([np.ones(n), x])
Q, _ = np.linalg.qr(X)
# centred, unit-length variance regressor
zc = (x - x.mean()) / np.linalg.norm(x - x.mean())
crit = stats.chi2.ppf(0.95, 1)
laws = {
    "normal": (lambda size: rng.normal(size=size), 3.0),
    "Laplace": (lambda size: rng.laplace(scale=np.sqrt(0.5), size=size), 6.0),
    "uniform": (lambda size: rng.uniform(-np.sqrt(3), np.sqrt(3), size=size),
                1.8),
}


def size_table(reps):
    table = {}
    for name, (draw, kappa) in laws.items():
        E = draw((reps, n))
        R = E - (E @ Q) @ Q.T                          # residuals
        R2 = R**2
        sig2 = R2.mean(axis=1, keepdims=True)
        # explained SS of e^2 on [1, x]
        ess = (R2 @ zc) ** 2
        bp = ess / (2 * sig2[:, 0] ** 2)
        koenker = n * ess / np.sum((R2 - sig2) ** 2, axis=1)
        limit_bp = stats.chi2.sf(crit / ((kappa - 1) / 2), 1)
        table[name] = (np.mean(bp > crit), np.mean(koenker > crit), limit_bp)
        print(f"{name:8s}: Breusch-Pagan rejects {table[name][0]:.3f} "
              f"(limit {limit_bp:.3f}), "
              f"Koenker rejects {table[name][1]:.3f}")
    return table


# a quick version (the book uses 40000)
quick = size_table(2000)
```

**White's test.** For the least squares standard errors, what matters is whether \( \sigma_i^2 \) is correlated with
the entries of \( \x_{(i)}\x_{(i)}\T \), the terms of \( \X\T\bSigma\X \) (White 1980). White's test is Koenker's statistic with
\( \bz_i \) the regressors, their squares and their cross products; with many regressors it has many degrees of freedom
and low power.

::: {#exm-het-engel-tests}
[Testing Engel's data]

With \( \bz_i=\text{income}_i \), the Breusch–Pagan statistic is \( 636.0 \) and Koenker's is \( 109.3 \), on one
degree of freedom. Their ratio is \( 5.820=(\hat\kappa-1)/2 \), since the residuals have sample kurtosis
\( 12.64 \). The heavy tails of the residuals inflate BP about sixfold, but the evidence is overwhelming either
way: Koenker's p-value is \( 1.4\times 10^{-25} \). White's test, with income and its square, gives
\( 181.1 \) on two degrees of freedom.
:::

```{.python .run #cell-engel-hetero-tests}
def explained_ss(v, Z):
    """Explained sum of squares of v on [1, Z], about the mean of v."""
    Z1 = np.column_stack([np.ones(len(v)), Z])
    coef, *_ = np.linalg.lstsq(Z1, v, rcond=None)
    return np.sum((Z1 @ coef - v.mean()) ** 2)


# maximum likelihood estimate under H0
sig2 = e @ e / n
u = e**2 / sig2
# Breusch-Pagan score statistic
bp = explained_ss(u, income) / 2
tss = np.sum((e**2 - sig2) ** 2)
# n R^2 from regressing e^2 on [1, income]
koenker = n * explained_ss(e**2, income) / tss
white = n * explained_ss(e**2, np.column_stack([income, income**2])) / tss
# sample kurtosis of the residuals
kappa = n * np.sum(e**4) / np.sum(e**2) ** 2
print(f"Breusch-Pagan {bp:.1f}, Koenker {koenker:.1f}, "
      f"White {white:.1f} (2 df)")
print(f"ratio BP/Koenker = {bp / koenker:.3f}, "
      f"(kurtosis - 1)/2 = {(kappa - 1) / 2:.3f}")

# Goldfeld-Quandt: sort, drop the middle fifth
order = np.argsort(income)
m = (n - 47) // 2
low, high = order[:m], order[-m:]
sse = []
for g in (low, high):
    b_g = np.linalg.lstsq(X[g], food[g], rcond=None)[0]
    sse.append(np.sum((food[g] - X[g] @ b_g) ** 2))
gq = (sse[1] / (m - 2)) / (sse[0] / (m - 2))
print(f"Goldfeld-Quandt F = {gq:.2f} on ({m - 2}, {m - 2}) df, "
      f"p = {stats.f.sf(gq, m - 2, m - 2):.1e}")
```

## An exact test: Goldfeld–Quandt

When the variance is suspected to change monotonically with one variable \( v \) that does not depend on the
responses, there is an exact test under normality. Sort the observations by \( v \), set aside a block in the
middle, and fit the model separately to the \( n_1 \) observations with the smallest values of \( v \) and the
\( n_2 \) with the largest.

::: {#prp-het-goldfeld-quandt}
[Goldfeld–Quandt]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \), let the two groups be chosen without reference to \( \Y \), and let
\( \X_1 \) and \( \X_2 \), their rows of \( \X \), have rank \( p<n_k \). If \( \text{SSE}_1 \) and \( \text{SSE}_2 \) are the
residual sums of squares of the separate fits, then
\[
F_{\text{GQ}}=\frac{\text{SSE}_2/(n_2-p)}{\text{SSE}_1/(n_1-p)}\sim F(n_2-p,\ n_1-p).
\]
:::

::: {.proof}
Each group is a normal linear model with error variance \( \sigma^2 \), so
\( \text{SSE}_k/\sigma^2\sim\chi^2(n_k-p) \) by @thm-opt-sampling(b). The two sums of squares are functions of
disjoint sets of independent errors, so they are independent, and the ratio has the stated law by
@def-qf-noncentral-f with \( \gamma=0 \).
:::

Dropping the middle block sharpens the contrast between the groups at the cost of degrees of freedom, and
how much to drop is a matter of judgement; here we omit a fifth of the data. For Engel's data, with the middle
\( 47 \) households omitted and \( 94 \) in each group, \( F_{\text{GQ}}=10.22 \) on \( (92,92) \) degrees of freedom,
with p-value \( 1.7\times 10^{-24} \). The test is exact only under normality. Like every comparison of variances, it
is sensitive to the tails of the error law. When the observations fall into groups, a more robust
comparison is Levene's test (Levene 1960), the one-way analysis of variance of the absolute residuals
\( \lvert\hat{\varepsilon}_i\rvert \), or its version with deviations from group medians (Brown and Forsythe 1974).
Box (1953) showed how badly the classical normal-theory test of equal group variances, Bartlett's, is
affected by kurtosis.

## Exercises

### A. Check your understanding

::: {#exr-het-residual-moment}
[A1]

Derive @eq-het-residual-second-moment. Show that under constant variance \( \E\hat{\varepsilon}_i^2=\sigma^2(1-h_{ii}) \),
and give an example in which \( \sigma_i^2 \) is the largest variance but \( \E\hat{\varepsilon}_i^2 \) is the smallest.
:::

### B. Practice

::: {#exr-het-random-coefficients}
[B1]

Let \( y_i=\x_{(i)}\T(\bbeta+\bu_i)+\varepsilon_i \), where \( \bu_i \) has mean \( \bzero \) and covariance \( \boldsymbol{\Omega} \),
\( \varepsilon_i \) has variance \( \sigma^2 \), and all are independent. Show that \( \Var(y_i) \) is a linear function of
the squares and cross products of the regressors. Which test of this section is the score-type test for
\( \boldsymbol{\Omega}=\mathbf{0} \)?
:::

::: {.solution}
\( y_i-\x_{(i)}\T\bbeta=\x_{(i)}\T\bu_i+\varepsilon_i \) has variance
\( \sigma^2+\x_{(i)}\T\boldsymbol{\Omega}\x_{(i)}=\sigma^2+\sum_{j,k}\omega_{jk}x_{ij}x_{ik} \), linear in the products
\( x_{ij}x_{ik} \) with coefficients \( \omega_{jk} \). This is @eq-het-variance-function with \( h(t)=t \) and \( \bz_i \) the
distinct squares and cross products (the squares of the intercept column merge with the intercept). The score
test is the Breusch–Pagan test with these \( \bz_i \), whose studentized form is White's test.
:::

::: {#exr-het-gq-ordering}
[B2]

Under the assumptions of @prp-het-goldfeld-quandt, suppose the observations are sorted by their fitted
values \( \hat Y_i \) instead of by a fixed variable. Show that \( F_{\text{GQ}} \) still has the \( F(n_2-p,n_1-p) \)
distribution. Show by an example that sorting by \( \lvert\hat{\varepsilon}_i\rvert \) destroys the result.
:::

::: {.solution}
Fix a partition, with selection matrices \( \bS_1,\bS_2 \), and let \( \M_k \) project onto \( \C(\X_k) \). The group
residuals \( (\I-\M_k)\bS_k\be \) and the fitted values \( \M\be \) are jointly normal with cross-covariance
\( \sigma^2(\I-\M_k)\bS_k\M=\sigma^2(\I-\M_k)\X_k(\X\T\X)^{-1}\X\T=\mathbf{0} \), since \( \bS_k\X=\X_k \). So, for every fixed
partition, \( (\text{SSE}_1,\text{SSE}_2) \) is independent of \( \M\Y \), and it has the law
of @prp-het-goldfeld-quandt. Sorting by fitted values makes the partition a function of \( \M\Y \), so given
\( \M\Y \) the statistic has the \( F \) law, and therefore also unconditionally. Sorting by
\( \lvert\hat{\varepsilon}_i\rvert \) puts the largest residuals in the second group: with \( p=1 \), \( \X=\bone \) and
\( n=4 \), no middle block, the second group contains the two observations farthest from \( \bar y \), and
\( F_{\text{GQ}} \) exceeds one with probability one.
:::

### C. Going deeper

::: {#exr-het-levene}
[C1]

Suppose the observations fall into \( g \) groups and \( \bz_i \) consists of \( g-1 \) group indicators. Show that
Koenker's statistic is \( n \) times the \( R^2 \) of a one-way analysis of variance of the squared residuals, and
compare it with Levene's test.
:::

