# Negative binomial models

The Poisson model spends its parameters on the mean and none on the variance.
Whether that is enough is an empirical question, and the answer is usually no.

## Where overdispersion comes from

Call a count model **overdispersed** relative to the Poisson if
\( \Var(Y_i)>\E(Y_i) \). Leave aside the wrong mean model, which @exm-cnt-grouping
ruled out first and which no dispersion parameter repairs. Two mechanisms remain:
**unobserved heterogeneity**, the units differing in ways the regressors do not
capture, which the calculation below shows always inflates the variance; and
**dependence between events**, one accident making the next more likely, so that
the count is not Poisson even conditionally. This section models the first.

::: {#lem-cnt-gamma-mixture}
[The negative binomial as a gamma mixture of Poissons]

Let \( \Lambda \) have the gamma distribution with shape \( \kappa>0 \) and rate
\( \kappa/\mu \), so that \( \E\Lambda=\mu \) and \( \Var\Lambda=\mu^{2}/\kappa \),
and let \( Y\mid\Lambda=\lambda \) be Poisson with mean \( \lambda \). Then:

::: {.enumerate options="label=(\alph*)"}
1. more generally, for any mixing distribution with
   \( \E\Lambda=\mu \) and \( \Var\Lambda=\tau^{2} \),
   \( \E Y=\mu \) and \( \Var Y=\mu+\tau^{2} \);

2. with the gamma mixing distribution above, the marginal mass function of
   \( Y \) is
   \[
   \Pr(Y=y)=\frac{\Gamma(y+\kappa)}{\Gamma(\kappa)\,y!}
   \Bigl(\frac{\kappa}{\kappa+\mu}\Bigr)^{\kappa}
   \Bigl(\frac{\mu}{\kappa+\mu}\Bigr)^{y},
   \qquad y=0,1,2,\dots
   \]{#eq-cnt-nb-pmf}

   with \( \E Y=\mu \) and \( \Var Y=\mu+\mu^{2}/\kappa \);

3. as \( \kappa\to\infty \) with \( \mu \) fixed, @eq-cnt-nb-pmf converges to the
   Poisson mass function with mean \( \mu \).
:::

:::

::: {.proof}
(a) \( \E Y=\E\Lambda=\mu \), and by the law of total variance, using
\( \Var(Y\mid\Lambda)=\Lambda \),
\( \Var Y=\E\Lambda+\Var\Lambda=\mu+\tau^{2} \).

(b) The gamma density is
\( f(\lambda)=\{(\kappa/\mu)^{\kappa}/\Gamma(\kappa)\}\lambda^{\kappa-1}e^{-\kappa\lambda/\mu} \).
Therefore
\[
\begin{aligned}
\Pr(Y=y)&=\int_0^{\infty}\frac{e^{-\lambda}\lambda^{y}}{y!}\,
\frac{(\kappa/\mu)^{\kappa}}{\Gamma(\kappa)}\lambda^{\kappa-1}
e^{-\kappa\lambda/\mu}\,d\lambda\\
&=\frac{(\kappa/\mu)^{\kappa}}{y!\,\Gamma(\kappa)}
\int_0^{\infty}\lambda^{y+\kappa-1}e^{-\lambda(1+\kappa/\mu)}\,d\lambda .
\end{aligned}
\]
The integral is \( \Gamma(y+\kappa)(1+\kappa/\mu)^{-(y+\kappa)} \), so
\[
\begin{aligned}
\Pr(Y=y)&=\frac{\Gamma(y+\kappa)}{\Gamma(\kappa)y!}
\Bigl(\frac{\kappa}{\mu}\Bigr)^{\kappa}
\Bigl(\frac{\mu}{\mu+\kappa}\Bigr)^{y+\kappa}\\
&=\frac{\Gamma(y+\kappa)}{\Gamma(\kappa)y!}
\Bigl(\frac{\kappa}{\kappa+\mu}\Bigr)^{\kappa}
\Bigl(\frac{\mu}{\kappa+\mu}\Bigr)^{y},
\end{aligned}
\]
which is @eq-cnt-nb-pmf. The moments follow from (a) with
\( \tau^{2}=\mu^{2}/\kappa \).

(c) Write \( \alpha=1/\kappa \) and use
\( \Gamma(y+\kappa)/\Gamma(\kappa)=\prod_{k=0}^{y-1}(\kappa+k)
=\kappa^{y}\prod_{k=0}^{y-1}(1+\alpha k) \). Then @eq-cnt-nb-pmf becomes
\[
\Pr(Y=y)=\frac{\mu^{y}}{y!}\Bigl\{\prod_{k=0}^{y-1}(1+\alpha k)\Bigr\}
(1+\alpha\mu)^{-1/\alpha-y} .
\]{#eq-cnt-nb-alpha}

As \( \alpha\downarrow0 \) the product tends to one and
\( (1+\alpha\mu)^{-1/\alpha}\to e^{-\mu} \), leaving \( e^{-\mu}\mu^{y}/y! \).
:::

Part (a) says the inflation is unavoidable: *any* unmodelled variation in the mean
adds to the variance, whatever its shape. The gamma is chosen in (b) for its
closed form, not because heterogeneity is really gamma distributed; the log-normal
gives no closed form, and is the natural choice once random effects enter the
picture ([Chapter 40](../ch40-glmm-gee/index.html)).

::: {#def-cnt-negbin}
[Negative binomial regression]

Let \( Y_1,\dots,Y_n \) be independent with the mass function @eq-cnt-nb-pmf, means
\( \mu_i \) and shapes \( \kappa_i \), and let \( \log\mu_i=\x_{(i)}\T\bbeta \).

::: {.enumerate options="label=(\alph*)"}
1. The **NB2** model takes \( \kappa_i=\kappa \) constant, so that
   \[
   \Var(Y_i)=\mu_i+\frac{\mu_i^{2}}{\kappa} .
   \]{#eq-cnt-nb2-variance}

   The variance is quadratic in the mean, and the coefficient of variation of
   \( \Lambda_i \) is the constant \( \kappa^{-1/2} \).

2. The **NB1** model takes \( \kappa_i=\alpha\mu_i \) proportional to the mean, so
   that
   \[
   \Var(Y_i)=\mu_i+\frac{\mu_i}{\alpha}=\phi\,\mu_i ,\qquad \phi=1+\alpha^{-1} .
   \]{#eq-cnt-nb1-variance}

   The variance is linear in the mean, with a constant inflation factor.
:::

The parameter \( \alpha=1/\kappa \) in the NB2 model is often called the
**dispersion parameter**; \( \alpha=0 \) is the Poisson.
:::

These are genuinely different models, not reparameterizations of one: NB2 says
units differ multiplicatively by a fixed relative amount, NB1 that the excess
variance is a fixed multiple of the mean.

## Estimation and inference

::: {#prp-cnt-negbin}
[Fitting the negative binomial]

In the NB2 model of @def-cnt-negbin with \( \kappa \) treated as a parameter:

::: {.enumerate options="label=(\alph*)"}
1. for fixed \( \kappa \) the distribution @eq-cnt-nb-pmf is an exponential
   dispersion family (@def-glm-edf) with natural parameter
   \( \theta=\log\{\mu/(\mu+\kappa)\} \), cumulant function
   \( b(\theta)=-\kappa\log(1-e^{\theta}) \), dispersion \( \phi=1 \) and variance
   function \( V(\mu)=\mu+\mu^{2}/\kappa \);

2. with the log link, the score for \( \bbeta \) is
   \[
   \frac{\partial\ell}{\partial\bbeta}
   =\sum_{i=1}^{n}\frac{y_i-\mu_i}{1+\mu_i/\kappa}\,\x_{(i)} ,
   \]{#eq-cnt-nb-score}

   so that Fisher scoring is iteratively reweighted least squares with weights
   \( w_i=\mu_i/(1+\mu_i/\kappa) \), and the expected information for \( \bbeta \)
   is \( \X\T\W\X \);

3. \( \bbeta \) and \( \kappa \) are orthogonal:
   \( \E\{\partial^{2}\ell/\partial\bbeta\,\partial\kappa\}=\bzero \). Hence the
   asymptotic distribution of \( \hbeta \) is the same whether \( \kappa \) is known
   or estimated, and the two can be updated alternately;

4. the score for \( \kappa \) is
   \( \sum_i\{\psi(y_i+\kappa)-\psi(\kappa)+\log\frac{\kappa}{\kappa+\mu_i}
   +1-\frac{\kappa+y_i}{\kappa+\mu_i}\} \), where \( \psi \) is the digamma
   function; it has no closed-form root and is solved numerically;

5. the Poisson model is the boundary case \( \alpha=1/\kappa=0 \). The likelihood
   ratio statistic \( L \) for \( H_0:\alpha=0 \) is therefore *not* asymptotically
   \( \chi^{2}(1) \); under the conditions of Self and Liang (1987) its limiting
   distribution is \( \tfrac12\chi^{2}(0)+\tfrac12\chi^{2}(1) \), so the
   \( p \)-value is \( \Pr\{\Normal(0,1)>\sqrt{L}\} \), half the naive one;

6. in the NB1 model of @def-cnt-negbin(b), with \( \kappa_i=\alpha\mu_i \) and a
   log link, the score for \( \bbeta \) at fixed \( \alpha \) is
   \[
   \frac{\partial\ell}{\partial\bbeta}=\alpha\sum_{i=1}^{n}\mu_i
   \Bigl\{\psi(y_i+\alpha\mu_i)-\psi(\alpha\mu_i)
   +\log\frac{\alpha}{1+\alpha}\Bigr\}\x_{(i)} .
   \]{#eq-cnt-nb1-score}

   The residual \( y_i-\mu_i \) does not appear, so this is not a weighted least
   squares equation at all. The quasi-Poisson model of [Chapter 38](../ch38-quasi-likelihood/index.html) shares
   the NB1
   variance function \( V(\mu)=\phi\mu \) but has the Poisson estimating equation
   \( \X\T(\y-\bmu)=\bzero \), so the two give different \( \hbeta \) from the same
   data (@exr-cnt-quasi-vs-nb): a variance function does not determine a fit.
   Overdispersion is the rule rather than the exception, and modelling it changes
   the coefficients little and the standard errors a great deal; which of the three
   devices to use turns on whether a likelihood is needed, as it is for AIC and for
   [Section 37.5](05-zero-inflation.html).
:::

:::

::: {.proof}
(a) Take logarithms in @eq-cnt-nb-pmf:
\[
\log\Pr(Y=y)=y\log\frac{\mu}{\mu+\kappa}+\kappa\log\frac{\kappa}{\kappa+\mu}
+\log\frac{\Gamma(y+\kappa)}{\Gamma(\kappa)y!} .
\]
With \( \theta=\log\{\mu/(\mu+\kappa)\} \), which is a strictly increasing
bijection from \( (0,\infty) \) onto \( (-\infty,0) \), we have
\( \mu=\kappa e^{\theta}/(1-e^{\theta}) \) and
\( \kappa/(\kappa+\mu)=1-e^{\theta} \), so the second term is
\( \kappa\log(1-e^{\theta})=-b(\theta) \). The last term does not involve
\( \mu \). Then
\( b'(\theta)=\kappa e^{\theta}/(1-e^{\theta})=\mu \) and
\( b''(\theta)=\kappa e^{\theta}/(1-e^{\theta})^{2}=\mu(1+\mu/\kappa) \), which by
@prp-glm-moments are the mean and the variance.

(b) By @thm-glm-score the score is
\( \sum_i\x_{(i)}(y_i-\mu_i)(\partial\mu_i/\partial\eta_i)/V(\mu_i) \). With the
log link \( \partial\mu_i/\partial\eta_i=\mu_i \), and
\( \mu_i/V(\mu_i)=\mu_i/\{\mu_i(1+\mu_i/\kappa)\}=1/(1+\mu_i/\kappa) \), giving
the score @eq-cnt-nb-score. The general IRLS weight of @thm-glm-irls is
\( (\partial\mu_i/\partial\eta_i)^{2}/V(\mu_i)=\mu_i/(1+\mu_i/\kappa) \).

(c) Differentiating @eq-cnt-nb-score, written as
\( \sum_i\x_{(i)}\kappa(y_i-\mu_i)/(\mu_i+\kappa) \), with respect to \( \kappa \)
at fixed \( \bbeta \) (so that \( \mu_i \) is fixed) gives
\[
\begin{aligned}
\frac{\partial^{2}\ell}{\partial\bbeta\,\partial\kappa}
&=\sum_i\x_{(i)}(y_i-\mu_i)\frac{\partial}{\partial\kappa}
\frac{\kappa}{\mu_i+\kappa}\\
&=\sum_i\x_{(i)}(y_i-\mu_i)\frac{\mu_i}{(\mu_i+\kappa)^{2}} .
\end{aligned}
\]
Every term is a constant times \( y_i-\mu_i \), which has mean zero, so the
expectation vanishes. The information matrix is therefore block diagonal, and the
leading block of its inverse is the inverse of the \( \bbeta \) block.

(d) Differentiate the logarithm of @eq-cnt-nb-pmf in \( \kappa \) at fixed
\( \mu \), using \( d\log\Gamma(x)/dx=\psi(x) \), and collect terms.

(f) With \( \kappa_i=\alpha\mu_i \) the two ratios in @eq-cnt-nb-pmf become
constants, \( \kappa_i/(\kappa_i+\mu_i)=\alpha/(1+\alpha) \) and
\( \mu_i/(\mu_i+\kappa_i)=1/(1+\alpha) \), so
\[
\ell=\sum_i\Bigl\{\log\frac{\Gamma(y_i+\alpha\mu_i)}{\Gamma(\alpha\mu_i)\,y_i!}
+\alpha\mu_i\log\frac{\alpha}{1+\alpha}-y_i\log(1+\alpha)\Bigr\} .
\]
Only \( \alpha\mu_i \) depends on \( \bbeta \); differentiating with
\( \partial\mu_i/\partial\bbeta=\mu_i\x_{(i)} \) gives @eq-cnt-nb1-score.

(e) By @lem-cnt-gamma-mixture(c) the Poisson is the limit \( \alpha\downarrow0 \),
and \( \alpha\ge0 \) is required for @eq-cnt-nb-pmf to be a distribution, so the
null value is on the boundary. The expansion behind the \( \chi^{2}(1) \) limit
assumes an interior maximum; here the unconstrained maximizer of the profile
likelihood falls below zero about half the time, and \( L=0 \) on that event.
@thm-mix-boundary computes the analogous mixture exactly in a balanced
variance-components problem; the general statement is Self and Liang's. The upper
tail of \( \tfrac12\chi^{2}(0)+\tfrac12\chi^{2}(1) \) at \( L \) is
\( \tfrac12\Pr\{\chi^{2}(1)>L\}=\Pr\{\Normal(0,1)>\sqrt{L}\} \).
:::

Part (c) makes the alternating algorithm honest as well as convenient: iterate the
weighted least squares step of (b) at the current \( \kappa \), then solve (d) at
the current \( \bbeta \), and the standard errors computed as though \( \kappa \)
were known are asymptotically correct.

## A score test for overdispersion

The likelihood ratio test in @prp-cnt-negbin(e) needs the negative binomial fit. A
score test needs only the Poisson fit, which makes it the natural diagnostic to
compute alongside any Poisson regression.

::: {#thm-cnt-dispersion-score}
[Score test against NB2 alternatives]

Let \( \hat\mu_i \) be the fitted means of a Poisson log-linear model with
\( \rank(\X)=p \), and consider the NB2 family @eq-cnt-nb-alpha indexed by
\( \alpha\ge0 \). Then the score for \( \alpha \) at \( \alpha=0 \) is
\[
\left.\frac{\partial\ell}{\partial\alpha}\right|_{\alpha=0}
=\tfrac12\sum_{i=1}^{n}\bigl\{(y_i-\mu_i)^{2}-y_i\bigr\},
\]{#eq-cnt-score-alpha}

its variance under the Poisson model is \( \tfrac12\sum_i\mu_i^{2} \), and its
covariance with the score for \( \bbeta \) is zero. Suppose in addition that, as
\( n\to\infty \) with \( p \) fixed, the true means stay in a fixed compact interval
inside \( (0,\infty) \) and no observation dominates, in the sense that
\[
\frac{\sum_{i=1}^{n}\mu_i^{4}}{\bigl(\sum_{i=1}^{n}\mu_i^{2}\bigr)^{2}}\longrightarrow0 .
\]{#eq-cnt-score-lyapunov}

Then
\[
S=\frac{\sum_{i=1}^{n}\bigl\{(y_i-\hat\mu_i)^{2}-y_i\bigr\}}
{\sqrt{2\sum_{i=1}^{n}\hat\mu_i^{2}}}
\]{#eq-cnt-score-test}

is asymptotically standard normal under the Poisson model, and the test rejects
for large positive \( S \).
:::

::: {.proof}
Take logarithms in @eq-cnt-nb-alpha:
\[
\ell_i(\alpha)=y_i\log\mu_i-\log y_i!+\sum_{k=0}^{y_i-1}\log(1+\alpha k)
-\Bigl(\frac1\alpha+y_i\Bigr)\log(1+\alpha\mu_i).
\]
Expand in \( \alpha \). The sum contributes
\( \alpha\sum_{k=0}^{y_i-1}k+O(\alpha^{2})=\alpha y_i(y_i-1)/2+O(\alpha^{2}) \).
For the last term, \( \log(1+\alpha\mu_i)=\alpha\mu_i-\alpha^{2}\mu_i^{2}/2+O(\alpha^{3}) \),
so
\[
\Bigl(\frac1\alpha+y_i\Bigr)\log(1+\alpha\mu_i)
=\mu_i+\alpha\Bigl(y_i\mu_i-\frac{\mu_i^{2}}{2}\Bigr)+O(\alpha^{2}).
\]
Hence
\[
\begin{aligned}
\ell_i(\alpha)&=y_i\log\mu_i-\mu_i-\log y_i!\\
&\qquad+\alpha\bigl\{y_i(y_i-1)/2-y_i\mu_i+\mu_i^{2}/2\bigr\}+O(\alpha^{2}).
\end{aligned}
\]

The leading term is the Poisson log-likelihood and the coefficient of \( \alpha \) is
\( \tfrac12\{y_i^{2}-y_i-2y_i\mu_i+\mu_i^{2}\}=\tfrac12\{(y_i-\mu_i)^{2}-y_i\} \),
which gives @eq-cnt-score-alpha.

For the variance, put \( Z_i=Y_i-\mu_i \) and
\( A_i=(Y_i-\mu_i)^{2}-Y_i=Z_i^{2}-Z_i-\mu_i \). For a Poisson variable,
\( \E Z_i=0 \), \( \E Z_i^{2}=\mu_i \), \( \E Z_i^{3}=\mu_i \) and
\( \E Z_i^{4}=\mu_i+3\mu_i^{2} \). Therefore \( \E A_i=\mu_i-0-\mu_i=0 \) and
\[
\begin{aligned}
\Var A_i&=\E\bigl(Z_i^{2}-Z_i-\mu_i\bigr)^{2}\\
&=\E Z_i^{4}-2\E Z_i^{3}-2\mu_i\E Z_i^{2}+\E Z_i^{2}+2\mu_i\E Z_i+\mu_i^{2}\\
&=2\mu_i^{2}.
\end{aligned}
\]
The \( A_i \) are independent, so the variance of @eq-cnt-score-alpha is
\( \tfrac14\sum_i2\mu_i^{2}=\tfrac12\sum_i\mu_i^{2} \).

For the covariance, differentiate @eq-cnt-score-alpha in \( \bbeta \), using
\( \partial\mu_i/\partial\bbeta=\mu_i\x_{(i)} \):
\[
\frac{\partial}{\partial\bbeta}\,
\tfrac12\bigl\{(y_i-\mu_i)^{2}-y_i\bigr\}=-\mu_i(y_i-\mu_i)\,\x_{(i)},
\]
whose expectation is zero. The information between \( \alpha \) and \( \bbeta \)
therefore vanishes at \( \alpha=0 \), so the usual adjustment for an estimated
nuisance parameter is zero and replacing \( \bbeta \) by \( \hbeta \) does not
change the asymptotic variance of the score to first order. Dividing by the
standard deviation, with \( \hat\mu_i \) for \( \mu_i \),
gives @eq-cnt-score-test. For the limit, the \( A_i \) are independent with mean zero and
\( \E A_i^{4}=O(\max(\mu_i,1)^{4}) \), a fourth-moment calculation like the one
above. Under the two hypotheses of the statement — means in a fixed compact interval
inside \( (0,\infty) \) and @eq-cnt-score-lyapunov — Lyapunov's condition holds with
\( \delta=2 \) and the standardized sum is asymptotically standard normal. Dean and
Lawless (1989) treat designs in which the means are not so confined.
:::

The statistic is one-sided by construction: a count *less* variable than the
Poisson is possible but not negative binomial, so a large negative \( S \) points
somewhere else (@exr-cnt-underdispersion).

## The physician visits again

::: {#exm-cnt-negbin}
[Negative binomial fits to the visit counts]

Refit the model of @exm-cnt-visits with the NB2 distribution, alternating the
weighted least squares step of @prp-cnt-negbin(b) with a one-dimensional solve for
\( \kappa \). The fit gives \( \hat\kappa=0.7637 \), that is
\( \hat\alpha=1.3094 \); at the overall mean of
\( 2.860 \) visits the model's variance is
\( 13.574 \) against the Poisson's \( 2.860 \).
The maximized log-likelihood rises from \( -62737.33 \) to
\( -43468.52 \), so the likelihood ratio statistic for
\( H_0:\alpha=0 \) is \( 38537.6 \); the boundary correction of
@prp-cnt-negbin(e) halves a \( p \)-value whose logarithm is about
\( -19275 \). The score statistic @eq-cnt-score-test, from the
Poisson fit alone, is \( 523.0 \).

What changes is inference. For the individual-deductible indicator the point
estimate moves only from \( -0.2648 \) to
\( -0.2475 \), but its standard error moves from
\( 0.0102 \) to \( 0.0218 \), a factor of
\( 2.125 \). Two other repairs agree: scaling the Poisson standard
errors by \( \sqrt{\hat\phi} \), with \( \hat\phi=6.319 \) the
Pearson statistic over its degrees of freedom, gives
\( 0.0257 \), and the sandwich estimator of @thm-het-sandwich gives
\( 0.0268 \). A Poisson analysis would have declared effects
significant on the strength of a variance assumption the data reject.

Fitting NB1 instead gives log-likelihood \( -43380.59 \) with
\( \hat\phi=4.766 \), a better fit than NB2 by AIC:
\( 86779.2 \) against \( 86955.0 \), with the
Poisson far behind at \( 125490.7 \).
[Figure 37.4.1](#fig-cnt-variance) shows why the comparison is close. Sorting the
person-years by fitted mean into twenty-five bins, the bin variances follow the
NB2 curve at the high end, while in the crowded low range, where most of the data
and most of the likelihood sit, NB1 is closer. AIC weighs the whole likelihood and
the figure the upper tail; when they disagree, choose by the use the model will be
put to.
:::

::: {when-format="html"}
![**Figure 37.4.1.** Sample variances of the physician-visit counts within
twenty-five bins of the fitted Poisson mean, against the three fitted variance
functions. The Poisson line is far below every bin.](variance_functions.svg){#fig-cnt-variance width=62%}
:::

::: {when-format="pdf"}
![Sample variances of the physician-visit counts within
twenty-five bins of the fitted Poisson mean, against the three fitted variance
functions. The Poisson line is far below every bin.](variance_functions.pdf){width=62%}
:::

```{.python .run #cell-negbin-poisson}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

data = sm.datasets.randhie.load_pandas().data
y = data["mdvis"].to_numpy(float)
names = ["intercept", "lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names[1:]])
n, p = X.shape

beta_pois = np.zeros(p)
beta_pois[0] = np.log(y.mean())
for _ in range(50):
    mu = np.exp(X @ beta_pois)
    XW = X * mu[:, None]
    beta_pois = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta_pois + (y - mu) / mu))
mu_pois = np.exp(X @ beta_pois)
info_pois = X.T @ (X * mu_pois[:, None])
se_pois = np.sqrt(np.diag(np.linalg.inv(info_pois)))
log_factorial = special.gammaln(y + 1)
loglik_pois = np.sum(y * np.log(mu_pois) - mu_pois - log_factorial)
print(f"Poisson log-likelihood {loglik_pois:.2f}")
```

```{.python .run #cell-negbin-nb2}
def nb_loglik(beta, kappa, X, y):
    """Negative binomial log-likelihood: mean exp(X beta), shape kappa (NB2 if scalar)."""
    mu = np.exp(X @ beta)
    return np.sum(special.gammaln(y + kappa) - special.gammaln(kappa)
                  - special.gammaln(y + 1)
                  + kappa * np.log(kappa / (kappa + mu))
                  + y * np.log(mu / (mu + kappa)))


def kappa_score(kappa, mu, y):
    """Derivative of the log-likelihood in the shape parameter."""
    return np.sum(special.digamma(y + kappa) - special.digamma(kappa)
                  + np.log(kappa / (kappa + mu)) + 1 - (kappa + y) / (kappa + mu))


def fit_nb2(X, y, tol=1e-11, maxit=200):
    """Alternate IRLS in beta with a one-dimensional search in kappa."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    kappa = 1.0
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        w = mu / (1 + mu / kappa)                  # NB2 working weights, log link
        XW = X * w[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        new_kappa = optimize.brentq(kappa_score, 1e-4, 1e4,
                                    args=(np.exp(X @ step), y), xtol=1e-12)
        done = max(np.abs(step - beta).max(), abs(new_kappa - kappa)) < tol
        beta, kappa = step, new_kappa
        if done:
            break
    return beta, kappa, nb_loglik(beta, kappa, X, y)


beta_nb2, kappa, loglik_nb2 = fit_nb2(X, y)
mu_nb2 = np.exp(X @ beta_nb2)
weights = mu_nb2 / (1 + mu_nb2 / kappa)            # IRLS weights for NB2 with a log link
se_nb2 = np.sqrt(np.diag(np.linalg.inv(X.T @ (X * weights[:, None]))))

print(f"kappa = {kappa:.4f}   log-likelihood {loglik_nb2:.2f}")
for name, b, s, sp in zip(names, beta_nb2, se_nb2, se_pois):
    print(f"{name:10s} {b:8.4f}  se {s:.4f}  (Poisson se {sp:.4f})")
```

```{.python .run #cell-negbin-tests}
lr = 2 * (loglik_nb2 - loglik_pois)
# the null puts 1/kappa on the boundary, so the null distribution of L is the mixture
# (1/2) chi2(0) + (1/2) chi2(1) and the p-value is the normal tail at sqrt(L)
log_p = stats.norm.logsf(np.sqrt(lr))

score = np.sum((y - mu_pois) ** 2 - y) / np.sqrt(2 * np.sum(mu_pois ** 2))
print(f"likelihood ratio {lr:.1f}, log p-value {log_p:.0f}")
print(f"score statistic {score:.1f}, log p-value {stats.norm.logsf(score):.0f}")
```

## Exercises

### A. Check your understanding

::: {#exr-cnt-variance-ratio}
[A1]

A negative binomial NB2 fit reports \( \hat\kappa=2 \). By what factor does the
model inflate the variance over the Poisson at a fitted mean of \( 1 \), of
\( 10 \), and of \( 100 \)? Repeat for an NB1 fit with \( \hat\phi=3 \), and
comment.
:::

::: {#exr-cnt-why-not-chisq}
[A2]

A colleague reports a likelihood ratio statistic of \( 2.9 \) for Poisson against
NB2 and, comparing with \( \chi^{2}(1) \), obtains \( p=0.089 \) and concludes
there is no overdispersion. What is the correct \( p \)-value, and does the
conclusion change at the five per cent level?
:::

### B. Practice

::: {#exr-cnt-score-hand}
[B1]

For an intercept-only Poisson model, \( \hat\mu_i=\bar y \) for every \( i \). Show
that the score statistic @eq-cnt-score-test reduces to
\( \sqrt{n/2}\,\{s^{2}_{n}/\bar y-1\} \) with
\( s^{2}_{n}=n^{-1}\sum_i(y_i-\bar y)^{2} \), so that it compares the sample
variance with the sample mean. Evaluate it for \( n=100 \), \( \bar y=3 \),
\( s^{2}_{n}=5 \).
:::

::: {.solution}
With \( \hat\mu_i=\bar y \) the numerator is
\( \sum_i(y_i-\bar y)^{2}-n\bar y=n(s^{2}_{n}-\bar y) \) and the denominator is
\( \sqrt{2n\bar y^{2}}=\bar y\sqrt{2n} \), so the ratio is
\( \sqrt{n/2}(s^{2}_{n}-\bar y)/\bar y \). For the numbers given this is
\( \sqrt{50}\times(5-3)/3 \), about \( 4.7 \): strong evidence of overdispersion.
:::

::: {#exr-cnt-nb-limit}
[B2]

Show that if \( Y \) is NB2 with mean \( \mu \) and shape \( \kappa \), then
\( \Pr(Y=0)=(1+\mu/\kappa)^{-\kappa} \), that this exceeds \( e^{-\mu} \) for every
finite \( \kappa \), and that it tends to one as \( \kappa\to0 \) with \( \mu \)
fixed. Interpret the last statement.
:::

::: {.solution}
Setting \( y=0 \) in @eq-cnt-nb-pmf leaves \( \{\kappa/(\kappa+\mu)\}^{\kappa} \),
which is \( (1+\mu/\kappa)^{-\kappa} \). Since \( \log(1+x)<x \) for \( x>0 \),
\( \kappa\log(1+\mu/\kappa)<\mu \), so the probability exceeds \( e^{-\mu} \). As
\( \kappa\to0 \), \( \kappa\log(1+\mu/\kappa)\to0 \) and the probability tends to
one: extreme heterogeneity concentrates the distribution on zero with occasional
very large counts, which is why the negative binomial can on its own account for a
great many zeros.
:::

::: {#exr-cnt-quasi-vs-nb}
[B3]

Show that the quasi-Poisson estimating equation for \( \bbeta \) is the *same* as
the Poisson one, so that \( \hbeta \) is unchanged, while the NB2 equation
@eq-cnt-nb-score downweights observations with large fitted means. Which model
would you expect to be more affected by a single very large count?
:::

::: {#exr-cnt-nb-deviance}
[B4]

Show that the deviance (@def-glm-deviance) of the NB2 model with \( \kappa \)
known is
\[
D=2\sum_{i=1}^{n}\Bigl\{y_i\log\frac{y_i}{\hat\mu_i}
-(y_i+\kappa)\log\frac{y_i+\kappa}{\hat\mu_i+\kappa}\Bigr\},
\]
and that it tends to the Poisson deviance @eq-cnt-deviance as
\( \kappa\to\infty \). Why does an analysis of deviance based on it need
\( \kappa \) to be held fixed across the models compared?
:::

::: {.solution}
Dropping the terms free of \( \mu \), the log-likelihood contribution of
observation \( i \) is
\( y_i\log\{\mu_i/(\mu_i+\kappa)\}+\kappa\log\{\kappa/(\mu_i+\kappa)\} \). The
saturated model puts \( \mu_i=y_i \), so twice the difference is
\[
2\Bigl\{y_i\log\frac{y_i}{\hat\mu_i}
+(y_i+\kappa)\log\frac{\hat\mu_i+\kappa}{y_i+\kappa}\Bigr\},
\]
which is the stated formula. For the limit write
\( (y_i+\kappa)\log\{1+(y_i-\hat\mu_i)/(\hat\mu_i+\kappa)\}
=(y_i-\hat\mu_i)+O(\kappa^{-1}) \). The deviance depends on \( \kappa \), so two
models fitted with different \( \hat\kappa \) have deviances on different scales
and their difference is not a likelihood ratio statistic; compare their maximized
log-likelihoods instead.
:::

### C. Going deeper

::: {#exr-cnt-underdispersion}
[C1]

Construct a count distribution with mean \( \mu \) and variance strictly less than
\( \mu \), and explain why no member of the negative binomial family has this
property. What mechanisms produce underdispersed counts in practice, and what
would you fit?
:::

::: {.solution}
The binomial with \( m \) trials and probability \( \pi \) has mean \( m\pi \) and
variance \( m\pi(1-\pi)<m\pi \), while by @lem-cnt-gamma-mixture(a) every Poisson
mixture has variance at least its mean. Underdispersion arises when events are
regular rather than random — inspections on a schedule, quotas — and the natural
models are the binomial, or a Conway–Maxwell–Poisson family with a dispersion
parameter free on both sides of one.
:::

::: {#exr-cnt-nb-canonical}
[C2]

@prp-cnt-negbin(a) identifies the canonical link of the NB2 family as
\( \theta=\log\{\mu/(\mu+\kappa)\} \). Work out the likelihood equations under this
link, show that they take the form \( \X\T(\y-\bmu)=\bzero \), and explain why the
log link is nevertheless the standard choice.
:::
