# Zero inflation

Many count data sets have far more zeros than a Poisson or negative binomial model
with the right mean can produce. The usual explanation is that zeros come from two
places: some units are not at risk at all — they do not smoke, do not drive, are
not eligible — and contribute a zero with certainty, while the rest contribute one
only when the count happens to be zero. This section writes that story down as a
model, does the same for its rival, and asks how the two can be told apart.

## Two ways to make extra zeros

::: {#def-cnt-zero-models}
[Zero-inflated and hurdle models]

Let \( f(y;\lambda) \) be a count distribution with parameter \( \lambda \), such
as the Poisson or the NB2 family of @def-cnt-negbin, and let
\( \bz_{(i)} \) and \( \x_{(i)} \) be rows of two model matrices, which may be
equal, overlapping or disjoint.

::: {.enumerate options="label=(\alph*)"}
1. The **zero-inflated** model gives \( Y_i \) the mass function
   \[
   \Pr(Y_i=y)=
   \begin{cases}
   \pi_i+(1-\pi_i)f(0;\lambda_i), & y=0,\\[2pt]
   (1-\pi_i)f(y;\lambda_i), & y=1,2,\dots,
   \end{cases}
   \]{#eq-cnt-zip}

   with \( \operatorname{logit}\pi_i=\bz_{(i)}\T\bgamma \) and
   \( \log\lambda_i=\x_{(i)}\T\bbeta \). Equivalently, a latent indicator
   \( Z_i \) is one with probability \( \pi_i \), and \( Y_i=0 \) when
   \( Z_i=1 \) while \( Y_i\sim f(\cdot;\lambda_i) \) when \( Z_i=0 \). Taking
   \( f \) Poisson gives the **zero-inflated Poisson** (ZIP), and \( f \) negative
   binomial the **zero-inflated negative binomial** (ZINB).

2. The **hurdle** model gives \( Y_i \) the mass function
   \[
   \Pr(Y_i=y)=
   \begin{cases}
   1-\rho_i, & y=0,\\[2pt]
   \rho_i\,\dfrac{f(y;\lambda_i)}{1-f(0;\lambda_i)}, & y=1,2,\dots,
   \end{cases}
   \]{#eq-cnt-hurdle}

   with \( \operatorname{logit}\rho_i=\bz_{(i)}\T\bgamma \) and
   \( \log\lambda_i=\x_{(i)}\T\bbeta \): a binary model decides whether the count
   is positive, and a **zero-truncated** version of \( f \) decides how large it
   is.
:::

:::

These are different stories. The zero-inflated model says there are two kinds of
unit and a zero does not reveal which; the hurdle model says there is one kind of
unit and two decisions, whether to cross the threshold and how far to go once
across. Which is right is a subject-matter question, and the data distinguish them
only indirectly.

::: {#prp-cnt-zero-inflation}
[Properties of the two models]

In the notation of @def-cnt-zero-models, write \( \mu_i=\E(Y_i) \) and
\( \mathbf{1}_i=\mathbf{1}\{y_i>0\} \). The zero-inflated log-likelihood is
\[
\begin{aligned}
&\sum_{i:y_i=0}\log\bigl\{\pi_i+(1-\pi_i)f(0;\lambda_i)\bigr\}\\
&\qquad+\sum_{i:y_i>0}\bigl[\log(1-\pi_i)+\log f(y_i;\lambda_i)\bigr],
\end{aligned}
\]

and the hurdle one is
\[
\begin{aligned}
&\sum_i\bigl[\mathbf{1}_i\log\rho_i+(1-\mathbf{1}_i)\log(1-\rho_i)\bigr]\\
&\qquad+\sum_{i:y_i>0}\bigl[\log f(y_i;\lambda_i)-\log\{1-f(0;\lambda_i)\}\bigr].
\end{aligned}
\]

Then:

::: {.enumerate options="label=(\alph*)"}
1. for the zero-inflated model with \( f \) Poisson,
   \[
   \mu_i=(1-\pi_i)\lambda_i,
   \qquad
   \Var(Y_i)=\mu_i\bigl(1+\pi_i\lambda_i\bigr),
   \]{#eq-cnt-zip-moments}

   so a zero-inflated Poisson is always overdispersed relative to its own mean
   unless \( \pi_i=0 \);

2. the hurdle log-likelihood separates,
   \( \ell(\bgamma,\bbeta,\cdot)=\ell_{\text{binary}}(\bgamma)+\ell_{\text{positive}}(\bbeta,\cdot) \),
   so the two parts can be fitted independently, and their estimates are
   asymptotically independent; the zero-inflated log-likelihood does not separate;

3. the hurdle model can represent **zero deflation**, \( \rho_i>1-f(0;\lambda_i) \),
   while @eq-cnt-zip needs \( \pi_i\ge0 \) and can only add zeros; and, for the
   hurdle model, if \( \bz_{(i)} \) contains an intercept then its fitted
   probabilities of a zero average to the observed proportion of zeros, whatever
   else is in the binary part. The zero-inflated model satisfies no such identity;

4. the zero-inflated model contains \( f \) itself as the limiting case
   \( \pi_i\equiv0 \), which lies outside the range where the usual expansion is
   valid, so a likelihood ratio test of \( f \) against its zero-inflated version
   is *not* asymptotically \( \chi^{2} \); nor is the Vuong statistic of
   @eq-cnt-vuong valid for such a comparison.
:::

:::

::: {.proof}
(a) With the latent representation,
\( \E(Y_i)=\E\{\E(Y_i\mid Z_i)\}=(1-\pi_i)\lambda_i=\mu_i \). Since
\( \E(Y_i\mid Z_i)=\lambda_i(1-Z_i) \) has variance
\( \lambda_i^{2}\pi_i(1-\pi_i) \), the law of total variance gives
\[
\begin{aligned}
\Var(Y_i)&=\E\{\Var(Y_i\mid Z_i)\}+\Var\{\E(Y_i\mid Z_i)\}\\
&=(1-\pi_i)\lambda_i+\lambda_i^{2}\pi_i(1-\pi_i).
\end{aligned}
\]
Factoring out \( (1-\pi_i)\lambda_i=\mu_i \) gives
\( \Var(Y_i)=\mu_i(1+\pi_i\lambda_i) \), which exceeds \( \mu_i \) whenever
\( \pi_i>0 \).

(b) With \( \mathbf{1}_i=\mathbf{1}\{y_i>0\} \), @eq-cnt-hurdle gives
\[
\begin{aligned}
\ell&=\sum_{i}\Bigl[\mathbf{1}_i\log\rho_i+(1-\mathbf{1}_i)\log(1-\rho_i)\Bigr]\\
&\qquad+\sum_{i:\,y_i>0}\Bigl[\log f(y_i;\lambda_i)-\log\{1-f(0;\lambda_i)\}\Bigr].
\end{aligned}
\]
The first sum involves only \( \bgamma \) and the second only \( \bbeta \) and any
shape parameter of \( f \), with unconstrained parameter spaces, so the joint
maximum is the pair of separate maxima and the information is block diagonal. For
@eq-cnt-zip the contribution of a zero is
\( \log\{\pi_i+(1-\pi_i)f(0;\lambda_i)\} \), which involves both parameter
vectors inside one logarithm and does not separate.

(c) In @eq-cnt-hurdle, \( \rho_i \) is a free probability, so
\( \Pr(Y_i=0)=1-\rho_i \) can take any value in \( (0,1) \); in @eq-cnt-zip,
\( \Pr(Y_i=0)=\pi_i+(1-\pi_i)f(0;\lambda_i)\ge f(0;\lambda_i) \) for
\( \pi_i\in[0,1] \). The first sum in the display in (b) is the log-likelihood of a
logistic regression of \( \mathbf{1}_i \) on \( \bz_{(i)} \), whose likelihood
equations are \( \sum_i(\mathbf{1}_i-\hat\rho_i)\bz_{(i)}=\bzero \) by
@thm-glm-score with the canonical link. If \( \bz_{(i)} \) contains an intercept,
one of those equations reads \( \sum_i\hat\rho_i=\sum_i\mathbf{1}_i \), so
\( n^{-1}\sum_i(1-\hat\rho_i) \) is the observed proportion of zeros. The
zero-inflated likelihood equation for the inflation intercept is instead
\( \sum_i\hat\pi_i=\sum_i\hat w_i \), with \( \hat w_i \) the posterior weights
@eq-cnt-em-estep, which says nothing about \( \sum_i\widehat{\Pr}(Y_i=0) \).

(d) Setting \( \pi_i\equiv0 \) in @eq-cnt-zip recovers \( f \), but how the
\( \chi^{2} \) argument fails depends on the parameterization. If inflation is a
single probability \( \pi \), the null sits on the boundary of \( [0,1] \) and
@prp-cnt-negbin(e) applies verbatim. If
\( \operatorname{logit}\pi_i=\bz_{(i)}\T\bgamma \), the null is reached only as the
intercept tends to \( -\infty \), and under it the remaining coefficients of
\( \bgamma \) are unidentified — a harder non-regularity, where no fixed mixture
of \( \chi^{2} \) laws is the limit. Either way the calibration fails. The Vuong
statistic assumes strictly non-nested models, so that the mean of its summands is
nonzero under the null; when they are nested that mean is zero and the limit is
not normal.
:::

Part (a) matters when choosing among the four combinations. A zero-inflated
Poisson is already overdispersed, so excess zeros with no other anomaly can be
described either by a ZIP or by a negative binomial; what a ZIP cannot do is
produce a long right tail without also producing extra zeros. Zero inflation on a
negative binomial addresses zeros and tail separately, which is why ZINB is the
workhorse.

## Fitting

The hurdle model needs no new machinery: fit a logistic regression to
\( \mathbf{1}\{y_i>0\} \), then maximize
\( \sum_{y_i>0}[\log f(y_i;\lambda_i)-\log\{1-f(0;\lambda_i)\}] \) numerically over
the positive observations. The zero-inflated model does not separate, but its
latent representation invites the EM algorithm, which is what Lambert used. The E
step computes the posterior probability that observation \( i \) came from the
degenerate component,
\[
w_i=\Pr(Z_i=1\mid Y_i=y_i)=
\begin{cases}
\dfrac{\pi_i}{\pi_i+(1-\pi_i)f(0;\lambda_i)}, & y_i=0,\\[8pt]
0, & y_i>0,
\end{cases}
\]{#eq-cnt-em-estep}

by Bayes's rule. The M step maximizes the expected complete-data log-likelihood,
which *does* separate: a logistic regression of the fractional response \( w_i \)
on \( \Z \), and a count regression of \( y_i \) on \( \X \) with prior weights
\( 1-w_i \). Both are weighted fits already in hand. The M step need only
*increase* that expected log-likelihood, not maximize it, so a single weighted
least squares step for \( \bbeta \) and a single Newton step for the shape suffice;
either way the observed-data log-likelihood increases at every
iteration (@exr-cnt-em-monotone).

## An example

::: {#exm-cnt-zeros}
[Five models for the physician visits]

A proportion \( 0.3124 \) of the \( 20190 \)
person-years record no physician visit; the Poisson model of @exm-cnt-visits
predicts \( 0.0809 \).

[**Table 37.5.1.** Five count models for the physician-visit data, all with the
same eight regressors; the zero-inflated models use them in both parts. The
observed proportion of zeros is \( 0.3124 \).]{#tab-cnt-zeros}

| model | parameters | log-likelihood | AIC | BIC | fitted \( \Pr(Y=0) \) |
|---|---|---|---|---|---|
| Poisson | \( 8 \) | \( -62737.3 \) | \( 125490.7 \) | \( 125554.0 \) | \( 0.0809 \) |
| NB2 | \( 9 \) | \( -43468.5 \) | \( 86955.0 \) | \( 87026.3 \) | \( 0.3200 \) |
| ZIP | \( 16 \) | \( -54978.6 \) | \( 109989.1 \) | \( 110115.8 \) | \( 0.3125 \) |
| ZINB | \( 17 \) | \( -43391.8 \) | \( 86817.5 \) | \( 86952.0 \) | \( 0.3255 \) |
| hurdle NB | \( 17 \) | \( -43244.6 \) | \( 86523.2 \) | \( 86657.7 \) | \( 0.3124 \) |

Three things stand out. First, zero inflation alone does not save the Poisson: the
ZIP reproduces the zeros almost exactly, at \( 0.3125 \), and still
trails the plain negative binomial by more than eleven thousand in log-likelihood,
because the trouble here is the whole shape of the distribution and not only its
zeros. Second, the negative binomial *already* matches the zeros, at
\( 0.3200 \) against \( 0.3124 \), as
@exr-cnt-nb-limit predicts: heterogeneity in the rate produces excess zeros as a
by-product, and adding zero inflation gains \( 153.5 \) in twice
the log-likelihood for eight parameters while slightly *overshooting*, at
\( 0.3255 \). Third, the hurdle model wins on both AIC and BIC, and
by construction reproduces the zero proportion exactly.

The hurdle model wins because of the positive counts, not the zeros. Splitting
each log-likelihood at \( y=0 \) shows it: on the zeros the hurdle model scores
\( -6966.1 \) against the ZINB's \( -6757.5 \),
a loss of \( 208.6 \); on the positive counts
\( -36278.5 \) against \( -36634.2 \), a gain
of \( 355.7 \), for a net \( 147.2 \). Above
two visits the two are indistinguishable in
[Figure 37.5.1](#fig-cnt-zeros); the difference is in how mass is arranged over
the small positive counts, and the hurdle model buys that freedom by giving the
truncated component its own shape, \( \hat\kappa=0.557 \) against
\( 0.838 \) in the ZINB. Estimated zero inflation in the ZINB
averages only \( 0.0475 \): it too spends its extra parameters
re-shaping the low counts, which is what a hurdle model does directly.
:::

::: {when-format="html"}
![**Figure 37.5.1.** Observed distribution of physician visits (bars) and the
fitted marginal distributions \( n^{-1}\sum_i\Pr(Y_i=k) \) under four models. The
Poisson is hopeless at zero and too heavy in the middle; the other three are close
to each other and to the data.](zero_inflation.svg){#fig-cnt-zeros width=90%}
:::

::: {when-format="pdf"}
![Observed distribution of physician visits (bars) and the
fitted marginal distributions \( n^{-1}\sum_i\Pr(Y_i=k) \) under four models. The
Poisson is hopeless at zero and too heavy in the middle; the other three are close
to each other and to the data.](zero_inflation.pdf){width=90%}
:::

```{.python .run #cell-zeros-pieces}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats

data = sm.datasets.randhie.load_pandas().data
y = data["mdvis"].to_numpy(float)
names = ["lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names])
n, p = X.shape


def poisson_irls(X, y, w=None, start=None, tol=1e-11, maxit=200):
    """Weighted Poisson regression with a log link; w are prior weights."""
    w = np.ones(len(y)) if w is None else w
    beta = np.zeros(X.shape[1]) if start is None else start.copy()
    if start is None:
        beta[0] = np.log(max((w * y).sum() / w.sum(), 1e-8))
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        XW = X * (w * mu)[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


def logistic_irls(X, r, tol=1e-11, maxit=200):
    """Logistic regression with a response r in [0, 1] (a fractional response is allowed)."""
    gamma = np.zeros(X.shape[1])
    for _ in range(maxit):
        pi = 1 / (1 + np.exp(-X @ gamma))
        v = np.clip(pi * (1 - pi), 1e-10, None)
        XW = X * v[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ gamma + (r - pi) / v))
        if np.max(np.abs(step - gamma)) < tol:
            return step
        gamma = step
    return gamma


def nb_logpmf(y, mu, kappa):
    """log P(Y = y) for the negative binomial with mean mu and shape kappa."""
    return (special.gammaln(y + kappa) - special.gammaln(kappa) - special.gammaln(y + 1)
            + kappa * np.log(kappa / (kappa + mu)) + y * np.log(mu / (mu + kappa)))


def poisson_logpmf(y, mu):
    return y * np.log(mu) - mu - special.gammaln(y + 1)


def nb_step(X, y, beta, kappa, w):
    """One IRLS step in beta and one Newton step in log kappa, at prior weights w."""
    mu = np.exp(X @ beta)
    XW = X * (w * mu / (1 + mu / kappa))[:, None]        # NB2 working weights, log link
    beta = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
    mu = np.exp(X @ beta)
    s = np.sum(w * (special.digamma(y + kappa) - special.digamma(kappa)
                    + np.log(kappa / (kappa + mu)) + 1 - (kappa + y) / (kappa + mu)))
    d = np.sum(w * (special.polygamma(1, y + kappa) - special.polygamma(1, kappa)
                    + 1 / kappa - 1 / (kappa + mu) - (mu - y) / (kappa + mu) ** 2))
    return beta, float(np.clip(np.exp(np.log(kappa) - s / (kappa * d)), 1e-4, 1e4))


def nb_fit(X, y, w=None, start=None, tol=1e-10, maxit=500):
    """Weighted NB2 maximum likelihood: iterate nb_step to convergence."""
    w = np.ones(len(y)) if w is None else w
    beta, kappa = (poisson_irls(X, y, w), 1.0) if start is None else start
    for _ in range(maxit):
        new_beta, new_kappa = nb_step(X, y, beta, kappa, w)
        done = max(np.abs(new_beta - beta).max(), abs(new_kappa - kappa)) < tol
        beta, kappa = new_beta, new_kappa
        if done:
            break
    return beta, kappa


beta_pois = poisson_irls(X, y)
loglik_pois = poisson_logpmf(y, np.exp(X @ beta_pois)).sum()
beta_nb, kappa_nb = nb_fit(X, y)
loglik_nb = nb_logpmf(y, np.exp(X @ beta_nb), kappa_nb).sum()
print(f"Poisson {loglik_pois:.1f}   NB2 {loglik_nb:.1f} (kappa {kappa_nb:.3f})")
```

```{.python .run #cell-zeros-em}
def zero_inflated_em(X, Z, y, negbin=False, tol=1e-10, maxit=2000):
    """EM for a zero-inflated Poisson or NB2 model; returns (gamma, beta, kappa, loglik).

    gamma indexes the logistic model for the degenerate component, beta the log-linear
    model for the count component. The latent indicator is z_i = 1 if observation i
    comes from the degenerate component at zero.
    """
    zero = y == 0
    gamma = np.zeros(Z.shape[1])
    beta, kappa = (poisson_irls(X, y), np.inf)
    if negbin:
        beta, kappa = nb_fit(X, y)
    previous = -np.inf
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        pi = 1 / (1 + np.exp(-Z @ gamma))
        log_count = nb_logpmf(y, mu, kappa) if negbin else poisson_logpmf(y, mu)
        # E step: the posterior probability that a zero came from the degenerate part
        w = np.where(zero, pi / (pi + (1 - pi) * np.exp(log_count)), 0.0)
        # M step: a fractional logistic regression and a weighted count regression
        gamma = logistic_irls(Z, w)
        if negbin:
            beta, kappa = nb_step(X, y, beta, kappa, 1 - w)   # one step is enough
        else:
            beta = poisson_irls(X, y, 1 - w, start=beta)
        pi = 1 / (1 + np.exp(-Z @ gamma))
        mu = np.exp(X @ beta)
        log_count = nb_logpmf(y, mu, kappa) if negbin else poisson_logpmf(y, mu)
        loglik = np.sum(np.where(zero,
                                 np.log(pi + (1 - pi) * np.exp(log_count)),
                                 np.log1p(-pi) + log_count))
        if loglik - previous < tol * abs(loglik):
            break
        previous = loglik
    return gamma, beta, kappa, loglik


gamma_zip, beta_zip, _, loglik_zip = zero_inflated_em(X, X, y)
gamma_zinb, beta_zinb, kappa_zinb, loglik_zinb = zero_inflated_em(X, X, y, negbin=True)
print(f"ZIP  log-likelihood {loglik_zip:.2f}")
print(f"ZINB log-likelihood {loglik_zinb:.2f}   kappa {kappa_zinb:.3f}")
```

```{.python .run #cell-zeros-hurdle}
def truncated_nb_loglik(theta, X, y):
    """Log-likelihood of the zero-truncated NB2 model, parameters (beta, log kappa)."""
    mu = np.exp(np.clip(X @ theta[:-1], -20.0, 20.0))
    k = np.exp(theta[-1])
    log_p0 = k * np.log(k / (k + mu))                 # log P(Y = 0) before truncation
    return np.sum(nb_logpmf(y, mu, k) - np.log(-np.expm1(log_p0)))


def hurdle_nb(X, y):
    """Logistic model for y > 0, zero-truncated NB2 for the positive counts."""
    gamma = logistic_irls(X, (y > 0).astype(float))
    pi = 1 / (1 + np.exp(-X @ gamma))
    loglik_zero = np.sum(np.where(y > 0, np.log(pi), np.log1p(-pi)))

    Xp, yp = X[y > 0], y[y > 0]
    beta0, kappa0 = nb_fit(Xp, yp)                    # an untruncated fit as the start
    out = optimize.minimize(lambda t: -truncated_nb_loglik(t, Xp, yp),
                            np.append(beta0, np.log(kappa0)), method="BFGS",
                            options={"gtol": 1e-5, "maxiter": 800})
    return gamma, out.x[:-1], np.exp(out.x[-1]), loglik_zero - out.fun


gamma_h, beta_h, kappa_h, loglik_hurdle = hurdle_nb(X, y)
print(f"hurdle log-likelihood {loglik_hurdle:.2f}   kappa {kappa_h:.3f}")
```

## Telling them apart

[Table 37.5.1](05-zero-inflation.html#tab-cnt-zeros) compares fits, not
mechanisms, and settles no question about which story is right. It is worth being
clear about what the data can and cannot decide.

**Identification.** When \( \bz_{(i)} \) and \( \x_{(i)} \) are the same, the two
parts compete to explain the same zeros and the likelihood is often very flat: a
larger \( \pi_i \) with a larger \( \lambda_i \) fits almost as well as a smaller
\( \pi_i \) with a smaller one. The parameters are identified in principle, since
the positive counts pin down \( \lambda_i \), but standard errors can be large,
the two parts' coefficients strongly correlated, and the fit sensitive to the
start. Restricting \( \bz_{(i)} \) to a few variables, or to an intercept, helps
and is easier to defend.

**The Vuong test.** For two non-nested models fitted to the same data,
Vuong (1989) proposed comparing the pointwise log-likelihood ratios
\( m_i=\log f_1(y_i;\hat\theta_1)-\log f_2(y_i;\hat\theta_2) \) with
\[
V=\frac{\sqrt{n}\,\bar m}{s_m},
\qquad
s_m^{2}=\frac1{n-1}\sum_{i=1}^{n}(m_i-\bar m)^{2},
\]{#eq-cnt-vuong}

which is asymptotically standard normal when the two models are equally close to
the truth, so that large \( |V| \) favours one or the other by its sign. For the
physician visits, ZINB against the hurdle model gives
\( V=-7.23 \): taken at face value it prefers the hurdle model,
agreeing with AIC.

The test is often misused, routinely applied to compare a count model with its own
zero-inflated version — exactly the nested case the theory excludes. Under the
null \( \bar m \) and \( s_m \) both tend to zero, the ratio has no normal limit,
and the reported \( p \)-value is meaningless. In @exm-cnt-zeros that comparison
gives \( V=5.90 \), a number with no calibration; Wilson (2015)
documents the practice. A further difficulty, quantified by Schennach and Wilhelm
(2017), is poor size even in the non-nested case when the two models are close.

```{.python .run #cell-zeros-vuong}
mu_p, mu_nb = np.exp(X @ beta_pois), np.exp(X @ beta_nb)
pi_zinb = 1 / (1 + np.exp(-X @ gamma_zinb))
mu_zinb = np.exp(X @ beta_zinb)
pi_h = 1 / (1 + np.exp(-X @ gamma_h))
mu_h = np.exp(X @ beta_h)
p0_h = np.exp(kappa_h * np.log(kappa_h / (kappa_h + mu_h)))


def vuong(log_f1, log_f2):
    """Vuong's statistic for model 1 against model 2 (not valid for nested models)."""
    m = log_f1 - log_f2
    return np.sqrt(len(m)) * m.mean() / m.std(ddof=1)


log_zinb = np.where(y == 0, np.log(pi_zinb + (1 - pi_zinb) * np.exp(nb_logpmf(0.0, mu_zinb, kappa_zinb))),
                    np.log1p(-pi_zinb) + nb_logpmf(y, mu_zinb, kappa_zinb))
log_hurdle = np.where(y == 0, np.log1p(-pi_h),
                      np.log(pi_h) + nb_logpmf(y, mu_h, kappa_h) - np.log1p(-p0_h))
log_nb = nb_logpmf(y, mu_nb, kappa_nb)

print(f"ZINB against NB2 (nested, so not valid):  V = {vuong(log_zinb, log_nb):6.2f}")
print(f"ZINB against the hurdle model:            V = {vuong(log_zinb, log_hurdle):6.2f}")
```

**What does distinguish them.** By @prp-cnt-zero-inflation(c) the two differ in
what they can say about the positive counts and about zero deflation, and those
are where to look: comparing fitted and observed frequencies at \( y=1,2,3 \), as
in [Figure 37.5.1](#fig-cnt-zeros), is more informative than any single statistic.

::: {.warning}
Excess zeros are a symptom, not a diagnosis. Before reaching for a two-part model,
check that the mean model is right, that the overdispersion is not heterogeneity a
negative binomial would absorb, and that the zeros are not an artefact of the
recording process: a question never asked, a period of non-observation, a value
truncated on entry. The physician-visit data make the point. Of the two models,
the hurdle likelihood factorizes and the zero-inflated one does not, which makes
the hurdle easier to fit and to interpret; the zero-inflated model has the better
story when there really are two kinds of unit.
:::

## Exercises

### A. Check your understanding

::: {#exr-cnt-zip-mean}
[A1]

A ZIP model has \( \pi=0.3 \) and \( \lambda=4 \). Compute \( \E Y \),
\( \Var Y \), \( \Pr(Y=0) \) and the ratio of variance to mean, and compare with a
Poisson having the same mean.
:::

### B. Practice

::: {#exr-cnt-truncated}
[B1]

Let \( Y \) be Poisson with mean \( \lambda \), truncated to \( \{1,2,\dots\} \).
Show that \( \E Y=\lambda/(1-e^{-\lambda}) \) and
\( \Var Y=\E Y\{1+\lambda-\E Y\} \), and deduce that the truncated distribution is
underdispersed relative to its own mean.
:::

::: {.solution}
Write \( c=1-e^{-\lambda} \). The \( y=0 \) term contributes nothing to either
moment, so \( \E Y=\lambda/c \) and \( \E\{Y(Y-1)\}=\lambda^{2}/c \), giving
\( \Var Y=\lambda^{2}/c+\lambda/c-\lambda^{2}/c^{2}=\E Y(1+\lambda-\E Y) \). Since
\( \E Y=\lambda/c>\lambda \), the bracket is less than one.
:::

::: {#exr-cnt-zinb-vs-nb}
[B2]

Derive the mean and variance of a ZINB with inflation probability \( \pi \), count
mean \( \lambda \) and shape \( \kappa \). Show that for a given mean and variance
there is a curve of \( (\pi,\lambda,\kappa) \) values giving the same first two
moments, and explain what this implies for the precision of the estimates.
:::

::: {.solution}
\( \E Y=(1-\pi)\lambda \) and, by the law of total variance,
\( \Var Y=(1-\pi)(\lambda+\lambda^{2}/\kappa)+\pi(1-\pi)\lambda^{2}
=\E Y\{1+\lambda(\pi+1/\kappa)\} \). The mean fixes \( (1-\pi)\lambda \) and the
variance \( \lambda(\pi+1/\kappa) \): two equations in three unknowns, so a
one-dimensional family shares the first two moments. Only higher moments separate
them, so the likelihood is flat along it and the estimates are strongly correlated
with large standard errors.
:::

::: {#exr-cnt-identification}
[B3]

Fit a ZINB to the physician-visit data using only an intercept in the inflation
part, and compare its AIC with the model of @exm-cnt-zeros that uses all the
regressors in both parts. Which does AIC prefer, and what does the comparison say
about the eight inflation coefficients?
:::

::: {#exr-cnt-em-monotone}
[B4]

Show that the EM step described around @eq-cnt-em-estep increases the
observed-data log-likelihood at every iteration, by writing the observed-data
log-likelihood as the expected complete-data log-likelihood minus an entropy term
and using Jensen's inequality.
:::

### C. Going deeper

::: {#exr-cnt-hurdle-nests}
[C1]

Take \( f \) Poisson, \( \bz_{(i)}=\x_{(i)} \), and a hurdle model whose binary
part uses the complementary log–log link,
\( \log\{-\log(1-\rho_i)\}=\bz_{(i)}\T\bgamma \), in place of the logit of
@def-cnt-zero-models(b). Show that this model contains the plain Poisson as the
submodel \( \bgamma=\bbeta \), so that a likelihood ratio test of the Poisson
against it *is* a standard interior test on \( p \) restrictions, unlike the
zero-inflated comparison of @prp-cnt-zero-inflation(d). Then show that the same is
not true if the binary part uses the logit.
:::

::: {.solution}
The hurdle model reduces to \( f \) when \( \rho_i=1-f(0;\lambda_i) \) for every
\( i \). For the Poisson that is \( 1-e^{-\lambda_i} \), so
\( -\log(1-\rho_i)=\lambda_i \) and the complementary log–log link makes the
constraint read \( \bz_{(i)}\T\bgamma=\log\lambda_i=\x_{(i)}\T\bbeta \), that is
\( \bgamma=\bbeta \): a linear restriction of \( p \) dimensions in the interior
of \( \Real^{2p} \), so the usual \( \chi^{2}(p) \) limit applies. With the logit
the constraint is \( \x_{(i)}\T\bgamma=\log(e^{\lambda_i}-1) \), and since
\( t\mapsto\log(e^{e^{t}}-1) \) is not affine, the vector on the right lies in
\( \C(\X) \) only in degenerate cases such as \( \X\bbeta \) constant. The Poisson
is then not a smooth submodel at all and no likelihood ratio test of one against
the other exists: a model "containing" another is not the same as nesting it.
:::

::: {#exr-cnt-marginal-effects}
[C2]

In a zero-inflated model \( \beta_j \) is a log rate ratio for the count component,
not for the response. Derive \( \partial\log\E(Y_i)/\partial x_{ij} \) when
\( x_j \) appears in both parts, and explain why reporting \( e^{\hat\beta_j} \) as
"the effect on the expected count" is wrong.
:::

::: {.solution}
\( \E(Y_i)=(1-\pi_i)\lambda_i \), so
\( \log\E(Y_i)=\log(1-\pi_i)+\x_{(i)}\T\bbeta \) and
\[
\frac{\partial\log\E(Y_i)}{\partial x_{ij}}
=\beta_j-\frac{1}{1-\pi_i}\frac{\partial\pi_i}{\partial x_{ij}}
=\beta_j-\pi_i\gamma_j ,
\]
using \( \partial\pi_i/\partial x_{ij}=\pi_i(1-\pi_i)\gamma_j \) for the logit
link. The effect on the expected response therefore mixes the two parts and
depends on \( \pi_i \), hence on all the covariates; \( e^{\hat\beta_j} \) is the
rate ratio only within the non-degenerate component, a conditional statement about
a latent class nobody observes.
:::
