# Overdispersion

The word has been used loosely so far. This section says where the extra variability comes
from, how to see it, what ignoring it costs, and which repair to make.

## Where it comes from

::: {#prp-ql-overdispersion}
[Sources of overdispersion, and the price of ignoring it]

::: {.enumerate options="label=(\alph*)"}
1. *Unobserved heterogeneity.* Let \( Y\mid\Lambda \) have mean \( \Lambda \) and variance
   \( V^{*}(\Lambda) \), and let \( \Lambda \) have mean \( \mu \) and variance
   \( \tau^2 \), with \( 0<\tau^2<\infty \) and \( \E\{V^{*}(\Lambda)\}<\infty \).
   Then \( \E Y=\mu \) and \( \Var Y=\E\{V^{*}(\Lambda)\}+\tau^2 \), which for a
   Poisson kernel \( V^{*}(\lambda)=\lambda \) is \( \mu+\tau^2>\mu \): *any*
   unmodelled variation in the mean overdisperses.

2. *Clustering.* If the response is an average of \( m \) units with pairwise correlation
   \( \rho>0 \), its variance is inflated by \( 1+\rho(m-1) \), by @lem-ql-exchangeable;
   the same factor applies to a sum of \( m \) correlated indicators.

3. *Contagion.* Let \( Y=\sum_{j=1}^NS_j \), where \( N \) is Poisson with mean
   \( \nu \) and the \( S_j \) are independent copies of a positive integer-valued
   \( S \) with \( \E S^2<\infty \), independent of \( N \): events arrive in clusters.
   Then \( \E Y=\nu\,\E S \) and \( \Var Y=\nu\,\E S^2 \), so
   \[
   \frac{\Var Y}{\E Y}=\frac{\E S^2}{\E S}\ \ge\ \E S\ \ge\ 1 ,
   \]{#eq-ql-contagion}

   with equality throughout if and only if \( S\equiv1 \), that is unless the events arrive
   singly.

4. *Neither of the above.* A wrong mean model, an omitted interaction or a few gross
   outliers all inflate \( X^2 \), and no dispersion parameter repairs them; for ungrouped
   binary data a large \( X^2 \) is always a statement about the mean model.

5. *The price.* Suppose the mean model is correct and the true variances are
   \( \sigma_i^2=\phi_0\,V(\mu_i^0)/w_i \) with \( \phi_0\ge1 \) unknown, while inference
   uses \( \phi=1 \). Then under (Q1)–(Q5) any statistic of @prp-glm-three-tests for a
   \( q \)-dimensional hypothesis converges to \( \phi_0\,\chi^2(q) \), so a nominal
   level-\( \alpha \) test has asymptotic size
   \[
   \Pr\Bigl\{\chi^2(q)>\chi^2_{1-\alpha}(q)\big/\phi_0\Bigr\}\ \ge\ \alpha ,
   \]{#eq-ql-size}

   increasing in \( \phi_0 \), and tending to \( 1 \) as \( q\to\infty \) for every
   fixed \( \phi_0>1 \). Standard errors are too small by the factor
   \( \phi_0^{1/2} \).
:::

:::

::: {.proof}
(a) The law of total variance gives \( \Var Y=\E\{V^{*}(\Lambda)\}+\tau^2 \), and
\( \E\Lambda=\mu \); this is @lem-cnt-gamma-mixture(a) with the mixing distribution left
arbitrary.

(b) is @lem-ql-exchangeable.

(c) Conditioning on \( N \) gives \( \E(Y\mid N)=N\,\E S \) and
\( \Var(Y\mid N)=N\Var S \), so \( \E Y=\nu\,\E S \) and
\[
\begin{aligned}
\Var Y&=\E(N)\Var S+\Var(N)(\E S)^2\\
&=\nu\bigl\{\Var S+(\E S)^2\bigr\}=\nu\,\E S^2,
\end{aligned}
\]
using \( \Var N=\E N=\nu \). Dividing gives \( \Var Y/\E Y=\E S^2/\E S \); the first
inequality in @eq-ql-contagion is \( \E S^2\ge(\E S)^2 \), equality only for degenerate
\( S \), and the second is \( \E S\ge1 \). Equality in both forces \( S\equiv1 \), and
then \( Y=N \) is Poisson.

(d) The binary claim is the warning of
[Section 38.1](01-variance-functions.html); the rest records what \( X^2 \) measures and
is not a mathematical assertion.

(e) By @thm-ql-asymptotics(a) with \( \B_n=\phi_0\A_n \) — which is @thm-ql-score(e) at
\( \sigma_i^2=\phi_0V_i/w_i \) with the working \( \phi=1 \) —
\( \A_n^{1/2}(\hbeta-\bbeta^0)\to\Normal_p(\bzero,\phi_0\I) \). The Wald statistic
computed as though the covariance were \( \A_n^{-1} \) is
\( \norm{(\I-\bP_n)\A_n^{1/2}(\hbeta-\bbeta^0)}^2 \) for a projection \( \bP_n \) of
rank \( p-q \), as in the proof of @prp-glm-three-tests, hence converges to \( \phi_0 \)
times a \( \chi^2(q) \) variable; the other two differ from it by \( o_p(1) \) by the
same proof, since @prp-ql-dispersion(c) extends the expansion to \( Q \). The
size @eq-ql-size follows and exceeds \( \alpha \) because \( \phi_0>1 \) lowers the
threshold; the limit as \( q\to\infty \) is @exr-ql-size-in-q.
:::

These sizes are worse than intuition suggests. At \( \phi_0=2 \) a nominal \( 5\% \) test
of one coefficient rejects with probability \( 0.166 \) and of six
coefficients with probability \( 0.391 \); at \( \phi_0=5 \) the
figures are \( 0.381 \) and \( 0.866 \). A dispersion
of five is not unusual — the strike data of @exm-ql-strikes-fit have
\( \hat\phi_P=39.89 \).

```{.python .run #cell-overdispersion-tests-size}
import numpy as np
from scipy import stats

def naive_size(phi, q):
    """Size of a nominal 5% chi-squared test when the statistic is inflated by phi."""
    return stats.chi2.sf(stats.chi2.ppf(0.95, q) / phi, q)

print("phi   q = 1    q = 3    q = 6")
for phi in (1.5, 2.0, 3.0, 5.0, 10.0):
    print(f"{phi:4.1f} " + "  ".join(f"{naive_size(phi, q):7.3f}" for q in (1, 3, 6)))
```

::: {#exm-ql-size}
[A test that rejects six times too often]

Counts were simulated with a constant mean \( 6.0 \), a variance
\( 4.0 \) times the mean and a true slope of zero, and a Poisson
log-linear model was fitted to each of \( 5000 \) data sets of size
\( 100 \). Referring the drop in deviance to \( \chi^2(1) \) rejects
at the nominal \( 5\% \) level in \( 0.322 \) of the replicates —
@eq-ql-size predicts \( 0.327 \) at \( \phi_0=4 \) — with the
statistic averaging \( 3.90 \) against a \( \chi^2(1) \) mean of
one. Dividing by \( \hat\phi_P \) and referring to \( F(1,98) \) rejects in
\( 0.051 \), and the sandwich Wald test @eq-ql-robust in
\( 0.058 \)
([Figure 38.4.1](04-overdispersion.html#fig-ql-tests)).
:::

::: {when-format="html"}
![**Figure 38.4.1.** Simulated null distributions of the drop in deviance for one
coefficient, on data whose variance is four times their mean: (a) the raw statistic
against \( \chi^2(1) \); (b) the statistic divided by the Pearson dispersion,
against \( F(1,98) \). The dashed line is the upper \( 5\% \) point.](overdispersion_tests.svg){#fig-ql-tests width=100%}
:::

::: {when-format="pdf"}
![Simulated null distributions of the drop in deviance for one
coefficient, on data whose variance is four times their mean: (a) the raw statistic
against \( \chi^2(1) \); (b) the statistic divided by the Pearson dispersion,
against \( F(1,98) \). The dashed line is the upper \( 5\% \) point.](overdispersion_tests.pdf){width=100%}
:::

```{.python .run #cell-overdispersion-tests-simulation}
def poisson_fit(X, y, steps=20):
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 0.5))
    for _ in range(steps):
        mu = np.exp(X @ beta)
        z = X @ beta + (y - mu) / mu
        beta = np.linalg.solve((X.T * mu) @ X, (X.T * mu) @ z)
    return beta, np.exp(X @ beta)

def poisson_deviance(y, mu):
    return 2 * np.sum(np.where(y > 0, y * np.log(np.maximum(y, 1e-300) / mu), 0.0) - (y - mu))

def null_statistics(n, phi, mu0, reps, rng):
    """Three statistics for the slope, on data whose true slope is zero."""
    x = np.linspace(-1.0, 1.0, n)
    X = np.column_stack([np.ones(n), x])
    shape = np.full(n, mu0 / (phi - 1.0))            # negative binomial with variance phi*mu
    out = np.empty((reps, 3))
    for r in range(reps):
        y = rng.negative_binomial(shape, shape / (shape + mu0)).astype(float)
        beta, mu = poisson_fit(X, y)
        beta0, mu_null = poisson_fit(X[:, :1], y)
        drop = poisson_deviance(y, mu_null) - poisson_deviance(y, mu)
        phi_hat = np.sum((y - mu) ** 2 / mu) / (n - 2)
        bread = np.linalg.inv((X.T * mu) @ X)
        sand = bread @ ((X * ((y - mu) ** 2)[:, None]).T @ X) @ bread
        out[r] = [drop, drop / phi_hat, beta[1] ** 2 / sand[1, 1]]
    return out

N, PHI, MU0 = 100, 4.0, 6.0
rng = np.random.default_rng(3841)
stat = null_statistics(N, PHI, MU0, 400, rng)
rates = [np.mean(stat[:, 0] > stats.chi2.ppf(0.95, 1)),
         np.mean(stat[:, 1] > stats.f.ppf(0.95, 1, N - 2)),
         np.mean(stat[:, 2] > stats.chi2.ppf(0.95, 1))]
print(f"rejection rates at the 5% level: naive {rates[0]:.3f},"
      f" quasi-F {rates[1]:.3f}, sandwich Wald {rates[2]:.3f}")
```

## Seeing it

Four diagnostics, in the order in which they should be used. **Check the mean model
first**: an omitted term can remove an apparent dispersion entirely, as @exm-cnt-grouping
showed, and only a mean model that has survived @def-glm-residuals and the plots of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) makes a dispersion parameter
mean what it says. **Compare \( \hat\phi_P=X^2/(n-p) \) with one**, remembering that it
is calibrated only where @thm-glm-deviance(b) applies, noisy for ungrouped counts with
small means, and meaningless for ungrouped binary data. **Use a score test**: for counts,
@thm-cnt-dispersion-score tests the Poisson against NB2 alternatives from the Poisson fit
alone, with a one-sided alternative. **Compare the two covariances**: a large discrepancy between
\( \hat\phi_P(\X\T\hat{\W}\X)^{-1} \) and the sandwich @eq-ql-robust says the *shape*
of the variance function is wrong and not merely its scale, while a small one, as
in @exm-ql-strikes-errors, is evidence that the shape is
right (@exr-ql-information-test). A plot of squared Pearson residuals against fitted means,
and against group sizes where they vary, shows the same in more detail.

## Choosing a repair

Four devices are in use, differing in what they assume and what they give.

| device | assumes | gives | does not give |
|---|---|---|---|
| quasi model with \( \hat\phi \) | mean and variance shape | \( \hbeta \), \( F \) tests, model-based errors | likelihood, AIC, predictions |
| quasi model with the sandwich | mean model only | \( \hbeta \), Wald tests | efficiency, small samples |
| NB2, NB1 or beta-binomial | a full distribution | likelihood, AIC, predictions | robustness to that choice |
| random effect | a full distribution and a mixing law | conditional coefficients, prediction of units | the marginal coefficients |

Putting a dispersion parameter into a quasi model changes the standard errors and leaves
\( \hbeta \) exactly where it was (@thm-ql-score(a)); changing to a different likelihood
with the same variance function can move it by more than a standard
error (@exm-ql-strikes-choice). The fourth device changes the estimand rather than the
estimate: a coefficient conditional on a unit's own random effect answers a different
question from one describing a population average, a distinction [Chapter 40](../ch40-glmm-gee/index.html) takes up and
does not treat as small.

::: {#exm-ql-strikes-choice}
[Four analyses of the strike data]

Continuing @exm-ql-strikes-fit and @exm-ql-strikes-errors, the production coefficient under
four treatments of the extra variability is

| fit | coefficient | standard error |
|---|---:|---:|
| Poisson, \( \phi=1 \) | \( -7.680 \) | \( 0.406 \) |
| quasi-Poisson, \( V=\phi\mu \) | \( -7.680 \) | \( 2.564 \) |
| NB1, \( \Var=\phi\mu \) | \( -4.934 \) | \( 1.983 \) |
| NB2, \( \Var=\mu+\mu^2/\kappa \) | \( -9.327 \) | \( 2.931 \) |

The NB1 fit has \( \hat\phi=42.78 \), close to the quasi-Poisson
\( \hat\phi_P=39.89 \), and *the same variance function*; yet its
coefficient differs by more than a standard error, because its estimating
equation @eq-cnt-nb1-score is not @eq-ql-equations — the warning of
[Section 38.1](01-variance-functions.html) in numbers. The NB2 fit, with
\( \hat\alpha=1/\hat\kappa=0.954 \), gives \( -9.327 \)
against the \( -9.353 \) of the quasi model with
\( V(\mu)=\mu^2 \): at these fitted means the quadratic term of the NB2 variance
dominates, so the two nearly coincide.

The tests agree on the conclusion and not on its strength. The drop in deviance is
\( 350.9 \), overwhelming as \( \chi^2(1) \); the quasi-\( F \)
statistic @eq-ql-quasi-f is \( 8.80 \) on \( (1,60) \) degrees of
freedom, \( p=0.0043 \), and the score test for overdispersion
@eq-cnt-score-test is \( 224.5 \). Unanticipated production predicts
shorter strikes, at about the evidence a \( t \) statistic of three provides; the Poisson
analysis overstates that by orders of magnitude.
:::

```{.python .run #cell-strikes-tests}
import numpy as np
import statsmodels.api as sm
from scipy import stats

data = sm.datasets.strikes.load_pandas().data
y = data["duration"].to_numpy(float)          # length of the strike, in days
x = data["iprod"].to_numpy(float)             # unanticipated industrial production
X = np.column_stack([np.ones(len(y)), x])
n, p = X.shape

def quasi_irls(X, y, V, steps=60):
    """Solve the quasi-score equations for a log link and variance function V."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(steps):
        eta = X @ beta
        mu = np.exp(eta)
        W = mu ** 2 / V(mu)                   # working weight (dmu/deta)^2 / V(mu)
        z = eta + (y - mu) / mu               # working response
        beta = np.linalg.solve((X.T * W) @ X, (X.T * W) @ z)
    return beta

beta_lin = quasi_irls(X, y, lambda m: m)              # V(mu) = mu
mu = np.exp(X @ beta_lin)
deviance = 2 * np.sum(y * np.log(y / mu) - (y - mu))          # quasi-deviance, V(mu) = mu
pearson = np.sum((y - mu) ** 2 / mu)                          # X^2 with V(mu) = mu
phi_pearson = pearson / (n - p)
beta_null = quasi_irls(X[:, :1], y, lambda m: m)
mu_null = np.exp(X[:, :1] @ beta_null)
deviance_null = 2 * np.sum(y * np.log(y / mu_null) - (y - mu_null))
naive_chisq = deviance_null - deviance                    # the Poisson likelihood ratio
quasi_f = (deviance_null - deviance) / phi_pearson        # referred to F(1, n - p)
print(f"naive chi-squared {naive_chisq:.1f}, p = {stats.chi2.sf(naive_chisq, 1):.2e}")
print(f"quasi F {quasi_f:.2f} on (1, {n - p}) df, p = {stats.f.sf(quasi_f, 1, n - p):.4f}")
```

::: {.idea}
The quasi-score does not know about the variance at all (@thm-ql-score(b)), so a wrong
variance costs efficiency and calibration but not consistency. Whenever a colleague reports
that "adding a dispersion parameter did not change the results", the coefficients were
being looked at and the standard errors were not.
:::

## Exercises

### A. Check your understanding

::: {#exr-ql-size-in-q}
[A1]

Show that for fixed \( \phi_0>1 \) the size @eq-ql-size tends to \( 1 \) as
\( q\to\infty \). (Use \( \chi^2(q)/q\to1 \) in probability, together with
\( \chi^2_{1-\alpha}(q)/q\to1 \), and compare with the event
\( \chi^2(q)/q>\chi^2_{1-\alpha}(q)/(q\phi_0) \).)
:::

::: {.solution}
By the law of large numbers \( \chi^2(q)/q\to1 \) in probability, and the same
concentration gives \( \chi^2_{1-\alpha}(q)/q=c_q\to1 \). The event in @eq-ql-size is
\( \{\chi^2(q)/q>c_q/\phi_0\} \), whose threshold tends to \( \phi_0^{-1}<1 \), so its
probability tends to \( 1 \); at \( q=1 \) it is the finite
\( \Pr\{\chi^2(1)>3.841/\phi_0\} \). The listing of this section shows the size
increasing with \( q \) over the range tabulated, but that is a computation, not a
proof.
:::

::: {#exr-ql-mean-model-first}
[A2]

A Poisson regression of accident counts on three covariates gives \( X^2/(n-p)=3.2 \). Give
three distinct explanations, and for each one say what you would compute next.
:::

### B. Practice

::: {#exr-ql-compound-numbers}
[B1]

In @prp-ql-overdispersion(c), let \( S \) be geometric on \( \{1,2,\dots\} \) with mean
\( 1/\theta \). Compute \( \Var Y/\E Y \) and show that it equals \( (2-\theta)/\theta \).
For which \( \theta \) does the compound Poisson have the variance-to-mean ratio of a
negative binomial with \( \kappa=\mu \)?
:::

::: {#exr-ql-binomial-detect}
[B2]

Grouped binomial data have group sizes \( m_i \) spread between \( 2 \) and \( 50 \). Show
that under the constant-dispersion specification the Pearson terms
\( m_i(y_i-\hat\mu_i)^2/\{\hat\mu_i(1-\hat\mu_i)\} \) have expectation \( \phi \) free of
\( m_i \), while under @eq-ql-exchangeable they have expectation \( 1+\rho(m_i-1) \).
Describe the plot this suggests and what each specification predicts for it.
:::

::: {.solution}
Under \( \Var(Y_i)=\phi\mu_i(1-\mu_i)/m_i \) the Pearson term at the true mean has
expectation \( \phi \); under @eq-ql-exchangeable the same calculation gives
\( 1+\rho(m_i-1) \), linear in \( m_i \) with slope \( \rho \). So plot squared Pearson
residuals against \( m_i \) with a smooth: flat supports a constant dispersion,
increasing the exchangeable form.
:::

::: {#exr-ql-nb1-vs-quasi}
[B3]

Both the NB1 model and the quasi-Poisson specification have \( \Var(Y_i)=\phi\mu_i \).
Using @eq-cnt-nb1-score and @eq-ql-equations, explain why a single very large count
influences the two fits differently, and predict which of the two coefficients
in @exm-ql-strikes-choice would move more if the longest strike were deleted.
:::

::: {.solution}
In the quasi-Poisson equation the log link cancels \( \mu_i \) against \( V(\mu_i) \), so
observation \( i \) enters through the raw residual \( y_i-\mu_i \): the longest strike,
a count of \( 216 \) against a fitted mean of
\( 58.99 \), contributes \( 157.0 \). The NB1
score @eq-cnt-nb1-score contains \( y_i \) only inside \( \psi(y_i+\alpha\mu_i) \), which
grows like \( \log y_i \), so the NB1 fit is far less sensitive to it and the
quasi-Poisson coefficient should move more when it is removed.
:::

### C. Going deeper

::: {#exr-ql-lognormal-shape}
[C1]

@prp-ql-overdispersion(a) says that unmodelled heterogeneity overdisperses. It says nothing about the
*shape* of the excess, and the shape decides the repair. Suppose that, given an unrecorded
\( Z_i\sim\Normal(0,1) \) independent of the covariate, \( Y_i \) is Poisson with mean
\( \Lambda_i=\exp(\beta_0+\beta_1x_i+\gamma Z_i) \), \( \gamma\ne0 \).

::: {.enumerate options="label=(\alph*)"}
1. Using @thm-mvn-mgf, show that \( \mu_i=\E(Y_i)=\exp(\beta_0+\gamma^2/2+\beta_1x_i) \): the marginal mean
   is still log-linear, with the same slope and a raised intercept. So the mean model is *not* at fault, and
   this is not case (d) of the proposition.

2. Show that \( \Var(Y_i)=\mu_i+(e^{\gamma^2}-1)\mu_i^2 \), so that
   \[
   \frac{\Var(Y_i)}{\E(Y_i)}=1+(e^{\gamma^2}-1)\mu_i ,
   \]
   which increases with the mean instead of being constant.

3. Which of the four devices in the table of this section is the right one, and what does the plot of
   squared Pearson residuals against \( \hat\mu_i \) look like under the true model and under the
   quasi-Poisson specification \( V(\mu)=\phi\mu \)?

4. Show that \( \hat\beta_1 \) from the quasi-Poisson fit is still consistent, and explain why no single
   \( \hat\phi \) makes the reported standard errors right across a design whose fitted means vary widely.
:::

:::

::: {.solution}
(a) By @thm-mvn-mgf,
\( \E(Y_i)=\E\Lambda_i=e^{\beta_0+\beta_1x_i}\E e^{\gamma Z_i}=e^{\beta_0+\beta_1x_i}e^{\gamma^2/2} \).
The slope of \( \log\mu_i \) in \( x_i \) is \( \beta_1 \), unchanged; only the intercept
moves, by \( \gamma^2/2 \).

(b) By the proposition with the Poisson kernel, \( \Var(Y_i)=\E\Lambda_i+\Var(\Lambda_i) \). Now
\( \E\Lambda_i^2=e^{2(\beta_0+\beta_1x_i)}\E e^{2\gamma Z_i}=e^{2(\beta_0+\beta_1x_i)}e^{2\gamma^2} \), while
\( \mu_i^2=e^{2(\beta_0+\beta_1x_i)}e^{\gamma^2} \), so
\( \Var(\Lambda_i)=\mu_i^2(e^{\gamma^2}-1) \). Add \( \mu_i \) and divide.

(c) The variance is \( \mu+\mu^2/\kappa \) with \( 1/\kappa=e^{\gamma^2}-1 \), which is exactly the NB2
variance function: the negative binomial with a lognormal mixing law replaced by a gamma one has the same
first two moments, and the third device in the table — a full distribution with that variance function — is
the right choice, as is a quasi model with \( V(\mu)=\mu+\mu^2/\kappa \). Under the true model the squared
Pearson residuals \( (y_i-\hat\mu_i)^2/\hat\mu_i \) have expectation \( 1+(e^{\gamma^2}-1)\hat\mu_i \), so a
smooth through them rises linearly in \( \hat\mu_i \); under the quasi-Poisson specification it should be
flat at \( \hat\phi \). This is the diagnostic of the "Seeing it" list, and it is also what makes the
sandwich covariance @eq-ql-robust and the model-based \( \hat\phi_P(\X\T\hat{\W}\X)^{-1} \) disagree, since
the shape and not merely the scale of \( V \) is wrong.

(d) The quasi-score with a log link and \( V(\mu)=\mu \) is \( \sum_i\x_{(i)}(y_i-\mu_i) \), whose mean is
zero at the true \( \bbeta \) because (a) says the mean model holds; consistency then follows from
@thm-ql-score(b) and @thm-ql-asymptotics(a), which never use the variance. The Pearson dispersion converges
to a weighted average of the ratios \( 1+(e^{\gamma^2}-1)\mu_i \), a number that is too small where the
fitted means are large and too large where they are small; the intervals for observations at the ends of
the design are wrong in opposite directions, and rescaling every variance by one constant cannot mend that.
The idea box of this section cuts the other way here: the coefficients are safe and the standard errors are
not, but the cure is a new variance function rather than a new number.
:::

::: {#exr-ql-no-scalar-dispersion}
[C2]

@prp-ql-overdispersion(b) attributes the inflation \( 1+\rho(m-1) \) to clustering. When the individual
units, and not their averages, are the data, a scalar dispersion parameter cannot represent it at all. Let
the data be \( Y_{ij} \), \( j=1,\dots,n_i \), \( i=1,\dots,m \), with means \( \mu_{ij} \), marginal
variance \( \sigma^2 \) throughout, exchangeable within-cluster correlation
\( \rho\in(0,1) \) (@eq-ql-exchangeable) and independence across clusters, and let
\( S=\sum_{ij}a_{ij}Y_{ij} \) for constants \( a_{ij} \).

::: {.enumerate options="label=(\alph*)"}
1. Show that
   \( \Var(S)=\sigma^2\sum_i\bigl\{(1-\rho)\sum_ja_{ij}^2+\rho\bigl(\sum_ja_{ij}\bigr)^2\bigr\} \).

2. Deduce that a *cluster-level* contrast, \( a_{ij}=c_i \) for all \( j \), has
   \( \Var(S)=\sigma^2\sum_ic_i^2n_i\{1+\rho(n_i-1)\} \): inflated by the design effect.

3. Deduce that a *within-cluster* contrast, \( \sum_ja_{ij}=0 \) for every \( i \), has
   \( \Var(S)=\sigma^2(1-\rho)\sum_{ij}a_{ij}^2 \): **deflated**, by the factor \( 1-\rho<1 \).

4. Show that the Pearson dispersion \( \hat\phi_P \) converges to \( 1 \) here, so it registers no
   overdispersion at all, and that no single \( \phi \) can correct both kinds of contrast even if one were
   supplied by hand. Which device in the table is then the right one?
:::

:::

::: {.solution}
(a) Clusters are independent, so \( \Var(S)=\sum_i\mathbf{a}_i\T\Cov(\Y_i)\mathbf{a}_i \) with
\( \Cov(\Y_i)=\sigma^2\{(1-\rho)\I+\rho\bone\bone\T\} \), and
\( \mathbf{a}\T\{(1-\rho)\I+\rho\bone\bone\T\}\mathbf{a}=(1-\rho)\norm{\mathbf{a}}^2+\rho(\bone\T\mathbf{a})^2 \).

(b) With \( a_{ij}=c_i \), \( \sum_ja_{ij}^2=c_i^2n_i \) and \( (\sum_ja_{ij})^2=c_i^2n_i^2 \), so the
bracket is \( c_i^2\{(1-\rho)n_i+\rho n_i^2\}=c_i^2n_i\{1+\rho(n_i-1)\} \).

(c) The second term vanishes for every \( i \), leaving \( \sigma^2(1-\rho)\sum_{ij}a_{ij}^2 \), which is
smaller than the independence value \( \sigma^2\sum_{ij}a_{ij}^2 \).

(d) The Pearson statistic uses only the marginal variances, which are correct by assumption:
\( \hat\phi_P=\sum_{ij}(y_{ij}-\hat\mu_{ij})^2/\{\sigma^2(N-p)\}\to1 \) by the law of large numbers over
clusters, whatever \( \rho \) is. Correlation is invisible to it, because it is a statement about pairs and
the Pearson statistic looks at one observation at a time. Even a \( \phi \) chosen with knowledge of
\( \rho \) fails: the cluster-level contrast needs the factor \( 1+\rho(n_i-1)>1 \) and the within-cluster
contrast the factor \( 1-\rho<1 \), and a single multiplier cannot be on both sides of one. (With unequal
\( n_i \) even the cluster-level factors differ among themselves.) So the right device is the second or the
fourth: a sandwich covariance with the *cluster* as the unit, which estimates \( \Var(S) \) directly from
the between-cluster variation, or an explicit model for the within-cluster covariance. The choice, and the
different estimands the two carry, is the business of
[Chapter 40](../ch40-glmm-gee/index.html).
:::
