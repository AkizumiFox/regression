# Sandwich covariance estimation

@thm-ql-score left an asymmetry: the estimating equations are unbiased as soon as the mean
model is right, while the identity that turns \( \A_n=\phi^{-1}\X\T\W\X \) into a
covariance needs the variance model too. When the latter fails, the covariance of
\( \hbeta \) becomes a product of three matrices instead of one, and the middle factor can
be estimated from the residuals even though its \( n \) ingredients cannot.

## The asymptotic distribution

Throughout \( \bbeta^0 \) is the true parameter, \( \eta_i^0=\x_{(i)}\T\bbeta^0 \),
\( \mu_i^0=h(\eta_i^0) \), and \( \sigma_i^2=\Var(Y_i) \) is *not* assumed equal to
\( \phi V(\mu_i^0)/w_i \). Write
\[
a(\eta)=\frac{h'(\eta)}{V\{h(\eta)\}},\qquad
\A_n=\frac1\phi\,\X\T\W\X=\frac1\phi\sum_{i=1}^nw_ia(\eta_i^0)h'(\eta_i^0)\,\x_{(i)}\x_{(i)}\T ,
\]
and let \( \B_n \) be the matrix @eq-ql-meat. The conditions are

- **(Q1)** the \( Y_i \) are independent with \( \E Y_i=\mu_i^0 \) for a \( \bbeta^0 \)
  interior to the admissible set, \( p \) is fixed and \( \rank(\X)=p \);
- **(Q2)** \( \norm{\x_{(i)}}\le C \) and \( 0<w_{\min}\le w_i\le w_{\max}<\infty \);
- **(Q3)** on a neighbourhood \( \mathcal N \) of \( \bbeta^0 \), \( h \) and \( a \) are
  continuously differentiable with \( h' \), \( a \), \( a' \) bounded and \( h' \),
  \( a \) bounded away from zero, uniformly in \( i \) and \( n \);
- **(Q4)** \( \sigma_i^2\ge\sigma_{\min}^2>0 \), \( \E(Y_i-\mu_i^0)^4\le K \), and
  \( n^{-1}\A_n\to\bar{\A} \), \( n^{-1}\B_n\to\bar{\B} \) with both limits positive
  definite;
- **(Q5)** with probability tending to one @eq-ql-equations has a root
  \( \hbeta\in\mathcal N \) with \( \hbeta\to\bbeta^0 \) in probability.

::: {#thm-ql-asymptotics}
[Consistency, normality and the sandwich]

Assume (Q1)–(Q5). Then:

::: {.enumerate options="label=(\alph*)"}
1. \[
   \B_n^{-1/2}\A_n\bigl(\hbeta-\bbeta^0\bigr)\ \xrightarrow{d}\ \Normal_p(\bzero,\I),
   \]{#eq-ql-sandwich-limit}

   so that \( \hbeta \) behaves as if its covariance were the **sandwich**
   \( \A_n^{-1}\B_n\A_n^{-1} \), which involves the distribution of \( \Y \) only
   through its first two moments.

2. If the variance specification @eq-ql-specification is correct, then \( \B_n=\A_n \)
   by @thm-ql-score(d), and @eq-ql-sandwich-limit becomes
   \( \A_n^{1/2}(\hbeta-\bbeta^0)\to\Normal_p(\bzero,\I) \): the **model-based** covariance
   \( \A_n^{-1}=\phi(\X\T\W\X)^{-1} \) is correct, and by @prp-ql-optimality no linear
   estimating function has a larger Godambe information — hence, whenever (Q1)–(Q5) hold
   for that estimating function too, none has a smaller asymptotic covariance.

3. Let \( \hat c_i=w_ih'(\hat\eta_i)/V(\hat\mu_i) \) and define the **empirical sandwich**
   \[
   \hat{\bSigma}_{\text{rob}}
   =\bigl(\X\T\hat{\W}\X\bigr)^{-1}
   \Bigl\{\sum_{i=1}^n\hat c_i^{\,2}\,(y_i-\hat\mu_i)^2\,\x_{(i)}\x_{(i)}\T\Bigr\}
   \bigl(\X\T\hat{\W}\X\bigr)^{-1},
   \]{#eq-ql-robust}

   in which \( \phi \) cancels. Then
   \( n\hat{\bSigma}_{\text{rob}}-n\A_n^{-1}\B_n\A_n^{-1}\to\bzero \) in probability, and
   @eq-ql-sandwich-limit holds with \( \hat{\bSigma}_{\text{rob}}^{-1/2}(\hbeta-\bbeta^0) \)
   in place of the left-hand side.
:::

:::

::: {.proof}
*Step 1: the quasi-score is asymptotically normal.* By @thm-ql-score(a) and (e),
\( \bU_n(\bbeta^0)=\phi^{-1}\sum_iw_ia(\eta_i^0)\,\xi_i\,\x_{(i)} \) with
\( \xi_i=Y_i-\mu_i^0 \) independent and mean zero and \( \Cov\{\bU_n(\bbeta^0)\}=\B_n \).
Fix a unit vector \( \bv \) and write \( \bv\T\B_n^{-1/2}\bU_n=\sum_id_{ni}\xi_i \) with
\( d_{ni}=\phi^{-1}w_ia(\eta_i^0)\,\bv\T\B_n^{-1/2}\x_{(i)} \), so that
\( \sum_id_{ni}^2\sigma_i^2=1 \). By (Q4), \( \sigma_i^2\le K^{1/2} \) and hence
\( \sum_id_{ni}^2\ge K^{-1/2} \), while (Q2)–(Q4) give
\( d_{ni}^2\le\phi^{-2}w_{\max}^2a_{\max}^2C^2\lambda_{\max}(\B_n^{-1})=O(n^{-1}) \),
because \( \lambda_{\min}(n^{-1}\B_n)\to\lambda_{\min}(\bar{\B})>0 \). So
\( \max_id_{ni}^2/\sum_id_{ni}^2\to0 \), which is @eq-het-no-dominant-weight, and the
fourth moments are bounded: @lem-het-clt(b) with \( \delta=2 \) and the Cramér–Wold device
give \( \B_n^{-1/2}\bU_n(\bbeta^0)\to\Normal_p(\bzero,\I) \).

*Step 2: the expansion.* By (Q3), \( \bbeta\mapsto\bU_n(\bbeta) \) is continuously
differentiable on \( \mathcal N \), so on \( \{\hbeta\in\mathcal N\} \),
\[
\bzero=\bU_n(\hbeta)=\bU_n(\bbeta^0)+\bar{\mathbf{J}}_n\,(\hbeta-\bbeta^0),
\]
where
\[
\bar{\mathbf{J}}_n=\int_0^1\frac{\partial\bU_n}{\partial\bbeta\T}
\Bigl(\bbeta^0+t(\hbeta-\bbeta^0)\Bigr)\,dt .
\]

*Step 3: the Jacobian.* By the proof of @thm-ql-score(c),
\[
\frac{\partial\bU_n}{\partial\bbeta\T}(\bbeta)
=\frac1\phi\sum_{i=1}^nw_i\bigl\{a'(\eta_i)\{y_i-\mu_i(\bbeta)\}
-a(\eta_i)h'(\eta_i)\bigr\}\x_{(i)}\x_{(i)}\T ,
\]
with \( \eta_i=\x_{(i)}\T\bbeta \); take \( \norm{\bbeta-\bbeta^0}\le\delta \), so that
every \( \eta_i \) moves by at most \( C\delta \). The second term is nonrandom and equals
\( -\A_n \) at \( \bbeta^0 \), and since \( ah' \) is uniformly continuous by (Q3) the
two differ, after division by \( n \), by at most
\( \phi^{-1}w_{\max}C^2\,\omega(C\delta) \) with \( \omega \) the modulus of continuity.
In the first term write \( y_i-\mu_i(\bbeta)=\xi_i+\{\mu_i^0-\mu_i(\bbeta)\} \): the
displacement is at most \( h'_{\max}C\delta \), contributing \( O(\delta) \) to
\( n^{-1} \) times the norm, while
\( n^{-1}\phi^{-1}\sum_iw_ia'(\eta_i)\xi_i\x_{(i)}\x_{(i)}\T \) has entries that
average independent mean-zero variables of bounded variance, hence \( o_p(1) \) by
Chebyshev. By (Q5) \( \delta \) may tend to zero, so
\( -n^{-1}\bar{\mathbf{J}}_n\to\bar{\A} \) in probability, with \( \bar{\A} \)
nonsingular.

*Step 4: assembling.* With probability tending to one \( \bar{\mathbf{J}}_n \) is
invertible and \( \hbeta-\bbeta^0=(-\bar{\mathbf{J}}_n)^{-1}\bU_n(\bbeta^0) \), so
\[
\B_n^{-1/2}\A_n(\hbeta-\bbeta^0)
=\Bigl(\B_n^{-1/2}\A_n(-\bar{\mathbf{J}}_n)^{-1}\A_n^{-1}\B_n^{1/2}\Bigr)\,
\B_n^{-1/2}\bU_n(\bbeta^0).
\]
By (Q4) and Step 3 the bracket tends in probability to \( \I \), since \( n^{-1}\A_n \),
\( n^{-1}\B_n \), their square roots and inverses all converge to fixed nonsingular
limits. Slutsky and Step 1 give @eq-ql-sandwich-limit, and (b) is the substitution
\( \B_n=\A_n \).

*Step 5: the empirical middle.* Write \( c_i=w_ia(\eta_i^0) \), so the middle of
@eq-ql-robust at the true parameter is \( \sum_ic_i^2\xi_i^2\x_{(i)}\x_{(i)}\T \) and
\( \phi^2\B_n=\sum_ic_i^2\sigma_i^2\x_{(i)}\x_{(i)}\T \). Each entry of
\( n^{-1}\sum_ic_i^2(\xi_i^2-\sigma_i^2)\x_{(i)}\x_{(i)}\T \) is an average of independent
mean-zero variables whose variances are bounded by \( c_{\max}^4C^4K \), so it is
\( o_p(1) \) by Chebyshev; this is Step 2 of the proof of @thm-het-sandwich. Replacing
\( \xi_i \) by \( y_i-\hat\mu_i \) and \( c_i \) by \( \hat c_i \) costs a further
\( o_p(1) \), both displacements being bounded uniformly in \( i \) by a constant times
\( \norm{\hbeta-\bbeta^0} \), so Step 3 of that proof applies verbatim. The same
reasoning gives
\( n^{-1}\X\T\hat{\W}\X\to\phi\bar{\A} \), and combining the three factors proves (c).
:::

::: {.remark}
[What (Q5) hides, and when it is free]

Assuming a consistent root assumes something real: the equations can have no solution, or
several. Under the canonical link of the variance function \( Q \) is strictly concave, and
the argument for @thm-glm-asymptotics gives existence and consistency without (Q5) as soon
as \( \lambda_{\min}(\A_n)\to\infty \); for a general link the conditions of Fahrmeir
and Kaufmann (1985), quoted in @prp-glm-noncanonical, do the same and need only
unbiasedness. Conditions (Q2)–(Q4) are the quasi-likelihood forms of the bounded-leverage
and moment conditions of @thm-het-sandwich, and fail in the same places.
:::

## Model-based or robust

The two covariance estimates answer different questions. Under the log link with
\( V(\mu)=\mu \) and \( w_i=1 \), where \( \W_{ii}=\mu_i \) and \( a(\eta_i)=1 \),
\[
\A_n^{-1}\B_n\A_n^{-1}
=\Bigl(\sum_i\mu_i\x_{(i)}\x_{(i)}\T\Bigr)^{-1}
\Bigl(\sum_i\sigma_i^2\,\x_{(i)}\x_{(i)}\T\Bigr)
\Bigl(\sum_i\mu_i\x_{(i)}\x_{(i)}\T\Bigr)^{-1},
\]
against the model-based \( \phi(\sum_i\mu_i\x_{(i)}\x_{(i)}\T)^{-1} \). Only the *shape*
of the dependence on \( \mu \) matters, since a constant factor is absorbed into
\( \phi \): a quasi-Poisson analysis is right about the standard errors whenever the true
variance is proportional to the mean, whatever the constant, and wrong for any other
shape.

::: {#exm-ql-coverage}
[Coverage when the variance function is right and when it is wrong]

A quasi-Poisson model \( \log\mu_i=\beta_0+\beta_1x_i \) was fitted to data simulated
with the same mean and two true variances, \( \Var(Y_i)=\mu_i+2\mu_i^2 \) (quadratic) and
\( \Var(Y_i)=9\mu_i \) (linear), over \( 3000 \) replicates at each of
seven sample sizes. [Figure 38.3.1](03-sandwich.html#fig-ql-coverage) records the coverage
of nominal \( 95\% \) intervals for \( \beta_1 \) from four standard errors: Poisson with
\( \phi=1 \), model-based scaled by \( \hat\phi \), the empirical
sandwich @eq-ql-robust, and a leverage-corrected sandwich described below.

Ignoring the dispersion is fatal in both panels and does not improve with \( n \):
\( 0.334 \) at \( n=20 \) and \( 0.323 \) at
\( n=1280 \) in the quadratic case, \( 0.513 \) and
\( 0.482 \) in the linear case. A linear true variance makes the
specification correct and the \( \hat\phi \)-scaled interval already good at
\( n=20 \), covering \( 0.919 \) against the sandwich's
\( 0.880 \); a quadratic one stops it improving, at
\( 0.914 \) for \( n=1280 \), while the sandwich climbs to
\( 0.941 \), buying correctness in the limit at the price of noise in
small samples, as in @exm-het-hc-coverage.
:::

::: {when-format="html"}
![**Figure 38.3.1.** Coverage of nominal \( 95\% \) intervals for the slope of a
quasi-Poisson fit with a correct mean model. In (a) the true variance is quadratic in the
mean, so the working variance function has the wrong shape; in (b) it is proportional to
the mean.](sandwich_coverage.svg){#fig-ql-coverage width=100%}
:::

::: {when-format="pdf"}
![Coverage of nominal \( 95\% \) intervals for the slope of a
quasi-Poisson fit with a correct mean model. In (a) the true variance is quadratic in the
mean, so the working variance function has the wrong shape; in (b) it is proportional to
the mean.](sandwich_coverage.pdf){width=100%}
:::

```{.python .run #cell-sandwich-coverage-coverage}
import numpy as np

BETA = np.array([1.5, 1.0])
HALF = 2.0            # the regressor runs over [-HALF, HALF]

def quasi_poisson(X, y, steps=25):
    """IRLS for the quasi-score equations with a log link and V(mu) = mu."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(max(y.mean(), 0.5))
    for _ in range(steps):
        mu = np.exp(X @ beta)
        z = X @ beta + (y - mu) / mu
        beta = np.linalg.solve((X.T * mu) @ X, (X.T * mu) @ z)
    return beta, np.exp(X @ beta)

def slope_variances(X, y):
    """The four estimates of Var(b1): naive, phi-scaled, sandwich, leverage-corrected."""
    beta, mu = quasi_poisson(X, y)
    n, p = X.shape
    bread = np.linalg.inv((X.T * mu) @ X)
    root = np.sqrt(mu)[:, None] * X
    hat = np.sum(root * (root @ bread), axis=1)
    r2 = (y - mu) ** 2
    phi = np.sum(r2 / mu) / (n - p)
    meat = lambda w: bread @ ((X * w[:, None]).T @ X) @ bread
    return beta[1], np.array([bread[1, 1], phi * bread[1, 1],
                              meat(r2)[1, 1], meat(r2 / (1 - hat) ** 2)[1, 1]])

def coverage(n, quadratic, reps, rng, kappa=0.5, phi=9.0):
    """Coverage of nominal 95% normal intervals for b1 over `reps` simulated data sets."""
    x = np.linspace(-HALF, HALF, n)
    X = np.column_stack([np.ones(n), x])
    mu = np.exp(X @ BETA)
    shape = np.full(n, kappa) if quadratic else mu / (phi - 1.0)
    hit = np.zeros(4)
    for _ in range(reps):
        y = rng.negative_binomial(shape, shape / (shape + mu)).astype(float)
        b1, var = slope_variances(X, y)
        hit += np.abs(b1 - BETA[1]) <= 1.96 * np.sqrt(var)
    return hit / reps

rng = np.random.default_rng(3821)
names = ["Poisson", "phi-scaled", "sandwich", "corrected"]
for n in (25, 100, 400):
    quad = coverage(n, True, 300, rng)
    lin = coverage(n, False, 300, rng)
    print(f"n = {n:4d}  quadratic variance " + "  ".join(
        f"{k} {v:.3f}" for k, v in zip(names, quad)))
    print(f"          linear variance    " + "  ".join(
        f"{k} {v:.3f}" for k, v in zip(names, lin)))
```

## The dispersion, the quasi-deviance and tests

The model-based covariance needs \( \phi \), and so does any \( F \) test; neither
\( \hbeta \) nor the sandwich @eq-ql-robust does. Since \( \phi \) cancels
from @eq-ql-equations and scales \( \A_n \) by \( \phi^{-1} \), it is orthogonal to
\( \bbeta \): estimating it multiplies the asymptotic covariance of \( \hbeta \) by a
consistent factor rather than inflating it.

::: {#prp-ql-dispersion}
[Estimating the dispersion, and quasi-deviance tests]

Assume @eq-ql-specification with the mean model correct.

::: {.enumerate options="label=(\alph*)"}
1. \( \E\bigl[w_i(Y_i-\mu_i^0)^2/V(\mu_i^0)\bigr]=\phi \) for every \( i \), so the
   **Pearson estimator** \( \hat\phi_P=X^2/(n-p) \) with
   \( X^2=\sum_iw_i(y_i-\hat\mu_i)^2/V(\hat\mu_i) \) — the statistic @eq-glm-pearson — is
   consistent for \( \phi \) under (Q1)–(Q5). The divisor \( n-p \) is an analogy with
   \( s^2 \) (@thm-lm-sigma2), exact only in the normal linear model.

2. The **deviance estimator** \( \hat\phi_D=D(\y,\hat{\bmu})/(n-p) \), with \( D \) the
   quasi-deviance @eq-ql-quasi-deviance, is *not* in general consistent for \( \phi \). The
   two agree only so far as a third-order remainder is negligible: large binomial groups
   and large Poisson means, not small counts. Only \( \hat\phi_P \) has expectation
   exactly \( \phi \) at the true mean.

3. For nested mean models \( \C(\X_0)\subseteq\C(\X) \) with
   \( \rank(\X)-\rank(\X_0)=q \), fitted with the same \( V \), the **quasi-\( F \)
   statistic** is
   \[
   F=\frac{\{D(\y,\tilde{\bmu})-D(\y,\hat{\bmu})\}/q}{\hat\phi_P}.
   \]{#eq-ql-quasi-f}

   Under the reduced model and (Q1)–(Q5), \( qF\to\chi^2(q) \), so referring \( F \) to
   \( F(q,n-p) \) is asymptotically correct and is the usual small-sample convention. If
   the variance specification is wrong the limit fails and the Wald statistic built
   from @eq-ql-robust should be used instead.

4. No convention makes \( -2Q(\hbeta) \) a log-likelihood, so @def-sel-aic-bic does not
   apply. The criterion \( \text{QAIC}=D(\y,\hat{\bmu})/\hat\phi+2p \), with
   \( \hat\phi \) fixed at its value in the largest model, is in common use, but it is a
   convention, not a theorem. Cross-validation (@def-sel-cv) needs no likelihood, and the
   analogue for correlated responses, the QIC of Chapter 40, rests on the same
   convention.
:::

:::

::: {.proof}
(a) Under @eq-ql-specification, \( \E(Y_i-\mu_i^0)^2=\phi V(\mu_i^0)/w_i \), which is the
identity. The \( T_i=w_i(Y_i-\mu_i^0)^2/V(\mu_i^0) \) are independent with mean
\( \phi \) and bounded variances, so \( n^{-1}\sum_iT_i\to\phi \) by Chebyshev; replacing
\( \mu_i^0 \) by \( \hat\mu_i \) costs \( o_p(1) \), by Step 5 of the proof
of @thm-ql-asymptotics with \( \x_{(i)}\x_{(i)}\T \) replaced by \( 1 \), and
\( n/(n-p)\to1 \).

(b) Each term of @eq-ql-quasi-deviance is
\( 2w_i\int_{\hat\mu_i}^{y_i}(y_i-t)V(t)^{-1}\,dt \). Expanding \( 1/V(t) \) about
\( t=y_i \) and integrating gives
\( w_i(y_i-\hat\mu_i)^2/V(y_i)+O(w_i\lvert y_i-\hat\mu_i\rvert^3) \), the calculation
recorded after @eq-glm-pearson: a Pearson term with \( V(y_i) \) in place of
\( V(\hat\mu_i) \). The two estimators therefore differ by the remainder and by the gap
between those denominators, both negligible when \( \lvert y_i-\hat\mu_i\rvert \) is
small relative to \( \mu_i \) and neither otherwise; there is no exact unbiasedness to
appeal to, since \( \E d_i\ne\phi \) in general.

(c) By @eq-ql-quasi-deviance the difference is
\( 2\phi\{Q(\hbeta)-Q(\tilde{\bbeta})\} \). The proof of @prp-glm-three-tests used the
log-likelihood only through the quadratic expansion @eq-glm-quadratic-expansion, whose
ingredients are asymptotic normality of the score and the identity between its variance and
minus the expected Hessian. Under
@eq-ql-specification, Step 1 of the proof of @thm-ql-asymptotics supplies the first and
@thm-ql-score(d) the second, with \( \A_n \) in both roles, so that argument goes through
word for word with \( Q \) in place of \( \ell \). Hence
\( \{D(\y,\tilde{\bmu})-D(\y,\hat{\bmu})\}/\phi\to\chi^2(q) \); part (a) lets
\( \hat\phi_P \) replace \( \phi \), and \( qF(q,n-p)\to\chi^2(q) \) makes the \( F \)
reference asymptotically equivalent.

(d) is a definition and a warning, not a claim.
:::

::: {#exm-ql-strikes-errors}
[Four standard errors for the strike data]

Continue @exm-ql-strikes-fit. The Pearson statistic is \( 2393.5 \)
on \( 60 \) degrees of freedom, so
\( \hat\phi_P=39.89 \), against
\( \hat\phi_D=37.22 \) from the quasi-deviance. The four standard
errors for the production coefficient \( -7.680 \) are

| covariance | standard error | \( t \) |
|---|---:|---:|
| Poisson, \( \phi=1 \) | \( 0.406 \) | \( -18.92 \) |
| model-based, \( \hat\phi_P \) | \( 2.564 \) | \( -3.00 \) |
| empirical sandwich | \( 2.573 \) | \( -2.98 \) |
| leverage-corrected sandwich | \( 2.896 \) | \( -2.65 \) |

The Poisson standard error is too small by a factor of about six. More interesting are the
second and third rows: the sandwich exceeds the model-based value by a factor of only
\( 1.004 \), which is evidence *for* the quasi-Poisson variance
function, so the model-based standard error is valid and slightly more efficient. The
fourth row divides each squared residual by \( (1-h_{ii})^2 \), with \( h_{ii} \) the
leverages of the weighted hat matrix of @def-glm-residuals; the largest is
\( 0.192 \), enough to matter at \( n=62 \).
:::

```{.python .run #cell-strikes-covariances}
import numpy as np
import statsmodels.api as sm

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
pearson = np.sum((y - mu) ** 2 / mu)                          # X^2 with V(mu) = mu
phi_pearson = pearson / (n - p)
W = mu                                        # log link with V(mu) = mu: W_ii = mu_i
bread = np.linalg.inv((X.T * W) @ X)          # (X' W X)^{-1}
root = np.sqrt(W)[:, None] * X
hat = np.sum(root * (root @ bread), axis=1)   # leverages of the weighted hat matrix
resid = y - mu

cov_model = bread                             # phi = 1: the Poisson standard errors
cov_quasi = phi_pearson * bread               # phi estimated from the same fit
cov_sand = bread @ ((X * resid[:, None] ** 2).T @ X) @ bread
cov_corr = bread @ ((X * (resid ** 2 / (1 - hat) ** 2)[:, None]).T @ X) @ bread

for name, C in [("Poisson", cov_model), ("quasi", cov_quasi),
                ("sandwich", cov_sand), ("corrected", cov_corr)]:
    se = np.sqrt(np.diag(C))
    print(f"{name:10s} slope {beta_lin[1]:8.3f}   se {se[1]:6.3f}   t {beta_lin[1] / se[1]:7.2f}")
print(f"largest leverage {hat.max():.4f}, sum {hat.sum():.4f}")
```

## Corrections, clusters and resampling

Three qualifications finish the section.

**Small samples.** The empirical sandwich is biased downwards, for the reason
@prp-het-hc-facts(a) identified in the linear model, and the remedy is the same: divide the
\( i \)th squared residual by \( 1-h_{ii} \) or by \( (1-h_{ii})^2 \), with \( h_{ii} \)
the diagonal of \( \hat{\W}^{1/2}\X(\X\T\hat{\W}\X)^{-1}\X\T\hat{\W}^{1/2} \), giving
the analogues of HC2 and HC3 of @def-het-hc. The second is the "corrected" line
of @exm-ql-coverage, better than the plain sandwich at \( n=80 \) in the quadratic panel,
\( 0.906 \) against \( 0.891 \); Kauermann and
Carroll (2001) give the leading bias term. The middle of @eq-ql-robust is itself an average
with very unequal weights, so it behaves like a variance estimate on far fewer than
\( n-p \) degrees of freedom (@exm-het-engel-leverage).

**Clusters.** If the observations come in \( G \) independent groups with arbitrary
within-group dependence, the middle of @eq-ql-robust becomes
\[
\sum_{g=1}^G\Bigl(\sum_{i\in g}\hat c_i(y_i-\hat\mu_i)\x_{(i)}\Bigr)
\Bigl(\sum_{i\in g}\hat c_i(y_i-\hat\mu_i)\x_{(i)}\Bigr)\T ,
\]
the cluster-robust form @eq-cls-cluster-sandwich with the quasi-score's coefficients in
place of the rows of \( \X \). Its consistency needs \( G\to\infty \), by Step 5 above
with the cluster in place of the observation. The quasi-score is then no longer optimal,
because the within-group correlations are ignored; Chapter 40 puts a generalized estimating
equation in its place.

**Resampling.** The bootstrap of [Chapter 23](../ch23-resampling-inference/index.html)
applies unchanged: whole observations (@prp-bs-case-linear), the wild
bootstrap (@def-bs-wild), or whole clusters (@exr-bs-wild-cluster) for small \( G \).

## Exercises

### A. Check your understanding

::: {#exr-ql-when-equal}
[A1]

Under a log link and \( V(\mu)=\mu \) with \( w_i\equiv1 \), suppose the true variances are
\( \sigma_i^2=\tau\mu_i \) for an unknown \( \tau>0 \). Show that the sandwich and the
model-based covariance have the same probability limit, and identify what \( \hat\phi_P \)
estimates. Repeat for \( \sigma_i^2=\tau\mu_i^2 \) and show that the two then differ unless
all \( \mu_i \) are equal.
:::

::: {#exr-ql-phi-cancels}
[A2]

Verify from @eq-ql-robust that the empirical sandwich does not change if \( V \) is replaced
by \( cV \) for a constant \( c>0 \), and explain why this is to be expected.
:::

### B. Practice

::: {#exr-ql-intercept-sandwich}
[B1]

Take an intercept-only quasi-Poisson model, \( \log\mu_i=\beta \) for all \( i \). Show that
\( \hat\mu_i=\bar y \), that the model-based variance of \( \hat\beta \) is
\( \hat\phi_P/(n\bar y) \), and that the sandwich variance is
\( \sum_i(y_i-\bar y)^2/(n\bar y)^2 \). Show that the two agree exactly when
\( \hat\phi_P \) is computed with divisor \( n \) rather than \( n-p \).
:::

::: {#exr-ql-efficiency-loss}
[B2]

Let \( \Var(Y_i)=\sigma_i^2 \) and consider the quasi-score built from a working variance
function \( V \). Using @prp-ql-optimality with the *true* variances, show that
\( \A_n^{-1}\B_n\A_n^{-1} \) — the asymptotic covariance of @thm-ql-asymptotics(a) — is
never smaller, in the nonnegative definite order, than the inverse Godambe information
\( \bigl(\sum_ih'(\eta_i)^2\sigma_i^{-2}\x_{(i)}\x_{(i)}\T\bigr)^{-1} \) of the best
linear estimating function, and that the two agree exactly when
\( \sigma_i^2\propto V(\mu_i)/w_i \).
:::

::: {.solution}
Repeat the proof of @prp-ql-optimality with \( \bSigma=\diag(\sigma_i^2) \) in place of
\( \diag\{\phi V(\mu_i)/w_i\} \); only its last display used the specification. The
quasi-score built from \( V \) is a linear estimating function with
\( \mathbf{a}_i=w_ih'(\eta_i)\x_{(i)}/\{\phi V(\mu_i)\} \), so by @thm-ql-score(c),(e)
its Godambe information is \( \A_n\B_n^{-1}\A_n \), at most the maximum
\( \mathbf{K}\T\mathbf{K}=\sum_ih'(\eta_i)^2\sigma_i^{-2}\x_{(i)}\x_{(i)}\T \);
inverting reverses the order. Equality needs coefficients proportional to
\( h'(\eta_i)\x_{(i)}/\sigma_i^2 \), that is
\( w_i/V(\mu_i)\propto1/\sigma_i^2 \).
:::

### C. Going deeper

::: {#exr-ql-information-test}
[C1]

A large discrepancy between \( \hat\phi_P(\X\T\hat{\W}\X)^{-1} \) and @eq-ql-robust is
evidence that the variance specification is wrong. Propose a test statistic from the
difference of the two matrices, say
\( \tr\{(\X\T\hat{\W}\X)\hat{\bSigma}_{\text{rob}}\}-p\hat\phi_P \), argue informally
that its mean is about zero when the specification is correct, and say what stands in the
way of a usable null distribution. (This is White's information matrix test;
compare @exr-het-white-test.)
:::
