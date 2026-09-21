# Diagnostics

A fitted generalized linear model can be wrong in three ways: the linear predictor may
omit a term or use a covariate on the wrong scale, the link may be wrong, or the variance
function may be wrong. A single number — the deviance, say — cannot tell these apart, and
for ungrouped binary data it cannot even tell that something is wrong (@prp-bin-fit(b)).
What follows is a set of plots and tests, each aimed at one of the three, and a worked case
in which they give different answers.

## Residuals in practice

The three residual types of @def-glm-residuals differ in what they standardize away. Their
common first-order behaviour is the generalized linear version
of @eq-res-residual-moments.

::: {#prp-prc-diagnostics}
[What the standard diagnostics measure]

Let a generalized linear model be fitted by maximum likelihood, let
\( \hat{\W}=\diag\{w_ih'(\hat\eta_i)^2/V(\hat\mu_i)\} \) be the working weights,
\( \bH=\hat{\W}^{1/2}\X(\X\T\hat{\W}\X)^{-1}\X\T\hat{\W}^{1/2} \) the weighted hat matrix
and \( h_{ii} \) its diagonal.

::: {.enumerate options="label=(\alph*)"}
1. *(Moments of the Pearson residuals.)* To first order in the expansion behind
   @thm-glm-asymptotics, \( \E(r_i^{P})\approx0 \) and
   \( \Var(r_i^{P})\approx1-h_{ii} \), and the standardized residual
   \( r_i^{P}/\sqrt{1-h_{ii}} \) has approximate variance one. Moreover
   \( \sum_ih_{ii}=p \), so \( h_{ii} \) averages \( p/n \).

2. *(Level of a simulated envelope.)* Let \( |r|_{(1)}\le\dots\le|r|_{(n)} \) be the
   ordered absolute residuals, and let \( B \) independent data sets be generated from the
   fitted model and refitted, giving vectors of ordered absolute residuals whose
   \( k \)th entries are \( R^{(1)}_k,\dots,R^{(B)}_k \). If the model is correct and
   the fitted parameters are the true ones, and the \( B+1 \) values
   \( |r|_{(k)},R^{(1)}_k,\dots,R^{(B)}_k \) are exchangeable with no ties, then
   \[
\Pr\bigl\{|r|_{(k)}<\min_bR^{(b)}_k\ \text{or}\ |r|_{(k)}>\max_bR^{(b)}_k\bigr\}
=\frac{2}{B+1}
\]{#eq-prc-envelope}

   for each fixed \( k \). The probability that *some* \( k \) leaves the band is larger
   than \( 2/(B+1) \), and no better than the union bound \( \min\{1,2n/(B+1)\} \)
   without further argument.

3. *(Influence.)* The one-step deletion results of @thm-res-deletion hold with
   \( \X \) replaced by \( \hat{\W}^{1/2}\X \) and \( \y \) by the working response
   \( \bz \) of @eq-glm-working-response, giving the Cook-type measure
   \( C_i=(r_i^{P})^2h_{ii}/\{p(1-h_{ii})^2\} \) of @def-res-cooks. The exact deletion
   identities do not hold, because removing a case changes the weights as well.

4. *(Checking the link.)* Let \( \hat\eta_i \) be the fitted linear predictor and refit the
   model with \( \hat\eta_i^2 \) added as a covariate. Under the fitted link the extra
   coefficient is zero, so the deviance drop is asymptotically \( \chi^2(1) \); a
   significant drop says the link is wrong, not which link is right.

5. *(Checking a covariate: two plots.)* The partial residual for covariate \( j \) on the
   link scale is \( u_{ij}=\hat\beta_jx_{ij}+(y_i-\hat\mu_i)/h'(\hat\eta_i) \). Its
   weighted least squares regression on \( [\bone,\x_j] \) with weights \( \hat{\W} \)
   has slope \( \hat\beta_j \) exactly, so a smooth of \( u_{ij} \) against
   \( x_{ij} \) that departs from a straight line indicates a wrong scale for
   \( x_j \), as in @prp-res-partial-residual. The added-variable plot is the other
   half of the pair: let \( \tilde{\bz} \) and \( \tilde{\x}_j \) be the residuals
   of the working response \( \bz \) and of \( \x_j \) after
   \( \hat{\W} \)-weighted regression on the remaining columns \( \X_{-j} \). Then
   the \( \hat{\W} \)-weighted regression of \( \tilde{\bz} \) on
   \( \tilde{\x}_j \) has slope exactly \( \hat\beta_j \), so curvature or a single
   dominating point in that plot is a statement about \( \hat\beta_j \) itself.

6. *(Checking the variance function.)* Suppose the assumed variance function is correct
   and each observation carries enough information for \( r_i^{P} \) to be approximately
   normal: grouped data with large group sizes, or counts with means that are not small.
   (This fails for ungrouped binary data, where \( r_i^{P} \) takes two values and the
   constant below is wrong.) Then
   \( \E|r_i^{P}|\approx\sqrt{2/\pi} \) for every \( i \), so a plot of \( |r_i^{P}| \)
   against \( \hat\eta_i \) should be level. A systematic trend indicates a variance
   function of the wrong shape; a level plot at the wrong height indicates a dispersion
   \( \phi\ne1 \).
:::
:::

::: {.proof}
(a) The expansion used in @thm-glm-asymptotics gives, to first order,
\( \hbeta-\bbeta\approx(\X\T\W\X)^{-1}\X\T\W\bD^{-1}(\y-\bmu) \) and hence
\( \hat\bmu-\bmu\approx\bD\X(\hbeta-\bbeta) \), where \( \bD=\diag(d\mu_i/d\eta_i) \). Put
\( \bu=\W^{1/2}\bD^{-1}(\y-\bmu) \), whose entries are the Pearson residuals evaluated at
the true parameter, with \( \Cov(\bu)=\I \) by the definition of \( \W \). Then
\[
\begin{aligned}
\W^{1/2}\bD^{-1}(\y-\hat\bmu)
&\approx\bu-\W^{1/2}\X(\X\T\W\X)^{-1}\X\T\W^{1/2}\bu\\
&=(\I-\bH)\bu ,
\end{aligned}
\]
so \( \E(r_i^{P})\approx0 \) and \( \Var(r_i^{P})\approx(\I-\bH)_{ii}=1-h_{ii} \). The
trace identity is @prp-proj-trace-rank applied to \( \bH \), which is a projection matrix.

(b) The event that \( |r|_{(k)} \) is the strict minimum or the strict maximum of the
\( B+1 \) exchangeable values has probability \( 2/(B+1) \), because each of the
\( B+1 \) values is equally likely to be the smallest and equally likely to be the largest,
and with no ties these two events are disjoint for \( B\ge1 \). The final sentence is the
usual multiplicity statement: the events for different \( k \) are highly dependent, so
the probability that some \( k \) leaves the band lies between \( 2/(B+1) \) and
\( \min\{1,2n/(B+1)\} \) by the union bound, and is not \( 2/(B+1) \).

(c) @thm-res-deletion is an algebraic identity about a least squares fit. At convergence
the generalized linear fit *is* a weighted least squares fit of \( \bz \) on \( \X \) with
weights \( \hat{\W} \) (@thm-glm-irls), so the identity applies to that problem, which is
what the substitution says. It is only one step of the algorithm, since refitting without
case \( i \) changes \( \hat\bmu \) and therefore \( \hat{\W} \).

(d) The enlarged model contains the fitted one at coefficient zero and adds one parameter,
so @thm-glm-deviance gives the \( \chi^2(1) \) limit under the smaller model. Note that
\( \hat\eta_i \) is itself estimated. Because \( \hat\eta^2 \) is a function of the
fitted values alone, the normal-theory analogue is @thm-cor-fitted-regressors, which makes
the corresponding \( F \) test *exact*; here the response is not normal, the limit is the
asymptotic one, and the effect of estimating \( \hat\eta_i \) is of smaller order under
the conditions of @thm-glm-asymptotics.

(e) The working residual is \( (y_i-\hat\mu_i)/h'(\hat\eta_i)=z_i-\hat\eta_i \)
by @eq-glm-working-response. The weighted least squares fit of \( \bz \) on \( \X \) has
residual vector \( \bz-\hat{\boldsymbol{\upeta}} \) orthogonal to \( \X \) in the
\( \hat{\W} \) inner product, so \( \bu_j=\hat\beta_j\x_j+(\bz-\hat{\boldsymbol{\upeta}}) \)
regressed on \( [\bone,\x_j] \) in that inner product returns \( \hat\beta_j \) exactly, by
@prp-res-partial-residual(a) applied in the \( \hat{\W} \) metric (@thm-proj-A-projection).
For the added-variable statement, at convergence the fit is the weighted least squares
regression of \( \bz \) on \( \X \) with weights \( \hat{\W} \) (@thm-glm-irls), which
is an ordinary least squares problem in the \( \hat{\W} \) inner product; @thm-proj-fwl
applied in that inner product gives the claim.

(f) By (a), \( r_i^{P} \) has approximate mean zero and variance \( 1-h_{ii}\approx1 \)
when the leverages are small. Its approximate normality is the hypothesis of the
statement, and it is a statement about the single observation \( Y_i \), not about
\( \hbeta \): @thm-glm-asymptotics concerns the estimate and does not deliver it. For \( Z\sim\Normal(0,1) \), \( \E|Z|=\sqrt{2/\pi} \). If instead
\( \Var(Y_i)=\phi V(\mu_i)/w_i \) with \( \phi\ne1 \), every \( r_i^{P} \) is inflated by
\( \sqrt\phi \) and the plot is level at \( \sqrt{2\phi/\pi} \).
:::

::: {.remark}
[Why absolute residuals]

A quantile plot of signed residuals asks two questions at once, and for a discrete response
the answer to the first is dominated by the discreteness. Plotting \( |r|_{(k)} \) against
\( \Phi^{-1}\{(k+n-1/8)/(2n+1/2)\} \), the expected order statistics of \( |Z| \), asks
only about magnitude, which is what an outlier or an inflated dispersion changes. The
envelope turns the picture into a parametric bootstrap test (@prp-bs-parametric), of only
approximate level because the data sets are generated at \( \hbeta \).
:::

## A model that fails, and why

::: {#exm-prc-tvnews}
[Days of television news]

The 1996 American National Election Studies survey asked its
\( 944 \) respondents on how many of the past seven days they had
watched the television news. Treating the answer as a binomial count out of seven and
fitting a logit model in age, education and income gives a deviance of
\( 4447.3 \) on \( 940 \) degrees of freedom and a
Pearson statistic of \( 3598.2 \), so
\( \hat\phi=3.828 \): the model is not close to fitting.

Panel (a) of [Figure 39.2.1](02-diagnostics.html#fig-prc-envelope) is the half-normal plot
with an envelope from \( B=19 \) simulated data sets, refitted each time. Every one of
the \( 944 \) ordered residuals lies above the band — not a little too
large in the tail but too large everywhere, the signature of a dispersion problem rather
than of outliers.

The targeted checks of @prp-prc-diagnostics say where the trouble is not. Adding
\( \hat\eta^2 \) drops the deviance by \( 0.31 \) on one degree of
freedom, \( p=0.580 \): the logit link is not the problem. Adding a
quadratic in age drops it by \( 0.03 \), and panel (a) of
[Figure 39.2.2](02-diagnostics.html#fig-prc-diagnostics) shows the partial residual for age
with a smooth that is mildly concave but close to the fitted line, so the scale of age is
not the problem either. The leverages are all tiny — the largest
\( 0.0144 \) against an average of
\( 0.0042 \) — and the largest Cook-type distance is
\( 0.0970 \).

What is left is the variance function, and panel (b) makes the case. The local means of
\( |r^{P}| \) in ten equal bins of the linear predictor run from
\( 1.25 \) to \( 1.91 \), against a reference
height of \( \sqrt{2/\pi}\approx0.80 \), and they drift upward: the regression of
\( \log(r^{P})^2 \) on the linear predictor has slope
\( 0.381 \) with \( t=4.21 \), so the excess variance is not even
constant. The binomial model treats a respondent's seven days as seven independent trials
with a common probability, and they are nothing of the kind:
\( 161 \) respondents watched on none of the seven days and
\( 288 \) on all seven. Habits are not coin flips.

The repair is a variance function with a second parameter. A beta-binomial model, in which
each respondent has a latent probability from a beta distribution with mean
\( \pi_i \) and precision \( s \), gives \( \hat s=1.072 \), a
within-respondent correlation of \( \rho=1/(1+s)=0.483 \), and the variance
\( 7\pi_i(1-\pi_i)(7+s)/(1+s) \), about \( 3.90 \) times the
binomial one, against the \( \hat\phi=3.828 \) estimated directly.
Panel (b) of [Figure 39.2.1](02-diagnostics.html#fig-prc-envelope) repeats the half-normal
plot under that model: only \( 0.2 \) per cent of the ordered residuals now
lie above the envelope, although \( 22.8 \) per cent still fall below it,
against a pointwise level of \( 5 \) per cent on that side — a residual misfit at the
bottom of the scale that the beta-binomial does not remove.

The consequences for inference are what
[Chapter 38](../ch38-quasi-likelihood/index.html) predicts. The age coefficient barely
moves, from \( 0.0424 \) to \( 0.0405 \); its
standard error moves from \( 0.0018 \) to
\( 0.0033 \), and the quasi-binomial standard error — the binomial
one multiplied by \( \sqrt{\hat\phi} \) — is \( 0.0034 \), which
agrees with the beta-binomial one to two figures. The education coefficient is the
casualty: \( 0.0233 \) with standard error
\( 0.0178 \) under the binomial, already unconvincing at
\( t=1.3 \), becomes \( -0.0013 \) with standard error
\( 0.0343 \) — an effect of zero with three times the uncertainty.
:::

::: {when-format="html"}
![**Figure 39.2.1.** Half-normal plots of the absolute Pearson residuals with an envelope
from 19 simulated data sets, refitted each time. (a) The binomial model: every point lies
above the envelope. (b) The beta-binomial model: the upper tail is contained, the lower one
is not.](halfnormal_envelope.svg){#fig-prc-envelope width=100%}
:::

::: {when-format="pdf"}
![Half-normal plots of the absolute Pearson residuals with an envelope
from 19 simulated data sets, refitted each time. (a) The binomial model: every point lies
above the envelope. (b) The beta-binomial model: the upper tail is contained, the lower one
is not.](halfnormal_envelope.pdf){width=100%}
:::

```{.python .run #cell-diagnostics-fit}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, special, stats
anes = sm.datasets.anes96.load_pandas().data
m = 7.0                                                   # days in a week
y = anes["TVnews"].to_numpy(float)                        # successes out of seven
X = np.column_stack([np.ones(len(y)), anes["age"], anes["educ"], anes["income"]])
n, p = X.shape

fit = sm.GLM(np.column_stack([y, m - y]), X, family=sm.families.Binomial()).fit()
mu = fit.fittedvalues * m                                 # fitted expected days
eta = fit.predict(which="linear")
pearson = (y - mu) / np.sqrt(mu * (1 - mu / m))           # r^P, with phi = 1

print(f"deviance {fit.deviance:8.1f} on {n - p} degrees of freedom")
print(f"Pearson  {np.sum(pearson**2):8.1f};  dispersion estimate "
      f"{np.sum(pearson**2) / (n - p):.3f}")
```

```{.python .run #cell-diagnostics-envelope}
B = 19                                                    # simulated data sets


def halfnormal_quantiles(k):
    """The horizontal axis: expected half-normal order statistics."""
    i = np.arange(1, k + 1)
    return stats.norm.ppf((i + k - 0.125) / (2 * k + 0.5))


def envelope(draw, refit_abs_resid, seed):
    """Smallest and largest of B ordered |residual| vectors from the fitted model."""
    rng = np.random.default_rng(seed)
    out = np.array([np.sort(refit_abs_resid(draw(rng))) for _ in range(B)])
    return out.min(axis=0), out.max(axis=0)


def binomial_abs_resid(ysim):
    f = sm.GLM(np.column_stack([ysim, m - ysim]), X, family=sm.families.Binomial()).fit()
    mu_s = f.fittedvalues * m
    return np.abs((ysim - mu_s) / np.sqrt(mu_s * (1 - mu_s / m)))


q = halfnormal_quantiles(n)
lo, hi = envelope(lambda rng: rng.binomial(int(m), fit.fittedvalues).astype(float),
                  binomial_abs_resid, seed=3901)
ordered = np.sort(np.abs(pearson))
outside, above = np.mean((ordered > hi) | (ordered < lo)), np.mean(ordered > hi)
print(f"binomial model: {100 * outside:.1f}% of the ordered residuals leave the envelope, "
      f"{100 * above:.1f}% above it")
```

```{.python .run #cell-diagnostics-checks}
# (i) is the link right?  Add the squared linear predictor and test it.
fit_link = sm.GLM(np.column_stack([y, m - y]), np.column_stack([X, eta**2]),
                  family=sm.families.Binomial()).fit()
lr_link = fit.deviance - fit_link.deviance

# (ii) is the age term right?  A partial residual, on the link scale.
working = (y - mu) / (mu * (1 - mu / m))                  # (y - mu) dEta/dMu
partial_age = fit.params[1] * anes["age"].to_numpy() + working

# (iii) is the variance function right?  Regress log (r^P)^2 on the linear predictor.
lin = sm.OLS(np.log(pearson**2 + 1e-12), sm.add_constant(eta)).fit()

print(f"link test: deviance drop {lr_link:.2f} on 1 d.f., p = {stats.chi2.sf(lr_link, 1):.3f}")
print(f"slope of log (r^P)^2 on the linear predictor: {lin.params[1]:.3f} "
      f"(t = {lin.tvalues[1]:.2f})")
```

::: {when-format="html"}
![**Figure 39.2.2.** Three targeted checks of the binomial fit. (a) Partial residual for
age, with a smooth (solid) and the fitted linear term (dotted). (b) Absolute Pearson
residuals against the linear predictor, with their means in ten equal bins and the
reference height \( \sqrt{2/\pi} \). (c) Cook-type distances. Only (b) shows a
problem.](glm_diagnostics.svg){#fig-prc-diagnostics width=100%}
:::

::: {when-format="pdf"}
![Three targeted checks of the binomial fit. (a) Partial residual for
age, with a smooth (solid) and the fitted linear term (dotted). (b) Absolute Pearson
residuals against the linear predictor, with their means in ten equal bins and the
reference height \( \sqrt{2/\pi} \). (c) Cook-type distances. Only (b) shows a
problem.](glm_diagnostics.pdf){width=100%}
:::

```{.python .run #cell-diagnostics-betabinom}
scale = np.array([1.0, X[:, 1].std(), X[:, 2].std(), X[:, 3].std()])   # keep BFGS conditioned


def betabinom_fit(ysim, start):
    """Maximize the beta-binomial likelihood: logit mean X beta, precision s."""
    def negll(theta):
        pi = special.expit((X / scale) @ theta[:p])
        s = np.exp(theta[p])
        v = -np.sum(special.betaln(ysim + s * pi, m - ysim + s * (1 - pi))
                    - special.betaln(s * pi, s * (1 - pi)))
        return v if np.isfinite(v) else 1e12

    theta = np.append(start[:p] * scale, start[p])
    for _ in range(3):                                    # restart until the gradient is flat
        opt = optimize.minimize(negll, theta, method="BFGS")
        theta = opt.x
    return np.append(theta[:p] / scale, theta[p]), opt


theta_bb, opt = betabinom_fit(y, np.append(fit.params, 0.1))
beta_bb, s_bb = theta_bb[:p], np.exp(theta_bb[p])
rho = 1.0 / (1.0 + s_bb)                                  # correlation within the week
print(f"beta-binomial: s = {s_bb:.3f}, within-respondent correlation {rho:.3f}")
for j, name in enumerate(["intercept", "age", "educ", "income"]):
    print(f"  {name:10s} binomial {fit.params[j]:8.4f}   beta-binomial {beta_bb[j]:8.4f}")
```

## Which repair

@exm-prc-tvnews had three candidate repairs and took the third; the choice among them is
the subject of @prp-ql-overdispersion. A **quasi-likelihood** fit keeps the mean model and
the shape of the variance function and lets \( \phi \)
float (@def-ql-variance-function): it changes no coefficient, because the estimating
equations are the same, and multiplies every standard error by \( \sqrt{\hat\phi} \). A
**mixture** model, such as the beta-binomial here or the negative binomial
of @def-cnt-negbin, names a mechanism and is a genuine likelihood, so it has an AIC and a
likelihood ratio test, and it generally does move the coefficients. A **random effect** in
the linear predictor is the third route, developed in
[Chapter 40](../ch40-glmm-gee/index.html); it changes the meaning of the coefficients as
well as their values (@prp-gmm-marginal).

::: {.idea}
Read diagnostics for a generalized linear model as three separate questions. Is the mean
right? Is the link right? Is the variance function right? A plot that fails is useful only
when it belongs to one of the three.
:::

## Exercises

### A. Check your understanding

::: {#exr-prc-envelope-level}
[A1]

@prp-prc-diagnostics(b) gives the pointwise level \( 2/(B+1) \). What is \( B \) for a
level of \( 5 \) per cent? If instead the envelope is taken to be the second smallest and
second largest of \( B \) simulated values, what is the pointwise level?
:::

::: {.solution}
\( 2/(B+1)=0.05 \) gives \( B=39 \). With the second order statistics as limits the
observed value must be among the two smallest or two largest of \( B+1 \), which has
probability \( 4/(B+1) \), so \( B=79 \) is needed: a less jagged envelope at the cost of
more simulation.
:::

::: {#exr-prc-link-test-reading}
[A2]

The link test of @prp-prc-diagnostics(d) came out non-significant in @exm-prc-tvnews. Is
that evidence that the logit link is correct? What would a significant result have shown?
:::

### B. Practice

::: {#exr-prc-standardized}
[B1]

Show that in a normal linear model with \( \phi=\sigma^2 \) and the identity link,
@prp-prc-diagnostics(a) reduces to the exact statement
\( \Cov(\he)=\sigma^2(\I-\M) \) of @eq-res-residual-moments, and that the standardization
by \( \sqrt{1-h_{ii}} \) is then exact rather than approximate.
:::

::: {.solution}
With the identity link and constant variance \( h'(\eta)=1 \), \( V(\mu)=1 \),
\( w_i=1 \), so \( \W=\bD=\I \) and \( \bH=\X(\X\T\X)^{-1}\X\T=\M \). The expansion in
the proof becomes an identity, \( \hat\bmu-\bmu=\M(\y-\bmu) \), whence
\( \Cov(\y-\hat\bmu)=\sigma^2(\I-\M) \) and \( r_i^{P}=(y_i-\hat\mu_i)/\sigma \) has
variance exactly \( 1-h_{ii} \).
:::

::: {#exr-prc-betabinom-variance}
[B2]

Let \( Y\mid\Pi\sim\text{Binomial}(m,\Pi) \) with \( \Pi\sim\text{Beta}(s\pi,s(1-\pi)) \).
Show that \( \E(Y)=m\pi \) and
\( \Var(Y)=m\pi(1-\pi)(m+s)/(1+s) \), and that the ratio to the binomial variance is
\( 1+(m-1)\rho \) with \( \rho=1/(1+s) \). Evaluate at \( m=7 \) and the fitted
\( \hat s=1.072 \).
:::

::: {.solution}
\( \E\Pi=\pi \) and \( \Var\Pi=\pi(1-\pi)/(1+s) \) for that beta distribution. By the
tower property \( \E Y=m\E\Pi=m\pi \), and
\[
\Var Y=\E\{m\Pi(1-\Pi)\}+\Var(m\Pi)
=m\{\pi-\E\Pi^2\}+m^2\Var\Pi .
\]
With \( \E\Pi^2=\Var\Pi+\pi^2 \) this is
\( m\pi(1-\pi)-m\Var\Pi+m^2\Var\Pi=m\pi(1-\pi)\{1+(m-1)/(1+s)\} \), which is the stated
formula since \( 1+(m-1)/(1+s)=(m+s)/(1+s) \). The ratio to \( m\pi(1-\pi) \) is
\( 1+(m-1)\rho \). At \( m=7 \) and \( s=1.072 \), \( \rho=0.483 \) and the ratio is
\( 1+6\times0.483=3.90 \), close to the \( \hat\phi=3.828 \) that
@exm-prc-tvnews estimates directly.
:::

### C. Going deeper

::: {#exr-prc-simultaneous-envelope}
[C1]

@prp-prc-diagnostics(b) says the simultaneous level of an envelope exceeds the pointwise
one. Construct a simultaneous band from the same \( B \) simulations, by computing for
each simulated data set the maximum over \( k \) of the standardized deviation
\( \{R^{(b)}_k-\bar R_k\}/\text{sd}_k \) and using its \( (1-\alpha) \) quantile. Show
that the resulting band has simultaneous level at most \( \alpha \) under the same
exchangeability assumption, and say what it costs.
:::
