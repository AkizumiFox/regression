# Poisson regression

Counts are nonnegative integers, their spread grows with their level, and a
regressor acts on them multiplicatively. No transformation satisfies all three at
once ([Chapter 22](../ch22-transformations/index.html)): the square root of a
Poisson count has nearly constant variance but a mean no longer linear in the
regressors, the logarithm a nearly linear mean but a variance that still
moves (@exr-tr-poisson-scales), and neither is defined at zero. The generalized linear
model separates the requirements: the link takes the mean, the distribution takes
the variance, and nothing is done to the data.

## The Poisson as an exponential dispersion family

Write the Poisson mass function in the form @def-glm-edf asks for. For
\( y=0,1,2,\dots \) and \( \mu>0 \),
\[
\Pr(Y=y)=\frac{e^{-\mu}\mu^{y}}{y!}
=\exp\bigl\{y\log\mu-\mu-\log y!\bigr\},
\]
so the natural parameter is \( \theta=\log\mu \), the cumulant function is
\( b(\theta)=e^{\theta} \), the dispersion is \( \phi=1 \) and the weight is
\( w=1 \). By @prp-glm-moments,
\[
\E(Y)=b'(\theta)=e^{\theta}=\mu,\qquad
\Var(Y)=\phi\,b''(\theta)=e^{\theta}=\mu,
\]
so the variance function is \( V(\mu)=\mu \). The Poisson has no free dispersion
parameter: once the mean is fixed, so is the variance. That is the most
consequential fact in this chapter, and Sections
[37.4](04-negative-binomial.html) and [37.5](05-zero-inflation.html) are about what
to do when the data disagree with it. The canonical link \( \theta=\log\mu \) also
makes the mean multiplicative and keeps it positive at every linear predictor: for
once convenience and interpretation point the same way.

::: {#def-cnt-poisson}
[The Poisson log-linear model]

Let \( Y_1,\dots,Y_n \) be independent, with \( Y_i \) Poisson distributed with mean
\( \mu_i>0 \), and let \( \x_{(i)} \) be the \( i \)th row of a model matrix \( \X \)
of full column rank \( p \). The **Poisson log-linear model** is the generalized
linear model (@def-glm-model) with this random component, the log link and the
linear predictor \( \eta_i=\x_{(i)}\T\bbeta \):
\[
\log\mu_i=\x_{(i)}\T\bbeta,
\qquad\text{equivalently}\qquad
\mu_i=\exp\bigl(\x_{(i)}\T\bbeta\bigr)
=e^{\beta_1x_{i1}}\cdots e^{\beta_px_{ip}} .
\]{#eq-cnt-model}

:::

The second form carries the interpretation. Increasing the \( j \)th regressor by
one unit, the others held anywhere, multiplies the mean by \( e^{\beta_j} \), the
**rate ratio** of that regressor. It does not depend on the levels of the others,
which is the log-scale version of "no interaction".

::: {.remark}
[The identity link]

Nothing forces the log link. The model \( \mu_i=\x_{(i)}\T\bbeta \) gives
regressors additive effects on the count, which is sometimes the question asked:
an intervention that prevents three events a year, not one that cuts events by ten
per cent. The price is a constrained parameter space, since
\( \x_{(i)}\T\bbeta>0 \) is needed for every \( i \) and the constraint can bind at
the maximum (@exr-cnt-identity-link).
:::

## The likelihood equations

::: {#thm-cnt-score}
[The Poisson likelihood equations]

In the model of @def-cnt-poisson the log-likelihood is
\[
\ell(\bbeta)=\sum_{i=1}^{n}\bigl\{y_i\,\x_{(i)}\T\bbeta-\exp(\x_{(i)}\T\bbeta)\bigr\}
-\sum_{i=1}^{n}\log y_i! ,
\]{#eq-cnt-loglik}

and, writing \( \bmu=\bmu(\bbeta) \) for the vector of means and
\( \W=\diag(\mu_1,\dots,\mu_n) \):

::: {.enumerate options="label=(\alph*)"}
1. the score is \( \partial\ell/\partial\bbeta=\X\T(\y-\bmu) \), so the likelihood
   equations are \( \X\T(\y-\hat{\bmu})=\bzero \);

2. the observed and expected information coincide and equal
   \( \X\T\W\X \);

3. \( \ell \) is concave, and strictly concave when \( \rank(\X)=p \), so a
   stationary point is the unique maximum.
:::

:::

::: {.proof}
The Poisson has \( \phi=1 \), \( w_i=1 \), \( V(\mu)=\mu \) and, under the canonical
link, \( \partial\mu_i/\partial\eta_i=\mu_i \). Substituting these in
@thm-glm-score(d) gives the score \( \X\T(\y-\bmu) \) and observed information equal
to expected information \( \X\T\diag\{w_iV(\mu_i)\}\X=\X\T\W\X \), which are (a) and
(b). For (c), that Hessian is nonnegative definite because
\( \mathbf{a}\T\X\T\W\X\mathbf{a}=\sum_i\mu_i(\x_{(i)}\T\mathbf{a})^2\ge0 \), and the
sum vanishes only if \( \X\mathbf{a}=\bzero \), since every \( \mu_i>0 \). Under full
column rank that forces \( \mathbf{a}=\bzero \), so \( \ell \) is strictly concave.
:::

The general score of @thm-glm-score carries a factor
\( (\partial\mu_i/\partial\eta_i)/V(\mu_i) \), here \( \mu_i/\mu_i=1 \); what is
left has the form of the normal equations \( \X\T(\y-\hY)=\bzero \), with the
fitted mean vector in place of the projection \( \M\y \). The residual
\( \y-\hat{\bmu} \) is orthogonal to every column of \( \X \), although
\( \hat{\bmu} \) is not a linear function of \( \y \) and there is no projection in
sight (@exr-lm-poisson-score).

::: {#thm-cnt-marginals}
[Fitted values reproduce observed margins]

Let \( \hat{\bmu} \) solve the likelihood equations of @thm-cnt-score. Then
\( \x_j\T\hat{\bmu}=\x_j\T\y \) for every column \( \x_j \) of \( \X \), and more
generally \( \mathbf{a}\T\hat{\bmu}=\mathbf{a}\T\y \) for every
\( \mathbf{a}\in\C(\X) \). In particular:

::: {.enumerate options="label=(\alph*)"}
1. if the model has an intercept, the fitted counts sum to the observed total,
   \( \sum_i\hat\mu_i=\sum_iy_i \);

2. if \( \X \) contains the indicator of a group \( G \), the fitted counts of that
   group sum to its observed total, \( \sum_{i\in G}\hat\mu_i=\sum_{i\in G}y_i \);

3. if \( \X \) contains a quantitative regressor \( \x_j \), the fitted counts
   reproduce the observed \( \x_j \)-weighted total.
:::

:::

::: {.proof}
The column identities are @cor-glm-marginals with \( w_i=1 \): the likelihood equations
say \( \X\T(\y-\hat{\bmu})=\bzero \), that is, \( \x_j\T(\y-\hat{\bmu})=0 \) for
each \( j \). The extension to \( \C(\X) \) is linearity: any
\( \mathbf{a}\in\C(\X) \) is \( \X\mathbf{c} \) for some \( \mathbf{c} \), and
\( \mathbf{a}\T(\y-\hat{\bmu})=\mathbf{c}\T\X\T(\y-\hat{\bmu})=0 \). Parts (a), (b)
and (c) are this statement for \( \mathbf{a}=\bone \), for the group indicator, and for
\( \x_j \).
:::

This is what makes log-linear models for contingency tables behave so tidily;
[Section 37.3](03-poisson-multinomial.html) uses it throughout.

## When the maximum likelihood estimate exists

Strict concavity says that a maximum, if attained, is unique; not that one is
attained. The log-likelihood @eq-cnt-loglik can increase along a ray for ever, and
then the fitting algorithm diverges. The condition below is the Poisson
counterpart of separation in logistic regression (@thm-bin-separation).

::: {#prp-cnt-existence}
[Existence of the maximum likelihood estimate]

Let \( \rank(\X)=p \). The log-likelihood @eq-cnt-loglik attains its maximum at a
finite \( \hbeta \) if and only if there is **no** vector \( \mathbf{d}\ne\bzero \)
with
\[
\x_{(i)}\T\mathbf{d}\le0\ \text{ for every } i,
\qquad\text{and}\qquad
\x_{(i)}\T\mathbf{d}=0\ \text{ whenever } y_i>0 .
\]{#eq-cnt-divergent}

:::

::: {.proof}
Write \( u_i=\x_{(i)}\T\mathbf{d} \) and, for a fixed \( \bbeta_0 \),
\[
\ell(\bbeta_0+t\mathbf{d})
=\sum_i y_i\eta_i^{0}+t\sum_i y_iu_i
-\sum_i e^{\eta_i^{0}}e^{tu_i}-\sum_i\log y_i! ,
\]
where \( \eta_i^{0}=\x_{(i)}\T\bbeta_0 \).
Suppose such a \( \mathbf{d} \) exists. Then \( \sum_iy_iu_i=0 \), because every
term with \( y_i>0 \) has \( u_i=0 \). Since every \( u_i\le0 \), each
\( e^{tu_i} \) is nonincreasing in \( t \), so \( \ell(\bbeta_0+t\mathbf{d}) \) is
nondecreasing in \( t \); and it is strictly increasing unless every \( u_i=0 \),
which full column rank excludes. So the supremum is approached only as
\( t\to\infty \) and is not attained.

Conversely, suppose no such \( \mathbf{d} \) exists, and fix a unit vector
\( \mathbf{d} \). Either some \( u_i>0 \), and
\( \ell(\bbeta_0+t\mathbf{d})\to-\infty \) because an exponential beats a linear
term; or every \( u_i\le0 \), in which case some \( i \) with \( y_i>0 \) has
\( u_i<0 \), so \( \sum_iy_iu_i<0 \) while \( -\sum_ie^{\eta_i^0}e^{tu_i}\le0 \),
and again \( \ell\to-\infty \). So there is a \( t_{\mathbf{d}} \) beyond which
\( \ell(\bbeta_0+t\mathbf{d})<\ell(\bbeta_0)-1 \). Continuity alone does not make
that threshold uniform in \( \mathbf{d} \); concavity does. By part (c) the
difference quotient \( \{\ell(\bbeta_0+t\mathbf{d})-\ell(\bbeta_0)\}/t \) is
nonincreasing in \( t \), so once it falls below \( -1/t_{\mathbf{d}} \) it stays
below, and by continuity it does so on a neighbourhood of \( \mathbf{d} \).
Finitely many such neighbourhoods cover the unit sphere, giving a single
\( t^{*} \) outside whose ball \( \ell<\ell(\bbeta_0)-1 \); and a continuous
function on that closed ball attains its maximum.
:::

The condition fails when a regressor is the indicator of a set of observations
whose counts are all zero: take \( \mathbf{d} \) to be minus that coordinate
direction. The coefficient runs to \( -\infty \) and the algorithm reports a huge
negative number with an enormous standard error. A ridge or lasso penalty
([Chapter 27](../ch27-shrinkage/index.html),
[Chapter 30](../ch30-regularization-boosting/index.html)) makes the penalized
log-likelihood coercive and restores a finite maximizer.

## Fitting

By @thm-glm-irls, Fisher scoring is iteratively reweighted least squares: at the
current \( \bbeta \) form the working response and weights
\[
z_i=\eta_i+\frac{y_i-\mu_i}{\mu_i},
\qquad
w_i=\mu_i ,
\]{#eq-cnt-working}

and set \( \bbeta\leftarrow(\X\T\W\X)^{-1}\X\T\W\bz \). Because the expected and
observed information agree (@thm-cnt-score(b)), Fisher scoring and Newton–Raphson
coincide here, and by concavity the iteration converges from any start whenever a
maximum exists. The same weighted least squares step appears in @prp-res-irls and
@thm-het-wls; only the weights differ.
[Figure 37.1.1](#fig-cnt-scoring) draws the log-likelihood surface of a
two-regressor fit with the path from a poor start: the contours are closed and
nested, as concavity requires, and the first step already lands near the
maximum.

::: {when-format="html"}
![**Figure 37.1.1.** The log-likelihood of a Poisson log-linear model with an
intercept and one regressor, for the physician-visit data used below in this
section, with the Fisher scoring path from the
start \( \bbeta=(\log\bar y,0)\T \). Contours are drawn at fixed drops below the
maximum.](scoring_path.svg){#fig-cnt-scoring width=62%}
:::

::: {when-format="pdf"}
![The log-likelihood of a Poisson log-linear model with an
intercept and one regressor, for the physician-visit data used below in this
section, with the Fisher scoring path from the
start \( \bbeta=(\log\bar y,0)\T \). Contours are drawn at fixed drops below the
maximum.](scoring_path.pdf){width=62%}
:::

## Interpretation and inference

The estimate \( \hbeta \) is consistent and asymptotically normal under the
conditions of @thm-glm-asymptotics, with asymptotic covariance
\( (\X\T\W\X)^{-1} \) at the estimate. Exponentiating the ends of the Wald interval
\( \hat\beta_j\pm z_{1-\alpha/2}\,\mathrm{se}(\hat\beta_j) \) gives an interval for
the rate ratio \( e^{\beta_j} \), better than a delta method applied to
\( e^{\hat\beta_j} \): the log scale is where the normal approximation is accurate,
and the transformed interval stays positive.

::: {#exm-cnt-visits}
[Physician visits in a health insurance experiment]

The RAND Health Insurance Experiment sample distributed with `statsmodels`
records \( 20190 \) person-years. The response is the number of outpatient
physician visits in the year, with mean \( 2.860 \) and a variance
\( 7.093 \) times as large; @exm-lm-counts used the same data to show
how badly a normal linear model describes it.

Fit @eq-cnt-model with an intercept, the logarithm of one plus the coinsurance
rate (`lncoins`), an indicator of an individual-deductible plan (`idp`), a
physical-limitation index, a chronic-disease score, and three indicators of
self-rated health (good, fair, poor, against excellent). Fisher scoring converges
in \( 7 \) steps. The coefficient of `lncoins` is
\( -0.0712 \), standard error \( 0.0022 \): a rate ratio of
\( 0.9313 \) per unit, Wald interval
\( (0.9272,\,0.9354) \). Moving across its whole range of
\( 4.615 \), from free care to the highest coinsurance plan,
multiplies expected visits by \( 0.720 \). The disease score is
not a count of diseases but a continuous index running to
\( 58.6 \); each unit of it multiplies expected visits by
\( 1.0353 \), so moving from its lower quartile
\( 6.90 \) to its upper quartile \( 13.73 \)
multiplies them by \( 1.267 \). Poor self-rated health multiplies
them by \( 1.2441 \) against excellent health. The fitted counts
sum to the observed total, and within each health group to that group's total, as
@thm-cnt-marginals requires.
:::

```{.python .run #cell-poisson-visits-fit}
import numpy as np
import statsmodels.api as sm

data = sm.datasets.randhie.load_pandas().data        # 20190 person-years, public domain
y = data["mdvis"].to_numpy(float)                    # outpatient physician visits in the year
names = ["intercept", "lncoins", "idp", "physlm", "disea", "hlthg", "hlthf", "hlthp"]
X = np.column_stack([np.ones(len(y))] + [data[v].to_numpy(float) for v in names[1:]])


def poisson_fit(X, y, offset=0.0, tol=1e-11, maxit=50):
    """Fisher scoring for log mu = X beta + offset; returns beta and its iterates."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    path = [beta]
    for _ in range(maxit):
        mu = np.exp(X @ beta + offset)               # current fitted means
        z = X @ beta + (y - mu) / mu                 # working response
        XW = X * mu[:, None]                         # weights w_i = mu_i
        step = np.linalg.solve(X.T @ XW, XW.T @ z)
        converged = np.max(np.abs(step - beta)) < tol
        beta = step
        path.append(beta)
        if converged:
            break
    return beta, np.array(path)


beta, path = poisson_fit(X, y)
mu = np.exp(X @ beta)
se = np.sqrt(np.diag(np.linalg.inv(X.T @ (X * mu[:, None]))))

print("largest entry of X^T (y - mu):", np.abs(X.T @ (y - mu)).max())
for name, b, s in zip(names, beta, se):
    print(f"{name:10s} {b:8.4f} ({s:.4f})   rate ratio {np.exp(b):.4f}")
```

The deviance of the Poisson fit (@def-glm-deviance) is
\[
D=2\sum_{i=1}^{n}\Bigl\{y_i\log\frac{y_i}{\hat\mu_i}-(y_i-\hat\mu_i)\Bigr\},
\]{#eq-cnt-deviance}

with the convention \( 0\log 0=0 \). If the model has an intercept the second term
sums to zero by @thm-cnt-marginals(a), and the deviance reduces to
\( 2\sum_iy_i\log(y_i/\hat\mu_i) \), the form usually written \( G^2 \); because
\( \phi=1 \) the deviance and the scaled deviance coincide.

Differences of deviances between nested models behave in the usual way: by
@thm-glm-deviance the drop is asymptotically \( \chi^2 \) on the difference in
dimension, and it is the likelihood ratio statistic of @prp-glm-three-tests. The
deviance of a *single* model is a different matter: its \( \chi^2 \) calibration
needs the fitted means to grow, which happens for grouped data with large cells
and not for ungrouped records with small counts. In @exm-cnt-visits the deviance
is \( 84569.7 \) on \( 20182 \) degrees of freedom
and the Pearson statistic is \( 127526.5 \), ratios of
\( 4.190 \) and \( 6.319 \). No formal test is
available, but ratios of this size are a loud signal of misfit.

```{.python .run #cell-poisson-visits-deviance}
dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))
pearson = np.sum((y - mu) ** 2 / mu)
df = len(y) - X.shape[1]
print(f"deviance {dev:.1f} on {df} df   (ratio {dev / df:.3f})")
print(f"Pearson  {pearson:.1f} on {df} df   (ratio {pearson / df:.3f})")
print("fitted total", mu.sum(), " observed total", y.sum())
```

Both statistics are sums of squared residuals: the Pearson residual
\( r^{P}_i=(y_i-\hat\mu_i)/\sqrt{\hat\mu_i} \) and the deviance residual
\( r^{D}_i=\operatorname{sign}(y_i-\hat\mu_i)\sqrt{d_i} \), with \( d_i \) the
\( i \)th term of @eq-cnt-deviance. These are the count analogues of
@def-res-residuals; both are badly skewed when \( \hat\mu_i \) is
small (@exr-cnt-anscombe), so a normal probability plot of them tells one little.
Leverage and influence transfer from
[Chapter 20](../ch20-residuals-leverage-influence/index.html) through the weights:
the hat matrix of the last weighted least squares step of @prp-res-irls plays the
role of \( \M \), and a one-step Cook's distance (@def-res-cooks) built from it is
the standard influence measure for a generalized linear model.

::: {.idea}
The Poisson log-linear model puts the linear predictor on the log scale, so its
coefficients are log rate ratios, and it ties the variance to the mean, so it has
no spare parameter for extra variability. The rest of the chapter follows from
taking one or other of those decisions seriously. A deviance far larger than its
degrees of freedom has three possible causes — a wrong mean model, a wrong
variance function, or dependent observations — and reaching for an overdispersion
parameter before checking the mean model is a common and expensive mistake,
because an inflated variance hides exactly the lack of fit that would have told
you what to add.
:::

## Exercises

### A. Check your understanding

::: {#exr-cnt-units}
[A1]

A Poisson log-linear model for the number of insurance claims has coefficient
\( 0.182 \) on an indicator of urban residence. State the effect on the expected
number of claims in words, and give the multiplicative effect of the indicator on
the *odds* that a policy has at least one claim under the same model. Why is the
second quantity not simply \( e^{0.182} \)?
:::

::: {#exr-cnt-intercept-only}
[A2]

Show directly from @thm-cnt-score that in the model with an intercept only, the
maximum likelihood estimate of the common mean is \( \bar y \), and that the
deviance equals \( 2\sum_iy_i\log(y_i/\bar y) \).
:::

::: {.solution}
With \( \X=\bone \) the likelihood equation is \( \sum_i(y_i-e^{\beta})=0 \), so
\( e^{\hat\beta}=\bar y \). Substituting into @eq-cnt-deviance and using
\( \sum_i(y_i-\bar y)=0 \) leaves \( 2\sum_iy_i\log(y_i/\bar y) \).
:::

### B. Practice

::: {#exr-cnt-one-way}
[B1]

Counts \( y_{ij} \), \( j=1,\dots,n_i \), are observed in \( c \) groups, and the
model is \( \log\mu_{ij}=\beta_i \), one parameter per group. Show that
\( \hat\mu_{ij}=\bar y_{i\cdot} \), that
\( \widehat{\Var}(\hat\beta_i)=1/(n_i\bar y_{i\cdot}) \), and that a Wald interval
for \( \mu_h/\mu_i \) is
\( \exp\bigl\{\hat\beta_h-\hat\beta_i\pm z\sqrt{1/(n_h\bar y_{h\cdot})+1/(n_i\bar y_{i\cdot})}\bigr\} \).
:::

::: {.solution}
@thm-cnt-marginals(b) applied to the group indicators gives
\( n_i\hat\mu_i=\sum_jy_{ij} \). The information \( \X\T\W\X \) is diagonal with
entries \( n_i\mu_i \), so the estimates are independent with variances
\( 1/(n_i\mu_i) \); the variance of their difference is the sum, and exponentiating
a Wald interval for \( \beta_h-\beta_i \) gives the stated interval.
:::

::: {#exr-cnt-one-way-lrt}
[B2]

In the setting of @exr-cnt-one-way, show that the likelihood ratio statistic for
\( H_0:\mu_1=\dots=\mu_c \) is
\( 2\sum_in_i\bar y_{i\cdot}\log(\bar y_{i\cdot}/\bar y) \), where \( \bar y \) is
the overall mean, and state the conditions under which its \( \chi^{2}(c-1) \)
calibration is reasonable.
:::

::: {#exr-cnt-anscombe}
[B3]

For a Poisson fit the Pearson residual is
\( r_i^{P}=(y_i-\hat\mu_i)/\sqrt{\hat\mu_i} \). Show that it is strongly skewed
when \( \mu_i \) is small, and that the Anscombe residual
\( \tfrac32(y_i^{2/3}-\hat\mu_i^{2/3})/\hat\mu_i^{1/6} \) has variance about one,
by applying the delta method to \( Y^{2/3} \).
:::

### C. Going deeper

::: {#exr-cnt-identity-link}
[C1]

Consider a Poisson model with the identity link, \( \mu_i=\x_{(i)}\T\bbeta \).
Show that the log-likelihood is concave on the open set where every
\( \mu_i>0 \), but that the maximum can occur on its boundary, where some
\( \hat\mu_i=0 \); what happens then to the standard errors an algorithm that
ignores the constraint reports? Repeat for the additive *rate* model
\( \mu_i=t_i\,\x_{(i)}\T\bbeta \) of [Section 37.2](02-offsets-and-rates.html),
whose coefficients are rate *differences*: show that its likelihood equations are
\( \sum_it_i(y_i-\mu_i)\x_{(i)}/\mu_i=\bzero \), and say why an epidemiologist
asking how many events an intervention prevents wants it.
:::

::: {#exr-cnt-existence-check}
[C2]

Give a data set with \( n=4 \), \( p=2 \) for which the condition of
@prp-cnt-existence fails. Then add a ridge penalty \( \lambda\norm{\bbeta}^2 \) (@def-shr-ridge)
and show that the penalized log-likelihood has a finite maximizer
for every \( \lambda>0 \).
:::

::: {.solution}
Take \( \X \) with intercept and a regressor \( x=(0,0,1,1)\T \) and counts
\( \y=(3,5,0,0)\T \). Then \( \mathbf{d}=(0,-1)\T \) has
\( \x_{(i)}\T\mathbf{d}\le0 \) everywhere and \( =0 \) at the two observations with
\( y_i>0 \), so the condition fails: the slope runs to \( -\infty \) while the
intercept converges to \( \log4 \). Adding \( -\lambda\norm{\bbeta}^2 \) makes the
objective tend to \( -\infty \) in every direction, since the log-likelihood grows
at most linearly along any ray on which it grows at all, so a finite maximizer
exists; it is unique by strict concavity of the penalized objective.
:::
