# A panorama of regression models to come

The linear model of this chapter makes three kinds of claim: about the *mean* of the response
(a linear combination of known columns), about its *variability* (a common variance, no
correlation) and, in the normal version, about its *distribution*. Every model in the rest of
the book relaxes one or more of these claims, or keeps them and asks harder questions about
them. This section surveys the route.

## Three components of a regression model

It helps to see any regression model as three separate choices.

1. **A distribution for the response given the regressors**: normal, binomial, Poisson,
   gamma, or something less specified, such as "any distribution with this mean and variance".

2. **A predictor** that combines the regressors: the linear predictor
   \( \eta_i=\x_{(i)}\T\bbeta \), or a sum of smooth functions, or a sum that includes random
   terms.

3. **A link** between the predictor and the mean, \( g\bigl(\E(Y_i)\bigr)=\eta_i \), where \( g \) is
   a known monotone function.

The normal linear model chooses the normal distribution with constant variance, the linear
predictor, and the identity link \( g(\mu)=\mu \). Changing one choice at a time generates
the models below.

## The ladder

| Model | What changes | Where |
|---|---|---|
| Linear model, second moments | the starting point: geometry, optimality, estimability, computation | Chapters 5–10 |
| Normal linear model | exact tests, intervals, multiple comparisons, designed experiments | Chapters 11–18 |
| Departures from the assumptions | diagnostics, transformations, resampling, measurement error, causation | Chapters 19–25 |
| Many regressors | collinearity, shrinkage, \( p>n \), selection, boosting | Chapters 26–30 |
| General covariance and random effects | \( \Cov(\be)=\sigma^2\V \); \( \Y=\X\bbeta+\bm Z\bm b+\be \) | Chapters 31–33 |
| Generalized linear models | non-normal distribution and a link | Chapters 34–41 |
| Nonparametric and additive models | smooth functions in the predictor | Chapters 42–44 |
| Distributional regression | regressors act on more than the mean | Chapter 45 |

**The normal linear model and its designs.** With normal errors, the least squares
estimator has an exact normal distribution and the residual sum of squares an exact
chi-squared distribution (@thm-opt-sampling). This gives \( t \) and \( F \) tests and exact
confidence regions (Chapters 11 and 12), and procedures for many comparisons at once
(Chapter 13). Analysis of variance and designed experiments (Chapters 15–18) are linear models
whose columns are indicator variables, often not of full rank; their theory rests on
[Chapter 6](../ch06-projections/index.html) and [Chapter 8](../ch08-estimability/index.html).

**When the assumptions fail, or there are many regressors.** Part V asks what happens when
the mean is misspecified, the variance is not constant, the errors are correlated or not
normal, the regressors are measured with error, or the question is causal; some remedies,
such as robust standard errors based on @eq-lm-sandwich (Chapter 21), stay inside the linear
model. Part VI treats large \( p \) and nearly dependent columns, where ridge regression,
principal components and the lasso trade a little bias for less variance. Regression trees,
random forests and neural networks also estimate \( m \), judged by prediction error rather than
by the meaning of coefficients; the book does not develop them, but Chapters 29 and 30 compare
predictors and treat boosting.

**General covariance and mixed models.** The general Gauss–Markov model replaces
\( \sigma^2\I \) by \( \sigma^2\V \) (Chapter 31), and generalized least squares is least squares
in a different inner product ([Section 6.9](../ch06-projections/09-inner-products.html)). Linear
mixed models write the response as
\[
\Y=\X\bbeta+\bm Z\bm b+\be,
\]
where \( \bm b \) is a vector of random effects shared by related cases, such as repeated
measurements on one person or pupils in one school (Chapters 32–33). The RAND data of
@exm-lm-health-coding, with the same people observed in several years, are of this kind.

**Generalized linear models.** For a binary response the mean is a probability, and a
linear predictor can leave \( [0,1] \). For a count the mean must be positive and the variance
grows with the mean. Generalized linear models keep the linear predictor, choose a
distribution from an exponential family, and connect the two by a link: the logit
\( \log\{\mu/(1-\mu)\} \) for binary responses (Chapter 35), the logarithm for counts
(Chapter 37). Coefficients act additively on the scale of the link, so in a Poisson model
with log link \( e^{\beta_j} \) is the factor by which the mean changes per unit of \( x_j \).
Chapters 34–41 develop the theory, overdispersion, correlated responses and missing data.

**Nonparametric and additive models.** Polynomials, as in @exm-lm-co2, force a global
shape on a curve. Splines and smoothers let the data choose the shape locally
(Chapters 42–43), and additive models replace each term \( \beta_jx_j \) by an unknown smooth
function \( f_j(x_j) \) (Chapter 44). Much of this is still penalized least squares on a model
matrix of basis functions.

**Beyond the mean.** Quantile regression models a chosen quantile of \( Y \), such as the
median or the tenth percentile, as a linear function of the regressors, and so describes how
the whole conditional distribution shifts. Distributional regression lets the regressors act
on the location, the scale and the shape of the response distribution through separate
predictors, turning process (b) of [Figure 5.1.1](01-what-a-model-claims.html#fig-lm-conditional)
into a model rather than a violation (Chapter 45).

**The Bayesian thread.** A prior distribution on \( (\bbeta,\sigma^2) \) turns the likelihood into a
posterior. For the normal linear model with a conjugate prior the posterior is available in
closed form (@thm-opt-bayes-conjugate), and it reappears for intervals (Chapter 12) and for
generalized linear models (Chapter 39).

## A first look beyond the linear model

::: {#exm-lm-counts}
[Physician visits are counts]

The number of physician visits in the RAND data is a count. It is zero for
\( 31.2 \) percent of person-years, and its variance is \( 7.09 \)
times its mean. Group the person-years into ten groups by the chronic-disease score, and
compute the mean and variance of visits in each ([Figure 5.8.1](#fig-lm-counts)). The
variance rises with the mean, and in every group it is between \( 4.8 \) and
\( 8.8 \) times the mean. (L2) fails badly.

Two models for the mean are fitted to the disease score \( x \). The linear model gives
\( 1.268+0.1416\,x \) visits. The Poisson regression with log link
(Chapter 37) gives \( \exp(0.549+0.0408\,x) \): each unit of the index multiplies the
mean by \( e^{0.0408} \), \( 4.17 \) percent more visits. Over the bulk of the data the
two curves are close, and both follow the group means. They differ in shape at high values
of the index, where the log-linear mean grows ever faster and the linear one at a constant
rate, and where the data are too sparse to decide between them. They differ more in what they
assume about variability.
The linear model assumes a constant variance; the Poisson model assumes variance equal to the
mean. The data contradict both. Chapter 38 shows how to keep the log-linear mean and let the
variance be a multiple of it.

The Poisson fit is not computed by ordinary least squares, but its estimating equations
have the same form as the normal equations: \( \X\T(\y-\hat{\bmu})=\bzero \), where
\( \hat{\bmu} \) is the vector of fitted means. The residuals are again orthogonal to every
column of \( \X \). The listing checks this.
:::

::: {when-format="html"}
![**Figure 5.8.1.** Annual physician visits in the RAND Health Insurance Experiment,
grouped into ten groups by a chronic-disease score. (a) Group means with a linear model and a
Poisson regression with log link. (b) Group variances against group means; the dotted line is
the Poisson relation variance = mean.](panorama_counts.svg){#fig-lm-counts width=100%}
:::

::: {when-format="pdf"}
![Annual physician visits in the RAND Health Insurance Experiment,
grouped into ten groups by a chronic-disease score. (a) Group means with a linear model and a
Poisson regression with log link. (b) Group variances against group means; the dotted line is
the Poisson relation variance = mean.](panorama_counts.pdf){width=100%}
:::

```{.python .run #cell-panorama-counts-counts}
import numpy as np
import pandas as pd
import statsmodels.api as sm
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
x = rand["disea"].to_numpy()
X = sm.add_constant(x)

bins = pd.qcut(x, 10, duplicates="drop")        # ten groups by disease score
groups = pd.DataFrame({"y": y, "x": x}).groupby(bins, observed=True)
table = groups.agg(x=("x", "mean"), mean=("y", "mean"), var=("y", "var"))
print(table.round(2).to_string(index=False))
print(f"share of zeros {np.mean(y == 0):.3f}; variance/mean overall {y.var() / y.mean():.2f}")

linear = sm.OLS(y, X).fit()
poisson = sm.GLM(y, X, family=sm.families.Poisson()).fit()
print("linear model:    ", np.round(linear.params, 4))
print("Poisson, log link:", np.round(poisson.params, 4))
mu = poisson.fittedvalues
print("X^T (y - mu_hat) =", X.T @ (y - mu))
```

## Exercises

### A. Check your understanding

::: {#exr-lm-which-model}
[A1]

For each response, name the rung of the ladder you would start from, and say which of the
three components (distribution, predictor, link) departs from the normal linear model:
(a) whether a loan is repaid; (b) blood pressure measured monthly on the same patients;
(c) the number of insurance claims per policy; (d) the \( 90 \)th percentile of birth weight
as a function of the mother's age; (e) crop yield in a randomized block experiment.
:::

### B. Practice

::: {#exr-lm-linear-probability}
[B1]

Let \( Y_i \) be binary with \( \Pr(Y_i=1)=\x_{(i)}\T\bbeta \) (the *linear probability model*).
Show that (L1) holds but (L2) fails, with \( \Var(Y_i)=\mu_i(1-\mu_i) \). Show that least squares is
still unbiased, find its covariance matrix from @prp-lm-misspecified(b), and explain why fitted
values outside \( [0,1] \) are possible.
:::

::: {.solution}
\( \E(Y_i)=\Pr(Y_i=1)=\mu_i=\x_{(i)}\T\bbeta \), which is (L1). \( \Var(Y_i)=\mu_i-\mu_i^2 \), which
varies with \( i \) unless all the \( \mu_i \) are equal or symmetric about \( 1/2 \), so (L2) fails.
Unbiasedness needs only (L1) (@thm-lm-moments(a)). With independent cases,
\( \bSigma=\diag\bigl(\mu_i(1-\mu_i)\bigr) \) and @eq-lm-sandwich gives
\( \Cov(\hbeta)=(\X\T\X)^{-1}\X\T\bSigma\X(\X\T\X)^{-1} \). Nothing in least squares constrains
\( \x\T\hbeta \) to \( [0,1] \): the estimate is a linear function of \( \y \), and at design points far
from the centre of the data the fitted line can cross \( 0 \) or \( 1 \) even when every true
probability is inside. The logit link of Chapter 35 removes the problem.
:::

::: {#exr-lm-poisson-score}
[B2]

For independent \( Y_i\sim\text{Poisson}(\mu_i) \) with \( \log\mu_i=\x_{(i)}\T\bbeta \), write the
log-likelihood and show that its gradient is \( \X\T(\y-\bmu) \). Deduce that when \( \X \) contains an
intercept the fitted means sum to \( \sum_iy_i \), and that for the intercept-only model the
fitted mean is \( \bar{y} \).
:::

### C. Going deeper

::: {#exr-lm-scale-model}
[C1]

In process (b) of @exm-lm-two-processes the error standard deviation is \( 0.15+0.25x \). Suppose
instead \( \log\Var(Y_i)=\gamma_0+\gamma_1x_i \) with unknown \( \gamma \). Propose a two-step
procedure: fit the mean by least squares, then regress \( \log\hat{\varepsilon}_i^2 \) on \( x_i \). Show that,
when the errors are normal and \( n\to\infty \) with \( p \) fixed, the slope of the second
regression estimates \( \gamma_1 \) consistently, while its intercept estimates
\( \gamma_0+\E(\log\chi^2_1) \). (Argue first with the errors \( \varepsilon_i \) in place of the
residuals, then explain why the difference does not matter in the limit.) Compute this bias,
\( -\gamma_{\mathrm E}-\log2\approx-1.27 \), where \( \gamma_{\mathrm E} \) is Euler's constant. (Chapter 45 treats models of this
kind by likelihood.)
:::
