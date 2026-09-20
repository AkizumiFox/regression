# Estimated covariance: feasible generalized least squares

Nothing in the previous three sections is usable until \( \V \) is known, and \( \V \)
is almost never known. The practical procedure is to choose a small family
\( \V(\boldsymbol{\uptheta}) \), estimate \( \boldsymbol{\uptheta} \), and proceed as if the estimate were the
truth. That costs asymptotically nothing and, in small samples, a great deal —
the damage falling on the standard errors rather than on the estimate.

## The two-step estimator

::: {#def-ggm-feasible}
[Feasible generalized least squares]

Let \( \{\V(\boldsymbol{\uptheta}):\boldsymbol{\uptheta}\in\Theta\} \) be a family of positive definite \( n\times n \)
matrices, \( \Theta\subseteq\Real^k \), and let \( \hat{\boldsymbol{\uptheta}} \) be an estimator of
\( \boldsymbol{\uptheta} \) computed from the data. The **feasible generalized least squares**
estimator is
\[
\hbeta_{F}=\bigl(\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X\bigr)^{-1}\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\Y ,
\]{#eq-ggm-fgls}

for \( \rank(\X)=p \), with reported covariance matrix
\( s_F^2\bigl(\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X\bigr)^{-1} \), where \( s_F^2 \) is
@eq-ggm-sse evaluated at \( \V(\hat{\boldsymbol{\uptheta}}) \).
:::

The families that matter are small. Equicorrelation inside clusters,
\( \V(\rho)=(1-\rho)\I+\rho\Z\Z\T \), and first-order autoregressive
errors (@def-het-ar1) both have \( k=1 \); a variance function
\( \Var(\varepsilon_i)=\sigma^2\exp(\bz_i\T\boldsymbol{\upgamma}) \) has \( k=\dim\boldsymbol{\upgamma} \), the
diagonal case treated in
[Section 21.3](../ch21-nonnormality-heteroscedasticity-serial/03-weighted-least-squares.html);
the linear mixed models of
[Chapter 32](../ch32-linear-mixed-models/index.html) give
\( \V(\boldsymbol{\uptheta})=\Z\G(\boldsymbol{\uptheta})\Z\T+\R(\boldsymbol{\uptheta}) \) with a handful of variance components.
In every case \( k \) stays fixed while \( n \) grows; a family with \( k \) of the
order of \( n \) cannot be estimated at all, which is why "estimate \( \V \) freely"
is not an option.

There are three routes to \( \hat{\boldsymbol{\uptheta}} \). The **moment** route matches a few
moments of the ordinary least squares residuals to the model, as in the analysis
of variance estimate of an intraclass correlation used below. The **likelihood**
route maximizes the profile log-likelihood of @exr-ggm-mle, taken up in
[Chapter 32](../ch32-linear-mixed-models/index.html) with the restricted version
that corrects its downward bias. The **iterative** route alternates between
\( \bbeta \) and \( \boldsymbol{\uptheta} \), as in the Cochrane–Orcutt procedure of
@thm-het-ar1(d). Only consistency is needed below.

One case needs no asymptotics. If \( \C(\V(\boldsymbol{\uptheta})\X)\subseteq\C(\X) \) for every
\( \boldsymbol{\uptheta} \), the best estimator is known even though the covariance is
not (@thm-ggm-ols-blue, @exr-ggm-gls-unknown-theta): in a balanced layout with
random block effects, only the standard errors need \( \hat{\boldsymbol{\uptheta}} \).

## Large-sample equivalence

Write \( \V=\V(\boldsymbol{\uptheta}_0) \) for the true covariance, \( \hbeta_{\V} \) for the
infeasible generalized least squares estimator that uses it, and
\( \mathbf{D}(\boldsymbol{\uptheta})=\V(\boldsymbol{\uptheta})^{-1}-\V(\boldsymbol{\uptheta}_0)^{-1} \). The whole question is whether
replacing \( \boldsymbol{\uptheta}_0 \) by \( \hat{\boldsymbol{\uptheta}} \) disturbs the two quantities
\( n^{-1}\X\T\V^{-1}\X \) and \( n^{-1/2}\X\T\V^{-1}\be \) that determine \( \hbeta_{\V} \).
The first is easy, because it does not involve the errors. The second is not,
because \( \hat{\boldsymbol{\uptheta}} \) is computed from the same \( \be \).

::: {#thm-ggm-feasible}
[Feasible generalized least squares]

Assume @eq-ggm-model with \( \rank(\X)=p \) fixed, and consider a sequence of models
indexed by \( n \). Suppose:

- (G1) there are constants \( 0<c\le C<\infty \) and a neighbourhood \( B \) of
  \( \boldsymbol{\uptheta}_0 \) such that every eigenvalue of \( \V(\boldsymbol{\uptheta}) \) lies in \( [c,C] \) for all
  \( \boldsymbol{\uptheta}\in B \) and all \( n \);
- (G2) \( \norm{\mathbf{D}(\boldsymbol{\uptheta})}_2\le L\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0} \) for \( \boldsymbol{\uptheta}\in B \), all \( n \),
  with \( \norm{\cdot}_2 \) the spectral norm;
- (G3) the rows of \( \X \) satisfy \( \norm{\x_{(i)}}\le K \), and
  \( n^{-1}\X\T\V^{-1}\X\to\A \) for a positive definite \( \A \);
- (G4) \( \hat{\boldsymbol{\uptheta}}\to\boldsymbol{\uptheta}_0 \) in probability.

Then:

::: {.enumerate options="label=(\alph*)"}
1. \( n^{-1}\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X\to\A \) in probability, so the reported
   covariance factor is consistent up to the scale \( \sigma^2 \).

2. If in addition \( \sqrt n(\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}_0)=O_p(1) \) and

   - (G5) \( \displaystyle\sup_{\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0}\le\delta}\ \bigl\lVert n^{-1/2}\X\T\mathbf{D}(\boldsymbol{\uptheta})\be\bigr\rVert=O_p(\delta) \)
     as \( \delta\to0 \), uniformly in \( n \),

   then \( \sqrt n(\hbeta_{F}-\hbeta_{\V})\to\bzero \) in probability. Consequently
   \( \hbeta_F \) has whatever limiting distribution \( \hbeta_{\V} \) has, and in
   particular \( \sqrt n(\hbeta_F-\bbeta)\to\Normal_p(\bzero,\sigma^2\A^{-1}) \)
   whenever \( \sqrt n(\hbeta_{\V}-\bbeta) \) does.
:::
:::

::: {.proof}
(a) On the event \( \hat{\boldsymbol{\uptheta}}\in B \), whose probability tends to one by (G4),
\[
\begin{aligned}
\bigl\lVert n^{-1}\X\T\mathbf{D}(\hat{\boldsymbol{\uptheta}})\X\bigr\rVert_2
&\le n^{-1}\norm{\mathbf{D}(\hat{\boldsymbol{\uptheta}})}_2\norm{\X}_2^2\\
&\le n^{-1}L\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}_0}\cdot nK^2
=LK^2\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}_0},
\end{aligned}
\]
using \( \norm{\X}_2^2\le\norm{\X}_F^2=\sum_i\norm{\x_{(i)}}^2\le nK^2 \). By (G4) this
tends to zero in probability, and adding \( n^{-1}\X\T\V^{-1}\X\to\A \) gives (a).

(b) Write \( \A_n(\boldsymbol{\uptheta})=n^{-1}\X\T\V(\boldsymbol{\uptheta})^{-1}\X \) and
\( \bb_n(\boldsymbol{\uptheta})=n^{-1/2}\X\T\V(\boldsymbol{\uptheta})^{-1}\be \), so that
\[
\begin{aligned}
\sqrt n(\hbeta_F-\bbeta)&=\A_n(\hat{\boldsymbol{\uptheta}})^{-1}\bb_n(\hat{\boldsymbol{\uptheta}}),\\
\sqrt n(\hbeta_{\V}-\bbeta)&=\A_n(\boldsymbol{\uptheta}_0)^{-1}\bb_n(\boldsymbol{\uptheta}_0).
\end{aligned}
\]
By (G1) and (G3), \( \E\norm{\bb_n(\boldsymbol{\uptheta}_0)}^2=\sigma^2\tr\A_n(\boldsymbol{\uptheta}_0)\le\sigma^2pK^2/c \)
is bounded, so \( \bb_n(\boldsymbol{\uptheta}_0)=O_p(1) \). By (G5) with
\( \delta_n=\norm{\hat{\boldsymbol{\uptheta}}-\boldsymbol{\uptheta}_0}=O_p(n^{-1/2}) \),
\[
\bb_n(\hat{\boldsymbol{\uptheta}})-\bb_n(\boldsymbol{\uptheta}_0)=n^{-1/2}\X\T\mathbf{D}(\hat{\boldsymbol{\uptheta}})\be=O_p(n^{-1/2})=o_p(1),
\]
so \( \bb_n(\hat{\boldsymbol{\uptheta}})=\bb_n(\boldsymbol{\uptheta}_0)+o_p(1) \) and both are \( O_p(1) \). By (a),
\( \A_n(\hat{\boldsymbol{\uptheta}})\) and \( \A_n(\boldsymbol{\uptheta}_0) \) both converge in probability to the
invertible \( \A \), and matrix inversion is continuous there, so
\( \A_n(\hat{\boldsymbol{\uptheta}})^{-1}-\A_n(\boldsymbol{\uptheta}_0)^{-1}=o_p(1) \). The difference of the two
displays is therefore
\[
\begin{aligned}
&\bigl\{\A_n(\hat{\boldsymbol{\uptheta}})^{-1}-\A_n(\boldsymbol{\uptheta}_0)^{-1}\bigr\}\bb_n(\hat{\boldsymbol{\uptheta}})\\
&\qquad{}+\A_n(\boldsymbol{\uptheta}_0)^{-1}\bigl\{\bb_n(\hat{\boldsymbol{\uptheta}})-\bb_n(\boldsymbol{\uptheta}_0)\bigr\},
\end{aligned}
\]
a sum of products of \( o_p(1) \) and \( O_p(1) \) terms, hence \( o_p(1) \). The last
claim is Slutsky's theorem.
:::

Condition (G5) is the substance of the theorem, and stating it as a hypothesis
rather than proving it is the honest course: it separates "the weights are
random" from "the weights might as well be fixed", and it is not automatic. The
naive bound \( \norm{n^{-1/2}\X\T\mathbf{D}\be}\le n^{-1/2}LK\sqrt n\,\delta\norm{\be} \) is
of order \( \delta\sqrt n \), far too crude, because it ignores the cancellation
among the \( n \) terms. The next lemma supplies that cancellation for one fixed
\( \boldsymbol{\uptheta} \) in the case that matters most, independent clusters of bounded size.

::: {#lem-ggm-drift}
[The drift term for independent clusters]

Suppose the observations fall into \( G \) independent clusters of size at most
\( \bar m \), that \( \V(\boldsymbol{\uptheta}) \) is block diagonal with one block per cluster for
every \( \boldsymbol{\uptheta} \), and that (G1)–(G3) hold. Then for every fixed \( \boldsymbol{\uptheta}\in B \),
\[
\E\bigl\lVert n^{-1/2}\X\T\mathbf{D}(\boldsymbol{\uptheta})\be\bigr\rVert^2\le \sigma^2\bar m\,K^2L^2C\,\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0}^2 .
\]
In particular \( n^{-1/2}\X\T\mathbf{D}(\boldsymbol{\uptheta})\be=O_p(\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0}) \).
:::

::: {.proof}
Write \( \X_k \), \( \be_k \) and \( \mathbf{D}_k \) for the blocks belonging to cluster \( k \), so
that \( \X\T\mathbf{D}\be=\sum_k\X_k\T\mathbf{D}_k\be_k \). The summands are independent with mean
zero, so their covariance matrices add and
\[
\E\bigl\lVert n^{-1/2}\X\T\mathbf{D}\be\bigr\rVert^2
=\frac1n\sum_k\E\norm{\X_k\T\mathbf{D}_k\be_k}^2
\le\frac1n\sum_k\norm{\X_k}_2^2\norm{\mathbf{D}_k}_2^2\,\E\norm{\be_k}^2 .
\]
Now \( \norm{\X_k}_2^2\le m_kK^2 \), \( \norm{\mathbf{D}_k}_2\le\norm{\mathbf{D}}_2\le L\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0} \)
by (G2), and \( \E\norm{\be_k}^2=\sigma^2\tr\V_{kk}\le\sigma^2m_kC \) by (G1). Since
\( m_k\le\bar m \) and \( \sum_km_k=n \), the sum is at most
\( \sigma^2\bar mK^2L^2C\norm{\boldsymbol{\uptheta}-\boldsymbol{\uptheta}_0}^2 \). Chebyshev's inequality
turns the moment bound into the stochastic order.
:::

::: {.remark}
@lem-ggm-drift gives the right order for each fixed \( \boldsymbol{\uptheta} \); (G5) asks for it
*uniformly* over a shrinking ball, which is what allows \( \boldsymbol{\uptheta} \) to be replaced
by the data-dependent \( \hat{\boldsymbol{\uptheta}} \). That passage is a stochastic
equicontinuity argument we do not carry out: with \( \boldsymbol{\uptheta} \) of fixed dimension,
a chaining bound over the ball converts the pointwise bound into the uniform one
under mild smoothness of \( \boldsymbol{\uptheta}\mapsto\V(\boldsymbol{\uptheta})^{-1} \). Amemiya (1985, ch. 6)
proves the conclusion of @thm-ggm-feasible(b) for several concrete families, and
Carroll and Ruppert (1988, ch. 3) for smooth variance functions; (G5) is verified
there, not assumed away.
:::

## What goes wrong in small samples

@thm-ggm-feasible says that \( \hbeta_F \) and \( \hbeta_{\V} \) are
indistinguishable in large samples. It says nothing about finite \( n \), where
three things go wrong.

**The estimator is no longer linear.** \( \hbeta_F \) is a ratio of random
quantities, so none of the exact theory of @thm-ggm-inference applies. A symmetry
argument rescues unbiasedness when \( \hat{\boldsymbol{\uptheta}} \) is an even function of the
residuals and the errors are symmetric, exactly as in @thm-het-wls(c), whose
argument goes through verbatim for a non-diagonal \( \V(\hat{\boldsymbol{\uptheta}})^{-1} \).

**The variance exceeds the generalized least squares variance.** Estimating
\( \boldsymbol{\uptheta} \) adds noise, with an excess of order \( 1/G \) in the cluster setting;
Kackar and Harville (1984) give the standard expansion of it.

**The reported standard error ignores that noise.** The matrix
\( s_F^2(\X\T\V(\hat{\boldsymbol{\uptheta}})^{-1}\X)^{-1} \) is the covariance \( \hbeta_{\V} \)
*would* have had if \( \hat{\boldsymbol{\uptheta}} \) were the truth: it omits the extra
variability, and is evaluated at a random \( \hat{\boldsymbol{\uptheta}} \) drawn towards values
that make the fit look good. Tests are liberal and intervals too short.

::: {#exm-ggm-fgls-simulation}
[How many clusters are enough?]

Simulate clusters of sizes \( 2, 4, 8 \) repeated \( 1,2,4,8,16,32 \) times, with a
cluster-level regressor \( z \), a within-cluster regressor \( w \), and equicorrelation
\( \rho=0.6 \) inside clusters. Estimate \( \rho \) by the analysis of variance moment
estimator from the ordinary least squares residuals and form \( \hbeta_F \).
[Figure 31.4.1](#fig-ggm-fgls) shows \( 20{,}000 \) replications for each size.

The efficiency story is almost trivial. Ordinary least squares attains about
\( 0.895 \) of the optimal precision for the coefficient of \( z \) — its variance
is \( 1.12 \) times the optimal one — at every number of clusters, while feasible
generalized least squares already reaches \( 0.987 \) with \( 3 \)
clusters and \( 1.000 \) with \( 96 \). Estimating one smooth
parameter costs almost nothing in efficiency, in sharp contrast with
@exm-het-replicates, where a separate variance was estimated in each of ten
groups.

The inference story is the opposite. The nominal \( 95\% \) interval based on
\( s_F^2(\X\T\V(\hat\rho)^{-1}\X)^{-1} \) covers the truth in
\( 0.728 \) of samples with \( 3 \) clusters, \( 0.916 \) with
\( 12 \), and \( 0.944 \) with \( 96 \), while the interval that uses the
true \( \rho \) covers \( 0.950 \) at three clusters and stays within simulation
error of \( 95\% \) at every size, as @thm-ggm-inference guarantees. The
cluster-robust interval built on \( \hbeta_F \) is uniformly a little better than
the model-based one, but no cure when clusters are few — \( 0.764 \) at
\( 3 \) — and the two become indistinguishable once there are several dozen,
neither reaching the nominal level even at \( 96 \), where both stand at
\( 0.944 \).
:::

::: {when-format="html"}
![**Figure 31.4.1.** Feasible generalized least squares with equicorrelated clusters of
sizes \( 2,4,8 \) and \( \rho=0.6 \), for the coefficient of a cluster-level regressor
(20,000 replications each). (a) Efficiency: the variance of generalized least squares
with \( \rho \) known divided by the estimator's own variance. Feasible generalized
least squares recovers essentially all of it at once. (b) Coverage of nominal 95%
intervals: both intervals need dozens of clusters, and the cluster-robust one is no
cure.](feasible_gls.svg){#fig-ggm-fgls width=100%}
:::

::: {when-format="pdf"}
![Feasible generalized least squares with equicorrelated clusters of
sizes \( 2,4,8 \) and \( \rho=0.6 \), for the coefficient of a cluster-level regressor
(20,000 replications each). (a) Efficiency: the variance of generalized least squares
with \( \rho \) known divided by the estimator's own variance. Feasible generalized
least squares recovers essentially all of it at once. (b) Coverage of nominal 95%
intervals: both intervals need dozens of clusters, and the cluster-robust one is no
cure.](feasible_gls.pdf){width=100%}
:::

::: {.idea}
Estimating a covariance parameter buys back the efficiency almost at once and
the correct standard error only slowly: efficiency depends on
\( \hat{\boldsymbol{\uptheta}} \) through a function whose derivative vanishes at the optimum,
so a mediocre \( \hat{\boldsymbol{\uptheta}} \) still gives nearly optimal weights, while the
reported variance depends on \( \hat{\boldsymbol{\uptheta}} \) to first order and inherits all its
error.
:::

## Which repair to use

*Ordinary least squares with robust standard errors* assumes nothing about the
shape of \( \V \), loses the few per cent bounded by @prp-ggm-ols-loss, and needs
many clusters. *Feasible generalized least squares* recovers that efficiency, but
only as far as \( \V(\boldsymbol{\uptheta}) \) is right; a wrong shape leaves the estimate
unbiased — by the argument of @thm-het-wls(b), applied with
\( \mathbf{W}=\V(\boldsymbol{\uptheta}_1)^{-1} \) for the limiting \( \boldsymbol{\uptheta}_1 \) — and the reported
covariance wrong in an unknown direction. *Both*, a fitted \( \V(\hat{\boldsymbol{\uptheta}}) \)
with a robust standard error, is the modern default for clustered data, and
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) develops it for
longitudinal designs.

None of the three helps when the independent units are few: with three clusters
\( \hat\rho \) is hopeless and a sandwich with three terms is not a covariance
estimate. The remedy is a better small-sample approximation — the Satterthwaite
and Kenward–Roger degrees of freedom of
[Chapter 32](../ch32-linear-mixed-models/index.html), or the bias-corrected
sandwiches and cluster bootstrap of
[Chapter 23](../ch23-resampling-inference/index.html) — or more clusters.

::: {#exm-ggm-nino-feasible}
[The El Niño trend without assuming \( \rho \)]

The analysis of @exm-ggm-elnino assumed \( \rho=0.5 \). The moment estimator from
the ordinary least squares residuals gives \( \hat\rho=0.610 \), and the
feasible generalized least squares fit then reports a trend of
\( 0.1346 \) degrees per decade with standard error \( 0.0618 \). For
comparison, the ordinary least squares trend, \( 0.1349 \), carries the usual
standard error \( 0.0223 \) and a cluster-robust standard error of
\( 0.0584 \) with years as clusters. The two honest routes agree closely and
the naive one is out by a factor of nearly three. With \( 61 \) clusters both
are in the range where @exm-ggm-fgls-simulation found them reliable.
:::

```{.python .run #cell-feasible-gls-helpers}
import numpy as np


def quasi_demean(A, Z, sizes, rho):
    """Whitening for equicorrelation rho: rows of A are vectors in R^n."""
    c = 1.0 - np.sqrt((1 - rho) / (1 - rho + sizes * rho))
    return (A - (A @ Z / sizes * c) @ Z.T) / np.sqrt(1 - rho)


def fit(X, Y, Z, sizes, rho):
    """GLS coefficients and the model-based covariance factor, for many data sets at once."""
    Xw = quasi_demean(X.T, Z, sizes, rho).T
    Yw = quasi_demean(Y, Z, sizes, rho)
    A = np.linalg.inv(Xw.T @ Xw)
    B = Yw @ Xw @ A
    R = Yw - B @ Xw.T
    s2 = np.sum(R * R, axis=1) / (X.shape[0] - X.shape[1])
    return B, A, s2


def moment_rho(R, Z, sizes):
    """ANOVA estimate of the intraclass correlation from residual vectors R."""
    G, n = len(sizes), R.shape[1]
    means = R @ Z / sizes
    between = np.sum(means**2 * sizes, axis=1) / G
    within = (np.sum(R * R, axis=1) - np.sum(means**2 * sizes, axis=1)) / (n - G)
    mbar = (n - np.sum(sizes**2) / n) / (G - 1)
    est = (between - within) / (between + (mbar - 1) * within)
    return np.clip(est, 0.0, 0.95)
```

```{.python .run #cell-feasible-gls-elnino}
import numpy as np
import statsmodels.api as sm

frame = sm.datasets.elnino.load_pandas().data
temps = frame.iloc[:, 1:].to_numpy()
G, m = temps.shape
n = G * m
y = temps.ravel()
month = np.tile(np.arange(1, m + 1), G)
time = (np.repeat(frame["YEAR"].to_numpy(), m) - 1980.0) + (month - 6.5) / 12.0
cols = [np.ones(n), time / 10.0]
for k in (1, 2, 3):
    cols += [np.cos(2 * np.pi * k * month / m), np.sin(2 * np.pi * k * month / m)]
X = np.column_stack(cols)
Zn = np.zeros((n, G))
Zn[np.arange(n), np.repeat(np.arange(G), m)] = 1.0
sizes_n = np.full(G, m)

b_ols = np.linalg.solve(X.T @ X, X.T @ y)
rho_hat = float(moment_rho((y - X @ b_ols)[None, :], Zn, sizes_n)[0])
b_f, A_f, s2_f = fit(X, y[None, :], Zn, sizes_n, rho_hat)
se_f = np.sqrt(s2_f[0] * A_f[1, 1])
print(f"El Nino: rho_hat = {rho_hat:.3f}, trend {b_f[0, 1]:.4f} (s.e. {se_f:.4f})")
```

## Exercises

### A. Check your understanding

::: {#exr-ggm-feasible-not-linear}
[A1]

Explain why \( \hbeta_F \) is not a linear function of \( \Y \), and why this alone
makes @cor-opt-aitken inapplicable to it. Which part of @thm-ggm-inference could
still hold exactly if \( \hat{\boldsymbol{\uptheta}} \) happened to be independent of \( \Y \)?
:::

::: {#exr-ggm-known-blue-family}
[A2]

Give a family \( \V(\boldsymbol{\uptheta}) \) with \( k=2 \) for which no estimation of \( \boldsymbol{\uptheta} \) is
needed to obtain the best estimator of \( \X\bbeta \), and say what \( \boldsymbol{\uptheta} \) is still
needed for.
:::

### B. Practice

::: {#exr-ggm-moment-rho}
[B1]

For balanced clusters of size \( m \) with equicorrelation \( \rho \), let \( B \) and \( W \)
be the between- and within-cluster mean squares of the ordinary least squares
residuals. Show that the moment estimator
\( \hat\rho=(B-W)/\{B+(m-1)W\} \) is the solution of matching \( \E B \) and \( \E W \) to
their model values when the model is fitted without regressors, and explain why
\( \hat\rho \) can be negative. What does truncating it at zero do to
@thm-ggm-feasible?
:::

::: {.solution}
For a one-way layout with no regressors, \( \E W=\sigma^2(1-\rho) \) and
\( \E B=\sigma^2\{1+(m-1)\rho\} \) (this is @exr-ss-random-oneway). Solving
\( B=\sigma^2\{1+(m-1)\rho\} \), \( W=\sigma^2(1-\rho) \) for \( \rho \) gives the stated
ratio. Since \( B \) and \( W \) are independent chi-squared multiples, \( B<W \) happens
with positive probability whenever \( \rho \) is small, so \( \hat\rho<0 \) is possible;
[Chapter 32](../ch32-linear-mixed-models/index.html) discusses negative variance-component
estimates at length. Truncation leaves \( \hat\rho \) consistent whenever
\( \rho>0 \), which is all (G4) requires; at \( \rho=0 \) the truncated estimator is
still consistent and still satisfies \( \sqrt n\,\hat\rho=O_p(1) \), so
@thm-ggm-feasible(b) still applies to \( \hbeta_F \), but \( \hat\rho \) itself is no
longer asymptotically normal and the usual Wald theory for it fails.
:::

::: {#exr-ggm-two-stage-weights}
[B2]

Suppose \( \V(\theta)=\diag(1+\theta z_1,\dots,1+\theta z_n) \) with known
\( z_i\ge0 \) and \( \theta\ge0 \). Verify (G1) and (G2) on a neighbourhood of an
interior \( \theta_0 \), giving explicit constants, and show that (G3) holds when the
\( z_i \) and the rows of \( \X \) are bounded.
:::

::: {#exr-ggm-efficiency-flat}
[B3]

In the balanced cluster model of @exr-ggm-two-groups, all clusters of size \( m \),
consider estimating \( \mu \) with weights built from a wrong value \( \rho_1 \) when
the truth is \( \rho_0 \). Show that the resulting estimator is \( \bar Y \) for every
\( \rho_1 \), and hence that the efficiency is exactly \( 1 \). Contrast with the
unbalanced case and explain the flat curve in
[Figure 31.4.1](#fig-ggm-fgls)(a).
:::

::: {.solution}
With equal sizes, @eq-ggm-cluster-weights gives \( w_k=m/(1-\rho_1+m\rho_1) \), the
same for every \( k \), so the weighted mean of the cluster means is \( \bar Y \)
whatever \( \rho_1 \) is: the estimator does not depend on the weights at all. With
unequal sizes the weights do depend on \( \rho_1 \), but only through the smooth
ratio \( m_k/(1-\rho_1+m_k\rho_1) \), and the variance of the resulting estimator is
stationary at \( \rho_1=\rho_0 \) because \( \rho_0 \) minimizes it. A first-order error
in \( \hat\rho \) therefore produces only a second-order loss of efficiency, which
is why panel (a) is flat.
:::

### C. Going deeper

::: {#exr-ggm-drift-sharp}
[C1]

Show that the bound of @lem-ggm-drift cannot be improved in order: construct a
family with \( k=1 \), independent clusters of size \( 1 \), and a design for which
\( \E\norm{n^{-1/2}\X\T\mathbf{D}(\theta)\be}^2 \) is bounded below by a positive multiple
of \( (\theta-\theta_0)^2 \).
:::

::: {#exr-ggm-plugin-bias}
[C2]

Let \( \hat{\boldsymbol{\uptheta}} \) be independent of \( \Y \) with \( \E\hat{\boldsymbol{\uptheta}}=\boldsymbol{\uptheta}_0 \), and write
\( g(\boldsymbol{\uptheta})=\mathbf{a}\T(\X\T\V(\boldsymbol{\uptheta})^{-1}\X)^{-1}\mathbf{a} \) for the reported variance of
\( \mathbf{a}\T\hbeta_F \) at scale \( \sigma^2=1 \), and \( h(\boldsymbol{\uptheta}) \) for its true variance.
Show that \( h(\boldsymbol{\uptheta})\ge g(\boldsymbol{\uptheta}_0) \) for every \( \boldsymbol{\uptheta} \), with equality at
\( \boldsymbol{\uptheta}_0 \), and deduce that \( \E g(\hat{\boldsymbol{\uptheta}}) \) understates
\( \E h(\hat{\boldsymbol{\uptheta}}) \) whenever \( g \) is concave near \( \boldsymbol{\uptheta}_0 \). Why does this
not settle the question in general?
:::

::: {.solution}
For fixed \( \boldsymbol{\uptheta} \) the statistic \( \mathbf{a}\T\hbeta_{\V(\boldsymbol{\uptheta})} \) is a linear unbiased
estimator of \( \mathbf{a}\T\bbeta \), so @cor-opt-aitken gives
\( h(\boldsymbol{\uptheta})\ge h(\boldsymbol{\uptheta}_0)=g(\boldsymbol{\uptheta}_0) \), with equality at \( \boldsymbol{\uptheta}_0 \). Conditionally
on \( \hat{\boldsymbol{\uptheta}} \) the report is \( g(\hat{\boldsymbol{\uptheta}}) \) and the truth is
\( h(\hat{\boldsymbol{\uptheta}})\ge g(\boldsymbol{\uptheta}_0) \); concavity of \( g \) and Jensen's inequality give
\( \E g(\hat{\boldsymbol{\uptheta}})\le g(\boldsymbol{\uptheta}_0) \). In general \( g \) need not be concave and
\( \hat{\boldsymbol{\uptheta}} \) is not independent of \( \Y \); that dependence is what
@exm-ggm-fgls-simulation measures.
:::
