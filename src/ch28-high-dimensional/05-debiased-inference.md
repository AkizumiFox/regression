# Post-selection and debiased inference

This section turns to confidence intervals and tests for a single coefficient \( \beta_j \) when \( p>n \). Two
tempting approaches fail: least squares intervals for the model the lasso selected, as if it had been
fixed in advance, and intervals around the lasso estimate itself, whose shrinkage bias is at least as large
as its standard error. The fix is a one-step correction of the lasso, which turns out to be the
Frisch–Waugh–Lovell theorem with a lasso in place of a projection.

## Why intervals after selection fail

The least squares interval of @thm-ci-estimable-interval has coverage \( 1-\alpha \) for a model fixed before
the data are seen. After selection, the interval for \( \beta_j \) is computed only in the samples where \( x_j \) was
selected, which are those in which its estimated effect happened to be large. Conditional on selection the
estimate is biased away from zero, and the interval is centred in the wrong place. This is the winner's curse of
[Chapter 29](../ch29-model-selection/index.html) (@prp-sel-selection-bias), and the multiplicity problem
of [Chapter 13](../ch13-multiplicity/index.html) in another form. With the lasso and \( p>n \) it is
severe, because the selection is among hundreds or thousands of candidates.

Leeb and Pötscher (2006) proved that this cannot be repaired by estimating the distribution of the
post-selection estimator, which depends on the unknown coefficients in a way that no estimator can track
uniformly. One can change the target (the selected model's coefficients, conditional on the selection)
or the estimator. This section changes the estimator.

The lasso estimate itself is no better as a centre for an interval. On the event of @thm-hd-support,
@eq-hd-oracle-lasso shows that it equals least squares on \( S \) minus \( \lambda\hat{\bSigma}_{SS}^{-1}\operatorname{sign}(\bbeta_S) \).
The shift is of order \( \lambda\asymp\sigma\sqrt{\log p/n} \), which is larger than the standard error \( \sigma/\sqrt n \)
of the oracle estimate by the factor \( \sqrt{\log p} \). An interval of the usual width centred at the lasso estimate
misses the truth more and more often as \( p \) grows.

## The debiased lasso

Write \( \hat{\bSigma}=\X\T\X/n \), and let \( \hat{\boldsymbol{\Theta}} \) be a \( p\times p \) matrix computed from \( \X \) alone, meant to
act like an inverse of \( \hat{\bSigma} \) (which does not exist when \( p>n \)). Its rows are \( \hat{\boldsymbol{\uptheta}}_1\T,\dots,\hat{\boldsymbol{\uptheta}}_p\T \).

::: {#def-hd-debiased}
[Debiased lasso]

Given a lasso estimate \( \hat{\bbeta}_\lambda \) and \( \hat{\boldsymbol{\Theta}} \), the **debiased** (or **desparsified**) lasso is
\[
\hat{\mathbf{b}}=\hat{\bbeta}_\lambda+\frac1n\,\hat{\boldsymbol{\Theta}}\X\T(\Y-\X\hat{\bbeta}_\lambda).
\]{#eq-hd-debiased}

Its \( j \)th entry uses only the row \( \hat{\boldsymbol{\uptheta}}_j \). The quantity
\( \mu_j=\norm{\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j-\mathbf{e}_j}_\infty \) measures how far \( \hat{\boldsymbol{\uptheta}}_j \) is from being the \( j \)th row of
an inverse of \( \hat{\bSigma} \).
:::

The correction adds to the lasso a multiple of the residual correlations \( \X\T(\Y-\X\hat{\bbeta}_\lambda)/n \). If
\( \hat{\bSigma} \) were invertible and \( \hat{\boldsymbol{\Theta}}=\hat{\bSigma}^{-1} \), \( \hat{\mathbf{b}} \) would be the least squares estimate (@exr-hd-debias-low-dim).

::: {#thm-hd-debiased}
[Debiased lasso]

In the setting of [Section 28.3](03-lasso-bounds.html), let \( \hat{\boldsymbol{\Theta}} \) depend on \( \X \) only, and let
\( \hat{\mathbf{b}} \) be the debiased lasso.

::: {.enumerate options="label=(\alph*)"}
1. *(Decomposition.)* \( \sqrt n(\hat{\mathbf{b}}-\bbeta)=\mathbf{w}+\mathbf{r} \), where
   \( \mathbf{w}=\hat{\boldsymbol{\Theta}}\X\T\be/\sqrt n\sim\Normal_p(\bzero,\sigma^2\hat{\boldsymbol{\Theta}}\hat{\bSigma}\hat{\boldsymbol{\Theta}}\T) \) exactly, and
   \( \mathbf{r}=-\sqrt n(\hat{\boldsymbol{\Theta}}\hat{\bSigma}-\I)(\hat{\bbeta}_\lambda-\bbeta) \).

2. *(Remainder.)* \( |r_j|\le\sqrt n\,\mu_j\norm{\hat{\bbeta}_\lambda-\bbeta}_1 \). If \( \lambda\ge2\lambda_0 \), \( \phi^2(S)>0 \) and
   \( \mathcal{T} \) occurs, then \( |r_j|\le4\sqrt n\,\mu_j\lambda s/\phi^2(S) \).

3. *(Variance.)* If \( \mu_j<1 \), then \( \omega_j^2=\hat{\boldsymbol{\uptheta}}_j\T\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j\ge(1-\mu_j)^2 \).

4. *(Asymptotic normality.)* Consider a sequence of such problems with \( n\to\infty \), in which \( p \), \( s \), \( \X \),
   \( \bbeta \) and \( \hat{\boldsymbol{\Theta}} \) may change with \( n \) and \( \sigma \) is fixed. Suppose \( \lambda=2\lambda_0 \) with \( \delta=\delta_n\to0 \)
   and \( \log(1/\delta_n)=O(\log p) \); \( \phi^2(S)\ge\phi_0^2>0 \); \( \mu_j\le C\sqrt{\log p/n} \) for a constant \( C \);
   \( s\log p/\sqrt n\to0 \); and \( \hat\sigma/\sigma\to1 \) in probability. Then
   \[
T_j=\frac{\sqrt n\,(\hat b_j-\beta_j)}{\hat\sigma\,\omega_j}\ \to\ \Normal(0,1)\quad\text{in distribution},
\]{#eq-hd-debiased-t}

   and the interval \( \hat b_j\pm z_{\alpha/2}\hat\sigma\omega_j/\sqrt n \) has coverage tending to \( 1-\alpha \).
:::

:::

::: {.proof}
(a) Substitute \( \Y=\X\bbeta+\be \) in @eq-hd-debiased and write \( \boldsymbol{\Delta}=\hat{\bbeta}_\lambda-\bbeta \):
\[
\hat{\mathbf{b}}-\bbeta=\boldsymbol{\Delta}+\frac1n\hat{\boldsymbol{\Theta}}\X\T(\be-\X\boldsymbol{\Delta})
=\frac1n\hat{\boldsymbol{\Theta}}\X\T\be-(\hat{\boldsymbol{\Theta}}\hat{\bSigma}-\I)\boldsymbol{\Delta} .
\]
Multiply by \( \sqrt n \). Since \( \hat{\boldsymbol{\Theta}} \) is a fixed matrix, \( \mathbf{w} \) is a linear function of \( \be \), hence normal with
mean zero and covariance \( \hat{\boldsymbol{\Theta}}\X\T(\sigma^2\I)\X\hat{\boldsymbol{\Theta}}\T/n=\sigma^2\hat{\boldsymbol{\Theta}}\hat{\bSigma}\hat{\boldsymbol{\Theta}}\T \) (@thm-mvn-linear).

(b) The \( j \)th row of \( \hat{\boldsymbol{\Theta}}\hat{\bSigma}-\I \) is \( (\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j-\mathbf{e}_j)\T \), because \( \hat{\bSigma} \) is symmetric. By
Hölder's inequality \( |r_j|\le\sqrt n\,\mu_j\norm{\boldsymbol{\Delta}}_1 \). The second bound follows from @thm-hd-lasso-estimation(a).

(c) The columns of \( \X \) have \( \norm{\x_j}^2=n \), so \( \mathbf{e}_j\T\hat{\bSigma}\mathbf{e}_j=1 \). The \( j \)th entry of \( \hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j \)
differs from \( 1 \) by at most \( \mu_j \), so \( \mathbf{e}_j\T\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j\ge1-\mu_j>0 \). The Cauchy–Schwarz inequality
for the nonnegative definite form \( \hat{\bSigma} \) (@prp-mat-cauchy-schwarz applied to \( \hat{\bSigma}^{1/2}\mathbf{e}_j \) and
\( \hat{\bSigma}^{1/2}\hat{\boldsymbol{\uptheta}}_j \)) gives \( (1-\mu_j)^2\le(\mathbf{e}_j\T\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j)^2\le(\mathbf{e}_j\T\hat{\bSigma}\mathbf{e}_j)\,\omega_j^2=\omega_j^2 \).

(d) Write \( T_j=(\sigma/\hat\sigma)\{w_j/(\sigma\omega_j)+r_j/(\sigma\omega_j)\} \). By (a), \( w_j/(\sigma\omega_j) \) has exactly the
\( \Normal(0,1) \) law. For the remainder, \( \log(2p/\delta_n)=\log2+\log p+\log(1/\delta_n)\le C'\log p \) for a constant
\( C' \) and large \( n \), so \( \lambda=2\lambda_0\le2\sigma\sqrt{2C'\log p/n} \). On \( \mathcal{T} \), which has probability at least
\( 1-\delta_n\to1 \), part (b) gives
\[
\begin{aligned}
|r_j|&\le4\sqrt n\cdot C\sqrt{\frac{\log p}n}\cdot2\sigma\sqrt{\frac{2C'\log p}n}\cdot\frac{s}{\phi_0^2}\\
&=\frac{8\sqrt{2C'}\,C\sigma}{\phi_0^2}\cdot\frac{s\log p}{\sqrt n}\ \to\ 0 .
\end{aligned}
\]
By (c) and \( \mu_j\to0 \), \( \omega_j\ge\frac12 \) for large \( n \). Hence \( r_j/(\sigma\omega_j)\to0 \) in probability, and
Slutsky's lemma (with \( \sigma/\hat\sigma\to1 \)) gives @eq-hd-debiased-t. The coverage statement is
\( \Pr(|T_j|\le z_{\alpha/2})\to1-\alpha \).
:::

The theorem separates what the design must deliver (a \( \hat{\boldsymbol{\Theta}} \) with small \( \mu_j \), and \( \phi^2(S) \) bounded
below) from what the lasso delivers (an \( \ell_1 \) error of order \( s\sqrt{\log p/n} \)). The condition
\( s\log p/\sqrt n\to0 \) is stronger than the \( s\log p/n\to0 \) that suffices for prediction: the remainder must be
small compared with the standard error \( 1/\sqrt n \).

## Constructing the approximate inverse

The remaining question is how to find a \( \hat{\boldsymbol{\Theta}} \) with small \( \mu_j \) when \( \hat{\bSigma} \) is singular. The
*node-wise lasso* regresses each column on all the others.

::: {#lem-hd-nodewise}
[Node-wise lasso]

For a fixed \( j \), let \( \hat{\bgamma}_j\in\Real^{p-1} \) be a lasso estimate, with penalty \( \lambda_j \), for the regression of
\( \x_j \) on the other columns \( \X_{-j} \). Put \( \hat{\bz}_j=\x_j-\X_{-j}\hat{\bgamma}_j \) and \( \hat\tau_j^2=\x_j\T\hat{\bz}_j/n \), and let
\( \hat{\boldsymbol{\uptheta}}_j \) have entry \( 1/\hat\tau_j^2 \) in position \( j \) and \( -\hat{\bgamma}_j/\hat\tau_j^2 \) elsewhere. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hat\tau_j^2=\norm{\hat{\bz}_j}^2/n+\lambda_j\norm{\hat{\bgamma}_j}_1\ge0 \);

2. if \( \hat\tau_j^2>0 \), then \( (\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j)_j=1 \) and \( \mu_j=\norm{\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j-\mathbf{e}_j}_\infty\le\lambda_j/\hat\tau_j^2 \);

3. the \( j \)th debiased estimate is
   \[
\hat b_j=\frac{\hat{\bz}_j\T(\Y-\X_{-j}\hat{\bbeta}_{\lambda,-j})}{\hat{\bz}_j\T\x_j}.
\]{#eq-hd-debiased-fwl}

:::

:::

::: {.proof}
The optimality conditions @eq-hd-optimality for the node-wise regression say that
\[
\X_{-j}\T\hat{\bz}_j/n=\lambda_j\hat{\boldsymbol{\upkappa}},
\]
where \( \hat{\boldsymbol{\upkappa}} \) has entries of size at most \( 1 \) and equal to
\( \operatorname{sign}(\hat\gamma_{j,k}) \) where \( \hat\gamma_{j,k}\ne0 \). Hence
\( \hat{\bgamma}_j\T\X_{-j}\T\hat{\bz}_j/n=\lambda_j\norm{\hat{\bgamma}_j}_1 \), and since \( \x_j=\hat{\bz}_j+\X_{-j}\hat{\bgamma}_j \),
\[
\hat\tau_j^2=\frac{\x_j\T\hat{\bz}_j}{n}=\frac{\norm{\hat{\bz}_j}^2}n+\lambda_j\norm{\hat{\bgamma}_j}_1 ,
\]
which is (a). For (b), \( \X\hat{\boldsymbol{\uptheta}}_j=(\x_j-\X_{-j}\hat{\bgamma}_j)/\hat\tau_j^2=\hat{\bz}_j/\hat\tau_j^2 \), so the \( k \)th entry of
\( \hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j \) is \( \x_k\T\hat{\bz}_j/(n\hat\tau_j^2) \). For \( k=j \) this is \( 1 \) by the definition of \( \hat\tau_j^2 \). For
\( k\ne j \) it is \( \lambda_j\hat\kappa_k/\hat\tau_j^2 \), of size at most \( \lambda_j/\hat\tau_j^2 \). For (c), \( \hat b_j \) is
\( \hat\beta_{\lambda,j}+\hat{\bz}_j\T(\Y-\X\hat{\bbeta}_\lambda)/(n\hat\tau_j^2) \). Put \( n\hat\tau_j^2=\hat{\bz}_j\T\x_j \) and
\( \X\hat{\bbeta}_\lambda=\x_j\hat\beta_{\lambda,j}+\X_{-j}\hat{\bbeta}_{\lambda,-j} \); the terms in \( \hat\beta_{\lambda,j} \) cancel.
:::

Part (c) shows what the debiased lasso is. By the Frisch–Waugh–Lovell theorem (@thm-proj-fwl), the least
squares coefficient of \( \x_j \) comes from regressing \( \Y \) on the part of \( \x_j \) orthogonal to the other columns,
which is zero when \( p>n \). The node-wise lasso replaces it by a lasso residual \( \hat{\bz}_j \), nearly orthogonal to
every other column by (b), and removes the other regressors' effect approximately through
\( \X_{-j}\hat{\bbeta}_{\lambda,-j} \). Javanmard and Montanari (2014) instead choose \( \hat{\boldsymbol{\uptheta}}_j \) to minimize
\( \hat{\boldsymbol{\uptheta}}\T\hat{\bSigma}\hat{\boldsymbol{\uptheta}} \) subject to \( \norm{\hat{\bSigma}\hat{\boldsymbol{\uptheta}}-\mathbf{e}_j}_\infty\le\mu \).

::: {.remark}
[Sketch: when the assumptions of part (d) hold]

This is a sketch; the full arguments are in van de Geer, Bühlmann, Ritov and Dezeure (2014) and Javanmard and
Montanari (2014). Suppose the rows of \( \X \) are independent \( \Normal_p(\bzero,\bSigma) \), the eigenvalues of \( \bSigma \)
are bounded away from \( 0 \) and \( \infty \), and row \( j \) of \( \bSigma^{-1} \) has at most \( s_j \) nonzero entries. The
population regression of \( x_j \) on the others then has \( s_j \) nonzero coefficients and residual variance
\( \tau_j^2=1/(\bSigma^{-1})_{jj} \). By the lemma \( \mu_j\le\lambda_j/\hat\tau_j^2 \) *deterministically*, so with
\( \lambda_j\asymp\sqrt{\log p/n} \) it remains to show \( \hat\tau_j^2\approx\tau_j^2 \). This follows from
[Section 28.3](03-lasso-bounds.html) applied to the node-wise regression, with the random-design bound of
[Section 28.2](02-restricted-eigenvalues.html), when \( s_j\log p/n\to0 \); the scaled lasso (Sun and Zhang 2012)
supplies \( \hat\sigma \). Then \( \omega_j^2\to(\bSigma^{-1})_{jj} \), which van de Geer and coauthors show to be the
smallest asymptotic variance of a regular estimator.
:::

::: {.warning}
The theorem covers one prespecified coefficient, or a fixed finite set. Coefficients chosen *because* they
look large bring back the selection problem, and simultaneous intervals need a multiplicity correction such
as Bonferroni (@thm-mc-bonferroni). With moderate \( n \) and a truth that is not very sparse, the remainder
\( r_j \) can be comparable to the noise, and coverage suffers.
:::

## A simulation

::: {#exm-hd-debiased-coverage}
[Coverage of three intervals]

Take \( n=200 \), \( p=400 \), a fixed Gaussian design with correlations \( 0.5^{|j-k|} \), \( \sigma=1 \), and
\( \beta_1=1 \), \( \beta_6=0.15 \), \( \beta_{11}=-1 \), \( \beta_{16}=0.15 \), all other coefficients zero. The targets are a strong
coefficient \( \beta_1 \), a weak one \( \beta_6 \), and \( \beta_2=0 \), whose regressor is correlated \( 0.5 \) with \( x_1 \). In
each of \( 1000 \) data sets we fit the lasso with \( \lambda=\sigma\sqrt{2\log p/n}=0.245 \), below \( 2\lambda_0 \) for the
reason seen in @exm-hd-lasso-bounds, and compute three 95% intervals:

- the **naive** interval: least squares on the selected variables and its \( t \) interval, reported only when the
  target was selected;
- the **debiased** interval \( \hat b_j\pm1.96\,\hat\sigma\omega_j/\sqrt n \), with node-wise lasso rows
  (\( \lambda_j=0.024 \)) and \( \hat\sigma^2 \) the lasso residual sum of squares divided by \( n \) minus the number
  of selected variables;
- the **oracle** \( t \) interval from least squares on the true support plus the target.

The results are in [Figure 28.5.1](#fig-hd-debiased). For the null coefficient \( \beta_2 \), the lasso selected \( x_2 \) in
a fraction \( 0.069 \) of the data sets, and in those the naive interval covered zero only in a fraction
\( 0.174 \). For the weak coefficient, selected in a fraction \( 0.329 \), the naive coverage was
\( 0.903 \), and for the strong one \( 0.915 \). The debiased intervals covered \( \beta_1 \), \( \beta_6 \) and
\( \beta_2 \) in fractions \( 0.960 \), \( 0.946 \) and \( 0.940 \) respectively, close to the oracle's
\( 0.943 \), \( 0.956 \) and \( 0.958 \) for the same three coefficients.

On the strong coefficient the lasso's average error is \( -0.210 \), close to \( -\lambda \), and the
debiased estimate's is \( 0.014 \). The standardized estimates \( T_1 \) have mean \( 0.14 \) and
standard deviation \( 0.95 \) (panel (a)). Two costs remain.

The debiased intervals are wider than the oracle's (\( 0.393 \) against \( 0.288 \) for \( \beta_1 \)): part of
this is the price of not knowing the support, and part is tuning. The node-wise penalty
\( \lambda_j=0.024 \), a tenth of \( \sqrt{2\log p/n} \), keeps the bound \( \mu_j\le\lambda_j/\hat\tau_j^2 \) of
@lem-hd-nodewise, and so the remainder, small. But it lets the node-wise fits explain much of each column:
\( \hat\tau_j^2 \) (\( 0.28 \), \( 0.25 \) and \( 0.27 \) for \( \beta_1 \), \( \beta_6 \) and \( \beta_2 \)) is well below the
population value \( 1/(\bSigma^{-1})_{jj} \), which is \( \frac34 \) for \( j=1 \) and \( \frac35 \) otherwise, and that
allows a larger \( \omega_j \) (@exr-hd-nodewise-bound).

The second cost is a remainder bias: \( 0.049 \) for \( \beta_2 \) and \( 0.046 \) for the weak \( \beta_6 \), about
a third of its value. This is the term \( r_j/\sqrt n \) of @thm-hd-debiased, at most
\( \mu_j\norm{\hat{\bbeta}_\lambda-\bbeta}_1 \), with \( \mu_j \) between \( 0.088 \) and \( 0.099 \) and an average
lasso \( \ell_1 \) error of \( 0.74 \). At \( n=200 \) that product is not negligible.
:::

::: {when-format="html"}
![**Figure 28.5.1.** Inference for single coefficients with \( n=200 \), \( p=400 \), \( 1000 \) data sets.
(a) Histogram of the standardized debiased estimates of \( \beta_1 \) and the \( \Normal(0,1) \) density.
(b) Coverage of nominal 95% intervals: naive least squares after lasso selection (computed only when the
target is selected), the debiased lasso, and the oracle.](debiased.svg){#fig-hd-debiased width=100%}
:::

::: {when-format="pdf"}
![Inference for single coefficients with \( n=200 \), \( p=400 \), \( 1000 \) data sets.
(a) Histogram of the standardized debiased estimates of \( \beta_1 \) and the \( \Normal(0,1) \) density.
(b) Coverage of nominal 95% intervals: naive least squares after lasso selection (computed only when the
target is selected), the debiased lasso, and the oracle.](debiased.pdf){width=100%}
:::

```{.python .run #cell-debiased-small}
import numpy as np
from scipy import stats

def soft(z, t):
    """Soft thresholding S(z, t) = sign(z) max(|z| - t, 0)."""
    return np.sign(z) * np.maximum(np.abs(z) - t, 0.0)

def lasso_cd(X, Y, lam, B=None, tol=1e-9, max_sweeps=2000):
    """Minimize ||y - X b||^2 / (2n) + lam ||b||_1 for every column y of Y, by coordinate descent."""
    n, p = X.shape
    Y = Y.reshape(n, -1)
    B = np.zeros((p, Y.shape[1])) if B is None else B.copy()
    R = Y - X @ B                                    # residuals, one column per response
    scale = np.sum(X ** 2, axis=0) / n
    active = np.arange(p)
    for sweep in range(max_sweeps):
        change = 0.0
        for j in active:
            z = X[:, j] @ R / n + scale[j] * B[j]
            new = soft(z, lam) / scale[j]
            d = new - B[j]
            if np.any(d != 0):
                R -= np.outer(X[:, j], d)
                B[j] = new
                change = max(change, np.max(np.abs(d)))
        if change < tol:
            if len(active) == p:                     # converged over all coordinates
                break
            active = np.arange(p)                    # check every coordinate once more
        else:
            active = np.flatnonzero(np.any(B != 0, axis=1)) if sweep % 5 else np.arange(p)
    return B

def nodewise_row(X, j, lam_j):
    """Row theta_j of Theta from the lasso of x_j on the other columns, and tau_j^2."""
    n, p = X.shape
    others = np.delete(np.arange(p), j)
    gamma = lasso_cd(X[:, others], X[:, j], lam_j)[:, 0]
    z = X[:, j] - X[:, others] @ gamma               # lasso residual of x_j on the rest
    tau2 = X[:, j] @ z / n
    theta = np.zeros(p)
    theta[j] = 1.0
    theta[others] = -gamma
    return theta / tau2, tau2

def debiased_interval(X, Y, B, j, theta, level=0.95):
    """Debiased estimate of beta_j and its normal-theory interval, one per column of Y."""
    n = X.shape[0]
    R = Y - X @ B
    b = B[j] + theta @ X.T @ R / n                    # one-step correction
    sigma_hat = np.sqrt(np.sum(R ** 2, axis=0) / (n - np.sum(B != 0, axis=0)))
    G = X.T @ X / n
    se = sigma_hat * np.sqrt(theta @ G @ theta / n)
    z = stats.norm.ppf((1 + level) / 2)
    return b, b - z * se, b + z * se

def toeplitz_design(rng, n, p, rho):
    idx = np.arange(p)
    L = np.linalg.cholesky(rho ** np.abs(idx[:, None] - idx[None, :]))
    X = rng.normal(size=(n, p)) @ L.T
    return X / np.sqrt(np.sum(X ** 2, axis=0) / n)

rng = np.random.default_rng(2850)
n, p, sigma = 200, 400, 1.0
X = toeplitz_design(rng, n, p, rho=0.5)
beta = np.zeros(p)
beta[[0, 5, 10, 15]] = [1.0, 0.15, -1.0, 0.15]
lam = sigma * np.sqrt(2 * np.log(p) / n)
targets = {"strong": 0, "weak": 5, "null": 1}

Ys = (X @ beta)[:, None] + sigma * rng.normal(size=(n, 200))
Bs = lasso_cd(X, Ys, lam)
for name, j in targets.items():
    theta, tau2 = nodewise_row(X, j, 0.1 * np.sqrt(2 * np.log(p) / n))
    b, lo, hi = debiased_interval(X, Ys, Bs, j, theta)
    print(f"{name:6s}: lasso bias {Bs[j].mean() - beta[j]:+.3f}, debiased bias {b.mean() - beta[j]:+.3f},"
          f" coverage {np.mean((lo <= beta[j]) & (beta[j] <= hi)):.3f}")
```

## Other routes to valid inference

Other approaches answer different questions. **Sample splitting** (Wasserman and Roeder 2009) selects on one
half of the data and fits on the other (@exr-hd-sample-splitting). **Selective inference** (Lee, Sun, Sun and
Taylor 2016) conditions on the lasso having chosen the reported model, a union of polyhedra in \( \Y \).
**Simultaneous inference over all models** (Berk, Brown, Buja, Zhang and Zhao 2013) widens the intervals, as
Scheffé's method does. These target the coefficients of the *selected* model, which generally differ from
\( \beta_j \). **Double selection** (Belloni, Chernozhukov and Hansen 2014) targets \( \beta_j \): it keeps the controls that
predict the response or \( x_j \) and refits least squares, the idea of @eq-hd-debiased-fwl.

## Exercises

### A. Check your understanding

::: {#exr-hd-debias-low-dim}
[A1]

Suppose \( p<n \), \( \X \) has full column rank and \( \hat{\boldsymbol{\Theta}}=\hat{\bSigma}^{-1} \). Show that the debiased lasso @eq-hd-debiased
equals the least squares estimate, whatever \( \hat{\bbeta}_\lambda \) is, and that \( \mu_j=0 \) and \( \mathbf{r}=\bzero \) in @thm-hd-debiased.
:::

::: {.solution}
\( \hat{\mathbf{b}}=\hat{\bbeta}_\lambda+(\X\T\X)^{-1}\X\T(\Y-\X\hat{\bbeta}_\lambda)=\hat{\bbeta}_\lambda+\hbeta-\hat{\bbeta}_\lambda=\hbeta \), where
\( \hbeta=(\X\T\X)^{-1}\X\T\Y \). Then \( \hat{\bSigma}\hat{\boldsymbol{\Theta}}=\I \), so every \( \mu_j=0 \) and \( \mathbf{r}=\bzero \), and part (a) of the theorem
reduces to \( \sqrt n(\hbeta-\bbeta)\sim\Normal_p(\bzero,\sigma^2\hat{\bSigma}^{-1}) \) (@thm-opt-sampling).
:::

::: {#exr-hd-lasso-interval}
[A2]

On the event of @thm-hd-support, with \( \hat{\bSigma}_{SS}=\I \), show that the interval
\( \hat\beta_{\lambda,j}\pm1.96\,\sigma/\sqrt n \) for \( j\in S \) covers \( \beta_j \) iff
\( |\x_j\T\be/\sqrt n\,\sigma^{-1}-\sqrt n\lambda\operatorname{sign}(\beta_j)/\sigma|\le1.96 \). If \( \lambda=2\lambda_0 \), what happens to its
coverage as \( p\to\infty \)?
:::

### B. Practice

::: {#exr-hd-sample-splitting}
[B1]

*(Sample splitting.)* Split the observations into halves \( I_1 \) and \( I_2 \) chosen independently of the data.
Select a set \( \hat S \) of columns using only the data in \( I_1 \), and suppose \( \X_{I_2,\hat S} \) has full column rank and
\( |\hat S|<|I_2| \). Show that if \( \hat S\supseteq S \), then conditionally on the data in \( I_1 \), the usual \( t \) interval for
\( \beta_j \) (\( j\in\hat S \)) from the least squares fit of \( \Y_{I_2} \) on \( \X_{I_2,\hat S} \) has exact coverage \( 1-\alpha \). What does
the fit estimate when \( \hat S \) misses an active variable?
:::

::: {.solution}
The errors in \( I_2 \) are independent of those in \( I_1 \), so conditionally on the data in \( I_1 \) the set \( \hat S \) is fixed
and \( \Y_{I_2}\sim\Normal(\X_{I_2}\bbeta,\sigma^2\I) \). If \( \hat S\supseteq S \), then \( \X_{I_2}\bbeta=\X_{I_2,\hat S}\bbeta_{\hat S} \), so the
fit on \( \X_{I_2,\hat S} \) is a correctly specified normal linear model with coefficients \( \bbeta_{\hat S} \), and the \( t \) interval
of @thm-ci-estimable-interval has exact conditional coverage. Averaging over the data in \( I_1 \), coverage is exact on
the event \( \hat S\supseteq S \). If \( \hat S \) misses an active variable, the fit estimates the projected coefficients
\( (\X_{I_2,\hat S}\T\X_{I_2,\hat S})^{-1}\X_{I_2,\hat S}\T\X_{I_2}\bbeta \) (@prp-lm-misspecified), which differ from \( \bbeta_{\hat S} \) when the
missed variables are correlated with the selected ones, and the residual mean square overestimates \( \sigma^2 \).
:::

::: {#exr-hd-nodewise-bound}
[B2]

In @lem-hd-nodewise, show that \( \omega_j^2=\hat{\boldsymbol{\uptheta}}_j\T\hat{\bSigma}\hat{\boldsymbol{\uptheta}}_j=\norm{\hat{\bz}_j}^2/(n\hat\tau_j^4) \), and hence
\( \omega_j^2\le1/\hat\tau_j^2 \). Interpret: a small \( \hat\tau_j^2 \) (a column that the others explain well) means a
large variance for \( \hat b_j \). Compare with the variance inflation factor (@def-col-vif).
:::

::: {.solution}
\( \X\hat{\boldsymbol{\uptheta}}_j=\hat{\bz}_j/\hat\tau_j^2 \), so \( \omega_j^2=\norm{\X\hat{\boldsymbol{\uptheta}}_j}^2/n=\norm{\hat{\bz}_j}^2/(n\hat\tau_j^4) \). By
@lem-hd-nodewise(a), \( \norm{\hat{\bz}_j}^2/n\le\hat\tau_j^2 \), so \( \omega_j^2\le1/\hat\tau_j^2 \). With least squares in place
of the lasso (\( p<n \), \( \lambda_j=0 \)), \( \hat\tau_j^2=\norm{(\I-\M_{(j)})\x_j}^2/n \), and \( \omega_j^2=1/\hat\tau_j^2 \) is the
variance factor \( n/\{S_{jj}(1-R_j^2)\} \) of @eq-proj-vif-preview for a column with \( \norm{\x_j}^2=n \) (and an intercept among
the other columns). The node-wise lasso residual plays the role of the FWL residual, and \( 1/\hat\tau_j^2 \) the role of
the variance inflation.
:::
