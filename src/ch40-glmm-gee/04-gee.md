# Generalized estimating equations

Suppose the question is the marginal one and the cluster structure a nuisance. Then a
generalized linear mixed model asks for more than is needed: a distribution for the random
effects, a conditional distribution for the responses, and an intractable integral.
Everything the question needs is in the *mean*, and
[Chapter 38](../ch38-quasi-likelihood/index.html) showed that a mean and a variance function
suffice to build a usable estimating equation without a likelihood.

Carrying that across to clustered data gives the method of Liang and Zeger (1986), the most
widely used in this chapter and the most widely misunderstood. It estimates a marginal
coefficient consistently even when the correlation structure it assumes is wrong. It does
not estimate a conditional coefficient at all; it produces no likelihood, so most of the
model-comparison machinery of the book does not apply; and its robustness to a wrong
correlation does not extend to data that go missing for the wrong reasons.

## The estimating equation

Clusters \( i=1,\dots,m \) are independent, with response vectors \( \Y_i \) of length
\( n_i \) and model matrices \( \X_i \) of size \( n_i\times p \). The model constrains the
marginal mean only:
\[
\E(\Y_i)=\bmu_i(\bbeta),\qquad
\mu_{ij}(\bbeta)=h\bigl(\x_{(ij)}\T\bbeta\bigr),
\]{#eq-gmm-gee-mean}

with \( h \) the inverse of a link as in @def-glm-link. Nothing is assumed about the joint
distribution of \( \Y_i \), and in particular nothing about its covariance — although a
guess at that covariance is used, deliberately.

::: {#def-gmm-gee}
[Generalized estimating equations]

Let \( V \) be a variance function and \( \phi \) a dispersion, and write
\[
\begin{aligned}
\A^{V}_i(\bbeta)&=\diag\{V(\mu_{i1}),\dots,V(\mu_{in_i})\},\\
\bD_i(\bbeta)&=\frac{\partial\bmu_i}{\partial\bbeta\T}
=\diag\{h'(\eta_{ij})\}\,\X_i .
\end{aligned}
\]
Let \( \R_i(\boldsymbol{\upalpha}) \) be a **working correlation matrix**: a correlation
matrix for the \( n_i \) observations of cluster \( i \), depending on a finite parameter
\( \boldsymbol{\upalpha} \), chosen by the analyst and *not* assumed correct. The **working
covariance** is
\[
\V_i(\bbeta,\boldsymbol{\upalpha},\phi)
=\phi\,(\A^{V}_i)^{1/2}\,\R_i(\boldsymbol{\upalpha})\,(\A^{V}_i)^{1/2},
\]{#eq-gmm-working-covariance}

and the **generalized estimating equation** is
\[
\bU_m(\bbeta;\boldsymbol{\upalpha},\phi)
=\sum_{i=1}^{m}\bD_i(\bbeta)\T\,\V_i^{-1}\,\bigl\{\Y_i-\bmu_i(\bbeta)\bigr\}=\bzero .
\]{#eq-gmm-gee}

In practice \( \boldsymbol{\upalpha} \) and \( \phi \) are replaced by moment estimates
computed from the Pearson residuals
\( r_{ij}=(y_{ij}-\hat\mu_{ij})/\sqrt{V(\hat\mu_{ij})} \) at the current \( \hat\bbeta \),
and @eq-gmm-gee is re-solved, until nothing changes.
:::

Within this chapter \( \R_i \) denotes the working correlation matrix, not the error
covariance \( \R \) of Part VII.

::: {.remark}
[Three familiar special cases]

With \( \R_i=\I \) and one observation per cluster, @eq-gmm-gee is exactly the likelihood
equation @eq-glm-score, with \( \W\bD^{-1} \) written out. With \( \R_i=\I \) and several
observations per cluster it is still that equation, applied as if the observations were
independent: the same \( \hat\bbeta \), different standard errors. And with a general
\( \R_i \) but one cluster, it is the quasi-score of
[Chapter 38](../ch38-quasi-likelihood/index.html) (@thm-ql-score) with a full
covariance matrix — which for the identity link and \( V\equiv1 \) is the generalized least
squares of @prp-cls-selection.
:::

Solving @eq-gmm-gee is Fisher scoring. Write \( \bS_i=\diag\{h'(\eta_{ij})/\sqrt{V(\mu_{ij})}\} \);
then \( \bD_i\T\V_i^{-1}=\phi^{-1}\X_i\T\bS_i\R_i^{-1}(\A^{V}_i)^{-1/2} \), so the update is
\[
\hbeta\ \leftarrow\ \hbeta+
\Bigl(\sum_i\X_i\T\bS_i\R_i^{-1}\bS_i\X_i\Bigr)^{-1}
\sum_i\X_i\T\bS_i\R_i^{-1}\br_i ,
\]{#eq-gmm-gee-update}

in which \( \phi \) cancels. It is the iteratively reweighted least squares of @thm-glm-irls
with \( \R_i^{-1} \) inserted between the weights, and it converges just as
quickly.

## What the estimate does and does not have

::: {#thm-gmm-gee}
[Consistency, asymptotic normality and the sandwich]

Let the clusters be independent, let \( \bbeta^{0} \) be the true value in @eq-gmm-gee-mean,
and let \( \bSigma_i^{0}=\Cov(\Y_i) \) be arbitrary and unknown. Assume

- (E1) \( n_i\le n_{\max} \) and \( \norm{\x_{(ij)}}\le C \) for all \( i,j \);
- (E2) \( h \) is twice continuously differentiable with \( h' \) bounded, \( V \) is
  continuous and bounded away from \( 0 \) on the relevant range, and \( \bbeta^{0} \) is
  interior to a compact parameter set;
- (E3) the eigenvalues of \( \R_i(\boldsymbol{\upalpha}) \) are bounded above and away from
  zero, uniformly in \( i \) and in \( \boldsymbol{\upalpha} \) near its limit;
- (E4) \( \E\norm{\Y_i}^{4}\le K \) for all \( i \);
- (E5) \( \hat{\boldsymbol{\upalpha}} \) and \( \hat\phi \) are \( \sqrt m \)-consistent for
  some limits \( \boldsymbol{\upalpha}^{*},\phi^{*} \), which need not make
  \( \V_i \) equal to \( \bSigma_i^{0} \);
- (E6) with \( \V_i \) evaluated at \( (\bbeta^{0},\boldsymbol{\upalpha}^{*},\phi^{*}) \),
  \[
  \begin{aligned}
  \A_m&=\frac1m\sum_i\bD_i\T\V_i^{-1}\bD_i\to\A,\\
  \B_m&=\frac1m\sum_i\bD_i\T\V_i^{-1}\bSigma_i^{0}\V_i^{-1}\bD_i\to\B,
  \end{aligned}
  \]
  both positive definite.

Then:

::: {.enumerate options="label=(\alph*)"}
1. **Unbiasedness.** For every fixed \( \boldsymbol{\upalpha} \) and \( \phi \),
   \( \E\{\bU_m(\bbeta^{0};\boldsymbol{\upalpha},\phi)\}=\bzero \). The estimating function
   is unbiased whatever the working correlation is.

2. **Consistency and asymptotic normality.** With probability tending to one, @eq-gmm-gee
   has a root \( \hbeta_m \) with \( \hbeta_m\to\bbeta^{0} \) in probability,
   and
   \[
   \sqrt m\,(\hbeta_m-\bbeta^{0})\ \xrightarrow{d}\ \Normal_p\bigl(\bzero,\ \A^{-1}\B\A^{-1}\bigr).
   \]{#eq-gmm-gee-asymptotics}

3. **The sandwich estimator.** Let \( \br_i=\Y_i-\bmu_i(\hbeta_m) \) and
   \[
   \begin{aligned}
   \hat{\V}_{\text{R}}&=\hat\A_m^{-1}
   \Bigl(\sum_i\hat\bD_i\T\hat\V_i^{-1}\br_i\br_i\T\hat\V_i^{-1}\hat\bD_i\Bigr)
   \hat\A_m^{-1},\\
   \hat\A_m&=\sum_i\hat\bD_i\T\hat\V_i^{-1}\hat\bD_i ,
   \end{aligned}
   \]{#eq-gmm-sandwich}

   with hats denoting evaluation at \( \hbeta_m,\hat{\boldsymbol{\upalpha}},\hat\phi \);
   \( \hat\A_m \) is \( m \) times the sample version of the \( \A_m \) of (E6).
   Then \( m\,\hat{\V}_{\text{R}}\to\A^{-1}\B\A^{-1} \) in probability, and the Wald
   statistic of a hypothesis \( \bLambda\T\bbeta \) of rank \( q \), formed with
   \( \hat{\V}_{\text{R}}\), tends to \( \chi^{2}(q) \).

4. **When the model-based covariance is valid.** If the working covariance is correct,
   \( \V_i=\bSigma_i^{0} \), then \( \B=\A \) and @eq-gmm-gee-asymptotics reduces to
   \( \A^{-1} \), which the **model-based** estimator
   \( \hat{\V}_{\text{M}}=(\sum_i\hat\bD_i\T\hat\V_i^{-1}\hat\bD_i)^{-1} \) estimates
   consistently. Otherwise \( \hat{\V}_{\text{M}} \) is inconsistent, in either direction.
:::

:::

::: {.proof}
(a) At \( \bbeta=\bbeta^{0} \) the matrices \( \bD_i \) and \( \V_i \) are deterministic —
they depend on the data only through \( \bbeta \), which is held fixed — so
\( \E\bU_m=\sum_i\bD_i\T\V_i^{-1}\E\{\Y_i-\bmu_i(\bbeta^{0})\}=\bzero \) by @eq-gmm-gee-mean.
Nothing about \( \R_i \) enters, and that is the whole source of the method's robustness.

(b) and (c). Expand \( \bU_m(\bbeta;\boldsymbol{\upalpha}) \) about \( \bbeta^{0} \) at
\( \boldsymbol{\upalpha}=\hat{\boldsymbol{\upalpha}} \). Since \( \bU_m \) is
vector-valued, no single intermediate point serves, so we use the integral form of Taylor's
theorem:
\[
\begin{aligned}
\bzero&=\bU_m(\hbeta_m;\hat{\boldsymbol{\upalpha}})
=\bU_m(\bbeta^{0};\hat{\boldsymbol{\upalpha}})
+\bar{\mathbf{M}}_m\,(\hbeta_m-\bbeta^{0}),\\
\bar{\mathbf{M}}_m&=\int_0^{1}\frac{\partial\bU_m}{\partial\bbeta\T}
\bigl(\bbeta^{0}+t(\hbeta_m-\bbeta^{0})\bigr)\,dt .
\end{aligned}
\]
Three facts give the result.

*The derivative.* \( \partial\bU_m/\partial\bbeta\T
=-\sum_i\bD_i\T\V_i^{-1}\bD_i+\bT_m \), where \( \bT_m \) collects the terms in which the
derivative falls on \( \bD_i\T\V_i^{-1} \) and therefore still multiplies
\( \Y_i-\bmu_i \). Each summand of \( \bT_m \) has mean zero at \( \bbeta^{0} \) and bounded
variance by (E1), (E2) and (E4), so \( m^{-1}\bT_m\to\mathbf{0} \) in probability by the weak
law for independent summands; and \( m^{-1}\sum_i\bD_i\T\V_i^{-1}\bD_i\to\A \) by (E6).
Hence \( m^{-1}\partial\bU_m/\partial\bbeta\T\to-\A \), uniformly for \( \bbeta \) in a
neighbourhood of \( \bbeta^{0} \) under (E2); the convergence being uniform over that
neighbourhood, \( m^{-1}\bar{\mathbf{M}}_m\to-\A \) as well, once \( \hbeta_m \) is in
it.

*The score.* \( m^{-1/2}\bU_m(\bbeta^{0};\boldsymbol{\upalpha}^{*}) \) is
\( m^{-1/2} \) times a sum of \( m \) independent mean-zero vectors
\( \mathbf{s}_i=\bD_i\T\V_i^{-1}(\Y_i-\bmu_i) \) with
\( \sum_i\Cov(\mathbf{s}_i)=m\B_m \). Conditions (E1)–(E4) bound \( \E\norm{\mathbf{s}_i}^{4} \)
uniformly, which gives the Lyapunov condition, so the Lindeberg–Feller theorem applies and
\( m^{-1/2}\bU_m\to\Normal_p(\bzero,\B) \).

*Estimating \( \boldsymbol{\upalpha} \) does not matter.* Differentiating @eq-gmm-gee with
respect to \( \boldsymbol{\upalpha} \) leaves the factor \( \Y_i-\bmu_i \) untouched, so
\( \partial\bU_m/\partial\boldsymbol{\upalpha}\T \) at \( \bbeta^{0} \) is again a sum of
mean-zero terms and \( m^{-1}\partial\bU_m/\partial\boldsymbol{\upalpha}\T\to\mathbf{0} \).
Therefore
\( m^{-1/2}\bU_m(\bbeta^{0};\hat{\boldsymbol{\upalpha}})
=m^{-1/2}\bU_m(\bbeta^{0};\boldsymbol{\upalpha}^{*})
+\{m^{-1}\partial\bU_m/\partial\boldsymbol{\upalpha}\T\}\sqrt m(\hat{\boldsymbol{\upalpha}}
-\boldsymbol{\upalpha}^{*})
=m^{-1/2}\bU_m(\bbeta^{0};\boldsymbol{\upalpha}^{*})+o_p(1) \) by (E5).

Combining, \( \sqrt m(\hbeta_m-\bbeta^{0})=\A^{-1}m^{-1/2}\bU_m(\bbeta^{0})+o_p(1) \), which
is @eq-gmm-gee-asymptotics. For (c), \( \br_i\br_i\T \) has mean
\( \bSigma_i^{0}+O(m^{-1}) \) at \( \hbeta_m \), and the argument of @thm-het-sandwich with
the cluster in place of the observation — the version already used for @prp-cls-selection(b)
— gives \( m^{-1}\sum_i\hat\bD_i\T\hat\V_i^{-1}\br_i\br_i\T\hat\V_i^{-1}\hat\bD_i\to\B \).
The Wald limit follows from @eq-gmm-gee-asymptotics and the continuous mapping theorem.

(d) If \( \V_i=\bSigma_i^{0} \) then
\( \bD_i\T\V_i^{-1}\bSigma_i^{0}\V_i^{-1}\bD_i=\bD_i\T\V_i^{-1}\bD_i \) for every \( i \), so
\( \B_m=\A_m \).
:::

::: {.remark}
[What is proved and what is assumed]

The existence of a consistent root and the uniform control of the derivative are the two
places where the argument leans on standard machinery rather than its own steam; both go
exactly as for @thm-glm-asymptotics, using the concavity-free argument for estimating
equations in Liang and Zeger (1986). Condition (E5) is a real assumption: moment estimators
of \( \boldsymbol{\upalpha} \) from Pearson residuals are \( \sqrt m \)-consistent under
(E1)–(E4), but an unstructured \( \R \) estimated from few clusters is not. Note the
notation: \( \bU_m \) is the estimating *function*, as in
[Chapter 38](../ch38-quasi-likelihood/index.html), while a hatted \( \V \) with a letter
subscript is a *covariance estimator* for \( \hbeta \), as in
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html); \( \V_i \), with a
cluster subscript, remains the working covariance of cluster \( i \).
:::

::: {.warning}
[The mean model must hold for the whole cluster]

@eq-gmm-gee-mean constrains \( \E(\Y_i) \), the mean of the whole cluster vector, and part
(a) uses it in that form; what one usually has in mind is the weaker
\( \E(Y_{ij}\mid\x_{(ij)})=h(\x_{(ij)}\T\bbeta) \). A *non-diagonal* \( \R_i \) mixes
observations across occasions, so unbiasedness then needs the **full-covariate conditional
mean** condition \( \E(Y_{ij}\mid\X_i)=h(\x_{(ij)}\T\bbeta) \): the mean at occasion
\( j \) must not depend on the covariates at other occasions. With time-varying covariates
and feedback — this year's investment changing next year's market value, a dose adjusted in
response to last month's reading — that fails even though the occasion-by-occasion mean
model is correct, and a non-diagonal working correlation is then *inconsistent* while
independence GEE is not (Pepe and Anderson 1994; Diggle, Heagerty, Liang and Zeger 2002,
chapter 12).
:::

::: {.idea}
The asymptotics are in \( m \), the number of clusters; the total sample size is
irrelevant. A trial with ten thousand patients in six hospitals has six independent pieces
of information about anything that varies between hospitals, and the sandwich cannot
manufacture more — the lesson of
[Section 33.2](../ch33-clustered-longitudinal-splitplot/02-clustered-data.html).
:::

## Travel-mode choice

::: {#exm-gmm-modechoice}
[Four options per traveller]

Each of \( 210 \) travellers between Australian cities chose one of four
modes — air, train, bus, car — and the data record, for every traveller and every mode, the
terminal waiting time and a generalized cost. Take \( Y_{ij}=1 \) if traveller \( i \)
chose mode \( j \), fit a logistic model with mode indicators, waiting time and cost, and
treat each traveller as a cluster of \( n_i=4 \). The sample deliberately over-represents
the rarer modes, so by @thm-bin-casecontrol the mode indicators and the intercept are
shifted by the sampling design and are not population quantities; the within-traveller
coefficients are unaffected, and standard errors are the point here in any case.

The four records of a traveller are anything but independent: exactly one is a \( 1 \). An
independence working correlation gives the same coefficients as ordinary logistic
regression, so the question is only what the standard errors should be. The ratio of the
sandwich standard error to the model-based one ranges from
\( 0.694 \) to \( 1.558 \) across the six
coefficients: for terminal waiting time, whose coefficient is
\( -1.0693 \) per ten minutes, the sandwich standard error is
\( 1.558 \) times larger (\( 0.1698 \)
against \( 0.1090 \)); for generalized cost, whose coefficient is
\( -0.0098 \) per ten dollars, it is
\( 0.694 \) times as large
(\( 0.0138 \) against \( 0.0199 \)).

The sandwich is not a device for making standard errors bigger. Negative within-cluster
correlation makes *within-cluster contrasts* better determined than independence would
suggest, and cost is such a contrast; waiting time is confounded with the mode indicators,
which the traveller's single choice constrains, and is therefore worse determined. Two
further observations. The Pearson dispersion estimate is
\( 5.155 \), which looks like extreme overdispersion and means nothing:
Bernoulli responses carry no information about a dispersion parameter (@prp-gmm-moments(b)
with \( N_i=1 \)). And the exchangeable working correlation estimate is
\( -0.3015 \), close to the smallest value
\( -1/(n_i-1)=-1/3 \) an exchangeable matrix of size four can have: there the smallest
eigenvalue of \( \R \) is \( 1+3\hat\alpha=0.0956 \). An almost
singular structure is a bad working covariance, and independence is the sensible choice.
:::

```{.python .run #cell-gee-modechoice-gee}
import numpy as np

def expit(z):
    return 1.0 / (1.0 + np.exp(-z))

def working_correlation(kind, r, p, scale_fixed=False):
    """Moment estimate of the working correlation from the Pearson residuals r (m by n)."""
    m, n = r.shape
    phi = 1.0 if scale_fixed else np.sum(r**2) / (m * n - p)
    if kind == "independence":
        return np.eye(n), 0.0, phi
    if kind == "exchangeable":
        off = (np.sum(r.sum(1) ** 2) - np.sum(r**2)) / 2
        a = off / (phi * (m * n * (n - 1) / 2 - p))
        return np.eye(n) + a * (1 - np.eye(n)), a, phi
    if kind == "ar1":
        a = np.sum(r[:, :-1] * r[:, 1:]) / (phi * (m * (n - 1) - p))
        lag = np.abs(np.arange(n)[:, None] - np.arange(n))
        return a**lag, a, phi
    if kind == "unstructured":
        R = (r.T @ r) / (phi * m)
        R = R / np.sqrt(np.outer(np.diag(R), np.diag(R)))
        return R, np.nan, phi
    raise ValueError(kind)

def gee(y, X, kind="independence", family="binomial", scale_fixed=False, maxit=200):
    """Solve sum_i D_i' V_i^{-1} (y_i - mu_i) = 0 for balanced clusters y (m, n), X (m, n, p).

    With S_i = diag(mu'_ij / sqrt(V(mu_ij))) and Pearson residuals r_i, the Fisher-scoring
    step is beta += (sum_i X_i'S_i R^{-1} S_i X_i)^{-1} sum_i X_i'S_i R^{-1} r_i."""
    m, n, p = X.shape
    beta, R, alpha, phi = np.zeros(p), np.eye(n), 0.0, 1.0
    for _ in range(maxit):
        eta = X @ beta
        if family == "binomial":
            mu = expit(eta)
            dmu, V = mu * (1 - mu), mu * (1 - mu)
        else:                                            # log link, V(mu) = mu^2 (gamma)
            mu = np.exp(eta)
            dmu, V = mu, mu**2
        r = (y - mu) / np.sqrt(V)
        R, alpha, phi = working_correlation(kind, r, p, scale_fixed)
        Rinv = np.linalg.inv(R)
        S = dmu / np.sqrt(V)
        SX = S[..., None] * X                            # (m, n, p)
        A = np.einsum("mjp,jk,mkq->pq", SX, Rinv, SX)
        g = np.einsum("mjp,jk,mk->p", SX, Rinv, r)
        step = np.linalg.solve(A, g)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-11:
            break
    u = np.einsum("mjp,jk,mk->mp", SX, Rinv, r)          # cluster contributions
    Ainv = np.linalg.inv(A)
    robust = Ainv @ (u.T @ u) @ Ainv                     # the sandwich covariance
    naive = phi * Ainv                                   # the model-based covariance
    return dict(beta=beta, naive=naive, robust=robust, alpha=alpha, phi=phi, R=R,
                A=A, Ainv=Ainv, S=S, r=r, Rinv=Rinv, X=X)
```

```{.python .run #cell-gee-modechoice-data}
import pandas as pd
import statsmodels.api as sm

d = sm.datasets.modechoice.load_pandas().data.sort_values(["individual", "mode"])
m, n = d["individual"].nunique(), 4                       # 210 travellers, 4 modes each
y = d["choice"].to_numpy().reshape(m, n)
D = pd.get_dummies(d["mode"].astype(int), drop_first=True).to_numpy(dtype=float)
cols = [np.ones(m * n)] + [D[:, j] for j in range(3)] + [d["ttme"] / 10, d["gc"] / 10]
X = np.stack([np.asarray(c, dtype=float).reshape(m, n) for c in cols], axis=2)

# a Bernoulli response carries no information about a dispersion parameter, so phi = 1
fit = gee(y, X, "independence", scale_fixed=True)
se_naive, se_robust = np.sqrt(np.diag(fit["naive"])), np.sqrt(np.diag(fit["robust"]))
names = ["intercept", "train", "bus", "car", "ttme/10", "gc/10"]
for k, name in enumerate(names):
    print(f"{name:10s} {fit['beta'][k]:8.4f}  naive {se_naive[k]:.4f}  robust {se_robust[k]:.4f}")
print("exchangeable alpha:", round(gee(y, X, 'exchangeable', scale_fixed=True)['alpha'], 4))
```

## Small samples, and a simulation

The sandwich @eq-gmm-sandwich replaces \( \bSigma_i^{0} \) by \( \br_i\br_i\T \), and the
residuals are too small because \( \hbeta_m \) was fitted to them; the bias is of order
\( p/m \) and downward, exactly as for @thm-het-sandwich and @exr-cls-crve-bias. Two
standard repairs, both used in @exm-gmm-coverage:

- **Bias-corrected residuals.** Mancl and DeRouen (2001) replace \( \br_i \) by
  \( (\I-\bH_i)^{-1}\br_i \), with
  \( \bH_i=\bD_i(\sum_k\bD_k\T\V_k^{-1}\bD_k)^{-1}\bD_i\T\V_i^{-1} \) cluster \( i \)'s
  leverage block — the clustered analogue of the HC3 correction of @def-het-hc. Kauermann and
  Carroll (2001) use \( (\I-\bH_i)^{-1/2} \), which corrects less.
- **A \( t \) reference.** Replace the normal quantile by \( t_{m-p} \). With
  \( m=15 \) and \( p=3 \) this widens a \( 95\% \) interval by
  \( t_{12,0.975}/z_{0.975}=2.179/1.960 \), about eleven per cent.

::: {#exm-gmm-coverage}
[Coverage when the working correlation is wrong]

Generate clusters of \( n_i=8 \) binary responses from a random-intercept logistic model with
\( \tau=1.5 \), one covariate varying within the cluster and one fixed for the whole cluster
(a treatment assigned to clusters), and fit the marginal model by GEE with an
*independence* working correlation — wrong by construction. The target is the marginal
cluster-level coefficient, which @eq-gmm-marginal-parameter puts at
\( 0.5748 \) when the conditional coefficient is \( 0.8 \). Over
\( 600 \) replications at \( m=240 \) the mean estimate is
\( 0.5681 \): consistent for the marginal coefficient,
and not for the conditional one.

The three intervals behave quite differently. The model-based interval covers
\( 0.757 \) of the time and does not improve with \( m \): it is wrong,
not merely imprecise. The sandwich covers \( 0.913 \) at
\( m=15 \) and \( 0.947 \) at \( m=240 \), and the
bias-corrected sandwich with a \( t \) reference covers
\( 0.973 \) at \( m=15 \) — slightly conservative, which at
that sample size is the right side to err on.
:::

::: {when-format="html"}
![**Figure 40.4.1.** GEE with a wrong working correlation. (a) Coverage of nominal
\( 95\% \) intervals for the marginal coefficient of a cluster-level covariate, against the
number of clusters, for three covariance estimators. (b) The sampling distribution of
\( \hat\beta_2 \) at \( m=240 \), against the marginal and conditional
values.](gee_coverage.svg){#fig-gmm-coverage width=100%}
:::

::: {when-format="pdf"}
![GEE with a wrong working correlation. (a) Coverage of nominal
\( 95\% \) intervals for the marginal coefficient of a cluster-level covariate, against the
number of clusters, for three covariance estimators. (b) The sampling distribution of
\( \hat\beta_2 \) at \( m=240 \), against the marginal and conditional
values.](gee_coverage.pdf){width=100%}
:::

## What GEE does not give

*No likelihood.* @eq-gmm-gee is not the derivative of anything in general, so there is no
deviance, no likelihood ratio test, and no AIC or BIC. Nested mean models are compared by
Wald tests built on @eq-gmm-sandwich or by score-type tests built on @eq-gmm-gee under the
null; selection uses QIC, the quasi-likelihood analogue of AIC
([Section 40.5](05-working-correlation.html)). Reporting a "deviance" from GEE output is a
category error.

*No random effects and no distribution.* There is nothing to predict, so if the object of
interest is a cluster — which hospital is doing badly, which firm is unusually profitable —
GEE has no answer and @prp-gmm-prediction has no analogue; and because only the mean is
modelled, the fit cannot be simulated from, cannot give a prediction interval, and cannot be
checked by an envelope or posterior predictive method.

*Few clusters.* Everything above is asymptotic in \( m \). Below about forty clusters the
sandwich needs correcting, and below about fifteen no correction reliably saves it; a
likelihood-based method with a covariance model from @def-cls-covariance is then the better
bet.

## Missing data

The robustness of @eq-gmm-gee has a sharp limit. If some responses are not observed, write
\( R_{ij}=1 \) when \( Y_{ij} \) is recorded, so that the estimating function is really
\( \sum_i\bD_i\T\V_i^{-1}\bR_i(\Y_i-\bmu_i) \) with \( \bR_i=\diag(R_{ij}) \). Part (a)
of @thm-gmm-gee needs \( \E\{\bR_i(\Y_i-\bmu_i)\}=\bzero \), a genuine restriction on the
missingness.

If the data are **missing completely at random**, \( \bR_i \) is independent of
\( \Y_i \), the expectation factorizes, and everything above holds unchanged; the same is
true if missingness depends only on covariates that are in the model and fully observed,
because the mean model is then still correct conditionally on being observed. But if
missingness depends on *previously observed responses* — a patient who was doing badly last
visit is more likely to drop out — then \( \bR_i \) and \( \Y_i-\bmu_i \) are dependent,
the estimating function is biased and GEE is inconsistent. This is the **missing at random**
case of @def-mis-mechanisms, and it is exactly the case in which a likelihood-based analysis
remains valid.

That asymmetry is the strongest practical argument for the mixed-model route in longitudinal
studies with dropout, and it is not widely enough known. The repair is to weight each
cluster's contribution by the inverse of its estimated probability of being observed — the
weighted GEE of Robins, Rotnitzky and Zhao (1995), which trades the assumption about the
correlation for one about the missingness model.
[Chapter 41](../ch41-missing-data/index.html) develops the mechanisms, the ignorability
theorem and the alternatives.

## Exercises

### A. Check your understanding

::: {#exr-gmm-independence-same}
[A1]

Show from @eq-gmm-gee that with \( \R_i=\I \) the estimate \( \hbeta \) does not depend on
\( \phi \) or on the cluster structure at all, and equals the ordinary generalized linear
model estimate. Which quantities *do* depend on the clustering?
:::

::: {#exr-gmm-what-target}
[A2]

A trial randomizes clinics to a new protocol and measures a binary outcome on patients within
clinics. The investigators want to report "the effect of the protocol on a patient's odds of
recovery". Which method should they use, and what would a GEE analysis report instead?
:::

### B. Practice

::: {#exr-gmm-two-clusters}
[B1]

Take \( m \) clusters of size \( 2 \), the identity link, \( V\equiv1 \), \( \phi=1 \), and
an exchangeable working correlation with parameter \( \alpha \). Show that @eq-gmm-gee
reduces to weighted least squares with weight matrix \( \R^{-1} \),
\( \R=(1-\alpha)\I+\alpha\mathbf{J} \), in each cluster, and that when
\( \x_{(i1)}=\x_{(i2)} \) for every \( i \) the estimate does not depend on \( \alpha \) at
all. Interpret.
:::

::: {.solution}
The estimating equation is \( \sum_i\X_i\T\R^{-1}(\y_i-\X_i\bbeta)=\bzero \), which is the
normal equation of generalized least squares with covariance \( \R \). If the two rows of
\( \X_i \) are equal, \( \X_i=\bone\x_{(i)}\T \) and
\( \X_i\T\R^{-1}\X_i=(\bone\T\R^{-1}\bone)\x_{(i)}\x_{(i)}\T \),
\( \X_i\T\R^{-1}\y_i=\x_{(i)}\bone\T\R^{-1}\y_i \). For an exchangeable \( \R \),
\( \R^{-1}\bone=\{1/(1+\alpha)\}\bone \), so both sides carry the same scalar factor
\( 1/(1+\alpha) \) and it cancels. A cluster-constant covariate sees only the cluster mean,
and reweighting a mean by an exchangeable covariance does nothing.
:::

::: {#exr-gmm-unbiased-family}
[B2]

Let \( \mathbf{G}_i \) be any \( n_i\times p \) matrix depending on \( \bbeta \) but not on
\( \Y_i \). Show that \( \sum_i\mathbf{G}_i\T(\Y_i-\bmu_i) \) is an unbiased estimating
function and that @eq-gmm-gee is the member with
\( \mathbf{G}_i=\V_i^{-1}\bD_i \). Which member is optimal, in the sense of @prp-ql-optimality,
and why is it not available?
:::

::: {.solution}
Unbiasedness is immediate from \( \E(\Y_i-\bmu_i)=\bzero \). The Godambe-optimal choice
within this linear class is \( \mathbf{G}_i=(\bSigma_i^{0})^{-1}\bD_i \), the *true*
covariance in place of the working one; @prp-ql-optimality is the same statement with
\( n_i=1 \). It is unavailable because \( \bSigma_i^{0} \) is unknown, which is precisely
what makes a working correlation necessary.
:::

### C. Going deeper

::: {#exr-gmm-efficiency-bound}
[C1]

Show that if \( \V_i=\bSigma_i^{0} \) then
\( \A^{-1}\B\A^{-1}=\A^{-1} \) is the smallest asymptotic covariance attainable in the class
of @exr-gmm-unbiased-family, in the nonnegative-definite ordering. (Hint: for any
\( \mathbf{G}_i \), compare
\( (\sum\mathbf{G}_i\T\bD_i)^{-1}(\sum\mathbf{G}_i\T\bSigma_i^{0}\mathbf{G}_i)
(\sum\bD_i\T\mathbf{G}_i)^{-1} \) with \( \A^{-1} \) using a Cauchy–Schwarz argument in the
inner product \(
\langle\mathbf{a},\mathbf{b}\rangle=\sum_i\mathbf{a}_i\T\bSigma_i^{0}\mathbf{b}_i \).)
:::
