# Support recovery

Good prediction does not require the right variables: the bounds of [Section 28.3](03-lasso-bounds.html)
allow false variables with small coefficients and missed true ones. Often, though, the variables are the
object of interest. This section asks when the lasso selects *exactly* the right variables with the right
signs. The answer needs a condition on the design much stronger than those of
[Section 28.2](02-restricted-eigenvalues.html), and it fails in simple, realistic examples.

## Uniqueness

The question "which variables does the lasso select?" needs an answer that does not depend on which
minimizer is reported.

::: {#lem-hd-lasso-unique}
[Uniqueness of the lasso fit]

::: {.enumerate options="label=(\alph*)"}
1. All lasso estimates at a given \( \lambda \) have the same fitted values \( \X\hat{\bbeta}_\lambda \) and the same \( \ell_1 \)
   norm. Hence they share the vector \( \hat{\bz}=\X\T(\Y-\X\hat{\bbeta}_\lambda)/(n\lambda) \) of @eq-hd-optimality.

2. Let \( \hat{\bbeta}_\lambda \) be a lasso estimate with support \( \hat S \). If \( |\hat z_j|<1 \) for every \( j\notin\hat S \)
   and \( \X_{\hat S} \) has full column rank, then \( \hat{\bbeta}_\lambda \) is the only lasso estimate.
:::

:::

::: {.proof}
(a) The first statement is @thm-shr-lasso-kkt(c).
The vector \( \hat{\bz} \) depends on \( \hat{\bbeta}_\lambda \) only through \( \X\hat{\bbeta}_\lambda \).

(b) Let \( \bb \) be any lasso estimate. By (a) its optimality conditions @eq-hd-optimality hold with the same
\( \hat{\bz} \). If \( b_j\ne0 \) then \( \hat z_j=\operatorname{sign}(b_j) \), so \( |\hat z_j|=1 \) and \( j\in\hat S \). Thus \( \bb \) is
supported on \( \hat S \), and \( \X_{\hat S}\bb_{\hat S}=\X\bb=\X\hat{\bbeta}_\lambda=\X_{\hat S}\hat{\bbeta}_{\lambda,\hat S} \). Full column rank
of \( \X_{\hat S} \) gives \( \bb=\hat{\bbeta}_\lambda \).
:::

## The irrepresentable condition

When will the lasso choose exactly \( S \)? Suppose for a moment that it does, with the correct signs
\( \bz_S=\operatorname{sign}(\bbeta_S) \). Then @eq-hd-optimality on \( S \) determines the estimate:
\( \X_S\T(\Y-\X_S\hat{\bbeta}_S)/n=\lambda\bz_S \), so
\[
\begin{aligned}
\hat{\bbeta}_S&=(\X_S\T\X_S)^{-1}(\X_S\T\Y-n\lambda\bz_S)\\
&=\bbeta_S+(\X_S\T\X_S)^{-1}\X_S\T\be-\lambda\hat{\bSigma}_{SS}^{-1}\bz_S ,
\end{aligned}
\]{#eq-hd-oracle-lasso}

where \( \hat{\bSigma}_{SS}=\X_S\T\X_S/n \). The estimate is least squares on the true support, shifted towards zero by
\( \lambda\hat{\bSigma}_{SS}^{-1}\bz_S \). For each \( k\notin S \), the optimality conditions also require
\( |\x_k\T(\Y-\X_S\hat{\bbeta}_S)|\le n\lambda \). Since \( \Y-\X_S\hat{\bbeta}_S=(\I-\M_S)\be+n\lambda\X_S(\X_S\T\X_S)^{-1}\bz_S \),
where \( \M_S \) projects onto \( \C(\X_S) \), this reads
\[
\Bigl|\,\underbrace{\x_k\T\X_S(\X_S\T\X_S)^{-1}\bz_S}_{a_k}+\frac{\x_k\T(\I-\M_S)\be}{n\lambda}\,\Bigr|\le1 .
\]{#eq-hd-dual-feasibility}

The first term does not involve the noise at all. The vector \( (\X_S\T\X_S)^{-1}\X_S\T\x_k \) holds the least
squares coefficients of \( \x_k \) on the active columns, and \( a_k \) is their sum weighted by the true signs.
If an inactive column is well explained by the active ones, with coefficients that line up with the signs of
\( \bbeta_S \), then \( |a_k| \) is large and the condition is hard to meet.

::: {#def-hd-irrepresentable}
[Irrepresentable and beta-min conditions]

Let \( \X_S \) have full column rank. The design satisfies the **irrepresentable condition** for \( \bbeta \)
with margin \( \eta\in(0,1] \) if
\[
\max_{k\notin S}\ \bigl|\x_k\T\X_S(\X_S\T\X_S)^{-1}\operatorname{sign}(\bbeta_S)\bigr|\le1-\eta .
\]
A **beta-min condition** is a lower bound on \( \beta_{\min}=\min_{j\in S}|\beta_j| \).
:::

The condition depends on the design and on the signs of the true coefficients, but not on their
sizes. It holds with \( \eta=1 \) when the inactive columns are orthogonal to the active ones. It is much
stronger than the compatibility condition. van de Geer and Bühlmann (2009) compare the two and show
that a uniform version of the irrepresentable condition implies compatibility; the example below shows
that the converse fails.

## Sufficient conditions

::: {#thm-hd-support}
[Sign recovery by the lasso]

Suppose \( \lambda_{\min}(\hat{\bSigma}_{SS})\ge C_{\min}>0 \) and the irrepresentable condition holds with margin
\( \eta \). Let \( \lambda>0 \), put \( g=\lambda\norm{\hat{\bSigma}_{SS}^{-1}\operatorname{sign}(\bbeta_S)}_\infty \), and suppose
\( \beta_{\min}>g \). Consider the events
\[
\mathcal{A}=\Bigl\{\max_{k\notin S}\frac{|\x_k\T(\I-\M_S)\be|}{n}<\eta\lambda\Bigr\},
\]
\[
\mathcal{B}=\Bigl\{\bigl\|(\X_S\T\X_S)^{-1}\X_S\T\be\bigr\|_\infty<\beta_{\min}-g\Bigr\}.
\]

::: {.enumerate options="label=(\alph*)"}
1. On \( \mathcal{A}\cap\mathcal{B} \) the lasso estimate is unique and \( \operatorname{sign}(\hat{\bbeta}_\lambda)=\operatorname{sign}(\bbeta) \). In
   particular its support is exactly \( S \).

2. The probabilities of the complements satisfy
   \[
   \Pr(\mathcal{A}^c)\le2(p-s)\exp\Bigl\{-\frac{n\eta^2\lambda^2}{2\sigma^2}\Bigr\},\qquad
   \Pr(\mathcal{B}^c)\le2s\exp\Bigl\{-\frac{nC_{\min}(\beta_{\min}-g)^2}{2\sigma^2}\Bigr\}.
   \]

3. Consequently, if
   \[
   \lambda=\frac{\sigma}{\eta}\sqrt{\frac{2\log\{4(p-s)/\delta\}}{n}}\quad\text{and}\quad
   \beta_{\min}\ge g+\sigma\sqrt{\frac{2\log(4s/\delta)}{nC_{\min}}},
   \]
   the lasso recovers the signs of \( \bbeta \) with probability at least \( 1-\delta \).
:::

:::

::: {.proof}
(a) Construct a candidate and check that it is optimal. This is the *primal–dual witness* argument.
Let \( \bz_S=\operatorname{sign}(\bbeta_S) \), let \( \check{\bbeta} \) be zero off \( S \) with \( \check{\bbeta}_S \) given by the right side of
@eq-hd-oracle-lasso, and let \( \check{\bz}=\X\T(\Y-\X\check{\bbeta})/(n\lambda) \).

*On \( S \).* By construction \( \X_S\T(\Y-\X_S\check{\bbeta}_S)=n\lambda\bz_S \), so \( \check{\bz}_S=\bz_S \). By
@eq-hd-oracle-lasso, for \( j\in S \),
\( |\check\beta_j-\beta_j|\le\|(\X_S\T\X_S)^{-1}\X_S\T\be\|_\infty+g<\beta_{\min}\le|\beta_j| \) on \( \mathcal{B} \). So
\( \check\beta_j \) is nonzero with the sign of \( \beta_j \), and \( \check z_j=\operatorname{sign}(\check\beta_j) \).

*Off \( S \).* Since \( \X\bbeta=\X_S\bbeta_S\in\C(\X_S) \), \( (\I-\M_S)\Y=(\I-\M_S)\be \), and
\[
\Y-\X_S\check{\bbeta}_S=\Y-\M_S\Y+n\lambda\X_S(\X_S\T\X_S)^{-1}\bz_S
=(\I-\M_S)\be+n\lambda\X_S(\X_S\T\X_S)^{-1}\bz_S .
\]
So for \( k\notin S \), \( \check z_k=a_k+\x_k\T(\I-\M_S)\be/(n\lambda) \), with \( |a_k|\le1-\eta \) by the irrepresentable
condition. On \( \mathcal{A} \), \( |\check z_k|<1-\eta+\eta=1 \).

So \( \check{\bbeta} \) and \( \check{\bz} \) satisfy @eq-hd-optimality, and \( \check{\bbeta} \) is a lasso estimate, because the
conditions are sufficient (@thm-shr-lasso-kkt(b)). Its support is \( S \), \( \X_S \) has full column rank because
\( C_{\min}>0 \), and \( |\check z_k|<1 \) off \( S \). By @lem-hd-lasso-unique(b) it is the only lasso estimate.

(b) \( \x_k\T(\I-\M_S)\be/n \) is normal with mean zero and variance
\( \sigma^2\norm{(\I-\M_S)\x_k}^2/n^2\le\sigma^2\norm{\x_k}^2/n^2=\sigma^2/n \). The \( j \)th entry of
\( (\X_S\T\X_S)^{-1}\X_S\T\be \) is normal with mean zero and variance
\( \sigma^2[(\X_S\T\X_S)^{-1}]_{jj}\le\sigma^2/(nC_{\min}) \), by @thm-mat-extremal-rayleigh applied to
\( (\X_S\T\X_S)^{-1} \), whose largest eigenvalue is \( 1/(nC_{\min}) \) at most. Apply @lem-hd-gaussian-max with
\( t=\eta\lambda \) and \( t=\beta_{\min}-g \). The events in the lemma use \( > \), and the events here use \( \ge \),
which have the same probability for normal variables with positive variance (and probability zero for
those with zero variance).

(c) With the stated \( \lambda \), the first bound in (b) equals \( \delta/2 \). With the stated \( \beta_{\min} \), the second
is at most \( \delta/2 \).
:::

The three hypotheses ask that the active columns be well conditioned (\( C_{\min} \)), that the inactive ones not
be representable by them (\( \eta \)), and that the coefficients stand out from noise and shrinkage (\( \beta_{\min} \)).
Since \( g\le\lambda\sqrt s/C_{\min} \), the last asks in the worst case for every nonzero coefficient to exceed a
multiple of \( \sigma\sqrt{s\log p/n} \). The theorem in this form is due to Wainwright (2009).

## The condition is also necessary

::: {#prp-hd-irrepresentable-necessary}
[Failure of the irrepresentable condition]

Let \( \X_S \) have full column rank, and suppose that for some \( k\notin S \)
\[
|a_k|=\bigl|\x_k\T\X_S(\X_S\T\X_S)^{-1}\operatorname{sign}(\bbeta_S)\bigr|>1 .
\]
Then the probability that for some \( \lambda>0 \) some lasso estimate has
\( \operatorname{sign}(\hat{\bbeta}_\lambda)=\operatorname{sign}(\bbeta) \) is at most \( \tfrac12 \).
:::

::: {.proof}
Suppose a lasso estimate at some \( \lambda>0 \) has the signs of \( \bbeta \). Its support is \( S \), so by the calculation
leading to @eq-hd-dual-feasibility, \( |a_k+W_k/(n\lambda)|\le1 \), where \( W_k=\x_k\T(\I-\M_S)\be \). Since
\( |a_k|>1 \), this forces \( \operatorname{sign}(a_k)\,W_k/(n\lambda)\le1-|a_k|<0 \), so \( \operatorname{sign}(a_k)W_k<0 \). The event in question
is therefore contained in \( \{\operatorname{sign}(a_k)W_k<0\} \), which does not depend on \( \lambda \). \( W_k \) is normal with mean
zero, so this event has probability \( \frac12 \) if \( (\I-\M_S)\x_k\ne\bzero \), and probability zero otherwise.
:::

The argument is due to Zhao and Yu (2006), who showed in addition that the condition with
\( \le1 \) is necessary for sign consistency. The bound \( \frac12 \) does not improve with more data. Sampling more
observations cannot repair a design in which an inactive variable imitates the active ones with
coefficients aligned with their signs. In the noiseless limit the lasso *never* recovers the signs (@exr-hd-noiseless).

::: {#exm-hd-irrepresentable}
[An inactive variable that looks active]

Let \( x_1,x_2 \) be independent standard normal regressors, let \( x_3 \) have correlation \( r \) with each of them,
and let \( x_4,\dots,x_{500} \) be independent noise variables, with \( Y=0.5x_1+0.5x_2+\varepsilon \) and \( \sigma=1 \). The
population coefficients of \( x_3 \) on \( (x_1,x_2) \) are \( (r,r) \), so \( a_3=2r \) (@exr-hd-irrepresentable-example). With
\( r=0.3 \) the condition holds with margin \( 0.4 \); with the moderate correlation \( r=0.6 \) it fails. Both
population designs have \( \lambda_{\min}(\bSigma)>0 \), so both are benign for prediction.

[Figure 28.4.1](#fig-hd-support-recovery)(a) shows, for each sample size, the fraction of \( 200 \) data sets
(fixed design, new errors) for which *some* penalty on a fine grid gives exactly the right signs. When the
condition holds, this rises from \( 0.12 \) at \( n=50 \) to \( 0.66 \) at \( n=100 \) and
\( 1.000 \) at \( n=800 \). When it fails, it rises only slowly, from \( 0.03 \) at
\( n=100 \) to \( 0.155 \) at \( n=800 \), and it can never exceed the bound \( \frac12 \) of
@prp-hd-irrepresentable-necessary, however large \( n \) is. The sample irrepresentable
quantities at \( n=800 \) are \( 0.66 \) and \( 1.21 \).

Panel (b) shows why. Along the lasso path for one data set with \( r=0.6 \) and \( n=800 \), the inactive \( x_3 \) is
the first variable to enter, because its inner product \( \x_3\T\X\bbeta/n \) with the signal (the mean of
the \( \x_3\T\Y/n \) that decides which variable enters first) is \( 0.62 \) in this sample (\( 0.6 \) in the
population) and exceeds that of \( x_1 \) and \( x_2 \) (\( 0.5 \) each). It never leaves. The lasso
finds \( x_1 \) and \( x_2 \) later, but always together with \( x_3 \).
:::

::: {when-format="html"}
![**Figure 28.4.1.** Support recovery with \( p=500 \), two active regressors and an inactive regressor
correlated \( r \) with both. (a) Fraction of data sets in which some penalty gives exactly the right signs,
against \( n \). (b) The lasso path of one data set with \( r=0.6 \), \( n=800 \): the inactive \( x_3 \) enters first
and stays; grey lines are the noise variables.](support_recovery.svg){#fig-hd-support-recovery width=100%}
:::

::: {when-format="pdf"}
![Support recovery with \( p=500 \), two active regressors and an inactive regressor
correlated \( r \) with both. (a) Fraction of data sets in which some penalty gives exactly the right signs,
against \( n \). (b) The lasso path of one data set with \( r=0.6 \), \( n=800 \): the inactive \( x_3 \) enters first
and stays; grey lines are the noise variables.](support_recovery.pdf){width=100%}
:::

```{.python .run #cell-support-recovery-small}
import numpy as np

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

def make_design(rng, n, p, r):
    """x_1, x_2 independent; x_3 correlated r with each; x_4, ..., x_p independent noise."""
    Z = rng.normal(size=(n, p))
    Z[:, 2] = r * Z[:, 0] + r * Z[:, 1] + np.sqrt(1 - 2 * r ** 2) * Z[:, 2]
    return Z / np.sqrt(np.sum(Z ** 2, axis=0) / n)  # ||x_j||^2 = n

def irrepresentable(X, S, signs):
    """max over j outside S of |x_j^T X_S (X_S^T X_S)^{-1} sign(beta_S)|."""
    A = X[:, S]
    w = A @ np.linalg.solve(A.T @ A, signs)
    out = np.delete(X.T @ w, S)
    return np.max(np.abs(out))

def recovers_signs(X, Y, beta, lams):
    """For each column of Y: does some penalty in lams give exactly sign(beta)?"""
    ok = np.zeros(Y.shape[1], dtype=bool)
    B = None
    for lam in lams:                                 # from large to small, warm starts
        B = lasso_cd(X, Y, lam, B)
        ok |= np.all(np.sign(B) == np.sign(beta)[:, None], axis=0)
    return ok

rng = np.random.default_rng(2840)
n, p, sigma = 100, 200, 1.0
beta = np.zeros(p)
beta[:2] = 0.5
S = np.array([0, 1])
lams = np.geomspace(1.0, 0.1, 30)
for r in [0.3, 0.6]:
    X = make_design(rng, n, p, r)
    Y = (X @ beta)[:, None] + sigma * rng.normal(size=(n, 100))
    print(f"r = {r}: irrepresentable quantity {irrepresentable(X, S, np.ones(2)):.2f},"
          f" sign recovery in {recovers_signs(X, Y, beta, lams).mean():.2f} of 100 data sets")
```

## Screening instead of selection

A weaker goal, often enough in practice, is **screening**: the selected set should *contain* \( S \), with false
variables left for a second stage to remove. It needs only the restricted eigenvalue and a beta-min condition.

::: {#prp-hd-screening}
[Variable screening]

Suppose \( \lambda\ge2\lambda_0 \), the event \( \mathcal{T} \) of @eq-hd-lambda0 occurs, \( \kappa^2(S)>0 \), and
\( \beta_{\min}>3\lambda\sqrt s/\kappa^2(S) \). Then every lasso estimate has \( \hat S\supseteq S \).
:::

::: {.proof}
If \( j\in S \) and \( \hat\beta_{\lambda,j}=0 \), then \( |\hat\beta_{\lambda,j}-\beta_j|=|\beta_j|\ge\beta_{\min} \). By
@thm-hd-lasso-estimation(b), \( |\hat\beta_{\lambda,j}-\beta_j|\le\norm{\hat{\bbeta}_\lambda-\bbeta}\le3\lambda\sqrt s/\kappa^2(S)<\beta_{\min} \),
a contradiction.
:::

In @exm-hd-irrepresentable with \( r=0.6 \), screening succeeds where selection fails: \( x_1 \) and \( x_2 \) are
found, together with \( x_3 \). Second stages build on this. The **thresholded lasso** discards small coefficients and
needs no irrepresentable condition (@exr-hd-thresholded). The **adaptive lasso** of Zou (2006) refits with
penalty weights \( 1/|\hat\beta_j| \) from a first stage, and **stability selection** (Meinshausen and Bühlmann 2010)
keeps the variables selected in most subsamples. All need beta-min conditions; none can find a coefficient
below the noise level.

::: {.warning}
The lasso's selected set is a random set that depends on small features of the data, especially
among correlated regressors. Reporting it as "the important variables" without an assessment of its
stability invites overinterpretation. [Section 28.5](05-debiased-inference.html) shows what goes wrong
if the selected set is then treated as fixed.
:::

## Exercises

### A. Check your understanding

::: {#exr-hd-irrepresentable-example}
[A1]

In the population version of @exm-hd-irrepresentable, with \( \bSigma_{SS}=\I_2 \) and
\( \bSigma_{3S}=(r,r) \), compute \( \bSigma_{3S}\bSigma_{SS}^{-1}\operatorname{sign}(\bbeta_S) \) for \( \bbeta_S=(0.5,0.5) \) and for
\( \bbeta_S=(0.5,-0.5) \). For which signs does the irrepresentable condition fail?
:::

::: {.solution}
\( \bSigma_{3S}\bSigma_{SS}^{-1}=(r,r) \). With signs \( (1,1) \) the quantity is \( 2r \), which exceeds \( 1 \) for \( r>\frac12 \). With
signs \( (1,-1) \) it is \( 0 \), so the condition holds with \( \eta=1 \). The trouble arises only when \( x_3 \) imitates the
*combination* \( x_1+x_2 \) that appears in the signal.
:::

### B. Practice

::: {#exr-hd-equicorrelated-irrepresentable}
[B1]

Let \( \bSigma=(1-\rho)\I+\rho\bone\bone\T \) with \( 0\le\rho<1 \). Show that
\( \bSigma_{kS}\bSigma_{SS}^{-1}\bz_S=\rho\,(\bone\T\bz_S)/\{1+(s-1)\rho\} \) for \( k\notin S \), and conclude that the population
irrepresentable condition holds for every sign pattern, with margin at least \( (1-\rho)/\{1+(s-1)\rho\} \).
:::

::: {#exr-hd-thresholded}
[B2]

*(Thresholded lasso.)* Under the assumptions of @thm-hd-lasso-estimation(b), let \( B=3\lambda\sqrt s/\kappa^2(S) \)
and suppose \( \beta_{\min}>2B \). Show that \( \{j:|\hat\beta_{\lambda,j}|>B\}=S \), and that the signs of the retained
coefficients are correct. No irrepresentable condition is needed.
:::

::: {.solution}
Every coordinate satisfies \( |\hat\beta_{\lambda,j}-\beta_j|\le\norm{\hat{\bbeta}_\lambda-\bbeta}\le B \). For \( j\notin S \), \( |\hat\beta_{\lambda,j}|\le B \),
so \( j \) is discarded. For \( j\in S \), \( |\hat\beta_{\lambda,j}|\ge|\beta_j|-B>2B-B=B \), so \( j \) is kept, and
\( \hat\beta_{\lambda,j} \) lies within \( B<|\beta_j| \) of \( \beta_j \), so it has the same sign. The price is the stronger beta-min
condition \( \beta_{\min}>6\lambda\sqrt s/\kappa^2(S) \) and a threshold that depends on the unknown \( \kappa \) and \( s \).
:::

### C. Going deeper

::: {#exr-hd-noiseless}
[C1]

Let \( \sigma=0 \), so \( \Y=\X\bbeta \), and let \( \X_S \) have full column rank. Show that if
\( \max_{k\notin S}|a_k|<1 \), then every \( \lambda \) with \( 0<\lambda<\beta_{\min}/\norm{\hat{\bSigma}_{SS}^{-1}\operatorname{sign}(\bbeta_S)}_\infty \)
gives a unique lasso estimate with the correct signs, and that if \( |a_k|>1 \) for some \( k\notin S \), then no \( \lambda>0 \) does.
:::

::: {.solution}
With \( \be=\bzero \), the events \( \mathcal{A} \) and \( \mathcal{B} \) of @thm-hd-support hold whenever \( \eta>0 \) and \( \beta_{\min}>g \), which
is the stated range of \( \lambda \), with \( \eta=1-\max_k|a_k| \). Part (a) of the theorem applies, and its proof uses no
probability. Conversely, if some lasso estimate had the correct signs, @eq-hd-dual-feasibility with \( \be=\bzero \) would
require \( |a_k|\le1 \) for all \( k\notin S \). So with perfect data, the irrepresentable condition (in its weak form) is
necessary and, with strict inequality, sufficient.
:::
