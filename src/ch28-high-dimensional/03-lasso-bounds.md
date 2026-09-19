# Prediction and estimation bounds for the lasso

This section proves that the lasso works when \( p>n \), in the form of explicit inequalities valid for every
\( n \) and \( p \). The proofs use three ingredients: the lasso minimizes its criterion, the largest of \( p \) normal
variables is small, and the conditions of [Section 28.2](02-restricted-eigenvalues.html) hold.

## Setting

Throughout this section and the next two, \( \Y=\X\bbeta+\be \) with \( \X \) a fixed \( n\times p \) matrix
whose columns satisfy \( \norm{\x_j}^2=n \), \( \be\sim\Normal_n(\bzero,\sigma^2\I) \), and \( p \) possibly much
larger than \( n \). The support of \( \bbeta \) is \( S \), with \( s=|S| \). A **lasso estimate** is any
\[
\hat{\bbeta}_\lambda\in\argmin_{\bb\in\Real^p}\Bigl\{\frac1{2n}\norm{\Y-\X\bb}^2+\lambda\norm{\bb}_1\Bigr\},
\qquad \lambda>0 .
\]{#eq-hd-lasso}

This is the lasso of @def-shr-lasso with the penalty measured per observation. Other scalings of the
criterion only rename \( \lambda \). By @thm-shr-lasso-kkt, a minimizer exists,
all minimizers share the fitted values
\( \X\hat{\bbeta}_\lambda \) and the norm \( \norm{\hat{\bbeta}_\lambda}_1 \), and the optimality conditions of
@def-shr-lasso are necessary and sufficient. When \( p>n \) the minimizer need not be unique, and
everything below holds for every minimizer. In the present scaling they say that
\( \bb \) is a lasso estimate iff
\[
\begin{gathered}
\frac1n\X\T(\Y-\X\bb)=\lambda\bz\quad\text{with}\\
z_j=\operatorname{sign}(b_j)\ \text{if}\ b_j\ne0,\qquad |z_j|\le1\ \text{if}\ b_j=0 .
\end{gathered}
\]{#eq-hd-optimality}

In words, every column's correlation with the residual is at most \( \lambda \) in absolute value, and exactly
\( \lambda \), with the sign of the coefficient, for every column in the fit. The error of the estimate is
written \( \boldsymbol{\Delta}=\hat{\bbeta}_\lambda-\bbeta \).

The benchmark is the **oracle** estimator, least squares on the columns in \( S \), which could be computed
if \( S \) were known. If \( \X_S \) has full column rank and \( \M_S \) is the projection onto \( \C(\X_S) \), its fitted
values are \( \M_S\Y=\X\bbeta+\M_S\be \), so its in-sample prediction error is
\( \E\norm{\M_S\be}^2/n=\sigma^2\tr(\M_S)/n=\sigma^2s/n \) (@prp-proj-trace-rank). No estimator
that does not know \( S \) can be expected to do this well. The question is how much is lost.

## The basic inequality

::: {#lem-hd-basic}
[Basic inequality]

Every lasso estimate satisfies
\[
\frac1{2n}\norm{\X\boldsymbol{\Delta}}^2\le\frac1n\be\T\X\boldsymbol{\Delta}+\lambda\bigl(\norm{\bbeta}_1-\norm{\hat{\bbeta}_\lambda}_1\bigr).
\]{#eq-hd-basic}

:::

::: {.proof}
The criterion at \( \hat{\bbeta}_\lambda \) is at most its value at \( \bbeta \):
\[
\frac1{2n}\norm{\Y-\X\hat{\bbeta}_\lambda}^2+\lambda\norm{\hat{\bbeta}_\lambda}_1\le\frac1{2n}\norm{\Y-\X\bbeta}^2+\lambda\norm{\bbeta}_1 .
\]
Since \( \Y-\X\bbeta=\be \) and \( \Y-\X\hat{\bbeta}_\lambda=\be-\X\boldsymbol{\Delta} \), the left side contains
\( \norm{\be-\X\boldsymbol{\Delta}}^2=\norm{\be}^2-2\be\T\X\boldsymbol{\Delta}+\norm{\X\boldsymbol{\Delta}}^2 \). Cancelling \( \norm{\be}^2/(2n) \) and
rearranging gives @eq-hd-basic.
:::

The inequality compares a quadratic in \( \boldsymbol{\Delta} \) on the left with a linear noise term and a
difference of \( \ell_1 \) norms on the right. By Hölder's inequality,
\( |\be\T\X\boldsymbol{\Delta}|/n\le\norm{\X\T\be/n}_\infty\norm{\boldsymbol{\Delta}}_1 \). So the noise enters only through the largest
correlation between the errors and a column of \( \X \),
\[
\norm{\X\T\be/n}_\infty=\max_{j\le p}\frac{|\x_j\T\be|}{n}.
\]
The strategy is to choose \( \lambda \) large enough to dominate this quantity with high probability.

## The noise event

::: {#lem-hd-gaussian-max}
[Maximum of normal variables]

Let \( W_1,\dots,W_m \) be normal random variables, not necessarily independent, with mean zero and
variances at most \( v \). For every \( t>0 \),
\[
\Pr\Bigl(\max_{j\le m}|W_j|>t\Bigr)\le2m\exp\Bigl(-\frac{t^2}{2v}\Bigr).
\]{#eq-hd-gaussian-max}

In particular \( \max_j|W_j|\le\sqrt{2v\log(2m/\delta)} \) with probability at least \( 1-\delta \).
:::

::: {.proof}
Let \( W\sim\Normal(0,v_j) \) with \( 0<v_j\le v \). For \( u>0 \), the function \( e^{uw} \) is increasing, so
Markov's inequality (\( \Pr(U\ge a)\le\E U/a \) for \( U\ge0 \)) and the normal moment generating function (@thm-mvn-mgf)
give
\[
\Pr(W>t)=\Pr(e^{uW}\ge e^{ut})\le e^{-ut}\,\E e^{uW}=\exp\bigl(\tfrac12u^2v_j-ut\bigr).
\]
The exponent is smallest at \( u=t/v_j \), where it equals \( -t^2/(2v_j)\le-t^2/(2v) \). By symmetry
\( \Pr(|W|>t)\le2e^{-t^2/(2v)} \), and this holds trivially when \( v_j=0 \). The union bound (@thm-mc-bonferroni(a))
over \( j=1,\dots,m \) gives @eq-hd-gaussian-max. Setting the right side equal to
\( \delta \) and solving for \( t \) gives the last statement.
:::

Each \( \x_j\T\be/n \) is normal with mean zero and variance \( \sigma^2\norm{\x_j}^2/n^2=\sigma^2/n \). For
\( 0<\delta<1 \) put
\[
\lambda_0=\sigma\sqrt{\frac{2\log(2p/\delta)}{n}},\qquad
\mathcal{T}=\Bigl\{\norm{\X\T\be/n}_\infty\le\lambda_0\Bigr\}.
\]{#eq-hd-lambda0}

By @lem-hd-gaussian-max, \( \Pr(\mathcal{T})\ge1-\delta \), whatever the correlations among the columns. The
threshold grows only like \( \sqrt{\log p} \), which is why \( p\gg n \) is not hopeless.

## Prediction bounds

::: {#thm-hd-lasso-prediction}
[Prediction bounds for the lasso]

Let \( \lambda\ge2\lambda_0 \), and suppose the event \( \mathcal{T} \) of @eq-hd-lambda0 occurs. Then every lasso
estimate satisfies:

::: {.enumerate options="label=(\alph*)"}
1. *(Cone.)* \( \dfrac1n\norm{\X\boldsymbol{\Delta}}^2+\lambda\norm{\boldsymbol{\Delta}_{S^c}}_1\le3\lambda\norm{\boldsymbol{\Delta}_S}_1 \). In
   particular \( \boldsymbol{\Delta}\in\mathcal{K}(S) \).

2. *(Slow rate.)* \( \dfrac1n\norm{\X\boldsymbol{\Delta}}^2\le3\lambda\norm{\bbeta}_1 \), with no condition on \( \X \).

3. *(Fast rate.)* If \( \phi^2(S)>0 \), then \( \dfrac1n\norm{\X\boldsymbol{\Delta}}^2\le\dfrac{9\lambda^2s}{\phi^2(S)} \).
:::

Consequently, with \( \lambda=2\lambda_0 \), with probability at least \( 1-\delta \),
\[
\begin{aligned}
\frac1n\norm{\X(\hat{\bbeta}_\lambda-\bbeta)}^2\le
\min\Bigl\{&6\sigma\norm{\bbeta}_1\sqrt{\frac{2\log(2p/\delta)}{n}},\\
&\frac{72\,\sigma^2}{\phi^2(S)}\,\frac{s\log(2p/\delta)}{n}\Bigr\}.
\end{aligned}
\]{#eq-hd-rates}

:::

::: {.proof}
On \( \mathcal{T} \), Hölder's inequality gives \( |\be\T\X\boldsymbol{\Delta}|/n\le\lambda_0\norm{\boldsymbol{\Delta}}_1\le(\lambda/2)\norm{\boldsymbol{\Delta}}_1 \).

(a) Because \( \bbeta \) vanishes off \( S \),
\( \norm{\hat{\bbeta}_\lambda}_1=\norm{\bbeta_S+\boldsymbol{\Delta}_S}_1+\norm{\boldsymbol{\Delta}_{S^c}}_1\ge\norm{\bbeta}_1-\norm{\boldsymbol{\Delta}_S}_1+\norm{\boldsymbol{\Delta}_{S^c}}_1 \).
Substituting both bounds in @eq-hd-basic, and writing \( \norm{\boldsymbol{\Delta}}_1=\norm{\boldsymbol{\Delta}_S}_1+\norm{\boldsymbol{\Delta}_{S^c}}_1 \),
\[
\begin{aligned}
\frac1{2n}\norm{\X\boldsymbol{\Delta}}^2&\le\frac\lambda2\bigl(\norm{\boldsymbol{\Delta}_S}_1+\norm{\boldsymbol{\Delta}_{S^c}}_1\bigr)
+\lambda\bigl(\norm{\boldsymbol{\Delta}_S}_1-\norm{\boldsymbol{\Delta}_{S^c}}_1\bigr)\\
&=\frac{3\lambda}2\norm{\boldsymbol{\Delta}_S}_1-\frac\lambda2\norm{\boldsymbol{\Delta}_{S^c}}_1 .
\end{aligned}
\]
Multiplying by \( 2 \) gives the inequality. Its left side is at least \( \lambda\norm{\boldsymbol{\Delta}_{S^c}}_1 \), so
\( \norm{\boldsymbol{\Delta}_{S^c}}_1\le3\norm{\boldsymbol{\Delta}_S}_1 \).

(b) Use instead \( \norm{\boldsymbol{\Delta}}_1\le\norm{\hat{\bbeta}_\lambda}_1+\norm{\bbeta}_1 \) in @eq-hd-basic:
\[
\begin{aligned}
\frac1{2n}\norm{\X\boldsymbol{\Delta}}^2&\le\frac\lambda2\bigl(\norm{\hat{\bbeta}_\lambda}_1+\norm{\bbeta}_1\bigr)
+\lambda\norm{\bbeta}_1-\lambda\norm{\hat{\bbeta}_\lambda}_1\\
&\le\frac{3\lambda}2\norm{\bbeta}_1 .
\end{aligned}
\]

(c) Write \( a=\norm{\X\boldsymbol{\Delta}}/\sqrt n \). By (a), \( \boldsymbol{\Delta}\in\mathcal{K}(S) \), so @def-hd-re gives
\( \phi^2(S)\le s\,a^2/\norm{\boldsymbol{\Delta}_S}_1^2 \), that is, \( \norm{\boldsymbol{\Delta}_S}_1\le\sqrt s\,a/\phi(S) \). (If
\( \boldsymbol{\Delta}=\bzero \) there is nothing to prove.) Then (a) gives \( a^2\le3\lambda\norm{\boldsymbol{\Delta}_S}_1\le3\lambda\sqrt s\,a/\phi(S) \),
so \( a\le3\lambda\sqrt s/\phi(S) \).

For @eq-hd-rates, substitute \( \lambda=2\lambda_0 \) in (b) and (c), and use \( \Pr(\mathcal{T})\ge1-\delta \).
:::

The **slow rate** needs nothing of \( \X \): the lasso predicts consistently whenever
\( \norm{\bbeta}_1\sqrt{\log p/n}\to0 \), which allows \( p \) to grow exponentially in \( n \). Greenshtein and Ritov (2004)
called this *persistence*. The **fast rate**, of order \( \sigma^2s\log p/n \), loses against the oracle's \( \sigma^2s/n \)
a factor \( \log p \) and the constant \( 72/\phi^2(S) \). The logarithm cannot be removed: no estimator that does not
know \( S \) has worst-case error below a multiple of \( \sigma^2s\log(p/s)/n \) under mild conditions on the design
(Raskutti, Wainwright and Yu 2011; not proved here). The constants are not sharp.

::: {.idea}
The proof has two separate parts. A *probabilistic* part shows that the noise is small in the
\( \ell_\infty \) norm dual to the penalty, at the level \( \lambda_0\asymp\sigma\sqrt{\log p/n} \). A *deterministic*
part says that on this event the estimate's error lies in a cone and is controlled by the geometry of
\( \X \) on that cone. Changing the noise distribution changes only the first part. Changing the penalty
changes the norm and the cone.
:::

## Estimation bounds

Prediction error measures \( \boldsymbol{\Delta} \) through \( \X \). Bounds on the coefficients themselves need the
design to separate the directions of the cone, which is what the restricted eigenvalue measures.

::: {#thm-hd-lasso-estimation}
[Estimation bounds for the lasso]

Let \( \lambda\ge2\lambda_0 \), and suppose \( \mathcal{T} \) occurs. Then every lasso estimate satisfies:

::: {.enumerate options="label=(\alph*)"}
1. if \( \phi^2(S)>0 \), then \( \norm{\hat{\bbeta}_\lambda-\bbeta}_1\le\dfrac{4\lambda s}{\phi^2(S)} \);

2. if \( \kappa^2(S)>0 \), then \( \norm{\hat{\bbeta}_\lambda-\bbeta}\le\dfrac{3\lambda\sqrt s}{\kappa^2(S)} \).
:::

With \( \lambda=2\lambda_0 \), both hold with probability at least \( 1-\delta \), and then
\( \norm{\hat{\bbeta}_\lambda-\bbeta}_1\le8\lambda_0s/\phi^2(S) \) and \( \norm{\hat{\bbeta}_\lambda-\bbeta}\le6\lambda_0\sqrt s/\kappa^2(S) \).
:::

::: {.proof}
Keep the notation \( a=\norm{\X\boldsymbol{\Delta}}/\sqrt n \) of the previous proof.

(a) Add \( \lambda\norm{\boldsymbol{\Delta}_S}_1 \) to both sides of @thm-hd-lasso-prediction(a), and use
\( \norm{\boldsymbol{\Delta}_S}_1\le\sqrt s\,a/\phi(S) \) and \( 2xy\le x^2+y^2 \):
\[
a^2+\lambda\norm{\boldsymbol{\Delta}}_1\le4\lambda\norm{\boldsymbol{\Delta}_S}_1\le2\cdot a\cdot\frac{2\lambda\sqrt s}{\phi(S)}
\le a^2+\frac{4\lambda^2s}{\phi^2(S)} .
\]
Cancel \( a^2 \) and divide by \( \lambda \).

(b) By @thm-hd-lasso-prediction(a) and @eq-hd-l1-l2, and since \( \boldsymbol{\Delta}\in\mathcal{K}(S) \), by @def-hd-re,
\[
\kappa^2(S)\norm{\boldsymbol{\Delta}}^2\le a^2\le3\lambda\norm{\boldsymbol{\Delta}_S}_1\le3\lambda\sqrt s\norm{\boldsymbol{\Delta}} .
\]
Dividing by \( \norm{\boldsymbol{\Delta}} \) (when it is not zero) gives the bound.
:::

With \( \lambda\asymp\sigma\sqrt{\log p/n} \), the \( \ell_2 \) error is of order \( \sigma\sqrt{s\log p/n} \). The oracle
has \( \ell_2 \) error of order \( \sigma\sqrt{s/n} \). Again only a logarithm is lost. Bickel, Ritov and Tsybakov
(2009) proved bounds of this form for the lasso and the Dantzig selector together. The compatibility
condition was introduced by van de Geer (see van de Geer and Bühlmann 2009), and Bühlmann and van de Geer (2011, chapter 6) present the prediction
bound in the form given here.

::: {.remark}
[What the theorems assume, and what they do not]

The theorems are *deterministic* on \( \mathcal{T} \); normality entered only through \( \Pr(\mathcal{T})\ge1-\delta \), so errors
with light (sub-Gaussian) tails need the same \( \lambda_0 \) up to a constant. For a random design they apply
conditionally on \( \X \). The choice \( \lambda_0 \) needs \( \sigma \); the *square-root lasso* of Belloni, Chernozhukov and
Wang (2011), which minimizes \( \norm{\Y-\X\bb}/\sqrt n+\lambda\norm{\bb}_1 \), and the *scaled lasso* of Sun and Zhang
(2012) avoid this. In practice \( \lambda \) is chosen by cross-validation
([Chapter 29](../ch29-model-selection/index.html)), without these guarantees.
:::

## The bounds in numbers

::: {#exm-hd-lasso-bounds}
[How tight are the bounds?]

Take the independent Gaussian design of @exm-hd-compatibility-values with \( n=200 \) and \( p=1000 \), five
nonzero coefficients equal to \( \pm1 \), \( \sigma=1 \) and \( \delta=0.05 \). Then
\( \lambda_0=0.3255 \). In \( 200 \) simulated data sets the event \( \mathcal{T} \) occurred in a
fraction \( 0.990 \), more than the guaranteed \( 0.95 \), because the Chernoff bound and the union
bound are both conservative. The cone property of @thm-hd-lasso-prediction(a) held in a fraction
\( 1.000 \) of the data sets.

At the theoretical choice \( \lambda=2\lambda_0=0.651 \) the average prediction error
\( \norm{\X\boldsymbol{\Delta}}^2/n \) is \( 2.3643 \), against the oracle's
\( \sigma^2s/n=0.025 \). The slow-rate bound \( 3\lambda\norm{\bbeta}_1 \) is
\( 9.77 \). With the compatibility constant \( 0.00447 \) of this design, the fast-rate bound
is \( 4269.8 \) (computed from the unrounded values), which is useless here. The same is true of the \( \ell_1 \) bound
\( 2914.8 \), against an actual average \( \ell_1 \) error of \( 3.594 \).

Almost all of the error at \( \lambda=2\lambda_0 \) is shrinkage. The lasso pulls each selected coefficient
towards zero by about \( \lambda \), and \( s\lambda^2=2.12 \) is close to the observed error. Panel (a)
of [Figure 28.3.1](#fig-hd-lasso-bounds) shows the whole path. The smallest average error on the path is
\( 0.1893 \), at \( \lambda\approx0.39\lambda_0 \), well below the theoretical choice. Below that
the lasso admits many false variables and the error rises again, to \( 0.587 \) at
\( \lambda=0.1\lambda_0 \).

Panel (b) holds \( n=200 \) and \( s=5 \) fixed and lets \( p \) grow from \( 100 \) to \( 5000 \) (the first \( p \)
columns of one larger design), with \( \lambda=2\lambda_0 \) for each \( p \). The error grows from
\( 1.818 \) to \( 2.665 \), a factor of \( 1.466 \), close to the ratio
\( 1.472 \) of the values of \( \log(2p/\delta) \). The error tracks the \( \lambda_0^2 \) in the bounds,
which is logarithmic in \( p \). The theorem gives only upper bounds, so this agreement is an observation
about the example, not a consequence of the theorem.
:::

::: {when-format="html"}
![**Figure 28.3.1.** The lasso with \( n=200 \), \( s=5 \), \( \sigma=1 \). (a) Average prediction error
\( \norm{\X(\hat{\bbeta}_\lambda-\bbeta)}^2/n \) and average number of selected variables (divided by \( 100 \))
along the path for \( p=1000 \), against \( \lambda/\lambda_0 \). The vertical line is the theoretical choice
\( \lambda=2\lambda_0 \), the dotted line the oracle error. (b) Prediction error at \( \lambda=2\lambda_0 \) as \( p \)
grows (logarithmic axis), with two standard errors.](lasso_bounds.svg){#fig-hd-lasso-bounds width=100%}
:::

::: {when-format="pdf"}
![The lasso with \( n=200 \), \( s=5 \), \( \sigma=1 \). (a) Average prediction error
\( \norm{\X(\hat{\bbeta}_\lambda-\bbeta)}^2/n \) and average number of selected variables (divided by \( 100 \))
along the path for \( p=1000 \), against \( \lambda/\lambda_0 \). The vertical line is the theoretical choice
\( \lambda=2\lambda_0 \), the dotted line the oracle error. (b) Prediction error at \( \lambda=2\lambda_0 \) as \( p \)
grows (logarithmic axis), with two standard errors.](lasso_bounds.pdf){width=100%}
:::

The bounds are worst-case statements over designs and noise vectors, and for a particular design they
can be very loose. What they get right is the *form* of the dependence on \( s \), \( p \) and \( n \), with one
tuning parameter of order \( \sigma\sqrt{\log p/n} \) that works for all sparse \( \bbeta \). They also expose the cost of
shrinkage, an error of order \( s\lambda^2 \) even when the right variables are selected, which motivates the next
two sections.

```{.python .run #cell-lasso-bounds-small}
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

rng = np.random.default_rng(2832)
n, p, s, sigma, delta = 100, 300, 5, 1.0, 0.05
X = rng.normal(size=(n, p))
X /= np.sqrt(np.sum(X ** 2, axis=0) / n)           # columns scaled to ||x_j||^2 = n
beta = np.zeros(p)
beta[:s] = [1.0, -1.0, 1.0, -1.0, 1.0]
lam0 = sigma * np.sqrt(2 * np.log(2 * p / delta) / n)
E = sigma * rng.normal(size=(n, 50))
Y = (X @ beta)[:, None] + E
print("lambda_0 =", round(lam0, 3))
print("frequency of the event:", np.mean(np.max(np.abs(X.T @ E) / n, axis=0) <= lam0))
for ratio in [2.0, 1.0, 0.5]:
    B = lasso_cd(X, Y, ratio * lam0)
    err = np.sum((X @ (B - beta[:, None])) ** 2, axis=0) / n
    print(f"lambda = {ratio} lambda_0: mean prediction error {err.mean():.3f},"
          f" mean number selected {np.mean(np.sum(B != 0, axis=0)):.1f}")
print("oracle sigma^2 s / n =", sigma ** 2 * s / n, "  slow-rate bound at 2 lambda_0:", 6 * lam0 * s)
```

## Exercises

### A. Check your understanding

::: {#exr-hd-lambda-max}
[A1]

By @thm-shr-lasso-kkt(d), \( \bzero \) is a lasso estimate iff
\( \lambda\ge\lambda_{\max}=\norm{\X\T\Y}_\infty/n \) (check this from @eq-hd-optimality). Show that it is then the
only lasso estimate. Why does the theory above not rule out \( 2\lambda_0>\lambda_{\max} \), and what
happens then?
:::

::: {.solution}
All lasso estimates have the same \( \ell_1 \) norm (@thm-shr-lasso-kkt(c)), and that of
\( \bzero \) is zero; so every lasso estimate is \( \bzero \). If \( 2\lambda_0>\lambda_{\max} \) and
\( \lambda\ge\lambda_{\max} \), the lasso returns \( \bzero \), and the theorems still hold: the error is \( -\bbeta \), and
the bounds then say that \( \bbeta \) is small compared with the noise level. They are not violated, only
uninformative.
:::

::: {#exr-hd-basic-general}
[A2]

Show that for every \( \bb\in\Real^p \) a lasso estimate satisfies
\[
\frac1{2n}\norm{\X(\hat{\bbeta}_\lambda-\bbeta)}^2\le\frac1{2n}\norm{\X(\bb-\bbeta)}^2
+\frac1n\be\T\X(\hat{\bbeta}_\lambda-\bb)+\lambda\bigl(\norm{\bb}_1-\norm{\hat{\bbeta}_\lambda}_1\bigr),
\]
and that \( \bb=\bbeta \) gives @lem-hd-basic. Deduce that on \( \mathcal{T} \), with \( \lambda\ge2\lambda_0 \),
\[
\frac1n\norm{\X(\hat{\bbeta}_\lambda-\bbeta)}^2\le\min_{\bb\in\Real^p}\Bigl\{\frac1n\norm{\X(\bb-\bbeta)}^2+3\lambda\norm{\bb}_1\Bigr\}.
\]
This oracle inequality suits an approximately sparse \( \bbeta \): \( \bb \) can keep its large entries and drop
the small ones.
:::

::: {.solution}
Compare the criterion @eq-hd-lasso at \( \hat{\bbeta}_\lambda \) and at \( \bb \), expanding
\( \norm{\X(\bbeta-\bb)+\be}^2=\norm{\X(\bb-\bbeta)}^2-2\be\T\X(\bb-\bbeta)+\norm{\be}^2 \) and likewise at
\( \hat{\bbeta}_\lambda \); the \( \norm{\be}^2 \) terms cancel. On \( \mathcal{T} \), Hölder's inequality gives
\( |\be\T\X(\hat{\bbeta}_\lambda-\bb)|/n\le\lambda_0\norm{\hat{\bbeta}_\lambda-\bb}_1\le(\lambda/2)(\norm{\hat{\bbeta}_\lambda}_1+\norm{\bb}_1) \).
So the right side is at most \( \norm{\X(\bb-\bbeta)}^2/(2n)+\frac32\lambda\norm{\bb}_1-\frac12\lambda\norm{\hat{\bbeta}_\lambda}_1 \). Drop the
last term, multiply by \( 2 \), and minimize over \( \bb \). At \( \bb=\bbeta \) this is the slow rate of @thm-hd-lasso-prediction(b).
An analogue with \( \lambda^2s/\phi^2 \) in place of \( 3\lambda\norm{\bb}_1 \),
under a compatibility condition, is in Bühlmann and van de Geer (2011, chapter 6).
:::

### B. Practice

::: {#exr-hd-orthonormal-lasso}
[B1]

Let \( p\le n \) and \( \hat{\bSigma}=\I \). Using @thm-shr-lasso-orthonormal, the lasso estimate is
\( \hat\beta_{\lambda,j}=\operatorname{sign}(z_j)(|z_j|-\lambda)_+ \) with \( \bz=\X\T\Y/n=\bbeta+\X\T\be/n \). On \( \mathcal{T} \), with
\( \lambda\ge2\lambda_0 \), show that \( \hat\beta_{\lambda,j}=0 \) for \( j\notin S \), that \( |\hat\beta_{\lambda,j}-\beta_j|\le\lambda+\lambda_0 \) for
\( j\in S \), and hence \( \norm{\X\boldsymbol{\Delta}}^2/n\le\frac94\lambda^2s \). Show also that if every \( |\beta_j| \), \( j\in S \), exceeds
\( \lambda+\lambda_0 \), then \( \norm{\X\boldsymbol{\Delta}}^2/n\ge\frac14\lambda^2s \). So the order \( \lambda^2s \) of @thm-hd-lasso-prediction(c)
is attained.
:::

::: {.solution}
With \( \hat{\bSigma}=\I \), \( \norm{\X\boldsymbol{\Delta}}^2/n=\norm{\boldsymbol{\Delta}}^2 \). For \( j\notin S \), \( |z_j|=|\x_j\T\be|/n\le\lambda_0<\lambda \), so the
estimate is zero. For \( j\in S \), soft thresholding moves \( z_j \) by at most \( \lambda \), and \( |z_j-\beta_j|\le\lambda_0 \), so
\( |\hat\beta_{\lambda,j}-\beta_j|\le\lambda+\lambda_0\le\frac32\lambda \). Summing squares over \( S \) gives the upper bound. If
\( |\beta_j|>\lambda+\lambda_0 \), then \( |z_j|>\lambda \), so \( \hat\beta_{\lambda,j}=z_j-\lambda\operatorname{sign}(z_j) \) with the sign of \( \beta_j \), and
\( |\hat\beta_{\lambda,j}-\beta_j|\ge\lambda-\lambda_0\ge\lambda/2 \). Summing gives the lower bound. The error of order
\( \lambda^2s \) is pure shrinkage.
:::

### C. Going deeper

::: {#exr-hd-dantzig}
[C1]

*(The Dantzig selector.)* Candès and Tao (2007) proposed to estimate \( \bbeta \) by a solution \( \tilde{\bbeta} \)
of: minimize \( \norm{\bb}_1 \) subject to \( \norm{\X\T(\Y-\X\bb)}_\infty/n\le\lambda_0 \). On \( \mathcal{T} \), show that \( \bbeta \) is
feasible, that \( \boldsymbol{\Delta}=\tilde{\bbeta}-\bbeta \) satisfies \( \norm{\boldsymbol{\Delta}_{S^c}}_1\le\norm{\boldsymbol{\Delta}_S}_1 \), and that
\( \norm{\X\boldsymbol{\Delta}}^2/n\le16\lambda_0^2s/\phi^2(S) \).
:::

::: {.solution}
On \( \mathcal{T} \), \( \norm{\X\T(\Y-\X\bbeta)}_\infty/n=\norm{\X\T\be}_\infty/n\le\lambda_0 \), so \( \bbeta \) is feasible and
\( \norm{\tilde{\bbeta}}_1\le\norm{\bbeta}_1 \). The argument of @thm-hd-lasso-prediction(a) gives
\( \norm{\bbeta}_1\ge\norm{\tilde{\bbeta}}_1\ge\norm{\bbeta}_1-\norm{\boldsymbol{\Delta}_S}_1+\norm{\boldsymbol{\Delta}_{S^c}}_1 \), which is the cone
property with constant \( 1 \). Since both \( \tilde{\bbeta} \) and \( \bbeta \) are feasible,
\( \norm{\X\T\X\boldsymbol{\Delta}}_\infty/n\le2\lambda_0 \), and by Hölder
\( \norm{\X\boldsymbol{\Delta}}^2/n=\boldsymbol{\Delta}\T\X\T\X\boldsymbol{\Delta}/n\le2\lambda_0\norm{\boldsymbol{\Delta}}_1\le4\lambda_0\norm{\boldsymbol{\Delta}_S}_1 \). The cone with
constant \( 1 \) lies inside \( \mathcal{K}(S) \), so \( \norm{\boldsymbol{\Delta}_S}_1\le\sqrt s\norm{\X\boldsymbol{\Delta}}/(\sqrt n\,\phi(S)) \). Combining,
\( \norm{\X\boldsymbol{\Delta}}/\sqrt n\le4\lambda_0\sqrt s/\phi(S) \), and squaring gives the bound.
:::
