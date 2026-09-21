# Choosing the smoothing parameter

Four smoothers have appeared: the local constant and local linear estimators of
[Section 43.1](01-kernels.html), the penalized spline of
[Section 43.3](03-penalized-splines.html) and the smoothing spline of
[Section 43.4](04-smoothing-splines.html). Each is linear in \( \y \), each has one
dial, and none has said how to set it. The rules of this section approximate a
risk that cannot be observed: they agree when the sample is large and the errors
independent, and fail in recognisably different ways when it is not.

## Every smoother here is a linear operator

::: {#prp-smo-linear-smoother}
[Linear smoothers]

Call \( \hat{\mathbf f}=\bS\y \) a **linear smoother** with **smoother matrix**
\( \bS \) when \( \bS \) does not depend on \( \y \). The estimators
of @def-smo-kernel, @def-smo-pspline and @thm-smo-smoothing-spline are linear
smoothers, with \( \bS \) given by @eq-smo-equivalent, @eq-smo-smoother
and @eq-smo-reinsch-fit. Suppose \( \E(\Y)=\mathbf f \) and \( \Cov(\Y)=\sigma^2\I \).
Then:

::: {.enumerate options="label=(\alph*)"}
1. *(Moments.)* \( \E(\hat{\mathbf f})=\bS\mathbf f \),
   \( \Cov(\hat{\mathbf f})=\sigma^2\bS\bS\T \), and
   \[
   \E\norm{\hat{\mathbf f}-\mathbf f}^2
   =\norm{(\I-\bS)\mathbf f}^2+\sigma^2\tr(\bS\bS\T) .
   \]{#eq-smo-risk}

2. *(Residuals.)*
   \( \E\norm{\y-\hat{\mathbf f}}^2
   =\norm{(\I-\bS)\mathbf f}^2+\sigma^2\{n-2\tr(\bS)+\tr(\bS\bS\T)\} \), so
   \[
   \hat\sigma^2=\frac{\norm{\y-\hat{\mathbf f}}^2}{n-2\tr(\bS)+\tr(\bS\bS\T)}
   \]{#eq-smo-sigma-hat}

   is unbiased for \( \sigma^2 \) when the bias term is negligible.

3. *(Degrees of freedom.)* If \( \bS \) is symmetric with eigenvalues in
   \( [0,1] \) — as it is for a penalized or smoothing spline, by
   @prp-smo-pspline-ridge(b) — then
   \( \tr(\bS\bS\T)\le\tr(\bS)\le2\tr(\bS)-\tr(\bS\bS\T) \), with equality
   throughout iff \( \bS \) is an orthogonal projection, in which case all three
   equal its rank.

4. *(Shrinkage.)* For a penalized spline, @eq-smo-smoother writes \( \bS_\lambda \)
   in the Demmler–Reinsch basis as a diagonal shrinkage
   \( \bu_j\T\y\mapsto\bu_j\T\y/(1+\lambda e_j) \): smooth directions
   (\( e_j \) small) pass through, rough ones are damped. The smoothing spline has
   the same form with \( e_j \) the eigenvalues of the roughness matrix relative to
   the identity.

5. *(Pointwise bands.)* \( \hat f(x_i)\pm z_{1-\alpha/2}\hat\sigma\sqrt{(\bS\bS\T)_{ii}} \)
   is an interval for \( \E\{\hat f(x_i)\}=(\bS\mathbf f)_i \), not for
   \( f(x_i) \). Its coverage of \( f(x_i) \) is \( 1-\alpha \) only to the extent
   that the bias \( \{(\bS-\I)\mathbf f\}_i \) is small next to
   \( \sigma\sqrt{(\bS\bS\T)_{ii}} \).
:::

:::

::: {.proof}
Linearity is @prp-smo-equivalent-kernel for the kernel estimators,
@prp-smo-pspline-ridge(a) for the penalized spline and @eq-smo-reinsch-fit for the
smoothing spline; in each case the matrix is built from the design, the knots and
the tuning constant alone.

(a) \( \E(\bS\Y)=\bS\mathbf f \) and \( \Cov(\bS\Y)=\bS\Cov(\Y)\bS\T \). For the
risk, \( \hat{\mathbf f}-\mathbf f=(\bS-\I)\mathbf f+\bS\be \), the two terms are
uncorrelated, and \( \E\norm{\bS\be}^2=\sigma^2\tr(\bS\bS\T) \).

(b) \( \y-\hat{\mathbf f}=(\I-\bS)\mathbf f+(\I-\bS)\be \), and
\[
\begin{aligned}
\E\norm{(\I-\bS)\be}^2&=\sigma^2\tr\{(\I-\bS)(\I-\bS)\T\}\\
&=\sigma^2\{n-2\tr(\bS)+\tr(\bS\bS\T)\} .
\end{aligned}
\]

(c) With eigenvalues \( s_j\in[0,1] \) we have \( \tr(\bS)=\sum_js_j \) and
\( \tr(\bS\bS\T)=\sum_js_j^2 \), while \( s_j^2\le s_j\le2s_j-s_j^2 \), the last
inequality because \( (1-s_j)^2\ge0 \). Equality in either needs every \( s_j\in\{0,1\} \), which for a
symmetric matrix means \( \bS \) is idempotent (@thm-proj-sym-idem); its trace is
then its rank (@prp-proj-trace-rank).

(d) is @eq-smo-smoother; for the smoothing spline apply @lem-smo-diagonalize with
\( \B=\I \) and the roughness matrix of @prp-smo-reinsch(b) in place of \( \bP \).

(e) By (a), \( \hat f(x_i) \) has mean \( (\bS\mathbf f)_i \) and standard
deviation \( \sigma\{(\bS\bS\T)_{ii}\}^{1/2} \), so under normal errors the
interval covers that mean with probability \( 1-\alpha \) if \( \hat\sigma \) is replaced by
\( \sigma \). Writing \( f(x_i)=(\bS\mathbf f)_i-\{(\bS-\I)\mathbf f\}_i \) shows
that its coverage of \( f(x_i) \) is the normal probability of an interval
displaced by the bias.
:::

::: {.warning}
Part (e) is the most important caveat in this chapter, and the one most often
ignored. A smoother is *deliberately* biased; that is what the penalty buys. The
interval above treats the bias as zero, so where \( f \) curves most — a peak, a
trough, the boundary — it is centred in the wrong place and undercovers, while on
flat stretches it overcovers. [Section 43.6](06-mixed-model-and-bayes.html) gives
a band that repairs the average of this, and says what "average" means.
:::

## Three rules

The quantity one would like to minimize is the risk @eq-smo-risk, which involves
the unknown \( \mathbf f \). Each rule estimates something close to it.

::: {#thm-smo-lambda}
[Cross-validation, generalized cross-validation and AIC for smoothers]

Let \( \bS_\lambda \) be the smoother matrix, \( \hat{\mathbf f}=\bS_\lambda\y \),
\( \mathrm{df}=\tr(\bS_\lambda) \) and \( s_{ii} \) the diagonal entries of
\( \bS_\lambda \) — the leverage of @prp-proj-leverage and
[Section 20.2](../ch20-residuals-leverage-influence/02-leverage.html), read now
for a smoother instead of a projection.

::: {.enumerate options="label=(\alph*)"}
1. *(The leave-one-out shortcut is exact.)* For a penalized spline and for a
   smoothing spline,
   \[
   \mathrm{CV}(\lambda)=\frac1n\sum_{i=1}^n\{y_i-\hat f^{(-i)}(x_i)\}^2
   =\frac1n\sum_{i=1}^n\Bigl(\frac{y_i-\hat f(x_i)}{1-s_{ii}}\Bigr)^2 ,
   \]{#eq-smo-cv}

   where \( \hat f^{(-i)} \) is the fit computed without case \( i \) at the same
   \( \lambda \). The same identity holds for a local polynomial estimator, with
   \( s_{ii}=\ell_i(x_i) \).

2. *(Generalized cross-validation.)* Replacing each \( s_{ii} \) by the average
   \( \mathrm{df}/n \) gives @def-sel-gcv,
   \[
   \mathrm{GCV}(\lambda)=\frac{n^{-1}\norm{\y-\hat{\mathbf f}}^2}{(1-\mathrm{df}/n)^2} ,
   \]{#eq-smo-gcv}

   which needs only the trace, and satisfies
   \( \mathrm{GCV}=n^{-1}\norm{\y-\hat{\mathbf f}}^2\,(1+2\,\mathrm{df}/n)+O(\mathrm{df}^2/n^2) \).

3. *(Information criteria.)* Under normal errors with \( \sigma^2 \) profiled out,
   @prp-sel-aic-bic gives
   \( \mathrm{AIC}(\lambda)=n\log(\norm{\y-\hat{\mathbf f}}^2/n)+2\,\mathrm{df} \).
   Its small-sample correction (Hurvich, Simonoff and Tsai, 1998) is
   \[
   \mathrm{AIC_c}(\lambda)=\log\Bigl(\frac{\norm{\y-\hat{\mathbf f}}^2}{n}\Bigr)
   +1+\frac{2(\mathrm{df}+1)}{n-\mathrm{df}-2} ,
   \]{#eq-smo-aicc}

   which penalizes large \( \mathrm{df} \) much more heavily and has a pole at
   \( \mathrm{df}=n-2 \).

4. *(What they estimate.)* If \( \Cov(\Y)=\sigma^2\I \), then
   \( \E\{\mathrm{CV}(\lambda)\} \) is, up to \( O(n^{-1}) \) terms, the average
   prediction error at a fresh response, which exceeds \( n^{-1} \) times the
   risk @eq-smo-risk by \( \sigma^2 \). A minimizer of any of the three is therefore aiming at
   @eq-smo-risk, not at \( \lambda \) itself, and nothing here estimates a "true"
   smoothing parameter, since no such thing is assumed to exist.
:::

:::

::: {.proof}
(a) A penalized spline minimizes \( \norm{\y-\B\bgamma}^2+\bgamma\T\boldsymbol{\Omega}\bgamma \)
with \( \boldsymbol{\Omega}=\lambda\bP \), and a smoothing spline minimizes
\( \norm{\y-\bz}^2+\bz\T(\lambda\bP)\bz \), which is the same form with
\( \B=\I \). Both have \( \boldsymbol{\Omega} \) nonnegative definite and the
deleted cross-product matrices nonsingular, so @thm-sel-loocv applies verbatim
and gives @eq-smo-cv. For a local polynomial estimator, fix \( i \) and let
\( \tilde{\y} \) be \( \y \) with \( y_i \) replaced by
\( \hat f^{(-i)}(x_i) \). The weighted least squares problem @eq-smo-localpoly at
\( x=x_i \) with data \( \tilde{\y} \) has the same minimizer as the problem with
case \( i \) deleted. Indeed the row of \( \bU_x \) belonging to case \( i \) is
\( (1,0,\dots,0) \), so the objective with \( \tilde{\y} \) is the deleted objective
plus \( K_h(0)(\tilde y_i-b_0)^2 \); the deleted minimizer makes the first term as
small as possible and the second term zero, so it minimizes the sum. Hence
\( \hat f^{(-i)}(x_i)=\sum_j\ell_j(x_i)\tilde y_j
=\hat f(x_i)+\ell_i(x_i)\{\hat f^{(-i)}(x_i)-y_i\} \); solving for
\( y_i-\hat f^{(-i)}(x_i) \) gives \( \{y_i-\hat f(x_i)\}/\{1-\ell_i(x_i)\} \).

(b) The substitution is the definition of GCV; the expansion is
\( (1-t)^{-2}=1+2t+O(t^2) \) with \( t=\mathrm{df}/n \).

(c) With a normal likelihood and \( \sigma^2 \) replaced by its maximum
likelihood estimate \( \norm{\y-\hat{\mathbf f}}^2/n \), the maximized log
likelihood is \( -\tfrac n2\log(\norm{\y-\hat{\mathbf f}}^2/n)+\text{const} \), and
@prp-sel-aic-bic counts \( \mathrm{df}+1 \) parameters, the \( +1 \) being
\( \sigma^2 \); dropping constants gives AIC. @eq-smo-aicc is the stated
correction, which reduces to AIC\( /n \) plus constants when
\( \mathrm{df}/n\to0 \).

(d) is @prp-sel-cv-target and @def-sel-prediction-error, which do not use the
form of the estimator.
:::

::: {.idea}
All three rules do the same arithmetic: residual sum of squares, inflated by a
factor that grows with the degrees of freedom — \( (1-s_{ii})^{-2} \) case by
case, \( (1-\mathrm{df}/n)^{-2} \) on average, \( \exp\{2\,\mathrm{df}/n\} \)
approximately for AIC. That is why they usually agree, and why they part company
at large \( \mathrm{df} \). Mallows' \( C_p \) (@def-sel-cp) inflates additively,
\( \norm{\y-\hat{\mathbf f}}^2+2\hat\sigma^2\tr(\bS_\lambda) \), but needs a
\( \hat\sigma^2 \) from outside — usually @eq-smo-sigma-hat at a deliberately
undersmoothed fit — and inherits whatever that choice got wrong.
:::

## What they do in practice

::: {#exm-smo-criteria}
[The three criteria on one data set, and over 200 of them]

Take \( n=120 \) equally spaced design points, the
mean function of @exm-smo-bandwidth,
\( \sigma=0.25 \), and a cubic P-spline with
\( K=25 \) knots and a second-order penalty. On one data
set the shortcut @eq-smo-cv agrees with \( 120 \) actual refits to
\( 4.6\times 10^{-15} \), and the three criteria choose
\( 9.83 \), \( 10.29 \) and
\( 9.40 \) degrees of freedom.

Over \( 200 \) replicates the three are hard to tell apart.
Each has median selected df \( 7.82 \), and the
average \( n^{-1}\norm{\hat{\mathbf f}-\mathbf f}^2 \) is
\( 0.00448 \) for CV,
\( 0.00449 \) for GCV and
\( 0.00429 \) for AIC\( _c \). The differences show in
the upper tail: CV and GCV choose more than \( 15 \) degrees of freedom in
\( 2.5 \) per cent of the samples, AIC\( _c \) in
\( 0.5 \) per cent. The criteria are flat to the
right of the minimum
([Figure 43.5.1](05-choosing-lambda.html#fig-smo-criteria)(a)), so an occasional
sample sends the minimum far out, and steep to the left, so gross oversmoothing
is rare.
:::

::: {.warning}
[Three small-sample failures]

**A flat criterion.** Near the minimum, GCV can change by a fraction of a per
cent over a factor of two in \( \lambda \); the choice is then arbitrary, and
reporting the fit at a single \( \lambda \) hides that. Plot the criterion.

**Interpolation.** Since \( (1-\mathrm{df}/n)^{-2} \) becomes large only as
\( \mathrm{df}\to n \), a small sample can give GCV a *local* minimum at the
interpolating end. AIC\( _c \) cannot, because of the pole in @eq-smo-aicc; this
is its main practical advantage.

**Correlated errors.** All three estimate prediction error at a *new,
independent* response. If the errors are serially correlated, a wiggly fit
predicts its neighbours well without describing \( f \) at all, and every
criterion rewards it (@exm-smo-correlated).
:::

::: {#exm-smo-correlated}
[What correlated errors do]

Repeat the simulation with errors from a stationary AR(1) process of the same
marginal standard deviation and autocorrelation
\( \rho=0.7 \) (@def-het-ar1). The median selected
degrees of freedom rise from \( 7.82 \) to
\( 25.40 \) for CV,
\( 24.53 \) for GCV and
\( 23.50 \) for AIC\( _c \), and the average squared
error rises from the values just quoted to
\( 0.04209 \),
\( 0.04145 \) and
\( 0.03982 \) — roughly a ninefold loss. The fitted
curve in [Figure 43.5.1](05-choosing-lambda.html#fig-smo-criteria)(b) uses
\( 27 \) degrees of freedom and traces every excursion of
the error process.

There is no repair inside this section: a smooth mean and smooth-looking noise
are not separately identified without an assumption about one of them. The honest
responses are to model the correlation (Part VII), to choose \( \lambda \) by a
criterion that accounts for it — leaving out whole blocks rather than single
points, or REML in the mixed-model form of
[Section 43.6](06-mixed-model-and-bayes.html) — or to say that the wiggles are not
interpretable. Opsomer, Wang and Yang (2001) survey what is known.
:::

::: {when-format="html"}
![**Figure 43.5.1.** (a) The three criteria on one data set, each divided by its own
minimum and plotted against effective degrees of freedom; dotted lines mark the
three choices. All are steep on the left and flat on the right. (b) One data set
with AR(1) errors, \( \rho=0.7 \), and the fit
chosen by GCV.](smoothing_criteria.svg){#fig-smo-criteria width=100%}
:::

::: {when-format="pdf"}
![(a) The three criteria on one data set, each divided by its own
minimum and plotted against effective degrees of freedom; dotted lines mark the
three choices. All are steep on the left and flat on the right. (b) One data set
with AR(1) errors, \( \rho=0.7 \), and the fit
chosen by GCV.](smoothing_criteria.pdf){width=100%}
:::

```{.python .run #cell-choose-lambda-criteria}
import numpy as np


def f_true(x):
    """The mean function of Section 43.1."""
    return 2 * x + np.exp(-25 * (x - 0.35) ** 2) - 0.75 * np.exp(-50 * (x - 0.80) ** 2)


def uniform_knots(a, b, K, d):
    """K interior knots equally spaced in (a, b), continued d further steps at each end."""
    delta = (b - a) / (K + 1)
    return a + delta * np.arange(-d, K + 2 + d)


def bspline_basis(x, kn, d):
    """The B-splines of degree d on the knot vector kn, by the de Boor recurrence."""
    x = np.atleast_1d(np.asarray(x, float))
    B = np.array([(x >= kn[j]) & (x < kn[j + 1]) for j in range(len(kn) - 1)], float).T
    last = np.max(np.nonzero(kn[:-1] < kn[1:])[0])
    B[x == kn[-1], last] = 1.0
    for deg in range(1, d + 1):
        C = np.zeros((len(x), B.shape[1] - 1))
        for j in range(C.shape[1]):
            if kn[j + deg] > kn[j]:
                C[:, j] += (x - kn[j]) / (kn[j + deg] - kn[j]) * B[:, j]
            if kn[j + deg + 1] > kn[j + 1]:
                C[:, j] += (kn[j + deg + 1] - x) / (kn[j + deg + 1] - kn[j + 1]) * B[:, j + 1]
        B = C
    return B


def smoother(B, P, lam):
    """The smoother matrix of the penalized fit with basis B and penalty matrix P."""
    return B @ np.linalg.solve(B.T @ B + lam * P, B.T)


def criteria(B, P, y, lam):
    """Leave-one-out CV by the leverage shortcut, GCV, and the corrected AIC."""
    n = len(y)
    S = smoother(B, P, lam)
    resid = y - S @ y
    df = np.trace(S)
    cv = np.mean((resid / (1 - np.diag(S))) ** 2)
    gcv = np.mean(resid**2) / (1 - df / n) ** 2
    aicc = np.log(np.mean(resid**2)) + 1 + 2 * (df + 1) / (n - df - 2)
    return df, cv, gcv, aicc


n, sigma, K, d, k = 120, 0.25, 25, 3, 2
x = (np.arange(1, n + 1) - 0.5) / n
kn = uniform_knots(0.0, 1.0, K, d)
B = bspline_basis(x, kn, d)
Dk = np.diff(np.eye(B.shape[1]), n=k, axis=0)
P = Dk.T @ Dk
rng = np.random.default_rng(4343)
y = f_true(x) + sigma * rng.normal(size=n)

lams = np.geomspace(1e-9, 1e3, 121)
table = np.array([criteria(B, P, y, lam) for lam in lams])
dfs, cvs, gcvs, aiccs = table.T
pick = {"CV": lams[np.argmin(cvs)], "GCV": lams[np.argmin(gcvs)], "AICc": lams[np.argmin(aiccs)]}
df_pick = {name: criteria(B, P, y, lam)[0] for name, lam in pick.items()}
print("degrees of freedom chosen:", {k_: round(v, 2) for k_, v in df_pick.items()})
```

## After the choice

::: {.remark}
[Inference after selecting the smoothing parameter]

Every interval, test and standard error in this chapter is computed *as if*
\( \lambda \) were fixed in advance. It was not: it was chosen by minimizing a
criterion computed from the same \( \y \). This is @prp-sel-selection-bias in a
continuous guise. The residual sum of squares at the selected \( \lambda \) is
the smallest of many, so \( \hat\sigma^2 \) from @eq-smo-sigma-hat is biased
downwards and a pointwise band does not attain its nominal coverage.

The effect is modest when the criterion has a well-determined minimum, which is
why the practice survives, and large when it does not. The defences are those of
[Chapter 29](../ch29-model-selection/index.html): hold out data, bootstrap the
whole procedure including the selection of \( \lambda \)
([Chapter 23](../ch23-resampling-inference/index.html)), or report the fit at
several \( \lambda \) and let the reader see the dependence. What is not a
defence is quoting a \( p \)-value from a fit whose flexibility was chosen to make
that \( p \)-value small.
:::


## Exercises

### A. Check your understanding

::: {#exr-smo-df-three}
[A1]

For the projection onto a \( p \)-dimensional space, compute \( \tr(\bS) \),
\( \tr(\bS\bS\T) \) and \( 2\tr(\bS)-\tr(\bS\bS\T) \). For a P-spline at
\( \mathrm{df}=8 \) with \( m=24 \), which of the three is largest?
Using @prp-smo-linear-smoother(c), say which one the residual degrees of freedom
\( n-2\tr(\bS)+\tr(\bS\bS\T) \) subtracts.
:::

### B. Practice

::: {#exr-smo-shortcut-check}
[B1]

Using the cell above, verify @eq-smo-cv against \( n \) refits at three values of
\( \lambda \), including one so large that \( \mathrm{df}\approx2 \) and one so
small that \( \mathrm{df}\approx m \). Then do the same for a *local linear*
estimator, using \( \ell_i(x_i) \) from
[Section 43.1](01-kernels.html).
:::

::: {#exr-smo-criteria-agree}
[B2]

Show that \( \log\mathrm{GCV}=\log(n^{-1}\norm{\y-\hat{\mathbf f}}^2)
-2\log(1-\mathrm{df}/n) \) and that
\( -2\log(1-t)=2t+t^2+O(t^3) \). Compare with @eq-smo-aicc expanded for small
\( \mathrm{df}/n \), and explain why GCV and AIC almost always choose the same
\( \lambda \) when \( n \) is large and differ when \( \mathrm{df} \) is an
appreciable fraction of \( n \).
:::

::: {.solution}
The first identity is immediate from @eq-smo-gcv. Expanding,
\( \log\mathrm{GCV}=\log(\mathrm{RSS}/n)+2\,\mathrm{df}/n+(\mathrm{df}/n)^2+\dots \),
while
\( \mathrm{AIC_c}=\log(\mathrm{RSS}/n)+1+2(\mathrm{df}+1)/(n-\mathrm{df}-2)
=\log(\mathrm{RSS}/n)+1+2\,\mathrm{df}/n+2\,\mathrm{df}^2/n^2+\dots \). To first
order in \( \mathrm{df}/n \) the two differ by a constant, so their minimizers
agree; the second-order terms differ by a factor of two, and
\( \mathrm{AIC_c} \) has a pole at \( \mathrm{df}=n-2 \) where GCV has one only at
\( \mathrm{df}=n \). Hence they separate exactly when \( \mathrm{df}/n \) is not
small.
:::

### C. Going deeper

::: {#exr-smo-gcv-invariance}
[C1]

GCV was introduced by Craven and Wahba (1979) as a version of cross-validation
invariant to rotation. Show that if \( \boldsymbol{\Gamma} \) is orthogonal and one replaces
\( \y \) by \( \boldsymbol{\Gamma}\y \) and \( \bS_\lambda \) by
\( \boldsymbol{\Gamma}\bS_\lambda\boldsymbol{\Gamma}\T \), then GCV is unchanged while CV is not. Why is
this a virtue for ridge regression (where such rotations are natural) and less
obviously one for a smoother (where they are not)?
:::

::: {.solution}
Under the substitution, the residual vector becomes
\( \boldsymbol{\Gamma}(\I-\bS_\lambda)\y \), whose squared length is unchanged, and
\( \tr(\boldsymbol{\Gamma}\bS_\lambda\boldsymbol{\Gamma}\T)=\tr(\bS_\lambda) \); so @eq-smo-gcv is invariant,
while CV depends on the individual \( s_{ii} \), which are not. In ridge
regression the canonical coordinates of @thm-shr-ridge are as natural as the
original ones, so invariance is attractive; for a smoother the coordinates are
cases in a fixed order along the \( x \)-axis and the \( s_{ii} \) carry real
information — small in dense regions, large in sparse ones — which replacing them
by their average throws away. Hence CV and GCV differ most in designs with
isolated points.
:::

::: {#exr-smo-flat-criterion}
[C2]

For the simulation of @exm-smo-criteria, record for each replicate the *set* of
\( \lambda \) within one per cent of the GCV minimum, and report the median ratio
of the largest to the smallest df in that set. Then plot, for one replicate, the
fits at both ends of the interval. Does the choice within the flat region matter
for what one would say about \( f \)?
:::
