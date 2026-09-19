# Transforming regressors

Transforming a regressor changes only the mean function. The errors are untouched, and there is no
Jacobian, because the data being explained stay the same. What is at stake is bias. If
\( \E(Y)=\X_2\bbeta_2+\beta_1\,g(\x_1) \) and the analyst fits \( \x_1 \) itself, the missing part
\( \beta_1\{g(\x_1)-\x_1\} \) acts as an omitted regressor (@prp-lm-misspecified(a) and @thm-dep-omitted). It shows
up as curvature in the residuals, and tests against a larger model detect it (@prp-cor-augmented-test and @thm-cor-fitted-regressors). This section finds the right \( g \), first by looking, then by estimating a power.

## Component-plus-residual plots

The component-plus-residual plot of [Section 20.1](../ch20-residuals-leverage-influence/01-kinds-of-residuals.html)
plots the partial residuals \( \he+\hat{\beta}_1\x_1 \) against \( \x_1 \). Its least squares line has slope
\( \hat{\beta}_1 \) and residuals \( \he \), as in the added-variable plot, and that section warns that it can mislead
when \( \x_1 \) is related nonlinearly to the other regressors. The next result measures how far.

::: {#prp-tr-cpr}
[What the plot shows]

Let \( \X=[\x_1,\X_2] \) have full column rank, with \( \bone\in\C(\X_2) \), and let \( \hbeta \) and \( \he \) be the
coefficients and residuals of the least squares fit of \( \y \) on \( \X \). Let \( \M \) and \( \M_2 \) project onto
\( \C(\X) \) and \( \C(\X_2) \), and \( \tilde{\x}_1=(\I-\M_2)\x_1 \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \E(\Y)=\mathbf{g}+\X_2\bbeta_2 \), where \( \mathbf{g} \) has entries \( g(x_{i1}) \) for some function \( g \), then
   \[
\E(\he+\hat{\beta}_1\x_1)=\mathbf{g}-\M_2(\mathbf{g}-b\,\x_1),\qquad b=\E(\hat{\beta}_1)=\frac{\tilde{\x}_1\T\mathbf{g}}{\norm{\tilde{\x}_1}^2}.
\]{#eq-tr-cpr-mean}

2. In particular, if \( \C(\X_2) \) is spanned by \( \bone \) alone, the expected partial residuals are
   \( \mathbf{g} \) shifted by a constant, so the plot shows \( g \) exactly apart from its level.
:::

:::

::: {.proof}
(a) Since \( \X_2\bbeta_2\in\C(\X) \), \( \E(\he)=(\I-\M)(\mathbf{g}+\X_2\bbeta_2)=(\I-\M)\mathbf{g} \). By @thm-proj-fwl(a),
\( \hat{\beta}_1=\tilde{\x}_1\T\Y/\norm{\tilde{\x}_1}^2 \), whose mean is \( \tilde{\x}_1\T\mathbf{g}/\norm{\tilde{\x}_1}^2=b \)
because \( \tilde{\x}_1\perp\C(\X_2) \). By @lem-proj-fwl-split(c),
\( \M=\M_2+\tilde{\x}_1\tilde{\x}_1\T/\norm{\tilde{\x}_1}^2 \), so \( (\I-\M)\mathbf{g}=\mathbf{g}-\M_2\mathbf{g}-b\tilde{\x}_1 \). Adding
\( b\x_1=b\tilde{\x}_1+b\M_2\x_1 \) gives \( \mathbf{g}-\M_2\mathbf{g}+b\M_2\x_1 \), which is @eq-tr-cpr-mean.
(b) If \( \C(\X_2)=\spn\{\bone\} \), then \( \M_2\mathbf{v}=\bar{v}\bone \) for every \( \mathbf{v} \).
:::

The distortion \( \M_2(\mathbf{g}-b\x_1) \) is the part of the *nonlinearity* \( \mathbf{g}-b\x_1 \) that the other
regressors reproduce linearly. Correlation with \( \x_1 \) is not in itself harmful; a regressor that resembles
the curvature is. Cook (1993) showed for random regressors that the plot recovers \( g \), up to a linear term,
when the conditional mean of the other regressors given \( x_1 \) is linear in \( x_1 \), and proposed CERES plots
for other cases. The **augmented partial residual plot** of Mallows (1986) adds \( \x_1^2 \) to the fit and plots
\( \he+\hat{\beta}_1\x_1+\hat{\beta}_{11}\x_1^2 \). Its expectation is \( \mathbf{g} \) minus the part of
\( \mathbf{g}-b_1\x_1-b_{11}\x_1^2 \) that the other regressors reproduce, small when \( g \) is nearly quadratic (@exr-tr-augmented).

::: {#exm-tr-cpr-simulated}
[When the plot misleads]

The script simulates \( 120 \) cases with \( \E(Y)=3\log x_1+x_2 \), \( x_1 \) uniform on \( [1,10] \), and measures how
much of the curvature of \( g \) the expected plot keeps (the coefficient of the nonlinear part of \( \mathbf{g} \) in
that of @eq-tr-cpr-mean; \( 1 \) is faithful). In panel (a) of [Figure 22.3.1](#fig-tr-cpr), \( x_2 \) is linear in
\( x_1 \) plus noise, with correlation \( 0.95 \), and the plot keeps \( 0.99 \). In panel (b), \( x_2 \) is a noisy
square root of \( x_1 \), with correlation \( 0.99 \); it resembles the logarithm and takes over much of the
curvature, and the plot keeps only \( 0.45 \), its slope even changing sign. The augmented plot, panel (c),
restores \( 0.88 \).
:::

::: {when-format="html"}
![**Figure 22.3.1.** Component-plus-residual plots for simulated data with \( \E(Y)=3\log x_1+x_2 \). Points:
partial residuals. Dashed: the true \( g \), centred. Small squares: the expected partial residuals @eq-tr-cpr-mean, which scatter because they depend on the simulated \( x_2 \). (a) \( x_2 \) linearly related to \( x_1 \). (b) \( x_2 \) close to \( \sqrt{x_1} \): the plot hides
the curvature. (c) The augmented plot for the data of (b).](cpr_plots.svg){#fig-tr-cpr width=100%}
:::

::: {when-format="pdf"}
![Component-plus-residual plots for simulated data with \( \E(Y)=3\log x_1+x_2 \). Points:
partial residuals. Dashed: the true \( g \), centred. Small squares: the expected partial residuals @eq-tr-cpr-mean, which scatter because they depend on the simulated \( x_2 \). (a) \( x_2 \) linearly related to \( x_1 \). (b) \( x_2 \) close to \( \sqrt{x_1} \): the plot hides
the curvature. (c) The augmented plot for the data of (b).](cpr_plots.pdf){width=100%}
:::

```{.python .run #cell-cpr-plots-cpr}
import numpy as np

rng = np.random.default_rng(2204)
n = 120
x1 = rng.uniform(1, 10, n)
g = 3 * np.log(x1)                                        # the true shape in x1
standardize = lambda v: (v - v.mean()) / v.std()
x2_lin = standardize(x1) + rng.normal(0, 0.3, n)          # linear in x1, plus noise
x2_sqrt = standardize(np.sqrt(x1)) + rng.normal(0, 0.1, n)  # nonlinear in x1, plus noise


def lstsq(X, v):
    coef, *_ = np.linalg.lstsq(X, v, rcond=None)
    return coef


def partial_residuals(y, x2, augmented=False):
    """e + b1 x1 (+ b11 x1^2 if augmented), from the fit of y on 1, x1 (, x1^2), x2."""
    X = np.column_stack([np.ones(len(y)), x1] + ([x1**2] if augmented else []) + [x2])
    coef = lstsq(X, y)
    e = y - X @ coef
    return e + coef[1] * x1 + (coef[2] * x1**2 if augmented else 0)


def curvature_retained(v):
    """Coefficient of the nonlinear part of g in the nonlinear part of v (1 = shape kept)."""
    Z = np.column_stack([np.ones(n), x1])
    nonlin = lambda u: u - Z @ lstsq(Z, u)
    return nonlin(v) @ nonlin(g) / (nonlin(g) @ nonlin(g))


noise = rng.normal(0, 0.6, n)
results = {}
for name, x2, aug in [("linear", x2_lin, False), ("square root", x2_sqrt, False),
                      ("square root, augmented", x2_sqrt, True)]:
    mean = g + x2                                         # E(Y); partial residuals are linear in y
    results[name] = (partial_residuals(mean + noise, x2, aug), partial_residuals(mean, x2, aug))
    print(f"{name:23s} corr(x1, x2) {np.corrcoef(x1, x2)[0, 1]:.2f}   curvature retained "
          f"{curvature_retained(results[name][1]):.2f}")
```

## Power transformations: Box and Tidwell

When the plot suggests a power, it can be estimated. Box and Tidwell (1962) considered the mean function
\[
\E(\Y)=\X_2\bbeta_2+\beta_1\,\x^{(\alpha)},\qquad \x^{(\alpha)}=(x_1^{(\alpha)},\dots,x_n^{(\alpha)})\T,
\]{#eq-tr-box-tidwell}

for a positive regressor \( \x \), with \( \bone\in\C(\X_2) \) and @eq-tr-box-cox applied to the regressor. For fixed
\( \alpha \) this is a linear model, and \( \hat{\alpha} \) minimizes the residual sum of squares \( \text{SSE}(\alpha) \) of the fit of
\( \y \) on \( [\X_2,\x^{(\alpha)}] \). Under normal errors it is the maximum likelihood estimate, and
\( \{\alpha:n\log\{\text{SSE}(\alpha)/\text{SSE}(\hat{\alpha})\}\le\chi^2_{1,0.05}\} \) is an asymptotic \( 95\% \) interval under the
conditions of Wilks's theorem, which include \( \beta_1\ne0 \) (@exr-tr-no-slope). The derivative of \( x^{(\alpha)} \) with respect to \( \alpha \) at
\( \alpha_0\ne0 \) is \( x^{\alpha_0}\log x/\alpha_0-x^{(\alpha_0)}/\alpha_0 \), and at \( \alpha_0=0 \) it is \( \tfrac12(\log x)^2 \) (@exr-tr-tidwell-derivative). Since \( x^{(\alpha_0)} \) and \( \bone \) are already in the model, the part that matters
is the **constructed variable**
\[
\bw=\begin{cases}\bigl(x_i^{\alpha_0}\log x_i/\alpha_0\bigr)_{i=1}^n,&\alpha_0\ne0,\\[0.6ex]
\bigl(\tfrac12(\log x_i)^2\bigr)_{i=1}^n,&\alpha_0=0.\end{cases}
\]{#eq-tr-power-constructed}

For \( \alpha_0=1 \) it is \( x\log x \).

::: {#thm-tr-box-tidwell}
[Box–Tidwell]

Let \( \X_0=[\X_2,\x^{(\alpha_0)}] \) have full column rank \( p \), with \( \bone\in\C(\X_2) \), let \( \bw \) be the
constructed variable @eq-tr-power-constructed, and suppose \( [\X_0,\bw] \) has full column rank \( p+1<n \). Let
\( \hat{\beta}_1 \) be the coefficient of \( \x^{(\alpha_0)} \) and \( \he_0 \) the residual vector in the least squares fit of
\( \y \) on \( \X_0 \), and let \( \hat{\gamma} \) be the coefficient of \( \bw \) in the fit of \( \y \) on \( [\X_0,\bw] \).

::: {.enumerate options="label=(\alph*)"}
1. **Exact test.** If \( \Y\sim\Normal_n(\X_0\bbeta,\sigma^2\I) \), which is @eq-tr-box-tidwell with \( \alpha=\alpha_0 \), then the
   \( t \) statistic of \( \hat{\gamma} \) has the \( t(n-p-1) \) distribution.

2. **The step.** If \( \hat{\beta}_1\ne0 \), one Gauss–Newton step for minimizing \( \norm{\y-\X_2\bbeta_2-\beta_1\x^{(\alpha)}}^2 \),
   started from the least squares coefficients at \( \alpha_0 \), moves \( \alpha \) to
   \[
\alpha_1=\alpha_0+\hat{\gamma}/\hat{\beta}_1 .
\]{#eq-tr-power-step}

3. **Stationarity.** \( \text{SSE}(\alpha) \) is differentiable at \( \alpha_0 \) with
   \( \text{SSE}'(\alpha_0)=-2\hat{\beta}_1\he_0\T\bw=-2\hat{\beta}_1\hat{\gamma}\,\norm{(\I-\M_0)\bw}^2 \), where \( \M_0 \) projects onto
   \( \C(\X_0) \). So if \( \hat{\beta}_1\ne0 \), then \( \hat{\gamma}=0 \) iff \( \alpha_0 \) is a stationary point of \( \text{SSE} \), and every
   fixed point of the iteration @eq-tr-power-step is a stationary point.
:::

:::

::: {.proof}
Write \( \mathbf{d}=\partial\x^{(\alpha)}/\partial\alpha \) at \( \alpha_0 \). By the derivative computed above,
\( \mathbf{d}-\bw \) is a multiple of \( \x^{(\alpha_0)} \), so \( \mathbf{d}-\bw\in\C(\X_0) \) and \( (\I-\M_0)\mathbf{d}=(\I-\M_0)\bw=:\tilde{\bw} \),
which is nonzero because \( [\X_0,\bw] \) has full rank. By @thm-proj-fwl(a),
\( \hat{\gamma}=\tilde{\bw}\T\y/\norm{\tilde{\bw}}^2 \), and \( \tilde{\bw}\T\y=\bw\T(\I-\M_0)\y=\bw\T\he_0 \).

(a) The column \( \bw \) depends on \( \x \) alone, not on \( \Y \). The \( t \) test for \( \hat{\gamma} \) is therefore the \( F \)
test of \( \C(\X_0) \) against the fixed larger space \( \C([\X_0,\bw]) \), with \( t^2=F \) and a symmetric null
distribution. By @prp-cor-augmented-test(a), \( F\sim F(1,n-p-1) \) when \( \E(\Y)\in\C(\X_0) \), so \( t\sim t(n-p-1) \).

(b) The mean function \( \bmu(\bbeta_2,\beta_1,\alpha)=\X_2\bbeta_2+\beta_1\x^{(\alpha)} \) has derivative matrix
\( \mathbf{J}=[\X_2,\x^{(\alpha_0)},\hat{\beta}_1\mathbf{d}] \) at the starting point, where the residual is \( \he_0 \). The
Gauss–Newton step is the least squares coefficient vector of \( \he_0 \) on \( \mathbf{J} \), and its last entry is the
change in \( \alpha \). By @thm-proj-fwl(a) that entry is
\[
\frac{\{(\I-\M_0)\hat{\beta}_1\mathbf{d}\}\T\he_0}{\norm{(\I-\M_0)\hat{\beta}_1\mathbf{d}}^2}
=\frac{1}{\hat{\beta}_1}\,\frac{\tilde{\bw}\T\he_0}{\norm{\tilde{\bw}}^2}
=\frac{\hat{\gamma}}{\hat{\beta}_1},
\]
using \( \tilde{\bw}\T\he_0=\tilde{\bw}\T\y \), which holds because \( \tilde{\bw}\perp\C(\X_0) \).

(c) Let \( \X(\alpha)=[\X_2,\x^{(\alpha)}] \), which has full column rank near \( \alpha_0 \), and let \( \hbeta(\alpha) \) be its least
squares coefficients, a differentiable function of \( \alpha \). Then \( \text{SSE}(\alpha)=\norm{\y-\X(\alpha)\hbeta(\alpha)}^2 \), and
\[
\text{SSE}'(\alpha_0)=-2\he_0\T\bigl\{\hat{\beta}_1\mathbf{d}+\X_0\hbeta'(\alpha_0)\bigr\}=-2\hat{\beta}_1\he_0\T\mathbf{d}
=-2\hat{\beta}_1\he_0\T\bw ,
\]
because \( \he_0\perp\C(\X_0) \) and \( \mathbf{d}-\bw\in\C(\X_0) \). Finally \( \he_0\T\bw=\hat{\gamma}\norm{\tilde{\bw}}^2 \) from the
first paragraph.
:::

Part (a) is exact because this constructed variable is built from the fixed regressor, not from the
response as in @prp-tr-constructed. Part (b) identifies the Box–Tidwell iteration with Gauss–Newton, fast
near a well-determined minimum and liable to oscillate otherwise, so in one dimension it is safer to
minimize \( \text{SSE}(\alpha) \) directly; by part (c) the two agree at convergence. Several powers can be stepped
at once (@exr-tr-several-powers).

::: {#exm-tr-engel-box-tidwell}
[Transforming income]

With \( \log(\text{food}) \) as the response and \( \text{income}^{(\alpha)} \) as the regressor, starting from \( \alpha_0=1 \), the Box–Tidwell iteration gives \( -0.053 \) after one step and
\( -0.182 \) after two, where it stays. The direct minimization agrees: \( \hat{\alpha}=-0.182 \), with \( 95\% \)
likelihood ratio interval \( [-0.334,\ -0.036] \). The constructed-variable \( t \) statistic is \( -16.58 \) at
\( \alpha_0=1 \), so untransformed income is decisively rejected, and \( -2.47 \) at \( \alpha_0=0 \), beyond the
critical value \( 1.97 \), so the logarithm is rejected at the \( 5\% \) level as well. As a check on
part (a), \( 4000 \) data sets simulated from the log–log fit give a rejection rate of \( 0.047 \) for this
test.

The rejection rests largely on one household. The richest, with income \( 4958 \) and food expenditure
\( 1827 \), sits far to the right in panel (c) of
[Figure 22.1.1](01-transforming-the-response.html#fig-tr-engel-scales), below the line. Without it,
\( \hat{\alpha}=-0.119 \), and the likelihood ratio statistic for \( \alpha=0 \) falls from \( 5.98 \) to \( 1.77 \). A negative
power says food expenditure flattens among the richest households, where the data are few; Cook's distance (@def-res-cooks)
pursues this.
:::

```{.python .run #cell-box-tidwell-tidwell}
import numpy as np
import statsmodels.api as sm

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
y = np.log(food)


def bc(x, a):
    """Box-Cox transform (x^a - 1)/a, with log x at a = 0."""
    return np.log(x) if abs(a) < 1e-12 else np.expm1(a * np.log(x)) / a


def fit(y, x, a):
    return sm.OLS(y, np.column_stack([np.ones(len(y)), bc(x, a)])).fit()


def box_tidwell_step(y, x, a0):
    """One Box-Tidwell (Gauss-Newton) step: a1 = a0 + gamma/beta1."""
    beta1 = fit(y, x, a0).params[1]                    # slope of x^(a0) in the null fit
    w = x**a0 * np.log(x) / a0 if a0 != 0 else np.log(x) ** 2 / 2   # d bc(x, a)/da, up to C(X0)
    aug = sm.OLS(y, np.column_stack([np.ones(len(y)), bc(x, a0), w])).fit()
    return a0 + aug.params[2] / beta1, aug.tvalues[2]


a, path = 1.0, [1.0]
for _ in range(8):
    a, t = box_tidwell_step(y, income, a)
    path.append(a)
print("Box-Tidwell iterates from alpha = 1:", " ".join(f"{v:.4f}" for v in path))
```

## Choosing both scales together

The best power for the response depends on the scale of the regressor, and conversely, so the two can be
estimated together. With response \( \mathbf{z}^{(\lambda)} \) and model matrix \( [\bone,\x^{(\alpha)}] \) the proof of
@thm-tr-box-cox goes through unchanged, since the Jacobian involves only \( \y \): the joint profile
log-likelihood is \( -\frac n2\log\text{SSE}_z(\lambda,\alpha) \) plus a constant.

::: {#exm-tr-engel-joint}
[Two powers for Engel's data]

[Figure 22.3.2](#fig-tr-joint) shows the joint profile. The maximum is at \( (\hat{\lambda},\hat{\alpha})=(-0.25,\ -0.39) \). The raw scale
\( (1,1) \) has likelihood ratio statistic \( 186.5 \), and the best pair with untransformed income,
\( (0.63,1) \) from @exm-tr-engel-box-cox, has \( 167.8 \): whatever is done to the response, income needs
transforming. The log–log pair \( (0,0) \) has statistic \( 11.2 \), outside the \( 95\% \) region, whose cut-off is
\( 5.99 \) for two parameters. Without the richest household the statistic falls to \( 6.55 \), still outside.

The region is a long tilted ellipse, since a lower response power can be partly compensated by a lower
income power (@exr-tr-power-law-ridge). The log–log model lies outside the region (\( 11.2 \) against \( 5.99 \)),
largely because of the richest household. It is still the natural working choice for interpretation, since its
slope is an elasticity, with the caveat that the data reject a constant elasticity and that it may fall among
the richest households (@prp-tr-elasticity).
:::

::: {when-format="html"}
![**Figure 22.3.2.** Joint profile log-likelihood of the response power \( \lambda \) and the income power
\( \alpha \) for Engel's data. The thick contour bounds the \( 95\% \) likelihood ratio region, the thin contours are
\( 10 \), \( 20 \) and \( 40 \) units of log-likelihood below the maximum (the dot). Squares mark the log–log and
raw models.](joint_profile.svg){#fig-tr-joint width=62%}
:::

::: {when-format="pdf"}
![Joint profile log-likelihood of the response power \( \lambda \) and the income power
\( \alpha \) for Engel's data. The thick contour bounds the \( 95\% \) likelihood ratio region, the thin contours are
\( 10 \), \( 20 \) and \( 40 \) units of log-likelihood below the maximum (the dot). Squares mark the log–log and
raw models.](joint_profile.pdf){width=62%}
:::

## Polynomials and other alternatives

Powers suit positive regressors with monotone, smoothly bending effects. Otherwise the alternative is to add
columns rather than change one: a polynomial or a spline. These buy flexibility with interpretability, since
no single coefficient measures the effect of \( x \) (@def-lm-coefficient-interpretation, third qualification),
and high-degree polynomials behave badly at the edges of the data. Royston and Altman's (1994) *fractional
polynomials* use one or two terms from a small set of powers including the logarithm. Chapter 42 treats
polynomial and spline regression.

## Exercises

### A. Check your understanding

::: {#exr-tr-tidwell-hand}
[A1]

At \( \alpha_0=1 \) a fit gives \( \hat{\beta}_1=2.5 \), and adding \( x\log x \) gives it the coefficient \( \hat{\gamma}=-1.5 \). What is
the next Box–Tidwell value of \( \alpha \)? What does the sign of \( \hat{\gamma} \) say about the curvature of the relation?
:::

### B. Practice

::: {#exr-tr-tidwell-derivative}
[B1]

Show that \( \partial x^{(\alpha)}/\partial\alpha=x^{\alpha}\log x/\alpha-x^{(\alpha)}/\alpha \) for \( \alpha\ne0 \), and that its limit as \( \alpha\to0 \)
is \( \tfrac12(\log x)^2 \). Conclude that in a model containing \( \bone \) and \( \x^{(\alpha_0)} \) the constructed variable
may be taken to be @eq-tr-power-constructed.
:::

::: {.solution}
Differentiate \( (x^{\alpha}-1)/\alpha \): the quotient rule gives
\( \{\alpha x^{\alpha}\log x-(x^{\alpha}-1)\}/\alpha^2=x^{\alpha}\log x/\alpha-x^{(\alpha)}/\alpha \). As \( \alpha\to0 \), expand
\( x^{\alpha}=1+\alpha\log x+\tfrac12\alpha^2(\log x)^2+O(\alpha^3) \): the numerator is
\( \alpha\log x+\alpha^2(\log x)^2-\alpha\log x-\tfrac12\alpha^2(\log x)^2+O(\alpha^3)=\tfrac12\alpha^2(\log x)^2+O(\alpha^3) \), and dividing by
\( \alpha^2 \) gives the limit. The term \( -\x^{(\alpha_0)}/\alpha_0 \) lies in the model space, and removing it changes neither the
residualized column \( \tilde{\bw} \) nor \( \hat{\gamma} \).
:::

::: {#exr-tr-power-law-ridge}
[B2]

Suppose \( y=c\,x^{\beta} \) exactly, with \( c,\beta>0 \). Show that for \( \lambda\ne0 \), \( y^{(\lambda)} \) is an affine function of
\( x^{(\alpha)} \) iff \( \alpha=\beta\lambda \), and that for \( \lambda=0 \) it is affine in \( x^{(\alpha)} \) iff \( \alpha=0 \). Explain why
data close to a power law give a joint likelihood for \( (\lambda,\alpha) \) with a ridge, and why the slope of the
ridge need not be exactly \( \beta \) when the variance also depends on \( \lambda \).
:::

::: {.solution}
For \( \lambda\ne0 \), \( y^{\lambda}=c^{\lambda}x^{\beta\lambda} \), so \( y^{(\lambda)}=(c^{\lambda}x^{\beta\lambda}-1)/\lambda \), an affine function of
\( x^{\beta\lambda} \). An affine function of \( x^{\beta\lambda} \) is affine in \( x^{\alpha} \) (equivalently in \( x^{(\alpha)} \)) for all
\( x>0 \) iff \( \alpha=\beta\lambda \), since distinct powers, and a power and the logarithm, are affinely independent
functions on \( (0,\infty) \). For \( \lambda=0 \), \( \log y=\log c+\beta\log x \) is affine in \( x^{(\alpha)} \) iff \( \alpha=0 \), which is again
\( \alpha=\beta\lambda \). So every pair on the line \( \alpha=\beta\lambda \) makes the mean exactly linear, and nearby data make
the residual sum of squares small along it. The likelihood also rewards a constant variance and
normal errors, which depend on \( \lambda \) alone, so the maximum lies somewhere along the tilted valley,
and its exact direction mixes the two requirements.
:::

::: {#exr-tr-augmented}
[B3]

Let \( \X=[\x_1,\x_1^2,\X_2] \), with \( \bone\in\C(\X_2) \) and full column rank, and let \( \E(\Y)=\mathbf{g}+\X_2\bbeta_2 \).
Show that the augmented partial residuals \( \he+\hat{\beta}_1\x_1+\hat{\beta}_{11}\x_1^2 \) have expectation
\( \mathbf{g}-\M_2(\mathbf{g}-b_1\x_1-b_{11}\x_1^2) \), where \( (b_1,b_{11}) \) are the expectations of \( (\hat{\beta}_1,\hat{\beta}_{11}) \).
Conclude that the plot is exact, up to a vertical shift, when \( g \) is a quadratic.
:::

::: {.solution}
Let \( \mathbf{Q}=[\x_1,\x_1^2] \) and \( \tilde{\mathbf{Q}}=(\I-\M_2)\mathbf{Q} \). As in the proof of @prp-tr-cpr,
\( \M=\M_2+\tilde{\mathbf{Q}}(\tilde{\mathbf{Q}}\T\tilde{\mathbf{Q}})^{-1}\tilde{\mathbf{Q}}\T \) and the coefficient vector \( \mathbf{b} \) of \( \mathbf{Q} \) has
expectation \( (\tilde{\mathbf{Q}}\T\tilde{\mathbf{Q}})^{-1}\tilde{\mathbf{Q}}\T\mathbf{g} \). Then
\( \E(\he)=(\I-\M)\mathbf{g}=\mathbf{g}-\M_2\mathbf{g}-\tilde{\mathbf{Q}}\mathbf{b} \), and adding
\( \mathbf{Q}\mathbf{b}=\tilde{\mathbf{Q}}\mathbf{b}+\M_2\mathbf{Q}\mathbf{b} \) gives \( \mathbf{g}-\M_2(\mathbf{g}-\mathbf{Q}\mathbf{b}) \). If
\( \mathbf{g}=c_0\bone+c_1\x_1+c_2\x_1^2 \), then \( \mathbf{g} \) lies in \( \C(\X) \), \( \mathbf{b}=(c_1,c_2)\T \), and
\( \mathbf{g}-\mathbf{Q}\mathbf{b}=c_0\bone \), whose projection is \( c_0\bone \).
:::

::: {#exr-tr-cpr-quadratic}
[B4]

In the setting of @prp-tr-cpr, let \( \X_2=[\bone,\mathbf{q}] \) with \( q_i=x_{i1}^2 \), so the second regressor is an exact
quadratic in the first, and write \( \bbeta_2=(\beta_{20},\delta)\T \).

::: {.enumerate options="label=(\alph*)"}
1. If \( g(x)=c\,x^2 \), show that \( b=0 \) and that the expected partial residuals are \( \mathbf{0} \): the plot is flat,
   although \( \E(Y) \) curves in \( x_1 \) when \( \delta\ne-c \).

2. If \( g(x)=c\,x \), show that the plot is exact.

3. Explain why no plot against \( x_1 \) can tell which of the two regressors carries the curvature in (a).
:::

:::

::: {.solution}
(a) Here \( \mathbf{g}=c\,\mathbf{q}\in\C(\X_2) \), so \( \tilde{\x}_1\T\mathbf{g}=0 \) and \( b=0 \); then
\( \M_2(\mathbf{g}-b\x_1)=\M_2\mathbf{g}=\mathbf{g} \) and @eq-tr-cpr-mean is \( \mathbf{0} \). The whole of \( g \) is the nonlinearity
\( \mathbf{g}-b\x_1 \), and \( \mathbf{q} \) reproduces all of it. The mean \( \E(Y)=\beta_{20}+(c+\delta)x_1^2 \) curves unless
\( c+\delta=0 \).
(b) Now \( b=c\,\tilde{\x}_1\T\x_1/\norm{\tilde{\x}_1}^2=c \), because \( \tilde{\x}_1\T\x_1=\norm{\tilde{\x}_1}^2 \), so
\( \mathbf{g}-b\x_1=\mathbf{0} \) and the expected partial residuals are \( \mathbf{g} \).
(c) The models \( \E(Y)=\beta_{20}+c\,x_1^2+\delta x_1^2 \) with the same \( c+\delta \) have the same mean for every case:
the split between \( g \) and \( \X_2\bbeta_2 \) is not identified from the data, so no display of them can recover it.
:::

### C. Going deeper

::: {#exr-tr-no-slope}
[C1]

Suppose \( \beta_1=0 \) in @eq-tr-box-tidwell. Show that the model does not depend on \( \alpha \), so \( \alpha \) is not
identified. Show that the test of @thm-tr-box-tidwell(a) still has exact size, and argue that its power
against changes in \( \alpha \) tends to zero as \( \beta_1\to0 \). What happens to the Box–Tidwell step, and to the likelihood
ratio interval for \( \alpha \)?
:::

::: {#exr-tr-several-powers}
[C2]

Let \( \E(\Y)=\X_3\bbeta_3+\beta_1\x_1^{(\alpha_1)}+\beta_2\x_2^{(\alpha_2)} \). Show that one Gauss–Newton step from
\( (\alpha_1^0,\alpha_2^0) \) and the least squares coefficients there changes \( \alpha_j \) by \( \hat{\gamma}_j/\hat{\beta}_j \), where
\( \hat{\gamma}_1,\hat{\gamma}_2 \) are the coefficients of the two constructed variables in a single regression of \( \y \) on
the model matrix augmented by both.
:::
