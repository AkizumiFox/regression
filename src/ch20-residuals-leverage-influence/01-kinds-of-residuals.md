# Kinds of residuals

The model makes its assumptions about the errors, which are never observed. The residuals are
their visible stand-ins, imperfect in a precise way: they are the projection of the errors onto a
subspace of dimension \( n-p \), so they are correlated, have unequal variances, and are shrunk most
at the cases that matter most. This section defines rescaled residuals, finds their
laws, and asks what residual plots reveal.

## The setting

Throughout the chapter,
\[
\Y=\X\bbeta+\be,\qquad \E(\be)=\bzero,\qquad \Cov(\be)=\sigma^2\I,
\]
with \( \X \) of size \( n\times p \) and full column rank \( p<n \). Results that need normal errors say so
and then assume @eq-opt-normal-model. The hat matrix is \( \M=\X(\X\T\X)^{-1}\X\T \), and \( h_{ij} \) denotes its
\( (i,j) \) entry, so that \( h_{ii} \) is the leverage of case \( i \) (@def-proj-leverage). The \( i \)th row of
\( \X \), written as a column, is \( \x_{(i)} \), and \( \vect{e}_i \) is the \( i \)th coordinate vector of \( \Real^n \).
The residual vector is \( \he=(\I-\M)\Y \), with entries \( \hat{\varepsilon}_i \), and
\( s^2=\norm{\he}^2/(n-p) \).

Everything in the chapter depends on \( \X \) only through \( \C(\X) \). If \( \X \) has rank \( r<p \), replace \( p \) by
\( r \) in every degree of freedom and every formula that involves \( \M \) alone. Formulas that involve
\( (\X\T\X)^{-1} \) hold for any full-rank reparameterization (@thm-proj-reparam).

Two facts from earlier chapters drive everything that follows. First, by @prp-proj-fit-residual(e),
\[
\he=(\I-\M)\be ,
\]{#eq-res-residual-error}

so the residuals are a fixed linear image of the errors. Second, by
@prp-lm-fit-moments(b), \( \Cov(\he)=\sigma^2(\I-\M) \). So \( \Var(\hat{\varepsilon}_i)=\sigma^2(1-h_{ii}) \) and, for \( i\ne j \),
\[
\operatorname{Corr}(\hat{\varepsilon}_i,\hat{\varepsilon}_j)=\frac{-h_{ij}}{\sqrt{(1-h_{ii})(1-h_{jj})}}.
\]{#eq-res-residual-moments}

A case with leverage near one has a residual with variance near zero: the fit passes close to it
whatever its response. Raw residuals are therefore not comparable across cases.

## Four residuals

::: {#def-res-residuals}
[Standardized and studentized residuals]

Let \( h_{ii}<1 \). For case \( i \):

::: {.enumerate options="label=(\alph*)"}
1. the **raw residual** is \( \hat{\varepsilon}_i=y_i-\x_{(i)}\T\hbeta \);

2. the **standardized residual** is \( z_i=\hat{\varepsilon}_i/\bigl(\sigma\sqrt{1-h_{ii}}\bigr) \), which needs the true \( \sigma \);

3. the **internally studentized residual** is
   \[
   r_i=\frac{\hat{\varepsilon}_i}{s\sqrt{1-h_{ii}}} ;
   \]

4. the **externally studentized residual** is
   \[
   t_i=\frac{\hat{\varepsilon}_i}{s_{(i)}\sqrt{1-h_{ii}}},
   \]
   where \( s_{(i)}^2=\text{SSE}_{(i)}/(n-p-1) \) is the residual mean square of the least squares fit
   to the \( n-1 \) cases other than \( i \).
:::

:::

Under normal errors \( z_i \) is exactly \( \Normal(0,1) \), but it cannot be computed. The studentized
residuals replace \( \sigma \) by an estimate: \( s \), which includes case \( i \), or \( s_{(i)} \), which does not, so that a
gross error in \( y_i \) cannot inflate the yardstick it is measured against.
[Section 20.3](03-deletion.html) shows that \( s_{(i)} \) costs nothing extra and that \( t_i \) has an exact \( t \) law.

::: {.remark}
[Names]

Terminology varies. Some authors call \( r_i \), or even \( \hat{\varepsilon}_i/s \), the standardized residual, and
\( t_i \) is also called the studentized deleted residual or RSTUDENT. R's `rstandard` and `rstudent` return
\( r_i \) and \( t_i \), as do `resid_studentized_internal` and `resid_studentized_external` in statsmodels.
:::

## The law of the internally studentized residual

The residual \( r_i \) is not \( t \)-distributed, because \( \hat{\varepsilon}_i \) is part of \( s \). Its law is a rescaled beta.

::: {#prp-res-internal}
[Law of the internally studentized residual]

Assume @eq-opt-normal-model, \( n-p\ge2 \) and \( h_{ii}<1 \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( r_i^2/(n-p)\sim\mathrm{Beta}\bigl(\tfrac12,\tfrac{n-p-1}{2}\bigr) \);

2. \( \lvert r_i\rvert\le\sqrt{n-p} \), \( \E(r_i)=0 \) and \( \Var(r_i)=1 \);

3. the law of \( r_i \) is the same for every case and does not depend on \( \X \), \( \bbeta \) or \( \sigma^2 \).
:::

:::

::: {.proof}
Put \( \bu=(\I-\M)\vect{e}_i/\sqrt{1-h_{ii}} \). Since \( \norm{(\I-\M)\vect{e}_i}^2=\vect{e}_i\T(\I-\M)\vect{e}_i=1-h_{ii} \),
\( \bu \) is a unit vector in \( \C(\X)\perpc \), and by @eq-res-residual-error
\( \hat{\varepsilon}_i=\vect{e}_i\T(\I-\M)\be=\sqrt{1-h_{ii}}\,\bu\T\be \). Let \( \bP_1=\bu\bu\T \) and
\( \bP_2=\I-\M-\bu\bu\T \). Because \( (\I-\M)\bu=\bu \), both are symmetric and idempotent, \( \bP_1\bP_2=\mathbf{0} \),
and their ranks are \( 1 \) and \( n-p-1 \). They split the residual sum of squares:
\( \text{SSE}=\be\T(\I-\M)\be=\norm{\bP_1\be}^2+\norm{\bP_2\be}^2 \). By @thm-qf-orthogonal-projections, applied to
\( \be\sim\Normal_n(\bzero,\sigma^2\I) \), the variables \( U=(\bu\T\be)^2/\sigma^2\sim\chi^2(1) \) and
\( V=\norm{\bP_2\be}^2/\sigma^2\sim\chi^2(n-p-1) \) are independent. Now
\[
\frac{r_i^2}{n-p}=\frac{\hat{\varepsilon}_i^2}{(1-h_{ii})\,\text{SSE}}=\frac{U}{U+V},
\]
which is \( \mathrm{Beta}\bigl(\tfrac12,\tfrac{n-p-1}2\bigr) \) by @exr-qf-beta. This is (a). A beta variable lies in
\( [0,1] \), which gives the bound in (b). Replacing \( \be \) by \( -\be \) changes the sign of \( r_i \) and not its law,
so \( \E(r_i)=0 \), and \( \E(r_i^2)=(n-p)\cdot\tfrac12\big/\tfrac{n-p}{2}=1 \). Part (c) holds because the law in
(a), together with the symmetry, determines the law of \( r_i \).
:::

The bound in (b) bites when \( n-p \) is small: with \( n-p=5 \), no \( \lvert r_i\rvert \) can exceed
\( \sqrt5\approx2.236 \), however wild the case, because a gross error inflates \( s \) along with \( \hat{\varepsilon}_i \).
The \( r_i \) are not independent either; their correlations are those of @eq-res-residual-moments (@exr-res-internal-correlation).

## What a residual plot can show

A residual plot shows residuals against the fitted values, a regressor, an omitted variable or the
order of collection, to reveal structure the model has missed. Least squares removes some structure
automatically, and it is worth knowing exactly which.

::: {#prp-res-plots}
[What least squares leaves in the residuals]

Suppose \( \bone\in\C(\X) \).

::: {.enumerate options="label=(\alph*)"}
1. The residual vector is orthogonal to \( \bone \), to every column of \( \X \) and to \( \hY \). Hence the least
   squares line through a plot of the residuals against the fitted values, or against any regressor in the
   model, is the horizontal axis.

2. The least squares line through a plot of the residuals against the responses \( y_i \) has slope \( 1-R^2 \).

3. If in fact \( \E(\Y)=\X\bbeta+\boldsymbol{\updelta} \) for some \( \boldsymbol{\updelta}\in\Real^n \), and \( \Cov(\Y)=\sigma^2\I \), then
   \( \E(\he)=(\I-\M)\boldsymbol{\updelta} \) and \( \Cov(\he)=\sigma^2(\I-\M) \).
:::

:::

::: {.proof}
(a) is @prp-proj-fit-residual(a) and (b). A line fitted with an intercept to points \( (a_i,\hat{\varepsilon}_i) \) has
slope \( \sum_i(a_i-\bar a)\hat{\varepsilon}_i/\sum_i(a_i-\bar a)^2 \), and the numerator is
\( \mathbf{a}\T\he-\bar a\,\bone\T\he=0 \) when \( \mathbf{a} \) is \( \hY \) or a column of \( \X \).
(b) With \( \mathbf{a}=\y \) the numerator is \( \y\T\he=(\hY+\he)\T\he=\norm{\he}^2=\text{SSE} \), and the denominator is
\( \sum_i(y_i-\bar y)^2=\text{SST} \). The slope is \( \text{SSE}/\text{SST}=1-R^2 \).
(c) \( \he=(\I-\M)\Y \), and @thm-rv-linear gives mean \( (\I-\M)(\X\bbeta+\boldsymbol{\updelta})=(\I-\M)\boldsymbol{\updelta} \) and covariance
\( \sigma^2(\I-\M)^2=\sigma^2(\I-\M) \).
:::

By (a), a plot against fitted values or regressors has had its linear trend removed, so any remaining
trend is curved, a sign of a misspecified mean. By (b), a plot against the responses always slopes
upward, even when the model is right, so it is not used.

Part (c) is the sobering one. The residuals are centred at \( (\I-\M)\boldsymbol{\updelta} \), the part of the error in the
mean orthogonal to the model space. The part of \( \boldsymbol{\updelta} \) inside \( \C(\X) \) leaves no trace, yet it is exactly
this part, \( (\X\T\X)^{-1}\X\T\boldsymbol{\updelta} \), that biases the coefficients (@prp-lm-misspecified(a)). An omitted regressor
highly correlated with the included ones biases their coefficients (@thm-dep-omitted) while leaving the residuals almost
clean. Residual plots can show that the *fitted values* are wrong, never that the *coefficients* are right.

The standard plots are these. Against **fitted values**: curvature suggests a misspecified mean, a fan shape variance that
grows with the mean, isolated points outliers. Against **each regressor**: curvature points to that regressor. Against a
variable **not in the model**: better replaced by the added-variable plot of [Section 6.6](../ch06-projections/06-fwl.html),
whose slope is the coefficient the variable would receive (@exr-res-plot-omitted). Against **collection order**: runs of equal
sign suggest serial correlation, tested in [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) (@thm-het-durbin-watson). A **normal quantile plot**: curvature suggests skewed or heavy-tailed errors (@prp-het-normal-check).
Studentized residuals are preferable to raw ones throughout, since raw residuals make high-leverage cases look better than they are.

::: {#exm-res-patterns}
[Four residual plots]

Each panel of [Figure 20.1.1](#fig-res-patterns) comes from \( 80 \) simulated cases with \( x \) uniform on
\( (0,10) \), fitted by a straight line. In (a) the model is right and the residuals form a structureless band. In (b)
the mean has a quadratic term, and the residuals show the parabola the line could not absorb. In (c) the error standard
deviation grows with \( x \), and the band opens into a fan. In (d) one response has an error of \( 8\sigma \); it stands out,
and the line has been pulled slightly towards it. In every panel the residuals are orthogonal to the fitted values, as
@prp-res-plots(a) requires; the listing prints the inner product.
:::

::: {when-format="html"}
![**Figure 20.1.1.** Residuals against fitted values for four straight-line fits: (a) correct model, (b) curved
mean, (c) variance growing with the mean, (d) one gross error.](residual_patterns.svg){#fig-res-patterns width=100%}
:::

::: {when-format="pdf"}
![Residuals against fitted values for four straight-line fits: (a) correct model, (b) curved
mean, (c) variance growing with the mean, (d) one gross error.](residual_patterns.pdf){width=100%}
:::

```{.python .run #cell-residual-plots-patterns}
import numpy as np

def ls_fit(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta, X @ beta, y - X @ beta

rng = np.random.default_rng(7)
n = 80
x = rng.uniform(0, 10, n)
eps = rng.normal(size=n)
X = np.column_stack([np.ones(n), x])
responses = {
    "(a) model correct": 1 + 0.5 * x + eps,
    "(b) curved mean": 1 + 0.5 * x + 0.08 * (x - 5) ** 2 + 0.6 * eps,
    "(c) variance grows": 1 + 0.5 * x + (0.1 + 0.25 * x) * eps,
    "(d) one gross error": 1 + 0.5 * x + eps + 8.0 * (np.arange(n) == 17),
}
for name, y in responses.items():
    beta, fitted, resid = ls_fit(X, y)
    print(f"{name:22s} slope {beta[1]:.3f}   residuals . fitted = {resid @ fitted:.1e}")
```

## Partial residuals

A plot of residuals against a regressor \( \x_j \) in the model has a flat least squares line (@prp-res-plots(a)). The **component-plus-residual plot**, or partial residual plot, adds the fitted
linear component back. It plots
\[
\hat{\varepsilon}_i+\hat{\beta}_jx_{ij}\quad\text{against}\quad x_{ij},\qquad i=1,\dots,n.
\]
Compare it with the added-variable plot of [Section 6.6](../ch06-projections/06-fwl.html), which plots
\( \tilde{\y}=(\I-\M_{(j)})\y \) against \( \tilde{\x}_j=(\I-\M_{(j)})\x_j \), where \( \M_{(j)} \) projects onto the span of the
columns other than \( \x_j \).

::: {#prp-res-partial-residual}
[Two plots with the same line]

Let \( \bone\in\C(\X) \), and let \( \x_j \) be a column of \( \X \) that is not a multiple of \( \bone \).

::: {.enumerate options="label=(\alph*)"}
1. The least squares line of \( \he+\hat{\beta}_j\x_j \) on \( [\bone,\x_j] \) has intercept \( 0 \), slope \( \hat{\beta}_j \), and
   residual vector \( \he \).

2. The least squares line through the origin of \( \tilde{\y} \) on \( \tilde{\x}_j \) has slope \( \hat{\beta}_j \) and residual vector \( \he \).
:::

:::

::: {.proof}
(a) The vector \( \he \) is orthogonal to \( \bone \) and to \( \x_j \) by @prp-res-plots(a). So
\( \he+\hat{\beta}_j\x_j=\hat{\beta}_j\x_j+\he \) is the sum of a vector in \( \C([\bone,\x_j]) \) and a vector orthogonal to that
space. By @thm-proj-ls-projection the first is the fitted vector and the second the residual.
(b) is @thm-proj-fwl(a) and (b) with \( \X_1=\x_j \).
:::

So the two plots share their line and their scatter about it, and differ only in the horizontal axis. The
added-variable plot uses the part of \( \x_j \) not explained by the other regressors: its spread is the
information for \( \hat{\beta}_j \), and cases far out along it pull on \( \hat{\beta}_j \) (partial leverage,
[Section 20.2](02-leverage.html)). The component-plus-residual plot uses \( \x_j \) on its own scale, the scale on which
a transformation would be chosen, so it shows the *form* of the dependence better. It can mislead when \( \x_j \) is
related nonlinearly to the other regressors.

::: {#exm-res-curvature}
[Curvature seen two ways]

In \( 120 \) simulated cases, \( x_1 \) is uniform on \( (0,4) \), \( x_2 \) is correlated with it (sample correlation
\( 0.87 \)), and \( \E(Y)=1+x_2+0.7(x_1-2)^2 \). The fitted model is linear in both, with coefficient \( -0.144 \) for \( x_1 \), the
slope of both panels of [Figure 20.1.2](#fig-res-partial-residual). A quadratic fitted to the component-plus-residual plot (b)
has leading coefficient \( 0.671 \), close to the true \( 0.7 \). In the added-variable plot (a) the curve is smeared, because
its horizontal axis mixes \( x_1 \) with \( x_2 \), and the corresponding coefficient is \( 0.309 \).
:::

::: {when-format="html"}
![**Figure 20.1.2.** A curved effect of \( x_1 \), fitted linearly: (a) added-variable plot, (b) component-plus-residual
plot. Same slope, same residuals (@prp-res-partial-residual).](partial_residual.svg){#fig-res-partial-residual width=100%}
:::

::: {when-format="pdf"}
![A curved effect of \( x_1 \), fitted linearly: (a) added-variable plot, (b) component-plus-residual
plot. Same slope, same residuals (@prp-res-partial-residual).](partial_residual.pdf){width=100%}
:::

```{.python .run #cell-residual-plots-partial}
rng = np.random.default_rng(11)
n2 = 120
x1 = rng.uniform(0, 4, n2)
x2 = 0.8 * x1 + rng.normal(0, 0.6, n2)             # correlated with x1
y2 = 1 + x2 + 0.7 * (x1 - 2) ** 2 + rng.normal(0, 0.4, n2)
X2 = np.column_stack([np.ones(n2), x1, x2])        # the fitted model is linear in x1
beta2, fitted2, resid2 = ls_fit(X2, y2)

component = resid2 + beta2[1] * x1                 # component-plus-residual for x1
Z = np.column_stack([np.ones(n2), x2])             # everything except x1
x1_tilde = x1 - Z @ np.linalg.lstsq(Z, x1, rcond=None)[0]
y_tilde = y2 - Z @ np.linalg.lstsq(Z, y2, rcond=None)[0]   # added-variable coordinates
print("coefficient of x1:", beta2[1])
print("slope in component-plus-residual plot:", ls_fit(np.column_stack([np.ones(n2), x1]), component)[0][1])
print("slope in added-variable plot:", (x1_tilde @ y_tilde) / (x1_tilde @ x1_tilde))
```

## Exercises

### A. Check your understanding

::: {#exr-res-small-df}
[A1]

A straight line is fitted to \( n=3 \) points. Show that every internally studentized residual satisfies
\( \lvert r_i\rvert=1 \) (unless all three residuals are zero), whatever the data. What does this say about using
\( r_i \) to judge a fit with one residual degree of freedom?
:::

::: {.solution}
Here \( n-p=1 \), so \( \C(\X)\perpc \) is spanned by one unit vector \( \bv \), and \( \he=c\bv \) for a scalar \( c \). Then
\( \hat{\varepsilon}_i=cv_i \), \( s^2=c^2 \) and \( 1-h_{ii}=v_i^2 \) (the \( i \)th diagonal entry of \( \I-\M=\bv\bv\T \)).
So \( r_i^2=c^2v_i^2/(c^2v_i^2)=1 \) whenever \( v_i\ne0 \). With one degree of freedom for error the studentized residuals
carry no information at all about which case fits worst.
:::

::: {#exr-res-weighted-sum}
[A2]

Show that \( \sum_i(1-h_{ii})\,r_i^2=n-p \) for any data. Deduce that if every leverage is at most \( h_{\max}<1 \), then
at most \( (n-p)/\bigl((1-h_{\max})c^2\bigr) \) cases can have \( \lvert r_i\rvert>c \).
:::


### B. Practice

::: {#exr-res-plot-omitted}
[B1]

Let \( \bone\in\C(\X) \), let \( \bz\notin\C(\X) \) be a candidate regressor, let \( \hat{\gamma} \) be its coefficient when it is
added to the model, and let \( R_z^2 \) be the \( R^2 \) of the regression of \( \bz \) on \( \X \). Show that the least
squares slope of the plot of \( \he \) against \( \bz \) is \( \hat{\gamma}\,(1-R_z^2) \). Why is the added-variable plot the
better tool for judging whether \( \bz \) belongs in the model?
:::

::: {.solution}
The slope is \( \sum_i(z_i-\bar z)\hat{\varepsilon}_i/S_{zz}=\bz\T\he/S_{zz} \), since \( \bone\T\he=0 \). Put
\( \tilde{\bz}=(\I-\M)\bz \). Then \( \bz\T\he=\bz\T(\I-\M)\y=\tilde{\bz}\T\y \), and by @thm-proj-fwl(a),
\( \hat{\gamma}=\tilde{\bz}\T\y/\norm{\tilde{\bz}}^2 \). So \( \bz\T\he=\hat{\gamma}\norm{\tilde{\bz}}^2 \), and
\( \norm{\tilde{\bz}}^2=S_{zz}(1-R_z^2) \) is the residual sum of squares of \( \bz \) on \( \X \) (@exr-proj-vif). The slope is
\( \hat{\gamma}(1-R_z^2) \). When \( \bz \) is highly correlated with the columns of \( \X \), the plot against \( \bz \) shrinks the
effect towards zero, while the added-variable plot shows \( \hat{\gamma} \) itself, against the horizontal spread
\( \norm{\tilde{\bz}} \) that actually carries information about it.
:::

::: {#exr-res-misspecified-quadratic}
[B2]

Simple linear regression is fitted to \( x_i=i \), \( i=1,\dots,n \), but \( \E(Y_i)=\beta_0+\beta_1x_i+\gamma x_i^2 \). Use
@prp-res-plots(c) to show that the expected residuals lie on a parabola in \( i \), and find its vertex.
:::


### C. Going deeper

::: {#exr-res-internal-correlation}
[C1]

Assume @eq-opt-normal-model with \( n-p\ge2 \) and leverages below one. Show that
\( \E(r_ir_j)=-h_{ij}/\sqrt{(1-h_{ii})(1-h_{jj})} \) for \( i\ne j \). *Hint:* write \( r_i=\sqrt{n-p}\,\bu_i\T\mathbf{W} \) with
\( \mathbf{W}=(\I-\M)\be/\norm{(\I-\M)\be} \), and show that \( \E(\mathbf{W}\mathbf{W}\T)=(\I-\M)/(n-p) \).
:::

::: {.solution}
With \( \bu_i \) as in the proof of @prp-res-internal, \( r_i=\hat{\varepsilon}_i/(s\sqrt{1-h_{ii}})=\sqrt{n-p}\,\bu_i\T\mathbf{W} \),
because \( \bu_i\T(\I-\M)\be=\bu_i\T\be \). For an orthogonal \( \Q \) with \( \Q(\I-\M)=(\I-\M)\Q \) that maps \( \C(\X)\perpc \)
to itself, \( \Q\be \) has the same law as \( \be \), so \( \Q\mathbf{W} \) has the same law as \( \mathbf{W} \) and
\( \mathbf{S}=\E(\mathbf{W}\mathbf{W}\T) \) satisfies \( \Q\mathbf{S}\Q\T=\mathbf{S} \). Since \( \mathbf{W}\in\C(\X)\perpc \), \( \mathbf{S} \) vanishes on
\( \C(\X) \); restricted to \( \C(\X)\perpc \) it commutes with every rotation of that space and is therefore a multiple of the
identity there. Its trace is \( \E\norm{\mathbf{W}}^2=1 \), so \( \mathbf{S}=(\I-\M)/(n-p) \). Hence
\( \E(r_ir_j)=(n-p)\,\bu_i\T\mathbf{S}\bu_j=\bu_i\T\bu_j=\vect{e}_i\T(\I-\M)\vect{e}_j/\sqrt{(1-h_{ii})(1-h_{jj})} \),
which is the stated value. Since \( \E(r_i)=0 \) and \( \Var(r_i)=1 \), this is also the correlation.
:::
