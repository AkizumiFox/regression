# Varying-coefficient models

The additive model of @def-add-model has no interactions: whatever \( z_2 \) does, the
contribution of \( z_1 \) is the same curve. Throwing that restriction away brings back the
curse of dimensionality, so relax it in one direction at a time. The most useful way keeps a
*linear* effect and lets its slope move.

## The model

::: {#def-add-varying}
[Varying-coefficient model]

Let \( u_1,\dots,u_n \) be the values of an **interaction variable** and
\( z_1,\dots,z_n \) those of an **effect modifier**. A **varying-coefficient term** is
\[
f(z_i)\,u_i,\qquad \sum_{i=1}^{n}f(z_i)=0 ,
\]{#eq-add-varying}

and the **varying-coefficient model** is a structured additive
model (@def-add-model) in which one or more of the terms have this form. Combined with a
parametric column \( \bu=(u_1,\dots,u_n)\T \) carrying a coefficient \( \beta_u \), the
total effect of \( \bu \) on the response is
\[
\E(y_i)=\dots+\{\beta_u+f(z_i)\}\,u_i ,
\]{#eq-add-varying-slope}

so that \( \beta_u+f(z) \) is the **coefficient function**: the slope in \( u \) at effect
modifier value \( z \). With \( u_i\equiv1 \) the term is an ordinary smooth main effect,
and with \( f \) constrained to be constant it is an ordinary linear effect.
:::

Three readings of the same object:

- **An interaction**, the smooth analogue of the product term \( \beta uz \) of
  [Section 5.7](../ch05-model-and-least-squares/07-interpreting-coefficients.html). Where a
  product term says the slope in \( u \) is *linear* in \( z \), this one says only that
  it is smooth.
- **A factor-by-curve model**, when \( \bu \) indicates a group: each group gets its own
  deviation curve in \( z \).
- **A time-varying coefficient**, when \( z \) is calendar time: the relationship has
  drifted, and the model does not say how.

## Estimation is the same computation

Varying coefficients cost nothing new because the design of the term is a row-scaled
basis.

::: {#prp-add-varying-fit}
[A varying-coefficient term is an ordinary penalized term]

Let \( B_1,\dots,B_d \) be a basis for the effect modifier and \( \mathbf{K}_0 \) the
penalty that the corresponding smooth main effect would use. Put
\( (\Z)_{il}=u_iB_l(z_i) \).

::: {.enumerate options="label=(\alph*)"}
1. With \( f=\sum_l\gamma_lB_l \), the vector of the term's contributions is
   \( \{f(z_i)u_i\}_i=\Z\bgamma \), and the roughness of \( f \) is
   \( \bgamma\T\mathbf{K}_0\bgamma \) as before. The term is therefore the pair
   \( (\Z,\mathbf{K}_0) \), and everything in
   [Sections 44.1](01-additive-models.html) and [44.2](02-backfitting.html) applies
   unchanged.

2. The constraint @eq-add-varying is \( \mathbf{c}\T\bgamma=0 \) with
   \( \mathbf{c}=\bigl(\sum_iB_1(z_i),\dots,\sum_iB_d(z_i)\bigr)\T \) — the *unweighted*
   column sums of the basis, not of \( \Z \) — and is absorbed by
   @prp-add-identifiability(c). Under it, \( \beta_u \) in @eq-add-varying-slope is the
   average of the coefficient function over the design points,
   \( \beta_u=n^{-1}\sum_i\{\beta_u+f(z_i)\} \).

3. If \( u_i\ne0 \) for at least \( d \) design points with distinct \( z_i \) and the
   basis is linearly independent there, \( \Z \) has full column rank. Observations with
   \( u_i=0 \) contribute nothing to the term: the coefficient function is estimated only
   where the interaction variable is nonzero.
:::

:::

::: {.proof}
(a) The \( i \)th entry of \( \Z\bgamma \) is \( \sum_l\gamma_lu_iB_l(z_i)=u_if(z_i) \).
The penalty is a statement about \( f \), so it is unchanged by multiplying the term by
\( u_i \).
(b) \( \sum_if(z_i)=\sum_l\gamma_l\sum_iB_l(z_i)=\mathbf{c}\T\bgamma \), which is the
stated constraint; @prp-add-identifiability(c) applies verbatim with \( \mathbf{c} \) in
place of \( \Z\T\bone \). The averaging statement is then
\( n^{-1}\sum_i\{\beta_u+f(z_i)\}=\beta_u+n^{-1}\mathbf{c}\T\bgamma=\beta_u \).
(c) Suppose \( \Z\bgamma=\bzero \). Then \( u_if(z_i)=0 \) for all \( i \), so
\( f(z_i)=0 \) at the \( d \) design points with \( u_i\ne0 \); linear independence of the
basis at those points forces \( \bgamma=\bzero \). The last sentence is the same
computation read forwards: rows with \( u_i=0 \) are zero rows of \( \Z \).
:::

::: {.warning}
Constraining the *weighted* sum \( \sum_iu_if(z_i)=0 \) instead — what centring the
columns of \( \Z \) does — also identifies the model, but \( \beta_u \) is then not the
average slope, and if \( \bu \) has been centred it means nothing at all, since
\( \sum_iu_i=0 \) leaves the constraint silent about the level of \( f \). Use (b).
:::

## An example, and how far it can be pushed

::: {#exm-add-statecrime}
[Does poverty matter more in urban states?]

@exm-proj-fwl-crime fitted the murder rate of the \( 50 \) US states in 2009
on the poverty rate, the single-parent percentage and the urban percentage, reporting one
coefficient for poverty. Here that coefficient may vary smoothly with how urban the state
is:
\[
\text{murder}_i=\beta_0+\beta_1\,\text{single}_i
+\{\beta_2+f(\text{urban}_i)\}\,(\text{poverty}_i-\overline{\text{poverty}})
+f_0(\text{urban}_i)+\varepsilon_i ,
\]
with P-spline terms of \( 8 \) basis functions for \( f \) and \( f_0 \) and
smoothing parameters chosen by the restricted-likelihood updates of
[Section 44.5](05-structured-additive.html). The fit uses
\( 6.83 \) effective degrees of freedom in all, of which
\( 1.82 \) go to the smooth main effect of urbanization and
\( 2.01 \) to the coefficient function.

The average poverty slope is \( 0.250 \). Keeping the smooth main effect of urbanization
but forcing the poverty slope to be constant gives \( 0.231 \), and
@exm-proj-fwl-crime, in which urbanization entered linearly as well, gave
\( 0.2538 \). The estimated
coefficient function rises from \( -0.034 \) at \( 20\% \) urban to
\( 0.364 \) at \( 60\% \) and falls back to \( 0.152 \)
at \( 90\% \): on its face, a story in which poverty is associated with homicide in
middling states and not in the most rural ones.

That story does not survive its standard errors. The band of \( \pm2 \) standard errors in
[Figure 44.3.1](03-varying-coefficients.html#fig-add-varying) contains the constant
coefficient over the whole range, and the contrast between \( 60\% \) and \( 20\% \)
urban is \( 0.398 \) with a standard error of \( 0.234 \), a
ratio of \( 1.70 \). Fifty states and a residual standard deviation of
\( 1.418 \) murders per \( 100{,}000 \) cannot resolve a coefficient
function; they can only say that one is not needed.
:::

::: {when-format="html"}
![**Figure 44.3.1.** A varying coefficient on 50 US states. (a) the coefficient of poverty against
the urban percentage, two standard errors wide; dashed, the constant slope from the model
that keeps the smooth urban effect. (b) fitted murder rate against poverty at three levels
of urbanization.](varying_coefficient.svg){#fig-add-varying width=100%}
:::

::: {when-format="pdf"}
![A varying coefficient on 50 US states. (a) the coefficient of poverty against
the urban percentage, two standard errors wide; dashed, the constant slope from the model
that keeps the smooth urban effect. (b) fitted murder rate against poverty at three levels
of urbanization.](varying_coefficient.pdf){width=100%}
:::

Fitting takes no new code, because @prp-add-varying-fit(a) has made this an ordinary term:
the design is the P-spline basis in the urban percentage scaled row by row by the centred
poverty column, the penalty is the same second-difference matrix, and the same solver
handles it. The band is read off the joint covariance of \( (\beta_2,\bgamma) \) with
\( \beta_2+f(z) \) evaluated on a grid; the code is
`code/ch44/varying_coefficient.py`.

::: {.remark}
[What a coefficient function means]

@exm-add-statecrime is an observational comparison across states, and
[Chapter 25](../ch25-causal-interpretation/index.html) applies to it word for word. A
varying coefficient describes how an *association* changes with the modifier. Reading it as
"raising poverty by one point in a \( 60\% \)-urban state would raise the murder rate by
\( 0.36 \)" requires everything @thm-cau-adjustment requires, and requires it separately at
each level of the modifier, since the confounding structure may differ across it. The extra
flexibility buys a richer description, not weaker causal assumptions.
:::

## Where varying coefficients take you

Two extensions use nothing beyond @prp-add-varying-fit. With \( z \) time and \( \bu \)
a covariate observed repeatedly, \( f(t_i)u_i \) is a coefficient drifting over the study;
for series data the band must then come from
[Part VII](../ch31-general-gauss-markov/index.html) rather than from the independence
assumption of @def-add-model (@exr-add-time-varying). With \( z \) a location and a spatial
term from [Section 44.4](04-spatial-effects.html) in place of the P-spline, the coefficient
varies over a map, which is disease mapping: the same code with a different penalty.

When neither covariate deserves to be the linear one, a genuine bivariate term
\( f(z_1,z_2) \) is needed. The standard construction is a **tensor-product** basis
\( B_{lm}(z_1,z_2)=B_l^{(1)}(z_1)B_m^{(2)}(z_2) \) with the penalty
\[
\mathbf{K}=\lambda_1\bigl(\I\otimes\mathbf{K}^{(1)}\bigr)
+\lambda_2\bigl(\mathbf{K}^{(2)}\otimes\I\bigr),
\]{#eq-add-tensor}

penalizing roughness along each coordinate separately, and so invariant to the units of the
two covariates. Everything in this chapter goes through for it, with two smoothing
parameters instead of one. The book does not develop the construction further; how to scale
the two penalties, how to separate the interaction from its main effects and how to display
the surface belong to a book on smoothing, and Wood (2017) is where to find them.
@exr-add-tensor builds the smallest case by hand.

## Exercises

### A. Check your understanding

::: {#exr-add-varying-special}
[A1]

What does a varying-coefficient term become when (i) \( u_i\equiv1 \); (ii) \( f \) is
forced to be constant; (iii) \( f \) is forced to be linear in \( z \)? Which of the three
has a familiar name from [Part II](../ch05-model-and-least-squares/index.html)?
:::

### B. Practice

::: {#exr-add-time-varying}
[B1]

A varying-coefficient model is fitted to a time series, but the errors follow
\( \Cov(\be)=\sigma^2\V \) with \( \V\ne\I \). Show that the penalized estimate
\( \hat{\bgamma}=\A^{-1}\Z\T\y \) is still unbiased for the penalized target
\( \A^{-1}\Z\T\E(\Y) \), but that its covariance is
\( \sigma^2\A^{-1}\Z\T\V\Z\A^{-1} \), not \( \sigma^2\A^{-1}\Z\T\Z\A^{-1} \). If \( \V \)
has positive autocorrelations, in which direction is the naive band wrong?
:::

::: {.solution}
Linearity gives both at once:
\( \E(\A^{-1}\Z\T\Y)=\A^{-1}\Z\T\E(\Y) \) and
\( \Cov(\A^{-1}\Z\T\Y)=\A^{-1}\Z\T\Cov(\Y)\Z\A^{-1} \) (@thm-rv-linear). Under positive
autocorrelation \( \Z\T\V\Z-\Z\T\Z \) is typically nonnegative definite for smooth
columns, so the true variance exceeds the naive one and the band is too narrow, as
@thm-dep-covariance found for a linear model. The repair is
[Chapter 31](../ch31-general-gauss-markov/index.html), or a random-effect term for the
serial structure as in
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html).
:::

::: {#exr-add-varying-concurvity}
[B2]

Show that a model containing both a smooth main effect \( f_0(z) \) and a
varying-coefficient term \( f(z)u \) has exact concurvity (@def-add-concurvity-index) as
soon as \( u_i \) is itself a function of \( z_i \), say \( u_i=h(z_i) \) with \( h \) in
the span of the basis. What does this say about the model of @exm-add-statecrime?
:::

::: {.solution}
If \( u_i=h(z_i) \) with \( h \) in the span, the varying term contributes
\( f(z_i)h(z_i) \), a function of \( z_i \) alone; whenever \( fh \) is again in the span,
the main-effect term realizes it too, and subtracting gives a nontrivial relation. In
@exm-add-statecrime the poverty rate is *not* a function of the urban percentage, so the
concurvity is partial rather than exact; but the two are correlated across states, so the
terms compete, which is part of why the bands are wide.
:::

### C. Going deeper

::: {#exr-add-tensor}
[C1]

Build the smallest tensor-product term by hand: take \( d_1=d_2=3 \) basis functions in
each coordinate and first-difference penalties \( \mathbf{K}^{(1)}=\mathbf{K}^{(2)} \).
Write out the \( 9\times9 \) penalty of @eq-add-tensor for \( \lambda_1=\lambda_2=1 \) and
find its null space. Why does that null space make the term identifiable only after the
main effects have been constrained away?
:::
