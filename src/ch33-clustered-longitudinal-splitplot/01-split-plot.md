# Split-plot designs

Some treatments cannot be applied to small units: an irrigation rig waters a strip
of field, not a square metre of it; a teaching method is delivered to a class, not
to a pupil. An experiment combining such a treatment with one applied as finely as
one likes gives the coarse treatment to large units and the fine treatment to
subdivisions of them — a **split-plot** design, named after the field trials it
came from. Two independent randomizations mean two error terms, so that an
analysis ignoring the distinction gets both wrong, in opposite directions.

## The design and its randomization

Take \( r \) blocks. Each is divided into \( a \) **whole plots**, carrying the
\( a \) levels of a factor \( A \) in random order: a randomized complete block
design (@def-dsn-rcbd). Each whole plot is then divided into \( b \) **subplots**,
carrying the \( b \) levels of a factor \( C \) in random order. The two
randomizations are independent, and the second is repeated separately inside every
whole plot. [Figure 33.1.1](#fig-cls-layout) shows one such layout.

::: {when-format="html"}
![**Figure 33.1.1.** A split-plot layout with \( r=4 \) blocks, \( a=3 \) irrigation
regimes assigned at random within each block, and \( b=4 \) varieties assigned at
random within each whole plot.](split_plot_layout.svg){#fig-cls-layout width=100%}
:::

::: {when-format="pdf"}
![A split-plot layout with \( r=4 \) blocks, \( a=3 \) irrigation
regimes assigned at random within each block, and \( b=4 \) varieties assigned at
random within each whole plot.](split_plot_layout.pdf){width=100%}
:::

Two things vary from subplot to subplot: whatever makes whole plots differ, shared
by all \( b \) subplots of a whole plot, and whatever makes subplots differ within
a whole plot, which is not. So the model carries two error terms.

::: {#def-cls-splitplot}
[The balanced split-plot model]

Index blocks by \( i=1,\dots,r \), whole-plot treatments by \( j=1,\dots,a \) and
subplot treatments by \( k=1,\dots,b \), and write \( y_{ijk} \) for the response on
the subplot of block \( i \) receiving whole-plot treatment \( j \) and subplot
treatment \( k \). The **split-plot model** is
\[
y_{ijk}=\mu+\xi_i+\omega_j+\eta_{ij}+\tau_k+\gamma_{jk}+e_{ijk},
\]{#eq-cls-splitplot-model}

where \( \xi_i \), \( \omega_j \), \( \tau_k \) and \( \gamma_{jk} \) are fixed effects, the
\( \eta_{ij} \) are uncorrelated **whole-plot errors** with mean zero and variance
\( \sigma_w^2 \), the \( e_{ijk} \) are uncorrelated **subplot errors** with mean zero and
variance \( \sigma_s^2 \), and the two sets are uncorrelated with each other. The
**normal split-plot model** takes both sets to be normal.
:::

The whole plot \( (i,j) \) is a cluster of \( b \) observations that share \( \eta_{ij} \).
Writing \( \Z_w \) for the \( n\times ra \) matrix of whole-plot indicators, with
\( n=rab \), the combined error \( \eta_{ij}+e_{ijk} \) has
\[
\Cov(\Y)=\V=\sigma_s^2\I+\sigma_w^2\Z_w\Z_w\T
 =\sigma_s^2\I+b\,\sigma_w^2\bP_W ,
\]{#eq-cls-splitplot-cov}

where \( \bP_W=b^{-1}\Z_w\Z_w\T \) replaces every observation by its whole-plot
mean (@exm-proj-oneway-M). Every pair of subplots in a whole plot has correlation
\( \rho=\sigma_w^2/(\sigma_w^2+\sigma_s^2) \), pairs in different whole plots
none. The next section studies this structure in general; here balance makes
everything explicit.

Where does the equicorrelation come from? In a randomized experiment, random
assignment plus unit–treatment additivity *produces* the model rather than
assuming it ([Section 15.1](../ch15-anova-subspaces/01-oneway-model.html)), and
here the argument runs twice: the second randomization is performed anew in every
whole plot, making the subplots of a whole plot exchangeable — exactly the
covariance @eq-cls-splitplot-cov. When the subplot factor is not randomized, as it is not
when that factor is time ([Section 33.3](03-repeated-measures.html)), the
justification disappears and with it the equicorrelation.

## Two strata

Order the observations with \( k \) fastest, then \( j \), then \( i \), and use the
averaging and centring matrices \( \bar{\mathbf{J}}_l \), \( \mathbf{C}_l \) of
[Section 16.1](../ch16-multiway-layouts/01-additive.html). Define
\[
\begin{aligned}
\bP_0&=\bar{\mathbf{J}}_r\otimes\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b, &
\bP_B&=\mathbf{C}_r\otimes\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b, &
\bP_A&=\bar{\mathbf{J}}_r\otimes\mathbf{C}_a\otimes\bar{\mathbf{J}}_b,\\
\bP_{AB}&=\mathbf{C}_r\otimes\mathbf{C}_a\otimes\bar{\mathbf{J}}_b, &
\bP_C&=\bar{\mathbf{J}}_r\otimes\bar{\mathbf{J}}_a\otimes\mathbf{C}_b, &
\bP_{AC}&=\bar{\mathbf{J}}_r\otimes\mathbf{C}_a\otimes\mathbf{C}_b,
\end{aligned}
\]{#eq-cls-splitplot-projections}

and \( \bP_E=\mathbf{C}_r\otimes\I_a\otimes\mathbf{C}_b \), the terms of the balanced
three-way layout with one observation per cell (@thm-tw-balanced): \( \bP_{AB} \)
is the block-by-\( A \) interaction, and \( \bP_E \) collects the block-by-\( C \)
and three-factor interactions. The split-plot analysis rests on the fact that
these seven projections fall into two groups according to how \( \V \) acts on
them.

::: {#lem-cls-strata}
[The two error strata]

Let \( \V \) be as in @eq-cls-splitplot-cov and put
\( \lambda_w=\sigma_s^2+b\,\sigma_w^2 \) and \( \lambda_s=\sigma_s^2 \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \bP_W=\I_r\otimes\I_a\otimes\bar{\mathbf{J}}_b \), and
   \( \V=\lambda_w\bP_W+\lambda_s(\I-\bP_W) \), so \( \V \) has the two eigenvalues
   \( \lambda_w \) and \( \lambda_s \) with eigenspaces \( \C(\bP_W) \) and
   \( \C(\bP_W)\perpc \), of dimensions \( ra \) and \( ra(b-1) \).

2. The six matrices of @eq-cls-splitplot-projections together with \( \bP_E \) are
   mutually orthogonal projections summing to \( \I \), with ranks
   \( 1,\ r-1,\ a-1,\ (r-1)(a-1),\ b-1,\ (a-1)(b-1) \) and \( (r-1)a(b-1) \).

3. \( \bP_0,\bP_B,\bP_A,\bP_{AB} \) have their ranges inside \( \C(\bP_W) \), and
   \( \bP_C,\bP_{AC},\bP_E \) have their ranges inside \( \C(\bP_W)\perpc \). Hence
   \( \V\bP=\lambda_w\bP \) for the first four and \( \V\bP=\lambda_s\bP \) for the last
   three.
:::

:::

::: {.proof}
(a) The whole plot \( (i,j) \) consists of the \( b \) observations with those first two
indices, so averaging within whole plots is \( \I_r\otimes\I_a\otimes\bar{\mathbf{J}}_b \),
which is symmetric and idempotent of rank \( ra \) by @lem-tw-kron(b). Substituting
\( \bP_W \) in @eq-cls-splitplot-cov gives
\( \V=\sigma_s^2(\I-\bP_W)+(\sigma_s^2+b\sigma_w^2)\bP_W \).

(b) Ranks, symmetry and idempotence are @lem-tw-kron(b), noting
\( \bP_E=\mathbf{C}_r\otimes\bar{\mathbf{J}}_a\otimes\mathbf{C}_b
 +\mathbf{C}_r\otimes\mathbf{C}_a\otimes\mathbf{C}_b \). Distinct terms differ in at
least one slot, so they are orthogonal by @lem-tw-kron(c), and they are the
\( 2^3 \) products of @lem-tw-kron(d) regrouped, so they sum to \( \I \).

(c) A product \( \mathbf{E}_1\otimes\mathbf{E}_2\otimes\bar{\mathbf{J}}_b \) satisfies
\( \bP_W\bP=\bP \) by the mixed-product rule and
\( \bar{\mathbf{J}}_b\bar{\mathbf{J}}_b=\bar{\mathbf{J}}_b \); a product
\( \mathbf{E}_1\otimes\mathbf{E}_2\otimes\mathbf{C}_b \) satisfies \( \bP_W\bP=\bzero \)
because \( \bar{\mathbf{J}}_b\mathbf{C}_b=\bzero \). The first four matrices are of the
first kind and the last three of the second. Then
\( \V\bP=\lambda_w\bP_W\bP+\lambda_s(\bP-\bP_W\bP) \) gives the two cases.
:::

Part (c) is the whole idea. The **whole-plot stratum** \( \C(\bP_W) \) carries the
whole-plot error intact, so its variance is \( \lambda_w \); in the **subplot
stratum** \( \C(\bP_W)\perpc \) that error is differenced away and only
\( \sigma_s^2 \) remains.

::: {.idea}
A split-plot design is two experiments sharing one field: a whole-plot experiment
with \( ra \) units and error variance \( \sigma_s^2+b\sigma_w^2 \), and a subplot
experiment with \( rab \) units and error variance \( \sigma_s^2 \). Effects
constant within whole plots belong to the first and effects averaging to zero
within every whole plot to the second — the interaction among them, which is why
such a design studies an interaction well and the coarse main effect badly.
:::

## The analysis

::: {#thm-cls-splitplot}
[The balanced split-plot design]

Assume @eq-cls-splitplot-model with \( r,a,b\ge2 \), and write the cell means
\( \mu_{ijk}=\E(y_{ijk}) \). Let \( \M \) be the orthogonal projection onto the model
space \( \C(\X)=\C(\bP_0+\bP_B+\bP_A+\bP_C+\bP_{AC}) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\V\X)\subseteq\C(\X) \), so least squares is generalized least squares and
   every least squares estimate of an estimable function is its best linear
   unbiased estimator, whatever \( \sigma_w^2 \) and \( \sigma_s^2 \) are.

2. The sums of squares are those of the balanced three-way layout, with
   \( \text{SS}_A=rb\sum_j(\bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot})^2 \),
   \( \text{SS}_C=ra\sum_k(\bar{y}_{\cdot\cdot k}-\bar{y}_{\cdot\cdot\cdot})^2 \),
   \( \text{SS}_{AC}=r\sum_{j,k}(\bar{y}_{\cdot jk}-\bar{y}_{\cdot j\cdot}
       -\bar{y}_{\cdot\cdot k}+\bar{y}_{\cdot\cdot\cdot})^2 \), and the two error sums
   of squares \( \text{SS}_{AB}=\norm{\bP_{AB}\y}^2 \) and
   \( \text{SSE}=\norm{\bP_E\y}^2 \).

3. The expected mean squares (@thm-ss-expected-mean-squares) are
   \[
\begin{aligned}
\E\,\text{MS}_B&=\lambda_w+\frac{ab}{r-1}\sum_i(\xi_i-\bar{\xi})^2, &
\E\,\text{MS}_{AB}&=\lambda_w=\sigma_s^2+b\,\sigma_w^2,\\
\E\,\text{MS}_A&=\lambda_w+\frac{rb}{a-1}\sum_j(\alpha_j-\bar{\alpha})^2, &
\E\,\text{MSE}&=\lambda_s=\sigma_s^2,\\
\E\,\text{MS}_C&=\lambda_s+\frac{ra}{b-1}\sum_k(\kappa_k-\bar{\kappa})^2, &
\E\,\text{MS}_{AC}&=\lambda_s+\frac{r}{(a-1)(b-1)}\sum_{j,k}\delta_{jk}^2,
\end{aligned}
\]{#eq-cls-splitplot-ems}

   where \( \alpha_j=\omega_j+\bar{\gamma}_{j\cdot} \), \( \kappa_k=\tau_k+\bar{\gamma}_{\cdot k} \)
   and \( \delta_{jk}=\gamma_{jk}-\bar{\gamma}_{j\cdot}-\bar{\gamma}_{\cdot k}+\bar{\gamma} \)
   are the equal-weight effects of @def-tw-interaction computed from
   \( \bar{\mu}_{\cdot jk} \).

4. Under normality the six sums of squares are independent,
   \( \text{SS}_{AB}/\lambda_w\sim\chi^2((r-1)(a-1)) \),
   \( \text{SSE}/\lambda_s\sim\chi^2((r-1)a(b-1)) \), and
   \[
\frac{\text{MS}_A}{\text{MS}_{AB}}\sim F\Bigl(a-1,(r-1)(a-1),\frac{rb}{\lambda_w}\sum_j(\alpha_j-\bar{\alpha})^2\Bigr),
\]
   while \( \text{MS}_C/\text{MSE} \) and \( \text{MS}_{AC}/\text{MSE} \) have noncentral
   \( F \) distributions on \( (b-1,(r-1)a(b-1)) \) and
   \( ((a-1)(b-1),(r-1)a(b-1)) \) degrees of freedom with noncentralities
   \( ra\sum_k(\kappa_k-\bar{\kappa})^2/\lambda_s \) and
   \( r\sum_{j,k}\delta_{jk}^2/\lambda_s \). So the whole-plot factor is tested against
   \( \text{MS}_{AB} \) and the subplot factor and the interaction against
   \( \text{MSE} \).

5. A contrast \( \sum_jc_j\alpha_j \) is estimated by \( \sum_jc_j\bar{y}_{\cdot j\cdot} \)
   with variance \( \lambda_w\sum_jc_j^2/(rb) \); a contrast \( \sum_kd_k\kappa_k \) is
   estimated by \( \sum_kd_k\bar{y}_{\cdot\cdot k} \) with variance
   \( \lambda_s\sum_kd_k^2/(ra) \). Both are best linear unbiased, and the two
   estimators are uncorrelated.

6. Fix a subplot treatment \( k \) and a contrast \( \sum_jc_j=0 \). The simple
   effect \( \sum_jc_j\bar{\mu}_{\cdot jk} \) is estimated by
   \( \sum_jc_j\bar{y}_{\cdot jk} \), best linear unbiased, with variance
   \[
\Bigl(\sum_jc_j^2\Bigr)\frac{\lambda_w+(b-1)\lambda_s}{rb}
 =\Bigl(\sum_jc_j^2\Bigr)\frac{\sigma_w^2+\sigma_s^2}{r},
\]{#eq-cls-splitplot-simple}

   which belongs to neither stratum: no single mean square estimates it, and the
   unbiased estimate is the combination
   \( (\sum_jc_j^2)\{\text{MS}_{AB}+(b-1)\text{MSE}\}/(rb) \).
:::

:::

::: {.proof}
(a) Every \( \bv\in\C(\X) \) splits as \( \bv=\bv_1+\bv_2 \) with
\( \bv_1\in\C(\bP_0+\bP_B+\bP_A) \) and \( \bv_2\in\C(\bP_C+\bP_{AC}) \), because the five
projections are mutually orthogonal and sum to \( \M \). By @lem-cls-strata(c),
\( \V\bv=\lambda_w\bv_1+\lambda_s\bv_2 \), again in \( \C(\X) \). So
\( \C(\V\X)\subseteq\C(\X) \), and @thm-proj-kruskal says the ordinary and the
\( \V^{-1} \) projection onto \( \C(\X) \) coincide. The estimator is therefore both the
least squares and the generalized least squares estimator, which is the best linear
unbiased estimator by @cor-opt-aitken. (This is the balanced case
of @thm-ggm-ols-blue.)

(b) Each projection averages over the slots carrying \( \bar{\mathbf{J}} \) and centres
in the slots carrying \( \mathbf{C} \). For instance
\( \bP_A=\bar{\mathbf{J}}_r\otimes\I_a\otimes\bar{\mathbf{J}}_b
 -\bar{\mathbf{J}}_r\otimes\bar{\mathbf{J}}_a\otimes\bar{\mathbf{J}}_b \) puts
\( \bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot\cdot\cdot} \) in position \( (i,j,k) \), a value
repeated \( rb \) times, giving \( \text{SS}_A \); the same expansion gives
\( \text{SS}_C \) and \( \text{SS}_{AC} \), as in @thm-tw-balanced(b). The error sums
of squares are the squared lengths of the remaining components.

(c) Each \( \bP \) in the list satisfies \( \V\bP=\lambda\bP \) with \( \lambda \) the
eigenvalue of its stratum, so
\( \E(\Y\T\bP\Y)=\tr(\bP\V)+\bmu\T\bP\bmu=\lambda\,\rank(\bP)+\norm{\bP\bmu}^2 \)
by @thm-rv-quadform-mean, where \( \bmu \) is the mean vector. Dividing by the rank
gives \( \E(\text{MS})=\lambda+\norm{\bP\bmu}^2/\rank(\bP) \). It remains to evaluate
\( \bP\bmu \). Since \( \bmu \) has entries
\( \mu+\xi_i+\omega_j+\tau_k+\gamma_{jk} \), the projection \( \bP_A \) returns
\( \bar{\mu}_{\cdot j\cdot}-\bar{\mu}_{\cdot\cdot\cdot}=\alpha_j-\bar{\alpha} \) on each of
the \( rb \) observations with whole-plot treatment \( j \); likewise \( \bP_C \) returns
\( \kappa_k-\bar{\kappa} \) on each of \( ra \) observations, \( \bP_{AC} \) returns
\( \delta_{jk} \) on each of \( r \), and \( \bP_B \) returns \( \xi_i-\bar{\xi} \) on each of
\( ab \). Finally \( \bP_{AB}\bmu=\bzero \) and \( \bP_E\bmu=\bzero \), because \( \bmu \) has
no block-by-treatment terms.

(d) Write \( \V^{1/2}=\lambda_w^{1/2}\bP_W+\lambda_s^{1/2}(\I-\bP_W) \), a symmetric
positive definite square root (@thm-mat-square-root), and put
\( \bz=\V^{-1/2}\Y\sim\Normal_n(\V^{-1/2}\bmu,\I) \) (@prp-rv-whitening). For any
\( \bP \) in the list, @lem-cls-strata(c) gives \( \V^{1/2}\bP\V^{1/2}=\lambda\bP \), so
\( \Y\T\bP\Y=\lambda\,\bz\T\bP\bz \). The projections are mutually orthogonal, so
@thm-qf-orthogonal-projections applied to \( \bz \) makes the quadratic forms
\( \bz\T\bP\bz \) independent with noncentral chi-squared distributions on
\( \rank(\bP) \) degrees of freedom and noncentrality
\( \norm{\bP\V^{-1/2}\bmu}^2=\norm{\bP\bmu}^2/\lambda \), the last equality because
\( \V^{-1/2} \) acts on \( \C(\bP) \) as multiplication by \( \lambda^{-1/2} \). Forming
ratios and cancelling \( \lambda \) between numerator and denominator *within the same
stratum* gives the stated \( F \) distributions (@def-qf-noncentral-f). Ratios across
strata would not cancel, which is precisely why \( \text{MS}_A/\text{MSE} \) is not an
\( F \) statistic.

(e) The coefficient vector of \( \sum_jc_j\bar{y}_{\cdot j\cdot} \) is
\( \boldsymbol{\uprho}=(rb)^{-1}\bone_r\otimes\mathbf{c}\otimes\bone_b\in\C(\bP_A) \),
since \( \mathbf{c}\perp\bone_a \), so by @lem-cls-strata(c)
\( \Var(\boldsymbol{\uprho}\T\Y)=\lambda_w\norm{\boldsymbol{\uprho}}^2
 =\lambda_w\sum_jc_j^2/(rb) \). The same computation with
\( \boldsymbol{\uprho}_C=(ra)^{-1}\bone_r\otimes\bone_a\otimes\mathbf{d}\in\C(\bP_C) \)
gives the second variance, and the two are uncorrelated because
\( \boldsymbol{\uprho}\T\boldsymbol{\uprho}_C=0 \). Both are least squares estimates
of estimable functions, hence best linear unbiased by (a).

(f) Split \( \vect{e}_k=b^{-1}\bone_b+(\vect{e}_k-b^{-1}\bone_b) \). The coefficient
vector \( r^{-1}\bone_r\otimes\mathbf{c}\otimes\vect{e}_k \) splits accordingly into
\( \boldsymbol{\uprho}_1\in\C(\bP_A) \) and
\( \boldsymbol{\uprho}_2\in\C(\bP_{AC}) \), of squared lengths
\( \sum_jc_j^2/(rb) \) and \( (b-1)\sum_jc_j^2/(rb) \). They lie in different
strata and are orthogonal, so the variance is
\( \lambda_w\norm{\boldsymbol{\uprho}_1}^2+\lambda_s\norm{\boldsymbol{\uprho}_2}^2 \),
which is @eq-cls-splitplot-simple; the estimate of it is unbiased
by @eq-cls-splitplot-ems, and the estimator is best linear unbiased by (a).
:::

Dropping \( (\xi\omega)_{ij} \) from a full three-way model is what turns the
block-by-\( A \) interaction into the error term for \( A \), on only
\( (r-1)(a-1) \) degrees of freedom: the formal content of "whole plots are the
experimental units for \( A \)". And \( \text{MS}_A \) tests a contrast in
\( \alpha_j=\omega_j+\bar{\gamma}_{j\cdot} \), the whole-plot effect averaged over the
subplot treatments used, not in \( \omega_j \)
alone ([Section 16.2](../ch16-multiway-layouts/02-interaction.html)).

Part (f) is the case most often got wrong. "The interaction is tested against
\( \text{MSE} \)" tempts one to compare two whole-plot treatments *within a single
subplot treatment* against \( \text{MSE} \) too; but that comparison carries the
whole-plot error intact, and @eq-cls-splitplot-simple prices it between the two
strata. Its estimated variance combines two mean squares, so no exact \( t \) or
\( F \) test exists, and Satterthwaite's rule (@prp-mix-inference) gives it
\[
\nu=\frac{\{\text{MS}_{AB}+(b-1)\text{MSE}\}^2}
   {\text{MS}_{AB}^2/\{(r-1)(a-1)\}+(b-1)^2\text{MSE}^2/\{(r-1)a(b-1)\}}
\]{#eq-cls-splitplot-satterthwaite}

approximate degrees of freedom, which lie between the two exact values.

::: {#exm-cls-irrigation}
[An irrigation and variety trial]

Four fields are each divided into three whole plots, and the regimes flood,
sprinkler and drip are assigned at random within each field; each whole plot is
divided into four subplots carrying four barley varieties in random order, giving
\( n=48 \) yields in tonnes per hectare, simulated from @eq-cls-splitplot-model
with \( \sigma_w=0.35 \) and \( \sigma_s=0.25 \) t/ha. The
analysis of variance is

| Source | df | SS | MS | \( F \) | \( p \) |
|:---|---:|---:|---:|---:|---:|
| blocks | 3 | 9.471 | 3.157 | | |
| irrigation \( A \) | 2 | 9.215 | 4.608 | 10.637 | 0.0106 |
| whole-plot error | 6 | 2.599 | 0.433 | | |
| variety \( C \) | 3 | 1.755 | 0.585 | 9.442 | 0.0002 |
| irrigation × variety | 6 | 0.958 | 0.160 | 2.577 | 0.0419 |
| subplot error | 27 | 1.673 | 0.062 | | |

The two error mean squares differ by a factor of seven. Solving
@eq-cls-splitplot-ems gives \( \hat{\sigma}_s=\sqrt{\text{MSE}}=0.2489 \) and
\( \hat{\sigma}_w=\sqrt{(\text{MS}_{AB}-\text{MSE})/b}=0.3046 \), close to the
generating values.

The irrigation means are \( 4.8075 \), \( 5.8000 \) and
\( 5.6575 \) t/ha; a difference of two has standard error
\( \sqrt{2\,\text{MS}_{AB}/(rb)}=0.2327 \) on \( 6 \) degrees of freedom, so
drip beats flood by \( 0.8500 \) t/ha with \( t=3.653 \). The variety
means are \( 5.134 \), \( 5.364 \), \( 5.574 \) and
\( 5.614 \), and a difference of two has standard error
\( \sqrt{2\,\text{MSE}/(ra)}=0.1016 \) on \( 27 \) degrees of freedom — more
than twice as precise, from the same \( 48 \) observations. Comparing drip with
flood *for one variety* is neither: by @eq-cls-splitplot-simple its standard error
is \( \sqrt{\{\text{MS}_{AB}+3\,\text{MSE}\}/8}=0.2782 \), on
\( 11.77 \) Satterthwaite degrees of freedom.
:::

```{.python .run #cell-split-plot-simulate}
import numpy as np

r, a, b = 4, 3, 4                                   # blocks, irrigation regimes, varieties
n = r * a * b
block_eff = np.array([-0.30, 0.10, 0.05, 0.15])     # field-to-field fertility
irrig_eff = np.array([0.00, 0.45, 0.70])            # flood, sprinkler, drip
variety_eff = np.array([-0.30, -0.10, 0.15, 0.25])
inter = np.array([[0.10, -0.05, -0.10, 0.05],       # irrigation x variety
                  [-0.05, 0.10, 0.00, -0.05],
                  [-0.05, -0.05, 0.10, 0.00]])
sigma_w, sigma_s = 0.35, 0.25                       # whole-plot and subplot errors

rng = np.random.default_rng(3301)
eta = rng.normal(0, sigma_w, (r, a))                # one error per whole plot
e = rng.normal(0, sigma_s, (r, a, b))               # one error per subplot
mean = (5.0 + block_eff[:, None, None] + irrig_eff[None, :, None]
        + variety_eff[None, None, :] + inter[None, :, :])
Y = np.round(mean + eta[:, :, None] + e, 2)         # yields, t/ha
y = Y.reshape(-1)                                   # stacked block, then regime, then variety
```

The analysis is seven inner products, each projection a Kronecker product as in
[Section 16.3](../ch16-multiway-layouts/03-balanced.html), with two denominators
in place of one.

```{.python .run #cell-split-plot-anova}
from scipy import stats

def Jbar(k):
    return np.full((k, k), 1.0 / k)

def Cen(k):
    return np.eye(k) - Jbar(k)

def kron3(E1, E2, E3):
    return np.kron(np.kron(E1, E2), E3)

P = {"blocks": kron3(Cen(r), Jbar(a), Jbar(b)),
     "irrigation": kron3(Jbar(r), Cen(a), Jbar(b)),
     "whole-plot error": kron3(Cen(r), Cen(a), Jbar(b)),
     "varieties": kron3(Jbar(r), Jbar(a), Cen(b)),
     "irrigation x variety": kron3(Jbar(r), Cen(a), Cen(b)),
     "subplot error": kron3(Cen(r), np.eye(a), Cen(b))}
df = {k: int(round(np.trace(M))) for k, M in P.items()}
SS = {k: y @ M @ y for k, M in P.items()}
MS = {k: SS[k] / df[k] for k in P}

F_irrigation = MS["irrigation"] / MS["whole-plot error"]
F_variety = MS["varieties"] / MS["subplot error"]
F_inter = MS["irrigation x variety"] / MS["subplot error"]
for name, F, d in [("irrigation", F_irrigation, df["whole-plot error"]),
                   ("varieties", F_variety, df["subplot error"]),
                   ("irrigation x variety", F_inter, df["subplot error"])]:
    p = stats.f.sf(F, df[name], d)
    print(f"{name:22s} F({df[name]},{d}) = {F:7.3f}   p = {p:.4f}")
```

## What pooling costs

Pool the two error sums of squares into one residual mean square on \( 33 \)
degrees of freedom — in @exm-cls-irrigation, \( 0.129 \) — and by
@eq-cls-splitplot-ems the consequences are opposite in the two strata. The
whole-plot denominator falls towards \( \lambda_s \), inflating the irrigation
\( F \) to \( 35.591 \) on \( (2,33) \) degrees of freedom,
\( p=5.78\times 10^{-9} \) against \( 0.0106 \), and shrinking the standard
error of an irrigation difference from \( 0.2327 \) to \( 0.1272 \), a
factor of \( 1.829 \) too small; the subplot denominator exceeds
\( \lambda_s \), and that of a variety difference rises from
\( 0.1016 \) to \( 0.1469 \), a factor of \( 1.445 \) too large.
Pooling manufactures significance for the treatment the design studies badly and
destroys it for the one it studies well.

::: {.warning}
The two error terms are a property of the *randomization*, not of the data. A
split-plot experiment analysed as a completely randomized three-factor experiment
produces a plausible table with one residual line, and nothing in the numbers
announces the mistake. Wherever one factor is applied to groups of units and
another within them, look for two errors.
:::

## The same analysis as a mixed model

@eq-cls-splitplot-model is a linear mixed model (@def-mix-model) with
\( \Z=\Z_w \), \( \G=\sigma_w^2\I_{ra} \) and \( \R=\sigma_s^2\I_n \), whose
marginal covariance \( \Z\G\Z\T+\R \) is @eq-cls-splitplot-cov. In the balanced
case, fitting it by restricted maximum likelihood (@def-mix-reml) reproduces the
analysis above exactly, provided \( \text{MS}_{AB}\ge\text{MSE} \): the restricted
likelihood depends on the data only through the two mean squares, and is maximized
at the moment estimates of @exm-cls-irrigation (@exr-cls-splitplot-reml).
Satterthwaite's rule (@prp-mix-inference) then returns exactly \( (r-1)(a-1) \)
degrees of freedom for a whole-plot contrast, and
@eq-cls-splitplot-satterthwaite for a simple effect.

The route earns its keep when the balance fails: a missing subplot, a whole plot
lost to flooding, a covariate measured at the subplot level. The strata are then
no longer orthogonal, no mean square estimates \( \lambda_w \) or \( \lambda_s \)
on its own, and Henderson's equations (@thm-mix-henderson) do the computing. Two
variations change nothing essential: the whole-plot design may be completely
randomized or a Latin square (@def-dsn-latin-square), which changes only the
whole-plot degrees of freedom, and splitting can be repeated, giving three
strata (@exr-cls-split-split).

## Exercises

### A. Check your understanding

::: {#exr-cls-strata-membership}
[A1]

For the design of @exm-cls-irrigation, say which stratum each of the following
belongs to, and on how many degrees of freedom its error is estimated: the
difference between two block means; between two variety means; between two
irrigation means; between two varieties within the drip regime; the difference in
the effect of drip versus flood between two varieties; and the difference between
drip and flood for the first variety alone.
:::

::: {.solution}
The block and irrigation differences are constant within whole plots, so they lie
in the whole-plot stratum and use \( \text{MS}_{AB} \) on \( 6 \) degrees of freedom.
The next three are contrasts summing to zero within each whole plot, so they lie
in the subplot stratum and use \( \text{MSE} \) on \( 27 \) degrees of freedom. The
last lies in neither: it is @thm-cls-splitplot(f), with estimated variance
\( \{\text{MS}_{AB}+3\,\text{MSE}\}/8 \) on about \( 11.8 \) Satterthwaite
degrees of freedom.
:::

::: {#exr-cls-more-subplots}
[A2]

With \( r \), \( a \), \( \sigma_w^2 \) and \( \sigma_s^2 \) fixed, what happens to the
variance of a whole-plot contrast, and to that of a subplot contrast, as \( b \)
grows? Which kind of comparison does adding subplots help?
:::

### B. Practice

::: {#exr-cls-splitplot-kruskal}
[B1]

Verify @thm-cls-splitplot(a) directly from @thm-proj-kruskal by computing
\( \V\X \) for the model matrix
\( \X=[\bone,\mathbf{B},\mathbf{A},\mathbf{C},\mathbf{AC}] \) of indicator columns.
Then show that adding the block-by-\( A \) columns leaves the condition true but
destroys the analysis, and say why.
:::

::: {#exr-cls-splitplot-reml}
[B2]

Show that for the balanced split-plot model the restricted likelihood of
\( (\sigma_w^2,\sigma_s^2) \) depends on the data only through \( \text{SS}_{AB} \) and
\( \text{SSE} \), and that it is maximized at \( \hat{\lambda}_w=\text{MS}_{AB} \),
\( \hat{\lambda}_s=\text{MSE} \) provided \( \text{MS}_{AB}\ge\text{MSE} \). *Hint:*
the error contrasts span \( \C(\bP_{AB})\dirsum\C(\bP_E) \), and \( \V \) acts as a
scalar on each.
:::

::: {.solution}
The error contrasts are the \( \mathbf{a}\T\y \) with
\( \mathbf{a}\in\C(\X)\perpc=\C(\bP_{AB})\dirsum\C(\bP_E) \). Take orthonormal bases
of the two summands, of sizes \( f_w=(r-1)(a-1) \) and \( f_s=(r-1)a(b-1) \). By
@lem-cls-strata(c) the transformed data are independent, \( \Normal(0,\lambda_w) \)
for the first basis and \( \Normal(0,\lambda_s) \) for the second. By @def-mix-reml
the restricted log likelihood is therefore
\( -\tfrac12\{f_w\log\lambda_w+\text{SS}_{AB}/\lambda_w
 +f_s\log\lambda_s+\text{SSE}/\lambda_s\} \) up to a constant, a sum of two
ordinary normal-variance likelihoods maximized at the two mean squares. Since
\( \sigma_s^2=\lambda_s \) and \( \sigma_w^2=(\lambda_w-\lambda_s)/b \), the constraint
\( \sigma_w^2\ge0 \) is \( \text{MS}_{AB}\ge\text{MSE} \).
:::

### C. Going deeper

::: {#exr-cls-split-split}
[C1]

A split-split-plot design divides each subplot into \( c \) sub-subplots and
randomizes a third factor \( D \) among them, adding an error \( \zeta_{ijk} \) of
variance \( \sigma_u^2 \) shared by the sub-subplots of a subplot. Write down the
covariance matrix in the style of @eq-cls-splitplot-cov, find its three
eigenvalues and eigenspaces, and state which effects are tested against which of
the three error mean squares.
:::
