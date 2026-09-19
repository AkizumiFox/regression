# The one-way layout as a linear model

Earlier chapters used the one-way layout as an example of rank deficiency ([Chapter 6](../ch06-projections/index.html)),
estimability ([Chapter 8](../ch08-estimability/index.html)), orthogonal decomposition
([Chapter 9](../ch09-sums-of-squares/index.html)) and simultaneous inference ([Chapter 13](../ch13-multiplicity/index.html)). Here it is
the subject. This section fixes the model, recalls its projection, and settles what its parameters mean and where they come from.

## A designed experiment

An agronomist wants to know how the yield of a wheat variety responds to nitrogen fertilizer. Thirty
field plots of equal size are available. Five application rates are chosen, \( 0 \), \( 40 \), \( 80 \), \( 120 \)
and \( 160 \) kilograms of nitrogen per hectare, and each rate is assigned to six plots chosen at random
from the thirty. At harvest the grain yield of each plot is recorded in tonnes per hectare. This is a
**completely randomized design** with one **factor** (nitrogen) at five **levels**, also called
**treatments**, with six **replicates** of each. The data used throughout this chapter were simulated
from a smooth response curve with plot-to-plot standard deviation \( 0.5 \) t/ha, which the analysis
does not know about.

| Rate (kg N/ha) | Yields (t/ha) | Mean | Standard deviation |
|:---:|:---|:---:|:---:|
| 0 | 3.69, 4.10, 3.41, 3.67, 4.54, 4.62 | \( 4.005 \) | \( 0.498 \) |
| 40 | 5.57, 5.21, 4.98, 4.99, 5.77, 4.90 | \( 5.237 \) | \( 0.357 \) |
| 80 | 5.86, 5.80, 6.64, 5.99, 5.72, 6.33 | \( 6.057 \) | \( 0.357 \) |
| 120 | 6.59, 6.85, 6.74, 6.32, 6.39, 6.20 | \( 6.515 \) | \( 0.254 \) |
| 160 | 6.60, 5.14, 7.82, 7.19, 6.61, 6.55 | \( 6.652 \) | \( 0.889 \) |

[Figure 15.1.1](#fig-aov-nitrogen-data) shows the thirty yields. Mean yield rises with the rate, steeply at first
and then more slowly, and the plots within a rate scatter by roughly half a tonne. The later sections compare particular
treatments ([Section 15.2](02-contrasts.html)), ask whether the rates differ at all ([Section 15.3](03-anova-table.html)), and
describe the shape of the response ([Section 15.4](04-orthogonal-contrasts.html)).

::: {when-format="html"}
![**Figure 15.1.1.** The nitrogen trial: yields of thirty plots at five nitrogen rates, six plots per rate
(slightly spread horizontally), with the rate means and the grand mean.](nitrogen_data.svg){#fig-aov-nitrogen-data width=66%}
:::

::: {when-format="pdf"}
![The nitrogen trial: yields of thirty plots at five nitrogen rates, six plots per rate
(slightly spread horizontally), with the rate means and the grand mean.](nitrogen_data.pdf){width=66%}
:::

## The model in two forms

Index the groups by \( k=1,\dots,g \) and the observations within group \( k \) by \( j=1,\dots,n_k \). The
observations may be listed in any order in the vector \( \Y \); what matters is which group each belongs to.

::: {#def-aov-oneway}
[The one-way layout]

A **one-way layout** has \( g\ge2 \) groups with \( n_k\ge1 \) observations in group \( k \) and
\( n=\sum_kn_k>g \) observations in all. Its **cell-means form** is
\[
Y_{kj}=\mu_k+\varepsilon_{kj},\qquad\text{that is,}\qquad \Y=\Z\bmu+\be ,
\]
where \( \Z \) is the \( n\times g \) matrix whose column \( k \) is the indicator \( \bz_k \) of group \( k \) and
\( \bmu=(\mu_1,\dots,\mu_g)\T \). Its **effects form** is
\[
Y_{kj}=\mu+\alpha_k+\varepsilon_{kj},\qquad\text{that is,}\qquad \Y=[\bone,\Z]\begin{pmatrix}\mu\\\boldsymbol{\upalpha}\end{pmatrix}+\be .
\]
In both, the errors have \( \E(\be)=\bzero \) and \( \Cov(\be)=\sigma^2\I \). The **normal one-way model** adds
\( \be\sim\Normal_n(\bzero,\sigma^2\I) \). The layout is **balanced** if \( n_1=\dots=n_g=m \), say. The treatment
parameters \( \mu_k \) (or \( \alpha_k \)) are **fixed effects**: unknown constants attached to the \( g \) treatments
actually studied.
:::

The two forms describe the same mean vectors: since \( \bz_1+\dots+\bz_g=\bone \), the column space of \( [\bone,\Z] \) is \( \C(\Z) \), the
vectors constant within each group (@exm-proj-oneway-rank). The effects matrix has \( g+1 \) columns but rank \( g \), so its parameters are
not identified. Fitted values, residuals, sums of squares and tests depend only on \( \C(\Z) \), so the choice of form decides which
parameters one talks about, never the results of the analysis.

A word on notation. Group means are written \( \bar{y}_k \) (the literature also writes \( \bar{y}_{k\cdot} \)),
the grand mean is \( \bar{y}=n^{-1}\sum_{k,j}y_{kj} \), and \( \bar{\mu}=n^{-1}\sum_kn_k\mu_k \) is the *weighted*
average of the group means, which is \( \E(\bar{Y}) \). The unweighted average \( g^{-1}\sum_k\mu_k \) is a different
quantity unless the layout is balanced.

## The model space and its projection

The orthogonal projection onto \( \C(\Z) \) was computed in @exm-proj-oneway-M:
\[
\M=\sum_{k=1}^g\frac{\bz_k\bz_k\T}{n_k},\qquad (\M\y)_{kj}=\bar{y}_k ,
\]{#eq-aov-M}

so the fitted value of every observation is its group mean and the residual is its deviation from that
mean. Whatever the parameterization and whatever generalized inverse is used to solve the normal equations,
this is the fit (@thm-proj-M-formula). In the cell-means form the least squares estimate is unique,
\( \hat{\bmu}=(\bar{y}_1,\dots,\bar{y}_g)\T \), because \( \Z\T\Z=\diag(n_1,\dots,n_g) \) is invertible.

The residual sum of squares pools the variation inside the groups:
\[
\text{SSE}=\norm{(\I-\M)\y}^2=\sum_{k=1}^g\sum_{j=1}^{n_k}(y_{kj}-\bar{y}_k)^2=\sum_{k=1}^g(n_k-1)s_k^2 ,
\]{#eq-aov-sse}

where \( s_k^2 \) is the sample variance of group \( k \) (taken as zero when \( n_k=1 \)). The error space
\( \C(\Z)\perpc \) has dimension \( n-g \), so the unbiased estimator of \( \sigma^2 \) (@thm-lm-sigma2) is
\[
s^2=\frac{\text{SSE}}{n-g}=\sum_{k=1}^g\frac{n_k-1}{n-g}\,s_k^2 ,
\]
a weighted average of the group variances with weights proportional to their degrees of freedom. In the nitrogen trial
\( s=0.5211 \) t/ha on \( 25 \) degrees of freedom, and since the design is balanced \( s^2 \) is the plain average
of the five group variances.

```{.python .run #cell-nitrogen-trial-fit}
import numpy as np

rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)

Z = (group[:, None] == np.arange(g)).astype(float)  # cell-means model matrix
n_k = Z.sum(axis=0)
means = Z.T @ y / n_k                               # least squares: the group means
fitted = Z @ means                                  # M y
resid = y - fitted                                  # (I - M) y
s2 = resid @ resid / (n - g)                        # pooled within-group variance
print("group means:", means.round(3))
print(f"s = {np.sqrt(s2):.4f} on {n - g} degrees of freedom")
```

## What the parameters mean

The group means \( \mu_k \) are estimable, and so is every linear function of them. The effects form adds one
parameter that the data cannot determine. By @thm-est-oneway-contrasts, a function
\( \lambda_0\mu+\sum_k\lambda_k\alpha_k \) is estimable iff \( \lambda_0=\sum_k\lambda_k \); in particular no single
\( \alpha_k \) and not \( \mu \) itself is estimable, while every difference \( \alpha_k-\alpha_l \) is.

Software nevertheless reports values of \( \mu \) and \( \alpha_k \), by adding a side condition,
one nonestimable linear restriction that picks a single solution of the normal equations (@thm-est-side-conditions).
Under the restriction, each parameter becomes a particular estimable function of the
group means (@prp-est-side-meaning), and it is worth knowing which:

| Side condition | \( \mu \) becomes | \( \alpha_k \) becomes |
|:---|:---|:---|
| \( \sum_k\alpha_k=0 \) | \( g^{-1}\sum_l\mu_l \), the unweighted average | \( \mu_k-g^{-1}\sum_l\mu_l \) |
| \( \sum_kn_k\alpha_k=0 \) | \( \bar{\mu}=n^{-1}\sum_ln_l\mu_l \), the weighted average | \( \mu_k-\bar{\mu} \) |
| \( \alpha_1=0 \) | \( \mu_1 \), the reference group | \( \mu_k-\mu_1 \) |

Each row solves the \( g \) equations \( \mu+\alpha_k=\mu_k \) together with the side condition, and the estimates put \( \bar{y}_k \)
for \( \mu_k \) (@exr-aov-side-conditions). The second gives the classical \( \hat{\mu}=\bar{y} \), \( \hat{\alpha}_k=\bar{y}_k-\bar{y} \). The
first two agree in a balanced layout; in an unbalanced one they define different "overall means", and a reported \( \hat{\alpha}_k \) means
different things under each (@exm-est-five-conditions). Contrasts are unaffected by the choice.

::: {.warning}
A side condition is a convention for printing, not a statement about the treatments. "The effect of treatment 3 is
\( 0.4 \)" means nothing until the side condition is named. Report group means and contrasts.
:::

## Where the model comes from

In a randomized experiment much of the one-way model can be derived rather than assumed. Suppose each plot \( i \) has a yield \( u_i \) it would give without nitrogen, reflecting its
soil and position, and suppose rate \( k \) adds a fixed amount \( \tau_k \), the same on every plot: this is **unit–treatment
additivity**. Random assignment then makes the set of plots receiving rate \( k \) a simple random sample of the thirty.
The observed yields in group \( k \) are \( \tau_k \) plus a sample of the \( u_i \), so each has expectation
\( \tau_k+\bar{u} \), with \( \bar{u} \) the average of all thirty plot values. The expectations therefore have exactly the
cell-means structure, with \( \mu_k=\tau_k+\bar{u} \), whatever the pattern of fertility in the field. Had the six plots
given the highest rate been chosen for convenience, say all at one end of the field, a fertility gradient would be
confounded with the rate, and no analysis could separate the two.

The derived errors are not exactly those of @def-aov-oneway: sampling without replacement makes them slightly negatively
correlated (@exr-aov-randomization-mean), and nothing makes them normal. But under the hypothesis that no treatment has any
effect, every rearrangement of the yields among the treatments was equally likely, so the **randomization distribution** of
the \( F \) statistic of [Section 15.3](03-anova-table.html) can be computed exactly by permuting the data, and for moderate
designs the normal-theory \( F \) distribution approximates it closely. Chapter 23 develops permutation tests. Here the normal
model is the working model, and randomization is the reason for trusting it.

## Fixed and random effects

The five nitrogen rates were chosen deliberately, and the conclusions are about them: the treatment parameters are fixed
constants, the setting of this chapter. In other studies the levels are themselves a sample. A laboratory measuring an enzyme
activity in batches of reagent from a production line wants to know how much batches vary, not how batch 3 compares with
batch 5. The natural model then takes the effects \( \alpha_k \) to be random, uncorrelated with mean zero and variance \( \sigma_a^2 \), and independent of
the errors. The mean vector is \( \mu\bone \) and the covariance matrix is \( \sigma^2\I+\sigma_a^2\Z\Z\T \), so the model is no
longer a linear model with covariance \( \sigma^2\I \).

The projections and sums of squares are the same in both models, but their expectations differ
([Section 15.3](03-anova-table.html)), and so does the target: contrasts among the treatments studied, or the variance
components \( \sigma_a^2 \) and \( \sigma^2 \) (@exr-ss-random-oneway). Chapter 32 develops random and mixed effects.

## Exercises

### A. Check your understanding

::: {#exr-aov-small-layout}
[A1]

For a one-way layout with group sizes \( 2,1,3 \), write down \( \Z \), \( \M \) and \( \I-\M \) with the observations listed
group by group. What is the residual of the observation in group 2, and why? How many degrees of freedom does the
error space have, and which groups contribute to them?
:::

::: {#exr-aov-all-contrasts-zero}
[A2]

Show that \( \mu_1=\dots=\mu_g \) iff \( \sum_kc_k\mu_k=0 \) for every contrast \( \mathbf{c} \) (\( \sum_kc_k=0 \)), and iff
\( \mu_k-\mu_g=0 \) for \( k=1,\dots,g-1 \).
:::

::: {.solution}
If all means are equal to \( \mu \), then \( \sum_kc_k\mu_k=\mu\sum_kc_k=0 \). Conversely, the differences
\( \mu_k-\mu_g \) are contrasts, and if they vanish the means are equal. So the three conditions are equivalent. The
differences \( \mu_k-\mu_g \) span the space of contrasts, which has dimension \( g-1 \).
:::

### B. Practice

::: {#exr-aov-side-conditions}
[B1]

Derive the three rows of the table in "What the parameters mean". Then show that under \( \sum_kn_k\alpha_k=0 \) the least
squares estimates are \( \hat{\mu}=\bar{y} \) and \( \hat{\alpha}_k=\bar{y}_k-\bar{y} \), and that
\( \sum_kn_k\hat{\alpha}_k^2=\sum_kn_k(\bar{y}_k-\bar{y})^2 \).
:::

::: {.solution}
A solution of the normal equations (@thm-proj-normal-equations) reproduces the fitted values, so \( \hat{\mu}+\hat{\alpha}_k=\bar{y}_k \) for every \( k \).
With \( \sum_k\alpha_k=0 \), summing over \( k \) gives \( g\hat{\mu}=\sum_k\bar{y}_k \). With \( \sum_kn_k\alpha_k=0 \), multiplying by \( n_k \)
and summing gives \( n\hat{\mu}=\sum_kn_k\bar{y}_k=n\bar{y} \). With \( \alpha_1=0 \), the first equation gives \( \hat{\mu}=\bar{y}_1 \).
In each case \( \hat{\alpha}_k=\bar{y}_k-\hat{\mu} \). The same algebra with \( \mu_k \) in place of \( \bar{y}_k \) gives the parameter
column. The last identity is immediate from \( \hat{\alpha}_k=\bar{y}_k-\bar{y} \).
:::

::: {#exr-aov-two-groups}
[B2]

For \( g=2 \), show that the pooled variance of the one-way model is the pooled variance of the two-sample \( t \) test, and
that the \( t \) statistic for \( \mu_1-\mu_2 \) is
\( (\bar{y}_1-\bar{y}_2)/\bigl(s\sqrt{1/n_1+1/n_2}\bigr) \). ([Section 15.3](03-anova-table.html) shows that its square is the
one-way \( F \) statistic.)
:::

::: {#exr-aov-randomization-mean}
[B3]

In the setting of "Where the model comes from", suppose \( N \) units with values \( u_1,\dots,u_N \) are divided at random
into groups of fixed sizes \( n_1,\dots,n_g \), and unit \( i \) in group \( k \) yields \( u_i+\tau_k \). Show that
\( \E(\bar{Y}_k)=\tau_k+\bar{u} \) and that
\( \Cov(Y_{kj},Y_{k'j'})=-S_u^2/N \) for two different units, where
\( S_u^2=\sum_i(u_i-\bar{u})^2/(N-1) \). *Hint:* the sum of the \( N \) deviations \( u_i-\bar{u} \) is zero.
:::

::: {.solution}
Under random assignment the unit placed in position \( (k,j) \) is equally likely to be any of the \( N \), so
\( \E(Y_{kj})=\tau_k+\bar{u} \), and so is the average \( \E(\bar{Y}_k) \). For two different positions, the pair of units is
equally likely to be any ordered pair of distinct units, so
\( \Cov(Y_{kj},Y_{k'j'})=\E\{(u_I-\bar{u})(u_J-\bar{u})\} \) with \( (I,J) \) a random ordered pair of distinct indices. Since
\( \sum_i(u_i-\bar{u})=0 \),
\[
\sum_{i\ne l}(u_i-\bar{u})(u_l-\bar{u})=\Bigl(\sum_i(u_i-\bar{u})\Bigr)^2-\sum_i(u_i-\bar{u})^2=-(N-1)S_u^2 ,
\]
and dividing by the \( N(N-1) \) ordered pairs gives \( -S_u^2/N \). The variance of a single observation is
\( (N-1)S_u^2/N \), so the correlation between two observations is \( -1/(N-1) \), small when \( N \) is large.
:::

### C. Going deeper

::: {#exr-aov-averaging-group}
[C1]

Let \( \mathcal G \) be the set of \( n\times n \) permutation matrices that permute observations only within groups. Show that
\[
\M=\frac1{\lvert\mathcal G\rvert}\sum_{\boldsymbol{\Pi}\in\mathcal G}\boldsymbol{\Pi} ,
\]
so the one-way projection is the average over the symmetries of the design. Deduce that \( \C(\Z) \) is exactly the set of
vectors fixed by every \( \boldsymbol{\Pi}\in\mathcal G \).
:::

::: {.solution}
Take observations \( a \) and \( b \). If they are in different groups, no \( \boldsymbol{\Pi}\in\mathcal G \) sends \( b \) to \( a \), so entry
\( (a,b) \) of the average is \( 0=m_{ab} \). If both are in group \( k \), the permutations of group \( k \) that send \( b \) to \( a \)
are a fraction \( 1/n_k \) of all permutations of group \( k \), and the other groups are permuted freely, so entry \( (a,b) \) is
\( 1/n_k=m_{ab} \). A vector fixed by every \( \boldsymbol{\Pi} \) is constant within groups, because a transposition of two members of a group
exchanges their entries. Conversely, a vector constant within groups is fixed by each \( \boldsymbol{\Pi} \), hence by their average \( \M \), so
it lies in \( \C(\M)=\C(\Z) \).
:::

::: {#exr-aov-random-effects-mean}
[C2]

In the random-effects model of "Fixed and random effects", show that \( \Var(\bar{Y})=\sigma^2/n+\sigma_a^2\sum_kn_k^2/n^2 \).
For fixed \( n \) and \( g \), which group sizes make the grand mean most precise? Contrast this with the fixed-effects model, where
\( \Var(\bar{Y})=\sigma^2/n \) whatever the group sizes.
:::
