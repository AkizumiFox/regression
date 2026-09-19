# Analysis of covariance as a projection problem

In many experiments something about each unit is known before the treatments are applied: a
pretest score, the weight of an animal at the start of a feeding trial, the yield of a plot in the
previous season, a patient's blood pressure at enrolment. Such a measurement is a **covariate**.
If it predicts the response, much of the variation that the one-way analysis
([Chapter 15](../ch15-anova-subspaces/index.html)) calls error is not error at all. Part of it is the
covariate at work, and it can be removed.

The **analysis of covariance** (ANCOVA) adds covariates to a factor model.
[Section 8.6](../ch08-estimability/06-several-factors.html) treated one factor and one covariate (@prp-est-adjusted). This
section states the result for any design part, rank deficient or not, and any number of covariates, and then asks what a
covariate does in a randomized experiment.

## The model

Write the mean of the response as the sum of a *design part* and a *covariate part*,
\[
\E(\Y)=\X\bbeta+\Z\bgamma ,\qquad \Cov(\Y)=\sigma^2\I_n .
\]{#eq-dsn-ancova-model}

Here \( \X \) is \( n\times p \) with rank \( r \), built from the indicator columns of the treatment
factor and of any blocking factors. It is usually overparameterized. The \( n\times q \) matrix
\( \Z \) holds the covariate values, and \( \bgamma\in\Real^q \) their coefficients. Let \( \M \) be the
orthogonal projection onto \( \C(\X) \). Everything below is written with the three
matrices
\[
\mathbf{E}_{zz}=\Z\T(\I-\M)\Z,\qquad \mathbf{E}_{zy}=\Z\T(\I-\M)\y,\qquad E_{yy}=\y\T(\I-\M)\y .
\]{#eq-dsn-within}

The letter \( E \) stands for *error* or *within*: in the one-way layout \( (\I-\M)\y \) holds the deviations from the
group means (@exm-proj-oneway-M), so these are within-group sums of squares and products. We assume throughout that
\[
\rank[\X,\Z]=\rank(\X)+q,\quad\text{equivalently}\quad (\I-\M)\Z\ \text{has full column rank } q .
\]{#eq-dsn-ancova-rank}

The two forms are equivalent by @lem-ss-residualized, since
\( \C([\X,\Z])=\C(\X)\dirsum\C((\I-\M)\Z) \). The condition fails exactly when some combination of
the covariates lies in the design space. With one factor and one covariate, it fails when the
covariate is constant within every group. Then the covariate cannot be separated from the group
effects (@exr-est-ancova-identifiable).

## The main theorem

::: {#thm-dsn-ancova}
[Analysis of covariance]

Assume @eq-dsn-ancova-model and @eq-dsn-ancova-rank. Let \( \M_W \) be the orthogonal projection onto
\( \C([\X,\Z]) \).

::: {.enumerate options="label=(\alph*)"}
1. \( \M_W=\M+(\I-\M)\Z\,\mathbf{E}_{zz}^{-1}\Z\T(\I-\M) \), and the two terms are orthogonal projections
   with orthogonal ranges.

2. \( \bgamma \) is estimable, with the unique least squares estimate
   \( \hat{\bgamma}=\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy} \).

3. A function \( \blambda\T\bbeta \) is estimable in @eq-dsn-ancova-model iff it is estimable in the
   model \( \E(\Y)=\X\bbeta \) without covariates. If \( \blambda=\X\T\boldsymbol{\uprho} \), its least squares
   estimate is
   \[
   \blambda\T\hbeta=\boldsymbol{\uprho}\T\M(\y-\Z\hat{\bgamma}) ,
   \]{#eq-dsn-adjusted-estimate}

   so the fitted design part is \( \X\hbeta=\M(\y-\Z\hat{\bgamma}) \): the covariance-adjusted response
   \( \y-\Z\hat{\bgamma} \) analysed by the formulas of the model without covariates.

4. The residual sum of squares is
   \[
   \text{SSE}=E_{yy}-\mathbf{E}_{zy}\T\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy},
   \]{#eq-dsn-ancova-sse}

   on \( n-r-q \) degrees of freedom.

5. \( \Cov(\hat{\bgamma})=\sigma^2\mathbf{E}_{zz}^{-1} \), \( \Cov(\blambda\T\hbeta,\hat{\bgamma})=-\sigma^2\boldsymbol{\uprho}\T\M\Z\,\mathbf{E}_{zz}^{-1} \), and
   \[
   \Var(\blambda\T\hbeta)=\sigma^2\bigl(\boldsymbol{\uprho}\T\M\boldsymbol{\uprho}+\boldsymbol{\uprho}\T\M\Z\,\mathbf{E}_{zz}^{-1}\Z\T\M\boldsymbol{\uprho}\bigr).
   \]{#eq-dsn-ancova-var}

6. Let \( \C(\X_0)\subseteq\C(\X) \) with rank \( r_0<r \) and projection \( \Mo \), and define
   \( \mathbf{T}_{zz},\mathbf{T}_{zy},T_{yy} \) as in @eq-dsn-within with \( \Mo \) in place of \( \M \). Then @eq-dsn-ancova-rank holds for \( \X_0 \) as well, and the hypothesis that the design part of the mean lies in
   \( \C(\X_0) \), with \( \bgamma \) unrestricted, has the **adjusted sum of squares**
   \[
   \text{SS}_H=\bigl(T_{yy}-\mathbf{T}_{zy}\T\mathbf{T}_{zz}^{-1}\mathbf{T}_{zy}\bigr)-\bigl(E_{yy}-\mathbf{E}_{zy}\T\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy}\bigr)
   \]{#eq-dsn-adjusted-ss}

   on \( r-r_0 \) degrees of freedom. If \( \Y \) is normal,
   \[
   F=\frac{\text{SS}_H/(r-r_0)}{\text{SSE}/(n-r-q)}\sim F(r-r_0,\,n-r-q,\,\gamma),
   \]
   where the noncentrality \( \gamma \), a scalar not to be confused with the coefficient vector \( \bgamma \), is given by
   \( \sigma^2\gamma \) equal to the squared distance from \( \X\bbeta+\Z\bgamma \) to \( \C([\X_0,\Z]) \).
:::

:::

::: {.proof}
(a) By @lem-ss-residualized, \( \C([\X,\Z])=\C(\X)\dirsum\C(\tilde{\Z}) \) with \( \tilde{\Z}=(\I-\M)\Z \),
and the summands are orthogonal. \( \tilde{\Z} \) has full column rank by @eq-dsn-ancova-rank, so its
projection is \( \tilde{\Z}(\tilde{\Z}\T\tilde{\Z})^{-1}\tilde{\Z}\T \) (@thm-proj-M-formula), and
\( \tilde{\Z}\T\tilde{\Z}=\mathbf{E}_{zz} \) because \( \I-\M \) is symmetric and idempotent. @thm-proj-sum adds the
two projections.

(b) and (c). Put \( \hat{\bgamma}=\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy} \). By (a),
\[
\M_W\y=\M\y+(\I-\M)\Z\hat{\bgamma}=\M(\y-\Z\hat{\bgamma})+\Z\hat{\bgamma}.
\]
The first term lies in \( \C(\X) \), so it equals \( \X\bb \) for some \( \bb \), and \( (\bb,\hat{\bgamma}) \) is a
least squares estimate (@thm-proj-ls-projection). If \( (\bb',\bgamma') \) is another one, then
\( \X(\bb-\bb')=\Z(\bgamma'-\hat{\bgamma}) \). Applying \( \I-\M \) gives
\( \tilde{\Z}(\bgamma'-\hat{\bgamma})=\bzero \), and so \( \bgamma'=\hat{\bgamma} \). The coefficient vector of the
covariates is therefore the same in every least squares estimate, and so \( \bgamma \) is estimable (@thm-est-characterization(g)). The same argument shows that \( \X\hbeta=\M(\y-\Z\hat{\bgamma}) \) for
every least squares estimate.

For estimability of \( \blambda\T\bbeta \): in the covariance model this function has coefficient vector
\( (\blambda\T,\bzero\T) \), which is estimable iff \( \mathbf{a}\T[\X,\Z]=(\blambda\T,\bzero\T) \) for some
\( \mathbf{a} \) (@thm-est-estimable-identifiable). Such an \( \mathbf{a} \) gives \( \blambda\in\C(\X\T) \).
Conversely, if \( \blambda=\X\T\boldsymbol{\uprho} \), put
\[
\mathbf{a}=\M\boldsymbol{\uprho}-\tilde{\Z}\,\mathbf{E}_{zz}^{-1}\Z\T\M\boldsymbol{\uprho} .
\]
Then \( \mathbf{a}\T\X=\boldsymbol{\uprho}\T\M\X-\boldsymbol{\uprho}\T\M\Z\mathbf{E}_{zz}^{-1}\tilde{\Z}\T\X=\blambda\T \), because
\( \tilde{\Z}\T\X=\Z\T(\I-\M)\X=\bzero \), and
\( \mathbf{a}\T\Z=\boldsymbol{\uprho}\T\M\Z-\boldsymbol{\uprho}\T\M\Z\mathbf{E}_{zz}^{-1}\tilde{\Z}\T\Z=\bzero\T \), because
\( \tilde{\Z}\T\Z=\mathbf{E}_{zz} \). So \( \blambda\T\bbeta \) is estimable. Since \( \mathbf{a} \) lies in \( \C([\X,\Z]) \), its least squares
estimate is \( \mathbf{a}\T\M_W\y=\mathbf{a}\T\y \) (@prp-est-ls-estimator(a)):
\[
\mathbf{a}\T\y=\boldsymbol{\uprho}\T\M\y-\boldsymbol{\uprho}\T\M\Z\mathbf{E}_{zz}^{-1}\mathbf{E}_{zy}=\boldsymbol{\uprho}\T\M(\y-\Z\hat{\bgamma}).
\]

(d) \( \text{SSE}=\norm{(\I-\M_W)\y}^2=\y\T(\I-\M)\y-\y\T\tilde{\Z}\mathbf{E}_{zz}^{-1}\tilde{\Z}\T\y \) by (a), and
\( \tilde{\Z}\T\y=\mathbf{E}_{zy} \). The rank of \( \M_W \) is \( r+q \).

(e) \( \hat{\bgamma}=\mathbf{E}_{zz}^{-1}\tilde{\Z}\T\Y \) has covariance
\( \sigma^2\mathbf{E}_{zz}^{-1}\tilde{\Z}\T\tilde{\Z}\mathbf{E}_{zz}^{-1}=\sigma^2\mathbf{E}_{zz}^{-1} \). The estimate of
\( \blambda\T\bbeta \) is \( \mathbf{a}\T\Y \), and \( \M\boldsymbol{\uprho} \) is orthogonal to \( \tilde{\Z} \). So
\( \Var(\mathbf{a}\T\Y)=\sigma^2\norm{\mathbf{a}}^2 \) is the sum of the two squared lengths in @eq-dsn-ancova-var, and \( \Cov(\mathbf{a}\T\Y,\hat{\bgamma})=\sigma^2\mathbf{a}\T\tilde{\Z}\mathbf{E}_{zz}^{-1}=-\sigma^2\boldsymbol{\uprho}\T\M\Z\mathbf{E}_{zz}^{-1} \).

(f) If \( (\I-\Mo)\Z\mathbf{c}=\bzero \), then \( \Z\mathbf{c}\in\C(\X_0)\subseteq\C(\X) \), so \( (\I-\M)\Z\mathbf{c}=\bzero \) and
\( \mathbf{c}=\bzero \). So @eq-dsn-ancova-rank holds for \( \X_0 \), and (d) applied to \( [\X_0,\Z] \) gives the
first bracket in @eq-dsn-adjusted-ss as the residual sum of squares of the reduced model. The rank of the
reduced model is \( r_0+q \), and \( \C([\X_0,\Z])\subseteq\C([\X,\Z]) \). The distribution is @thm-glh-f-test.
:::

The proof is the Frisch–Waugh–Lovell theorem (@thm-proj-fwl) with the covariates as the residualized columns, stated
through \( \M \) so that the design part may be rank deficient.

::: {.idea}
Analysis of covariance is two ordinary analyses joined by one projection. The within-group regression of
\( \y \) on the covariates gives \( \hat{\bgamma} \), and the analysis of variance of the adjusted response
\( \y-\Z\hat{\bgamma} \) gives the treatment estimates. Their variances @eq-dsn-ancova-var carry an extra term for the
uncertainty in \( \hat{\bgamma} \). Tests must compare models, as in @eq-dsn-adjusted-ss: the between-groups sum of squares of
the adjusted response is not the numerator of any \( F \) statistic, because the reduced model must be allowed to
re-estimate \( \bgamma \).
:::

## The one-way layout with one covariate

For a single factor with \( g \) levels and a single covariate \( x \), the matrices of @eq-dsn-within are scalars:
\[
E_{xx}=\sum_{k}\sum_{i\in k}(x_i-\bar{x}_k)^2,\qquad
E_{xy}=\sum_{k}\sum_{i\in k}(x_i-\bar{x}_k)(y_i-\bar{y}_k),
\]
and \( E_{yy} \) is the within-groups sum of squares of \( y \). With \( \X_0=\bone \), the \( T \) quantities are the
same sums around the grand means. @thm-dsn-ancova gives:

::: {.enumerate options="label=(\roman*)"}
1. slope \( \hat{\beta}=E_{xy}/E_{xx} \);

2. \( \text{SSE}=E_{yy}-E_{xy}^2/E_{xx} \) on \( n-g-1 \) degrees of freedom;

3. adjusted treatment sum of squares
   \( \text{SS}_{\text{trt}}=\bigl(T_{yy}-T_{xy}^2/T_{xx}\bigr)-\bigl(E_{yy}-E_{xy}^2/E_{xx}\bigr) \) on \( g-1 \) degrees of freedom;

4. level means adjusted to a common covariate value, \( \bar{y}_k-\hat{\beta}(\bar{x}_k-x_0) \), which differ by
   \( (\bar{y}_k-\bar{y}_l)-\hat{\beta}(\bar{x}_k-\bar{x}_l) \), as in @prp-est-adjusted;

5. the test of \( \beta=0 \), with sum of squares \( E_{xy}^2/E_{xx} \) on one degree of freedom (@cor-ss-single-column).
:::

::: {#exm-dsn-tutoring}
[Four study formats with a pretest]

In a synthetic randomized experiment, \( n=40 \) students are assigned at random, ten to each of four
study formats for an introductory statistics unit. A pretest score \( x \) is recorded before the
assignment, and the response \( y \) is the score on the final examination. The data were generated from a
parallel lines model. Randomization balanced the pretest only on average: the four groups have pretest means \( 62.1 \), \( 57.2 \),
\( 59.9 \) and \( 60.5 \), and final means \( 72.86 \), \( 75.15 \), \( 69.09 \) and \( 78.09 \).
[Figure 18.1.1](#fig-dsn-ancova-lines) shows the data, and the sums of squares and products are:

| | \( xx \) | \( xy \) | \( yy \) | df |
|---|---|---|---|---|
| between formats | \( 124.9 \) | \( -31.4 \) | \( 432.94 \) | 3 |
| within formats | \( 3643.9 \) | \( 3159.2 \) | \( 4146.41 \) | 36 |
| total | \( 3768.8 \) | \( 3127.8 \) | \( 4579.35 \) | 39 |

Without the covariate, the one-way analysis has error mean square \( 4146.41/36=115.2 \) and treatment
\( F=1.25 \) on 3 and 36 degrees of freedom (\( p=0.305 \)). There is no evidence that the formats
differ.

With the covariate, the within-groups slope is \( \hat{\beta}=3159.2/3643.9=0.8670 \) points on the final per
point on the pretest. The residual sum of squares falls to
\[
\text{SSE}=4146.41-2738.94=1407.47
\]
on 35 degrees of freedom, so the error mean square falls from \( 115.2 \) to \( 40.2 \).
The residual standard deviation falls from \( 10.73 \) to \( 6.34 \). The reduced model, one line
for all forty students, has residual sum of squares \( 4579.35-2595.83=1983.52 \). The adjusted
treatment sum of squares is therefore \( 1983.52-1407.47=576.05 \), and
\[
F=\frac{576.05/3}{1407.47/35}=4.77
\]
on 3 and 35 degrees of freedom, with \( p=0.0068 \): clear differences between the formats. The test of
\( \beta=0 \) has \( F=68.1 \). The adjusted treatment sum of squares, \( 576.05 \), is even *larger* than the
unadjusted \( 432.94 \): format 2 had the weakest students on entry and still finished second, and adjustment gives it
credit for that.
:::

::: {when-format="html"}
![**Figure 18.1.1.** The tutoring experiment: final score against pretest score for the four formats, with the
fitted parallel lines of the covariance model. Open symbols mark the raw group means
\( (\bar{x}_k,\bar{y}_k) \), and filled symbols on the dotted vertical line mark the adjusted means at the overall pretest mean
\( \bar{x} \). The lines of formats 2 (dashed, drawn on top) and 4 almost coincide.](ancova_lines.svg){#fig-dsn-ancova-lines width=72%}
:::

::: {when-format="pdf"}
![The tutoring experiment: final score against pretest score for the four formats, with the
fitted parallel lines of the covariance model. Open symbols mark the raw group means
\( (\bar{x}_k,\bar{y}_k) \), and filled symbols on the dotted vertical line mark the adjusted means at the overall pretest mean
\( \bar{x} \). The lines of formats 2 (dashed, drawn on top) and 4 almost coincide.](ancova_lines.pdf){width=72%}
:::

```{.python .run #cell-tutoring-ancova}
import numpy as np
from scipy import stats

rng = np.random.default_rng(84)
g, m = 4, 10                                          # four formats, ten students each
n = g * m
group = rng.permutation(np.repeat(np.arange(g), m))   # the randomization
x = np.round(rng.normal(60, 10, n))                   # pretest, measured before assignment
mu_true = np.array([70.0, 74.0, 71.0, 77.0])
y = np.round(mu_true[group] + 0.8 * (x - 60) + rng.normal(0, 6, n), 1)

Z = np.eye(g)[group]                                  # indicator matrix of the formats
n_k = Z.sum(axis=0)
x_bar, y_bar = Z.T @ x / n_k, Z.T @ y / n_k           # group means

def within(v):
    """(I - M) v: deviations from the group means (M projects onto C(Z))."""
    return v - (Z.T @ v / n_k)[group]

E_xx, E_xy, E_yy = within(x) @ within(x), within(x) @ within(y), within(y) @ within(y)
beta = E_xy / E_xx                                    # slope from within-group variation only
sse = E_yy - E_xy ** 2 / E_xx                         # ANCOVA residual sum of squares
nu = n - g - 1

T_xx = np.sum((x - x.mean()) ** 2)                    # the same with I - M0, M0 = mean only
T_xy = np.sum((x - x.mean()) * (y - y.mean()))
T_yy = np.sum((y - y.mean()) ** 2)
sse0 = T_yy - T_xy ** 2 / T_xx                        # reduced model: one line for everyone
ss_trt = sse0 - sse                                   # adjusted treatment sum of squares
F_adj = (ss_trt / (g - 1)) / (sse / nu)
F_unadj = ((T_yy - E_yy) / (g - 1)) / (E_yy / (n - g))
print(f"slope {beta:.4f};  SSE: ANOVA {E_yy:.2f} on {n - g} df, ANCOVA {sse:.2f} on {nu} df")
print(f"treatments: unadjusted F = {F_unadj:.2f} (p = {stats.f.sf(F_unadj, g - 1, n - g):.3f}),"
      f" adjusted F = {F_adj:.2f} (p = {stats.f.sf(F_adj, g - 1, nu):.4f})")
```

## What a covariate does in a randomized experiment

In a randomized experiment the groups differ on the covariate only by chance. Adjustment still helps, for two reasons.

**Precision.**  The error variance of the covariance model is the variance of \( y \) *given* \( x \), which is smaller
whenever the two are correlated. The price is the estimated slope, which adds the second term in @eq-dsn-ancova-var. The
next result balances gain against price when the covariate is random, as it is in practice.

::: {#prp-dsn-precision}
[The precision of adjustment]

Consider a one-way layout with \( g \) groups of \( m \) units, \( N=gm \). Suppose the covariate is measured before
randomization, so that the \( x_i \) are independent \( \Normal(\xi,\tau^2) \) whatever their group, and that given all the
\( x \)'s the \( y_i \) are independent and normal, with \( \E(y_i\mid x_i)=\mu_k+\beta x_i \) for unit \( i \) in group \( k \)
and \( \Var(y_i\mid x_i)=\sigma^2 \). Let \( \rho \) be the
within-group correlation of \( x \) and \( y \). For two groups \( k\ne l \), the unadjusted estimate
\( \bar{y}_k-\bar{y}_l \) and the adjusted estimate \( (\bar{y}_k-\bar{y}_l)-\hat{\beta}(\bar{x}_k-\bar{x}_l) \) are both unbiased
for \( \mu_k-\mu_l \). If \( N-g>2 \),
\[
\frac{\Var(\text{adjusted})}{\Var(\text{unadjusted})}=(1-\rho^2)\Bigl(1+\frac{1}{N-g-2}\Bigr).
\]{#eq-dsn-precision}

:::

::: {.proof}
Write \( \bar{y}_k=\mu_k+\beta\bar{x}_k+\bar{\varepsilon}_k \), where \( \varepsilon_i=y_i-\mu_k-\beta x_i \) has mean zero and
variance \( \sigma^2 \) and is independent of all the \( x \)'s. Unconditionally \( \bar{y}_k-\bar{y}_l \) has mean
\( \mu_k-\mu_l \), because \( \E\bar{x}_k=\E\bar{x}_l \), and variance
\( (2/m)(\beta^2\tau^2+\sigma^2) \). The within-group variance of \( y \) is \( \beta^2\tau^2+\sigma^2 \), and
\( \sigma^2=(1-\rho^2)(\beta^2\tau^2+\sigma^2) \).

Given all the \( x \)'s, the covariance model holds with fixed covariates. By @thm-dsn-ancova(e), with
\( \boldsymbol{\uprho} \) the difference of the two scaled group indicators (so \( \M\boldsymbol{\uprho}=\boldsymbol{\uprho} \) and
\( \boldsymbol{\uprho}\T\Z=\bar{x}_k-\bar{x}_l \)), the adjusted estimate is conditionally unbiased with conditional variance
\( \sigma^2\{2/m+(\bar{x}_k-\bar{x}_l)^2/E_{xx}\} \). Its unconditional mean is therefore \( \mu_k-\mu_l \), and its
unconditional variance is the mean of the conditional variance. The group means of \( x \) are independent of
the within-group sum of squares \( E_{xx} \) (@cor-qf-sample-variance, applied group by group), with
\( \bar{x}_k-\bar{x}_l\sim\Normal(0,2\tau^2/m) \) and \( E_{xx}/\tau^2\sim\chi^2(N-g) \). For \( V\sim\chi^2(\nu) \) with
\( \nu>2 \), \( \E(1/V)=1/(\nu-2) \) (@exr-dsn-inverse-chisq). Hence
\[
\E\frac{(\bar{x}_k-\bar{x}_l)^2}{E_{xx}}=\frac{2\tau^2}{m}\cdot\frac{1}{\tau^2(N-g-2)},
\qquad
\Var(\text{adjusted})=\frac{2\sigma^2}{m}\Bigl(1+\frac1{N-g-2}\Bigr).
\]
Dividing by the unadjusted variance gives @eq-dsn-precision.
:::

The factor \( 1-\rho^2 \) is what adjustment would achieve if the slope were known. The factor
\( 1+1/(N-g-2) \) is the cost of estimating it. For the design of @exm-dsn-tutoring the cost is
\( 1.0294 \). Adjustment pays whenever \( \rho^2>1/(N-g-1) \), that is, when \( |\rho| \) exceeds
\( 0.169 \). At \( \rho=0.8 \), the correlation used to generate the tutoring data, the adjusted comparison has
\( 0.371 \) of the variance of the unadjusted one (simulation: \( 0.372 \);
[Figure 18.1.2](#fig-dsn-precision)). Without the covariate the experiment would need almost three times as many
students.

::: {when-format="html"}
![**Figure 18.1.2.** Variance of the adjusted difference of two group means relative to the unadjusted difference,
for four groups of ten, as a function of the within-group correlation between covariate and response. The
curve is @eq-dsn-precision, and the points are simulated from \( 100000 \) experiments each. The dashed curve is
the ratio when the slope is known.](precision_gain.svg){#fig-dsn-precision width=66%}
:::

::: {when-format="pdf"}
![Variance of the adjusted difference of two group means relative to the unadjusted difference,
for four groups of ten, as a function of the within-group correlation between covariate and response. The
curve is @eq-dsn-precision, and the points are simulated from \( 100000 \) experiments each. The dashed curve is
the ratio when the slope is known.](precision_gain.pdf){width=66%}
:::

**Chance imbalance.**  Randomization makes \( \bar{y}_k-\bar{y}_l \) unbiased *over repetitions of the experiment*.
In the experiment actually run, the proof above shows that, given the covariates,
\[
\E\bigl(\bar{y}_k-\bar{y}_l\mid x\text{'s}\bigr)=\mu_k-\mu_l+\beta(\bar{x}_k-\bar{x}_l),
\]
which is *not* \( \mu_k-\mu_l \) unless the covariate means agree. Once the pretests are known, the raw difference
has an error of predictable sign. In @exm-dsn-tutoring, format 1 started \( 4.9 \) points ahead of format 2 on the pretest, which predicts
\( 0.8670\times4.9\approx4.25 \) points of the final difference. The raw difference between formats 1 and 2
is \( -2.29 \), and the conditionally unbiased adjusted difference is \( -2.29-4.25=-6.54 \).

All this assumes the covariate was chosen in advance. Choosing, after seeing the data, the covariates that make the
effect look best is a form of multiple testing, which is why trial protocols name their covariates beforehand.

::: {.remark}
[Random covariates]

The theorem (@thm-dsn-ancova) treats \( \Z \) as fixed. If, given random covariates, the model holds with normal errors, every
\( F \) test and \( t \) interval of the theorem keeps its level unconditionally (@thm-cor-conditional, applied after
reparameterizing the design part to full rank, which changes neither the tests nor the estimable functions). Only the power
depends on the covariate values that occur.
:::

## Adjusting for a covariate that the treatment affects

::: {.warning}
A covariate must not be affected by the treatment. If it is measured after assignment and the treatment
changes it, adjusting for it removes the part of the treatment effect that acts through the covariate, and the
adjusted comparison no longer estimates the effect of assigning the treatment.
:::

If hours of study were used as a covariate, a format that works *by* getting students to study more would be
compared with the others at equal hours of study, which removes the very thing it does. In a synthetic randomized experiment with two
arms of \( 100 \), let the treatment raise a post-assignment variable \( w \) by \( 2 \) units and let
\( \E(y)=\theta T+0.75w \) with direct effect \( \theta=1 \). The effect of assigning the treatment is
\( 1+0.75\times2=2.50 \). Over \( 20000 \) replications, the unadjusted difference of means averages
\( 2.500 \), and the coefficient of the treatment after adjustment for \( w \) averages
\( 1.000 \), the direct effect, which is not what a decision about the treatment needs.

```{.python .run #cell-bad-control-mediator}
import numpy as np

rng = np.random.default_rng(2025)
n, reps = 200, 20_000
theta, delta, gamma = 1.0, 2.0, 0.75                  # direct effect; T -> w; w -> y
total = theta + gamma * delta                          # effect of assigning treatment
est_unadj, est_adj = np.empty(reps), np.empty(reps)
for r in range(reps):
    T = rng.permutation(np.repeat([0.0, 1.0], n // 2))    # randomized assignment
    w = rng.normal(size=n) + delta * T                    # measured after treatment
    y = theta * T + gamma * w + rng.normal(size=n)
    X = np.column_stack([np.ones(n), T, w])
    est_unadj[r] = y[T == 1].mean() - y[T == 0].mean()
    est_adj[r] = np.linalg.lstsq(X, y, rcond=None)[0][1]
print(f"total effect {total:.2f}: unadjusted mean {est_unadj.mean():.3f},"
      f" adjusted for w mean {est_adj.mean():.3f}")
```

[Chapter 25](../ch25-causal-interpretation/index.html) develops the language for this in observational data: mediators, on the causal path from treatment to
response, must not be adjusted for when the total effect is wanted. In a randomized experiment the safe rule is simple:
*adjust only for variables fixed before randomization.* By @prp-dsn-precision this costs little even when the covariate
turns out to be useless. A test of whether the groups are "balanced" on such a covariate is pointless, since its null
hypothesis is true by construction. What matters is whether the covariate predicts the response.

## Exercises

### A. Check your understanding

::: {#exr-dsn-constant-covariate}
[A1]

In the one-way layout with one covariate, suppose \( x \) is constant within each group. Show that @eq-dsn-ancova-rank fails, and describe the set of least squares estimates of \( \beta \). Is the adjusted
difference of two group means estimable?
:::

::: {.solution}
If \( x_i=c_k \) for every unit in group \( k \), then \( \x \) is a combination of the group indicators, so
\( \x\in\C(\X) \), \( (\I-\M)\x=\bzero \), and @eq-dsn-ancova-rank fails. The mean of group \( k \) is
\( \theta_k=\mu_k+\beta c_k \), and only the \( \theta_k \) are estimable. Every value of \( \beta \) belongs to some least
squares estimate: for any \( b \), the group parameters \( \bar{y}_k-bc_k \) together with \( b \) reproduce the fitted
group means. The adjusted difference \( \mu_k-\mu_l=(\theta_k-\theta_l)-\beta(c_k-c_l) \) is therefore estimable iff
\( c_k=c_l \). The data cannot separate a covariate that varies only between groups from the group effects.
:::

::: {#exr-dsn-slope-f}
[A2]

From the numbers in @exm-dsn-tutoring, compute the \( F \) statistic for \( \beta=0 \) and check that its square
root is the \( t \) statistic \( \hat{\beta}/\bigl(s/\sqrt{E_{xx}}\bigr) \).
:::

### B. Practice

::: {#exr-dsn-inverse-chisq}
[B1]

Let \( V\sim\chi^2(\nu) \) with \( \nu>2 \). Show that \( \E(1/V)=1/(\nu-2) \), and that \( \E(1/V)=\infty \) for
\( \nu\le2 \).
:::

::: {.solution}
The density of \( V \) is \( f_\nu(v)=v^{\nu/2-1}e^{-v/2}/\{2^{\nu/2}\Gamma(\nu/2)\} \). Then
\[
\E(1/V)=\int_0^\infty\frac{v^{\nu/2-2}e^{-v/2}}{2^{\nu/2}\Gamma(\nu/2)}\,dv
=\frac{2^{\nu/2-1}\Gamma(\nu/2-1)}{2^{\nu/2}\Gamma(\nu/2)}=\frac{1}{2(\nu/2-1)}=\frac1{\nu-2},
\]
using \( \Gamma(s+1)=s\Gamma(s) \). The integral converges at \( 0 \) iff \( \nu/2-2>-1 \), that is, iff \( \nu>2 \).
:::

::: {#exr-dsn-adjusted-ss-formula}
[B2]

In the one-way layout with one covariate, write \( B_{xx}=T_{xx}-E_{xx} \), and similarly \( B_{xy} \) and \( B_{yy} \), for
the between-groups sums of squares and products. Show that the adjusted treatment sum of squares equals
\[
B_{yy}-\Bigl(\frac{T_{xy}^2}{T_{xx}}-\frac{E_{xy}^2}{E_{xx}}\Bigr),
\]
and deduce that it exceeds \( B_{yy} \) iff \( E_{xy}^2/E_{xx}>T_{xy}^2/T_{xx} \). Show that this happens in @exm-dsn-tutoring, and explain it in terms of \( B_{xy} \).
:::

::: {.solution}
The formula is @eq-dsn-adjusted-ss with \( T_{yy}-E_{yy}=B_{yy} \). In the example
\( E_{xy}^2/E_{xx}=2738.94 \) and \( T_{xy}^2/T_{xx}=2595.83 \), so the adjusted sum of squares is larger, by
\( 143.11 \). The total regression uses both the within-group and the between-group variation of \( x \). Here the
between-group cross-product \( B_{xy}=-31.4 \) is negative: across the four formats, higher pretest means go with
*lower* final means, the opposite of the within-group relationship. The single line for all students therefore
explains less than the within-group lines do, and the difference is attributed to the formats.
:::

::: {#exr-dsn-general-var}
[B3]

In the one-way layout with \( q \) covariates, use @eq-dsn-ancova-var to show that the adjusted difference
\( (\bar{y}_k-\bar{y}_l)-(\bar{\bz}_k-\bar{\bz}_l)\T\hat{\bgamma} \) has variance
\( \sigma^2\{1/n_k+1/n_l+(\bar{\bz}_k-\bar{\bz}_l)\T\mathbf{E}_{zz}^{-1}(\bar{\bz}_k-\bar{\bz}_l)\} \).
:::

### C. Going deeper

::: {#exr-dsn-shift-invariance}
[C1]

In the one-way covariance model, show that replacing \( x \) by \( a+cx \) with \( c\ne0 \) leaves the fitted values,
SSE, the adjusted treatment sum of squares and the adjusted differences unchanged. Then replace \( x \) by
\( x+d_k \) on group \( k \), with group-specific constants \( d_k \). Show that \( \hat{\beta} \) and SSE are unchanged
but the adjusted differences and the adjusted treatment sum of squares change. Explain why a covariate
recorded with a group-specific calibration error is dangerous.
:::
