# Growth curves

A growth-curve model gives each subject a *curve of its own* — a straight line, a
polynomial in time — the population curve being their average. How is the average
curve estimated, and how is an individual's curve recovered from its own few
observations? The second answer is the best linear unbiased predictor: a
matrix-weighted average of the subject's own noisy curve and the precise
population curve, weighted by signal against signal-plus-noise — the shrinkage
of [Chapter 27](../ch27-shrinkage/index.html) by a different road.

## Two formulations

Let each of \( N \) subjects be observed at the same \( m \) times, and collect the
data in an \( m\times N \) matrix \( \Y \) whose column \( i \) is subject \( i \)'s
profile. Let \( \X \) be the \( m\times p \) matrix of a chosen basis in time, say
\( [\bone,\mathbf{t}] \) for a straight line.

The **Potthoff–Roy** formulation (Potthoff and Roy 1964) writes
\[
\E(\Y)=\X\,\boldsymbol{\Xi}\,\A\T,\qquad
\Cov(\text{column }i)=\bSigma,\ \text{columns independent},
\]{#eq-cls-potthoff-roy}

where \( \A \) is an \( N\times k \) matrix of between-subject covariates (a column
of ones if the subjects form one group, group indicators if not) and
\( \boldsymbol{\Xi} \) is the \( p\times k \) matrix of curve coefficients, with
\( \bSigma \) unstructured. It is a multivariate regression with a design on *both*
sides: \( \X \) acts on the times and \( \A \) on the subjects.

The **random-coefficient** formulation gives each subject its own coefficient
vector,
\[
\y_i=\X(\bbeta+\bu_i)+\be_i,\qquad
\bu_i\sim(\bzero,\bD)\ \text{independent of}\ \be_i\sim(\bzero,\sigma^2\I_m),
\]{#eq-cls-random-coefficient}

with \( \bu_1,\dots,\bu_N \) uncorrelated. This is the linear mixed
model (@def-mix-model) with \( \Z_i=\X \), \( \G=\bD \) and
\( \R=\sigma^2\I \), and the implied within-subject covariance is
\[
\bSigma=\X\bD\X\T+\sigma^2\I_m ,
\]{#eq-cls-rc-covariance}

a covariance with \( p(p+1)/2+1 \) parameters instead of the \( m(m+1)/2 \) of
@eq-cls-potthoff-roy — four, for a straight line, whatever \( m \) is. The two
describe the same *mean* structure: with one group of subjects, the
random-coefficient model is @eq-cls-potthoff-roy with \( \A=\bone_N \) and
\( \bSigma \) restricted to @eq-cls-rc-covariance. That restriction gives
\( \C(\bSigma\X)\subseteq\C(\X) \), so the population curve is ordinary least
squares and needs no covariance parameter (@prp-cls-growth(a) below); with
\( \bSigma \) unstructured the inclusion fails unless \( p=m \), and the efficient
estimate is a genuine generalized least squares fit.

## Estimation

::: {#prp-cls-growth}
[The balanced random-coefficient growth model]

Assume @eq-cls-random-coefficient with a common \( \X \) of full column rank
\( p\le m \), and \( \bD \) positive definite. Write \( \mathbf{K}=\X\T\X \),
\( \mathbf{C}=\sigma^2\mathbf{K}^{-1} \), \( \bS=\bD+\mathbf{C} \), and let
\( \bb_i=\mathbf{K}^{-1}\X\T\y_i \) be subject \( i \)'s own least squares coefficient vector.

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\bSigma\X)\subseteq\C(\X) \), so ordinary and generalized least squares
   coincide, and
   \[
\hbeta=\mathbf{K}^{-1}\X\T\bar{\y}=\bar{\bb},\qquad
\Cov(\hbeta)=\frac{\bD+\sigma^2\mathbf{K}^{-1}}{N}=\frac{\bS}{N},
\]{#eq-cls-growth-beta}

   exactly, for every \( N \). No covariance parameter is needed to estimate the
   population curve.

2. The best linear unbiased predictor of \( \bu_i \) (@thm-mix-blup, which is
   @thm-cor-blup applied coordinatewise with \( \x_0=\bzero \)) is
   \[
\hat{\bu}_i=\bD\X\T\bSigma^{-1}(\y_i-\X\hbeta)=\W(\bb_i-\hbeta),
\qquad \W=\bD\bS^{-1},
\]{#eq-cls-growth-blup}

   with \( \W=\bD(\bD+\mathbf{C})^{-1} \) and
   \( \Cov(\hat{\bu}_i-\bu_i)=\bD-(1-1/N)\W\bD \). The predicted curve of subject
   \( i \) is \( \X(\hbeta+\hat{\bu}_i) \), with
   \[
\Cov\bigl(\hbeta+\hat{\bu}_i-\bbeta-\bu_i\bigr)
 =\bigl(\bD^{-1}+\mathbf{C}^{-1}\bigr)^{-1}+(\I-\W)\frac{\mathbf{C}}{N}.
\]{#eq-cls-growth-pev}

3. The restricted maximum likelihood estimates of \( \sigma^2 \) and \( \bD \) are
   \[
\hat{\sigma}^2=\frac{1}{N(m-p)}\sum_{i=1}^N\norm{\y_i-\X\bb_i}^2,
\qquad
\hat{\bD}=\bS_b-\hat{\sigma}^2\mathbf{K}^{-1},
\]{#eq-cls-growth-reml}

   where \( \bS_b=(N-1)^{-1}\sum_i(\bb_i-\bar{\bb})(\bb_i-\bar{\bb})\T \), provided
   \( \hat{\bD} \) is nonnegative definite.
:::

:::

::: {.proof}
(a) \( \bSigma\X=\X\bD\mathbf{K}+\sigma^2\X=\X(\bD\mathbf{K}+\sigma^2\I) \), whose columns lie in
\( \C(\X) \). Stacking the subjects, the model matrix is \( \bone_N\otimes\X \) and the
covariance is \( \I_N\otimes\bSigma \), so
\( (\I_N\otimes\bSigma)(\bone_N\otimes\X)=\bone_N\otimes\bSigma\X \) has columns in
\( \C(\bone_N\otimes\X) \), and @thm-proj-kruskal gives the equality of the two fits
(this is @thm-ggm-ols-blue in the balanced case). The ordinary least squares
estimate from the stacked data is
\( (N\mathbf{K})^{-1}\sum_i\X\T\y_i=\mathbf{K}^{-1}\X\T\bar{\y}=\bar{\bb} \). Its covariance is
\( \mathbf{K}^{-1}\X\T(\bSigma/N)\X\mathbf{K}^{-1} \), and
\( \mathbf{K}^{-1}\X\T\bSigma\X\mathbf{K}^{-1}=\mathbf{K}^{-1}(\mathbf{K}\bD\mathbf{K}+\sigma^2\mathbf{K})\mathbf{K}^{-1}
 =\bD+\sigma^2\mathbf{K}^{-1} \).

(b) The first expression is @thm-mix-blup. For the second, @thm-mat-woodbury gives
\( \bSigma^{-1}=\sigma^{-2}\{\I-\X(\sigma^2\bD^{-1}+\mathbf{K})^{-1}\X\T\} \), so
\( \bD\X\T\bSigma^{-1}
 =\sigma^{-2}\bD\{\I-\mathbf{K}(\sigma^2\bD^{-1}+\mathbf{K})^{-1}\}\X\T \). The braced
matrix is \( \sigma^2\bD^{-1}(\sigma^2\bD^{-1}+\mathbf{K})^{-1} \), so
\( \bD\X\T\bSigma^{-1}=(\sigma^2\bD^{-1}+\mathbf{K})^{-1}\X\T \); applying this to
\( \y_i-\X\hbeta \) and inserting \( \mathbf{K}\mathbf{K}^{-1} \) gives
\( \hat{\bu}_i=(\sigma^2\bD^{-1}+\mathbf{K})^{-1}\mathbf{K}(\bb_i-\hbeta) \), and
\( (\sigma^2\bD^{-1}+\mathbf{K})^{-1}\mathbf{K}=(\sigma^2\I+\bD\mathbf{K})^{-1}\bD\mathbf{K}
 =\bD(\sigma^2\mathbf{K}^{-1}+\bD)^{-1}=\W \).

For the prediction error covariance, @thm-mix-blup gives
\( \Cov(\hat{\bu}-\bu)=\G-\G\Z\T\bP\Z\G \) with
\( \bP=\V^{-1}-\V^{-1}\X_*(\X_*\T\V^{-1}\X_*)^{-1}\X_*\T\V^{-1} \) for the stacked
model. Here \( \Z=\I_N\otimes\X \), \( \G=\I_N\otimes\bD \),
\( \X_*=\bone_N\otimes\X \) and \( \V=\I_N\otimes\bSigma \), so
\( \Z\T\bP\Z=(\I_N-N^{-1}\bone_N\bone_N\T)\otimes(\X\T\bSigma^{-1}\X) \) and the
\( i \)th diagonal block of \( \G-\G\Z\T\bP\Z\G \) is
\( \bD-(1-1/N)\bD\X\T\bSigma^{-1}\X\bD=\bD-(1-1/N)\W\bD \), using
\( \bD\X\T\bSigma^{-1}\X=\W \) from the first part.

For @eq-cls-growth-pev, write \( \bb_i=\bbeta+\bu_i+\mathbf{c}_i \) with
\( \mathbf{c}_i=\mathbf{K}^{-1}\X\T\be_i\sim(\bzero,\mathbf{C}) \) independent of \( \bu_i \), so
\( \hbeta=\bbeta+\bar{\bu}+\bar{\mathbf{c}} \) and
\( \hat{\bu}_i=\W(\bb_i-\hbeta) \). Then
\( \Cov(\hbeta,\hat{\bu}_i)=\{\Cov(\bar{\bb},\bb_i)-\Cov(\bar{\bb})\}\W\T=\bzero \) and
\( \Cov(\hbeta,\bu_i)=\bD/N \), so
\[
\Cov(\hbeta+\hat{\bu}_i-\bbeta-\bu_i)
=\frac{\bS}{N}+\bD-\Bigl(1-\frac1N\Bigr)\W\bD-\frac{2\bD}{N}.
\]
Now \( \bD-\W\bD=\bD-\bD\bS^{-1}\bD=\bD\bS^{-1}(\bS-\bD)=\W\mathbf{C} \) and
\( \I-\W=(\bS-\bD)\bS^{-1}=\mathbf{C}\bS^{-1} \), so the right-hand side is
\( \W\mathbf{C}+(\bS-2\bD+\bD\bS^{-1}\bD)/N
 =\W\mathbf{C}+\mathbf{C}\bS^{-1}\mathbf{C}/N=\W\mathbf{C}+(\I-\W)\mathbf{C}/N \).
Finally \( \W\mathbf{C}=\bD\bS^{-1}\mathbf{C}
 =\{\mathbf{C}^{-1}\bS\bD^{-1}\}^{-1}=(\mathbf{C}^{-1}+\bD^{-1})^{-1} \), which
is @eq-cls-growth-pev.

(c) An error contrast is a linear function of the data with mean zero for every
\( \bbeta \) (@def-mix-reml). Let \( \bU \) be \( m\times(m-p) \) with columns an
orthonormal basis of \( \C(\X)\perpc \), and consider the \( Nm-p \) functions
\[
\bw_i=\bU\T\y_i\ (i=1,\dots,N),\qquad \bb_i-\bar{\bb}\ (i=1,\dots,N).
\]
They are error contrasts, and a maximal set: the map \( \y_i\mapsto(\bb_i,\bw_i) \)
is a bijection, the \( \bw_i \) contribute \( N(m-p) \) independent functions and the
\( \bb_i-\bar{\bb} \) contribute \( (N-1)p \), totalling \( Nm-p \). Now
\( \bw_i=\bU\T\be_i\sim\Normal_{m-p}(\bzero,\sigma^2\I) \) independently, while
\( \bb_i\sim\Normal_p(\bbeta,\bS) \) independently, and
\( \Cov(\bb_i,\bw_i)=\sigma^2\mathbf{K}^{-1}\X\T\bU=\bzero \); the two blocks are
independent, so the restricted likelihood factors as
\[
L_R(\sigma^2,\bD)=L_1(\sigma^2)\cdot L_2(\bS),\qquad \bS=\bD+\sigma^2\mathbf{K}^{-1},
\]
with \( L_1 \) the likelihood of \( N(m-p) \) independent \( \Normal(0,\sigma^2) \)
variables and \( L_2 \) that of \( \bS \) from the contrasts \( \bb_i-\bar{\bb} \). Since
\( (\sigma^2,\bD)\mapsto(\sigma^2,\bS) \) is a bijection onto
\( \{\bS-\sigma^2\mathbf{K}^{-1}\succeq\bzero\} \), the unconstrained maximizers are
the restricted maximum likelihood estimates whenever they satisfy the constraint.
Maximizing \( L_1 \) gives
\( \hat{\sigma}^2=\sum_i\norm{\y_i-\X\bb_i}^2/\{N(m-p)\} \), and maximizing
\( L_2 \), the restricted likelihood of a multivariate normal mean model with design
\( \bone_N \), gives \( \hat{\bS}=\bS_b \) (@thm-mix-reml). Subtracting yields
\( \hat{\bD} \).
:::

Part (a) is a statement about the restriction @eq-cls-rc-covariance, not about
balance alone: the population curve is ordinary least squares on the *average
profile*, which is also the average of the subjects' own curves, and nothing has
to be estimated first. Potthoff and Roy's own proposal was not to estimate
\( \bSigma \) at all: transforming each profile by a positive definite weight
matrix reduces @eq-cls-potthoff-roy to an ordinary multivariate analysis of
variance, whose tests are exact for weights chosen without looking at the data.
The efficient estimate came later, with Rao (1965); estimating \( \bSigma \) from
the residual profiles and weighting by it is feasible generalized least
squares (@def-ggm-feasible).

::: {.idea}
\( \W=\bD(\bD+\mathbf{C})^{-1} \) compares two covariances: \( \bD \), how much
subjects really differ, and \( \mathbf{C}=\sigma^2(\X\T\X)^{-1} \), how noisily each
subject's own curve is measured. When subjects differ a lot relative to the noise
each keeps its own curve, \( \W\approx\I \); when they hardly differ every subject
is given the population curve, \( \W\approx\bzero \). Between those the predictor
borrows strength across subjects, coordinate by coordinate: an intercept may be
barely shrunk while a slope is shrunk by a quarter.
:::

Part (c) says that restricted maximum likelihood here is the obvious moment
estimator: pool the within-subject residual sums of squares for \( \sigma^2 \),
and from the spread of the subjects' own coefficient vectors subtract the part
that is estimation noise. The subtraction can produce a matrix that is not
nonnegative definite — the growth-curve version of the negative variance component
of @prp-mix-oneway — and the estimate then sits on the boundary, where the closed
form no longer applies.

## Whose curve

Two prediction problems must be kept apart. For a subject **in the data** the
predicted curve is \( \X(\hbeta+\hat{\bu}_i) \), with error covariance
@eq-cls-growth-pev: the first term is what remains even with \( \bbeta \) known,
the second the price of estimating it. For a subject **not in the data** nothing is
known but its membership of the population, so the best predictor is
\( \X\hbeta \) with error covariance \( \bD+\bS/N \) — far larger, as it should
be.

## An example

::: {#exm-cls-grunfeld-growth}
[The investment trajectories of eleven firms]

The panel of @exm-cls-grunfeld-covariance, seen a third time — after the random
firm effect of [Chapter 32](../ch32-linear-mixed-models/index.html) and the
covariance models of [Section 33.4](04-covariance-models.html) — now with a
straight line in \( t \) (decades from the middle of the record) whose intercept
*and slope* vary by firm. The within-firm residual standard deviation is
\( \hat{\sigma}=0.2797 \), and
\[
\hat{\bD}=\begin{pmatrix}2.3635&0.0455\\0.0455&0.0361\end{pmatrix},
\]
so firms' log-investment levels have standard deviation \( 1.5374 \) —
enormous, the difference between General Motors and Diamond Match — and their
slopes \( 0.1899 \) per decade, correlated \( 0.156 \) with
the level. Maximizing the restricted likelihood numerically returns the same three
values to four decimals, as @prp-cls-growth(c) predicts.

The population curve has intercept \( 3.9109 \) and slope
\( 0.5577 \) per decade, with standard errors \( 0.4639 \) and
\( 0.0659 \) from @eq-cls-growth-beta; the slope is \( 8.46 \)
standard errors from zero. The shrinkage matrix is
\[
\W=\begin{pmatrix}0.9983&0.0016\\0.0048&0.7496\end{pmatrix},
\]
so intercepts are left alone — with twenty observations per firm a firm's own
level is known almost exactly — while slopes are pulled a quarter of the way
towards the population slope. The firms' own slopes span \( 0.158 \) to
\( 0.943 \), the predicted slopes \( 0.248 \) to
\( 0.846 \), a standard deviation \( 0.755 \) times that of the
raw slopes; the largest adjustment is IBM's, from \( 0.943 \) to
\( 0.846 \). [Figure 33.5.1](#fig-cls-growth) shows the eleven firms
with all three lines.

For the last year of the record, \( t=0.95 \), the predicted value on a firm's own
curve has standard error \( 0.1102 \), against \( 1.6461 \)
for a firm outside the data — fifteen times larger, and simulation from the fitted
model confirms both, at \( 0.1103 \) and \( 1.6603 \). The
interval for a new firm, \( 1.214 \) to \( 7.667 \) on the log
scale, is nearly useless, and rightly so: without knowing a firm's size, its
investment cannot be predicted within three orders of magnitude.
:::

```{.python .run #cell-growth-curves-data}
import numpy as np
import statsmodels.api as sm

df = sm.datasets.grunfeld.load_pandas().data
firms = sorted(df["firm"].unique())
years = np.sort(df["year"].unique())
N, m = len(firms), len(years)
t = (years - years.mean()) / 10                      # decades from the middle of the record
X = np.column_stack([np.ones(m), t])                 # the same design for every firm
p = X.shape[1]
Y = np.column_stack([np.log(df[df["firm"] == f].sort_values("year")["invest"].to_numpy())
                     for f in firms])                # m x N: one column per firm
```

```{.python .run #cell-growth-curves-twostage}
XtXinv = np.linalg.inv(X.T @ X)
B = XtXinv @ X.T @ Y                                 # p x N: each firm's own least squares line
beta_hat = B.mean(axis=1)                            # the average of the firm lines
E = Y - X @ B                                        # within-firm residuals
sigma2 = np.sum(E**2) / (N * (m - p))                # pooled within-firm variance
Sb = np.cov(B, ddof=1)                               # spread of the firm lines
D_hat = Sb - sigma2 * XtXinv                         # ... minus what estimation noise adds
print("beta_hat", beta_hat.round(4), " sigma", round(np.sqrt(sigma2), 4))
print("D_hat\n", D_hat.round(4))
print("implied sd of firm intercepts and slopes:", np.sqrt(np.diag(D_hat)).round(4))
```

```{.python .run #cell-growth-curves-blup}
Sigma = X @ D_hat @ X.T + sigma2 * np.eye(m)
cov_beta = (D_hat + sigma2 * XtXinv) / N             # exact covariance of beta_hat
W = D_hat @ np.linalg.inv(D_hat + sigma2 * XtXinv)   # the shrinkage matrix
B_blup = beta_hat[:, None] + W @ (B - beta_hat[:, None])
u_direct = D_hat @ X.T @ np.linalg.solve(Sigma, Y - (X @ beta_hat)[:, None])
print("shrinkage matrix\n", W.round(4))
print("largest difference from the direct formula:",
      np.abs(B_blup - beta_hat[:, None] - u_direct).max())
```

```{.python .run #cell-growth-curves-prediction}
C = sigma2 * XtXinv                                  # noise in a firm's own line
V_known = np.linalg.inv(np.linalg.inv(D_hat) + np.linalg.inv(C)) \
    + (np.eye(p) - W) @ C / N                        # error of the predicted firm curve
V_new = D_hat + cov_beta                             # error for a firm not in the data
x_end = np.array([1.0, t[-1]])
se_curve_known = np.sqrt(x_end @ V_known @ x_end)
se_curve_new = np.sqrt(x_end @ V_new @ x_end)
pred_new = x_end @ beta_hat
print(f"1954 curve: se for a firm in the data {se_curve_known:.4f}, "
      f"for a new firm {se_curve_new:.4f}")
```

::: {when-format="html"}
![**Figure 33.5.1.** Log gross investment of eleven firms, 1935–1954, with each
firm's own least squares line (dotted), its predicted curve
\( \X(\hbeta+\hat{\bu}_i) \) (solid) and the population curve \( \X\hbeta \) (dashed).
The last panel plots own slope against predicted slope.](growth_curves.svg){#fig-cls-growth width=100%}
:::

::: {when-format="pdf"}
![Log gross investment of eleven firms, 1935–1954, with each
firm's own least squares line (dotted), its predicted curve
\( \X(\hbeta+\hat{\bu}_i) \) (solid) and the population curve \( \X\hbeta \) (dashed).
The last panel plots own slope against predicted slope.](growth_curves.pdf){width=100%}
:::

::: {.remark}
[Beyond the balanced straight line]

Three extensions matter in practice and none changes the ideas. Subjects observed
at different times, or a different number of times: @eq-cls-growth-beta fails, the
fit becomes genuine generalized least squares with
\( \bSigma_i=\X_i\bD\X_i\T+\sigma^2\I \), and Henderson's
equations (@thm-mix-henderson) do the computing. Between-subject covariates:
replace \( \bbeta \) by \( \boldsymbol{\Xi}\mathbf{a}_i \) as in the
formulation @eq-cls-potthoff-roy. Serial correlation on top of the random coefficients: take
\( \R \) autoregressive, as in
[Section 33.4](04-covariance-models.html). Beyond all of these lies the question
whether a polynomial in time is the right shape at all: [Chapter 42](../ch42-polynomials-piecewise/index.html)
and [Chapter 43](../ch43-smoothing/index.html).
:::

## Exercises

### A. Check your understanding

::: {#exr-cls-shrinkage-scalar}
[A1]

Take \( p=1 \), so \( \X=\bone_m \) and each subject has a single random level. Show
that \( \W \) reduces to \( \sigma_u^2/(\sigma_u^2+\sigma^2/m) \) and that the
predictor is the shrunken mean \( \hbeta+\W(\bar{y}_i-\hbeta) \). What happens as
\( m\to\infty \), and as \( \sigma_u^2\to0 \)?
:::

::: {.solution}
With \( \X=\bone_m \), \( \mathbf{K}=m \), \( \mathbf{C}=\sigma^2/m \) and
\( \bD=\sigma_u^2 \), so \( \W=\sigma_u^2/(\sigma_u^2+\sigma^2/m) \) and
\( \bb_i=\bar{y}_i \). As \( m\to\infty \) the noise in \( \bar{y}_i \) vanishes,
\( \W\to1 \) and no shrinkage occurs; as \( \sigma_u^2\to0 \), \( \W\to0 \) and every
subject is given the common mean, correctly, since the subjects are then
identical.
:::

::: {#exr-cls-blup-not-unbiased}
[A2]

The predicted curves are less spread out than the true ones. Is this a bias?
Compute \( \E(\hat{\bu}_i) \) and \( \Cov(\hat{\bu}_i) \) and compare the latter with
\( \bD \).
:::

### B. Practice

::: {#exr-cls-growth-quadratic}
[B1]

Suppose the fixed effects use a quadratic in time,
\( \X=[\bone,\mathbf{t},\mathbf{t}^2] \), but the random effects only a line,
\( \Z=[\bone,\mathbf{t}] \). Show that \( \C(\bSigma\X)\subseteq\C(\X) \) still holds,
and identify what changes in \( \Cov(\hbeta) \). What if instead the random effects
used the quadratic and the fixed effects only the line?
:::

::: {#exr-cls-growth-unbalanced}
[B2]

Let subject \( i \) be observed at its own times, giving \( \X_i \) of full column
rank \( p \), and let \( \bb_i=(\X_i\T\X_i)^{-1}\X_i\T\y_i \). Show that
\( \bb_i\sim(\bbeta,\bD+\sigma^2(\X_i\T\X_i)^{-1}) \) independently over \( i \), and
deduce that the generalized least squares estimate from the \( \bb_i \) is
\( \hbeta=(\sum_i\W_i)^{-1}\sum_i\W_i\bb_i \) with
\( \W_i=\{\bD+\sigma^2(\X_i\T\X_i)^{-1}\}^{-1} \), a precision-weighted average of
the subjects' own curves. Recover @eq-cls-growth-beta as the balanced case.
:::

::: {.solution}
\( \bb_i=\bbeta+\bu_i+(\X_i\T\X_i)^{-1}\X_i\T\be_i \) is a sum of independent terms
with covariances \( \bD \) and \( \sigma^2(\X_i\T\X_i)^{-1} \), and different subjects
are independent. The \( \bb_i \) are therefore independent with common mean
\( \bbeta \) and covariances \( \W_i^{-1} \), so @cor-opt-aitken applied to the stacked
\( \bb_i \) with design \( \bone_N\otimes\I_p \) gives the weighted average. When every
\( \X_i=\X \) the weights cancel, leaving \( \hbeta=\bar{\bb} \) with covariance
\( \W_1^{-1}/N \), which is @eq-cls-growth-beta.
:::

### C. Going deeper

::: {#exr-cls-growth-negative-d}
[C1]

Simulate data for which \( \hat{\bD} \) of @eq-cls-growth-reml is not nonnegative
definite. Describe what the restricted maximum likelihood estimate does instead,
relate it to the boundary problem of @prp-mix-inference, and explain why the
two-stage estimate is nevertheless unbiased for \( \bD \).
:::
