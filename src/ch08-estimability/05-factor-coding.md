# Indicator variables and the coding of a factor

Most regressors in practice are not measurements on a continuous scale. They record a
category: treatment arm, region, occupation, party, a yes-or-no answer. Such a variable is a
**factor**, and its possible values are its **levels**. A factor enters a linear model through
*indicator variables* (also called dummy variables), columns of zeros and ones that mark
membership of a level. There are many ways to build these columns, and they give quite
different coefficient tables. This section shows that they all describe the same fitted model.
Each choice of coding is a choice of which linear functions of the level means the coefficients
report. It gives a single formula that says exactly which ones.

## One binary regressor

Suppose each observation belongs to one of two groups, and let \( d_i=1 \) if observation \( i \) is in
group 1 and \( d_i=0 \) if it is in group 0. The model
\[
\E(y_i)=\beta_0+\beta_1d_i
\]
says that the mean is \( \mu_0=\beta_0 \) in group 0 and \( \mu_1=\beta_0+\beta_1 \) in group 1. So
\( \beta_1=\mu_1-\mu_0 \) is the difference between the two group means. The least squares fit sets
the fitted value in each group to that group's mean (the model space is the set of vectors constant
within groups, @exm-proj-oneway-M), so \( \hat{\beta}_0=\bar{y}_0 \) and
\( \hat{\beta}_1=\bar{y}_1-\bar{y}_0 \). The \( t \) statistic for \( \beta_1=0 \) is the pooled two-sample
\( t \) statistic (@exr-est-two-sample). The familiar two-sample comparison is a regression on one
indicator.

It is tempting to add a second indicator, \( 1-d_i \), for group 0. Then
\( \X=[\bone,\mathbf{d},\bone-\mathbf{d}] \) has three columns, but they satisfy
\( \bone-\mathbf{d}-(\bone-\mathbf{d})=\bzero \), so \( \rank(\X)=2 \). This is the **dummy variable trap**: an
intercept together with an indicator for every level produces an exactly collinear model
matrix. By [Section 8.1](01-identifiability.html), none of the three coefficients is then identifiable. The
trap is not a mistake in the model, which still describes two means, but in the coordinates.

## A factor with \( g \) levels

Let a factor have levels \( 1,\dots,g \), with \( n_k\ge1 \) observations at level \( k \). Its
**indicator matrix** is the \( n\times g \) matrix \( \Z=[\bz_1,\dots,\bz_g] \) with \( z_{ik}=1 \) if
observation \( i \) is at level \( k \) and \( 0 \) otherwise. Each row of \( \Z \) contains a single \( 1 \),
so
\[
\Z\bone_g=\bone_n,\qquad \Z\T\Z=\bD=\diag(n_1,\dots,n_g),\qquad \Z\T\y=(y_{1\cdot},\dots,y_{g\cdot})\T,
\]{#eq-est-indicator-facts}

where \( y_{k\cdot} \) is the total of the responses at level \( k \). The columns have disjoint supports, so
\( \Z \) has full column rank \( g \), and \( \C(\Z) \) is the set of vectors that are constant within levels.

The *cell-means model* \( \E(\Y)=\Z\bmu \), \( \bmu=(\mu_1,\dots,\mu_g)\T \), gives each level its own mean. It
has full rank, and by @eq-est-indicator-facts
\[
\hat{\bmu}=\bD^{-1}\Z\T\y=(\bar{y}_1,\dots,\bar{y}_g)\T ,
\]
the vector of level means. Any other model containing \( \C(\Z) \) and nothing more has the same fit.
The overparameterized model \( [\bone,\Z] \) of @exm-proj-oneway-rank is one example. The full-rank
codings below are others.

::: {#def-est-coding}
[Coding of a factor]

A **coding** of a factor with \( g \) levels is a nonsingular \( g\times g \) matrix \( \mathbf{K} \). The coded model
matrix is \( \Z\mathbf{K} \), and the model is \( \E(\Y)=\Z\mathbf{K}\bgamma \). When the first column of \( \mathbf{K} \) is \( \bone_g \),
so that \( \mathbf{K}=[\bone_g,\mathbf{C}] \), the coded matrix is \( [\bone_n,\Z\mathbf{C}] \): an intercept and the \( g-1 \) columns
\( \Z\mathbf{C} \). The \( g\times(g-1) \) matrix \( \mathbf{C} \) is then called the **coding matrix** of the factor.
:::

Row \( k \) of \( \mathbf{K} \) lists the values that the coded columns take for an observation at level \( k \). In the
software of [Section 8.4](04-side-conditions.html) the matrix \( \mathbf{C} \) is what one specifies (in R, the
“contrasts” attribute of a factor), and the program builds the columns \( \Z\mathbf{C} \) from it.

::: {#prp-est-coding}
[What a coding estimates]

Let \( \mathbf{K} \) be a coding and \( \hat{\bmu}=\bD^{-1}\Z\T\y \) the vector of level means. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\Z\mathbf{K})=\C(\Z) \), so every coding has the same fitted values, residuals and residual sum of squares;

2. the coefficients are \( \bgamma=\mathbf{K}^{-1}\bmu \), and their least squares estimates are
           \( \hat{\bgamma}=\mathbf{K}^{-1}\hat{\bmu} \);

3. \( [\bone_g,\mathbf{C}] \) is nonsingular iff the columns of \( \mathbf{C} \) are linearly independent and
           \( \bone_g\notin\C(\mathbf{C}) \);

4. to make the coefficients estimate prescribed linear functions \( \mathbf{L}\bmu \) of the level means, with \( \mathbf{L} \)
           a nonsingular \( g\times g \) matrix, use the coding \( \mathbf{K}=\mathbf{L}^{-1} \).
:::

:::

::: {.proof}
(a) The columns of \( \Z\mathbf{K} \) lie in \( \C(\Z) \), and \( \rank(\Z\mathbf{K})=\rank(\Z)=g \) because \( \mathbf{K} \) is nonsingular
(@prp-mat-rank-product(b)). A subspace of \( \C(\Z) \) of the same dimension is \( \C(\Z) \) itself. The fit depends
only on the column space (@thm-proj-reparam).
(b) The two models describe the same means iff \( \Z\mathbf{K}\bgamma=\Z\bmu \), and since \( \Z \) has full column rank this
means \( \mathbf{K}\bgamma=\bmu \). For the estimates, \( \Z\mathbf{K}\hat{\bgamma}=\M\y=\Z\hat{\bmu} \), and the same cancellation
gives \( \mathbf{K}\hat{\bgamma}=\hat{\bmu} \).
(c) A square matrix is nonsingular iff its columns are linearly independent. The columns \( \bone_g,\mathbf{c}_1,\dots,\mathbf{c}_{g-1} \) are
independent iff the \( \mathbf{c}_j \) are independent and \( \bone_g \) is not in their span.
(d) With \( \mathbf{K}=\mathbf{L}^{-1} \), (b) gives \( \bgamma=\mathbf{L}\bmu \).
:::

Part (b) is the whole story. Whatever the coding, the coefficients are fixed linear combinations of the
level means, and these combinations are the *rows of \( \mathbf{K}^{-1} \)*, not the columns of \( \mathbf{K} \). Each such
combination is estimable in the overparameterized model, since it is a function of the means. Part (d)
reverses the logic, and in practice it is the useful direction. Decide which questions the coefficients should
answer, write them as the rows of \( \mathbf{L} \), and invert.

## The standard codings

For a factor with \( g \) levels, the codings in common use are the following. For each we give the meaning of
the intercept \( \gamma_0 \) and of the other coefficients \( \gamma_1,\dots,\gamma_{g-1} \), that is, the rows of \( \mathbf{K}^{-1} \).

- **Reference (treatment) coding.** \( \mathbf{C}=\begin{psmallmatrix}\bzero\T\\\I_{g-1}\end{psmallmatrix} \): the columns are the
  indicators of levels \( 2,\dots,g \). Then \( \gamma_0=\mu_1 \), and \( \gamma_j=\mu_{j+1}-\mu_1 \) is the difference between level
  \( j+1 \) and the *reference level* 1. The level chosen as reference is arbitrary.
- **Sum-to-zero (deviation or effect) coding.** \( \mathbf{C}=\begin{psmallmatrix}\I_{g-1}\\-\bone\T\end{psmallmatrix} \). Level \( k<g \) has
  mean \( \gamma_0+\gamma_k \) and level \( g \) has mean \( \gamma_0-\sum_{j<g}\gamma_j \). Summing over the \( g \) levels gives
  \( \gamma_0=\bar{\mu}=g^{-1}\sum_k\mu_k \), the unweighted average of the level means. Then \( \gamma_k=\mu_k-\bar{\mu} \) is the
  deviation of level \( k \) from that average. The deviation of the last level is \( -\sum_{j<g}\gamma_j \).
- **Weighted effect coding.** \( \mathbf{C}=\begin{psmallmatrix}\I_{g-1}\\-(n_1,\dots,n_{g-1})/n_g\end{psmallmatrix} \). Now
  \( \sum_kn_k\mu_k=n\gamma_0 \), so \( \gamma_0=\sum_kn_k\mu_k/n \), the size-weighted average, estimated by the overall mean
  \( \bar{y} \). The coefficient \( \gamma_k=\mu_k-\gamma_0 \) is the deviation from it.
- **Successive differences.** With \( \mathbf{L} \) having first row \( g^{-1}\bone\T \) and rows \( \mathbf{e}_{j+1}\T-\mathbf{e}_j\T \), part (d) gives
  a coding in which \( \gamma_0=\bar{\mu} \) and \( \gamma_j=\mu_{j+1}-\mu_j \). This is natural for ordered levels.
- **Helmert and orthogonal polynomial codings**, whose columns are contrasts, are treated in
  [Section 8.7](07-contrasts.html).

The first three are the side conditions of @exm-est-five-conditions in another language. In the overparameterized
model \( \mu+\alpha_k \), reference coding is the side condition \( \alpha_1=0 \), sum-to-zero coding is \( \sum_k\alpha_k=0 \), and
weighted effect coding is \( \sum_kn_k\alpha_k=0 \). In each case the constrained \( (\mu,\alpha_2,\dots,\alpha_g) \) or
\( (\mu,\alpha_1,\dots,\alpha_{g-1}) \) equals \( \bgamma \) (@exr-est-coding-side).

::: {#exm-est-party-coding}
[Party identification and ideology]

An extract of the 1996 American National Election Study, distributed with statsmodels as a public-domain
data set, contains \( 944 \) respondents who placed themselves on a seven-point scale from
extremely liberal (1) to extremely conservative (7). It also records their party identification in seven levels from
strong Democrat to strong Republican. The level means and sizes are

| party identification | strong Dem | weak Dem | ind. Dem | independent | ind. Rep | weak Rep | strong Rep |
|---|---|---|---|---|---|---|---|
| \( n_k \) | 200 | 180 | 108 | 37 | 94 | 150 | 175 |
| \( \bar{y}_k \) | 3.335 | 3.617 | 3.685 | 4.054 | 4.957 | 5.033 | 5.691 |

Here are the coefficients of four codings of the same factor (\( \gamma_1,\dots,\gamma_6 \) refer to the columns of
\( \mathbf{C} \) in the order defined above; for the two deviation codings they belong to the first six levels):

| coding | \( \hat{\gamma}_0 \) | \( \hat{\gamma}_1 \) | \( \hat{\gamma}_2 \) | \( \hat{\gamma}_3 \) | \( \hat{\gamma}_4 \) | \( \hat{\gamma}_5 \) | \( \hat{\gamma}_6 \) |
|---|---|---|---|---|---|---|---|
| reference (strong Dem) | 3.335 | 0.282 | 0.350 | 0.719 | 1.622 | 1.698 | 2.356 |
| sum to zero | 4.339 | -1.004 | -0.722 | -0.654 | -0.285 | 0.618 | 0.694 |
| weighted effect | 4.325 | -0.990 | -0.709 | -0.640 | -0.271 | 0.632 | 0.708 |
| successive differences | 4.339 | 0.282 | 0.069 | 0.369 | 0.903 | 0.076 | 0.658 |

Each row is \( \mathbf{K}^{-1}\hat{\bmu} \) for its own \( \mathbf{K} \), as @prp-est-coding(b) says, and all four have residual sum of
squares \( 1178.3 \). The intercepts are the reference mean \( 3.335 \), the unweighted mean of the
seven means \( 4.339 \), the overall mean \( \bar{y}=4.325 \) and again \( 4.339 \).
[Figure 8.5.1](#fig-est-coding) shows the three reference points.
The successive differences show that the largest step in ideology, \( 0.903 \), falls between
pure independents and independents leaning Republican, and the step between weak and strong Republicans is
\( 0.658 \). The reference coding cannot show this directly. The two deviation codings nearly agree here
because the groups are of similar size, apart from the small group of pure independents.
:::

::: {when-format="html"}
![**Figure 8.5.1.** Mean self-placement by party identification in the 1996 ANES (point area proportional to
group size), with the three baselines used by the codings of @exm-est-party-coding: the reference level, the
unweighted mean of the seven means, and the overall mean.](coding_means.svg){#fig-est-coding width=70%}
:::

::: {when-format="pdf"}
![Mean self-placement by party identification in the 1996 ANES (point area proportional to
group size), with the three baselines used by the codings of @exm-est-party-coding: the reference level, the
unweighted mean of the seven means, and the overall mean.](coding_means.pdf){width=70%}
:::

```{.python .run #cell-factor-coding-data}
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].to_numpy()        # 1 = very liberal ... 7 = very conservative
level = anes["PID"].to_numpy().astype(int)    # 0 = strong Dem. ... 6 = strong Rep.
g = 7
Z = np.eye(g)[level]                          # n x 7 indicator (cell-means) matrix
n_k = Z.sum(axis=0)
cell_means = (Z.T @ y) / n_k
print("group sizes:", n_k.astype(int))
print("cell means: ", np.round(cell_means, 3))
```

```{.python .run #cell-factor-coding-codings}
def coding_matrices(n_k):
    """Each coding is X = Z K for a nonsingular g x g K; the coefficients are K^{-1} mu."""
    g = len(n_k)
    one = np.ones((g, 1))
    reference = np.vstack([np.zeros(g - 1), np.eye(g - 1)])            # level 0 is baseline
    deviation = np.vstack([np.eye(g - 1), -np.ones(g - 1)])            # sum-to-zero
    weighted = np.vstack([np.eye(g - 1), -n_k[:-1] / n_k[-1]])         # sum n_k * effect = 0
    successive = np.array([[(j + 1 - g if i <= j else j + 1) / g for j in range(g - 1)]
                           for i in range(g)])                         # backward differences
    return {"cell means": np.eye(g),
            "reference": np.hstack([one, reference]),
            "sum to zero": np.hstack([one, deviation]),
            "weighted effect": np.hstack([one, weighted]),
            "successive differences": np.hstack([one, successive])}

for name, K in coding_matrices(n_k).items():
    Xc = Z @ K
    gamma = np.linalg.lstsq(Xc, y, rcond=None)[0]
    print(f"{name:23s}", np.round(gamma, 3))
```

## Choosing a coding

Since every coding gives the same fit, the choice is about communication and, to a lesser extent, precision.
Four points are worth keeping in mind.

*Anything that depends on the fit is invariant.* This includes fitted values, residuals, \( R^2 \),
\( \hat{\sigma}^2 \), and the \( F \) test that all levels have the same mean. That test compares \( \C(\Z) \)
with \( \C(\bone) \), and the coding does not enter it (@thm-glh-f-test).

*Coefficient tests answer different questions.* The \( t \) test for a reference-coded coefficient tests
whether level \( k \) differs from the reference level. The \( t \) test for a sum-coded coefficient tests
whether level \( k \) differs from the average of all levels, a question that changes if a level is added. Neither
question is wrong, but they are different questions. A table of coefficients from an unfamiliar coding
should never be read as if it came from the familiar one.

*Precision depends on the reference.* In reference coding \( \Var(\hat{\gamma}_j)=\sigma^2(1/n_{j+1}+1/n_1) \). A
small reference group inflates the variance of every coefficient. In @exm-est-party-coding, comparing each
of the other five levels with the \( 37 \) pure independents instead of the \( 200 \) strong Democrats would
multiply the standard errors by factors between \( 1.55 \) and \( 1.76 \). A common default is to use the
largest group, or a natural control condition, as the reference.

*With other regressors, coefficients are adjusted comparisons.* When the model also contains other variables, a
reference-coded coefficient compares level \( k \) with the reference level at equal values of those variables. That is
the “holding the others fixed” reading of @def-lm-coefficient-interpretation, with all its caveats.
[Section 8.6](06-several-factors.html) gives an example in which adjustment changes the answer a great deal.

::: {.idea}
A coding is a choice of coordinates for the \( g \)-dimensional space of level means. The coefficients are
\( \mathbf{K}^{-1}\hat{\bmu} \). To get the comparisons you want, write them down first as the rows of a
nonsingular matrix \( \mathbf{L} \), then code with \( \mathbf{K}=\mathbf{L}^{-1} \).
:::

For an ordered factor, a useful instance of this principle is the *cumulative* coding, with columns
\( \mathbf{1}\{\text{level}\ge k\} \) for \( k=2,\dots,g \). Then \( \mathbf{K} \) is lower triangular with ones on and below the diagonal. Its
inverse takes successive differences, so the intercept is \( \mu_1 \) and the coefficient of the \( k \)th column is
\( \mu_k-\mu_{k-1} \) (@exr-est-cumulative). An observation's mean is built up as a sum of steps, which is
often how an ordered factor is best described.

## Exercises

### A. Check your understanding

::: {#exr-est-two-sample}
[A1]

In the model \( \E(y_i)=\beta_0+\beta_1d_i \), with \( n_0 \) observations in group 0 (level 1) and \( n_1 \) in group 1
(level 2), write down the coding matrix \( \mathbf{K} \) with \( \bmu=\mathbf{K}\bbeta \). Use @prp-est-coding(b), not an inverse
of \( \X\T\X \), to find \( \hbeta \). Then show that \( \Var(\hat{\beta}_1)=\sigma^2(1/n_0+1/n_1) \) and that
\( \hat{\beta}_1/\operatorname{se}(\hat{\beta}_1) \) is the pooled two-sample \( t \) statistic, where \( s^2=\text{SSE}/(n-2) \)
replaces \( \sigma^2 \).
:::

::: {.solution}
The two means are \( \mu_0=\beta_0 \) and \( \mu_1=\beta_0+\beta_1 \), so
\( \mathbf{K}=\begin{psmallmatrix}1&0\\1&1\end{psmallmatrix} \) and
\( \mathbf{K}^{-1}=\begin{psmallmatrix}1&0\\-1&1\end{psmallmatrix} \). By @prp-est-coding(b),
\( \hbeta=\mathbf{K}^{-1}\hat{\bmu}=(\bar{y}_0,\,\bar{y}_1-\bar{y}_0)\T \). The level means are the estimates of the full-rank
cell-means model, with \( \Cov(\hat{\bmu})=\sigma^2\bD^{-1}=\sigma^2\diag(1/n_0,1/n_1) \) (@thm-lm-moments), so
\( \Var(\hat{\beta}_1)=(-1,1)\,\sigma^2\bD^{-1}(-1,1)\T=\sigma^2(1/n_0+1/n_1) \). The residuals are deviations from the
group means (@prp-est-coding(a)), so \( s^2 \) is the pooled variance, and the ratio is
\( (\bar{y}_1-\bar{y}_0)/\{s\sqrt{1/n_0+1/n_1}\} \).
:::

::: {#exr-est-trap}
[A2]

For a factor with three levels, write down the null space of \( [\bone,\bz_1,\bz_2,\bz_3] \). Show that adding a
second factor's full set of indicators, \( [\bone,\Z_A,\Z_B] \), creates a null space of dimension at least two, and give
two independent null vectors.
:::

### B. Practice

::: {#exr-est-single-dummy}
[B1]

Let \( \X \) have full column rank, let \( \mathbf{e}_i \) be the \( i \)th coordinate vector with \( \mathbf{e}_i\notin\C(\X) \), and
fit \( [\X,\mathbf{e}_i] \) with coefficients \( (\bbeta,\delta) \). Write \( \hat{e}_i \) and \( h_{ii} \) for the residual and
the leverage (@def-proj-leverage) of observation \( i \) in the fit of \( \X \) alone.

::: {.enumerate options="label=(\alph*)"}
1. Use @thm-proj-fwl, with \( \mathbf{e}_i \) as the regressor of interest and \( \X \) as the one partialled out, to
   show that \( \hat{\delta}=\hat{e}_i/(1-h_{ii}) \), and explain why \( h_{ii}<1 \).

2. Show that the enlarged fit reproduces \( y_i \) exactly, and that its \( \hbeta \) is the estimate
   \( \hbeta_{(i)} \) without observation \( i \).

3. Deduce the deletion identity \( y_i-\x_{(i)}\T\hbeta_{(i)}=\hat{e}_i/(1-h_{ii}) \): the error in predicting
   \( y_i \) from the other observations is the ordinary residual inflated by \( 1/(1-h_{ii}) \).
:::

:::

::: {.solution}
(a) Apply @thm-proj-fwl(a) with \( \X_1=\mathbf{e}_i \) and \( \X_2=\X \), and let \( \M \) project onto \( \C(\X) \). Then
\( \hat{\delta} \) is the coefficient in the regression of \( \y \) on \( (\I-\M)\mathbf{e}_i \):
\( \hat{\delta}=\mathbf{e}_i\T(\I-\M)\y/\mathbf{e}_i\T(\I-\M)\mathbf{e}_i=\hat{e}_i/(1-h_{ii}) \). The denominator is
\( \norm{(\I-\M)\mathbf{e}_i}^2 \), which is positive because \( \mathbf{e}_i\notin\C(\X) \).
(b) The residual vector of the enlarged fit is orthogonal to every column, in particular to \( \mathbf{e}_i \), so its
\( i \)th entry is zero. For fixed \( \bbeta \), the residual sum of squares
\( \sum_{j\ne i}(y_j-\x_{(j)}\T\bbeta)^2+(y_i-\x_{(i)}\T\bbeta-\delta)^2 \) is minimized over \( \delta \) by making the last term
zero. What remains is the deletion criterion, and its minimizer \( \hbeta_{(i)} \) is unique: if \( \X\mathbf{b}=c\mathbf{e}_i \) with
\( \mathbf{b}\ne\bzero \), then \( c\ne0 \) by full rank, and \( \mathbf{e}_i\in\C(\X) \). So \( \hbeta=\hbeta_{(i)} \).
(c) By (b), \( \hat{\delta}=y_i-\x_{(i)}\T\hbeta_{(i)} \); compare with (a). The \( t \) statistic for \( \delta \) is the
externally studentized residual of Chapter 20.
:::

::: {#exr-est-cumulative}
[B2]

For an ordered factor with \( g \) levels, let the coded columns be the indicators of \( \{\text{level}\ge k\} \),
\( k=2,\dots,g \), together with an intercept. Show that \( \mathbf{K} \) is lower triangular with ones on and below the
diagonal, that \( \mathbf{K}^{-1} \) is the first-difference matrix, and hence that the coefficients are \( \mu_1 \) and
\( \mu_k-\mu_{k-1} \).
:::

::: {#exr-est-counts}
[B3]

On each of \( m \) days a warehouse runs \( N_{1i} \) machines of type 1 and \( N_{2i} \) of type 2, and meters its total
energy use \( V_i \). Suppose each running machine of type \( k \) uses an amount of energy with mean \( \theta_k \) and
variance \( \tau_k^2 \), independently of the other machines and the same on every day. Show that
\( \E(V_i)=\theta_1N_{1i}+\theta_2N_{2i} \), a linear model *without* intercept, and that \( (\theta_1,\theta_2) \) is
identifiable iff the vectors \( (N_{1i}) \) and \( (N_{2i}) \) are not proportional. Find \( \Cov(\mathbf{V}) \); why is it not a
multiple of \( \I \)? Why is it essential that \( \theta_k \) does not vary by day?
:::

### C. Going deeper

::: {#exr-est-coding-side}
[C1]

Show that fitting the reference-coded model is the same as solving the normal equations of the overparameterized
one-way model under the side condition \( \alpha_1=0 \), with \( \gamma_0=\mu \) and \( \gamma_j=\alpha_{j+1} \). State and prove the
corresponding facts for the sum-to-zero and weighted effect codings.
:::
