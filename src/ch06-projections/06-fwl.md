# The Frisch–Waugh–Lovell theorem

A coefficient in a multiple regression is often described as the effect of its
regressor “holding the other regressors fixed”. In observational data
nobody held anything fixed, so the phrase needs a precise meaning. The
Frisch–Waugh–Lovell (FWL) theorem supplies one. A coefficient in a multiple
regression equals the slope of a simple regression in which the other regressors
have been *projected out* of both the response and the regressor of interest.

## Statement and proof

Partition the model matrix and coefficients as
\[
\X=[\X_1,\X_2],\qquad \bbeta=\begin{pmatrix}\bbeta_1\\\bbeta_2\end{pmatrix},
\]
with \( \X_1 \) of size \( n\times p_1 \) and \( \X_2 \) of size \( n\times p_2 \), and assume
\( \X \) has full column rank. Let \( \M_2 \) project onto \( \C(\X_2) \), and write
\[
\tilde{\X}_1=(\I-\M_2)\X_1,\qquad \tilde{\y}=(\I-\M_2)\y
\]
for the parts of \( \X_1 \) and \( \y \) that \( \X_2 \) cannot explain linearly. Each
column of \( \tilde{\X}_1 \) is the residual vector from regressing the corresponding
column of \( \X_1 \) on \( \X_2 \).

::: {#lem-proj-fwl-split}
[Orthogonalizing a partitioned model space]

Under the assumptions above:

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\X)=\C(\X_2)\dirsum\C(\tilde{\X}_1) \), and the two summands are orthogonal;

2. \( \tilde{\X}_1 \) has full column rank;

3. \( \M=\M_2+\tilde{\M}_1 \), where \( \tilde{\M}_1=\tilde{\X}_1(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T \).
:::

:::

::: {.proof}
(a) The columns of \( \tilde{\X}_1=\X_1-\M_2\X_1 \) lie in \( \C(\X) \) and are orthogonal to
\( \C(\X_2) \). Conversely, for any \( \bb=(\bb_1\T,\bb_2\T)\T \),
\[
\X\bb=\X_1\bb_1+\X_2\bb_2=\tilde{\X}_1\bb_1+(\M_2\X_1\bb_1+\X_2\bb_2),
\]
with the second term in \( \C(\X_2) \). So \( \C(\X)\subseteq\C(\tilde{\X}_1)+\C(\X_2) \).
(b) If \( \tilde{\X}_1\mathbf{a}=\bzero \), then \( \X_1\mathbf{a}=\M_2\X_1\mathbf{a}=\X_2\mathbf{c} \) for some
\( \mathbf{c} \), so \( \X(\mathbf{a}\T,-\mathbf{c}\T)\T=\bzero \). Full column rank of \( \X \) forces
\( \mathbf{a}=\bzero \).
(c) Combine (a) with @thm-proj-sum and @thm-proj-M-formula.
:::

::: {#thm-proj-fwl}
[Frisch–Waugh–Lovell]

Let \( \hbeta=(\hbeta_1\T,\hbeta_2\T)\T \) be the least squares estimate in the full
model. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta_1=(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T\tilde{\y}
         =(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T\y \). So \( \hbeta_1 \) is the coefficient
            from regressing \( \tilde{\y} \) (or \( \y \) itself) on \( \tilde{\X}_1 \);

2. the residual vector from regressing \( \tilde{\y} \) on \( \tilde{\X}_1 \) equals the
           residual vector \( (\I-\M)\y \) of the full model;

3. \( \hbeta_2=(\X_2\T\X_2)^{-1}\X_2\T(\y-\X_1\hbeta_1) \).
:::

:::

::: {.proof}
Let \( \tilde{\bb}=(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T\y \). Since \( \I-\M_2 \) is
symmetric and idempotent,
\( \tilde{\X}_1\T\y=\X_1\T(\I-\M_2)\y=\tilde{\X}_1\T\tilde{\y} \), so the two expressions in
(a) agree. By @lem-proj-fwl-split(c),
\[
\M\y=\M_2\y+\tilde{\X}_1\tilde{\bb}
      =\M_2\y+\X_1\tilde{\bb}-\M_2\X_1\tilde{\bb}
      =\X_1\tilde{\bb}+\M_2(\y-\X_1\tilde{\bb}).
\]
The last term lies in \( \C(\X_2) \), so it equals \( \X_2\mathbf{c} \) with
\( \mathbf{c}=(\X_2\T\X_2)^{-1}\X_2\T(\y-\X_1\tilde{\bb}) \). Thus
\( \M\y=\X_1\tilde{\bb}+\X_2\mathbf{c} \). Because \( \X \) has full column rank, the coefficient
vector representing \( \M\y \) is unique, so \( \hbeta_1=\tilde{\bb} \) and \( \hbeta_2=\mathbf{c} \).
This proves (a) and (c). For (b),
\( \y-\M\y=(\I-\M_2)\y-\tilde{\M}_1\y=\tilde{\y}-\tilde{\X}_1\tilde{\bb} \), which is the
residual of the regression of \( \tilde{\y} \) on \( \tilde{\X}_1 \), because
\( \tilde{\X}_1\tilde{\bb}=\tilde{\M}_1\y=\tilde{\M}_1\tilde{\y} \).
:::

The proof uses nothing beyond the orthogonal decomposition of the model space.
Part (a) contains a small surprise. Residualizing \( \y \) is optional, because the
part of \( \y \) lying in \( \C(\X_2) \) is orthogonal to \( \tilde{\X}_1 \) and drops out of
\( \tilde{\X}_1\T\y \) anyway. Residualizing the regressor of interest is essential.

::: {.idea}
In a multiple regression, \( \hat{\beta}_j \) measures how \( \y \) moves with
*the part of \( \x_j \) that is orthogonal to the other regressors*. A regressor
that is almost a linear combination of the others has a tiny orthogonal part,
so its coefficient rests on very little information. That is collinearity,
seen geometrically.
:::

## Consequences

**Variances.**  The FWL formula writes \( \hbeta_1 \) as a linear function of
\( \y \). If \( \Cov(\Y)=\sigma^2\I \), then
\[
\Cov(\hbeta_1)=\sigma^2(\tilde{\X}_1\T\tilde{\X}_1)^{-1}
\tilde{\X}_1\T\tilde{\X}_1(\tilde{\X}_1\T\tilde{\X}_1)^{-1}
=\sigma^2(\tilde{\X}_1\T\tilde{\X}_1)^{-1}.
\]
Comparing with \( \Cov(\hbeta)=\sigma^2(\X\T\X)^{-1} \) (@thm-lm-moments)
identifies the leading block of the inverse of a partitioned matrix,
\[
\bigl[(\X\T\X)^{-1}\bigr]_{11}
=\bigl(\X_1\T(\I-\M_2)\X_1\bigr)^{-1}
=\bigl(\X_1\T\X_1-\X_1\T\X_2(\X_2\T\X_2)^{-1}\X_2\T\X_1\bigr)^{-1},
\]{#eq-proj-partitioned-inverse}

a Schur-complement identity that [Chapter 1](../ch01-matrix-algebra/index.html) proves by algebra.
For a single regressor \( \x_j \), with all other columns (including the intercept)
in \( \X_2 \),
\[
\Var(\hat{\beta}_j)=\frac{\sigma^2}{\norm{(\I-\M_{(j)})\x_j}^2}
=\frac{\sigma^2}{S_{jj}\,(1-R_j^2)},
\]{#eq-proj-vif-preview}

where \( S_{jj}=\sum_i(x_{ij}-\bar x_j)^2 \) and \( R_j^2 \) is the coefficient of
determination from regressing \( \x_j \) on the other regressors. The second equality
is @exr-proj-vif. The factor \( 1/(1-R_j^2) \) is the *variance inflation
factor* of [Chapter 26](../ch26-collinearity/index.html) (@def-col-vif).

**Centring.**  Take \( \X_2=\bone \). Then \( \M_2=n^{-1}\bone\bone\T \) and
\( \tilde{\X}_1 \), \( \tilde{\y} \) are the centred regressors and response. FWL says the
slopes in a model with an intercept are the coefficients of a no-intercept
regression of centred \( \y \) on centred regressors, and part (c) gives the
intercept \( \bar y-\bar{\x}\T\hbeta_1 \). @exm-proj-simple-nested was the case
\( p_1=1 \).

**Detrending and deseasonalizing.**  The original application was
economic time series. Frisch and Waugh (1933) showed that regressing
on a time trend together with other variables gives the same slope as first
detrending every series. Lovell (1963) proved the analogous
statement for seasonal dummy variables, and noted a trap: if only the
*response* is seasonally adjusted and the regressors are not, the result is
generally wrong. Part (a) of the theorem shows why. What must be residualized
is the regressor.

**A warning about standard errors.**  If the residualized regression of
\( \tilde{\y} \) on \( \tilde{\X}_1 \) is fitted with ordinary software, the coefficient and
residuals are right, but the software divides the residual sum of squares by
\( n-p_1 \) instead of \( n-p_1-p_2 \). It does not know that \( p_2 \) dimensions were
used up in forming \( \tilde{\y} \). The reported standard error is too small by the
factor \( \sqrt{(n-p_1)/(n-p)} \).

::: {#exm-proj-fwl-crime}
[Poverty and murder rates]

For the \( 50 \) US states in 2009, regress the murder rate (per 100,000)
on the poverty rate, the percentage of single-parent households, and the
percentage of the population living in urban areas, with an intercept. The
coefficient of poverty is \( 0.2538 \). On its own, a simple regression
of murder rate on poverty has slope \( 0.4706 \), nearly twice as large.
Poverty and single parenthood are correlated, and the simple regression gives
poverty credit for part of their shared association.

[Figure 6.6.1](06-fwl.html#fig-proj-fwl) shows both regressions. Panel (b) is an
**added-variable plot**: the residual of murder rate after regression on the
other two regressors and the intercept, against the corresponding residual of
poverty. By FWL its least squares slope is exactly the multiple-regression
coefficient. The listing confirms the equality to machine precision.

The standard error of the coefficient in the full model is
\( 0.0934 \). Fitting the residualized regression naively gives
\( 0.0905 \), smaller by the factor \( 1.0321 \),
which is \( \sqrt{(n-1)/(n-4)} \) as predicted.
:::

::: {when-format="html"}
![**Figure 6.6.1.** Murder rate and poverty for 50 US states. (a) The marginal
relationship. (b) The added-variable plot: both variables residualized on the
single-parent percentage, urban percentage and intercept. The slope in (b)
is the multiple-regression coefficient of poverty
(@thm-proj-fwl).](fwl_added_variable.svg){#fig-proj-fwl width=100%}
:::

::: {when-format="pdf"}
![Murder rate and poverty for 50 US states. (a) The marginal
relationship. (b) The added-variable plot: both variables residualized on the
single-parent percentage, urban percentage and intercept. The slope in (b)
is the multiple-regression coefficient of poverty
(@thm-proj-fwl).](fwl_added_variable.pdf){width=100%}
:::

```{.python .run #cell-fwl-fwl}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data
data = data.drop(index="District of Columbia")        # an extreme point; see the exercises of Section 6.6

y = data["murder"].to_numpy()
x1 = data["poverty"].to_numpy()                       # the coefficient we care about
X2 = np.column_stack([np.ones(len(y)), data["single"], data["urban"]])
X = np.column_stack([x1, X2])

def resid(v, Z):
    """Residual of v after projecting onto C(Z)."""
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef

beta_full, *_ = np.linalg.lstsq(X, y, rcond=None)     # regress y on everything

y_tilde = resid(y, X2)                                # (I - M2) y
x_tilde = resid(x1, X2)                               # (I - M2) x1
beta_fwl = (x_tilde @ y_tilde) / (x_tilde @ x_tilde)  # simple regression, no intercept

print(f"full regression:  {beta_full[0]:.6f}")
print(f"FWL regression:   {beta_fwl:.6f}")
```

::: {.remark}
[What FWL does not say]

FWL is an algebraic identity about fitted coefficients. It says what the number
\( \hat{\beta}_j \) *is*: an association with the unexplained part of \( \x_j \). It
does not say that changing \( x_j \) in the world, with the other regressors held
fixed, would change the response by \( \hat{\beta}_j \). That is a causal claim, and it
needs assumptions about how the data arose. [Chapter 25](../ch25-causal-interpretation/index.html) takes
up this question, including cases where adding a regressor to \( \X_2 \) makes the
coefficient *less* interpretable rather than more.
:::

## Fixed effects are a projection {#sec-proj-within}

One of the most widely used estimators in applied statistics is the FWL theorem in
disguise. Suppose units \( k=1,\dots,K \) are observed repeatedly: firms over years,
patients over visits, schools over cohorts. The analyst wants the effect of
time-varying regressors \( \X_1 \) while allowing each unit its own intercept:
\[
\E(\Y)=\X_1\bbeta_1+\mathbf{D}\boldsymbol{\upalpha} ,
\]
where \( \mathbf{D} \) is the \( n\times K \) matrix of unit indicators. Fitting this directly means a
regression with \( K \) extra columns, and \( K \) may be in the thousands.

FWL with \( \X_2=\mathbf{D} \) removes that burden. The projection onto \( \C(\mathbf{D}) \) replaces each
entry by its unit's average (@exm-proj-oneway-M), so \( \I-\M_D \) subtracts unit means.
Therefore \( \hbeta_1 \) is the least squares coefficient from regressing the
*within-unit deviations* \( y_{kt}-\bar{y}_{k\cdot} \) on
\( \x_{kt}-\bar{\x}_{k\cdot} \), with no indicators at all. This is the **within
estimator** of panel-data econometrics. Its standard error needs the degrees-of-freedom
correction described above, with \( p_2=K \), which is large.

::: {#exm-proj-grunfeld}
[Grunfeld's investment data]

Grunfeld's classic panel follows \( 11 \) US firms for
\( 20 \) years, giving \( n=220 \) observations of gross investment,
market value and capital stock. With a separate intercept for each firm, the
coefficients of value and capital are \( 0.1101 \) and \( 0.3100 \).
The same numbers come from demeaning every variable within firm and fitting two
columns. Ignoring the firms altogether (one common intercept) gives
\( 0.1145 \) and \( 0.2275 \) instead. That
comparison uses variation *between* firms, of the sort “large firms invest
more”, which the within estimator deliberately discards.
:::

```{.python .run #cell-within-estimator-within}
import numpy as np
import pandas as pd
import statsmodels.api as sm
df = sm.datasets.grunfeld.load_pandas().data
y = df["invest"].to_numpy()
X1 = df[["value", "capital"]].to_numpy()

# (i) least squares with one indicator column per firm ("dummy variables")
D = pd.get_dummies(df["firm"]).to_numpy(dtype=float)
beta_dummies, *_ = np.linalg.lstsq(np.column_stack([X1, D]), y, rcond=None)

# (ii) FWL: projecting onto C(D) replaces each value by its firm mean,
#      so (I - M_D) subtracts firm means -- the "within" transformation
def demean_by(v, groups):
    v = pd.DataFrame(v)
    return (v - v.groupby(groups.to_numpy()).transform("mean")).to_numpy()

X1_within = demean_by(X1, df["firm"])
y_within = demean_by(y, df["firm"]).ravel()
beta_within, *_ = np.linalg.lstsq(X1_within, y_within, rcond=None)

print("with firm indicators:", beta_dummies[:2])
print("within transformation:", beta_within)
```

Whether discarding between-unit variation is wise depends on why units differ. That
question is taken up for random-effects alternatives in
Chapter 32 and Chapter 33.

Added-variable plots are among the most useful diagnostics for multiple
regression. They show the evidence for a single coefficient in two dimensions,
together with any points that dominate it. They reappear in
[Section 20.1](../ch20-residuals-leverage-influence/01-kinds-of-residuals.html).

## Exercises

### A. Check your understanding

::: {#exr-proj-dc}
[A1]

Refit the model of @exm-proj-fwl-crime with the District of Columbia included. How
does the poverty coefficient change? Draw the added-variable plot and locate the District.
Relate what you see to its leverage (@exm-proj-leverage).
:::

### B. Practice

::: {#exr-proj-vif}
[B1]

Let \( \x_j \) be a column of \( \X \), let \( \X_{(j)} \) be \( \X \) with that column removed, and assume
\( \bone\in\C(\X_{(j)}) \). Show that
\( \norm{(\I-\M_{(j)})\x_j}^2=S_{jj}(1-R_j^2) \), where \( R_j^2 \) is the \( R^2 \) from regressing
\( \x_j \) on \( \X_{(j)} \). Deduce @eq-proj-vif-preview.
:::

::: {.solution}
Since \( \bone\in\C(\X_{(j)}) \), the vector
\( (\I-\M_{(j)})\x_j \) is the residual from regressing \( \x_j \) on \( \X_{(j)} \). Its squared length is that
regression's residual sum of squares, which equals \( \text{SST}(1-R_j^2) \) with
\( \text{SST}=S_{jj} \) by the definition of \( R^2 \). Substituting into the variance formula
from FWL gives @eq-proj-vif-preview.
:::

::: {#exr-proj-fwl-residuals}
[B2]

In @thm-proj-fwl, regressing \( \y \) (not \( \tilde{\y} \)) on \( \tilde{\X}_1 \) gives the
coefficient \( \hbeta_1 \). Show that its residual vector is \( (\I-\M)\y+\M_2\y \), so the
residuals are generally wrong even though the coefficient is right.
:::

::: {.solution}
The fitted vector of the regression of \( \y \)
on \( \tilde{\X}_1 \) is \( \tilde{\M}_1\y \). By @lem-proj-fwl-split,
\( \y-\tilde{\M}_1\y=\y-(\M-\M_2)\y=(\I-\M)\y+\M_2\y \). The extra term \( \M_2\y \) is zero
only if \( \y\perp\C(\X_2) \).
:::

::: {#exr-proj-partial-correlation}
[B3]

For a single regressor \( \x_j \), the *partial correlation* of \( \y \) and \( \x_j \) given the
other regressors is the sample correlation of \( \tilde{\y} \) and \( \tilde{\x}_j \) (both
residualized on \( \X_{(j)} \), which contains \( \bone \)). Let
\( t=\hat{\beta}_j/\text{se}(\hat{\beta}_j) \) with
\( \text{se}(\hat{\beta}_j)^2=\text{SSE}/\bigl[(n-p)\norm{\tilde{\x}_j}^2\bigr] \). Show that the
squared partial correlation equals \( t^2/(t^2+n-p) \).
:::

::: {.solution}
Write \( \tilde{\mathbf{x}}=\tilde{\x}_j \) and
\( \tilde{\mathbf{y}}=\tilde{\y} \). Both have mean zero because \( \bone\in\C(\X_{(j)}) \), so their
correlation is \( r=\tilde{\mathbf{x}}\T\tilde{\mathbf{y}}/(\norm{\tilde{\mathbf{x}}}\norm{\tilde{\mathbf{y}}}) \).
By FWL, \( \hat{\beta}_j=\tilde{\mathbf{x}}\T\tilde{\mathbf{y}}/\norm{\tilde{\mathbf{x}}}^2 \), and the full-model
residual sum of squares equals that of the residualized regression:
\( \text{SSE}=\norm{\tilde{\mathbf{y}}}^2-(\tilde{\mathbf{x}}\T\tilde{\mathbf{y}})^2/\norm{\tilde{\mathbf{x}}}^2
=\norm{\tilde{\mathbf{y}}}^2(1-r^2) \). Hence
\[
t^2=\frac{\hat{\beta}_j^2\,(n-p)\norm{\tilde{\mathbf{x}}}^2}{\text{SSE}}
=\frac{(n-p)\,r^2}{1-r^2},
\]
and solving gives \( r^2=t^2/(t^2+n-p) \).
:::

### C. Going deeper

::: {#exr-proj-fwl-rank-deficient}
[C1]

Drop the full-rank assumption in @thm-proj-fwl. Show that \( \bb_1 \) solves the normal
equations of the regression of \( \tilde{\y} \) on \( \tilde{\X}_1 \) iff there exists \( \bb_2 \) such
that \( (\bb_1\T,\bb_2\T)\T \) is a least squares estimate in the full model.
:::
