# Interpreting coefficients

A fitted regression is usually reported as a table of coefficients. The short answer, "the change in the mean response per unit
change in the regressor, holding the other regressors fixed", is a useful slogan and a
frequent source of error. This section makes the slogan precise within the model, shows
how the meaning of a coefficient depends on units, centring, transformations and the other
columns of the model, and marks where the slogan stops being true.

## What a coefficient measures

::: {#def-lm-coefficient-interpretation}
[Interpretation of a coefficient]

In a linear model \( \E(Y\mid\bm x)=\beta_0+\beta_1x_1+\dots+\beta_kx_k \), suppose the regressor
\( x_j \) enters the model only through its own column, so that no other column is a
function of \( x_j \). For any two regressor profiles \( \bm x \) and \( \bm x' \) that agree in every
regressor except the \( j \)th, and have \( x_j'=x_j+1 \),
\[
\E(Y\mid\bm x')-\E(Y\mid\bm x)=\beta_j .
\]{#eq-lm-partial-effect}

The coefficient \( \beta_j \) is the **partial effect** of \( x_j \) *in this model*: the difference in
mean response per unit difference in \( x_j \), between profiles that agree on all the other
regressors in the model. Its estimate \( \hat{\beta}_j \) is the least squares estimate of this difference.
:::

Three qualifications are built into the definition.

1. **It compares profiles, not interventions.** @eq-lm-partial-effect compares the mean response
   of two kinds of case that exist in the population the model describes. Whether changing
   \( x_j \) *in a given case* would change its response by \( \beta_j \) is a causal question,
   which the model alone cannot answer (Chapter 25).

2. **"The other regressors" means the ones in the model.** Adding or removing a
   regressor changes the comparison being made, and so changes \( \beta_j \), unless the
   regressors are orthogonal. @prp-lm-omitted gives the exact relation.

3. **It needs \( x_j \) to vary alone.** If \( x_j \) also appears in other columns, as in a
   polynomial or an interaction, there are no two profiles that differ only in the
   \( x_j \) column, and \( \beta_j \) alone does not measure the effect of \( x_j \).

The geometric content of "holding the others fixed" in a data set, where nobody held
anything fixed, is the Frisch–Waugh–Lovell theorem of
[Section 6.6](../ch06-projections/06-fwl.html): \( \hat{\beta}_j \) is the slope of the response on
the part of \( x_j \) that is not linearly explained by the other regressors.

## Units, rescaling and reparameterization

The coefficient \( \beta_j \) is measured in units of the response per unit of \( x_j \).
Changing units rescales it. More generally, any invertible linear recoding of the columns
changes the coefficients in a predictable way and leaves the fit untouched.

::: {#prp-lm-reparameterization}
[Reparameterization]

Let \( \X \) be \( n\times p \) of rank \( p \), let \( \bm K \) be a nonsingular \( p\times p \) matrix, and put
\( \bm Z=\X\bm K \). Then \( \bm Z \) has rank \( p \), and the least squares estimates from the two
model matrices are related by
\[
\hbeta_Z=\bm K^{-1}\hbeta_X .
\]
The fitted values, residuals, SSE and \( s^2 \) are the same for both, and the estimated
covariance matrices satisfy \( s^2(\bm Z\T\bm Z)^{-1}=\bm K^{-1}\bigl[s^2(\X\T\X)^{-1}\bigr](\bm K^{-1})\T \).
:::

::: {.proof}
\( \rank(\bm Z)=\rank(\X)=p \) by @prp-mat-rank-product(b). Then
\[
\begin{aligned}
\hbeta_Z&=(\bm K\T\X\T\X\bm K)^{-1}\bm K\T\X\T\y\\
&=\bm K^{-1}(\X\T\X)^{-1}(\bm K\T)^{-1}\bm K\T\X\T\y=\bm K^{-1}\hbeta_X,
\end{aligned}
\]
and \( \bm Z\hbeta_Z=\X\bm K\bm K^{-1}\hbeta_X=\X\hbeta_X \). Equal fitted values give equal residuals,
SSE and \( s^2 \). The covariance identity follows from \( (\bm Z\T\bm Z)^{-1}=\bm K^{-1}(\X\T\X)^{-1}(\bm K^{-1})\T \).
:::

With \( \bm K=\diag(1,c_1,\dots,c_k) \), measuring \( x_j \) in units \( 1/c_j \) times as large divides its
coefficient and its standard error by \( c_j \), and leaves everything else alone. Centring,
effect coding of a factor (@exr-lm-effect-coding) and many other changes of coding are also of the
form \( \X\bm K \). [Section 6.7](../ch06-projections/07-reparameterization.html) shows that what
really stays fixed is the column space of \( \X \).

**The intercept and centring.** The intercept \( \beta_0 \) is the mean response at the profile
where every regressor is zero. That profile may be impossible (zero income) or far outside the
data (an August sea temperature of \( 0 \)°C in @exm-lm-elnino). If the regressors are centred,
the intercept becomes the mean response at the average profile \( \bar{\bm x} \), and its
estimate is \( \bar{y} \) (@prp-lm-centred). Centring does not change the slopes. It makes the
intercept interpretable, and it makes \( \hat{\beta}_0 \) uncorrelated with the slopes
(@exr-lm-centred-uncorrelated). With interactions or polynomials, centring also changes
what the lower-order coefficients mean, as the next subsections show.

**Standardized coefficients.** Measuring every variable in units of its own sample
standard deviation gives the **standardized coefficients**
\( \hat{\beta}_j^{*}=\hat{\beta}_j\,s_{x_j}/s_y \): the difference in the response, in standard
deviations, per standard deviation of \( x_j \). They make coefficients unit-free and are
sometimes used to rank the "importance" of regressors. That use is fragile. A standard
deviation is a property of the sample, not of the relation: the same relation observed over a
wider range of \( x_j \) gives a larger standardized coefficient.

## Logarithms

When the response or a regressor is on a logarithmic scale, the coefficients describe
proportional rather than absolute differences.

- **Log–log**, \( \E(\log Y)=\beta_0+\beta_1\log x+\dots \): a \( 1 \) percent difference in \( x \)
  goes with approximately a \( \beta_1 \) percent difference in \( Y \) (exactly, a factor
  \( 1.01^{\beta_1} \) in the geometric mean). Economists call \( \beta_1 \) an **elasticity**.

- **Log–level**, \( \E(\log Y)=\beta_0+\beta_1x+\dots \): a unit difference in \( x \) goes with a
  factor \( e^{\beta_1} \) in the geometric mean of \( Y \), a difference of
  \( 100(e^{\beta_1}-1) \) percent, approximately \( 100\beta_1 \) percent when \( \beta_1 \) is small.

- **Level–log**, \( \E(Y)=\beta_0+\beta_1\log x+\dots \): a \( 1 \) percent difference in \( x \) goes
  with a difference of about \( \beta_1/100 \) in \( Y \).

The phrase "geometric mean" matters. \( \E(\log Y) \) is not \( \log\E(Y) \), so a model for the mean of
\( \log Y \) is not a model for the mean of \( Y \). Chapter 22 returns to back-transformation.

::: {#exm-lm-money-interpret}
[Reading the money demand coefficients]

In @exm-lm-money the response is \( \log \) real money and the first regressor is \( \log \) real
income, so \( \hat{\beta}_1=1.2958 \) is an income elasticity. Comparing
quarters with the same interest rates, \( 1 \) percent more real income goes with about
\( 1.3 \) percent more real money. The interest rates are measured as fractions, so the
bond rate coefficient \( -2.6163 \) is per unit of the fraction, that is, per
\( 100 \) percentage points, an absurdly large change. Rescaling the rate to percentage points
(@prp-lm-reparameterization with \( c=100 \)) divides the coefficient by \( 100 \): a bond rate one
percentage point higher, with income and the deposit rate the same, goes with money holdings
lower by the factor \( e^{-0.026163} \), a difference of \( -2.58 \) percent. The
deposit rate has a positive coefficient, which is plausible if deposits are part of the
money measured, since a higher deposit rate makes holding them more attractive. But the
coefficient is smaller than its standard error.

The standardized coefficients of income, bond rate and deposit rate are \( 0.619 \),
\( -0.529 \) and \( 0.057 \). The sample standard deviation of log income is
\( 0.0728 \), about \( 7 \) percent, and that of the bond rate is
\( 0.0308 \), about \( 3 \) percentage points; within those ranges of variation,
income and the bond rate move money holdings by similar amounts.
:::

```{.python .run #cell-money-demand-rescale}
import numpy as np
import statsmodels.api as sm
dk = sm.datasets.danish_data.load_pandas().data
y = dk["lrm"].to_numpy()
X = np.column_stack([np.ones(len(y)), dk["lry"], dk["ibo"], dk["ide"]])
XtX = X.T @ X
Xty = X.T @ y
beta_hat = np.linalg.solve(XtX, Xty)       # solve X^T X b = X^T y

X_pct = X.copy()
X_pct[:, 2:] *= 100                         # rates in percentage points
beta_pct = np.linalg.solve(X_pct.T @ X_pct, X_pct.T @ y)
print("rates as fractions:", np.round(beta_hat, 4))
print("rates in percent:  ", np.round(beta_pct, 5))
print("same fitted values:", np.allclose(X @ beta_hat, X_pct @ beta_pct))

sd = X[:, 1:].std(axis=0, ddof=1)
standardized = beta_hat[1:] * sd / y.std(ddof=1)
print("standardized slopes:", np.round(standardized, 3))
```

## Terms that move together

If \( x \) enters through both \( x \) and \( x^2 \), the mean response
\( \beta_0+\beta_1x+\beta_2x^2 \) changes at the rate
\[
\frac{\partial\E(Y\mid x)}{\partial x}=\beta_1+2\beta_2x ,
\]
which depends on \( x \). The coefficient \( \beta_1 \) alone is this rate at \( x=0 \), and
\( \beta_2 \) is half the rate at which the slope itself changes. Neither is "the effect of
\( x \)". In @exm-lm-co2 the growth rate of carbon dioxide was
\( 0.86 \) ppm per year in 1960 and \( 1.81 \) ppm per year in 2000. The
coefficient \( \hat{\beta}_1=1.3344 \) is the growth rate in 1980, the year chosen as the
origin of \( t_c \), and would have been different had the origin been chosen differently.

Interactions behave the same way. In @eq-lm-interaction the slope in \( x_1 \) is
\( \beta_1+\beta_3x_2 \). Reporting "the effect of \( x_1 \)" requires choosing a value of
\( x_2 \), and \( \beta_1 \) is the choice \( x_2=0 \). Centring \( x_2 \) at a representative value before
forming the product makes \( \beta_1 \) the slope at that value. The fitted surface does not
change (@exr-lm-interaction-centring).

## Indicator variables

With reference coding, the coefficient of an indicator is the difference in mean
response between its category and the reference category, among profiles that agree on the
other regressors. Changing the reference changes every indicator coefficient but not the
fit, since the reference coding for one reference category is a reparameterization of the
coding for any other.

::: {#exm-lm-health-adjusted}
[Self-rated health, with and without adjustment]

In @exm-lm-health-coding, people reporting poor health made \( 3.160 \) more
physician visits a year, on average, than those reporting excellent health. Add two
regressors: a chronic-disease score (the variable `disea`) and an indicator of a physical limitation. The
coefficient of poor health falls to \( 1.554 \), and that of fair health from
\( 1.058 \) to \( 0.202 \). The good-health coefficient changes sign, to
\( -0.081 \). The adjusted comparison is between people with the same disease score and
limitation status. Much of the difference in visits between health categories goes with
differences in measured disease, and the self-rating adds less once disease is taken into
account. Neither set of numbers is "the" effect of self-rated health. They answer different
questions.

The data record the same people in several years, so the cases are not independent and the
second-moment assumptions do not hold exactly; Chapter 33 treats such clustered data.
:::

```{.python .run #cell-rand-dummies-adjusted}
import numpy as np
import statsmodels.api as sm
rand = sm.datasets.randhie.load_pandas().data
y = rand["mdvis"].to_numpy(dtype=float)
D = rand[["hlthg", "hlthf", "hlthp"]].to_numpy(dtype=float)   # good, fair, poor
X = np.column_stack([np.ones(len(y)), D])                      # excellent = reference

Xa = np.column_stack([X, rand["disea"], rand["physlm"]])
beta_a = np.linalg.solve(Xa.T @ Xa, Xa.T @ y)
print("adjusted for chronic-disease score and physical limitation:")
print(np.round(beta_a, 4))
```

## Omitted regressors

How does a coefficient change when other regressors are added or removed? Partition the
model matrix as \( \X=[\X_1,\X_2] \), with \( \hbeta=(\hbeta_1\T,\hbeta_2\T)\T \) the estimate from the full
("long") regression, and let \( \tilde{\bbeta}_1=(\X_1\T\X_1)^{-1}\X_1\T\y \) be the estimate from the
"short" regression on \( \X_1 \) alone.

::: {#prp-lm-omitted}
[Omitted regressors]

Let \( \X=[\X_1,\X_2] \) have full rank, and let \( \bm\Pi=(\X_1\T\X_1)^{-1}\X_1\T\X_2 \) be the
coefficient matrix from regressing each column of \( \X_2 \) on \( \X_1 \). Then:

::: {.enumerate options="label=(\alph*)"}
1. for every data vector \( \y \), \( \tilde{\bbeta}_1=\hbeta_1+\bm\Pi\hbeta_2 \);

2. if \( \E(\Y)=\X_1\bbeta_1+\X_2\bbeta_2 \), then \( \E(\tilde{\bbeta}_1)=\bbeta_1+\bm\Pi\bbeta_2 \).
:::

In particular the short and long coefficients agree, for all \( \y \), iff \( \X_1\T\X_2=\bzero \).
:::

::: {.proof}
*(a)* Write \( \y=\X_1\hbeta_1+\X_2\hbeta_2+\he \), where \( \X\T\he=\bzero \) and so
\( \X_1\T\he=\bzero \). Then
\[
\begin{aligned}
\tilde{\bbeta}_1&=(\X_1\T\X_1)^{-1}\X_1\T\y\\
&=\hbeta_1+(\X_1\T\X_1)^{-1}\X_1\T\X_2\hbeta_2+(\X_1\T\X_1)^{-1}\X_1\T\he
=\hbeta_1+\bm\Pi\hbeta_2 .
\end{aligned}
\]
*(b)* Take expectations in (a), using \( \E(\hbeta)=\bbeta \) (@thm-lm-moments), or apply
@prp-lm-misspecified(a) with \( \bm\delta=\X_2\bbeta_2 \). For the last claim, if \( \X_1\T\X_2=\bzero \) then
\( \bm\Pi=\bzero \). Conversely, if \( \bm\Pi\ne\bzero \), take \( \y=\X_2\bm c \) with \( \bm\Pi\bm c\ne\bzero \): the
long regression gives \( \hbeta_2=\bm c \), and the two coefficients differ by \( \bm\Pi\bm c \).
:::

Part (b) is the **omitted-variable bias** formula. The short regression does not estimate
\( \bbeta_1 \) but \( \bbeta_1+\bm\Pi\bbeta_2 \). Each omitted regressor contributes its own
coefficient times its association with the included ones. The bias is zero if the omitted
regressors have no effect or are uncorrelated with the included ones, in the sense
\( \X_1\T\X_2=\bzero \); otherwise it is there, and it can have either sign. Part (a) says the same
thing about the numbers actually computed from a data set, with no model assumptions at all.

::: {#exm-lm-money-omitted}
[Income without interest rates]

Regressing log money on log income alone, with an intercept, gives an income elasticity of
\( 1.819 \), against \( 1.2958 \) when both interest rates are included.
@prp-lm-omitted(a) accounts for the whole difference. Regressing each rate on log income gives
slopes \( -0.2215 \) (bond rate) and \( -0.0917 \) (deposit rate): in these data,
quarters with higher income had lower interest rates, and the correlation of log income with the
bond rate is \( -0.52 \). The short-regression elasticity therefore absorbs
\[
(-0.2215)(-2.6163)+(-0.0917)(0.6186)=0.5229,
\]
mostly the association of higher income with lower bond rates, which in turn go with more
money. Which elasticity is the "right" one depends on the question. The long regression
compares quarters with equal interest rates. The short one compares quarters as they came,
with whatever rates accompanied a given income.
:::

```{.python .run #cell-money-demand-omitted}
X1 = X[:, :2]                               # intercept and log income
X2 = X[:, 2:]                               # the two interest rates
b_short = np.linalg.solve(X1.T @ X1, X1.T @ y)
A = np.linalg.solve(X1.T @ X1, X1.T @ X2)   # regress each rate on (1, lry)
print("short regression slope on income:", round(b_short[1], 4))
print("long slope + A @ (rate coefs):   ", round(beta_hat[1] + A[1] @ beta_hat[2:], 4))
print("auxiliary slopes of the rates on income:", np.round(A[1], 4))
```

::: {.warning}
[Adjustment is not a cure]

@prp-lm-omitted suggests a tempting strategy: include every variable that might be relevant, so
that nothing is omitted. Two problems remain. First, some variables should *not* be held fixed.
If a treatment works partly by changing an intermediate variable, adjusting for that variable
removes part of the effect one wants to measure. Adjusting for a common consequence of the
regressor and the response can create an association where there is none. Second, adding
regressors that are nearly collinear with \( x_j \) inflates the variance of \( \hat{\beta}_j \)
(Chapter 26). Which variables to adjust for is a question about how the data were generated, not
about the fit. Chapter 25 develops the tools for answering it.
:::

## Exercises

### A. Check your understanding

::: {#exr-lm-fahrenheit}
[A1]

In @exm-lm-elnino, suppose both temperatures are converted to degrees Fahrenheit,
\( F=32+1.8\,C \). What happens to the slope, the intercept, \( R^2 \) and \( s \)? Answer
using @prp-lm-reparameterization, noting that the response is also transformed.
:::

### B. Practice

::: {#exr-lm-standardized}
[B1]

Show that the standardized coefficients are the slopes obtained by regressing the
standardized response \( (y_i-\bar{y})/s_y \) on the standardized regressors
\( (x_{ij}-\bar{x}_j)/s_{x_j} \), without an intercept. For a single regressor, show that the
standardized slope is the correlation \( r \).
:::

::: {.solution}
Let \( \X_c \) be the centred regressors and \( \bm D=\diag(s_{x_1},\dots,s_{x_k}) \). The standardized
regressors are \( \X_c\bm D^{-1} \) and the standardized response is \( (\y-\bar{y}\bone)/s_y \). By
@prp-lm-centred the slopes of \( \y \) on \( \X_c \) are \( \hbeta_1 \). Dividing the response by \( s_y \)
divides them by \( s_y \), since least squares is linear in the response; replacing \( \X_c \) by
\( \X_c\bm D^{-1} \) multiplies them by \( \bm D \) (@prp-lm-reparameterization with \( \bm K=\bm D^{-1} \)). The
result is \( \bm D\hbeta_1/s_y \), whose entries are \( \hat{\beta}_js_{x_j}/s_y \). For one regressor this is
\( (S_{xy}/S_{xx})\sqrt{S_{xx}/S_{yy}}=r \), as in @eq-lm-toward-mean.
:::

::: {#exr-lm-centred-uncorrelated}
[B2]

In the centred parameterization \( \E(\Y)=\alpha\bone+\X_c\bbeta_1 \) of @prp-lm-centred, show that
\( \Cov(\hat{\alpha},\hbeta_1)=\bzero \), \( \Var(\hat{\alpha})=\sigma^2/n \) and
\( \Cov(\hbeta_1)=\sigma^2(\X_c\T\X_c)^{-1} \). Deduce that
\( \Var(\hat{\beta}_0)=\sigma^2\bigl(1/n+\bar{\x}\T(\X_c\T\X_c)^{-1}\bar{\x}\bigr) \), which generalizes @cor-lm-simple-moments.
:::

::: {.solution}
The matrix \( [\bone,\X_c] \) has full rank and \( [\bone,\X_c]\T[\bone,\X_c]=\diag(n,\X_c\T\X_c) \),
because \( \bone\T\X_c=\bzero\T \). By @thm-lm-moments the covariance of \( (\hat{\alpha},\hbeta_1\T)\T \) is
\( \sigma^2\diag\bigl(1/n,(\X_c\T\X_c)^{-1}\bigr) \), which gives all three claims. Since
\( \hat{\beta}_0=\hat{\alpha}-\bar{\x}\T\hbeta_1 \), @cor-lm-linear-functions gives
\( \Var(\hat{\beta}_0)=\sigma^2/n+\bar{\x}\T\sigma^2(\X_c\T\X_c)^{-1}\bar{\x} \), the cross term being zero.
With one regressor, \( \X_c\T\X_c=S_{xx} \).
:::

::: {#exr-lm-omitted-sign}
[B3]

In a model with an intercept and two regressors, let the short regression omit \( x_2 \). Show
that the bias of the slope on \( x_1 \) is \( \beta_2S_{12}/S_{11} \), where
\( S_{12}=\sum_i(x_{i1}-\bar{x}_1)(x_{i2}-\bar{x}_2) \) and \( S_{11}=\sum_i(x_{i1}-\bar{x}_1)^2 \). Tabulate
the sign of the bias according to the signs of \( \beta_2 \) and the correlation of the two regressors.
:::

### C. Going deeper

::: {#exr-lm-sign-reversal}
[C1]

Construct a data set with \( n=8 \) cases, two regressors and an intercept, in which the slope of
\( y \) on \( x_1 \) is positive in the simple regression and negative in the regression on both
\( x_1 \) and \( x_2 \). Explain the reversal with @prp-lm-omitted, and draw the scatterplot of \( y \)
against \( x_1 \) with the cases labelled by \( x_2 \).
:::
