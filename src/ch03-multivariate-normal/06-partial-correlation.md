# Partial correlation

Two variables can be strongly correlated only because both depend on a third.
@exm-mvn-trivariate showed a positive association between \( Y_1 \) and \( Y_2 \)
turning negative once \( Y_3 \) is held fixed. The conditional covariance matrix
\( \bSigma_{11\cdot2} \) measures association that remains after conditioning, and its
correlation form is the partial correlation.

## Definition and interpretations

::: {#def-mvn-partial-correlation}
[Partial correlation]

Let \( \Y=(\Y_1\T,\Y_2\T)\T \) be a random vector with covariance matrix \( \bSigma \),
partitioned as in [Section 3.5](05-conditional.html), and write
\( \bSigma_{11\cdot2}=\bSigma_{11}-\bSigma_{12}\bSigma_{22}\ginv\bSigma_{21}=(\sigma_{ij\cdot2}) \).
The **partial covariance** of \( Y_i \) and \( Y_j \) (two entries of \( \Y_1 \)) given \( \Y_2 \) is
\( \sigma_{ij\cdot2} \), and their **partial correlation** given \( \Y_2 \) is
\[
\rho_{ij\cdot2}=\frac{\sigma_{ij\cdot2}}{\sqrt{\sigma_{ii\cdot2}\,\sigma_{jj\cdot2}}},
\]
provided \( \sigma_{ii\cdot2} \) and \( \sigma_{jj\cdot2} \) are positive. When the
conditioning variables are \( Y_k,Y_l,\dots \) we also write \( \rho_{ij\cdot kl\dots} \). The
matrix of partial correlations is
\( \bD^{-1/2}\bSigma_{11\cdot2}\bD^{-1/2} \) with \( \bD=\diag(\bSigma_{11\cdot2}) \), the
correlation matrix (@def-rv-correlation) built from \( \bSigma_{11\cdot2} \).
:::

The definition uses only \( \bSigma \), so partial correlations exist for any random
vector with finite second moments. The next result gives them two meanings, one
that always applies and one that needs normality.

::: {#prp-mvn-partial-meaning}
In the setting of @def-mvn-partial-correlation:

::: {.enumerate options="label=(\alph*)"}
1. \( \bSigma_{11\cdot2} \) is the covariance matrix of the prediction errors
           \( \W=\Y_1-\bmu_1-\bSigma_{12}\bSigma_{22}\ginv(\Y_2-\bmu_2) \), and \( \W \) is uncorrelated
           with \( \Y_2 \). If \( \bSigma_{22} \) is positive definite, the entries of \( \W \) are the
           errors of the best linear predictors \( L(Y_i\mid\Y_2) \), so \( \rho_{ij\cdot2} \) is the
           correlation between \( Y_i-L(Y_i\mid\Y_2) \) and \( Y_j-L(Y_j\mid\Y_2) \).

2. If \( \Y \) is jointly normal, \( \rho_{ij\cdot2} \) is the correlation of \( Y_i \) and \( Y_j \)
           in their conditional distribution given \( \Y_2=\y_2 \), the same for every \( \y_2 \).
           Moreover \( \rho_{ij\cdot2}=0 \) iff \( Y_i \) and \( Y_j \) are conditionally independent
           given \( \Y_2 \).
:::

:::

::: {.proof}
(a) This is @thm-rv-blp, restated:
\( \Cov(\W,\Y_2)=\bzero \) and \( \Cov(\W)=\bSigma_{11\cdot2} \) are the moment calculations in the
proof of @thm-mvn-conditional(a), which use no normality. When \( \bSigma_{22} \) is
positive definite, the \( i \)th entry of \( \bmu_1+\bSigma_{12}\bSigma_{22}^{-1}(\Y_2-\bmu_2) \) is
\( L(Y_i\mid\Y_2) \) by @thm-proj-blp.
(b) By @thm-mvn-conditional, the conditional distribution of \( (Y_i,Y_j) \) given
\( \Y_2=\y_2 \) is bivariate normal with covariance matrix the corresponding \( 2\times2 \) block of
\( \bSigma_{11\cdot2} \). Its correlation is \( \rho_{ij\cdot2} \), and by
@thm-mvn-independence applied to this conditional normal distribution, the two are
independent iff the off-diagonal entry is zero.
:::

Part (a) is the definition to remember: *a partial correlation is the correlation of
two residuals*. Remove from each variable the part that is linearly predictable from the
conditioning variables, and correlate what is left. This is the population version of the
Frisch–Waugh–Lovell construction (@thm-proj-fwl).

## Computing partial correlations

Two routes avoid computing \( \bSigma_{11\cdot2} \) in full.

With a single conditioning variable, \( \bSigma_{11\cdot2} \) is a \( 2\times2 \) matrix with entries
\( \sigma_{ij}-\sigma_{i3}\sigma_{j3}/\sigma_{33} \), and in correlation units this is the familiar
three-variable formula of @exr-rv-partial-three,
\[
\rho_{12\cdot3}=\frac{\rho_{12}-\rho_{13}\rho_{23}}{\sqrt{(1-\rho_{13}^2)(1-\rho_{23}^2)}},
\qquad \lvert\rho_{13}\rvert<1,\ \lvert\rho_{23}\rvert<1 .
\]{#eq-mvn-partial-recursion}

The residual interpretation shows that the same formula can be applied one conditioning
variable at a time.

::: {#prp-mvn-partial-recursion}
[Recursion]

Let \( S \) be a set of indices and \( k\notin S\cup\{1,2\} \) one more. Suppose
\( \sigma_{11\cdot S} \), \( \sigma_{22\cdot S} \) and \( \sigma_{kk\cdot S} \) are positive and
\( \lvert\rho_{1k\cdot S}\rvert<1 \), \( \lvert\rho_{2k\cdot S}\rvert<1 \). Then
\[
\sigma_{12\cdot Sk}=\sigma_{12\cdot S}-\frac{\sigma_{1k\cdot S}\,\sigma_{2k\cdot S}}{\sigma_{kk\cdot S}},
\qquad
\rho_{12\cdot Sk}=\frac{\rho_{12\cdot S}-\rho_{1k\cdot S}\,\rho_{2k\cdot S}}
{\sqrt{(1-\rho_{1k\cdot S}^2)(1-\rho_{2k\cdot S}^2)}} .
\]
With \( S \) empty this is @eq-mvn-partial-recursion.
:::

::: {.proof}
No normality is needed. For \( i\in\{1,2,k\} \) let \( e_i \) be the error of the best linear
predictor of \( Y_i \) from \( \Y_S \), as in @prp-mvn-partial-meaning(a). Then each \( e_i \) has
mean zero, is uncorrelated with \( \Y_S \), and \( \Cov(e_i,e_j)=\sigma_{ij\cdot S} \). Put
\( r_i=e_i-c_ie_k \) with \( c_i=\sigma_{ik\cdot S}/\sigma_{kk\cdot S} \), for \( i=1,2 \).

We claim that \( r_i \) is the error of the best linear predictor of \( Y_i \) from \( (\Y_S,Y_k) \).
First, \( Y_i-r_i \) is an affine function of \( (\Y_S,Y_k) \), because \( e_i \) and \( e_k \) differ from
\( Y_i \) and \( Y_k \) by affine functions of \( \Y_S \). Second, \( r_i \) has mean zero and is uncorrelated
with \( \Y_S \) and with \( e_k \), by the choice of \( c_i \), hence with \( Y_k=e_k+(\text{affine in }\Y_S) \).
Now let \( r_i' \) be the error from @prp-mvn-partial-meaning(a) for predicting \( Y_i \) from
\( (\Y_S,Y_k) \), which has the same two properties. The difference \( d=r_i-r_i' \) is an affine
function of \( (\Y_S,Y_k) \) with mean zero and is uncorrelated with \( (\Y_S,Y_k) \), so
\( \Var(d)=\Cov(d,d)=0 \). Thus \( r_i=r_i' \) with probability one.

Therefore \( \sigma_{12\cdot Sk}=\Cov(r_1,r_2)=\sigma_{12\cdot S}-c_2\sigma_{1k\cdot S}-c_1\sigma_{2k\cdot S}
+c_1c_2\sigma_{kk\cdot S} \), which simplifies to the first formula, and similarly
\( \sigma_{ii\cdot Sk}=\sigma_{ii\cdot S}(1-\rho_{ik\cdot S}^2) \), which is positive by hypothesis.
Dividing \( \sigma_{12\cdot Sk} \) by \( \sqrt{\sigma_{11\cdot Sk}\sigma_{22\cdot Sk}} \) and then numerator
and denominator by \( \sqrt{\sigma_{11\cdot S}\sigma_{22\cdot S}} \) gives the second.
:::

::: {#prp-mvn-precision}
[Partial correlations from the precision matrix]

Let \( \bSigma \) be positive definite and \( \bm\Omega=\bSigma^{-1}=(\omega_{ij}) \). The partial
correlation of \( Y_i \) and \( Y_j \) given all the other entries of \( \Y \) is
\[
\rho_{ij\cdot\text{rest}}=-\frac{\omega_{ij}}{\sqrt{\omega_{ii}\,\omega_{jj}}},
\]
and the conditional variance of \( Y_i \) given all the others is \( 1/\omega_{ii} \).
:::

::: {.proof}
Permute \( \Y \) so that \( (Y_i,Y_j) \) comes first. By @thm-mvn-conditional(c) (a statement
about matrices, valid whether or not \( \Y \) is normal), the leading \( 2\times2 \) block of
\( \bm\Omega \) is \( \bSigma_{11\cdot2}^{-1} \). Inverting,
\[
\bSigma_{11\cdot2}=\begin{pmatrix}\omega_{ii}&\omega_{ij}\\\omega_{ij}&\omega_{jj}\end{pmatrix}^{-1}
=\frac{1}{\omega_{ii}\omega_{jj}-\omega_{ij}^2}
\begin{pmatrix}\omega_{jj}&-\omega_{ij}\\-\omega_{ij}&\omega_{ii}\end{pmatrix},
\]
whose correlation is \( -\omega_{ij}/\sqrt{\omega_{ii}\omega_{jj}} \). The same argument with a
\( 1\times1 \) leading block gives the conditional variance.
:::

@prp-mvn-precision has a striking consequence for normal vectors: *a zero entry
\( \omega_{ij} \) of the inverse covariance matrix means that \( Y_i \) and \( Y_j \) are conditionally
independent given all the other variables*. Zeros in \( \bSigma \) describe marginal
independence, zeros in \( \bSigma^{-1} \) describe conditional independence. The pattern of zeros in
\( \bSigma^{-1} \) defines a graph on the variables, which is the starting point of Gaussian
graphical models (Dempster 1972; Lauritzen 1996).

::: {#exm-mvn-trivariate-partial}
[The trivariate example, continued]

For \( \bSigma \) of @exm-mvn-trivariate, the listing gives the covariance of \( (Y_1,Y_2) \)
given \( Y_3 \) as
\( \begin{psmallmatrix}2&-1\\-1&1\end{psmallmatrix} \).
So \( \rho_{12}=0.289 \) while \( \rho_{12\cdot3}=-0.707 \).
The script confirms that @eq-mvn-partial-recursion, @prp-mvn-precision, and
the correlation of simulated residuals all give this value.
:::

## Sample partial correlations

For data, replace \( \bSigma \) by the sample covariance matrix. Because the empirical
distribution of the rows is a probability distribution in its own right
([Section 6.11](../ch06-projections/11-population.html)), all the identities above hold exactly for the sample
versions. In particular, for two variables \( \y,\x \) and conditioning variables collected with
an intercept in \( \bm Z \), the following numbers are *identical*:

::: {.enumerate options="label=(\arabic*)"}
1. the correlation of the residuals of \( \y \) and of \( \x \) after least squares regression
           on \( \bm Z \);

2. the recursion @eq-mvn-partial-recursion applied to sample correlations
           (one conditioning variable) or repeatedly (several);

3. \( -w_{12}/\sqrt{w_{11}w_{22}} \) computed from the inverse of the sample covariance or
           correlation matrix of \( (\y,\x,\text{conditioning variables}) \);

4. \( \operatorname{sign}(t)\sqrt{t^2/(t^2+n-p)} \), where \( t \) is the \( t \) statistic of \( \x \) in the
           regression of \( \y \) on \( [\bm Z,\x] \) with \( p \) columns (@exr-proj-partial-correlation).
:::

::: {#exm-mvn-longley}
[Employment and unemployment]

Longley's macroeconomic series record, among other things, total employment and the number
unemployed in the United States for the \( 16 \) years
\( 1947 \)–\( 1962 \). Their correlation is
\( 0.502 \): years with more people employed had more people unemployed. Both
series grew with the population. Their correlations with the calendar year are
\( 0.971 \) and \( 0.668 \). Adjusting for year, the partial correlation
is \( -0.829 \). Once the common trend is removed, a year with unemployment above
its trend is a year with employment below its trend, as one would expect.
[Figure 3.6.1](06-partial-correlation.html#fig-mvn-partial-correlation) shows both scatterplots. The listing computes the partial
correlation in all four ways above. They agree to ten decimal places. The slope in panel (b),
\( -0.996 \), is the coefficient of unemployment in the regression of employment on
unemployment and year, as @thm-proj-fwl requires. Its \( t \) statistic is
\( -5.34 \).

Removing a linear time trend before correlating two series is the problem that
Frisch and Waugh (1933) solved, and it is a special case of partial correlation.
:::

::: {when-format="html"}
![**Figure 3.6.1.** Employment and unemployment in the United States, 1947–1962, coloured by year.
(a) The raw series are positively correlated because both trend upwards. (b) After
removing the linear trend in year from each, the residuals are strongly negatively
correlated. Their correlation is the partial correlation given year, and the fitted slope
is the multiple-regression coefficient.](partial_correlation.svg){#fig-mvn-partial-correlation width=100%}
:::

::: {when-format="pdf"}
![Employment and unemployment in the United States, 1947–1962, coloured by year.
(a) The raw series are positively correlated because both trend upwards. (b) After
removing the linear trend in year from each, the residuals are strongly negatively
correlated. Their correlation is the partial correlation given year, and the fitted slope
is the multiple-regression coefficient.](partial_correlation.pdf){width=100%}
:::

```{.python .run #cell-partial-correlation-pcor}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.longley.load_pandas().data

emp = data["TOTEMP"].to_numpy() / 1000       # thousands -> millions of persons
unemp = data["UNEMP"].to_numpy() / 1000
year = data["YEAR"].to_numpy()
n = len(emp)
Z = np.column_stack([np.ones(n), year])      # the variables we adjust for

def resid(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef

# 1. correlation of residuals (detrended series)
r_resid = np.corrcoef(resid(emp, Z), resid(unemp, Z))[0, 1]

# 2. the one-variable recursion from ordinary correlations
R = np.corrcoef(np.vstack([emp, unemp, year]))
r12, r13, r23 = R[0, 1], R[0, 2], R[1, 2]
r_formula = (r12 - r13 * r23) / np.sqrt((1 - r13**2) * (1 - r23**2))

# 3. from the inverse of the correlation (or covariance) matrix
W = np.linalg.inv(R)
r_precision = -W[0, 1] / np.sqrt(W[0, 0] * W[1, 1])

# 4. from the t statistic of unemployment in the regression of emp on (1, unemp, year)
fit = sm.OLS(emp, np.column_stack([Z, unemp])).fit()
t = fit.tvalues[2]
r_t = np.sign(t) * np.sqrt(t**2 / (t**2 + n - 3))

print(f"marginal r = {r12:.4f}")
print(f"partial r:  {r_resid:.10f} {r_formula:.10f} {r_precision:.10f} {r_t:.10f}")
```

::: {.remark}
[Adjusting is not always right]

A partial correlation answers a precise question: how are the parts of two variables not
linearly explained by \( \Y_2 \) related? Whether that is the *right* question depends on
what \( \Y_2 \) is. Adjusting for a common cause, as the year stands in for population growth
here, removes a spurious association. Adjusting for a common effect can *create* an
association between variables that are independent. Chapter 25 treats
the difference. Under normality, the sampling distribution of a sample partial correlation
given \( q \) variables is that of an ordinary sample correlation from \( n-q \) observations, a
result used for tests in [Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html).
:::

## Exercises

### A. Check your understanding

::: {#exr-mvn-partial-zero}
[A1]

Show that \( \rho_{12\cdot3}=0 \) iff \( \rho_{12}=\rho_{13}\rho_{23} \). Find \( \rho_{12\cdot3} \) when
\( \rho_{12}=0 \) and \( \rho_{13}=\rho_{23}=0.6 \), and interpret.
:::

### B. Practice

::: {#exr-mvn-collider}
[B1]

Let \( X \), \( Z \) and \( e \) be independent, with \( X,Z\sim\Normal(0,1) \) and \( e\sim\Normal(0,\tau^2) \), and
let \( C=X+Z+e \). Show that \( \rho_{XZ}=0 \) but \( \rho_{XZ\cdot C}=-1/(1+\tau^2) \). What happens as
\( \tau\to0 \) and as \( \tau\to\infty \)?
:::

::: {.solution}
The covariance matrix of \( (X,Z,C) \) is
\( \begin{psmallmatrix}1&0&1\\0&1&1\\1&1&2+\tau^2\end{psmallmatrix} \). Conditioning on \( C \),
\[
\sigma_{XZ\cdot C}=0-\frac{1}{2+\tau^2},\qquad
\sigma_{XX\cdot C}=\sigma_{ZZ\cdot C}=1-\frac1{2+\tau^2}=\frac{1+\tau^2}{2+\tau^2},
\]
so \( \rho_{XZ\cdot C}=-1/(1+\tau^2) \). As \( \tau\to0 \) it tends to \( -1 \): if \( X+Z \) is known exactly,
a larger \( X \) forces a smaller \( Z \). As \( \tau\to\infty \) it tends to \( 0 \), since \( C \) carries little
information about \( X+Z \). Two independent causes become dependent once their common effect is
held fixed.
:::

::: {#exr-mvn-partial-product}
[B2]

Let \( \bSigma \) be positive definite. Show that the coefficient of \( Y_j \) in the best linear
predictor of \( Y_1 \) from all the other entries is \( -\omega_{1j}/\omega_{11} \). Deduce that
\( \rho_{1j\cdot\text{rest}}^2 \) equals the product of the coefficient of \( Y_j \) for predicting \( Y_1 \)
and the coefficient of \( Y_1 \) for predicting \( Y_j \), each from all the other variables.
:::

::: {.solution}
Put \( Y_1 \) first. By
@exr-mvn-precision-conditional with a scalar first block, the coefficient vector of the
best linear predictor of \( Y_1 \) from the rest is \( -\bm\Omega_{12}/\omega_{11} \), whose \( j \)th entry is
\( -\omega_{1j}/\omega_{11} \). Likewise the coefficient of \( Y_1 \) in the predictor of \( Y_j \) is
\( -\omega_{j1}/\omega_{jj} \). Their product is \( \omega_{1j}^2/(\omega_{11}\omega_{jj}) \), which is
\( \rho_{1j\cdot\text{rest}}^2 \) by @prp-mvn-precision.
:::

::: {#exr-mvn-sample-precision}
[B3]

Let \( \bm D \) be an \( n\times k \) data matrix whose columns are centred, and let \( \bS=\bm D\T\bm D \).
Show, using @thm-proj-fwl and the partitioned inverse, that
\( -w_{12}/\sqrt{w_{11}w_{22}} \), computed from \( \bS^{-1}=(w_{ij}) \), equals the correlation of the
residuals of columns \( 1 \) and \( 2 \) after regression on the remaining columns.
:::

::: {#exr-mvn-longley-gnp}
[B4]

Using the Longley data of @exm-mvn-longley, compute the partial correlation of employment and
unemployment given both year and GNP, in at least two of the four ways. Draw the corresponding
residual scatterplot. Does adjusting for GNP as well as year change the conclusion?
:::
