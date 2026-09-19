# Partial correlation

A partial correlation measures linear association after the linear effects of conditioning
variables have been removed from both variables (@def-mvn-partial-correlation, @prp-mvn-partial-meaning).
[Section 3.6](../ch03-multivariate-normal/06-partial-correlation.html) computed the sample version
and @thm-ss-partial-r2 identified its square with a partial \( R^2 \); the sampling distribution is still
missing. This section shows that a partial correlation from \( n \) cases with
\( q \) conditioning variables behaves exactly like an ordinary correlation from \( n-q \) cases, and that its test
is the \( t \) test of a regression coefficient.

## The sample partial correlation

::: {#def-cor-sample-partial}
[Sample partial correlation]

Let \( \y,\x\in\Real^n \), let \( \Z \) be an \( n\times q \) matrix, and let \( \bm D=[\bone,\Z] \) have rank \( q+1<n \). With
\( \M_D \) the projection onto \( \C(\bm D) \), put \( \tilde{\y}=(\I-\M_D)\y \) and \( \tilde{\x}=(\I-\M_D)\x \), and suppose both are
nonzero. The **sample partial correlation** of \( \y \) and \( \x \) given \( \Z \) is
\[
r_{yx\cdot Z}=\frac{\tilde{\y}\T\tilde{\x}}{\norm{\tilde{\y}}\,\norm{\tilde{\x}}} .
\]
:::

Both residual vectors have mean zero, because \( \bone\in\C(\bm D) \), so \( r_{yx\cdot Z} \) is the ordinary sample
correlation of the two residual vectors. With \( q=0 \) it is the ordinary correlation of \( \y \) and \( \x \).

## Its distribution

The key is that residualizing on \( q+1 \) columns leaves \( n-q-1 \) dimensions, exactly as centring
alone leaves \( n-1 \).

::: {#thm-cor-partial-law}
[Distribution of a sample partial correlation]

Let \( \Z \) be a random \( n\times q \) matrix with \( \rank[\bone,\Z]=q+1 \) with probability one, and \( n\ge q+3 \).
Suppose that, given \( \Z \), the rows \( (y_i,x_i) \) are independent bivariate normal vectors with means
\( \bm a+\B\bz_{(i)} \), affine in the \( i \)th row \( \bz_{(i)} \) of \( \Z \), and a common covariance matrix with correlation \( \rho \).
Then \( r_{yx\cdot Z} \) has the same distribution as the sample correlation of \( n-q \) independent bivariate
normal pairs with correlation \( \rho \), and it is independent of \( \Z \).

In particular, if the rows \( (y_i,x_i,\bz_{(i)}\T) \) are independent draws from a multivariate normal distribution
with a nonsingular covariance matrix, the conclusion holds with \( \rho=\rho_{yx\cdot Z} \), the population partial
correlation.
:::

::: {.proof}
Fix \( \Z \) and write \( \bm D=[\bone,\Z] \). Let \( \bm K \) be an \( (n-q-1)\times n \) matrix whose rows form an orthonormal
basis of \( \C(\bm D)\perpc \), so that \( \bm K\bm D=\bzero \) and \( \bm K\T\bm K=\I-\M_D \). Then
\( \tilde{\y}\T\tilde{\x}=\y\T(\I-\M_D)\x=(\bm K\y)\T(\bm K\x) \), and similarly for the squared lengths, so
\[
r_{yx\cdot Z}=\frac{(\bm K\y)\T(\bm K\x)}{\norm{\bm K\y}\,\norm{\bm K\x}} .
\]
Stack the pairs as the \( n\times2 \) matrix \( [\y,\x] \). Given \( \Z \), \( [\y,\x]=\bm D\bm A+\bm E \) for a \( (q+1)\times2 \) matrix of
coefficients \( \bm A \), where \( \bm E \) has independent \( \Normal_2(\bzero,\bm C) \) rows. Hence \( \bm K[\y,\x]=\bm K\bm E \). Writing
\( \bm E=\bm E_0\bm C^{1/2} \) with independent standard normal entries in \( \bm E_0 \), each column of \( \bm K\bm E_0 \) is
\( \Normal(\bzero,\bm K\bm K\T)=\Normal(\bzero,\I_{n-q-1}) \) (@thm-mvn-linear), so the \( n-q-1 \) rows of \( \bm K\bm E \) are independent
\( \Normal_2(\bzero,\bm C) \) vectors. Thus, given \( \Z \), \( r_{yx\cdot Z} \) is the *uncentred* correlation of \( m=n-q-1 \)
independent \( \Normal_2(\bzero,\bm C) \) pairs.

The same argument with \( q=0 \) and \( \Z \) absent shows that the ordinary sample correlation of \( m+1 \)
bivariate normal pairs is the uncentred correlation of \( m \) independent \( \Normal_2(\bzero,\bm C) \) pairs. So the
conditional law of \( r_{yx\cdot Z} \) is that of an ordinary correlation from \( m+1=n-q \) pairs, whatever \( \Z \) is,
and @lem-cor-conditioning gives the unconditional statement.

For the multivariate normal case, @thm-mvn-conditional says that given the conditioning variables the
pair is bivariate normal with mean affine in them and covariance \( \bSigma_{11\cdot2} \), whose correlation is
\( \rho_{yx\cdot Z} \) by @def-mvn-partial-correlation.
:::

This is the result promised in [Section 3.6](../ch03-multivariate-normal/06-partial-correlation.html):
each conditioning variable costs one case. Everything in [Section 14.4](04-correlation-coefficient.html)
transfers with \( n \) replaced by \( n-q \): the exact law @eq-cor-r-cdf, the \( t \) test with \( n-q-2 \) degrees of
freedom, and Fisher's interval with \( \sqrt{n-q-3} \) (Fisher 1924). The conditioning variables need not be normal.

::: {#exm-cor-partial-simulation}
[Twelve cases, two skewed conditioning variables]

Take \( n=12 \), \( q=2 \) conditioning variables drawn as squares of exponential variables, which are
very skewed, a pair whose conditional means depend linearly on them, and conditional correlation
\( \rho=0.6 \). Over \( 100000 \) samples, the \( 10\% \), \( 50\% \) and \( 90\% \) points of the sample partial
correlation are \( 0.257 \), \( 0.624 \) and \( 0.834 \). For ordinary correlations from
\( n-q=10 \) normal pairs they are \( 0.256 \), \( 0.625 \) and \( 0.835 \). A two-sample
Kolmogorov–Smirnov test finds no difference between the two sets of simulated values (\( p \)-value
\( 0.09 \)).
:::

```{.python .run #cell-partial-correlation-simulate}
import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

rng = np.random.default_rng(1451)
n, q, rho, reps = 12, 2, 0.6, 100_000
B = np.array([[1.0, -2.0], [0.5, 3.0]])              # dependence of the pair on Z
C = np.array([[1.0, rho], [rho, 1.0]])               # conditional covariance of the pair
Lc = np.linalg.cholesky(C)

partial = np.empty(reps)
for s in range(reps):
    Zc = rng.exponential(size=(n, q)) ** 2           # skewed conditioning variables
    pair = Zc @ B.T + rng.normal(size=(n, 2)) @ Lc.T
    D = np.column_stack([np.ones(n), Zc])
    coef, *_ = np.linalg.lstsq(D, pair, rcond=None)
    e = pair - D @ coef                              # residuals of both variables
    partial[s] = e[:, 0] @ e[:, 1] / np.sqrt((e[:, 0] @ e[:, 0]) * (e[:, 1] @ e[:, 1]))

# ordinary correlations from n - q pairs with correlation rho
x = rng.normal(size=(reps, n - q))
y = rho * x + np.sqrt(1 - rho ** 2) * rng.normal(size=(reps, n - q))
ordinary = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(x, y)])
print("quantiles, partial :", np.quantile(partial, [0.1, 0.5, 0.9]).round(3))
print("quantiles, ordinary:", np.quantile(ordinary, [0.1, 0.5, 0.9]).round(3))
print("two-sample KS p-value:", round(stats.ks_2samp(partial, ordinary).pvalue, 3))
```

## The test of a partial correlation is a coefficient test

In the population, a zero partial correlation and a zero regression coefficient are the same thing.

::: {#prp-cor-partial-coefficient}
[Partial correlation and the regression coefficient]

Let \( (Y,X,\bm Z) \) have finite second moments, with \( \Cov\bigl((\bm Z\T,X)\T\bigr) \) positive definite, and let
\( \beta_X \) be the coefficient of \( X \) in the best linear predictor of \( Y \) from \( (\bm Z,X) \). Then
\[
\beta_X=\frac{\sigma_{yx\cdot Z}}{\sigma_{xx\cdot Z}}=\rho_{yx\cdot Z}\sqrt{\frac{\sigma_{yy\cdot Z}}{\sigma_{xx\cdot Z}}} ,
\]
so \( \beta_X=0 \) iff \( \rho_{yx\cdot Z}=0 \).
:::

::: {.proof}
Let \( e_X=X-L(X\mid\bm Z) \) and \( e_Y=Y-L(Y\mid\bm Z) \) be the prediction errors of @prp-mvn-partial-meaning(a),
with \( \Var(e_X)=\sigma_{xx\cdot Z}>0 \) and \( \Cov(e_Y,e_X)=\sigma_{yx\cdot Z} \). The affine span of \( (\bm Z,X) \) is the orthogonal
sum, in \( L^2 \), of the affine span of \( \bm Z \) and the line spanned by \( e_X \), which is orthogonal to it. By the
population version of @thm-proj-sum, projecting \( Y \) onto the sum projects it onto each piece:
\( L(Y\mid\bm Z,X)=L(Y\mid\bm Z)+\bigl(\Cov(Y,e_X)/\Var(e_X)\bigr)e_X \). Since \( L(Y\mid\bm Z) \) is orthogonal to \( e_X \),
\( \Cov(Y,e_X)=\Cov(e_Y,e_X)=\sigma_{yx\cdot Z} \). The coefficient of \( X \) in this expression comes only from the last
term, because \( e_X=X-(\text{affine function of }\bm Z) \); it is \( \sigma_{yx\cdot Z}/\sigma_{xx\cdot Z} \).
:::

This is the Frisch–Waugh–Lovell theorem in the population, and the sample version is
@thm-proj-fwl itself. Combining the two with @thm-cor-partial-law:

::: {#thm-cor-partial-test}
[The partial correlation test is the coefficient \( t \) test]

Regress \( \y \) on the \( p=q+2 \) columns \( [\bone,\Z,\x] \), assumed of full rank, and let \( t \) be the \( t \) statistic of
the coefficient of \( \x \) for the value zero. Then:

::: {.enumerate options="label=(\alph*)"}
1. for every data set, \( t=\sqrt{n-p}\;r_{yx\cdot Z}\big/\sqrt{1-r_{yx\cdot Z}^2} \);

2. under the assumptions of @thm-cor-partial-law with \( \rho=0 \), \( t\sim t(n-p) \), whatever the distribution of
   \( \Z \);

3. when the rows are multivariate normal, the \( t \) test of the coefficient of \( \x \) and the test of
   \( \rho_{yx\cdot Z}=0 \) are the same test, and it has exact size.
:::

:::

::: {.proof}
(a) By @thm-ss-partial-r2(c) (see also @exr-proj-partial-correlation), \( r_{yx\cdot Z}^2=t^2/(t^2+n-p) \). By @thm-proj-fwl the coefficient of \( \x \) is
\( \tilde{\x}\T\tilde{\y}/\norm{\tilde{\x}}^2 \), which has the sign of \( r_{yx\cdot Z} \), so \( t \) and \( r_{yx\cdot Z} \) have the same sign.
Solving for \( t \) gives (a).

(b) By @thm-cor-partial-law, \( r_{yx\cdot Z} \) has the law of an ordinary correlation from \( n-q \) pairs with \( \rho=0 \),
and by @prp-cor-r-exact(c) the statistic in (a) is then \( t\bigl((n-q)-2\bigr)=t(n-p) \).

(c) Under multivariate normality the best linear predictor (@thm-rv-blp) is the conditional mean (@prp-mvn-best-predictor), so the regression coefficient of \( \x \) in the conditional model is \( \beta_X \) of
@prp-cor-partial-coefficient, which vanishes iff \( \rho_{yx\cdot Z}=0 \). By (a) the two test statistics are
monotone functions of each other, and by (b) or by @thm-cor-conditional the null distribution is exact.
:::

So the test of a single coefficient in [Chapter 11](../ch11-general-linear-hypothesis/index.html) (@thm-glh-t-test) and the test of a partial correlation are one and the same, under either reading
of the regressors. What the partial correlation adds is a scale-free measure of effect size, with its
own confidence interval, which a coefficient does not provide.

::: {#exm-cor-states-partial}
[Violent crime and urbanization]

Across the \( 50 \) states, the violent crime rate and the urbanization rate have correlation
\( 0.345 \). Given the percentage of single-parent households and the poverty rate (\( q=2 \)), their
partial correlation is \( r=0.1863 \). In the regression of violent crime on the three variables and an
intercept, the coefficient of urbanization is \( 1.388 \) with \( t=1.286 \) on \( 46 \) degrees of
freedom (\( p \)-value \( 0.2049 \)), and indeed \( \sqrt{46}\,r/\sqrt{1-r^2}=1.286 \). With \( n-q=48 \), the
Fisher interval for \( \rho_{yx\cdot Z} \) is \( (-0.103, 0.447) \) and the exact interval is
\( (-0.103, 0.443) \). Much of the marginal association is shared with single parenthood and poverty, and what remains
is compatible with none.
:::

```{.python .run #cell-partial-correlation-states}
def r_cdf(c, rho, n, m=4000):
    """Exact P(r <= c) for an ordinary sample correlation of n bivariate normal pairs."""
    if c < 0:
        return 1 - r_cdf(-c, -rho, n, m)
    t = c * np.sqrt(n - 2) / np.sqrt(1 - c ** 2)
    if rho == 0:
        return stats.t.cdf(t, n - 2)
    V = stats.chi2.ppf((np.arange(m) + 0.5) / m, n - 1)
    return np.mean(stats.nct.cdf(t, n - 2, rho / np.sqrt(1 - rho ** 2) * np.sqrt(V)))

data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["violent"].to_numpy()
x = data["urban"].to_numpy()
Z = np.column_stack([np.ones(len(y)), data["single"], data["poverty"]])
n_s, q_s = len(y), Z.shape[1] - 1                    # q = 2 conditioning variables

def resid(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef

ey, ex = resid(y, Z), resid(x, Z)
r_p = ey @ ex / np.sqrt((ey @ ey) * (ex @ ex))       # sample partial correlation
fit = sm.OLS(y, np.column_stack([Z, x])).fit()       # p = 4 columns
t = fit.tvalues[-1]
df = n_s - 4
print(f"partial r = {r_p:.4f},  t = {t:.3f},  sqrt(df) r/sqrt(1-r^2) = "
      f"{np.sqrt(df) * r_p / np.sqrt(1 - r_p**2):.3f}")
zq = stats.norm.ppf(0.975)
lo_z, hi_z = np.tanh(np.arctanh(r_p) + np.array([-1, 1]) * zq / np.sqrt(n_s - q_s - 3))
lo_x = optimize.brentq(lambda rh: r_cdf(r_p, rh, n_s - q_s) - 0.975, -0.99, 0.99)
hi_x = optimize.brentq(lambda rh: r_cdf(r_p, rh, n_s - q_s) - 0.025, -0.99, 0.99)
print(f"Fisher z interval ({lo_z:.3f}, {hi_z:.3f}); exact ({lo_x:.3f}, {hi_x:.3f})")
```

::: {.remark}
[Which variables to condition on]

Conditioning can shrink an association, as in the example, reverse its sign, as in @exm-mvn-longley
and Simpson's paradox, or create one where none exists, when the conditioning variable is a common
effect of the two (@exr-cor-collider). The distribution theory is silent on the choice of conditioning
set; Chapter 25 takes it up.
:::

## Exercises

### A. Check your understanding

::: {#exr-cor-partial-numbers}
[A1]

From \( n=40 \) cases, the partial correlation of \( y \) and \( x \) given \( q=3 \) other variables is \( -0.30 \). Compute
the \( t \) statistic, its degrees of freedom, and the \( 95\% \) Fisher interval for the population partial
correlation.
:::

::: {#exr-cor-partial-df}
[A2]

Explain, in terms of dimensions of subspaces, why the \( t \) test of a partial correlation given \( q \) variables
has \( n-q-2 \) degrees of freedom, and why this equals \( n-p \) for the corresponding regression.
:::

### B. Practice

::: {#exr-cor-partial-r2-beta}
[B1]

Let \( \X=[\X_1,\X_2] \) with \( \bone\in\C(\X_1) \), \( \rank\X_1=r_1 \), \( \rank\X=r \) almost surely, and suppose that, given the regressors, \( \Y \)
follows the normal linear model with mean \( \X_1\bbeta_1 \). Show that the partial
\( R^2 \) of @def-ss-partial-r2 has the \( \mathrm{Beta}\bigl(\tfrac{r-r_1}2,\tfrac{n-r}2\bigr) \) distribution, whatever the
distribution of the regressors.
:::

::: {.solution}
By @thm-ss-partial-r2(b), \( R^2_{\X\mid\X_1}=qF/(qF+n-r) \) with \( q=r-r_1 \), where \( F \) is the \( F \) ratio for the
added block. Given the regressors, the null hypothesis holds, so \( F\sim F(q,n-r) \) for every design (@thm-cor-conditional(b)), and by @lem-cor-conditioning unconditionally. Writing \( F=(U/q)/(V/(n-r)) \) with independent
chi-squared variables, \( R^2_{\X\mid\X_1}=U/(U+V) \), which is \( \mathrm{Beta}(q/2,(n-r)/2) \) by @exr-qf-beta.
:::

::: {#exr-cor-semipartial}
[B2]

The *semipartial* correlation of \( \y \) and \( \x \) given \( \Z \) is the correlation between \( \y \) itself and the
residual \( \tilde{\x}=(\I-\M_D)\x \). Show that its square equals \( R^2_{\text{full}}-R^2_{\text{reduced}} \), the increase in
\( R^2 \) when \( \x \) is added to \( [\bone,\Z] \). How does it relate to the squared partial correlation?
:::

::: {.solution}
By @lem-ss-residualized, the extra sum of squares of \( \x \) is \( (\tilde{\x}\T\y)^2/\norm{\tilde{\x}}^2 \). Dividing by
\( \text{SST}=\norm{(\I-\bP_1)\y}^2 \) gives \( R^2_{\text{full}}-R^2_{\text{reduced}} \). Since \( \tilde{\x} \) has mean zero,
\( \tilde{\x}\T\y=\tilde{\x}\T(\I-\bP_1)\y \), so the quotient is the squared correlation of \( \y \) and \( \tilde{\x} \). The squared partial
correlation divides the same extra sum of squares by \( \text{SSE}_{\text{reduced}} \) instead of \( \text{SST} \), so it
equals the squared semipartial correlation divided by \( 1-R^2_{\text{reduced}} \), and is never smaller.
:::

### C. Going deeper

::: {#exr-cor-collider}
[C1]

Let \( X \) and \( Y \) be independent \( \Normal(0,1) \) and \( C=X+Y+U \) with \( U\sim\Normal(0,1) \) independent of both.
Compute \( \rho_{XY\cdot C} \). With \( n=200 \) cases, what does the \( t \) test of the coefficient of \( X \) in the
regression of \( Y \) on \( C \) and \( X \) tend to conclude? Is the conclusion about the relation between \( X \) and \( Y \)
wrong, or is the question wrong?
:::

