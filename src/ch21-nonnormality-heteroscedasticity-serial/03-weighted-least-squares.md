# Weighted least squares

Once heteroscedasticity has been found, one can transform the response
([Chapter 22](../ch22-transformations/index.html)), keep least squares and correct its standard errors
([Section 21.4](04-sandwich-estimators.html)), or give the noisier observations less weight. This section is
about weighting. With known variances it is optimal and its theory is already in place. With estimated
variances the gain is real but the theory is asymptotic, and a careless implementation can do worse than no
weighting at all.

## Known weights

Suppose \( \Var(\varepsilon_i)=\sigma^2/w_i \), with the \( w_i>0 \) known and \( \sigma^2 \) unknown, and the errors
uncorrelated. With \( \mathbf{W}=\diag(w_1,\dots,w_n) \) this is \( \Cov(\be)=\sigma^2\mathbf{W}^{-1} \). The generalized least
squares estimate of @def-proj-gls with \( \V=\mathbf{W}^{-1} \) minimizes
\[
\sum_{i=1}^nw_i\,(y_i-\x_{(i)}\T\bb)^2 ,
\qquad\text{so}\qquad
\hbeta_W=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\y .
\]{#eq-het-wls}

This is the **weighted least squares** (WLS) estimate. [Section 6.9](../ch06-projections/09-inner-products.html)
described its geometry: it is ordinary least squares after each row of \( \X \) and \( \y \) is multiplied by
\( \sqrt{w_i} \), which turns the errors into \( \sqrt{w_i}\,\varepsilon_i \) with common variance \( \sigma^2 \). Aitken's
theorem (@cor-opt-aitken) makes it the best linear unbiased estimator. Everything in Part III then applies
to the rescaled data. The error variance is estimated by
\[
s_w^2=\frac{1}{n-p}\sum_{i=1}^nw_i\,(y_i-\x_{(i)}\T\hbeta_W)^2 ,
\]
the covariance of \( \hbeta_W \) by \( s_w^2(\X\T\mathbf{W}\X)^{-1} \), and under normality the \( t \) and \( F \) procedures are
exact. The residuals to plot are the rescaled ones, \( \sqrt{w_i}\,(y_i-\x_{(i)}\T\hbeta_W) \), studentized with the
leverages of the weighted fit. The raw residuals \( y_i-\x_{(i)}\T\hbeta_W \) should still fan out, since their
variances still differ.

Known weights arise when \( y_i \) is the mean of \( m_i \) observations (\( w_i=m_i \)), a total (\( w_i=1/m_i \)), or a
measurement from an instrument of known precision. The weights are then fixed by the design, not by the data.

## Estimating the weights

More often the variances are unknown, and the weights must come from the data. There are three common routes.

1. **Replicates.** If there are several observations at each distinct row of \( \X \), or in each of a few
   groups known to differ in precision, the sample variance within each group estimates its \( \sigma_g^2 \).
2. **A variance function in other variables.** Model \( \sigma_i^2=\sigma^2\exp(\bz_i\T\boldsymbol{\upgamma}) \), as in
   @eq-het-variance-function, and estimate \( \boldsymbol{\upgamma} \) by regressing \( \log\hat{\varepsilon}_i^2 \) on an intercept
   and \( \bz_i \). The intercept of that regression is biased, because \( \E\log\chi^2(1)\ne0 \), but weights are needed
   only up to a common factor, so the intercept does not matter.
3. **A variance function of the mean.** Model \( \sigma_i=\sigma\lvert\mu_i\rvert^\theta \) with
   \( \mu_i=\x_{(i)}\T\bbeta \), estimate \( \theta \) from the slope of \( \log\lvert\hat{\varepsilon}_i\rvert \) on
   \( \log\lvert\hat Y_i\rvert \), and weight by \( \lvert\hat Y_i\rvert^{-2\theta} \). Since the weights depend on the
   fit, the estimate is usually iterated.

The resulting estimate \( \hbeta_{\hat W} \) is called **feasible** (or estimated) weighted least squares.

::: {#exm-het-engel-fwls}
[Engel's data, reweighted]

Regressing \( \log\hat{\varepsilon}_i^2 \) on an intercept and \( \log\text{income}_i \) gives slope
\( \hat\gamma=1.90 \), so the fitted variance is proportional to \( \text{income}^{1.90} \), close to the square law
suggested by [Figure 21.2.1](02-detecting-heteroscedasticity.html#fig-het-engel). Weighting by
\( \text{income}^{-1.90} \) changes the slope from \( 0.4852 \) to \( 0.5711 \). The model-based standard error
\( s_w\{(\X\T\hat{\mathbf{W}}\X)^{-1}\}_{22}^{1/2} \) of the weighted slope is \( 0.0149 \). The sandwich standard error of
[Section 21.4](04-sandwich-estimators.html), which does not trust the variance model, is \( 0.0172 \) for the
weighted slope and \( 0.0664 \) for the unweighted one. So weighting has cut the standard error by a factor of
nearly four.

The variance function is not sharply determined by these data. Re-estimating \( \gamma \) from the residuals of
the weighted fit, and repeating until it settles (\( 10 \) iterations), moves \( \hat\gamma \) to \( 2.52 \) and the
slope to \( 0.5875 \), about one standard error away. Reporting the weighted estimate with the sandwich standard
error protects against this uncertainty in the weights.
:::

```{.python .run #cell-feasible-wls-engel}
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = sm.datasets.engel.load_pandas().data
income, food = df["income"].to_numpy(), df["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])


def wls(X, y, w):
    """Weighted least squares: coefficients, model-based and HC3 covariances."""
    sw = np.sqrt(w)
    Q, R = np.linalg.qr(X * sw[:, None])              # OLS on the whitened data
    b = np.linalg.solve(R, Q.T @ (y * sw))
    ew = (y - X @ b) * sw                             # whitened residuals
    # leverages in the weighted geometry
    hw = np.sum(Q**2, axis=1)
    Rinv = np.linalg.inv(R)
    bread = Rinv @ Rinv.T                             # (X'WX)^{-1}
    s2w = ew @ ew / (len(y) - X.shape[1])
    meat = (Q * (ew / (1 - hw))[:, None]).T @ (Q * (ew / (1 - hw))[:, None])
    return b, s2w * bread, Rinv @ meat @ Rinv.T


b_ols, V_ols, V_ols_hc3 = wls(X, food, np.ones(n))
e = food - X @ b_ols
Z = np.column_stack([np.ones(n), np.log(income)])
# variance ~ income^gamma
gamma = np.linalg.lstsq(Z, np.log(e**2), rcond=None)[0][1]
b_fw, V_fw, V_fw_hc3 = wls(X, food, income**-gamma)
print(f"estimated power gamma = {gamma:.2f}")
print(f"OLS  slope {b_ols[1]:.4f}  model se {np.sqrt(V_ols[1, 1]):.4f}  "
      f"HC3 se {np.sqrt(V_ols_hc3[1, 1]):.4f}")
print(f"FWLS slope {b_fw[1]:.4f}  model se {np.sqrt(V_fw[1, 1]):.4f}  "
      f"HC3 se {np.sqrt(V_fw_hc3[1, 1]):.4f}")
```

## What feasible weighting achieves

The weights in \( \hbeta_{\hat W} \) are random and depend on the same data as the fit, so none of the finite-sample
results for known weights carries over directly. The next theorem collects what can be said in general. Part (b)
shows that weighting with the wrong weights is harmless for bias. Part (c) is an exact finite-sample result,
and part (d) is the large-sample justification in the setting of replicated groups.

::: {#thm-het-wls}
[Weighted least squares with known and estimated weights]

Let \( \Y=\X\bbeta+\be \) with \( \E\be=\bzero \) and \( \rank(\X)=p \).

::: {.enumerate options="label=(\alph*)"}
1. *(Known weights.)* If \( \Cov(\be)=\sigma^2\mathbf{W}^{-1} \) with \( \mathbf{W} \) known, then \( \hbeta_W \) is the best linear
   unbiased estimator of \( \bbeta \), \( \Cov(\hbeta_W)=\sigma^2(\X\T\mathbf{W}\X)^{-1} \) and \( \E s_w^2=\sigma^2 \). If \( \be \) is
   normal, \( \hbeta_W\sim\Normal_p(\bbeta,\sigma^2(\X\T\mathbf{W}\X)^{-1}) \) independently of
   \( (n-p)s_w^2/\sigma^2\sim\chi^2(n-p) \).

2. *(Wrong weights.)* If \( \Cov(\be)=\bSigma \) is arbitrary and \( \mathbf{W} \) is any fixed positive diagonal matrix, then
   \( \hbeta_W \) is unbiased with
   \[
\Cov(\hbeta_W)=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\bSigma\mathbf{W}\X(\X\T\mathbf{W}\X)^{-1}.
\]

3. *(Estimated weights, symmetric errors.)* Let \( \hat{\mathbf{W}}=\mathbf{W}(\he) \) be a positive diagonal matrix that
   depends on the data only through the ordinary least squares residuals, with \( \mathbf{W}(-\mathbf{e})=\mathbf{W}(\mathbf{e}) \). If
   \( -\be \) has the same distribution as \( \be \) and \( \E\norm{\hbeta_{\hat W}}<\infty \), then
   \( \E\hbeta_{\hat W}=\bbeta \).

4. *(Estimated weights from groups.)* Suppose the observations form \( G \) groups, fixed in advance, that the
   errors are independent, and that within group \( g \) they are identically distributed with variance
   \( \sigma_g^2>0 \). Let \( n_g\to\infty \) for each \( g \) with \( n^{-1}\X_g\T\X_g\to\A_g \), where \( \X_g \) holds the rows
   of group \( g \), and let \( \A_w=\sum_g\A_g/\sigma_g^2 \) be positive definite. Let \( \mathbf{W} \) have entries
   \( 1/\sigma_g^2 \), and \( \hat{\mathbf{W}} \) have entries \( 1/\hat{\sigma}_g^2 \), where \( \hat{\sigma}_g^2=\norm{\he_g}^2/n_g \) is the
   mean squared ordinary least squares residual in group \( g \). Then
   \[
\sqrt n\,(\hbeta_{\hat W}-\hbeta_W)\to\bzero\quad\text{and}\quad n(\X\T\hat{\mathbf{W}}\X)^{-1}\to\A_w^{-1}
\]
   in probability,
   while \( n\Cov(\hbeta_W)=n(\X\T\mathbf{W}\X)^{-1}\to\A_w^{-1} \).
:::

:::

::: {.proof}
(a) Multiply @eq-het-wls through by \( \mathbf{W}^{1/2} \): the model \( \mathbf{W}^{1/2}\Y=\mathbf{W}^{1/2}\X\bbeta+\mathbf{W}^{1/2}\be \)
satisfies \( \Cov(\mathbf{W}^{1/2}\be)=\sigma^2\I \), and \( \hbeta_W \) is its least squares estimate. The claims are
@cor-opt-aitken, @thm-lm-moments, @thm-lm-sigma2 and @thm-opt-sampling applied to the rescaled model.

(b) \( \hbeta_W-\bbeta=(\X\T\mathbf{W}\X)^{-1}\X\T\mathbf{W}\be \) is a fixed linear map of \( \be \); apply @thm-rv-linear.

(c) The residuals \( \he=(\I-\M)\be \) do not involve \( \bbeta \), and they change sign with \( \be \). So
\( \hbeta_{\hat W}-\bbeta=\mathbf{g}(\be) \) with \( \mathbf{g}(\be)=(\X\T\mathbf{W}(\he)\X)^{-1}\X\T\mathbf{W}(\he)\be \), and
\( \mathbf{g}(-\be)=-\mathbf{g}(\be) \) because \( \mathbf{W}(-\he)=\mathbf{W}(\he) \). If \( -\be \) has the law of \( \be \), then
\( \mathbf{g}(\be) \) has the law of \( \mathbf{g}(-\be)=-\mathbf{g}(\be) \), and its mean, which exists, equals its negative.

(d) *Consistency of the group variances.* Write \( \mathbf{a}=\M\be \), so \( \he_g=\be_g-\mathbf{a}_g \). By the weak law,
\( \norm{\be_g}^2/n_g\to\sigma_g^2 \). Also \( \E\norm{\mathbf{a}}^2=\tr(\M\bSigma)\le p\max_g\sigma_g^2 \), so
\( \norm{\mathbf{a}_g}^2/n_g\le\norm{\mathbf{a}}^2/n_g\to0 \) in probability. The triangle inequality
\( \bigl\lvert\norm{\he_g}-\norm{\be_g}\bigr\rvert\le\norm{\mathbf{a}_g} \) gives \( \hat{\sigma}_g^2\to\sigma_g^2 \), hence
\( \hat w_g=1/\hat{\sigma}_g^2\to w_g=1/\sigma_g^2 \).

*Comparison.* Put \( \A_n(\mathbf{v})=n^{-1}\sum_gv_g\X_g\T\X_g \) for weights \( \mathbf{v}=(v_1,\dots,v_G) \), and
\( \bb_g=n^{-1/2}\X_g\T\be_g \). Then
\[
\begin{aligned}
\sqrt n(\hbeta_{\hat W}-\bbeta)&=\A_n(\hat{\mathbf{w}})^{-1}\sum_g\hat w_g\bb_g,\\
\sqrt n(\hbeta_W-\bbeta)&=\A_n(\mathbf{w})^{-1}\sum_gw_g\bb_g .
\end{aligned}
\]
Each \( \bb_g \) has mean zero and \( \E\norm{\bb_g}^2=\sigma_g^2\tr(n^{-1}\X_g\T\X_g) \), which is bounded, so
\( \bb_g=O_p(1) \). Both \( \A_n(\hat{\mathbf{w}}) \) and \( \A_n(\mathbf{w}) \) converge in probability to \( \A_w \), which is
invertible, and matrix inversion is continuous there. The difference of the two displays is
\[
\bigl\{\A_n(\hat{\mathbf{w}})^{-1}-\A_n(\mathbf{w})^{-1}\bigr\}\sum_g\hat w_g\bb_g+\A_n(\mathbf{w})^{-1}\sum_g(\hat w_g-w_g)\bb_g,
\]
a sum of products of \( o_p(1) \) and \( O_p(1) \) terms, so it is \( o_p(1) \). Finally
\( n(\X\T\hat{\mathbf{W}}\X)^{-1}=\A_n(\hat{\mathbf{w}})^{-1}\to\A_w^{-1} \), and \( n(\X\T\mathbf{W}\X)^{-1}=\A_n(\mathbf{w})^{-1}\to\A_w^{-1} \).
:::

Part (c) is exact but covers only weights computed from the least squares residuals, not those computed from
fitted values (@exr-het-fitted-weights). Part (d) says that with enough replicates estimating the weights costs
nothing in large samples. The scaled errors \( \sqrt n(\hbeta_{\hat W}-\bbeta) \) and \( \sqrt n(\hbeta_W-\bbeta) \)
differ by \( o_p(1) \), so whatever limiting distribution the optimal estimator has, the feasible one has the
same; it is normal when, for instance, the rows of \( \X \) are bounded, by @lem-het-clt(a) applied within each
group. The model-based covariance estimate is consistent as well. Carroll and Ruppert (1988) prove the same
conclusion for smooth parametric variance functions, including functions of the mean, under conditions
we do not reproduce. What (d) does not say is how many replicates are "enough". A simulation shows that the
answer can be surprisingly many.

::: {#exm-het-replicates}
[Weights from a few replicates]

Take ten design points \( x=1,\dots,10 \), \( m \) replicates at each, a straight-line mean, and normal errors
with standard deviation proportional to \( x \). [Figure 21.3.1](#fig-het-fwls) compares four slope estimates by
their variance relative to weighted least squares with the true weights. Ordinary least squares has
\( 2.44 \) times the optimal variance when \( m=2 \), and about the same for every \( m \). Weights \( 1/s_g^2 \) from the
within-group sample variances are a disaster for small groups: with \( m=2 \) their variance is \( 5.44 \) times
the optimum, twice that of the unweighted fit, and with \( m=3 \) it is still \( 2.67 \), worse than no
weighting. They approach the optimum only slowly (\( 1.22 \) for \( m=10 \), \( 1.05 \) for \( m=30 \)). Worse, the
model-based interval, which treats the estimated weights as known, covers the true slope in only
\( 0.47 \) of the data sets for \( m=2 \), \( 0.65 \) for \( m=3 \) and \( 0.80 \) for \( m=5 \). Weights from a fitted variance
function, the regression of \( \log s_g^2 \) on \( \log x_g \), are far better: \( 1.41 \) times the optimum for
\( m=2 \), \( 1.12 \) for \( m=3 \) and \( 1.04 \) for \( m=5 \).
:::

::: {when-format="html"}
![**Figure 21.3.1.** Feasible weighted least squares with \( m \) replicates at each of ten design points and
standard deviation proportional to \( x \) (\( 20{,}000 \) simulated data sets for each \( m \)). (a) Variance of the
slope estimate relative to weighted least squares with the true weights. (b) Coverage of the nominal 95%
interval computed as if the within-group weights were known.](feasible_wls.svg){#fig-het-fwls width=100%}
:::

::: {when-format="pdf"}
![Feasible weighted least squares with \( m \) replicates at each of ten design points and
standard deviation proportional to \( x \) (\( 20{,}000 \) simulated data sets for each \( m \)). (a) Variance of the
slope estimate relative to weighted least squares with the true weights. (b) Coverage of the nominal 95%
interval computed as if the within-group weights were known.](feasible_wls.pdf){width=100%}
:::

```{.python .run #cell-feasible-wls-replicates}
rng = np.random.default_rng(2106)
# ten design points, sd proportional to x
xg = np.arange(1.0, 11.0)


def replicate_study(m, reps):
    """Variance of four slope estimators relative to WLS with true weights."""
    x = np.repeat(xg, m)
    Xr = np.column_stack([np.ones(len(x)), x])
    Y = 1 + 0.5 * x + rng.normal(size=(reps, len(x))) * x

    def slopes(W):
        A = np.einsum("ri,ij,ik->rjk", W, Xr, Xr)     # X'WX for each replicate
        c = np.einsum("ri,ij,ri->rj", W, Xr, Y)
        return np.linalg.solve(A, c[..., None])[..., 0], A

    b_o, _ = slopes(np.ones_like(Y))
    b_w, _ = slopes(np.tile(1 / x**2, (reps, 1)))
    # within-group variances
    s2g = Y.reshape(reps, len(xg), m).var(axis=2, ddof=1)
    Wg = np.repeat(1 / s2g, m, axis=1)
    b_g, A_g = slopes(Wg)
    # smooth model: log s2 on log x
    Zg = np.column_stack([np.ones(len(xg)), np.log(xg)])
    coef = np.linalg.lstsq(Zg, np.log(s2g).T, rcond=None)[0]
    b_m, _ = slopes(np.repeat(np.exp(-(Zg @ coef).T), m, axis=1))
    var = np.array([np.var(b[:, 1]) for b in (b_o, b_w, b_g, b_m)])
    # nominal interval for b_g
    res = Y - b_g @ Xr.T
    se = np.sqrt((Wg * res**2).sum(axis=1) / (len(x) - 2)
                 * np.linalg.inv(A_g)[:, 1, 1])
    q = stats.t.ppf(0.975, len(x) - 2)
    cover = np.mean(np.abs(b_g[:, 1] - 0.5) <= q * se)
    return var / var[1], cover


for m in (2, 3, 5, 10, 30):
    rel, cover = replicate_study(m, 2000)
    print(f"m = {m:2d}: variance relative to WLS  OLS {rel[0]:.2f}  "
          f"group weights {rel[2]:.2f}  fitted weights {rel[3]:.2f};  "
          f"coverage {cover:.2f}")
```

The failure of the group weights has a simple cause. With \( m=2 \) each \( s_g^2 \) is \( \sigma_g^2\chi^2(1) \), which is
often tiny, and \( 1/s_g^2 \) then gives one group an enormous weight on the strength of two nearly equal
responses. In fact \( \E(1/s_g^2)=\infty \) (@exr-het-inverse-chisq). A fitted variance function pools the
information from all groups into a few parameters and cannot be fooled by one lucky pair. The general lessons are
to estimate a smooth variance function whenever possible, and to report feasible weighted estimates with
standard errors that do not depend on the variance model being right. Part (b) of the theorem shows what those
are: the sandwich form with \( \bSigma \) estimated from the residuals, which is the subject of the next section.

## Exercises

### A. Check your understanding

::: {#exr-het-wls-residuals}
[A1]

Under @thm-het-wls(a), show that \( \E\sum_iw_i(Y_i-\x_{(i)}\T\hbeta_W)^2=\sigma^2(n-p) \), and that the raw residual
\( Y_i-\x_{(i)}\T\hbeta_W \) has variance \( \sigma^2(1-g_{ii})/w_i \), where \( g_{ii} \) is the \( i \)th leverage of the
weighted fit.
:::

::: {#exr-het-averages}
[A2]

A laboratory reports only the mean \( \bar y_i \) of \( m_i \) replicate measurements at each of \( k \) settings \( \x_{(i)} \),
where the individual measurements have common variance \( \sigma^2 \). Which weights should be used? What does
\( s_w^2 \) then estimate, and why is it not the replicate variance \( \sigma^2 \) unless the model for the means is correct?
:::

### B. Practice

::: {#exr-het-origin-efficiency}
[B1]

In the model \( y_i=\beta x_i+\varepsilon_i \) with \( x_i>0 \) and \( \Var(\varepsilon_i)=\sigma^2x_i^2 \), show that the weighted
estimate is the mean of the ratios \( y_i/x_i \), and that the efficiency of ordinary least squares relative to
it is \( (\sum_ix_i^2)^2/(n\sum_ix_i^4)\le1 \). Evaluate the efficiency for \( x=1,\dots,10 \).
:::

::: {.solution}
With \( w_i=x_i^{-2} \), \( \hat\beta_W=\sum_iw_ix_iy_i/\sum_iw_ix_i^2=n^{-1}\sum_iy_i/x_i \), with variance \( \sigma^2/n \).
Ordinary least squares gives \( \sum_ix_iy_i/\sum_ix_i^2 \), with variance
\( \sigma^2\sum_ix_i^4/(\sum_ix_i^2)^2 \). The ratio is \( (\sum_ix_i^2)^2/(n\sum_ix_i^4) \), at most one by the
Cauchy–Schwarz inequality applied to the vectors \( (x_i^2) \) and \( (1) \). For \( x=1,\dots,10 \),
\( \sum x_i^2=385 \) and \( \sum x_i^4=25{,}333 \), so the efficiency is \( 148{,}225/253{,}330=0.585 \).
:::

::: {#exr-het-duplicated-rows}
[B2]

Some programs implement integer weights \( w_i \) by repeating row \( i \) of the data \( w_i \) times and running ordinary
least squares. Show that this gives the correct \( \hbeta_W \) but divides the weighted residual sum of squares by
\( \sum_iw_i-p \) instead of \( n-p \). Which standard errors are then wrong, and in which direction?
:::

::: {#exr-het-fitted-weights}
[B3]

Weights \( \hat w_i=\lvert\hat Y_i\rvert^{-2\theta} \) depend on the fitted values. Show that they are not functions
of \( \he \) alone, and find where the symmetry argument of @thm-het-wls(c) fails. What happens for a single
regressor through the origin, \( Y_i=\beta x_i+\varepsilon_i \) with every \( x_i\ne0 \)?
:::

::: {.solution}
\( \hat{\Y}=\X\bbeta+\M\be \) depends on \( \bbeta \), and replacing \( \be \) by \( -\be \) turns it into \( \X\bbeta-\M\be \),
which does not have the same absolute entries. So \( \hat{\mathbf{W}} \) is not an even function of \( \be \), and
\( \mathbf{g}(-\be)\ne-\mathbf{g}(\be) \) in general. The feasible estimate then need not be unbiased in finite
samples, and the argument of @thm-het-wls(c) gives no guarantee.

The single regressor through the origin is an exception. There \( \hat Y_i=\hat\beta x_i \), so
\( \hat w_i=\lvert\hat\beta\rvert^{-2\theta}\lvert x_i\rvert^{-2\theta} \) whenever \( \hat\beta\ne0 \), an event of
probability one for continuously distributed errors. The random factor \( \lvert\hat\beta\rvert^{-2\theta} \) is common
to all the weights and cancels from \( \hat\beta_{\hat W}=\sum_i\hat w_ix_iY_i/\sum_i\hat w_ix_i^2 \). So the feasible
estimate equals weighted least squares with the fixed weights \( \lvert x_i\rvert^{-2\theta} \). It is unbiased by
@thm-het-wls(b), and if \( \Var(\varepsilon_i) \) is proportional to \( \lvert x_i\rvert^{2\theta} \) it is the best linear
unbiased estimator by @thm-het-wls(a); for \( \theta=1 \) this is the ratio-type case of variance proportional to
\( x_i^2 \). With two or more regressors the ratios \( \hat w_i/\hat w_j=\lvert\hat Y_j/\hat Y_i\rvert^{2\theta} \)
depend on the data, and the argument fails again.
:::

### C. Going deeper

::: {#exr-het-inverse-chisq}
[C1]

Show that \( \E(1/U)=\infty \) if \( U\sim\chi^2(1) \), and that \( \E(1/U)=1/(k-2) \) if \( U\sim\chi^2(k) \) with \( k\ge3 \).
Explain what this implies for weights \( 1/s_g^2 \) computed from \( m \) replicates.
:::

::: {.solution}
The density of \( \chi^2(k) \) is proportional to \( u^{k/2-1}e^{-u/2} \), so \( \E(1/U) \) is proportional to
\( \int_0^\infty u^{k/2-2}e^{-u/2}\,du \), which diverges at zero when \( k/2-2\le-1 \), that is, for \( k\le2 \). For
\( k\ge3 \) the integral is a gamma integral and gives \( \E(1/U)=\Gamma(k/2-1)/\{2\Gamma(k/2)\}=1/(k-2) \). With \( m \)
normal replicates, \( (m-1)s_g^2/\sigma_g^2\sim\chi^2(m-1) \), so \( 1/s_g^2 \) has infinite mean for \( m\le3 \) and
mean \( (m-1)/\{(m-3)\sigma_g^2\} \) for \( m\ge4 \). The group weights are both inflated and wildly variable in
small groups.
:::
