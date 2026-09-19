# Collinearity and prediction

By @thm-col-variance(d), collinearity leaves the total variance of the fitted values, the predictions at
the design points, untouched. At a new point \( \x_0 \), collinearity damages a prediction exactly to the
extent that the point asks for information along a weak direction. @exm-dep-along-across showed this
contrast on simulated data; this section proves it in general and measures it on real data.

## The variance of a prediction

Keep the notation of [Section 26.1](01-variance.html): \( \X\T\X=\sum_\ell\lambda_\ell\bv_\ell\bv_\ell\T \), and let
\[
m_\ell=\frac{\lambda_\ell}{n}=\frac1n\sum_{i=1}^n\bigl(\x_{(i)}\T\bv_\ell\bigr)^2
\]
be the mean square of the design rows along the direction \( \bv_\ell \). A small \( m_\ell \) means that the rows
hardly vary along \( \bv_\ell \): the data cloud is thin in that direction. The estimate of the mean response at
\( \x_0 \) is \( \hat{Y}_0=\x_0\T\hbeta \).

::: {#prp-col-prediction}
[Prediction under collinearity]

Let \( \X \) have rank \( p \) and \( \Cov(\Y)=\sigma^2\I \).

::: {.enumerate options="label=(\alph*)"}
1. \( \displaystyle\Var(\hat{Y}_0)=\frac{\sigma^2}{n}\sum_{\ell=1}^p\frac{(\x_0\T\bv_\ell)^2}{m_\ell} \).

2. If \( \x_0=\X\T\bw \) for some \( \bw\in\Real^n \), then \( \Var(\hat{Y}_0)=\sigma^2\norm{\M\bw}^2\le\sigma^2\norm{\bw}^2 \). If \( \x_0 \) is a
   convex combination of the rows of \( \X \), then \( \Var(\hat{Y}_0)\le\sigma^2\max_ih_{ii}\le\sigma^2 \). Neither bound involves
   the eigenvalues.

3. Let \( \x_0 \) be random, independent of \( \Y \), with \( \E(\x_0\x_0\T)=\bSigma_0 \). Then
   \[
   \E\bigl\{\Var(\hat{Y}_0\mid\x_0)\bigr\}=\frac{\sigma^2}{n}\sum_{\ell=1}^p\frac{\bv_\ell\T\bSigma_0\bv_\ell}{m_\ell} .
   \]
   If \( \x_0 \) is drawn uniformly from the rows of \( \X \), this equals \( \sigma^2p/n \), whatever the eigenvalues. For every
   \( \bSigma_0 \) it is at least \( (\sigma^2/n)\,\bv_p\T\bSigma_0\bv_p/m_p \).
:::

:::

::: {.proof}
(a) This is @prp-dep-collinear-directions(a) with \( \mathbf{a}=\x_0 \) and \( \lambda_\ell=nm_\ell \). The second expression for
\( m_\ell \) is \( \norm{\X\bv_\ell}^2/n \).

(b) By @thm-proj-M-formula,
\[
\Var(\hat{Y}_0)=\sigma^2\bw\T\X(\X\T\X)^{-1}\X\T\bw=\sigma^2\bw\T\M\bw=\sigma^2\norm{\M\bw}^2,
\]
and \( \norm{\M\bw}\le\norm{\bw} \) (@prp-proj-trace-rank). If \( \x_0=\sum_iw_i\x_{(i)} \) with
\( w_i\ge0 \) and \( \sum_iw_i=1 \), then \( \x_0=\X\T\bw \), and by the triangle inequality and @prp-proj-leverage(a),
\[
\norm{\M\bw}\le\sum_iw_i\norm{\M\mathbf{e}_i}=\sum_iw_i\sqrt{h_{ii}}\le\max_i\sqrt{h_{ii}} ,
\]
and \( h_{ii}\le1 \).

(c) Conditionally on \( \x_0 \), (a) holds. Taking expectations,
\( \E(\x_0\T(\X\T\X)^{-1}\x_0)=\tr\bigl((\X\T\X)^{-1}\bSigma_0\bigr)=\sum_\ell\bv_\ell\T\bSigma_0\bv_\ell/\lambda_\ell \). If \( \x_0 \) is uniform on the rows,
\( \bSigma_0=\X\T\X/n \), and the trace is \( \tr(\I_p)/n=p/n \). The lower bound keeps the term \( \ell=p \), since all
terms are nonnegative.
:::

Part (b) answers the first question of @exr-dep-convex-hull: inside the convex hull of the data, a
prediction is never worse than the worst fitted value, however collinear the regressors. Part (a) says
what "the pattern of the data" means: a point outside the hull is still predicted well if its coordinates
along the weak directions are small compared with \( \sqrt{m_\ell} \), and a point inside the range of every
regressor separately can be far off the pattern (@exm-ci-hidden-extrapolation). By part (c), if future
regressors are generated the way the design was, the average prediction variance is \( \sigma^2p/n \); if a
relation between the regressors breaks down, the term \( \bv_p\T\bSigma_0\bv_p/m_p \) can be enormous.

::: {.warning}
[Variance is not the whole risk]

@prp-col-prediction is about variance under a correct model. Off the pattern of the data the model
itself has not been checked, so its bias there is unknown
([Section 12.3](../ch12-intervals-and-bands/03-prediction.html)). A small prediction variance is not a
guarantee.
:::

::: {#exm-col-macro-prediction}
[Predicting consumption along and across the pattern]

In the consumption regression of @exm-col-macro-vif, \( s=58.1 \) (billions of 2005 dollars). At 1985Q1, a
design point with population \( 237.9 \) million, the estimated mean has standard deviation \( 8.5 \), below the
design average \( s\sqrt{p/n}=10.8 \). Midway between the first quarters of 1969 and 1999 the standard
deviation is \( 7.8 \), and by @prp-col-prediction(b) no point in the convex hull of the data can exceed
\( s\sqrt{\max_ih_{ii}}=23.3 \).

Now keep every regressor at its 1985 value except population, which is lowered by ten per cent, to
\( 214.1 \) million. Population did take that value, in 1974, and income took its 1985 value, but never
together. The standard deviation of the estimated mean becomes \( 36.8 \), \( 4.3 \) times as large. Decomposing
the variance over the canonical directions of the scaled design, as in the proof of
@thm-col-decomposition, the direction with condition index \( 270.9 \), the income–population
near-dependency of @exm-col-macro-bkw, contributes \( 93.9 \) per cent of the variance at the new point and
\( 0.2 \) per cent at the 1985 quarter itself.

[Figure 26.4.1](#fig-col-prediction) shows \( 2000 \) parametric replicates: responses generated from the fitted
model with normal errors of standard deviation \( s \), each refitted. The income and population coefficients
are strongly negatively correlated: their exact correlation is \( -0.842 \) (\( -0.832 \) in the replicates). The prediction at 1985Q1 barely moves while the one
off the pattern spreads widely.
:::

::: {when-format="html"}
![**Figure 26.4.1.** Parametric replicates for @exm-col-macro-prediction. (a) Income and population
coefficients. (b) Estimated mean, centred, at 1985Q1 and at the same point with population ten per cent lower.](prediction_along_across.svg){#fig-col-prediction width=100%}
:::

::: {when-format="pdf"}
![Parametric replicates for @exm-col-macro-prediction. (a) Income and population
coefficients. (b) Estimated mean, centred, at 1985Q1 and at the same point with population ten per cent lower.](prediction_along_across.pdf){width=100%}
:::

```{.python .run #cell-prediction-setup}
import numpy as np
import statsmodels.api as sm
macro = sm.datasets.macrodata.load_pandas().data
names = ["realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
X = sm.add_constant(macro[names].to_numpy())
y = macro["realcons"].to_numpy()
n, p = X.shape
fit = sm.OLS(y, X).fit()
s = np.sqrt(fit.scale)                                # estimate of sigma
G = np.linalg.inv(X.T @ X)

x_on = X[104]                                         # 1985Q1, a design point
x_mix = 0.5 * (X[40] + X[160])                        # midway between 1969Q1 and 1999Q1
x_off = x_on.copy()
x_off[2] *= 0.90                                      # 1985Q1 with population 10% lower
for label, x0 in [("1985Q1", x_on), ("midway", x_mix), ("off pattern", x_off)]:
    print(f"{label:13s} fitted {x0 @ fit.params:8.1f}   sd {s * np.sqrt(x0 @ G @ x0):6.1f}")
```

```{.python .run #cell-prediction-decompose}
d = np.linalg.norm(X, axis=0)
U, mu, Vt = np.linalg.svd(X / d, full_matrices=False) # scaled design, as in Section 26.3

def shares(x0):
    """Contributions of the canonical directions to Var(x0' beta_hat), as fractions."""
    c = (Vt @ (x0 / d)) ** 2 / mu**2                  # (x0_scaled' v_k)^2 / mu_k^2
    return c / c.sum()

np.set_printoptions(precision=3, suppress=True)
print("condition indices:", mu[0] / mu)
print("shares, 1985Q1     :", shares(x_on))
print("shares, off pattern :", shares(x_off))
```

```{.python .run #cell-prediction-simulate}
rng = np.random.default_rng(2604)
reps = 2000
Ystar = X @ fit.params + s * rng.normal(size=(reps, n))      # parametric replicates, same design
Bstar = np.linalg.solve(X.T @ X, X.T @ Ystar.T).T
print("sd of income, population coefficients:", Bstar[:, 1:3].std(axis=0))
print("correlation of the two:", np.corrcoef(Bstar[:, 1], Bstar[:, 2])[0, 1])
print("sd of predictions, 1985Q1 / off:", (Bstar @ x_on).std(), (Bstar @ x_off).std())
```

## Coefficients against predictions

A coefficient \( \beta_j \) is the difference between the mean responses at two points that differ by one unit in
\( x_j \) and agree in every other regressor; with collinear data at least one of them lies off the pattern. A
regression built for forecasting, where the relations among the regressors are expected to persist, can
ignore collinearity; one built to attribute effects to individual regressors cannot.

## Exercises

### A. Check your understanding

::: {#exr-col-three-points}
[A1]

For the design of @exm-col-two-regressors, rank the following new points by the variance of the fitted
mean, without computing: the centroid; \( (x_1,x_2)=(3,3) \); \( (x_1,x_2)=(0.5,-0.5) \); the average of the two design
points with the largest \( x_1 \). Which bound of @prp-col-prediction applies to each? Compare
[Figure 19.6.1](../ch19-theory-of-departures/06-collinearity.html#fig-dep-collinearity).
:::

### B. Practice

::: {#exr-col-min-norm}
[B1]

Show that the bound \( \Var(\hat{Y}_0)\le\sigma^2\norm{\bw}^2 \) of @prp-col-prediction(b) is attained iff \( \bw\in\C(\X) \), and that
\( \min\{\norm{\bw}^2:\X\T\bw=\x_0\} \) equals \( \Var(\hat{Y}_0)/\sigma^2 \). Compare with @prp-ci-new-leverage(a).
:::

::: {.solution}
\( \norm{\M\bw}^2=\norm{\bw}^2-\norm{(\I-\M)\bw}^2 \), so equality holds iff \( (\I-\M)\bw=\bzero \), that is, \( \bw\in\C(\X) \). Every
\( \bw \) with \( \X\T\bw=\x_0 \) has the same \( \M\bw \), because two such vectors differ by an element of
\( \Null(\X\T)=\C(\X)\perpc \). The minimum norm is therefore attained at \( \bw=\M\bw \) and equals
\( \norm{\M\bw}^2=\Var(\hat{Y}_0)/\sigma^2 \). This is @prp-ci-new-leverage(a) with \( h_0=\Var(\hat{Y}_0)/\sigma^2 \): the leverage of a new
point is the smallest squared norm of a set of weights that builds \( \x_0 \) from the rows of the design.
:::

::: {#exr-col-along-trend}
[B2]

Let the model have an intercept and centred regressors, so that the rows of \( \X \) are \( (1,\bz_{(i)}\T) \) with
\( \sum_i\bz_{(i)}=\bzero \). For a design row \( i \) and \( t>1 \), let \( \x_0=(1,t\bz_{(i)}\T)\T \), an extrapolation from the
centroid along the ray through a design point. Show that \( \Var(\hat{Y}_0)=\sigma^2\bigl(1/n+t^2(h_{ii}-1/n)\bigr) \). Why does
collinearity not enter, although \( \x_0 \) is outside the convex hull of the data when \( t \) is large?
:::

::: {.solution}
With centred regressors \( \X\T\X=\diag(n,\Z\T\Z) \), where \( \Z \) has rows \( \bz_{(i)}\T \), so
\( \Var(\hat{Y}_0)=\sigma^2\bigl(1/n+t^2\bz_{(i)}\T(\Z\T\Z)^{-1}\bz_{(i)}\bigr) \), and
\( h_{ii}=1/n+\bz_{(i)}\T(\Z\T\Z)^{-1}\bz_{(i)} \) by @prp-proj-leverage-mahalanobis. The slope part of \( \x_0 \) is \( t \) times that of a
row, so its coordinate along every eigenvector is \( t \) times the row's, and in each term of @prp-col-prediction(a)
it is penalized relative to the data exactly as much as the row is, times \( t^2 \). Collinearity could enter only
through a component along a weak direction that is large compared with the spread of the data there, and the
row itself has no such component.
:::

### C. Going deeper

::: {#exr-col-structural-change}
[C1]

In the design of @exm-col-two-regressors, suppose future regressor values satisfy
\( x_2-x_1=\tau\zeta \) with \( \zeta \) independent of \( x_1 \), \( \E\zeta=0 \), \( \Var\zeta=1 \), and \( x_1 \) distributed as in the sample.
Using @prp-col-prediction(c), show that the average prediction variance is approximately
\( (\sigma^2/n)\bigl(2+\tau^2/(2m_2)\bigr) \), where \( m_2=\lambda_2/n \) is the mean square of the design along the weak
slope direction. Check that this is \( \sigma^2p/n \) when \( \tau^2 \) equals the within-sample mean square of
\( x_2-x_1 \), and discuss what happens when \( \tau \) is several times larger. What would a sensible forecaster do?
:::
