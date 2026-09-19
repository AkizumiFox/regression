# A single bad observation

Some failures affect a single case, such as a misplaced decimal point. How much one bad case matters depends on its leverage and on how far its response is from the rest of
the data. This section makes that exact; [Chapter 20](../ch20-residuals-leverage-influence/index.html) turns it into diagnostics.

## The effect of contaminating one response

::: {#thm-dep-outlier}
[One contaminated response]

Let \( \X \) have full column rank. For a data vector \( \y \), a case \( k \) and a shift \( \Delta \), let
\( \y_\Delta=\y+\Delta\mathbf{e}_k \), where \( \mathbf{e}_k \) is the \( k \)th standard basis vector.

::: {.enumerate options="label=(\alph*)"}
1. The estimate moves by \( \hbeta(\y_\Delta)-\hbeta(\y)=\Delta(\X\T\X)^{-1}\x_{(k)} \). The fitted vector moves by
   \( \Delta\M\mathbf{e}_k \), a vector of length \( \lvert\Delta\rvert\sqrt{h_{kk}} \). The fitted value of case \( k \) moves by
   \( \Delta h_{kk} \), and its residual by only \( \Delta(1-h_{kk}) \).

2. If \( \E(\Y)=\X\bbeta \) and \( \Cov(\Y)=\sigma^2\I \), then \( \hbeta(\Y_\Delta) \) has mean
   \( \bbeta+\Delta(\X\T\X)^{-1}\x_{(k)} \), and the residual sum of squares has mean
   \( \sigma^2(n-p)+\Delta^2(1-h_{kk}) \).

3. Suppose \( h_{kk}<1 \), and let \( \hbeta_{(k)} \) be the estimate without case \( k \). Then for every \( \y \),
   \[
   \norm{\X\bigl(\hbeta-\hbeta_{(k)}\bigr)}^2=\frac{h_{kk}}{1-h_{kk}}\cdot\frac{\hat{\varepsilon}_k^2}{1-h_{kk}} .
   \]{#eq-dep-influence}

4. Under the conditions of (c), adding \( \mathbf{e}_k \) to the model matrix as an extra column gives the
   coefficients \( \hbeta_{(k)} \) for \( \X \) and \( \hat{\delta}=y_k-\x_{(k)}\T\hbeta_{(k)}=\hat{\varepsilon}_k/(1-h_{kk}) \) for the new column.
:::

:::

::: {.proof}
(a) \( \hbeta(\y_\Delta)-\hbeta(\y)=\Delta(\X\T\X)^{-1}\X\T\mathbf{e}_k=\Delta(\X\T\X)^{-1}\x_{(k)} \). Multiplying by \( \X \)
gives \( \Delta\M\mathbf{e}_k \), with squared length \( \Delta^2\mathbf{e}_k\T\M\mathbf{e}_k=\Delta^2h_{kk} \), since \( \M \) is symmetric
and idempotent. Its \( k \)th entry is \( \Delta h_{kk} \), and the residual vector moves by
\( \Delta(\I-\M)\mathbf{e}_k \), whose \( k \)th entry is \( \Delta(1-h_{kk}) \).
(b) The mean follows from (a) and @thm-lm-moments(a). The residual vector \( (\I-\M)\Y_\Delta \) has mean
\( \Delta(\I-\M)\mathbf{e}_k \), so @eq-rv-rss-mean gives
\( \sigma^2(n-p)+\Delta^2\norm{(\I-\M)\mathbf{e}_k}^2=\sigma^2(n-p)+\Delta^2(1-h_{kk}) \).
(c) By @prp-cmp-loo(a), \( \hbeta-\hbeta_{(k)}=(\X\T\X)^{-1}\x_{(k)}\hat{\varepsilon}_k/(1-h_{kk}) \), and
\( \norm{\X(\X\T\X)^{-1}\x_{(k)}}^2=\x_{(k)}\T(\X\T\X)^{-1}\x_{(k)}=h_{kk} \).
(d) With the extra column, the criterion is
\( \sum_{i\ne k}(y_i-\x_{(i)}\T\bb)^2+(y_k-\x_{(k)}\T\bb-\delta)^2 \). For any \( \bb \) the choice
\( \delta=y_k-\x_{(k)}\T\bb \) makes the last term zero, so the minimizing \( \bb \) minimizes the first sum alone,
which is \( \hbeta_{(k)} \). The augmented matrix has full rank because \( h_{kk}<1 \) means \( \mathbf{e}_k\notin\C(\X) \) (@prp-proj-leverage). Then \( \hat{\delta}=y_k-\x_{(k)}\T\hbeta_{(k)} \), which is \( \hat{\varepsilon}_k/(1-h_{kk}) \) by @prp-cmp-loo(b).
:::

Part (a), which needs no model, says that a high-leverage case *hides* its own contamination: of a shift
\( \Delta \), the fraction \( h_{kk} \) goes into its fitted value and only \( 1-h_{kk} \) shows in its residual. Raw
residuals therefore miss the cases that matter most, and [Chapter 20](../ch20-residuals-leverage-influence/index.html) rescales them.

Part (c) makes "influence is leverage times residual" exact: a leverage factor \( h_{kk}/(1-h_{kk}) \) times the
squared *deleted* residual \( \hat{\varepsilon}_k^2/(1-h_{kk}) \). Either factor alone is harmless. Divided by \( p\,s^2 \),
@eq-dep-influence is Cook's distance (@def-res-cooks). Part (d) is the *mean-shift outlier model*: giving case \( k \) its
own parameter is the same as deleting it. Part (b) adds a second cost: \( s^2 \) is inflated by about
\( \Delta^2(1-h_{kk})/(n-p) \), which weakens every test, even for coefficients the bad case hardly moves.

::: {#exm-dep-stackloss}
[The last day of the stack loss data]

In the stack loss data of @exm-opt-stackloss-mle, \( n=21 \) days of a plant oxidizing ammonia, regress the
stack loss on air flow, water temperature and acid concentration with an intercept. The average leverage is
\( p/n=0.190 \). By @eq-dep-influence, the case whose deletion moves the fitted vector most is day
\( 21 \), with leverage \( 0.285 \) and residual \( -7.238 \). Deleting it moves the fitted
vector by a distance of \( 5.396 \). The next largest move, for day \( 1 \), is
\( 2.543 \), and the median over all days is \( 1.228 \). The largest leverage, \( 0.412 \) for day
\( 17 \), comes with a residual of only \( -1.520 \), and that day has little influence.

Without day \( 21 \), the air-flow coefficient rises from \( 0.716 \) to \( 0.889 \), the water
temperature coefficient falls from \( 1.295 \) to \( 0.817 \), and \( s \) drops from
\( 3.243 \) to \( 2.569 \). The mean-shift estimate for day \( 21 \) is
\( \hat{\delta}=-10.116 \). The listing computes every
deletion effect from the single fit (@thm-dep-outlier(c)).

Adding \( \Delta=10 \) to the response of the lowest-leverage day, day \( 5 \) with
\( h_{55}=0.052 \), moves the fitted vector by \( 2.285 \) and its own residual by
\( 9.478 \). The same shift on day \( 17 \) moves the fitted vector by \( 6.420 \) and its
residual by only \( 5.879 \).
:::

```{.python .run #cell-outliers-stackloss}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.stackloss.load_pandas().data
y = data["STACKLOSS"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP", "ACIDCONC"]].to_numpy()])
n, p = X.shape

G = np.linalg.inv(X.T @ X)
beta_hat = G @ X.T @ y
e_hat = y - X @ beta_hat
h = np.einsum("ij,jk,ik->i", X, G, X)             # leverages

# the change in beta_hat from deleting case k, without refitting
change = (G @ X.T) * (e_hat / (1 - h))            # column k is beta_hat - beta_hat_(k)
shift = h * e_hat**2 / (1 - h) ** 2               # ||X (beta_hat - beta_hat_(k))||^2
k = int(np.argmax(shift))
print("most influential case:", k + 1, " leverage", round(h[k], 3), " residual", round(e_hat[k], 3))
print("beta_hat           ", np.round(beta_hat, 3))
print("beta_hat without it", np.round(beta_hat - change[:, k], 3))
print("mean-shift estimate e_k/(1-h_kk):", round(e_hat[k] / (1 - h[k]), 3))
```

```{.python .run #cell-outliers-contaminate}
Delta = 10.0
for j in [int(np.argmin(h)), int(np.argmax(h))]:
    y_c = y.copy()
    y_c[j] += Delta                                        # contaminate one response
    b_c = G @ X.T @ y_c
    print(f"case {j + 1:2d}: h = {h[j]:.3f}, ||X(b_c - b)|| = {np.linalg.norm(X @ (b_c - beta_hat)):.3f},"
          f" |Delta| sqrt(h) = {Delta * np.sqrt(h[j]):.3f}")
```

::: {when-format="html"}
![**Figure 19.4.1.** A straight line through \( 15 \) points. (a) The same upward shift of \( 6 \) applied to a
point in the middle of the design barely moves the line; applied to the high-leverage point on the right, it
tilts the line, and that point's residual stays small. (b) A *bad leverage point*: one case keeps its response
while its \( x \) value is moved to the right, and the fitted slope (dotted: its original value) is dragged towards
zero.](outliers.svg){#fig-dep-outliers width=100%}
:::

::: {when-format="pdf"}
![A straight line through \( 15 \) points. (a) The same upward shift of \( 6 \) applied to a
point in the middle of the design barely moves the line; applied to the high-leverage point on the right, it
tilts the line, and that point's residual stays small. (b) A *bad leverage point*: one case keeps its response
while its \( x \) value is moved to the right, and the fitted slope (dotted: its original value) is dragged towards
zero.](outliers.pdf){width=100%}
:::

## Breakdown

By @thm-dep-outlier(a) the effect of one case grows without bound with \( \Delta \). Contaminated regressors can be
as bad: in panel (b) of [Figure 19.4.1](#fig-dep-outliers) a single case moved far out in \( x \) takes control of the
slope. The following definition, from Donoho and Huber (1983), makes this precise.

::: {#def-dep-breakdown}
[Breakdown point]

Let \( T \) be an estimator that maps a data set \( Z \) of \( n \) cases \( (\x_{(i)},y_i) \) to a vector in \( \Real^p \). For
\( m\le n \), let \( b(m;T,Z) \) be the supremum of \( \norm{T(Z')-T(Z)} \) over all data sets \( Z' \) obtained from \( Z \) by
replacing any \( m \) of its cases with arbitrary values, for which \( T(Z') \) is defined. The **finite-sample
breakdown point** of \( T \) at \( Z \) is
\[
\varepsilon^*(T,Z)=\min\Bigl\{\frac mn:\ b(m;T,Z)=\infty\Bigr\}.
\]
:::

::: {#prp-dep-breakdown}
[Least squares breaks down at one case]

For least squares and every data set with a full-rank model matrix, \( \varepsilon^*=1/n \).
:::

::: {.proof}
Since \( \X \) has full rank, some row \( \x_{(k)} \) is nonzero, and then \( (\X\T\X)^{-1}\x_{(k)}\ne\bzero \). Replace case \( k \)
by \( (\x_{(k)},y_k+\Delta) \). The model matrix is unchanged, so the estimate is defined, and by
@thm-dep-outlier(a) it moves by \( \lvert\Delta\rvert\norm{(\X\T\X)^{-1}\x_{(k)}} \), which is unbounded in \( \Delta \). So
\( b(1)=\infty \). No smaller fraction than \( 1/n \) is possible.
:::

So least squares offers no protection against gross errors. The sample median, by contrast, has breakdown point close to one half (@exr-dep-median-breakdown). High-breakdown
regression estimators exist too, at a cost in efficiency discussed with robust fitting in [Section 20.6](../ch20-residuals-leverage-influence/06-robust.html). Hampel (1974) introduced the related *influence function*, the
effect of an infinitesimal contamination at a point \( (\x_0,y_0) \). For least squares it is unbounded in both
\( \x_0 \) and \( y_0 \) (@exr-dep-influence-function).

::: {.idea}
High leverage makes a case both more influential and better hidden, and one case can move the estimate
arbitrarily far.
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-median-breakdown}
[A1]

Show that the sample median of \( n \) numbers has breakdown point \( \lfloor(n+1)/2\rfloor/n \), in the sense of
@def-dep-breakdown with \( p=1 \) and no regressors.
:::

::: {#exr-dep-other-residuals}
[A2]

In @thm-dep-outlier(a), show that the residual of every other case \( i \) moves by \( -\Delta m_{ik} \), where \( m_{ik} \) is
an entry of \( \M \). Which cases' residuals are most affected by a contaminated high-leverage case?
:::

### B. Practice

::: {#exr-dep-bad-leverage}
[B1]

In a straight-line regression, move the \( x \) value of case \( k \) to \( t \) while keeping \( y_k \) fixed. Show that the
fitted slope tends to \( 0 \) as \( t\to\infty \). If instead \( y_k=c\,t+d \) moves with \( t \), show that the slope tends to \( c \).
:::

::: {.solution}
Let \( \bar{x}_t \) and \( \bar{y} \) be the means with case \( k \) at \( t \). The slope is
\( \sum_i(x_i-\bar{x}_t)(y_i-\bar{y})/\sum_i(x_i-\bar{x}_t)^2 \). As \( t\to\infty \), \( t-\bar{x}_t=(1-1/n)t+O(1) \) while the other
deviations \( x_i-\bar{x}_t \) are \( -t/n+O(1) \). The denominator is \( (1-1/n)t^2+O(t) \). With \( y_k \) fixed, the numerator is
\( O(t) \), since each term is a deviation of order \( t \) times a bounded \( y_i-\bar{y} \), so the slope is \( O(1/t) \). With
\( y_k=ct+d \), put \( r_i=y_i-cx_i \), so that \( r_k=d \) and all the \( r_i \) stay bounded. Since \( y_i=cx_i+r_i \) and least
squares is linear in the response, the slope is \( c \) plus the slope of \( r \) on \( x \), which is \( O(1/t) \) as before. The far case dictates the slope.
:::

::: {#exr-dep-contamination-power}
[B2]

Suppose case \( k \) is contaminated by \( \Delta \), and consider the \( t \) statistic for a coefficient \( \beta_j \) whose estimate is
unaffected, \( \bigl[(\X\T\X)^{-1}\x_{(k)}\bigr]_j=0 \). Using @thm-dep-outlier(b), explain how the contamination still
changes the test, and show that its effect on \( \E(s^2) \) is largest for low-leverage cases.
:::

### C. Going deeper

::: {#exr-dep-influence-function}
[C1]

Let \( (\x,Y) \) have a joint distribution \( F \) with \( \bSigma=\E(\x\x\T) \) positive definite, and let
\( \bbeta^*(F)=\bSigma^{-1}\E(\x Y) \) be the population least squares coefficient. For the mixture
\( F_\epsilon=(1-\epsilon)F+\epsilon\,\delta_{(\x_0,y_0)} \), which puts mass \( \epsilon \) at a single point, show that
\[
\frac{d}{d\epsilon}\bbeta^*(F_\epsilon)\Big|_{\epsilon=0}=\bSigma^{-1}\x_0\bigl(y_0-\x_0\T\bbeta^*(F)\bigr),
\]
and compare with @thm-dep-outlier(c) for a sample of size \( n \).
:::

::: {.solution}
Under \( F_\epsilon \), \( \bSigma_\epsilon=(1-\epsilon)\bSigma+\epsilon\x_0\x_0\T \) and
\( \mathbf{m}_\epsilon=(1-\epsilon)\E(\x Y)+\epsilon\x_0y_0 \), and \( \bSigma_\epsilon\bbeta^*(F_\epsilon)=\mathbf{m}_\epsilon \). Differentiate at
\( \epsilon=0 \), writing \( \mathbf{d} \) for the derivative: \( (\x_0\x_0\T-\bSigma)\bbeta^*+\bSigma\mathbf{d}=\x_0y_0-\E(\x Y) \). Since \( \bSigma\bbeta^*=\E(\x Y) \), this gives
\( \bSigma\mathbf{d}=\x_0(y_0-\x_0\T\bbeta^*) \). The influence is the product of a leverage-like factor \( \bSigma^{-1}\x_0 \) and a
residual, unbounded in both. In a sample, replacing \( F \) by the empirical distribution and \( \epsilon \) by \( 1/n \) gives
\( n^{-1}(\X\T\X/n)^{-1}\x_{(k)}\hat{\varepsilon}_k \), which is @prp-cmp-loo(a) without the factor \( 1/(1-h_{kk}) \).
:::
