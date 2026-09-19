# Deletion diagnostics

A case cannot fairly judge itself, since its own response pulled the fit towards it. So set the case aside, fit the other \( n-1 \) cases, and ask how well that fit predicts it and how far the fit moved.
[Chapter 10](../ch10-computation/index.html) showed that one regression gives all \( n \) answers (@prp-cmp-loo). This
section collects the formulas and proves that the externally studentized residual has an exact \( t \) law, which turns
it into a test for an outlier.

## Leave-one-out formulas

Write \( \hbeta_{(i)} \) for the least squares estimate from the cases other than \( i \), \( \hY_{(i)}=\X\hbeta_{(i)} \) for the
fitted values it gives at *all* \( n \) cases (including case \( i \)), \( \text{SSE}_{(i)} \) for its residual sum of squares
and \( s_{(i)}^2=\text{SSE}_{(i)}/(n-p-1) \).

::: {#thm-res-deletion}
[Leave-one-out formulas]

Let \( \X \) have full column rank, and let \( \X_{(i)} \) be \( \X \) without row \( i \).

::: {.enumerate options="label=(\alph*)"}
1. \( \X_{(i)} \) has full column rank iff \( h_{ii}<1 \). Assume this from now on.

2. \( \hbeta-\hbeta_{(i)}=(\X\T\X)^{-1}\x_{(i)}\,\hat{\varepsilon}_i/(1-h_{ii}) \).

3. \( \hY-\hY_{(i)}=\X(\X\T\X)^{-1}\x_{(i)}\,\hat{\varepsilon}_i/(1-h_{ii}) \); that is, the fitted value of case \( j \) moves by
   \( h_{ji}\hat{\varepsilon}_i/(1-h_{ii}) \). In particular case \( i \)'s own fitted value moves by \( h_{ii}\hat{\varepsilon}_i/(1-h_{ii}) \).

4. The **prediction residual** of case \( i \) is \( y_i-\x_{(i)}\T\hbeta_{(i)}=\hat{\varepsilon}_i/(1-h_{ii}) \).

5. \( (n-p-1)\,s_{(i)}^2=(n-p)\,s^2-\hat{\varepsilon}_i^2/(1-h_{ii}) \).
:::

:::

::: {.proof}
(a) \( \X_{(i)}\T\X_{(i)}=\X\T\X-\x_{(i)}\x_{(i)}\T \), whose determinant is \( (1-h_{ii})\det(\X\T\X) \) by @exm-mat-deletion.
So it is nonsingular iff \( h_{ii}\ne1 \), and \( h_{ii}\le1 \) always. Parts (b), (d) and (e) are @prp-cmp-loo(a), (b)
and (c), the last divided through by the degrees of freedom. Part (c) is (b) multiplied by \( \X \); its \( j \)th entry
is \( \x_{(j)}\T(\X\T\X)^{-1}\x_{(i)}\hat{\varepsilon}_i/(1-h_{ii}) \), and \( \x_{(j)}\T(\X\T\X)^{-1}\x_{(i)}=h_{ji} \).
:::

Every quantity is a residual times a function of the design: a case with zero residual changes nothing when deleted,
and high leverage magnifies a residual by \( 1/(1-h_{ii}) \). When \( h_{ii}=1 \) the case determines one direction of
\( \C(\X) \) by itself (@prp-proj-leverage(d)), and without it a parameter is not estimable; no diagnostic can check it.

## Prediction residuals and PRESS

The prediction residual \( \hat{\varepsilon}_i/(1-h_{ii}) \) is the error in predicting \( y_i \) by a fit that never saw it.
Its variance is \( \sigma^2/(1-h_{ii}) \), and dividing it by an estimate of \( \sigma/\sqrt{1-h_{ii}} \) gives back \( r_i \) (with
\( s \)) or \( t_i \) (with \( s_{(i)} \)). The sum of squares of the prediction residuals is the PRESS statistic of Allen (1974),
\[
\text{PRESS}=\sum_{i=1}^n\Bigl(\frac{\hat{\varepsilon}_i}{1-h_{ii}}\Bigr)^2 ,
\]{#eq-res-press}

a leave-one-out estimate of prediction error that Chapter 29 uses for choosing between models.

::: {#prp-res-press}
[Expected PRESS]

Under the second-moment model of [Section 20.1](01-kinds-of-residuals.html),
\[
\E(\text{PRESS})=\sigma^2\sum_{i=1}^n\frac1{1-h_{ii}}\ \ge\ \sigma^2\,\frac{n^2}{n-p},
\]
with equality iff all leverages equal \( p/n \). By comparison, the expected sum of squared errors in predicting
new responses at the same \( n \) design points is \( \sigma^2(n+p) \) (@prp-cor-optimism(a)), which is smaller than
\( \sigma^2n^2/(n-p) \) by \( \sigma^2p^2/(n-p) \).
:::

::: {.proof}
\( \E(\hat{\varepsilon}_i^2)=\Var(\hat{\varepsilon}_i)=\sigma^2(1-h_{ii}) \), and dividing by \( (1-h_{ii})^2 \) and summing gives the
first identity. The function \( u\mapsto1/(1-u) \) is strictly convex on \( [0,1) \), so by Jensen's inequality
\( n^{-1}\sum_i1/(1-h_{ii})\ge1/(1-\bar h) \) with \( \bar h=p/n \), with equality iff all \( h_{ii} \) are equal. This gives the
bound. Finally \( n^2-(n+p)(n-p)=p^2 \).
:::

So PRESS overcorrects the optimism of the in-sample error (@prp-cor-optimism) slightly, because each case is predicted
from a design that lacks it. It also weights cases by \( (1-h_{ii})^{-2} \): at \( h_{ii}=0.5 \), four times the squared raw
residual. For the state regression with all \( 51 \) jurisdictions, \( \text{SSE}=142.53 \) and \( \text{PRESS}=230.52 \) (@exm-cmp-dc), and the District of Columbia alone contributes \( 9.066^2 \), a share of \( 0.36 \). For this design
\( \E(\text{PRESS})/(n\sigma^2)=1.096 \), against the bound \( n/(n-p)=1.085 \) and \( (n+p)/n=1.078 \): the District's residual, not
the unequal leverages, inflates PRESS.

## The externally studentized residual

As a test statistic, \( r_i \) is bounded by \( \sqrt{n-p} \) and has a nonstandard law (@prp-res-internal). The externally
studentized residual repairs both, and its law is clearest from the alternative it is meant to detect.

The **mean-shift outlier model** for case \( i \) is
\[
\Y=\X\bbeta+\Delta\,\vect{e}_i+\be,\qquad \be\sim\Normal_n(\bzero,\sigma^2\I),
\]{#eq-res-mean-shift}

in which case \( i \) alone has its mean displaced by an unknown \( \Delta \), as by a misrecorded digit. Testing \( \Delta=0 \)
tests whether case \( i \) is an outlier.

::: {#thm-res-external-t}
[Law of the externally studentized residual]

Assume @eq-res-mean-shift with \( h_{ii}<1 \) and \( n-p\ge2 \).

::: {.enumerate options="label=(\alph*)"}
1. In the model with columns \( [\X,\vect{e}_i] \), the least squares estimate of the shift is
   \( \hat{\Delta}=\hat{\varepsilon}_i/(1-h_{ii})=y_i-\x_{(i)}\T\hbeta_{(i)} \), the estimate of \( \bbeta \) is \( \hbeta_{(i)} \), the residual mean square is
   \( s_{(i)}^2 \), and the \( t \) statistic for \( \Delta=0 \) is \( t_i \).

2. \( t_i\sim t\bigl(n-p-1,\ \Delta\sqrt{1-h_{ii}}/\sigma\bigr) \). In particular, if \( \Delta=0 \), then \( t_i\sim t(n-p-1) \).

3. \( t_i=r_i\sqrt{(n-p-1)/(n-p-r_i^2)} \), an increasing function of \( r_i \) that maps
   \( (-\sqrt{n-p},\sqrt{n-p}) \) onto the whole real line.

4. *(Bonferroni outlier test.)* If no case is shifted, that is under @eq-opt-normal-model, and \( h_{jj}<1 \) for every \( j \), then
   \[
   \Pr\Bigl\{\max_j\lvert t_j\rvert>t_{n-p-1,\,\alpha/(2n)}\Bigr\}\le\alpha .
   \]
   Equivalently, the adjusted \( p \)-value \( \min\bigl(1,\,n\,p_{(1)}\bigr) \), where \( p_{(1)} \) is the smallest of the
   single-case two-sided \( p \)-values from (b), is valid.
:::

:::

::: {.proof}
(a) Since \( h_{ii}<1 \), \( \vect{e}_i\notin\C(\X) \) (@prp-proj-leverage(d)), so \( [\X,\vect{e}_i] \) has full column rank \( p+1 \). By
@exr-est-single-dummy, \( \hat{\Delta}=\hat{\varepsilon}_i/(1-h_{ii}) \), the enlarged fit reproduces \( y_i \) exactly, and its estimate of
\( \bbeta \) is \( \hbeta_{(i)} \). Its residuals are therefore zero at case \( i \) and equal to the residuals of the deletion
fit elsewhere, so its residual sum of squares is \( \text{SSE}_{(i)} \) on \( n-p-1 \) degrees of freedom, and its residual
mean square is \( s_{(i)}^2 \). By @eq-proj-partitioned-inverse, the diagonal entry of
\( \bigl([\X,\vect{e}_i]\T[\X,\vect{e}_i]\bigr)^{-1} \) belonging to \( \vect{e}_i \) is \( 1/\norm{(\I-\M)\vect{e}_i}^2=1/(1-h_{ii}) \). So the standard
error of \( \hat{\Delta} \) is \( s_{(i)}/\sqrt{1-h_{ii}} \), and
\[
\frac{\hat{\Delta}}{s_{(i)}/\sqrt{1-h_{ii}}}=\frac{\hat{\varepsilon}_i}{(1-h_{ii})}\cdot\frac{\sqrt{1-h_{ii}}}{s_{(i)}}=t_i .
\]
(b) The enlarged model is a normal linear model of full rank \( p+1<n \), and \( \Delta \) is one of its coefficients. By
@thm-glh-t-test(a) the \( t \) statistic for \( \Delta=0 \) has the \( t(n-p-1,\delta) \) law with
\( \delta=\Delta/\bigl(\sigma\sqrt{1/(1-h_{ii})}\bigr)=\Delta\sqrt{1-h_{ii}}/\sigma \).
(c) By @thm-res-deletion(e) and \( \hat{\varepsilon}_i^2/(1-h_{ii})=r_i^2s^2 \),
\( (n-p-1)s_{(i)}^2=s^2(n-p-r_i^2) \). Hence \( t_i^2=r_i^2s^2/s_{(i)}^2=r_i^2(n-p-1)/(n-p-r_i^2) \), and \( t_i \) has the sign of
\( r_i \). The map \( r\mapsto r\sqrt{(n-p-1)/(n-p-r^2)} \) is increasing on \( (-\sqrt{n-p},\sqrt{n-p}) \) and tends to
\( \pm\infty \) at the ends.
(d) By (b) with \( \Delta=0 \), each event \( A_j=\{\lvert t_j\rvert>t_{n-p-1,\alpha/(2n)}\} \) has probability exactly \( \alpha/n \).
Apply @thm-mc-bonferroni(a) to their union. The \( p \)-value form is @thm-mc-bonferroni(c) with every \( \alpha_j=\alpha/n \).
:::

By (c), \( r_i \) and \( t_i \) carry the same information, but \( t_i \) is unbounded and has a standard law, so it is the one to
report. @exr-res-direct-t gives a second proof of (b). Part (b) also shows how leverage hides outliers: the noncentrality is
\( \Delta\sqrt{1-h_{ii}}/\sigma \), because a shift at a high-leverage case is partly absorbed by the fit. At the state design,
a shift of \( 4\sigma \) at the District of Columbia (\( h_{ii}=0.506 \)) gives noncentrality \( 2.81 \), and the Bonferroni test
of (d) detects it with probability \( 0.258 \); at Missouri, the case of lowest leverage (\( h_{ii}=0.021 \)), the
noncentrality is \( 3.96 \) and the probability \( 0.666 \).

## Testing for an outlier

A case suspected in advance can be tested with @thm-res-external-t(b). Usually the suspect is chosen *because* its
residual is largest, and the question is whether \( \max_i\lvert t_i\rvert \) is surprising under a model with no outliers. The
\( t_i \) are dependent, but the Bonferroni test of @thm-res-external-t(d) needs no joint law, and it is nearly exact
because two cases rarely exceed the high cut-off together. At the state design (\( n=51 \), \( p=4 \)), the simulated
familywise level of the \( 5\% \) test is \( 0.0497 \) (standard error \( 0.0005 \), from \( 200000 \) normal samples).

```{.python .run #cell-bonferroni-level-level}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data
X = np.column_stack([np.ones(len(data)), data[["poverty", "single", "urban"]]])

def max_abs_t(X, Y):
    """max_i |t_i| for each column of Y (one simulated response vector per column)."""
    n, p = X.shape
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)[:, None]
    E = Y - Q @ (Q.T @ Y)
    sse = np.sum(E ** 2, axis=0)
    s2_del = (sse - E ** 2 / (1 - h)) / (n - p - 1)
    return np.max(np.abs(E) / np.sqrt(s2_del * (1 - h)), axis=0)

n, p = X.shape
alpha = 0.05
crit = stats.t.isf(alpha / (2 * n), n - p - 1)
rng = np.random.default_rng(51)
reps = 4000
T = max_abs_t(X, rng.normal(size=(n, reps)))          # no outliers: errors N(0, 1)
print(f"Bonferroni bound {alpha}, simulated familywise level {np.mean(T > crit):.4f}")
```

::: {#exm-res-dc-outlier}
[Is the District of Columbia an outlier?]

For the murder-rate regression with all \( 51 \) jurisdictions, the District of Columbia has raw residual \( 4.475 \),
leverage \( 0.506 \), \( r_i=3.658 \) and \( t_i=4.278 \) (as in @exm-cmp-dc); deleting it lowers \( s=1.741 \) to
\( s_{(i)}=1.489 \). Its single-case \( p \)-value is \( 9.4\times 10^{-5} \), and the Bonferroni-adjusted value \( 0.0048 \)
(critical value \( 3.522 \)). The next case, Mississippi, has \( t_i=-2.877 \) and adjusted \( p \)-value \( 0.310 \). Among the
\( 50 \) states alone the most extreme is Louisiana, \( t_i=3.299 \), adjusted \( p \)-value \( 0.095 \); such sequential testing is
not covered by @thm-res-external-t(d), and [Section 20.5](05-masking.html) shows it can fail badly.

The test says that the District's murder rate is incompatible with the linear relation that holds for the states. It
does not say the number is wrong or that the case should be dropped. The District is a city, with regressors outside
the range of the states; a defensible analysis restricts the population to states and says so.
:::

```{.python .run #cell-state-diagnostics-diagnostics}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.statecrime.load_pandas().data           # District of Columbia kept
y = data["murder"].to_numpy()
X = np.column_stack([np.ones(len(y)), data[["poverty", "single", "urban"]]])
names = np.array(data.index)

n, p = X.shape
Q, R = np.linalg.qr(X)
beta = np.linalg.solve(R, Q.T @ y)
e = y - X @ beta                                   # raw residuals
h = np.sum(Q ** 2, axis=1)                         # leverages
s2 = e @ e / (n - p)
r = e / np.sqrt(s2 * (1 - h))                      # internally studentized
s2_del = (e @ e - e ** 2 / (1 - h)) / (n - p - 1)  # s^2 without case i
t = e / np.sqrt(s2_del * (1 - h))                  # externally studentized

p_single = 2 * stats.t.sf(np.abs(t), n - p - 1)    # one case, chosen in advance
p_bonf = np.minimum(1, n * p_single)               # Bonferroni over all n cases
crit = stats.t.isf(0.05 / (2 * n), n - p - 1)      # 5% Bonferroni critical value
for i in np.argsort(-np.abs(t))[:3]:
    print(f"{names[i]:22s} h = {h[i]:.3f}  r = {r[i]:6.3f}  t = {t[i]:6.3f}"
          f"  Bonferroni p = {p_bonf[i]:.4f}")
print(f"critical value |t| > {crit:.3f}")
```

::: {.warning}
An outlier test tests *the model* at one case, including normality, constant variance and the mean function. A large
\( \lvert t_i\rvert \) can come from a gross error, from heavy-tailed errors, from variance that grows where case \( i \) lies, or from a
wrong mean function there, and the test cannot tell these apart
([Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html), [Chapter 22](../ch22-transformations/index.html)).
:::

## Exercises

### A. Check your understanding

::: {#exr-res-dc-conversion}
[A1]

Use @thm-res-external-t(c) to compute the District of Columbia's \( t_i \) from its \( r_i=3.658 \), with \( n=51 \) and \( p=4 \).
:::

::: {.solution}
\( t_i=3.658\sqrt{46/(47-3.658^2)}=3.658\sqrt{46/33.62}=3.658\times1.1697\approx4.279 \), which agrees with
the recorded \( 4.278 \) to the rounding of \( r_i \).
:::

::: {#exr-res-half-leverage}
[A2]

A case has leverage \( 0.5 \) and residual \( 1 \). By how much does its own fitted value change when it is deleted, and
what is its prediction residual? What is the prediction residual of a case with leverage \( 0.9 \) and the same raw residual?
:::


### B. Practice

::: {#exr-res-direct-t}
[B1]

Prove @thm-res-external-t(b) with \( \Delta=0 \) directly. Show that \( y_i-\x_{(i)}\T\hbeta_{(i)} \) is independent of \( s_{(i)}^2 \), that its
variance is \( \sigma^2\bigl(1+\x_{(i)}\T(\X_{(i)}\T\X_{(i)})^{-1}\x_{(i)}\bigr) \), and that this equals \( \sigma^2/(1-h_{ii}) \).
:::

::: {.solution}
The fit without case \( i \) is a normal linear model for the other \( n-1 \) cases, and \( y_i \) is independent of those cases.
So \( y_i-\x_{(i)}\T\hbeta_{(i)} \) is the error of predicting a new observation at \( \x_{(i)} \), and by @thm-ci-prediction-interval(a)
it is normal with mean zero and variance \( \sigma^2(1+g) \), \( g=\x_{(i)}\T(\X_{(i)}\T\X_{(i)})^{-1}\x_{(i)} \), independent of
\( s_{(i)}^2 \), with \( (y_i-\x_{(i)}\T\hbeta_{(i)})/(s_{(i)}\sqrt{1+g})\sim t(n-p-1) \). By @exm-mat-deletion,
\( (\X_{(i)}\T\X_{(i)})^{-1}=(\X\T\X)^{-1}+(\X\T\X)^{-1}\x_{(i)}\x_{(i)}\T(\X\T\X)^{-1}/(1-h_{ii}) \), so
\( g=h_{ii}+h_{ii}^2/(1-h_{ii})=h_{ii}/(1-h_{ii}) \) and \( 1+g=1/(1-h_{ii}) \). Using
@thm-res-deletion(d), the statistic is \( \bigl[\hat{\varepsilon}_i/(1-h_{ii})\bigr]\sqrt{1-h_{ii}}/s_{(i)}=t_i \).
:::

::: {#exr-res-deleted-variance}
[B2]

Show that \( s_{(i)}<s \) iff \( r_i^2>1 \). Interpret.
:::


::: {#exr-res-outlier-power}
[B3]

In the setting of @thm-res-external-t, the shift at case \( i \) is detected by the single-case test with probability
\( \Pr\{\lvert T\rvert>c\} \), \( T\sim t(n-p-1,\Delta\sqrt{1-h_{ii}}/\sigma) \). Show that for fixed \( \Delta \) and \( c \) this probability
decreases as \( h_{ii} \) increases. *Hint:* \( T^2 \) has a noncentral \( F \) law; use @thm-qf-f-power.
:::

### C. Going deeper

::: {#exr-res-heavy-tails-outlier}
[C1]

Suppose the errors are independent with a \( t \) distribution on \( 3 \) degrees of freedom and the design is that of the
state data. By simulation, estimate the probability that the \( 5\% \) Bonferroni outlier test rejects. Explain the result in
terms of the warning at the end of this section.
:::
