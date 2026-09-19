# Heterogeneous slopes

The covariance model of [Section 18.1](01-ancova.html) assumes a common slope, so each treatment comparison is a
vertical shift between parallel lines, one number. If the slopes differ, the treatment effect depends on the covariate,
and "the" adjusted difference no longer exists. This section tests parallelism and says what to report when it fails.

## The separate lines model and the test of parallelism

Keep one factor with \( g \) levels and one covariate, and allow each level its own line:
\[
\E(y_i)=\mu_k+\beta_kx_i\qquad\text{for unit } i \text{ in group } k .
\]
This is the separate lines model @eq-est-separate of [Section 8.6](../ch08-estimability/06-several-factors.html). The
parallel lines model is the special case \( \beta_1=\dots=\beta_g \). For group \( k \) write
\[
E_{xx,k}=\sum_{i\in k}(x_i-\bar{x}_k)^2,\qquad E_{xy,k}=\sum_{i\in k}(x_i-\bar{x}_k)(y_i-\bar{y}_k),
\]
so that the pooled within-group quantities of [Section 18.1](01-ancova.html) are \( E_{xx}=\sum_kE_{xx,k} \) and
\( E_{xy}=\sum_kE_{xy,k} \).

::: {#thm-dsn-slopes}
[The test of parallel slopes]

Assume \( x \) varies within every group, so that every \( E_{xx,k}>0 \), and \( n>2g \).

::: {.enumerate options="label=(\alph*)"}
1. In the separate lines model, \( \hat{\beta}_k=E_{xy,k}/E_{xx,k} \), these estimates are uncorrelated with variances
   \( \sigma^2/E_{xx,k} \), and the residual sum of squares is
   \[
   \text{SSE}_{\text{sep}}=\sum_k\Bigl(E_{yy,k}-\frac{E_{xy,k}^2}{E_{xx,k}}\Bigr)
   \]
   on \( n-2g \) degrees of freedom.

2. The pooled slope of the parallel lines model is the weighted mean
   \( \hat{\beta}=\sum_kE_{xx,k}\hat{\beta}_k/\sum_kE_{xx,k} \), and the extra sum of squares for the separate slopes is
   \[
   \text{SS}_{\text{slopes}}=\sum_kE_{xx,k}\,\bigl(\hat{\beta}_k-\hat{\beta}\bigr)^2
   =\sum_k\frac{E_{xy,k}^2}{E_{xx,k}}-\frac{E_{xy}^2}{E_{xx}}
   \]{#eq-dsn-slopes-ss}

   on \( g-1 \) degrees of freedom.

3. If \( \Y \) is normal,
   \[
   F=\frac{\text{SS}_{\text{slopes}}/(g-1)}{\text{SSE}_{\text{sep}}/(n-2g)}\sim F\bigl(g-1,\,n-2g,\,\gamma\bigr),
   \]
   with
   \[
   \gamma=\frac1{\sigma^2}\sum_kE_{xx,k}\,(\beta_k-\bar{\beta})^2 ,
   \]
   where \( \bar{\beta}=\sum_kE_{xx,k}\beta_k/\sum_kE_{xx,k} \). The test that rejects for large \( F \) is the \( F \) test of
   parallelism.
:::

:::

::: {.proof}
Let \( \M \) project onto the span of the group indicators, and let \( \bw_k \) be the vector that equals \( x_i-\bar{x}_k \) for
units of group \( k \) and zero elsewhere. The \( k \)th slope column of the separate lines model is \( \x \) restricted to
group \( k \), and its residual after projection onto \( \C(\M) \) is \( \bw_k \). By @lem-ss-residualized the separate lines
model space is \( \C(\M)\dirsum\spn\{\bw_1,\dots,\bw_g\} \), and the parallel lines space is
\( \C(\M)\dirsum\spn\{\bw\} \) with \( \bw=\sum_k\bw_k=(\I-\M)\x \). The vectors \( \bw_k \) have disjoint supports, so they are
mutually orthogonal, with \( \norm{\bw_k}^2=E_{xx,k} \) and \( \bw_k\T\y=E_{xy,k} \).

(a) The projection of \( \y \) onto \( \spn\{\bw_k\} \) is \( \sum_k(E_{xy,k}/E_{xx,k})\bw_k \). Reading off the coefficient of
the \( k \)th slope column gives \( \hat{\beta}_k \), which is \( \bw_k\T\Y/\norm{\bw_k}^2 \). Hence
\( \Cov(\hat{\beta}_k,\hat{\beta}_l)=\sigma^2\bw_k\T\bw_l/(E_{xx,k}E_{xx,l}) \), which is \( \sigma^2/E_{xx,k} \) for \( k=l \) and zero
otherwise. The residual sum of squares is
\( \norm{(\I-\M)\y}^2-\sum_kE_{xy,k}^2/E_{xx,k} \), and \( \norm{(\I-\M)\y}^2=\sum_kE_{yy,k} \). The model has rank \( 2g \).

(b) The pooled slope is \( \bw\T\y/\norm{\bw}^2=E_{xy}/E_{xx} \), and \( E_{xy}=\sum_kE_{xx,k}\hat{\beta}_k \). By @thm-proj-nested the extra sum of squares is the squared length of the projection onto the orthogonal complement of
\( \bw \) in \( \spn\{\bw_k\} \), which equals \( \sum_kE_{xy,k}^2/E_{xx,k}-E_{xy}^2/E_{xx} \). Writing
\( E_{xy,k}=E_{xx,k}\hat{\beta}_k \) and expanding
\( \sum_kE_{xx,k}(\hat{\beta}_k-\hat{\beta})^2=\sum_kE_{xx,k}\hat{\beta}_k^2-\hat{\beta}^2E_{xx} \) gives the first form.

(c) By @thm-glh-f-test, with \( \sigma^2\gamma \) equal to the extra sum of squares evaluated at the mean vector. The mean
vector has \( \bw_k\T\bmu=\beta_kE_{xx,k} \), because \( \bw_k \) is orthogonal to the group indicators. Substituting into (b) with \( \hat{\beta}_k \) replaced by \( \beta_k \) gives \( \gamma \).
:::

The test is a weighted test of homogeneity of \( g \) independent slope estimates, with weights equal to their
precisions \( E_{xx,k}/\sigma^2 \).

::: {#exm-dsn-tutoring-slopes}
[Parallelism in the tutoring experiment]

For the four formats of @exm-dsn-tutoring, the separate slopes are \( 0.703 \), \( 0.823 \),
\( 1.035 \) and \( 0.992 \), and their weighted mean is the pooled slope \( 0.8670 \). The separate lines model
has residual sum of squares \( 1347.47 \) on 32 degrees of freedom, against \( 1407.47 \) on 35 for parallel
lines. The test gives \( F=0.47 \) on 3 and 32 degrees of freedom, \( p=0.70 \), and the covariance analysis
stands.
:::

```{.python .run #cell-tutoring-parallel}
import numpy as np
from scipy import stats

rng = np.random.default_rng(84)
g, m = 4, 10                                          # four formats, ten students each
n = g * m
group = rng.permutation(np.repeat(np.arange(g), m))   # the randomization
x = np.round(rng.normal(60, 10, n))                   # pretest, measured before assignment
mu_true = np.array([70.0, 74.0, 71.0, 77.0])
y = np.round(mu_true[group] + 0.8 * (x - 60) + rng.normal(0, 6, n), 1)
Z = np.eye(g)[group]                                  # indicator matrix of the formats
n_k = Z.sum(axis=0)

def within(v):
    """(I - M) v: deviations from the group means (M projects onto C(Z))."""
    return v - (Z.T @ v / n_k)[group]

E_xx, E_xy, E_yy = within(x) @ within(x), within(x) @ within(y), within(y) @ within(y)
sse = E_yy - E_xy ** 2 / E_xx                         # ANCOVA residual sum of squares
X_sep = np.column_stack([Z, Z * x[:, None]])          # a separate line for each format
E_xx_k = np.array([within(x)[group == k] @ within(x)[group == k] for k in range(g)])
E_xy_k = np.array([within(x)[group == k] @ within(y)[group == k] for k in range(g)])
slopes = E_xy_k / E_xx_k
sse_sep = sse - (np.sum(E_xy_k ** 2 / E_xx_k) - E_xy ** 2 / E_xx)
F_par = ((sse - sse_sep) / (g - 1)) / (sse_sep / (n - 2 * g))
print("separate slopes", slopes.round(3),
      f" F = {F_par:.2f}, p = {stats.f.sf(F_par, g - 1, n - 2 * g):.3f}")
```

With ten students per group the test has little power, so the case for parallel lines rests mostly on the
science. When a treatment is meant to help some units more than others, such as a remedial programme aimed at weak
students, the separate lines model should be the starting point rather than a check.

## When the slopes differ

If the slopes differ, the difference between treatments \( k \) and \( l \) at covariate value \( x \) is the line
\[
\delta(x)=(\mu_k-\mu_l)+(\beta_k-\beta_l)\,x ,
\]
estimated by \( \hat{\delta}(x)=(\hat{\mu}_k-\hat{\mu}_l)+(\hat{\beta}_k-\hat{\beta}_l)x \). The two lines are independent simple
regressions, so by @cor-lm-simple-moments
\[
\Var\bigl(\hat{\delta}(x)\bigr)=\sigma^2v(x),\qquad
v(x)=\frac1{n_k}+\frac{(x-\bar{x}_k)^2}{E_{xx,k}}+\frac1{n_l}+\frac{(x-\bar{x}_l)^2}{E_{xx,l}} .
\]{#eq-dsn-delta-var}

Honest reports include:

::: {.enumerate options="label=(\roman*)"}
1. Report \( \hat{\delta}(x) \) with its standard error at a few covariate values chosen for their meaning, for instance the
   quartiles of \( x \), or a clinical threshold. Each is a \( t \) interval (@thm-ci-estimable-interval) with \( n-2g \) degrees of
   freedom.

2. Draw \( \hat{\delta}(x) \) as a function with a confidence band. The pointwise band uses \( t_{n-2g,\alpha/2} \). A band that
   covers the whole line \( \delta(\cdot) \) with probability exactly \( 1-\alpha \) replaces this multiplier by
   \( \sqrt{2F_\alpha(2,n-2g)} \), because \( \{\delta(x):x\in\Real\} \) lies in a two-dimensional space of estimable functions.
   This is the argument of @thm-mc-scheffe, as used for the Working–Hotelling band (@thm-ci-working-hotelling).

3. Give the set of covariate values at which the treatments differ significantly. This is the **Johnson–Neyman
   region** (@prp-dsn-johnson-neyman below).
:::

Reporting \( \hat{\delta}(\bar{x}) \) alone as *the* treatment effect hides the slope, which is the whole point.

::: {.remark}
[The average effect in a randomized experiment]

There is one exception. In a randomized experiment, \( \delta \) at the population mean of the covariate is the *average*
treatment effect, because \( \delta \) is linear. Its estimate \( \hat{\delta}(\bar{x}) \) is the treatment coefficient in the
regression with a treatment-by-covariate interaction and the covariate centred at \( \bar{x} \) (@exr-dsn-interacted). This
"interacted" adjustment is asymptotically never less precise than the unadjusted difference, even when the model is wrong.
[Chapter 25](../ch25-causal-interpretation/index.html) returns to it (@thm-cau-adjustment).
:::

## The Johnson–Neyman region

::: {#prp-dsn-johnson-neyman}
[Johnson–Neyman region]

In the separate lines model with normal errors, let \( s^2=\text{SSE}_{\text{sep}}/(n-2g) \) and fix a multiplier \( c>0 \). Write
\( \Delta_\mu=\hat{\mu}_k-\hat{\mu}_l \) and \( \Delta_\beta=\hat{\beta}_k-\hat{\beta}_l \). The set of \( x \) with
\( |\hat{\delta}(x)|>c\,s\sqrt{v(x)} \) is \( \{x:Ax^2+Bx+C>0\} \), where
\[
\begin{aligned}
A&=\Delta_\beta^2-c^2s^2\Bigl(\frac1{E_{xx,k}}+\frac1{E_{xx,l}}\Bigr),\qquad
B=2\Delta_\mu\Delta_\beta+2c^2s^2\Bigl(\frac{\bar{x}_k}{E_{xx,k}}+\frac{\bar{x}_l}{E_{xx,l}}\Bigr),\\
C&=\Delta_\mu^2-c^2s^2\Bigl(\frac1{n_k}+\frac1{n_l}+\frac{\bar{x}_k^2}{E_{xx,k}}+\frac{\bar{x}_l^2}{E_{xx,l}}\Bigr).
\end{aligned}
\]
If \( A>0 \), which happens iff the \( t \) test of \( \beta_k=\beta_l \) with multiplier \( c \) rejects, the region is
the whole line with at most an interval removed: \( \delta \) is significantly positive on one side and significantly
negative on the other. If \( A<0 \), the region is a bounded interval or empty. With \( c=t_{n-2g,\alpha/2} \) each point of
the region is significant at level \( \alpha \). With \( c=\sqrt{2F_\alpha(2,n-2g)} \) the statement "\( \delta(x)\neq0 \) for every
\( x \) in the region" holds with probability at least \( 1-\alpha \).
:::

::: {.proof}
Square both sides of \( |\hat{\delta}(x)|>c\,s\sqrt{v(x)} \) and expand \( \hat{\delta}(x)^2=\Delta_\mu^2+2\Delta_\mu\Delta_\beta x+\Delta_\beta^2x^2 \) and
\( v(x) \) from @eq-dsn-delta-var as polynomials in \( x \). Collecting terms gives \( Ax^2+Bx+C \). The coefficient \( A \) is
positive iff \( \Delta_\beta^2/\{s^2(1/E_{xx,k}+1/E_{xx,l})\}>c^2 \), which is the squared \( t \) statistic for \( \beta_k-\beta_l \) exceeding
\( c^2 \). A quadratic with \( A>0 \) is positive outside the interval between its real roots (everywhere if it has none), and with
\( A<0 \) it is positive between its real roots (nowhere if it has none). The last statement follows from the simultaneous band of
item (ii) above: on the event that the band covers \( \delta(\cdot) \) everywhere, \( \delta(x) \) cannot be zero at any point where the
band excludes zero.
:::

::: {#exm-dsn-reading}
[Two versions of a reading programme]

In a second synthetic randomized experiment, \( 48 \) children with baseline reading fluency between \( 42 \) and
\( 118 \) words per minute are assigned, 24 each, to the current and a new version of a reading programme. The
response is the gain in fluency. The new version was designed for weak readers, and the fitted lines of
[Figure 18.3.1](#fig-dsn-johnson-neyman)(a) show why that matters. The slope is \( 0.0156 \) for the current version and
\( -0.1615 \) for the new one. The test of parallelism has \( F=11.2 \) on 1 and 44 degrees of freedom,
\( p=0.0017 \). The pooled slope of the parallel lines model would be \( -0.0638 \), a compromise that describes neither
group.

The estimated lines cross at \( x=107.7 \). At a baseline of 60 words per minute the new version adds
\( 8.46 \) words per minute (standard error \( 1.56 \)). At 100 it adds \( 1.37 \) (standard error
\( 1.46 \)), which is not significant. With \( c=t_{44,0.025} \), the Johnson–Neyman region (@prp-dsn-johnson-neyman) is \( x<93.5 \) together with
\( x>152.0 \). The second part lies far beyond the largest baseline in the data, and it is an artefact of extending two
fitted lines. The useful conclusion is that the new version is better for children who start below about 93 words per
minute, which covers \( 33 \) of the \( 48 \) children in the study. With the simultaneous multiplier
\( \sqrt{2F_{0.05}(2,44)} \), the region shrinks slightly, to \( x<90.8 \) (and \( x>195.6 \)).
:::

::: {when-format="html"}
![**Figure 18.3.1.** The reading programme. (a) Gain against baseline fluency with the separate least squares lines. (b) The
estimated difference \( \hat{\delta}(x) \) (new minus current) with a pointwise 95% band (dashed) and a simultaneous 95% band
(shaded). The dotted line marks the lower end of the Johnson–Neyman region, \( x=93.5 \); to its left the
difference is significantly positive.](johnson_neyman.svg){#fig-dsn-johnson-neyman width=100%}
:::

::: {when-format="pdf"}
![The reading programme. (a) Gain against baseline fluency with the separate least squares lines. (b) The
estimated difference \( \hat{\delta}(x) \) (new minus current) with a pointwise 95% band (dashed) and a simultaneous 95% band
(shaded). The dotted line marks the lower end of the Johnson–Neyman region, \( x=93.5 \); to its left the
difference is significantly positive.](johnson_neyman.pdf){width=100%}
:::

```{.python .run #cell-slopes-jn}
import numpy as np
from scipy import stats

rng = np.random.default_rng(31)
n_per = 24
arm = rng.permutation(np.repeat([0, 1], n_per))      # 0 = current version, 1 = new version
x = np.round(rng.uniform(40, 120, 2 * n_per))        # baseline fluency
y = np.where(arm == 0, 14 + 0.06 * x, 36 - 0.16 * x) + rng.normal(0, 4, 2 * n_per)
y = np.round(y, 1)
n = len(y)

def line(k):
    """Least squares line within arm k, and the arm's summaries."""
    xs, ys = x[arm == k], y[arm == k]
    Sxx = np.sum((xs - xs.mean()) ** 2)
    b = np.sum((xs - xs.mean()) * (ys - ys.mean())) / Sxx
    a = ys.mean() - b * xs.mean()
    sse = np.sum((ys - a - b * xs) ** 2)
    return a, b, xs.mean(), Sxx, len(xs), sse

(a0, b0, xb0, S0, n0, e0), (a1, b1, xb1, S1, n1, e1) = line(0), line(1)
nu = n - 4
s2 = (e0 + e1) / nu                                   # separate-lines model

def significant_region(crit2):
    """Roots of d(x)^2 - crit2 * s2 * v(x) = 0, a quadratic in x."""
    A = (b1 - b0) ** 2 - crit2 * s2 * (1 / S0 + 1 / S1)
    B = 2 * (a1 - a0) * (b1 - b0) + 2 * crit2 * s2 * (xb0 / S0 + xb1 / S1)
    C = (a1 - a0) ** 2 - crit2 * s2 * (1 / n0 + 1 / n1 + xb0 ** 2 / S0 + xb1 ** 2 / S1)
    return np.sort(np.roots([A, B, C]).real)

t2 = stats.t.ppf(0.975, nu) ** 2                      # pointwise, 5% level
w2 = 2 * stats.f.ppf(0.95, 2, nu)                     # simultaneous over all x (Scheffe, q = 2)
x_cross = -(a1 - a0) / (b1 - b0)
print(f"lines cross at x = {x_cross:.1f}")
print("pointwise 5% region: difference significant outside", significant_region(t2).round(1))
print("simultaneous region: difference significant outside", significant_region(w2).round(1))
```

The region is only as good as the linearity of each line, it extrapolates outside the data, and it marks
significance, not practical importance: with enough data every value except the crossing point becomes significant.

## Several covariates, several factors

With \( q \) covariates, the test of a common \( \bgamma \) against \( \bgamma_1,\dots,\bgamma_g \) has \( (g-1)q \) and
\( n-g(q+1) \) degrees of freedom. With two factors, slope terms can be built term by term like the means model of
[Chapter 16](../ch16-multiway-layouts/index.html). Each such test compares nested column spaces (@thm-glh-f-test).

## Exercises

### A. Check your understanding

::: {#exr-dsn-two-slopes}
[A1]

For \( g=2 \), show that @eq-dsn-slopes-ss reduces to
\( \text{SS}_{\text{slopes}}=\dfrac{E_{xx,1}E_{xx,2}}{E_{xx,1}+E_{xx,2}}\,(\hat{\beta}_1-\hat{\beta}_2)^2 \), and that the \( F \) statistic is the
square of the two-sample \( t \) statistic for the difference of slopes.
:::

::: {#exr-dsn-slopes-crossing}
[A2]

In @exm-dsn-reading, the estimated lines cross at \( x=107.7 \). Explain why the Johnson–Neyman region never contains the
crossing point, whatever the multiplier.
:::

### B. Practice

::: {#exr-dsn-jn-bounded}
[B1]

Give a numerical example of two groups for which \( A<0 \) and the Johnson–Neyman region of @prp-dsn-johnson-neyman is a bounded interval. Interpret the
region.
:::

::: {#exr-dsn-slopes-noncentrality}
[B2]

In the tutoring experiment, suppose the true slopes were \( 0.8,0.8,0.8,1.2 \) and \( \sigma=6 \). Using the \( E_{xx,k} \) of the
example (their sum is \( 3643.9 \); take them equal for the calculation), compute the noncentrality of the test of parallelism and its
power at level \( 0.05 \) (@thm-glh-power).
:::

::: {#exr-dsn-band-derivation}
[B3]

Prove the claim in item (ii): with \( c=\sqrt{2F_\alpha(2,n-2g)} \), the band \( \hat{\delta}(x)\pm c\,s\sqrt{v(x)} \) contains \( \delta(x) \)
for all \( x \) simultaneously with probability at least \( 1-\alpha \). Is the probability exactly \( 1-\alpha \)?
:::

::: {.solution}
The functions \( \delta(x)=\blambda_x\T\bbeta \) for \( x\in\Real \) have \( \blambda_x \) in a two-dimensional space \( \mathcal L \) spanned by the
coefficient vectors of \( \mu_k-\mu_l \) and \( \beta_k-\beta_l \), both estimable. By @thm-mc-scheffe(b) with \( q=2 \), with probability
exactly \( 1-\alpha \) *every* \( \blambda\in\mathcal L \) satisfies \( |\blambda\T\hbeta-\blambda\T\bbeta|\le c\,\text{se}(\blambda\T\hbeta) \). The
functions \( \blambda_x \) form a subset of \( \mathcal L \), so the probability for them is at least \( 1-\alpha \). It is in fact exactly
\( 1-\alpha \). In the coordinates \( (\mu_k-\mu_l,\beta_k-\beta_l) \) we have \( \blambda_x=(1,x) \), and the studentized ratio
\( (\blambda\T\hbeta-\blambda\T\bbeta)^2/\text{se}(\blambda\T\hbeta)^2 \) is unchanged when \( \blambda \) is rescaled and is continuous on
\( \mathcal L\setminus\{\bzero\} \). Every direction of \( \mathcal L \) except that of \( (0,1) \) is the direction of some \( \blambda_x \), and
\( \blambda_x/\norm{\blambda_x}\to\pm(0,1) \) as \( x\to\pm\infty \). So the supremum of the ratio over \( x\in\Real \) equals its supremum over
\( \mathcal L\setminus\{\bzero\} \), and the two events coincide. This is the argument of @thm-ci-working-hotelling(b) for a single line.
:::

### C. Going deeper

::: {#exr-dsn-interacted}
[C1]

In a two-group randomized experiment, show that the treatment coefficient in the least squares fit of
\( y \) on \( 1 \), \( T \), \( x-\bar{x} \) and \( T(x-\bar{x}) \), where \( \bar{x} \) is the overall mean, equals
\( \bigl(\bar{y}_1-\hat{\beta}_1(\bar{x}_1-\bar{x})\bigr)-\bigl(\bar{y}_0-\hat{\beta}_0(\bar{x}_0-\bar{x})\bigr) \), with \( \hat{\beta}_k \) the within-group slopes.
Compare with the adjusted difference of the parallel lines model.
:::

::: {.solution}
The model with an interaction is the separate lines model, so its fitted lines are the two within-group least squares lines,
\( \bar{y}_k+\hat{\beta}_k(x-\bar{x}_k) \). In the given parameterization the line for group \( k \) is
\( (\beta_0+\beta_TT)+(\beta_x+\beta_{Tx}T)(x-\bar{x}) \), so the coefficient of \( T \) is the difference of the two lines at
\( x=\bar{x} \), which is the stated expression. The parallel lines model replaces both \( \hat{\beta}_k \) by the pooled slope. The two
estimates agree when the slopes are equal or when \( \bar{x}_0=\bar{x}_1=\bar{x} \).
:::
