# Checking normality

Normality turns the first two moments of least squares into exact distributions. How much of Part III
needs it? [Chapter 19](../ch19-theory-of-departures/index.html)
answered the question for coefficients. When no observation carries much leverage, \( \hbeta \) is
approximately normal whatever the error law, and the \( t \) and \( F \) tests keep their nominal levels in
large samples (@thm-dep-nonnormal). Two kinds of inference are not protected in this way. The interval for
\( \sigma^2 \) of [Section 12.1](../ch12-intervals-and-bands/01-estimable-intervals.html) depends on the
fourth moment of the errors (@prp-opt-var-sse). A prediction interval is a statement about a single new
error, and no averaging makes a single error normal. So the error distribution is worth checking mainly when
the analysis will predict individual responses, when the sample is small, or when a few observations may be
of a different kind from the rest.

The errors are never observed, so the first task is to say how well the residuals stand in for them.

## Residuals are not the errors

If \( \Y=\X\bbeta+\be \) and \( \M \) is the hat matrix, the residual vector is
\( \he=(\I-\M)\be \) (@prp-proj-fit-residual). Write \( m_{ij} \) for the entries of \( \M \), so that \( h_{ii}=m_{ii} \) are the
leverages, and \( c_{ij} \) for the entries of \( \I-\M \).
Each residual is a weighted sum of all \( n \) errors,
\[
\hat{\varepsilon}_i=(1-h_{ii})\varepsilon_i-\sum_{j\ne i}m_{ij}\varepsilon_j ,
\]
so even when the errors are independent and identically distributed the residuals are neither. They are
correlated, and their variances \( \sigma^2(1-h_{ii}) \) differ.

::: {#prp-het-normal-check}
[What residuals reveal about the errors]

Let \( \varepsilon_1,\dots,\varepsilon_n \) be independent and identically distributed with mean \( 0 \) and
variance \( \sigma^2 \), and let \( \rank(\X)=p \).

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\varepsilon}_i-\varepsilon_i=-(\M\be)_i \), which has mean \( 0 \) and variance \( \sigma^2h_{ii} \). The
   total \( \E\norm{\he-\be}^2 \) equals \( \sigma^2p \).

2. Let \( \hat{\varepsilon}_{(1)}\le\dots\le\hat{\varepsilon}_{(n)} \) and
   \( \varepsilon_{(1)}\le\dots\le\varepsilon_{(n)} \) be the ordered residuals and ordered errors. Then
   \[
\frac1n\sum_{i=1}^n\bigl(\hat{\varepsilon}_{(i)}-\varepsilon_{(i)}\bigr)^2\le\frac1n\norm{\M\be}^2,
\]
   and the right-hand side has expectation \( \sigma^2p/n \).

3. Suppose the errors have skewness \( \gamma_1=\E\varepsilon^3/\sigma^3 \) and excess kurtosis
   \( \gamma_2=\E\varepsilon^4/\sigma^4-3 \). If \( h_{ii}<1 \), the standardized residual
   \( \hat{\varepsilon}_i/(\sigma\sqrt{1-h_{ii}}) \) has skewness \( a_i\gamma_1 \) and excess kurtosis
   \( b_i\gamma_2 \), where
   \[
a_i=\frac{\sum_jc_{ij}^3}{(1-h_{ii})^{3/2}},\qquad b_i=\frac{\sum_jc_{ij}^4}{(1-h_{ii})^2},
\]
   with \( \lvert a_i\rvert\le1 \) and \( (1-h_{ii})^2\le b_i\le1 \).
:::

:::

::: {.proof}
(a) \( \he-\be=-\M\be \), with mean zero and covariance \( \sigma^2\M\M\T=\sigma^2\M \) (@thm-rv-linear), whose
diagonal is \( \sigma^2h_{ii} \) and whose trace is \( \sigma^2p \).

(b) Let \( \mathbf{u},\mathbf{v}\in\Real^n \). Pairing the sorted values gives
\( \sum_iu_{(i)}v_{(i)}\ge\sum_iu_iv_{\pi(i)} \) for every permutation \( \pi \): if a pairing matches
\( u_{(i)}<u_{(j)} \) with \( v \)-values \( b>b' \) respectively, exchanging the two partners changes the sum by
\( (u_{(j)}-u_{(i)})(b-b')\ge0 \), and finitely many such exchanges reach the sorted pairing. Since sorting
does not change \( \norm{\mathbf{u}} \) or \( \norm{\mathbf{v}} \), it follows that
\( \sum_i(u_{(i)}-v_{(i)})^2\le\norm{\mathbf{u}-\mathbf{v}}^2 \). Take \( \mathbf{u}=\he \) and \( \mathbf{v}=\be \)
and use (a).

(c) Row \( i \) of \( \I-\M \) has \( \sum_jc_{ij}^2=(\I-\M)_{ii}=1-h_{ii} \) because \( \I-\M \) is symmetric and
idempotent. Expand \( \hat{\varepsilon}_i=\sum_jc_{ij}\varepsilon_j \). By independence and
\( \E\varepsilon_j=0 \), every cross term with an index appearing once vanishes, so \( \E\hat{\varepsilon}_i^3=\E\varepsilon^3\sum_jc_{ij}^3 \) and
\[
\begin{aligned}
\E\hat{\varepsilon}_i^4&=\E\varepsilon^4\sum_jc_{ij}^4+3\sigma^4\sum_{j\ne k}c_{ij}^2c_{ik}^2\\
&=(\E\varepsilon^4-3\sigma^4)\sum_jc_{ij}^4+3\sigma^4(1-h_{ii})^2 .
\end{aligned}
\]
Dividing by \( \sigma^3(1-h_{ii})^{3/2} \) and \( \sigma^4(1-h_{ii})^2 \) gives the skewness and the excess
kurtosis. For the bounds, \( \max_jc_{ij}^2\le\sum_jc_{ij}^2=1-h_{ii} \), so
\( \lvert\sum_jc_{ij}^3\rvert\le(1-h_{ii})^{1/2}(1-h_{ii}) \) and \( \sum_jc_{ij}^4\le(1-h_{ii})^2 \), while
\( \sum_jc_{ij}^4\ge c_{ii}^4=(1-h_{ii})^4 \).
:::

Parts (a) and (b) say that the residuals are good proxies for the errors when \( p/n \) is small. Part (b) is
the statement that matters for plots, because a normal probability plot shows sorted values: on average the
sorted residuals lie within \( \sigma\sqrt{p/n} \) of the sorted errors, in root mean square. Part (c)
describes what happens when \( p/n \) is not small. Each residual is an average of many errors, and averaging
pushes a distribution towards normality, so residuals look *more* normal than the errors that produced
them. This is sometimes called **supernormality**. The factor \( b_i \) measures it exactly.

::: {#exm-het-supernormality}
[How much nonnormality survives]

Take \( n=30 \) observations, an intercept, and \( p-1 \) columns of standard normal regressors. Averaged over
the observations, the share of the errors' excess kurtosis that survives in the standardized residuals is
\( 0.873 \) for \( p=2 \), \( 0.467 \) for \( p=10 \) and \( 0.162 \) for \( p=20 \). The share of skewness
is \( 0.898 \), \( 0.532 \) and \( 0.162 \). With centred exponential errors (skewness \( 2 \), excess
kurtosis \( 6 \)) and \( p=10 \), the residual with the largest leverage, \( h_{ii}=0.57 \), has predicted excess
kurtosis \( 1.23 \), against \( 1.26 \) in \( 200{,}000 \) simulated data sets. For the ordered values, the
average squared distance between sorted residuals and sorted errors was \( 0.133 \), while the bound in (b)
averaged \( 0.330 \), close to its expectation \( p/n=1/3 \).
:::

```{.python .run #cell-residual-errors-design}
import numpy as np
from scipy import stats

rng = np.random.default_rng(2101)
n = 30
designs = {}
for p in (2, 10, 20):
    X = np.column_stack([np.ones(n), rng.normal(size=(n, p - 1))])
    Q, _ = np.linalg.qr(X)
    M = Q @ Q.T                                      # the hat matrix
    designs[p] = (M, np.diag(M).copy())

for p, (M, h) in designs.items():
    # residual i = (row i of C) @ errors
    C = np.eye(n) - M
    # excess kurtosis of residual i / that of an error
    kurt_kept = (C**4).sum(axis=1) / (1 - h) ** 2
    # skewness of residual i / that of an error
    skew_kept = (C**3).sum(axis=1) / (1 - h) ** 1.5
    print(f"p = {p:2d}: max leverage {h.max():.2f}, "
          f"average share of skewness kept {skew_kept.mean():.3f}, "
          f"of excess kurtosis kept {kurt_kept.mean():.3f}")
```

## Normal probability plots

A **normal probability plot** graphs the sorted values of a sample against the corresponding quantiles of the
standard normal distribution. If \( U_1,\dots,U_n \) are independent \( \Normal(\mu,\tau^2) \), the \( i \)th
smallest is close to \( \mu+\tau q_i \), where \( q_i \) approximates the expected \( i \)th order statistic of a
standard normal sample. A common choice, due to Blom (1958), is
\[
q_i=\Phi^{-1}\Bigl(\frac{i-3/8}{n+1/4}\Bigr).
\]
Normal data give points near a straight line with intercept \( \mu \) and slope \( \tau \). Heavy tails bend
both ends away from the line, into an S shape. Skewness bends the plot into a curve that is convex for a
long right tail, and isolated points at either end suggest outliers.

For residuals, the unequal variances \( \sigma^2(1-h_{ii}) \) should be removed first. The natural
choice is the externally studentized residual \( t_i \), which divides \( \hat{\varepsilon}_i \) by
\( s_{(i)}\sqrt{1-h_{ii}} \), where \( s_{(i)}^2 \) is the error variance estimate without observation
\( i \) (@def-res-residuals). Under normal errors every \( t_i \) has the \( t(n-p-1) \)
distribution (@thm-res-external-t), so the \( t_i \) share a common marginal law and their sorted values can be compared
with one reference. They are not independent, and that makes it hard to judge by eye how much wiggle is
normal.

A simulated **envelope** answers that. Under normal errors the vector \( (t_1,\dots,t_n) \) is a function of
\( \he/\norm{\he}=(\I-\M)\be/\norm{(\I-\M)\be} \) alone, so its joint distribution depends on the design but not
on \( \bbeta \) or \( \sigma \) (@exr-het-envelope-free). We can therefore simulate it exactly: generate
normal errors, compute the sorted studentized residuals, repeat, and plot pointwise percentiles
(Atkinson 1981). The envelope is pointwise, so a few points outside it are expected even for normal data.
A run of consecutive points outside it, or an end of the plot leaving it, is the signal.

::: {#exm-het-normal-plots}
[Three error laws]

[Figure 21.1.1](#fig-het-normal-plots) shows plots for one simulated data set of \( n=60 \) observations with an
intercept and two regressors (\( p=3 \)), with three error laws of variance one. The normal case stays inside
its envelope except for \( 3 \) scattered points, and the Shapiro–Wilk test described below gives \( p=0.519 \).
Errors from a \( t(3) \) distribution give the S shape of heavy tails, with the largest residual far above
the envelope (Shapiro–Wilk \( p<0.001 \)). Centred exponential errors give a curved plot whose left end is too
short and whose right end is too long (\( p=0.001 \)).
:::

::: {when-format="html"}
![**Figure 21.1.1.** Normal probability plots of externally studentized residuals, \( n=60 \), \( p=3 \), with
pointwise 95% envelopes from 1999 data sets simulated with normal errors (shaded) and their pointwise
median (line). Points outside the envelope are drawn in the second colour. (a) Normal errors. (b) Errors
from \( t(3) \). (c) Centred exponential errors.](normal_plots.svg){#fig-het-normal-plots width=100%}
:::

::: {when-format="pdf"}
![Normal probability plots of externally studentized residuals, \( n=60 \), \( p=3 \), with
pointwise 95% envelopes from 1999 data sets simulated with normal errors (shaded) and their pointwise
median (line). Points outside the envelope are drawn in the second colour. (a) Normal errors. (b) Errors
from \( t(3) \). (c) Centred exponential errors.](normal_plots.pdf){width=100%}
:::

```{.python .run #cell-normal-plots-setup}
import numpy as np
from scipy import stats

rng = np.random.default_rng(2104)
n, p = 60, 3
X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
Q, _ = np.linalg.qr(X)
h = np.sum(Q**2, axis=1)                          # leverages


def ext_studentized(y):
    """Externally studentized residuals t_i = e_i / (s_(i) sqrt(1 - h_ii))."""
    e = y - Q @ (Q.T @ y)
    s2 = e @ e / (n - p)
    # squared internally studentized residuals
    r2 = e**2 / (s2 * (1 - h))
    return np.sign(e) * np.sqrt(r2 * (n - p - 1) / (n - p - r2))


# plotting positions
blom = stats.norm.ppf((np.arange(1, n + 1) - 0.375) / (n + 0.25))
sims = np.sort([ext_studentized(rng.normal(size=n)) for _ in range(1999)],
               axis=1)
# pointwise envelope
lower, upper = np.percentile(sims, [2.5, 97.5], axis=0)
```

```{.python .run #cell-normal-plots-plots}
laws = {
    "normal": rng.normal(size=n),
    "t(3)": rng.standard_t(3, size=n) / np.sqrt(3.0),
    "centred exponential": rng.exponential(size=n) - 1.0,
}
beta = np.array([1.0, 2.0, -1.0])
for name, eps in laws.items():
    t = np.sort(ext_studentized(X @ beta + eps))
    outside = np.sum((t < lower) | (t > upper))
    print(f"{name:20s}: {outside} of {n} points outside the envelope, "
          f"Shapiro-Wilk p = {stats.shapiro(t).pvalue:.3f}")
```

## Formal tests

Many tests of normality exist. Four kinds are in common use for residuals.

- The **Shapiro–Wilk** statistic (Shapiro and Wilk 1965) is
  \( W=\bigl(\sum_ia_iu_{(i)}\bigr)^2/\sum_i(u_i-\bar u)^2 \). The weights \( a_i \) come from the means and
  the covariance matrix of standard normal order statistics, and they make the numerator the square of the
  best linear unbiased estimate of the standard deviation from the ordered sample, up to a constant. So
  \( W \) compares two estimates of scale, one that is efficient only under normality and one that is always
  valid. It lies in \( (0,1] \), and small values reject. Its null distribution has no closed form and is
  computed numerically.
- The **Shapiro–Francia** statistic (Shapiro and Francia 1972) replaces the weights by the plotting
  positions, which makes it the squared correlation between the sorted sample and the \( q_i \). It is a
  one-number summary of how straight the normal probability plot is.
- **Moment tests** use the sample skewness \( \sqrt{b_1} \) and kurtosis \( b_2 \). The Jarque–Bera statistic
  \( n\{b_1/6+(b_2-3)^2/24\} \) is approximately \( \chi^2(2) \) for large normal samples (Jarque and Bera
  1980). These tests are directed at the two departures that matter most for inference.
- **Distance tests** compare the empirical and fitted normal distribution functions: Kolmogorov–Smirnov with
  estimated parameters (Lilliefors 1967) and Anderson–Darling (Anderson and Darling 1952). They are usually
  less powerful than Shapiro–Wilk.

Applied to residuals, none of these has its nominal null distribution, because the residuals are not a
sample. The effect on the level is small. In the design of @exm-het-supernormality with \( p=20 \) and
normal errors, the 5% Shapiro–Wilk test applied to the raw residuals rejected in \( 0.045 \) of \( 4000 \)
simulated data sets, and applied to the internally studentized residuals in \( 0.035 \). The effect on power is
not small. Supernormality hides the nonnormality that the test is looking for.

::: {#exm-het-sw-power}
[Power lost to supernormality]

In the same designs (\( n=30 \)), the 5% Shapiro–Wilk test applied to the unobservable centred exponential
errors rejects in \( 0.96 \) to \( 0.97 \) of simulated data sets. Applied to the residuals it rejects in \( 0.92 \) for
\( p=2 \), \( 0.45 \) for \( p=10 \) and only \( 0.07 \) for \( p=20 \). For errors from \( t(5) \), rescaled to
variance one, the rates are \( 0.24 \) to \( 0.25 \) for the errors and \( 0.23 \), \( 0.14 \) and \( 0.05 \) for the
residuals.
:::

```{.python .run #cell-residual-errors-power}
def errors(kind, size, rng):
    """Mean-zero, variance-one errors."""
    if kind == "exp":
        return rng.exponential(size=size) - 1.0
    return rng.standard_t(5, size=size) / np.sqrt(5 / 3)


def power_table(reps, rng):
    """Rejection rates of the 5% Shapiro-Wilk test on errors and on residuals."""
    table = {}
    for kind in ("exp", "t5"):
        for p, (M, h) in designs.items():
            E = errors(kind, (reps, n), rng)
            # residual vectors (M is symmetric)
            R = E - E @ M
            rej_err = np.mean([stats.shapiro(e).pvalue < 0.05 for e in E])
            rej_res = np.mean([stats.shapiro(r).pvalue < 0.05 for r in R])
            table[kind, p] = (rej_err, rej_res)
            print(f"{kind:3s} p = {p:2d}: Shapiro-Wilk rejects "
                  f"errors {rej_err:.2f}, "
                  f"residuals {rej_res:.2f}")
    return table


# a quick version (the book uses 4000 replications)
quick = power_table(400, rng)
```

A test of normality also answers the wrong question. In small samples it passes error laws that would
distort a prediction interval badly; in large ones it rejects departures too small to matter. The useful
question is whether the departure is large enough to affect the inference at hand. For coefficient tests
without high-leverage points it rarely is. For prediction intervals it can be, as the next result shows.

## What nonnormality does to prediction intervals

Under normal errors the interval \( \hat{Y}_0\pm t_{n-p,\alpha/2}\,s\sqrt{1+h_0} \) covers a new response
\( Y_0=\x_0\T\bbeta+\varepsilon_0 \) with probability exactly \( 1-\alpha \), where
\( h_0=\x_0\T(\X\T\X)^{-1}\x_0 \) (@thm-ci-prediction-interval). Without normality its coverage does not tend to \( 1-\alpha \) as the sample
grows. By @prp-ci-pi-limit, when \( h_0\to0 \) and the errors are independent and identically distributed
with a continuous distribution function \( F \), the coverage tends to
\( \Pr(\lvert\varepsilon_0\rvert\le z\sigma)=F(z\sigma)-F(-z\sigma) \), where \( z=\Phi^{-1}(1-\alpha/2) \), and
the probabilities of missing above and below tend to \( 1-F(z\sigma) \) and \( F(-z\sigma) \).

The limit is the probability that a single error falls within \( z \) standard deviations of zero, and it
depends on the whole shape of \( F \). [Section 12.3](../ch12-intervals-and-bands/03-prediction.html)
evaluated it for centred exponential and uniform errors at the 80% and 95% levels.
[Table 21.1.1](01-checking-normality.html#tab-het-prediction) extends that comparison to heavy-tailed laws and
to the 99% level, and splits the misses between the two tails.

[**Table 21.1.1.** Limiting coverage of the usual prediction interval (@prp-ci-pi-limit) when the errors are
not normal, with the limiting probabilities of missing below and above in parentheses. All five laws have
mean zero and variance one.]{#tab-het-prediction}

| error law | nominal 80% | nominal 95% | nominal 99% |
|---|---|---|---|
| normal | 0.800 (0.100, 0.100) | 0.950 (0.025, 0.025) | 0.990 (0.005, 0.005) |
| \( t(5) \) | 0.841 (0.079, 0.079) | 0.947 (0.026, 0.026) | 0.979 (0.010, 0.010) |
| \( t(3) \) | 0.887 (0.057, 0.057) | 0.957 (0.021, 0.021) | 0.979 (0.010, 0.010) |
| uniform | 0.740 (0.130, 0.130) | 1.000 (0.000, 0.000) | 1.000 (0.000, 0.000) |
| centred exponential | 0.898 (0.000, 0.102) | 0.948 (0.000, 0.052) | 0.972 (0.000, 0.028) |

Heavy tails make a 99% interval too short and an 80% interval too long. The 95% level lies near the
crossing point, so checking a 95% interval alone would miss the problem. For the skewed law every miss is
above the interval, since a centred exponential error is never below \( -\sigma \); at the 99% level it misses
above \( 0.028 \) of the time, more than five times the nominal \( 0.005 \) in that tail.

```{.python .run #cell-prediction-coverage-limits}
import numpy as np
from scipy import stats

laws = {
    "normal": stats.norm(),
    "t(5)": stats.t(5, scale=np.sqrt(3 / 5)),
    "t(3)": stats.t(3, scale=np.sqrt(1 / 3)),
    "uniform": stats.uniform(loc=-np.sqrt(3), scale=2 * np.sqrt(3)),
    "centred exponential": stats.expon(loc=-1.0),
}
for name, law in laws.items():
    assert abs(law.var() - 1) < 1e-9 and abs(law.mean()) < 1e-9
    row = []
    for level in (0.80, 0.95, 0.99):
        z = stats.norm.ppf(0.5 + level / 2)
        below, above = law.cdf(-z), law.sf(z)            # miss low, miss high
        row.append(f"{1 - below - above:.3f} ({below:.3f}, {above:.3f})")
    print(f"{name:20s}", "  ".join(row))
```

The remedy is to use the residuals for what they are good at, estimating the error distribution:
empirical quantiles of the residuals in place of \( \pm z\sigma \) give intervals whose coverage tends to the
nominal level, and [Chapter 23](../ch23-resampling-inference/index.html) refines this by resampling. A
transformation of the response ([Chapter 22](../ch22-transformations/index.html)) can remove skewness at its
source.

## Exercises

### A. Check your understanding

::: {#exr-het-residual-variance}
[A1]

In simple regression with \( n=5 \) and \( x=(1,2,3,4,10) \), compute the leverages and the variances of the five
residuals when \( \sigma^2=1 \). Which residual is least like its error, in the sense of
@prp-het-normal-check(a)?
:::

### B. Practice

::: {#exr-het-envelope-free}
[B1]

Assume \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \). Show that the internally and externally studentized residuals
are functions of \( \mathbf{g}=\he/\norm{\he} \), and that the law of \( \mathbf{g} \) does not depend on \( \bbeta \) or
\( \sigma \). Conclude that the simulated envelope of this section is exact, whatever \( \bbeta \) and \( \sigma \)
are used in the simulation.
:::

::: {.solution}
\( \he=(\I-\M)\be \), and \( \be/\sigma\sim\Normal_n(\bzero,\I) \) whatever \( \bbeta \) and \( \sigma \) are, so
\( \mathbf{g}=(\I-\M)(\be/\sigma)/\norm{(\I-\M)(\be/\sigma)} \) has a fixed law. The internally studentized
residual is \( r_i=\hat{\varepsilon}_i/(s\sqrt{1-h_{ii}}) \) with \( s=\norm{\he}/\sqrt{n-p} \), so
\( r_i=\sqrt{n-p}\,g_i/\sqrt{1-h_{ii}} \). The externally studentized residual is the monotone function
\( t_i=r_i\sqrt{(n-p-1)/(n-p-r_i^2)} \) of \( r_i \) (@def-res-residuals). Hence the whole vector of sorted
\( t_i \) has a law that depends only on \( \M \), and it can be simulated with any \( \bbeta \) and \( \sigma \).
:::

::: {#exr-het-exponential-coverage}
[B2]

For centred exponential errors, \( \varepsilon=E-1 \) with \( E \) standard exponential, show that the limit
in @prp-ci-pi-limit at the 95% level is \( 1-e^{-2.96} \), and find the probabilities of missing below and
above.
:::

::: {.solution}
Here \( \sigma=1 \) and \( z=1.96 \). The event \( \lvert E-1\rvert\le1.96 \) is \( E\le2.96 \), since \( E\ge0>-0.96 \)
always. Its probability is \( 1-e^{-2.96}=1-0.0518=0.948 \). The interval never misses below, and it
misses above with probability \( e^{-2.96}=0.0518 \).
:::

::: {#exr-het-oneway-kurtosis}
[B3]

In a one-way layout in which every group has \( m \) observations, show that the kurtosis factor of
@prp-het-normal-check(c) is the same for every residual and equals
\[
b=\frac{(1-1/m)^4+(m-1)/m^4}{(1-1/m)^2}.
\]
Evaluate it for \( m=2 \) and \( m=5 \), and explain why residual plots from small groups say little about the
tails of the errors.
:::

### C. Going deeper

::: {#exr-het-ecdf}
[C1]

Under the assumptions of @prp-het-normal-check, suppose \( p/n\to0 \) along a sequence of designs. Let
\( \hat F_n \) and \( F_n \) be the empirical distribution functions of the residuals and of the errors. Show that
\( \hat F_n(x)-F_n(x)\to0 \) in probability at every continuity point \( x \) of \( F \).
:::

