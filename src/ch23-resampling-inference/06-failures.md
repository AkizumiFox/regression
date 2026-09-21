# When resampling fails

The consistency theorems of this chapter assume exchangeable resampled quantities, a smooth statistic, a fixed number of
parameters, a finite error variance, no dominant case, independence, and a sample large enough for limits to be relevant. When any of these fails, the bootstrap can
fail silently, producing a distribution that looks reasonable and is wrong.

## A statistic that is not smooth: the maximum

The cleanest failure is not specific to regression. Let \( Y_1,\dots,Y_n \) be independent and uniform on
\( [0,\theta] \), and estimate \( \theta \) by the sample maximum \( Y_{(n)} \).

::: {#prp-bs-maximum}
[The bootstrap of a maximum]

::: {.enumerate options="label=(\alph*)"}
1. \( n(\theta-Y_{(n)})/\theta \) converges in distribution to the standard exponential law.

2. For every sample with distinct values, the case bootstrap gives
   \( {\Pr}_*(Y^*_{(n)}=Y_{(n)})=1-(1-1/n)^n\to1-e^{-1} \).
:::

So the bootstrap distribution of \( n(Y_{(n)}-Y^*_{(n)})/Y_{(n)} \) has an atom at zero of mass at least
\( 1-e^{-1} \) for every sample, and it cannot approach the continuous exponential law.
:::

::: {.proof}
(a) For \( 0\le t\le n \), \( \Pr\bigl(n(\theta-Y_{(n)})/\theta>t\bigr)=\Pr(\text{all }Y_i<\theta(1-t/n))=(1-t/n)^n\to e^{-t} \).
(b) \( Y^*_{(n)}=Y_{(n)} \) unless none of the \( n \) draws picks the largest observation, which happens with probability
\( (1-1/n)^n \). And \( (1-1/n)^n=\exp\bigl(n\log(1-1/n)\bigr)\to e^{-1} \), from below, since \( n\log(1-1/n)<-1 \).
:::

With \( n=50 \) the atom is \( 1-0.98^{50}=0.636 \), and a simulation with \( 20000 \) resamples finds
\( 0.639 \). The bootstrap needs the statistic to depend smoothly on the empirical distribution, and extremes do not.
The \( m \)-out-of-\( n \) bootstrap, with resamples of size \( m \) and \( m/n\to0 \), repairs many such failures at a price in
efficiency (Bickel, Götze and van Zwet 1997). In regression the analogue is any statistic driven by the most extreme residual,
such as the largest studentized residual of the outlier test (@thm-res-external-t).

```{.python .run #cell-failures-maximum}
import numpy as np
rng = np.random.default_rng(2310)
n = 50
u = rng.uniform(size=n)
B = 20_000
boot_max = u[rng.integers(0, n, size=(B, n))].max(axis=1)
atom = np.mean(boot_max == u.max())
print(f"bootstrap P*(max* = max) = {atom:.3f}; 1 - (1 - 1/n)^n = {1 - (1 - 1 / n) ** n:.3f}")
```

## Heavy tails

Without a finite error variance the sample mean, and any least squares coefficient, is not asymptotically normal. Athreya (1987)
showed that the bootstrap of the mean then converges not to the correct stable limit but to a *random* limit that differs from
sample to sample: resampling cannot learn a tail that shows up only through a few enormous values.

[Figure 23.6.1](#fig-bs-failures)(b) shows the consequence for coverage. The data are samples from a Pareto law with
tail index \( 1.5 \), which has a finite mean and an infinite variance, and the target is the mean. The nominal \( 95\% \)
\( t \) interval covers \( 0.665 \) of the time at \( n=25 \) and \( 0.758 \) at
\( n=1600 \); the bootstrap percentile interval covers \( 0.673 \) and \( 0.780 \). Neither improves between
\( n=400 \) and \( n=1600 \). The same holds for regression coefficients, which is one reason to prefer estimators less
sensitive to extreme errors, such as the least absolute deviations fit of [Chapter 45](../ch45-quantile-gamlss/index.html) (@def-qnt-quantile), to resampling least squares.

## One dominant case

When one case has leverage near one, @eq-bs-no-dominant fails and the case bootstrap splits into resamples that contain the case
and resamples, about a third of them, that do not.

::: {#exm-bs-leverage-failure}
[A slope determined by one case]

Take \( 15 \) cases, fourteen with \( x \) spread over \( [0,1] \) and one at \( x=4 \), with leverage
\( 0.902 \), and a response that puts the far case a little below the line of the others. The least squares slope is
\( 0.254 \). Without the far case it would be \( 0.437 \). In \( 20000 \) case-bootstrap resamples,
a fraction \( 0.356 \) omit the far case, close to \( (14/15)^{15}=0.355 \).
[Figure 23.6.1](#fig-bs-failures)(a) shows the result. The resamples that contain the far case give slopes tightly
concentrated near the full-data value, with standard deviation \( 0.037 \). Those that omit it give slopes spread
around the fourteen-case value, with standard deviation \( 0.274 \). The bootstrap distribution mixes two answers to two
questions, and neither its spread nor its percentiles describe the uncertainty of \( \hat{\beta}_1 \).
:::

::: {when-format="html"}
![**Figure 23.6.1.** (a) Case-bootstrap slopes when one of 15 cases has leverage 0.902, split by whether the resample contains it.
(b) Coverage of nominal 95% intervals for a Pareto mean with infinite variance.](failures.svg){#fig-bs-failures width=100%}
:::

::: {when-format="pdf"}
![(a) Case-bootstrap slopes when one of 15 cases has leverage 0.902, split by whether the resample contains it.
(b) Coverage of nominal 95% intervals for a Pareto mean with infinite variance.](failures.pdf){width=100%}
:::

```{.python .run #cell-failures-leverage}
import numpy as np
rng = np.random.default_rng(2313)
n_l = 15
x_l = np.append(np.linspace(0, 1, n_l - 1), 4.0)    # one case far from the rest
X_l = np.column_stack([np.ones(n_l), x_l])
h_l = np.sum((X_l @ np.linalg.inv(X_l.T @ X_l)) * X_l, axis=1)
y_l = 1 + 0.5 * x_l + 0.3 * rng.normal(size=n_l)
y_l[-1] -= 1.0                                     # and not quite on the line of the others
idx = rng.integers(0, n_l, size=(20_000, n_l))
xs, ys = x_l[idx], y_l[idx]
xc = xs - xs.mean(axis=1, keepdims=True)
b_case = np.sum(xc * ys, axis=1) / np.sum(xc ** 2, axis=1)
absent = ~np.any(idx == n_l - 1, axis=1)
print(f"leverage of the far case {h_l[-1]:.3f}; resamples without it: {absent.mean():.3f}")
```

The residual and wild bootstraps keep the far case, but inherit the problem (@exm-bs-engel-wild): its residual has variance
\( \sigma^2(1-h_{ii}) \), nearly zero, so the data say almost nothing about its error.

## Many parameters

With \( p \) proportional to \( n \), the plain residual bootstrap underestimates variances by the non-vanishing factor \( (n-p)/n \) of
@prp-bs-residual-moments(b). The case bootstrap errs the other way: a resample has only about \( 0.632n \) distinct cases (@exr-bs-distinct), and when \( p/n \) is not small its design is much worse conditioned than the original.

A simulation with \( n=60 \), \( p=30 \), a Gaussian design and normal errors shows both effects. The median ratio of bootstrap
variance to true variance for one coefficient, over \( 200 \) data sets, is \( 0.47 \) for the plain residual bootstrap,
\( 0.97 \) after rescaling the residuals by \( \sqrt{n/(n-p)} \), and \( 3.38 \) for the case bootstrap. El Karoui and
Purdom (2018) study these failures, and Mammen (1993) gives conditions under which the residual and wild bootstraps remain
valid as \( p \) grows. [Chapter 28](../ch28-high-dimensional/index.html) takes up regression with more parameters than cases.

```{.python .run #cell-failures-dimension}
import numpy as np

def variance_ratios(n, p, reps, B, rng):
    """Bootstrap variance of the first coefficient divided by its true variance (errors N(0, 1))."""
    X = rng.normal(size=(n, p))                   # fixed design
    XtX_inv = np.linalg.inv(X.T @ X)
    true_var = XtX_inv[0, 0]
    H = X @ XtX_inv @ X.T
    out = {"residual": [], "rescaled": [], "case": []}
    for _ in range(reps):
        y = rng.normal(size=n)                    # beta = 0
        e = y - H @ y
        fit = H @ y
        for name, r in [("residual", e), ("rescaled", e * np.sqrt(n / (n - p)))]:
            Ys = fit + r[rng.integers(0, n, size=(B, n))]
            out[name].append(np.var(Ys @ (X @ XtX_inv)[:, 0]) / true_var)
        bs = []
        for _ in range(B):
            i = rng.integers(0, n, size=n)
            bs.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0][0])
        out["case"].append(np.var(bs) / true_var)
    return {k: float(np.median(v)) for k, v in out.items()}

rng = np.random.default_rng(2314)
print(variance_ratios(60, 30, reps=20, B=100, rng=rng))
```

## Dependence

Resampling independently destroys serial correlation. With positively autocorrelated errors and a smooth regressor, the true
variance of a slope far exceeds the independent-errors formula (@thm-dep-covariance), which every bootstrap of this chapter
reproduces. Take a straight line in time with
\( n=100 \) equally spaced points and AR(1) errors with parameter \( 0.7 \). Over \( 2000 \) series the slope has
standard deviation \( 1.130 \), while the residual bootstrap's average standard error is \( 0.449 \), a ratio of
\( 0.40 \).

The **moving-block bootstrap** (Künsch 1989) resamples blocks of consecutive residuals and so keeps the dependence within blocks.
With blocks of length \( 10 \) its average standard error is \( 0.818 \), a ratio of \( 0.72 \): better, still short,
because dependence between blocks is lost (Lahiri 2003). When the dependence has a known form it is better to model it, as in
@thm-het-ar1, and resample the estimated innovations.

```{.python .run #cell-failures-dependence}
import numpy as np

def dependence_study(n, phi, reps, B, block, rng):
    """True sd of the slope, and the average iid and moving-block residual-bootstrap standard errors."""
    x = np.linspace(0, 1, n)
    X = np.column_stack([np.ones(n), x])
    A = np.linalg.solve(X.T @ X, X.T)
    slopes, se_iid, se_block = [], [], []
    starts_max = n - block + 1
    for _ in range(reps):
        eps = np.empty(n)
        eps[0] = rng.normal() / np.sqrt(1 - phi ** 2)
        for i in range(1, n):
            eps[i] = phi * eps[i - 1] + rng.normal()   # AR(1) errors
        y = 1 + 2 * x + eps
        b = A @ y
        e = y - X @ b
        slopes.append(b[1])
        Ys = X @ b + e[rng.integers(0, n, size=(B, n))]
        se_iid.append(np.std(Ys @ A[1]))
        starts = rng.integers(0, starts_max, size=(B, n // block))
        idx = (starts[:, :, None] + np.arange(block)).reshape(B, -1)
        Ys = X @ b + e[idx]
        se_block.append(np.std(Ys @ A[1]))
    return np.std(slopes), np.mean(se_iid), np.mean(se_block)

rng = np.random.default_rng(2316)
print(dependence_study(100, 0.7, reps=100, B=200, block=10, rng=rng))
```

## Small samples

A bootstrap distribution computed from eight residuals is an estimate from eight numbers. Take a straight line with \( n=8 \) equally spaced points and normal errors, the most favourable case. In
\( 4000 \) simulated data sets, the classical \( t \) interval covers \( 0.953 \), as it must. The residual
bootstrap percentile interval covers only \( 0.857 \). The reason is @prp-bs-residual-moments(b) combined with
@prp-bs-parametric(b). The percentile interval is roughly \( \hat{\beta}_1\pm1.96\sqrt{(n-2)/n}\,s\norm{\mathbf{a}} \), so its
coverage is about
\[
\Pr\bigl(\lvert t_6\rvert\le1.96\sqrt{6/8}\bigr)=\Pr\bigl(\lvert t_6\rvert\le1.697\bigr)=0.859 .
\]
The studentized residual bootstrap covers \( 0.952 \). In small samples, studentize, and compare with the classical interval.

```{.python .run #cell-failures-small}
import numpy as np
from scipy import stats

def small_n_coverage(n, reps, B, rng):
    """Coverage of 95% intervals for a slope with normal errors: t, residual percentile, bootstrap-t."""
    x = np.linspace(0, 1, n)
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    A = XtX_inv @ X.T
    c = np.sqrt(XtX_inv[1, 1])
    hits = {"t": 0, "percentile": 0, "studentized": 0}
    for _ in range(reps):
        y = rng.normal(size=n)                    # true slope 0
        b = A[1] @ y
        e = y - X @ (A @ y)
        s = np.sqrt(e @ e / (n - 2))
        hits["t"] += abs(b) <= stats.t.ppf(0.975, n - 2) * s * c
        Ys = X @ (A @ y) + e[rng.integers(0, n, size=(B, n))]
        bs = Ys @ A[1]
        ss = np.sqrt(np.sum((Ys - (Ys @ A.T) @ X.T) ** 2, axis=1) / (n - 2))
        lo, hi = np.quantile(bs, [0.025, 0.975])
        hits["percentile"] += lo <= 0 <= hi
        tl, th = np.quantile((bs - b) / (ss * c), [0.025, 0.975])
        hits["studentized"] += b - th * s * c <= 0 <= b - tl * s * c
    return {k: v / reps for k, v in hits.items()}

rng = np.random.default_rng(2318)
print(small_n_coverage(8, reps=300, B=499, rng=rng))
```

::: {.idea}
Before trusting a bootstrap, check each assumption listed at the start of this section, and studentize the root. A permutation
test needs only exchangeability under the null hypothesis, and is exact or not at all.
:::

## Exercises

### A. Check your understanding

::: {#exr-bs-atom-n}
[A1]

Compute \( 1-(1-1/n)^n \) for \( n=10 \) and \( n=1000 \). Does the failure of @prp-bs-maximum become less severe with more data?
:::

::: {.solution}
The values are \( 0.651 \) and \( 0.632 \). The atom decreases slightly toward its limit \( 1-e^{-1} \) but never
vanishes, so the failure does not go away.
:::

::: {#exr-bs-block-mean}
[A2]

In the moving-block bootstrap with blocks of length \( \ell \), explain why the bootstrap variance of a sample mean tends to be
too small when \( \ell \) is much shorter than the range of dependence, and why a very long block makes the estimate noisy.
:::

### B. Practice

::: {#exr-bs-leverage-absent}
[B1]

In @exm-bs-leverage-failure, show that the probability that the far case appears at least twice in a resample is
\( 1-(1-1/n)^n-(1-1/n)^{n-1} \), and evaluate it for \( n=15 \). What does a resample with the far case drawn three times estimate?
:::

::: {.solution}
The number of times the case is drawn is binomial with \( n \) trials and probability \( 1/n \). It is zero with probability \( (1-1/n)^n \)
and one with probability \( n(1/n)(1-1/n)^{n-1}=(1-1/n)^{n-1} \), which gives the formula. For \( n=15 \) this is
\( 1-0.3553-0.3806=0.2641 \). A resample with the case drawn three times fits a line in which that case has triple weight, a
weighted least squares estimate that pulls even harder toward the far case, and not a draw from anything like the sampling
distribution of \( \hat{\beta}_1 \).
:::

### C. Going deeper

::: {#exr-bs-m-out-of-n}
[C1]

For the uniform maximum of @prp-bs-maximum, show that the \( m \)-out-of-\( n \) bootstrap with \( m\to\infty \) and \( m/n\to0 \) gives
\( {\Pr}_*(Y^*_{(m)}=Y_{(n)})=1-(1-1/n)^m\to0 \), and that the bootstrap distribution of \( m(Y_{(n)}-Y^*_{(m)})/Y_{(n)} \) approaches the
standard exponential law.
:::
