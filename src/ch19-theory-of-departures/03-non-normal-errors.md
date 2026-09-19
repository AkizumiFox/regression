# Non-normal errors

Normality entered Part III at one point only: it turned the moments of \( \hbeta \) and \( s^2 \) into exact
distributions (@thm-opt-sampling). Everything that needs only two moments holds for any error distribution
with finite variance. Inference about the mean structure survives in large samples, provided no single observation
dominates the fit. Inference about \( \sigma^2 \) does not survive at any sample size, and neither does prediction
of a single new response.

Throughout, the errors \( \varepsilon_1,\dots,\varepsilon_n \) are independent and identically distributed with mean
\( 0 \) and variance \( \sigma^2\in(0,\infty) \), from a distribution that does not change with \( n \). Where fourth moments
are needed, the *excess kurtosis* is \( \kappa=\mu_4/\sigma^4-3 \), where \( \mu_4=\E(\varepsilon_i^4) \). It is zero for normal errors,
\( 6 \) for centred exponential errors, \( -1.2 \) for uniform errors, and at least \( -2 \) for every distribution, since
\( \mu_4\ge\sigma^4 \) by Jensen's inequality.

## Asymptotic normality of least squares

Large-sample statements concern a sequence of models \( \Y_n=\X_n\bbeta+\be_n \), \( n=1,2,\dots \), with \( p \)
fixed and \( \X_n \) of rank \( p \). We write \( h^*_n=\max_ih_{ii} \) for the largest leverage of \( \X_n \). Besides the
three limit theorems recalled in [Section 14.4](../ch14-correlation-lack-of-fit-prediction/04-correlation-coefficient.html)
(the central limit theorem, Slutsky's lemma and the delta method), we need two more standard results
of probability theory.

::: {.remark}
[Two more limit theorems]

(i) *Lindeberg–Feller.* For each \( n \) let \( Z_{n1},\dots,Z_{nn} \) be independent with mean zero and
\( \sum_i\Var(Z_{ni})=1 \). If \( \sum_i\E\bigl[Z_{ni}^2\,\mathbf 1\{\lvert Z_{ni}\rvert>\eta\}\bigr]\to0 \) for every \( \eta>0 \), then
\( \sum_iZ_{ni}\to\Normal(0,1) \) in distribution (Billingsley 1995, section 27).

(ii) *Cramér–Wold and continuous mapping.* Random vectors \( \mathbf{T}_n \) in \( \Real^q \) converge in distribution to \( \mathbf{T} \) iff \( \mathbf{c}\T\mathbf{T}_n\to\mathbf{c}\T\mathbf{T} \) for every
\( \mathbf{c}\in\Real^q \), and then \( g(\mathbf{T}_n)\to g(\mathbf{T}) \) for every continuous \( g \) (van der Vaart 1998, chapter 2).
:::

::: {#lem-dep-weighted-clt}
[Weighted sums of errors]

Let \( \mathbf{w}_n=(w_{n1},\dots,w_{nn})\T \) be constant vectors with \( \norm{\mathbf{w}_n}=1 \) and
\( m_n=\max_i\lvert w_{ni}\rvert\to0 \). Then \( \sigma^{-1}\sum_iw_{ni}\varepsilon_i\to\Normal(0,1) \) in distribution.
:::

::: {.proof}
Put \( Z_{ni}=w_{ni}\varepsilon_i/\sigma \), so that \( \sum_i\Var(Z_{ni})=\sum_iw_{ni}^2=1 \). With \( U=\varepsilon_1/\sigma \), for every \( \eta>0 \),
\[
\begin{aligned}
\sum_i\E\bigl[Z_{ni}^2\mathbf 1\{\lvert Z_{ni}\rvert>\eta\}\bigr]
&=\sum_iw_{ni}^2\,\E\bigl[U^2\mathbf 1\{\lvert w_{ni}U\rvert>\eta\}\bigr]\\
&\le\E\bigl[U^2\mathbf 1\{\lvert U\rvert>\eta/m_n\}\bigr],
\end{aligned}
\]
because \( \lvert w_{ni}\rvert\le m_n \) and the weights \( w_{ni}^2 \) sum to one. Since \( \E U^2=1 \) and \( \eta/m_n\to\infty \), the
right side tends to zero by dominated convergence. The Lindeberg–Feller theorem applies.
:::

::: {#thm-dep-nonnormal}
[Least squares with non-normal errors]

In the sequence of models above, suppose \( h^*_n\to0 \). Let \( s^2=\norm{(\I-\M)\Y_n}^2/(n-p) \).

::: {.enumerate options="label=(\alph*)"}
1. For every sequence of nonzero vectors \( \mathbf{a}_n\in\Real^p \),
   \[
   \frac{\mathbf{a}_n\T\hbeta-\mathbf{a}_n\T\bbeta}{\sigma\sqrt{\mathbf{a}_n\T(\X_n\T\X_n)^{-1}\mathbf{a}_n}}\ \longrightarrow\ \Normal(0,1)
   \quad\text{in distribution.}
   \]

2. \( s^2\to\sigma^2 \) in probability. This part does not need \( h^*_n\to0 \).

3. The \( t \) statistic, with \( s \) in place of \( \sigma \) in (a), tends to \( \Normal(0,1) \), and the nominal \( 1-\alpha \)
   interval \( \mathbf{a}_n\T\hbeta\pm t_{n-p,\alpha/2}\,s\sqrt{\mathbf{a}_n\T(\X_n\T\X_n)^{-1}\mathbf{a}_n} \) has coverage tending to \( 1-\alpha \).

4. Let \( \bLambda \) be a fixed \( p\times q \) matrix of rank \( q \), and suppose \( \bLambda\T\bbeta=\mathbf{d} \). The statistic
   \[
   F_n=\frac{(\bLambda\T\hbeta-\mathbf{d})\T\bigl[\bLambda\T(\X_n\T\X_n)^{-1}\bLambda\bigr]^{-1}(\bLambda\T\hbeta-\mathbf{d})}{q\,s^2}
   \]
   satisfies \( qF_n\to\chi^2(q) \) in distribution, and \( \Pr\{F_n>F_\alpha(q,n-p)\}\to\alpha \).
:::

:::

::: {.proof}
(a) Put \( \mathbf{G}_n=(\X_n\T\X_n)^{-1} \) and \( \tau_n^2=\mathbf{a}_n\T\mathbf{G}_n\mathbf{a}_n \). Since
\( \hbeta-\bbeta=\mathbf{G}_n\X_n\T\be_n \), the numerator is \( \sum_ic_{ni}\varepsilon_i \) with \( c_{ni}=\mathbf{a}_n\T\mathbf{G}_n\x_{(i)} \), and
\( \sum_ic_{ni}^2=\mathbf{a}_n\T\mathbf{G}_n\X_n\T\X_n\mathbf{G}_n\mathbf{a}_n=\tau_n^2 \). The Cauchy–Schwarz inequality in the inner product
\( \inner{\bu}{\bv}=\bu\T\mathbf{G}_n\bv \) (@prp-mat-cauchy-schwarz applied to \( \mathbf{G}_n^{1/2}\mathbf{a}_n \) and \( \mathbf{G}_n^{1/2}\x_{(i)} \)) gives
\[
c_{ni}^2\le(\mathbf{a}_n\T\mathbf{G}_n\mathbf{a}_n)(\x_{(i)}\T\mathbf{G}_n\x_{(i)})=\tau_n^2h_{ii}.
\]
So \( w_{ni}=c_{ni}/\tau_n \) has unit length and \( \max_i\lvert w_{ni}\rvert\le\sqrt{h^*_n}\to0 \). Apply @lem-dep-weighted-clt.

(b) \( (n-p)s^2=\be_n\T\be_n-\be_n\T\M\be_n \), since \( (\I-\M)\X_n\bbeta=\bzero \). By the weak law of large numbers,
\( \be_n\T\be_n/n\to\sigma^2 \) in probability. The second term is nonnegative with
\( \E(\be_n\T\M\be_n)=\sigma^2\tr\M=\sigma^2p \) (@thm-rv-quadform-mean), so by Markov's inequality
\( \Pr(\be_n\T\M\be_n/n>\delta)\le\sigma^2p/(n\delta)\to0 \). Since \( n/(n-p)\to1 \), \( s^2\to\sigma^2 \). Under a finite fourth moment this is
@cor-lm-s2-consistent; the argument here needs only a finite variance.

(c) The \( t \) statistic is the ratio in (a) multiplied by \( \sigma/s\to1 \), so Slutsky's lemma gives the limit. The
quantile \( t_{n-p,\alpha/2} \) tends to \( z_{\alpha/2} \), because the \( t(\nu) \) distribution tends to \( \Normal(0,1) \) and the normal
distribution function is continuous and strictly increasing. Hence
\( \lvert T_n\rvert-t_{n-p,\alpha/2}\to\lvert Z\rvert-z_{\alpha/2} \), and the coverage tends to
\( \Pr(\lvert Z\rvert\le z_{\alpha/2})=1-\alpha \).

(d) Let \( \W_n=\bLambda\T\mathbf{G}_n\bLambda \), which is positive definite, and
\( \mathbf{Z}_n=\W_n^{-1/2}\bLambda\T(\hbeta-\bbeta)/\sigma \). Under the hypothesis the numerator of \( F_n \) is
\( \sigma^2\norm{\mathbf{Z}_n}^2 \). For a unit vector \( \mathbf{c}\in\Real^q \), \( \mathbf{c}\T\mathbf{Z}_n=\mathbf{a}_n\T(\hbeta-\bbeta)/\sigma \) with
\( \mathbf{a}_n=\bLambda\W_n^{-1/2}\mathbf{c} \), and \( \mathbf{a}_n\T\mathbf{G}_n\mathbf{a}_n=\mathbf{c}\T\W_n^{-1/2}\W_n\W_n^{-1/2}\mathbf{c}=1 \). By (a),
\( \mathbf{c}\T\mathbf{Z}_n\to\Normal(0,1) \), which is the law of \( \mathbf{c}\T\mathbf{Z} \) for \( \mathbf{Z}\sim\Normal_q(\bzero,\I) \). By Cramér–Wold,
\( \mathbf{Z}_n\to\mathbf{Z} \), and by continuous mapping \( \norm{\mathbf{Z}_n}^2\to\chi^2(q) \). Then
\( qF_n=\norm{\mathbf{Z}_n}^2\sigma^2/s^2\to\chi^2(q) \) by (b) and Slutsky's lemma. Finally, \( qF(q,\nu) \) is
\( \chi^2(q) \) divided by an independent \( \chi^2(\nu)/\nu \), which tends to one, so \( qF_\alpha(q,n-p) \) tends to the upper
\( \alpha \) point of \( \chi^2(q) \), and the rejection probability tends to \( \alpha \) as in (c).
:::

Part (a) allows \( \mathbf{a}_n \) to depend on the design. The condition \( h^*_n\to0 \) says that no observation keeps a non-negligible
share of the fit. For a straight line at \( x=1,\dots,n \) the largest leverage is below \( 4/n \) (@exr-dep-spread-leverage). The condition is not a technicality, as the next example shows.

::: {#exm-dep-indicator}
[An observation with its own parameter]

Fit an intercept together with an indicator of the first observation, \( \X_n=[\bone,\mathbf{e}_1] \), so that
\( h_{11}=1 \) for every \( n \). The least squares estimates are \( \hat{\mu}=\bar{Y}_{(1)} \), the mean of observations
\( 2,\dots,n \), and \( \hat{\delta}=Y_1-\bar{Y}_{(1)} \). Hence
\[
\hat{\delta}-\delta=\varepsilon_1-\bar{\varepsilon}_{(1)},\qquad \Var(\hat{\delta})=\sigma^2\Bigl(1+\frac1{n-1}\Bigr),
\]
and since \( \bar{\varepsilon}_{(1)}\to0 \) in probability, \( \hat{\delta}-\delta\to\varepsilon_1 \) in distribution. The limit is the error
distribution itself. No amount of other data can average away the error of one observation.
:::

## How far from normal in a finite sample

@thm-dep-nonnormal is a limit statement. A finite-sample measure of the approach comes from the third and
fourth cumulants of the numerator, which leverage controls directly.

::: {#prp-dep-cumulants}
[Skewness and kurtosis of a linear estimate]

Let the errors have finite fourth moment, skewness \( \gamma_1=\mu_3/\sigma^3 \) with \( \mu_3=\E(\varepsilon_i^3) \), and excess
kurtosis \( \kappa \). For \( \mathbf{a}\ne\bzero \) let \( Z=(\mathbf{a}\T\hbeta-\mathbf{a}\T\bbeta)/\sqrt{\Var(\mathbf{a}\T\hbeta)} \), and let
\( w_i=c_i/\norm{\mathbf{c}} \) with \( \mathbf{c}=\X(\X\T\X)^{-1}\mathbf{a} \). Then
\[
\E(Z^3)=\gamma_1\sum_iw_i^3,\qquad \E(Z^4)-3=\kappa\sum_iw_i^4,
\]
and \( \bigl\lvert\sum_iw_i^3\bigr\rvert\le\sqrt{h^*} \), \( \sum_iw_i^4\le h^* \), where \( h^*=\max_ih_{ii} \).
:::

::: {.proof}
\( Z=\sigma^{-1}\sum_iw_i\varepsilon_i \) with \( \sum_iw_i^2=1 \). Expanding the cube, every term with some index appearing
exactly once has mean zero, by independence and \( \E\varepsilon_i=0 \). Only the terms \( w_i^3\varepsilon_i^3 \) remain, so
\( \E(Z^3)=\sigma^{-3}\mu_3\sum_iw_i^3 \). In the fourth power the surviving terms are \( w_i^4\varepsilon_i^4 \) and, for
\( i\ne j \), the \( 3 \) pairings of \( w_i^2\varepsilon_i^2w_j^2\varepsilon_j^2 \). So
\[
\begin{aligned}
\E(Z^4)&=\frac{\mu_4}{\sigma^4}\sum_iw_i^4+3\sum_{i\ne j}w_i^2w_j^2\\
&=(\kappa+3)\sum_iw_i^4+3\Bigl(1-\sum_iw_i^4\Bigr),
\end{aligned}
\]
using \( \sum_{i\ne j}w_i^2w_j^2=(\sum_iw_i^2)^2-\sum_iw_i^4 \). As in the proof of @thm-dep-nonnormal(a),
\( w_i^2\le h_{ii}\le h^* \). Hence \( \lvert\sum_iw_i^3\rvert\le\max_i\lvert w_i\rvert\sum_iw_i^2\le\sqrt{h^*} \) and
\( \sum_iw_i^4\le\max_iw_i^2\le h^* \).
:::

So skewness is averaged away more slowly than heavy tails, although with symmetric weights it can cancel
exactly (@exr-dep-contrast-skewness). The denominator \( s \) of the \( t \) statistic adds terms of the same order,
so the proposition describes the leading effect only.

::: {#exm-dep-coverage}
[Coverage of the t interval]

Fit a straight line with normal, centred exponential and scaled \( t(3) \) errors, each with variance one, in two
designs. In the *spread* design \( x=1,\dots,n \), and \( h^*_n \) is \( 0.345 \), \( 0.096 \) and
\( 0.025 \) at \( n=10 \), \( 40 \) and \( 160 \). In the *far-point* design the last value is moved to
\( x=4n \), and \( h^*_n \) is \( 0.954 \), \( 0.800 \) and \( 0.485 \). With
\( 100000 \) data sets per case, the coverages of the nominal \( 95\% \) interval for the slope are:

| Design | Errors | \( n=10 \) | \( n=40 \) | \( n=160 \) |
|---|---|---|---|---|
| spread | normal | 0.950 | 0.951 | 0.950 |
| spread | exponential | 0.955 | 0.953 | 0.950 |
| spread | \( t(3) \) | 0.955 | 0.952 | 0.951 |
| far point | normal | 0.951 | 0.950 | 0.951 |
| far point | exponential | 0.933 | 0.948 | 0.950 |
| far point | \( t(3) \) | 0.931 | 0.944 | 0.953 |

The two-sided coverages hide unbalanced tails.
With exponential errors and the far point, the interval misses above the slope in
\( 0.067 \) of the data sets at \( n=10 \) and below it in only \( 0.001 \). At \( n=160 \)
the two proportions are \( 0.038 \) and \( 0.012 \), still far from
\( 0.025 \) each, as the slowly falling leverage predicts. A one-sided test in this design has
nearly three times its nominal size at \( n=10 \). Panel (a) of [Figure 19.3.1](#fig-dep-nonnormal) shows the skewed
\( t \) statistic at \( n=40 \). The \( t(3) \) errors have no fourth moment, but @thm-dep-nonnormal still applies.
:::

```{.python .run #cell-nonnormal-laws}
import numpy as np
from scipy import stats

def errors(law, size, rng):
    """Errors with mean 0 and variance 1."""
    if law == "normal":
        return rng.normal(size=size)
    if law == "exponential":                       # skewed: skewness 2, excess kurtosis 6
        return rng.exponential(size=size) - 1.0
    if law == "t3":                                # heavy tails: variance 3 before scaling
        return rng.standard_t(3, size=size) / np.sqrt(3.0)
    if law == "uniform":                           # short tails: excess kurtosis -1.2
        return rng.uniform(-np.sqrt(3), np.sqrt(3), size=size)

def design(kind, n):
    x = np.arange(1.0, n + 1)
    if kind == "far point":
        x[-1] = 4.0 * n                            # one point far to the right
    return np.column_stack([np.ones(n), x])
```

```{.python .run #cell-nonnormal-slope}
def slope_coverage(kind, law, n, reps, rng):
    """Proportion of nominal 95% t intervals for the slope that cover it."""
    X = design(kind, n)
    G = np.linalg.inv(X.T @ X)
    A = G @ X.T
    E = errors(law, (reps, n), rng)
    B = E @ A.T                                    # beta_hat - beta
    s2 = np.sum((E - B @ X.T) ** 2, axis=1) / (n - 2)
    T = B[:, 1] / np.sqrt(s2 * G[1, 1])
    tq = stats.t.ppf(0.975, n - 2)
    return np.mean(np.abs(T) <= tq), np.mean(T > tq), np.mean(T < -tq)
```

```{.python .run #cell-nonnormal-small}
rng = np.random.default_rng(1903)
for kind in ["spread", "far point"]:
    for law in ["normal", "exponential", "t3"]:
        cov, up, lo = slope_coverage(kind, law, 10, 20_000, rng)
        print(f"{kind:9s} {law:11s} n=10: coverage {cov:.3f} (misses above {up:.3f}, below {lo:.3f})")
```

## Inference about the error variance

The interval for \( \sigma^2 \) (@prp-ci-sigma) rests on \( \text{SSE} \), a sum of *squared* errors, whose variance
depends on the kurtosis (@prp-lm-var-s2, @prp-opt-var-sse). A central limit theorem applies, but with the wrong variance.

::: {#prp-dep-sigma-interval}
[The interval for sigma squared]

In the sequence of models above, let the errors have finite fourth moment with \( \kappa>-2 \). Then
\( \sqrt n\,(s^2/\sigma^2-1)\to\Normal(0,2+\kappa) \), and the coverage of the nominal \( 1-\alpha \) interval of
@prp-ci-sigma tends to
\[
2\,\Phi\Bigl(z_{\alpha/2}\sqrt{2/(2+\kappa)}\Bigr)-1,
\]{#eq-dep-sigma-limit}

where \( \Phi \) is the standard normal distribution function. The limit equals \( 1-\alpha \) iff \( \kappa=0 \).
:::

::: {.proof}
As in the proof of @thm-dep-nonnormal(b),
\[
\begin{aligned}
\sqrt n\,(s^2-\sigma^2)={}&\frac{n}{n-p}\sqrt n\Bigl(\frac{\be_n\T\be_n}{n}-\sigma^2\Bigr)\\
&+\frac{\sqrt n\,p\,\sigma^2}{n-p}-\frac{\sqrt n}{n-p}\,\be_n\T\M\be_n .
\end{aligned}
\]
The \( \varepsilon_i^2 \) are independent with mean \( \sigma^2 \) and variance \( \mu_4-\sigma^4=\sigma^4(2+\kappa) \), so by the central limit
theorem the first term tends to \( \Normal\bigl(0,\sigma^4(2+\kappa)\bigr) \). The second term tends to zero, and so
does the third in probability, because \( \E(\be_n\T\M\be_n)=\sigma^2p \) and Markov's inequality applies to
\( \be_n\T\M\be_n/\sqrt n \). Slutsky's lemma gives the first claim.

Write \( \nu=n-p \) and let \( c_\nu^- \) and \( c_\nu^+ \) be the lower and upper \( \alpha/2 \) points of \( \chi^2(\nu) \). The interval
covers \( \sigma^2 \) iff \( c_\nu^-\le\nu s^2/\sigma^2\le c_\nu^+ \), that is, iff
\[
\frac{c_\nu^--\nu}{\sqrt\nu}\ \le\ \sqrt\nu\Bigl(\frac{s^2}{\sigma^2}-1\Bigr)\ \le\ \frac{c_\nu^+-\nu}{\sqrt\nu}.
\]
A \( \chi^2(\nu) \) variable is a sum of \( \nu \) independent squared standard normals, each with variance \( 2 \), so by the
central limit theorem \( (c_\nu^\pm-\nu)/\sqrt{2\nu}\to\pm z_{\alpha/2} \), and the bounds tend to \( \pm z_{\alpha/2}\sqrt2 \).
The middle term tends to \( \Normal(0,2+\kappa) \) because \( \nu/n\to1 \). The coverage therefore tends to
\( \Pr\bigl(\lvert Z\rvert\sqrt{2+\kappa}\le z_{\alpha/2}\sqrt2\bigr) \), which is @eq-dep-sigma-limit.
:::

For centred exponential errors, \( \kappa=6 \), the nominal \( 95\% \) interval has limiting coverage
\( 2\Phi(1.96/2)-1=0.673 \). For uniform errors, \( \kappa=-1.2 \), it tends to
\( 0.998 \). Panel (b) of [Figure 19.3.1](#fig-dep-nonnormal) shows simulated coverages with the spread design.
For exponential errors they are \( 0.786 \) at \( n=10 \), \( 0.717 \) at
\( n=40 \) and \( 0.672 \) at \( n=2560 \): more data make it *worse*. This answers
@exr-opt-sigma-interval and explains the failure reported in [Section 12.1](../ch12-intervals-and-bands/01-estimable-intervals.html).

::: {when-format="html"}
![**Figure 19.3.1.** (a) The \( t \) statistic for the slope in the far-point design with \( n=40 \) and
centred exponential errors (\( 100000 \) data sets), with the \( t(38) \) density. (b) Coverage of the nominal
\( 95\% \) interval for \( \sigma^2 \) in the spread design, against \( n \), with the limits @eq-dep-sigma-limit dashed.](nonnormal.svg){#fig-dep-nonnormal width=100%}
:::

::: {when-format="pdf"}
![(a) The \( t \) statistic for the slope in the far-point design with \( n=40 \) and
centred exponential errors (\( 100000 \) data sets), with the \( t(38) \) density. (b) Coverage of the nominal
\( 95\% \) interval for \( \sigma^2 \) in the spread design, against \( n \), with the limits @eq-dep-sigma-limit dashed.](nonnormal.pdf){width=100%}
:::

```{.python .run #cell-nonnormal-sigma}
def sigma2_coverage(law, n, reps, rng):
    X = design("spread", n)
    E = errors(law, (reps, n), rng)
    B = np.linalg.solve(X.T @ X, X.T @ E.T).T
    sse = np.sum((E - B @ X.T) ** 2, axis=1)     # sigma^2 = 1
    lo, hi = stats.chi2.ppf([0.025, 0.975], n - 2)
    return np.mean((sse >= lo) & (sse <= hi))

def limit(kappa, z=stats.norm.ppf(0.975)):
    """Limiting coverage of the nominal 95% chi-squared interval."""
    return 2 * stats.norm.cdf(z * np.sqrt(2 / (2 + kappa))) - 1

print("limits:", {k: round(limit(k), 3) for k in [0, 6, -1.2]})
```

Prediction intervals fail for the same reason: they must cover a single new error, which no averaging acts on (@prp-ci-pi-limit).

## Why the F test is protected and the variance-ratio test is not

The \( F \) statistic is a ratio of two quadratic forms, each with a kurtosis-dependent variance. The effects can
cancel, as Box and Watson (1962) and Atiqullah (1962) found. Since \( \log F \) is approximately the difference of
the two mean squares divided by \( \sigma^2 \), it suffices to study that difference.

::: {#prp-dep-balance}
[Kurtosis and a ratio of mean squares]

Let \( \bP_1 \) and \( \bP_2 \) be symmetric idempotent with \( \bP_1\bP_2=\mathbf{0} \), ranks \( f_1 \) and \( f_2 \), and
\( \bP_1\X\bbeta=\bP_2\X\bbeta=\bzero \). Let the errors have finite fourth moment and excess kurtosis \( \kappa \), put
\( S_j^2=\Y\T\bP_j\Y/f_j \), and let \( \mathbf{p}_j \) be the vector of diagonal entries of \( \bP_j \). Then
\( L=(S_1^2-S_2^2)/\sigma^2 \) has mean zero and
\[
\Var(L)=\frac2{f_1}+\frac2{f_2}+\kappa\,\norm{\mathbf{d}}^2,\qquad \mathbf{d}=\frac{\mathbf{p}_1}{f_1}-\frac{\mathbf{p}_2}{f_2} .
\]
In particular, if both \( \bP_1 \) and \( \bP_2 \) have constant diagonals, then \( \mathbf{d}=\bzero \) and \( \Var(L) \) is the same as
under normality.
:::

::: {.proof}
Since \( \bP_j\Y=\bP_j\be \), \( L=\be\T\A\be/\sigma^2 \) with \( \A=\bP_1/f_1-\bP_2/f_2 \), whose diagonal is \( \mathbf{d} \). By
@thm-rv-quadform-mean, \( \E(L)=\tr\A=1-1=0 \). By @thm-rv-quadform-variance with zero means,
\( \Var(\be\T\A\be)=(\mu_4-3\sigma^4)\mathbf{d}\T\mathbf{d}+2\sigma^4\tr(\A^2) \), and
\( \tr(\A^2)=\tr\bP_1/f_1^2+\tr\bP_2/f_2^2-2\tr(\bP_1\bP_2)/(f_1f_2)=1/f_1+1/f_2 \). Divide by \( \sigma^4 \). If the diagonals
are constant, \( \mathbf{p}_j=(f_j/n)\bone \) because \( \tr\bP_j=f_j \), and \( \mathbf{d}=\bzero \).
:::

For the \( F \) test of a linear hypothesis, \( \bP_1=\M-\Mo \) and \( \bP_2=\I-\M \). Designs in which both have
constant diagonals are called *quadratically balanced*. Every balanced one-way or complete factorial layout is
an example. In such designs the kurtosis of the errors does not affect the \( F \) statistic to first order. The
variance-ratio test for comparing the error variances of two independent samples is the opposite case. With
\( n_1 \) and \( n_2 \) observations, each sample centred at its own mean, \( \mathbf{d} \) has entries \( 1/n_1 \) on the first
sample and \( -1/n_2 \) on the second (@exr-dep-variance-ratio), so
\[
\Var(L)=\frac2{f_1}+\frac2{f_2}+\kappa\Bigl(\frac1{n_1}+\frac1{n_2}\Bigr),
\]
and the kurtosis term is as large as the normal-theory term. The test is not robust at any sample size. For
exponential errors and \( n_1=n_2=20 \), the nominal \( 5\% \) test rejected equal variances in
\( 0.268 \) of \( 100000 \) simulated pairs of samples. Its large-sample limit is \( 0.327 \) (@exr-dep-variance-ratio). The \( F \) test for equal means in a balanced one-way layout with three groups of eight had size
\( 0.042 \), and with unbalanced groups of \( 2 \), \( 4 \) and \( 18 \) it had size
\( 0.058 \). In the unbalanced layout \( \Var(L) \) is \( 1.870 \) instead of
\( 1.095 \), which overstates the effect on the size: in the tail, skewness and the curvature of the logarithm
also matter.

::: {.idea}
Inference about \( \bbeta \) is robust to the *shape* of the errors when no observation has large leverage.
Inference about \( \sigma^2 \), variance ratios or single future responses is not, at any sample size.
:::

## Exercises

### A. Check your understanding

::: {#exr-dep-spread-leverage}
[A1]

For a straight line at \( x=1,\dots,n \), show that the largest leverage is \( 1/n+3(n-1)/\bigl(n(n+1)\bigr)<4/n \).
:::

::: {#exr-dep-far-leverage}
[A2]

In the far-point design of @exm-dep-coverage, show that the leverage of the far point is approximately
\( 147/(n+147) \) for large \( n \). Why does the table still show nearly nominal two-sided coverage at \( n=160 \),
when the leverage is almost one half?
:::

### B. Practice

::: {#exr-dep-contrast-skewness}
[B1]

In a one-way layout with \( m \) observations per group and skewed errors, compute \( \sum_iw_i^3 \) of @prp-dep-cumulants for
(i) a single group mean and (ii) the difference of two group means. What does this say about one-sample and
two-sample \( t \) tests with skewed data?
:::

::: {.solution}
(i) The mean of group \( 1 \) has \( c_i=1/m \) on its \( m \) observations, so \( w_i=m^{-1/2} \) there and
\( \sum_iw_i^3=m\cdot m^{-3/2}=m^{-1/2} \). (ii) The difference has \( c_i=\pm1/m \) on the two groups, so \( w_i=\pm(2m)^{-1/2} \), and the cubes
cancel: \( \sum_iw_i^3=0 \). The numerator of a one-sample \( t \) statistic has skewness \( \gamma_1/\sqrt m \). That of a
two-sample comparison with equal group sizes has none, and its non-normality starts with the kurtosis term, of order \( 1/m \).
Equal-sized two-sample comparisons are therefore much better protected against skewness.
:::

::: {#exr-dep-variance-ratio}
[B2]

Two independent samples of size \( m \) are each centred at their own mean. Let \( S_1^2 \) and \( S_2^2 \) be the sample
variances. (i) Show that the vector \( \mathbf{d} \) of @prp-dep-balance has entries \( 1/m \) and \( -1/m \). (ii) Using the first part
of @prp-dep-sigma-interval for each sample, show that when the two variances are equal the two-sided
variance-ratio test of nominal size \( \alpha \) has size tending to \( 2\bigl(1-\Phi(z_{\alpha/2}\sqrt{2/(2+\kappa)})\bigr) \) as \( m\to\infty \).
:::

::: {.solution}
(i) \( \bP_1 \) is the centring matrix on the first sample and zero elsewhere, so its diagonal is \( 1-1/m=f_1/m \) there, and
\( \mathbf{p}_1/f_1 \) has entries \( 1/m \) on the first sample. Likewise \( \mathbf{p}_2/f_2 \) has entries \( 1/m \) on the second, which enter
\( \mathbf{d} \) with a minus sign. (ii) By @prp-dep-sigma-interval (with \( p=1 \)), \( \sqrt m(S_j^2/\sigma^2-1)\to\Normal(0,2+\kappa) \) for
\( j=1,2 \), independently. By the delta method, \( \sqrt m\log(S_1^2/S_2^2)\to\Normal\bigl(0,2(2+\kappa)\bigr) \). Under normality the same
statistic tends to \( \Normal(0,4) \), and the critical values of the \( F(m-1,m-1) \) distribution satisfy
\( \sqrt m\log F_{\alpha/2}\to2z_{\alpha/2} \). The rejection probability therefore tends to
\( \Pr\bigl(\lvert Z\rvert\sqrt{2(2+\kappa)}>2z_{\alpha/2}\bigr)=2\bigl(1-\Phi(z_{\alpha/2}\sqrt{2/(2+\kappa)})\bigr) \). For
\( \kappa=6 \) and \( \alpha=0.05 \) this is \( 0.327 \).
:::

### C. Going deeper

::: {#exr-dep-hetero-clt}
[C1]

Let the errors be independent with mean zero, variances \( \sigma_i^2\in[m_0,M_0] \) with \( m_0>0 \), and
\( \sup_i\E\lvert\varepsilon_i\rvert^3<\infty \). Suppose \( h^*_n\to0 \). Show that for every nonzero \( \mathbf{a}_n \),
\( (\mathbf{a}_n\T\hbeta-\mathbf{a}_n\T\bbeta)/\sqrt{\Var(\mathbf{a}_n\T\hbeta)}\to\Normal(0,1) \), where the variance is the sandwich @eq-lm-sandwich. Explain why this does not make the usual \( t \) statistic asymptotically standard normal.
:::

::: {.solution}
With \( c_{ni} \) as in the proof of @thm-dep-nonnormal, \( \Var(\mathbf{a}_n\T\hbeta)=\sum_ic_{ni}^2\sigma_i^2\ge m_0\tau_n^2 \). Put
\( Z_{ni}=c_{ni}\varepsilon_i/\sqrt{\Var(\mathbf{a}_n\T\hbeta)} \). Lyapunov's bound
\( \E\bigl[Z^2\mathbf 1\{\lvert Z\rvert>\eta\}\bigr]\le\E\lvert Z\rvert^3/\eta \) gives the Lindeberg sum at most
\[
\begin{aligned}
&\frac{\sup_i\E\lvert\varepsilon_i\rvert^3\sum_i\lvert c_{ni}\rvert^3}{\eta\,(m_0\tau_n^2)^{3/2}}\\
&\quad\le\frac{\sup_i\E\lvert\varepsilon_i\rvert^3}{\eta\,m_0^{3/2}}\max_i\frac{\lvert c_{ni}\rvert}{\tau_n}\\
&\quad\le\frac{\sup_i\E\lvert\varepsilon_i\rvert^3}{\eta\,m_0^{3/2}}\sqrt{h^*_n}\to0,
\end{aligned}
\]
using \( \sum_i\lvert c_{ni}\rvert^3\le\max_i\lvert c_{ni}\rvert\,\tau_n^2 \). The usual \( t \) statistic divides by
\( s\,\tau_n \), and \( s^2\tau_n^2 \) estimates \( \bar{\sigma}^2\tau_n^2 \) with \( \bar\sigma^2 \) an average variance, not
\( \sum_ic_{ni}^2\sigma_i^2 \). The ratio of the two is the quantity of @exr-dep-hetero-origin, and it need not tend to one.
[Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html) replaces \( s^2\tau_n^2 \) by a consistent estimate of the sandwich.
:::
