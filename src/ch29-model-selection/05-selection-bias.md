# Selection bias and honest inference

Everything in Part III assumed that the model was fixed before the data were seen. After a search, the reported
model, its coefficients, its error estimate and its test statistics are the ones that looked best, and each is biased:
effects look larger, errors smaller, and intervals cover less often than they claim. This section proves the direction
of the bias, measures it, and describes three remedies.

## The winner's curse

::: {#prp-sel-selection-bias}
[What selection does to estimates and intervals]

::: {.enumerate options="label=(\alph*)"}
1. *(Selected maxima are biased upwards.)* Let \( Z_1,\dots,Z_m \) have finite means \( \theta_1,\dots,\theta_m \), and let
   \( \hat{J} \) be an index at which \( Z_{\hat{J}}=\max_jZ_j \). Then
   \[
   \E(Z_{\hat{J}})\ \ge\ \max_j\theta_j\ \ge\ \E(\theta_{\hat{J}}),
   \]
   so \( \E(Z_{\hat{J}}-\theta_{\hat{J}})\ge0 \). If the \( Z_j \) are independent \( \Normal(\theta,1) \) with a common mean and
   \( m\ge2 \), the bias is \( \E\max_jZ_j-\theta>0 \).

2. *(Selected minima of risk estimates are optimistic.)* Let \( \hat{R}_k \) be an unbiased estimate of \( R_k \) for each
   candidate \( k \) in a finite set, and let \( \hat{k} \) minimize \( \hat{R}_k \). Then
   \[
   \E(\hat{R}_{\hat{k}})\ \le\ \min_kR_k\ \le\ \E(R_{\hat{k}}).
   \]

3. *(Selected intervals under-cover.)* Let \( Z\sim\Normal(\theta,1) \), and report the interval \( Z\pm z_{\alpha/2} \) only when
   \( \lvert Z\rvert>c \). If \( c\ge z_{\alpha/2} \), the probability that the reported interval covers \( \theta \), given that it is
   reported, is \( 0 \) at \( \theta=0 \). It tends to \( 1-\alpha \) as \( \lvert\theta\rvert\to\infty \).
:::

:::

::: {.proof}
(a) \( \max_jZ_j\ge Z_i \) for every \( i \), so \( \E\max_jZ_j\ge\theta_i \) for every \( i \), which is the first inequality. The second
holds pointwise: \( \theta_{\hat{J}}\le\max_j\theta_j \). In the normal case, \( \max_jZ_j-\theta=\max_j(Z_j-\theta) \) has the
distribution of the largest of \( m\ge2 \) independent standard normals, which exceeds \( Z_1-\theta \) with positive probability
and is never smaller, so its mean exceeds \( \E(Z_1-\theta)=0 \).
(b) \( \hat{R}_{\hat{k}}\le\hat{R}_k \) for every \( k \); take expectations to get \( \E\hat{R}_{\hat{k}}\le R_k \) for every \( k \). The second
inequality holds pointwise.
(c) At \( \theta=0 \) the interval covers iff \( \lvert Z\rvert\le z_{\alpha/2} \), which is incompatible with \( \lvert Z\rvert>c\ge z_{\alpha/2} \).
As \( \lvert\theta\rvert\to\infty \), \( \Pr(\lvert Z\rvert>c)\to1 \), so the conditional coverage tends to the unconditional
\( \Pr(\lvert Z-\theta\rvert\le z_{\alpha/2})=1-\alpha \).
:::

Part (a) is the **winner's curse**: the largest of ten unbiased estimates of equal means overstates its mean by
\( 1.539 \) standard errors on average, and the coefficient and \( t \) statistic of a regressor that "won" a screening are
biased in the same way. Part (b) applies to every criterion of this chapter, since \( C_p \) and cross-validation are
unbiased for each fixed candidate (@prp-sel-cp(a) and @prp-sel-cv-target); the search degrees of freedom of
[Section 29.1](01-prediction-error.html) measure the same effect.

::: {#exm-sel-freedman}
[Freedman's paradox]

Freedman (1983) described the following procedure, which many analyses follow in some form. Fit \( n=100 \) cases on
\( 50 \) regressors, none of which is related to the response. Keep those with \( p \)-value below \( 0.25 \), refit on the
survivors, and report. In one simulated data set \( 15 \) regressors survive the first pass, \( 7 \) of them are significant at
\( 5\% \) in the refit, and the overall \( F \) test of the refit has \( p \)-value \( 0.0030 \). Over \( 2000 \) data sets, an average of
\( 12.3 \) regressors survive and \( 5.6 \) are significant, and the overall \( F \) test rejects at \( 5\% \) in a fraction
\( 0.94 \) of them. Pure noise produces a model that every standard test endorses.
:::

```{.python .run #cell-selection-bias-freedman}
import numpy as np
from scipy import stats
rng = np.random.default_rng(2913)
n, m = 100, 50


def ols(X, y):
    """Coefficients, t statistics (slopes), overall F p-value; X without intercept."""
    A = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    e = y - A @ b
    df = len(y) - A.shape[1]
    s2 = e @ e / df
    t = b / np.sqrt(s2 * np.diag(np.linalg.inv(A.T @ A)))
    ssr = np.sum((A @ b - y.mean()) ** 2)
    F = (ssr / X.shape[1]) / s2
    return t[1:], df, stats.f.sf(F, X.shape[1], df)


def screen_and_refit():
    X = rng.normal(size=(n, m))
    y = rng.normal(size=n)                                  # no regressor matters
    t, df, _ = ols(X, y)
    keep = np.abs(t) > stats.t.ppf(1 - 0.25 / 2, df)        # first pass: p < 0.25
    t2, df2, p_F = ols(X[:, keep], y)                       # second pass on the survivors
    return keep.sum(), np.sum(np.abs(t2) > stats.t.ppf(0.975, df2)), p_F


kept, signif, p_F = screen_and_refit()
print(f"one data set: {kept} kept, {signif} significant at 5%, overall F p-value {p_F:.4f}")
runs = np.array([screen_and_refit() for _ in range(500)])
print("over 500 data sets: mean kept", runs[:, 0].mean(), " mean significant", runs[:, 1].mean(),
      " Pr(F significant)", np.mean(runs[:, 2] < 0.05))
```

## Coverage after selection

The failure in @prp-sel-selection-bias(c) concerns intervals reported only after selection. A subtler failure affects
intervals for a coefficient that is always reported, when the *model* around it was chosen from the data.

::: {#exm-sel-pretest-coverage}
[Coverage after a pretest]

Take two standardized regressors with correlation \( r \), known \( \sigma \), and the question: what is \( \beta_1 \)? The
analyst first tests \( \beta_2=0 \) at level \( 5\% \). If the test rejects, she reports the usual \( 95\% \) interval for
\( \beta_1 \) from the model with both regressors; otherwise, the interval from the model with \( \x_1 \) alone. Each interval is
correct for its own model. Panel (b) of [Figure 29.5.1](#fig-sel-post-selection) shows the coverage of the reported
interval as a function of \( \delta=\beta_2/\text{sd}(\hat{\beta}_2) \), the standardized size of the coefficient being tested.
When \( \delta=0 \) the short model is correct and coverage is close to nominal. When \( \delta \) is large the test always keeps
\( \x_2 \) and coverage is nominal. In between, the test often drops \( \x_2 \) although the short slope is biased by
\( r\beta_2 \) (@thm-dep-omitted), and the coverage falls to a minimum of \( 0.873 \) (at \( \delta=2.2 \)) for \( r=0.5 \),
\( 0.599 \) (at \( \delta=1.8 \)) for \( r=0.8 \) and \( 0.234 \) (at \( \delta=1.2 \)) for \( r=0.95 \).

Panel (a) is the one-dimensional case of @prp-sel-selection-bias(c) with \( c=1.96 \): the conditional coverage is \( 0 \)
at \( \theta=0 \), \( 0.844 \) at \( \theta=1 \), \( 0.951 \) at \( \theta=2 \) and \( 0.971 \) at \( \theta=3 \), where it slightly
over-covers.
:::

::: {when-format="html"}
![**Figure 29.5.1.** (a) Coverage of \( Z\pm1.96 \), \( Z\sim\Normal(\theta,1) \), given that it is reported only when
\( \lvert Z\rvert>1.96 \). (b) Coverage for \( \beta_1 \) of the usual \( 95\% \) interval from the model chosen by a \( 5\% \) test of
\( \beta_2=0 \), for three correlations between the regressors (\( 200000 \) simulations per point).](post_selection.svg){#fig-sel-post-selection width=100%}
:::

::: {when-format="pdf"}
![(a) Coverage of \( Z\pm1.96 \), \( Z\sim\Normal(\theta,1) \), given that it is reported only when
\( \lvert Z\rvert>1.96 \). (b) Coverage for \( \beta_1 \) of the usual \( 95\% \) interval from the model chosen by a \( 5\% \) test of
\( \beta_2=0 \), for three correlations between the regressors (\( 200000 \) simulations per point).](post_selection.pdf){width=100%}
:::

```{.python .run #cell-selection-bias-coverage}
z = stats.norm.ppf(0.975)


def conditional_coverage(theta, c=z):
    """Pr(|Z - theta| <= z | |Z| > c) for Z ~ N(theta, 1): the reported interval."""
    lo, hi = theta - z, theta + z
    inside_sel = (stats.norm.cdf(max(hi, c) - theta) - stats.norm.cdf(max(lo, c) - theta)) \
        + (stats.norm.cdf(min(hi, -c) - theta) - stats.norm.cdf(min(lo, -c) - theta))
    return inside_sel / (stats.norm.sf(c - theta) + stats.norm.cdf(-c - theta))


def pretest_coverage(delta, r, reps=200_000):
    """Coverage of the usual 95% interval for beta_1 after a 5% pretest of beta_2 = 0.
    Two standardized regressors with correlation r, sigma known; delta = beta_2 / sd(beta_2 hat)."""
    sd2 = 1 / np.sqrt(1 - r**2)                             # sd of either long-model slope
    C = np.array([[1, -r], [-r, 1]]) / (1 - r**2)           # covariance of (b1, b2) in the long model
    b = rng.normal(size=(reps, 2)) @ np.linalg.cholesky(C).T + [0.0, delta * sd2]
    keep = np.abs(b[:, 1] / sd2) > z                        # the pretest keeps x2
    b1_short = b[:, 0] + r * b[:, 1]                        # short slope (sd 1)
    cover_long = np.abs(b[:, 0]) <= z * sd2                 # true beta_1 = 0
    cover_short = np.abs(b1_short) <= z
    return np.mean(np.where(keep, cover_long, cover_short))


for r, d_min in ((0.5, 2.2), (0.8, 1.8), (0.95, 1.2)):              # the minimizing delta of the figure
    print(f"r = {r}: coverage at delta = 0, {d_min}, 5:",
          ", ".join(f"{pretest_coverage(d, r, 50_000):.3f}" for d in (0.0, d_min, 5.0)))
```

The coverage depends on the unknown \( \delta \), and it is worst where the data are least able to tell whether \( \x_2 \)
belongs. No simple correction repairs this: Leeb and Pötscher (2005, 2006) showed that the distribution of a
post-selection estimator cannot be estimated consistently, uniformly in the parameters, so bootstrapping the whole
procedure does not rescue the naive interval either.

## Honest inference, I: splitting the sample

The simplest remedy separates the two uses of the data. Choose the model on one part of the sample, and compute
intervals from the other part, as if the model had been fixed in advance, which, for that part, it was (Cox 1975).

::: {#prp-sel-splitting}
[Sample splitting]

Let \( \Y=(\Y_1\T,\Y_2\T)\T \), where \( \Y_2\sim\Normal_{n_2}(\bmu_2,\sigma^2\I) \) is independent of \( \Y_1 \), and let \( \X_2 \) be the
fixed design of the second part. Let \( \hat{M} \) be any rule that chooses a set of columns from \( \Y_1 \) (and the design), with values in the set
of \( M \) for which \( \X_{2,M} \) has full column rank \( p_M<n_2 \). For each such \( M \), let
\( \bbeta_M=(\X_{2,M}\T\X_{2,M})^{-1}\X_{2,M}\T\bmu_2 \), and let \( I_{M,j} \) be the usual level \( 1-\alpha \) \( t \) interval for the \( j \)th
coefficient from the least squares fit of \( \Y_2 \) on \( \X_{2,M} \). Then, given \( \Y_1 \), \( I_{\hat{M},j} \) covers
\( \beta_{\hat{M},j} \) with probability at least \( 1-\alpha \), with equality if \( \bmu_2\in\C(\X_{2,\hat{M}}) \). The same holds
unconditionally.
:::

::: {.proof}
Given \( \Y_1 \), the model \( M=\hat{M} \) is fixed and \( \Y_2 \) still has its distribution, by independence. For fixed \( M \),
write \( \M_M \) for the projection onto \( \C(\X_{2,M}) \) and \( w=[(\X_{2,M}\T\X_{2,M})^{-1}]_{jj} \). The estimate
\( \hat{\beta}_{M,j} \) is normal with mean \( \beta_{M,j} \) and variance \( \sigma^2w \), so
\( Z=(\hat{\beta}_{M,j}-\beta_{M,j})/(\sigma\sqrt w)\sim\Normal(0,1) \). The residual sum of squares
\( \norm{(\I-\M_M)\Y_2}^2 \) is a function of \( (\I-\M_M)\Y_2 \), which is uncorrelated with \( \M_M\Y_2 \) and hence independent of it (@thm-mvn-independence), and \( W=\norm{(\I-\M_M)\Y_2}^2/\sigma^2\sim\chi^2(\nu,\gamma) \) with \( \nu=n_2-p_M \) and
\( \gamma=\norm{(\I-\M_M)\bmu_2}^2/\sigma^2 \) (@thm-qf-chisq(a)). The interval covers iff
\( \lvert Z\rvert\le t_{\alpha/2}(\nu)\sqrt{W/\nu} \), iff \( W>\nu Z^2/t_{\alpha/2}(\nu)^2 \). Conditioning on \( Z \) and using
@prp-qf-ncchisq-monotone, this probability is nondecreasing in \( \gamma \), and at \( \gamma=0 \) it is \( 1-\alpha \) by the
definition of the \( t \) distribution. Averaging over \( \Y_1 \) gives the unconditional statement.
:::

The target \( \beta_{M,j} \) is the coefficient of \( \x_j \) in the best approximation to the mean by the chosen columns,
which is what the fitted coefficient estimates, right model or wrong. The costs: selection sees only part of the data,
the intervals are wider, and the result depends on the random split.
Wasserman and Roeder (2009) and Rinaldo, Wasserman and G'Sell (2019) develop the idea for high-dimensional
and assumption-lean settings.

## Honest inference, II: covering every model

The second remedy keeps all the data and widens the intervals until they cover the coefficients of every model that
could have been chosen, simultaneously. Then it does not matter which one was.

::: {#prp-sel-simultaneous}
[Simultaneous inference over submodels]

Assume the normal model @eq-opt-normal-model with \( \rank(\X)=P<n \), let \( s^2 \) be the full-model estimate with
\( \nu=n-P \), and let \( \mathcal M \) be a collection of column subsets \( M \) with \( \X_M \) of full column rank. For \( M\in\mathcal M \) and
\( j\in M \) let \( \bv_{M,j}=\X_M(\X_M\T\X_M)^{-1}\mathbf{e}_j \), where \( \mathbf{e}_j \) picks out the coefficient of \( \x_j \), and
\[
\beta_{M,j}=\bv_{M,j}\T\X\bbeta,\qquad \hat{\beta}_{M,j}=\bv_{M,j}\T\Y,
\]
and for \( K>0 \) let
\[
I_{M,j}(K)=\bigl[\hat{\beta}_{M,j}\pm K\,s\norm{\bv_{M,j}}\bigr].
\]
Let \( N \) be the number of pairs \( (M,j) \) and \( q \) the dimension of the span of the \( \bv_{M,j} \).

::: {.enumerate options="label=(\alph*)"}
1. With \( K=\sqrt{qF_\alpha(q,\nu)} \), \( \Pr\bigl(\beta_{M,j}\in I_{M,j}(K)\text{ for all }M\in\mathcal M,\ j\in M\bigr)\ge1-\alpha \).

2. The same holds with \( K=t_{\alpha/(2N)}(\nu) \).

3. If \( K \) has the property in (a), then for every selection rule \( \hat{M} \) with values in \( \mathcal M \), a function of
   \( \Y \) that may also depend on the unknown \( \bbeta \) and \( \sigma \),
   \( \Pr\bigl(\beta_{\hat{M},j}\in I_{\hat{M},j}(K)\text{ for all }j\in\hat{M}\bigr)\ge1-\alpha \).

4. The smallest \( K \) with the property in (c) for every such rule, including rules that depend on the unknown
   parameters (the convention of Berk, Brown, Buja, Zhang and Zhao 2013), is \( K_{\text{PoSI}} \), the \( 1-\alpha \) quantile of
   \( \max_{M,j}\lvert\bv_{M,j}\T(\Y-\X\bbeta)\rvert/(s\norm{\bv_{M,j}}) \), whose distribution does not depend on
   \( \bbeta \) or \( \sigma \).
:::

:::

::: {.proof}
In each case \( \hat{\beta}_{M,j}-\beta_{M,j}=\bv_{M,j}\T(\Y-\X\bbeta) \) with \( \bv_{M,j}\in\C(\X) \).
(a) Let \( \mathcal V \) be the span of the \( \bv_{M,j} \), a subspace of \( \C(\X) \) of dimension \( q \). By @thm-mc-scheffe(a), with
probability exactly \( 1-\alpha \), \( \lvert\bv\T(\Y-\X\bbeta)\rvert\le\sqrt{qF_\alpha(q,\nu)}\,s\norm{\bv} \) for every \( \bv\in\mathcal V \),
in particular for every \( \bv_{M,j} \).
(b) Each \( (\hat{\beta}_{M,j}-\beta_{M,j})/(s\norm{\bv_{M,j}}) \) has the \( t(\nu) \) distribution, because \( \bv_{M,j}\T\Y \) is independent of
\( s^2 \) when \( \bv_{M,j}\in\C(\X) \) (@thm-opt-sampling). Apply @thm-mc-bonferroni(b) with \( \alpha_{M,j}=\alpha/N \).
(c) The event in (c) contains the simultaneous event of (a).
(d) The simultaneous event of (a) is \( \{\max_{M,j}T_{M,j}\le K\} \) with
\( T_{M,j}=\lvert\bv_{M,j}\T(\Y-\X\bbeta)\rvert/(s\norm{\bv_{M,j}}) \). Its distribution is free of \( \bbeta \) and \( \sigma \) because
\( \Y-\X\bbeta=\sigma\bz \) with \( \bz\sim\Normal_n(\bzero,\I) \) and \( s=\sigma\norm{(\I-\M)\bz}/\sqrt\nu \). So \( K_{\text{PoSI}} \) has the
property in (a) and hence in (c). If \( K<K_{\text{PoSI}} \), take the rule that selects a model containing a pair at which
\( T_{M,j} \) is largest, a rule that uses \( \bbeta \) and so is allowed in (c); its coverage is \( \Pr(\max T_{M,j}\le K)<1-\alpha \).
:::

This is the PoSI (post-selection inference) approach of Berk, Brown, Buja, Zhang and Zhao (2013): model selection
treated as a multiple-comparison problem of [Chapter 13](../ch13-multiplicity/index.html). The Scheffé bound never exceeds
\( \sqrt{PF_\alpha(P,\nu)} \), however many models are searched, and the price is paid whether or not the search was
aggressive.

::: {#exm-sel-state-posi}
[Honest intervals for the state regressions]

For the state data, with the intercept always in, the \( 32 \) candidate models contain \( N=80 \) slope coefficients. Every
\( \bv_{M,j} \) is orthogonal to \( \bone \), so they span at most \( q=5 \) dimensions, and \( \nu=44 \). The multipliers are
\( 2.015 \) for the naive \( t \) interval, \( 3.12 \) for PoSI (by simulation), \( 3.484 \) for Scheffé and \( 3.684 \) for
Bonferroni. In the model chosen by AIC, single parenthood and urbanization, the \( t \) statistics computed with the full-model
\( s \) are \( 5.46 \) and \( 1.83 \). Single parenthood clears every multiplier. Urbanization clears none, not even the naive one:
AIC kept it because its squared \( t \) exceeds about \( 2 \) (the threshold of [Section 29.2](02-cp-aic-bic.html)), which is a
decision about prediction, not evidence of an effect.
:::

```{.python .run #cell-selection-bias-posi}
import itertools
import statsmodels.api as sm
rng = np.random.default_rng(2914)
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
ys = np.log(data["violent"].to_numpy())
Xs = np.column_stack([np.ones(len(ys)), data[names].to_numpy()])
ns, P = Xs.shape
nu = ns - P

# unit vectors u with  coefficient of x_j in submodel M  =  (u' y) / ||v||,  for every M and j in M
U = []
for k in range(1, 6):
    for M in itertools.combinations(range(1, 6), k):
        XM = Xs[:, [0, *M]]
        V = XM @ np.linalg.inv(XM.T @ XM)                   # column i is v for coefficient i
        for i in range(1, len(M) + 1):
            U.append(V[:, i] / np.linalg.norm(V[:, i]))
U = np.array(U)                                              # 80 unit vectors, all in C(X) minus 1
Qb, _ = np.linalg.qr(Xs)
Qb = Qb[:, 1:]                                               # orthonormal basis of C(X) orthogonal to 1
W = U @ Qb                                                   # coordinates of the u's in that basis
sims = 100_000
Zp = rng.normal(size=(sims, P - 1))
s = np.sqrt(rng.chisquare(nu, size=sims) / nu)
K_posi = np.quantile(np.max(np.abs(Zp @ W.T), axis=1) / s, 0.95)
K_t = stats.t.ppf(0.975, nu)
K_bonf = stats.t.ppf(1 - 0.025 / len(U), nu)
K_sch = np.sqrt((P - 1) * stats.f.ppf(0.95, P - 1, nu))
print(f"{len(U)} coefficients; multipliers: t {K_t:.3f}, PoSI {K_posi:.3f}, "
      f"Bonferroni {K_bonf:.3f}, Scheffe {K_sch:.3f}")
```

::: {#exm-sel-three-ways}
[Three ways to report after forward selection]

Take \( n=100 \) cases of ten normal regressors with pairwise correlations \( 0.5 \), \( \sigma=1 \), choose a model by forward
selection with AIC, and report \( 95\% \) intervals for the chosen coefficients, with the projection coefficients
\( \beta_{\hat{M},j} \) as targets. Over \( 3000 \) data sets, the fraction of reported intervals that cover their targets is:

| | naive | split in half | Scheffé |
|---|---|---|---|
| no signal: coverage | 0.440 | 0.946 | 1.000 |
| no signal: mean half-width | 0.228 | 0.351 | 0.516 |
| three slopes of 0.3: coverage | 0.868 | 0.956 | 1.000 |
| three slopes of 0.3: mean half-width | 0.248 | 0.379 | 0.559 |

With no signal, fewer than half of the naive intervals cover: a coefficient is reported only because its estimate was
large. Splitting restores the nominal level at the cost of intervals about one and a half times as wide. The Scheffé
intervals, with multiplier \( \sqrt{10F_{0.05}(10,89)}=4.403 \), are valid for every model and every rule, and here they are
very conservative, because forward selection explores only a small part of the \( 5120 \) coefficients they protect.
:::

## Honest inference, III: conditioning on the selection

The third remedy uses all the data and asks for coverage *given* that the model was selected. In the one-dimensional
setting of @prp-sel-selection-bias(c) this can be done exactly.

::: {#exm-sel-truncated}
[A selective interval]

Let \( Z\sim\Normal(\theta,1) \) be reported only when \( \lvert Z\rvert>1.96 \), and suppose \( Z=2.5 \) is observed. Given selection, \( Z \)
has the normal distribution truncated to \( \lvert z\rvert>1.96 \). Let \( F_\theta \) be its distribution function. Then
\( F_\theta(Z) \) is uniform on \( (0,1) \) given selection, whatever \( \theta \), and \( F_\theta(z) \) decreases in \( \theta \) (@exr-sel-trunc-monotone). So \( \{\theta:0.025\le F_\theta(2.5)\le0.975\} \) is an interval with conditional coverage
exactly \( 95\% \) for every \( \theta \). Numerically it is \( [-0.42,4.33] \), against the naive \( [0.54,4.46] \). An estimate
just past the threshold is weak evidence once the threshold is taken into account, and the selective interval
includes zero.
:::

```{.python .run #cell-selection-bias-selective}
from scipy.optimize import brentq


def trunc_cdf(x, theta, c=z):
    """Distribution function at x of N(theta, 1) conditioned on |Z| > c."""
    den = stats.norm.sf(c - theta) + stats.norm.cdf(-c - theta)
    if x <= -c:
        num = stats.norm.cdf(x - theta)
    elif x < c:
        num = stats.norm.cdf(-c - theta)
    else:
        num = stats.norm.cdf(-c - theta) + stats.norm.cdf(x - theta) - stats.norm.cdf(c - theta)
    return num / den


z_obs = 2.5
lower = brentq(lambda th: trunc_cdf(z_obs, th) - 0.975, -10, z_obs)
upper = brentq(lambda th: trunc_cdf(z_obs, th) - 0.025, -10, 10)
print(f"naive interval [{z_obs - z:.2f}, {z_obs + z:.2f}], selective interval [{lower:.2f}, {upper:.2f}]")
```

Lee, Sun, Sun and Taylor (2016) extended the construction to selection events that are polyhedra
\( \{\A\Y\le\bb\} \), which covers the lasso at a fixed penalty, forward stepwise selection (Tibshirani, Taylor, Lockhart and
Tibshirani 2016) and thresholding: a linear function of \( \Y \), given the selection and the part of \( \Y \) orthogonal to it,
is a truncated normal. The intervals can be very long near the selection boundary. Selecting on part of the data but
using all of it for inference ("data carving", Fithian, Sun and Taylor 2014) shortens them. In high dimensions the
debiased lasso of @thm-hd-debiased gives intervals for the full model's coefficients instead.

None of the three routes is free. The cheapest honest inference is still to fix the model before seeing the data, and
to use selection for prediction, where its errors can be measured.

## Exercises

### A. Check your understanding

::: {#exr-sel-cond-cov}
[A1]

In @prp-sel-selection-bias(c), take \( c=1.5<z_{0.025} \). Compute the conditional coverage at \( \theta=0 \).
:::

### B. Practice

::: {#exr-sel-mills}
[B1]

In the orthonormal design of @prp-sel-orthonormal with \( \sigma=1 \) and a zero coefficient, show that
\( \E(\lvert Z_j\rvert\mid\lvert Z_j\rvert>c)=\phi(c)/(1-\Phi(c)) \). Evaluate it for \( c=\sqrt2 \) and \( c=1.96 \), and interpret it as the
average size of a selected null coefficient.
:::

::: {.solution}
By symmetry \( \E(\lvert Z\rvert\mid\lvert Z\rvert>c)=\E(Z\mid Z>c)=\int_c^\infty z\phi(z)\,dz/(1-\Phi(c))=\phi(c)/(1-\Phi(c)) \). For
\( c=\sqrt2 \) this is \( 0.1468/0.0786=1.87 \); for \( c=1.96 \) it is \( 0.0584/0.0250=2.34 \). A coefficient that is truly zero, once
selected, is reported with an estimate of about two standard errors on average: the winner's curse of
@prp-sel-selection-bias(a), in its purest form.
:::

::: {#exr-sel-posi-count}
[B2]

Verify that the \( 32 \) candidate models of @exm-sel-state-posi contain \( 80 \) slope coefficients, and that every
\( \bv_{M,j} \) for a slope is orthogonal to \( \bone \). Why does this let Scheffé's bound use \( q=5 \) instead of \( 6 \)?
:::

::: {.solution}
Each of the \( 5 \) regressors appears in \( 2^4=16 \) of the models, so there are \( 5\times16=80 \) pairs. The vector
\( \bv_{M,j}=\X_M(\X_M\T\X_M)^{-1}\mathbf{e}_j \) satisfies \( \X_M\T\bv_{M,j}=\mathbf{e}_j \), and the intercept column of \( \X_M \) is
\( \bone \), with \( \mathbf{e}_j \) zero in that position, so \( \bone\T\bv_{M,j}=0 \). All the \( \bv_{M,j} \) lie in \( \C(\X)\cap\bone\perpc \),
of dimension \( 5 \), and @thm-mc-scheffe(a) needs only a subspace containing them.
:::

### C. Going deeper

::: {#exr-sel-trunc-monotone}
[C1]

Let \( F_\theta \) be the distribution function of \( \Normal(\theta,1) \) truncated to a set \( A \) of positive probability. Show that for
\( \theta_1<\theta_2 \) the density ratio \( f_{\theta_2}(z)/f_{\theta_1}(z) \) is increasing in \( z \in A \), and deduce that
\( F_{\theta_2}(z)\le F_{\theta_1}(z) \) for every \( z \). Why does this make the set in @exm-sel-truncated an interval?
:::

::: {.solution}
On \( A \), \( f_\theta(z)=\phi(z-\theta)/\Pr_\theta(A)=\phi(z)e^{\theta z-\theta^2/2}/\Pr_\theta(A) \), so the ratio is a positive constant
times \( e^{(\theta_2-\theta_1)z} \), increasing in \( z \). Put \( g=f_{\theta_2}/f_{\theta_1} \) on \( A \). If \( g(z)\le1 \) then \( g\le1 \) on
\( A\cap(-\infty,z] \), so \( F_{\theta_2}(z)\le F_{\theta_1}(z) \). If \( g(z)>1 \) then \( g>1 \) on \( A\cap[z,\infty) \), so
\( 1-F_{\theta_2}(z)\ge1-F_{\theta_1}(z) \). Either way \( F_{\theta_2}(z)\le F_{\theta_1}(z) \). For fixed \( z \), \( \theta\mapsto F_\theta(z) \)
is therefore nonincreasing (and continuous), so \( \{\theta:0.025\le F_\theta(z)\le0.975\} \) is an interval.
:::

::: {#exr-sel-pretest-formula}
[C2]

In @exm-sel-pretest-coverage, let \( \hat{\beta}_1^{\mathrm{S}} \) be the short slope and \( \hat{\beta}_2 \) the long-model estimate of \( \beta_2 \). Use
@lem-dep-short-long to show that the two are independent, and write the coverage of the reported interval as
\[
\Pr(\lvert T_2\rvert\le z)\,\Pr\bigl(\lvert\hat{\beta}_1^{\mathrm{S}}-\beta_1\rvert\le z\,\text{sd}(\hat{\beta}_1^{\mathrm{S}})\bigr)
+\Pr\bigl(\lvert T_2\rvert>z,\ \lvert\hat{\beta}_1-\beta_1\rvert\le z\,\text{sd}(\hat{\beta}_1)\bigr),
\]
where \( T_2=\hat{\beta}_2/\text{sd}(\hat{\beta}_2) \). Explain why only the first term depends on the bias \( r\beta_2 \) of the short slope in
this simple way.
:::

