# Freedman–Lane and related schemes

With a continuous nuisance regressor in @eq-bs-partial-model, no permutation of the responses leaves the null distribution
unchanged. Several schemes permute something else, or permute the responses and hope. None is exact in general, but their
differences can be understood with the projection algebra of [Chapter 6](../ch06-projections/index.html).

## Five schemes

Keep the notation of @eq-bs-partial-model: \( \Z \) is \( n\times p_0 \) with \( \bone\in\C(\Z) \), \( \X_1 \) is \( n\times q \), the full
model matrix \( [\Z,\X_1] \) has rank \( p=p_0+q<n \), \( \Mo \) projects onto \( \C(\Z) \) and \( \M \) onto \( \C([\Z,\X_1]) \). The statistic
for \( H_0:\bbeta_1=\bzero \) is the \( F \) statistic of @thm-glh-f-test, written as a function of the response vector:
\[
F(\bu)=\frac{\norm{(\M-\Mo)\bu}^2/q}{\norm{(\I-\M)\bu}^2/(n-p)} .
\]{#eq-bs-F-function}

Since \( (\M-\Mo)\Z=\bzero \) and \( (\I-\M)\Z=\bzero \), \( F(\bu+\Z\mathbf{g})=F(\bu) \) for every \( \mathbf{g} \). Let
\[
\he_0=(\I-\Mo)\y,\qquad \he=(\I-\M)\y
\]
be the residuals of the reduced and the full model. By @thm-proj-fwl, \( \M-\Mo \) projects onto \( \C(\tilde{\X}_1) \) with
\( \tilde{\X}_1=(\I-\Mo)\X_1 \).

::: {#def-bs-schemes}
[Permutation schemes for a partial hypothesis]

For a permutation matrix \( \bP \), each scheme computes a permuted statistic; its p-value is @eq-bs-mc-p with \( T(\y) \) the
observed \( F(\y) \) and \( T(\bP_b\y) \) replaced by the scheme's permuted statistic.

::: {.enumerate options="label=(\alph*)"}
1. **Raw permutation** (Manly 1997): \( F(\bP\y) \), refitting the full model to permuted responses.

2. **Permuting the regressors of interest** (Draper and Stoneman 1966): the \( F \) statistic computed from \( \y \) with model matrix
   \( [\Z,\bP\X_1] \).

3. **Kennedy's scheme** (Kennedy 1995): regress \( \bP\he_0 \) on \( \tilde{\X}_1 \) alone and compute
   \( F_K(\bP\he_0) \), where
   \( F_K(\bu)=\dfrac{\norm{(\M-\Mo)\bu}^2/q}{\bigl(\norm{\bu}^2-\norm{(\M-\Mo)\bu}^2\bigr)/(n-p)} \).

4. **Freedman–Lane** (Freedman and Lane 1983): form \( \y^*=\Mo\y+\bP\he_0 \), the reduced-model fit plus permuted
   reduced-model residuals, and compute \( F(\y^*) \).

5. **ter Braak's scheme** (ter Braak 1992): form \( \y^*=\M\y+\bP\he \), the full-model fit plus permuted full-model residuals,
   and compute the \( F \) statistic for the hypothesis \( \bbeta_1=\hbeta_1 \), which is true in the permuted world.
:::

:::

## What can be proved

::: {#prp-bs-freedman-lane}
[The Freedman–Lane scheme]

In the setting above:

::: {.enumerate options="label=(\alph*)"}
1. \( F(\Mo\y+\bP\he_0)=F(\bP\he_0) \) for every \( \bP \), and \( F(\y)=F(\he_0) \). ter Braak's permuted statistic equals \( F(\bP\he) \).

2. Under \( H_0 \), \( \he_0=(\I-\Mo)\be \). The p-values of the Freedman–Lane, ter Braak, Kennedy and Draper–Stoneman schemes depend
   on \( \y \) only through \( (\I-\Mo)\y \), so under \( H_0 \) their distributions do not depend on \( \bgamma \). The raw permutation
   statistic involves \( \bgamma \).

3. Under \( H_0 \), if \( \bP\Z=\Z \) for every \( \bP \) in a group \( \mathcal{G} \) and \( \be \) is \( \mathcal{G} \)-invariant in distribution, then the
   Freedman–Lane test over \( \mathcal{G} \) is exact, and it coincides with the raw permutation test over \( \mathcal{G} \).

4. Suppose the errors are independent draws from one law with mean zero and variance \( \sigma^2>0 \), \( \bP \) is uniform
   over all permutations and independent of \( \be \), and \( p \) is fixed as \( n\to\infty \). Then under \( H_0 \),
   \( F(\bP\he_0)-F(\bP\be)\to0 \) in probability, and so does \( F(\bP\he)-F(\bP\be) \).
:::

:::

::: {.proof}
(a) \( \Mo\y\in\C(\Z) \), so the invariance of @eq-bs-F-function gives the first claim, and \( \y=\Mo\y+\he_0 \) gives the second. In
ter Braak's scheme, the \( F \) statistic for \( \bbeta_1=\hbeta_1 \) computed from \( \y^* \) is the \( F \) statistic for \( \bbeta_1=\bzero \)
computed from \( \y^*-\X_1\hbeta_1 \), whose reduced model is the set of means \( \X_1\hbeta_1+\C(\Z) \): the numerators agree by @thm-glh-general-f(b), and
the denominators agree because subtracting \( \X_1\hbeta_1\in\C(\X) \) leaves the residuals unchanged. Since
\( \y^*-\X_1\hbeta_1=\Z\hat{\bgamma}+\bP\he \) with \( \Z\hat{\bgamma}\in\C(\Z) \), the statistic is \( F(\bP\he) \).

(b) Under \( H_0 \), \( (\I-\Mo)\y=(\I-\Mo)\be \). Freedman–Lane and Kennedy use \( \y \) only through \( \he_0 \), and ter Braak through
\( \he=(\I-\M)\he_0 \). The Draper–Stoneman statistic is unchanged when a vector in \( \C(\Z) \) is added to \( \y \). In each case the observed statistic is \( F(\y)=F(\he_0) \) by (a). The raw permutation
statistic \( F(\bP\y)=F(\bP\Z\bgamma+\bP\be) \) involves \( \bgamma \) whenever \( \bP\Z\bgamma\notin\C(\Z) \), and its null
distribution in general depends on \( \bgamma \), as [Figure 23.5.1](#fig-bs-freedman-lane-size) shows; we do not prove this.

(c) \( \bP\Z=\Z \) implies \( \bP\Mo\bP\T=\Mo \), since \( \bP\Mo\bP\T \) is the orthogonal projection onto \( \C(\bP\Z)=\C(\Z) \). So
\( \bP\he_0=\bP(\I-\Mo)\be=(\I-\Mo)\bP\be \) under \( H_0 \), and \( F(\bP\he_0)=F(\bP\be) \) by invariance. Also \( F(\y)=F(\be) \). So the
Freedman–Lane p-value is the permutation p-value of the statistic \( F \) applied to \( \be \), which is exact by @thm-bs-permutation.
The raw permutation statistic is \( F(\Z\bgamma+\bP\be)=F(\bP\be) \) as well.

(d) Put \( \mathbf{v}=\Mo\be \), so \( \bP\he_0=\bP\be-\bP\mathbf{v} \). For a fixed vector \( \mathbf{v} \) and a uniform random permutation,
\( \E(\bP\mathbf{v})=\bar{v}\bone \) and \( \Cov(\bP\mathbf{v})=\frac{S_v}{n-1}\bigl(\I-n^{-1}\bone\bone\T\bigr) \), where
\( S_v=\sum_i(v_i-\bar{v})^2\le\norm{\mathbf{v}}^2 \) (@exr-bs-permutation-moments). Since \( (\M-\Mo)\bone=\bzero \),
\[
\begin{aligned}
\E\bigl[\norm{(\M-\Mo)\bP\mathbf{v}}^2\mid\be\bigr]&=\frac{S_v}{n-1}\,\tr(\M-\Mo)\\
&\le\frac{q\,\norm{\Mo\be}^2}{n-1},
\end{aligned}
\]
and \( \E\norm{\Mo\be}^2=\sigma^2p_0 \), so \( \norm{(\M-\Mo)\bP\mathbf{v}}\to0 \) in probability. The numerators of the two statistics
are \( \norm{\mathbf{N}}^2/q \) and \( \norm{\mathbf{N}-\boldsymbol{\updelta}}^2/q \) with \( \mathbf{N}=(\M-\Mo)\bP\be \) and
\( \boldsymbol{\updelta}=(\M-\Mo)\bP\mathbf{v} \). Here \( \bP\be \) has the distribution of \( \be \), so
\( \E\norm{\mathbf{N}}^2=q\sigma^2 \), and
\( \bigl\lvert\norm{\mathbf{N}-\boldsymbol{\updelta}}^2-\norm{\mathbf{N}}^2\bigr\rvert\le2\norm{\mathbf{N}}\norm{\boldsymbol{\updelta}}+\norm{\boldsymbol{\updelta}}^2\to0 \)
in probability. For the denominators,
\( \bigl\lvert\norm{(\I-\M)\bP\he_0}-\norm{(\I-\M)\bP\be}\bigr\rvert\le\norm{\mathbf{v}}=O_p(1) \), while
\( \norm{(\I-\M)\bP\be}^2/(n-p) \) has the distribution of \( s^2 \) and tends to \( \sigma^2 \) in probability (@cor-lm-s2-consistent needs
fourth moments; the argument of @lem-bs-residual-tails needs only the second). Both norms are \( O_p(\sqrt n) \), so the difference of
their squares is at most \( \norm{\mathbf{v}} \) times \( O_p(\sqrt n) \). Dividing by \( n-p \), the two denominators differ by a quantity that
tends to zero, and both tend to \( \sigma^2 \). The claim for \( F(\bP\he) \) is the same argument with
\( \mathbf{v}=\M\be \), \( \E\norm{\M\be}^2=\sigma^2p \).
:::

Part (d) compares Freedman–Lane with an **oracle** that permutes the true errors. The oracle is exact: \( F(\y)=F(\be) \) under \( H_0 \),
and it compares \( F(\be) \) with \( F(\bP_b\be) \), a permutation test for the exchangeable \( \be \). Freedman–Lane replaces \( \bP_b\be \) by
\( \bP_b\he_0 \), which changes the statistics by vanishing amounts. To conclude that its size tends to \( \alpha \), one needs also that the
oracle statistics have a continuous joint limit, so that vanishing perturbations rarely change their ranks. That holds under
conditions like @eq-bs-no-dominant by a permutation central limit theorem that we do not prove; Anderson and Robinson (2001) give the
argument.

Part (b) is the practical point: four of the schemes see the data only after \( \Z \) has been projected out, while for large
\( \Z\bgamma \) the raw permutation mostly permutes \( \Z\bgamma \), not the errors. Kennedy's scheme can be pinned down exactly.

::: {#prp-bs-kennedy}
[Kennedy's scheme is anticonservative relative to Freedman–Lane]

For every permutation \( \bP \), \( F_K(\bP\he_0)\le F(\bP\he_0) \), with equality when \( \bP=\I \). Consequently, for any set of permutations
\( \bP_1,\dots,\bP_B \), the Kennedy p-value is at most the Freedman–Lane p-value.
:::

::: {.proof}
The two statistics have the same numerator. Write \( \bu=\bP\he_0 \). Since \( \I=\Mo+(\M-\Mo)+(\I-\M) \) with orthogonal ranges,
\[
\begin{aligned}
\norm{(\I-\M)\bu}^2&=\norm{\bu}^2-\norm{(\M-\Mo)\bu}^2-\norm{\Mo\bu}^2\\
&\le\norm{\bu}^2-\norm{(\M-\Mo)\bu}^2,
\end{aligned}
\]
so Kennedy's denominator is at least as large. For \( \bP=\I \), \( \Mo\he_0=\bzero \), and the two agree. Each permuted Kennedy
statistic is therefore at most the corresponding Freedman–Lane one, while the observed statistics coincide, so fewer permuted values
reach the observed one.
:::

The discrepancy \( \norm{\Mo\bP\he_0}^2 \), the part of the permuted residuals that \( \Z \) could have explained, has expectation about
\( (p_0-1)\sigma^2 \), small beside \( (n-p)\sigma^2 \) for large \( n \); so Kennedy's scheme is liberal only in small samples.

## Stack loss again, and a simulation

::: {#exm-bs-stackloss-freedman-lane}
[Acid concentration given air flow and water temperature]

In the full stack loss regression, test whether acid concentration matters once air flow and water temperature are in the
model. The \( F \) statistic is \( 0.947 \) on \( 1 \) and \( 17 \) degrees of freedom, with p-value
\( 0.344 \). With \( 19999 \) permutations each, the schemes give: raw permutation \( 0.346 \),
Draper–Stoneman \( 0.343 \), Kennedy \( 0.315 \), Freedman–Lane \( 0.345 \) and ter Braak
\( 0.345 \). Only Kennedy's differs noticeably, in the direction @prp-bs-kennedy predicts. The marginal association of
@exm-bs-stackloss-permutation disappears once the correlated regressors are in the model.
:::

```{.python .run #cell-freedman-lane-stackloss}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.stackloss.load_pandas().data

def proj(Z):
    """Orthogonal projection onto C(Z), from a QR factorization."""
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T

def perm_rows(v, B, rng):
    """B random permutations of the vector v, one per row."""
    return v[np.argsort(rng.random((B, len(v))), axis=1)]

def perm_pvalues(y, x, Z, B, rng):
    """p-values for H0: beta = 0 in E(y) = Z gamma + x beta, by six methods."""
    n, p = len(y), Z.shape[1] + 1
    MZ = proj(Z)
    xt = x - MZ @ x                               # x residualized on Z
    yt = y - MZ @ y                               # reduced-model residuals
    e_full = yt - xt * (xt @ yt) / (xt @ xt)      # full-model residuals (FWL)

    def F(Y):
        """F statistic for beta = 0, for each row of Y (model matrix [Z, x] fixed)."""
        num = (Y @ xt) ** 2 / (xt @ xt)
        rss_Z = np.sum(Y ** 2, axis=-1) - np.sum((Y @ MZ) * Y, axis=-1)
        return num / ((rss_Z - num) / (n - p))

    F_obs = F(y)
    out = {"F test": stats.f.sf(F_obs, 1, n - p)}

    def pval(F_star, F0=F_obs):
        return (1 + np.sum(F_star >= F0)) / (len(F_star) + 1)

    out["Manly"] = pval(F(perm_rows(y, B, rng)))                    # permute y
    Xs = perm_rows(x, B, rng)                                       # permute x (Draper-Stoneman)
    xts = Xs - Xs @ MZ
    num = (xts @ y) ** 2 / np.sum(xts ** 2, axis=1)
    out["Draper-Stoneman"] = pval(num / ((yt @ yt - num) / (n - p)))
    Yk = perm_rows(yt, B, rng)                                      # Kennedy: no Z in the refit
    num = (Yk @ xt) ** 2 / (xt @ xt)
    out["Kennedy"] = pval(num / ((np.sum(Yk ** 2, axis=1) - num) / (n - p)))
    out["Freedman-Lane"] = pval(F(perm_rows(yt, B, rng)))           # permute reduced residuals
    out["ter Braak"] = pval(F(perm_rows(e_full, B, rng)))           # permute full residuals
    return out, F_obs

rng = np.random.default_rng(2307)
y = data["STACKLOSS"].to_numpy()
x = data["ACIDCONC"].to_numpy()
Z = np.column_stack([np.ones(len(y)), data[["AIRFLOW", "WATERTEMP"]]])
pvals, F_obs = perm_pvalues(y, x, Z, B=19_999, rng=rng)
print(f"F = {F_obs:.3f}")
for m, pv in pvals.items():
    print(f"{m:16s} p = {pv:.4f}")
```

Real data with a p-value near \( 0.34 \) cannot separate the schemes. The simulation in
[Figure 23.5.1](#fig-bs-freedman-lane-size) is designed to. It has \( n=12 \) cases, normal errors, a nuisance
regressor \( z \) correlated with the regressor of interest \( x \), and one high-leverage value of \( x \). The null hypothesis holds, and the
nuisance coefficient \( \gamma \) ranges from \( 0 \) to \( 30 \). Each point is the rejection rate of a nominal \( 5\% \) test in
\( 5000 \) data sets, with \( 199 \) permutations each, so the Monte Carlo standard error is about
\( 0.0031 \).

::: {when-format="html"}
![**Figure 23.5.1.** Null rejection rates of nominal 5% tests of one coefficient, \( n=12 \), against the nuisance coefficient
\( \gamma \). The band is two Monte Carlo standard errors around 0.05.](freedman_lane_size.svg){#fig-bs-freedman-lane-size width=85%}
:::

::: {when-format="pdf"}
![Null rejection rates of nominal 5% tests of one coefficient, \( n=12 \), against the nuisance coefficient
\( \gamma \). The band is two Monte Carlo standard errors around 0.05.](freedman_lane_size.pdf){width=85%}
:::

The same random numbers are used for every \( \gamma \), so by @prp-bs-freedman-lane(b) five of the six tests make identical
decisions at every \( \gamma \), and their lines are flat. The \( F \) test rejects in \( 0.048 \) of the data sets, Freedman–Lane
in \( 0.054 \), ter Braak and Draper–Stoneman in \( 0.046 \), and Kennedy in \( 0.066 \), the liberal bias of
@prp-bs-kennedy made visible by the small sample. The raw permutation test is exact at \( \gamma=0 \) (\( 0.048 \)), but its
size moves with \( \gamma \), reaching \( 0.064 \) at \( \gamma=30 \).

```{.python .run #cell-freedman-lane-size}
def size_study(gamma, reps, B, seed, n=12, alpha=0.05):
    """Rejection rates at level alpha when beta = 0 and the nuisance coefficient is gamma."""
    rng = np.random.default_rng(seed)
    z = np.linspace(-1, 1, n)
    x = 0.8 * z + 0.6 * np.random.default_rng(1).normal(size=n)   # fixed, correlated with z
    x[0] = 6.0                                    # one high-leverage value of x
    Z = np.column_stack([np.ones(n), z])
    rejections = {}
    for _ in range(reps):
        y = 1.0 + gamma * z + rng.normal(size=n)  # the null hypothesis beta = 0 holds
        for m, pv in perm_pvalues(y, x, Z, B, rng)[0].items():
            rejections[m] = rejections.get(m, 0) + (pv <= alpha)
    return {m: r / reps for m, r in rejections.items()}

for gamma in (0.0, 30.0):
    print(gamma, size_study(gamma, reps=300, B=99, seed=2308))
```

As Anderson and Legendre (1999) concluded from a much larger simulation, Freedman–Lane is the sensible default: exact whenever
an exact test exists, blind to the nuisance coefficients, and close to the oracle in large samples. ter Braak's scheme behaves similarly, and since its reference distribution does not need the null hypothesis, it
suits the inversion of tests into intervals. Winkler and coauthors (2014) review the schemes for the general linear model.

## Exercises

### A. Check your understanding

::: {#exr-bs-fl-intercept-only}
[A1]

Show that when \( \Z=\bone \) the Freedman–Lane, Kennedy and raw permutation tests all coincide with the permutation test of
[Section 23.4](04-permutation-tests.html).
:::

::: {.solution}
With \( \Z=\bone \), every permutation fixes \( \Z \), so by @prp-bs-freedman-lane(c) Freedman–Lane and the raw permutation coincide and
are exact. Kennedy's scheme differs only through \( \norm{\Mo\bP\he_0}^2=n\,\overline{(\bP\he_0)}^2 \), which is zero because
\( \he_0=\y-\bar{y}\bone \) sums to zero and permuting does not change the sum.
:::

::: {#exr-bs-terbraak-identity}
[A2]

In ter Braak's scheme, what is the permuted statistic for \( \bP=\I \)? Explain why the observed statistic must still be counted in the
numerator and denominator of @eq-bs-mc-p.
:::

### B. Practice

::: {#exr-bs-permutation-moments}
[B1]

Let \( \mathbf{v}\in\Real^n \) be fixed and \( \bP \) a uniformly random permutation matrix. Show that \( \E(\bP\mathbf{v})=\bar{v}\bone \) and
\( \Cov(\bP\mathbf{v})=\frac{S_v}{n-1}\bigl(\I-n^{-1}\bone\bone\T\bigr) \) with \( S_v=\sum_i(v_i-\bar{v})^2 \). (The computation is
the one in the proof of @prp-dsn-randomization, there for permutations within a block.)
:::

### C. Going deeper

::: {#exr-bs-fl-heteroscedastic}
[C1]

Suppose the errors are independent with variances that depend on \( \X_1 \). Explain why the schemes of @def-bs-schemes need not be
asymptotically valid for \( H_0:\bbeta_1=\bzero \), and where the argument of @prp-bs-freedman-lane(d) breaks down. Propose a studentized
version in the spirit of [Section 23.4](04-permutation-tests.html).
:::
