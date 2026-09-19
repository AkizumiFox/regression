# Permutation tests for linear hypotheses

The bootstrap replaces an unknown distribution by an estimate. A permutation test instead looks for transformations of the
data that leave their distribution unchanged under the null hypothesis, and compares the observed statistic with its values on
the transformed data. When such transformations exist the test is exact at every sample size, whatever the error law. The idea
is Fisher's (1935), developed by Pitman (1937); in designed experiments the transformations are the rearrangements the
randomization could have produced ([Section 18.4](../ch18-covariance-and-design/04-randomized-blocks.html)).

## The general construction

Let \( \mathcal{G} \) be a finite group of \( n\times n \) permutation matrices. A random vector \( \Y \) is
**\( \mathcal{G} \)-invariant in distribution** if \( \bP\Y \) has the distribution of \( \Y \) for every \( \bP\in\mathcal{G} \). For the
group of all permutations this is exchangeability.

::: {#def-bs-permutation-test}
[Permutation test]

Let \( T:\Real^n\to\Real \) be a statistic, large values of which count against the null hypothesis.

::: {.enumerate options="label=(\alph*)"}
1. The **permutation p-value** is
   \( p(\y)=\lvert\mathcal{G}\rvert^{-1}\,\#\{\bP\in\mathcal{G}:T(\bP\y)\ge T(\y)\} \).

2. The **Monte Carlo permutation p-value** draws \( \bP_1,\dots,\bP_B \) independently and uniformly from \( \mathcal{G} \),
   independently of \( \Y \), and is
   \[
p_B(\y)=\frac{1+\#\{b:T(\bP_b\y)\ge T(\y)\}}{B+1} .
\]{#eq-bs-mc-p}

:::

The test of level \( \alpha \) rejects when the p-value is at most \( \alpha \).
:::

The \( 1 \) in @eq-bs-mc-p counts the observed arrangement among those compared; without it the p-value can be zero and the
test is not exact (Phipson and Smyth 2010). Exactness rests on a counting fact.

::: {#lem-bs-rank}
[Ranks of exchangeable values]

For real numbers \( t_0,\dots,t_m \), let \( R_j=\#\{i:t_i\ge t_j\} \). Then \( \#\{j:R_j\le k\}\le k \) for every \( k \). Consequently, if
\( (T_0,\dots,T_m) \) is an exchangeable random vector and \( R_0=\#\{i:T_i\ge T_0\} \), then \( \Pr(R_0\le k)\le k/(m+1) \).
:::

::: {.proof}
Let \( S=\{j:R_j\le k\} \), and suppose it is not empty. Choose \( j^*\in S \) with the smallest \( t_{j^*} \). Every \( j\in S \) has
\( t_j\ge t_{j^*} \), so \( R_{j^*}\ge\lvert S\rvert \), and \( \lvert S\rvert\le R_{j^*}\le k \). For the second claim, exchangeability gives
\( \Pr(R_j\le k)=\Pr(R_0\le k) \) for every \( j \), so
\( (m+1)\Pr(R_0\le k)=\E\,\#\{j:R_j\le k\}\le k \).
:::

::: {#thm-bs-permutation}
[Exactness of permutation tests]

If \( \Y \) is \( \mathcal{G} \)-invariant in distribution, then for every statistic \( T \), every \( \alpha\in[0,1] \) and every \( B\ge1 \),
\[
\Pr\bigl(p(\Y)\le\alpha\bigr)\le\alpha,\qquad \Pr\bigl(p_B(\Y)\le\alpha\bigr)\le\alpha .
\]
:::

::: {.proof}
Let \( \boldsymbol{\Gamma} \) be uniform on \( \mathcal{G} \), independent of everything else. Invariance gives
\( \boldsymbol{\Gamma}\Y\overset{d}{=}\Y \), since this holds conditionally on each value of \( \boldsymbol{\Gamma} \).

*Full group.* \( p(\Y) \) has the distribution of \( p(\boldsymbol{\Gamma}\Y) \). As \( \bP \) runs over \( \mathcal{G} \) so does
\( \bP\boldsymbol{\Gamma} \), so
\( p(\boldsymbol{\Gamma}\Y)=\lvert\mathcal{G}\rvert^{-1}\#\{\bP\in\mathcal{G}:T(\bP\Y)\ge T(\boldsymbol{\Gamma}\Y)\} \). Given \( \Y \), list the
values \( t_{\bP}=T(\bP\Y) \), \( \bP\in\mathcal{G} \). Then \( \lvert\mathcal{G}\rvert\,p(\boldsymbol{\Gamma}\Y)=R_{\boldsymbol{\Gamma}} \) in the notation
of @lem-bs-rank, with \( \boldsymbol{\Gamma} \) uniform over the \( \lvert\mathcal{G}\rvert \) indices. By the first part of the lemma,
\( \Pr\bigl(R_{\boldsymbol{\Gamma}}\le k\mid\Y\bigr)\le k/\lvert\mathcal{G}\rvert \). With \( k=\lfloor\alpha\lvert\mathcal{G}\rvert\rfloor \) this gives
\( \Pr(p\le\alpha\mid\Y)\le\alpha \), and we take expectations.

*Monte Carlo.* Put \( \bP_0=\I \) and \( T_b=T(\bP_b\Y) \), so that \( (B+1)\,p_B=R_0 \). We show \( (T_0,\dots,T_B) \) is exchangeable. Since
\( \Y \) is independent of \( (\bP_1,\dots,\bP_B) \) and \( \boldsymbol{\Gamma}\Y\overset{d}{=}\Y \) with \( \boldsymbol{\Gamma}\Y \) also independent of them,
\( (T_0,\dots,T_B) \) has the distribution of
\( \bigl(T(\boldsymbol{\Gamma}\Y),T(\bP_1\boldsymbol{\Gamma}\Y),\dots,T(\bP_B\boldsymbol{\Gamma}\Y)\bigr) \). The map
\( (\mathbf{G}_0,\mathbf{G}_1,\dots,\mathbf{G}_B)\mapsto(\mathbf{G}_0,\mathbf{G}_1\mathbf{G}_0,\dots,\mathbf{G}_B\mathbf{G}_0) \) is a bijection of
\( \mathcal{G}^{B+1} \), so \( \boldsymbol{\Gamma},\bP_1\boldsymbol{\Gamma},\dots,\bP_B\boldsymbol{\Gamma} \) are independent and uniform on \( \mathcal{G} \),
and independent of \( \Y \). Given \( \Y \), the \( B+1 \) values are therefore independent and identically distributed, hence
exchangeable, and so they are unconditionally. The second part of @lem-bs-rank with \( k=\lfloor\alpha(B+1)\rfloor \) completes the
proof.
:::

Nothing about \( T \) was used: the statistic affects power, never validity. The full-group test has size exactly \( \alpha \) when
\( \alpha\lvert\mathcal{G}\rvert \) is an integer and the values \( T(\bP\Y) \) are almost surely distinct. In the Monte Carlo test a
draw can repeat the identity or another draw, and the resulting ties make \( R_0 \) not exactly uniform, so its size is at most
\( \alpha \) and approaches \( \alpha \) as \( \lvert\mathcal{G}\rvert \) grows, when \( \alpha(B+1) \) is an integer. And the theorem says nothing about \( \bbeta \). The
hypothesis it tests is the invariance, and translating a regression hypothesis into an invariance is the whole difficulty.

## Simple regression

In the straight-line model with fixed \( x_i \) and independent, identically distributed errors, \( \beta_1=0 \) says
\( \Y=\beta_0\bone+\be \), which is exchangeable, so a permutation test of \( \beta_1=0 \) is exact for any error law. The usual
statistics give the same test.

::: {#prp-bs-simple-equivalence}
[Equivalent statistics for the slope]

With \( \mathcal{G} \) the group of all permutations, the statistics \( \lvert\hat{\beta}_1\rvert \), \( \lvert r\rvert \), \( \lvert t\rvert \)
and \( F=t^2 \) for the simple regression of \( \y \) on \( \x \) give the same permutation p-value, for every data set and, in the
Monte Carlo version, for every draw of permutations.
:::

::: {.proof}
Permuting \( \y \) leaves \( \bar{y} \) and \( S_{yy}=\sum_i(y_i-\bar{y})^2 \) unchanged, and \( S_{xx} \) does not involve \( \y \). Now
\( \hat{\beta}_1=S_{xy}/S_{xx} \) and \( r=S_{xy}/\sqrt{S_{xx}S_{yy}} \), so \( \lvert\hat{\beta}_1\rvert \) and \( \lvert r\rvert \) are
the same increasing function of \( \lvert S_{xy}\rvert \) for every rearrangement. And
\( t^2=(n-2)r^2/(1-r^2) \), which increases with \( r^2 \) on \( [0,1) \). So all four order the arrangements \( \bP\y \) in the same way,
and the counts in @def-bs-permutation-test agree.
:::

::: {#exm-bs-stackloss-permutation}
[Acid concentration and stack loss]

Brownlee's stack loss data (@exm-opt-stackloss-mle) record \( 21 \) days of operation of a plant oxidizing ammonia. Regress
stack loss on the acid concentration of the absorbing liquid alone. The slope is \( 0.7590 \), with
\( t=1.901 \) and a two-sided \( t \)-test p-value of \( 0.0725 \). The \( 5.11\times 10^{19} \) arrangements are too many to enumerate; Monte Carlo with
\( B=19999 \) permutations gives the p-value \( 0.0698 \), with Monte Carlo standard error
\( 0.0018 \). [Figure 23.4.1](#fig-bs-stackloss-permutation) shows the permutation distribution of \( t \) with the
\( t(19) \) density. The agreement is close, but the permutation p-value assumes nothing about the error law.
:::

::: {when-format="html"}
![**Figure 23.4.1.** Permutation distribution of \( t \) for stack loss on acid concentration, with the \( t(19) \) density and the
observed \( \pm t \).](stackloss_permutation.svg){#fig-bs-stackloss-permutation width=80%}
:::

::: {when-format="pdf"}
![Permutation distribution of \( t \) for stack loss on acid concentration, with the \( t(19) \) density and the
observed \( \pm t \).](stackloss_permutation.pdf){width=80%}
:::

```{.python .run #cell-permutation-stackloss}
import numpy as np
import statsmodels.api as sm
from scipy import stats
data = sm.datasets.stackloss.load_pandas().data

rng = np.random.default_rng(2304)
y = data["STACKLOSS"].to_numpy()
x = data["ACIDCONC"].to_numpy()
n = len(y)
xc = x - x.mean()

def slope_stats(Y):
    """Slope, correlation and t statistic for each row of Y, regressed on the fixed x."""
    Yc = Y - Y.mean(axis=-1, keepdims=True)
    b = Yc @ xc / (xc @ xc)
    r = Yc @ xc / np.sqrt((xc @ xc) * np.sum(Yc ** 2, axis=-1))
    t = r * np.sqrt(n - 2) / np.sqrt(1 - r ** 2)
    return b, r, t

b_obs, r_obs, t_obs = slope_stats(y)
p_t = 2 * stats.t.sf(abs(t_obs), n - 2)

B = 19_999
Yperm = y[np.argsort(rng.random((B, n)), axis=1)]  # B random permutations of y
b_perm, r_perm, t_perm = slope_stats(Yperm)
p_perm = (1 + np.sum(np.abs(b_perm) >= abs(b_obs))) / (B + 1)
print(f"slope {b_obs:.4f}, t = {t_obs:.3f}, t-test p = {p_t:.4f}, permutation p = {p_perm:.4f}")
```

## What the permutation test tests

Without identically distributed errors, invariance and \( \beta_1=0 \) come apart. If independent cases have
\( \E(Y\mid X)=0 \) but \( \Var(Y\mid X) \) depending on \( X \), the slope is zero but the responses are not exchangeable given the
\( x \)'s. The permutation test tests independence, and rejects too often as a test of zero slope. By the central limit theorem for the sample covariance, when \( \Cov(X,Y)=0 \),
\( \sqrt n\,r\to\Normal(0,\kappa^2) \) with
\[
\kappa^2=\frac{\E\bigl[(X-\mu_X)^2(Y-\mu_Y)^2\bigr]}{\Var(X)\Var(Y)} ,
\]
which equals \( 1 \) under independence but not in general. The permutation distribution of \( \sqrt n\,r \) tends to \( \Normal(0,1) \)
under mild conditions whatever the dependence, because permuting creates independence; we do not prove this permutation central
limit theorem (see Lehmann and Romano 2005, chapter 15). Let \( Y=X\varepsilon \), \( X \) and \( \varepsilon \) independent standard normals. Then \( \Var(Y)=1 \) and
\( \E(X^2Y^2)=\E X^4\,\E\varepsilon^2=3 \), so \( \kappa^2=3 \). A nominal \( 5\% \) permutation test based on \( \lvert r\rvert \), or by
@prp-bs-simple-equivalence on \( \lvert\hat{\beta}_1\rvert \), has limiting size
\[
\Pr\bigl(\lvert\Normal(0,3)\rvert>1.96\bigr)=2\bigl(1-\Phi(1.96/\sqrt3)\bigr)=2\bigl(1-\Phi(1.1316)\bigr)=0.258 .
\]
A simulation with \( 4000 \) data sets and \( 399 \) permutations each gives \( 0.236 \) at
\( n=20 \) and \( 0.248 \) at \( n=200 \).

The remedy (DiCiccio and Romano 2017) is to studentize: \( \hat{\beta}_1/\widehat{\text{se}}_{\text{HC0}} \) has a standard normal
limit under zero slope whether or not \( X \) and \( Y \) are independent (@exr-bs-studentized-limit), and so does its permutation
distribution. The test is exact under independence and asymptotically valid for zero slope. In the simulation its size is \( 0.117 \) at \( n=20 \) and
\( 0.056 \) at \( n=200 \): much better, though at \( n=20 \) still far from exact.

```{.python .run #cell-permutation-hetero}
import numpy as np

def perm_sizes(n, reps, B, rng, alpha=0.05):
    """Rejection rates of permutation tests of zero slope when y = x * eps (x, eps iid normal)."""
    reject = {"raw": 0, "studentized": 0}
    for _ in range(reps):
        x = rng.normal(size=n)
        y = x * rng.normal(size=n)                # E(y | x) = 0, Var(y | x) = x^2
        Y = np.vstack([y, y[np.argsort(rng.random((B, n)), axis=1)]])   # row 0: the data
        xc = x - x.mean()
        Yc = Y - Y.mean(axis=1, keepdims=True)
        b = Yc @ xc / (xc @ xc)
        e = Yc - b[:, None] * xc
        se_hc0 = np.sqrt((e ** 2) @ xc ** 2) / (xc @ xc)
        for name, T in [("raw", np.abs(b)), ("studentized", np.abs(b / se_hc0))]:
            p = np.mean(T >= T[0])                # (1 + #{T_b >= T_0}) / (B + 1)
            reject[name] += p <= alpha
    return {k: v / reps for k, v in reject.items()}

rng = np.random.default_rng(2305)
for n_sim in (20, 200):
    print(n_sim, perm_sizes(n_sim, reps=400, B=199, rng=rng))
```

## Nuisance regressors

Most regression hypotheses concern some coefficients in the presence of others:
\[
\Y=\Z\bgamma+\X_1\bbeta_1+\be,\qquad H_0:\bbeta_1=\bzero,
\]{#eq-bs-partial-model}

with \( \Z \) of size \( n\times p_0 \) containing the nuisance regressors (including the intercept) and \( \X_1 \) of size \( n\times q \).
Under \( H_0 \), \( \Y=\Z\bgamma+\be \), which is not exchangeable unless \( \Z\bgamma \) is a multiple of \( \bone \). There is an exact test
only when some group of permutations leaves the nuisance part alone.

::: {#prp-bs-invariant-nuisance}
[Permutations that fix the nuisance regressors]

Suppose \( \be \) is \( \mathcal{G} \)-invariant in distribution, and \( \bP\Z=\Z \) for every \( \bP\in\mathcal{G} \). Then under \( H_0 \) of
@eq-bs-partial-model the response \( \Y \) is \( \mathcal{G} \)-invariant in distribution, and the permutation test over \( \mathcal{G} \) with
any statistic is exact, whatever the value of \( \bgamma \).
:::

::: {.proof}
\( \bP\Y=\bP\Z\bgamma+\bP\be=\Z\bgamma+\bP\be\overset{d}{=}\Z\bgamma+\be=\Y \). Apply @thm-bs-permutation.
:::

The condition \( \bP\Z=\Z \) says that \( \bP \) permutes only cases with identical rows of \( \Z \). For block indicators, \( \mathcal{G} \)
is the group of permutations within blocks, and the proposition gives, under errors exchangeable within blocks, the same test as the design-based
randomization test of a randomized complete block design
([Section 18.4](../ch18-covariance-and-design/04-randomized-blocks.html), @def-dsn-rcbd), which needs no error model at all. When \( \Z \) contains a continuous covariate with distinct values only the identity qualifies, and no exact
permutation test exists. The schemes of the next section give up exactness for generality.

## Exercises

### A. Check your understanding

::: {#exr-bs-smallest-p}
[A1]

With \( B=99 \) random permutations, what is the smallest possible value of @eq-bs-mc-p? For which \( \alpha \) between \( 0.01 \) and
\( 0.10 \) is the size of the Monte Carlo test exactly \( \alpha \) when there are no ties (ignoring the chance of drawing the
same arrangement twice)?
:::

::: {.solution}
The smallest value is \( 1/100 \), reached when no permuted statistic reaches the observed one. When there are no ties
(ignoring the chance of drawing the same arrangement twice),
\( R_0 \) is uniform on \( \{1,\dots,100\} \), so \( \Pr(p_B\le\alpha)=\lfloor100\alpha\rfloor/100 \), which equals \( \alpha \) exactly when
\( 100\alpha \) is an integer: \( \alpha=0.01,0.02,\dots,0.10 \).
:::

### B. Practice

::: {#exr-bs-exact-enumeration-perm}
[B1]

For \( x=(1,2,3,4) \) and \( y=(1.2,\ 0.8,\ 2.9,\ 3.1) \), enumerate all \( 24 \) arrangements of \( y \) and compute the exact two-sided
permutation p-value for the slope. Compare it with the \( t \)-test p-value, and explain why the permutation p-value cannot be
smaller than \( 2/24 \) for any data with \( n=4 \) and distinct values.
:::

::: {.solution}
By @prp-bs-simple-equivalence it suffices to rank the arrangements by \( \lvert S_{xy}\rvert \), where
\( S_{xy}=\sum_i(x_i-\bar{x})y_{\pi(i)} \) with centred \( x \) equal to \( (-1.5,-0.5,0.5,1.5) \). For the data,
\( S_{xy}=-1.8-0.4+1.45+4.65=3.9 \). The largest values of \( S_{xy} \) come from putting the two largest responses at
\( x=3,4 \) and the two smallest at \( x=1,2 \); the four such arrangements give \( 3.7 \), \( 3.9 \), \( 4.1 \) and \( 4.3 \), and every
other arrangement gives at most \( 2.6 \). Reversing an arrangement negates \( S_{xy} \). So \( 6 \) of the \( 24 \) arrangements have
\( \lvert S_{xy}\rvert\ge3.9 \), and the p-value is \( 6/24=0.25 \). The \( t \)-test gives \( 0.139 \): with four points the
permutation distribution is coarse. The observed arrangement and its reversal always count, so \( p\ge2/24 \).
:::

::: {#exr-bs-kappa}
[B2]

Compute \( \kappa^2 \) when \( Y=\varepsilon\sqrt{1+\lambda X^2} \) with \( X \) and \( \varepsilon \) independent standard normals, and the
limiting size of the nominal \( 5\% \) permutation test of zero slope as a function of \( \lambda \). What happens as
\( \lambda\to\infty \)?
:::

### C. Going deeper

::: {#exr-bs-studentized-limit}
[C1]

For independent pairs with \( \Cov(X,Y)=0 \) and finite fourth moments, show that
\( \hat{\beta}_1/\widehat{\text{se}}_{\text{HC0}}\to\Normal(0,1) \) in distribution. (Use the central limit theorem for
\( n^{-1/2}\sum_i(X_i-\mu_X)(Y_i-\mu_Y) \) and Slutsky's lemma.)
:::
