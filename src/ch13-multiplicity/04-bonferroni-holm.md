# Bonferroni and Holm

Scheffé's and Tukey's methods need an exact distribution for the worst error in a family. The
methods of this section need only the level of each separate test and a probability inequality, so
they apply to any finite list of hypotheses in any model, whatever the dependence. Holm's step-down
procedure recovers part of the resulting conservatism at no cost in generality.

## The Bonferroni inequality

::: {#thm-mc-bonferroni}
[Bonferroni]

::: {.enumerate options="label=(\alph*)"}
1. For any events \( A_1,\dots,A_m \), \( \Pr\bigl(\bigcup_jA_j\bigr)\le\sum_j\Pr(A_j) \).

2. Let \( \alpha_1,\dots,\alpha_m>0 \) with \( \sum_j\alpha_j=\alpha \). If \( [L_j,U_j] \) is a confidence
   interval for \( \psi_j \) with coverage at least \( 1-\alpha_j \), for each \( j \), then the \( m \)
   intervals have simultaneous coverage at least \( 1-\alpha \).

3. Let \( p_j \) be a \( p \)-value for \( H_j \), valid in the sense that \( \Pr(p_j\le u)\le u \) for
   \( 0\le u\le1 \) whenever \( H_j \) is true. The procedure that rejects \( H_j \) when
   \( p_j\le\alpha_j \) has \( \text{PFER}\le\sum_{j\in I_0}\alpha_j\le\alpha \), and hence controls the familywise
   error rate strongly at level \( \alpha \).
:::

No assumption is made about the dependence among the events, intervals or \( p \)-values.
:::

::: {.proof}
*(a)* Pointwise, the indicator of a union is at most the sum of the indicators:
\( 1_{\bigcup_jA_j}\le\sum_j1_{A_j} \), since the right side is at least \( 1 \) wherever the left side is.
Take expectations. *(b)* Apply (a) to the events \( A_j=\{\psi_j\notin[L_j,U_j]\} \).
*(c)* The number of false rejections is \( V=\sum_{j\in I_0}1_{\{p_j\le\alpha_j\}} \), so
\( \E(V)=\sum_{j\in I_0}\Pr(p_j\le\alpha_j)\le\sum_{j\in I_0}\alpha_j \). Strong control follows from @prp-mc-error-rate-order.
:::

Part (c) bounds even the expected *number* of false rejections, which is why the procedure is
conservative when false rejections tend to occur together.

In the linear model, the **Bonferroni intervals** for \( k \) estimable functions
\( \blambda_1\T\bbeta,\dots,\blambda_k\T\bbeta \) chosen in advance are
\[
\blambda_j\T\hbeta\pm t_{\nu,\alpha/(2k)}\,\text{se}(\blambda_j\T\hbeta),\qquad j=1,\dots,k ,
\]{#eq-mc-bonferroni-interval}

each the \( t \) interval of @thm-ci-estimable-interval at level \( 1-\alpha/k \). The equal split is
conventional; any split fixed before the data are seen is valid.

::: {.warning}
The Bonferroni guarantee is for a list of \( k \) hypotheses fixed in advance. If the analyst
looks at the estimates and then chooses which \( k \) to test, the hypotheses actually tested are a
data-dependent selection from a larger family, and the relevant \( k \) is the size of that larger
family. For data-driven choices within a subspace, use Scheffé's method.
:::

## Independent tests: Šidák's improvement

If the tests of the true hypotheses are independent and each has exact level \( a \), then by @eq-mc-independent-rate
the familywise error rate is \( 1-(1-a)^{m_0} \). Choosing
\[
a=1-(1-\alpha)^{1/m}
\]{#eq-mc-sidak-level}

makes this at most \( 1-(1-a)^m=\alpha \), with equality when all hypotheses are true. This
**Šidák level** is slightly larger than \( \alpha/m \), because \( (1-\alpha/m)^m>1-\alpha \) for \( m\ge2 \)
(Bernoulli's inequality), so the Šidák
procedure is slightly less conservative than Bonferroni's.

In a linear model, however, \( t \) statistics are never independent: they share \( s \). The
following result shows that the sharing works in the right direction.

::: {#lem-mc-chebyshev}
[Chebyshev's association inequality]

Let \( U \) be a real random variable and \( f_1,\dots,f_k \) nonnegative, nondecreasing functions with
\( \E f_j(U)^k<\infty \). Then \( \E\bigl[\prod_jf_j(U)\bigr]\ge\prod_j\E f_j(U) \).
:::

::: {.proof}
For two functions \( f,h \), let \( U' \) be an independent copy of \( U \). Since \( f \) and \( h \) are both
nondecreasing, \( (f(U)-f(U'))(h(U)-h(U'))\ge0 \). Taking expectations gives
\( 2\E[f(U)h(U)]-2\E f(U)\,\E h(U)\ge0 \). For \( k \) functions, argue by induction: the product
\( f_1\cdots f_{k-1} \) of nonnegative nondecreasing functions is nonnegative and nondecreasing, so
\( \E[f_1\cdots f_k(U)]\ge\E[f_1\cdots f_{k-1}(U)]\,\E f_k(U) \), and the induction hypothesis
applies to the first factor. The moment condition makes every expectation finite, by Hölder's
inequality.
:::

::: {#prp-mc-sidak}
[Šidák intervals]

Let \( Z_1,\dots,Z_k \) be independent \( \Normal(0,1) \) variables and \( V\sim\chi^2(\nu) \) independent of
them, and put \( T_j=Z_j/\sqrt{V/\nu} \). Then for every \( c>0 \),
\[
\Pr\bigl(\lvert T_j\rvert\le c\ \text{for all }j\bigr)\ge\prod_{j=1}^k\Pr\bigl(\lvert T_j\rvert\le c\bigr).
\]
Consequently, if \( \blambda_1\T\hbeta,\dots,\blambda_k\T\hbeta \) are uncorrelated, the intervals
\( \blambda_j\T\hbeta\pm c\,\text{se}(\blambda_j\T\hbeta) \) with \( c \) the upper
\( \tfrac12\bigl(1-(1-\alpha)^{1/k}\bigr) \) point of \( t(\nu) \) have simultaneous coverage at least
\( 1-\alpha \).
:::

::: {.proof}
Conditionally on \( V=v \), the events \( \{\lvert Z_j\rvert\le c\sqrt{v/\nu}\} \) are independent, each with
probability \( f(v)=2\Phi(c\sqrt{v/\nu})-1 \). So
\( \Pr(\lvert T_j\rvert\le c\ \forall j)=\E\bigl[f(V)^k\bigr] \), while \( \Pr(\lvert T_j\rvert\le c)=\E f(V) \). The
function \( f \) is nonnegative, bounded and nondecreasing, and @lem-mc-chebyshev gives
\( \E[f(V)^k]\ge(\E f(V))^k \). For the intervals, the estimates are jointly normal (@thm-opt-sampling),
so uncorrelated estimates are independent (@thm-mvn-independence), and they
are independent of \( s \). The standardized errors are therefore of the form \( T_j \), and each
\( \Pr(\lvert T_j\rvert\le c)=(1-\alpha)^{1/k} \).
:::

The proposition covers orthogonal contrasts (@prp-est-contrast-ss(c)) and regressions with orthogonal
columns. Šidák (1967) proved the corresponding inequality
\( \Pr(\lvert Z_j\rvert\le c_j\ \forall j)\ge\prod_j\Pr(\lvert Z_j\rvert\le c_j) \) for *any* centred jointly normal
\( (Z_1,\dots,Z_k) \); with the conditioning argument above it makes Šidák's intervals valid for any
finite family of estimable functions. We use that theorem without proof. The gain over Bonferroni is
small: for \( k=10 \) and \( \nu=29 \) the multipliers are \( 3.0289 \) and \( 3.0380 \), and with independent
numerators and a shared \( s \) the exact coverage of the nominal 95% intervals, computed by integrating
over \( s \), is \( 0.9520 \) and \( 0.9530 \).

## The closure principle

Bonferroni's procedure keeps the threshold \( \alpha/m \) even after rejections, although a rejected
false hypothesis can no longer be falsely rejected. The closure principle exploits this. For each
nonempty \( I\subseteq\{1,\dots,m\} \) let \( H_I=\bigcap_{i\in I}H_i \) be the **intersection hypothesis**, and let
\( \phi_I \) be a level-\( \alpha \) **local test** of it.

::: {#thm-mc-closed-testing}
[Closed testing]

The procedure that rejects \( H_j \) iff \( \phi_I \) rejects \( H_I \) for every \( I \) containing \( j \)
controls the familywise error rate strongly at level \( \alpha \).
:::

::: {.proof}
Fix the parameters and suppose \( I_0\ne\emptyset \); otherwise there is nothing to prove. The
intersection hypothesis \( H_{I_0} \) is true, so \( \Pr(\phi_{I_0}\text{ rejects})\le\alpha \). If some true
\( H_j \) is rejected, then \( j\in I_0 \), so by the definition of the procedure \( \phi_{I_0} \) rejects.
Hence \( \Pr(V\ge1)\le\Pr(\phi_{I_0}\text{ rejects})\le\alpha \).
:::

Rejecting one hypothesis requires rejecting every intersection containing it, in general
\( 2^{m-1} \) tests. Holm's procedure is a case where this collapses to \( m \) steps.

## Holm's step-down procedure

Order the \( p \)-values as \( p_{(1)}\le p_{(2)}\le\dots\le p_{(m)} \), breaking ties arbitrarily, and
let \( H_{(i)} \) be the hypothesis with the \( i \)th smallest \( p \)-value.

::: {.algorithm}
**Holm's procedure at level \( \alpha \).** For \( i=1,2,\dots,m \) in turn: if
\( p_{(i)}\le\alpha/(m-i+1) \), reject \( H_{(i)} \) and continue; otherwise stop, and reject none of
\( H_{(i)},\dots,H_{(m)} \).
:::

The first threshold is Bonferroni's, \( \alpha/m \). Each rejection relaxes the next threshold: to
\( \alpha/(m-1) \), then \( \alpha/(m-2) \), and so on, up to \( \alpha \) for the last hypothesis.

::: {#thm-mc-holm}
[Holm]

Let each \( p_j \) be valid for \( H_j \), as in @thm-mc-bonferroni(c), with any dependence among them. Then:

::: {.enumerate options="label=(\alph*)"}
1. Holm's procedure controls the familywise error rate strongly at level \( \alpha \);

2. it rejects every hypothesis that Bonferroni's procedure (all \( \alpha_j=\alpha/m \)) rejects, and
   possibly more.
:::

:::

::: {.proof}
*(a)* Suppose \( m_0\ge1 \), and let \( j^* \) be the rank, in the ordering of all \( m \) \( p \)-values, of the
smallest \( p \)-value belonging to a true hypothesis, so that \( p_{(j^*)}=\min_{i\in I_0}p_i \). The
hypotheses ranked before \( j^* \) are all false, and there are only \( m-m_0 \) false hypotheses, so
\( j^*\le m-m_0+1 \) and \( \alpha/(m-j^*+1)\le\alpha/m_0 \). Holm's procedure rejects an initial segment
\( H_{(1)},\dots,H_{(i)} \) of the ordering. If that segment contains any true hypothesis, it contains
\( H_{(j^*)} \), so \( p_{(j^*)}\le\alpha/(m-j^*+1)\le\alpha/m_0 \). Therefore
\[
\Pr(V\ge1)\le\Pr\Bigl(\min_{i\in I_0}p_i\le\frac{\alpha}{m_0}\Bigr)\le\sum_{i\in I_0}\Pr\Bigl(p_i\le\frac{\alpha}{m_0}\Bigr)\le\alpha,
\]
by @thm-mc-bonferroni(a) and the validity of the \( p \)-values.
*(b)* If Bonferroni rejects \( H_{(i)} \), then \( p_{(j)}\le p_{(i)}\le\alpha/m\le\alpha/(m-j+1) \) for every
\( j\le i \), so Holm's procedure does not stop before step \( i \) and rejects \( H_{(i)} \).
:::

Bonferroni's \( \alpha/m \) pays for \( m \) potential errors when only \( m_0 \) are possible. Holm's procedure is closed testing with the Bonferroni local tests (@exr-mc-holm-closed).
It gives no simple simultaneous intervals, since the threshold for a hypothesis
depends on the other \( p \)-values; its output is a list of rejections, or **adjusted \( p \)-values**
\[
\tilde p_{(i)}=\max_{j\le i}\,\min\bigl\{1,\ (m-j+1)\,p_{(j)}\bigr\},
\]{#eq-mc-holm-adjusted}

with \( H_{(i)} \) rejected at level \( \alpha \) iff \( \tilde p_{(i)}\le\alpha \) (@exr-mc-holm-adjusted).

::: {#exm-mc-holm-education}
[Holm for the education pairs]

For the \( 21 \) pairwise comparisons of @exm-mc-tukey-education the smallest \( p \)-value is
\( 1.1\times 10^{-5} \). Bonferroni's threshold is \( 0.05/21=0.00238 \), and \( 8 \) of the \( p \)-values fall below it.
Holm's procedure starts at the same threshold but relaxes it after each rejection. It rejects \( 10 \)
hypotheses: the tenth smallest \( p \)-value is \( 0.0037 \), below its threshold \( 0.05/12 \), and the
eleventh, \( 0.0126 \), exceeds \( 0.05/11=0.0045 \), so the procedure stops.
:::

```{.python .run #cell-education-stepwise}
import itertools
import numpy as np
import statsmodels.api as sm
from scipy import stats
anes = sm.datasets.anes96.load_pandas().data
y = anes["ClinLR"].to_numpy()
level = anes["educ"].to_numpy().astype(int) - 1       # 0, ..., 6
g = 7
n_k = np.bincount(level).astype(float)
means = np.bincount(level, weights=y) / n_k
n = len(y)
nu = n - g                                            # error degrees of freedom
s2 = np.sum((y - means[level]) ** 2) / nu             # pooled variance estimate
s = np.sqrt(s2)
alpha = 0.05
pairs = list(itertools.combinations(range(g), 2))
diff = np.array([means[k] - means[l] for k, l in pairs])
se_d = np.array([s * np.sqrt(1 / n_k[k] + 1 / n_k[l]) for k, l in pairs])
t = diff / se_d
pval = 2 * stats.t.sf(np.abs(t), nu)

def holm(p, alpha):
    """Holm's step-down procedure: reject the i smallest p-values, i maximal with
    p_(j) <= alpha / (m - j + 1) for all j <= i."""
    m = len(p)
    order = np.argsort(p)
    reject = np.zeros(m, dtype=bool)
    for j, idx in enumerate(order):                   # j = 0, 1, ...
        if p[idx] > alpha / (m - j):
            break
        reject[idx] = True
    return reject

def benjamini_hochberg(p, q):
    """Reject the k smallest p-values, k the largest index with p_(k) <= k q / m."""
    m = len(p)
    order = np.argsort(p)
    below = np.nonzero(p[order] <= q * np.arange(1, m + 1) / m)[0]
    reject = np.zeros(m, dtype=bool)
    if below.size:
        reject[order[: below[-1] + 1]] = True
    return reject

print("Holm rejects", holm(pval, alpha).sum(), "  BH (q = 0.05) rejects",
      benjamini_hochberg(pval, alpha).sum())
```

## How much does Holm gain?

If all hypotheses are true, Holm and Bonferroni make a false rejection in the same data sets. The more
hypotheses are false and rejected early, the larger the later thresholds and the gain.

::: {#exm-mc-holm-power}
[Ten regression coefficients]

A regression has an intercept and ten centred, mutually orthogonal regressors, with \( n=40 \). The ten
\( t \) statistics for the slopes are then \( T_j=Z_j/\sqrt{V/29} \) with independent
\( Z_j\sim\Normal(\delta_j,1) \) and a shared \( V\sim\chi^2(29) \) (@thm-opt-sampling). Let \( k \) of the slopes be
nonzero, each with \( \delta_j=3.5 \), and test all ten at familywise level \( 0.05 \).
[Figure 13.4.1](#fig-mc-holm-power) shows the result of \( 20000 \) simulated data sets for each \( k \). All
three procedures keep the familywise error rate at or below \( 0.05 \), up to simulation error. The average power of Bonferroni's
procedure, the proportion of nonzero slopes it detects, hardly depends on \( k \): it is \( 0.675 \) with
one nonzero slope and \( 0.676 \) with five or nine. Šidák's is larger by a few thousandths. Holm's power is \( 0.675 \) with one nonzero slope, \( 0.718 \) with five and
\( 0.794 \) with nine.
:::

::: {when-format="html"}
![**Figure 13.4.1.** Ten orthogonal regression coefficients, \( k \) of them nonzero, tested
at familywise level \( 0.05 \): (a) familywise error rate, (b) average power, the share of nonzero
coefficients detected.](holm_power.svg){#fig-mc-holm-power width=100%}
:::

::: {when-format="pdf"}
![Ten orthogonal regression coefficients, \( k \) of them nonzero, tested
at familywise level \( 0.05 \): (a) familywise error rate, (b) average power, the share of nonzero
coefficients detected.](holm_power.pdf){width=100%}
:::

```{.python .run #cell-holm-power-simulate}
import numpy as np
from scipy import stats

alpha, m, nu, delta, reps = 0.05, 10, 29, 3.5, 20_000

def simulate(k, reps, seed):
    """FWER and average power of three procedures when k of m coefficients are nonzero."""
    rng = np.random.default_rng(seed)
    shift = np.r_[np.full(k, delta), np.zeros(m - k)]
    T = (rng.normal(size=(reps, m)) + shift) / np.sqrt(rng.chisquare(nu, (reps, 1)) / nu)
    p = 2 * stats.t.sf(np.abs(T), nu)
    rules = {
        "Bonferroni": p <= alpha / m,
        "Sidak": p <= 1 - (1 - alpha) ** (1 / m),
    }
    # Holm: sort each row, compare p_(j) with alpha / (m - j + 1), stop at the first failure
    order = np.argsort(p, axis=1)
    ps = np.take_along_axis(p, order, axis=1)
    passed = np.cumprod(ps <= alpha / (m - np.arange(m)), axis=1).astype(bool)
    holm = np.zeros_like(passed)
    np.put_along_axis(holm, order, passed, axis=1)
    rules["Holm"] = holm
    out = {}
    for name, rej in rules.items():
        fwer = np.mean(rej[:, k:].any(axis=1)) if k < m else np.nan
        power = rej[:, :k].mean() if k > 0 else np.nan
        out[name] = (fwer, power, rej)
    return out

res = simulate(5, reps, seed=5)
for name, (fwer, power, _) in res.items():
    print(f"{name:10s} FWER = {fwer:.3f}   average power = {power:.3f}")
```

There is no reason to use Bonferroni's *tests* rather than Holm's, though Bonferroni's *intervals*
remain useful. With independent \( p \)-values, Hochberg (1988) proposed the
step-*up* version: reject \( H_{(1)},\dots,H_{(i)} \) for the largest \( i \) with \( p_{(i)}\le\alpha/(m-i+1) \). It
rejects everything Holm's procedure does (@exr-mc-hochberg-holm), but its validity rests on an
inequality of Simes (1986) that holds under independence and some positive dependence, not in
general. The step-up idea returns in the next section, controlling a different error rate.

## Exercises

### A. Check your understanding

::: {#exr-mc-five-pvalues}
[A1]

Five hypotheses have \( p \)-values \( 0.001,\ 0.011,\ 0.012,\ 0.04,\ 0.2 \). At familywise level \( 0.05 \),
which are rejected by Bonferroni's procedure and which by Holm's? Compute Holm's adjusted
\( p \)-values @eq-mc-holm-adjusted.
:::

::: {.solution}
Bonferroni's threshold is \( 0.01 \), so only the first is rejected. Holm compares
\( 0.001\le0.05/5 \), \( 0.011\le0.05/4 \), \( 0.012\le0.05/3 \), and then \( 0.04>0.05/2 \), so it rejects the first
three and stops. The adjusted \( p \)-values are \( 5\times0.001=0.005 \), \( \max(0.005,4\times0.011)=0.044 \),
\( \max(0.044,3\times0.012)=0.044 \), \( \max(0.044,2\times0.04)=0.08 \) and \( \max(0.08,0.2)=0.2 \). The three
below \( 0.05 \) are the rejected hypotheses.
:::

### B. Practice

::: {#exr-mc-holm-adjusted}
[B1]

Show that Holm's procedure rejects \( H_{(i)} \) iff the adjusted \( p \)-value @eq-mc-holm-adjusted
satisfies \( \tilde p_{(i)}\le\alpha \), and that the adjusted \( p \)-values are nondecreasing in \( i \).
:::

::: {#exr-mc-holm-pfer}
[B2]

All \( m\ge2 \) hypotheses are true and the \( p \)-values are independent and uniform. Show that the
per-family error rate of Bonferroni's procedure is exactly \( \alpha \), while that of Holm's procedure
is larger than \( \alpha \). So Holm's procedure controls the familywise rate but not the per-family
rate.
:::

::: {#exr-mc-holm-closed}
[B3]

Show that Holm's procedure is closed testing (@thm-mc-closed-testing) with the Bonferroni local
tests \( \phi_I \): reject \( H_I \) when \( \min_{i\in I}p_i\le\alpha/\lvert I\rvert \).
:::

::: {.solution}
Suppose Holm rejects \( H_{(1)},\dots,H_{(i)} \) and stops at step \( i+1 \). Take \( j\le i \) and any \( I\ni(j) \),
writing \( (j) \) for the index of \( H_{(j)} \). If \( I \) contains some \( (l) \) with \( l\le j \), let \( l \) be the smallest
such; then \( I \) contains no \( (1),\dots,(l-1) \), so \( \lvert I\rvert\le m-l+1 \) and
\( \min_Ip\le p_{(l)}\le\alpha/(m-l+1)\le\alpha/\lvert I\rvert \): \( \phi_I \) rejects. So closed testing rejects
\( H_{(j)} \). Conversely, the set \( I=\{(i+1),\dots,(m)\} \) has \( m-i \) elements and minimum \( p_{(i+1)}>\alpha/(m-i) \),
so \( \phi_I \) does not reject, and no hypothesis in \( I \) is rejected by closed testing. The two procedures
therefore reject the same hypotheses.
:::

### C. Going deeper

::: {#exr-mc-hochberg-holm}
[C1]

Show that Hochberg's step-up procedure rejects every hypothesis rejected by Holm's step-down
procedure with the same thresholds. Give an example with \( m=2 \) where Hochberg rejects both
hypotheses and Holm rejects neither.
:::

