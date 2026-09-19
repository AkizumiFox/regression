# False discovery rate control

With hundreds of hypotheses, as when a regression screens many candidate predictors, familywise
control becomes self-defeating: the threshold \( \alpha/m \) shrinks with \( m \), and power disappears
exactly when there is most to find. The false discovery rate of @def-mc-error-rates asks instead that,
on average, only a small *proportion* of the rejections be false. This section proves that the
procedure of Benjamini and Hochberg controls it for independent \( p \)-values, and describes what is
known under dependence.

## The Benjamini–Hochberg procedure

As before, \( p_{(1)}\le\dots\le p_{(m)} \) are the ordered \( p \)-values, \( H_{(i)} \) the corresponding
hypotheses, and \( q\in(0,1) \) is the target false discovery rate.

::: {.algorithm}
**Benjamini–Hochberg (BH) procedure at level \( q \).** Let
\[
\hat k=\max\Bigl\{k\in\{1,\dots,m\}:p_{(k)}\le\frac{kq}{m}\Bigr\},
\]
and reject \( H_{(1)},\dots,H_{(\hat k)} \). If no \( k \) qualifies, reject nothing.
:::

The thresholds rise linearly from Bonferroni's \( q/m \) to \( q \). The procedure is a **step-up**
procedure: it rejects everything up to the *last* \( p \)-value on or below the line, including smaller
\( p \)-values that lie above the line at their own rank ([Figure 13.5.1](#fig-mc-bh-stepup)).

Why is the line the right one? Suppose \( k \) hypotheses are rejected, all with \( p \)-values at most
\( t=kq/m \). The \( m_0 \) true hypotheses have uniform \( p \)-values, so about \( m_0t\le mt \) of them fall
below \( t \) by chance, and the proportion of false ones among the \( k \) rejections is about
\( mt/k=q \) at most. The BH procedure chooses the largest threshold at which this estimate of the
false discovery proportion is at most \( q \). The theorem makes the heuristic exact.

Two simple descriptions of the rejected set are useful. Write
\( N(t)=\#\{j:p_j\le t\} \) for the number of \( p \)-values at most \( t \). Since \( p_{(k)}\le t \) iff
\( N(t)\ge k \),
\[
\hat k=\max\Bigl\{k:N\bigl(kq/m\bigr)\ge k\Bigr\}.
\]{#eq-mc-bh-count}

And the rejected hypotheses are exactly those with \( p_j\le\hat kq/m \). Indeed, every rejected
\( p \)-value is at most \( p_{(\hat k)}\le\hat kq/m \). Conversely, if some \( p_j \) with \( j \) not rejected were at most
\( \hat kq/m \), then \( N(\hat kq/m)\ge\hat k+1 \), hence \( N\bigl((\hat k+1)q/m\bigr)\ge\hat k+1 \), and \( \hat k+1 \) would
qualify, contradicting the maximality of \( \hat k \).

## Control under independence

::: {#lem-mc-bh-leave-one-out}
[Replacing one \( p \)-value by zero]

Fix \( i \), and let \( \hat k^{(i)} \) be the number of rejections of the BH procedure applied to the
\( p \)-values with \( p_i \) replaced by \( 0 \). Then \( \hat k^{(i)} \) is a function of \( (p_j)_{j\ne i} \) alone, and for
each \( k=1,\dots,m \),
\[
\bigl\{H_i\text{ is rejected and }\hat k=k\bigr\}=\bigl\{p_i\le kq/m\text{ and }\hat k^{(i)}=k\bigr\}.
\]
:::

::: {.proof}
The first claim is clear. Let \( N^{(i)} \) count the modified \( p \)-values. Then \( N^{(i)}(t)=N(t) \)
for \( t\ge p_i \), because at such \( t \) the \( p \)-value \( p_i \) is counted in both. By the second
description above, \( H_i \) is rejected and \( \hat k=k \) iff \( p_i\le kq/m \) and \( \hat k=k \). Assume
\( p_i\le kq/m \). By @eq-mc-bh-count, \( \hat k=k \) iff \( N(kq/m)\ge k \) and \( N(k'q/m)<k' \) for every
\( k'>k \). These conditions involve \( N \) only at points \( t\ge kq/m\ge p_i \), where \( N \) and \( N^{(i)} \)
agree, so they hold iff the same conditions hold for \( N^{(i)} \), that is, iff \( \hat k^{(i)}=k \).
:::

::: {#thm-mc-bh}
[Benjamini–Hochberg]

Suppose that for each true hypothesis \( H_i \), \( i\in I_0 \), the \( p \)-value \( p_i \) is independent of the
other \( m-1 \) \( p \)-values and satisfies \( \Pr(p_i\le u)\le u \) for \( 0\le u\le1 \). Then the BH procedure
at level \( q \) has
\[
\text{FDR}\le\frac{m_0}{m}\,q\le q ,
\]
with equality in the first inequality when every \( p_i \), \( i\in I_0 \), is uniform on \( (0,1) \).
:::

::: {.proof}
The false discovery proportion is \( Q=\sum_{i\in I_0}1_{\{H_i\text{ rejected}\}}/\max(\hat k,1) \). A rejection
forces \( \hat k\ge1 \), so
\[
\text{FDR}=\sum_{i\in I_0}\sum_{k=1}^m\frac1k\Pr\bigl(H_i\text{ rejected},\ \hat k=k\bigr)
=\sum_{i\in I_0}\sum_{k=1}^m\frac1k\Pr\Bigl(p_i\le\frac{kq}{m}\Bigr)\Pr\bigl(\hat k^{(i)}=k\bigr),
\]
by @lem-mc-bh-leave-one-out and the independence of \( p_i \) from \( \hat k^{(i)} \). Now
\( \Pr(p_i\le kq/m)\le kq/m \), so the inner sum is at most
\( (q/m)\sum_{k=1}^m\Pr(\hat k^{(i)}=k)=q/m \), because \( \hat k^{(i)}\ge1 \) always (the replaced
\( p \)-value \( 0 \) lies below the first threshold). Summing over the \( m_0 \) true hypotheses gives
\( \text{FDR}\le m_0q/m \). With uniform null \( p \)-values every inequality is an equality.
:::

The false hypotheses' \( p \)-values may depend on each other in any way. With uniform null
\( p \)-values the false discovery rate is exactly \( \pi_0q \), \( \pi_0=m_0/m \). With \( m=100 \) independent \( z \)-tests, \( 90 \) of them
null and the rest shifted by \( 3 \), \( 20000 \) simulated repetitions estimate it at \( 0.0439 \), against
the exact \( 0.0450 \).

When all hypotheses are true the false discovery rate is the familywise error rate (@prp-mc-error-rate-order),
so BH is a level-\( q \) test of the global null. It rejects everything
Bonferroni's and Holm's procedures at level \( q \) reject (@exr-mc-bh-holm). But it does *not* control
the familywise error rate strongly: once many false hypotheses are rejected, the thresholds for the
rest approach \( q \).

## Dependent tests

In a regression the \( t \) statistics share \( s \), and their numerators are correlated unless \( \X \) has
orthogonal columns. Two results, stated here without proof, cover much of what is needed.

**Positive dependence.** Benjamini and Yekutieli (2001) showed that \( \text{FDR}\le m_0q/m \) remains true
if the \( p \)-values are *positively regression dependent on each one from the subset of true
hypotheses* (PRDS): for every true \( H_i \) and every increasing set \( D \) of \( p \)-value vectors,
\( \Pr\bigl((p_1,\dots,p_m)\in D\mid p_i=u\bigr) \) is nondecreasing in \( u \). One-sided tests based on jointly
normal statistics with nonnegative correlations are an example; two-sided tests of correlated
statistics need not be.

**Arbitrary dependence.** They also showed that under *any* dependence, BH at level \( q/c_m \) with
\( c_m=\sum_{i=1}^m1/i\approx\log m+0.577 \) has false discovery rate at most \( m_0q/m \) (@exr-mc-by). For
\( m=100 \), \( c_m=5.187 \), so a target of \( 0.05 \) becomes an effective level of \( 0.0096 \).

BH is routinely applied to regression \( t \) tests, and simulations, including the one below, find its
false discovery rate below \( q \). That is evidence, not a theorem; when a guarantee is essential, the
Benjamini–Yekutieli version is the safe choice.

## Adjusted p-values and q-values

As for Holm's procedure, the BH decisions for all levels at once are summarized by **adjusted
\( p \)-values**
\[
\tilde p_{(i)}=\min_{j\ge i}\,\min\Bigl\{1,\ \frac{m\,p_{(j)}}{j}\Bigr\},
\]{#eq-mc-bh-adjusted}

with \( H_{(i)} \) rejected at level \( q \) iff \( \tilde p_{(i)}\le q \) (@exr-mc-bh-adjusted). Since BH
controls the rate at \( \pi_0q \), it is conservative when many hypotheses are false. **Adaptive**
procedures estimate \( \pi_0 \), for instance by twice the proportion of \( p \)-values above one half, and
run BH at level \( q/\hat{\pi}_0 \). Storey (2002) developed this approach and defined the **\( q \)-value** of a hypothesis, the
smallest false discovery rate at which it would be rejected; with \( \hat{\pi}_0=1 \) the \( q \)-values are the
adjusted values @eq-mc-bh-adjusted. The gain is largest when \( \pi_0 \) is well below one.

## Screening regression coefficients

::: {#exm-mc-bh-screening}
[One hundred candidate regressors]

A simulated regression has \( n=250 \) observations, an intercept and \( 100 \) regressors with
independent standard normal entries, leaving \( \nu=149 \) error degrees of freedom. The first \( 10 \)
slopes are nonzero, each chosen so that the noncentrality of its \( t \) statistic is \( 3 \); the other
\( 90 \) are zero. The \( 100 \) two-sided \( t \) tests are run on one data set at level \( 0.05 \):

| procedure | rejections | of which false |
|---|---|---|
| unadjusted, \( p\le0.05 \) | \( 9 \) | \( 2 \) |
| Bonferroni | \( 1 \) | \( 0 \) |
| Holm | \( 1 \) | \( 0 \) |
| BH, \( q=0.05 \) | \( 7 \) | \( 1 \) |

Bonferroni's threshold is \( 0.05/100=0.0005 \), and only one \( p \)-value is below it. The BH line
reaches \( 7\times0.05/100=0.0035 \) at rank \( 7 \), and the seventh smallest \( p \)-value is \( 0.0034 \)
([Figure 13.5.1](#fig-mc-bh-stepup)). Six of the seven BH discoveries are real.
:::

::: {when-format="html"}
![**Figure 13.5.1.** The twenty smallest of \( 100 \) ordered \( p \)-values (log scale), with the BH
line \( kq/m \) and Bonferroni's threshold \( q/m \), \( q=0.05 \). Filled points are nonzero coefficients;
BH rejects the seven left of the vertical line.](bh_stepup.svg){#fig-mc-bh-stepup width=80%}
:::

::: {when-format="pdf"}
![The twenty smallest of \( 100 \) ordered \( p \)-values (log scale), with the BH
line \( kq/m \) and Bonferroni's threshold \( q/m \), \( q=0.05 \). Filled points are nonzero coefficients;
BH rejects the seven left of the vertical line.](bh_stepup.pdf){width=80%}
:::

```{.python .run #cell-bh-screening-design}
import numpy as np
from scipy import stats

q = 0.05

rng = np.random.default_rng(20260919)
n, m = 250, 100
X = np.column_stack([np.ones(n), rng.normal(size=(n, m))])
Qx, Rx = np.linalg.qr(X)
nu = n - m - 1                                        # 149 error degrees of freedom
c_diag = np.diag(np.linalg.inv(Rx.T @ Rx))[1:]        # [(X'X)^{-1}]_jj for the slopes

def t_statistics(Y):
    """t statistics of the 100 slopes for each column of Y."""
    coef = np.linalg.solve(Rx, Qx.T @ Y)
    resid = Y - X @ coef
    s2 = np.sum(resid ** 2, axis=0) / nu
    return coef[1:] / np.sqrt(np.outer(c_diag, s2))
```

```{.python .run #cell-bh-screening-procedures}
def bh(p, q):
    """Benjamini-Hochberg step-up at level q; returns a boolean rejection vector."""
    m = len(p)
    order = np.argsort(p)
    below = np.nonzero(p[order] <= q * np.arange(1, m + 1) / m)[0]
    reject = np.zeros(m, dtype=bool)
    if below.size:
        reject[order[: below[-1] + 1]] = True
    return reject

def holm(p, alpha):
    m = len(p)
    order = np.argsort(p)
    ok = np.cumprod(p[order] <= alpha / (m - np.arange(m))).astype(bool)
    reject = np.zeros(m, dtype=bool)
    reject[order] = ok
    return reject

procedures = {
    "unadjusted": lambda p: p <= q,
    "Bonferroni": lambda p: p <= q / m,
    "Holm": lambda p: holm(p, q),
    "BH": lambda p: bh(p, q),
}
```

```{.python .run #cell-bh-screening-one}
m1 = 10                                               # the first 10 slopes are nonzero
beta = np.zeros(m + 1)
beta[1:m1 + 1] = 3.0 * np.sqrt(c_diag[:m1])           # noncentrality 3 when sigma = 1
y = X @ beta + rng.normal(size=n)
p = 2 * stats.t.sf(np.abs(t_statistics(y[:, None])[:, 0]), nu)
for name, rule in procedures.items():
    rej = rule(p)
    print(f"{name:10s} rejects {rej.sum():3d}, of which {rej[m1:].sum()} are false discoveries")
```

A simulation with the same design, \( m_1=0,5,10,25,50 \) nonzero slopes and \( 2000 \) responses for
each, shows the average behaviour ([Figure 13.5.2](#fig-mc-fdr-simulation)). With \( m_1=10 \), the false
discovery rate is \( 0.327 \) without adjustment, \( 0.014 \) for Bonferroni and Holm and \( 0.046 \) for BH,
and the average power is \( 0.845 \), \( 0.289 \), \( 0.291 \) and \( 0.476 \). With \( m_1=50 \), BH's power rises
to \( 0.731 \) and its false discovery rate falls to \( 0.024 \), in line with \( \pi_0q \), while Bonferroni's
and Holm's power stay near \( 0.3 \). BH pays with familywise errors: at least one of its rejections is
false in \( 0.238 \) of the data sets with \( m_1=10 \) and \( 0.570 \) with \( m_1=50 \). With all slopes zero, its
false discovery rate and familywise error rate coincide, at \( 0.050 \).

::: {when-format="html"}
![**Figure 13.5.2.** False discovery rate (a, log scale) and average power (b) of four
procedures when \( m_1 \) of \( 100 \) regression slopes are nonzero, \( 2000 \) simulated data sets per
point, nominal level \( 0.05 \). The pale line in (a) marks \( 0.05 \).](fdr_simulation.svg){#fig-mc-fdr-simulation width=100%}
:::

::: {when-format="pdf"}
![False discovery rate (a, log scale) and average power (b) of four
procedures when \( m_1 \) of \( 100 \) regression slopes are nonzero, \( 2000 \) simulated data sets per
point, nominal level \( 0.05 \). The pale line in (a) marks \( 0.05 \).](fdr_simulation.pdf){width=100%}
:::

For the \( 21 \) education comparisons of @exm-mc-holm-education, BH at \( q=0.05 \) rejects \( 13 \)
hypotheses: the thirteenth smallest \( p \)-value is \( 0.0303\le13\times0.05/21=0.0310 \), and no later one
comes back below the line. The next section discusses how to choose between the two error rates.

## Exercises

### A. Check your understanding

::: {#exr-mc-bh-five}
[A1]

Apply the BH procedure at \( q=0.05 \) to the five \( p \)-values of @exr-mc-five-pvalues, and compute the
adjusted \( p \)-values @eq-mc-bh-adjusted.
:::

::: {.solution}
The thresholds are \( 0.01,0.02,0.03,0.04,0.05 \). From the top: \( 0.2>0.05 \), but \( 0.04\le0.04 \), so
\( \hat k=4 \) and the first four hypotheses are rejected, one more than Holm's procedure rejects. The
adjusted values are \( \min(5\times0.001,\dots)=0.005 \), then \( \min(5\times0.011/2,\ 5\times0.012/3,\ 5\times0.04/4,\ 1)=0.02 \),
\( 0.02 \), \( 0.05 \) and \( 0.2 \).
:::

::: {#exr-mc-fdr-not-fdp}
[A2]

A report says: “Using BH at \( q=0.1 \) we made \( 40 \) discoveries, so at most \( 4 \) of them are false.”
What is wrong with this statement? What would be a correct one?
:::

### B. Practice

::: {#exr-mc-bh-adjusted}
[B1]

Prove that BH at level \( q \) rejects \( H_{(i)} \) iff the adjusted value @eq-mc-bh-adjusted satisfies
\( \tilde p_{(i)}\le q \).
:::

::: {#exr-mc-bh-holm}
[B2]

Show that the BH procedure at level \( q \) rejects every hypothesis that Holm's procedure at level \( q \)
rejects. *Hint:* \( q/(m-i+1)\le iq/m \) for \( 1\le i\le m \).
:::

::: {.solution}
The hint holds because \( i(m-i+1)\ge m \) for \( 1\le i\le m \): the product of two positive integers summing
to \( m+1 \) is at least \( m \). If Holm rejects \( H_{(1)},\dots,H_{(i)} \), then \( p_{(i)}\le q/(m-i+1)\le iq/m \), so
\( \hat k\ge i \).
:::

::: {#exr-mc-bh-fwer}
[B3]

With \( m_0=2 \) true hypotheses with independent uniform \( p \)-values and \( m-2 \) false hypotheses whose
\( p \)-values are \( 0 \), find the familywise error rate of BH at level \( q \). Compare it with the false
discovery rate given by @thm-mc-bh.
:::

### C. Going deeper

::: {#exr-mc-by}
[C1]

Prove the Benjamini–Yekutieli bound: under arbitrary dependence, BH at level \( q/c_m \) has false
discovery rate at most \( m_0q/m \). *Hint:* write the false discovery rate as
\( \sum_{i\in I_0}\sum_k\Pr(p_i\in(\tfrac{(j-1)q'}{m},\tfrac{jq'}{m}],\ \hat k=k)/k \) summed over \( j\le k \), with
\( q'=q/c_m \), and exchange the order of summation.
:::
