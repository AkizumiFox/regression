# The multiple testing problem

Seven groups give \( 21 \) pairs, and a two-sided \( t \) test at level \( 0.05 \) for each pair is a
correct test. But an analyst who runs all \( 21 \) and reports the significant ones is not using a
procedure with error rate \( 0.05 \): even if all seven means are equal, about one pair is expected
to be significant, and the reader cannot tell it from a real difference. This section defines the
error rates of a family of inferences, measures by simulation how large they become without
adjustment, and separates strong control from a weak form that protects only against the
situation in which nothing is going on.

## A family of hypotheses

Let \( H_1,\dots,H_m \) be null hypotheses about the parameters of a model, tested on the same
data. The model's true parameter makes some of them true and the others false. Write
\( I_0\subseteq\{1,\dots,m\} \) for the set of true hypotheses and \( m_0=\lvert I_0\rvert \) for
their number. A **multiple testing procedure** maps the data to a set
\( \mathcal R\subseteq\{1,\dots,m\} \) of rejected hypotheses. The outcome can be summarized
by three counts:

| | not rejected | rejected | total |
|---|---|---|---|
| true hypotheses | \( m_0-V \) | \( V \) | \( m_0 \) |
| false hypotheses | \( m-m_0-S \) | \( S \) | \( m-m_0 \) |
| total | \( m-R \) | \( R \) | \( m \) |

Here \( R=\lvert\mathcal R\rvert \) is the number of rejections, \( V=\lvert\mathcal R\cap I_0\rvert \)
the number of **false rejections** (Type I errors), and \( S=R-V \) the number of correct ones.
The counts \( R \), \( V \) and \( S \) are random, and \( m_0 \) is unknown. The analyst sees
\( R \) but not how it splits into \( V \) and \( S \).

With one true hypothesis, a level-\( \alpha \) test makes \( \Pr(V=1)\le\alpha \). For a family there are several reasonable
ways to summarize the distribution of \( V \), and they lead to different procedures.

::: {#def-mc-error-rates}
[Error rates of a family]

For a multiple testing procedure and a given value of the parameters:

::: {.enumerate options="label=(\alph*)"}
1. the **familywise error rate** is \( \text{FWER}=\Pr(V\ge1) \), the probability of at
   least one false rejection;

2. the **per-family error rate** is \( \text{PFER}=\E(V) \), the expected number of false
   rejections, and the **per-comparison error rate** is \( \text{PCER}=\E(V)/m \);

3. the **false discovery proportion** is \( Q=V/\max(R,1) \), which is \( V/R \) when something is
   rejected and \( 0 \) otherwise, and the **false discovery rate** is its expectation,
   \( \text{FDR}=\E(Q) \).
:::

A procedure **controls** one of these rates at level \( \alpha \) if the rate is at most
\( \alpha \) for every value of the parameters.
:::

The familywise error rate treats one false claim as a failure of the whole family; it suits
conclusions that will each be acted on. The false discovery rate accepts some errors in a long list
of discoveries and bounds their expected share. The unadjusted procedure, testing each hypothesis
at level \( \alpha \), controls only the per-comparison rate.

::: {#prp-mc-error-rate-order}
[How the error rates compare]

For every procedure and every value of the parameters,
\[
\text{PCER}\le\text{FWER}\le\text{PFER},\qquad \text{FDR}\le\text{FWER}.
\]
If all the hypotheses are true (\( m_0=m \)), then \( \text{FDR}=\text{FWER} \). Consequently a
procedure that controls the per-family error rate at \( \alpha \) controls the familywise
rate at \( \alpha \), and a procedure that controls the familywise rate controls the false
discovery rate.
:::

::: {.proof}
All four rates are expectations of functions of \( V \) and \( R \), so it is enough to compare
the functions pointwise. Since \( 0\le V\le m \), we have \( V/m\le1_{\{V\ge1\}}\le V \), which
gives the first chain after taking expectations. Since \( V\le R \), the ratio
\( V/\max(R,1) \) is at most \( 1 \), and it is \( 0 \) when \( V=0 \); so
\( V/\max(R,1)\le1_{\{V\ge1\}} \), which gives \( \text{FDR}\le\text{FWER} \). If every hypothesis is
true, every rejection is false, so \( V=R \) and \( V/\max(R,1)=1_{\{V\ge1\}} \) exactly.
:::

In particular, a procedure that controls the false discovery rate is a valid level-\( \alpha \) test
of the hypothesis that *all* the \( H_j \) are true, but not in general a procedure with familywise
control ([Section 13.5](05-fdr.html)).

## How often something is significant by chance

If the \( m \) tests are independent and each has exact level \( \alpha \), then when all
hypotheses are true,
\[
\text{FWER}=1-\Pr(\text{no rejection})=1-(1-\alpha)^m .
\]{#eq-mc-independent-rate}

With \( \alpha=0.05 \) this is \( 0.401 \) for ten tests and \( 0.994 \) for a hundred. In a
linear model the tests are rarely independent. They share the estimate \( s \) of \( \sigma \), and
comparisons involving a common group share that group's mean. The next example measures the
effect of that dependence for the most common family of all.

::: {#exm-mc-chance}
[All pairs of equal means]

Take \( g \) groups of \( 5 \) observations each, with all group means equal, and compare every pair
with an unadjusted two-sided \( t \) test at level \( 0.05 \), using the pooled variance estimate on
\( g(5-1) \) degrees of freedom. The listing simulates \( 20000 \) data sets for each \( g \).
For \( g=3,5,7,10 \) (so \( 3,10,21,45 \) pairs), the proportion of data sets with at least one
significant pair is
\[
0.116,\qquad 0.265,\qquad 0.412,\qquad 0.589 .
\]
Independent tests would give \( 0.143 \), \( 0.401 \), \( 0.659 \) and \( 0.901 \). The average
number of significant pairs is \( 0.15 \), \( 0.50 \), \( 1.08 \) and \( 2.23 \), in each case
\( 0.05 \) times the number of pairs up to simulation error.
[Figure 13.1.1](#fig-mc-chance) plots the familywise rate against the number of tests.
:::

::: {when-format="html"}
![**Figure 13.1.1.** Probability of at least one false rejection among \( m \) unadjusted
level-\( 0.05 \) tests: all pairs of \( g=2,\dots,12 \) equal means (dots, simulated), independent tests
(curve, @eq-mc-independent-rate) and the Bonferroni bound \( m\alpha \) (dashed).](chance_significance.svg){#fig-mc-chance width=80%}
:::

::: {when-format="pdf"}
![Probability of at least one false rejection among \( m \) unadjusted
level-\( 0.05 \) tests: all pairs of \( g=2,\dots,12 \) equal means (dots, simulated), independent tests
(curve, @eq-mc-independent-rate) and the Bonferroni bound \( m\alpha \) (dashed).](chance_significance.pdf){width=80%}
:::

```{.python .run #cell-chance-significance-simulate}
import itertools
import numpy as np
from scipy import stats

alpha, reps, m_rep = 0.05, 20_000, 5                  # 5 observations per group

def any_significant(g, m_rep, reps, alpha):
    """Share of data sets with at least one unadjusted pairwise rejection, and mean count."""
    rng = np.random.default_rng([20260919, g])        # one fixed stream per g
    nu = g * (m_rep - 1)                              # error degrees of freedom
    y = rng.normal(size=(reps, g, m_rep))             # all group means equal
    means = y.mean(axis=2)
    s2 = ((y - means[:, :, None]) ** 2).sum(axis=(1, 2)) / nu
    se = np.sqrt(s2 * 2 / m_rep)                      # standard error of a difference
    crit = stats.t.ppf(1 - alpha / 2, nu)
    count = np.zeros(reps)
    for k, l in itertools.combinations(range(g), 2):
        count += np.abs(means[:, k] - means[:, l]) / se > crit
    return np.mean(count > 0), np.mean(count)

for g in (3, 5, 10):
    pairs = g * (g - 1) // 2
    fwer, mean_count = any_significant(g, m_rep, reps, alpha)
    print(f"g = {g:2d}, {pairs:2d} pairs: P(some rejection) = {fwer:.3f}, "
          f"mean number = {mean_count:.2f} (alpha x pairs = {alpha * pairs:.2f}), "
          f"independent tests would give {1 - (1 - alpha) ** pairs:.3f}")
```

The expected number of false rejections is additive whatever the dependence, which is why the
Bonferroni bound of [Section 13.4](04-bonferroni-holm.html) needs no assumptions. The probability of
at least one grows more slowly than for independent tests, because false rejections cluster when
\( s \) is small or one mean is extreme. Still, with ten groups some pair is “significant” in well over
half of all experiments in which nothing differs.

## Simultaneous intervals and the tests they induce

Let \( \psi_j \), \( j\in J \), be real parameters, where the index set \( J \) may be infinite. Intervals
\( [L_j,U_j] \) computed from the data are **simultaneous confidence intervals** at level \( 1-\alpha \) if
\[
\Pr\bigl(L_j\le\psi_j\le U_j\ \text{for every }j\in J\bigr)\ge1-\alpha
\]
for every value of the parameters. Scheffé's and Tukey's methods and the Bonferroni intervals all
produce them, and they yield tests with the strongest guarantees.

::: {#prp-mc-intervals-to-tests}
[Tests from simultaneous intervals]

Let \( [L_j,U_j] \), \( j\in J \), be simultaneous confidence intervals at level \( 1-\alpha \),
and for given constants \( d_j \) reject \( H_j:\psi_j=d_j \) exactly when
\( d_j\notin[L_j,U_j] \). Then:

::: {.enumerate options="label=(\alph*)"}
1. the probability of rejecting at least one true \( H_j \) is at most \( \alpha \), for every value
   of the parameters;

2. if, whenever \( H_j \) is rejected, one also declares \( \psi_j>d_j \) or \( \psi_j<d_j \)
   according to the side of \( d_j \) on which the interval lies, the probability of at least
   one false rejection *or* wrong declaration of sign is at most \( \alpha \).
:::

:::

::: {.proof}
*(a)* If a true \( H_j \) is rejected, then \( \psi_j=d_j\notin[L_j,U_j] \), so some interval fails
to cover. That event has probability at most \( \alpha \). *(b)* A wrong declaration, say
\( \psi_j>d_j \) when in fact \( \psi_j\le d_j \), happens only if \( L_j>d_j\ge\psi_j \); again
an interval fails to cover. Both kinds of error are contained in the event that some
interval misses, whose probability is at most \( \alpha \).
:::

Part (b) matters because few parameters are exactly equal to their hypothesized values; the real
question is which way each difference goes. Some good procedures, such as Holm's, give tests but no
simple companion intervals.

## Weak control is not enough

The error rates depend on which hypotheses are true. The strongest requirement bounds the rate in
every configuration; a much weaker one only when every hypothesis is true.

::: {#def-mc-strong-control}
[Strong and weak control]

A procedure controls the familywise error rate **strongly** at level \( \alpha \) if
\( \Pr(V\ge1)\le\alpha \) for every value of the parameters, and **weakly** if
\( \Pr(V\ge1)\le\alpha \) whenever all of \( H_1,\dots,H_m \) are true.
:::

Weak control is easy to obtain. Let \( H_0=\bigcap_jH_j \) be the **global null hypothesis**
that every \( H_j \) holds.

::: {#prp-mc-gatekeeper}
[A gate gives weak control]

Let \( \phi_0 \) be a level-\( \alpha \) test of the global null \( H_0 \). A procedure that rejects no
\( H_j \) unless \( \phi_0 \) rejects \( H_0 \) controls the familywise error rate weakly at level
\( \alpha \).
:::

::: {.proof}
When all \( H_j \) are true, \( V\ge1 \) requires a rejection by \( \phi_0 \), which has probability
at most \( \alpha \) because \( H_0 \) is true.
:::

The best-known gated procedure is Fisher's **protected least significant difference** (LSD): in a
one-way layout, test equality of all \( g \) means with the \( F \) test at level \( \alpha \), and only if it
rejects, compare each pair with an unadjusted \( t \) test at level \( \alpha \). Weak control is all it has.

::: {#exm-mc-protected-lsd}
[The protected LSD without strong control]

With \( g=10 \) groups of \( 5 \) observations and all means equal, the simulation below gives a
familywise error rate of \( 0.049 \), as @prp-mc-gatekeeper promises. Now move one mean far
from the other nine, which remain equal. The \( F \) test rejects in essentially every data set,
the gate is always open, and the \( 36 \) pairwise comparisons among the nine equal means are
made without adjustment. The familywise error rate rises to \( 0.539 \).
:::

```{.python .run #cell-lsd-protected}
import itertools
import numpy as np
from scipy import stats

alpha, reps, m_rep = 0.05, 20_000, 5

def protected_lsd_fwer(mu, m_rep, reps, alpha, seed):
    """Share of data sets in which the protected LSD rejects some true equality."""
    rng = np.random.default_rng(seed)
    g = len(mu)
    nu = g * (m_rep - 1)
    y = np.asarray(mu)[None, :, None] + rng.normal(size=(reps, g, m_rep))
    means = y.mean(axis=2)
    s2 = ((y - means[:, :, None]) ** 2).sum(axis=(1, 2)) / nu
    between = m_rep * ((means - means.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    gate = between / (g - 1) / s2 > stats.f.ppf(1 - alpha, g - 1, nu)   # the F test
    crit = stats.t.ppf(1 - alpha / 2, nu) * np.sqrt(2 * s2 / m_rep)
    false = np.zeros(reps, dtype=bool)
    for k, l in itertools.combinations(range(g), 2):
        if mu[k] == mu[l]:                            # a true null hypothesis
            false |= gate & (np.abs(means[:, k] - means[:, l]) > crit)
    return false.mean()

all_equal = protected_lsd_fwer([0.0] * 10, m_rep, reps, alpha, seed=1)
one_apart = protected_lsd_fwer([0.0] * 9 + [10.0], m_rep, reps, alpha, seed=2)
print(f"10 groups, all means equal:        FWER = {all_equal:.3f}")
print(f"10 groups, one mean far from rest: FWER = {one_apart:.3f}")
```

Weak control protects only against the one configuration in which there is nothing to find. Every
procedure recommended in this chapter controls its error rate strongly.

## Exercises

### A. Check your understanding

::: {#exr-mc-independent-levels}
[A1]

Ten and twenty independent tests of true hypotheses are made, each at level \( 0.01 \). Find the
familywise error rate in each case. What common level for twenty independent tests gives a
familywise rate of exactly \( 0.05 \)?
:::

::: {.solution}
By @eq-mc-independent-rate the rates are \( 1-0.99^{10}=0.0956 \) and
\( 1-0.99^{20}=0.1821 \). The level \( a \) with \( 1-(1-a)^{20}=0.05 \) is \( a=1-0.95^{1/20}=0.00256 \),
slightly more than \( 0.05/20=0.00250 \). This is the Šidák level of [Section 13.4](04-bonferroni-holm.html).
:::

::: {#exr-mc-fdr-small}
[A2]

Give a configuration of true and false hypotheses, and a procedure, for which the familywise error
rate is close to one while the false discovery rate is close to zero.
:::

### B. Practice

::: {#exr-mc-lsd-three}
[B1]

Show that with \( g=3 \) groups the protected LSD controls the familywise error rate strongly.
*Hint:* consider separately the configurations with three, one and no true equalities among the
three pairwise hypotheses.
:::

::: {.solution}
The three hypotheses are \( \mu_1=\mu_2 \), \( \mu_1=\mu_3 \), \( \mu_2=\mu_3 \). If two of them are
true, the third is true as well, so the possible numbers of true hypotheses are \( 3 \), \( 1 \)
and \( 0 \). With three, all means are equal and weak control (@prp-mc-gatekeeper) applies.
With none, no false rejection is possible. With exactly one, say \( \mu_1=\mu_2\ne\mu_3 \), a false
rejection needs the unadjusted \( t \) test of \( \mu_1=\mu_2 \) to reject, an event of
probability \( \alpha \) whether or not the gate opens. So the familywise rate is at most
\( \alpha \) in every configuration. With four or more groups, a configuration can hold two
disjoint true equalities, and the argument fails.
:::

::: {#exr-mc-unadjusted-fdr}
[B2]

Let \( m_0 \) true hypotheses have independent uniform \( p \)-values, and let \( m_1 \) false
hypotheses have \( p \)-values equal to \( 0 \). The unadjusted procedure rejects when
\( p_j\le\alpha \). Find its familywise error rate, show that its false discovery rate is at most
\( m_0\alpha/m_1 \), and conclude that the false discovery rate tends to \( 0 \) as \( m_1\to\infty \)
with \( m_0 \) fixed.
:::

### C. Going deeper

::: {#exr-mc-lsd-limit}
[C1]

In a one-way layout with \( g\ge3 \) groups of \( m \) observations and known \( \sigma \), run the protected
LSD with \( z \) tests: the gate is the level-\( \alpha \) test that rejects when
\( \sum_km(\bar{y}_k-\bar{y})^2/\sigma^2 \) exceeds the upper \( \alpha \) point of \( \chi^2(g-1) \), and a pair is
declared different when \( \lvert\bar{y}_k-\bar{y}_l\rvert>z_{\alpha/2}\,\sigma\sqrt{2/m} \), where \( z_{\alpha/2} \) is the
upper \( \alpha/2 \) point of \( \Normal(0,1) \). Let
\( \mu_1=\dots=\mu_{g-1}=0 \) and \( \mu_g=\Delta \). Show that as \( \Delta\to\infty \) the familywise error rate
tends to \( \Pr(W_{g-1}>\sqrt2\,z_{\alpha/2}) \), where \( W_{g-1} \) is the range of \( g-1 \) independent
\( \Normal(0,1) \) variables. Show that this limit is at least \( 1-(1-\alpha)^{\lfloor(g-1)/2\rfloor} \) and
that it tends to \( 1 \) as \( g\to\infty \).
:::

::: {.solution}
With \( Z_k=\sqrt m\,\bar{y}_k/\sigma \), \( k<g \), independent \( \Normal(0,1) \) and free of \( \Delta \), some pair
\( k,l<g \) is declared different iff \( W_{g-1}=\max_{k<g}Z_k-\min_{k<g}Z_k>\sqrt2\,z_{\alpha/2} \), so the rate
is \( \Pr(G\cap\{W_{g-1}>\sqrt2\,z_{\alpha/2}\}) \), with \( G \) the event that the gate opens. The gate
statistic is at least \( m(\bar{y}_g-\bar{y})^2/\sigma^2 \), where \( \bar{y}_g-\bar{y} \) is \( \tfrac{g-1}{g}\Delta \) plus a
variable free of \( \Delta \); so \( \Pr(G)\to1 \), giving the limit. The disjoint pairs
\( (1,2),(3,4),\dots \) give \( \lfloor(g-1)/2\rfloor \) independent statistics
\( \lvert Z_{2i-1}-Z_{2i}\rvert/\sqrt2 \), each exceeding \( z_{\alpha/2} \) with probability \( \alpha \), and
\( W_{g-1} \) is at least \( \sqrt2 \) times each of them; hence the lower bound, which tends to \( 1 \).
:::
