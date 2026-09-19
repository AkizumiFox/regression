# Choosing a procedure

None of the procedures of this chapter is best in general, because they answer different
questions. This section compares them on one data set, collects the comparison into a guide, and
explains why several procedures common in software are not developed here.

## Decide the family first

Every guarantee in this chapter is for a *family* of inferences, and is only as good as the honesty
with which the family is defined. A short list of
functions written down before the data arrive is a **planned** family. All \( g(g-1)/2 \) pairwise
differences form a family with its own exact method. A function chosen *because* the data make it
look interesting belongs to the family of everything the analyst might have chosen, usually a
whole subspace; this is **post hoc** inference, or data snooping. And a long list screened for leads
worth following up calls for a different error rate altogether.

::: {.warning}
The family must be fixed before the data are examined, and reported. Choosing the adjustment after
seeing the results, for example reporting Bonferroni's intervals for the three comparisons that
turned out largest, reproduces the problem that adjustment was meant to solve. So does splitting
one analysis into several “families” so that each is small.
:::

## One data set, five sets of intervals

The education example of @exm-mc-scheffe-education and @exm-mc-tukey-education makes the
differences concrete. The table gives, for each method, the multiplier of the standard error for
all \( 21 \) pairwise differences, the resulting half-width for two pairs with very different
precision, and the number of pairs declared different at familywise level \( 0.05 \) (false
discovery rate \( 0.05 \) for BH). The pair of levels 1 and 2 involves the two smallest groups, with
standard error \( 0.420 \); the pair of levels 3 and 6 involves two large groups, with standard error
\( 0.124 \).

| method | multiplier | ratio to \( t \) | half-width, 1 vs 2 | half-width, 3 vs 6 | pairs declared different |
|---|---|---|---|---|---|
| unadjusted \( t \) | \( 1.962 \) | \( 1.00 \) | \( 0.824 \) | \( 0.244 \) | \( 14 \) |
| Tukey–Kramer | \( 2.955 \) | \( 1.51 \) | \( 1.241 \) | \( 0.368 \) | \( 9 \) |
| Šidák | \( 3.039 \) | \( 1.55 \) | \( 1.276 \) | \( 0.378 \) | \( 8 \) |
| Bonferroni | \( 3.046 \) | \( 1.55 \) | \( 1.279 \) | \( 0.379 \) | \( 8 \) |
| Scheffé | \( 3.557 \) | \( 1.81 \) | \( 1.494 \) | \( 0.442 \) | \( 6 \) |
| Holm | | | | | \( 10 \) |
| Benjamini–Hochberg | | | | | \( 13 \) |

The ordering is the one the theory predicts: the method designed for pairs is shortest, the general
finite-family methods are a few percent longer, and the method covering every contrast is much
longer. Holm's procedure rejects more by stepping down, and BH by controlling a weaker error rate.
The Šidák row relies on Šidák's inequality for correlated normal variables, stated without proof in
[Section 13.4](04-bonferroni-holm.html).

```{.python .run #cell-education-table}
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
c_sch = np.sqrt((g - 1) * stats.f.ppf(1 - alpha, g - 1, nu))   # Scheffe multiplier, q = g - 1
pairs = list(itertools.combinations(range(g), 2))
m = len(pairs)                                        # 21 comparisons
diff = np.array([means[k] - means[l] for k, l in pairs])
se_d = np.array([s * np.sqrt(1 / n_k[k] + 1 / n_k[l]) for k, l in pairs])
t = diff / se_d
mult = {                                              # half-width = multiplier x se
    "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
    "Bonferroni": stats.t.ppf(1 - alpha / (2 * m), nu),
    "Sidak": stats.t.ppf((1 + (1 - alpha) ** (1 / m)) / 2, nu),
    "Tukey-Kramer": stats.studentized_range.ppf(1 - alpha, g, nu) / np.sqrt(2),
    "Scheffe": c_sch,
}
for name, c in mult.items():
    print(f"{name:13s} multiplier {c:.3f}   pairs declared different: {np.sum(np.abs(t) > c)}")
```

## A decision guide

The choice can be put as a short sequence of questions.

::: {.enumerate options="label=(\arabic*)"}
1. *Is there one question?* Then use the ordinary \( t \) or \( F \) procedure of
   [Chapter 11](../ch11-general-linear-hypothesis/index.html) or
   [Chapter 12](../ch12-intervals-and-bands/index.html). No adjustment is needed.

2. *Will the functions be chosen after seeing the data?* If they lie in a known subspace, use
   Scheffé's intervals for that subspace, choosing the smallest subspace that contains every function
   one might pick ([Section 13.2](02-scheffe.html)).

3. *Is the family all pairs of means?* Use Tukey's intervals if the design is balanced and the
   means are independent, Tukey–Kramer if the one-way layout is unbalanced, and the exact
   multivariate \( t \) computation otherwise.

4. *Is the family a short list fixed in advance?* Use Bonferroni's (or Šidák's) intervals, or Holm's
   tests if intervals are not needed, or the exact multivariate \( t \) value if software is at hand.

5. *Is the family long, and the goal a list of leads for further study?* Control the false
   discovery rate with the Benjamini–Hochberg procedure, or with the Benjamini–Yekutieli version
   if dependence among the tests is strong and of unknown form.
:::

Where intervals are available they are preferable to tests alone, because they report the size
of each difference as well as its sign (@prp-mc-intervals-to-tests). Adjusted \( p \)-values, @eq-mc-holm-adjusted
and @eq-mc-bh-adjusted, let readers apply a different level. And a report
should say how many comparisons were made in all, not only which were significant.

## Planned versus post hoc

A post hoc comparison is chosen because it looks large, and its \( t \) statistic has the
distribution of a maximum. In the education example the most extreme contrast has \( t=6.88 \). Under
equal means its square is not an \( F(1,937) \) variable but \( 6 \) times an \( F(6,937) \) variable (@cor-mc-scheffe-f),
and treating it as planned would overstate the evidence by exactly the amount
Scheffé's multiplier corrects. The same happens at a larger scale when regressors are selected by
their \( t \) statistics; Chapter 29 takes up inference after model selection.

## Procedures not developed in this book

Several procedures widely available in software are deliberately left out. They control the
familywise error rate only weakly, or not at all.

**Fisher's protected LSD** controls the familywise error rate only weakly, as
@exm-mc-protected-lsd and @exr-mc-lsd-limit show, except with three groups (@exr-mc-lsd-three).

**The Newman–Keuls procedure** orders the means and tests the range of each set of \( w \) adjacent
means against \( q_\alpha(w,\nu)\,s/\sqrt m \), stepping down from \( w=g \) and never testing inside a set whose
range was not significant. Under the global null it rejects exactly when Tukey's test does, so its
control is weak; with means in well-separated groups, each small group is tested nearly unprotected.

::: {#exm-mc-newman-keuls}
[Newman–Keuls without strong control]

Take \( g=10 \) groups of \( 5 \) observations, with the true means equal in pairs and the five pairs far
apart, and \( \nu=40 \). The larger ranges are always significant, so the procedure reaches every pair and
tests it with \( q_\alpha(2,\nu) \), the critical value of a single unadjusted comparison. In \( 10000 \)
simulated data sets, the Newman–Keuls procedure declares some pair of equal means different in a
proportion \( 0.219 \) of them. Five independent level-\( 0.05 \) tests would err with probability
\( 0.226 \); the shared \( s \) lowers the rate slightly, as @prp-mc-sidak implies. Tukey's intervals err
in \( 0.011 \) of the data sets.
:::

```{.python .run #cell-newman-keuls-newman-keuls}
import numpy as np
from scipy import stats

alpha, g, m_rep, reps = 0.05, 10, 5, 10_000
nu = g * (m_rep - 1)
mu = np.repeat(np.arange(5) * 20.0, 2)                # pairs 0,0, 20,20, 40,40, ...

def newman_keuls(means, s, m_rep, q_crit):
    """Pairs (k, l) declared different by the Newman-Keuls step-down procedure;
    q_crit[w] is the upper-alpha point of the studentized range of w means."""
    g = len(means)
    order = np.argsort(means)
    ybar = means[order]
    se = s / np.sqrt(m_rep)
    kept = []                                         # windows found homogeneous
    different = set()
    for width in range(g, 1, -1):                     # ranges of width, width-1, ..., 2
        crit = q_crit[width] * se
        for i in range(g - width + 1):
            j = i + width - 1
            if any(a <= i and j <= b for a, b in kept):
                continue                              # inside a homogeneous window
            if ybar[j] - ybar[i] > crit:
                different.add((min(order[i], order[j]), max(order[i], order[j])))
            else:
                kept.append((i, j))
    return different

q_crit = {w: stats.studentized_range.ppf(1 - alpha, w, nu) for w in range(2, g + 1)}
rng = np.random.default_rng(20260919)
q_tukey = stats.studentized_range.ppf(1 - alpha, g, nu)
nk_false = tukey_false = 0
for r in range(reps):
    y = mu[:, None] + rng.normal(size=(g, m_rep))
    means = y.mean(axis=1)
    s = np.sqrt(np.sum((y - means[:, None]) ** 2) / nu)
    true_pairs = [(2 * i, 2 * i + 1) for i in range(5)]
    nk = newman_keuls(means, s, m_rep, q_crit)
    nk_false += any(pair in nk for pair in true_pairs)
    tukey_false += any(abs(means[a] - means[b]) > q_tukey * s / np.sqrt(m_rep)
                       for a, b in true_pairs)
print(f"P(some false difference): Newman-Keuls {nk_false / reps:.3f}, "
      f"Tukey {tukey_false / reps:.3f}")
```

With \( 2J \) means in \( J \) separated pairs, the familywise error rate of Newman–Keuls approaches
\( 1-(1-\alpha)^J \) as the degrees of freedom grow (@exr-mc-nk-limit). No bound below one holds
uniformly in \( g \).

**Duncan's multiple range test.** Duncan's procedure has the structure of Newman–Keuls but raises
the level used for a set of \( w \) means to \( 1-(1-\alpha)^{w-1} \). Its first step, the range of all \( g \)
means, is therefore a test at level \( 1-(1-\alpha)^{g-1} \), and if it rejects, at least the most
extreme pair is declared different. Under the global null its familywise error rate is exactly
\( 1-(1-\alpha)^{g-1} \), which is \( 0.370 \) for \( g=10 \) and \( \alpha=0.05 \). It does not control the rate
even weakly.

These procedures buy power over Tukey's by relaxing the error rate when some hypotheses are
false. If that trade is wanted, it is better made openly, by controlling the false discovery rate,
whose meaning is clear. Closed testing (@thm-mc-closed-testing) with range or \( F \) local tests gives
step-down procedures that are more powerful than Tukey's and still control the familywise error rate
strongly; the books cited in the notes describe them.

A Bayesian response to multiplicity is a hierarchical model that treats the group means as draws
from a common distribution and shrinks the estimated differences toward each other (@exr-mc-bayes-shrink). Chapter 32 develops random-effects models
of this kind, and the conjugate analysis of @thm-opt-bayes-conjugate is their simplest case.

## Exercises

### A. Check your understanding

::: {#exr-mc-choose}
[A1]

Choose a procedure for each situation, and say what family its guarantee refers to:
*(a)* a trial with three doses and a placebo, where the protocol specifies the three comparisons
with placebo; *(b)* eight fertilizers compared in a balanced field trial, with interest in every
pair; *(c)* \( 5000 \) genes screened for association with a trait, to select candidates for a second
study; *(d)* a contrast suggested by a plot of the group means.
:::

### B. Practice

::: {#exr-mc-duncan}
[B1]

Show that under the global null the familywise error rate of Duncan's procedure is
\( 1-(1-\alpha)^{g-1} \) in a balanced one-way layout.
:::

::: {#exr-mc-nk-limit}
[B2]

In @exm-mc-newman-keuls with \( J \) pairs, let the separation between pairs tend to infinity. Show
that the familywise error rate of Newman–Keuls converges to
\( 1-\Pr\bigl(\lvert T_j\rvert\le q_\alpha(2,\nu)/\sqrt2,\ j=1,\dots,J\bigr) \), where the \( T_j \) are \( t \)
statistics with independent numerators and a common \( s \). Using @prp-mc-sidak, show that the limit is
at most \( 1-(1-\alpha)^J \), and that it tends to this value as \( \nu\to\infty \).
:::

### C. Going deeper

::: {#exr-mc-bayes-shrink}
[C1]

In the balanced one-way layout with known \( \sigma^2 \), give the group means the prior
\( \mu_k\iid\Normal(\mu_0,\tau^2) \), with \( \mu_0 \) and \( \tau^2 \) known. Find the posterior distribution of
\( \mu_k-\mu_l \) and compare the width and centre of its central \( 95\% \) posterior interval with the
unadjusted confidence interval. Explain in what sense the shrinkage addresses the multiplicity of
comparisons, and in what sense it does not.
:::
