# Tukey's studentized range

For the family of all pairwise differences \( \mu_k-\mu_l \), Tukey saw that the largest standardized
error is the range of the estimated means divided by an estimate of their common standard
deviation. When the means are estimated independently and with equal precision, as in a balanced
one-way layout, this ratio is pivotal, and its upper \( \alpha \) point gives intervals for all pairs
with simultaneous coverage exactly \( 1-\alpha \), shorter than Scheffé's or Bonferroni's.

## The studentized range distribution

::: {#def-mc-studentized-range}
[Studentized range]

Let \( Z_1,\dots,Z_g \) be independent \( \Normal(0,1) \) variables, \( g\ge2 \), and let
\( V\sim\chi^2(\nu) \) be independent of them. The law of
\[
Q=\frac{\max_kZ_k-\min_kZ_k}{\sqrt{V/\nu}}
\]
is the **studentized range distribution** with parameters \( g \) and \( \nu \). Its upper
\( \alpha \) point is written \( q_\alpha(g,\nu) \).
:::

The numerator is the range of a normal sample of size \( g \). Dividing by \( \sqrt{V/\nu} \)
replaces the known standard deviation \( 1 \) by an independent estimate of it with \( \nu \)
degrees of freedom, just as a \( t \) variable is a studentized normal. For \( g=2 \) the range is
\( \lvert Z_1-Z_2\rvert \), a \( \Normal(0,2) \) variable in absolute value, so
\( Q/\sqrt2 \) is the absolute value of a \( t(\nu) \) variable and
\[
q_\alpha(2,\nu)=\sqrt2\,t_{\nu,\alpha/2}.
\]{#eq-mc-q-two}

For larger \( g \) the distribution function is a double integral: conditioning on the smallest
\( Z_k \) gives \( \Pr(\max_kZ_k-\min_kZ_k\le w)=g\int\varphi(z)[\Phi(z+w)-\Phi(z)]^{g-1}dz \), and one further
integration over \( V \) gives the distribution of \( Q \). SciPy evaluates it as
`scipy.stats.studentized_range`. For \( g=5 \) and \( \nu=20 \) it gives \( q_{0.05}(5,20)=4.232 \), and
\( 200000 \) simulated values of \( Q \) put the upper \( 5\% \) point at \( 4.230 \).

```{.python .run #cell-tukey-sim-quantile}
import itertools
import numpy as np
from scipy import stats

alpha, reps = 0.05, 40_000

g, nu = 5, 20
rng = np.random.default_rng(20260919)
Z = rng.normal(size=(200_000, g))
V = rng.chisquare(nu, size=200_000)
Q = (Z.max(axis=1) - Z.min(axis=1)) / np.sqrt(V / nu)   # studentized range
q_sim = np.quantile(Q, 1 - alpha)
q_exact = stats.studentized_range.ppf(1 - alpha, g, nu)
print(f"upper 5% point of the studentized range, g = {g}, nu = {nu}:"
      f" simulated {q_sim:.3f}, scipy {q_exact:.3f}")
```

## Tukey's intervals for a balanced layout

Consider the one-way model with \( g \) levels and the same number \( m \) of observations at every
level, so \( n=gm \), and assume normal errors. The level means \( \bar{y}_1,\dots,\bar{y}_g \) are the least
squares estimates of \( \mu_1,\dots,\mu_g \), the pooled variance is \( s^2=\text{SSE}/\nu \) with
\( \nu=g(m-1) \), and the standard error of every difference \( \bar{y}_k-\bar{y}_l \) is
\( s\sqrt{2/m} \).

::: {#lem-mc-range-contrast}
[A contrast is bounded by the range]

If \( \sum_kc_k=0 \), then for all real \( x_1,\dots,x_g \),
\[
\Bigl\lvert\sum_kc_kx_k\Bigr\rvert\le\frac12\sum_k\lvert c_k\rvert\,\bigl(\max_kx_k-\min_kx_k\bigr),
\]
with equality when \( \mathbf{c} \) is a multiple of a pairwise difference \( \mathbf{e}_k-\mathbf{e}_l \) and \( x_k \), \( x_l \) are the
largest and smallest of the \( x \)'s.
:::

::: {.proof}
Let \( a=\tfrac12(\max_kx_k+\min_kx_k) \). Because \( \sum_kc_k=0 \),
\( \sum_kc_kx_k=\sum_kc_k(x_k-a) \), and \( \lvert x_k-a\rvert\le\tfrac12(\max_kx_k-\min_kx_k) \) for every \( k \).
The triangle inequality gives the bound. For \( \mathbf{c}=\mathbf{e}_k-\mathbf{e}_l \), the left side is
\( \lvert x_k-x_l\rvert \) and the right side is the range, which agree when \( x_k,x_l \) are the extremes.
:::

::: {#thm-mc-tukey}
[Tukey's simultaneous intervals]

In the balanced one-way model with normal errors, \( g\ge2 \) and \( \nu=g(m-1)\ge1 \), let
\( c_T=q_\alpha(g,\nu) \). Then:

::: {.enumerate options="label=(\alph*)"}
1. with probability exactly \( 1-\alpha \),
   \[
   \mu_k-\mu_l\in\Bigl[\bar{y}_k-\bar{y}_l-c_T\frac{s}{\sqrt m},\ \ \bar{y}_k-\bar{y}_l+c_T\frac{s}{\sqrt m}\Bigr]
   \quad\text{for all }k\ne l ;
   \]{#eq-mc-tukey}

2. on the same event, for every contrast \( \psi=\sum_kc_k\mu_k \),
   \[
   \Bigl\lvert\sum_kc_k\bar{y}_k-\psi\Bigr\rvert\le c_T\frac{s}{\sqrt m}\cdot\frac12\sum_k\lvert c_k\rvert ,
   \]
   so these intervals for all contrasts also have simultaneous coverage exactly \( 1-\alpha \).
:::

:::

::: {.proof}
By @thm-opt-sampling, the means \( \bar{y}_k \) are independent \( \Normal(\mu_k,\sigma^2/m) \), and they are
independent of \( s^2 \), with \( \nu s^2/\sigma^2\sim\chi^2(\nu) \). Put
\( Z_k=\sqrt m(\bar{y}_k-\mu_k)/\sigma \) and \( V=\nu s^2/\sigma^2 \). These satisfy the assumptions of @def-mc-studentized-range,
whatever the \( \mu_k \) are. Now
\[
\max_{k\ne l}\bigl\lvert(\bar{y}_k-\bar{y}_l)-(\mu_k-\mu_l)\bigr\rvert
=\max_{k,l}\bigl[(\bar{y}_k-\mu_k)-(\bar{y}_l-\mu_l)\bigr]
=\frac{\sigma}{\sqrt m}\bigl(\max_kZ_k-\min_kZ_k\bigr),
\]
so all intervals @eq-mc-tukey cover iff \( (\max_kZ_k-\min_kZ_k)/\sqrt{V/\nu}\le c_T \), an event of
probability \( 1-\alpha \), since the studentized range has a continuous distribution. For (b), apply @lem-mc-range-contrast
with \( x_k=\bar{y}_k-\mu_k \): on that event the range of the \( x_k \) is at most
\( c_Ts/\sqrt m \). The pairwise differences belong to the family of contrasts, so the coverage of the
enlarged family cannot exceed that of the pairs, and it equals \( 1-\alpha \).
:::

::: {#prp-mc-tukey-shortest}
[Tukey's intervals are the shortest for all pairs]

In the balanced one-way model, suppose intervals \( \bar{y}_k-\bar{y}_l\pm c\,s\sqrt{2/m} \) with a common
multiplier \( c \) have simultaneous coverage at least \( 1-\alpha \) for all pairs. Then
\( c\ge q_\alpha(g,\nu)/\sqrt2 \). In particular the Bonferroni, Šidák and Scheffé intervals for all
pairs are at least as long as Tukey's.
:::

::: {.proof}
By the proof of @thm-mc-tukey, the simultaneous coverage of such intervals is
\( \Pr\bigl(Q\le c\sqrt2\bigr) \), where \( Q \) has the studentized range distribution. This is at least
\( 1-\alpha \) only if \( c\sqrt2\ge q_\alpha(g,\nu) \). The other three methods give intervals of this form with
coverage at least \( 1-\alpha \) (@thm-mc-scheffe and [Section 13.4](04-bonferroni-holm.html)).
:::

[Figure 13.3.1](#fig-mc-multipliers) compares the multipliers of the standard error of a
difference for \( \nu=30 \). With \( g=5 \) groups they are \( 2.901 \) for Tukey, \( 3.021 \) for Šidák,
\( 3.030 \) for Bonferroni and \( 3.280 \) for Scheffé. With \( g=10 \) they are \( 3.411 \), \( 3.597 \),
\( 3.607 \) and \( 4.461 \). The unadjusted multiplier is \( 2.042 \) throughout. For \( g=2 \) all methods
coincide, by @eq-mc-q-two. Scheffé's method falls further behind as \( g \) grows, because it covers
a space of dimension \( g-1 \) while the family of pairs is finite.

::: {when-format="html"}
![**Figure 13.3.1.** Multipliers of the standard error for \( 95\% \) simultaneous
intervals for all pairwise differences of \( g \) balanced means, \( \nu=30 \) error degrees of
freedom. Tukey's multiplier is \( q_{0.05}(g,30)/\sqrt2 \). The Šidák and Bonferroni curves nearly
coincide.](multipliers.svg){#fig-mc-multipliers width=80%}
:::

::: {when-format="pdf"}
![Multipliers of the standard error for \( 95\% \) simultaneous
intervals for all pairwise differences of \( g \) balanced means, \( \nu=30 \) error degrees of
freedom. Tukey's multiplier is \( q_{0.05}(g,30)/\sqrt2 \). The Šidák and Bonferroni curves nearly
coincide.](multipliers.pdf){width=80%}
:::

```{.python .run #cell-multipliers-multipliers}
import numpy as np
from scipy import stats

alpha, nu = 0.05, 30

def pairwise_multipliers(g, nu, alpha=0.05):
    m = g * (g - 1) // 2                              # number of pairs
    return {
        "unadjusted": stats.t.ppf(1 - alpha / 2, nu),
        "Tukey": stats.studentized_range.ppf(1 - alpha, g, nu) / np.sqrt(2),
        "Sidak": stats.t.ppf((1 + (1 - alpha) ** (1 / m)) / 2, nu),
        "Bonferroni": stats.t.ppf(1 - alpha / (2 * m), nu),
        "Scheffe": np.sqrt((g - 1) * stats.f.ppf(1 - alpha, g - 1, nu)),
    }

for g in (2, 3, 5, 10):
    row = pairwise_multipliers(g, nu)
    print(f"g = {g:2d}: " + "  ".join(f"{k} {v:.3f}" for k, v in row.items()))
```

For contrasts other than pairs the comparison can reverse. Part (b) of @thm-mc-tukey gives the
half-width \( q_\alpha(g,\nu)\,(s/\sqrt m)\,\tfrac12\sum_k\lvert c_k\rvert \), while Scheffé's is
\( \sqrt{(g-1)F_\alpha(g-1,\nu)}\,(s/\sqrt m)\,\bigl(\sum_kc_k^2\bigr)^{1/2} \). With \( g=6 \) and \( \nu=30 \), in
units of \( s/\sqrt m \), a pairwise difference has half-width \( 4.301 \) by Tukey's method and
\( 5.033 \) by Scheffé's; for \( c=\tfrac13(1,1,1,-1,-1,-1) \) the half-widths are \( 4.301 \) and \( 2.906 \).
Tukey's bound is sharp for pairs and crude for contrasts spread over many means.

## Unbalanced layouts: the Tukey–Kramer intervals

With unequal group sizes \( n_k \) the maximum standardized difference is no longer a studentized
range. Tukey and Kramer proposed keeping the multiplier and using the standard error of each
difference:
\[
\bar{y}_k-\bar{y}_l\;\pm\;\frac{q_\alpha(g,\nu)}{\sqrt2}\;s\sqrt{\frac1{n_k}+\frac1{n_l}},\qquad
\nu=n-g .
\]{#eq-mc-tukey-kramer}

Tukey conjectured, and Hayter (1984) proved, that the simultaneous coverage of these intervals is
never less than \( 1-\alpha \), whatever the group sizes. The proof is beyond this chapter. How
conservative the intervals are is a separate question; the simulation below estimates the coverage
of 95% intervals in three designs from \( 40000 \) data sets each (standard error about \( 0.0011 \)):

| design | coverage |
|---|---|
| balanced, five groups of 5 | \( 0.9511 \) |
| the education example (sizes 13 to 248) | \( 0.9562 \) |
| four groups of 2 and four of 40 | \( 0.9614 \) |

The balanced design is exact; the unbalanced ones overcover only modestly, even with very different
group sizes.

```{.python .run #cell-tukey-sim-coverage}
import itertools
import numpy as np
from scipy import stats

alpha, reps = 0.05, 40_000

def tukey_kramer_coverage(n_k, reps, alpha, seed):
    """Share of data sets in which every Tukey-Kramer interval covers its difference."""
    rng = np.random.default_rng(seed)
    n_k = np.asarray(n_k)
    g, nu = len(n_k), n_k.sum() - len(n_k)
    q = stats.studentized_range.ppf(1 - alpha, g, nu)
    # group means (true means 0) and an independent pooled variance, as in the model
    means = rng.normal(size=(reps, g)) / np.sqrt(n_k)
    s = np.sqrt(rng.chisquare(nu, size=reps) / nu)
    ok = np.ones(reps, dtype=bool)
    for k, l in itertools.combinations(range(g), 2):
        half = q / np.sqrt(2) * s * np.sqrt(1 / n_k[k] + 1 / n_k[l])
        ok &= np.abs(means[:, k] - means[:, l]) <= half
    return ok.mean()

balanced = tukey_kramer_coverage([5] * 5, reps, alpha, seed=1)
education = tukey_kramer_coverage([13, 52, 248, 187, 90, 227, 127], reps, alpha, seed=2)
lopsided = tukey_kramer_coverage([2, 2, 2, 2, 40, 40, 40, 40], reps, alpha, seed=3)
print(f"balanced (5 x 5):           {balanced:.4f}")
print(f"education design:           {education:.4f}")
print(f"four groups of 2, four of 40: {lopsided:.4f}")
```

::: {#exm-mc-tukey-education}
[All pairs of education levels]

In @exm-mc-scheffe-education there are \( g=7 \) education levels, so \( 21 \) pairs, and
\( \nu=937 \). The Tukey–Kramer multiplier is \( q_{0.05}(7,937)/\sqrt2=4.179/\sqrt2=2.955 \). It
declares \( 9 \) of the \( 21 \) differences nonzero, against \( 14 \) for unadjusted \( t \) intervals
([Section 13.6](06-choosing.html) compares all the methods). [Figure 13.3.2](#fig-mc-tukey-education) shows the intervals. The
significant differences all separate a lower education level from a higher one, and the three most
educated levels do not differ from one another. The smallest group, with only \( 13 \)
respondents, has the highest mean, but its intervals are so wide that only its differences from
levels 5, 6 and 7 are established.
:::

::: {when-format="html"}
![**Figure 13.3.2.** Differences in mean placement of Clinton between education levels, with
95% Tukey–Kramer simultaneous intervals (thin lines) and unadjusted \( t \) intervals (thick pale
bands). Differences whose Tukey–Kramer interval excludes zero are drawn in blue.](tukey_education.svg){#fig-mc-tukey-education width=70%}
:::

::: {when-format="pdf"}
![Differences in mean placement of Clinton between education levels, with
95% Tukey–Kramer simultaneous intervals (thin lines) and unadjusted \( t \) intervals (thick pale
bands). Differences whose Tukey–Kramer interval excludes zero are drawn in blue.](tukey_education.pdf){width=70%}
:::

```{.python .run #cell-education-pairwise}
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
pval = 2 * stats.t.sf(np.abs(t), nu)

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

## Beyond the one-way layout

The proof of @thm-mc-tukey used only that the \( g \) estimates are independent, with a common
variance \( \sigma^2c \) for a known \( c \), and independent of \( s^2 \). That also holds for the levels of
one factor in a balanced multi-way design (@exr-mc-tukey-two-way), the subject of Chapter 15 and
Chapter 16. For correlated estimates, such as regression coefficients or adjusted means, neither
form is exact; the exact critical value for any planned family then comes from a multivariate \( t \)
distribution (@exr-mc-multivariate-t), which standard software computes.

## Exercises

### A. Check your understanding

::: {#exr-mc-q-two}
[A1]

Prove @eq-mc-q-two, and deduce that Tukey's intervals for \( g=2 \) are the ordinary \( t \) intervals.
:::

### B. Practice

::: {#exr-mc-tukey-two-way}
[B1]

In a two-way additive model \( \E(y_{ijk})=\mu+\alpha_i+\beta_j \) with \( a \) rows, \( b \) columns and \( N \)
observations in every cell, show that the column means \( \bar{y}_{\cdot j\cdot} \) are independent, each with
variance \( \sigma^2/(aN) \), and independent of the residual mean square. Deduce Tukey intervals for
all differences \( \beta_j-\beta_{j'} \), and state their degrees of freedom.
:::

::: {.solution}
Each \( \bar{y}_{\cdot j\cdot} \) averages the \( aN \) observations of column \( j \), and different columns use
disjoint sets of independent observations, so the column means are independent with variance
\( \sigma^2/(aN) \). Their expectations are \( \mu+\bar{\alpha}+\beta_j \), with \( \bar{\alpha} \) the average row
effect, so differences estimate \( \beta_j-\beta_{j'} \). The column means are functions of the fitted
values (the fitted value of cell \( (i,j) \) is \( \bar{y}_{i\cdot\cdot}+\bar{y}_{\cdot j\cdot}-\bar{y} \), and averaging it
over \( i \) gives \( \bar{y}_{\cdot j\cdot} \)), so they are independent of \( s^2 \) by @thm-opt-sampling. The
residual degrees of freedom are \( \nu=abN-a-b+1 \), and the intervals are
\( \bar{y}_{\cdot j\cdot}-\bar{y}_{\cdot j'\cdot}\pm q_\alpha(b,\nu)\,s/\sqrt{aN} \), exact for all pairs of columns.
:::

::: {#exr-mc-control}
[B2]

Treatments \( 1,\dots,g-1 \) are to be compared with a control, level \( 0 \), in a balanced layout. Show
that the exact simultaneous multiplier for the \( g-1 \) differences \( \mu_k-\mu_0 \) is the upper \( \alpha \)
point of \( \max_{k\ge1}\lvert Z_k-Z_0\rvert/\sqrt{V/\nu} \), and that it is smaller than \( q_\alpha(g,\nu) \).
:::

::: {.solution}
As in the proof of @thm-mc-tukey, the largest standardized error among the differences with the
control is \( (\sigma/\sqrt m)\max_{k\ge1}\lvert Z_k-Z_0\rvert \), and the statistic
\( D=\max_{k\ge1}\lvert Z_k-Z_0\rvert/\sqrt{V/\nu} \) is pivotal. Let \( Q \) be the studentized range of
\( Z_0,\dots,Z_{g-1} \), so that \( D\le Q \). The event \( \{D\le q_\alpha(g,\nu)<Q\} \) has positive
probability: it contains an open set of values of \( (Z_0,\dots,Z_{g-1},V) \), for instance those with
\( Z_0 \) near the middle of the others and the range just above the critical value, and these variables
have a positive joint density. Hence \( \Pr(D>q_\alpha(g,\nu))=\alpha-\Pr(D\le q_\alpha(g,\nu)<Q)<\alpha \),
and since \( D \) has a continuous distribution its upper \( \alpha \) point is smaller than
\( q_\alpha(g,\nu) \). This is Dunnett's many-to-one procedure.
:::

::: {#exr-mc-range-test}
[B3]

In the balanced one-way layout, the **studentized range test** rejects \( \mu_1=\dots=\mu_g \) when
\( Q_{\text{obs}}=(\max_k\bar{y}_k-\min_k\bar{y}_k)/(s/\sqrt m)>q_\alpha(g,\nu) \). Show that it has exact level
\( \alpha \), and that it rejects iff some Tukey interval @eq-mc-tukey excludes \( 0 \). Show that the
\( F \) test rejects iff \( \sum_k(\bar{y}_k-\bar{y})^2>(g-1)F_\alpha(g-1,\nu)\,s^2/m \), and describe, for \( g=3 \),
configurations of the means that one test rejects and the other does not.
:::

::: {.solution}
Under the hypothesis \( Q_{\text{obs}} \) has the studentized range distribution (proof of @thm-mc-tukey with
all \( \mu_k \) equal), so the level is exactly \( \alpha \), and \( Q_{\text{obs}}>q_\alpha(g,\nu) \) says that the
extreme pair's interval excludes \( 0 \); every other pair has a smaller difference. The \( F \) statistic is
\( m\sum_k(\bar{y}_k-\bar{y})^2/\{(g-1)s^2\} \). For \( g=3 \) and a fixed range \( 2a \), the sum of squares
is \( 2a^2 \) for means \( (-a,0,a) \) and \( 8a^2/3 \) for \( (-a,-a,a) \). In the limit \( \nu\to\infty \), in units of
\( s/\sqrt m \), the range test rejects when \( 2a>3.314 \) and the \( F \) test when the sum of squares exceeds
\( 5.991 \): so \( (-a,0,a) \) with \( 1.657<a<1.731 \) is rejected only by the range test, and \( (-a,-a,a) \) with
\( 1.499<a<1.657 \) only by the \( F \) test. The range test is the natural gate for pairs, the \( F \) test for
contrasts in general.
:::

### C. Going deeper

::: {#exr-mc-multivariate-t}
[C1]

Let \( \blambda_1,\dots,\blambda_k\in\C(\X\T) \) be chosen in advance, with \( t \) statistics
\( T_j=(\blambda_j\T\hbeta-\blambda_j\T\bbeta)/\text{se}(\blambda_j\T\hbeta) \). Show that
\( (T_1,\dots,T_k)=\mathbf{Z}/(s/\sigma) \), where \( \mathbf{Z}\sim\Normal_k(\bzero,\R) \) is independent of \( s \), and
\( \R \) is the correlation matrix with entries
\( \blambda_i\T\G\blambda_j/\sqrt{\blambda_i\T\G\blambda_i\cdot\blambda_j\T\G\blambda_j} \). Deduce that the
constant \( c \) with \( \Pr(\max_j\lvert T_j\rvert\le c)=1-\alpha \) depends only on \( \R \) and \( \nu \), and
that it is at most the Bonferroni multiplier and at most Scheffé's multiplier for the span of the
\( \blambda_j \).
:::
