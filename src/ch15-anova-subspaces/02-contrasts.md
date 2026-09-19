# Estimable contrasts

Most questions about the treatments of a one-way experiment are contrasts. Does fertilizer help at all? Is the highest rate
worth its cost? Which pairs of treatments differ? Earlier chapters supplied the theory: estimability
([Chapter 8](../ch08-estimability/index.html)), the \( t \) test and interval
([Chapter 11](../ch11-general-linear-hypothesis/index.html), [Chapter 12](../ch12-intervals-and-bands/index.html)) and simultaneous
methods ([Chapter 13](../ch13-multiplicity/index.html)). This section collects them for the one-way analysis, adds a geometric fact about power,
and applies them to the nitrogen trial.

## Contrasts are the treatment space

Recall that a contrast is a linear function \( \psi=\sum_kc_k\mu_k \) with \( \sum_kc_k=0 \) and \( \mathbf{c}\ne\bzero \) (@def-est-contrast).
In the effects form it equals \( \sum_kc_k\alpha_k \), and the contrasts are exactly the estimable
functions of the effects alone (@thm-est-oneway-contrasts). Their estimates are \( \hat{\psi}=\sum_kc_k\bar{y}_k \), whatever
side condition or coding was used to fit the model.

Behind this lies a simple geometry. Let
\[
\mathcal T=\C(\Z)\cap\bone\perpc
\]
be the **treatment space**: the vectors that are constant within groups and sum to zero. Its dimension is \( g-1 \) (@thm-proj-nested).
A vector \( \bu\in\mathcal T \) takes some value \( v_k \) on every observation of group \( k \), with
\( \sum_kn_kv_k=0 \). Put \( c_k=n_kv_k \). Then \( \mathbf{c} \) is a contrast and \( \bu \) is the vector
\( \bu_{\mathbf{c}}=\Z\bD^{-1}\mathbf{c} \) of @eq-est-contrast-vector, with \( \bD=\diag(n_1,\dots,n_g) \). Conversely every
\( \bu_{\mathbf{c}} \) lies in \( \mathcal T \) (@prp-est-contrast-ss). So contrasts and nonzero vectors of the treatment space
correspond one to one, and \( \hat{\psi}=\bu_{\mathbf{c}}\T\y \). *A contrast is a direction in the treatment space, and its
estimate is the component of the data along that direction.*

Equal treatment means say that the mean vector \( \boldsymbol{\uptheta}=\E(\Y) \) has no component in \( \mathcal T \), that is, that every
contrast vanishes (@exr-aov-all-contrasts-zero). The overall test of [Section 15.3](03-anova-table.html) looks at all of
\( \mathcal T \); a single contrast looks along one line of it.

## The test and interval for one contrast

For a contrast \( \mathbf{c} \) write
\[
v_{\mathbf{c}}=\sum_{k=1}^g\frac{c_k^2}{n_k},
\]
so that \( \Var(\hat{\psi})=\sigma^2v_{\mathbf{c}} \) (@eq-est-contrast-var).

::: {#prp-aov-contrast-ss}
[Sum of squares and test for a contrast]

In the normal one-way model, let \( \psi=\sum_kc_k\mu_k \) be a contrast, \( \hat{\psi}=\sum_kc_k\bar{Y}_k \), and
\( \psi_0\in\Real \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\psi}\sim\Normal(\psi,\sigma^2v_{\mathbf{c}}) \), and \( \hat{\psi} \) is independent of \( s^2 \);

2. the statistic
   \[
T=\frac{\hat{\psi}-\psi_0}{s\sqrt{v_{\mathbf{c}}}}\sim t(n-g,\delta),\qquad \delta=\frac{\psi-\psi_0}{\sigma\sqrt{v_{\mathbf{c}}}},
\]
   so the test of \( \psi=\psi_0 \) rejects when \( \lvert T\rvert>t_{n-g,\alpha/2} \), and
   \( \hat{\psi}\pm t_{n-g,\alpha/2}\,s\sqrt{v_{\mathbf{c}}} \) is an exact \( 1-\alpha \) confidence interval for \( \psi \);

3. the **sum of squares of the contrast**
   \[
\text{SS}(\mathbf{c})=\frac{\hat{\psi}^2}{v_{\mathbf{c}}}=\norm{\bP_{\mathbf{c}}\Y}^2 ,
\]
   where \( \bP_{\mathbf{c}} \) projects onto the line spanned by \( \bu_{\mathbf{c}} \), has one degree of freedom, and
   \( \text{SS}(\mathbf{c})/s^2=T^2 \) for \( \psi_0=0 \) has the \( F(1,n-g,\delta^2) \) distribution. Replacing \( \mathbf{c} \) by
   \( a\mathbf{c} \) with \( a\ne0 \) changes neither \( \text{SS}(\mathbf{c}) \) nor the test;

4. \( \text{SS}(\mathbf{c})\le\sum_kn_k(\bar{Y}_k-\bar{Y})^2 \), with equality when \( c_k=n_k(\bar{Y}_k-\bar{Y}) \) for every \( k \).
:::

:::

::: {.proof}
Use the cell-means form, with model matrix \( \Z \) of full rank \( g \), so that \( (\Z\T\Z)^{-1}=\bD^{-1} \) and
\( \mathbf{c}\T\bD^{-1}\mathbf{c}=v_{\mathbf{c}} \). The function \( \mathbf{c}\T\bmu \) is estimable with estimate \( \mathbf{c}\T\hat{\bmu}=\hat{\psi} \).
(a) and (b) are @thm-glh-t-test(a) and @cor-opt-t with \( \blambda=\mathbf{c} \), \( r=g \), and the interval is @thm-ci-estimable-interval.
The independence of \( \hat{\psi} \) and \( s^2 \) is part of @thm-opt-sampling. (c) The identity
\( \hat{\psi}^2/v_{\mathbf{c}}=\norm{\bP_{\mathbf{c}}\Y}^2 \) is @prp-est-contrast-ss(b), and \( T^2 \) is the \( F \) statistic of the
hypothesis \( \psi=0 \) by @thm-glh-t-test(b), with noncentrality \( \delta^2 \). Scaling \( \mathbf{c} \) by \( a \) multiplies
\( \hat{\psi}^2 \) and \( v_{\mathbf{c}} \) by \( a^2 \), and leaves the line spanned by \( \bu_{\mathbf{c}} \) unchanged. (d) is @exr-est-max-contrast.
:::

By part (c) a contrast's sum of squares is a projection onto a line of \( \mathcal T \), and by part (d) no line captures more than
the whole space. The line that captures all of it is chosen by the data, which is why its \( t \) statistic needs Scheffé's multiplier (@thm-mc-scheffe).


## Planned contrasts in the nitrogen trial

Before the trial the agronomist wrote down two questions. First, does nitrogen raise yield at all? The natural contrast
compares the average of the four fertilized rates with the control,
\( \psi_1=\tfrac14(\mu_2+\mu_3+\mu_4+\mu_5)-\mu_1 \). Second, is the highest rate better than the moderate one?
That is \( \psi_2=\mu_5-\mu_3 \), the difference between \( 160 \) and \( 80 \) kg N/ha.

::: {#exm-aov-nitrogen-planned}
[Two planned contrasts]

The pooled standard deviation is \( s=0.5211 \) on \( 25 \) degrees of freedom. For \( \psi_1 \),
\( v_{\mathbf{c}}=\tfrac16(1+4\cdot\tfrac1{16}) \), the estimate is \( 2.110 \) t/ha with standard error \( 0.238 \), so
\( t=8.87 \) and the \( 95\% \) interval is \( [1.620,\ 2.600] \). Fertilizer raises yield by roughly two tonnes per
hectare. Its sum of squares is
\( 21.370 \).

For \( \psi_2 \) the estimate is \( 0.595 \) t/ha with standard error \( 0.301 \), \( t=1.98 \), two-sided
\( p=0.059 \), and the interval \( [-0.025,\ 1.215] \). The data are consistent with anything from no gain to a gain of more
than a tonne: with six plots per rate, half a tonne is only about \( 1.7 \) standard errors of a pairwise difference.
:::

```{.python .run #cell-nitrogen-contrasts-contrast}
import numpy as np
from scipy import stats

rates = np.array([0, 40, 80, 120, 160])            # kg N per hectare
g, m = len(rates), 6                                # 5 rates, 6 plots per rate
true_means = 6.8 - 2.7 * np.exp(-0.018 * rates)     # unknown to the analyst
rng = np.random.default_rng(1515)
group = np.repeat(np.arange(g), m)                  # plot i received rate group[i]
y = np.round(true_means[group] + rng.normal(0, 0.5, g * m), 2)   # yield, t/ha
n = len(y)

n_k = np.bincount(group).astype(float)
means = np.bincount(group, weights=y) / n_k
nu = n - g
s = np.sqrt(np.sum((y - means[group]) ** 2) / nu)

def contrast(c, psi0=0.0, mult=None):
    """Estimate, standard error, t statistic, SS and interval for sum_k c_k mu_k."""
    est = c @ means
    v = np.sum(c ** 2 / n_k)                        # Var(estimate) / sigma^2
    se = s * np.sqrt(v)
    mult = stats.t.ppf(0.975, nu) if mult is None else mult
    return est, se, (est - psi0) / se, est ** 2 / v, (est - mult * se, est + mult * se)

c_fert = np.array([-1, 0.25, 0.25, 0.25, 0.25])     # fertilized (average) minus control
c_top = np.array([0, 0, -1, 0, 1])                  # 160 against 80 kg N/ha
for name, c in (("fertilized - control", c_fert), ("160 - 80", c_top)):
    est, se, t, ss, (lo, hi) = contrast(c)
    print(f"{name:21s} {est:6.3f}  se {se:.3f}  t = {t:5.2f}  SS = {ss:6.3f}"
          f"  95% interval [{lo:.3f}, {hi:.3f}]")
```

Each was a single planned question, so each gets an unadjusted \( t \) interval; joint coverage for the pair would use Bonferroni's
\( t_{25,0.0125} \). Writing the questions down in advance fixes a small family ([Section 13.6](../ch13-multiplicity/06-choosing.html)).

## Families of comparisons

Three families come up in almost every one-way analysis, and [Chapter 13](../ch13-multiplicity/index.html) supplies a method for
each. The multipliers below multiply the standard error of a difference of two means, \( s\sqrt{2/m}=0.3009 \) here.

**Comparisons with a control.** The four differences \( \mu_k-\mu_1 \), \( k=2,\dots,5 \), form a planned family of
size four. Bonferroni's method (@thm-mc-bonferroni) uses the multiplier \( t_{25,0.05/8}=2.692 \) in place of the
unadjusted \( 2.060 \):

| Rate against control | Difference | Bonferroni \( 95\% \) interval |
|:---:|:---:|:---:|
| 40 | \( 1.232 \) | \( [0.422,\ 2.041] \) |
| 80 | \( 2.052 \) | \( [1.242,\ 2.861] \) |
| 120 | \( 2.510 \) | \( [1.700,\ 3.320] \) |
| 160 | \( 2.647 \) | \( [1.837,\ 3.456] \) |

Every fertilized rate beats the control. The four estimates share \( \bar{y}_1 \), so they are
positively correlated, with correlation \( \tfrac12 \) in a balanced design (@exr-aov-dunnett-correlation). Dunnett (1955)
computed the exact joint distribution of the four \( t \) statistics under this correlation; his multiplier here is
\( 2.606 \) (@exr-mc-control characterizes it), a few percent shorter than Bonferroni's. We do not derive its tables.

**All pairs.** Tukey's method (@thm-mc-tukey) covers all ten differences \( \mu_k-\mu_l \) with exact joint coverage in a
balanced design. With \( q_{0.05}(5,25)=4.153 \) the half-width is \( q\,s/\sqrt m=0.884 \), and six of the ten pairs are
declared different. The four that are not are \( 40 \) and \( 80 \), \( 80 \) and \( 120 \), \( 80 \) and \( 160 \), and
\( 120 \) and \( 160 \). The conclusions are not transitive: \( 40 \) is not distinguished from \( 80 \), nor \( 80 \) from \( 160 \),
yet \( 40 \) and \( 160 \) are declared different, since \( 6.652-5.237=1.415 \) exceeds \( 0.884 \). "Not declared different" means "not enough evidence of a difference", which is not an equivalence relation, so displays that
sort treatments into groups of "equal" means should be read with care, especially in unbalanced designs.

**Contrasts suggested by the data.** Looking at [Figure 15.1.1](01-oneway-model.html#fig-aov-nitrogen-data), one might
ask whether the response has levelled off, by comparing the average of the two highest rates with \( 80 \) kg N/ha,
\( \tfrac12(\mu_4+\mu_5)-\mu_3 \). The estimate is \( 0.527 \) with standard error \( 0.261 \). The unadjusted interval
\( [-0.010,\ 1.063] \) is not legitimate, because the contrast was chosen after seeing the means. Scheffé's interval (@eq-mc-scheffe-contrasts),
with multiplier \( \sqrt{4F_{0.05}(4,25)}=3.322 \), is \( [-0.339,\ 1.392] \). It is valid however the
contrast was found, and it says the question is open.

```{.python .run #cell-nitrogen-contrasts-families}
import itertools

k_ctrl = g - 1                                      # four comparisons with the control
bonf = stats.t.ppf(1 - 0.05 / (2 * k_ctrl), nu)     # Bonferroni multiplier
se_pair = s * np.sqrt(2 / m)
for k in range(1, g):
    d = means[k] - means[0]
    print(f"{rates[k]:3d} - 0: {d:6.3f}  Bonferroni interval"
          f" [{d - bonf * se_pair:.3f}, {d + bonf * se_pair:.3f}]")

q = stats.studentized_range.ppf(0.95, g, nu)        # q_{0.05}(g, nu)
hsd = q * s / np.sqrt(m)                            # Tukey half-width, all pairs
pairs = list(itertools.combinations(range(g), 2))
different = [(int(rates[k]), int(rates[l])) for k, l in pairs if abs(means[k] - means[l]) > hsd]
print(f"Tukey half-width {hsd:.3f}; pairs declared different: {len(different)} of {len(pairs)}")
print("not declared different:",
      [(int(rates[k]), int(rates[l])) for k, l in pairs if abs(means[k] - means[l]) <= hsd])
```

## Why planned contrasts are powerful

A single contrast is tested with one numerator degree of freedom, the overall test with \( g-1 \). When the scientific
question is sharp, the contrast is the better test, and the reason is geometric. Write \( \bP_{\mathcal T} \) for the projection
onto the treatment space and
\[
\gamma_{\mathcal T}=\frac{\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2}{\sigma^2}=\frac1{\sigma^2}\sum_{k=1}^gn_k(\mu_k-\bar{\mu})^2 ,
\]
the squared length of the treatment component of the mean in units of \( \sigma^2 \). ([Section 15.3](03-anova-table.html) shows
that it is the noncentrality of the overall \( F \) test.)

::: {#prp-aov-contrast-angle}
[The noncentrality of a contrast]

Let \( \mathbf{c} \) be a contrast with vector \( \bu_{\mathbf{c}}\in\mathcal T \), and suppose \( \bP_{\mathcal T}\boldsymbol{\uptheta}\ne\bzero \). The
noncentrality of the test of \( \psi=0 \) in @prp-aov-contrast-ss is
\[
\delta^2=\frac{\psi^2}{\sigma^2v_{\mathbf{c}}}=\gamma_{\mathcal T}\cos^2\varphi ,
\]
where \( \varphi \) is the angle between \( \bu_{\mathbf{c}} \) and \( \bP_{\mathcal T}\boldsymbol{\uptheta} \). In particular \( \delta^2\le\gamma_{\mathcal T} \),
with equality iff \( \bu_{\mathbf{c}} \) is parallel to \( \bP_{\mathcal T}\boldsymbol{\uptheta} \), that is, iff \( c_k\propto n_k(\mu_k-\bar{\mu}) \).
:::

::: {.proof}
By @prp-est-contrast-ss(a) applied to the mean vector, \( \psi=\bu_{\mathbf{c}}\T\boldsymbol{\uptheta} \), and \( v_{\mathbf{c}}=\norm{\bu_{\mathbf{c}}}^2 \) by the
computation in the proof of that result. Since \( \bu_{\mathbf{c}}\in\mathcal T \), \( \bu_{\mathbf{c}}\T\boldsymbol{\uptheta}=\bu_{\mathbf{c}}\T\bP_{\mathcal T}\boldsymbol{\uptheta} \). So
\[
\delta^2=\frac{(\bu_{\mathbf{c}}\T\bP_{\mathcal T}\boldsymbol{\uptheta})^2}{\sigma^2\norm{\bu_{\mathbf{c}}}^2}
=\frac{\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2}{\sigma^2}\cdot\frac{(\bu_{\mathbf{c}}\T\bP_{\mathcal T}\boldsymbol{\uptheta})^2}{\norm{\bu_{\mathbf{c}}}^2\norm{\bP_{\mathcal T}\boldsymbol{\uptheta}}^2}
=\gamma_{\mathcal T}\cos^2\varphi
\]
by @eq-proj-cosine. Equality in \( \cos^2\varphi\le1 \) holds iff the two vectors are parallel. The vector
\( \bP_{\mathcal T}\boldsymbol{\uptheta} \) has entry \( \mu_k-\bar{\mu} \) in group \( k \) (@thm-proj-nested), so it is \( \bu_{\mathbf{c}} \) for
\( c_k=n_k(\mu_k-\bar{\mu}) \).
:::

If the departure from equal means lies along the planned contrast, the contrast test gets all of the noncentrality on one
degree of freedom, while the overall test spreads the same noncentrality over \( g-1 \). For the nitrogen design, with \( 25 \)
error degrees of freedom, the powers at level \( 0.05 \) are:

| \( \gamma_{\mathcal T} \), all along the contrast | \( 4 \) | \( 8 \) | \( 12 \) |
|:---|:---:|:---:|:---:|
| contrast test, \( F(1,25) \) | \( 0.485 \) | \( 0.776 \) | \( 0.914 \) |
| overall test, \( F(4,25) \) | \( 0.269 \) | \( 0.516 \) | \( 0.714 \) |

The price is paid when the question is wrong: a contrast orthogonal to \( \bP_{\mathcal T}\boldsymbol{\uptheta} \) has power \( \alpha \) however large
the treatment differences, and a direction chosen at random captures on average only \( 1/(g-1) \) of \( \gamma_{\mathcal T} \) (@exr-aov-random-direction).

```{.python .run #cell-nitrogen-contrasts-power}
def power_f(gamma, q, nu, alpha=0.05):
    return stats.ncf.sf(stats.f.ppf(1 - alpha, q, nu), q, nu, gamma)

for gamma in (4.0, 8.0, 12.0):                      # all of it along the planned contrast
    print(f"gamma = {gamma:4.1f}: contrast test {power_f(gamma, 1, nu):.3f},"
          f" omnibus F test {power_f(gamma, g - 1, nu):.3f}")
```

## Exercises

### A. Check your understanding

::: {#exr-aov-contrast-se}
[A1]

In the nitrogen trial, compute by hand the estimate and standard error of
\( \tfrac12(\mu_4+\mu_5)-\tfrac12(\mu_1+\mu_2) \), the average of the two highest rates minus the average of the two lowest, from
the rate means of [Section 15.1](01-oneway-model.html) and \( s=0.5211 \).
:::

::: {.solution}
The estimate is \( \tfrac12(6.515+6.652)-\tfrac12(4.005+5.237)=6.5835-4.621=1.9625 \). With \( m=6 \),
\( v_{\mathbf{c}}=4\cdot\tfrac14/6=1/6 \), so the standard error is \( 0.5211/\sqrt6\approx0.213 \).
:::

::: {#exr-aov-contrast-scale}
[A2]

Show directly that the \( t \) statistic of @prp-aov-contrast-ss(b) with \( \psi_0=0 \) is unchanged when \( \mathbf{c} \) is multiplied by
\( a>0 \) and changes sign when \( a<0 \). What happens to the confidence interval?
:::

### B. Practice

::: {#exr-aov-balanced-contrast-power}
[B1]

In a balanced layout with \( m \) observations per group, show that the noncentrality of the test of a contrast is
\( \delta^2=m\psi^2/(\sigma^2\sum_kc_k^2) \). How many observations per group are needed for the test of
\( \mu_5-\mu_3 \) in the nitrogen trial to have power \( 0.8 \) at level \( 0.05 \) when the true difference is \( 0.5 \) t/ha and
\( \sigma=0.5 \)? (Use software for the noncentral \( F \), or the normal approximation.)
:::

::: {.solution}
Here \( v_{\mathbf{c}}=\sum_kc_k^2/m \), so \( \delta^2=\psi^2/(\sigma^2v_{\mathbf{c}})=m\psi^2/(\sigma^2\sum_kc_k^2) \). For \( \mu_5-\mu_3 \),
\( \sum_kc_k^2=2 \), \( \psi/\sigma=1 \), and \( \delta^2=m/2 \). The normal approximation asks for
\( \delta\approx z_{0.025}+z_{0.2}=1.96+0.84=2.80 \), so \( m\approx2\cdot2.80^2\approx15.7 \). The exact noncentral \( F \) with
\( 5(m-1) \) error degrees of freedom gives \( m=17 \), since the \( t \) quantile exceeds the normal one. The trial's
\( m=6 \) is far short.
:::

::: {#exr-aov-contrast-allocation}
[B2]

A single contrast \( \mathbf{c} \) matters above all. For a fixed total \( n \), show that \( v_{\mathbf{c}}=\sum_kc_k^2/n_k \) is minimized, over
positive real \( n_k \) with \( \sum_kn_k=n \), by \( n_k\propto\lvert c_k\rvert \), with minimum \( (\sum_k\lvert c_k\rvert)^2/n \). What
does this say about how to allocate plots between two treatments, and between a control and the average of four
treatments?
:::

::: {.solution}
By the Cauchy–Schwarz inequality,
\( \bigl(\sum_k\lvert c_k\rvert\bigr)^2=\bigl(\sum_k(\lvert c_k\rvert/\sqrt{n_k})\sqrt{n_k}\bigr)^2\le\sum_k(c_k^2/n_k)\sum_kn_k=nv_{\mathbf{c}} \),
with equality iff \( \lvert c_k\rvert/\sqrt{n_k}\propto\sqrt{n_k} \), that is, \( n_k\propto\lvert c_k\rvert \). Groups with \( c_k=0 \) get no
observations. For \( \mu_1-\mu_2 \), equal allocation is best. For \( \mu_1-\tfrac14(\mu_2+\dots+\mu_5) \), the control should get as many
plots as the four treatments together, \( n_1=n/2 \) and \( n_k=n/8 \) for the others.
:::

::: {#exr-aov-sqrt-rule}
[B3]

In a comparison of \( g-1 \) treatments with a control, suppose the control gets \( n_0 \) observations and each treatment gets
\( n_t \), with \( n_0+(g-1)n_t=n \) fixed. Show that \( \sum_{k}\Var(\bar{Y}_k-\bar{Y}_0) \) is minimized when
\( n_0=\sqrt{g-1}\,n_t \). For \( g=5 \) and \( n=30 \), what allocation does this suggest?
:::

::: {.solution}
The sum is \( \sigma^2(g-1)(1/n_0+1/n_t) \). Minimizing \( 1/n_0+1/n_t \) subject to \( n_0+(g-1)n_t=n \) with a Lagrange multiplier gives
\( 1/n_0^2=\lambda \) and \( 1/n_t^2=(g-1)\lambda \), so \( n_0/n_t=\sqrt{g-1} \). For \( g=5 \) the control should get twice as many
observations as each treatment: \( n_0=10 \) and \( n_t=5 \). Contrast @exr-aov-contrast-allocation: for the average of the four
comparisons the optimal control share is \( n/2 \), that is \( n_0=4n_t \). The best allocation depends on whether the four
comparisons or their average is the target.
:::

::: {#exr-aov-dunnett-correlation}
[B4]

With \( n_0 \) observations on the control and \( n_t \) on each treatment, show that the estimates \( \bar{Y}_k-\bar{Y}_0 \) and
\( \bar{Y}_l-\bar{Y}_0 \) of two treatment-versus-control differences have correlation \( n_t/(n_t+n_0) \). In a balanced design it is
\( \tfrac12 \).
:::

::: {.solution}
The covariance is \( \Var(\bar{Y}_0)=\sigma^2/n_0 \), since the treatment means are independent of each other and of \( \bar{Y}_0 \). Each
variance is \( \sigma^2(1/n_t+1/n_0) \). The ratio is \( (1/n_0)/(1/n_t+1/n_0)=n_t/(n_t+n_0) \).
:::

### C. Going deeper

::: {#exr-aov-random-direction}
[C1]

Let \( q=g-1\ge2 \) and let \( \bu \) be uniformly distributed on the unit sphere of the \( q \)-dimensional space \( \mathcal T \) (for
example \( \bu=\bw/\norm{\bw} \) with \( \bw\sim\Normal \) on \( \mathcal T \) with covariance \( \bP_{\mathcal T} \)). Show that
\( \E\cos^2\varphi=1/q \), where \( \varphi \) is the angle between \( \bu \) and a fixed nonzero \( \bv\in\mathcal T \). Interpret this
with @prp-aov-contrast-angle.
:::

::: {.solution}
Choose an orthonormal basis of \( \mathcal T \) whose first vector is \( \bv/\norm{\bv} \). In these coordinates \( \cos^2\varphi=u_1^2 \). By
symmetry the \( u_i^2 \) have a common mean, and they sum to \( 1 \), so \( \E u_1^2=1/q \). By @prp-aov-contrast-angle, a contrast whose
direction is unrelated to the truth captures on average a fraction \( 1/q \) of \( \gamma_{\mathcal T} \) on its one degree of freedom,
the same average share per degree of freedom as the overall test. Knowledge, not the count of degrees of freedom, is what makes a planned contrast
more powerful.
:::

::: {#exr-aov-contrast-family-rank}
[C2]

The Scheffé multiplier for a family of contrasts depends on the dimension of the space the family spans (@thm-mc-scheffe).
Show that the \( g-1 \) comparisons with a control span \( \mathcal T \), and that the family of the four
"successive differences" \( \mu_{k+1}-\mu_k \) does too. Why, then, are Bonferroni or Dunnett intervals shorter than Scheffé's
for these families?
:::
