# Random effects

A one-way layout (@def-aov-oneway) has \( g \) groups and \( g \) group means. In
Part IV those means were the object of the study: the five nitrogen rates were
chosen, and every conclusion was about them.
[Section 15.1](../ch15-anova-subspaces/01-oneway-model.html) previewed the other
possibility and left it here; this section carries the preview through. Consider a
study of a different shape. A reference material with a known lead content is sent to twelve
accredited laboratories, each of which analyses four aliquots of it. Nobody cares
whether laboratory 3 reads higher than laboratory 7. The question is how much
accredited laboratories disagree with one another, because that number, and not
the twelve particular deviations, is what a client using an unnamed laboratory
faces. The twelve are a sample from a population of laboratories, and the
conclusion must reach the population.

## When a factor is random

The test, met first in [Section 15.1](../ch15-anova-subspaces/01-oneway-model.html),
is about the *scope of the conclusion*. If the levels in the data exhaust
the levels you want to talk about, the factor is fixed; if they are a sample and the
conclusion should extend to levels you did not observe, it is random. Doses,
treatments, varieties and the two arms of a trial are almost always fixed;
laboratories, batches of reagent, subjects, litters, classrooms and interviewers are
usually random. Three things change when a factor is declared random.

::: {.enumerate options="label=(\arabic*)"}
1. **The parameters change.** The \( g \) group effects are replaced by a single
   number, their variance. A study with \( g=12 \) laboratories estimates one
   variance, not twelve means, and its precision is governed by \( g \), not by
   the total number of observations.

2. **Estimation becomes prediction.** The realized effect of laboratory 3 is a
   random variable, not a parameter. Asking for its value is a prediction
   problem, and [Section 32.3](03-blup.html) shows that the best answer is not the
   laboratory's own average.

3. **The covariance changes.** Two aliquots analysed by the same laboratory share
   its effect, so they are correlated. The model is no longer a linear model with
   \( \Cov(\Y)=\sigma^2\I \), and this is the reason the chapter belongs to
   Part VII.
:::

::: {.warning}
A variance is estimated from the number of *levels*, not the number of observations.
With \( g=3 \) laboratories there are two degrees of freedom for the
between-laboratory variance, and any interval for it is nearly useless however many
aliquots each laboratory analyses. A treated and a control arm should never be a
random factor: two levels carry no information about a population of levels, and
there is no population. When a factor has few levels but its levels really are a
sample, the honest options are to report the fixed-effect analysis alongside, or to
supply prior information ([Section 32.6](06-bayes.html)).
:::

## The one-way random effects model

::: {#def-mix-oneway}
[One-way random effects model]

Let \( g\ge2 \) groups contain \( n_k\ge1 \) observations each, with
\( n=\sum_kn_k \). The **one-way random effects model** is
\[
Y_{kj}=\mu+a_k+\varepsilon_{kj},\qquad k=1,\dots,g,\ \ j=1,\dots,n_k,
\]
where \( \mu \) is an unknown constant, the **group effects** \( a_1,\dots,a_g \)
are uncorrelated with mean \( 0 \) and variance \( \sigma_a^2\ge0 \), the errors
\( \varepsilon_{kj} \) are uncorrelated with mean \( 0 \) and variance
\( \sigma^2>0 \), and the two sets are uncorrelated with each other. In matrix
form, with \( \Z \) the \( n\times g \) matrix of group indicators,
\[
\Y=\mu\bone+\Z\mathbf{a}+\be .
\]
The **normal** one-way random effects model adds
\( \mathbf{a}\sim\Normal_g(\bzero,\sigma_a^2\I) \) independent of
\( \be\sim\Normal_n(\bzero,\sigma^2\I) \). The layout is **balanced** if
\( n_1=\dots=n_g=m \). The pair \( (\sigma_a^2,\sigma^2) \) are the **variance
components**.
:::

The model has three parameters whatever \( g \) is, where @def-aov-oneway had
\( g+1 \). That is the economy of random effects, and also their risk: it is bought
by assuming that the \( g \) effects are an uncorrelated sample from one
distribution.

::: {#prp-mix-oneway}
[The one-way random effects model]

In @def-mix-oneway, write \( \gamma=\sigma_a^2/\sigma^2 \) and let \( \M \) be the
projection onto \( \C(\Z) \), which replaces each entry by its group mean.

::: {.enumerate options="label=(\alph*)"}
1. \( \E(\Y)=\mu\bone \) and
   \( \Cov(\Y)=\sigma^2\I+\sigma_a^2\Z\Z\T \), a block diagonal matrix whose
   \( k \)th block is \( \sigma^2\I_{n_k}+\sigma_a^2\bone\bone\T \).

2. *(Intraclass correlation.)* Observations in different groups are uncorrelated,
   and two observations in the same group have correlation
   \[
   \rho=\frac{\sigma_a^2}{\sigma_a^2+\sigma^2}=\frac{\gamma}{1+\gamma}\in[0,1).
   \]

3. *(Moment estimates.)* Let \( \text{MSB} \) and \( \text{MSE} \) be the between
   and within mean squares of the one-way analysis of variance, on \( g-1 \) and
   \( n-g \) degrees of freedom. Then \( \E(\text{MSE})=\sigma^2 \) and
   \( \E(\text{MSB})=\sigma^2+m_0\sigma_a^2 \) with
   \[
   m_0=\frac{1}{g-1}\Bigl(n-\frac{\sum_kn_k^2}{n}\Bigr),
   \]{#eq-mix-m0}

   which equals \( m \) in a balanced layout. Hence
   \( \hat\sigma^2=\text{MSE} \) and
   \( \hat\sigma_a^2=(\text{MSB}-\text{MSE})/m_0 \) are unbiased.

4. *(Distributions, balanced and normal.)* If the layout is balanced and the model
   is normal, then with \( \lambda_1=\sigma^2+m\sigma_a^2 \),
   \[
   \frac{\text{SSB}}{\lambda_1}\sim\chi^2(g-1)
   \quad\text{independently of}\quad
   \frac{\text{SSE}}{\sigma^2}\sim\chi^2(n-g),
   \]
   and therefore \( \text{MSB}/\text{MSE}\sim(1+m\gamma)\,F(g-1,n-g) \).
:::

:::

::: {.proof}
(a) By @thm-rv-linear applied to \( \Y=\mu\bone+[\Z,\I](\mathbf{a}\T,\be\T)\T \),
\( \Cov(\Y)=\Z(\sigma_a^2\I_g)\Z\T+\sigma^2\I \). The matrix \( \Z\Z\T \) has
\( (i,i') \) entry \( 1 \) when observations \( i \) and \( i' \) are in the same
group and \( 0 \) otherwise, which is the block form stated.

(b) Read off (a): the covariance of two distinct observations in group \( k \) is
\( \sigma_a^2 \) and each variance is \( \sigma_a^2+\sigma^2 \).

(c) The two sums of squares are \( \text{SSB}=\Y\T(\M-\bP_1)\Y \) and
\( \text{SSE}=\Y\T(\I-\M)\Y \), where \( \bP_1=n^{-1}\bone\bone\T \). Both
projections annihilate \( \bone \), so @thm-rv-quadform-mean leaves only the trace
term. Since the columns of \( \Z\Z\T \) lie in \( \C(\Z) \) we have
\( \M\Z\Z\T=\Z\Z\T \), and \( (\I-\M)\Z=\bzero \), so
\[
\E(\text{SSE})=\sigma^2\tr(\I-\M)=\sigma^2(n-g),
\]
while \( \tr(\Z\Z\T)=n \) and
\( \tr(\bP_1\Z\Z\T)=n^{-1}\bone\T\Z\Z\T\bone=n^{-1}\sum_kn_k^2 \) give
\[
\E(\text{SSB})=\sigma^2(g-1)+\sigma_a^2\Bigl(n-\frac{\sum_kn_k^2}{n}\Bigr).
\]
Dividing by the degrees of freedom gives the two expectations, and unbiasedness of
the two estimates follows by linearity. In a balanced layout \( n=gm \) and
\( \sum_kn_k^2=gm^2 \), so \( m_0=\{gm-gm^2/(gm)\}/(g-1)=m(g-1)/(g-1)=m \).

(d) In the balanced case \( \Z\Z\T=m\M \), so by (a)
\[
\Cov(\Y)=\lambda_1\M+\sigma^2(\I-\M),
\]{#eq-mix-balanced-cov}

a spectral decomposition into the two orthogonal pieces \( \C(\Z) \) and
\( \C(\Z)\perpc \). Write \( \bSigma=\Cov(\Y) \). For the between projection
\( \A=\M-\bP_1 \), which satisfies \( \A\M=\A \) and \( \A(\I-\M)=\bzero \), we get
\( \A\bSigma=\lambda_1\A \), so \( \A\bSigma/\lambda_1=\A \) is idempotent of rank
\( g-1 \) and \( \A\bmu=\bzero \); @thm-qf-chisq gives
\( \text{SSB}/\lambda_1\sim\chi^2(g-1) \). For \( \B=\I-\M \) we get
\( \B\bSigma=\sigma^2\B \) and \( \text{SSE}/\sigma^2\sim\chi^2(n-g) \) in the same
way. Independence follows from @thm-qf-indep-quadratic, because
\( \A\bSigma\B=\lambda_1\A\B=\bzero \). The ratio of the two mean squares is then
\( (\lambda_1/\sigma^2) \) times an \( F(g-1,n-g) \) variable, and
\( \lambda_1/\sigma^2=1+m\gamma \).
:::

Part (c) is @exr-ss-random-oneway carried to unbalanced layouts, and part (d) is the
fact that makes the balanced case so pleasant: the same \( F \) ratio that tests
equality of fixed group means tests \( \sigma_a^2=0 \), and it does more, because its
distribution is known for every \( \gamma \), not only at \( \gamma=0 \). Inverting
it gives an exact interval.

::: {#cor-mix-icc-interval}
[Exact interval for the variance ratio]

In the balanced normal one-way random effects model, let \( F=\text{MSB}/\text{MSE} \)
and let \( F_q \) denote the \( q \) quantile of \( F(g-1,n-g) \). Then
\[
\Bigl[\frac{1}{m}\Bigl(\frac{F}{F_{1-\alpha/2}}-1\Bigr),\ \
      \frac{1}{m}\Bigl(\frac{F}{F_{\alpha/2}}-1\Bigr)\Bigr]
\]{#eq-mix-gamma-interval}

is an exact \( 1-\alpha \) confidence interval for \( \gamma=\sigma_a^2/\sigma^2 \),
and applying the increasing map \( \gamma\mapsto\gamma/(1+\gamma) \) to its endpoints
gives an exact \( 1-\alpha \) interval for the intraclass correlation \( \rho \).
Endpoints below zero are replaced by zero.
:::

::: {.proof}
By @prp-mix-oneway(d), \( F/(1+m\gamma)\sim F(g-1,n-g) \), so
\( \Pr\{F_{\alpha/2}\le F/(1+m\gamma)\le F_{1-\alpha/2}\}=1-\alpha \) whatever
\( \gamma \) is. Solving the two inequalities for \( \gamma \) gives
@eq-mix-gamma-interval, an event of the same probability. A strictly increasing
transformation of both endpoints keeps the coverage, and
\( \rho=\gamma/(1+\gamma) \) is strictly increasing on \( [0,\infty) \). Truncating
at zero can only raise the coverage, since \( \gamma\ge0 \).
:::

Notice what @eq-mix-gamma-interval is not: it is an interval for the *ratio*. No
exact interval of this kind exists for \( \sigma_a^2 \) alone, because
\( \hat\sigma_a^2 \) is a difference of two mean squares with different scale
factors and its distribution depends on both components separately. Approximate
intervals for \( \sigma_a^2 \) come either from Satterthwaite's method
([Section 32.5](05-inference.html)) or from the profile likelihood
([Section 32.4](04-likelihood.html)).

::: {#exm-mix-proficiency}
[A laboratory proficiency study]

Twelve laboratories each analyse four aliquots of one homogeneous batch of soil for
lead, in mg/kg; the data are simulated so that the truth is known, with
\( \mu=2.50 \), \( \sigma_a=0.20 \) and \( \sigma=0.25 \), so that the population
intraclass correlation is \( 0.390 \). The analysis of variance gives
\( \text{MSB}=0.2166 \) on \( 11 \) degrees of freedom and
\( \text{MSE}=0.0620 \) on \( 36 \), so
\[
\hat\sigma^2=0.0620,\qquad
\hat\sigma_a^2=\frac{0.2166-0.0620}{4}=0.0386 ,
\]
that is, a within-laboratory (repeatability) standard deviation of
\( 0.249 \) and a between-laboratory standard deviation of
\( 0.197 \). The estimated intraclass correlation is
\( 0.384 \): about two fifths of the variation in a single
measurement is attributable to which laboratory made it.

The \( F \) ratio is \( 3.49 \), and @eq-mix-gamma-interval gives
\( \gamma\in(0.118,\,2.439) \), that is
\( \rho\in(0.105,\,0.709) \). Twelve laboratories and
\( 48 \) measurements pin the intraclass correlation down only to
a factor of seven. That is not a defect of the method; it is what
\( 11 \) degrees of freedom buy.
:::

::: {when-format="html"}
![**Figure 32.1.1.** (a) The proficiency study: four aliquots in each of twelve
laboratories, with the laboratory means (short bars) and the overall mean (dashed).
(b) The probability that the moment estimate \( \hat\sigma_a^2 \) is negative, as a
function of the intraclass correlation, for \( g=6 \), \( 12 \) and \( 25 \)
laboratories with \( m=4 \) aliquots each.](lab_study.svg){#fig-mix-lab width=100%}
:::

::: {when-format="pdf"}
![(a) The proficiency study: four aliquots in each of twelve laboratories, with the
laboratory means (short bars) and the overall mean (dashed). (b) The probability
that the moment estimate \( \hat\sigma_a^2 \) is negative, as a function of the
intraclass correlation, for \( g=6 \), \( 12 \) and \( 25 \) laboratories with
\( m=4 \) aliquots each.](lab_study.pdf){width=100%}
:::

```{.python .run #cell-lab-study-data}
import numpy as np
from scipy import stats

g, m = 12, 4                                    # laboratories, aliquots each
n = g * m
mu, sigma_a, sigma = 2.50, 0.20, 0.25           # the values used to simulate

rng = np.random.default_rng(20320048)
lab_effect = sigma_a * rng.standard_normal(g)
y = mu + np.repeat(lab_effect, m) + sigma * rng.standard_normal(n)

Y = y.reshape(g, m)                             # one row per laboratory
lab_mean, grand_mean = Y.mean(axis=1), Y.mean()

ss_between = m * np.sum((lab_mean - grand_mean) ** 2)
ss_within = np.sum((Y - lab_mean[:, None]) ** 2)
ms_between = ss_between / (g - 1)
ms_within = ss_within / (n - g)

sigma2_hat = ms_within                          # E(MS within) = sigma^2
sigma2_a_hat = (ms_between - ms_within) / m     # E(MS between) = sigma^2 + m sigma_a^2
icc_hat = sigma2_a_hat / (sigma2_a_hat + sigma2_hat)

print(f"MS between {ms_between:.4f}   MS within {ms_within:.4f}")
print(f"sigma^2 {sigma2_hat:.4f}   sigma_a^2 {sigma2_a_hat:.4f}   ICC {icc_hat:.3f}")
```

```{.python .run #cell-lab-study-interval}
f_obs = ms_between / ms_within                  # ~ (1 + m gamma) F(g-1, n-g)
gamma_lo = (f_obs / stats.f.ppf(0.975, g - 1, n - g) - 1) / m
gamma_hi = (f_obs / stats.f.ppf(0.025, g - 1, n - g) - 1) / m
icc_lo, icc_hi = gamma_lo / (1 + gamma_lo), gamma_hi / (1 + gamma_hi)

print(f"F = {f_obs:.2f};  gamma in ({gamma_lo:.3f}, {gamma_hi:.3f})")
print(f"intraclass correlation in ({icc_lo:.3f}, {icc_hi:.3f})")
```

## Negative estimates

Nothing in @prp-mix-oneway(c) prevents \( \text{MSB} \) from falling below
\( \text{MSE} \), and then the unbiased estimate of a variance is negative. In the
balanced normal model the probability is exactly
\[
\Pr(\hat\sigma_a^2<0)=\Pr\Bigl(F(g-1,n-g)<\frac{1}{1+m\gamma}\Bigr),
\]{#eq-mix-negative}

by @prp-mix-oneway(d). Panel (b) of [Figure 32.1.1](01-random-effects.html#fig-mix-lab)
plots it. With \( g=6 \) laboratories, \( m=4 \) aliquots and an intraclass
correlation of \( 0.05 \) the probability is
\( 0.453 \); when \( \sigma_a^2 \) is exactly zero it is
\( 0.554 \), a little over one half, because
\( \Pr\{F(5,18)<1\} \) exceeds \( 1/2 \). In the study of @exm-mix-proficiency, with
twelve laboratories and a true correlation of \( 0.390 \), it is only
\( 0.0143 \).

A negative estimate is not a computational failure but evidence, and there are three
honest responses.

::: {.enumerate options="label=(\arabic*)"}
1. **Report it.** A meta-analysis that averages several studies needs the unbiased
   value, negative or not; truncating every study at zero and then averaging biases
   the average upwards.

2. **Truncate.** The estimate \( \max(\hat\sigma_a^2,0) \) is biased upwards but has
   smaller mean squared error for every parameter value (@exr-mix-truncation).
   Maximum likelihood and REML produce it automatically, since they maximize over
   \( \sigma_a^2\ge0 \).

3. **Question the model.** A persistently negative estimate can mean that
   observations within a group are genuinely *negatively* correlated, as happens when
   the group is a fixed resource the observations compete for: litter-mates sharing a
   food supply, or plots sharing water. The matrix
   \( \sigma^2\I+\sigma_a^2\Z\Z\T \) with \( \sigma_a^2<0 \) is a perfectly good
   covariance as long as \( \sigma^2+m\sigma_a^2>0 \), and fitting it as such is
   more honest than pretending a variance was estimated.
:::

## Unbalanced layouts, and beyond one way

Part (c) is the simplest instance of a general recipe: write down quadratic forms,
equate them to their expectations, and solve. Henderson (1953) organized this into
three methods for unbalanced layouts, of which the one used here is the first. The
recipe always gives unbiased estimates and never needs normality, but it is not
unique — different quadratic forms give different estimates, and in unbalanced
layouts none of them is best. That is why
[Section 32.4](04-likelihood.html) turns to the likelihood.

The same calculation runs through any balanced layout. In the two-stage nested design
of @exr-tw-random-classrooms, with \( a \) levels of \( A \), \( b \) of \( B \)
within each and \( m \) observations per cell, that exercise gave
\( \E\,\text{MS}_{B(A)}=\sigma^2+m\sigma_B^2 \), so the inner component is estimated
by \( (\text{MS}_{B(A)}-\text{MSE})/m \). Taking \( A \) random too, with effects of
variance \( \sigma_A^2 \), \( \text{MS}_A \) is \( bm \) times the sample variance of
\( a \) independent group means of variance
\( \sigma_A^2+\sigma_B^2/b+\sigma^2/(bm) \), so
\( \E\,\text{MS}_A=\sigma^2+m\sigma_B^2+bm\sigma_A^2 \) and the outer component is
estimated by \( (\text{MS}_A-\text{MS}_{B(A)})/(bm) \). Every mean square in a
balanced layout turns out to estimate \( \sigma^2 \) plus a positive combination of
the components nested inside its factor, a fact
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) establishes and uses
to read a split-plot analysis off its strata.

Finally, a warning about the grand mean. In the fixed-effects model
\( \Var(\bar Y)=\sigma^2/n \), so precision is bought with observations. Here it is
different: @exr-aov-random-effects-mean showed that
\[
\Var(\bar Y)=\frac{\sigma^2}{n}+\sigma_a^2\frac{\sum_kn_k^2}{n^2}
\ \xrightarrow[\text{balanced}]{}\ \frac{\sigma^2+m\sigma_a^2}{gm},
\]
which for fixed \( g \) is bounded below by \( \sigma_a^2/g \) however large \( m \)
grows. Measuring the same twelve laboratories a thousand times each does not improve
the estimate of the population mean beyond a limit set by \( g=12 \). This is the
commonest way of overstating precision in clustered data, and
[Chapter 33](../ch33-clustered-longitudinal-splitplot/index.html) quantifies it as the
design effect.

## Exercises

### A. Check your understanding

::: {#exr-mix-fixed-or-random}
[A1]

For each study, say whether the named factor is better treated as fixed or random, and
what the corresponding target of inference is. (a) Four fertilizer formulations
compared on a research farm. (b) Twenty operators sampled from a factory's workforce,
each measuring the same set of parts. (c) The two arms of a randomized trial.
(d) Thirty schools sampled from a district, in a study of a teaching method.
(e) Three furnaces, the only three the plant owns, in a study of their temperature
uniformity.
:::

### B. Practice

::: {#exr-mix-m0}
[B1]

Verify @eq-mix-m0 directly from the group means: show that
\( \E\{\sum_kn_k(\bar Y_k-\bar Y)^2\}=\sigma^2(g-1)+\sigma_a^2(n-\sum_kn_k^2/n) \) by
computing \( \Var(\bar Y_k) \), \( \Var(\bar Y) \) and \( \Cov(\bar Y_k,\bar Y) \).
Confirm that \( m_0=m \) when the layout is balanced, and that
\( m_0\le(n-\min_kn_k)/(g-1) \) in general.
:::

::: {.solution}
\( \bar Y_k=\mu+a_k+\bar\varepsilon_k \) has variance
\( \sigma_a^2+\sigma^2/n_k \), and \( \bar Y=\mu+\sum_kn_ka_k/n+\bar\varepsilon \) has
variance \( \sigma_a^2\sum_kn_k^2/n^2+\sigma^2/n \). Also
\( \Cov(\bar Y_k,\bar Y)=\sigma_a^2n_k/n+\sigma^2/n \). Hence
\( \E(\bar Y_k-\bar Y)^2=\sigma_a^2(1-2n_k/n+\sum_jn_j^2/n^2)+\sigma^2(1/n_k-1/n) \).
Multiplying by \( n_k \) and summing, the \( \sigma^2 \) terms give
\( \sigma^2(g-1) \) and the \( \sigma_a^2 \) terms give
\( \sigma_a^2(n-2\sum_kn_k^2/n+\sum_kn_k^2/n)=\sigma_a^2(n-\sum_kn_k^2/n) \). For the
balanced case, \( n=gm \) and \( \sum_kn_k^2=gm^2 \), so
\( m_0=(gm-m)/(g-1)=m \). In general \( \sum_kn_k^2\ge n\min_kn_k \), whence
\( m_0\le(n-\min_kn_k)/(g-1) \).
:::

::: {#exr-mix-variance-of-estimate}
[B2]

In the balanced normal model, use @prp-mix-oneway(d) to show that
\[
\Var(\hat\sigma_a^2)=\frac{2}{m^2}\Bigl\{\frac{(\sigma^2+m\sigma_a^2)^2}{g-1}
+\frac{\sigma^4}{n-g}\Bigr\}.
\]
With \( \sigma^2 \) and \( \sigma_a^2 \) fixed and \( n=gm \) held constant, is the
variance smaller with many small groups or few large ones?
:::

::: {.solution}
\( \text{MSB} \) and \( \text{MSE} \) are independent, with
\( \Var(\text{MSB})=\lambda_1^2\cdot 2/(g-1) \) and
\( \Var(\text{MSE})=\sigma^4\cdot 2/(n-g) \), because a \( \chi^2(r) \) variable has
variance \( 2r \). The formula follows from
\( \hat\sigma_a^2=(\text{MSB}-\text{MSE})/m \). Neither extreme answers the second
question. With \( g-1\approx n/m \) and \( n-g\approx n \),
\[
\Var(\hat\sigma_a^2)\approx\frac{2}{n}
\Bigl(\frac{\sigma^4}{m}+2\sigma^2\sigma_a^2+m\sigma_a^4
+\frac{\sigma^4}{m^2}\Bigr),
\]
a decreasing part plus an increasing one, minimized near
\( m=\sigma^2/\sigma_a^2 \), that is at \( m\gamma\approx1 \): many small groups win
only once each group is large enough for \( \text{MSB} \) to stand clear of
\( \text{MSE} \). In the exact formula with \( n=100 \) and \( \sigma^2=1 \), at
\( \sigma_a^2=0.01 \) the variance is
\( 0.0206 \), \( 0.0056 \),
\( 0.0029 \) and \( 0.0017 \) for
\( m=2,5,10,25 \) — four groups of twenty-five beat fifty pairs by a factor of twelve — while
at \( \sigma_a^2=1 \) it is \( 0.1018 \),
\( 0.1526 \) and \( 0.2691 \) for
\( m=2,5,10 \), and the small groups win. The extreme \( m=1 \) always fails, because
then the two components are not separately identified.
:::

::: {#exr-mix-negative-probability}
[B3]

Using @eq-mix-negative, compute \( \Pr(\hat\sigma_a^2<0) \) when \( \sigma_a^2=0 \),
for \( (g,m)=(5,4) \) and \( (g,m)=(25,4) \). Explain why the probability does not
tend to \( 0 \) as \( g\to\infty \) with \( \sigma_a^2=0 \) fixed, and why it does
tend to \( 0 \) for any fixed \( \sigma_a^2>0 \).
:::

### C. Going deeper

::: {#exr-mix-truncation}
[C1]

Let \( T \) be any estimator of \( \sigma_a^2\ge0 \) and \( T^+=\max(T,0) \). Show
that \( \E(T^+-\sigma_a^2)^2\le\E(T-\sigma_a^2)^2 \), with strict inequality whenever
\( \Pr(T<0)>0 \). Deduce that the truncated moment estimator dominates the unbiased
one in mean squared error, and explain why this does not make the unbiased estimator
useless.
:::

::: {.solution}
On \( \{T\ge0\} \) the two estimators agree. On \( \{T<0\} \),
\( |T^+-\sigma_a^2|=\sigma_a^2\le\sigma_a^2-T=|T-\sigma_a^2| \), because \( T<0 \) and
\( \sigma_a^2\ge0 \); the inequality is strict unless \( T=0 \) almost surely on that
event. Taking expectations gives the claim. The unbiased estimator remains useful
because unbiasedness is preserved by averaging: combining \( K \) independent studies,
\( K^{-1}\sum_kT_k \) is unbiased and consistent, whereas \( K^{-1}\sum_kT_k^+ \) is
biased upwards by an amount that does not vanish as \( K \) grows.
:::

::: {#exr-mix-two-way-random}
[C2]

Consider a balanced two-way layout with \( a \) levels of \( A \), \( b \) levels of
\( B \), \( m \) observations per cell, and *all* effects random:
\( Y_{ijk}=\mu+\alpha_i+\beta_j+(\alpha\beta)_{ij}+\varepsilon_{ijk} \) with the four
sets uncorrelated and variances
\( \sigma_A^2,\sigma_B^2,\sigma_{AB}^2,\sigma^2 \). Using the Kronecker projections of
[Section 16.3](../ch16-multiway-layouts/03-balanced.html) and @thm-rv-quadform-mean,
show that
\[
\E\,\text{MS}_A=\sigma^2+m\sigma_{AB}^2+bm\sigma_A^2,\qquad
\E\,\text{MS}_{AB}=\sigma^2+m\sigma_{AB}^2 .
\]
Which ratio tests \( \sigma_A^2=0 \)? Contrast this with the fixed-effects analysis,
where \( \text{MS}_A/\text{MSE} \) is the right ratio.
:::
