# Missingness mechanisms

A data set with holes carries two kinds of information: the values that were
recorded, and the pattern of which values were recorded at all. The second decides
whether any analysis of the first means what it appears to. This section gives the
pattern a probability model.

## Notation

Write the **complete data** as an array \( \mathbf{D} \) with \( n \) rows and
\( q \) columns, whose \( i \)th row \( \mathbf{d}_{(i)} \) holds the variables of
unit \( i \). Response and covariates sit in the same array: which is the response
is a question for the analysis, not for the mechanism. Alongside it put the
**missingness indicator** array \( \R \) of the same shape, with \( R_{ij}=1 \)
when \( D_{ij} \) is recorded; \( \R \) is free for this use here, as \( \bH \)
was in Chapters 5 and 10.

The array \( \R \) is always completely observed, which is what makes it data. Split
\( \mathbf{D} \) into recorded entries \( \mathbf{D}_{\mathrm{obs}} \) and unrecorded
ones \( \mathbf{D}_{\mathrm{mis}} \). The analyst holds
\( (\mathbf{D}_{\mathrm{obs}},\R) \); every model here describes \( (\mathbf{D},\R) \)
by the factorization
\[
f(\mathbf{D},\R\mid\boldsymbol{\uptheta},\boldsymbol{\uppsi})
 =f(\mathbf{D}\mid\boldsymbol{\uptheta})\;f(\R\mid\mathbf{D},\boldsymbol{\uppsi}) .
\]{#eq-mis-selection}

Here \( \boldsymbol{\uptheta} \) is the parameter the analyst cares about — a
regression coefficient lives in it — and the second factor is the **missingness
mechanism**. Conditioning it on the *complete* data is the essential move: it lets
the probability that a value is recorded depend on the value itself, which is exactly
the case that causes trouble. This is the **selection** factorization;
[Section 41.6](06-sensitivity.html) uses the other one.

## Three classes

::: {#def-mis-mechanisms}
[Missing completely at random, at random, not at random]

Let \( \mathbf{D} \) and \( \R \) have the joint distribution @eq-mis-selection, and
fix the realized values \( \R=\mathbf{r} \) and
\( \mathbf{D}_{\mathrm{obs}}=\mathbf{a} \). The mechanism is

::: {.enumerate options="label=(\alph*)"}
1. **missing completely at random** (MCAR) if
   \( f(\mathbf{r}\mid\mathbf{D},\boldsymbol{\uppsi})=f(\mathbf{r}\mid\boldsymbol{\uppsi}) \)
   for every \( \mathbf{D} \) and \( \boldsymbol{\uppsi} \): the pattern is
   independent of the data, recorded or not;

2. **missing at random** (MAR) if
   \( f(\mathbf{r}\mid\mathbf{D},\boldsymbol{\uppsi}) \) takes the same value for
   every \( \mathbf{D} \) whose observed part is \( \mathbf{a} \), and every
   \( \boldsymbol{\uppsi} \), so that it may be written
   \( f(\mathbf{r}\mid\mathbf{a},\boldsymbol{\uppsi}) \);

3. **missing not at random** (MNAR) otherwise.
:::

:::

MCAR implies MAR, and the two buy different things. Under MCAR the recorded units are
a random subsample, so *every* analysis of them estimates what it would have estimated
with no losses, at the price of a smaller sample. MAR says less: the recorded units
are a biased subsample, but the bias is a function of recorded variables, so a method
using those variables correctly can undo it.

::: {.warning}
The names have misled generations of readers. "Missing at random" does not mean the
holes are haphazard; it means they are haphazard *within each configuration of the
recorded values* — a strong assumption in modest language. "Missing completely at
random" is the haphazard one.
:::

Definition @def-mis-mechanisms states MAR at the realized \( (\mathbf{r},\mathbf{a}) \);
requiring it at every pattern is stronger, and the distinction is not pedantry. A
Bayesian, conditioning on the data in hand, needs only the realized version; a
frequentist, averaging over data sets that did not occur, needs the stronger one
(Seaman, Galati, Jackson and Carlin, 2013). This chapter states the realized version
and flags where the stronger one is wanted.

## What each class looks like

::: {#exm-mis-survey}
[A survey that asks about income]

A survey records age \( x \) for everyone and income \( y \) for those who answer.
Three stories. *Completely at random:* a printing fault left the income question off
one batch of forms, assembled without regard to anything, so
\( \Pr(R=0\mid x,y) \) is constant. *At random:* younger respondents skip the
question more often, and among respondents of a given age, answering has nothing to
do with earnings, so \( \Pr(R=0\mid x,y)=\Pr(R=0\mid x) \). *Not at random:* high
earners skip the question at every age, so \( \Pr(R=0\mid x,y) \) depends on
\( y \) even given \( x \). The three are mutually exclusive descriptions of the
world, and @prp-mis-untestable says the observed data cannot tell them apart.
:::

So "the mechanism" is a property not of a variable but of the pair (variable, recorded
variables). If the survey also recorded occupation, and high earners skip only because
they hold occupations whose members skip, the third story becomes MAR once occupation
joins \( \mathbf{D}_{\mathrm{obs}} \). Enlarging the recorded data can turn MNAR into
MAR; it can never do the reverse.

::: {#exm-mis-co2}
[A record with real gaps]

The weekly averages of the continuous Mauna Loa record of atmospheric carbon dioxide
run from March 1958 to December 2001 and cover \( 2284 \) weeks, of which
\( 59 \) carry no value — \( 2.6\% \) of the
record, in \( 22 \) runs, the longest of
\( 18 \) consecutive weeks.

They are not spread evenly. Before 1970, \( 8.6\% \) of weeks are
missing; from 1970 on, \( 0.4\% \), and a logistic regression of
\( 1-R_i \) on calendar time multiplies the odds of a gap by
\( 0.192 \) per decade. The mechanism depends on a variable
recorded for every week, so the record is not MCAR.

Whether it is MAR given time is a question about the instrument, not about the
numbers. The published rule is that a day contributes only if the analyser ran
steadily for at least six hours, and a week with no such day contributes nothing. So
the question is whether the analyser's failure to deliver six steady hours is related
to the concentration it would have read — plausibly not, if the failures are
electrical or logistical; plausibly so, if unsettled air both disturbs the instrument
and moves the reading. Nothing in
[Figure 41.1.1](01-mechanisms.html#fig-mis-co2) can settle that.
:::

::: {when-format="html"}
![**Figure 41.1.1.** The weekly Mauna Loa record. (a) The series, with a vertical
rule at each missing week. (b) The fraction of weeks missing in each calendar year.
The mechanism depends strongly on an observed variable, calendar time, and so is not
missing completely at random.](co2_gaps.svg){#fig-mis-co2 width=100%}
:::

::: {when-format="pdf"}
![The weekly Mauna Loa record. (a) The series, with a vertical
rule at each missing week. (b) The fraction of weeks missing in each calendar year.
The mechanism depends strongly on an observed variable, calendar time, and so is not
missing completely at random.](co2_gaps.pdf){width=100%}
:::

```{.python .run #cell-mechanisms-pattern}
import numpy as np
import statsmodels.api as sm

co2 = sm.datasets.co2.load_pandas().data["co2"]       # one value per week, some absent
t = co2.index.year + (co2.index.dayofyear - 0.5) / 365.25
r = (~co2.isna()).to_numpy().astype(float)            # 1 = recorded, 0 = missing
n, n_mis = len(r), int((1 - r).sum())
print(f"{n} weeks, {n_mis} missing ({100 * n_mis / n:.1f}%)")

era = (t >= 1970).astype(float)                       # one observed covariate: is it 1970 or later?
rate_early = 1 - r[era == 0].mean()
rate_late = 1 - r[era == 1].mean()
print(f"missing before 1970: {100 * rate_early:.1f}%   from 1970 on: {100 * rate_late:.1f}%")
```

One feature of \( \R \) matters for computation rather than validity. A pattern is
**monotone** when the columns can be ordered so that \( R_{ij}=0 \) forces
\( R_{ik}=0 \) for \( k>j \): once a unit drops out it never returns, as in
longitudinal dropout. Such patterns are far easier to handle, the joint density
factoring into a chain of conditionals each fitted on a nested
subsample (@exr-mis-monotone).

## Why the class cannot be read off the data

::: {#prp-mis-untestable}
[The observed data do not identify the mechanism]

Let \( \mathbf{d}=(d_1,d_2) \) be a unit's record, with \( d_1 \) always recorded
and \( d_2 \) recorded exactly when \( R=1 \). The distribution of the observed data
is determined by the three pieces \( f(d_1) \), \( \Pr(R=1\mid d_1) \) and
\( f(d_2\mid d_1,R=1) \), and determines nothing else. Consequently:

::: {.enumerate options="label=(\alph*)"}
1. for every conditional density \( g(d_2\mid d_1) \) there is exactly one joint
   distribution of \( (\mathbf{d},R) \) with the given observed-data distribution
   and with \( f(d_2\mid d_1,R=0)=g(d_2\mid d_1) \);

2. one of these choices, \( g=f(\cdot\mid d_1,R=1) \), is MAR; every other is MNAR;
   and all assign the observed data the same likelihood.
:::

:::

::: {.proof}
The observed data are \( (d_1,R) \), with \( d_2 \) when \( R=1 \). Their density
is \( f(d_1)\Pr(R=0\mid d_1) \) on \( R=0 \) and
\( f(d_1)\Pr(R=1\mid d_1)f(d_2\mid d_1,R=1) \) on \( R=1 \); only the three listed
pieces enter.

(a) Given the three pieces and a candidate \( g \), define
\[
\begin{aligned}
f(d_1,d_2,R=1)&=f(d_1)\Pr(R=1\mid d_1)\,f(d_2\mid d_1,R=1),\\
f(d_1,d_2,R=0)&=f(d_1)\Pr(R=0\mid d_1)\,g(d_2\mid d_1).
\end{aligned}
\]
This integrates to one, has the given observed-data distribution, and has
conditional density \( g \) for \( d_2 \) given \( (d_1,R=0) \). Conversely any joint
distribution with those properties has these two densities, so it is unique.

(b) Abbreviate \( p_r=\Pr(R=r\mid d_1) \) and \( f_1=f(d_2\mid d_1,R=1) \). Under the
constructed distribution,
\[
\Pr(R=1\mid d_1,d_2)=\frac{p_1f_1}{p_1f_1+p_0\,g(d_2\mid d_1)} ,
\]
which is free of \( d_2 \) for all \( d_1 \) precisely when \( g \) equals
\( f(\cdot\mid d_1,R=1) \). Every member of the family assigns the observed data the
same likelihood by construction, so no statistic computed from those data can prefer
one to another.
:::

::: {.idea}
The missing values are missing. No cleverness recovers the conditional distribution
of what was not recorded given what was; a choice of that distribution *is* the
assumption, and different choices fit the observed data identically. MAR is one such
choice: convenient, often defensible, never verifiable.
:::

Two consequences. A test of MCAR against MAR is possible, that comparison involving
recorded variables only, and Little (1988) gives one; a test of MAR against MNAR is
not. And a sensitivity analysis is not an optional extra but the only way to report
what the data do not determine ([Section 41.6](06-sensitivity.html)).

Finally, the classes of @def-mis-mechanisms describe the world; whether a procedure
works is a property of the *pair*. Under MCAR nearly every method is valid, including
dropping incomplete records. Under MAR with distinct parameters, likelihood and
Bayesian inference with a correct model for \( \mathbf{D} \) are valid
by @thm-mis-ignorable, while dropping incomplete records usually is
not (@prp-mis-complete-case). And under a mechanism depending on a regression's
covariates but not on its response, dropping incomplete records is valid *for that
regression*, even though the recorded units are a badly biased sample for any other
purpose. That third case is the most useful and least well known, and the next
section starts there.

## Exercises

### A. Check your understanding

::: {#exr-mis-classify}
[A1]

Classify each mechanism as MCAR, MAR or MNAR with respect to the recorded variables
named.

::: {.enumerate options="label=(\alph*)"}
1. A blood assay is run on a random one-third of the samples chosen by a random
   number generator; age, sex and diagnosis are recorded for all.

2. In a weight-loss trial, participants who have gained weight since the last visit
   stop attending. Weight at every previous visit is recorded.

3. A scale reads "over range" and records nothing whenever the true weight exceeds
   \( 150 \) kg.

4. A questionnaire item is missing for everyone interviewed by one interviewer,
   whose identity is recorded.
:::

:::

::: {.solution}
(a) MCAR. (b) MNAR for the current visit's weight, since whether it is recorded
depends on that weight itself; it would be MAR if attendance depended only on the
previous, recorded visits. (c) MNAR, and of an extreme kind: missingness is a
deterministic function of the unrecorded value. (d) MAR given the interviewer, and
not MCAR; without the interviewer's identity recorded it would be MNAR *given the
recorded variables*, which is what @def-mis-mechanisms measures.
:::

::: {#exr-mis-mcar-implies-mar}
[A2]

Show that MCAR implies MAR, and explain why no enlargement of
\( \mathbf{D}_{\mathrm{obs}} \) can turn a MAR mechanism into an MNAR one.
:::

### B. Practice

::: {#exr-mis-monotone}
[B1]

Let \( \mathbf{d}=(d_1,\dots,d_q) \) have a monotone pattern, so that the recorded
entries are \( d_1,\dots,d_K \) for some \( K\in\{0,\dots,q\} \). Write the joint
density of \( (\mathbf{d},K) \) using the chain rule
\( f(\mathbf{d})=\prod_jf(d_j\mid d_1,\dots,d_{j-1}) \), and show that MAR is
equivalent to \( \Pr(K=k\mid\mathbf{d}) \) depending on \( \mathbf{d} \) only through
\( d_1,\dots,d_k \).
:::

::: {.solution}
Definition @def-mis-mechanisms asks that \( \Pr(K=k\mid\mathbf{d}) \) be the same for
all \( \mathbf{d} \) agreeing in their first \( k \) coordinates, which is the stated
condition. Multiplying the chain-rule factorization by it gives
\[
f(\mathbf{d},K=k)=\Bigl\{\prod_{j\le k}f(d_j\mid d_{<j})\Bigr\}\Pr(K=k\mid d_1,\dots,d_k)
\Bigl\{\prod_{j>k}f(d_j\mid d_{<j})\Bigr\},
\]
and integrating out \( d_{k+1},\dots,d_q \) removes the last brace entirely. That is
the factorization [Section 41.3](03-likelihood.html) exploits.
:::

::: {#exr-mis-little-test}
[B2]

Two variables \( d_1 \) (always recorded) and \( d_2 \) (recorded when \( R=1 \)).
Propose a statistic with a known null distribution under MCAR and power against "the
mechanism depends on \( d_1 \)", and explain why no statistic can have power against
"the mechanism depends on \( d_2 \) given \( d_1 \)".
:::

::: {#exr-mis-two-patterns}
[B3]

Let \( \mathbf{d}=(d_1,d_2) \) be bivariate normal with mean \( (\mu_1,\mu_2) \) and
covariance \( \bSigma \), and let \( d_2 \) be recorded with probability
\( \Phi(a+bd_1) \). Show that the mechanism is MAR, compute \( \E(d_2\mid R=1) \),
and confirm that it differs from \( \mu_2 \) unless \( b=0 \) or
\( \sigma_{12}=0 \). Which class does the mechanism fall in if the probability is
instead \( \Phi(a+bd_2) \)?
:::

::: {.solution}
The probability depends only on the recorded \( d_1 \), so the mechanism is MAR.
Writing \( d_2=\mu_2+\beta(d_1-\mu_1)+e \) with \( \beta=\sigma_{12}/\sigma_{11} \)
and \( e \) independent of \( d_1 \),
\( \E(d_2\mid R=1)=\mu_2+\beta\{\E(d_1\mid R=1)-\mu_1\} \), and
\( \E(d_1\mid R=1)\ne\mu_1 \) whenever \( b\ne0 \), since the selection probability
is monotone in \( d_1 \). With \( \Phi(a+bd_2) \) the mechanism depends on the
possibly unrecorded \( d_2 \), so it is MNAR for \( b\ne0 \).
:::

### C. Going deeper

::: {#exr-mis-realized-mar}
[C1]

Let the units be independent, each with \( d_1\in\{0,1\} \) always recorded, and let
\( d_2 \) be recorded with probability \( p_0 \) when \( d_1=0 \) and with
probability \( q(d_2) \) when \( d_1=1 \), where \( q \) is not constant. Show that
@def-mis-mechanisms(b) fails at some patterns, but that on any sample in which no
unit has \( d_1=1 \) the realized-value condition holds exactly. What does this say
about a frequentist argument, which averages over samples that did not occur,
compared with a Bayesian one, which conditions on the sample in hand?
:::
