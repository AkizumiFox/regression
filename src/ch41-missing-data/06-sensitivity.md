# Sensitivity analysis

Everything in the last three sections rested on missingness at random, and
@prp-mis-untestable says the data cannot check it. That is a reason neither to
despair nor to pretend, but to report how much the conclusion depends on the
assumption — which is what a sensitivity analysis does.

## Two factorizations

The selection factorization @eq-mis-selection splits the joint law into a model for
the data and a model for who responds. The other order gives the **pattern-mixture**
factorization
\[
f(\mathbf{D},\R\mid\boldsymbol{\uptheta}_1,\boldsymbol{\uptheta}_2)
=f(\mathbf{D}\mid\R,\boldsymbol{\uptheta}_1)\,f(\R\mid\boldsymbol{\uptheta}_2),
\]{#eq-mis-pattern-mixture}

a mixture, over patterns, of a different distribution of the data within each
pattern. The two describe the same object and differ only in what is easy to say.

::: {#prp-mis-sensitivity}
[Selection and pattern-mixture models]

Consider a single variable \( v \) recorded when \( R=1 \), alongside always
recorded variables \( \bu \).

::: {.enumerate options="label=(\alph*)"}
1. Every selection model \( f(v\mid\bu)\Pr(R=1\mid v,\bu) \) corresponds to the
   pattern-mixture model with
   \[
   \begin{aligned}
   f(v\mid\bu,R=r)&\propto f(v\mid\bu)\Pr(R=r\mid v,\bu),\\
   \Pr(R=r\mid\bu)&=\int f(v\mid\bu)\Pr(R=r\mid v,\bu)\,dv ,
   \end{aligned}
   \]
   and conversely, so neither formulation is more general.

2. The observed data identify \( f(v\mid\bu,R=1) \) and \( \Pr(R=1\mid\bu) \), and
   nothing about \( f(v\mid\bu,R=0) \). The marginal
   \[
   \begin{aligned}
   f(v\mid\bu)=\ &\Pr(R=1\mid\bu)\,f(v\mid\bu,R=1)\\
   +\ &\Pr(R=0\mid\bu)\,f(v\mid\bu,R=0)
   \end{aligned}
   \]{#eq-mis-pm-marginal}

   is therefore identified only after the analyst supplies the last factor. Missing
   at random is the supply \( f(v\mid\bu,R=0)=f(v\mid\bu,R=1) \).

3. Let the mechanism be \( \Pr(R=1\mid v,\bu)=h\{a(\bu)+cv\} \), with \( h \)
   continuous, strictly increasing and valued in \( (0,1) \) and \( a(\cdot) \)
   an unrestricted function. Then \( c \) is not identified: for every \( c \) at
   which the tilt @eq-mis-tilt below is integrable and its integral covers the
   required value, there are an \( a(\cdot) \) and an \( f(v\mid\bu,R=0) \)
   reproducing the observed-data distribution exactly. Any standard error reported
   for \( c \) comes from restricting \( a(\cdot) \) — to \( a+\bu\T\mathbf{b} \),
   say — and from the assumed forms of \( h \) and \( f(v\mid\bu) \), not from the
   data.
:::

:::

::: {.proof}
(a) Both expressions are the joint density \( f(v,R=r\mid\bu) \) written in the two
orders, by the definition of conditional probability; the displayed formulas are
Bayes's rule in each direction, and the correspondence is one-to-one wherever the
denominators are positive.

(b) The observed data consist of \( (\bu,R) \) for all units and \( v \) for those
with \( R=1 \). Its distribution is determined by the law of \( \bu \),
\( \Pr(R=1\mid\bu) \) and \( f(v\mid\bu,R=1) \), and by nothing else, exactly as
in @prp-mis-untestable. Equation @eq-mis-pm-marginal is the law of total probability.
Setting the unidentified factor equal to the identified one is the MAR choice,
since then \( \Pr(R=1\mid v,\bu) \) computed by Bayes's rule is free of \( v \).

(c) Fix \( c \), write \( \pi(\bu)=\Pr(R=1\mid\bu) \) and
\( f_1(v\mid\bu)=f(v\mid\bu,R=1) \) for the two identified pieces, and abbreviate
\( h_c(v)=h\{a(\bu)+cv\} \). Bayes's rule applied to (a) forces the unrecorded
conditional to be the recorded one *tilted* by the odds of not responding:
\[
f(v\mid\bu,R=0)=\frac{f_1(v\mid\bu)}{K(\bu)}\cdot\frac{1-h_c(v)}{h_c(v)},\qquad
K(\bu)=\int f_1\,\frac{1-h_c}{h_c}\,dv .
\]{#eq-mis-tilt}

Take that as the definition and assemble \( f(v\mid\bu) \) from @eq-mis-pm-marginal.
Then
\( f(v\mid\bu)h_c(v)=\pi f_1h_c+\{(1-\pi)/K\}f_1(1-h_c) \), which collapses to
\( \pi(\bu)f_1(v\mid\bu) \) — free of \( h_c \) — exactly when
\( K(\bu)=\{1-\pi(\bu)\}/\pi(\bu) \). Integrating over \( v \) then returns
\( \Pr(R=1\mid\bu)=\pi(\bu) \), and dividing returns
\( f(v\mid\bu,R=1)=f_1(v\mid\bu) \): both identified pieces are reproduced, for this
\( c \).

It remains to solve \( K(\bu)=\{1-\pi(\bu)\}/\pi(\bu) \) for the intercept
\( a(\bu) \), one scalar equation per \( \bu \). Since \( (1-h_c)/h_c \) decreases
pointwise in \( a(\bu) \), \( K \) is continuous and strictly decreasing, with
\( K\to0 \) as \( a(\bu)\to\infty \); a solution exists whenever the target lies in
its range. Both provisos bite. For \( h=\Phi \) and \( f_1 \) normal with standard
deviation \( s \), the integrand's logarithm behaves like
\( (c^2-s^{-2})v^2/2 \) as \( v\to-\infty \), so \( K \) is finite exactly when
\( |c|s<1 \), and then runs over all of \( (0,\infty) \); for larger \( |c| \) no
such pattern-mixture model exists and the construction is unavailable. Within the
range, the observed-data likelihood profiled over the remaining parameters is
constant in \( c \), so \( c \) is not identified. Restricting \( a(\bu) \) to
\( a+\bu\T\mathbf{b} \) leaves one equation per \( \bu \) and only
\( 1+\dim\bu \) unknowns; the profile is then no longer flat, and the curvature —
the reported standard error — comes from that restriction, not from the data.
:::

Part (c) is worth dwelling on. Fit a Heckman-style selection model — a normal outcome
equation and a probit selection equation, as in Heckman (1979) — and the software
returns an estimate of the selection coefficient with a standard error; part (c) says
that number is an artefact of the assumed functional forms. With an **exclusion
restriction**, a variable entering the selection equation and not the outcome
equation, the parameter becomes identified by more than functional form; without one
it does not. Reporting such a fit as though the data had spoken is the central
dishonesty of MNAR modelling, and the reason this section prefers to display the
untestable part rather than estimate it.

::: {.idea}
Under MNAR nothing is estimated that was not assumed. The useful question is not
"what is the answer?" but "how far must the assumption move before the answer
changes?"
:::

## Delta-adjustment and tipping points

The pattern-mixture form suggests an honest device. Impute under MAR, then shift
the imputed values by a fixed amount \( \delta \) before analysing them: that is,
assume
\[
\E(v\mid\bu,R=0)=\E(v\mid\bu,R=1)-\delta .
\]{#eq-mis-delta}

Here \( \delta \) is a *sensitivity parameter*: not estimated, but set. Zero is MAR,
and a positive \( \delta \) says the unrecorded values are worse than their recorded
counterparts with the same covariates. Because @eq-mis-delta constrains only the
unidentified factor of @eq-mis-pm-marginal, every \( \delta \) fits the observed
records equally well — in the example below all give the observed-data log-likelihood
\( -373.62 \) — so the choice is a statement about the world,
not about the data. Repeating the analysis over a grid of \( \delta \) and reporting
the curve is a **tipping-point analysis**: the value \( \delta^{*} \) at which the
conclusion reverses is the smallest departure from MAR that would overturn it.

::: {#exm-mis-tipping}
[How far must the assumption move?]

A two-arm study of \( 400 \) units records a baseline
measurement for everyone and an outcome for those who stay. Dropout is likelier in the
treated arm and at high baseline values, and \( 147 \) outcomes
are lost — \( 83 \) treated and
\( 64 \) control. Multiple imputation under MAR gives a
treatment effect of \( 0.475 \) with \( 95\% \)
interval \( (0.199,\ 0.750) \).

Now suppose the treated units who dropped out would have done worse than the MAR
imputation says, by \( \delta \). The estimate falls steadily and the lower limit
reaches zero at \( \delta^{*}=0.482 \), which is
\( 0.452 \) of a residual standard deviation. The finding
therefore survives the assumption that treated dropouts were up to about half a
standard deviation worse than comparable stayers, and not more. Whether that is
reassuring is a question for whoever knows why people dropped out.
:::

::: {when-format="html"}
![**Figure 41.6.1.** Estimated treatment effect and its \( 95\% \)
interval as the imputed outcomes of the treated arm are shifted down by
\( \delta \). Every \( \delta \) fits the recorded outcomes equally well. The
tipping point is where the lower limit reaches
zero.](tipping_point.svg){#fig-mis-tipping width=68%}
:::

::: {when-format="pdf"}
![Estimated treatment effect and its \( 95\% \)
interval as the imputed outcomes of the treated arm are shifted down by
\( \delta \). Every \( \delta \) fits the recorded outcomes equally well. The
tipping point is where the lower limit reaches
zero.](tipping_point.pdf){width=68%}
:::

```{.python .run #cell-sensitivity-data}
import numpy as np
from scipy import stats

rng = np.random.default_rng(4108)
n, tau_true = 400, 0.6
arm = np.repeat([0, 1], n // 2)                       # 0 = control, 1 = treated
base = rng.normal(size=n)                             # a baseline measurement, never missing
y = 0.5 + tau_true * arm + 0.6 * base + rng.normal(size=n)
eta = -1.0 + 0.7 * arm + 0.8 * base                   # dropout: depends on arm and baseline
seen = rng.uniform(size=n) > 1 / (1 + np.exp(-eta))
print(f"{n - seen.sum()} of {n} outcomes lost:"
      f" {int((~seen & (arm == 1)).sum())} treated, {int((~seen & (arm == 0)).sum())} control")
```

```{.python .run #cell-sensitivity-tipping}
Z = np.column_stack([np.ones(n), arm, base])


def mi_estimate(delta, M, rng):
    """Multiple imputation under MAR, then shift the treated arm's imputations by -delta."""
    Zo, k = Z[seen], int(seen.sum())
    g, *_ = np.linalg.lstsq(Zo, y[seen], rcond=None)
    sse = float(np.sum((y[seen] - Zo @ g) ** 2))
    qs, us = [], []
    for _ in range(M):
        s2 = sse / rng.chisquare(k - Z.shape[1])
        g_star = rng.multivariate_normal(g, s2 * np.linalg.inv(Zo.T @ Zo))
        draw = Z @ g_star + np.sqrt(s2) * rng.normal(size=n) - delta * arm
        y_full = np.where(seen, y, draw)
        b, *_ = np.linalg.lstsq(Z, y_full, rcond=None)
        resid = y_full - Z @ b
        v = (resid @ resid / (n - 3)) * np.linalg.inv(Z.T @ Z)[1, 1]
        qs.append(b[1])
        us.append(v)
    M_ = len(qs)
    q_bar, u_bar, b_var = np.mean(qs), np.mean(us), np.var(qs, ddof=1)
    total = u_bar + (1 + 1 / M_) * b_var
    gamma = (1 + 1 / M_) * b_var / total
    return q_bar, total, (M_ - 1) / gamma**2


deltas = np.linspace(0.0, 1.6, 33)
est, lo, hi = [], [], []
for d in deltas:
    q, t, nu = mi_estimate(d, 50, np.random.default_rng(4109))
    half = stats.t.ppf(0.975, nu) * np.sqrt(t)
    est.append(q)
    lo.append(q - half)
    hi.append(q + half)
est, lo, hi = map(np.array, (est, lo, hi))
tipping = float(np.interp(0.0, -lo, deltas))          # where the lower limit reaches zero
print(f"delta = 0: effect {est[0]:.3f}, 95% interval ({lo[0]:.3f}, {hi[0]:.3f})")
print(f"tipping point: delta = {tipping:.3f}")
```

Variations abound. The shift may be applied to one arm or both, on the response scale
or a link scale, and may depend on covariates or grow with time since dropout; in
trials it is often replaced by **jump to reference**, in which a treated dropout is
imputed as though from the control arm, or **copy increments in reference**, in which
only the subsequent changes are borrowed. Each is a particular choice of the
unidentified factor in @eq-mis-pm-marginal.

## Honest reporting

The technical content of this section is thin; the discipline is not. A report should
say how much is missing and where, by variable and by pattern; why, as far as the
reasons recorded at collection time reveal, those being the only evidence about the
mechanism that will ever exist; what was assumed, naming the variables MAR is
conditional on, since "missing values were imputed" states no assumption at all; what
the method was, including whether the imputation model contained the response; and
how far the conclusion travels, as a sensitivity curve or at least a tipping point.
That last item is where regulators have converged: the National Research Council
(2010) report recommends designing to prevent missingness and pre-specifying both the
primary analysis and the sensitivity analyses, rather than choosing them after seeing
the result.

::: {.warning}
A sensitivity analysis chosen after seeing the primary result is not a sensitivity
analysis; it is a selection procedure, with the properties
[Chapter 29](../ch29-model-selection/index.html) describes (@prp-sel-selection-bias).
Fix the grid of \( \delta \) and the summary to be reported before the analysis is
run.
:::

## Exercises

### A. Check your understanding

::: {#exr-mis-two-factorizations}
[A1]

Write the bivariate normal of @exr-mis-two-patterns in both factorizations,
selection and pattern-mixture, and identify which parameters are estimable from the
observed data in each. Why is the pattern-mixture form better suited to a
sensitivity analysis?
:::

::: {#exr-mis-tipping-units}
[A2]

In @exm-mis-tipping the tipping point is \( 0.482 \) on the
outcome scale. Express it as a multiple of the estimated treatment effect and as a
multiple of the residual standard deviation, and say which you would put in a
report and why.
:::

### B. Practice

::: {#exr-mis-delta-bias}
[B1]

For the model of @exm-mis-tipping, show analytically that the treatment effect
estimated after a \( \delta \)-shift of the treated arm's imputations is
\( \hat\tau(\delta)=\hat\tau(0)-\delta\,p_1 \), where \( p_1 \) is the fraction of
treated units whose outcome is missing, provided the analysis model contains the
arm indicator. Check the slope against the figure.
:::

::: {.solution}
With the arm indicator in the model and a completed data set, the estimated effect
is a difference of arm means adjusted for the baseline. Shifting the imputed
outcomes of the treated arm down by \( \delta \) lowers that arm's adjusted mean by
\( \delta \) times the fraction of its units that were imputed, and leaves the
control arm alone; the baseline adjustment is unaffected because the shift is
constant within the arm. Hence \( \hat\tau(\delta)=\hat\tau(0)-\delta p_1 \). With
\( 83 \) of \( 200 \) treated units missing,
\( p_1\approx0.42 \), and the slope of the curve in
[Figure 41.6.1](06-sensitivity.html#fig-mis-tipping) matches.
:::

### C. Going deeper

::: {#exr-mis-ccmv}
[C1]

For a monotone pattern with three time points, state the **complete-case missing
value** restriction — the unidentified conditional for a unit that dropped out
after time \( j \) equals the corresponding conditional among units seen at all
three times — and the **neighbouring-case** alternative. Show that both are
particular choices of the unidentified factor in @eq-mis-pm-marginal, and that the
MAR choice lies between them in a sense you should make precise.
:::

::: {#exr-mis-worst-case}
[C2]

Suppose the response is bounded in \( [0,1] \). Compute the **worst-case bounds**
for its mean: impute every missing value by \( 0 \), then by \( 1 \). Show that the
width of the resulting interval is the fraction missing, that it requires no
assumption at all, and that it is usually too wide to be useful. Discuss where such
bounds nevertheless earn their place.
:::

::: {.solution}
If a fraction \( \lambda \) of values is missing, the mean of the completed data
ranges over an interval of width \( \lambda \) as the imputations range over
\( [0,1] \), and every value in that interval is attained by some choice, so the
bounds are sharp and assumption-free — they are the Manski-style bounds for a
bounded outcome. They are useless when \( \lambda \) is large compared with the
effect of interest, which is the usual case. They earn their place as a first line
in a report: if the worst-case bound already excludes the null, no assumption about
the mechanism is needed at all, and that is worth knowing before any modelling
begins.
:::
