# Marginal versus conditional models

A clinician asks what a treatment would do to *this* patient's odds of recovery; a health
minister asks what it would do to the recovery rate of the *whole population*. In a linear
model the two questions have the same answer, because averaging a straight line over
patients returns the same straight line. As soon as the link is nonlinear they do not, and
no amount of data removes the difference: it is a property of the model, not of the sample.
This section sets up the two targets and proves that the population-averaged coefficient of
a logistic model is strictly smaller in magnitude than the individual one.

## Two coefficients, one model

Throughout this chapter the data are **clustered**: units \( i=1,\dots,m \) each carry
\( n_i \) observations \( Y_{i1},\dots,Y_{in_i} \), with \( n=\sum_i n_i \) — a patient
observed repeatedly, a school with its pupils, a firm across years. The simplest model that
makes observations within a cluster dependent gives each cluster its own level \( U_i \),
with the link \( g \) of [Chapter 34](../ch34-exponential-families-glm/index.html):
\[
g\bigl\{\E(Y_{ij}\mid U_i)\bigr\}=\alpha+\beta x_{ij}+U_i,
\qquad U_i\sim\Normal(0,\tau^2),
\]{#eq-gmm-random-intercept}

with \( Y_{i1},\dots,Y_{in_i} \) conditionally independent given \( U_i \). Here
\( \beta \) is a **conditional**, or *subject-specific*, coefficient: it compares two
covariate values *within one cluster*, with that cluster's level held wherever it happens to
be. Write it \( \beta^{C} \) when the distinction matters.

The other target concerns the response averaged over clusters. Write \( h=g^{-1} \) and
\[
m(x)=\E\bigl\{h(\alpha+\beta x+U)\bigr\},\qquad U\sim\Normal(0,\tau^2),
\]{#eq-gmm-marginal-curve}

for the **marginal mean function**: the mean response at \( x \) in a cluster drawn at
random. The **marginal**, or *population-averaged*, coefficients are the ones a model with
the same link would report if fitted to the population: the solution
\( (\alpha^{M},\beta^{M}) \) of
\[
\E\left[\begin{pmatrix}1\\X\end{pmatrix}
\Bigl\{m(X)-h\bigl(\alpha^{M}+\beta^{M}X\bigr)\Bigr\}\right]=\bzero ,
\]{#eq-gmm-marginal-parameter}

the expectation being over the covariate distribution, and assumed to have a unique
solution — which it does whenever \( \E\{(1,X)\T(1,X)\} \) is nonsingular and the marginal
mean function stays away from the ends of the range of \( h \). This is the population
version of the likelihood equations @eq-glm-score with the dispersion and the prior weights
dropped; for a canonical link the remaining factor \( h'(\eta)/V(\mu) \) is identically one,
so for the logit, log and identity links of this section \( (\alpha^{M},\beta^{M}) \) is
exactly what an ordinary generalized linear model fitted to the clustered data estimates
consistently — a point [Section 40.4](04-gee.html) returns to. For a binary covariate
@eq-gmm-marginal-parameter is solved by matching the two means exactly,
\( h(\alpha^{M})=m(0) \) and \( h(\alpha^{M}+\beta^{M})=m(1) \), so that \( \beta^{M} \) is
the log odds ratio between the two marginal probabilities.

::: {.warning}
Neither coefficient is "the" right one and neither is an approximation to the other. They
are answers to different questions, and a paper that reports one while discussing the other
has made a substantive error, not a numerical one.
:::

## The attenuation inequality

::: {#prp-gmm-marginal}
[Marginal and conditional coefficients]

Let \( U \) have mean zero and let \( h \) be the inverse link of
the model @eq-gmm-random-intercept.

::: {.enumerate options="label=(\alph*)"}
1. **Identity link.** If \( h(\eta)=\eta \) then \( m(x)=\alpha+\beta x \), so
   \( \alpha^{M}=\alpha \) and \( \beta^{M}=\beta \) for every distribution of \( U \)
   with mean zero. The two targets coincide.

2. **Log link.** If \( h(\eta)=e^{\eta} \) and \( \E e^{U}<\infty \), then
   \( m(x)=\E(e^{U})\,e^{\alpha+\beta x} \), so \( \beta^{M}=\beta \) and
   \( \alpha^{M}=\alpha+\log\E e^{U} \). For \( U\sim\Normal(0,\tau^2) \) the shift is
   \( \tau^{2}/2 \). Only the intercept moves.

3. **Probit link.** If \( h=\Phi \) and \( U\sim\Normal(0,\tau^{2}) \) then
   \[
   m(x)=\Phi\!\left(\frac{\alpha+\beta x}{\sqrt{1+\tau^{2}}}\right),
   \]{#eq-gmm-probit-exact}

   so the marginal mean function is again a probit, with
   \( \alpha^{M}=\alpha/\sqrt{1+\tau^{2}} \) and
   \( \beta^{M}=\beta/\sqrt{1+\tau^{2}} \).

4. **Logit link.** Let \( h(\eta)=e^{\eta}/(1+e^{\eta}) \), let \( U \) be any random
   variable with \( 0<\Var\{h(\theta+U)\}<\infty \) for every \( \theta \), and put
   \( P_\theta=h(\theta+U) \) and \( F(\theta)=g\bigl\{\E P_\theta\bigr\} \), the
   marginal linear predictor. Then
   \[
   F'(\theta)=1-\frac{\Var(P_\theta)}{\E(P_\theta)\,\E(1-P_\theta)}\in(0,1).
   \]{#eq-gmm-attenuation-derivative}

   Consequently, if \( \beta\ne0 \), then for any \( x_1<x_2 \) the marginal log odds
   ratio is strictly attenuated:
   \[
   0<\operatorname{sign}(\beta)\bigl\{g\,m(x_2)-g\,m(x_1)\bigr\}
   <\lvert\beta\rvert\,(x_2-x_1).
   \]
   In particular a binary covariate has \( 0<\lvert\beta^{M}\rvert<\lvert\beta^{C}\rvert \).
:::

:::

::: {.proof}
(a) \( \E(\alpha+\beta x+U)=\alpha+\beta x \), and this is exactly @eq-gmm-marginal-parameter
with \( h \) the identity.

(b) \( \E e^{\alpha+\beta x+U}=e^{\alpha+\beta x}\E e^{U} \). For a normal \( U \), @thm-mvn-mgf
gives \( \E e^{U}=e^{\tau^{2}/2} \). The right-hand side is again
of the form \( \exp(\alpha^{M}+\beta^{M}x) \), so it solves @eq-gmm-marginal-parameter
exactly.

(c) Let \( Z\sim\Normal(0,1) \) be independent of \( U \) and write
\( \eta=\alpha+\beta x \). Then
\[
\E\Phi(\eta+U)=\E\Pr(Z\le\eta+U\mid U)=\Pr(Z-U\le\eta),
\]
and \( Z-U\sim\Normal(0,1+\tau^{2}) \), which is @eq-gmm-probit-exact. Again the result is
of the assumed form, so it solves @eq-gmm-marginal-parameter.

(d) Write \( P_\theta=h(\theta+U) \) and \( \bar p(\theta)=\E P_\theta \). The logistic
inverse link satisfies \( h'=h(1-h) \), so differentiating under the expectation (legitimate
because \( 0\le h'\le 1/4 \), so the difference quotients are dominated) gives
\( \bar p'(\theta)=\E\{P_\theta(1-P_\theta)\} \). Since
\( F(\theta)=\log\bar p(\theta)-\log\{1-\bar p(\theta)\} \),
\[
F'(\theta)=\frac{\bar p'(\theta)}{\bar p(\theta)\{1-\bar p(\theta)\}}
=\frac{\E(P_\theta)-\E(P_\theta^{2})}{\E(P_\theta)\,\E(1-P_\theta)} .
\]
The denominator is \( \E(P_\theta)-\{\E(P_\theta)\}^{2} \), so numerator minus denominator
equals \( \{\E(P_\theta)\}^{2}-\E(P_\theta^{2})=-\Var(P_\theta) \), which gives @eq-gmm-attenuation-derivative.
The ratio is positive because \( P_\theta\in(0,1) \) makes
\( \E\{P_\theta(1-P_\theta)\}>0 \), and it is strictly less than one because
\( \Var(P_\theta)>0 \). Integrating \( F' \) from \( \alpha+\beta x_1 \) to
\( \alpha+\beta x_2 \), an interval of positive length because \( \beta\ne0 \), gives the
stated inequality; for a binary covariate the left-hand side is \( \beta^{M} \) by the
remark after @eq-gmm-marginal-parameter.
:::

Part (d) needs no normality of \( U \), no particular covariate and no asymptotics; it is
an identity about the logistic function. Its factor \( F'(\theta) \) is one minus the ratio
of the variance of the individual probabilities to the variance a single Bernoulli variable
with the average probability would have, so the more the clusters differ, the flatter the
population curve.

::: {.idea}
A nonlinear link and a random effect do not commute. Averaging \( h(\eta+U) \) over \( U \)
gives a *different function* of \( \eta \), and for the logit and probit it is a flatter
one. The linear and log links are the exceptions, because for them averaging changes at
most an additive or multiplicative constant.
:::

## How large is the difference?

The probit factor in @eq-gmm-probit-exact is exact; the logit has no closed form, but the
two links are so nearly proportional that a good approximation comes from pretending it
does. Matching a standard logistic distribution to a normal one scaled by
\( c=16\sqrt3/(15\pi)=0.5881 \) and applying @eq-gmm-probit-exact gives the classical
factor
\[
\beta^{M}\approx\frac{\beta^{C}}{\sqrt{1+c^{2}\tau^{2}}},
\qquad c=\frac{16\sqrt3}{15\pi},
\]{#eq-gmm-logit-factor}

due to Zeger, Liang and Albert (1988). It is an approximation, not a theorem.

::: {#exm-gmm-attenuation}
[Attenuation in a logistic random-intercept model]

Take \( \alpha=0 \), \( \beta^{C}=1 \) and a covariate \( X\sim\Normal(0,1) \), and solve
@eq-gmm-marginal-parameter numerically. With \( \tau=0.5 \) the marginal slope is
\( 0.9502 \); with \( \tau=1 \) it is \( 0.8379 \);
with \( \tau=2 \) it is \( 0.6142 \). The approximation @eq-gmm-logit-factor
gives \( 0.8620 \) and
\( 0.6477 \) in the last two cases, too large by about three and five per
cent. The probit factor at \( \tau=1 \) is \( 0.7071 \), far
smaller — but that is a different model, with a different meaning for \( \beta^{C} \). The
covariate distribution matters little: at \( \tau=1 \) a binary covariate gives a marginal
log odds ratio of \( 0.8318 \).
[Figure 40.1.1](01-marginal-versus-conditional.html#fig-gmm-attenuation)(a) shows why the
approximation works at all — the average of the conditional curves is almost exactly a
flatter logistic curve — and panel (b) tracks \( \beta^{M}/\beta^{C} \) as the clusters
become more heterogeneous.
:::

::: {when-format="html"}
![**Figure 40.1.1.** Attenuation in a random-intercept binary model. (a) With
\( \tau=1.5 \): five conditional logistic curves (grey), their average, and the logistic
curve solving @eq-gmm-marginal-parameter. (b) The ratio \( \beta^{M}/\beta^{C} \) against
\( \tau \) for the logit (exact and approximate) and for the
probit.](attenuation.svg){#fig-gmm-attenuation width=100%}
:::

::: {when-format="pdf"}
![Attenuation in a random-intercept binary model. (a) With
\( \tau=1.5 \): five conditional logistic curves (grey), their average, and the logistic
curve solving @eq-gmm-marginal-parameter. (b) The ratio \( \beta^{M}/\beta^{C} \) against
\( \tau \) for the logit (exact and approximate) and for the
probit.](attenuation.pdf){width=100%}
:::

```{.python .run #cell-attenuation-setup}
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import root
from scipy.stats import norm

def expit(z):
    return 1.0 / (1.0 + np.exp(-z))

# nodes and weights that integrate g(U) over U ~ N(0, 1): sum_k w_k g(z_k)
z_node, z_w = hermegauss(80)
z_w = z_w / np.sqrt(2 * np.pi)

def marginal_curve(eta, tau, link="logit"):
    """E[H(eta + U)] for U ~ N(0, tau^2), by Gauss-Hermite quadrature."""
    H = expit if link == "logit" else norm.cdf
    grid = np.atleast_1d(eta)[..., None] + tau * z_node
    return np.sum(z_w * H(grid), axis=-1)

def beta_marginal(alpha, beta, tau, link="logit"):
    """Population-averaged coefficients: the root of E[(1, x){m(x) - H(a + b x)}] = 0
    when x ~ N(0, 1), with m the marginal curve of the conditional model."""
    H = expit if link == "logit" else norm.cdf
    m = marginal_curve(alpha + beta * z_node, tau, link)

    def score(par):
        resid = m - H(par[0] + par[1] * z_node)
        return [np.sum(z_w * resid), np.sum(z_w * z_node * resid)]

    sol = root(score, [alpha, beta], tol=1e-13)
    return sol.x

C_LOGIT = 16 * np.sqrt(3) / (15 * np.pi)

def approx_factor(tau, link="logit"):
    c = C_LOGIT if link == "logit" else 1.0
    return 1.0 / np.sqrt(1.0 + c**2 * tau**2)

for t in [0.5, 1.0, 2.0]:
    print(f"tau = {t}:  logit {beta_marginal(0.0, 1.0, t)[1]:.4f}"
          f"   approximation {approx_factor(t):.4f}"
          f"   probit {1 / np.sqrt(1 + t**2):.4f}")
```

## Attenuation is not bias

Attenuation is not bias: both coefficients are correct answers to different questions. The
phenomenon is the **non-collapsibility** of the odds ratio, met in
[Section 35.5](../ch35-binary-responses/05-case-control.html) as @exm-bin-noncollapsible —
averaging \( 2\times2 \) tables that all share one odds ratio produces a table whose odds
ratio is closer to one, with no confounding whatever. Risk differences and risk ratios do
not behave this way, which is parts (a) and (b).

It follows that the conditional coefficient depends on what else is in the model: adding a
covariate that explains part of the cluster heterogeneity lowers \( \tau^{2} \) and, by
@eq-gmm-attenuation-derivative, moves the conditional and marginal coefficients of the
*other* covariates closer together. That is not evidence of confounding. Nor is either
coefficient a causal quantity by itself: randomization (@thm-cau-randomized) licenses a
causal reading of the marginal contrast between treated and untreated, but does not turn the
marginal log odds ratio into the conditional one, and does not license a conditional reading
of \( \beta^{C} \) unless the random effect really is a property of the cluster that the
intervention leaves alone.

Which target, then? The conditional coefficient answers questions about a change *within a
unit* — a dose for a patient, a policy for a school — and is comparable across studies with
different amounts of unmeasured heterogeneity, provided the random-effect distribution is
right. The marginal coefficient answers questions about a population total or rate, and is
safer when the cluster structure is a nuisance, because it assumes nothing about the cluster
effects. For the identity link the distinction disappears, which is why
[Chapter 32](../ch32-linear-mixed-models/index.html) never raised it: in @def-mix-model the
marginal model @eq-mix-marginal has the same \( \X\bbeta \) as the conditional one.

## Exercises

### A. Check your understanding

::: {#exr-gmm-two-targets}
[A1]

A random-intercept logistic model for whether a pupil passes an examination has
\( \hat\beta^{C}=0.60 \) for one extra hour of study per week and
\( \hat\tau=1.2 \). Use @eq-gmm-logit-factor to say roughly what a
population-averaged logistic model would report, and state in one sentence each what the
two numbers mean.
:::

::: {.solution}
The factor is \( (1+0.5881^{2}\times1.2^{2})^{-1/2} \). Now
\( 0.5881^{2}=0.3459 \) and \( 1.2^{2}=1.44 \), so the factor is
\( (1+0.4981)^{-1/2}=1.4981^{-1/2}\approx0.817 \), giving
\( \hat\beta^{M}\approx0.49 \). The first says: for a given pupil, an extra hour of study
multiplies the odds of passing by \( e^{0.60}\approx1.82 \). The second says: comparing two
randomly chosen pupils from the population who differ by an hour of study, the odds ratio is
about \( e^{0.49}\approx1.63 \).
:::

::: {#exr-gmm-log-link-check}
[A2]

A Poisson random-intercept model with a log link is fitted and reports
\( \hat\alpha=-1.2 \), \( \hat\beta=0.35 \), \( \hat\tau^{2}=0.8 \). What does a
population-averaged log-linear model estimate? Which of the three numbers changes?
:::

### B. Practice

::: {#exr-gmm-attenuation-bound}
[B1]

Use @eq-gmm-attenuation-derivative to show that in a logistic random-intercept model with
\( \Pr(Y=1\mid U)=h(U) \) and \( U \) symmetric about \( 0 \), the attenuation factor at
\( \theta=0 \) equals \( 1-4\Var\{h(U)\} \). Evaluate it for
\( U\sim\Normal(0,1) \), for which \( \Var\{h(U)\}=0.0434 \), and compare with
the ratio \( 0.8379 \) of @exm-gmm-attenuation.
:::

::: {.solution}
By symmetry \( \E h(U)=1/2 \), so the denominator of @eq-gmm-attenuation-derivative
is \( 1/4 \) and \( F'(0)=1-4\Var\{h(U)\} \). With
\( \Var\{h(U)\}=0.0434 \) this is \( 1-4\times0.0434=0.8264 \). This is the
*local* attenuation at \( \theta=0 \), and it is the strongest: as
\( \lvert\theta\rvert \) grows, \( P_\theta \) concentrates near \( 0 \) or \( 1 \),
its variance vanishes faster than \( \E P_\theta\E(1-P_\theta) \) and
\( F'\to1 \). The global ratio \( \beta^{M}/\beta^{C}=0.8379 \) of @exm-gmm-attenuation
is accordingly a little larger, because @eq-gmm-marginal-parameter
averages \( F' \) over the whole range of \( X \).
:::

::: {#exr-gmm-probit-general}
[B2]

Extend @prp-gmm-marginal(c) to a probit model with a random *slope*:
\( \Pr(Y=1\mid U_0,U_1)=\Phi(\alpha+\beta x+U_0+U_1 x) \) with
\( (U_0,U_1)\sim\Normal_2(\bzero,\G) \). Show that the marginal curve is
\( \Phi\{(\alpha+\beta x)/s(x)\} \) with
\( s(x)^2=1+(1,x)\G(1,x)\T \), and explain why the marginal model is now not a probit model
in \( x \) at all.
:::

::: {.solution}
Exactly as in the proof of (c), with \( Z-U_0-U_1x \) normal with mean zero and variance
\( 1+g_{00}+2g_{01}x+g_{11}x^{2}=s(x)^{2} \). The marginal linear predictor
\( (\alpha+\beta x)/s(x) \) is a ratio of a linear function to the square root of a
quadratic, which is linear in \( x \) only if \( g_{11}=g_{01}=0 \). With a random slope the
attenuation is different at different \( x \), so no single marginal coefficient describes
it.
:::

### C. Going deeper

::: {#exr-gmm-two-point}
[C1]

@eq-gmm-logit-factor is an approximation built for a normal random effect. This exercise computes the
attenuation exactly for a two-point one and shows how different the answer can be. Take a binary covariate
and the model \( \operatorname{logit}\Pr(Y=1\mid x,U)=\beta x+U \), where \( U=\pm c \) with probability
\( \tfrac12 \) each, \( c>0 \) and \( \beta>0 \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( m(0)=\tfrac12 \) exactly, so that \( \alpha^{M}=0 \) and \( \beta^{M} \) is the marginal log
   odds at \( x=1 \).

2. Show that
   \[
   \beta^{M}=\log\Bigl\{1+\frac{e^{2\beta}-1}{1+e^{\beta}\cosh c}\Bigr\}
   =\beta+\log\frac{\cosh c+e^{\beta}}{1+e^{\beta}\cosh c} .
   \]

3. Check that \( \beta^{M}=\beta \) at \( c=0 \), that \( \beta^{M} \) decreases in \( c \), and that
   \( 0<\beta^{M}<\beta \) for \( c>0 \), as @prp-gmm-marginal(d) requires.

4. Show that \( \beta^{M}\to0 \) as \( c\to\infty \), and that
   \( \beta^{M}\sim4\sinh(\beta)\,e^{-c} \). Since \( \Var(U)=c^{2} \), compare with the prediction of
   @eq-gmm-logit-factor at \( \tau=c \). What does the comparison say about that formula?
:::

:::

::: {.solution}
(a) \( m(0)=\tfrac12\{h(c)+h(-c)\}=\tfrac12 \), because \( h(-\eta)=1-h(\eta) \) for the logistic \( h \).
So \( g\,m(0)=0 \), and by the remark after @eq-gmm-marginal-parameter a binary covariate matches the two
marginal means exactly, giving \( \alpha^{M}=0 \) and \( \beta^{M}=g\,m(1) \).

(b) Write \( b=e^{\beta} \) and \( t=e^{c} \). Then
\( h(\beta+c)=bt/(1+bt) \) and \( h(\beta-c)=b/(t+b) \), so
\[
h(\beta+c)+h(\beta-c)=\frac{b\{t(t+b)+(1+bt)\}}{(1+bt)(t+b)}=\frac{b(t^{2}+1+2bt)}{t(1+b^{2})+b(1+t^{2})} .
\]
Hence \( 2m(1) \) is that expression, and
\( 2\{1-m(1)\}=\{2t(1+b^{2})+2b(1+t^{2})-b(t^{2}+1+2bt)\}/\{t(1+b^{2})+b(1+t^{2})\} \), whose numerator
simplifies to \( 2t+b+bt^{2} \). The marginal odds are therefore
\[
\frac{m(1)}{1-m(1)}=\frac{b(t^{2}+1+2bt)}{b(1+t^{2})+2t}
=\frac{b(\cosh c+b)}{1+b\cosh c},
\]
after dividing numerator and denominator by \( 2t \) and using \( (t^{2}+1)/(2t)=\cosh c \). Taking
logarithms gives the second form; and
\( b(\cosh c+b)/(1+b\cosh c)=1+(b^{2}-1)/(1+b\cosh c) \) gives the first.

(c) At \( c=0 \), \( \cosh c=1 \) and the first form is \( \log\{1+(b^{2}-1)/(1+b)\}=\log b=\beta \). The
first form is visibly decreasing in \( \cosh c \), hence in \( c \), because \( b^{2}>1 \); and it stays
positive for the same reason. This is @prp-gmm-marginal(d) in closed form.

(d) As \( c\to\infty \), \( (b^{2}-1)/(1+b\cosh c)\to0 \), so \( \beta^{M}\to0 \); and
\( \log(1+u)\sim u \) with \( \cosh c\sim e^{c}/2 \) gives
\( \beta^{M}\sim(b^{2}-1)/(b\,e^{c}/2)=2(b-b^{-1})e^{-c}=4\sinh(\beta)e^{-c} \). The approximation
@eq-gmm-logit-factor, applied at \( \tau=c \), gives
\( \beta/\sqrt{1+c^{2}\cdot0.5881^{2}} \), which falls only like \( 1/c \). At \( c=6 \), say, the exact
attenuation is of order \( e^{-6} \) while the formula still reports a third of \( \beta \). The formula is
not a general fact about attenuation: it is a calculation for a normal \( U \), where the tails of the
random effect keep some clusters at intermediate probabilities, and it fails for a \( U \) that pushes every
cluster to the flat ends of the logistic curve. Attenuation depends on the whole distribution of \( U \),
not on \( \tau \) alone.
:::

::: {#exr-gmm-slope-reversal}
[C2]

With a random *slope* the marginal coefficient need not even share the sign of the conditional one. Extend
@exr-gmm-probit-general to the logit: let \( x \) be binary and
\[
\operatorname{logit}\Pr(Y=1\mid x,U_1)=\alpha+\beta x+U_1x,\qquad
U_1=\pm d\ \text{with probability }\tfrac12 ,
\]
and take \( \alpha=\log3 \), \( \beta=\log2 \), \( d=\log6 \).

::: {.enumerate options="label=(\alph*)"}
1. Show that the two clusters have conditional log odds ratios \( \beta+d=\log12 \) and
   \( \beta-d=\log\tfrac13 \), whose mean is \( \beta=\log2>0 \).

2. Show that \( m(0)=\tfrac{111}{148} \) and \( m(1)=\tfrac{109}{148} \), so that the marginal odds ratio is
   \( \tfrac{109}{117}<1 \) and \( \beta^{M}<0 \).

3. Show that the marginal risk difference is \( -\tfrac1{74} \) and equals the average of the two
   cluster-specific risk differences \( \tfrac{33}{148} \) and \( -\tfrac{37}{148} \), while the marginal
   log odds ratio is *not* the average of the two conditional ones.

4. Explain why @prp-gmm-marginal(d) does not apply, and why "attenuation" is the wrong word for what has
   happened.
:::

:::

::: {.solution}
(a) Within a cluster, \( U_1 \) is fixed, so the log odds at \( x=1 \) minus the log odds at \( x=0 \) is
\( \beta+U_1 \), namely \( \log2+\log6=\log12 \) or \( \log2-\log6=\log\tfrac13 \). The mean over clusters
is \( \beta \) because \( \E U_1=0 \).

(b) At \( x=0 \) the random slope drops out, so \( m(0)=h(\alpha)=3/4=111/148 \). At \( x=1 \) the linear
predictors are \( \alpha+\beta+d=\log36 \) and \( \alpha+\beta-d=\log1=0 \), so
\( m(1)=\tfrac12(\tfrac{36}{37}+\tfrac12)=\tfrac12\cdot\tfrac{109}{74}=\tfrac{109}{148} \). The marginal
odds are \( (109/148)/(39/148)=109/39 \) at \( x=1 \) and \( 3 \) at \( x=0 \), so the marginal odds ratio is
\( 109/117 \), less than one.

(c) \( m(1)-m(0)=-2/148=-1/74 \). The cluster-specific risk differences are
\( \tfrac{36}{37}-\tfrac34=\tfrac{144-111}{148}=\tfrac{33}{148} \) and
\( \tfrac12-\tfrac34=-\tfrac{37}{148} \), whose average is \( (33-37)/296=-1/74 \). Equality holds because
\( m(1)=\E h(\alpha+\beta+U_1) \) is itself an average, so risk differences pass through it; log odds ratios
do not, and \( \log(109/117)\ne\tfrac12\{\log12+\log\tfrac13\}=\log2 \).

(d) Part (d) of the proposition concerns \( h(\theta+U) \), a random *level* added to a common linear
predictor, and its proof differentiates \( \theta\mapsto g\,\E h(\theta+U) \). Here the random term
multiplies \( x \), so the two conditional curves have different slopes and the argument does not start. The
conclusion fails too: \( \beta^{M} \) and \( \beta^{C} \) have opposite signs. Calling the difference
attenuation would suggest that the marginal coefficient is a shrunken version of the conditional one, and
here it is not a version of it at all — \( \beta^{C} \) is an average on the log odds scale, dominated by
the cluster in which the treatment helps a great deal, while \( m(1) \) is an average of probabilities, in
which the cluster where the treatment *hurts* costs more than the other gains. Two averages on two scales,
of a quantity that varies across clusters, and they need not agree even in sign; the warning after
@eq-gmm-marginal-parameter applies with full force.
:::
