# Probit and complementary log–log links

Nothing forces a binary model to use the logit. Any function mapping the real line
onto \( (0,1) \) can serve as the response function \( \pi=h(\eta) \), and its inverse
is then the link. The candidates that survive in practice come with a story about how
the binary outcome was produced, because such a story says what the coefficients mean.
There are two such stories, and they give the probit and the complementary log–log
link.

## Thresholds

::: {#def-bin-latent}
[Threshold model]

Let \( \varepsilon_1,\dots,\varepsilon_n \) be independent with common continuous
distribution function \( F \), and let
\[
Y_i=\begin{cases}1,&\varepsilon_i\le\x_{(i)}\T\bbeta,\\ 0,&\text{otherwise.}\end{cases}
\]
The variable \( \x_{(i)}\T\bbeta-\varepsilon_i \) is the **latent response**, and \( Y_i \)
records the sign of it.
:::

Immediately \( \pi_i=\Pr(\varepsilon_i\le\x_{(i)}\T\bbeta)=F(\x_{(i)}\T\bbeta) \), so a
threshold model *is* a binary regression with response function \( F \) and link
\( F^{-1} \), and the content lies in which \( F \) the situation suggests.

::: {#prp-bin-links}
[Three links and where they come from]

In @def-bin-latent:

::: {.enumerate options="label=(\alph*)"}
1. if \( F \) is the standard logistic distribution function \( F(u)=e^{u}/(1+e^{u}) \),
   the model is @eq-bin-logistic; if \( F=\Phi \), it is the **probit model**
   \( \Phi^{-1}(\pi_i)=\x_{(i)}\T\bbeta \). Both \( F \) are symmetric about zero, so both
   links satisfy \( g(\pi)=-g(1-\pi) \);

2. if \( F(u)=1-\exp(-e^{u}) \), the distribution function of the smallest extreme
   value law, the model is the **complementary log–log** model
   \( \log\{-\log(1-\pi_i)\}=\x_{(i)}\T\bbeta \). This \( F \) is not symmetric: \( \pi \) leaves
   \( 0 \) slowly and approaches \( 1 \) abruptly. Reflecting it, that is, applying the
   complementary log–log link to \( 1-\pi \) and changing the sign of \( \bbeta \), gives the
   **log–log** model \( -\log\{-\log\pi_i\}=\x_{(i)}\T\bbeta \), with the two tails exchanged;

3. suppose a positive lifetime \( T \) has hazard function
   \( \lambda(t\mid\x)=\lambda_0(t)\exp(\x\T\bgamma) \), and that \( Y \) records whether the
   event has occurred by a fixed time \( t_0 \). Then \( Y \) follows a complementary
   log–log model with linear predictor \( \log\Lambda_0(t_0)+\x\T\bgamma \), where
   \( \Lambda_0 \) is the integrated baseline hazard. Moreover, if \( g \) is a continuous
   strictly increasing link with the property that for every \( c>0 \) replacing the hazard
   \( \lambda \) by \( c\lambda \) changes \( g(\pi) \) by a quantity \( \kappa(c) \) depending
   on \( c \) alone, and not on \( \Lambda_0(t_0) \) or \( \x \), then \( g \) is an increasing
   affine function of the complementary log–log link.
:::

:::

::: {.proof}
(a) and (b) are the identity \( \pi=F(\eta) \) read backwards, together with the
computation \( F^{-1}(\pi)=\log\{\pi/(1-\pi)\} \), \( \Phi^{-1}(\pi) \) and
\( \log\{-\log(1-\pi)\} \) for the three \( F \). For the reflection in (b), if
\( 1-\pi=1-\exp(-e^{-\eta}) \) then \( \log\{-\log\pi\}=-\eta \).

(c) The survivor function of a hazard model is
\( \Pr(T>t_0\mid\x)=\exp\{-\int_0^{t_0}\lambda(u\mid\x)\,du\}
=\exp\{-\Lambda_0(t_0)e^{\x\T\bgamma}\} \). Hence
\[
1-\pi=\exp\bigl(-e^{\log\Lambda_0(t_0)+\x\T\bgamma}\bigr),
\]
and taking \( \log(-\log\cdot) \) of both sides gives the claim. For uniqueness, write
\( H=\Lambda_0(t_0)e^{\x\T\bgamma} \) for the cumulative hazard and
\( G(u)=g\bigl(1-e^{-e^{u}}\bigr) \) with \( u=\log H \), so that \( u \) ranges over the
whole line as \( c \) and the baseline vary. The hypothesis reads
\( G(u+\log c)=G(u)+\kappa(c) \) for all real \( u \) and all \( c>0 \), so that
\( \kappa(v)=G(u+v)-G(u) \) does not depend on \( u \), and hence
\( \kappa(v+v')=\kappa(v)+\kappa(v') \) for all real \( v,v' \). Continuity of \( g \) makes
\( \kappa \) continuous, and the only continuous solutions of that additive Cauchy equation
are the linear ones, so \( \kappa(v)=av \); monotonicity gives \( a>0 \). Then
\( G(u)=G(0)+au \). Undoing the substitution,
\( g(\pi)=a\log\{-\log(1-\pi)\}+b \).
:::

Part (c) is why the complementary log–log link is standard for binary data derived
from times: a "yes" by a fixed follow-up date, a machine failed within a warranty
period. Its coefficients are log hazard ratios, whereas logistic coefficients for the
same data would be log odds ratios depending on
\( t_0 \) (@exr-bin-cloglog-interval).

## What the latent scale can and cannot tell you

::: {#prp-bin-latent-scale}
[Scale is not identified]

In @def-bin-latent, replace \( \varepsilon_i \) by \( \sigma\varepsilon_i+\tau \) with
\( \sigma>0 \), and \( \bbeta \) by \( \sigma\bbeta+\tau\mathbf{e}_1 \), where \( \mathbf{e}_1 \)
picks out an intercept column. The observed responses have exactly the same
distribution. Consequently \( \bbeta \) is identified only up to the location and scale of
the latent errors, while ratios \( \beta_j/\beta_k \) of coefficients of non-intercept
covariates are identified.
:::

::: {.proof}
The event \( \{\sigma\varepsilon_i+\tau\le\x_{(i)}\T(\sigma\bbeta+\tau\mathbf{e}_1)\} \)
is \( \{\varepsilon_i\le\x_{(i)}\T\bbeta\} \), because the intercept column is
identically \( 1 \) so \( \x_{(i)}\T\mathbf{e}_1=1 \) and the \( \tau \) cancels, after which
dividing by \( \sigma \) removes it. The two models therefore give the same
\( \pi_i \) for every \( i \). Ratios of non-intercept coefficients are unchanged by the
substitution.
:::

This is why the probit fixes \( \Var(\varepsilon)=1 \) and the logistic model fixes
\( \Var(\varepsilon)=\pi^2/3 \): a convention must be adopted. It is also why probit
and logit coefficients are not comparable until one is rescaled. Matching the standard
deviations multiplies a probit coefficient by \( \pi/\sqrt3=1.8138 \);
matching the slope at \( \pi=1/2 \) instead, where the logistic density is \( 1/4 \)
and the normal density is \( \phi(0) \), multiplies it by
\( 4\phi(0)=1.5958 \). The conventions differ by about fourteen per
cent, a fair measure of how meaningless the comparison is in the tails.

::: {.warning}
The latent-variable story tempts one into calling the probit coefficients "the real
effects, on the underlying continuous variable". They are not: the latent variable is
not observed and its scale is a convention that @prp-bin-latent-scale says the data
cannot pin down. What the story buys is a rationale for a particular \( F \).
:::

## How much does the choice matter?

Near \( \pi=1/2 \) all reasonable links agree, because smooth response curves look
alike near their steepest point. They separate in the tails: the logistic tails are
exponential, the normal tails decay like \( e^{-\eta^2/2} \), and the complementary
log–log is not symmetric at all.

::: {#exm-bin-three-links}
[The survey under three links]

Fitting the model of @exm-bin-anes with the three links gives the following
coefficients of party identification.

| link | \( \hat\beta \) | rescaled by \( 1.8138 \) | maximized log-likelihood |
|---|---|---|---|
| logit | 1.219 | — | \( -264.14 \) |
| probit | 0.673 | 1.221 | \( -265.96 \) |
| complementary log–log | 0.801 | — | \( -273.48 \) |

The rescaled probit coefficient reproduces the logistic one to two decimals, and the
largest difference between the two sets of fitted probabilities is
\( 0.0305 \). The complementary log–log fit is a different animal: its
fitted probabilities differ from the logistic ones by as much as
\( 0.1596 \), and with equal numbers of parameters AIC prefers the logit
by \( 18.67 \). Here the choice between symmetric links is immaterial and
the choice of an asymmetric one is not
([Figure 35.2.1](02-probit-cloglog.html#fig-bin-links)).
:::

::: {when-format="html"}
![**Figure 35.2.1.** (a) The three response curves, each rescaled to pass through
\( (0,1/2) \) with slope \( 1/4 \); the symmetric ones are close everywhere. (b)
Fitted probabilities for the survey of @exm-bin-anes: the probit differs from the
logit by at most \( 0.03 \), the complementary log–log by five times as
much.](link_comparison.svg){#fig-bin-links width=100%}
:::

::: {when-format="pdf"}
![(a) The three response curves, each rescaled to pass through
\( (0,1/2) \) with slope \( 1/4 \); the symmetric ones are close everywhere. (b)
Fitted probabilities for the survey of @exm-bin-anes: the probit differs from the
logit by at most \( 0.03 \), the complementary log–log by five times as
much.](link_comparison.pdf){width=100%}
:::

```{.python .run #cell-links-fits}
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
X = pd.DataFrame({"const": 1.0, "PID": anes["PID"], "age": anes["age"] / 10,
                  "educ": anes["educ"], "income": anes["income"]})
y = anes["vote"]

links = {"logit": sm.families.links.Logit(), "probit": sm.families.links.Probit(),
         "cloglog": sm.families.links.CLogLog()}
fits = {k: sm.GLM(y, X, family=sm.families.Binomial(link=v)).fit() for k, v in links.items()}

table = pd.DataFrame({k: f.params for k, f in fits.items()})
table["probit x 1.814"] = table["probit"] * np.pi / np.sqrt(3)
print(table.round(3).to_string())
print("maximized log-likelihoods:",
      {k: round(f.llf, 2) for k, f in fits.items()})
```

The choice therefore matters when the interesting probabilities are extreme, where
the links disagree and few observations settle the disagreement; when a mechanism
names the link, as in @prp-bin-links(c); and under the retrospective sampling of
[Section 35.5](05-case-control.html), where the logit is the only link whose
coefficients survive.

## Exercises

### A. Check your understanding

::: {#exr-bin-link-symmetry}
[A1]

Show that the logit and probit links satisfy \( g(\pi)=-g(1-\pi) \) while the
complementary log–log link does not. Deduce that for the first two, exchanging the
labels \( 0 \) and \( 1 \) merely reverses the sign of every coefficient, and say what
the exchange does under the complementary log–log link.
:::

::: {.solution}
For the logit, \( \log\{(1-\pi)/\pi\}=-\log\{\pi/(1-\pi)\} \); for the probit,
\( \Phi^{-1}(1-\pi)=-\Phi^{-1}(\pi) \). Relabelling therefore turns \( \eta \) into
\( -\eta \), so \( \bbeta \) becomes \( -\bbeta \) and the fit is the same model. For the
complementary log–log link, relabelling produces the *log–log* model of
@prp-bin-links(b): a genuinely different fit.
:::

### B. Practice

::: {#exr-bin-gumbel}
[B1]

Verify that \( F(u)=1-\exp(-e^{u}) \) is a distribution function, and that if
\( \varepsilon \) has this law then \( \E\varepsilon=-\gamma \) and
\( \Var\varepsilon=\pi^2/6 \), where \( \gamma=0.5772 \) is Euler's constant. Deduce
the factor that makes a complementary log–log coefficient's latent scale comparable
with a logistic one.
:::

::: {.solution}
\( F \) is continuous and increasing with limits \( 0 \) and \( 1 \). Substituting
\( v=e^{u} \) turns \( -\varepsilon \) into a standard Gumbel variable, whose mean is
\( \gamma \) and whose variance is \( \pi^2/6 \); reflecting gives the stated moments. The
standard deviation is \( \pi/\sqrt6=1.2825 \) against
\( \pi/\sqrt3=1.8138 \) for the logistic, so the multiplier is
\( \sqrt2 \). The location differs too, which shifts the intercept, so only the slopes
are comparable.
:::

::: {#exr-bin-cloglog-interval}
[B2]

Follow-up is divided into intervals \( (t_{k-1},t_k] \). Let \( Y_{ik}=1 \) if subject
\( i \) has the event in interval \( k \), given that it has not had it before. Show that a
proportional hazards model for the underlying continuous time implies
\( \log\{-\log(1-\Pr(Y_{ik}=1))\}=\alpha_k+\x_{(i)}\T\bgamma \) with
\( \alpha_k=\log\{\Lambda_0(t_k)-\Lambda_0(t_{k-1})\} \), and explain how to fit it with
ordinary binary regression software.
:::

### C. Going deeper

::: {#exr-bin-link-family}
[C1]

Consider the asymmetric family of Aranda-Ordaz (1981),
\( \pi=1-(1+\alpha e^{\eta})^{-1/\alpha} \) for \( \alpha>0 \). Show that
\( \alpha=1 \) gives the logistic response function and \( \alpha\to0 \) the complementary
log–log one. Explain how to test \( \alpha=1 \) with grouped data, and why the test has
little power with ungrouped data.
:::

