# Case–control sampling

A rare disease affects one person in a thousand, so a prospective study hoping to see
a hundred cases must follow a hundred thousand people, most of whom carry almost none
of the weight \( \pi_i(1-\pi_i) \). The economy is to go looking for the cases. In a
**retrospective**, or **case–control**, design the response decides who enters the
sample and the covariates are observed afterwards. This inverts the conditioning, and
one might expect nothing to be estimable. In fact, for the logit link and no other,
almost everything is.

## Only the intercept moves

::: {#thm-bin-casecontrol}
[Retrospective sampling under a logistic model]

Let \( (\,\X,Y) \) follow \( \operatorname{logit}\Pr(Y=1\mid\X=\x)=\alpha+\x\T\bbeta \) in the
population, and let \( S\in\{0,1\} \) indicate inclusion in the sample, with
\[
\Pr(S=1\mid Y=y,\X=\x)=\rho_y\in(0,1],\qquad y=0,1,
\]
depending on \( y \) but not on \( \x \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \operatorname{logit}\Pr(Y=1\mid\X=\x,S=1)
      =\alpha+\log(\rho_1/\rho_0)+\x\T\bbeta \). The slopes are unchanged; only the
   intercept moves, and it moves by a constant that does not involve \( \x \);

2. let \( F \) be a continuous strictly increasing distribution function, suppose
   \( \Pr(Y=1\mid\x)=F(\alpha+\x\T\bbeta) \), and suppose that for every ratio
   \( r=\rho_1/\rho_0 \) in some open interval containing \( 1 \) there is a constant \( c(r) \)
   with \( \Pr(Y=1\mid\x,S=1)=F\{\alpha+c(r)+\x\T\bbeta\} \) for all \( \x \). If the linear
   predictor \( \alpha+\x\T\bbeta \) ranges over an interval \( J \) with nonempty interior
   as \( \x \) varies, then \( F(u)=\{1+e^{-(au+b)}\}^{-1} \) for every \( u\in J \) and some
   \( a>0 \), \( b \): wherever the population identifies it, the link is logistic;

3. if \( \rho_1 \) and \( \rho_0 \) are known, or if the population prevalence
   \( \bar\pi=\Pr(Y=1) \) is known and the sample contains \( n_1 \) cases and \( n_0 \)
   controls, then \( \alpha \) itself is estimable: in the second case
   \( \hat\alpha=\hat\alpha_S-\log\{(n_1/n_0)(1-\bar\pi)/\bar\pi\} \), where \( \hat\alpha_S \)
   is the intercept fitted to the sample.
:::

:::

::: {.proof}
(a) By Bayes' theorem, writing \( \pi=\Pr(Y=1\mid\x) \),
\[
\Pr(Y=1\mid\x,S=1)=\frac{\rho_1\pi}{\rho_1\pi+\rho_0(1-\pi)},
\]
so that
\[
\frac{\Pr(Y=1\mid\x,S=1)}{\Pr(Y=0\mid\x,S=1)}
=\frac{\rho_1}{\rho_0}\cdot\frac{\pi}{1-\pi}.
\]
Taking logarithms and using \( \operatorname{logit}\pi=\alpha+\x\T\bbeta \) gives the
claim. The step needs \( \rho_y \) to be free of \( \x \), and nothing else.

(b) Put \( G=\operatorname{logit}\circ F \), a continuous strictly increasing function, and
write \( \eta=\alpha+\x\T\bbeta \). The computation in (a) used only Bayes' theorem, so it
holds for any \( F \):
\( \operatorname{logit}\Pr(Y=1\mid\x,S=1)=\log r+G(\eta) \). The hypothesis says this
equals \( G\{\eta+c(r)\} \). Hence
\[
G\{\eta+c(r)\}-G(\eta)=\log r\qquad\text{for all }\eta\text{ in an interval }J .
\]
Since \( G \) is strictly increasing, \( c \) is strictly increasing in \( r \) with
\( c(1)=0 \), and as \( r \) runs over an open interval containing \( 1 \), \( v=c(r) \) runs
over an open interval \( V \) containing \( 0 \). The displayed identity says that
\( \kappa(v)=G(\eta+v)-G(\eta) \) is the same for every \( \eta\in J \) with
\( \eta+v\in J \). Taking \( v,v' \) and \( \eta \) with \( v,v',v+v'\in V \) and
\( \eta,\eta+v,\eta+v+v' \) all in \( J \), which the nonempty interior of \( J \) permits,
gives \( \kappa(v+v')=\kappa(v)+\kappa(v') \); \( \kappa \) is continuous, so
\( \kappa(v)=av \) near \( 0 \), with \( a>0 \) by monotonicity. Fix \( \eta_0 \) in the
interior of \( J \). Then \( G(\eta)=G(\eta_0)+a(\eta-\eta_0) \) for \( \eta\in J \) within
\( V \) of \( \eta_0 \), and chaining such steps along the interval \( J \) extends the
identity to all of \( J \). Inverting, \( F(u)=\operatorname{logit}^{-1}(au+b) \) on \( J \):
nothing pins \( F \) down outside \( J \), and nothing can, since the sampling scheme says
nothing about linear predictors the population never produces.

(c) By (a) the fitted intercept estimates \( \alpha+\log(\rho_1/\rho_0) \), so knowing the
ratio recovers \( \alpha \). If instead \( \bar\pi \) is known and the population has \( N \)
members, then \( \rho_1=n_1/(N\bar\pi) \) and \( \rho_0=n_0/\{N(1-\bar\pi)\} \) in
expectation, so \( \rho_1/\rho_0=(n_1/n_0)(1-\bar\pi)/\bar\pi \).
:::

Part (a) is why the odds ratio dominates epidemiology. It is not that odds ratios are
more natural than risk ratios — they are less natural — but that they are the only
effect measure a retrospective study can estimate without knowing how the sample was
drawn. Part (b) says this belongs to the logit alone.

::: {.remark}
[Why the prospective likelihood is legitimate]

Part (a) says the *model* transfers, not that maximizing the ordinary prospective
likelihood on retrospective data is a correct procedure: the actual likelihood of such
a sample is \( \prod_if(\x_{(i)}\mid y_i) \), which involves the unknown population
distribution of \( \X \). Prentice and Pyke (1979) settled this. Treating that
distribution nonparametrically and profiling it out of the retrospective likelihood
yields exactly the prospective logistic likelihood, so the ordinary fit is the
semiparametric maximum likelihood estimate of \( \bbeta \) and the usual inverse
information is its correct asymptotic covariance. Their argument profiles an
infinite-dimensional nuisance parameter and is not reproduced here.
:::

::: {#exm-bin-case-control}
[A simulated retrospective study]

Take a population of \( 200000 \) with a continuous exposure
\( x_1\sim\Normal(0,1) \), an independent binary covariate \( x_2 \), and
\( \operatorname{logit}\pi=-5+1.2x_1+0.8x_2 \). The prevalence is
\( 0.0178 \), so there are \( 3565 \) cases, and
fitting the whole population gives \( \hat\alpha=-4.959 \),
\( \hat\beta_1=1.202 \), \( \hat\beta_2=0.748 \).

Now take every case and four random controls per case, and fit the same model to those
records as if they had been collected prospectively. Averaged over
\( 400 \) control samples the slopes are recovered, at
\( \hat\beta_1=1.201 \) and \( \hat\beta_2=0.755 \), while
the intercept averages \( -2.339 \) against the predicted
\( \hat\alpha+\log(\rho_1/\rho_0)=-4.959+2.623=-2.337 \).

With data generated from a *probit* model and fitted with a probit link, the
population slope is \( 0.547 \) while the retrospective samples average
\( 0.719 \), too large by a factor of \( 1.31 \)
([Figure 35.5.1](05-case-control.html#fig-bin-cc)).
:::

::: {when-format="html"}
![**Figure 35.5.1.** Retrospective estimates of a slope from \( 400 \)
case–control samples, with the whole-population estimate marked. (a) Logistic data
under the logit link: centred on the population value. (b) Probit data under the
probit link: systematically too large.](case_control.svg){#fig-bin-cc width=100%}
:::

::: {when-format="pdf"}
![Retrospective estimates of a slope from \( 400 \)
case–control samples, with the whole-population estimate marked. (a) Logistic data
under the logit link: centred on the population value. (b) Probit data under the
probit link: systematically too large.](case_control.pdf){width=100%}
:::

```{.python .run #cell-case-control-population}
import numpy as np
import statsmodels.api as sm

rng = np.random.default_rng(20)
N = 200_000
alpha, beta1, beta2 = -5.0, 1.2, 0.8
x1 = rng.normal(size=N)                                # a continuous exposure
x2 = (rng.random(N) < 0.3).astype(float)               # a binary covariate
Xpop = np.column_stack([np.ones(N), x1, x2])
ypop = (rng.random(N) < 1 / (1 + np.exp(-(Xpop @ [alpha, beta1, beta2])))).astype(float)

prospective = sm.GLM(ypop, Xpop, family=sm.families.Binomial()).fit()
cases = np.flatnonzero(ypop == 1)
controls = np.flatnonzero(ypop == 0)
n1 = len(cases)

print("prevalence", ypop.mean().round(4), " cases", n1)
print("whole population:", prospective.params.round(3))
```

```{.python .run #cell-case-control-retrospective}
def case_control(rng, ratio=4, link=None):
    """Take every case and `ratio` controls per case; fit the same model prospectively."""
    keep = np.r_[cases, rng.choice(controls, ratio * n1, replace=False)]
    family = sm.families.Binomial(link=link) if link else sm.families.Binomial()
    return sm.GLM(ypop[keep], Xpop[keep], family=family).fit().params


rng2 = np.random.default_rng(7)
draws = np.array([case_control(rng2) for _ in range(400)])
shift = np.log(1.0 / (4 * n1 / len(controls)))          # log(rho_1 / rho_0)

print("mean case-control estimate:", draws.mean(0).round(3))
print("predicted intercept:", round(prospective.params[0] + shift, 3))
print("slopes in the population:", prospective.params[1:].round(3))
```

## Matched designs

Controls are often matched to cases on variables that matter but are not of interest:
age, sex, neighbourhood. A matched set is one case and \( M \) controls sharing those
variables, and the model gives each set its own intercept:
\[
\operatorname{logit}\Pr(Y_{kj}=1\mid\x_{kj})=\alpha_k+\x_{kj}\T\bbeta,
\qquad k=1,\dots,K,\ \ j=0,1,\dots,M .
\]
Estimating the \( K \) intercepts along with \( \bbeta \) is the incidental-parameter
situation of @exm-opt-neyman-scott, in which maximum likelihood is not consistent.
Conditioning removes them.

::: {#prp-bin-conditional}
[Conditional logistic regression]

In the model above, suppose exactly one member of each set is a case, and condition on
that fact. Then:

::: {.enumerate options="label=(\alph*)"}
1. the conditional probability that member \( j^* \) of set \( k \) is the case, given the
   covariates of the set and that exactly one member is, is
   \[
\frac{\exp(\x_{kj^*}\T\bbeta)}{\sum_{j=0}^{M}\exp(\x_{kj}\T\bbeta)} ,
\]{#eq-bin-conditional}

   free of \( \alpha_k \);

2. the product of these over \( k \) is the **conditional likelihood**; maximizing it gives
   a consistent, asymptotically normal estimate of \( \bbeta \) as \( K\to\infty \) with \( M \)
   fixed, and the ordinary inverse-information standard errors apply to it;

3. for \( M=1 \) (matched pairs), let \( \mathbf{d}_k=\x_{k1}-\x_{k0} \) be the case-minus-control
   difference in a discordant pair. The conditional log-likelihood is
   \( -\sum_k\log\{1+\exp(-\bbeta\T\mathbf{d}_k)\} \), while the profile log-likelihood
   obtained by maximizing the *unconditional* likelihood over the \( \alpha_k \) is
   \( -2\sum_k\log\{1+\exp(-\tfrac12\bbeta\T\mathbf{d}_k)\} \). Consequently the
   unconditional estimate is exactly twice the conditional one, and converges to
   \( 2\bbeta \).
:::

:::

::: {.proof}
(a) With \( \pi_{kj}=\operatorname{logit}^{-1}(\alpha_k+\x_{kj}\T\bbeta) \), the probability
that \( j^* \) alone is a case is \( \pi_{kj^*}\prod_{j\ne j^*}(1-\pi_{kj}) \). Dividing by
the sum of the same expression over all choices, and cancelling
\( \prod_j(1-\pi_{kj}) \) from numerator and denominator, leaves
\( \{\pi_{kj^*}/(1-\pi_{kj^*})\}/\sum_j\{\pi_{kj}/(1-\pi_{kj})\} \). Each odds is
\( e^{\alpha_k}e^{\x_{kj}\T\bbeta} \), and \( e^{\alpha_k} \) cancels, giving @eq-bin-conditional.

(b) The sets are independent and @eq-bin-conditional is a genuine probability model
for a fixed parameter \( \bbeta \) of fixed dimension, so the standard likelihood
asymptotics apply as \( K\to\infty \). (It is the likelihood of a multinomial logit model,
which [Chapter 36](../ch36-multinomial-ordinal/index.html) studies in its own right.)

(c) Part (a) with \( M=1 \) gives
\( \{1+\exp(-\bbeta\T\mathbf{d}_k)\}^{-1} \), whose logarithm summed over \( k \) is the
stated conditional log-likelihood. For the unconditional calculation, write
\( u=\alpha_k+\tfrac12\bbeta\T(\x_{k1}+\x_{k0}) \) and
\( c=\tfrac12\bbeta\T\mathbf{d}_k \), so the two linear predictors are \( u+c \) and
\( u-c \). The pair's log-likelihood is
\( (u+c)-\log(1+e^{u+c})-\log(1+e^{u-c}) \); differentiating in \( u \) gives
\( 1-\operatorname{logit}^{-1}(u+c)-\operatorname{logit}^{-1}(u-c)=0 \), that is,
\( \operatorname{logit}^{-1}(u+c)=\operatorname{logit}^{-1}(c-u) \), whose unique
solution is \( u=0 \). Substituting \( u=0 \),
\[
c-\log(1+e^{c})-\log(1+e^{-c})=-2\log(1+e^{-c}),
\]
which is the stated profile contribution with \( c=\tfrac12\bbeta\T\mathbf{d}_k \). A
concordant pair contributes a supremum of \( 0 \), approached as \( \alpha_k\to\pm\infty \)
and free of \( \bbeta \). The profile log-likelihood is therefore twice the conditional
log-likelihood evaluated at \( \bbeta/2 \), so its maximizer is twice the conditional
maximizer; since the latter is consistent for \( \bbeta \), the former converges to
\( 2\bbeta \).
:::

::: {#exm-bin-matched}
[Doubling the effect]

Simulate \( 3000 \) matched pairs with one covariate, stratum intercepts
drawn from a normal distribution and \( \beta=1 \), keeping only the discordant pairs.
The conditional estimate is \( 0.942 \) with standard error
\( 0.037 \). Profiling the unconditional likelihood numerically over
the three thousand intercepts gives \( 1.884 \), twice that to three
decimals as @prp-bin-conditional(c) requires, and about twenty-four standard errors
from the truth.
:::

```{.python .run #cell-case-control-matched}
from scipy.optimize import minimize


def matched_pairs(K=3000, beta=1.0, seed=11):
    """K strata, each with one case and one control, and a stratum-specific intercept."""
    rng = np.random.default_rng(seed)
    a = rng.normal(-1.0, 2.0, K)                        # the nuisance intercepts
    xc, xk = np.empty(K), np.empty(K)
    for k in range(K):
        while True:                                     # sample the stratum until it is discordant
            x = rng.normal(size=2)
            pr = 1 / (1 + np.exp(-(a[k] + beta * x)))
            yy = rng.random(2) < pr
            if yy[0] != yy[1]:
                xc[k], xk[k] = (x[0], x[1]) if yy[0] else (x[1], x[0])
                break
    return xc, xk                                       # exposure of the case, of the control


x_case, x_control = matched_pairs()
# the conditional likelihood for 1:1 matching depends on beta only through the difference
d = x_case - x_control
cond = minimize(lambda b: np.sum(np.log1p(np.exp(-b[0] * d))), [0.0]).x[0]
print("conditional estimate:", round(cond, 3))
```

Matching buys precision when the matching variables are strong confounders, and costs
all information about them, since they are absorbed into the \( \alpha_k \).
Overmatching, matching on a variable on the causal path from exposure to outcome,
removes part of the effect one is measuring: a bad control (@thm-cau-bad-controls).

## Two warnings about odds ratios

::: {#exm-bin-noncollapsible}
[The odds ratio does not collapse]

Generate a population of \( 200000 \) with two *independent*
covariates, \( x_1 \) standard normal and \( x_2 \) Bernoulli with probability
\( 0.3 \), and \( \operatorname{logit}\pi=1.0\,x_1+4.0\,x_2 \), giving a prevalence of
\( 0.638 \). Fitting both covariates and an intercept recovers the coefficient
of \( x_1 \) as \( 1.006 \). Fitting an intercept and \( x_1 \) alone, and
thereby averaging over \( x_2 \), gives \( 0.745 \), a quarter smaller.
There is no confounding here: the two covariates are independent by construction, and
the marginal and conditional odds ratios simply differ.
:::

The odds ratio is not *collapsible*: the best logistic approximation to an average of
logistic curves is flatter than they are, by an amount growing with the variance of
the omitted term (@exr-bin-noncollapsible-normal). Conditional and marginal odds ratios
therefore answer different questions, and adding a covariate moves the other
coefficients away from zero even when it is independent of them, so the reasoning
"the coefficient changed, therefore the new variable is a confounder" is invalid here.
Risk differences and risk ratios do collapse.

The second warning is causal. @thm-bin-casecontrol requires selection to depend on
\( Y \) alone; if recruitment also depends on the exposure, say because the cases come
from a hospital serving an exposed population, then \( S \) is a collider on a path
between \( \X \) and \( Y \) and conditioning on \( S=1 \) opens
it (@prp-cau-selection-bias). And the conditional odds ratio it estimates is the causal
effect only if the covariates satisfy the backdoor criterion
(@def-cau-backdoor-criterion and @thm-cau-backdoor), which no sampling scheme can
supply.

::: {.idea}
Retrospective sampling is the one design in which the choice of link is forced:
@thm-bin-casecontrol makes the odds ratio estimable under a scheme that destroys every
other effect measure. The price is that odds ratios do not collapse, do not
approximate risk ratios unless the outcome is rare, and are causal only under the
usual assumptions.
:::

## Exercises

### A. Check your understanding

::: {#exr-bin-prevalence}
[A1]

A case–control study has \( 200 \) cases and \( 600 \) controls, and the disease affects
\( 2\% \) of the population. The fitted intercept is \( -0.8 \). Estimate the population
intercept, and the probability of disease at the reference covariate values.
:::

::: {.solution}
By @thm-bin-casecontrol(c),
\( \log\{(200/600)(0.98/0.02)\}=\log(16.333)=2.793 \), so
\( \hat\alpha=-0.8-2.793=-3.593 \) and the probability is \( 1/(1+e^{3.593})=0.0268 \),
close to the population prevalence as it should be.
:::

### B. Practice

::: {#exr-bin-selection-on-x}
[B1]

Suppose \( \Pr(S=1\mid Y=y,\X=\x)=\rho_y\,g(\x) \) for some \( g \) with values in
\( (0,1] \). Show that @thm-bin-casecontrol(a) still holds. Then suppose instead that
cases are sampled at rate \( \rho_1e^{\x\T\bgamma} \) while controls are sampled at a
rate free of \( \x \); find the bias in the estimated slopes.
:::

::: {#exr-bin-matched-2x2}
[B2]

For matched pairs with a single binary exposure, let \( n_{10} \) be the number of pairs
with an exposed case and unexposed control and \( n_{01} \) the reverse. Show that the
conditional estimate of the log odds ratio is \( \log(n_{10}/n_{01}) \), and that the
conditional test of \( \beta=0 \) is a binomial test of \( 1/2 \) in \( n_{10}+n_{01} \)
trials.
:::

::: {.solution}
With a binary exposure \( d_k\in\{-1,0,1\} \), and pairs with \( d_k=0 \) contribute a
constant. The conditional log-likelihood is
\( -n_{10}\log(1+e^{-\beta})-n_{01}\log(1+e^{\beta}) \); setting the derivative to zero
gives \( e^{\beta}=n_{10}/n_{01} \). Writing \( p=e^{\beta}/(1+e^{\beta}) \), the conditional
likelihood is that of \( n_{10} \) successes in \( n_{10}+n_{01} \) Bernoulli trials, and
\( \beta=0 \) is \( p=1/2 \). The resulting test is McNemar's, its exact version the
binomial test.
:::

::: {#exr-bin-efficiency}
[B3]

With \( n_1 \) cases fixed and \( n_0=Mn_1 \) controls, show that for a single binary
exposure and a rare disease the asymptotic variance of the log odds ratio is
proportional to \( 1+1/M \), and deduce that there is little gain beyond \( M=4 \) or
\( 5 \).
:::

### C. Going deeper

::: {#exr-bin-noncollapsible-normal}
[C1]

Let \( \operatorname{logit}\Pr(Y=1\mid X,Z)=\beta X+\gamma Z \) with \( Z \) independent
of \( X \). Show that \( \Pr(Y=1\mid X) \) is not of the form
\( \operatorname{logit}^{-1}(\beta^*X) \), and explain why the best logistic
approximation has \( |\beta^*|<|\beta| \). What happens as \( \gamma\to0 \)?
:::

::: {.solution}
\( \Pr(Y=1\mid X=x)=\E_Z\operatorname{logit}^{-1}(\beta x+\gamma Z) \) is a mixture of
horizontally shifted logistic curves, and a nondegenerate such mixture is not itself
logistic: its derivative at its median is strictly smaller than that of a single
component with the same slope, by concavity of the response function above the median
and convexity below. The marginal curve is flatter, so the best logistic fit has a
smaller slope in modulus. As \( \gamma\to0 \), \( \beta^*\to\beta \); the attenuation
grows with \( \gamma^2\Var Z \), which is why the effect in @exm-bin-noncollapsible is
large.
:::
