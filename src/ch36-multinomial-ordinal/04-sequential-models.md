# Sequential models

The cumulative models of the last two sections treat an ordered response as a coarsened
continuum. Some ordered responses are not like that. A student reaches a master's degree
only by passing through a bachelor's; an unemployment spell reaches its second year only
by surviving its first. There the categories are *stages*, reached in order, and each
stage is a decision taken by those who got that far. A model built on that picture asks,
for each stage, a question conditional on having reached it.

::: {#def-mlt-sequential}
[Sequential model]

With \( \Y_i\sim\text{Mult}(m_i,\boldsymbol{\uppi}_i) \) as in @def-mlt-baseline, let
\[
\delta_{ir}=P(Y_i=r\mid Y_i\ge r,\ \x_{(i)}),\qquad r=1,\dots,c-1,
\]{#eq-mlt-continuation}

be the **continuation ratios**, or discrete stopping probabilities. The **sequential
model** with link \( F^{-1} \) is
\[
\delta_{ir}=F\bigl(\theta_r-\x_{(i)}\T\bbeta_r\bigr),\qquad r=1,\dots,c-1,
\]{#eq-mlt-sequential}

for a continuous, strictly increasing distribution function \( F \). With \( F=\Lambda \) and
common slopes \( \bbeta_r=\bbeta \) it is the **continuation-ratio logit model**. The
category probabilities are
\[
\pi_{ir}=\delta_{ir}\prod_{s<r}(1-\delta_{is})\ \ (r<c),
\qquad
\pi_{ic}=\prod_{s<c}(1-\delta_{is}).
\]{#eq-mlt-sequential-probs}

:::

Two features are visible at once. The sequential model is *not* symmetric in the
direction of the scale: reversing the categories reverses which conditional
probabilities are modelled. And the cutpoints are completely free, since
@eq-mlt-sequential-probs delivers a probability distribution for *any* \( \theta_r \) and
\( \bbeta_r \). There is no analogue of the ordering constraint @eq-mlt-ordering, so
cut-specific slopes cost nothing in coherence here, where in a cumulative model they
cost the model's existence.

## The likelihood falls apart

::: {#thm-mlt-factorization}
[Factorization of the sequential likelihood]

Let \( R_{ir}=\sum_{s\ge r}y_{is} \) be the number of subjects at covariate pattern \( i \)
still "at risk" at stage \( r \), so that \( R_{i1}=m_i \) and \( R_{ir}=y_{ir}+R_{i,r+1} \).
Under @def-mlt-sequential,
\[
\prod_{r=1}^{c}\pi_{ir}^{\,y_{ir}}
=\prod_{r=1}^{c-1}\delta_{ir}^{\,y_{ir}}\,(1-\delta_{ir})^{\,R_{i,r+1}} .
\]{#eq-mlt-factorization}

Consequently:

::: {.enumerate options="label=(\alph*)"}
1. the log-likelihood is \( \ell=\sum_{r=1}^{c-1}\ell_r \), where
      \( \ell_r(\theta_r,\bbeta_r)=\sum_i\{y_{ir}\log\delta_{ir}+R_{i,r+1}\log(1-\delta_{ir})\} \)
      is the log-likelihood of a binary regression with \( R_{ir} \) trials and \( y_{ir} \)
      successes;

2. if the parameters \( (\theta_r,\bbeta_r) \) of different stages are unrestricted and
      distinct, the maximum likelihood estimate is obtained by \( c-1 \) *separate* binary
      regressions, one on the subjects at risk at each stage; the information matrix is block
      diagonal, and the estimates for different stages are asymptotically independent;

3. \( Y_{ir}\mid R_{ir}\sim\text{Bin}(R_{ir},\delta_{ir}) \), and these conditional
      distributions determine the multinomial law.
:::

:::

::: {.proof}
Count exponents on the two sides of @eq-mlt-factorization. On the left, substituting
@eq-mlt-sequential-probs, the factor \( \delta_{ir} \) appears only in \( \pi_{ir} \) for
\( r<c \), with exponent \( y_{ir} \). The factor \( (1-\delta_{is}) \) appears in \( \pi_{ir} \)
for every \( r>s \) and in \( \pi_{ic} \), so its total exponent is
\( \sum_{r>s}y_{ir}=R_{i,s+1} \). That is exactly @eq-mlt-factorization.

(a) Take logarithms of @eq-mlt-factorization and sum over \( i \); by @eq-mlt-sequential
the \( r \)th term involves only \( (\theta_r,\bbeta_r) \), and since
\( R_{ir}=y_{ir}+R_{i,r+1} \), the \( r \)th term is the binomial log-likelihood with
\( R_{ir} \) trials, \( y_{ir} \) successes and success probability \( \delta_{ir} \), up to a
term free of the parameters.

(b) A sum of functions of disjoint groups of parameters is maximized by maximizing each
separately, and its second derivative matrix is block diagonal with the \( r \)th block
\( \partial^2\ell_r/\partial(\theta_r,\bbeta_r)^{\otimes 2} \). Taking expectations, the
information is block diagonal, so by @thm-glm-asymptotics the joint limiting
distribution of the stacked estimate is normal with block diagonal covariance: the
stages are asymptotically independent.

(c) Write the multinomial mass function as
\( P(\Y_i=\y_i)=h(\y_i)\prod_r\pi_{ir}^{y_{ir}} \) and use @eq-mlt-factorization. The
combinatorial factor \( h \) likewise factors as
\( \prod_{r<c}\binom{R_{ir}}{y_{ir}} \), because choosing which subjects stop at each stage
in turn is the same as choosing the whole allocation. So the mass function is the
product of \( c-1 \) binomial mass functions with the stated parameters.
:::

::: {.idea}
A sequential model is barely a multinomial model at all: it is \( c-1 \) binary
regressions with a shared parameterization. Separation (@thm-bin-separation), the
links (@prp-bin-links), residuals (@prp-bin-fit) and case-control sampling
of @thm-bin-casecontrol apply stage by stage, with no new theory.
:::

::: {#cor-mlt-expanded}
[One regression on expanded data]

Let the slopes be common, \( \bbeta_r=\bbeta \), and the data ungrouped (\( m_i=1 \)). Build an
expanded data set with one row for each pair \( (i,r) \) with \( r\le\min(y_i,c-1) \), response
\( \bone\{y_i=r\} \), and regressors \( (\mathbf{e}_r,-\x_{(i)}) \), where \( \mathbf{e}_r \) is the
\( r \)th coordinate vector of length \( c-1 \). Then an ordinary binary regression with link
\( F^{-1} \) on the expanded data has, as a function of \( (\boldsymbol{\uptheta},\bbeta) \), the
log-likelihood of the sequential model; the two fits, their standard errors and their
tests coincide.
:::

::: {.proof}
With \( m_i=1 \), \( y_{ir}=\bone\{y_i=r\} \) and \( R_{i,r+1}=\bone\{y_i>r\} \). The \( (i,r) \)
term of @eq-mlt-factorization is \( \delta_{ir} \) when \( y_i=r \), \( 1-\delta_{ir} \) when
\( y_i>r \), and \( 1 \) when \( y_i<r \), so the rows with \( r>\min(y_i,c-1) \) contribute nothing
and may be dropped. The surviving rows contribute exactly the Bernoulli likelihood with
success probability \( \delta_{ir}=F(\theta_r-\x_{(i)}\T\bbeta) \), which is the model with
linear predictor \( (\mathbf{e}_r,-\x_{(i)})\T(\boldsymbol{\uptheta}\T,\bbeta\T)\T \).
:::

The expanded data set has \( \sum_i\min(y_i,c-1) \) rows, one per stage each subject
reached. That rows from the same subject are not independent looks alarming and is not a
problem: @thm-mlt-factorization says the likelihood is *exactly* the product over rows —
a factorization, not an approximation.

## The complementary log-log link joins the two families

The cumulative and sequential families are genuinely different — except at one link,
where they are the same.

::: {#prp-mlt-cloglog}
[Sequential and cumulative complementary log-log models agree]

Take \( F(z)=1-\exp(-e^{z}) \) and common slopes. If the sequential model
@eq-mlt-sequential holds with cutpoints \( \tilde{\theta}_1,\dots,\tilde{\theta}_{c-1} \),
then the cumulative complementary log-log model @eq-mlt-cumulative holds with the same
\( \bbeta \) and cutpoints
\[
\theta_r=\log\Bigl(\sum_{s\le r}e^{\tilde{\theta}_s}\Bigr),\qquad r=1,\dots,c-1,
\]{#eq-mlt-laara}

which are automatically increasing. The map \( \tilde{\boldsymbol{\uptheta}}\mapsto
\boldsymbol{\uptheta} \) is a bijection onto the increasing sequences, so the two models are
the same family, and the maximum likelihood fits agree.
:::

::: {.proof}
Since \( \{Y_i\ge r\}=\{Y_i>r-1\} \),
\[
P(Y_i>r\mid Y_i>r-1)=1-\delta_{ir}=\exp\bigl\{-e^{\tilde{\theta}_r}e^{-\x_{(i)}\T\bbeta}\bigr\}.
\]
Multiplying these for \( s=1,\dots,r \) telescopes to
\( P(Y_i>r)=\exp\{-H_re^{-\x_{(i)}\T\bbeta}\} \) with \( H_r=\sum_{s\le r}e^{\tilde{\theta}_s} \),
so \( \log\{-\log P(Y_i>r)\}=\log H_r-\x_{(i)}\T\bbeta \), which is @eq-mlt-cumulative with
\( G(z)=1-\exp(-e^{z}) \) and \( \theta_r=\log H_r \). The \( H_r \) are strictly increasing
because \( e^{\tilde{\theta}_s}>0 \), and \( \tilde{\theta}_r=\log(H_r-H_{r-1}) \) recovers the
sequential cutpoints from any increasing \( \boldsymbol{\uptheta} \), with \( H_0=0 \); so the map
is a bijection. Two parameterizations of the same set of distributions have the same
maximized likelihood.
:::

This is the observation of Läärä and Matthews (1985), and it explains why the
complementary log-log link is the natural one for stage data: by @eq-mlt-grouped-ph the
same model is a grouped proportional hazards model, with \( \delta_{ir} \) the discrete
hazard of stopping at stage \( r \). No such identity holds for the logit link, and
[Figure 36.4.1](04-sequential-models.html#fig-mlt-seq) shows how far apart the two families
then are.

::: {when-format="html"}
![**Figure 36.4.1.** Each family is simple in its own coordinates and curved in the
other's. (a) Cumulative logits implied by a sequential logit model with common slopes.
(b) Continuation-ratio logits implied by a proportional odds model. Parallel straight
lines in one panel would mean the other family contained this
one.](sequential_vs_cumulative.svg){#fig-mlt-seq width=100%}
:::

::: {when-format="pdf"}
![Each family is simple in its own coordinates and curved in the other's. (a) Cumulative
logits implied by a sequential logit model with common slopes. (b) Continuation-ratio
logits implied by a proportional odds model. Parallel straight lines in one panel would
mean the other family contained this one.](sequential_vs_cumulative.pdf){width=100%}
:::

## The two families on the same data

::: {#exm-mlt-sequential-anes}
[Self-placement, fitted sequentially]

The self-placement data of @exm-mlt-selfplacement are not a stage process — nobody
becomes extremely conservative by first passing through "extremely liberal" — but the
sequential model is still a parameterization of the seven category probabilities, and
comparing the fits shows what each family buys.

The at-risk counts are 944, 928, 825,
678, 422 and 252, and the expanded data set has
4049 rows. One binary logistic regression on it, as in @cor-mlt-expanded,
gives cutpoints \( (-3.9705,-1.9870,-1.4252,-0.3784,
-0.2744,1.9824) \) and slopes
\( 0.0454 \) (age, per decade),
\( -0.0866 \) (education) and
\( 0.0174 \) (income), with standard errors
\( 0.0270 \), \( 0.0303 \) and
\( 0.0079 \). The signs agree with the proportional odds fit and the
magnitudes are smaller, as one should expect: these are effects on a single stage
transition, not on the whole scale. The multinomial log-likelihood is
\( -1618.2981 \), which the script confirms equals the
log-likelihood the binary regression reports — @thm-mlt-factorization to machine
precision.

Akaike's criterion is 3254.60 for the sequential model against
3245.71 for proportional odds (@def-sel-aic-bic): with nine parameters
each, the cumulative model describes these data better, as a scale of opinion rather
than a ladder of stages should. Letting every stage have its own slopes raises the
log-likelihood to \( -1578.3794 \) and lowers the criterion to
3204.76, just below the unrestricted cumulative model's
3208.29. Once the common-slope restriction goes, the two families have
the same \( (c-1)(1+p) \) parameters and fit these data about equally well.
:::

```{.python .run #cell-sequential-data}
import numpy as np
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
y = anes["selfLR"].astype(int).to_numpy()
X = np.column_stack([anes["age"] / 10.0, anes["educ"], anes["income"]])
n, p = X.shape
c = 7


def seq_probs(theta, beta, X, F=lambda z: 1.0 / (1.0 + np.exp(-z))):
    """Category probabilities of the sequential model with P(Y = r | Y >= r) = F(theta_r - x'beta)."""
    delta = F(theta[None, :] - (X @ beta)[:, None])
    surv = np.column_stack([np.ones(len(X)), np.cumprod(1.0 - delta, axis=1)])
    return np.column_stack([delta * surv[:, :-1], surv[:, -1]])


def seq_loglik(theta, beta, X, y, F=lambda z: 1.0 / (1.0 + np.exp(-z))):
    P = seq_probs(theta, beta, X, F)
    return float(np.sum(np.log(P[np.arange(len(y)), y - 1])))
```

With cut-specific slopes the fit is @thm-mlt-factorization(b): one logistic regression
per stage, on the subjects still at risk.

```{.python .run #cell-sequential-separate}
# unrestricted sequential model: one binary logistic fit per stage, on those still at risk
stage_beta, stage_ll = [], 0.0
for r in range(1, c):
    at_risk = y >= r
    Z = np.column_stack([np.ones(at_risk.sum()), -X[at_risk]])
    fit_r = sm.Logit((y[at_risk] == r).astype(float), Z).fit(disp=0)
    stage_beta.append(np.asarray(fit_r.params))
    stage_ll += fit_r.llf
    print(f"stage {r}: at risk {at_risk.sum():4d}, stop {int((y == r).sum()):4d}, "
          f"theta {fit_r.params[0]:7.4f}, beta {np.round(fit_r.params[1:], 4)}")

print(f"sum of the {c - 1} binary log-likelihoods: {stage_ll:.4f}")
```

With a common slope the fit is @cor-mlt-expanded: one logistic regression on the
expanded data.

```{.python .run #cell-sequential-expanded}
# common slope: one binary logistic regression on the expanded data set
rows_X, rows_y = [], []
for i in range(n):
    for r in range(1, min(y[i], c - 1) + 1):        # stages this respondent reached
        stage = np.zeros(c - 1)
        stage[r - 1] = 1.0
        rows_X.append(np.concatenate([stage, -X[i]]))
        rows_y.append(float(y[i] == r))
Xe, ye = np.array(rows_X), np.array(rows_y)

fit_seq = sm.Logit(ye, Xe).fit(disp=0)
theta_seq, beta_seq = fit_seq.params[:c - 1], fit_seq.params[c - 1:]
se_seq = fit_seq.bse[c - 1:]
ll_seq = seq_loglik(theta_seq, beta_seq, X, y)

print(f"expanded data: {len(ye)} rows from {n} respondents")
print("cutpoints", np.round(theta_seq, 4))
print("slopes   ", np.round(beta_seq, 4), " se", np.round(se_seq, 4))
print(f"binary log-likelihood {fit_seq.llf:.4f} = multinomial log-likelihood {ll_seq:.4f}")
```

Changing the link to the complementary log-log turns the same expanded regression into a
cumulative model, by @prp-mlt-cloglog. The script computes the cumulative cutpoints from
@eq-mlt-laara and confirms that the two ways of assembling the category probabilities
agree to within \( 3\times 10^{-16} \).

```{.python .run #cell-sequential-cloglog}
# Laara and Matthews: with a complementary log-log link the sequential model IS a cumulative model
cll = lambda z: 1.0 - np.exp(-np.exp(z))
fit_cll = sm.GLM(ye, Xe, family=sm.families.Binomial(sm.families.links.CLogLog())).fit()
theta_t, beta_c = fit_cll.params[:c - 1], fit_cll.params[c - 1:]
theta_cum = np.log(np.cumsum(np.exp(theta_t)))                   # the cumulative model's cutpoints

P_step = seq_probs(theta_t, beta_c, X, F=cll)                    # built stage by stage
gamma_cum = cll(theta_cum[None, :] - (X @ beta_c)[:, None])      # built from cumulative probabilities
P_cum = np.diff(np.column_stack([np.zeros(n), gamma_cum, np.ones(n)]), axis=1)

print("sequential cutpoints ", np.round(theta_t, 4))
print("cumulative cutpoints ", np.round(theta_cum, 4))
print("largest difference between the two sets of fitted probabilities:",
      float(np.max(np.abs(P_step - P_cum))))
```

::: {#exm-mlt-sequential-simulation}
[When the mechanism really is sequential]

Generate 4000 observations from a four-stage sequential logit model with
one standard normal regressor, cutpoints \( (-0.6,0,0.5) \) and slope \( \beta=1 \); the
category counts come out as \( 1487 \), \( 1094 \), \( 693 \) and \( 726 \). The sequential fit
recovers the slope as \( 0.9761 \), standard error
\( 0.0309 \). The proportional odds fit of the same data, with the same
number of parameters, returns \( 1.1352 \), which is not estimating
\( \beta \) at all, and Akaike's criterion is 9444.56 against
9515.66. A gap of seventy is not a close call.

The moral is not that one family is better — reverse the simulation and the verdict
reverses with it — but that the choice is a modelling decision about *how the response
came about*, on which the data have plenty to say.
:::

```{.python .run #cell-sequential-simulate}
# a mechanism that really is sequential, fitted by both models
def simulate_sequential(n, theta, beta, rng):
    X = rng.normal(size=(n, len(beta)))
    delta = 1.0 / (1.0 + np.exp(-(theta[None, :] - (X @ beta)[:, None])))
    stop = rng.random(delta.shape) < delta
    y = np.where(stop.any(axis=1), stop.argmax(axis=1) + 1, len(theta) + 1)
    return X, y


rng = np.random.default_rng(20360436)
theta0, beta0 = np.array([-0.6, 0.0, 0.5]), np.array([1.0])
Xs, ys = simulate_sequential(4000, theta0, beta0, rng)
print("simulated counts:", np.bincount(ys)[1:])
```

## Choosing a family

Three families are now available: the cumulative model (@def-mlt-cumulative), the
sequential model (@def-mlt-sequential) and, through @exr-mlt-adjacent,
adjacent-category logits. They are non-nested, have the same dimension in their
common-slope forms, and often fit real data about equally well, so a criterion such as
AIC will decline to choose. The mechanism should choose.

- **Cumulative** when the categories are a graded measurement of something continuous:
  severity, agreement, quality. Its slope is invariant to how the scale was cut
  (@thm-mlt-proportional-odds), which makes results comparable across studies.
- **Sequential** when the categories are stages that must be passed through, especially
  for durations and attainments. Its coefficients describe one transition at a time, and
  @thm-mlt-factorization makes cut-specific slopes cheap.
- **Adjacent-category** when local comparisons are the meaningful ones: category
  \( r+1 \) rather than \( r \), among those in one of the two.

None of them addresses overdispersion relative to the multinomial.
[Chapter 38](../ch38-quasi-likelihood/index.html) treats quasi-likelihood and overdispersion, and
[Chapter 40](../ch40-glmm-gee/index.html) the random-effects models that repeated categorical
measurements need.

## Exercises

### A. Check your understanding

::: {#exr-mlt-seq-count}
[A1]

How many parameters has the sequential model with common slopes? With stage-specific
slopes? How many rows has the expanded data set of @cor-mlt-expanded when every subject
falls in the last category?
:::

::: {.solution}
Common slopes: \( (c-1)+p \). Stage-specific: \( (c-1)(1+p) \). If every \( y_i=c \) then each
subject contributes \( c-1 \) rows, giving \( n(c-1) \); every response in the expanded data is
zero, and the fit does not exist — a separation of the kind in @thm-bin-separation.
:::

::: {#exr-mlt-seq-reverse}
[A2]

Show by a three-category example that reversing the order of the categories changes the
sequential model, but not the cumulative one (@exr-mlt-reverse).
:::

### B. Practice

::: {#exr-mlt-seq-hazard}
[B1]

Interpret \( \delta_{ir} \) as a discrete hazard and \( \prod_{s\le r}(1-\delta_{is}) \) as a
survivor function. Show that under @prp-mlt-cloglog the model says
\( P(Y_i>r\mid\x_{(i)})=S_0(r)^{\exp(-\x_{(i)}\T\bbeta)} \) for a baseline survivor
function \( S_0 \), and identify \( S_0 \).
:::

::: {.solution}
From the proof of @prp-mlt-cloglog, \( P(Y_i>r)=\exp\{-H_re^{-\x\T\bbeta}\} \), so with
\( S_0(r)=e^{-H_r} \) we get \( P(Y_i>r)=S_0(r)^{\exp(-\x\T\bbeta)} \). Here \( S_0 \) is the
survivor function at \( \x=\bzero \) and \( H_r=\sum_{s\le r}e^{\tilde{\theta}_s} \) is the
baseline cumulative hazard, which is exactly the grouped-data statement @eq-mlt-grouped-ph.
:::

### C. Going deeper

::: {#exr-mlt-seq-discrete-survival}
[C1]

Let \( T_i \) be a discrete survival time taking values \( 1,2,\dots \), possibly
right-censored at a known time \( c_i \). Write the likelihood in terms of the discrete
hazards \( \lambda_{ir}=P(T_i=r\mid T_i\ge r) \), and show that it is again the expanded
binary likelihood of @cor-mlt-expanded, now with subject-specific numbers of rows. What
does censoring correspond to in @thm-mlt-factorization?
:::

::: {.solution}
The contribution of an uncensored subject is
\( \lambda_{i,t_i}\prod_{s<t_i}(1-\lambda_{is}) \) and of a censored one
\( \prod_{s\le c_i}(1-\lambda_{is}) \). Both are products of Bernoulli terms over the
periods the subject was at risk, so the expanded binary regression fits them, with a
censored subject contributing only zero responses. In @thm-mlt-factorization censoring is
the case of a subject who is removed from the risk set without stopping: it reduces
\( R_{ir} \) for later \( r \) without adding to any \( y_{ir} \). This is why continuation-ratio
models and discrete-time survival models are the same subject.
:::

::: {#exr-mlt-seq-collapse}
[C2]

Does the sequential model have a collapsing invariance like
@thm-mlt-proportional-odds(b)? Show that merging the *last* two categories leaves the
model and the remaining parameters intact, but merging any other adjacent pair does not.
:::
