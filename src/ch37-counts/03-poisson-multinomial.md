# Poisson and multinomial

A contingency table is a table of counts, so a Poisson log-linear model can be
fitted to it. But the counts usually come from a survey with a fixed sample size,
or from strata with fixed totals, and then they are not independent Poisson
variables at all: they are multinomial. The two descriptions lead to the same fit.
That is why software for Poisson regression is also software for contingency
tables, and it is the bridge to the baseline-category logit models of
[Chapter 36](../ch36-multinomial-ordinal/index.html).

## Conditioning on the total

::: {#thm-cnt-poisson-multinomial}
[Poisson counts conditioned on their total]

::: {.enumerate options="label=(\alph*)"}
1. Let \( Y_1,\dots,Y_c \) be independent, \( Y_j \) Poisson with mean
   \( \mu_j>0 \), and \( N=\sum_jY_j \). Then \( N \) is Poisson with mean
   \( \mu_{+}=\sum_j\mu_j \), and conditionally on \( N=n \) the vector
   \( (Y_1,\dots,Y_c) \) is multinomial with \( n \) trials and probabilities
   \( \pi_j=\mu_j/\mu_{+} \).

2. Suppose the cells are indexed by a stratum \( s=1,\dots,S \) and a category
   \( j=1,\dots,c \), the counts \( Y_{sj} \) are independent Poisson with means
   \( \mu_{sj} \), and the log-linear model is
   \[
   \log\mu_{sj}=\alpha_s+\mathbf{z}_{sj}\T\bgamma ,
   \]{#eq-cnt-strata}

   with \( \boldsymbol{\upalpha}\in\Real^{S} \) unrestricted: that is, the model
   contains a free parameter for every stratum. Then the Poisson log-likelihood
   splits as
   \[
   \ell(\boldsymbol{\upalpha},\bgamma)=\ell_1(\mu_{1+},\dots,\mu_{S+})
   +\ell_2(\bgamma)+\text{constant},
   \]
   where \( \ell_1 \) is the log-likelihood of \( S \) independent Poisson totals
   and \( \ell_2 \) is the log-likelihood of \( S \) independent multinomial
   samples with category probabilities
   \[
   \pi_{sj}(\bgamma)
   =\frac{\exp(\mathbf{z}_{sj}\T\bgamma)}{\sum_{l=1}^{c}\exp(\mathbf{z}_{sl}\T\bgamma)} ,
   \]{#eq-cnt-multinomial-pi}

   and the two pieces have no parameters in common.

3. Consequently the maximum likelihood estimate \( \hat{\bgamma} \) from the
   Poisson fit is the maximum likelihood estimate from the multinomial fit; the
   fitted cell counts are \( \hat\mu_{sj}=n_s\hat\pi_{sj} \) with
   \( n_s=\sum_jy_{sj} \); the observed information for \( \bgamma \) is the same
   in the two models; and deviance differences between models that all contain
   the stratum parameters agree. The stratum parameters
   \( \hat\alpha_s \) themselves have no multinomial counterpart: they are
   determined by \( \hat\mu_{s+}=n_s \).
:::

:::

::: {.proof}
(a) Sums of independent Poisson variables are Poisson with the sum of the means
(convolve the mass functions, or multiply moment generating functions), so
\( \Pr(N=n)=e^{-\mu_{+}}\mu_{+}^{n}/n! \). For nonnegative integers \( y_j \) with
\( \sum_jy_j=n \),
\[
\begin{aligned}
\Pr\bigl(Y_1=y_1,\dots,Y_c=y_c\mid N=n\bigr)
&=\frac{\prod_{j}e^{-\mu_j}\mu_j^{y_j}/y_j!}{e^{-\mu_{+}}\mu_{+}^{n}/n!}\\
&=\frac{n!}{\prod_jy_j!}\prod_{j}\pi_j^{y_j},
\end{aligned}
\]
because \( \sum_j\mu_j=\mu_{+} \) makes the exponential factors cancel and
\( \prod_j\mu_j^{y_j}=\mu_{+}^{n}\prod_j\pi_j^{y_j} \). That is the multinomial
mass function.

(b) Write \( \mu_{s+}=\sum_j\mu_{sj} \) and \( \pi_{sj}=\mu_{sj}/\mu_{s+} \). Under
@eq-cnt-strata, \( \mu_{s+}=e^{\alpha_s}\sum_l\exp(\mathbf{z}_{sl}\T\bgamma) \), and
the ratio \( \pi_{sj} \) is @eq-cnt-multinomial-pi, free of
\( \boldsymbol{\upalpha} \). Substituting \( \mu_{sj}=\mu_{s+}\pi_{sj} \) in the
Poisson log-likelihood,
\[
\begin{aligned}
\ell&=\sum_{s,j}\bigl\{y_{sj}\log\mu_{sj}-\mu_{sj}\bigr\}-\sum_{s,j}\log y_{sj}!\\
&=\sum_{s}\bigl\{n_s\log\mu_{s+}-\mu_{s+}-\log n_s!\bigr\}\\
&\qquad+\sum_{s}\Bigl\{\log n_s!-\sum_j\log y_{sj}!+\sum_jy_{sj}\log\pi_{sj}\Bigr\},
\end{aligned}
\]
which is the stated split, the first sum being \( \ell_1 \) and the second
\( \ell_2 \) plus the multinomial coefficients. For fixed \( \bgamma \), the map
\( \boldsymbol{\upalpha}\mapsto(\mu_{1+},\dots,\mu_{S+}) \) is a bijection from
\( \Real^{S} \) onto \( (0,\infty)^{S} \), so the parameters
\( (\mu_{1+},\dots,\mu_{S+}) \) and \( \bgamma \) vary independently.

(c) The two pieces share no parameters and their ranges are unconstrained, so the
joint maximum is attained by maximizing each separately: \( \ell_1 \) gives
\( \hat\mu_{s+}=n_s \) and \( \ell_2 \) is the multinomial problem, giving
\( \hat{\bgamma} \) and \( \hat\mu_{sj}=n_s\hat\pi_{sj} \). The information for
\( \bgamma \) is minus the second derivative of \( \ell_2 \), and differences of
maximized log-likelihoods between models sharing the stratum parameters lose the
common \( \ell_1 \).
:::

Part (c) is the practical statement: *fit the table as independent Poisson counts,
and provided the model contains a free parameter for every total the design fixed,
you get the multinomial answer.* The converse matters as much: if the design fixed
the row totals and the model omits the row parameters, the two fits differ and the
Poisson one answers a question the design cannot. The proof used \( s \) only to
group cells, so any partition into blocks with fixed totals works the same way.

## Independence in a two-way table

Take \( S=r \) rows and \( c \) columns, and let \( \mathbf{z}_{ij} \) consist of the
column indicators alone, so that @eq-cnt-strata reads
\[
\log\mu_{ij}=\lambda+\lambda^{R}_i+\lambda^{C}_j ,
\]{#eq-cnt-independence}

the log-linear model with row and column main effects and no interaction. The
parameters are redundant by two, and any of the usual side conditions removes the
redundancy without changing the fit (@thm-proj-reparam).

::: {#cor-cnt-independence}
[The independence model]

Under @eq-cnt-independence the conditional distribution of the column given the
row does not depend on the row, so the two classifications are independent. The
maximum likelihood fitted values are
\[
\hat\mu_{ij}=\frac{y_{i+}\,y_{+j}}{n},
\qquad
y_{i+}=\sum_jy_{ij},\quad y_{+j}=\sum_iy_{ij},\quad n=y_{++},
\]{#eq-cnt-independence-fit}

the deviance is \( G^{2}=2\sum_{i,j}y_{ij}\log(y_{ij}/\hat\mu_{ij}) \), the Pearson
statistic is \( X^{2}=\sum_{i,j}(y_{ij}-\hat\mu_{ij})^{2}/\hat\mu_{ij} \), and both
have \( (r-1)(c-1) \) residual degrees of freedom — and, when the expected counts are
large, both are asymptotically \( \chi^{2} \) on that many degrees of
freedom (@thm-glm-deviance).
:::

::: {.proof}
With \( \mathbf{z}_{ij} \) free of \( i \), @eq-cnt-multinomial-pi gives
\( \pi_{ij}=\pi_j \), which is independence. By @thm-cnt-marginals(b) applied to
the row and to the column indicators, \( \hat\mu_{i+}=y_{i+} \) and
\( \hat\mu_{+j}=y_{+j} \). Under @eq-cnt-independence,
\( \hat\mu_{ij}=\hat\mu_{i+}\hat\mu_{+j}/\hat\mu_{++} \), because
\( \hat\mu_{ij}=e^{\lambda}a_ib_j \) factorizes and summing gives
\( \hat\mu_{i+}\hat\mu_{+j}=e^{2\lambda}a_ib_j\sum_lb_l\sum_ka_k
=\hat\mu_{ij}\hat\mu_{++} \). Substituting the matched margins gives
the fit @eq-cnt-independence-fit. The model has \( 1+(r-1)+(c-1) \) free parameters
against \( rc \) cells, leaving \( rc-r-c+1=(r-1)(c-1) \). The deviance formula is
@eq-cnt-deviance with the second term vanishing by @thm-cnt-marginals(a).
:::

The fitted values @eq-cnt-independence-fit are the expected counts of the
elementary chi-squared test for independence, available since Pearson (1900); what
the log-linear framework adds is everything around it.

## Association is interaction on the log scale

Adding an interaction term to @eq-cnt-independence gives the saturated model
\[
\log\mu_{ij}=\lambda+\lambda^{R}_i+\lambda^{C}_j+\lambda^{RC}_{ij},
\]{#eq-cnt-saturated}

with \( rc \) free parameters and \( \hat\mu_{ij}=y_{ij} \).

::: {#prp-cnt-interaction}
[Interaction parameters are log odds ratios]

In @eq-cnt-saturated, for any rows \( i\ne k \) and columns \( j\ne l \),
\[
\lambda^{RC}_{ij}+\lambda^{RC}_{kl}-\lambda^{RC}_{il}-\lambda^{RC}_{kj}
=\log\frac{\mu_{ij}\mu_{kl}}{\mu_{il}\mu_{kj}} ,
\]
which does not depend on the side conditions (@thm-est-side-conditions) used to
identify the parameters.
Under the multinomial reading of @thm-cnt-poisson-multinomial the right-hand side
is the log odds ratio
\( \log\{\pi_{ij}\pi_{kl}/(\pi_{il}\pi_{kj})\} \). The model @eq-cnt-independence
holds if and only if every such quantity is zero.
:::

::: {.proof}
In the four terms of @eq-cnt-saturated the constant, the row effects and the
column effects each appear twice with a plus and twice with a minus sign, so they
cancel, leaving the stated combination equal to
\( \log\mu_{ij}+\log\mu_{kl}-\log\mu_{il}-\log\mu_{kj} \). The cancellation uses
nothing about the parameterization, so any two parameterizations of the same
fitted means agree. Dividing every \( \mu \) by \( \mu_{++} \) leaves the ratio
unchanged, giving the multinomial form; and all such quantities vanish exactly
when \( \log\mu_{ij} \) is additive in \( i \) and \( j \), which
is @eq-cnt-independence.
:::

The combination in @prp-cnt-interaction is the tetrad difference of
@def-tw-interaction applied to \( \log\mu_{ij} \) instead of to a cell mean, so
the apparatus of [Chapter 16](../ch16-multiway-layouts/index.html) transfers:
interaction contrasts survive reparameterization (@prp-tw-interaction-equiv) and
main effects do not. The difference is that there the scale was chosen for
convenience and could be changed by a transformation (@prp-tw-removable), whereas
here the log scale is forced by the link and "no interaction" means
"independence".

::: {#exm-cnt-party}
[Education and party identification]

The 1996 American National Election Studies sample distributed with `statsmodels`
records \( 944 \) respondents. Cross-classify them by education in
four levels (no diploma, high school graduate, some college, degree) and by party
identification collapsed to three, with column totals \( 380 \),
\( 239 \) and \( 325 \).

Fitting the independence model @cor-cnt-independence gives
\( G^{2}=20.541 \) and \( X^{2}=19.338 \) on six
degrees of freedom. The smallest expected count is \( 16.46 \), so
the \( \chi^{2} \) calibration is safe, and the \( p \)-value is
\( 0.0022 \): independence is rejected. The largest Pearson residual
\( (y_{ij}-\hat\mu_{ij})/\sqrt{\hat\mu_{ij}} \) is \( -2.62 \), at
Republicans without a diploma, where \( 10 \) were observed against
\( 22.38 \) expected; the largest positive one,
\( 2.31 \), is at Democrats in the same row.

The alternative need not be saturated. Education is ordered, so give it scores
\( u_i=1,2,3,4 \) and fit
\( \log\mu_{ij}=\lambda+\lambda^{R}_i+\lambda^{C}_j+u_i\delta_j \), spending two
parameters instead of six on the association. Its deviance is
\( 7.367 \) on four degrees of freedom, \( p \)-value
\( 0.1177 \): this model fits, and the drop from independence is
\( 13.174 \) on two degrees of freedom,
\( p=0.00138 \).

Now read the same three fits as multinomial logit models for party given
education (@def-mlt-baseline), with Democrat as the baseline: independence becomes the
intercept-only logit model, the score model a logit model linear in \( u_i \), the
saturated model a logit model with a free effect per education level. Each pair
gives identical fitted counts, and the association parameters \( \delta_j \) are
exactly the education slopes of the logit model,
\( 0.1513 \) for Independent against Democrat and
\( 0.2748 \) for Republican against Democrat: each step up the
education scale multiplies the odds of Republican rather than Democrat by
\( 1.3163 \). In the saturated model the interaction contrast
comparing the degree row with the no-diploma row and the Republican column with
the Democrat column is \( 1.3658 \), so those odds are
\( 3.919 \) times as large among respondents with a degree as
among those without a diploma.
:::

```{.python .run #cell-contingency-table}
import numpy as np
import pandas as pd
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
party = pd.cut(anes["PID"], [-0.5, 1.5, 4.5, 6.5],
               labels=["Democrat", "Independent", "Republican"])
education = pd.cut(anes["educ"], [0.5, 2.5, 3.5, 4.5, 7.5],
                   labels=["no diploma", "high school", "some college", "degree"])
table = pd.crosstab(education, party)
print(table)
```

```{.python .run #cell-contingency-loglinear}
def poisson_mle(X, y, tol=1e-12, maxit=80):
    """Fisher scoring for log mu = X beta."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.mean())
    for _ in range(maxit):
        mu = np.exp(X @ beta)
        XW = X * mu[:, None]
        step = np.linalg.solve(X.T @ XW, XW.T @ (X @ beta + (y - mu) / mu))
        if np.max(np.abs(step - beta)) < tol:
            return step
        beta = step
    return beta


def dummies(levels, k):
    """Indicator columns for levels 1, ..., k-1 of a factor taking the values 0..k-1."""
    return np.column_stack([(levels == j).astype(float) for j in range(1, k)])


counts = table.to_numpy(float)
r, c = counts.shape
score = np.arange(1.0, r + 1)                 # education score u_i = 1, 2, 3, 4
rows, cols = np.divmod(np.arange(counts.size), c)
y = counts.ravel()
main = np.column_stack([np.ones(counts.size), dummies(rows, r), dummies(cols, c)])
linear = np.column_stack([main, score[rows][:, None] * dummies(cols, c)])
saturated = np.column_stack([main] + [dummies(rows, r)[:, [i]] * dummies(cols, c)[:, [j]]
                                      for i in range(r - 1) for j in range(c - 1)])

fits = {}
for name, X in [("independence", main), ("linear", linear), ("saturated", saturated)]:
    beta = poisson_mle(X, y)
    fits[name] = (beta, np.exp(X @ beta).reshape(counts.shape))

print("independence fit:\n", np.round(fits["independence"][1], 2))
print("row totals match:", np.allclose(fits["independence"][1].sum(1), counts.sum(1)))
print("column totals match:", np.allclose(fits["independence"][1].sum(0), counts.sum(0)))
```

```{.python .run #cell-contingency-tests}
from scipy import stats


def deviance(observed, fitted):
    return 2 * np.sum(observed * np.log(observed / fitted))


for name, dfree in [("independence", 6), ("linear", 4)]:
    dev = deviance(counts, fits[name][1])
    pear = np.sum((counts - fits[name][1]) ** 2 / fits[name][1])
    print(f"{name:13s} G2 = {dev:6.3f}  X2 = {pear:6.3f}  df = {dfree}  "
          f"p = {stats.chi2.sf(dev, dfree):.4f}")
```

The next listing fits the three multinomial logit models directly, by Fisher
scoring on the multinomial likelihood, and confirms that each reproduces the
fitted counts of its log-linear partner.

```{.python .run #cell-contingency-multinomial}
def multinomial_logit(counts, design, tol=1e-12, maxit=100):
    """Baseline-category logits: log(pi_ik / pi_i0) = design_i . gamma_k, k = 1..c-1."""
    r, c = counts.shape
    q = design.shape[1]
    gamma = np.zeros((c - 1) * q)
    total = counts.sum(1)
    for _ in range(maxit):
        eta = np.column_stack([np.zeros(r)] +
                              [design @ gamma[k * q:(k + 1) * q] for k in range(c - 1)])
        pi = np.exp(eta) / np.exp(eta).sum(1, keepdims=True)
        mu = total[:, None] * pi
        u = np.concatenate([design.T @ (counts[:, k + 1] - mu[:, k + 1])
                            for k in range(c - 1)])
        info = np.zeros((len(gamma), len(gamma)))
        for k in range(c - 1):
            for l in range(c - 1):
                w = total * pi[:, k + 1] * ((k == l) - pi[:, l + 1])
                info[k * q:(k + 1) * q, l * q:(l + 1) * q] = design.T @ (design * w[:, None])
        step = np.linalg.solve(info, u)
        gamma = gamma + step
        if np.max(np.abs(step)) < tol:
            break
    return gamma, mu


designs = {"independence": np.ones((r, 1)),
           "linear": np.column_stack([np.ones(r), score]),
           "saturated": np.column_stack([np.ones(r), dummies(np.arange(r), r)])}
for name, design in designs.items():
    gamma, mu = multinomial_logit(counts, design)
    print(f"{name:13s} largest difference in fitted counts: "
          f"{np.abs(mu - fits[name][1]).max():.2e}")
```

## More than two classifications

With three classifications the same ideas give a family of models, one for each
independence structure. Writing the saturated model as
\[
\begin{aligned}
\log\mu_{ijk}&=\lambda+\lambda^{A}_i+\lambda^{B}_j+\lambda^{C}_k\\
&\qquad+\lambda^{AB}_{ij}+\lambda^{AC}_{ik}+\lambda^{BC}_{jk}+\lambda^{ABC}_{ijk},
\end{aligned}
\]
the standard submodels delete terms, always a higher-order term before any
lower-order term it contains (a **hierarchical** model). Dropping all two-factor
terms gives mutual independence; keeping \( \lambda^{AC} \) and
\( \lambda^{BC} \) makes \( A \) and \( B \) conditionally independent given
\( C \); dropping only \( \lambda^{ABC} \) gives **homogeneous association**, in
which the \( AB \) odds ratios are the same at every level of \( C \)
(@exr-cnt-three-way, @exr-cnt-conditional).

@thm-cnt-marginals makes the fitting transparent: the likelihood equations equate
observed and fitted counts in exactly the margins spanned by the model's terms.
For mutual independence those are the one-way margins and the fit is a product of
proportions; for homogeneous association they are all three two-way margins, which
has no closed form and is what iterative fitting is for.

::: {.remark}
[Graphical models]

Draw a vertex for each classification and an edge for each two-factor term the
model retains. For the hierarchical models whose structure that graph captures,
conditional independence given a separating set can be read off the picture: these
are the **graphical log-linear models** of Darroch, Lauritzen and Speed (1980),
the discrete ancestors of graphical models generally. They also say when
collapsing over a variable leaves an association unchanged, which is never safe to
assume: the \( AB \) association survives collapsing over \( C \) only when
\( C \) is conditionally independent of \( A \), or of \( B \), given the other
(Bishop, Fienberg and Holland, 1975). Otherwise marginal and conditional
associations may differ in size or in sign, as
[Section 14.5](../ch14-correlation-lack-of-fit-prediction/05-partial-correlation.html)
found for partial correlations.
:::

## Exercises

### A. Check your understanding

::: {#exr-cnt-df-count}
[A1]

A three-way table is \( 3\times4\times2 \). How many cells does it have, how many
free parameters does the homogeneous association model have, and what are the
residual degrees of freedom?
:::

::: {#exr-cnt-poisson-sum}
[A2]

Verify from @thm-cnt-poisson-multinomial(a) that the conditional distribution of
\( Y_1 \) given \( N=n \) is binomial with \( n \) trials and probability
\( \mu_1/\mu_{+} \), and use it to explain why a Poisson model for two counts is
the same as a binomial model for one of them given the total.
:::

### B. Practice

::: {#exr-cnt-independence-hand}
[B1]

For the \( 2\times2 \) table with entries \( y_{11},y_{12},y_{21},y_{22} \), compute
the independence fit @eq-cnt-independence-fit by hand and show that
\( y_{ij}-\hat\mu_{ij}=\pm(y_{11}y_{22}-y_{12}y_{21})/n \) with alternating signs.
Deduce that \( X^{2}=n(y_{11}y_{22}-y_{12}y_{21})^{2}/(y_{1+}y_{2+}y_{+1}y_{+2}) \).
:::

::: {#exr-cnt-conditional}
[B2]

For the conditional independence model
\( \log\mu_{ijk}=\lambda+\lambda^{A}_i+\lambda^{B}_j+\lambda^{C}_k+\lambda^{AC}_{ik}+\lambda^{BC}_{jk} \),
show that the likelihood equations equate the \( AC \) and \( BC \) margins, that
\( \hat\mu_{ijk}=y_{i+k}y_{+jk}/y_{++k} \), and that the residual degrees of freedom
are \( (r-1)(c-1)\ell \) for an \( r\times c\times\ell \) table.
:::

::: {.solution}
The terms in the model span the indicators of the \( AC \) and \( BC \) cells, so
@thm-cnt-marginals equates those margins. The model says
\( \hat\mu_{ijk}=f_{ik}g_{jk} \), so within each level \( k \) it is the
independence model of @cor-cnt-independence applied to the \( k \)th layer, whose
fit is \( y_{i+k}y_{+jk}/y_{++k} \). Each layer contributes \( (r-1)(c-1) \)
residual degrees of freedom, and there are \( \ell \) layers.
:::

::: {#exr-cnt-logit-equivalence}
[B3]

In a \( 2\times c\times\ell \) table with \( A \) binary, show that the homogeneous
association model implies that the logit of \( \Pr(A=1\mid B=j,C=k) \) is additive
in \( j \) and \( k \), and identify the logistic regression that corresponds to it.
Which log-linear terms determine the logistic intercept, which the slopes, and
which do not appear at all?
:::

::: {.solution}
The logit is \( \log(\mu_{1jk}/\mu_{2jk}) \). Substituting the model, the terms
\( \lambda \), \( \lambda^{B}_j \), \( \lambda^{C}_k \) and \( \lambda^{BC}_{jk} \)
cancel, leaving
\( (\lambda^{A}_1-\lambda^{A}_2)+(\lambda^{AB}_{1j}-\lambda^{AB}_{2j})
+(\lambda^{AC}_{1k}-\lambda^{AC}_{2k}) \), which is additive: it is the
main-effects logistic regression of \( A \) on \( B \) and \( C \). The \( A \) main
effect gives the intercept and the \( AB \), \( AC \) terms the slopes;
\( \lambda^{BC} \), the association between the explanatory factors, does not
appear.
:::

### C. Going deeper

::: {#exr-cnt-uniform-association}
[C1]

The model of @exm-cnt-party with row scores,
\( \log\mu_{ij}=\lambda+\lambda^{R}_i+\lambda^{C}_j+u_i\delta_j \), is called a
row-effects model. Show that its local log odds ratios
\( \log\{\mu_{ij}\mu_{i+1,j+1}/(\mu_{i,j+1}\mu_{i+1,j})\} \) equal
\( (u_{i+1}-u_i)(\delta_{j+1}-\delta_j) \). Then show that giving *both*
classifications scores, with association term \( \theta u_iv_j \), makes every
local log odds ratio \( \theta(u_{i+1}-u_i)(v_{j+1}-v_j) \) — the **uniform
association** model — and count its parameters.
:::

::: {#exr-cnt-three-way}
[C2]

For a three-way table, show that the likelihood equations of the mutual
independence model are \( \hat\mu_{i++}=y_{i++} \), \( \hat\mu_{+j+}=y_{+j+} \) and
\( \hat\mu_{++k}=y_{++k} \), and solve them to get
\( \hat\mu_{ijk}=y_{i++}y_{+j+}y_{++k}/n^{2} \). Explain why the homogeneous
association model has no such closed form.
:::
