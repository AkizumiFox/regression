# Orthogonal designs

[Section 9.3](03-sequential-partial.html) showed how much the sums of squares for a term can
depend on the other terms. This section asks when they do not. The answer is geometric: after
the mean is removed, the terms must span mutually orthogonal subspaces. Balanced experiments
achieve this by design and observational data essentially never do, but orthogonality
can be manufactured inside a term by a change of basis, as with orthogonal polynomials and
orthogonal contrasts.

## When the order does not matter

Let the model matrix be \( \X=[\bone,\X_1,\dots,\X_k] \) with an intercept and \( k \) terms. Write
\( \bP_0=n^{-1}\bone\bone\T \), and let
\[
\tilde{\X}_j=(\I-\bP_0)\X_j,\qquad \mathcal V_j=\C(\tilde{\X}_j),
\]
be the centred columns of term \( j \) and the subspace they span. By @lem-ss-residualized with
\( \X_1=\bone \), \( \C(\bone,\X_j)=\spn(\bone)\dirsum\mathcal V_j \), so
\( \text{SS}(\X_j\mid\bone)=\norm{\bP_{\mathcal V_j}\y}^2 \), the sum of squares of the term
fitted on its own after the mean. Call the design **orthogonal** (for these terms) if the
subspaces \( \mathcal V_1,\dots,\mathcal V_k \) are mutually orthogonal, that is, if
\( \tilde{\X}_i\T\tilde{\X}_j=\bzero \) for \( i\neq j \). For numerical regressors this says that
every pair of regressors from different terms has sample correlation zero.

::: {#thm-ss-orthogonal-design}
[Orthogonal designs]

For a model \( \X=[\bone,\X_1,\dots,\X_k] \), the following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. the centred subspaces \( \mathcal V_1,\dots,\mathcal V_k \) are mutually orthogonal;

2. \( \M=\bP_0+\bP_{\mathcal V_1}+\dots+\bP_{\mathcal V_k} \);

3. for every \( j \) and every \( \y \), the marginal and partial sums of squares agree:
   \( \text{SS}(\X_j\mid\bone)=\text{SS}(\X_j\mid\X_{(-j)}) \);

4. for every \( \y \), the marginal sums of squares add up to the regression sum of squares:
   \( \sum_j\text{SS}(\X_j\mid\bone)=\text{SS}(\X_1,\dots,\X_k\mid\bone) \).
:::

When they hold, \( \text{SS}(\X_j\mid\bone,\X_S)=\norm{\bP_{\mathcal V_j}\y}^2 \) for every set \( S \) of
other terms. In particular the sequential sum of squares of \( \X_j \) is the same in every order,
and it equals the partial and the Type II sums of squares of \( \X_j \).
:::

::: {.proof}
*(a) implies (b) and the final statement.* Let \( S \) be any set of terms. The subspaces
\( \spn(\bone) \) and \( \mathcal V_l \), \( l\in S \), are mutually orthogonal, because each
\( \mathcal V_l\subseteq\bone\perpc \), and they sum to \( \C(\bone,\X_S) \), because each
\( \C(\bone,\X_l)=\spn(\bone)\dirsum\mathcal V_l \). By @thm-proj-sum the projection onto
\( \C(\bone,\X_S) \) is \( \bP_0+\sum_{l\in S}\bP_{\mathcal V_l} \). With \( S \) all terms this is (b).
With \( S \) and \( S\cup\{j\} \), the difference of the two projections is \( \bP_{\mathcal V_j} \), which
gives \( \text{SS}(\X_j\mid\bone,\X_S)=\norm{\bP_{\mathcal V_j}\y}^2 \). Sequential, partial and Type II
sums of squares of \( \X_j \) are all of this form.

*(b) implies (d)*: \( \text{SS}(\X_1,\dots,\X_k\mid\bone)=\y\T(\M-\bP_0)\y=\sum_j\y\T\bP_{\mathcal V_j}\y \).

*(d) implies (a).* The two sides of (d) are the quadratic forms of the symmetric matrices
\( \sum_j\bP_{\mathcal V_j} \) and \( \M-\bP_0 \). Two symmetric matrices with the same quadratic form
are equal, because \( 2\bu\T\A\bv=(\bu+\bv)\T\A(\bu+\bv)-\bu\T\A\bu-\bv\T\A\bv \). So
\( \sum_j\bP_{\mathcal V_j}=\M-\bP_0 \), which is idempotent (@thm-proj-nested). By the algebraic
form of Cochran's theorem (@thm-qf-cochran-algebra), projections whose sum is idempotent
satisfy \( \bP_{\mathcal V_i}\bP_{\mathcal V_j}=\bzero \) for \( i\neq j \), which is (a).

*(a) implies (c)* is part of the final statement.

*(c) implies (a).* As before, equality of the quadratic forms gives
\( \bP_{\mathcal V_j}=\M-\M_{(-j)} \). The right-hand side projects onto
\( \C(\X)\ominus\C(\X_{(-j)}) \) (@thm-proj-nested), which is orthogonal to \( \C(\X_{(-j)}) \). The
space \( \C(\X_{(-j)}) \) contains \( \mathcal V_i \) for every \( i\neq j \), since
\( \tilde{\X}_i=\X_i-\bP_0\X_i \) has columns in \( \C(\bone,\X_i) \). So \( \mathcal V_j\perp\mathcal V_i \).
:::

::: {.idea}
Orthogonality is a property of the *centred* subspaces. The raw columns of two indicator
variables are never orthogonal to the intercept, and a regressor and its square are rarely
orthogonal. What matters is whether the parts of the terms that vary around their means point
in perpendicular directions.
:::

Orthogonality also stabilizes the coefficients.

::: {#cor-ss-orthogonal-coefficients}
[Coefficients in orthogonal designs]

Suppose \( \X=[\bone,\X_1,\dots,\X_k] \) has full column rank and satisfies the conditions of @thm-ss-orthogonal-design.
Then the least squares coefficients of the columns of \( \X_j \) are the
same in every model that contains \( \bone \) and \( \X_j \), whichever other terms it contains. Only
the intercept changes.
:::

::: {.proof}
In the model with terms \( S\ni j \), the fitted vector is
\( \bP_0\y+\sum_{l\in S}\bP_{\mathcal V_l}\y \) by the proof of the theorem, and
\( \bP_{\mathcal V_l}\y=\tilde{\X}_l\bb_l \) with \( \bb_l=(\tilde{\X}_l\T\tilde{\X}_l)^{-1}\tilde{\X}_l\T\y \),
which does not depend on \( S \). Since \( \tilde{\X}_l\bb_l=\X_l\bb_l-\bP_0\X_l\bb_l \), the fitted
vector equals \( \sum_{l\in S}\X_l\bb_l \) plus a multiple of \( \bone \). The coefficients are unique
by full rank, so the coefficient vector of \( \X_j \) is \( \bb_j \) in every such model.
:::

With uncorrelated regressors, dropping one leaves the others' estimates as they were, and each
regressor has an order-free sum of squares. This is one reason factorial experiments are laid
out on orthogonal arrays (Chapters [16](../ch16-multiway-layouts/index.html) and
[18](../ch18-covariance-and-design/index.html)).

## Balanced and proportional layouts

For two factors the orthogonality condition reduces to a condition on the cell counts.

::: {#thm-ss-proportional}
[Proportional frequencies]

In a two-way layout with factor \( A \) at \( a \) levels, factor \( B \) at \( b \) levels, cell counts
\( n_{ij} \), and all row totals \( n_{i\cdot} \) and column totals \( n_{\cdot j} \) positive, the centred
subspaces of \( A \) and \( B \) are orthogonal iff
\[
n_{ij}=\frac{n_{i\cdot}\,n_{\cdot j}}{n}\qquad\text{for all }i,j .
\]{#eq-ss-proportional}

:::

::: {.proof}
A vector in the centred subspace of \( A \) takes a value \( \alpha_i \) on every observation in row
\( i \), with \( \sum_in_{i\cdot}\alpha_i=0 \). A vector in that of \( B \) takes the value \( \beta_j \) in
column \( j \), with \( \sum_jn_{\cdot j}\beta_j=0 \). Their inner product is
\( \sum_{i,j}n_{ij}\alpha_i\beta_j=\boldsymbol{\upalpha}\T\mathbf{N}\boldsymbol{\upbeta} \), with \( \mathbf{N}=(n_{ij}) \). Write
\( \br=(n_{1\cdot},\dots,n_{a\cdot})\T \) and \( \mathbf{c}=(n_{\cdot1},\dots,n_{\cdot b})\T \). The subspaces are
orthogonal iff \( \boldsymbol{\upalpha}\T\mathbf{N}\boldsymbol{\upbeta}=0 \) whenever \( \boldsymbol{\upalpha}\T\br=0 \) and \( \mathbf{c}\T\boldsymbol{\upbeta}=0 \).

If @eq-ss-proportional holds, \( \mathbf{N}=\br\mathbf{c}\T/n \) and
\( \boldsymbol{\upalpha}\T\mathbf{N}\boldsymbol{\upbeta}=(\boldsymbol{\upalpha}\T\br)(\mathbf{c}\T\boldsymbol{\upbeta})/n=0 \). Conversely, let
\( \mathbf{E}=\mathbf{N}-\br\mathbf{c}\T/n \). Its row sums are \( \br-\br(\mathbf{c}\T\bone)/n=\bzero \) and its column sums are
\( \mathbf{c}\T-(\bone\T\br)\mathbf{c}\T/n=\bzero\T \), since \( \bone\T\br=\mathbf{c}\T\bone=n \). So \( \boldsymbol{\upalpha}\T\mathbf{E}\boldsymbol{\upbeta} \)
does not change when \( \boldsymbol{\upalpha} \) or \( \boldsymbol{\upbeta} \) is shifted by a multiple of \( \bone \). Given any
\( \boldsymbol{\upalpha} \) and \( \boldsymbol{\upbeta} \), shift them to \( \boldsymbol{\upalpha}'=\boldsymbol{\upalpha}-(\boldsymbol{\upalpha}\T\br/n)\bone \) and
\( \boldsymbol{\upbeta}'=\boldsymbol{\upbeta}-(\mathbf{c}\T\boldsymbol{\upbeta}/n)\bone \), which satisfy \( {\boldsymbol{\upalpha}'}\T\br=0 \) and
\( \mathbf{c}\T\boldsymbol{\upbeta}'=0 \). Then
\( \boldsymbol{\upalpha}\T\mathbf{E}\boldsymbol{\upbeta}={\boldsymbol{\upalpha}'}\T\mathbf{E}\boldsymbol{\upbeta}'={\boldsymbol{\upalpha}'}\T\mathbf{N}\boldsymbol{\upbeta}'-0=0 \) by assumption. As
\( \boldsymbol{\upalpha} \) and \( \boldsymbol{\upbeta} \) are arbitrary, \( \mathbf{E}=\bzero \).
:::

Equal replication, \( n_{ij}=m \), is the common case (@exm-proj-two-factor,
@exr-proj-two-factor-interaction). If the interaction subspace is defined as the orthogonal
complement of the additive model within the cell-constant vectors, then under @eq-ss-proportional
the subspaces for the mean, \( A \), \( B \) and \( AB \) form a decomposition whatever the order.

Orthogonality settles Type I against Type II, but not Type III. Under proportional frequencies the
Type I and Type II sums of squares of \( A \) coincide and address equality of row means weighted by
the column proportions \( n_{\cdot j}/n \) (@exr-ss-weights-coincide), while Type III weights the
columns equally. They agree with equal replication or without interaction.

::: {#exm-ss-grunfeld-balanced}
[Firms and years]

Grunfeld's investment panel (@exm-proj-grunfeld) records \( 11 \) firms in each of
\( 20 \) years, one observation per firm and year. As a two-way layout with factors
firm and year, it has \( n_{ij}=1 \), so @eq-ss-proportional holds, and the additive model's
sums of squares for the logarithm of investment are order-free:

| Source | df | Sum of squares |
|---|---|---|
| firm | 10 | 473.482 |
| year | 19 | 26.749 |
| residual | 190 | 14.667 |
| corrected total | 219 | 514.898 |

Firms differ enormously in the scale of their investment, and years matter much less.
Adding the logarithm of market value as a covariate destroys the orthogonality, because value
varies with both firm and year. The sum of squares for firm then depends on what else is in the
model: it is \( 90.701 \) after year and value, and \( 88.443 \) after value
alone. Most of the variation between firms is variation in size, which value already
measures. This is the setting of analysis of covariance ([Chapter 18](../ch18-covariance-and-design/index.html), @thm-dsn-ancova).
:::

```{.python .run #cell-orthogonal-design-grunfeld}
import numpy as np
import pandas as pd
import statsmodels.api as sm

def proj(*blocks):
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]
    return U @ U.T

g = sm.datasets.grunfeld.load_pandas().data
y = np.log(g["invest"].to_numpy())
n = len(y)
one = np.ones((n, 1))
F = pd.get_dummies(g["firm"]).to_numpy(float)          # 11 firms
T = pd.get_dummies(g["year"]).to_numpy(float)          # 20 years
P0, PF, PT, PFT = proj(one), proj(one, F), proj(one, T), proj(one, F, T)

print("centred spaces orthogonal:", np.abs((PF - P0) @ (PT - P0)).max() < 1e-10)
print(f"firm:  first {y @ (PF - P0) @ y:8.3f}   after year {y @ (PFT - PT) @ y:8.3f}")
print(f"year:  first {y @ (PT - P0) @ y:8.3f}   after firm {y @ (PFT - PF) @ y:8.3f}")
```

```{.python .run #cell-orthogonal-design-covariate}
v = np.log(g["value"].to_numpy())[:, None]             # log market value
PFv, PTv, PFTv = proj(one, F, v), proj(one, T, v), proj(one, F, T, v)
print(f"with log value: firm after year, value {y @ (PFTv - PTv) @ y:8.3f}"
      f"   firm after value {y @ (PFv - proj(one, v)) @ y:8.3f}")
```

## Orthogonal polynomials

A polynomial regression in a single variable \( x \), with distinct values
\( x_1,\dots,x_n \), has a built-in chain of models,
\[
\mathcal P_0\subseteq\mathcal P_1\subseteq\dots\subseteq\mathcal P_d,
\qquad \mathcal P_k=\spn(\bone,\x,\x^2,\dots,\x^k),
\]
where \( \x^k \) is the vector with entries \( x_i^k \). The decomposition of interest is that of the chain (@thm-ss-decomposition(e)). Its \( k \)th subspace,
\( \mathcal P_{k-1}\perpc\cap\mathcal P_k \), is a line. Applying Gram–Schmidt (@prp-proj-gram-schmidt) to \( \bone,\x,\x^2,\dots \) produces
unit vectors \( \mathbf{q}_0,\mathbf{q}_1,\dots,\mathbf{q}_d \) with
\( \spn(\mathbf{q}_0,\dots,\mathbf{q}_k)=\mathcal P_k \) for every \( k \). So \( \mathbf{q}_k \) spans that line. The
vectors \( \mathbf{q}_k \), evaluated at the design points, are the **orthogonal polynomials** of the
design. Three consequences follow at once.

1. The sequential sum of squares for degree \( k \), the drop in residual sum of squares from degree
   \( k-1 \) to degree \( k \), is \( (\mathbf{q}_k\T\y)^2 \), with one degree of freedom.

2. The columns \( \mathbf{q}_1,\dots,\mathbf{q}_d \) are orthogonal and centred, so as separate terms they form
   an orthogonal design (@thm-ss-orthogonal-design), and their sums of squares do not depend on
   the order of entry.

3. The coefficient of \( \mathbf{q}_k \) is \( \mathbf{q}_k\T\y \) in every fit that includes it
   (@cor-ss-orthogonal-coefficients). Raising the degree never changes the lower coefficients.

For equally spaced \( x \) the orthogonal polynomials have integer multiples that are widely
tabulated. For three points they are \( (-1,0,1) \) and \( (1,-2,1) \). For four points they are
\( (-3,-1,1,3) \), \( (1,-1,-1,1) \) and \( (-1,3,-3,1) \). They serve as contrasts for a quantitative
factor with equally spaced levels, splitting the treatment sum of squares into linear, quadratic
and higher components. Their numerical merits and the recurrence that generates them belong to
[Chapter 42](../ch42-polynomials-piecewise/index.html) (@thm-ply-recurrence).

::: {#exm-ss-nile-polynomials}
[A trend in the Nile]

The annual flow of the Nile at Aswan was recorded for the \( 100 \) years 1871–1970. How
much of its variation is a smooth trend? The listing builds the orthonormal polynomials
\( \mathbf{q}_0,\dots,\mathbf{q}_4 \) for the years by a QR factorization of the matrix of powers, and reads
off the sequential sums of squares (in units of \( (10^8\,\text{m}^3)^2 \)):

| Degree | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Sequential sum of squares | 613893 | 309415 | 1894 | 134142 |

The corrected total is \( 2835157 \). The linear term accounts for a fraction
\( 0.217 \) of it, and the first three terms together for \( 0.326 \). Against
the residual mean square of the cubic fit, \( 19895 \) on \( 96 \) degrees of freedom,
the quadratic component has \( F=15.55 \) and the cubic only \( 0.10 \).
The quartic sum of squares is large again, a sign that the pattern is not a gentle
polynomial at all. The flow dropped abruptly around the end of the nineteenth century
(Cobb 1978), and a low-degree polynomial can only approximate a step
([Figure 9.4.1](#fig-ss-nile)).

With the raw powers of the centred and scaled year, the sums of squares depend on the order:
the linear column entered after the quadratic and cubic ones receives \( 74751 \) instead of
\( 613893 \). The raw linear coefficient is \( -134.4 \) in the linear and
quadratic fits, but \( -117.3 \) in the cubic and quartic ones. It survives the
quadratic term because the years are symmetric about their mean, which makes odd and even
centred powers orthogonal (@exr-ss-symmetric-design), but not the cubic term. The coefficient of
\( \mathbf{q}_1 \) is \( -783.5 \) in every fit.
:::

::: {when-format="html"}
![**Figure 9.4.1.** (a) The orthonormal polynomials \( \mathbf{q}_1,\mathbf{q}_2,\mathbf{q}_3 \) for the years
1871–1970. (b) The Nile flow with the nested least squares fits of degrees 1, 2 and 3. The
quadratic and cubic fits nearly coincide, because the cubic sum of squares is
tiny.](orthogonal_poly.svg){#fig-ss-nile width=100%}
:::

::: {when-format="pdf"}
![(a) The orthonormal polynomials \( \mathbf{q}_1,\mathbf{q}_2,\mathbf{q}_3 \) for the years
1871–1970. (b) The Nile flow with the nested least squares fits of degrees 1, 2 and 3. The
quadratic and cubic fits nearly coincide, because the cubic sum of squares is
tiny.](orthogonal_poly.pdf){width=100%}
:::

```{.python .run #cell-orthogonal-poly-basis}
import numpy as np
import statsmodels.api as sm
nile = sm.datasets.nile.load_pandas().data
year = nile["year"].to_numpy()
y = nile["volume"].to_numpy()

t = year - year.mean()                              # centring only helps rounding
V = np.vander(t / t.max(), 5, increasing=True)      # 1, t, t^2, t^3, t^4 (scaled)
Q, R = np.linalg.qr(V)                              # Gram-Schmidt in matrix form
Q = Q * np.sign(np.diag(R))                         # make each q_k agree in sign with t^k
coef = Q.T @ y                                      # coordinates of y on q_0, ..., q_4
seq_ss = coef[1:] ** 2                              # sequential SS: linear, ..., quartic
for k, s in enumerate(seq_ss, start=1):
    print(f"degree {k}: sequential SS = {s:10.0f}")
print(f"residual after degree 4: {y @ y - np.sum(coef ** 2):10.0f}")
```

## Single-degree-of-freedom decompositions

Orthogonal polynomials split a sum of squares into one-dimensional pieces along directions
chosen for their meaning. The same can be done for any hypothesis, by choosing single
estimable functions whose estimates are uncorrelated.

::: {#thm-ss-single-df}
[Orthogonal single-degree-of-freedom sums of squares]

Let \( \blambda_1,\dots,\blambda_s \) be nonzero vectors with \( \blambda_k\T\bbeta \) estimable, so that
\( \blambda_k=\X\T\boldsymbol{\uprho}_k \), and let \( \G \) be a generalized inverse of \( \X\T\X \).

::: {.enumerate options="label=(\alph*)"}
1. The sum of squares for \( \blambda_k\T\bbeta=0 \) is
   \[
\text{SS}_k=\frac{(\blambda_k\T\hbeta)^2}{\blambda_k\T\G\blambda_k}=\norm{\bP_{\mathbf{v}_k}\y}^2,
   \qquad \mathbf{v}_k=\M\boldsymbol{\uprho}_k ,
\]
   the squared length of the projection of \( \y \) onto the line spanned by \( \mathbf{v}_k \).

2. \( \mathbf{v}_h\T\mathbf{v}_k=\blambda_h\T\G\blambda_k \). If \( \blambda_h\T\G\blambda_k=0 \) for all \( h\neq k \),
   the lines are mutually orthogonal, and \( \sum_k\text{SS}_k \) is the sum of squares for the joint
   hypothesis \( \blambda_1\T\bbeta=\dots=\blambda_s\T\bbeta=0 \), on \( s \) degrees of freedom.

3. Conversely, every orthonormal basis \( \bu_1,\dots,\bu_s \) of a subspace of \( \C(\X) \) arises in this
   way, with \( \blambda_k=\X\T\bu_k \).
:::

:::

::: {.proof}
(a) is @thm-ss-hypothesis with \( q=1 \) and \( \mathbf{d}=0 \): the test space is spanned by
\( \M\boldsymbol{\uprho}_k=\mathbf{v}_k \), and the projection of \( \y \) onto it has squared length
\( (\mathbf{v}_k\T\y)^2/\mathbf{v}_k\T\mathbf{v}_k \), where \( \mathbf{v}_k\T\y=\boldsymbol{\uprho}_k\T\M\y=\blambda_k\T\hbeta \) and
\( \mathbf{v}_k\T\mathbf{v}_k=\blambda_k\T\G\blambda_k \). (b) \( \mathbf{v}_h\T\mathbf{v}_k=\boldsymbol{\uprho}_h\T\M\boldsymbol{\uprho}_k=\boldsymbol{\uprho}_h\T\X\G\X\T\boldsymbol{\uprho}_k=\blambda_h\T\G\blambda_k \)
by @thm-proj-M-formula. Orthogonal nonzero vectors are independent, so
\( \bLambda=[\blambda_1,\dots,\blambda_s] \) has rank \( s \) (because \( \X\T\mathbf{V}=\bLambda \) with
\( \mathbf{V}=[\mathbf{v}_1,\dots,\mathbf{v}_s] \), and \( \X\T \) is one-to-one on \( \C(\X) \)). The test space of the joint
hypothesis is \( \C(\M\bT)=\C(\mathbf{V}) \), with \( \bT=[\boldsymbol{\uprho}_1,\dots,\boldsymbol{\uprho}_s] \), by @thm-ss-hypothesis(b), and the projection onto a span of
orthogonal lines is the sum of the projections onto the lines (@thm-proj-sum). (c) With
\( \boldsymbol{\uprho}_k=\bu_k \), \( \mathbf{v}_k=\M\bu_k=\bu_k \), and \( \blambda_h\T\G\blambda_k=\bu_h\T\bu_k \).
:::

In the one-way layout with group sizes \( n_1,\dots,n_g \), take the generalized inverse of @exr-ss-oneway-contrast.
A contrast \( \sum_ic_i\alpha_i \) with \( \sum_ic_i=0 \) (@def-est-contrast)
is estimated by \( \sum_ic_i\bar{y}_i \), and \( \blambda\T\G\blambda=\sum_ic_i^2/n_i \). So two
contrasts \( \mathbf{c} \) and \( \mathbf{c}' \) give orthogonal sums of squares iff
\[
\sum_{i=1}^g\frac{c_ic_i'}{n_i}=0 ,
\]
and \( g-1 \) such contrasts split the between-groups sum of squares into single degrees of
freedom. With equal group sizes this is ordinary orthogonality, \( \sum_ic_ic_i'=0 \); with unequal
sizes, contrasts that look orthogonal need not be.

::: {#exm-ss-education-contrasts}
[Splitting the education sum of squares]

Consider the \( 944 \) self-placements of @exm-ss-anes-types by education alone. The group sizes
are \( 313 \), \( 277 \) and \( 354 \), and the group means are \( 4.476 \),
\( 4.357 \) and \( 4.167 \). The between-groups sum of squares is
\( 16.306 \) on two degrees of freedom. The contrast of college against graduate,
\( \mathbf{c}_2=(0,1,-1) \), has sum of squares \( 5.653 \). The contrast of high school against
the other two, weighted by their sizes, \( \mathbf{c}_1=(1,-0.439,-0.561) \), is orthogonal
to \( \mathbf{c}_2 \) in the sense above, and has sum of squares \( 10.652 \). The two add up to
the between-groups sum of squares, up to rounding. The more obvious contrast \( (1,-\tfrac12,-\tfrac12) \),
high school against the plain average of the others, has sum of squares \( 9.534 \),
and together with \( \mathbf{c}_2 \) it gives \( 15.187 \): with unequal group sizes it is not
orthogonal to \( \mathbf{c}_2 \).
:::

```{.python .run #cell-single-df-contrasts}
import numpy as np
import pandas as pd
import statsmodels.api as sm
d = sm.datasets.anes96.load_pandas().data
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"])
y = d["selfLR"].to_numpy()

counts = educ.value_counts(sort=False).to_numpy().astype(float)
means = d.groupby(educ, observed=True)["selfLR"].mean().to_numpy()
grand = y.mean()
ss_treat = np.sum(counts * (means - grand) ** 2)          # SS(educ | 1), 2 df

def contrast_ss(c):
    """(c' ybar)^2 / sum(c_i^2 / n_i): the SS for the hypothesis c' mu = 0."""
    return (c @ means) ** 2 / np.sum(c ** 2 / counts)

def orthogonal(c, h):
    return np.isclose(np.sum(c * h / counts), 0.0)

n2, n3 = counts[1], counts[2]
c1 = np.array([1.0, -n2 / (n2 + n3), -n3 / (n2 + n3)])   # HS against the rest, weighted
c2 = np.array([0.0, 1.0, -1.0])                          # College against Graduate
c1_plain = np.array([1.0, -0.5, -0.5])                   # HS against the plain average
print("orthogonal:", orthogonal(c1, c2), orthogonal(c1_plain, c2))
print(f"treatment SS {ss_treat:.3f} = {contrast_ss(c1):.3f} + {contrast_ss(c2):.3f}")
print(f"with the plain contrast: {contrast_ss(c1_plain):.3f} + {contrast_ss(c2):.3f}")
```

Single-degree-of-freedom sums of squares reappear in [Chapter 13](../ch13-multiplicity/index.html), where many of them are tested at
once, and in [Chapter 15](../ch15-anova-subspaces/index.html), where orthogonal contrasts organize the one-way
analysis of variance (@thm-aov-orthogonal-contrasts).

## Exercises

### A. Check your understanding

::: {#exr-ss-four-point}
[A1]

With \( x=1,2,3,4 \) and \( \y=(3,5,4,8)\T \), use the orthogonal polynomial contrasts for four equally
spaced points to compute the linear, quadratic and cubic sums of squares, and check that they add
up to the corrected total sum of squares.
:::

::: {#exr-ss-small-layouts}
[A2]

Which of the \( 2\times2 \) count tables \( \begin{psmallmatrix}2&4\\3&6\end{psmallmatrix} \),
\( \begin{psmallmatrix}2&4\\4&2\end{psmallmatrix} \) and \( \begin{psmallmatrix}1&3\\3&9\end{psmallmatrix} \)
give orthogonal centred factor spaces? For the one that does not, compute the inner product of the
two centred indicator vectors of the first levels.
:::

::: {#exr-ss-uncorrelated-regressors}
[A3]

Show that if the regressors \( \x_1,\dots,\x_k \) of a model with an intercept are pairwise
uncorrelated, then \( R^2=\sum_jr_{yj}^2 \), where \( r_{yj} \) is the sample correlation of \( \y \) and
\( \x_j \).
:::

### B. Practice

::: {#exr-ss-balanced-types}
[B1]

In a two-way layout with interaction and \( n_{ij}=m\ge1 \) for all cells, show that the Type I sums
of squares of \( A \) in either order, and its Type II and Type III sums of squares, are all equal to
\( bm\sum_i(\bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot})^2 \).
:::

::: {.solution}
Equal counts satisfy @eq-ss-proportional, so by @thm-ss-orthogonal-design the Type I sums of
squares of \( A \) in both orders and its Type II sum of squares equal \( \text{SS}(A\mid\mu)=\norm{(\bP_A-\bP_0)\y}^2 \).
Each row has \( bm \) observations, and \( \bP_A-\bP_0 \) replaces an observation in row \( i \) by
\( \bar{y}_{i\cdot\cdot}-\bar{y}_{\cdot\cdot\cdot} \), which gives the formula. For Type III, the
hypothesis of @thm-ss-two-way-hypotheses(b) is that of (a), because with equal counts weighted and
unweighted row means coincide. A sum of squares for a hypothesis depends only on the hypothesis
(it is \( \text{SSE}_0-\text{SSE} \)), so the Type III sum of squares equals \( \text{SS}(A\mid\mu) \) too.
:::

::: {#exr-ss-proportional-not-equal}
[B2]

Take counts \( n_{11}=1 \), \( n_{12}=2 \), \( n_{21}=2 \), \( n_{22}=4 \). Check @eq-ss-proportional. Find cell
means for which the Type III hypothesis for \( A \) holds but the Type I (and Type II) hypothesis
fails.
:::

::: {#exr-ss-one-way-helmert}
[B3]

In a one-way layout with \( g \) groups of equal size \( m \), show that the Helmert contrasts
\( \mathbf{c}_k=(1,\dots,1,-(k-1),0,\dots,0) \), with \( k-1 \) leading ones, \( k=2,\dots,g \), are mutually
orthogonal, and write the between-groups sum of squares as a sum of \( g-1 \) single-degree-of-freedom
sums of squares. What goes wrong with unequal group sizes?
:::

### C. Going deeper

::: {#exr-ss-symmetric-design}
[C1]

Suppose the design points are symmetric about their mean: after sorting,
\( x_i-\bar{x}=-(x_{n+1-i}-\bar{x}) \). Show that every odd power of the centred variable is orthogonal
to every even power (including the zeroth), that \( \mathbf{q}_k \) is an even or odd vector according to
the parity of \( k \), and explain the pattern of raw coefficients in @exm-ss-nile-polynomials.
:::

::: {.solution}
Let \( \mathbf{t} \) have entries \( t_i=x_i-\bar{x} \) and let \( \mathbf{J} \) be the reversal permutation, so that
\( \mathbf{J}\mathbf{t}=-\mathbf{t} \). A vector \( \bu \) is even if \( \mathbf{J}\bu=\bu \) and odd if \( \mathbf{J}\bu=-\bu \). Powers
\( \mathbf{t}^k \) have the parity of \( k \), and an even vector \( \bu \) and an odd vector \( \bv \) are
orthogonal because \( \bu\T\bv=(\mathbf{J}\bu)\T(\mathbf{J}\bv)=-\bu\T\bv \). Gram–Schmidt on
\( \bone,\mathbf{t},\mathbf{t}^2,\dots \) subtracts from \( \mathbf{t}^k \) only its projections on the earlier
\( \mathbf{q}_j \) of the same parity, by induction, so \( \mathbf{q}_k \) has the parity of \( k \). Hence the
column \( \mathbf{t} \) is orthogonal to the column \( \mathbf{t}^2 \) and to \( \bone \), and adding the quadratic
term does not change the linear coefficient (@cor-ss-orthogonal-coefficients, applied to the terms
\( \mathbf{t} \) and \( \mathbf{t}^2-\bar{\mathbf{t}^2}\bone \)). The cubic column is odd and not orthogonal to \( \mathbf{t} \),
so it changes the linear coefficient. The quartic column is even. By @thm-ss-adding(b), adding it
changes the cubic fit's coefficients by \( -\bL\hat{\gamma} \) with
\( \bL=(\X\T\X)^{-1}\X\T\mathbf{t}^4 \), where \( \X=[\bone,\mathbf{t},\mathbf{t}^2,\mathbf{t}^3] \). Inner products
between odd and even columns vanish, so \( \X\T\X \) and its inverse do not link the odd and even
columns, and the entries of \( \X\T\mathbf{t}^4 \) for the odd columns are zero. Hence the rows of \( \bL \) for
\( \mathbf{t} \) and \( \mathbf{t}^3 \) vanish, and the linear coefficient is unchanged.
:::

