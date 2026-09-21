# Quantile regression

Least squares answers a question about a loss function: which constant minimizes
expected squared error? Change the loss and the answer changes. The median minimizes
expected absolute error, and one modification of absolute error — weighting the two
sides differently — picks out any quantile at all. Koenker and Bassett (1978) built a
regression method on that observation.

## The check loss

::: {#def-qnt-quantile}
[Check loss and quantile regression]

For \( \tau\in(0,1) \) the **check function** is
\[
\rho_{\tau}(u)=u\bigl\{\tau-1\{u<0\}\bigr\}
=\begin{cases}\tau u,&u\ge0,\\(\tau-1)u,&u<0,\end{cases}
\]{#eq-qnt-check}

a convex, piecewise linear, nonnegative function with a kink at the origin, equal to
\( \lvert u\rvert/2 \) when \( \tau=1/2 \). The **linear quantile regression model** at
level \( \tau \) asserts
\[
Q_Y(\tau\mid\x)=\x\T\bbeta_{\tau},
\]
and given data \( (y_i,\x_{(i)}) \), \( i=1,\dots,n \), a \( \tau \)**-regression
quantile** is any minimizer
\[
\hbeta_{\tau}\in\argmin_{\bb\in\Real^p}\;
   R_{\tau}(\bb)=\sum_{i=1}^n\rho_{\tau}\bigl(y_i-\x_{(i)}\T\bb\bigr).
\]{#eq-qnt-sample}

:::

The subscript \( \tau \) is there because the coefficients may differ from level to
level: the model is a separate statement for each \( \tau \), not one model with a
nuisance parameter. Nothing is assumed about the shape of the conditional
distribution, not even that it has a density.

## What the check loss estimates

::: {#thm-qnt-check}
[The check loss elicits the quantile]

Let \( Y \) be a random variable with distribution function \( F \) and
\( \E\lvert Y\rvert<\infty \), and let \( G(q)=\E\rho_{\tau}(Y-q) \).

::: {.enumerate options="label=(\alph*)"}
1. **Population characterization.** \( G \) is finite and convex, and for \( q<q' \)
   \[
   G(q')-G(q)=\int_{q}^{q'}\{F(s)-\tau\}\,ds .
   \]{#eq-qnt-difference}

   Consequently the set of minimizers of \( G \) is exactly the set of
   \( \tau \)-quantiles of \( F \), that is
   \( \{q:F(q^{-})\le\tau\le F(q)\} \).

2. **Linear program.** Writing \( u_i=(y_i-\x_{(i)}\T\bb)_{+} \) and
   \( v_i=(\x_{(i)}\T\bb-y_i)_{+} \), problem @eq-qnt-sample is equivalent to
   \[
   \min_{\bb,\bu,\bv}\;\bigl\{\tau\bone\T\bu+(1-\tau)\bone\T\bv\bigr\}
   \quad\text{subject to}\quad \X\bb+\bu-\bv=\y,\ \ \bu\ge\bzero,\ \bv\ge\bzero,
   \]{#eq-qnt-lp}

   a linear program in \( 2n+2p \) nonnegative variables with \( n \) equality
   constraints.

3. **Equivariance of the estimate.** Let \( \hbeta_{\tau}(\y,\X) \) denote any
   solution of @eq-qnt-sample. For \( a>0 \), \( \bgamma\in\Real^p \) and nonsingular
   \( \A \),
   \[
   \begin{gathered}
   \hbeta_{\tau}(a\y,\X)=a\,\hbeta_{\tau}(\y,\X),\qquad
   \hbeta_{\tau}(-\y,\X)=-\hbeta_{1-\tau}(\y,\X),\\
   \hbeta_{\tau}(\y+\X\bgamma,\X)=\hbeta_{\tau}(\y,\X)+\bgamma,\\
   \hbeta_{\tau}(\y,\X\A)=\A^{-1}\hbeta_{\tau}(\y,\X).
   \end{gathered}
   \]

4. **Equivariance under monotone maps.** For strictly increasing left-continuous
   \( h \), \( Q_{h(Y)}(\tau\mid\x)=h\{Q_Y(\tau\mid\x)\} \).
:::

:::

::: {.proof}
(a) Split the loss at the kink: \( \rho_{\tau}(u)=\tau u_{+}+(1-\tau)u_{-} \), where
\( u_{-}=(-u)_{+} \), so
\( G(q)=\tau\E(Y-q)_{+}+(1-\tau)\E(q-Y)_{+} \), finite because both terms are bounded
by \( \E\lvert Y\rvert+\lvert q\rvert \). By Fubini's theorem,
\[
\E(Y-q)_{+}=\int_{q}^{\infty}\{1-F(s)\}\,ds,\qquad
\E(q-Y)_{+}=\int_{-\infty}^{q}F(s)\,ds .
\]
Subtracting the two expressions for \( q' \) and \( q \) gives
\[
\begin{gathered}
G(q')-G(q)=-\tau\int_{q}^{q'}\{1-F(s)\}\,ds+(1-\tau)\int_{q}^{q'}F(s)\,ds\\
=\int_{q}^{q'}\{F(s)-\tau\}\,ds,
\end{gathered}
\]
which is @eq-qnt-difference; convexity follows because the integrand is
nondecreasing. Let \( q_{\tau}=\inf\{q:F(q)\ge\tau\} \). If \( q<q_{\tau} \) then
\( F(s)<\tau \) for every \( s<q_{\tau} \), so
\( G(q)-G(q_{\tau})=\int_{q}^{q_{\tau}}\{\tau-F(s)\}\,ds>0 \) and \( q \) is not a
minimizer. If \( q>q_{\tau} \) then \( F(s)\ge F(q_{\tau})\ge\tau \) for \( s\ge q_{\tau} \)
by right continuity, so \( G(q)\ge G(q_{\tau}) \), with equality if and only if
\( F\equiv\tau \) almost everywhere on \( (q_{\tau},q) \), that is \( F(q^{-})=\tau \).
Hence \( q_{\tau} \) is a minimizer and the minimizers are exactly the points with
\( F(q^{-})\le\tau\le F(q) \).

(b) For any \( \bb \) the choices \( u_i=(y_i-\x_{(i)}\T\bb)_{+} \) and
\( v_i=(\x_{(i)}\T\bb-y_i)_{+} \) are feasible with objective value
\( R_{\tau}(\bb) \), so the linear program's optimum is at most
\( \min_{\bb}R_{\tau}(\bb) \). Conversely, if \( (\bb,\bu,\bv) \) is feasible then
\( u_i-v_i=y_i-\x_{(i)}\T\bb \), and for a fixed difference the sum
\( \tau u_i+(1-\tau)v_i \) is minimized by \( \min(u_i,v_i)=0 \), which gives exactly
\( \rho_{\tau}(y_i-\x_{(i)}\T\bb) \). The two problems therefore have the same value
and the same optimal \( \bb \).

(c) For \( a>0 \), \( \rho_{\tau}(au)=a\rho_{\tau}(u) \), so
\( R_{\tau}(\bb;a\y)=aR_{\tau}(\bb/a;\y) \) and the minimizers correspond under
\( \bb\mapsto a\bb \). A direct check on @eq-qnt-check gives
\( \rho_{\tau}(-u)=\rho_{1-\tau}(u) \), whence
\( R_{\tau}(\bb;-\y)=R_{1-\tau}(-\bb;\y) \) and the second identity. For the shift,
\( y_i+\x_{(i)}\T\bgamma-\x_{(i)}\T\bb=y_i-\x_{(i)}\T(\bb-\bgamma) \), so the objective
is the old one at \( \bb-\bgamma \). For the reparameterization, the residuals of
\( \X\A \) at \( \bb \) are those of \( \X \) at \( \A\bb \).

(d) This is @lem-qnt-monotone applied conditionally on \( \x \).
:::

Part (a) is the whole justification of the method, and it needs remarkably little: no
density, no moments beyond the first, no continuity. Part (d) is a property least
squares does not have: the median of \( \log Y \) is the logarithm of the median of
\( Y \), while the mean of \( \log Y \) is not the logarithm of the mean, the
retransformation problem of
[Section 22.5](../ch22-transformations/05-interpreting.html). A quantile regression
fitted on the log scale reads back without a smearing correction.

::: {.remark}
[The influence of one observation]

Away from the kink \( \rho_{\tau}'(u)=\tau-1\{u<0\} \) takes only the values
\( \tau \) and \( \tau-1 \), so an observation enters the estimating equation *only
through the sign of its residual*: moving one response from \( +10 \) to
\( +10^{6} \) changes nothing. This is the bounded-influence property that made the
median robust in
[Section 20.6](../ch20-residuals-leverage-influence/06-robust.html). It buys nothing
against a bad *design* point, whose \( \x_{(i)} \) enters unbounded.
:::

::: {.remark}
[Without moments]

Minimizing \( G_0(q)=\E\{\rho_{\tau}(Y-q)-\rho_{\tau}(Y)\} \) instead removes the
moment condition: since
\( \lvert\rho_{\tau}(u-v)-\rho_{\tau}(u)\rvert\le\max(\tau,1-\tau)\lvert v\rvert \),
this is finite whatever the tails of \( Y \), and differs from \( G \) by a constant
whenever \( G \) is finite, so (a) holds verbatim.
:::

## Solving the problem

::: {#prp-qnt-counts}
[Interpolation and counting]

Let \( \rank(\X)=p \) and let \( \hbeta_{\tau} \) minimize @eq-qnt-sample.

::: {.enumerate options="label=(\alph*)"}
1. There is a minimizer whose residual vector vanishes at \( p \) observations whose
   rows of \( \X \) are linearly independent: the fit passes exactly through \( p \)
   of the data points.

2. If \( \bone\in\C(\X) \), then for *any* minimizer, with
   \( n_{-},n_{0},n_{+} \) the numbers of negative, zero and positive residuals,
   \[
   n_{-}\;\le\;\tau n\;\le\;n_{-}+n_{0}.
   \]{#eq-qnt-counting}

:::

:::

::: {.proof}
(a) Let \( \bb \) be any minimizer, \( H=\{i:y_i=\x_{(i)}\T\bb\} \) the interpolated
set, and suppose \( \{\x_{(i)}:i\in H\} \) spans a space of dimension less than
\( p \). Choose \( \mathbf{d}\ne\bzero \) orthogonal to all of them. For small
\( t \) the residuals in \( H \) stay zero and no other residual changes sign, so
\[
R_{\tau}(\bb+t\mathbf{d})
=\sum_{i\notin H}\bigl(r_i-t\,\x_{(i)}\T\mathbf{d}\bigr)
   \bigl\{\tau-1\{r_i<0\}\bigr\},\qquad r_i=y_i-\x_{(i)}\T\bb,
\]
which is affine in \( t \). An affine function minimized at \( t=0 \) is constant, so
\( R_{\tau} \) stays minimal as \( t \) grows, until some residual outside \( H \)
reaches zero. Because \( \rank(\X)=p \) and \( \mathbf{d}\ne\bzero \), some
\( \x_{(i)}\T\mathbf{d}\ne0 \), so such a \( t \) exists and is finite; there the
interpolated set is strictly larger and the objective still minimal. Repeating raises
the rank of \( \{\x_{(i)}:i\in H\} \) to \( p \).

(b) Let \( \bone=\X\mathbf{c} \) and consider \( \phi(t)=R_{\tau}(\hbeta_{\tau}+t\mathbf{c}) \),
which shifts every fitted value by \( t \) and so every residual by \( -t \). For
\( t>0 \) small, an observation with \( r_i>0 \) contributes \( \tau(r_i-t) \), one with
\( r_i\le0 \) contributes \( (\tau-1)(r_i-t) \); hence
\( \phi'(0^{+})=-\tau n_{+}+(1-\tau)(n_{-}+n_{0}) \). Similarly
\( \phi'(0^{-})=-\tau(n_{+}+n_{0})+(1-\tau)n_{-} \). Minimality forces
\( \phi'(0^{+})\ge0 \) and \( \phi'(0^{-})\le0 \). The first gives
\( \tau n_{+}\le(1-\tau)(n_{-}+n_{0}) \); adding
\( \tau(n_{-}+n_{0}) \) to both sides and using \( n_{+}+n_{-}+n_{0}=n \) yields
\( \tau n\le n_{-}+n_{0} \). The second gives
\( \tau(n_{+}+n_{0})\ge(1-\tau)n_{-} \), and adding \( \tau n_{-} \) to both sides
yields \( \tau n\ge n_{-} \).
:::

Part (b) says the fitted surface has about \( \tau n \) observations below it, which
is what the name promises; part (a) says it is pinned to \( p \) of the data points,
as the median of an odd number of points is one of the points. Together they make the
fit a combinatorial object: finitely many candidate hyperplanes, one per choice of
\( p \) rows, and the optimum is one of them. Hence the simplex method of Koenker and
d'Orey (1987) rather than differentiation, or the interior-point route of Portnoy and
Koenker (1997), which makes quantile regression competitive in cost with least
squares for large \( n \).

::: {#exm-qnt-lp}
[Seven quantile lines for Engel's budgets]

Fitting @eq-qnt-lp at seven levels gives slopes rising steadily with \( \tau \):

| \( \tau \) | 0.05 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.95 |
|---|---|---|---|---|---|---|---|
| slope | 0.3434 | 0.4018 | 0.4741 | 0.5602 | 0.6440 | 0.6863 | 0.7091 |
| intercept | 124.9 | 110.1 | 95.5 | 81.5 | 62.4 | 67.4 | 64.1 |

The least squares slope, \( 0.4852 \), sits between the quartiles. A household at
the tenth percentile of its income group spends about \( 0.4018 \) of each extra
franc on food, one at the ninetieth \( 0.6863 \); at an income of \( 1000 \) francs
the fitted median budget share is \( 0.642 \) and the tenth-percentile share
\( 0.512 \). Every fit interpolates \( 2 \) observations and
satisfies @eq-qnt-counting, as the listings verify.
:::

```{.python .run #cell-quantile-engel-lp}
import numpy as np
import statsmodels.api as sm
from scipy.optimize import linprog
data = sm.datasets.engel.load_pandas().data
TAUS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

income = data["income"].to_numpy()
food = data["foodexp"].to_numpy()
n = len(food)
X = np.column_stack([np.ones(n), income])

def check_loss(u, tau):
    """The check function rho_tau applied elementwise."""
    return u * (tau - (u < 0))

def qreg(X, y, tau):
    """Quantile regression by linear programming.

    Variables are (b+, b-, u, v), all nonnegative, with beta = b+ - b- and
    u - v = y - X beta the positive and negative parts of the residual.
    """
    n, p = X.shape
    cost = np.concatenate([np.zeros(2 * p), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    out = linprog(cost, A_eq=A, b_eq=y, bounds=(0, None), method="highs")
    return out.x[:p] - out.x[p:2 * p]

beta = {tau: qreg(X, food, tau) for tau in TAUS}
for tau in (0.10, 0.50, 0.90):
    print(f"tau = {tau:.2f}:  intercept {beta[tau][0]:8.2f}   slope {beta[tau][1]:.4f}")
```

```{.python .run #cell-quantile-engel-counts}
for tau in TAUS:
    r = food - X @ beta[tau]
    n_zero = np.sum(np.abs(r) < 1e-8)
    n_neg = np.sum(r < -1e-8)
    assert n_zero == X.shape[1]                      # p residuals vanish exactly
    assert n_neg <= tau * n <= n_neg + n_zero        # the counting identity
print("at every tau: 2 residuals are zero and #negative <= tau*n <= #negative + 2")
```

::: {when-format="html"}
![**Figure 45.2.1.** (a) Seven fitted quantile lines for Engel's budgets, least
squares dashed. (b) The fitted slope against \( \tau \).](engel_quantiles.svg){#fig-qnt-engel-quantiles width=100%}
:::

::: {when-format="pdf"}
![(a) Seven fitted quantile lines for Engel's budgets, least
squares dashed. (b) The fitted slope against \( \tau \).](engel_quantiles.pdf){width=100%}
:::

Koenker and Machado (1999) proposed the analogue of \( R^2 \): with
\( \hat V(\tau) \) the minimized check loss and \( \tilde V(\tau) \) its value in the
intercept-only model, \( R_1(\tau)=1-\hat V(\tau)/\tilde V(\tau) \). For Engel's data
it rises from \( 0.474 \) at \( \tau=0.05 \) to \( 0.805 \) at \( \tau=0.95 \):
income explains far more about the upper tail than about the lower. Like \( R^2 \) it
says nothing about whether the linear specification is right.

## Reading a coefficient

The coefficient \( \beta_{\tau j} \) is the derivative of the \( \tau \)-quantile of
the conditional distribution with respect to \( x_j \). It is *not* the effect of
\( x_j \) on the household that happens to sit at rank \( \tau \).

::: {.remark}
[Quantiles of distributions, not of people]

Suppose incomes rise by one franc for every household. The \( 0.9 \)-quantile of the
new conditional distribution exceeds the old by \( \beta_{0.9} \); but the household
at rank \( 0.9 \) before need not be there after, and nothing in the data says what
happened to it. Reading \( \beta_{\tau} \) as an individual effect assumes the change
preserves ranks — the kind of untestable structure
[Chapter 25](../ch25-causal-interpretation/index.html) insists on making
explicit (@def-cau-potential-outcomes). Without it, quantile regression describes how
a *distribution* moves, which is usually the question anyway.
:::

Two cautions carry over: extrapolation is as unsafe as in
[Section 12.3](../ch12-intervals-and-bands/03-prediction.html), and
[Section 45.3](03-inference.html) shows it is the first place quantile curves
misbehave; and linearity of \( Q_Y(\tau\mid\x) \) is an assumption, testable in the
spirit of @thm-cor-lack-of-fit.

## Exercises

### A. Check your understanding

::: {#exr-qnt-count-example}
[A1]

A quantile regression with an intercept and one regressor is fitted at
\( \tau=0.25 \) to \( n=40 \) observations. How many residuals can be negative, and how
many are exactly zero? Use @prp-qnt-counts.
:::

::: {.solution}
Generic data give \( n_0=p=2 \), and @eq-qnt-counting reads
\( n_{-}\le10\le n_{-}+2 \), so \( n_{-}\in\{8,9,10\} \).
:::

### B. Practice

::: {#exr-qnt-sample-quantile}
[B1]

Take \( \X=\bone \), so that the problem is to fit one number. Show that the
minimizers of \( \sum_i\rho_{\tau}(y_i-b) \) are exactly the \( \tau \)-quantiles of the
empirical distribution of \( y_1,\dots,y_n \), and that for \( \tau n \) not an integer
the minimizer is the unique order statistic \( y_{(\lceil\tau n\rceil)} \).
:::

::: {.solution}
Apply @thm-qnt-check(a) to the empirical distribution \( F_n \), which is legitimate
because \( F_n \) is a distribution function with finite mean. The minimizers are
\( \{b:F_n(b^{-})\le\tau\le F_n(b)\} \). If \( \tau n\notin\mathbb{Z} \), no order
statistic has \( F_n \) equal to \( \tau \), and the condition picks out the single
point \( y_{(\lceil\tau n\rceil)} \); if \( \tau n \) is an integer, the whole interval
\( [y_{(\tau n)},y_{(\tau n+1)}] \) minimizes.
:::

::: {#exr-qnt-lp-dual}
[B2]

Write down the dual of the linear program @eq-qnt-lp and show that it is
\[
\max_{\mathbf{d}}\;\y\T\mathbf{d}\quad\text{subject to}\quad
\X\T\mathbf{d}=\bzero,\ \ \mathbf{d}\in[\tau-1,\tau]^n,
\]
or, after the shift \( \mathbf{a}=\mathbf{d}+(1-\tau)\bone \),
\( \max\{\y\T\mathbf{a}:\X\T\mathbf{a}=(1-\tau)\X\T\bone,\ \mathbf{a}\in[0,1]^n\} \).
Show that at an optimum \( a_i=1 \) for a positive residual and \( a_i=0 \) for a
negative one, so that \( \mathbf{a} \) records the "signs" of the residuals, and
explain why these are the raw material of the rank tests of
[Section 45.3](03-inference.html).
:::

### C. Going deeper

::: {#exr-qnt-uniqueness}
[C1]

Show that the set of minimizers of \( R_{\tau} \) is a convex polyhedron, and that if
it contains two distinct points then \( R_{\tau} \) is constant on the segment joining
them. Deduce that a unique solution is the generic case: non-uniqueness requires a
coincidence in the data, such as \( p+1 \) observations lying on one hyperplane.
:::
