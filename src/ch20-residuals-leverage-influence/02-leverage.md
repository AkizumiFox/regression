# Leverage and the hat matrix

[Section 6.8](../ch06-projections/08-leverage.html) proved the basic properties of the leverage \( h_{ii} \): it lies in
\( [0,1] \), the leverages sum to \( p \), \( h_{ii}\ge1/n \) with an intercept, and \( h_{ii} \) is an affine function of the
Mahalanobis distance of \( \x_{(i)} \) from the centroid of the regressors (@prp-proj-leverage and @prp-proj-leverage-mahalanobis). This section asks when a leverage deserves attention, how designed and observational data differ,
and how to measure the pull of a case on one coefficient.

The fitted value of case \( i \) is
\[
\hat y_i=h_{ii}y_i+\sum_{j\ne i}h_{ij}y_j ,
\]
so \( h_{ii} \) is the weight of a case in its own fitted value. Since \( \sum_jh_{ij}^2=h_{ii} \) (@prp-proj-leverage(a)), a
case with \( h_{ii} \) near one leaves little weight for the others: its fitted value is essentially its own response, and
the fit carries no information about whether \( y_i \) is right. Leverage is *potential* influence, a property of the
design alone; [Section 20.4](04-influence.html) combines it with the residual.

## How large is large?

The average leverage is \( p/n \). The most widely used rule (Hoaglin and Welsch 1978) calls a case of high
leverage if \( h_{ii}>2p/n \); a stricter variant uses \( 3p/n \). To see what these rules mean, compute how often they
fire under a reference design with multivariate normal regressors.

::: {#prp-res-leverage-beta}
[Leverage under a Gaussian design]

Let \( \X=[\bone,\mathbf{Z}] \), where the rows of the \( n\times k \) matrix \( \mathbf{Z} \) are independent \( \Normal_k(\bmu,\bSigma) \)
vectors with \( \bSigma \) positive definite, and let \( n\ge k+2 \). Then \( \X \) has full column rank \( p=k+1 \) with
probability one, and for every \( i \),
\[
\frac{h_{ii}-1/n}{1-1/n}\sim\mathrm{Beta}\Bigl(\frac k2,\frac{n-1-k}2\Bigr).
\]
:::

::: {.proof}
*Reduction.* Let \( \mathbf{C}=\I-n^{-1}\bone\bone\T \). By @lem-proj-fwl-split with \( \X_2=\bone \), \( \M=n^{-1}\bone\bone\T+\bP \),
where \( \bP \) is the orthogonal projection onto \( \C(\mathbf{C}\mathbf{Z}) \), a subspace of \( \bone\perpc \). So
\( h_{ii}-1/n=\vect{e}_i\T\bP\vect{e}_i \). Replacing \( \mathbf{Z} \) by \( (\mathbf{Z}-\bone\bmu\T)\bSigma^{-1/2} \) does not change \( \C(\mathbf{C}\mathbf{Z}) \),
because \( \mathbf{C}\bone=\bzero \) and \( \bSigma^{-1/2} \) is nonsingular. So we may assume that the entries of \( \mathbf{Z} \) are
independent \( \Normal(0,1) \). Then \( \mathbf{C}\mathbf{Z} \) has rank \( k \) with probability one (@exr-res-gaussian-rank), and so
does \( \X \).

*Invariance.* Let \( \bv=\mathbf{C}\vect{e}_i/\norm{\mathbf{C}\vect{e}_i} \), a unit vector in \( \bone\perpc \), with
\( \norm{\mathbf{C}\vect{e}_i}^2=1-1/n \). Since \( \bP=\bP\mathbf{C}=\mathbf{C}\bP \),
\[
h_{ii}-\frac1n=(\mathbf{C}\vect{e}_i)\T\bP(\mathbf{C}\vect{e}_i)=\Bigl(1-\frac1n\Bigr)\bv\T\bP\bv .
\]
Write \( \bP(\mathbf{Z}) \) to show the dependence on \( \mathbf{Z} \). Let \( \bw \) be any other unit vector in \( \bone\perpc \), and let
\( \mathbf{H}=\I-2\mathbf{a}\mathbf{a}\T/\mathbf{a}\T\mathbf{a} \) with \( \mathbf{a}=\bv-\bw \) (or \( \mathbf{H}=\I \) if \( \bw=\bv \)). This is a symmetric orthogonal
matrix, and \( \mathbf{H}\bv=\bw \) because \( \norm{\bv}=\norm{\bw} \). It fixes \( \bone \), because \( \bone\perp\mathbf{a} \), so \( \mathbf{H} \) commutes with \( \mathbf{C} \) and maps
\( \C(\mathbf{C}\mathbf{Z}) \) onto \( \C(\mathbf{C}\mathbf{H}\mathbf{Z}) \). Hence \( \bP(\mathbf{H}\mathbf{Z})=\mathbf{H}\bP(\mathbf{Z})\mathbf{H} \), and
\( \bw\T\bP(\mathbf{Z})\bw=\bv\T\mathbf{H}\bP(\mathbf{Z})\mathbf{H}\bv=\bv\T\bP(\mathbf{H}\mathbf{Z})\bv \). Each column of \( \mathbf{H}\mathbf{Z} \) is
\( \Normal_n(\bzero,\mathbf{H}\mathbf{H}\T)=\Normal_n(\bzero,\I) \), and the columns are independent, so \( \mathbf{H}\mathbf{Z} \) has the law of \( \mathbf{Z} \).
Therefore the law of \( \bv\T\bP(\mathbf{Z})\bv \) is the same for every unit vector \( \bv\in\bone\perpc \).

*Averaging.* Let \( \mathbf{g}\sim\Normal_n(\bzero,\mathbf{C}) \) be independent of \( \mathbf{Z} \), and \( \mathbf{U}=\mathbf{g}/\norm{\mathbf{g}} \), a random unit
vector in \( \bone\perpc \). By the previous paragraph, applied conditionally on \( \mathbf{U} \), the variable
\( \mathbf{U}\T\bP(\mathbf{Z})\mathbf{U} \) has the same law as \( \bv\T\bP(\mathbf{Z})\bv \). Now condition on \( \mathbf{Z} \) instead. Then \( \bP \) is a fixed
projection of rank \( k \) onto a subspace of \( \bone\perpc \), and \( \mathbf{C}-\bP \) is a projection of rank \( n-1-k \) orthogonal to it,
so by @thm-qf-orthogonal-projections \( \norm{\bP\mathbf{g}}^2\sim\chi^2(k) \) and \( \norm{(\mathbf{C}-\bP)\mathbf{g}}^2\sim\chi^2(n-1-k) \)
are independent. Since \( \mathbf{g}=\mathbf{C}\mathbf{g} \),
\[
\mathbf{U}\T\bP\mathbf{U}=\frac{\norm{\bP\mathbf{g}}^2}{\norm{\bP\mathbf{g}}^2+\norm{(\mathbf{C}-\bP)\mathbf{g}}^2}
\sim\mathrm{Beta}\Bigl(\frac k2,\frac{n-1-k}2\Bigr)
\]
by @exr-qf-beta. This conditional law does not depend on \( \mathbf{Z} \), so it is also the unconditional law.
:::

The mean of the beta law is \( k/(n-1) \), which gives \( \E(h_{ii})=1/n+(1-1/n)k/(n-1)=p/n \), as it must. As
\( n\to\infty \) with \( k \) fixed, \( n(h_{ii}-1/n) \) converges in law to \( \chi^2(k) \) (@exr-res-leverage-limit), so
\[
\Pr\{h_{ii}>c\,p/n\}\approx\Pr\{\chi^2(k)>c(k+1)-1\}.
\]
With \( k=3 \), the rule \( 2p/n \) fires for \( 0.072 \) of the cases in large samples (\( 0.067 \) when \( n=50 \)), and
the rule \( 3p/n \) for \( 0.012 \). With \( k=10 \) and \( n=500 \), the rule \( 2p/n \) fires for only \( 0.020 \). "Twice the
average" flags a different share of well-behaved cases for every \( p \): it is a convention, not a test.
[Figure 20.2.1](#fig-res-leverage-design)(a) checks the law against \( 4000 \) simulated designs.

```{.python .run #cell-leverage-design-beta}
import numpy as np
from scipy import stats

def leverages(X):
    Q, _ = np.linalg.qr(X)
    return np.sum(Q ** 2, axis=1)

rng = np.random.default_rng(2020)
n, k = 30, 3                                       # intercept plus k Gaussian regressors
p = k + 1
reps = 4000
h1 = np.empty(reps)
for r in range(reps):
    Z = rng.normal(size=(n, k)) @ np.array([[1.0, 0.0, 0.0], [0.8, 0.6, 0.0], [0.0, 2.0, 1.0]])
    X = np.column_stack([np.ones(n), Z + 5.0])     # correlated, shifted regressors
    h1[r] = leverages(X)[0]                        # leverage of case 1
B = (h1 - 1 / n) / (1 - 1 / n)                     # should be Beta(k/2, (n-1-k)/2)
law = stats.beta(k / 2, (n - 1 - k) / 2)
print("mean of h_11:", h1.mean(), " theory p/n =", p / n)
print("Kolmogorov-Smirnov p-value against the Beta law:", stats.kstest(B, law.cdf).pvalue)
```

## Designed and observational data

In a designed experiment the analyst chooses the rows of \( \X \). The largest leverage is at least the average
\( p/n \), with equality iff all are equal (@exr-res-equal-leverage), and many designs achieve it. If the columns are
orthogonal with entries \( \pm1 \), as in two-level factorials, then \( \X\T\X=n\I \) and \( h_{ii}=\norm{\x_{(i)}}^2/n=p/n \): a
\( 2^3 \) factorial run twice (\( n=16 \)) with intercept, main effects and one interaction (\( p=5 \)) gives every case
leverage \( 5/16=0.3125 \). A balanced one-way layout gives every case \( 1/n_k=p/n \) (@exr-proj-leverage-monotone). Kiefer and Wolfowitz (1960) showed
that, over a region of possible design points, the designs maximizing \( \det(\X\T\X) \) are exactly those minimizing the
largest leverage \( \x\T(\X\T\X)^{-1}\x \) over the region, with minimum \( p/n \) (for designs viewed as probability
measures). For D-optimal designs, good design and small maximum leverage are the same thing.

In observational data the rows are drawn, not chosen, and with light enough tails no case dominates in large samples.

::: {#prp-res-max-leverage}
[Leverage vanishes in large samples]

Let the rows \( \x_{(1)},\x_{(2)},\dots \) be independent copies of a random vector \( \x \) with \( \E\norm{\x}^2<\infty \) and
\( \E(\x\x\T) \) positive definite, and let \( h_{ii}^{(n)} \) be the leverages of the first \( n \) rows. Then
\( \max_{i\le n}h_{ii}^{(n)}\to0 \) in probability.
:::

::: {.proof}
Write \( \mathbf{S}_n=n^{-1}\X\T\X=n^{-1}\sum_{i\le n}\x_{(i)}\x_{(i)}\T \) and let \( \lambda_n \) be its smallest eigenvalue. By
the strong law of large numbers \( \mathbf{S}_n\to\E(\x\x\T) \) with probability one, and eigenvalues are continuous
functions of a symmetric matrix, so \( \lambda_n\to\lambda>0 \), the smallest eigenvalue of \( \E(\x\x\T) \). Hence with probability one \( \X\T\X \) is nonsingular
for all large \( n \), and the leverages are defined from then on. By
@thm-mat-extremal-rayleigh applied to \( (\X\T\X)^{-1} \), whose largest eigenvalue is \( 1/(n\lambda_n) \),
\[
\max_{i\le n}h_{ii}^{(n)}=\max_{i\le n}\x_{(i)}\T(\X\T\X)^{-1}\x_{(i)}\le\frac{\max_{i\le n}\norm{\x_{(i)}}^2}{n\lambda_n}.
\]
It remains to show that \( \max_{i\le n}\norm{\x_{(i)}}^2/n\to0 \) in probability. For \( \epsilon>0 \), by the union bound
and Markov's inequality,
\[
\Pr\Bigl\{\max_{i\le n}\norm{\x_{(i)}}^2>\epsilon n\Bigr\}\le n\Pr\{\norm{\x}^2>\epsilon n\}
\le\frac1\epsilon\,\E\bigl[\norm{\x}^2\,1\{\norm{\x}^2>\epsilon n\}\bigr]\to0,
\]
by dominated convergence, since \( \E\norm{\x}^2<\infty \).
:::

The condition matters for inference with non-normal errors. Huber (1973) showed that every fitted value
\( \x\T\hbeta \) is asymptotically normal, for every error law with finite variance, iff the largest leverage tends to zero,
and @thm-dep-nonnormal of [Chapter 19](../ch19-theory-of-departures/index.html) rests on a condition of this Lindeberg type.
A case that keeps a fixed share of the leverage keeps a fixed share of its own error in its fitted value, and no
averaging makes that error normal.

[Figure 20.2.1](#fig-res-leverage-design)(b) shows the median, over \( 300 \) replications, of the largest leverage with two
independent regressors. For normal regressors it falls from \( 0.322 \) at \( n=25 \) to \( 0.0030 \) at \( n=6400 \). For
\( t_3 \) regressors, heavy-tailed but with finite variance, it falls much more slowly, to \( 0.066 \). For Cauchy regressors,
which have no variance, it does not fall at all (\( 0.857 \) at \( n=25 \), \( 0.831 \) at \( n=6400 \)): one case keeps most of
the pull however many are added.

::: {when-format="html"}
![**Figure 20.2.1.** (a) Leverage of one case in \( 4000 \) Gaussian designs (\( n=30 \), \( p=4 \)), the beta density of
@prp-res-leverage-beta, and the cut-offs \( 2p/n \) (dashed) and \( 3p/n \) (dotted). (b) Median largest leverage for normal, \( t_3 \)
and Cauchy regressors.](leverage_design.svg){#fig-res-leverage-design width=100%}
:::

::: {when-format="pdf"}
![(a) Leverage of one case in \( 4000 \) Gaussian designs (\( n=30 \), \( p=4 \)), the beta density of
@prp-res-leverage-beta, and the cut-offs \( 2p/n \) (dashed) and \( 3p/n \) (dotted). (b) Median largest leverage for normal, \( t_3 \)
and Cauchy regressors.](leverage_design.pdf){width=100%}
:::

::: {#exm-res-engel-leverage}
[Leverage in Engel's household budgets]

Engel's 1857 survey records the annual income and food expenditure of \( 235 \) Belgian working-class households.
Income is strongly skewed (skewness \( 2.78 \)): the median is \( 884 \) and the largest, row \( 138 \) of the file, is
\( 4958 \). In the straight-line regression of food expenditure on income the average leverage is \( 2/235=0.0085 \), and
the richest household has leverage \( 0.255 \), about \( 30 \) times the average; the next largest is \( 0.058 \). Since
\( h_{ii}-1/n=(x_i-\bar x)^2/S_{xx} \), this one household supplies a quarter (\( 0.251 \)) of the spread of income about its
mean, and the fitted line will pass near it whatever it spends on food. On the logarithmic scale its leverage is
\( 0.070 \), still the largest but no longer in a class of its own. Choosing the scale
([Chapter 22](../ch22-transformations/index.html)) decides which cases the fit depends on.
:::

```{.python .run #cell-engel-leverage-leverage}
import numpy as np
import statsmodels.api as sm
engel = sm.datasets.engel.load_pandas().data

income = engel["income"].to_numpy()
food = engel["foodexp"].to_numpy()
n = len(income)
for label, x in [("income", income), ("log income", np.log(income))]:
    X = np.column_stack([np.ones(n), x])
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)
    top = np.argsort(h)[::-1][:3]
    print(f"{label:11s} average p/n = {2 / n:.4f};  largest:",
          ", ".join(f"row {i + 1}: h = {h[i]:.3f}" for i in top))
```

Leverage also depends on the model: adding a column never decreases a leverage (@exr-proj-leverage-monotone), and a
case ordinary for a line can be extreme for a quadratic (@exr-res-quadratic-leverage). A new point inside the range of
every regressor can have high leverage because its *combination* of values is unusual (@exm-ci-hidden-extrapolation).

## Partial leverage

The leverage \( h_{ii} \) measures the pull of case \( i \) on the whole fitted surface. The Frisch–Waugh–Lovell theorem
isolates one coefficient. Let \( \M_{(j)} \) project onto the span of the columns other than \( \x_j \), and let
\( \tilde{\x}_j=(\I-\M_{(j)})\x_j \).

::: {#def-res-partial-leverage}
[Partial leverage]

The **partial leverage** of case \( i \) for the coefficient of \( \x_j \) is
\[
\ell_{ij}=\frac{\tilde x_{ij}^2}{\norm{\tilde{\x}_j}^2}.
\]
:::

::: {#prp-res-partial-leverage}
[Properties of partial leverage]

Let \( h_{ii}^{(j)} \) be the leverages of the model without \( \x_j \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( h_{ii}=h_{ii}^{(j)}+\ell_{ij} \) for every \( i \);

2. \( 0\le\ell_{ij}\le1 \) and \( \sum_i\ell_{ij}=1 \);

3. \( \ell_{ij} \) is the leverage of case \( i \) in the regression through the origin of \( \tilde{\y}=(\I-\M_{(j)})\y \) on
   \( \tilde{\x}_j \), the regression drawn in the added-variable plot for \( \x_j \);

4. \( \hat{\beta}_j=\sum_i\bigl(\tilde x_{ij}/\norm{\tilde{\x}_j}^2\bigr)y_i \), so the weight of \( y_i \) in \( \hat{\beta}_j \) has
   square \( \ell_{ij}/\norm{\tilde{\x}_j}^2=\ell_{ij}\Var(\hat{\beta}_j)/\sigma^2 \).
:::

:::

::: {.proof}
By @lem-proj-fwl-split(c) with \( \X_1=\x_j \),
\( \M=\M_{(j)}+\tilde{\x}_j\tilde{\x}_j\T/\norm{\tilde{\x}_j}^2 \). The diagonal entries give (a), and the second term is the hat
matrix of the regression in (c). It is a rank-one projection, which gives (b). Part (d) is @thm-proj-fwl(a),
together with \( \Var(\hat{\beta}_j)=\sigma^2/\norm{\tilde{\x}_j}^2 \) (@eq-proj-vif-preview).
:::

A case can have high leverage and negligible partial leverage for the coefficient of interest, when its unusual
position lies along the other regressors. Partial leverages average \( 1/n \), so a few times \( 1/n \) is already large.

::: {#exm-res-state-partial}
[Where the District of Columbia's leverage comes from]

In the murder-rate regression for the \( 51 \) jurisdictions, with poverty, single-parent households and urbanization
as regressors, the District of Columbia has leverage \( 0.506 \) (@exm-proj-leverage) but partial leverage for poverty
only \( 0.032 \): its leverage comes from its extreme single-parent and urban percentages, not from its poverty rate. The
largest partial leverage for poverty is Alaska's, \( 0.141 \). [Section 20.4](04-influence.html) shows that the District
nevertheless moves the poverty coefficient more than any other case, because influence also depends on the residual.
:::

## Exercises

### A. Check your understanding

::: {#exr-res-origin-leverage}
[A1]

For regression through the origin on one regressor, show that \( h_{ii}=x_i^2/\sum_kx_k^2 \). Give data with \( n=10 \)
in which one leverage exceeds \( 0.99 \).
:::

::: {.solution}
\( \M=\x\x\T/\x\T\x \), so \( h_{ii}=x_i^2/\x\T\x \). With \( x_1=100 \) and \( x_2=\dots=x_{10}=1 \),
\( h_{11}=10^4/(10^4+9)>0.99 \).
:::

::: {#exr-res-affine-invariance}
[A2]

Show that with an intercept the leverages are unchanged when the regressors are replaced by \( \mathbf{a}+\B\x \) with \( \B \)
nonsingular, but not under a nonlinear transformation (@exm-res-engel-leverage). Use @prp-proj-leverage-mahalanobis.
:::

### B. Practice

::: {#exr-res-leverage-limit}
[B1]

In the setting of @prp-res-leverage-beta, show that \( n(h_{ii}-1/n) \) converges in law to \( \chi^2(k) \) as \( n\to\infty \).
:::

::: {.solution}
By the proof of @prp-res-leverage-beta, \( h_{ii}-1/n \) has the law of \( (1-1/n)\,U/(U+V) \) with \( U\sim\chi^2(k) \) and
\( V\sim\chi^2(n-1-k) \) independent. So \( n(h_{ii}-1/n) \) has the law of
\( (n-1)\,U/(U+V)=U\big/\bigl[(U+V)/(n-1)\bigr] \). By the law of large numbers \( V/(n-1)\to1 \) in probability, and
\( U/(n-1)\to0 \), so the denominator tends to \( 1 \) in probability. Slutsky's theorem gives convergence in law to
\( U\sim\chi^2(k) \).
:::

::: {#exr-res-equal-leverage}
[B2]

Show that \( \max_ih_{ii}\ge p/n \), with equality iff all leverages equal \( p/n \). Show that if \( \X \) has orthogonal
columns and every entry is \( \pm1 \), then all leverages equal \( p/n \).
:::


::: {#exr-res-quadratic-leverage}
[B3]

Take \( x=-3,-2,-1,0,1,2,3 \). Compute the leverage of the case at \( x=0 \) for the straight-line model and for the
quadratic model \( \E(Y)=\beta_0+\beta_1x+\beta_2x^2 \). Explain the difference.
:::

::: {.solution}
For the line, \( h=1/7+0^2/28=1/7 \). For the quadratic, apply @prp-res-partial-leverage(a) with \( \x_j \) the column of
squares: \( h=1/7+\ell \), where \( \ell \) is the partial leverage of the quadratic term. The residual of \( x^2 \) on
\( (1,x) \) is \( x^2-4 \), because \( x^2 \) has mean \( 4 \) and \( \sum_x x(x^2-4)=0 \) by symmetry. Its values are
\( 5,0,-3,-4,-3,0,5 \), with sum of squares \( 84 \). So \( \ell=16/84=4/21 \) and \( h=1/7+4/21=1/3 \). The case at \( x=0 \)
is the only one near the centre, and in a quadratic model the centre is where the curvature is pinned down.
:::

### C. Going deeper

::: {#exr-res-gaussian-rank}
[C1]

Complete the proof of @prp-res-leverage-beta by showing that if \( \mathbf{Z} \) is \( n\times k \) with independent \( \Normal(0,1) \)
entries and \( n\ge k+2 \), then \( \mathbf{C}\mathbf{Z} \) has rank \( k \) with probability one. *Hint:* add columns one at a time, and use
that a nonzero normal vector with a nonsingular covariance restricted to a subspace falls in a fixed proper subspace
with probability zero.
:::
