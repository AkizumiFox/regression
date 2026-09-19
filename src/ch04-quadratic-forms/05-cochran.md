# The Cochran theorem

An analysis of variance starts from an identity
\[
\Y\T\Y=\Y\T\A_1\Y+\dots+\Y\T\A_k\Y,\qquad \A_1+\dots+\A_k=\I,
\]
that splits the total sum of squares into pieces attributed to sources of variation.
@thm-qf-orthogonal-projections says the pieces are independent chi-squared variables
when the \( \A_i \) are mutually orthogonal projections. Checking that directly means checking
\( k \) idempotence conditions and \( \binom k2 \) orthogonality conditions. Cochran's theorem shows
that one count of ranks does all of it.

## The algebra

The matrix algebra behind this was done in [Chapter 1](../ch01-matrix-algebra/index.html). For symmetric
matrices, any two of the three properties “each \( \A_i \) is idempotent”,
“\( \A_i\A_j=\mathbf{0} \) for \( i\ne j \)” and “the sum is idempotent” imply the third, and
symmetric matrices that add to \( \I \) with ranks adding to \( n \) are idempotent and mutually
orthogonal (@thm-mat-idempotent-sum). The next theorem adds the rank count for a sum that
is idempotent but not necessarily \( \I \), which is what nested models need.

::: {#thm-qf-cochran-algebra}
[Cochran, algebraic form]

Let \( \A_1,\dots,\A_k \) be symmetric \( n\times n \) matrices and \( \A=\sum_i\A_i \). Consider

::: {.enumerate options="label=(\alph*)"}
1. each \( \A_i \) is idempotent;

2. \( \A_i\A_j=\mathbf{0} \) for all \( i\ne j \);

3. \( \A \) is idempotent;

4. \( \rank(\A)=\sum_i\rank(\A_i) \).
:::

Then any two of (a), (b), (c) imply all four, and (c) together with (d) implies all four.
:::

::: {.proof}
*Any two of (a), (b), (c).* By @thm-mat-idempotent-sum, any two of them imply the
third. When all three hold, rank equals trace for symmetric idempotent matrices
(@thm-mat-idempotent), so \( \rank(\A)=\tr(\A)=\sum_i\tr(\A_i)=\sum_i\rank(\A_i) \), which is (d).

*(c) and (d) imply (a) and (b).* Let \( r=\rank(\A) \) and \( r_i=\rank(\A_i) \). Every
\( \A\bv=\sum_i\A_i\bv \) lies in \( \C(\A_1)+\dots+\C(\A_k) \), whose dimension is at most
\( \sum_ir_i=r=\dim\C(\A) \). So the two spaces coincide, and in particular
\( \C(\A_i)\subseteq\C(\A) \). By @thm-mat-idempotent(a), \( \A=\Q\Q\T \) with \( \Q \) of size
\( n\times r \) and \( \Q\T\Q=\I_r \). Since \( \A \) is symmetric and idempotent it fixes every vector of
its column space, so \( \A\A_i=\A_i \), and transposing, \( \A_i\A=\A_i \). Hence
\( \A_i=\A\A_i\A=\Q\B_i\Q\T \) with \( \B_i=\Q\T\A_i\Q \). The \( r\times r \) matrices \( \B_i \) are symmetric,
satisfy \( \sum_i\B_i=\Q\T\A\Q=\I_r \), and have \( \rank(\B_i)=r_i \) because \( \Q \) has orthonormal
columns. By the last statement of @thm-mat-idempotent-sum, \( \B_i^2=\B_i \) and
\( \B_i\B_j=\mathbf{0} \) for \( i\ne j \), and multiplying by \( \Q \) on the left and \( \Q\T \) on the right
transfers both identities to the \( \A_i \).
:::

## The distributional statement

::: {#thm-qf-cochran}
[Cochran]

Let \( \Y\sim\Normal_n(\bmu,\sigma^2\I) \), and let \( \A_1,\dots,\A_k \) be symmetric with ranks
\( r_1,\dots,r_k \), whose sum \( \A=\sum_i\A_i \) is idempotent of rank \( r \) (for example \( \A=\I \),
\( r=n \)). The following are equivalent:

::: {.enumerate options="label=(\roman*)"}
1. the forms \( \Y\T\A_i\Y/\sigma^2 \) are mutually independent, and
           \( \Y\T\A_i\Y/\sigma^2\sim\chi^2(r_i,\bmu\T\A_i\bmu/\sigma^2) \) for each \( i \);

2. \( r_1+\dots+r_k=r \);

3. each \( \A_i \) is idempotent;

4. \( \A_i\A_j=\mathbf{0} \) for all \( i\ne j \).
:::

:::

::: {.proof}
Condition (c) of @thm-qf-cochran-algebra holds by assumption, so that theorem shows
that (ii), (iii) and (iv) are equivalent, and each implies the other two. If (iii) and (iv)
hold, @thm-qf-orthogonal-projections gives (i). If (i) holds, apply
@thm-qf-chisq-identity to \( \sigma^{-1}\Y\sim\Normal_n(\sigma^{-1}\bmu,\I) \). Each \( \A_i \)
with \( r_i\ge1 \) is idempotent, and an \( \A_i \) of rank zero is \( \mathbf{0} \). This is (iii).
:::

In words: *if the degrees of freedom add up, the sums of squares are independent
chi-squared variables.* Independence in (i) is not needed to reach (iii). The marginal
chi-squared laws alone force the whole structure.

::: {#exm-qf-oneway-cochran}
[The one-way layout]

Let observations fall into \( g \) groups of sizes \( n_1,\dots,n_g \), with \( Y_{\ell j} \) independent
\( \Normal(\mu_\ell,\sigma^2) \). Let \( \M \) replace each observation by its group mean, and
\( \bP_1=n^{-1}\bone\bone\T \) by the grand mean. Then
\[
\I=\bP_1+(\M-\bP_1)+(\I-\M)
\]
splits \( \norm{\Y}^2 \) into \( n\bar{Y}^2 \), the between-group sum of squares
\( \sum_\ell n_\ell(\bar{Y}_\ell-\bar{Y})^2 \) and the within-group sum of squares
\( \sum_{\ell,j}(Y_{\ell j}-\bar{Y}_\ell)^2 \). The three matrices are symmetric with ranks \( 1 \),
\( g-1 \) and \( n-g \), which add to \( n \). Cochran's theorem needs nothing more. The three sums of
squares, divided by \( \sigma^2 \), are independent noncentral chi-squared variables. Their
noncentralities are \( n\bar{\mu}^2/\sigma^2 \), \( \sum_\ell n_\ell(\mu_\ell-\bar{\mu})^2/\sigma^2 \) with
\( \bar{\mu}=\sum_\ell n_\ell\mu_\ell/n \), and \( 0 \). The last is zero because \( \M\bmu=\bmu \). The
ratio of the between and within mean squares is \( F(g-1,n-g,\gamma) \) by
@thm-qf-nested-f, and it is central iff all group means are equal.

The script uses \( g=3 \) groups of sizes \( 3, 5, 4 \), means
\( 2, 3, 1.5 \), and \( \sigma=1 \), and
draws \( 200{,}000 \) data sets. [Table 4.5.1](05-cochran.html#tab-qf-cochran) compares the theoretical means
\( r_i+\gamma_i \) and variances \( 2r_i+4\gamma_i \) with the simulated ones. The largest absolute
correlation between the three sums of squares is \( 0.0013 \). As a check that
goes beyond correlation, the proportion of data sets in which both the between and the within
sums of squares exceed their medians is \( 0.2497 \), against \( 0.25 \) for independent
variables.
:::

```{.python .run #cell-cochran-cochran}
import numpy as np
from scipy import stats
rng = np.random.default_rng(5150)

sizes, means, sigma = [3, 5, 4], [2.0, 3.0, 1.5], 1.0
groups = np.repeat(np.arange(3), sizes)
n = len(groups)
Z = (groups[:, None] == np.arange(3)).astype(float)    # group indicators
M = Z @ np.diag(1 / np.array(sizes)) @ Z.T             # group-mean projection
J = np.ones((n, n)) / n
pieces = [J, M - J, np.eye(n) - M]
ranks = [np.linalg.matrix_rank(P) for P in pieces]
print("ranks", ranks, "sum", sum(ranks), "n", n)

mu = np.array(means)[groups]
gammas = [mu @ P @ mu / sigma**2 for P in pieces]
Y = mu + sigma * rng.standard_normal((200_000, n))
Q = np.stack([np.einsum("ij,jk,ik->i", Y, P, Y) for P in pieces], axis=1) / sigma**2
print(np.corrcoef(Q, rowvar=False).round(3))
```

[**Table 4.5.1.** Cochran decomposition for a one-way layout with \( n=12 \): theory
against \( 200{,}000 \) simulated data sets.]{#tab-qf-cochran}

| Sum of squares | \( r_i \) | \( \gamma_i \) | mean | simulated | variance | simulated |
|---|---|---|---|---|---|---|
| grand mean | 1 | 60.750 | 61.750 | 61.763 | 245.00 | 244.95 |
| between groups | 2 | 5.250 | 7.250 | 7.230 | 25.00 | 24.85 |
| within groups | 9 | 0.000 | 9.000 | 8.994 | 18.00 | 17.96 |

::: {.remark}
[What the theorem does not cover]

Two assumptions carry the weight. The first is normality. Without it,
@thm-qf-mean-var(d) fails, and with independent non-normal errors the between and
within sums of squares are in general correlated, with a correlation driven by the fourth
cumulant (@exr-qf-kurtosis). The second is the spherical covariance \( \sigma^2\I \). With
\( \Cov(\Y)=\sigma^2\V \), \( \V \) positive definite, apply the theorem to \( \V^{-1/2}\Y \), whose forms
have matrices \( \V^{1/2}\A_i\V^{1/2} \). Their ranks are those of the \( \A_i \), so the rank count is
unchanged. What changes is the condition on the sum: \( \V^{1/2}\A\V^{1/2} \) is idempotent iff
\( \A\V\A=\A \) (for example \( \sum_i\A_i=\V^{-1} \)), and conditions (iii) and (iv) become
\( \A_i\V\A_i=\A_i \) and \( \A_i\V\A_j=\mathbf{0} \). These are the generalized least squares decompositions
of Chapter 31.
:::

Most analysis of variance tables in this book come from a chain of nested model spaces.
The successive differences of their projections are mutually orthogonal projections whose
ranks add up, so Cochran's theorem applies without further work
(@exr-qf-sequential, @thm-ss-decomposition).

## Exercises

### A. Check your understanding

::: {#exr-qf-cochran-fails}
[A1]

Let \( \A_1=\diag(2,0) \) and \( \A_2=\diag(-1,1) \), so that \( \A_1+\A_2=\I_2 \). Which hypotheses of
@thm-qf-cochran fail, and what are the laws of the two forms when \( \Y\sim\Normal_2(\bzero,\I) \)?
:::

::: {.solution}
The sum \( \I_2 \) is idempotent, but every one of the conditions
(ii)–(iv) fails. The ranks are \( 1 \) and \( 2 \), which add to \( 3\ne2 \). Neither matrix is idempotent, since
\( \A_1^2=\diag(4,0) \) and \( \A_2^2=\I_2 \). And \( \A_1\A_2=\diag(-2,0)\ne\mathbf{0} \). The forms are \( 2Y_1^2 \), which is
twice a \( \chi^2(1) \) variable, and \( Y_2^2-Y_1^2 \), the difference of two independent \( \chi^2(1) \) variables,
which takes negative values. Neither is chi-squared. They are dependent, because
\[
\Cov(2Y_1^2,\,Y_2^2-Y_1^2)=2\Cov(Y_1^2,Y_2^2)-2\Var(Y_1^2)=0-2\cdot2=-4\ne0 .
\]
:::

### B. Practice

::: {#exr-qf-twoway}
[B1]

A two-way layout has one observation \( Y_{ij} \) for each of \( a \) rows and \( b \) columns. Write
\( \norm{\Y}^2 \) as the sum of the grand-mean, row, column and residual sums of squares, identify the four
matrices, and use @thm-qf-cochran to show that under
\( Y_{ij}\sim\Normal(\mu+\alpha_i+\beta_j,\sigma^2) \) independently the four sums of squares are
independent, with degrees of freedom \( 1 \), \( a-1 \), \( b-1 \) and \( (a-1)(b-1) \).
:::

::: {#exr-qf-sequential}
[B2]

Let \( \bP_1,\dots,\bP_m \) be symmetric idempotent with
\( \C(\bP_1)\subset\C(\bP_2)\subset\dots\subset\C(\bP_m) \), and put \( \bP_0=\mathbf{0} \), \( \bP_{m+1}=\I \). Show
that the matrices \( \bP_i-\bP_{i-1} \), \( i=1,\dots,m+1 \), satisfy the conditions of
@thm-qf-cochran, so that the sequential sums of squares
\( \norm{(\bP_i-\bP_{i-1})\Y}^2 \) are independent noncentral chi-squared variables when
\( \Y\sim\Normal_n(\bmu,\sigma^2\I) \).
:::

### C. Going deeper

::: {#exr-qf-hogg-craig}
[C1]

Let \( \Y\sim\Normal_n(\bzero,\I) \), and \( Q=Q_1+Q_2 \) with \( Q=\Y\T\A\Y\sim\chi^2(r) \),
\( Q_1=\Y\T\A_1\Y\sim\chi^2(r_1) \) and \( Q_2=\Y\T\A_2\Y \) for a nonnegative definite \( \A_2 \). Show that
\( Q_1 \) and \( Q_2 \) are independent and \( Q_2\sim\chi^2(r-r_1) \).
:::

::: {.solution}
By @thm-qf-chisq-identity, \( \A \) and \( \A_1 \) are idempotent,
with ranks \( r \) and \( r_1 \). If \( \A\bv=\bzero \), then
\( 0=\bv\T\A_1\bv+\bv\T\A_2\bv \) with both terms nonnegative, so \( \norm{\A_1\bv}^2=\bv\T\A_1\bv=0 \).
Thus \( \Null(\A)\subseteq\Null(\A_1) \), hence \( \C(\A_1)\subseteq\C(\A) \) and \( \A\A_1=\A_1=\A_1\A \). Then
\( \A_2^2=(\A-\A_1)^2=\A-2\A_1+\A_1=\A_2 \), so \( \A_2 \) is idempotent, and \( \A_1\A_2=\A_1\A-\A_1=\mathbf{0} \). By
@thm-qf-orthogonal-projections, \( Q_1 \) and \( Q_2 \) are independent and
\( Q_2\sim\chi^2(\tr\A_2)=\chi^2(r-r_1) \).
:::
