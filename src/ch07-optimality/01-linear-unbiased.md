# Linear unbiased estimation

To call an estimator the best, we have to say what it is compared with and by what
criterion. The comparison can't be made against every estimator. The constant estimator
\( \tilde{\beta}_1\equiv 0.5 \) of a slope has variance zero, and it is perfect when the true
slope happens to be \( 0.5 \) and useless otherwise. No estimator has smaller variance, yet nobody
would call it good. Some restriction on the class of competitors is unavoidable. The classical
restriction has two parts. The estimator must be *unbiased*, which rules out estimators that
buy precision by ignoring the data. And it must be *linear* in the response, which makes its
variance computable from the first two moments of the errors alone. This section describes the
class these two conditions carve out. The description is geometric, and the Gauss–Markov
theorem of [Section 7.2](02-gauss-markov.html) will follow from it in a few lines.

## The setting

Throughout this chapter the response satisfies the linear model of
[Chapter 5](../ch05-model-and-least-squares/index.html) (@def-lm-linear-model),
\[
\Y=\X\bbeta+\be,\qquad \E(\be)=\bzero,\qquad \Cov(\be)=\sigma^2\I,
\]{#eq-opt-model}

with \( \X \) a fixed \( n\times p \) matrix of rank \( r\le p \), \( \bbeta\in\Real^p \) unknown and
\( \sigma^2>0 \) unknown. These are the **second-moment assumptions**. The errors have mean
zero, a common variance, and are uncorrelated. Nothing is assumed about their distribution
beyond that. From [Section 7.3](03-maximum-likelihood.html) on we add normality, and we will
say so each time. As in [Chapter 6](../ch06-projections/index.html), \( \M \) denotes the orthogonal
projection onto \( \C(\X) \), \( \hbeta \) is any least squares estimate, and
\( \text{SSE}=\norm{(\I-\M)\Y}^2 \) is the residual sum of squares.

The quantities we want to estimate are linear functions \( \blambda\T\bbeta \), for fixed
\( \blambda\in\Real^p \). This covers a single coefficient (\( \blambda \) a coordinate vector), a
difference of two coefficients, and the mean response \( \x_0\T\bbeta \) at a new point \( \x_0 \).
Several such functions at once are collected as \( \bLambda\T\bbeta \), where \( \bLambda \) is a
\( p\times q \) matrix.

## Linear and unbiased

::: {#def-opt-blue}
[Linear unbiased and best linear unbiased estimators]

Let \( \blambda\in\Real^p \).

::: {.enumerate options="label=(\alph*)"}
1. A **linear estimator** is a statistic of the form \( c+\mathbf{a}\T\Y \), with a constant
   \( c\in\Real \) and a constant vector \( \mathbf{a}\in\Real^n \) (constant means not depending on
   \( \Y \) or on the parameters).

2. A linear estimator is a **linear unbiased estimator** (LUE) of \( \blambda\T\bbeta \) if
   \( \E(c+\mathbf{a}\T\Y)=\blambda\T\bbeta \) for every \( \bbeta\in\Real^p \) and every \( \sigma^2>0 \).

3. An LUE \( c+\mathbf{a}\T\Y \) of \( \blambda\T\bbeta \) is a **best linear unbiased estimator**
   (BLUE) if \( \Var(c+\mathbf{a}\T\Y)\le\Var(d+\bb\T\Y) \) for every LUE \( d+\bb\T\Y \) of
   \( \blambda\T\bbeta \), for every \( \bbeta \) and every \( \sigma^2>0 \).

4. For \( q \) functions \( \bLambda\T\bbeta \), a **linear estimator** is \( \mathbf{c}+\A\Y \) with a
   constant \( q\times n \) matrix \( \A \) and \( \mathbf{c}\in\Real^q \). It is an LUE of \( \bLambda\T\bbeta \)
   if each coordinate is an LUE of the corresponding coordinate of \( \bLambda\T\bbeta \), and a BLUE
   if \( \Cov(\B\Y+\mathbf{d})-\Cov(\A\Y+\mathbf{c}) \) is nonnegative definite for every LUE
   \( \B\Y+\mathbf{d} \) (@def-mat-nnd).
:::
:::

The words "for every \( \bbeta \)" in (b) carry all the weight. An estimator that is unbiased at
one parameter value only, like the constant \( 0.5 \) when \( \beta_1=0.5 \), does not qualify. The
first result (@prp-opt-lue) turns unbiasedness into linear equations.

::: {#prp-opt-lue}
[Unbiasedness as a linear equation]

The linear estimator \( c+\mathbf{a}\T\Y \) is an LUE of \( \blambda\T\bbeta \) iff \( c=0 \) and
\( \X\T\mathbf{a}=\blambda \). Consequently \( \blambda\T\bbeta \) has an LUE iff
\( \blambda\in\C(\X\T) \).
:::

::: {.proof}
By @thm-rv-linear, \( \E(c+\mathbf{a}\T\Y)=c+\mathbf{a}\T\X\bbeta \). This equals \( \blambda\T\bbeta \)
for all \( \bbeta \) iff the affine function \( \bbeta\mapsto c+(\X\T\mathbf{a}-\blambda)\T\bbeta \) vanishes
identically. Putting \( \bbeta=\bzero \) gives \( c=0 \). Putting \( \bbeta \) equal to each coordinate vector
in turn then gives \( \X\T\mathbf{a}=\blambda \). The converse is immediate. An LUE therefore exists iff
the equation \( \X\T\mathbf{a}=\blambda \) has a solution, that is, iff \( \blambda\in\C(\X\T) \).
:::

So the constant \( c \) is always zero, and from now on we drop it. A linear function
\( \blambda\T\bbeta \) with \( \blambda\in\C(\X\T) \) is called **estimable**. For these functions the
data carry linear information. For the others no linear combination of the observations is centred
on the target. [Chapter 8](../ch08-estimability/index.html) takes this as the definition of
estimability (@def-est-estimable) and develops its many equivalent forms. For now we need only
@prp-opt-lue and one comparison. @thm-proj-invariant-functions of [Chapter 6](../ch06-projections/index.html) found the same
condition \( \blambda\in\C(\X\T) \) as the answer to a different question: when is
\( \blambda\T\hbeta \) the same for every least squares estimate \( \hbeta \)? The functions that least squares
estimates unambiguously are exactly the functions that have linear unbiased estimators at all.
When \( \X \) has full column rank, \( \C(\X\T)=\Real^p \), so every linear function of \( \bbeta \) is
estimable.

## The geometry of the unbiased class

Fix an estimable \( \blambda\T\bbeta \) and write \( \blambda=\X\T\boldsymbol{\uprho} \) for some
\( \boldsymbol{\uprho}\in\Real^n \). The coefficient vectors of its LUEs form the set
\[
\mathcal A_{\blambda}=\{\mathbf{a}\in\Real^n:\X\T\mathbf{a}=\blambda\}.
\]
This set is the key object of the chapter's first half.

::: {#prp-opt-lue-set}
[The set of unbiased coefficient vectors]

Let \( \blambda=\X\T\boldsymbol{\uprho} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \mathcal A_{\blambda}=\M\boldsymbol{\uprho}+\C(\X)\perpc \), an affine subspace of \( \Real^n \) of
   dimension \( n-r \), parallel to \( \C(\X)\perpc \);

2. every \( \mathbf{a}\in\mathcal A_{\blambda} \) has the same projection \( \M\mathbf{a}=\M\boldsymbol{\uprho} \), and
   \( \M\boldsymbol{\uprho} \) is the only element of \( \mathcal A_{\blambda} \) that lies in \( \C(\X) \);

3. \( \M\boldsymbol{\uprho}=\X(\X\T\X)\ginv\blambda \) for every generalized inverse, and
   \( (\M\boldsymbol{\uprho})\T\Y=\blambda\T\hbeta \) for every least squares estimate \( \hbeta \).
:::
:::

::: {.proof}
(a) By @prp-opt-lue, \( \mathbf{a}\in\mathcal A_{\blambda} \) means \( \X\T\mathbf{a}=\blambda \). Now \( \X\T\mathbf{a}=\X\T\boldsymbol{\uprho} \) iff \( \X\T(\mathbf{a}-\boldsymbol{\uprho})=\bzero \), iff
\( \mathbf{a}-\boldsymbol{\uprho}\in\Null(\X\T)=\C(\X)\perpc \) (@lem-proj-null-colspace). So
\( \mathcal A_{\blambda}=\boldsymbol{\uprho}+\C(\X)\perpc \). Since \( \boldsymbol{\uprho}-\M\boldsymbol{\uprho}=(\I-\M)\boldsymbol{\uprho}\in\C(\X)\perpc \),
the same set is \( \M\boldsymbol{\uprho}+\C(\X)\perpc \). Its dimension is
\( \dim\C(\X)\perpc=n-r \).

(b) If \( \mathbf{a}=\M\boldsymbol{\uprho}+\bv \) with \( \bv\in\C(\X)\perpc \), then \( \M\mathbf{a}=\M\boldsymbol{\uprho}+\M\bv=\M\boldsymbol{\uprho} \).
If moreover \( \mathbf{a}\in\C(\X) \), then \( \mathbf{a}=\M\mathbf{a}=\M\boldsymbol{\uprho} \).

(c) By @thm-proj-M-formula, \( \M\boldsymbol{\uprho}=\X(\X\T\X)\ginv\X\T\boldsymbol{\uprho}=\X(\X\T\X)\ginv\blambda \).
By @thm-proj-invariant-functions, \( \blambda\T\hbeta=\boldsymbol{\uprho}\T\M\Y=(\M\boldsymbol{\uprho})\T\Y \).
:::

Part (c) says that least squares is itself one of the linear unbiased estimators. Its
coefficient vector is \( \mathbf{a}_*=\M\boldsymbol{\uprho} \), the one member of the class that lies in the model space.
It does not depend on which \( \boldsymbol{\uprho} \) was used to represent \( \blambda \), because every such
\( \boldsymbol{\uprho} \) has the same projection. That \( \mathbf{a}_*\T\Y \) is unbiased can also be seen directly:
\( \E(\mathbf{a}_*\T\Y)=\boldsymbol{\uprho}\T\M\X\bbeta=\boldsymbol{\uprho}\T\X\bbeta=\blambda\T\bbeta \).

The variance of any linear estimator follows from @thm-rv-linear:
\[
\Var(\mathbf{a}\T\Y)=\mathbf{a}\T(\sigma^2\I)\mathbf{a}=\sigma^2\norm{\mathbf{a}}^2 .
\]{#eq-opt-var-norm}

Under the second-moment assumptions, then, the variance of a linear estimator is \( \sigma^2 \)
times the squared length of its coefficient vector. Finding the BLUE means finding the shortest
vector in the affine set \( \mathcal A_{\blambda} \).

::: {.idea}
The linear unbiased estimators of \( \blambda\T\bbeta \) correspond to the points of a flat
\( \mathcal A_{\blambda}=\mathbf{a}_*+\C(\X)\perpc \) in \( \Real^n \). The variance of each is
\( \sigma^2 \) times its squared distance from the origin. The question "which LUE is best?"
is a nearest-point problem, the same kind of problem that least squares solves in
[Chapter 6](../ch06-projections/index.html), now posed in the space of coefficient vectors.
:::

A second description of the class is useful too. Every \( \mathbf{a}\in\mathcal A_{\blambda} \) is
\( \mathbf{a}_*+(\I-\M)\bw \) for some \( \bw \), so every LUE has the form
\[
\mathbf{a}\T\Y=\blambda\T\hbeta+\bw\T(\I-\M)\Y=\blambda\T\hbeta+\bw\T\he .
\]{#eq-opt-lue-residual}

An LUE is the least squares estimate plus a linear function of the residuals. The linear functions
of the residuals are exactly the linear unbiased estimators of zero (@exr-opt-unbiased-zero).
Adding one to an estimator keeps it unbiased and can only add noise. [Section 7.2](02-gauss-markov.html)
makes "can only add noise" precise.

## Three estimators of a slope

::: {#exm-opt-three-slopes}
[Three linear unbiased estimators of a slope]

A straight line \( \E(Y_i)=\beta_0+\beta_1x_i \) is fitted to \( n=12 \) observations at the
regressor values \( x=1,2,2.5,3,4,5,5.5,7,8,9,10.5,12 \), for which
\( S_{xx}=\sum_i(x_i-\bar{x})^2=137.229 \). Here \( \X=[\bone,\x] \) has full rank and the target is
\( \blambda=(0,1)\T \). Three natural estimators of the slope are linear in \( \Y \).

::: {.enumerate options="label=(\roman*)"}
1. *Least squares*, \( \hat{\beta}_1=\sum_i(x_i-\bar{x})Y_i/S_{xx} \) (@thm-lm-simple-ls), with
   \( a_i=(x_i-\bar{x})/S_{xx} \).

2. *End points*, the slope of the line through the first and last observations,
   \( (Y_{12}-Y_1)/(x_{12}-x_1) \).

3. *Group means*, which splits the data into the six smallest and six largest \( x \) values and
   takes \( (\bar{Y}_{\text{high}}-\bar{Y}_{\text{low}})/(\bar{x}_{\text{high}}-\bar{x}_{\text{low}}) \).
:::

Each coefficient vector satisfies \( \bone\T\mathbf{a}=0 \) and \( \x\T\mathbf{a}=1 \), which is the condition
\( \X\T\mathbf{a}=\blambda \) of @prp-opt-lue. So all three are unbiased, whatever the distribution of
the errors. By @eq-opt-var-norm their variances are \( \sigma^2\norm{\mathbf{a}}^2 \):
\[
0.00729\,\sigma^2,\qquad 0.01653\,\sigma^2,\qquad 0.01008\,\sigma^2 .
\]
Least squares has the smallest, and the other two have efficiencies \( 0.441 \) and
\( 0.723 \) relative to it. The listing also confirms what @prp-opt-lue-set(b) predicts: the
three coefficient vectors differ, but their projections onto \( \C(\X) \) all equal the least
squares vector.
:::

```{.python .run #cell-gauss-markov-design}
import numpy as np

x = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 7.0, 8.0, 9.0, 10.5, 12.0])
n = len(x)
X = np.column_stack([np.ones(n), x])
lam = np.array([0.0, 1.0])                      # target: the slope beta_1
M = X @ np.linalg.solve(X.T @ X, X.T)           # projection onto C(X)

# coefficient vectors a of three linear estimators a^T y of the slope
a_ls = X @ np.linalg.solve(X.T @ X, lam)        # least squares: a = X (X^T X)^{-1} lambda
a_end = np.zeros(n)
a_end[[0, -1]] = [-1.0, 1.0]
a_end /= x[-1] - x[0]                           # slope through the two end points
low, high = np.arange(n) < n // 2, np.arange(n) >= n // 2
a_grp = (high / high.sum() - low / low.sum()) / (x[high].mean() - x[low].mean())  # group means

for name, a in [("least squares", a_ls), ("end points", a_end), ("group means", a_grp)]:
    print(f"{name:14s} X^T a = {X.T @ a}   Var/sigma^2 = {a @ a:.5f}   Ma = a_ls: {np.allclose(M @ a, a_ls)}")
```

The end-point estimator throws away ten of the twelve observations, so its poor showing is no
surprise. The group-means estimator uses every observation, and it is still clearly worse than least
squares. It weights all observations in a half equally, while least squares weights each
observation in proportion to its distance from \( \bar{x} \). Estimators of this grouping type were
proposed by Wald (1940) for regressions in which \( x \) itself is measured with error.
Chapter 24 returns to that problem.

## Why these restrictions?

Linearity and unbiasedness are restrictions of convenience as well as of principle, and it is
worth being clear about what each one does.

**Unbiasedness** removes estimators that are good only near one parameter value. It is not
sacred. [Section 7.2](02-gauss-markov.html) shows biased estimators with smaller mean squared
error, and [Section 7.6](06-bayes-conjugate.html) derives a whole family of them from prior
information. Unbiasedness can also be too demanding: a nonestimable \( \blambda\T\bbeta \) has no
linear unbiased estimator, and under normality no unbiased estimator of any kind
(@exr-opt-no-unbiased).

**Linearity** is what makes the second-moment assumptions enough. The variance of a nonlinear
estimator, such as a median or a trimmed mean, depends on the whole error distribution. Its
variance can't be compared with that of least squares until that distribution is specified.
For linear estimators @eq-opt-var-norm needs only \( \sigma^2 \). The restriction therefore has a
cost: optimality among linear estimators says nothing about nonlinear ones. [Section 7.2](02-gauss-markov.html)
gives an error distribution under which a simple nonlinear estimator is far better than least
squares. [Section 7.4](04-minimum-variance.html) shows that under normality this cannot happen.

## Exercises

### A. Check your understanding

::: {#exr-opt-location}
[A1]

In the location model \( Y_i=\mu+\varepsilon_i \) (\( \X=\bone \)), show that the LUEs of \( \mu \) are
the weighted averages \( \sum_ia_iY_i \) with \( \sum_ia_i=1 \). Identify \( \mathcal A_{\blambda} \) and
\( \mathbf{a}_* \) and find the LUE of smallest variance directly.
:::

::: {#exr-opt-oneway-lue}
[A2]

In the one-way layout \( \E(Y_{ij})=\mu+\alpha_i \) with three groups and an intercept
(@exm-proj-ginverse-numeric), show that \( \alpha_1 \) has no linear unbiased estimator. Write down two
different LUEs of \( \alpha_1-\alpha_2 \), and verify that their coefficient vectors have the same
projection onto \( \C(\X) \).
:::

### B. Practice

::: {#exr-opt-unbiased-zero}
[B1]

A linear estimator \( \mathbf{a}\T\Y \) is an **unbiased estimator of zero** if \( \E(\mathbf{a}\T\Y)=0 \) for every
\( \bbeta \). Show that this holds iff \( \mathbf{a}\in\C(\X)\perpc \), iff \( \mathbf{a}\T\Y \) is a linear function of
the residual vector \( \he \). Show that every unbiased estimator of zero is uncorrelated with
\( \blambda\T\hbeta \) for every estimable \( \blambda\T\bbeta \).
:::

::: {.solution}
By @prp-opt-lue with \( \blambda=\bzero \), \( \mathbf{a}\T\Y \) is unbiased for zero iff
\( \X\T\mathbf{a}=\bzero \), iff \( \mathbf{a}\in\C(\X)\perpc \). If so, \( \mathbf{a}=(\I-\M)\mathbf{a} \) and
\( \mathbf{a}\T\Y=\mathbf{a}\T(\I-\M)\Y=\mathbf{a}\T\he \). Conversely, \( \bw\T\he=((\I-\M)\bw)\T\Y \) has coefficient
vector in \( \C(\X)\perpc \). For the last claim, write \( \blambda\T\hbeta=\boldsymbol{\uprho}\T\M\Y \). By
@thm-rv-linear,
\( \Cov(\boldsymbol{\uprho}\T\M\Y,\mathbf{a}\T\Y)=\sigma^2\boldsymbol{\uprho}\T\M\mathbf{a}=0 \), because \( \M\mathbf{a}=\bzero \).
:::

::: {#exr-opt-lue-matrix}
[B2]

Let \( \bLambda \) be \( p\times q \) with \( \C(\bLambda)\subseteq\C(\X\T) \). Show that \( \A\Y \) is an LUE of
\( \bLambda\T\bbeta \) iff \( \A\X=\bLambda\T \), and that the set of such \( \A \) is
\( \{\bLambda\T(\X\T\X)\ginv\X\T+\mathbf{W}(\I-\M):\mathbf{W}\in\Real^{q\times n}\} \).
:::

### C. Going deeper

::: {#exr-opt-no-unbiased}
[C1]

Assume \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) and \( \blambda\notin\C(\X\T) \). Show that no statistic
\( \delta(\Y) \), linear or not, satisfies \( \E\,\delta(\Y)=\blambda\T\bbeta \) for all \( \bbeta \). *Hint:*
the distribution of \( \Y \) depends on \( \bbeta \) only through \( \X\bbeta \).
:::

::: {.solution}
Since \( \blambda\notin\C(\X\T)=\Null(\X)\perpc \), there is \( \bv\in\Null(\X) \) with
\( \blambda\T\bv\ne0 \). The parameters \( \bbeta \) and \( \bbeta+\bv \) give the same mean
\( \X\bbeta \), hence the same distribution of \( \Y \), hence the same value of
\( \E\,\delta(\Y) \). An unbiased \( \delta \) would need
\( \blambda\T\bbeta=\blambda\T(\bbeta+\bv) \), which is false. The argument uses only that the law of
\( \Y \) is a function of \( \X\bbeta \) and \( \sigma^2 \). [Chapter 8](../ch08-estimability/index.html) calls this *nonidentifiability* (@def-est-identifiable).
:::

::: {#exr-opt-grouping}
[C2]

For the group-means estimator of @exm-opt-three-slopes with general \( x_1\le\dot{s}\le x_n \)
(\( n=2m \)), find its variance in terms of \( \bar{x}_{\text{high}}-\bar{x}_{\text{low}} \). Show that its
efficiency relative to least squares is \( (\bar{x}_{\text{high}}-\bar{x}_{\text{low}})^2\,n/(4S_{xx}) \), and
that this is \( 1 \) iff the \( x \) values take only two distinct values, with half of the observations at each.
:::
