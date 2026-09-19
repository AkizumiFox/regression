# Estimable functions

Identifiability asks whether a quantity is determined by the distribution of the data.
Estimation asks for more. We want a statistic whose distribution is centred on the quantity,
however the parameter is set. For linear functions of \( \bbeta \), the natural candidates are linear
functions of \( \Y \), and the natural centring requirement is unbiasedness. This leads to the
classical notion of an estimable function.

## Linear unbiased estimators

Here and below the model is \( \E(\Y)=\X\bbeta \) with \( \X \) an \( n\times p \) matrix of rank \( r \);
no assumption on \( \Cov(\Y) \) is needed until variances are computed.

::: {#def-est-estimable}
[Estimable function]

A linear function \( \blambda\T\bbeta \), \( \blambda\in\Real^p \), is **estimable** if there is a vector
\( \bm a\in\Real^n \) such that
\[
\E(\bm a\T\Y)=\blambda\T\bbeta\qquad\text{for every }\bbeta\in\Real^p .
\]
A vector \( \bLambda\T\bbeta \) of linear functions, with \( \bLambda \) a \( p\times q \) matrix, is
estimable if each of its \( q \) components is. Equivalently, there is an \( n\times q \) matrix \( \A \)
with \( \E(\A\T\Y)=\bLambda\T\bbeta \) for all \( \bbeta \).
:::

The statistic \( \bm a\T\Y \) is a **linear unbiased estimator** of \( \blambda\T\bbeta \). The
requirement “for every \( \bbeta \)” is essential. For a single fixed \( \bbeta \), almost any
statistic can be made unbiased by adding a constant. What the definition demands is one
statistic that is centred correctly whatever the true coefficients are. [Chapter 7](../ch07-optimality/index.html)
asked which linear unbiased estimator is best (@def-opt-blue). Here we ask only whether there
is one.

One might allow an *affine* estimator \( a_0+\bm a\T\Y \). This gains nothing, as
[Chapter 7](../ch07-optimality/index.html) already showed.

::: {#lem-est-affine-unbiased}
The statistic \( a_0+\bm a\T\Y \) is unbiased for \( \blambda\T\bbeta \), for every \( \bbeta \), iff
\( a_0=0 \) and \( \X\T\bm a=\blambda \).
:::

::: {.proof}
This is @prp-opt-lue, proved in [Chapter 7](../ch07-optimality/index.html) by evaluating the
expectation \( a_0+\bm a\T\X\bbeta \) at \( \bbeta=\bzero \) and at the coordinate vectors.
:::

So the linear unbiased estimators of \( \blambda\T\bbeta \) correspond exactly to the solutions
\( \bm a \) of the linear system \( \X\T\bm a=\blambda \). There is at least one iff the system is
consistent, and when it is, the solutions form the affine subspace \( \bm a_0+\Null(\X\T) \)
(@thm-mat-consistency).

## Estimable is the same as identifiable

::: {#thm-est-estimable-identifiable}
For \( \blambda\in\Real^p \) the following are equivalent:

::: {.enumerate options="label=(\alph*)"}
1. \( \blambda\T\bbeta \) is estimable;

2. \( \blambda\T\bbeta \) is identifiable;

3. \( \blambda\in\C(\X\T) \), that is, \( \blambda\T=\bm\rho\T\X \) for some \( \bm\rho\in\Real^n \).
:::

For a \( p\times q \) matrix \( \bLambda \), \( \bLambda\T\bbeta \) is estimable iff \( \bLambda=\X\T\bm P \) for
some \( n\times q \) matrix \( \bm P \), that is, iff \( \C(\bLambda)\subseteq\C(\X\T) \).
:::

::: {.proof}
(a)\( \Leftrightarrow \)(c) is @lem-est-affine-unbiased: an unbiased \( \bm a\T\Y \) exists iff
\( \X\T\bm a=\blambda \) has a solution, which is the statement \( \blambda\in\C(\X\T) \).
(b)\( \Leftrightarrow \)(c) is @prp-est-linear-identifiable. The matrix version applies the vector
version to each column of \( \bLambda \) and collects the vectors \( \bm\rho \) as the columns of \( \bm P \).
:::

The theorem connects two different ideas. Identifiability is about what the model could
reveal with unlimited data, and it involves no estimator. Estimability is about a concrete
property of a concrete statistic. For linear functions they are the same thing, and both reduce
to membership of \( \blambda \) in the row space of \( \X \). The row space is spanned by the rows
\( \x_{(1)}\T,\dots,\x_{(n)}\T \), so the theorem says:

*A linear function of \( \bbeta \) is estimable iff it is a linear combination of the expected
responses \( \E(y_i)=\x_{(i)}\T\bbeta \).*

This reading is the most useful one in practice. To show that \( \blambda\T\bbeta \) is estimable,
write it as a combination of the means of observations that were actually made. To show that it
is not, find a direction \( \bv \) with \( \X\bv=\bzero \) and \( \blambda\T\bv\ne0 \).

::: {#cor-est-estimable-space}
[The estimable space]

The set of \( \blambda \) for which \( \blambda\T\bbeta \) is estimable is the subspace \( \C(\X\T) \) of
dimension \( r=\rank(\X) \). In particular:

::: {.enumerate options="label=(\alph*)"}
1. every expected response \( \x_{(i)}\T\bbeta \), and hence the whole mean vector \( \X\bbeta \), is
           estimable;

2. every entry of \( \X\T\X\bbeta \) is estimable;

3. at most \( r \) estimable functions can be linearly independent, and a set of \( r \) linearly
           independent ones determines all the others as linear combinations;

4. if \( r=p \), every linear function of \( \bbeta \) is estimable.
:::

:::

::: {.proof}
The first sentence is @thm-est-estimable-identifiable. For (a) take \( \bm\rho=\bm e_i \), the
\( i \)th coordinate vector. For (b), \( \X\T\X\bbeta=\X\T(\X\bbeta) \) and the rows of \( \X\T\X \) lie
in \( \C(\X\T) \). Part (c) is the dimension count, and (d) is the case \( \C(\X\T)=\Real^p \).
:::

## The mean at a new point

A frequent use of a fitted model is to estimate the mean response at a setting
\( \x_0 \) of the regressors, which may or may not have been observed. The quantity
\( \x_0\T\bbeta \) is a linear function of \( \bbeta \), so the theorem applies. It is
estimable iff \( \x_0\in\C(\X\T) \). In a full-rank model this is automatic, and the only
danger in predicting at a new point is extrapolation, meaning that the model may not hold
there. In a rank-deficient model there is a sharper limit. Some settings are
*algebraically* out of reach. The data contain no information about the mean there,
even if the model is exactly true.

::: {#exm-est-fixed-ratio}
[A fertilizer trial with a fixed ratio]

Ten plots receive nitrogen at \( 0,20,40,60,80 \) kg/ha (two plots each), and phosphorus is always
applied at half the nitrogen rate, as in a commercial blend. The model
\[
\E(y)=\beta_0+\beta_1N+\beta_2P
\]
is natural agronomically, but in this design \( \x_P=\tfrac12\x_N \), so \( \bv=(0,1,-2)\T \) spans
\( \Null(\X) \) and \( \rank(\X)=2 \). A linear function is estimable iff its coefficient vector
satisfies \( \lambda_1-2\lambda_2=0 \). So:

- \( \beta_0 \), the mean yield without fertilizer, is estimable. Plots with \( N=P=0 \) were observed.
- The mean yield at \( (N,P)=(40,20) \) is estimable, since \( 40-2\cdot20=0 \), and so is the
  mean at any dose of the blend, observed or not, such as \( (50,25) \).
- The mean yield at \( (40,0) \), nitrogen alone, is not estimable, since \( 40-0\neq0 \).
- Neither \( \beta_1 \) nor \( \beta_2 \) is estimable, but \( \beta_1+\tfrac12\beta_2 \), the gain per
  kilogram of nitrogen *when it comes with its phosphorus*, is.

The listing computes three least squares solutions: the minimum-norm one, the one that
omits \( P \), and the one that omits \( N \). They give intercepts \( 2.9579 \),
\( 2.9579 \) and \( 2.9579 \), all equal as they must be. The pairs of slopes
\( (0.0296, 0.0148) \), \( (0.0370, 0.0000) \) and \( (0.0000, 0.0740) \) all
disagree. At the blend dose \( (40,20) \) all three predict \( 4.439 \) t/ha. At the
nitrogen-only dose \( (40,0) \) they predict \( 4.143 \), \( 4.439 \) and
\( 2.958 \). The third value is just \( \hat{\beta}_0 \), because that solution
puts the whole fertilizer effect on phosphorus. Each of these is a legitimate least squares
answer, and none of them estimates anything.

The two kinds of prediction differ in kind, not in degree. Nothing in this experiment says
what nitrogen does without phosphorus. A larger experiment of the same design would not
help. What is needed is at least one plot off the line \( P=N/2 \) (@exr-est-new-row).
:::

```{.python .run #cell-fixed-ratio-design}
import numpy as np

N = np.array([0, 0, 20, 20, 40, 40, 60, 60, 80, 80], dtype=float)   # kg/ha nitrogen
P = N / 2                                                           # phosphorus, always half
X = np.column_stack([np.ones(N.size), N, P])
rng = np.random.default_rng(80)
y = 3.0 + 0.020 * N + 0.030 * P + rng.normal(scale=0.25, size=N.size)   # yield, t/ha
print("rank(X) =", np.linalg.matrix_rank(X), "  null vector (0, 1, -2):", X @ [0, 1, -2])
```

```{.python .run #cell-fixed-ratio-solutions}
b_min = np.linalg.pinv(X) @ y                                  # minimum-norm solution
b_noP = np.append(np.linalg.lstsq(X[:, :2], y, rcond=None)[0], 0.0)     # drop P
b_noN = np.insert(np.linalg.lstsq(X[:, [0, 2]], y, rcond=None)[0], 1, 0.0)  # drop N
new_points = {"N=40, P=20 (on the line)": [1, 40, 20],
              "N=40, P=0  (off the line)": [1, 40, 0]}
for name, x0 in new_points.items():
    print(name, [round(float(np.dot(x0, b)), 3) for b in (b_min, b_noP, b_noN)])
```

## Many unbiased estimators, one least squares estimator

An estimable function usually has many linear unbiased estimators. In
@exm-est-fixed-ratio, the mean at the blend dose \( (40,20) \) is estimated without bias by the
average of the two plots that received exactly that dose, \( \bm a_1\T\Y \) with
\( \bm a_1 \) putting weight \( \tfrac12 \) on each. It is also estimated without bias by the value of
any least squares fit at that dose. If \( \Cov(\Y)=\sigma^2\I \), the first estimator has
variance \( \sigma^2\norm{\bm a_1}^2=0.500\,\sigma^2 \), while the least squares
estimator has variance \( 0.100\,\sigma^2 \), five times smaller, because it uses
all ten plots.

The least squares estimator can be written in the same form. If \( \bm a \) is any solution of
\( \X\T\bm a=\blambda \), then @thm-proj-invariant-functions gives
\[
\blambda\T\hbeta=\bm a\T\X\hbeta=\bm a\T\M\y=(\M\bm a)\T\y ,
\]{#eq-est-ls-as-linear}

for every least squares solution \( \hbeta \). So the least squares estimator is the linear
unbiased estimator whose coefficient vector is the projection \( \M\bm a \) of any one coefficient
vector onto \( \C(\X) \). It is unbiased, because \( \X\T\M\bm a=\X\T\bm a=\blambda \).
Every solution of \( \X\T\bm a=\blambda \) has the form \( \M\bm a+\bw \) with
\( \bw\in\Null(\X\T)=\C(\X)\perpc \), and Pythagoras gives
\( \norm{\M\bm a+\bw}^2=\norm{\M\bm a}^2+\norm{\bw}^2 \). This is the
Gauss–Markov theorem in its simplest form: among linear unbiased estimators of an estimable
function, least squares has the smallest variance when \( \Cov(\Y)=\sigma^2\I \).
[Chapter 7](../ch07-optimality/index.html) states and proves it in general (@thm-opt-gauss-markov).

## What estimability does not depend on

Three remarks keep the concept in its place.

*It does not depend on the covariance.* Estimability involves only \( \E(\Y)=\X\bbeta \), so
it is the same whether the errors are uncorrelated, correlated, heteroscedastic or normal,
provided \( \Cov(\Y) \) is free of \( \bbeta \). Generalized least squares in Part VII has the same
estimable functions as ordinary least squares.

*It is a property of the design, not the response.* Whether \( \blambda\T\bbeta \) is
estimable can be decided before any response is measured, from \( \X \) alone. This makes
estimability a planning tool. If a comparison matters, the design must contain the rows that
make it estimable.

*It is a statement about linear functions.* A turning point \( -\beta_1/(2\beta_2) \) or a quotient of
two estimable quantities is not “estimable” in the technical sense, and usually has no unbiased
estimator of any kind. Yet when it is identifiable it is a perfectly sensible target, to be estimated by
other means. The line that must not be crossed is identifiability, not linearity. In this chapter
“nonestimable” is used only for linear functions, and for them it means the same as “not identifiable”.

::: {.warning}
Software fitted to a rank-deficient model reports a number for every coefficient and every
requested combination, estimable or not. For a nonestimable function the number is an artefact
of the particular solution the software chose, as the three predictions at \( (40,0) \) show. It
comes with a standard error that is equally arbitrary ([Section 8.3](03-characterizations.html)).
:::

## Exercises

### A. Check your understanding

::: {#exr-est-fixed-ratio-list}
[A1]

In @exm-est-fixed-ratio, decide which of the following are estimable: \( \beta_1-2\beta_2 \);
\( 2\beta_1+\beta_2 \); the mean yield at \( (N,P)=(100,50) \); the mean yield at \( (0,10) \); the
difference between the mean yields at \( (60,30) \) and \( (20,10) \). For the estimable ones give a
vector \( \bm a \) with \( \X\T\bm a=\blambda \).
:::

::: {#exr-est-two-points}
[A2]

In simple linear regression \( \E(y_i)=\beta_0+\beta_1x_i \), show that \( \beta_1 \) is estimable iff
the \( x_i \) are not all equal. If all \( x_i=c \), which linear functions are estimable, and what is
their least squares estimator?
:::

### B. Practice

::: {#exr-est-ls-unique-in-C}
[B1]

Let \( \blambda\T\bbeta \) be estimable and let \( \bm a^*=\M\bm a \) as in @eq-est-ls-as-linear.
Show that \( \bm a^* \) does not depend on which solution \( \bm a \) of \( \X\T\bm a=\blambda \) is used,
that \( \bm a^*=\X\G\blambda \) for any generalized inverse \( \G \) of \( \X\T\X \), and that
\( \norm{\bm a^*}^2=\blambda\T\G\blambda \).
:::

::: {.solution}
Two solutions differ by a vector of
\( \Null(\X\T)=\C(\X)\perpc \), which \( \M \) sends to \( \bzero \), so \( \M\bm a \) is the same for both. With
\( \M=\X\G\X\T \) (@thm-proj-M-formula), \( \M\bm a=\X\G\X\T\bm a=\X\G\blambda \). Finally,
\( \norm{\M\bm a}^2=\bm a\T\M\bm a=\bm a\T\X\G\X\T\bm a=\blambda\T\G\blambda \).
:::

::: {#exr-est-vector-estimable}
[B2]

Show that \( \bLambda\T\bbeta \) is estimable iff \( \bLambda\T\bm b_1=\bLambda\T\bm b_2 \) whenever
\( \X\bm b_1=\X\bm b_2 \), and that in that case \( \bLambda\T\hbeta=\bm P\T\M\y \) for any \( \bm P \) with
\( \X\T\bm P=\bLambda \). Show also that \( \bm P\T\M \) does not depend on the choice of \( \bm P \).
:::

::: {#exr-est-design-planning}
[B3]

A study will compare three treatments, \( A \), \( B \) and \( C \), in the additive model
\( \E(y)=\mu+\tau_{\text{trt}}+\delta_{\text{site}} \). Two sites are available, and each site can
host two of the three treatments. Which allocations make every difference \( \tau_s-\tau_t \)
estimable? Show that putting \( \{A,B\} \) at both sites leaves \( \tau_A-\tau_C \) nonestimable,
while \( \{A,B\} \) at one site and \( \{B,C\} \) at the other makes all three differences estimable.
:::

### C. Going deeper

::: {#exr-est-nonlinear-unbiased}
[C1]

Suppose the distribution of \( \be \) does not depend on \( \bbeta \). Show that if some statistic
\( f(\Y) \), linear or not, satisfies \( \E f(\Y)=g(\bbeta) \) for every \( \bbeta \), then \( g \) is
identifiable. Conclude that no statistic whatever is unbiased for a nonestimable linear function.
:::

::: {.solution}
The distribution of \( \Y=\X\bbeta+\be \) depends on \( \bbeta \) only through \( \X\bbeta \).
If \( \X\bbeta_1=\X\bbeta_2 \), then \( \Y \) has the same distribution under both, so
\( g(\bbeta_1)=\E_{\bbeta_1}f(\Y)=\E_{\bbeta_2}f(\Y)=g(\bbeta_2) \). A linear function with an unbiased
estimator is therefore identifiable, hence estimable by @thm-est-estimable-identifiable.
:::

