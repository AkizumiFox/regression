# Linear transformations and marginals

## Affine maps

::: {#thm-mvn-linear}
[Affine transformations]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \), let \( \A \) be any \( m\times n \) matrix and
\( \bb\in\Real^m \). Then
\[
\A\Y+\bb\sim\Normal_m\bigl(\A\bmu+\bb,\ \A\bSigma\A\T\bigr).
\]
:::

::: {.proof}
For \( \mathbf{c}\in\Real^m \), \( \mathbf{c}\T(\A\Y+\bb)=(\A\T\mathbf{c})\T\Y+\mathbf{c}\T\bb \) is a linear
combination of \( \Y \) plus a constant, which is normal by @def-mvn and
@lem-mvn-univariate(c). So \( \A\Y+\bb \) is normal, and its mean and
covariance are given by @thm-rv-linear.
:::

The proof is one line because the definition was chosen to make it one line. There is
no rank condition on \( \A \). If \( \A \) has more rows than \( \rank(\bSigma) \), the result
is a singular normal vector, and that is allowed. With a density-based
definition, the same theorem needs \( \A \) to have full row rank and \( \bSigma \) to be
positive definite, and the most useful cases fall outside it. The moment
generating function gives a second proof:
\( \E e^{\mathbf{t}\T(\A\Y+\bb)}=e^{\mathbf{t}\T\bb}M_{\Y}(\A\T\mathbf{t})
=\exp\{\mathbf{t}\T(\A\bmu+\bb)+\tfrac12\mathbf{t}\T\A\bSigma\A\T\mathbf{t}\} \).

::: {#cor-mvn-marginals}
[Marginal distributions]

Every subvector of a normal vector is normal. If \( \Y\sim\Normal_n(\bmu,\bSigma) \) and
\( \Y_1 \) consists of the entries with indices in a set \( I \), then
\( \Y_1\sim\Normal(\bmu_I,\bSigma_{II}) \), where \( \bmu_I \) and \( \bSigma_{II} \) keep the
entries, rows and columns indexed by \( I \). In particular \( Y_i\sim\Normal(\mu_i,\sigma_{ii}) \).
:::

::: {.proof}
\( \Y_1=\A\Y \), where the rows of \( \A \) are the coordinate vectors \( \mathbf{e}_i\T \),
\( i\in I \). Then \( \A\bmu=\bmu_I \) and \( \A\bSigma\A\T=\bSigma_{II} \).
:::

::: {#cor-mvn-standardize}
[Standardization and rotation]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \).

::: {.enumerate options="label=(\alph*)"}
1. If \( \bSigma \) is positive definite and \( \bSigma=\A\A\T \) with \( \A \) square, then
           \( \A^{-1}(\Y-\bmu)\sim\Normal_n(\bzero,\I) \). Its entries are independent
           standard normal variables.

2. If \( \bSigma=\sigma^2\I \) and \( \Q \) is an \( n\times n \) orthogonal matrix, then
           \( \Q\T\Y\sim\Normal_n(\Q\T\bmu,\sigma^2\I) \).
:::

:::

::: {.proof}
Both follow from @thm-mvn-linear: \( \A^{-1}\A\A\T(\A^{-1})\T=\I \) and
\( \Q\T(\sigma^2\I)\Q=\sigma^2\I \). Independence in (a) was shown in the proof of @prp-mvn-mahalanobis.
:::

Part (b) is the probabilistic basis of the canonical form of the linear model.
Spherical normal errors look the same in every orthonormal coordinate system, so
one may rotate to coordinates adapted to \( \C(\X) \), in which fitted values and
residuals occupy separate independent coordinates. [Chapter 4](../ch04-quadratic-forms/index.html) uses
this rotation to find the distributions of sums of squares.

## The rank of a normal vector

@thm-mvn-existence writes a normal vector as an affine function of
independent standard normals *in distribution*. When \( \bSigma \) is singular the
representation can be made exact, with as few standard normals as the rank.

::: {#prp-mvn-rank-representation}
Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) with \( \rank(\bSigma)=r\ge1 \), and let
\( \bSigma=\Q\bLambda\Q\T \), where \( \Q \) is \( n\times r \) with orthonormal columns and
\( \bLambda=\diag(\lambda_1,\dots,\lambda_r) \) holds the positive eigenvalues. Put
\( \W=\bLambda^{-1/2}\Q\T(\Y-\bmu) \). Then \( \W\sim\Normal_r(\bzero,\I_r) \), and
\[
\Y=\bmu+\Q\bLambda^{1/2}\W\qquad\text{with probability one.}
\]
:::

::: {.proof}
By @thm-mvn-linear, \( \W \) is normal with covariance
\( \bLambda^{-1/2}\Q\T\Q\bLambda\Q\T\Q\bLambda^{-1/2}=\I_r \). Next,
\( \Q\bLambda^{1/2}\W=\Q\Q\T(\Y-\bmu) \), and \( \Q\Q\T \) is the symmetric idempotent
matrix with range \( \C(\bSigma) \) (the orthogonal projection onto it, in the
language of [Chapter 6](../ch06-projections/index.html)). Since \( \Y-\bmu\in\C(\bSigma) \) with probability
one (@thm-rv-cov-nnd), \( \Q\Q\T(\Y-\bmu)=\Y-\bmu \) with probability one.
:::

So a normal vector of rank \( r \) is exactly \( r \) independent standard normals placed
in \( \Real^n \) by an injective affine map. @exr-rv-singular-whitening gives the
same construction for any random vector, with uncorrelated entries in place of independent
normal ones. Singular normal vectors are
lower-dimensional normal vectors in disguise.

::: {#exm-mvn-ls-estimator}
[Least squares under normal errors]

Let \( \Y\sim\Normal_n(\X\bbeta,\sigma^2\I) \) with \( \X \) of full column rank \( p \), and let
\( \M=\X(\X\T\X)^{-1}\X\T \), the symmetric idempotent matrix of rank \( p \) with
\( \M\X=\X \) ([Chapter 5](../ch05-model-and-least-squares/index.html) and [Chapter 6](../ch06-projections/index.html)). By
@thm-mvn-linear:

::: {.enumerate options="label=(\alph*)"}
1. \( \hbeta=(\X\T\X)^{-1}\X\T\Y\sim\Normal_p\bigl(\bbeta,\sigma^2(\X\T\X)^{-1}\bigr) \),
           a nonsingular normal vector;

2. \( \blambda\T\hbeta\sim\Normal(\blambda\T\bbeta,\sigma^2\blambda\T(\X\T\X)^{-1}\blambda) \)
           for every \( \blambda\in\Real^p \);

3. the fitted values \( \hY=\M\Y\sim\Normal_n(\X\bbeta,\sigma^2\M) \), singular of rank \( p \);

4. the residuals \( \he=(\I-\M)\Y\sim\Normal_n(\bzero,\sigma^2(\I-\M)) \), singular of
           rank \( n-p \).
:::

No density for \( \hY \) or \( \he \) was needed, and none exists.
:::

## Independent normal vectors

::: {#lem-mvn-stack}
If \( \Y_1\sim\Normal_{n_1}(\bmu_1,\bSigma_1) \) and \( \Y_2\sim\Normal_{n_2}(\bmu_2,\bSigma_2) \)
are independent, then \( (\Y_1\T,\Y_2\T)\T \) is jointly normal with mean
\( (\bmu_1\T,\bmu_2\T)\T \) and block-diagonal covariance \( \diag(\bSigma_1,\bSigma_2) \).
Consequently, if also \( n_1=n_2 \), then \( \Y_1+\Y_2\sim\Normal(\bmu_1+\bmu_2,\bSigma_1+\bSigma_2) \).
:::

::: {.proof}
\( \mathbf{a}_1\T\Y_1+\mathbf{a}_2\T\Y_2 \) is a sum of independent normal variables, which is
normal by @lem-mvn-univariate(c). Independence makes the cross-covariance
zero. The sum is \( [\I,\I] \) applied to the stacked vector.
:::

Independence is essential in this lemma. Two normal vectors defined on the same
probability space need not be jointly normal, as we now see.

## Normal marginals are not enough {#sec-mvn-normal-marginals}

@cor-mvn-marginals cannot be reversed. A vector whose entries are each normal
need not be jointly normal. Two examples make the point, one extreme and one
realistic.

::: {#exm-mvn-sign-flip}
[A random sign]

Let \( Z\sim\Normal(0,1) \) and let \( S \) be independent of \( Z \) with
\( \Pr(S=1)=\Pr(S=-1)=\tfrac12 \). Put \( W=SZ \). By symmetry of \( Z \), \( W\sim\Normal(0,1) \).
Also \( \Cov(Z,W)=\E(S)\E(Z^2)=0 \). But \( Z+W=(1+S)Z \) equals zero with probability
\( \tfrac12 \) and is otherwise \( \Normal(0,4) \). A variable with an atom at zero and a
continuous part is not normal, so \( (Z,W) \) is not jointly normal. And \( Z \) and \( W \) are
far from independent: \( \lvert W\rvert=\lvert Z\rvert \) always.
:::

::: {#exm-mvn-mixture}
[A mixture of two correlations]

Fix \( 0<\rho<1 \). With probability \( \tfrac12 \) draw \( (Y_1,Y_2) \) from the bivariate
normal with unit variances and correlation \( \rho \), and otherwise from the one with
correlation \( -\rho \). Both components have \( \Normal(0,1) \) margins, so the mixture
does too, exactly. The covariance is \( \tfrac12\rho+\tfrac12(-\rho)=0 \). The pair is
nevertheless dependent: since \( \E(Y_1^2Y_2^2)=1+2r^2 \) for a standard bivariate normal
with correlation \( r \), the squares have correlation \( \rho^2 \), not zero. The
direction \( \mathbf{a}=(1,1)\T/\sqrt2 \) exposes the failure of joint normality.
\( S=(Y_1+Y_2)/\sqrt2 \) is an equal mixture of \( \Normal(0,1+\rho) \) and
\( \Normal(0,1-\rho) \). It has variance \( 1 \) but kurtosis
\( \E S^4=\tfrac32\{(1+\rho)^2+(1-\rho)^2\}=3(1+\rho^2) \), where a normal variable
would have \( 3 \).

[Figure 3.3.1](03-linear.html#fig-mvn-normal-marginals) shows a simulation with \( \rho=0.8 \)
and \( 200{,}000 \) draws. Kolmogorov–Smirnov tests of normality for the
two margins give \( p \)-values \( 0.84 \) and \( 0.29 \). The
sample correlation is \( -0.0033 \), the correlation of the squares is
\( 0.640 \) (theory \( 0.64 \)), and the kurtosis of
\( S \) is \( 4.88 \) (theory \( 4.92 \)). The
Kolmogorov–Smirnov distance of \( S \) from \( \Normal(0,1) \) is \( 0.067 \),
which with this many draws rejects normality at any conventional level.
:::

::: {when-format="html"}
![**Figure 3.3.1.** An equal mixture of bivariate normals with correlations \( \pm0.8 \).
(a) Both margins are exactly \( \Normal(0,1) \) and the correlation is zero, but
the scatter is an X, not an ellipse. (b) The linear combination
\( (Y_1+Y_2)/\sqrt2 \) has variance one but is sharply peaked and heavy-tailed
compared with \( \Normal(0,1) \). By @def-mvn, one nonnormal direction is
enough to rule out joint normality.](normal_marginals.svg){#fig-mvn-normal-marginals width=100%}
:::

::: {when-format="pdf"}
![An equal mixture of bivariate normals with correlations \( \pm0.8 \).
(a) Both margins are exactly \( \Normal(0,1) \) and the correlation is zero, but
the scatter is an X, not an ellipse. (b) The linear combination
\( (Y_1+Y_2)/\sqrt2 \) has variance one but is sharply peaked and heavy-tailed
compared with \( \Normal(0,1) \). By @def-mvn, one nonnormal direction is
enough to rule out joint normality.](normal_marginals.pdf){width=100%}
:::

```{.python .run #cell-normal-marginals-mixture}
import numpy as np
from scipy import integrate, stats
rng = np.random.default_rng(302)
rho, N = 0.8, 200_000

sign = rng.choice([-1.0, 1.0], size=N)          # which component each draw uses
Z1, Z2 = rng.standard_normal((2, N))
Y1 = Z1
Y2 = sign * rho * Z1 + np.sqrt(1 - rho**2) * Z2  # corr(Y1, Y2 | sign) = sign * rho

S = (Y1 + Y2) / np.sqrt(2)                     # a direction that is not normal
print("marginal KS p-values:", stats.kstest(Y1, "norm").pvalue, stats.kstest(Y2, "norm").pvalue)
print("correlation:", np.corrcoef(Y1, Y2)[0, 1])
print("corr of squares:", np.corrcoef(Y1**2, Y2**2)[0, 1])
print("kurtosis of S:", stats.kurtosis(S, fisher=False), " theory:", 3 * (1 + rho**2))
print("KS p-value of S:", stats.kstest(S, "norm").pvalue)
```

::: {.idea}
Normality of a vector is a property of *every direction*, not of the
coordinate axes. Normal margins with zero correlation can hide strong
dependence. The theorems of the next two sections, which turn zero covariance into
independence and give linear conditional means, need *joint* normality.
They are not about variables that merely happen to be normal one at a time.
:::

## Exercises

### A. Check your understanding

::: {#exr-mvn-linear-numeric}
[A1]

Let \( \Y\sim\Normal_3(\bmu,\bSigma) \) with \( \bmu=(0,1,2)\T \) and
\( \bSigma=\begin{psmallmatrix}2&1&0\\1&2&1\\0&1&2\end{psmallmatrix} \). Find the joint distribution of
\( U=Y_1+Y_2+Y_3 \), \( V=Y_1-Y_3 \) and \( W=Y_1-2Y_2+Y_3 \). Which of the pairs are independent?
Is \( V \) independent of \( (U,W) \)?
:::

::: {#exr-mvn-residual-leverage}
[A2]

In @exm-mvn-ls-estimator, show that the \( i \)th residual is
\( \Normal\bigl(0,\sigma^2(1-m_{ii})\bigr) \), where \( m_{ii} \) is the \( i \)th diagonal entry of \( \M \),
and that it is zero with probability one iff \( m_{ii}=1 \).
:::

### B. Practice

::: {#exr-mvn-rotation}
[B1]

Let \( \Y\sim\Normal_n(\mu\bone,\sigma^2\I) \) and let \( \Q \) be an orthogonal matrix whose first row
is \( n^{-1/2}\bone\T \). Show that the entries of \( \mathbf{V}=\Q\Y \) are independent, with
\( V_1\sim\Normal(\sqrt n\mu,\sigma^2) \) and \( V_i\sim\Normal(0,\sigma^2) \) for \( i\ge2 \). Show that
\( \sum_i(Y_i-\bar{Y})^2=\sum_{i\ge2}V_i^2 \), and conclude again that \( \bar{Y} \) and \( S^2 \) are
independent.
:::

::: {#exr-mvn-truncated-sign}
[B2]

Let \( Z\sim\Normal(0,1) \) and, for \( c>0 \), let \( W=Z \) if \( \lvert Z\rvert\le c \) and \( W=-Z \)
otherwise. Show that \( W\sim\Normal(0,1) \), that there is a \( c \) with \( \Cov(Z,W)=0 \), and that
\( (Z,W) \) is not jointly normal for any \( c \).
:::

::: {.solution}
Because \( Z \) and \( -Z \) have the same distribution and
the event \( \{\lvert Z\rvert>c\} \) is unchanged by the sign change,
\( \Pr(W\le w)=\Pr(Z\le w,\lvert Z\rvert\le c)+\Pr(-Z\le w,\lvert Z\rvert>c)=\Pr(Z\le w) \).
Next, \( \Cov(Z,W)=\E(Z^2;\lvert Z\rvert\le c)-\E(Z^2;\lvert Z\rvert>c)=2\E(Z^2;\lvert Z\rvert\le c)-1 \),
which is continuous and increasing in \( c \), tends to \( -1 \) as \( c\to0 \) and to \( 1 \) as \( c\to\infty \),
and so vanishes at some \( c \). Finally \( Z+W=2Z\,\mathbf{1}\{\lvert Z\rvert\le c\} \) equals \( 0 \) with
probability \( \Pr(\lvert Z\rvert>c)\in(0,1) \) but is not constant. A normal variable is either
constant or has a density, so \( Z+W \) is not normal.
:::

::: {#exr-mvn-sign-flip-directions}
[B3]

For the pair \( (Z,W) \) of @exm-mvn-sign-flip, find every \( \mathbf{a}\ne\bzero \) for which
\( a_1Z+a_2W \) is normal. *Hint:* an equal mixture of \( \Normal(0,v_1) \) and \( \Normal(0,v_2) \)
has fourth moment \( \tfrac32(v_1^2+v_2^2) \). Conclude that \( (Z,W) \) is normal along exactly two
lines through the origin, the coordinate axes.
:::

### C. Going deeper

::: {#exr-mvn-four-directions}
[C1]

@exm-mvn-mixture is normal along the two coordinate directions only. Normality along more
directions is still not enough. Let \( \varphi \) be the \( \Normal(0,1) \) density, fix
\( 0<\varepsilon\le e^2/(8\pi) \), and put
\[
f(x,y)=\varphi(x)\varphi(y)+\varepsilon\,xy(x^2-y^2)\,e^{-(x^2+y^2)} .
\]

::: {.enumerate options="label=(\alph*)"}
1. Show that \( f \) is a probability density on \( \Real^2 \). *Hint:* in polar coordinates
           \( xy(x^2-y^2)=\tfrac14r^4\sin4\theta \), and \( r^4e^{-r^2/2} \) is largest at \( r=2 \).

2. Show that both margins are \( \Normal(0,1) \) and that \( \Cov(X,Y)=0 \).

3. Show that \( (X+Y)/\sqrt2 \) and \( (X-Y)/\sqrt2 \) are \( \Normal(0,1) \) as well. *Hint:* write the
           perturbation in the coordinates \( u=(x+y)/\sqrt2 \), \( v=(x-y)/\sqrt2 \).

4. Show that \( (X,Y) \) is nevertheless not jointly normal, by computing the sign of
           \( \E(X^3Y-XY^3) \). Which four directions of @def-mvn have been checked, and which
           direction would expose the failure?
:::

:::

::: {.solution}
(a) The perturbation is continuous and integrates to \( 0 \) by (b) below, so
only \( f\ge0 \) needs checking. In polar coordinates
\( \lvert xy(x^2-y^2)\rvert e^{-r^2}=\tfrac14r^4\lvert\sin4\theta\rvert e^{-r^2}\le\tfrac14r^4e^{-r^2} \),
while \( \varphi(x)\varphi(y)=(2\pi)^{-1}e^{-r^2/2} \). So it suffices that
\( \varepsilon\,\tfrac14r^4e^{-r^2/2}\le(2\pi)^{-1} \) for all \( r\ge0 \). The function \( r^4e^{-r^2/2} \)
has derivative \( (4r^3-r^5)e^{-r^2/2} \), so its maximum is at \( r=2 \), with value \( 16e^{-2} \). The
condition is \( 4\varepsilon e^{-2}\le(2\pi)^{-1} \), that is \( \varepsilon\le e^2/(8\pi) \).

(b) For fixed \( x \),
\( \int xy(x^2-y^2)e^{-(x^2+y^2)}dy=xe^{-x^2}\int y(x^2-y^2)e^{-y^2}dy=0 \), the integrand being
odd in \( y \). So the density of \( X \) is \( \varphi \), and by symmetry so is that of \( Y \). Integrating
once more shows that \( f \) integrates to \( 1 \). For the covariance,
\( \E(XY)=\varepsilon\iint x^2y^2(x^2-y^2)e^{-(x^2+y^2)}dx\,dy \), which is zero because the
integrand changes sign under the interchange \( (x,y)\mapsto(y,x) \).

(c) With \( x=(u+v)/\sqrt2 \) and \( y=(u-v)/\sqrt2 \) one has \( xy=(u^2-v^2)/2 \),
\( x^2-y^2=2uv \) and \( x^2+y^2=u^2+v^2 \), so \( xy(x^2-y^2)=uv(u^2-v^2) \). The map is orthogonal,
so the joint density of \( (U,V)=\bigl((X+Y)/\sqrt2,(X-Y)/\sqrt2\bigr) \) is \( f(u,v) \): the *same*
function. By (b) its margins are \( \Normal(0,1) \).

(d) Under the normal part \( \E(X^3Y-XY^3)=0 \), and the perturbation contributes
\( \varepsilon\iint\{xy(x^2-y^2)\}^2e^{-(x^2+y^2)}dx\,dy>0 \). If \( (X,Y) \) were jointly normal it
would be \( \Normal_2(\bzero,\I) \) by (b), hence would have independent entries by
@thm-mvn-independence, and then \( \E(X^3Y)=\E(X^3)\E(Y)=0=\E(XY^3) \). So it is not. We have
checked the directions \( \mathbf{e}_1,\mathbf{e}_2 \) and \( (1,\pm1)\T/\sqrt2 \); a direction such as
\( (\cos\pi/8,\sin\pi/8)\T \), where \( \sin4\theta \) does not vanish, gives a projection that is
not normal.
:::

::: {#exr-mvn-residual-subvector}
[C2]

In @exm-mvn-ls-estimator let \( I\subseteq\{1,\dots,n\} \) be a set of cases, let \( \he_I \) be the
corresponding subvector of residuals and let \( \X_{I^c} \) be the rows of \( \X \) outside \( I \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \Cov(\he_I)=\sigma^2(\I-\M)_{II} \) is singular if and only if some nonzero vector
           of \( \C(\X) \) vanishes outside \( I \).

2. Deduce that \( \he_I \) has a density on \( \Real^{\lvert I\rvert} \) if and only if
           \( \rank(\X_{I^c})=\rank(\X) \): the residuals of a set of cases are nondegenerate exactly
           when deleting those cases leaves a design of unchanged rank (@thm-mvn-density).

3. Recover @exr-mvn-residual-leverage as the case \( \lvert I\rvert=1 \), and give a design with
           \( n=4 \) in which two residuals are degenerate as a pair although each is nondegenerate
           alone.
:::

:::

::: {.solution}
(a) For \( \mathbf{a}\in\Real^{\lvert I\rvert} \) let \( \tilde{\mathbf{a}} \) be \( \mathbf{a} \)
extended by zeros. Then
\( \mathbf{a}\T(\I-\M)_{II}\mathbf{a}=\tilde{\mathbf{a}}\T(\I-\M)\tilde{\mathbf{a}}=\norm{(\I-\M)\tilde{\mathbf{a}}}^2 \),
because \( \I-\M \) is symmetric and idempotent. A nonnegative definite matrix is singular iff
its quadratic form vanishes somewhere off the origin, so \( (\I-\M)_{II} \) is singular iff there
is \( \mathbf{a}\ne\bzero \) with \( (\I-\M)\tilde{\mathbf{a}}=\bzero \), that is, with
\( \tilde{\mathbf{a}}\in\C(\X) \) and \( \tilde{\mathbf{a}} \) zero outside \( I \).

(b) Such a \( \tilde{\mathbf{a}}=\X\bb \) is nonzero with \( \X_{I^c}\bb=\bzero \), so
\( \Null(\X_{I^c}) \) is strictly larger than \( \Null(\X) \); conversely any
\( \bb\in\Null(\X_{I^c})\setminus\Null(\X) \) gives such a vector. Since
\( \Null(\X)\subseteq\Null(\X_{I^c}) \) always, the two null spaces differ iff their dimensions
do, which by @thm-mat-rank-nullity happens iff \( \rank(\X_{I^c})<\rank(\X) \). A normal vector
has a density iff its covariance matrix is nonsingular (@thm-mvn-density), which gives the
statement.

(c) With \( I=\{i\} \) the condition reads \( 1-m_{ii}=0 \), that is, \( \mathbf{e}_i\in\C(\X) \): the
design has a column combination concentrated on case \( i \). For a pair, take
\( \X=[\bone,\x] \) with \( \x=(0,0,1,1)\T \) and \( I=\{1,2\} \). No nonzero vector of \( \C(\X) \) is
supported on a single case, so each residual is nondegenerate; but
\( \bone-\x=(1,1,0,0)\T\in\C(\X) \) vanishes outside \( I \), so \( (\he_1,\he_2) \) is degenerate.
Indeed \( \he_1=-\he_2 \), the two cases sharing one group mean. Deleting both cases leaves the
design \( [\bone,\x] \) with \( \x \) constant, of rank \( 1<2 \).
:::
