# Conditional distributions

Throughout this section \( \Y\sim\Normal_n(\bmu,\bSigma) \) is split into
\( \Y_1\in\Real^{n_1} \) and \( \Y_2\in\Real^{n_2} \), with
\[
\bmu=\begin{pmatrix}\bmu_1\\\bmu_2\end{pmatrix},\qquad
\bSigma=\begin{pmatrix}\bSigma_{11}&\bSigma_{12}\\\bSigma_{21}&\bSigma_{22}\end{pmatrix},
\qquad \bSigma_{21}=\bSigma_{12}\T ,
\]
as in @thm-rv-partitioned. The question is what \( \Y_1 \) looks like once \( \Y_2 \)
has been observed. The answer is the probabilistic core of regression.

## Partial covariance with a generalized inverse

The natural answer involves \( \bSigma_{22}^{-1} \), but \( \bSigma_{22} \) may be singular,
for instance when \( \Y_2 \) contains a sum and its summands. A generalized inverse (@def-mat-ginverse)
removes the need for a separate treatment. The facts that make
this work were proved for arbitrary random vectors in [Section 2.5](../ch02-random-vectors/05-partitioned.html). We
collect them in the notation of this chapter.

::: {#lem-mvn-schur}
Let \( \bSigma \) be nonnegative definite and partitioned as above, and let \( \G \) be any
generalized inverse of \( \bSigma_{22} \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\bSigma_{21})\subseteq\C(\bSigma_{22}) \), so \( \bSigma_{21}=\bSigma_{22}\mathbf{H} \) for some
           matrix \( \mathbf{H} \);

2. \( \bSigma_{12}\G\bSigma_{22}=\bSigma_{12} \) and \( \bSigma_{22}\G\bSigma_{21}=\bSigma_{21} \);

3. \( \bSigma_{12}\G\bSigma_{21}=\mathbf{H}\T\bSigma_{22}\mathbf{H} \), which does not depend on the
           choice of \( \G \) and is nonnegative definite.
:::

The matrix \( \bSigma_{11\cdot2}=\bSigma_{11}-\bSigma_{12}\bSigma_{22}\ginv\bSigma_{21} \), the
**Schur complement** of \( \bSigma_{22} \) in \( \bSigma \), is therefore well defined. It is the
partial covariance matrix of \( \Y_1 \) given \( \Y_2 \) in the sense of @thm-rv-blp.
:::

::: {.proof}
By @thm-mvn-existence, \( \bSigma \) is the covariance matrix of some random vector, so
(a) and (b) are @thm-rv-partitioned(c). For (c),
\( \bSigma_{12}\G\bSigma_{21}=\mathbf{H}\T\bSigma_{22}\G\bSigma_{22}\mathbf{H}=\mathbf{H}\T\bSigma_{22}\mathbf{H} \).
:::

## The conditional distribution

::: {#thm-mvn-conditional}
[Conditional distribution]

Let \( \Y\sim\Normal_n(\bmu,\bSigma) \) be partitioned as above, let \( \G \) be a generalized
inverse of \( \bSigma_{22} \), and put \( \B=\bSigma_{12}\G \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \W=\Y_1-\bmu_1-\B(\Y_2-\bmu_2) \) is independent of \( \Y_2 \) and
           \( \W\sim\Normal_{n_1}(\bzero,\bSigma_{11\cdot2}) \);

2. the conditional distribution of \( \Y_1 \) given \( \Y_2=\y_2 \) is
     \[
\Y_1\mid\Y_2=\y_2\ \sim\ \Normal_{n_1}\Bigl(\bmu_1+\bSigma_{12}\bSigma_{22}\ginv(\y_2-\bmu_2),\
      \bSigma_{11}-\bSigma_{12}\bSigma_{22}\ginv\bSigma_{21}\Bigr),
\]{#eq-mvn-conditional}

           and for \( \y_2 \) in \( \bmu_2+\C(\bSigma_{22}) \), where \( \Y_2 \) lies with probability
           one, the conditional mean does not depend on the choice of \( \G \);

3. if \( \bSigma \) is positive definite, then \( \bSigma_{22}\ginv=\bSigma_{22}^{-1} \),
           \( \bSigma_{11\cdot2} \) is positive definite,
           \( \bSigma_{11\cdot2}^{-1} \) is the leading \( n_1\times n_1 \) block of \( \bSigma^{-1} \), and
           \( \det\bSigma=\det\bSigma_{22}\,\det\bSigma_{11\cdot2} \).
:::

:::

::: {.proof}
(a) \( \W=[\I,-\B](\Y-\bmu) \), so \( (\W\T,\Y_2\T)\T \) is an affine function of \( \Y \) and is
jointly normal (@thm-mvn-linear). The vector \( \W \) is the error of the best linear
predictor of \( \Y_1 \) from \( \Y_2 \), so by @thm-rv-blp(a,b) it has mean zero,
\( \Cov(\W,\Y_2)=\bzero \) and \( \Cov(\W)=\bSigma_{11\cdot2} \). These are moment calculations, valid
for any random vector. What normality adds is that zero covariance becomes independence
(@thm-mvn-independence), and that \( \W \), a linear function of \( \Y \), is normal.

(b) Write \( \Y_1=\bmu_1+\B(\Y_2-\bmu_2)+\W \). The first two terms are a function of
\( \Y_2 \), and \( \W \) is independent of \( \Y_2 \). For independent \( \W \) and \( \Y_2 \), the
conditional distribution of \( h(\W,\Y_2) \) given \( \Y_2=\y_2 \) is the distribution of
\( h(\W,\y_2) \) (Billingsley 1995). With
\( h(\bw,\y_2)=\bmu_1+\B(\y_2-\bmu_2)+\bw \), this is the normal distribution @eq-mvn-conditional.
If \( \y_2-\bmu_2=\bSigma_{22}\mathbf{k} \), then
\( \B(\y_2-\bmu_2)=\bSigma_{12}\G\bSigma_{22}\mathbf{k}=\bSigma_{12}\mathbf{k} \), which does not
involve \( \G \).

(c) If \( \bSigma \) is positive definite, so is its principal submatrix \( \bSigma_{22} \), and
its only generalized inverse is its inverse. By the block determinant formula
(@thm-mat-block-determinant), which needs only \( \bSigma_{22} \) nonsingular,
\( \det\bSigma=\det\bSigma_{22}\,\det\bSigma_{11\cdot2} \). The left side and \( \det\bSigma_{22} \) are
positive, so \( \bSigma_{11\cdot2} \) is nonsingular. It is nonnegative definite, being \( \Cov(\W) \)
by (a), hence positive definite. With both \( \bSigma_{22} \) and its Schur complement
nonsingular, the partitioned-inverse formula (@thm-mat-partitioned-inverse) gives
\( \bSigma_{11\cdot2}^{-1} \) as the leading block of \( \bSigma^{-1} \).

When \( \bSigma \) is positive definite, (b) can also be read off the densities. By (c),
\( \bSigma_{11\cdot2} \) is positive definite, so \( \W \) has a density. The map
\( (\y_1,\y_2)\mapsto(\bw,\y_2) \) has unit Jacobian, so
\( f_{\Y}(\y_1,\y_2)=f_{\W}\bigl(\y_1-\bmu_1-\B(\y_2-\bmu_2)\bigr)\,f_{\Y_2}(\y_2) \) by
independence, and dividing by \( f_{\Y_2}(\y_2) \) gives the conditional density.
:::

::: {.idea}
Every jointly normal vector splits as
\[
\Y_1=\underbrace{\bmu_1+\bSigma_{12}\bSigma_{22}\ginv(\Y_2-\bmu_2)}_{\text{function of }\Y_2}
\;+\;\underbrace{\W}_{\text{independent of }\Y_2}.
\]
The first term is a linear regression of \( \Y_1 \) on \( \Y_2 \). The second is a normal error
whose distribution is the same whatever value \( \Y_2 \) takes.
:::

## What the theorem says

Three features of @eq-mvn-conditional are special to the normal distribution,
and each is an assumption of the classical linear model.

::: {.enumerate options="label=(\roman*)"}
1. *The regression is linear.* \( \E(\Y_1\mid\Y_2=\y_2) \) is an affine function of
           \( \y_2 \), with coefficient matrix \( \bSigma_{12}\bSigma_{22}\ginv \).

2. *The errors are homoscedastic.* The conditional covariance
           \( \bSigma_{11\cdot2} \) is the same for every \( \y_2 \). Observing \( \Y_2 \) tells us
           where \( \Y_1 \) is centred, but not how spread out it is.

3. *Conditioning cannot increase uncertainty.*
           \( \bSigma_{11}-\bSigma_{11\cdot2}=\bSigma_{12}\bSigma_{22}\ginv\bSigma_{21} \) is
           nonnegative definite (@lem-mvn-schur(c)), so every linear combination
           \( \mathbf{a}\T\Y_1 \) has conditional variance at most its unconditional variance. The
           reduction is zero iff \( \bSigma_{12}=\bzero \).
:::

None of these holds in general. In the mixture of @exm-mvn-mixture, given
\( Y_1=y_1 \) the variable \( Y_2 \) is an equal mixture of \( \Normal(\rho y_1,1-\rho^2) \) and
\( \Normal(-\rho y_1,1-\rho^2) \). Its conditional mean is \( 0 \) for every \( y_1 \), but its
conditional variance is \( 1-\rho^2+\rho^2y_1^2 \), which grows with \( \lvert y_1\rvert \).
Linearity and constant variance are not consequences of normal margins or of zero
correlation. They are consequences of joint normality.

**One response, several predictors.**  The case used most in this book has
a scalar response \( Y \) and a vector \( \X \) of \( k \) predictors, jointly normal, with
\( \Var(Y)=\sigma_Y^2 \), \( \boldsymbol{\upsigma}_{XY}=\Cov(\X,Y) \) and \( \bSigma_{XX}=\Cov(\X) \) positive
definite. @thm-mvn-conditional says
\[
Y=\alpha+\bbeta\T\X+e,\qquad
\bbeta=\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY},\quad \alpha=\mu_Y-\bbeta\T\bmu_X,\quad
e\sim\Normal(0,\sigma^2)\text{ independent of }\X,
\]{#eq-mvn-random-x-model}

with \( \sigma^2=\sigma_Y^2-\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY}
=\sigma_Y^2(1-\rho_{Y\cdot X}^2) \), where
\( \rho_{Y\cdot X}^2=\boldsymbol{\upsigma}_{XY}\T\bSigma_{XX}^{-1}\boldsymbol{\upsigma}_{XY}/\sigma_Y^2 \) is the
population squared multiple correlation (@prp-rv-multiple-correlation). So if observations \( (Y_i,\X_i) \) are drawn
independently from a joint normal distribution, then conditionally on the predictors,
the normal linear model with independent homoscedastic errors holds exactly.
Inference that treats the predictors as fixed is valid for random predictors, as
[Chapter 14](../ch14-correlation-lack-of-fit-prediction/index.html) develops (@thm-cor-conditional).

**The bivariate case.**  With \( k=1 \), correlation \( \rho \) and standard deviations
\( \sigma_X,\sigma_Y \),
\[
\E(Y\mid X=x)=\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X),\qquad
\Var(Y\mid X=x)=\sigma_Y^2(1-\rho^2).
\]
In standard units, \( (\E(Y\mid X=x)-\mu_Y)/\sigma_Y=\rho\,(x-\mu_X)/\sigma_X \). The
predicted response is fewer standard deviations from its mean than the predictor is
from its mean, whenever \( \lvert\rho\rvert<1 \). This is **regression toward the mean**,
Galton's observation about the heights of parents and children
(Galton 1886), and the origin of the word “regression”. It
is also a warning. Regressing \( X \) on \( Y \) is a different line,
\( \E(X\mid Y=y)=\mu_X+\rho(\sigma_X/\sigma_Y)(y-\mu_Y) \), and when \( 0<\lvert\rho\rvert<1 \) neither
line is the major axis of the density contours.

::: {#exm-mvn-two-lines}
[Two regression lines]

Take \( \mu_X=1 \), \( \mu_Y=2 \), \( \sigma_X=1 \),
\( \sigma_Y=1.5 \) and \( \rho=0.7 \). The regression of \( Y \) on \( X \)
has slope \( 1.05 \) and residual standard deviation
\( 1.071 \). Drawn in the same \( (x,y) \) plane, the regression of \( X \) on \( Y \)
has slope \( 2.14 \), and the major axis of the ellipses lies between
them, with slope \( 1.76 \). [Figure 3.5.1](05-conditional.html#fig-mvn-conditional) shows the
three lines. The line \( \E(Y\mid X=x) \) passes through the points where each contour
ellipse has a vertical tangent (@exr-mvn-tangent): on a vertical line \( X=x \),
the density is highest at the conditional mean. A least squares line fitted to
\( 60 \) draws has slope \( 1.12 \) (standard error
\( 0.16 \)). It estimates the regression of \( Y \) on \( X \), not the major axis.
:::

::: {when-format="html"}
![**Figure 3.5.1.** Contours \( \Delta^2=1,4,9 \) of a bivariate normal density with \( \rho=0.7 \). The
conditional mean line \( \E(Y\mid X=x) \) passes through the points of vertical
tangency (dots). The line \( \E(X\mid Y=y) \) passes through the points of
horizontal tangency and is steeper, and the major axis lies between the two.
Least squares on a sample of \( 60 \) points estimates \( \E(Y\mid X=x) \).](conditional.svg){#fig-mvn-conditional width=66%}
:::

::: {when-format="pdf"}
![Contours \( \Delta^2=1,4,9 \) of a bivariate normal density with \( \rho=0.7 \). The
conditional mean line \( \E(Y\mid X=x) \) passes through the points of vertical
tangency (dots). The line \( \E(X\mid Y=y) \) passes through the points of
horizontal tangency and is steeper, and the major axis lies between the two.
Least squares on a sample of \( 60 \) points estimates \( \E(Y\mid X=x) \).](conditional.pdf){width=66%}
:::

::: {#exm-mvn-trivariate}
[Holding a third variable fixed]

Let \( \Y\sim\Normal_3(\bmu,\bSigma) \) with
\[
\bmu=\begin{pmatrix}1\\0\\-1\end{pmatrix},\qquad
\bSigma=\begin{pmatrix}4&1&2\\1&3&2\\2&2&2\end{pmatrix}.
\]
Conditioning \( Y_1 \) on \( Y_2 \) alone gives slope \( \sigma_{12}/\sigma_{22}=\tfrac{1}{3} \)
and conditional variance \( \tfrac{11}{3} \). Conditioning on \( (Y_2,Y_3) \), the
listing computes \( \B=\bSigma_{12}\bSigma_{22}^{-1} \) in exact arithmetic and finds
\[
\E(Y_1\mid Y_2=y_2,Y_3=y_3)=3-y_2+2y_3,
\qquad \Var(Y_1\mid Y_2,Y_3)=1.
\]
The coefficient of \( y_2 \) changes sign. Among outcomes with the same \( Y_3 \), larger
\( Y_2 \) goes with smaller \( Y_1 \). Across all outcomes, larger \( Y_2 \) goes with larger \( Y_1 \),
because both tend to be large when \( Y_3 \) is large. The conditional variance
\( 1 \) also equals \( \det\bSigma/\det\bSigma_{22}=2/2 \), as
@thm-mvn-conditional(c) predicts.
:::

```{.python .run #cell-conditional-exact}
from fractions import Fraction as F
import numpy as np
rng = np.random.default_rng(305)
mx, my, sx, sy, rho = 1.0, 2.0, 1.0, 1.5, 0.7
Sigma = np.array([[sx**2, rho * sx * sy], [rho * sx * sy, sy**2]])
slope_y_on_x = rho * sy / sx                    # E(Y | X = x)
slope_x_on_y = rho * sx / sy                    # E(X | Y = y), as dx/dy
cond_sd = sy * np.sqrt(1 - rho**2)
n = 60
XY = rng.multivariate_normal([mx, my], Sigma, size=n)
b_ls, a_ls = np.polyfit(XY[:, 0], XY[:, 1], 1)
Linv = np.linalg.inv(np.linalg.cholesky(Sigma))
theta = np.linspace(0, 2 * np.pi, 200_001)
circle = np.vstack([np.cos(theta), np.sin(theta)])
L = np.linalg.cholesky(Sigma)
big = rng.multivariate_normal([mx, my], Sigma, size=4_000_000)
x0 = 2.0
sl = big[np.abs(big[:, 0] - x0) < 0.01, 1]
resid = XY[:, 1] - a_ls - b_ls * XY[:, 0]
se = np.sqrt(resid @ resid / (n - 2) / np.sum((XY[:, 0] - XY[:, 0].mean())**2))
lam, U = np.linalg.eigh(Sigma)
slope_major = U[1, 1] / U[0, 1]
levels = [1, 4, 9]
gx, gy = np.meshgrid(np.linspace(-3, 5, 300), np.linspace(-4, 8, 300))
D = np.stack([gx - mx, gy - my], axis=-1)
D2 = np.einsum("...i,ij,...j->...", D, np.linalg.inv(Sigma), D)
xs = np.linspace(-2.2, 4.2, 2)
ys = np.linspace(-3.0, 7.0, 2)

mu = [F(1), F(0), F(-1)]
S = [[F(4), F(1), F(2)],
     [F(1), F(3), F(2)],
     [F(2), F(2), F(2)]]

def inv2(M):
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]

# Y1 given (Y2, Y3): B = sigma_12 Sigma_22^{-1}, variance sigma_11 - B sigma_21
S22inv = inv2([[S[1][1], S[1][2]], [S[2][1], S[2][2]]])
s12 = [S[0][1], S[0][2]]
B = [s12[0] * S22inv[0][j] + s12[1] * S22inv[1][j] for j in range(2)]
var_1_23 = S[0][0] - (B[0] * s12[0] + B[1] * s12[1])
const = mu[0] - B[0] * mu[1] - B[1] * mu[2]
print("E(Y1 | y2, y3) =", const, "+", B[0], "y2 +", B[1], "y3;  variance", var_1_23)

# (Y1, Y2) given Y3, and the partial correlation of Y1 and Y2 given Y3
C = [[S[i][j] - S[i][2] * S[2][j] / S[2][2] for j in range(2)] for i in range(2)]
print("Cov((Y1, Y2) | Y3) =", C)
```

## The conditional mean is the best predictor

[Section 6.11](../ch06-projections/11-population.html) defines the best linear predictor \( L(\Y_1\mid\Y_2) \) as a
projection, and shows that the conditional mean is the best predictor among *all*
square-integrable functions of \( \Y_2 \). For a scalar response, the gap between the two
is the approximation error of @prp-proj-conditional-expectation(c). Under joint
normality the gap is zero.

::: {#prp-mvn-best-predictor}
Let \( \Y \) be jointly normal and partitioned as above, and let
\( \mathbf{m}(\Y_2)=\bmu_1+\bSigma_{12}\bSigma_{22}\ginv(\Y_2-\bmu_2) \). For every function
\( \mathbf{g} \) with \( \E\norm{\mathbf{g}(\Y_2)}^2<\infty \),
\[
\E\norm{\Y_1-\mathbf{g}(\Y_2)}^2=\tr(\bSigma_{11\cdot2})+\E\norm{\mathbf{m}(\Y_2)-\mathbf{g}(\Y_2)}^2 .
\]
So \( \mathbf{m}(\Y_2) \) minimizes the mean squared prediction error over all predictors, it
equals the best linear predictor of @thm-proj-blp (entry by entry, when
\( \bSigma_{22} \) is positive definite), and the minimum is \( \tr(\bSigma_{11\cdot2}) \).
:::

::: {.proof}
\( \Y_1-\mathbf{g}(\Y_2)=\W+\bigl(\mathbf{m}(\Y_2)-\mathbf{g}(\Y_2)\bigr) \) with \( \W \) as in @thm-mvn-conditional.
The cross term
\( \E\bigl[\W\T(\mathbf{m}(\Y_2)-\mathbf{g}(\Y_2))\bigr] \) vanishes, because \( \W \) has mean zero and is
independent of \( \Y_2 \). And \( \E\norm{\W}^2=\tr\Cov(\W)=\tr(\bSigma_{11\cdot2}) \). The
minimizer \( \mathbf{m} \) is affine, so it is also the best affine predictor, which is the
best linear predictor.
:::

The proof is short because independence of \( \W \) and \( \Y_2 \) is much stronger than the
zero correlation that defines a projection. For nonnormal vectors only zero
correlation is available, and the best linear predictor is merely the best
*linear* predictor.

## Exercises

### A. Check your understanding

::: {#exr-mvn-conditional-numeric}
[A1]

For \( \bmu \) and \( \bSigma \) of @exm-mvn-trivariate, find the conditional distribution of
\( (Y_1,Y_2) \) given \( Y_3=y_3 \) and of \( Y_3 \) given \( (Y_1,Y_2)=(y_1,y_2) \). Verify that the
conditional variances are smaller than the unconditional ones.
:::

::: {#exr-mvn-precision-row}
[A2]

In the setting of @eq-mvn-random-x-model, let \( \boldsymbol{\Omega} \) be the inverse of
\( \Cov\bigl((Y,\X\T)\T\bigr) \). Show that the row of \( \boldsymbol{\Omega} \) belonging to \( Y \) is
\( (1,-\bbeta\T)/\sigma^2 \), where \( \sigma^2=\sigma_Y^2(1-\rho_{Y\cdot X}^2) \). So the regression
coefficients and the error variance can be read off one row of the precision matrix.
:::

::: {.solution}
Put \( Y \) first and apply
@exr-mvn-precision-conditional with the scalar block \( Y \) and \( \Y_2=\X \). The conditional
variance is \( \omega_{YY}^{-1}=\Var(Y\mid\X)=\sigma^2 \), and the coefficient vector of the
conditional mean is \( -\omega_{YY}^{-1}\boldsymbol{\Omega}_{YX}=\bbeta\T \). Hence \( \omega_{YY}=1/\sigma^2 \) and
\( \boldsymbol{\Omega}_{YX}=-\bbeta\T/\sigma^2 \). By @eq-mvn-random-x-model,
\( \sigma^2=\sigma_Y^2(1-\rho_{Y\cdot X}^2) \).
:::

### B. Practice

::: {#exr-mvn-conditional-to-joint}
[B1]

Suppose \( \Y_2\sim\Normal(\bmu_2,\bSigma_{22}) \) and that, given \( \Y_2=\y_2 \), \( \Y_1 \) is
\( \Normal(\mathbf{a}+\B\y_2,\bS) \) with \( \mathbf{a} \), \( \B \) and \( \bS \) not depending on \( \y_2 \). Show that
\( (\Y_1\T,\Y_2\T)\T \) is jointly normal and find its mean and covariance. Explain why
@exm-mvn-mixture does not contradict this.
:::

::: {.solution}
Condition on \( \Y_2 \):
\[
\E e^{\mathbf{t}_1\T\Y_1+\mathbf{t}_2\T\Y_2}
=\E\Bigl[e^{\mathbf{t}_2\T\Y_2}e^{\mathbf{t}_1\T(\mathbf{a}+\B\Y_2)+\frac12\mathbf{t}_1\T\bS\mathbf{t}_1}\Bigr]
=e^{\mathbf{t}_1\T\mathbf{a}+\frac12\mathbf{t}_1\T\bS\mathbf{t}_1}\,M_{\Y_2}(\mathbf{t}_2+\B\T\mathbf{t}_1).
\]
Inserting \( M_{\Y_2} \) from @thm-mvn-mgf, the exponent is linear plus
\( \tfrac12\mathbf{t}\T\mathbf{C}\mathbf{t} \) with
\( \mathbf{C}=\begin{psmallmatrix}\bS+\B\bSigma_{22}\B\T&\B\bSigma_{22}\\\bSigma_{22}\B\T&\bSigma_{22}\end{psmallmatrix} \),
so the pair is normal with mean \( \bigl((\mathbf{a}+\B\bmu_2)\T,\bmu_2\T\bigr)\T \) and covariance \( \mathbf{C} \). In
@exm-mvn-mixture, the conditional distribution of \( Y_2 \) given \( Y_1 \) is a mixture, not a
normal distribution, and its variance depends on \( y_1 \).
:::

::: {#exr-mvn-tangent}
[B2]

Show that for a bivariate normal density with \( \lvert\rho\rvert<1 \), the points of each contour
ellipse at which the tangent is vertical lie on the line \( y=\E(Y\mid X=x) \). Explain the
result in terms of the conditional density of \( Y \) given \( X=x \).
:::

::: {.solution}
Write the ellipse as \( \{\bmu+\bu:\bu\T\bSigma^{-1}\bu=c^2\} \). At
a point with vertical tangent, the normal vector \( \bSigma^{-1}\bu \) is horizontal, so
\( \bSigma^{-1}\bu=\lambda(1,0)\T \) and \( \bu=\lambda(\sigma_X^2,\rho\sigma_X\sigma_Y)\T \). The slope
\( u_2/u_1=\rho\sigma_Y/\sigma_X \) is that of the conditional mean line, and both such points
(\( \lambda>0 \) and \( \lambda<0 \)) lie on it. Along the vertical line \( X=x \), the joint density is
proportional to the conditional density of \( Y \), which is highest at \( \E(Y\mid X=x) \). The contour
that just touches the vertical line touches it at the point of highest density on that line.
:::

::: {#exr-mvn-random-walk}
[B3]

Let \( Z_1,\dots,Z_n \) be independent \( \Normal(0,1) \) and \( S_k=Z_1+\dots+Z_k \). Show that
\( (S_1,\dots,S_n) \) is jointly normal with \( \Cov(S_j,S_k)=\min(j,k) \). Find the conditional
distribution of \( S_j \) given \( S_n=s \), and comment on its variance as a function of \( j \).
:::

::: {.solution}
\( (S_1,\dots,S_n)\T=\bL\Z \) with \( \bL \) lower triangular and all
entries on and below the diagonal equal to one, so the vector is normal with covariance \( \bL\bL\T \),
whose \( (j,k) \) entry counts the common indices, \( \min(j,k) \). By @thm-mvn-conditional with
\( \Y_1=S_j \), \( \Y_2=S_n \),
\[
S_j\mid S_n=s\ \sim\ \Normal\Bigl(\frac jn\,s,\ j-\frac{j^2}{n}\Bigr)=\Normal\Bigl(\frac jn\,s,\ \frac{j(n-j)}{n}\Bigr).
\]
The conditional mean interpolates linearly between \( 0 \) and \( s \). The variance is zero at both ends,
where the value is known, and largest in the middle.
:::

::: {#exr-mvn-singular-conditioning}
[B4]

Let \( (Y,X_1,X_2) \) be normal with positive definite covariance, and let
\( \Y_2=(X_1,X_2,X_1+X_2)\T \). Show that \( \bSigma_{22} \) is singular. Using a generalized inverse,
show that the conditional distribution of \( Y \) given \( \Y_2 \) from @thm-mvn-conditional
coincides with its conditional distribution given \( (X_1,X_2) \).
:::

::: {#exr-mvn-precision-conditional}
[B5]

With \( \boldsymbol{\Omega}=\bSigma^{-1} \) partitioned like \( \bSigma \), show that
\( \E(\Y_1\mid\Y_2=\y_2)=\bmu_1-\boldsymbol{\Omega}_{11}^{-1}\boldsymbol{\Omega}_{12}(\y_2-\bmu_2) \) and
\( \Cov(\Y_1\mid\Y_2)=\boldsymbol{\Omega}_{11}^{-1} \).
:::

::: {.solution}
By @thm-mat-partitioned-inverse,
\( \boldsymbol{\Omega}_{11}=\bSigma_{11\cdot2}^{-1} \) and
\( \boldsymbol{\Omega}_{12}=-\bSigma_{11\cdot2}^{-1}\bSigma_{12}\bSigma_{22}^{-1} \). Hence
\( -\boldsymbol{\Omega}_{11}^{-1}\boldsymbol{\Omega}_{12}=\bSigma_{12}\bSigma_{22}^{-1} \) and
\( \boldsymbol{\Omega}_{11}^{-1}=\bSigma_{11\cdot2} \). Substitute into @eq-mvn-conditional.
:::

### C. Going deeper

::: {#exr-mvn-conditional-density}
[C1]

Assume \( \bSigma \) is positive definite. Prove @thm-mvn-conditional(b) by dividing the joint
density by the marginal density of \( \Y_2 \) directly, using the partitioned inverse and the
block determinant formula.
:::
