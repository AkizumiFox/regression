# The singular value decomposition and conditioning

The spectral theorem needs a symmetric matrix. A model matrix is rectangular, and the
decomposition that fits it is the singular value decomposition.

::: {#thm-mat-svd}
[Singular value decomposition]

Let \( \A \) be \( m\times n \) of rank \( r\ge1 \). There are an \( m\times r \) matrix \( \bU \) and an \( n\times r \)
matrix \( \V \), both with orthonormal columns, and numbers \( \sigma_1\ge\dots\ge\sigma_r>0 \), such
that
\[
\A=\bU\bD\V\T=\sum_{i=1}^r\sigma_i\bu_i\bv_i\T,\qquad \bD=\diag(\sigma_1,\dots,\sigma_r).
\]{#eq-mat-svd}

The **singular values** \( \sigma_i \) are the positive square roots of the nonzero
eigenvalues of \( \A\T\A \) (equivalently of \( \A\A\T \)). The columns of \( \V \) and \( \bU \) are
eigenvectors of \( \A\T\A \) and \( \A\A\T \), with \( \C(\V)=\C(\A\T) \) and \( \C(\bU)=\C(\A) \). Moreover
\( \A^+=\V\bD^{-1}\bU\T \) is the Moore–Penrose inverse of \( \A \).
:::

::: {.proof}
\( \A\T\A \) is nonnegative definite of rank \( r \) (@prp-mat-rank-product(c)). By the spectral
theorem its eigenvalues are \( \sigma_1^2\ge\dots\ge\sigma_r^2>0 \) and \( n-r \) zeros, with orthonormal
eigenvectors \( \V=[\bv_1,\dots,\bv_r] \) for the nonzero ones and \( \V_0 \) for the zeros. Define
\( \bU=\A\V\bD^{-1} \). Then
\( \bU\T\bU=\bD^{-1}\V\T\A\T\A\V\bD^{-1}=\bD^{-1}\bD^2\bD^{-1}=\I_r \). For each column \( \bv \) of \( \V_0 \),
\( \norm{\A\bv}^2=\bv\T\A\T\A\bv=0 \), so \( \A\V_0=\bzero \). Since \( \V\V\T+\V_0\V_0\T=\I \)
(@prp-mat-orthogonal(d)),
\[
\A=\A\V\V\T+\A\V_0\V_0\T=\bU\bD\V\T .
\]
Then \( \A\A\T\bU=\bU\bD\V\T\V\bD\bU\T\bU=\bU\bD^2 \), so the columns of \( \bU \) are eigenvectors of
\( \A\A\T \). \( \C(\bU)\subseteq\C(\A) \) by definition of \( \bU \), and both have dimension \( r \). The
statement for \( \V \) is the same argument applied to \( \A\T=\V\bD\bU\T \). The four Penrose
conditions for \( \V\bD^{-1}\bU\T \) follow from \( \bU\T\bU=\V\T\V=\I_r \). For instance,
\( \A\A^+=\bU\bU\T \) and \( \A^+\A=\V\V\T \) are symmetric.
:::

Completing \( \bU \) and \( \V \) to orthogonal matrices gives the *full* SVD
\( \A=\tilde{\bU}\tilde{\bD}\tilde{\V}\T \), with \( \tilde{\bD} \) of size \( m\times n \) carrying
\( \sigma_1,\dots,\sigma_r \) on its diagonal and zeros elsewhere. By convention,
\( \sigma_{r+1}=\dots=\sigma_{\min(m,n)}=0 \). Geometrically, \( \A \) maps the unit ball of
\( \Real^n \) onto a solid ellipsoid in \( \C(\A) \) with semi-axes \( \sigma_i\bu_i \). The input direction
\( \bv_i \) is stretched by the factor \( \sigma_i \) and turned into \( \bu_i \). The Moore–Penrose inverse
undoes the stretching on \( \C(\A) \) and sends \( \C(\A)\perpc \) to zero. The matrix \( \bU\bU\T \) is
the orthogonal projection onto \( \C(\A) \), whose theory is in [Chapter 6](../ch06-projections/index.html).

::: {#prp-mat-svd-norms}
[Norms and extremal properties]

Let \( \A \) be \( m\times n \) with singular values \( \sigma_1\ge\sigma_2\ge\cdots \), and let
\( \norm{\A}_2=\max_{\x\ne\bzero}\norm{\A\x}/\norm{\x} \) be the **spectral norm**.

::: {.enumerate options="label=(\alph*)"}
1. \( \norm{\A}_2=\sigma_1 \) and \( \norm{\A}_F^2=\sum_i\sigma_i^2 \).

2. If \( \A \) has full column rank \( n \), then \( \min_{\x\ne\bzero}\norm{\A\x}/\norm{\x}=\sigma_n \),
           attained at \( \bv_n \).

3. (Eckart–Young.) For \( k<r \), let \( \A_k=\sum_{i\le k}\sigma_i\bu_i\bv_i\T \). Then
           \( \norm{\A-\B}_2\ge\sigma_{k+1}=\norm{\A-\A_k}_2 \) for every \( \B \) with \( \rank(\B)\le k \).
:::

:::

::: {.proof}
(a) and (b): \( \norm{\A\x}^2/\norm{\x}^2=\x\T\A\T\A\x/\x\T\x \), and the eigenvalues of \( \A\T\A \)
are the \( \sigma_i^2 \) (padded with zeros), so @thm-mat-extremal-rayleigh applies. For the
Frobenius norm, \( \norm{\A}_F^2=\tr(\A\T\A)=\sum\sigma_i^2 \) by @prp-mat-eigen-basic(a).
(c) \( \A-\A_k=\sum_{i>k}\sigma_i\bu_i\bv_i\T \) has largest singular value \( \sigma_{k+1} \). If
\( \rank(\B)\le k \), then \( \dim\Null(\B)\ge n-k \), and since \( (n-k)+(k+1)>n \), a basis of \( \Null(\B) \) together
with \( \bv_1,\dots,\bv_{k+1} \) is dependent. So the two subspaces share a unit vector \( \x=\sum_{i\le k+1}c_i\bv_i \). Then
\( \norm{(\A-\B)\x}^2=\norm{\A\x}^2=\sum_{i\le k+1}\sigma_i^2c_i^2\ge\sigma_{k+1}^2 \).
:::

The Frobenius version of (c), \( \norm{\A-\B}_F^2\ge\sum_{i>k}\sigma_i^2 \), also holds
(Eckart and Young 1936). It is the basis of principal components.

## The condition number

::: {#def-mat-condition-number}
[Condition number]

The **condition number** of an \( m\times n \) matrix \( \A \) of full column rank is
\[
\kappa(\A)=\frac{\sigma_{\max}(\A)}{\sigma_{\min}(\A)}=\frac{\sigma_1}{\sigma_n}.
\]
For a nonsingular square matrix this equals \( \norm{\A}_2\norm{\A^{-1}}_2 \). A matrix
without full column rank has \( \kappa(\A)=\infty \) by convention.
:::

By @prp-mat-svd-norms, \( \kappa(\A) \) is the ratio of the largest to the smallest
stretching that \( \A \) applies to a unit vector. It is at least \( 1 \), is unchanged by scaling
\( \A \) or multiplying it on the left by a matrix with orthonormal columns, and equals \( 1 \) iff
\( \A\T\A=c\I \), that is, iff the columns are orthogonal with equal lengths. Two further facts
explain its role in regression. First, the singular values of \( \A\T\A \) are the \( \sigma_i^2 \),
so
\[
\kappa(\A\T\A)=\kappa(\A)^2 .
\]{#eq-mat-kappa-squared}

Second, \( \kappa \) bounds how much relative errors in the data can grow when a linear system
is solved. If \( \A \) is square and nonsingular, \( \A\x=\bb \) and \( \A(\x+\delta\x)=\bb+\delta\bb \), then
\( \delta\x=\A^{-1}\delta\bb \), and \( \norm{\bb}\le\sigma_1\norm{\x} \) gives
\[
\frac{\norm{\delta\x}}{\norm{\x}}\le\frac{\norm{\A^{-1}}_2\norm{\delta\bb}}{\norm{\bb}/\sigma_1}
=\kappa(\A)\,\frac{\norm{\delta\bb}}{\norm{\bb}} .
\]{#eq-mat-perturbation}

In floating point, relative errors of about \( 10^{-16} \) in the data can therefore become errors
of about \( \kappa\times10^{-16} \) in the solution. Solving the normal equations works with
\( \X\T\X \) and so, by @eq-mat-kappa-squared, squares the damage.
[Chapter 6](../ch06-projections/index.html) and [Chapter 10](../ch10-computation/index.html) take this up.

A small \( \sigma_{\min} \) has a direct statistical reading. By
@prp-mat-svd-norms(b), \( \norm{\X\bv_n}=\sigma_n \) with \( \norm{\bv_n}=1 \). The coefficients
\( \bv_n \) combine the columns of \( \X \) into a vector that is nearly zero, which is a
near-dependence among the regressors. This is the basis of the collinearity diagnostics of
[Chapter 26](../ch26-collinearity/index.html), such as the condition indices of @def-col-condition-index. The size of \( \kappa \) depends on the units of the columns, so
columns are usually rescaled before it is interpreted.

::: {#exm-mat-macro-svd}
[A macroeconomic design]

Take the \( 203 \) quarterly US observations of real GDP, real investment, real
government spending, real disposable income, population and unemployment, with an intercept,
as a model matrix with \( p=7 \) columns. As recorded, \( \kappa(\X)=5.10\times 10^{5} \).
Most of that is units: GDP is measured in billions and unemployment in percent. Rescaling each
column to unit length reduces it to \( 383 \), and centring the regressors before
rescaling reduces it further to \( 124 \). By @eq-mat-kappa-squared the
corresponding cross-product matrix has condition number \( 15386 \).

What remains is not an artefact of units. [Figure 1.8.1](08-svd.html#fig-mat-svd)(b) shows the right singular
vector for the smallest singular value, \( \sigma_{\min}=0.0175 \), of the centred and
scaled design. Its weight is concentrated on GDP (\( 0.78 \)) and disposable income
(\( -0.61 \)). The combination of the scaled columns with these coefficients
is a vector of length \( 0.0175 \), although every column has length one. The
two series move almost in lockstep: their correlation is \( 0.9990 \). No rescaling can separate their coefficients
well. Centring is not a neutral choice, however: it removes near-dependences that involve the
intercept, and Belsley et al. (1980) recommend scaling without centring for that
reason, so their diagnostics would use the unit-length figure \( 383 \).
[Section 26.3](../ch26-collinearity/03-condition-indices.html) returns to the choice. The listing computes the singular values, and the script also checks
@eq-mat-kappa-squared, the Penrose conditions for \( \V\bD^{-1}\bU\T \), and the
Eckart–Young bound.
:::

::: {when-format="html"}
![**Figure 1.8.1.** Singular value decomposition of a macroeconomic model matrix
(@exm-mat-macro-svd). (a) Singular values relative to the largest, for the columns as
recorded, rescaled to unit length, and centred then rescaled. The ratio of the first to the last
is the condition number. (b) The right singular vector for the smallest singular value of the
centred design. Its two large loadings reveal a near-dependence between GDP and disposable
income.](svd_condition.svg){#fig-mat-svd width=100%}
:::

::: {when-format="pdf"}
![Singular value decomposition of a macroeconomic model matrix
(@exm-mat-macro-svd). (a) Singular values relative to the largest, for the columns as
recorded, rescaled to unit length, and centred then rescaled. The ratio of the first to the last
is the condition number. (b) The right singular vector for the smallest singular value of the
centred design. Its two large loadings reveal a near-dependence between GDP and disposable
income.](svd_condition.pdf){width=100%}
:::

```{.python .run #cell-svd-condition-svd}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.macrodata.load_pandas().data
names = ["realgdp", "realinv", "realgovt", "realdpi", "pop", "unemp"]
labels = ["GDP", "investment", "government", "disp. income", "population", "unemployment"]
Z = data[names].to_numpy()
n = len(Z)

X = np.column_stack([np.ones(n), Z])                          # as recorded
X_unit = X / np.linalg.norm(X, axis=0)                        # unit-length columns
Xc = np.column_stack([np.ones(n), Z - Z.mean(axis=0)])      # regressors centred
X_std = Xc / np.linalg.norm(Xc, axis=0)                       # ... then length 1

for name, M in [("raw", X), ("unit length", X_unit), ("centred", X_std)]:
    U, s, Vt = np.linalg.svd(M, full_matrices=False)          # M = U diag(s) V'
    print(f"{name:12s} kappa = {s[0] / s[-1]:10.1f}   singular values / s_1:",
          np.array2string(s / s[0], precision=4))
```

## Exercises

### B. Practice

::: {#exr-mat-cond-two}
[B1]

Let \( \X=[\x_1,\x_2] \) with \( \norm{\x_1}=\norm{\x_2}=1 \) and \( \x_1\T\x_2=r \), \( \lvert r\rvert<1 \). Show
that \( \kappa(\X)=\sqrt{(1+\lvert r\rvert)/(1-\lvert r\rvert)} \). How large must \( \lvert r\rvert \) be for
\( \kappa(\X)>30 \)?
:::

::: {.solution}
\( \X\T\X=\begin{psmallmatrix}1&r\\r&1\end{psmallmatrix} \) has
eigenvalues \( 1\pm r \), with eigenvectors \( (1,\pm1)\T \). The singular values of \( \X \) are their square
roots, so \( \kappa(\X)^2=(1+\lvert r\rvert)/(1-\lvert r\rvert) \). For \( \kappa>30 \) we need
\( (1+\lvert r\rvert)/(1-\lvert r\rvert)>900 \), that is, \( \lvert r\rvert>899/901\approx0.998 \).
For two columns, even \( \kappa(\X)=100 \) requires \( \lvert r\rvert>9999/10001\approx0.9998 \).
:::

::: {#exr-mat-cond-scaling}
[B2]

Show that \( \kappa(\X)\ge\max_j\norm{\x_j}/\min_k\norm{\x_k} \), so that columns on very different
scales force a large condition number. Show by example that rescaling the columns can also
increase \( \kappa \).
:::
