# Computation

The formulas \( \M=\X(\X\T\X)\ginv\X\T \) and \( \hbeta=(\X\T\X)^{-1}\X\T\y \) are the right way
to *reason* about least squares and the wrong way to *compute* it. There are
two problems. \( \M \) is \( n\times n \), so forming it costs \( O(n^2) \) memory for a quantity
that is never needed in full. And forming \( \X\T\X \) squares the condition number of the
problem, which can throw away half of the available floating-point digits before any
equation is solved. This section shows how the geometry itself suggests the fix, which
is to build an orthonormal basis for \( \C(\X) \). [Chapter 10](../ch10-computation/index.html) develops the
numerical linear algebra in full.

## QR: Gram–Schmidt as a factorization

Apply Gram–Schmidt (@prp-proj-gram-schmidt) to the columns \( \x_1,\dots,\x_p \) of a
full-column-rank \( \X \). Each \( \x_j \) is a combination of \( \bu_1,\dots,\bu_j \) with
coefficients \( r_{ij}=\bu_i\T\x_j \) for \( i<j \) and \( r_{jj}=\norm{\bw_j}>0 \). In matrix form,
\[
\X=\Q\R,
\]
where \( \Q=[\bu_1,\dots,\bu_p] \) is \( n\times p \) with \( \Q\T\Q=\I_p \), and \( \R \) is \( p\times p \)
upper triangular with positive diagonal. This is the **thin QR factorization**. Since
\( \C(\Q)=\C(\X) \), @prp-proj-orthonormal-formula gives \( \M=\Q\Q\T \) at once, and the
rest of least squares follows:
\[
\hY=\Q(\Q\T\y),\qquad
\R\hbeta=\Q\T\y,\qquad
\text{SSE}=\norm{\y}^2-\norm{\Q\T\y}^2,\qquad
h_{ii}=\textstyle\sum_k q_{ik}^2 .
\]{#eq-proj-qr-ls}

The coefficient equation holds because the normal equations become
\( \R\T\R\bb=\R\T\Q\T\y \), and \( \R\T \) is nonsingular. It is solved by back substitution,
with no matrix inverse. Nothing of size \( n\times n \) is ever formed.

In floating point, Gram–Schmidt as written in @prp-proj-gram-schmidt gradually
loses orthogonality among the computed \( \bu_j \). Production software computes the same
factorization with Householder reflections, which are backward stable
(Golub and Van Loan 2013; Higham 2002). The geometry does not change: \( \Q \) is
an orthonormal basis for the model space.

## Conditioning

The sensitivity of least squares to rounding is governed by the **condition number**
\[
\kappa(\X)=\frac{\sigma_{\max}(\X)}{\sigma_{\min}(\X)},
\]
the ratio of the largest and smallest singular values of \( \X \) ([Chapter 1](../ch01-matrix-algebra/index.html)).
Geometrically, \( \kappa(\X) \) is large when some combination of columns of unit
coefficient norm is nearly zero. The columns are then close to linearly dependent, and
the coordinates of a point in \( \C(\X) \) are poorly determined, even though the point itself
is not. The singular values of \( \X\T\X \) are the squares of those of \( \X \), so
\[
\kappa(\X\T\X)=\kappa(\X)^2 .
\]

A useful rule of thumb from numerical analysis says that solving a linear system with
condition number \( \kappa \) in double precision loses about \( \log_{10}\kappa \) of the roughly
16 available significant digits. The normal equations solve a system whose matrix is
\( \X\T\X \), so they lose about \( 2\log_{10}\kappa(\X) \) digits. QR solves the triangular
system \( \R\bb=\Q\T\y \), and \( \kappa(\R)=\kappa(\X) \) because \( \Q \) has orthonormal
columns, so it loses about \( \log_{10}\kappa(\X) \). The precise perturbation theory has one
more subtlety. When the residual is large, the least squares *problem* itself has a
sensitivity term proportional to \( \kappa(\X)^2\norm{\he} \), and no algorithm can avoid it
(Higham 2002, ch. 20). QR's advantage is largest when the model fits
well.

::: {#exm-proj-qr-polynomial}
[Polynomial regression]

The columns \( 1,t,t^2,\dots,t^d \) evaluated at 200 equally spaced points in \( [0,1] \) become
nearly dependent as \( d \) grows. To isolate the numerical error, the listing sets
\( \y=\X\bbeta \) exactly for a random \( \bbeta \), so that the true least squares solution is
\( \bbeta \) itself. For \( d=8 \) the condition number is \( 6.8\times 10^{5} \). The normal
equations return \( \bbeta \) with relative error \( 1.0\times 10^{-5} \), while QR achieves
\( 1.7\times 10^{-11} \). At \( d=13 \), with
\( \kappa(\X)=4.1\times 10^{9} \), the normal-equation answer has relative error
\( 9.3\times 10^{-1} \) and not one correct digit, while QR still gets
\( 4.4\times 10^{-9} \). [Figure 6.10.1](10-computation.html#fig-proj-qr) plots the whole sequence against the
reference lines \( \kappa\epsilon \) and \( \kappa^2\epsilon \), where \( \epsilon\approx2.2\times10^{-16} \)
is machine precision.
:::

::: {when-format="html"}
![**Figure 6.10.1.** Accuracy of least squares coefficients for polynomial regressions of degree
\( 1 \) to \( 13 \), against the condition number of the model matrix. The error of the normal
equations tracks \( \kappa^2\epsilon \) until it reaches order one. The error of QR tracks
\( \kappa\epsilon \).](qr_vs_normal.svg){#fig-proj-qr width=62%}
:::

::: {when-format="pdf"}
![Accuracy of least squares coefficients for polynomial regressions of degree
\( 1 \) to \( 13 \), against the condition number of the model matrix. The error of the normal
equations tracks \( \kappa^2\epsilon \) until it reaches order one. The error of QR tracks
\( \kappa\epsilon \).](qr_vs_normal.pdf){width=62%}
:::

```{.python .run #cell-qr-vs-normal-solvers}
import numpy as np
import statsmodels.api as sm

def fit_normal_equations(X, y):
    return np.linalg.solve(X.T @ X, X.T @ y)

def fit_qr(X, y):
    Q, R = np.linalg.qr(X)                  # X = QR, Q has orthonormal columns
    return np.linalg.solve(R, Q.T @ y)      # triangular system R b = Q'y
```

```{.python .run #cell-qr-vs-normal-experiment}
rng = np.random.default_rng(10)
n = 200
t = np.linspace(0, 1, n)
rows = []

for degree in range(1, 14):
    X = np.vander(t, degree + 1, increasing=True)          # 1, t, t^2, ..., t^degree
    beta = rng.normal(size=degree + 1)
    y = X @ beta                                           # exact fit: the true answer is beta
    kappa = np.linalg.cond(X)
    err_ne = np.linalg.norm(fit_normal_equations(X, y) - beta) / np.linalg.norm(beta)
    err_qr = np.linalg.norm(fit_qr(X, y) - beta) / np.linalg.norm(beta)
    rows.append((degree, kappa, err_ne, err_qr))
    print(f"degree {degree:2d}  cond(X) = {kappa:9.2e}   normal eq. {err_ne:8.1e}   QR {err_qr:8.1e}")
```

The polynomial example is extreme by design, but ill-conditioning is common in real
data. The macroeconomic series of Longley's classic test problem (in
`statsmodels.datasets.longley`) give a model matrix with
\( \kappa(\X)\approx4.86\times 10^{9} \). Then
\( \kappa(\X\T\X)\approx2.36\times 10^{19} \), which is larger than \( 1/\epsilon \).
The rule of thumb predicts about \( 6 \) correct digits from
QR and about \( 0 \) from the normal equations. Part of that
condition number is only a matter of units: one column is a year, another is
population in thousands. Rescaling columns (@thm-proj-reparam) changes \( \kappa \)
without changing the fit (@exr-proj-longley-scaling). The part that remains
after rescaling reflects genuine near-dependence among the regressors, the subject of
Chapter 26.

Two practical lessons follow. First, never compute \( (\X\T\X)^{-1} \) just to obtain
\( \hbeta \). Solve with QR, or call a least squares routine that does. Second,
reparameterizations that leave \( \C(\X) \) unchanged, such as centring, scaling and
orthogonal polynomials, are statistically neutral but can be numerically decisive.

## Rank deficiency and the singular value decomposition

QR with positive diagonal needs full column rank. For rank-deficient or nearly
rank-deficient \( \X \), the most informative tool is the **singular value
decomposition**
\[
\X=\bU\bD\V\T=\bU_r\bD_r\V_r\T ,
\]
where \( \bU \) (\( n\times n \)) and \( \V \) (\( p\times p \)) are orthogonal, \( \bD \) is
\( n\times p \) with nonnegative diagonal entries \( \sigma_1\ge\dots\ge\sigma_r>0 \), and the
subscript \( r \) keeps only the first \( r \) columns. Then:

- \( \bU_r \) is an orthonormal basis for \( \C(\X) \), so \( \M=\bU_r\bU_r\T \);

- the last \( p-r \) columns of \( \V \) form an orthonormal basis for \( \Null(\X) \);

- the Moore–Penrose inverse is \( \X^+=\V_r\bD_r^{-1}\bU_r\T \), and \( \M=\X\X^+ \).

::: {#prp-proj-min-norm}
[Minimum-norm least squares]

\( \hbeta^+=\X^+\y \) is a least squares estimate, and it has the smallest Euclidean norm
among all least squares estimates.
:::

::: {.proof}
\( \X\X^+\y=\bU_r\bD_r\V_r\T\V_r\bD_r^{-1}\bU_r\T\y=\bU_r\bU_r\T\y=\M\y \), so \( \hbeta^+ \) is a
least squares estimate. It lies in \( \C(\V_r)=\C(\X\T)=\Null(\X)\perpc \). Every other
least squares estimate is \( \hbeta^++\bv \) with \( \bv\in\Null(\X) \)
(@thm-proj-normal-equations), and
\( \norm{\hbeta^++\bv}^2=\norm{\hbeta^+}^2+\norm{\bv}^2 \).
:::

In @exm-proj-ginverse-numeric, the Moore–Penrose row is this minimum-norm
solution. In floating point, “\( \sigma_j>0 \)” becomes “\( \sigma_j \) larger
than a tolerance”, and choosing the tolerance decides the numerical rank. That choice is
a statistical decision in disguise. [Chapter 10](../ch10-computation/index.html) discusses it, together with
pivoted QR (@thm-cmp-rank-revealing), which gives most of the SVD's rank-revealing ability at lower cost.

## Exercises

### A. Check your understanding

::: {#exr-proj-qr-identities}
[A1]

Prove the four identities in @eq-proj-qr-ls.
:::

### B. Practice

::: {#exr-proj-longley-scaling}
[B1]

Compute \( \kappa(\X) \) for the Longley data (with intercept) as given, after scaling every
non-intercept column to unit Euclidean norm, and after centring and scaling. Confirm that
the fitted values agree to the precision you expect in all three cases. Which part of the
original condition number is “just units”?
:::

::: {#exr-proj-gram-schmidt-stability}
[B2]

Implement classical Gram–Schmidt as in @prp-proj-gram-schmidt, and the
*modified* version that subtracts each projection from the working vector before
computing the next inner product. For the degree-10 polynomial design of
@exm-proj-qr-polynomial, compare \( \norm{\Q\T\Q-\I} \) for the two and for
`numpy.linalg.qr`.
:::

### C. Going deeper

::: {#exr-proj-cond-centering}
[C1]

Let \( \X=[\bone,\x] \) with \( \bar x=c \) and \( S_{xx}=s^2 n \). Find \( \kappa(\X) \) as a function of \( c \),
\( s \) and \( n \), and show that it grows without bound as \( c\to\infty \) with \( s \) fixed. Show
that centring the regressor, a reparameterization, reduces the condition number to
\( \max(1,s)/\min(1,s) \).
:::

::: {.solution}
Write \( \x=c\bone+\mathbf{d} \) with \( \mathbf{d}\perp\bone \) and
\( \norm{\mathbf{d}}^2=ns^2 \). Then
\( \X\T\X=n\begin{psmallmatrix}1&c\\c&c^2+s^2\end{psmallmatrix} \). The matrix in parentheses
has trace \( T=1+c^2+s^2 \) and determinant \( s^2 \), so its eigenvalues are
\( \lambda_{\pm}=\bigl(T\pm\sqrt{T^2-4s^2}\bigr)/2 \), with \( \lambda_+\lambda_-=s^2 \) and
\( \lambda_+\ge T/2 \). Hence
\[
\kappa(\X)^2=\frac{\lambda_+}{\lambda_-}=\frac{\lambda_+^2}{s^2}\;\ge\;\frac{T^2}{4s^2}
\;\to\;\infty\quad\text{as } c\to\infty .
\]
After centring, \( \X=[\bone,\mathbf{d}] \) has orthogonal columns of lengths \( \sqrt n \) and
\( s\sqrt n \), so its singular values are these lengths and
\( \kappa=\max(1,s)/\min(1,s) \).
:::
