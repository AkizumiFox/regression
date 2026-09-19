# Correlation matrices

Covariances carry units: the covariance of income in dollars with age in years is in
dollar-years. To compare the strength of linear association across pairs of variables,
each covariance is divided by the two standard deviations.

::: {#def-rv-correlation}
[Correlation matrix]

Let \( \Y \) have covariance matrix \( \bSigma=(\sigma_{ij}) \) with every \( \sigma_{ii}>0 \), and let
\( \bD_\sigma=\diag(\sigma_{11}^{1/2},\dots,\sigma_{pp}^{1/2}) \) be the diagonal matrix of standard
deviations. The **correlation matrix** of \( \Y \) is
\[
\R=\bD_\sigma^{-1}\bSigma\bD_\sigma^{-1},
\qquad\text{with entries}\qquad
\rho_{ij}=\frac{\sigma_{ij}}{\sqrt{\sigma_{ii}\sigma_{jj}}}=\operatorname{corr}(Y_i,Y_j).
\]
Equivalently \( \bSigma=\bD_\sigma\R\bD_\sigma \). For random vectors \( \bU \) and \( \V \) with
standard-deviation matrices \( \bD_U \) and \( \bD_V \), the cross-correlation matrix is
\( \bD_U^{-1}\Cov(\bU,\V)\bD_V^{-1} \).
:::

Some books write \( \mathbf{P}_\rho \) for the correlation matrix. We use \( \R \) in this chapter, and
\( \hat{\R} \) for its sample counterpart in [Section 2.7](07-sample.html).

::: {#prp-rv-correlation}
[Properties of correlation matrices]

Let \( \R \) be the correlation matrix of \( \Y \).

::: {.enumerate options="label=(\alph*)"}
1. \( \R \) is the covariance matrix of the standardized vector \( \Z=\bD_\sigma^{-1}(\Y-\bmu) \).
           Hence \( \R \) is nonnegative definite with unit diagonal, and \( \tr(\R)=p \).

2. \( \lvert\rho_{ij}\rvert\le1 \), with \( \lvert\rho_{ij}\rvert=1 \) iff \( Z_j=\rho_{ij}Z_i \) with
           probability one.

3. \( \rank(\R)=\rank(\bSigma) \). In particular, \( \R \) is positive definite iff \( \bSigma \) is.

4. If \( \Y^*=\bD\Y+\mathbf{c} \) with \( \bD=\diag(d_1,\dots,d_p) \) and every \( d_i\ne0 \), then
           \( \rho^*_{ij}=\operatorname{sign}(d_id_j)\,\rho_{ij} \). Positive rescaling and shifting
           leave \( \R \) unchanged.

5. Conversely, every nonnegative definite matrix with unit diagonal is the correlation
           matrix of some random vector.
:::

:::

::: {.proof}
(a) @thm-rv-linear gives \( \Cov(\Z)=\bD_\sigma^{-1}\bSigma\bD_\sigma^{-1}=\R \), and
@thm-rv-cov-nnd(a) applies. (b) The \( 2\times2 \) matrix of rows and columns \( i,j \) of \( \R \)
is the covariance matrix of \( (Z_i,Z_j) \), so it is nonnegative definite and its determinant
\( 1-\rho_{ij}^2 \) is nonnegative. If \( \rho_{ij}=\pm1 \), the vector \( (\rho_{ij},-1)\T \) is in the null
space of that \( 2\times2 \) matrix, so \( \rho_{ij}Z_i-Z_j=0 \) with probability one by
@thm-rv-cov-nnd(c). Conversely, \( Z_j=cZ_i \) with \( \Var(Z_i)=\Var(Z_j)=1 \) forces \( c=\pm1 \)
and \( \rho_{ij}=c \). (c) \( \bD_\sigma \) is nonsingular. (d) \( \Cov(\Y^*)=\bD\bSigma\bD \) has entries
\( d_id_j\sigma_{ij} \) and standard deviations \( \lvert d_i\rvert\sigma_{ii}^{1/2} \). (e) Such a
matrix is a covariance matrix by @thm-rv-cov-nnd(b), and its variances are one.
:::

Part (d) is the reason correlation is the natural scale-free summary. Part (e) is less
obvious than it looks. A symmetric matrix with unit diagonal and off-diagonal entries in
\( [-1,1] \) need not be a correlation matrix, because the correlations must be jointly
consistent.

::: {#prp-rv-three-correlations}
[Three correlations]

Let \( \lvert\rho_{12}\rvert<1 \) and \( \lvert\rho_{13}\rvert<1 \). The matrix
\[
\R=\begin{pmatrix}1&\rho_{12}&\rho_{13}\\\rho_{12}&1&\rho_{23}\\\rho_{13}&\rho_{23}&1\end{pmatrix}
\]
is a correlation matrix iff
\[
\rho_{12}\rho_{13}-\sqrt{(1-\rho_{12}^2)(1-\rho_{13}^2)}
\;\le\;\rho_{23}\;\le\;
\rho_{12}\rho_{13}+\sqrt{(1-\rho_{12}^2)(1-\rho_{13}^2)} .
\]
:::

::: {.proof}
By @prp-rv-correlation(e) we must decide when \( \R \) is nonnegative definite. Write
\( \mathbf{r}=(\rho_{12},\rho_{13})\T \) and let \( \R_2 \) be the lower right \( 2\times2 \) block. The
nonsingular matrix \( \bT=\begin{psmallmatrix}1&\bzero\T\\-\mathbf{r}&\I_2\end{psmallmatrix} \) gives
\[
\bT\R\bT\T=\begin{pmatrix}1&\bzero\T\\\bzero&\R_2-\mathbf{r}\mathbf{r}\T\end{pmatrix},
\]
and \( \mathbf{a}\T\R\mathbf{a}=\mathbf{b}\T(\bT\R\bT\T)\mathbf{b} \) with \( \mathbf{b}=(\bT\T)^{-1}\mathbf{a} \) ranges over all
vectors as \( \mathbf{a} \) does. So \( \R \) is nonnegative definite iff the Schur complement
\[
\R_2-\mathbf{r}\mathbf{r}\T=\begin{pmatrix}1-\rho_{12}^2&\rho_{23}-\rho_{12}\rho_{13}\\
\rho_{23}-\rho_{12}\rho_{13}&1-\rho_{13}^2\end{pmatrix}
\]
is. A symmetric \( 2\times2 \) matrix with positive diagonal is nonnegative definite iff its
determinant is nonnegative, because its eigenvalues have positive sum and their product is
the determinant. The determinant is
\( (1-\rho_{12}^2)(1-\rho_{13}^2)-(\rho_{23}-\rho_{12}\rho_{13})^2 \), and it is nonnegative exactly on
the stated interval.
:::

If \( Y_2 \) and \( Y_3 \) both correlate \( 0.8 \) with \( Y_1 \), the interval is \( [0.28,1] \). Two variables that
are strongly related to a third cannot be negatively correlated with each other, nor
even uncorrelated. The quantity
\( (\rho_{23}-\rho_{12}\rho_{13})/\sqrt{(1-\rho_{12}^2)(1-\rho_{13}^2)} \), which the proof forces into
\( [-1,1] \), is the correlation of \( Y_2 \) and \( Y_3 \) after the linear effect of \( Y_1 \) has been removed
from both. [Section 2.5](05-partitioned.html) develops this residualization, and
@def-mvn-partial-correlation names the result the partial correlation.

Correlation matrices do not transform simply under linear maps: the correlation matrix of
\( \A\Y \) must be recomputed from \( \A\bSigma\A\T \). Theory is therefore done with covariance
matrices, and correlation matrices are used for description.

## Exercises

### A. Check your understanding

::: {#exr-rv-equicorrelation}
[A1]

Show that \( (1-\rho)\I_p+\rho\bone\bone\T \) is a correlation matrix iff \( -1/(p-1)\le\rho\le1 \). Explain the
lower bound through the variance of \( \sum_iZ_i \).
:::

::: {.solution}
The matrix sends \( \bone \) to \( (1+(p-1)\rho)\bone \) and every
\( \bu\perp\bone \) to \( (1-\rho)\bu \). Its eigenvalues are \( 1+(p-1)\rho \) (once) and \( 1-\rho \) (\( p-1 \) times). By
@prp-rv-correlation(e) it is a correlation matrix iff both are nonnegative. If \( \Z \) has this
correlation matrix, \( \Var(\sum_iZ_i)=\bone\T\R\bone=p[1+(p-1)\rho] \), which must be nonnegative.
Many variables cannot all be strongly negatively correlated with one another, because their sum would
then have negative variance.
:::

### B. Practice

::: {#exr-rv-hadamard}
[B1]

Show that \( \det(\bSigma)\le\prod_i\sigma_{ii} \) for every covariance matrix, with equality for positive
definite \( \bSigma \) iff \( \bSigma \) is diagonal. *Hint:* factor out \( \bD_\sigma \) and use the eigenvalues
of \( \R \).
:::

::: {.solution}
If some \( \sigma_{ii}=0 \), the \( i \)th row of \( \bSigma \) vanishes
(apply @thm-rv-cov-nnd(c) with \( \mathbf{a} \) the \( i \)th coordinate vector) and both sides
are zero. Otherwise \( \det\bSigma=\det(\bD_\sigma)^2\det\R=\prod_i\sigma_{ii}\cdot\det\R \). The eigenvalues
\( \lambda_i\ge0 \) of \( \R \) sum to \( \tr\R=p \), so by the inequality of arithmetic and geometric means
\( \det\R=\prod_i\lambda_i\le(p^{-1}\sum_i\lambda_i)^p=1 \). Equality holds iff every \( \lambda_i=1 \), that is,
iff \( \R=\I \), that is, iff \( \bSigma \) is diagonal.
:::

### C. Going deeper

::: {#exr-rv-correlation-angles}
[C1]

(a) Show that a \( p\times p \) matrix is a correlation matrix iff it is the Gram matrix
\( (\bu_i\T\bu_j) \) of \( p \) unit vectors in \( \Real^p \). (b) For \( p=3 \), write
\( \rho_{ij}=\cos\vartheta_{ij} \) with \( \vartheta_{ij}\in[0,\pi] \). Show that
@prp-rv-three-correlations is equivalent to the spherical triangle inequalities
\( \lvert\vartheta_{12}-\vartheta_{13}\rvert\le\vartheta_{23}\le\min(\vartheta_{12}+\vartheta_{13},\,2\pi-\vartheta_{12}-\vartheta_{13}) \).
:::
