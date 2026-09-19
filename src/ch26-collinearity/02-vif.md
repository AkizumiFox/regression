# Variance inflation factors

Analysts usually want something coarser than the eigenvalues: one number per coefficient saying how
much its precision has suffered. [Chapter 6](../ch06-projections/index.html) supplied that number as a
by-product of the Frisch–Waugh–Lovell theorem. This section states what it measures and what it does
not, and extends it to a group of columns such as the indicators of a factor.

## Definition and meaning

Throughout, the model contains an intercept, \( \x_j \) is one of the other columns, and \( \X_{(j)} \)
is \( \X \) with \( \x_j \) removed, with projection \( \M_{(j)} \). Write
\( S_{jj}=\norm{\x_j-\bar x_j\bone}^2 \) for the centred sum of squares of \( \x_j \).

::: {#def-col-vif}
[Variance inflation factor]

Let \( R_j^2 \) be the coefficient of determination of the least squares regression of \( \x_j \) on the
columns of \( \X_{(j)} \), which include \( \bone \), and suppose \( R_j^2<1 \). The **variance inflation
factor** of \( \x_j \) is
\[
\text{VIF}_j=\frac1{1-R_j^2},
\]
and \( \text{TOL}_j=1-R_j^2=1/\text{VIF}_j \) is its **tolerance**.
:::

The identity @eq-proj-vif-preview, proved with the Frisch–Waugh–Lovell theorem, reads
\( \Var(\hat{\beta}_j)=\sigma^2\text{VIF}_j/S_{jj} \). The next proposition collects what follows from it.

::: {#prp-col-vif}
[What a variance inflation factor measures]

::: {.enumerate options="label=(\alph*)"}
1. \( \norm{(\I-\M_{(j)})\x_j}^2=S_{jj}/\text{VIF}_j \), and \( \text{VIF}_j \) is the ratio of \( \Var(\hat{\beta}_j) \) to
   \( \sigma^2/S_{jj} \), the variance \( \hat{\beta}_j \) would have if the centred \( \x_j \) were orthogonal to every other
   centred column. By @exr-dep-vif-one, \( \text{VIF}_j\ge1 \).

2. \( \text{VIF}_j \) depends on \( \X \) only through \( \x_j \) and \( \C(\X_{(j)}) \). It is unchanged when \( \x_j \) is
   replaced by \( a\x_j+b\bone \) with \( a\ne0 \), and when \( \X_{(j)} \) is replaced by any matrix with the same column
   space.

3. If \( \R \) is the correlation matrix of the non-intercept columns, then \( \text{VIF}_j=[\R^{-1}]_{jj} \).

4. The \( t \) interval for \( \beta_j \) at any level is \( \sqrt{\text{VIF}_j} \) times as long, on average, as in the orthogonal
   design of (a), and the noncentrality of the \( F \) test of \( \beta_j=0 \) (the square of the \( t \)
   statistic's) is divided by \( \text{VIF}_j \):
   \[
   \frac{\beta_j^2}{\Var(\hat{\beta}_j)}=\frac{\beta_j^2S_{jj}}{\sigma^2}\cdot\frac1{\text{VIF}_j}.
   \]
:::

:::

::: {.proof}
(a) Since \( \bone\in\C(\X_{(j)}) \), the vector \( (\I-\M_{(j)})\x_j \) is the residual of the regression of \( \x_j \) on \( \X_{(j)} \), and
its squared length is \( S_{jj}(1-R_j^2) \) by the definition of \( R^2 \) (@exr-proj-vif). The variance ratio is
@eq-proj-vif-preview divided by \( \sigma^2/S_{jj} \), which is the variance in the orthogonal case by
@prp-opt-orthogonal-design applied to the model with centred columns.

(b) By (a), \( 1/\text{VIF}_j=\norm{(\I-\M_{(j)})\x_j}^2/S_{jj} \), and \( \M_{(j)} \) depends only on \( \C(\X_{(j)}) \). Replacing
\( \x_j \) by \( a\x_j+b\bone \) multiplies the numerator by \( a^2 \), because \( (\I-\M_{(j)})\bone=\bzero \), and multiplies
\( S_{jj} \) by \( a^2 \).

(c) This is @exr-dep-vif-inverse.

(d) The \( t \) interval has half-length \( t_{\alpha/2}(n-p)\,s\sqrt{\text{VIF}_j/S_{jj}} \), and \( s \) has the same distribution
in both designs. The noncentrality is \( \beta_j^2/\Var(\hat{\beta}_j) \) with the variance from (a) (@thm-glh-t-test).
:::

By (b) the factor is scale-free: a property of the angle between the centred \( \x_j \) and the space of
the other centred columns, whose squared cosine is \( R_j^2 \) (@thm-proj-r2-cosine). By (a), with the spread per observation
fixed, \( \hat{\beta}_j \) is as precise as in an orthogonal design with \( n/\text{VIF}_j \) observations.

::: {.idea}
A variance inflation factor of \( 10 \) costs \( \hat{\beta}_j \) as much precision as discarding nine observations in ten
from an orthogonal design with the same spread.
:::

## Reading variance inflation factors

Rules of thumb flag a regressor when \( \text{VIF}_j \) exceeds \( 10 \), or \( 5 \), or \( 4 \). The thresholds are arbitrary,
and they measure the wrong thing. The factor compares the design with a hypothetical orthogonal
design of the same size; it says nothing about whether the actual standard error is small enough for
the purpose, which depends also on \( \sigma^2 \) and on \( S_{jj} \), hence on \( n \). As \( n \) grows with the regressors drawn
from a fixed population, \( \text{VIF}_j \) settles to a population value while \( S_{jj} \) grows like \( n \), so a factor that is
ruinous at \( n=20 \) can be irrelevant at \( n=20{,}000 \). O'Brien (2007) makes this case in detail. A large
factor is a reason to look at the interval, never a verdict by itself.

::: {#exm-col-macro-vif}
[Large factors, precise estimate]

The quarterly US macroeconomic series of @exm-mat-macro-svd run from 1959 to 2009, \( n=203 \).
Regress real personal consumption on real disposable income, population, the consumer price index,
the money stock M1, the three-month Treasury bill rate and the unemployment rate, with an intercept.
The variance inflation factors are

| Regressor | income | population | CPI | M1 | T-bill | unemployment |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| \( \text{VIF}_j \) | \( 154.5 \) | \( 183.7 \) | \( 140.4 \) | \( 97.0 \) | \( 2.9 \) | \( 1.3 \) |

By every rule of thumb, four regressors are badly collinear. Yet the coefficient of income is
\( 1.0659 \) with standard error \( 0.0210 \), a \( t \) statistic of \( 50.9 \): the standard error is
\( \sqrt{154.5}\approx12.4 \) times what an orthogonal design would give, and still small, because income
varies widely and consumption tracks it closely. (The errors are autocorrelated, so all the standard
errors are too small, [Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html); the
factors depend on \( \X \) alone.)
:::

```{.python .run #cell-vif-macro}
import numpy as np
import statsmodels.api as sm

def vifs(Z):
    """VIFs of the columns of Z (an intercept is always included): diagonal of R^{-1}."""
    return np.diag(np.linalg.inv(np.corrcoef(Z, rowvar=False)))

macro = sm.datasets.macrodata.load_pandas().data
names = ["realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
Z = macro[names].to_numpy()
fit = sm.OLS(macro["realcons"].to_numpy(), sm.add_constant(Z)).fit()
for name, v, b, se in zip(names, vifs(Z), fit.params[1:], fit.bse[1:]):
    print(f"{name:9s} VIF {v:7.1f}   estimate {b:9.4f}   se {se:7.4f}   t {b / se:6.1f}")
```

## Collinearity that a reparameterization removes

Some large factors are produced by the parameterization, not by the data. Polynomial and interaction
terms built from regressors whose origin lies far from the data are the standard case. Such
collinearity is sometimes called *nonessential*.

::: {#exm-col-quadratic}
[A quadratic far from the origin]

Take \( n=21 \) equally spaced values \( x_i=10,10.5,\dots,20 \) and the quadratic model
\( \E(Y)=\beta_0+\beta_1x+\beta_2x^2 \). Over this range \( x \) and \( x^2 \) are almost perfectly correlated:
\( R^2=0.9920 \), and both have \( \text{VIF}=124.6 \). Rewrite the model in the centred regressors
\( x-15 \) and \( (x-15)^2 \). The fitted values are identical, and both variance inflation factors are now exactly \( 1 \),
because the design is symmetric about \( 15 \) and an odd function of \( x-15 \) is orthogonal to an even one
once both are centred.

What changed is the question. In the raw coding \( \beta_1 \) is the slope of the mean function at \( x=0 \),
ten units below the data, with standard deviation \( 0.804\,\sigma \). In the centred coding the coefficient of
\( x-15 \) is the slope at \( x=15 \), with standard deviation \( 0.0721\,\sigma \), and the raw coding estimates the
same slope, \( \beta_1+30\beta_2 \), with exactly the same standard deviation. The large factors described
the difficulty of a question, the slope at the origin, that the analyst probably did not mean to ask (@exr-col-quadratic-vif
derives \( 124.6 \) by hand).
:::

```{.python .run #cell-vif-quadratic}
x = np.linspace(10, 20, 21)                           # a regressor far from its origin

raw = np.column_stack([x, x**2])
cen = np.column_stack([x - x.mean(), (x - x.mean()) ** 2])
print("VIFs, raw     :", np.round(vifs(raw), 1))
print("VIFs, centred :", np.round(vifs(cen), 3))

X_raw = np.column_stack([np.ones_like(x), raw])
X_cen = np.column_stack([np.ones_like(x), cen])
G_raw, G_cen = np.linalg.inv(X_raw.T @ X_raw), np.linalg.inv(X_cen.T @ X_cen)
a_mean = np.array([0, 1, 2 * x.mean()])               # slope at x = 15 in the raw coefficients
print("sd of the slope at 0 (raw b1)       :", round(np.sqrt(G_raw[1, 1]), 3))
print("sd of the slope at 15, raw coding   :", round(np.sqrt(a_mean @ G_raw @ a_mean), 4))
print("sd of the slope at 15, centred (b1) :", round(np.sqrt(G_cen[1, 1]), 4))
```

## Groups of columns

A factor enters the model through several indicator columns, a polynomial through several powers.
The individual coefficients of such a block depend on the coding (@def-est-coding), and so do their
variance inflation factors. What is intrinsic is the joint precision of the block, which the next
definition, due to Fox and Monette (1992), measures.

For the rest of the section, split the non-intercept columns of \( \X \) into \( \Z_1 \) (\( n\times q \)) and \( \Z_2 \)
(\( n\times(k-q) \)), let \( \tilde{\Z}_1 \) and \( \tilde{\Z}_2 \) be the centred versions, and let \( \tilde{\M}_2 \) project onto \( \C(\tilde{\Z}_2) \).
Assume \( \X=[\bone,\Z_1,\Z_2] \) has full column rank.

::: {#def-col-gvif}
[Generalized variance inflation factor]

The **generalized variance inflation factor** of the block \( \Z_1 \) is
\[
\text{GVIF}_1=\frac{\det\bigl(\tilde{\Z}_1\T\tilde{\Z}_1\bigr)}{\det\bigl(\tilde{\Z}_1\T(\I-\tilde{\M}_2)\tilde{\Z}_1\bigr)} .
\]
:::

The **canonical correlations** between the two blocks are defined through orthonormal bases
\( \Q_1 \) of \( \C(\tilde{\Z}_1) \) and \( \Q_2 \) of \( \C(\tilde{\Z}_2) \): they are the singular values
\( \rho_1\ge\rho_2\ge\dots\ge\rho_m \), \( m=\min(q,k-q) \), of the \( q\times(k-q) \) matrix \( \Q_1\T\Q_2 \). The largest,
\( \rho_1 \), is the largest correlation between a linear combination of the centred columns of \( \Z_1 \) and a
linear combination of the centred columns of \( \Z_2 \) (@exr-col-canonical-max).

::: {#prp-col-gvif}
[Properties of the generalized factor]

Let \( \bbeta_1 \) be the coefficients of \( \Z_1 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \det\Cov(\hbeta_1)=\text{GVIF}_1\cdot\det\Cov_\perp \), where \( \Cov_\perp=\sigma^2(\tilde{\Z}_1\T\tilde{\Z}_1)^{-1} \) is the covariance
   \( \hbeta_1 \) would have if \( \C(\tilde{\Z}_1) \) were orthogonal to \( \C(\tilde{\Z}_2) \). So \( \text{GVIF}_1 \) is the squared ratio of the volumes of
   the joint confidence ellipsoids for \( \bbeta_1 \) in the two designs.

2. \( \text{GVIF}_1=\prod_{i=1}^m1/(1-\rho_i^2) \).

3. \( \text{GVIF}_1 \) is unchanged when \( \Z_1 \) is replaced by \( \Z_1\A+\bone\mathbf{f}\T \) and \( \Z_2 \) by
   \( \Z_2\B+\bone\mathbf{g}\T \), for nonsingular \( \A \) and \( \B \) and any vectors \( \mathbf{f} \) and \( \mathbf{g} \). In particular it
   does not depend on the coding of a factor. It is symmetric:
   \( \text{GVIF}_1=\text{GVIF}_2 \).

4. If \( q=1 \), \( \text{GVIF}_1 \) is the variance inflation factor of @def-col-vif. In general, with \( \R \) the correlation
   matrix of \( [\Z_1,\Z_2] \) partitioned conformably,
   \( \text{GVIF}_1=\det\R_{11}\det\R_{22}/\det\R \).
:::

:::

::: {.proof}
(a) Since \( \C(\bone,\Z_2)=\C(\bone)\dirsum\C(\tilde{\Z}_2) \) with orthogonal summands, the projection onto \( \C(\bone,\Z_2) \)
is \( n^{-1}\bone\bone\T+\tilde{\M}_2 \) (@thm-proj-sum), and the residual of \( \Z_1 \) after it is
\( (\I-\tilde{\M}_2)\tilde{\Z}_1 \). By @thm-proj-fwl and @eq-proj-partitioned-inverse,
\( \Cov(\hbeta_1)=\sigma^2\bigl(\tilde{\Z}_1\T(\I-\tilde{\M}_2)\tilde{\Z}_1\bigr)^{-1} \). If \( \C(\tilde{\Z}_1)\perp\C(\tilde{\Z}_2) \), then
\( \tilde{\M}_2\tilde{\Z}_1=\mathbf{0} \) and the covariance is \( \Cov_\perp \). Take the ratio of determinants. By
@thm-ci-ellipsoid(b), the volume of the confidence ellipsoid is proportional to \( (\det\bSigma)^{1/2} \), where
\( \bSigma \) is the covariance matrix, with the same constant in both designs.

(b) Centring annihilates \( \bone\mathbf{f}\T \) and commutes with multiplication on the right, so replacing \( \Z_1 \)
by \( \Z_1\A+\bone\mathbf{f}\T \) replaces \( \tilde{\Z}_1 \) by \( \tilde{\Z}_1\A \), and both determinants in @def-col-gvif are
multiplied by \( \det(\A)^2 \). For the same reason \( \Z_2\mapsto\Z_2\B+\bone\mathbf{g}\T \) replaces \( \tilde{\Z}_2 \) by
\( \tilde{\Z}_2\B \), which leaves \( \C(\tilde{\Z}_2) \), and so \( \tilde{\M}_2 \), unchanged. We may therefore take
\( \tilde{\Z}_1=\Q_1 \). Then \( \tilde{\M}_2=\Q_2\Q_2\T \), and with \( \mathbf{C}=\Q_1\T\Q_2 \),
\[
\text{GVIF}_1=\frac{\det(\Q_1\T\Q_1)}{\det(\Q_1\T\Q_1-\mathbf{C}\mathbf{C}\T)}=\frac1{\det(\I_q-\mathbf{C}\mathbf{C}\T)} .
\]
The eigenvalues of \( \mathbf{C}\mathbf{C}\T \) are \( \rho_1^2,\dots,\rho_m^2 \) and \( q-m \) zeros, so the determinant is
\( \prod_i(1-\rho_i^2) \).

(c) The invariance was shown in the proof of (b). Symmetry follows from (b), because \( \Q_1\T\Q_2 \) and
\( \Q_2\T\Q_1 \) have the same singular values.

(d) For \( q=1 \) the ratio is \( S_{jj}/\norm{(\I-\M_{(j)})\x_j}^2=\text{VIF}_j \) by @prp-col-vif(a). For the determinant form, let
\( \bD \) be the diagonal matrix with \( \R=\bD[\tilde{\Z}_1,\tilde{\Z}_2]\T[\tilde{\Z}_1,\tilde{\Z}_2]\bD \), with diagonal blocks \( \bD_1,\bD_2 \). By
@thm-mat-block-determinant, \( \det\R=\det\R_{22}\det(\R_{11}-\R_{12}\R_{22}^{-1}\R_{21}) \), and the Schur complement is
\( \bD_1\tilde{\Z}_1\T(\I-\tilde{\M}_2)\tilde{\Z}_1\bD_1 \). Since \( \R_{11}=\bD_1\tilde{\Z}_1\T\tilde{\Z}_1\bD_1 \), the factors \( \det(\bD_1)^2 \) cancel in
\( \det\R_{11}\det\R_{22}/\det\R \), which leaves \( \text{GVIF}_1 \).
:::

A generalized factor inflates a \( q \)-dimensional volume, so for comparison with single columns Fox and
Monette suggest \( \text{GVIF}_1^{1/(2q)} \), the geometric mean of the factors by which the axes of the confidence
ellipsoid lengthen; for \( q=1 \) it is \( \sqrt{\text{VIF}_j} \).

::: {#exm-col-grunfeld}
[A factor that explains a regressor]

Grunfeld's panel (@exm-proj-grunfeld) has \( 220 \) firm-years of gross investment, market value
and capital stock for \( 11 \) firms. Regress investment on value, capital and the firm, entered as
\( 10 \) indicator columns with the first firm as reference. The generalized factors are

| Term | columns \( q \) | \( \text{GVIF} \) | \( \text{GVIF}^{1/(2q)} \) |
|:---|:---:|:---:|:---:|
| value | \( 1 \) | \( 18.32 \) | \( 4.280 \) |
| capital | \( 1 \) | \( 2.04 \) | |
| firm | \( 10 \) | \( 23.61 \) | \( 1.171 \) |

The canonical correlations between the firm block and the pair (value, capital) are \( 0.9694 \) and
\( 0.5459 \), and \( 1/\bigl((1-0.9694^2)(1-0.5459^2)\bigr)\approx23.6 \) as @prp-col-gvif(b) requires. Market value
is largely a characteristic of the firm. With the firm indicators in the model, the coefficient of value is
identified by movements within firms (the within estimator of
[Section 6.6](../ch06-projections/06-fwl.html#sec-proj-within)), and the part of value that neither firm nor
capital explains is \( 1/18.32 \), about \( 5.5 \) per cent, of its variation.

The ordinary factors of the ten indicator columns range from \( 1.8 \) to \( 16.8 \); with the last firm as
reference the largest is \( 12.7 \). The generalized factors are the same under every coding, including
sum-to-zero coding (@prp-col-gvif(c)).
:::

```{.python .run #cell-vif-gvif}
import pandas as pd
grun = sm.datasets.grunfeld.load_pandas().data
D = pd.get_dummies(grun["firm"], drop_first=True).to_numpy(dtype=float)   # 10 indicator columns
Zg = np.column_stack([grun["value"], grun["capital"], D])

def gvif(Z, block):
    """det(R_11) det(R_22) / det(R) for the columns in block against the rest."""
    R = np.corrcoef(Z, rowvar=False)
    rest = [k for k in range(Z.shape[1]) if k not in block]
    return (np.linalg.det(R[np.ix_(block, block)]) * np.linalg.det(R[np.ix_(rest, rest)])
            / np.linalg.det(R))

firm = list(range(2, Zg.shape[1]))
print("GVIF value  :", round(gvif(Zg, [0]), 2))
print("GVIF capital:", round(gvif(Zg, [1]), 2))
print("GVIF firm   :", round(gvif(Zg, firm), 2), " per dimension:", round(gvif(Zg, firm) ** (1 / 20), 3))
print("VIFs of the firm indicators:", np.round(vifs(Zg)[2:], 1))
```

## Exercises

### A. Check your understanding

::: {#exr-col-vif-thresholds}
[A1]

Which values of \( R_j^2 \) correspond to variance inflation factors of \( 4 \), \( 5 \) and \( 10 \)? By what factor does
each multiply the length of a confidence interval for \( \beta_j \), relative to the orthogonal design of
@prp-col-vif(a)?
:::

::: {#exr-col-vif-two}
[A2]

For a model with an intercept and two regressors with sample correlation \( r \), show that
\( \text{VIF}_1=\text{VIF}_2=1/(1-r^2) \). What does @prp-col-vif(a) say about the design of @exm-col-two-regressors?
:::

### B. Practice

::: {#exr-col-vif-monotone}
[B1]

Show that appending a column \( \bw\notin\C(\X) \) to \( \X \) cannot decrease the variance inflation factor of any
column already in the model, and that it leaves \( \text{VIF}_j \) unchanged iff the residual of \( \bw \) on \( \X_{(j)} \) is
orthogonal to the residual of \( \x_j \) on \( \X_{(j)} \). Is \( \Var(\hat{\beta}_j) \) also nondecreasing?
:::

::: {.solution}
Appending a column enlarges \( \C(\X_{(j)}) \), so \( \norm{(\I-\M_{(j)})\x_j}^2 \) cannot increase, and by @prp-col-vif(a)
\( \text{VIF}_j \) cannot decrease. Let \( \tilde{\bw}\ne\bzero \) be the residual of \( \bw \) on \( \X_{(j)} \). The enlarged
space is \( \C(\X_{(j)})\dirsum\C(\tilde{\bw}) \) with orthogonal summands, so the squared residual of \( \x_j \) drops by
\( (\tilde{\bw}\T\x_j)^2/\norm{\tilde{\bw}}^2 \), and \( \tilde{\bw}\T\x_j=\tilde{\bw}\T(\I-\M_{(j)})\x_j \). So \( \text{VIF}_j \) is unchanged iff
the two residuals are orthogonal. Since \( S_{jj} \) does not change, \( \Var(\hat{\beta}_j)=\sigma^2\text{VIF}_j/S_{jj} \) is
nondecreasing too. This is the variance half of @thm-dep-overfit.
:::

::: {#exr-col-quadratic-vif}
[B2]

Let \( x_i=m+u_i \), where the \( u_i \) are symmetric about zero. Write \( v \) for the mean of the \( u_i^2 \) and \( w \) for the
variance of the \( u_i^2 \) (both with divisor \( n \)). Show that the variance inflation factor of \( x \) and of \( x^2 \) in the
quadratic model is
\[
\text{VIF}=1+\frac{4m^2v}{w}.
\]
Check the value in @exm-col-quadratic.
:::

::: {.solution}
With one other regressor, \( \text{VIF}=1/(1-r^2) \), where \( r \) is the correlation of \( x \) and \( x^2 \) (@exr-col-vif-two).
Since \( x^2=m^2+2mu+u^2 \), the covariance of \( x \) and \( x^2 \) is that of \( u \) and \( 2mu+u^2 \), which is \( 2mv \) because
the third moment of \( u \) vanishes. The variance of \( x^2 \) is \( 4m^2v+w \), since \( u \) and \( u^2 \) are uncorrelated. So
\( r^2=4m^2v^2/\bigl(v(4m^2v+w)\bigr)=4m^2v/(4m^2v+w) \) and \( 1/(1-r^2)=1+4m^2v/w \). In the example, \( m=15 \) and \( u \) runs
over \( -5,-4.5,\dots,5 \): \( v=192.5/21\approx9.167 \), the mean of \( u^4 \) is \( 3166.6/21\approx150.79 \), so
\( w\approx150.79-84.03=66.76 \) and \( \text{VIF}\approx1+900\cdot9.167/66.76\approx124.6 \).
:::

::: {#exr-col-gvif-member}
[B3]

Let \( \x_j \) be a column of the block \( \Z_1 \). Show that \( \text{GVIF}_1\ge1/(1-R^2_{j\mid2}) \), where \( R^2_{j\mid2} \) is the
coefficient of determination of \( \x_j \) regressed on \( \bone \) and \( \Z_2 \) alone. Why is there no such inequality
between \( \text{GVIF}_1 \) and the ordinary \( \text{VIF}_j \) computed in the full model?
:::

::: {.solution}
\( R_{j\mid2} \) is the correlation between the centred \( \x_j \) and its projection on \( \C(\tilde{\Z}_2) \), which is a correlation
between a vector \( \Q_1\mathbf{s} \) of \( \C(\tilde{\Z}_1) \) and a vector \( \Q_2\mathbf{t} \) of \( \C(\tilde{\Z}_2) \), with \( \mathbf{s} \) and
\( \mathbf{t} \) unit vectors. That correlation is \( \mathbf{s}\T\Q_1\T\Q_2\mathbf{t}\le\rho_1 \) (@prp-mat-svd-norms; this is the
first part of @exr-col-canonical-max). Hence
\( \text{GVIF}_1\ge1/(1-\rho_1^2)\ge1/(1-R_{j\mid2}^2) \) by @prp-col-gvif(b). The ordinary \( \text{VIF}_j \) regresses \( \x_j \) also on the
other columns of its own block, which can explain it well even when \( \Z_2 \) cannot. The indicators of a
factor with a small reference level are an example: each is nearly determined by the others.
:::

### C. Going deeper

::: {#exr-col-canonical-max}
[C1]

Show that \( \rho_1=\max\{\cos\angle(\bu,\bv):\bu\in\C(\tilde{\Z}_1),\ \bv\in\C(\tilde{\Z}_2),\ \bu,\bv\ne\bzero\} \). More generally,
show that the generalized eigenvalues of \( \Cov(\hbeta_1) \) relative to \( \Cov_\perp \), the roots \( \theta \) of
\( \det(\Cov(\hbeta_1)-\theta\Cov_\perp)=0 \), are \( 1/(1-\rho_i^2) \), \( i\le m \), together with \( q-m \) ones. Deduce that for every
\( \mathbf{a}\ne\bzero \) the ratio \( \Var(\mathbf{a}\T\hbeta_1)/\mathbf{a}\T\Cov_\perp\mathbf{a} \) lies between \( 1 \) and \( 1/(1-\rho_1^2) \).
:::

::: {.solution}
Write \( \bu=\Q_1\mathbf{s} \) and \( \bv=\Q_2\mathbf{t} \) with unit \( \mathbf{s},\mathbf{t} \). Then
\( \cos\angle(\bu,\bv)=\mathbf{s}\T\Q_1\T\Q_2\mathbf{t} \), whose maximum over unit vectors is the largest singular value of
\( \Q_1\T\Q_2 \) (@prp-mat-svd-norms). For the second part, take \( \tilde{\Z}_1=\Q_1 \) by @prp-col-gvif(c), which changes both
matrices by the same congruence and so leaves the generalized eigenvalues unchanged.
Then \( \Cov_\perp=\sigma^2\I \) and \( \Cov(\hbeta_1)=\sigma^2(\I-\mathbf{C}\mathbf{C}\T)^{-1} \), whose eigenvalues are the stated ones. The
ratio in the last part is a ratio of quadratic forms, bounded by @cor-mat-generalized-rayleigh, so it lies
between the smallest and the largest generalized eigenvalue.
:::
