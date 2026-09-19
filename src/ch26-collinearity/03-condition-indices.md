# Condition indices and variance decomposition

A variance inflation factor says how much one coefficient has suffered, but not which other
regressors are responsible, how many near-dependencies there are, or whether the intercept is
involved. The eigenvalue analysis of [Section 26.1](01-variance.html) answers those questions but depends on
the units of the columns. This section fixes a scaling and derives the variance-decomposition
proportions of Belsley, Kuh and Welsch (1980).

## Scaling the columns

The canonical directions of \( \X\T\X \) change when a column is rescaled (@exr-col-units-direction), and a
merely short column produces a small eigenvalue without any near-dependence (@exr-col-max-coefficient).
The diagnostics are therefore computed for the design with every column, including \( \bone \), scaled to unit
length:
\[
\tilde{\X}=\X\bD_0,\qquad \bD_0=\diag\bigl(1/\norm{\x_1},\dots,1/\norm{\x_p}\bigr).
\]
By @prp-cmp-van-der-sluis this is within a factor \( \sqrt p \) of the best diagonal scaling for the condition
number. The coefficients of \( \tilde{\X} \) are \( \tilde{\beta}_j=\norm{\x_j}\beta_j \).

::: {#def-col-condition-index}
[Condition indices]

Let \( \tilde{\X} \) have singular values \( \mu_1\ge\mu_2\ge\dots\ge\mu_p>0 \). The **condition indices** of \( \X \) are
\[
\eta_k=\frac{\mu_1}{\mu_k},\qquad k=1,\dots,p .
\]
The largest, \( \eta_p=\kappa(\tilde{\X}) \), is the **scaled condition number**.
:::

Whether to centre before scaling is a real choice. Let \( \R \) be the correlation matrix of the non-intercept columns, the cross-product matrix of
the centred and scaled design, and write \( s_j^2=S_{jj}/n \) and \( \text{CV}_j=s_j/\lvert\bar x_j\rvert \) for the
**coefficient of variation** of a column with nonzero mean.

::: {#prp-col-scaling}
[Scaled condition numbers]

Let \( \X \) contain the column \( \bone \) and have full column rank, and let \( \R \) be nonsingular.

::: {.enumerate options="label=(\alph*)"}
1. \( \sum_k\mu_k^2=p \), so \( 1\le\mu_1^2\le p \). Every diagonal entry of \( (\tilde{\X}\T\tilde{\X})^{-1} \) is at most
   \( 1/\mu_p^2\le\eta_p^2 \).

2. \( \lambda_{\min}(\R)\ge\mu_p^2 \). Hence \( \eta_p^2\ge1/\lambda_{\min}(\R)\ge\max_j\text{VIF}_j \): every near-dependency
   among the centred regressors is at least as severe in the uncentred scaled design.

3. If a non-intercept column \( \x_j \) has \( \bar x_j\ne0 \), then \( \mu_p^2\le\text{CV}_j^2/(\text{CV}_j^2+2) \), so that
   \( \eta_p\ge\sqrt{(\text{CV}_j^2+2)/\text{CV}_j^2}>\sqrt2/\text{CV}_j \).
:::

:::

::: {.proof}
(a) \( \sum_k\mu_k^2=\norm{\tilde{\X}}_F^2=p \) (@prp-mat-svd-norms), since each column has unit length. The largest of
\( p \) numbers with sum \( p \) is between \( 1 \) and \( p \). The \( j \)th diagonal entry of \( (\tilde{\X}\T\tilde{\X})^{-1} \) is a
Rayleigh quotient of a matrix whose largest eigenvalue is \( 1/\mu_p^2 \) (@thm-mat-extremal-rayleigh), and
\( 1/\mu_p^2=\eta_p^2/\mu_1^2\le\eta_p^2 \).

(b) Number the columns so that \( \x_0=\bone \) and \( \x_1,\dots,\x_k \) are the others. The centred and scaled
columns are \( \bz_j=(\x_j-\bar x_j\bone)/(\sqrt n\,s_j) \), so that \( \R=\Z\T\Z \) with \( \Z=[\bz_1,\dots,\bz_k] \). Since
\( \bone=\sqrt n\,\tilde{\x}_0 \) and \( \x_j=\norm{\x_j}\tilde{\x}_j \), where \( \tilde{\x}_j \) are the columns of \( \tilde{\X} \),
\[
\Z\mathbf{c}=\sum_{j=1}^kc_j\frac{\norm{\x_j}\tilde{\x}_j-\sqrt n\,\bar x_j\tilde{\x}_0}{\sqrt n\,s_j}=\tilde{\X}\bb,
\]
where
\[
b_j=\frac{c_j\norm{\x_j}}{\sqrt n\,s_j}\quad(j\ge1),\qquad b_0=-\sum_{j=1}^k\frac{c_j\bar x_j}{s_j}.
\]
Because \( \norm{\x_j}^2=n(s_j^2+\bar x_j^2)\ge ns_j^2 \), \( \lvert b_j\rvert\ge\lvert c_j\rvert \), so \( \norm{\bb}\ge\norm{\mathbf{c}} \). For
a unit vector \( \mathbf{c} \), \( \mathbf{c}\T\R\mathbf{c}=\norm{\tilde{\X}\bb}^2\ge\mu_p^2\norm{\bb}^2\ge\mu_p^2 \) by
@thm-mat-extremal-rayleigh; minimizing over \( \mathbf{c} \) gives \( \lambda_{\min}(\R)\ge\mu_p^2 \). Then
\( \eta_p^2\ge1/\mu_p^2\ge1/\lambda_{\min}(\R) \) by (a), and \( 1/\lambda_{\min}(\R)\ge\text{VIF}_j \) by @prp-dep-collinear-directions(c).

(c) Take \( \bb \) with \( b_j=\norm{\x_j} \), \( b_0=-\sqrt n\,\bar x_j \) and all other entries zero. Then
\( \tilde{\X}\bb=\x_j-\bar x_j\bone \), with squared length \( ns_j^2 \), while
\( \norm{\bb}^2=n(s_j^2+\bar x_j^2)+n\bar x_j^2 \). So
\( \mu_p^2\le ns_j^2/\bigl(n(s_j^2+2\bar x_j^2)\bigr)=\text{CV}_j^2/(\text{CV}_j^2+2) \), and \( \eta_p\ge1/\mu_p \) by (a).
:::

Part (b) is the formal content of the argument of Belsley (1984) for not centring: centring can only
hide near-dependencies, never create them. Part (c) shows what it hides. A regressor that varies
little relative to its mean is nearly a multiple of the column of ones, and that near-dependency with
the intercept is invisible in \( \R \) and in every variance inflation factor. The design of
@exm-mat-macro-svd shows the size of the effect: scaled to unit length, its condition number is \( 383 \);
centred and scaled, \( 124 \).

Whether the hidden dependency matters depends on the question, as in @exm-col-quadratic. It inflates
the variance of the intercept, the mean response at the origin, and of the coefficients whose columns
take part in it. If the origin is far from the data, it concerns a question nobody asked; if the origin is
meaningful, the uncentred diagnostic is the honest one. This chapter uses the uncentred scaled design, as
Belsley, Kuh and Welsch do, and reports the centred condition number alongside.

## The variance-decomposition proportions

Let the singular value decomposition of the scaled design be
\[
\tilde{\X}=\sum_{k=1}^p\mu_k\bu_k\mathbf{t}_k\T,
\]
with orthonormal right singular vectors \( \mathbf{t}_k=(t_{1k},\dots,t_{pk})\T \). The \( p\times p \) matrix with
columns \( \mathbf{t}_k \) is orthogonal, so its rows also have unit length: \( \sum_kt_{jk}^2=1 \) for every \( j \).

::: {#thm-col-decomposition}
[Variance-decomposition proportions]

Let \( \X \) have rank \( p \), \( \Cov(\Y)=\sigma^2\I \), and let \( \phi_{jk}=t_{jk}^2/\mu_k^2 \).

::: {.enumerate options="label=(\alph*)"}
1. \( \Var(\hat{\beta}_j)=\dfrac{\sigma^2}{\norm{\x_j}^2}\sum_{k=1}^p\phi_{jk} \). The **variance-decomposition
   proportions**
   \[
   \pi_{kj}=\frac{\phi_{jk}}{\sum_{l=1}^p\phi_{jl}}
   \]
   satisfy \( \sum_k\pi_{kj}=1 \), do not depend on \( \sigma^2 \), and vanish iff \( t_{jk}=0 \). The term \( \pi_{kj} \) is the share of
   \( \Var(\hat{\beta}_j) \) contributed by the canonical coordinate with condition index \( \eta_k \).

2. For each \( k \), \( \norm{\sum_jt_{jk}\tilde{\x}_j}=\mu_k \) while \( \norm{\mathbf{t}_k}=1 \): the vector \( \mathbf{t}_k \) holds the coefficients of a
   near-dependency among the scaled columns, of strength \( \mu_k \).

3. Let \( S \) be a nonempty proper subset of \( \{1,\dots,p\} \) such that \( \mu_k\le\mu_S \) for \( k\in S \) and \( \mu_k\ge\mu_L \)
   for \( k\notin S \), with \( \mu_S<\mu_L \). Let \( w_j=\sum_{k\in S}t_{jk}^2 \) measure the participation of column \( j \) in the
   near-dependencies indexed by \( S \). Then
   \[
   \sum_{k\in S}\pi_{kj}\ \ge\ \frac{w_j}{w_j+(1-w_j)\,\mu_S^2/\mu_L^2},
   \qquad
   \Var(\hat{\beta}_j)\ \ge\ \frac{\sigma^2w_j}{\norm{\x_j}^2\mu_S^2}.
   \]
:::

:::

::: {.proof}
(a) Since \( \X\bbeta=\tilde{\X}\bD_0^{-1}\bbeta \), the least squares coefficients of \( \tilde{\X} \) are
\( \hat{\tilde{\beta}}_j=\norm{\x_j}\hat{\beta}_j \), and by @thm-col-variance(a) applied to \( \tilde{\X} \), whose cross-product matrix
has eigenvalues \( \mu_k^2 \) and eigenvectors \( \mathbf{t}_k \), \( \Cov(\hat{\tilde{\bbeta}})=\sigma^2\sum_k\mathbf{t}_k\mathbf{t}_k\T/\mu_k^2 \). Its
\( j \)th diagonal entry is \( \sigma^2\sum_k\phi_{jk} \). The statements about the proportions follow from the definition, since
\( \mu_k>0 \).

(b) \( \sum_jt_{jk}\tilde{\x}_j=\tilde{\X}\mathbf{t}_k=\mu_k\bu_k \).

(c) Put \( A=\sum_{k\in S}\phi_{jk} \) and \( B=\sum_{k\notin S}\phi_{jk} \). Then \( A\ge w_j/\mu_S^2 \) and, since
\( \sum_{k\notin S}t_{jk}^2=1-w_j \), \( B\le(1-w_j)/\mu_L^2 \). The sum of the proportions over \( S \) is \( A/(A+B) \), which
increases with \( A \) and decreases with \( B \), so it is at least
\( (w_j/\mu_S^2)/\bigl(w_j/\mu_S^2+(1-w_j)/\mu_L^2\bigr) \), which is the stated bound. The variance bound is
\( \Var(\hat{\beta}_j)\ge\sigma^2A/\norm{\x_j}^2 \) together with \( A\ge w_j/\mu_S^2 \).
:::

Part (c) is why the proportions locate near-dependencies. If the condition indices in \( S \) are at least
ten times those outside it, \( \mu_S^2/\mu_L^2\le0.01 \), then any column with even a tenth of its unit row inside
\( S \), \( w_j\ge0.1 \), has at least \( 0.1/(0.1+0.9\times0.01)\approx0.92 \) of its variance attributed to \( S \), and its
variance is inflated accordingly. A column that does not take part has proportions near zero there.
When several small singular values are close together, the individual \( \mathbf{t}_k \)
are poorly determined, and a column's proportions can be split between their rows in unpredictable
ratios. Part (c) bounds the *sum* over the group, which is stable, and that is what should be read (@exr-col-tied-proportions).

::: {.idea}
Belsley, Kuh and Welsch's reading of the table: scale the columns to unit length; flag condition
indices of about \( 30 \) or more (indices near \( 10 \) signal weak dependencies); in each flagged row, or
group of rows with similar indices, look for two or more coefficients with large proportions,
conventionally above \( 0.5 \). Those columns form a near-dependency. Confirm it by regressing one of them
on the others.
:::

::: {#exm-col-macro-bkw}
[A consumption function]

Return to the regression of @exm-col-macro-vif, with an intercept and six regressors. The scaled
design has condition indices and proportions

| Index | const. | income | pop. | CPI | M1 | T-bill | unemp. |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| \( 1.0 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.001 \) | \( 0.001 \) |
| \( 3.4 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.001 \) | \( 0.053 \) | \( 0.007 \) |
| \( 7.7 \) | \( 0.000 \) | \( 0.000 \) | \( 0.000 \) | \( 0.001 \) | \( 0.001 \) | \( 0.347 \) | \( 0.090 \) |
| \( 13.4 \) | \( 0.002 \) | \( 0.001 \) | \( 0.000 \) | \( 0.001 \) | \( 0.002 \) | \( 0.001 \) | \( 0.688 \) |
| \( 38.9 \) | \( 0.005 \) | \( 0.149 \) | \( 0.000 \) | \( 0.001 \) | \( 0.290 \) | \( 0.049 \) | \( 0.036 \) |
| \( 74.3 \) | \( 0.007 \) | \( 0.130 \) | \( 0.000 \) | \( 0.995 \) | \( 0.569 \) | \( 0.241 \) | \( 0.116 \) |
| \( 270.9 \) | \( 0.985 \) | \( 0.720 \) | \( 0.999 \) | \( 0.001 \) | \( 0.138 \) | \( 0.307 \) | \( 0.060 \) |

[Figure 26.3.1](#fig-col-proportions) shows the table graphically. Three indices exceed \( 30 \), and the table
separates them.

- *Index \( 270.9 \): intercept, income and population* (proportions \( 0.985 \), \( 0.720 \), \( 0.999 \)). The singular
  vector has entries \( -0.522 \), \( -0.236 \) and \( 0.817 \) on these columns and small ones elsewhere: population,
  which grew steadily with income, is nearly an affine function of it (\( R^2=0.9885 \)). Its small coefficient of
  variation, \( 0.156 \), alone forces an index of at least \( 9.1 \) by @prp-col-scaling(c); the shared trend does the rest.
- *Index \( 74.3 \): the price level and money.* CPI has proportion \( 0.995 \) and M1 \( 0.569 \). Regressing CPI on income
  and M1 gives \( R^2=0.9868 \).
- *Index \( 38.9 \): a weaker dependency involving M1 and income.* Taken with the previous row, M1 has
  \( 0.290+0.569=0.859 \) of its variance in the two rows, and income has \( 0.149+0.130+0.720=0.999 \) in the
  three rows above \( 30 \). Income takes part in all three near-dependencies.

A proportion is a share, not a size: the \( 0.307 \) of the T-bill rate at index \( 270.9 \) is a share of a
variance that collinearity has barely inflated (its variance inflation factor is \( 2.9 \)). Centred and scaled,
the design has condition number \( 37.9 \) against \( 270.9 \) uncentred: centring removes the intercept from the
largest near-dependency, and what is left is the relation between income and population that the
variance inflation factors \( 154.5 \) and \( 183.7 \) already showed.
:::

::: {when-format="html"}
![**Figure 26.3.1.** Variance-decomposition proportions of @exm-col-macro-bkw: rows are condition
indices, columns coefficients, dot area proportional to \( \pi_{kj} \). Highlighted: \( \pi_{kj}\ge0.5 \) with index at least \( 30 \).](variance_proportions.svg){#fig-col-proportions width=85%}
:::

::: {when-format="pdf"}
![Variance-decomposition proportions of @exm-col-macro-bkw: rows are condition
indices, columns coefficients, dot area proportional to \( \pi_{kj} \). Highlighted: \( \pi_{kj}\ge0.5 \) with index at least \( 30 \).](variance_proportions.pdf){width=85%}
:::

```{.python .run #cell-condition-indices-bkw}
import numpy as np
import statsmodels.api as sm
macro = sm.datasets.macrodata.load_pandas().data
names = ["const", "realdpi", "pop", "cpi", "m1", "tbilrate", "unemp"]
X = sm.add_constant(macro[names[1:]].to_numpy())
X_s = X / np.linalg.norm(X, axis=0)                   # unit-length columns, not centred

U, mu, Vt = np.linalg.svd(X_s, full_matrices=False)   # X_s = U diag(mu) V'
eta = mu[0] / mu                                      # condition indices
phi = Vt.T ** 2 / mu**2                               # phi[j, k] = v_jk^2 / mu_k^2
prop = (phi / phi.sum(axis=1, keepdims=True)).T       # row k: proportions for index eta_k

print("index   " + " ".join(f"{nm[:7]:>7s}" for nm in names))
for k in range(len(mu)):
    print(f"{eta[k]:7.1f} " + " ".join(f"{p:7.3f}" for p in prop[k]))
```

```{.python .run #cell-condition-indices-centred}
Z = X[:, 1:] - X[:, 1:].mean(axis=0)
Z_s = Z / np.linalg.norm(Z, axis=0)                   # centred, then unit length
mu_c = np.linalg.svd(Z_s, compute_uv=False)
print("scaled condition number, not centred:", round(eta[-1], 1))
print("scaled condition number, centred    :", round(mu_c[0] / mu_c[-1], 1))
cv = X[:, 1:].std(axis=0) / X[:, 1:].mean(axis=0)
print("coefficients of variation:", dict(zip(names[1:], np.round(cv, 3))))
```

```{.python .run #cell-condition-indices-auxiliary}
def r_squared(target, others):
    fit = sm.OLS(macro[target], sm.add_constant(macro[others])).fit()
    return fit.rsquared

print("pop on realdpi     : R^2 =", round(r_squared("pop", ["realdpi"]), 4))
print("cpi on realdpi, m1 : R^2 =", round(r_squared("cpi", ["realdpi", "m1"]), 4))
```

The script behind the example checks the bound of @thm-col-decomposition(c) for the groups formed by the
last one, two and three indices, and the inequalities of @prp-col-scaling.

## Single cases

Collinearity diagnostics can be dominated by a few cases: one far out along a near-dependency can
create it, one lying across it can hide it. The deletion formulas of
[Chapter 20](../ch20-residuals-leverage-influence/index.html) show that only cases of high leverage can
have much effect.

::: {#prp-col-deletion}
[Deleting a case]

Let \( \X \) have rank \( p \), let \( \x_{(i)} \) be its \( i \)th row with leverage \( h_{ii}<1 \), let \( \X_{[i]} \) be \( \X \) without that
row (@thm-res-deletion writes \( \X_{(i)} \); square brackets here avoid a clash with \( \X_{(j)} \), the matrix without
column \( j \)), and write \( \G=(\X\T\X)^{-1} \), \( \G_{[i]}=(\X_{[i]}\T\X_{[i]})^{-1} \).

::: {.enumerate options="label=(\alph*)"}
1. \( \mathbf{a}\T\G_{[i]}\mathbf{a}=\mathbf{a}\T\G\mathbf{a}+(\mathbf{a}\T\G\x_{(i)})^2/(1-h_{ii}) \) for every \( \mathbf{a} \).

2. \( \mathbf{a}\T\G\mathbf{a}\le\mathbf{a}\T\G_{[i]}\mathbf{a}\le\mathbf{a}\T\G\mathbf{a}/(1-h_{ii}) \): deleting case \( i \) multiplies the variance of every
   linear function of \( \hbeta \) by a factor between \( 1 \) and \( 1/(1-h_{ii}) \).

3. \( (1-h_{ii})\lambda_{\min}(\X\T\X)\le\lambda_{\min}(\X_{[i]}\T\X_{[i]})\le\lambda_{\min}(\X\T\X) \).
:::

:::

::: {.proof}
(a) This is the rank-one deletion formula \( \G_{[i]}=\G+\G\x_{(i)}\x_{(i)}\T\G/(1-h_{ii}) \) of @exm-mat-deletion, on which
@thm-res-deletion rests, sandwiched between \( \mathbf{a}\T \) and \( \mathbf{a} \).
(b) The lower bound is clear from (a). For the upper, apply the Cauchy–Schwarz inequality to the inner
product \( \langle\mathbf{a},\bb\rangle=\mathbf{a}\T\G\bb \), which is legitimate because \( \G \) is positive definite:
\( (\mathbf{a}\T\G\x_{(i)})^2\le(\mathbf{a}\T\G\mathbf{a})\,h_{ii} \). Then
\( \mathbf{a}\T\G_{[i]}\mathbf{a}\le\mathbf{a}\T\G\mathbf{a}\,(1+h_{ii}/(1-h_{ii})) \).
(c) For a positive definite matrix, \( 1/\lambda_{\min} \) is the largest eigenvalue of the inverse, which is
\( \max_{\norm{\mathbf{a}}=1}\mathbf{a}\T(\cdot)^{-1}\mathbf{a} \) (@thm-mat-extremal-rayleigh). Taking the maximum in (b) over unit vectors gives
the two inequalities.
:::

So removing a case of small leverage changes the smallest eigenvalue by at most the factor
\( 1-h_{ii} \), and no variance by more than \( 1/(1-h_{ii}) \). In the consumption regression
the largest leverage is \( 0.161 \), at 2009Q3, so by (c) no single quarter can reduce the smallest
eigenvalue of the scaled cross-product matrix by more than the factor \( 1.192 \), with the column scaling of the full data held fixed. The largest
reduction actually achieved by deleting one of the \( 203 \) quarters is by the factor \( 1.055 \). The near-dependencies of
this design are carried by the whole sample, not by a few quarters.

Adding a case never increases a variance, by (b) read the other way round, but it can increase every
variance inflation factor and the condition number (@exr-col-add-point).

## Exercises

### A. Check your understanding

::: {#exr-col-read-table}
[A1]

In the table of @exm-col-macro-bkw, why is the intercept's variance almost entirely attributed to the
largest index, although the intercept is not an economic variable? Which coefficients would you expect
to change most if a few more years of data were added in which population grew faster than income?
:::

### B. Practice

::: {#exr-col-cv-two}
[B1]

Let \( \X=[\bone,\x] \) with \( \bar x\ne0 \). Show that the scaled design has
\( \mu_{1,2}^2=1\pm\cos\theta \), where \( \cos\theta=1/\sqrt{1+\text{CV}^2} \) is the cosine of the angle between \( \bone \) and \( \x \),
and hence \( \eta_2^2=(1+\cos\theta)/(1-\cos\theta) \). Compare with the bound of @prp-col-scaling(c). What is
the variance inflation factor of \( \x \)?
:::

::: {.solution}
The scaled columns are unit vectors at angle \( \theta \), so \( \tilde{\X}\T\tilde{\X}=\begin{psmallmatrix}1&c\\c&1\end{psmallmatrix} \) with
\( c=\cos\theta \), whose eigenvalues are \( 1\pm c \) (@exr-mat-cond-two). Now
\( c=\bone\T\x/(\sqrt n\norm{\x})=n\bar x/\sqrt{n\cdot n(s^2+\bar x^2)}=1/\sqrt{1+\text{CV}^2} \) (taking \( \bar x>0 \)). So
\( \eta_2^2=(1+c)/(1-c) \). The bound of the proposition gives \( \mu_2^2\le\text{CV}^2/(\text{CV}^2+2) \). For small CV,
\( 1-c\approx\text{CV}^2/2 \), so \( \mu_2^2=1-c\approx\text{CV}^2/2 \), and the bound is attained to first order. With a single regressor
\( R_j^2=0 \) and \( \text{VIF}=1 \): the variance inflation factor cannot see this near-dependency at all.
:::

::: {#exr-col-add-point}
[B2]

Take \( n=20 \) cases with \( (x_1,x_2) \) at the points of a \( 4\times5 \) grid in \( [-1,1]^2 \), so that the regressors are
uncorrelated, and add a twenty-first case at \( (t,t) \). Show that the correlation of the two regressors tends to
one as \( t\to\infty \), so that both variance inflation factors grow without bound, while no variance of any
linear function of \( \hbeta \) increases. Reconcile the two statements.
:::

::: {.solution}
By @prp-col-deletion(b) read backwards, adding a case can only decrease \( \mathbf{a}\T\G\mathbf{a} \) for every \( \mathbf{a} \). The sample
covariance of \( x_1,x_2 \) in the enlarged data is \( t^2(1-1/21)/21 \) (the grid contributes nothing, being centred and
uncorrelated), and the sample variance of \( x_j \) is \( (S_j+t^2(1-1/21))/21 \), where \( S_j \) is the grid's sum of
squares for \( x_j \). As \( t\to\infty \)
the correlation tends to one and \( \text{VIF}_j=1/(1-r^2)\to\infty \). But \( \Var(\hat{\beta}_j)=\sigma^2\text{VIF}_j/S_{jj} \), and
\( S_{jj} \) grows like \( t^2 \) as fast as \( \text{VIF}_j \) does, so the variances converge to finite limits no larger than
before. The variance inflation factor compares with an orthogonal design with the *new* spread, which the
added point has made enormous; the point added information along \( (1,1) \) and none across it.
:::

### C. Going deeper

::: {#exr-col-tied-proportions}
[C1]

Suppose \( \mu_k=\mu \) for all \( k\in S \). Show that the singular vectors \( \mathbf{t}_k \), \( k\in S \), are determined only up to an
orthogonal transformation of their span, that the individual proportions \( \pi_{kj} \), \( k\in S \), depend on the choice,
and that \( \sum_{k\in S}\pi_{kj} \) does not.
:::

::: {.solution}
If \( \mathbf{T}_S \) has columns \( \mathbf{t}_k \), \( k\in S \), and \( \mathbf{O} \) is orthogonal, then
\( \tilde{\X}\T\tilde{\X}\,\mathbf{T}_S\mathbf{O}=\mu^2\mathbf{T}_S\mathbf{O} \), so the columns of \( \mathbf{T}_S\mathbf{O} \) are equally valid singular vectors, and
the individual \( t_{jk}^2 \) change with \( \mathbf{O} \). But
\( \sum_{k\in S}\phi_{jk}=\mu^{-2}\sum_{k\in S}t_{jk}^2=\mu^{-2}\,\mathbf{e}_j\T\mathbf{T}_S\mathbf{T}_S\T\mathbf{e}_j \), and \( \mathbf{T}_S\mathbf{T}_S\T \) is the projection
onto the span, which does not depend on \( \mathbf{O} \). The denominator \( \sum_l\phi_{jl} \) is a diagonal entry of
\( (\tilde{\X}\T\tilde{\X})^{-1} \) and does not depend on the choice either.
:::
