# Leverage as geometry

The projection \( \M \) is an \( n\times n \) matrix whose entries have a direct
interpretation. Since \( \hat y_i=\sum_j m_{ij}y_j \), the entry \( m_{ij} \) is the weight
observation \( j \) receives in the fitted value for observation \( i \). So
\[
m_{ij}=\frac{\partial\hat y_i}{\partial y_j} .
\]
The diagonal entries measure how strongly each observation pulls the fitted
surface towards itself.

::: {#def-proj-leverage}
[Leverage]

The **leverage** of observation \( i \) is \( h_{ii}=m_{ii} \), the \( i \)th diagonal entry
of \( \M \). In this context \( \M \) is often called the **hat matrix** and written \( \bH \),
because it “puts the hat on” \( \y \).
:::

Leverage is a property of the design, not of the response. It depends on \( \y \) not at
all, and on \( \X \) only through \( \C(\X) \), so it is unchanged by reparameterization.

::: {#prp-proj-leverage}
[Properties of leverage]

Let \( \rank(\X)=r \), and let \( \vect{e}_i \) denote the \( i \)th standard basis vector. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( h_{ii}=\norm{\M\vect{e}_i}^2=\sum_j m_{ij}^2 \), the squared length of the projection
           of \( \vect{e}_i \) onto \( \C(\X) \);

2. \( 0\le h_{ii}\le1 \), and \( \sum_i h_{ii}=r \), so the average leverage is \( r/n \);

3. if \( \bone\in\C(\X) \), then \( h_{ii}\ge1/n \);

4. \( h_{ii}=1 \) iff \( \vect{e}_i\in\C(\X) \). In that case \( \hat y_i=y_i \) for every \( \y \), the \( i \)th
           residual is always zero, and \( m_{ij}=0 \) for \( j\ne i \);

5. if \( \X=\Q\R \) with \( \Q \) having orthonormal columns spanning \( \C(\X) \), then
           \( h_{ii} \) is the squared length of the \( i \)th row of \( \Q \).
:::

:::

::: {.proof}
(a) \( h_{ii}=\vect{e}_i\T\M\vect{e}_i=\vect{e}_i\T\M\T\M\vect{e}_i \), and \( \M\vect{e}_i \) is the \( i \)th column of \( \M \).
(b) By @prp-proj-trace-rank(d), \( 0\le\norm{\M\vect{e}_i}^2\le\norm{\vect{e}_i}^2=1 \). The sum is
\( \tr(\M)=r \).
(c) With \( \Mo=n^{-1}\bone\bone\T \), the difference \( \M-\Mo \) is a projection
(@thm-proj-nested), hence nonnegative definite, so
\( h_{ii}\ge(\Mo)_{ii}=1/n \).
(d) Equality holds in \( \norm{\M\vect{e}_i}\le\norm{\vect{e}_i} \) iff \( (\I-\M)\vect{e}_i=\bzero \), that is, iff
\( \vect{e}_i\in\C(\X) \). Then \( \M\vect{e}_i=\vect{e}_i \), so column \( i \) of \( \M \) is \( \vect{e}_i \), and by symmetry so
is row \( i \). Hence \( \hat y_i=y_i \).
(e) \( \M=\Q\Q\T \), so \( m_{ii}=\sum_k q_{ik}^2 \).
:::

Part (d) describes the extreme case. An observation with leverage one has a
coordinate direction of its own in the model space. This happens, for example, when it is
the only member of a group in a one-way layout, or when a column of \( \X \) is an indicator
for that single observation. The model then fits the observation exactly, whatever its
value, and the data give no check on it. High leverage is a weaker form of the same
thing. If \( h_{ii} \) is close to one, then \( \hat y_i\approx y_i \), and a
gross error in \( y_i \) barely shows up in its residual.

This shows up in the residual variances. If \( \Cov(\be)=\sigma^2\I \), then by
@prp-proj-fit-residual(e)
\[
\Var(\hat{\varepsilon}_i)=\sigma^2(1-h_{ii}),
\]
so residuals at high-leverage points are *less* variable than others, even
though those are exactly the points where an error would do the most damage.
Chapter 20 standardizes residuals to correct for this.

## Leverage as distance from the centre

When the model has an intercept, leverage has a second description in variable
space: distance from the centroid of the regressors. Write
\( \X=[\bone,\X_1] \), let \( \bar{\x} \) be the vector of column means of \( \X_1 \), and let
\( \tilde{\X}_1=(\I-n^{-1}\bone\bone\T)\X_1 \) be the centred regressors, assumed to have full column rank.

::: {#prp-proj-leverage-mahalanobis}
Under these assumptions,
\[
h_{ii}=\frac1n+(\x_{(i)}-\bar{\x})\T\bigl(\tilde{\X}_1\T\tilde{\X}_1\bigr)^{-1}(\x_{(i)}-\bar{\x}),
\]
where \( \x_{(i)} \) is the \( i \)th row of \( \X_1 \), written as a column vector.
:::

::: {.proof}
By @lem-proj-fwl-split with \( \X_2=\bone \),
\( \M=n^{-1}\bone\bone\T+\tilde{\X}_1(\tilde{\X}_1\T\tilde{\X}_1)^{-1}\tilde{\X}_1\T \). The \( i \)th row
of \( \tilde{\X}_1 \) is \( (\x_{(i)}-\bar{\x})\T \).
:::

The quadratic form is \( (n-1)^{-1} \) times the squared Mahalanobis distance of \( \x_{(i)} \)
from the centroid, measured in the metric of the sample covariance matrix of the
regressors. Leverage is therefore high for points that are unusual *jointly*,
even when no single coordinate is extreme. A state with moderate poverty and
moderate single parenthood, but an unusual *combination* of the two, can have
high leverage. For simple regression the formula reduces to
\( h_{ii}=n^{-1}+(x_i-\bar x)^2/S_{xx} \).

::: {#exm-proj-leverage}
[High-leverage states]

Return to the state data, now with all \( 51 \) jurisdictions including the
District of Columbia, and the model of @exm-proj-fwl-crime. The average
leverage is \( p/n=0.078 \), and the smallest is
\( 0.021 \), just above the bound \( 1/n \). The three largest are

| Jurisdiction | \( h_{ii} \) | multiple of average |
|---|---|---|
| District of Columbia | 0.506 | 6.5 |
| Mississippi | 0.212 | 2.7 |
| Alaska | 0.177 | 2.3 |

The District of Columbia is entirely urban, with by far the highest proportion of
single-parent households. It sits at a corner of the regressor cloud that no state comes
near, and its leverage is about half the maximum possible. It was omitted from @exm-proj-fwl-crime.
@exr-proj-dc asks what happens to the poverty
coefficient when it is put back.
:::

```{.python .run #cell-leverage-leverage}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data          # District of Columbia kept here
X = np.column_stack([np.ones(len(data)), data["poverty"], data["single"], data["urban"]])

Q, _ = np.linalg.qr(X)                 # orthonormal basis for C(X)
h = np.sum(Q ** 2, axis=1)             # h_ii = m_ii = ||row i of Q||^2
n, p = X.shape
print("sum of leverages:", h.sum(), " (= p =", p, ")")
print("range:", h.min(), h.max(), " bounds:", 1 / n, 1)
top = np.argsort(h)[::-1][:3]
for i in top:
    print(f"{data.index[i]:22s} h = {h[i]:.3f}   ({h[i] / (p / n):.1f} x average)")
```

The listing computes leverages from a QR factorization rather than by forming \( \M \).
For large \( n \) the \( n\times n \) matrix \( \M \) should never be formed.
[Section 6.10](10-computation.html) explains why.

## Exercises

### A. Check your understanding

::: {#exr-proj-simple-leverage}
[A1]

For simple regression, derive \( h_{ii}=n^{-1}+(x_i-\bar x)^2/S_{xx} \) directly from
@exm-proj-simple-nested, and show that the largest leverage belongs to the
observation with \( x_i \) farthest from \( \bar x \).
:::

### B. Practice

::: {#exr-proj-offdiag}
[B1]

Show that \( m_{ij}^2\le m_{ii}(1-m_{ii}) \) for \( i\neq j \), and hence \( \lvert m_{ij}\rvert\le1/2 \).
When is the bound attained?
:::

::: {.solution}
From @prp-proj-leverage(a),
\( m_{ii}=\sum_k m_{ik}^2\ge m_{ii}^2+m_{ij}^2 \). So \( m_{ij}^2\le m_{ii}(1-m_{ii})\le1/4 \).
Equality needs \( m_{ii}=1/2 \) and \( m_{ij} \) the only other nonzero entry of row \( i \). The
intercept-only model with \( n=2 \), where every entry of \( \M \) is \( 1/2 \), attains it.
:::

::: {#exr-proj-leverage-monotone}
[B2]

Show that adding a column to \( \X \) cannot decrease any leverage. Show that the leverages in
a one-way layout are \( 1/n_k \) for members of group \( k \), and explain why an observation
forming a group on its own has leverage one.
:::
