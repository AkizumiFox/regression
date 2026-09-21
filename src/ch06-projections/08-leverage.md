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
[Chapter 20](../ch20-residuals-leverage-influence/index.html) standardizes residuals to correct for this (@def-res-residuals).

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

### C. Going deeper

::: {#exr-proj-leverage-one-rank}
[C1]

@prp-proj-leverage(d) says that \( h_{ii}=1 \) iff \( \vect{e}_i\in\C(\X) \). Here is what that means
for the design. Write \( \X_{(i)} \) for \( \X \) with row \( i \) deleted, and let
\( \bD \) be the \( (n-1)\times n \) matrix that deletes the \( i \)th coordinate of a vector, so that
\( \C(\X_{(i)})=\bD\,\C(\X) \) and \( \Null(\bD)=\spn(\vect{e}_i) \).

::: {.enumerate options="label=(\alph*)"}
1. Show that \( \rank(\X_{(i)})=\rank(\X)-\dim\bigl(\C(\X)\cap\spn(\vect{e}_i)\bigr) \), and deduce
   that \( h_{ii}=1 \) iff deleting case \( i \) lowers the rank of the design.

2. Deduce that at most \( \rank(\X) \) cases can have leverage one, and that if \( k \) of them do,
   the remaining leverages sum to \( \rank(\X)-k \).

3. Show that \( \bD\M\bD\T \) is the orthogonal projection onto \( \C(\X_{(i)}) \) when
   \( h_{ii}=1 \), and conclude that deleting such a case changes none of the other fitted
   values. In what sense does the data set then contain no information about \( y_i \)?
:::
:::

::: {.solution}
(a) The map \( \bD \) restricted to \( \C(\X) \) is linear with image \( \C(\X_{(i)}) \) and kernel
\( \C(\X)\cap\Null(\bD)=\C(\X)\cap\spn(\vect{e}_i) \), so the rank–nullity
theorem (@eq-proj-rank-nullity) applied to that restriction gives the identity. The intersection is
\( \spn(\vect{e}_i) \) if \( \vect{e}_i\in\C(\X) \) and \( \{\bzero\} \) otherwise, so the rank falls by one
exactly when \( \vect{e}_i\in\C(\X) \), that is, when \( h_{ii}=1 \).

(b) Leverages are nonnegative and sum to \( r=\rank(\X) \) by @prp-proj-leverage(b), so \( k \)
entries equal to one force \( k\le r \) and leave \( r-k \) for the rest.

(c) Suppose \( h_{ii}=1 \), so \( \M\vect{e}_i=\vect{e}_i \) and, by @prp-proj-leverage(d), row \( i \) and
column \( i \) of \( \M \) are \( \vect{e}_i\T \) and \( \vect{e}_i \). The matrix \( \bD\M\bD\T \) is symmetric, and
since \( \bD\T\bD=\I-\vect{e}_i\vect{e}_i\T \),
\[
(\bD\M\bD\T)^2=\bD\M(\I-\vect{e}_i\vect{e}_i\T)\M\bD\T=\bD(\M-\vect{e}_i\vect{e}_i\T)\bD\T=\bD\M\bD\T ,
\]
because \( \bD\vect{e}_i=\bzero \). Its column space is \( \bD\C(\M)=\C(\X_{(i)}) \), so it is the
projection of @thm-proj-sym-idem. The reduced fit is therefore
\( \bD\M\bD\T\y_{(i)}=\bD\M(\I-\vect{e}_i\vect{e}_i\T)\y=\bD\M\y-y_i\,\bD\vect{e}_i=\bD\M\y \),
whose \( j \)th entry is \( \hat y_j \): the other fitted values do not move, and they never
involved \( y_i \) in the first place, since \( m_{ji}=0 \). The response \( y_i \) influences no
fitted value but its own, which it reproduces exactly: the model has a free coordinate reserved
for case \( i \), and the data can neither check it nor borrow strength for it.
:::

::: {#exr-proj-leverage-replication}
[C2]

Replication bounds leverage. Suppose \( m \) cases share the same row of \( \X \), say
\( \x_{(i)}=\x_{(j)} \) for all \( i,j\in D \) with \( \lvert D\rvert=m \).

::: {.enumerate options="label=(\alph*)"}
1. Show that rows \( i \) and \( j \) of \( \M \) are equal for \( i,j\in D \), and hence that
   \( m_{ij}=h_{ii}=h_{jj} \) for all \( i,j\in D \).

2. Deduce \( h_{ii}\le1/m \), with equality iff \( m_{il}=0 \) for every \( l\notin D \). Check the
   bound against the one-way layout of @exr-proj-leverage-monotone.

3. A design point of high leverage is to be measured again. Explain what (b) says about
   repeating it \( m \) times, and why this does not contradict @exr-proj-leverage-monotone.
:::
:::

::: {.solution}
(a) By @thm-proj-ls-projection, \( \M\y=\X\hbeta \) for any least squares estimate, so the \( i \)th
fitted value is \( \x_{(i)}\T\hbeta \). If \( \x_{(i)}=\x_{(j)} \), then \( \hat y_i=\hat y_j \) for every
\( \y \), which forces rows \( i \) and \( j \) of \( \M \) to be equal. Taking the \( k \)th entry with
\( k\in D \) and using symmetry, \( m_{ik}=m_{kk}=h_{kk} \), and also \( m_{ik}=m_{ki}=h_{ii} \), so all
these numbers coincide.

(b) By @prp-proj-leverage(a),
\( h_{ii}=\sum_lm_{il}^2\ge\sum_{l\in D}m_{il}^2=m\,h_{ii}^2 \), and \( h_{ii}>0 \) unless the row is
\( \bzero \); dividing gives \( h_{ii}\le1/m \). Equality needs the discarded terms to vanish, that
is, \( m_{il}=0 \) for \( l\notin D \). In a one-way layout the \( m \) members of a group are exactly
such a set, and the projection averages within groups, so \( m_{il}=0 \) outside the group and the
bound is attained: \( h_{ii}=1/m \).

(c) Replication is the one design change that caps leverage: with \( m \) copies of a point, none
of them can have leverage above \( 1/m \), whatever the rest of the design, so an isolated design
point stops being able to dictate its own fitted value. There is no conflict with
@exr-proj-leverage-monotone, which adds a *column* to a fixed set of cases; here we add
*cases*, and the total \( \sum_ih_{ii}=\rank(\X) \) is redistributed over more of them.
:::

::: {#exr-proj-leverage-add-case}
[C3]

Let \( \X \) have full column rank \( p \), write \( \A=\X\T\X \), and append one more case with row
\( \x\T \), giving \( \X_+ \) with \( \X_+\T\X_+=\A+\x\x\T \). Put \( d=\x\T\A^{-1}\x \).

::: {.enumerate options="label=(\alph*)"}
1. Using @thm-mat-woodbury, show that the new case has leverage \( d/(1+d) \), which is always
   less than one.

2. Show that no old leverage increases: the new leverage of case \( i \) is
   \( h_{ii}-(\x_{(i)}\T\A^{-1}\x)^2/(1+d) \). When is it unchanged?

3. Reconcile (b) with \( \sum_ih_{ii}=p \) in both designs.
:::
:::

::: {.solution}
(a) By @eq-mat-sherman-morrison,
\( (\A+\x\x\T)^{-1}=\A^{-1}-\A^{-1}\x\x\T\A^{-1}/(1+d) \), and \( 1+d>0 \) because \( \A^{-1} \) is
positive definite. So the new case has leverage
\( \x\T(\A+\x\x\T)^{-1}\x=d-d^2/(1+d)=d/(1+d)<1 \).

(b) Applying the same formula to \( \x_{(i)} \) gives
\( h_{ii}^{+}=h_{ii}-(\x_{(i)}\T\A^{-1}\x)^2/(1+d)\le h_{ii} \), with equality iff
\( \x_{(i)}\T\A^{-1}\x=0 \), that is, iff the two rows are orthogonal in the metric of \( \A^{-1} \).

(c) The old leverages fall by \( \sum_i(\x_{(i)}\T\A^{-1}\x)^2/(1+d)=\x\T\A^{-1}\A\A^{-1}\x/(1+d)=d/(1+d) \),
which is exactly the leverage of the new case. The total stays at \( p \): one more case shares
the same \( p \) units of leverage.
:::
