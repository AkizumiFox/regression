# Masking and swamping

Every diagnostic so far removes *one* case. But bad cases often come in groups: a batch of records in the wrong
units, or a sensor that failed for an afternoon. **Masking** is when a bad case escapes
detection because others hide it: deleting it alone leaves its companions to hold the fit in place. **Swamping** is when
good cases are flagged because the bad ones have dragged the fit away from them. This section deletes sets of cases and
proves how a cluster masks itself and swamps the rest.

## Deleting a set of cases

Let \( D\subset\{1,\dots,n\} \) be a set of \( d \) cases. Write \( \X_D \) for the \( d\times p \) matrix of their rows,
\( \X_{(D)} \) for \( \X \) with those rows removed, \( \he_D \) for their residuals in the full fit, and
\[
\bH_D=\X_D(\X\T\X)^{-1}\X_D\T
\]
for the \( d\times d \) block of the hat matrix that belongs to them. Let \( \hbeta_{(D)} \) and \( \text{SSE}_{(D)} \) come from the least
squares fit to the cases not in \( D \), and let \( \mathbf{E}_D \) be the \( n\times d \) matrix whose columns are the coordinate vectors
\( \vect{e}_i \), \( i\in D \).

::: {#thm-res-group-deletion}
[Deleting a set of cases]

Let \( \X \) and \( \X_{(D)} \) have full column rank. Then \( \I-\bH_D \) is nonsingular, and

::: {.enumerate options="label=(\alph*)"}
1. \( \y_D-\X_D\hbeta_{(D)}=(\I-\bH_D)^{-1}\he_D \);

2. \( \hbeta-\hbeta_{(D)}=(\X\T\X)^{-1}\X_D\T(\I-\bH_D)^{-1}\he_D \);

3. \( \text{SSE}_{(D)}=\text{SSE}-\he_D\T(\I-\bH_D)^{-1}\he_D \);

4. if \( n-p-d\ge1 \), then under @eq-opt-normal-model
   \[
   F_D=\frac{\he_D\T(\I-\bH_D)^{-1}\he_D/d}{\text{SSE}_{(D)}/(n-p-d)}\sim F(d,\ n-p-d),
   \]
   and \( F_D \) is the \( F \) statistic for \( \bgamma=\bzero \) in the mean-shift model \( \E(\Y)=\X\bbeta+\mathbf{E}_D\bgamma \).
:::

:::

::: {.proof}
Let \( \mathbf{Z}=[\X,\mathbf{E}_D] \). If \( \X\bb+\mathbf{E}_D\mathbf{c}=\bzero \), the rows outside \( D \) give \( \X_{(D)}\bb=\bzero \), so \( \bb=\bzero \), and then
\( \mathbf{c}=\bzero \). So \( \mathbf{Z} \) has full column rank. By @lem-proj-fwl-split(b) with \( \X_1=\mathbf{E}_D \), the matrix
\( (\I-\M)\mathbf{E}_D \) has full column rank, so its Gram matrix \( \mathbf{E}_D\T(\I-\M)\mathbf{E}_D=\I-\bH_D \) is nonsingular.

Fit the mean-shift model. By @thm-proj-fwl(a) with \( \X_1=\mathbf{E}_D \) and \( \X_2=\X \), the estimated shifts are
\[
\hat{\bgamma}=\bigl(\mathbf{E}_D\T(\I-\M)\mathbf{E}_D\bigr)^{-1}\mathbf{E}_D\T(\I-\M)\y=(\I-\bH_D)^{-1}\he_D .
\]
On the other hand, for any \( \bb \) the residual sum of squares \( \sum_{i\notin D}(y_i-\x_{(i)}\T\bb)^2+\norm{\y_D-\X_D\bb-\mathbf{c}}^2 \) is
minimized over \( \mathbf{c} \) by \( \mathbf{c}=\y_D-\X_D\bb \), which removes the second term. What remains is minimized by
\( \bb=\hbeta_{(D)} \), uniquely because \( \X_{(D)} \) has full column rank. So the mean-shift fit has coefficients
\( \hbeta_{(D)} \) and \( \hat{\bgamma}=\y_D-\X_D\hbeta_{(D)} \), and its residual sum of squares is \( \text{SSE}_{(D)} \). Comparing the two
expressions for \( \hat{\bgamma} \) gives (a). By @thm-proj-fwl(c), the coefficients of \( \X \) in the mean-shift fit are
\( (\X\T\X)^{-1}\X\T(\y-\mathbf{E}_D\hat{\bgamma})=\hbeta-(\X\T\X)^{-1}\X_D\T\hat{\bgamma} \), which gives (b). By @lem-proj-fwl-split(c) the
projection onto \( \C(\mathbf{Z}) \) is \( \M \) plus the projection onto \( \C\bigl((\I-\M)\mathbf{E}_D\bigr) \), so the fall in the residual sum of
squares from the fit of \( \X \) to the fit of \( \mathbf{Z} \) is
\[
\y\T(\I-\M)\mathbf{E}_D(\I-\bH_D)^{-1}\mathbf{E}_D\T(\I-\M)\y=\he_D\T(\I-\bH_D)^{-1}\he_D ,
\]
which gives (c). Finally \( \C(\X)\subseteq\C(\mathbf{Z}) \) with dimensions \( p \) and \( p+d \), and \( F_D \) is the statistic @eq-glh-F for this
pair of nested models, so (d) is @thm-glh-f-test(d).
:::

With \( d=1 \) this is @thm-res-deletion, and \( F_D=t_i^2 \) (@exr-res-group-single). Part (d) is exact only for a set chosen
in advance; a Bonferroni correction over all \( \binom nd \) subsets is valid but very conservative.

## How a cluster masks itself

The simplest mechanism is a group of identical bad cases: \( n_0 \) clean cases plus \( m \) copies of one bad case.

::: {#prp-res-masking}
[Masking and swamping by a cluster]

Let \( \X_0 \) (\( n_0\times p \), full column rank) and \( \y_0 \) be the clean cases, with least squares estimate \( \hbeta_0 \), residuals
\( \breve{\varepsilon}_j=y_j-\x_{(j)}\T\hbeta_0 \) and residual sum of squares \( \text{SSE}_0 \). Add \( m\ge1 \) cases, each with row
\( \x_\star\T \) and response \( y_\star \). Put \( \A=\X_0\T\X_0 \),
\[
h_\star=\x_\star\T\A^{-1}\x_\star,\qquad \delta=y_\star-\x_\star\T\hbeta_0 ,
\]
the leverage of \( \x_\star \) relative to the clean design and the distance of the cluster from the clean fit. In the least
squares fit to all \( n_0+m \) cases:

::: {.enumerate options="label=(\alph*)"}
1. each cluster case has leverage \( h=h_\star/(1+mh_\star)<1/m \);

2. each cluster case has residual \( \delta/(1+mh_\star) \);

3. deleting one cluster case gives it prediction residual \( \delta/\bigl(1+(m-1)h_\star\bigr) \) and changes the estimate by
   \( \hbeta-\hbeta_{(i)}=\A^{-1}\x_\star\,\delta\big/\bigl[(1+mh_\star)(1+(m-1)h_\star)\bigr] \);

4. deleting the whole cluster gives each member prediction residual \( \delta \), and
   \( \hbeta-\hbeta_{(D)}=m\bigl(1+(m-1)h_\star\bigr)\bigl(\hbeta-\hbeta_{(i)}\bigr) \);

5. clean case \( j \) has residual \( \breve{\varepsilon}_j-m\,\delta\,h_{j\star}/(1+mh_\star) \), where \( h_{j\star}=\x_{(j)}\T\A^{-1}\x_\star \), and
   \( \text{SSE}=\text{SSE}_0+m\delta^2/(1+mh_\star) \).
:::

:::

::: {.proof}
Let \( D \) be the cluster. Then \( \X\T\X=\A+m\x_\star\x_\star\T \), and
\( (\A+m\x_\star\x_\star\T)\A^{-1}\x_\star=(1+mh_\star)\x_\star \), so
\[
(\X\T\X)^{-1}\x_\star=\frac{\A^{-1}\x_\star}{1+mh_\star}.
\]{#eq-res-cluster-solve}

(a) \( h=\x_\star\T(\X\T\X)^{-1}\x_\star=h_\star/(1+mh_\star) \), and \( mh=mh_\star/(1+mh_\star)<1 \).
(b) and (d) Deleting \( D \) leaves the clean cases, so \( \hbeta_{(D)}=\hbeta_0 \) and every cluster case has prediction residual \( \delta \).
All rows of \( \X_D \) equal \( \x_\star\T \), so \( \bH_D=h\,\bone\bone\T \), and by @thm-res-group-deletion(a)
\( \he_D=(\I-h\bone\bone\T)\,\delta\bone=(1-mh)\,\delta\bone=\delta\bone/(1+mh_\star) \). By @thm-res-group-deletion(b) and @eq-res-cluster-solve,
\( \hbeta-\hbeta_0=(\X\T\X)^{-1}\X_D\T\,\delta\bone=m\delta(\X\T\X)^{-1}\x_\star=m\delta\A^{-1}\x_\star/(1+mh_\star) \).
(c) By @thm-res-deletion(d) and (b), with \( 1-h=\bigl(1+(m-1)h_\star\bigr)/(1+mh_\star) \), the prediction residual is
\( \delta/\bigl(1+(m-1)h_\star\bigr) \). By @thm-res-deletion(b) and @eq-res-cluster-solve,
\( \hbeta-\hbeta_{(i)}=(\X\T\X)^{-1}\x_\star\cdot\delta/\bigl(1+(m-1)h_\star\bigr) \), which is the stated vector. Comparing with the expression for
\( \hbeta-\hbeta_0 \) in the previous step gives the ratio in (d).
(e) \( y_j-\x_{(j)}\T\hbeta=\breve{\varepsilon}_j-\x_{(j)}\T(\hbeta-\hbeta_0) \), and \( \x_{(j)}\T(\hbeta-\hbeta_0)=m\delta h_{j\star}/(1+mh_\star) \). For
the sum of squares, @thm-res-group-deletion(c) gives
\( \text{SSE}=\text{SSE}_0+\he_D\T(\I-\bH_D)^{-1}\he_D \), and \( \he_D\T(\I-\bH_D)^{-1}\he_D=\delta^2\,\bone\T(\I-h\bone\bone\T)\bone=m\delta^2(1-mh) \),
using \( (\I-h\bone\bone\T)^{-1}\he_D=\delta\bone \) from part (b).
:::

Each part is a way of hiding. **Leverage is masked** (a): however remote the cluster, each copy has leverage below \( 1/m \); a single copy would have \( h_\star/(1+h_\star) \), near one.
**Residuals are masked** (b): each is only \( 1/(1+mh_\star) \) of the cluster's distance from the clean fit. **Deletion is
masked** (c, d): deleting one copy leaves \( m-1 \) to hold the fit, recovering only \( 1/\bigl(1+(m-1)h_\star\bigr) \) of \( \delta \).
**Good cases are swamped** (e): clean residuals shift in proportion to their cross-leverage with the cluster, and \( s^2 \) is
inflated.

::: {#exm-res-masking}
[Four identical bad records]

Thirty clean cases have \( x \) uniform on \( (0,10) \) and \( y=2+0.5x+\varepsilon \), \( \sigma=0.5 \). Four identical records are added at
\( x=20 \), \( y=2 \), where the clean line predicts about \( 12 \) ([Figure 20.5.1](#fig-res-masking)(a)). The slope drops from
\( 0.514 \) (clean cases) to \( 0.007 \), and \( s \) rises from \( 0.608 \) to \( 1.632 \). Here \( h_\star=1.126 \) (it can exceed one,
since \( x=20 \) is outside the clean design, as in @prp-ci-new-leverage) and \( \delta=-10.150 \). Each bad case has leverage
\( 0.205 \), below \( 1/m=0.25 \), where a single copy would have \( 0.530 \); residual \( -1.844 \); prediction residual \( -2.318 \);
and \( t_i=-1.28 \). Deleting the group changes \( \hbeta \) \( 17.51 \) times as much as deleting one member.

The single-case diagnostics fail as predicted. The Bonferroni test (critical value \( 3.49 \)) flags nothing, and the only case
with \( \lvert t_i\rvert>2 \) is a clean one: case \( 30 \), the clean case with largest \( x \) (\( 9.86 \)), has \( t_i=2.25 \), because the
dragged line passes below it. Cook's distance does a little better: each bad case has \( D_i=0.206 \), above \( 4/n=0.118 \) (no
clean case exceeds it; the largest is \( 0.100 \)) but far below the median of \( F(2,32) \). A single bad record alone would
have had \( t_i=-11.45 \) and \( D_i=13.46 \). Deleting the four together tells the truth: their prediction residuals are all
\( \delta \), and @thm-res-group-deletion(d) gives \( F=50.6 \) on \( 4 \) and \( 28 \) degrees of freedom, \( p \)-value \( 2.0\times 10^{-12} \).
The difficulty is not in testing the cluster but in finding it.
:::

::: {when-format="html"}
![**Figure 20.5.1.** Masking and swamping. (a) Clean cases, four identical bad records (diamond), and least squares lines
with and without them. (b) Externally studentized residuals (bad records last) and the \( 5\% \) Bonferroni cut-offs.](masking.svg){#fig-res-masking width=100%}
:::

::: {when-format="pdf"}
![Masking and swamping. (a) Clean cases, four identical bad records (diamond), and least squares lines
with and without them. (b) Externally studentized residuals (bad records last) and the \( 5\% \) Bonferroni cut-offs.](masking.pdf){width=100%}
:::

```{.python .run #cell-masking-single}
import numpy as np
from scipy import stats

def single_case(X, y):
    """Leverages, residuals, s^2, externally studentized residuals and Cook's distances."""
    n, p = X.shape
    Q, _ = np.linalg.qr(X)
    h = np.sum(Q ** 2, axis=1)
    e = y - Q @ (Q.T @ y)
    s2 = e @ e / (n - p)
    s2_del = (e @ e - e ** 2 / (1 - h)) / (n - p - 1)
    t = e / np.sqrt(s2_del * (1 - h))
    cook = e ** 2 * h / (p * s2 * (1 - h) ** 2)
    return h, e, s2, t, cook

rng = np.random.default_rng(20)
n0, m = 30, 4                                      # clean cases, size of the cluster
x_clean = np.sort(rng.uniform(0, 10, n0))
y_clean = 2 + 0.5 * x_clean + rng.normal(0, 0.5, n0)
x0, y0 = 20.0, 2.0                                 # four identical records
x = np.r_[x_clean, np.full(m, x0)]
y = np.r_[y_clean, np.full(m, y0)]
X = np.column_stack([np.ones(len(x)), x])
n, p = X.shape
cluster = np.arange(n0, n)

h, e, s2, t, cook = single_case(X, y)
crit = stats.t.isf(0.05 / (2 * n), n - p - 1)      # Bonferroni outlier test at 5%
print("slope with the cluster:", np.linalg.lstsq(X, y, rcond=None)[0][1])
print("cluster: h =", h[-1].round(3), " t =", t[-1].round(2), " Cook's D =", cook[-1].round(3))
i = int(np.argmax(np.abs(t)))
print(f"largest |t|: case {i + 1} at x = {x[i]:.2f}, t = {t[i]:.2f}; Bonferroni cut-off {crit:.2f}")
```

```{.python .run #cell-masking-group}
D = cluster                                        # delete the whole cluster at once
K = np.linalg.inv(X.T @ X)
H_D = X[D] @ K @ X[D].T
pred_D = np.linalg.solve(np.eye(m) - H_D, e[D])    # y_D minus its prediction from the rest
sse = e @ e
sse_D = sse - e[D] @ pred_D
F = (e[D] @ pred_D / m) / (sse_D / (n - p - m))
print("group prediction residuals:", pred_D.round(3))
print(f"F = {F:.1f} on {m} and {n - p - m} degrees of freedom,"
      f" p-value {stats.f.sf(F, m, n - p - m):.1e}")
```

## Finding the set

There are \( \binom nd \) candidate sets of each size \( d \), far too many to search. The natural shortcut, deleting the case
with the largest \( \lvert t_i\rvert \) and repeating, starts in the wrong place: in @exm-res-masking it would first remove case
\( 30 \), a clean case. Three ideas get around the combinatorics. First, **start from a fit the cluster cannot capture**, one
determined by about half the data (the next section); the cases it flags form a candidate set to examine and test, remembering
that it was chosen by looking. Second, **grow a clean subset**: the forward search of Atkinson and Riani (2000) adds cases
one at a time in order of their residuals from a robust start, and a cluster enters late and all at once, as a jump in the
monitored residuals and estimates; Hadi and Simonoff (1993) proposed a related procedure. Third, **unmask leverage
separately**: by @prp-res-masking(a) the hat matrix is fooled by clusters, since leverage is a Mahalanobis distance built
from the sample covariance, which a cluster inflates in its own direction. Rousseeuw and van Zomeren (1990) replaced the mean
and covariance by estimates that ignore up to half the cases.

::: {.idea}
Single-case diagnostics answer "what if this case alone were different?", the wrong question when bad cases come
in groups. Group deletion asks the right one, but only for a suspected group; finding the group needs a fit it cannot
capture.
:::

## Exercises

### A. Check your understanding

::: {#exr-res-masking-numbers}
[A1]

Check the leverage \( 0.205 \) and the residual \( -1.844 \) in @exm-res-masking from \( h_\star=1.126 \), \( \delta=-10.150 \) and \( m=4 \).
:::


::: {#exr-res-group-single}
[A2]

Show that @thm-res-group-deletion with \( D=\{i\} \) gives parts (b), (d) and (e) of @thm-res-deletion, and that \( F_D=t_i^2 \).
:::

::: {.solution}
With \( d=1 \), \( \bH_D=h_{ii} \) and \( \he_D=\hat{\varepsilon}_i \), so (a) is the prediction residual \( \hat{\varepsilon}_i/(1-h_{ii}) \), (b) is
@thm-res-deletion(b), and (c) is \( \text{SSE}_{(i)}=\text{SSE}-\hat{\varepsilon}_i^2/(1-h_{ii}) \). Then
\( F_D=\hat{\varepsilon}_i^2/\{(1-h_{ii})s_{(i)}^2\}=t_i^2 \).
:::

### B. Practice

::: {#exr-res-rank-one-inverse}
[B1]

Show that for \( \bH_D=h\bone\bone\T \) (\( d\times d \)) with \( dh\ne1 \), \( (\I-h\bone\bone\T)^{-1}=\I+h\bone\bone\T/(1-dh) \), and use it to recompute the
group term \( \he_D\T(\I-\bH_D)^{-1}\he_D \) in the proof of @prp-res-masking(e) from \( \he_D=\delta\bone/(1+mh_\star) \).
:::


::: {#exr-res-swamping-sign}
[B2]

In simple regression with clean mean \( \bar x_0 \) and \( S_0=\sum(x_j-\bar x_0)^2 \), show that
\( h_{j\star}=1/n_0+(x_j-\bar x_0)(x_\star-\bar x_0)/S_0 \). In @exm-res-masking, with \( x_\star \) far to the right and \( \delta<0 \), which clean
residuals does @prp-res-masking(e) push up, and which down? Explain why case \( 30 \) was flagged.
:::

::: {.solution}
Use the centred rows \( (1,x-\bar x_0) \), which span the same column space; then \( \A^{-1}=\operatorname{diag}(1/n_0,1/S_0) \),
so \( h_{j\star}=1/n_0+(x_j-\bar x_0)(x_\star-\bar x_0)/S_0 \). The shift in clean
residual \( j \) is \( -m\delta h_{j\star}/(1+mh_\star) \), which has the sign of \( h_{j\star} \) because \( \delta<0 \). Since \( x_\star-\bar x_0>0 \),
\( h_{j\star} \) grows with \( x_j \): positive on the right of the clean data and negative on the far left (once
\( (x_j-\bar x_0)(x_\star-\bar x_0)/S_0<-1/n_0 \)). Residuals on the right are pushed up and those on the far left down. Case \( 30 \) is the
clean case with the largest \( x \), so it has the largest upward shift, and its studentized residual crosses two.
:::

::: {#exr-res-partial-cluster}
[B3]

In the setting of @prp-res-masking, delete \( k<m \) of the cluster cases. Show that each deleted case has prediction residual
\( \delta/\bigl(1+(m-k)h_\star\bigr) \). What does this say about searching for a cluster of unknown size by deleting sets of size \( k \)?
:::


### C. Going deeper

::: {#exr-res-two-clusters}
[C1]

Place two clusters of \( m \) identical bad cases at the two ends of a straight-line design, one above the clean line and one below.
Describe, using @thm-res-group-deletion, what single-case deletion, deletion of one cluster, and deletion of both clusters report.
Is it possible for each cluster to mask the other?
:::

::: {#exr-res-group-cook}
[C2]

For a set \( D \) of cases, define \( D_D=(\hbeta_{(D)}-\hbeta)\T\X\T\X(\hbeta_{(D)}-\hbeta)/(p\,s^2) \). Using
@thm-res-group-deletion(b), show that
\( D_D=\he_D\T(\I-\bH_D)^{-1}\bH_D(\I-\bH_D)^{-1}\he_D/(p\,s^2) \).
:::
