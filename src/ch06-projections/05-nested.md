# Nested models and sequential projection

Statistical questions about linear models are almost always comparisons. Does a
regressor matter? Are several group means equal? Is a curve needed, or will a line do? Each
comparison sets a full model against a *reduced* model whose mean space is
a subspace of the full one. This section studies the geometry of such nested
pairs. Chapter 9 and Chapter 11 build the analysis of
variance and the \( F \) test on it.

## Two nested subspaces

Let \( \C(\X_0)\subseteq\C(\X) \), with ranks \( r_0\le r \). Write \( \Mo \) and \( \M \) for the
two projections.

::: {#thm-proj-nested}
[Nested projections]

If \( \C(\X_0)\subseteq\C(\X) \), then:

::: {.enumerate options="label=(\alph*)"}
1. \( \M\Mo=\Mo\M=\Mo \);

2. \( \M-\Mo \) is the orthogonal projection onto \( \C(\X_0)\perpc\cap\C(\X) \), the
           orthogonal complement of \( \C(\X_0) \) *within* \( \C(\X) \);

3. \( \rank(\M-\Mo)=r-r_0 \);

4. \( (\I-\Mo)-(\I-\M)=\M-\Mo \), and \( \Mo(\M-\Mo)=(\M-\Mo)(\I-\M)=\bzero \).
:::

:::

::: {.proof}
(a) Every column of \( \Mo \) lies in \( \C(\X_0)\subseteq\C(\X) \) and is left unchanged
by \( \M \), so \( \M\Mo=\Mo \). Transposing gives \( \Mo\M=\Mo \).
(b) \( \M-\Mo \) is symmetric, and by (a)
\( (\M-\Mo)^2=\M-\M\Mo-\Mo\M+\Mo=\M-\Mo \). By @thm-proj-sym-idem it projects
onto its column space. If \( \bv\in\C(\X) \) and \( \bv\perp\C(\X_0) \), then
\( (\M-\Mo)\bv=\bv-\bzero=\bv \), so \( \bv\in\C(\M-\Mo) \). Conversely, a vector
\( \bv=(\M-\Mo)\bu \) lies in \( \C(\X) \) and satisfies \( \Mo\bv=(\Mo-\Mo)\bu=\bzero \).
(c) \( \tr(\M-\Mo)=r-r_0 \) and @prp-proj-trace-rank.
(d) Direct multiplication using (a).
:::

The theorem gives a three-way orthogonal decomposition of \( \Real^n \),
\[
\Real^n=\underbrace{\C(\X_0)}_{\dim r_0}
\;\dirsum\;
\underbrace{\C(\X_0)\perpc\cap\C(\X)}_{\dim r-r_0}
\;\dirsum\;
\underbrace{\C(\X)\perpc}_{\dim n-r},
\]
and correspondingly \( \I=\Mo+(\M-\Mo)+(\I-\M) \), with the three terms mutually
orthogonal projections. Applying the decomposition to \( \y \) and using Pythagoras,
\[
\norm{\y}^2=\norm{\Mo\y}^2+\norm{(\M-\Mo)\y}^2+\norm{(\I-\M)\y}^2 .
\]{#eq-proj-three-way}

The middle term has a second reading that is the key to model comparison.
Since \( (\I-\Mo)=(\M-\Mo)+(\I-\M) \) is an orthogonal sum,
\[
\norm{(\M-\Mo)\y}^2
=\norm{(\I-\Mo)\y}^2-\norm{(\I-\M)\y}^2
=\text{SSE}_0-\text{SSE}.
\]{#eq-proj-extra-ss}

The *drop in residual sum of squares* from the reduced to the full model is
the squared length of a projection onto an \( (r-r_0) \)-dimensional space. That is
why it has \( r-r_0 \) degrees of freedom, and why it is independent of SSE under
normality: the two lie in orthogonal subspaces ([Chapter 4](../ch04-quadratic-forms/index.html)). It
also shows that SSE can never increase when the model space grows, with
equality iff \( \M\y=\Mo\y \).

\begin{center}
\begin{tikzpicture}[scale=1.05, >=Stealth, line cap=round]
\definecolor{ink}{HTML}{1B1F27}\definecolor{accent}{HTML}{2F5F96}\definecolor{thread}{HTML}{74509A}\definecolor{draft}{HTML}{A45C25}
\fill[accent!8] (-0.6,-0.5) -- (4.6,-0.5) -- (5.6,1.6) -- (0.4,1.6) -- cycle;
\node[accent!80!black, font=\small] at (5.1,0.2) {$\C(\X)$};
\draw[accent!60, line width=2.2pt, opacity=0.5] (-0.4,-0.25) -- (4.9,0.1);
\node[accent!80!black, font=\small, anchor=north] at (4.4,0.02) {$\C(\X_0)$};
\coordinate (O) at (0,0);
\coordinate (Y) at (3.2,3.1);
\coordinate (MY) at (3.2,1.0);
\coordinate (M0Y) at (3.28,-0.01);
\draw[->, thick, ink] (O) -- (Y) node[above right] {$\y$};
\draw[->, thick, accent] (O) -- (MY) node[pos=0.55, above left, font=\small] {$\M\y$};
\draw[->, thick, accent!60!black] (O) -- (M0Y) node[pos=0.5, below, font=\small] {$\Mo\y$};
\draw[->, thick, draft] (MY) -- (Y) node[midway, right, font=\small] {$(\I-\M)\y$};
\draw[->, thick, thread] (M0Y) -- (MY) node[midway, right, font=\small] {$(\M-\Mo)\y$};
\draw[dashed, black!40] (M0Y) -- (Y);
\draw[black!60] ($(MY)+(0,0.22)$) -- ++(0.2,0) -- ++(0,-0.22);
\draw[black!60] ($(M0Y)+(0,0.22)$) -- ++(-0.2,0.003) -- ++(0,-0.22);
\end{tikzpicture}
\end{center}

[**Figure 6.5.1.** Nested projections. The fitted vector of the reduced model, \( \Mo\y \), is also
the projection of \( \M\y \) onto \( \C(\X_0) \). The three coloured pieces
\( \Mo\y \), \( (\M-\Mo)\y \) and \( (\I-\M)\y \) are mutually orthogonal. The dashed segment is
the reduced model's residual, whose squared length is the sum of the squared
lengths of the two segments meeting at \( \M\y \).]{#fig-proj-nested}

Part (a) also says \( \Mo\y=\Mo(\M\y) \). *Fitting the reduced model to the data
gives the same result as fitting it to the full model's fitted values.* Projection
can be done in stages ([Figure 6.5.1](05-nested.html#fig-proj-nested)).

::: {#exm-proj-simple-nested}
[Intercept inside simple regression]

Let \( \X_0=\bone \) and \( \X=[\bone,\x] \) with \( \x \) not constant. Then
\( \C(\X_0)\perpc\cap\C(\X) \) is spanned by the centred regressor
\( \x-\bar x\bone=(\I-n^{-1}\bone\bone\T)\x \), which lies in \( \C(\X) \) and is
orthogonal to \( \bone \). By @exm-proj-line,
\[
\M-\Mo=\frac{(\x-\bar x\bone)(\x-\bar x\bone)\T}{S_{xx}},
\qquad S_{xx}=\sum_i(x_i-\bar x)^2,
\]
and so
\[
\M\y=\bar y\bone+\frac{S_{xy}}{S_{xx}}(\x-\bar x\bone),
\qquad S_{xy}=\sum_i(x_i-\bar x)(y_i-\bar y).
\]
Reading off coordinates, the slope is \( S_{xy}/S_{xx} \) and the intercept is
\( \bar y-\bar x\,S_{xy}/S_{xx} \). These are the familiar formulas, obtained without
solving any equations, and the decomposition @eq-proj-three-way becomes
\[
\sum_i y_i^2 = n\bar y^2+\frac{S_{xy}^2}{S_{xx}}+\text{SSE}.
\]
:::

## Chains of subspaces

The argument extends to any chain
\( \C(\X_1)\subseteq\C(\X_2)\subseteq\dots\subseteq\C(\X_k) \) with projections
\( \M_1,\dots,\M_k \). Setting \( \M_0=\bzero \) and \( \M_{k+1}=\I \), the differences
\( \M_j-\M_{j-1} \) for \( j=1,\dots,k+1 \) are mutually orthogonal projections summing
to \( \I \). The corresponding lengths
\[
\norm{(\M_j-\M_{j-1})\y}^2=\text{SSE}_{j-1}-\text{SSE}_j
\]
are the **sequential sums of squares** of Chapter 9. Their
values depend on the *order* in which the subspaces are nested. A regressor
added early can take credit that it would not get if added last. The
Frisch–Waugh–Lovell theorem of [Section 6.6](06-fwl.html) describes what the regressor
added last actually contributes.

::: {#exm-proj-sequential-order}
[Order matters]

For the state data of @exm-proj-fwl-crime, enter the three regressors after the
intercept in two different orders. The listing computes each squared length
\( \norm{(\M_j-\M_{j-1})\y}^2 \) directly from the projections.

| Order A |  | Order B |  |
|---|---|---|---|
| poverty | 102.49 | single-parent | 162.84 |
| single-parent | 80.04 | urban | 3.59 |
| urban | 0.27 | poverty | 16.37 |
| residual | 101.96 | residual | 101.96 |

Entered first, poverty accounts for \( 102.49 \) of the corrected total
sum of squares \( 284.76 \). Entered last, it accounts for only
\( 16.37 \). The residual sum of squares, and the total explained by
all three regressors together (\( 182.80 \)), do not depend on the order,
because the final subspace does not. The individual pieces do depend on the order,
because the intermediate subspaces differ. Poverty and single parenthood are
correlated, so whichever enters first takes the variation they share. The
piece for the *last* regressor is the only one with an order-free meaning. By
@thm-proj-fwl, it is the squared length of the projection of \( \y \) onto the part of
that regressor orthogonal to all the others.
:::

```{.python .run #cell-sequential-ss-sequential}
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
cols = {"poverty": data["poverty"], "single": data["single"], "urban": data["urban"]}

def projection(Z):
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T                       # fine for n = 50; never do this for large n

def sequential_ss(y, columns):
    """Squared lengths ||(M_j - M_{j-1}) y||^2 as columns enter one at a time."""
    n = len(y)
    Z = np.ones((n, 1))
    M_prev = projection(Z)
    pieces = {"mean": y @ M_prev @ y}
    for name, x in columns:
        Z = np.column_stack([Z, x])
        M = projection(Z)
        pieces[name] = y @ (M - M_prev) @ y
        M_prev = M
    pieces["residual"] = y @ (np.eye(n) - M_prev) @ y
    return pieces

order_a = sequential_ss(y, [(k, cols[k]) for k in ["poverty", "single", "urban"]])
order_b = sequential_ss(y, [(k, cols[k]) for k in ["single", "urban", "poverty"]])
for name in ["poverty", "single", "urban", "residual"]:
    print(f"{name:9s}  order A: {order_a[name]:8.2f}   order B: {order_b[name]:8.2f}")
```

## Two factors and when the pieces are orthogonal

A sum of projections is a projection exactly when the pieces are orthogonal
(@thm-proj-sum). In designed experiments, orthogonality comes from
balance, and when it fails, it fails in ways that later chapters have to face.

::: {#exm-proj-two-factor}
[Additive two-factor model, one observation per cell]

Let factor \( A \) have \( a \) levels and factor \( B \) have \( b \) levels, with one observation
\( y_{ij} \) for each of the \( n=ab \) combinations. Let \( \bP_A \) project onto vectors
that depend only on the \( A \) level, so \( \bP_A\y \) has entries \( \bar y_{i\cdot} \), and
define \( \bP_B \) similarly. Let \( \bP_0=n^{-1}\bone\bone\T \). Both
\( \C(\bP_A) \) and \( \C(\bP_B) \) contain \( \bone \), so \( \bP_A\bP_0=\bP_B\bP_0=\bP_0 \) by
@thm-proj-nested.

The crucial computation is \( \bP_A\bP_B \). Averaging over \( j \) with \( i \) fixed a
vector that depends only on \( j \) gives the overall average of that vector. So for
every \( \y \) the entries of \( \bP_A(\bP_B\y) \) are all equal to
\( b^{-1}\sum_j\bar y_{\cdot j}=\bar y_{\cdot\cdot} \). Therefore
\[
\bP_A\bP_B=\bP_0,
\]
and so
\( (\bP_A-\bP_0)(\bP_B-\bP_0)=\bP_0-\bP_0-\bP_0+\bP_0=\bzero \). The two
*centred* factor spaces are orthogonal. By @thm-proj-sum, the additive model
space has projection \( \bP_0+(\bP_A-\bP_0)+(\bP_B-\bP_0)=\bP_A+\bP_B-\bP_0 \), and
\[
\I=\bP_0+(\bP_A-\bP_0)+(\bP_B-\bP_0)+(\I-\bP_A-\bP_B+\bP_0)
\]
is an orthogonal decomposition with ranks \( 1 \), \( a-1 \), \( b-1 \) and \( (a-1)(b-1) \).
The residual projection sends \( y_{ij} \) to
\( y_{ij}-\bar y_{i\cdot}-\bar y_{\cdot j}+\bar y_{\cdot\cdot} \).

The key step used the fact that every \( A \) level appears with every \( B \) level
equally often. If cells are missing or unequally replicated, \( \bP_A\bP_B\neq\bP_0 \),
the centred factor spaces are no longer orthogonal, and the sum of squares for
\( A \) depends on whether \( B \) was fitted first. Chapter 17 is about
that situation.
:::

## Reduced models defined by constraints

Reduced models often come as *constraints* on the parameters of the full model,
for example \( \beta_2=\beta_3 \) or \( \beta_4=0 \), rather than as a smaller model matrix.
In mean space, a constraint \( \bLambda\T\bbeta=\bzero \) on an estimable function
\( \bLambda\T=\bm P\T\X \) restricts \( \bmu=\X\bbeta \) to the subspace
\[
\mathcal S_0=\C(\X)\cap\Null(\bm P\T)=\{\bmu\in\C(\X):\bm P\T\bmu=\bzero\}.
\]
To test the constraint we need the piece of \( \C(\X) \) that \( \mathcal S_0 \) removes.

::: {#thm-proj-constraint-space}
Let \( \M \) project onto \( \C(\X) \) and let \( \bm P \) be an \( n\times q \) matrix. If
\( \mathcal S_0=\C(\X)\cap\Null(\bm P\T) \), then
\[
\mathcal S_0\perpc\cap\C(\X)=\C(\M\bm P).
\]
Hence the projection onto \( \mathcal S_0 \) is \( \M-\bP_{\C(\M\bm P)} \), and
\( \dim\mathcal S_0=r-\rank(\M\bm P) \).
:::

::: {.proof}
By @lem-proj-complement-intersection and @lem-proj-null-colspace,
\( \mathcal S_0\perpc=\C(\X)\perpc+\Null(\bm P\T)\perpc=\C(\X)\perpc+\C(\bm P) \).
If \( \bv\in\mathcal S_0\perpc\cap\C(\X) \), write \( \bv=\bw+\bm P\bm a \) with
\( \bw\perp\C(\X) \). Applying \( \M \) gives \( \bv=\M\bv=\M\bm P\bm a\in\C(\M\bm P) \).
Conversely, \( \M\bm P\bm a \) lies in \( \C(\X) \) and equals
\( \bm P\bm a-(\I-\M)\bm P\bm a\in\C(\bm P)+\C(\X)\perpc=\mathcal S_0\perpc \). The
remaining claims follow from @thm-proj-nested with the roles
\( \C(\X_0)=\mathcal S_0 \).
:::

The subspace \( \C(\M\bm P) \) is the *test space* of the hypothesis. Its dimension
\( \rank(\M\bm P) \) counts how many independent restrictions the constraint really
places on the mean. It can be smaller than \( q \) if some rows of the constraint are
redundant. Chapter 11 uses exactly this space for the
numerator of the \( F \) statistic.

## Exercises

### A. Check your understanding

::: {#exr-proj-two-factor-numeric}
[A1]

Construct the projections \( \bP_A \), \( \bP_B \) and \( \bP_0 \) of @exm-proj-two-factor for
\( a=3 \), \( b=2 \) and verify \( \bP_A\bP_B=\bP_0 \) numerically. Delete one observation and
recompute. Show that the centred factor spaces are no longer orthogonal, and find the
largest entry of \( (\bP_A-\bP_0)(\bP_B-\bP_0) \).
:::

### B. Practice

::: {#exr-proj-extra-ss-distance}
[B1]

With \( \C(\X_0)\subseteq\C(\X) \), show that \( \text{SSE}_0-\text{SSE}=\norm{\hY-\hY_0}^2 \),
the squared distance between the two fitted vectors.
:::

::: {.solution}
\( \hY-\hY_0=(\M-\Mo)\y \), and
\( \norm{(\M-\Mo)\y}^2=\text{SSE}_0-\text{SSE} \) by @eq-proj-extra-ss.
:::

::: {#exr-proj-equal-coefficients}
[B2]

In the model \( \E(\Y)=\beta_0\bone+\beta_1\x_1+\beta_2\x_2+\beta_3\x_3 \) with full-rank
\( \X \), the hypothesis \( \beta_2=\beta_3 \) defines a reduced model with model matrix
\( \X_0=[\bone,\x_1,\x_2+\x_3] \). Show that the test space \( \C(\X_0)\perpc\cap\C(\X) \) is
spanned by \( (\I-\Mo)\x_2 \), and that \( (\I-\Mo)\x_2=-(\I-\Mo)\x_3 \). Verify that this agrees
with @thm-proj-constraint-space for a suitable \( \bm P \).
:::

::: {.solution}
The two model spaces differ in
dimension by one, so the test space is a line. The vector \( (\I-\Mo)\x_2 \) lies in \( \C(\X) \)
and is orthogonal to \( \C(\X_0) \). It is nonzero, because \( \x_2\in\C(\X_0) \) would make
the columns of \( \X \) dependent. Since \( \x_2+\x_3\in\C(\X_0) \),
\( (\I-\Mo)(\x_2+\x_3)=\bzero \). For @thm-proj-constraint-space, write the constraint as
\( \blambda\T\bbeta=0 \) with \( \blambda=(0,0,1,-1)\T \), and take
\( \bm P=\X(\X\T\X)^{-1}\blambda \), so that \( \bm P\T\X=\blambda\T \). Then \( \M\bm P=\bm P \).
Each column of \( \X_0 \) is \( \X\bm k \) with \( \blambda\T\bm k=0 \), so \( \bm P\perp\C(\X_0) \) and
\( \bm P \) spans the same line.
:::

### C. Going deeper

::: {#exr-proj-two-factor-interaction}
[C1]

Extend @exm-proj-two-factor to \( m\ge2 \) observations in every cell and a model with
interaction. Find five mutually orthogonal projections summing to \( \I \), corresponding to
the grand mean, factor \( A \), factor \( B \), the interaction, and within-cell error, and find
their ranks. Describe the action of each on \( \y \) in terms of cell, row, column and grand
means.
:::
