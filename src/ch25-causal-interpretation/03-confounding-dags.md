# Confounding and causal graphs

The functions and disturbances of a structural causal model are usually unknown; its parent sets, saying which
variables can influence which, are qualitative and best displayed as a graph. This section shows how the graph
constrains the joint distribution and identifies confounding as a property of paths. For linear models the connection is exact, and
the omitted-variable bias of @prp-lm-omitted becomes a sum over paths.

## Graphs

::: {#def-cau-dag}
[Directed acyclic graph]

A **directed graph** on nodes \( 1,\dots,m \) is a set of ordered pairs \( (l,j) \), \( l\ne j \), drawn as arrows
\( l\to j \). It is **acyclic** (a DAG) if no sequence of arrows leads from a node back to itself. In a DAG:

::: {.enumerate options="label=(\alph*)"}
1. \( l \) is a **parent** of \( j \), and \( j \) a **child** of \( l \), if \( l\to j \);

2. a **path** between \( l \) and \( j \) is a sequence of distinct nodes, starting at \( l \) and ending at \( j \), in which
   consecutive nodes are joined by an arrow in either direction;

3. a path is **directed** (from \( l \) to \( j \)) if every arrow on it points towards \( j \); then \( l \) is an
   **ancestor** of \( j \), and \( j \) a **descendant** of \( l \);

4. an interior node \( c \) of a path is a **collider** on the path if both arrows of the path that meet at \( c \) point
   into it, \( \to c\leftarrow \); every other interior node is a **non-collider** on the path.
:::

The **graph of a structural causal model** has an arrow \( l\to j \) iff \( l\in\mathrm{pa}(j) \). It is acyclic because
parents have smaller indices.
:::

Being a collider is a property of a node *on a given path*: in [Figure 25.3.1](#fig-cau-structures)(c), \( C \) is a
collider on \( X\to C\leftarrow Y \) but not on a path that continues along an arrow \( C\to D \).

## What the graph says about the distribution

Assume that for each \( j \) and each value \( \bv \) of the
parents, the variable \( f_j(\bv,U_j) \) has a density \( p_j(\cdot\mid\bv) \) with respect to Lebesgue or counting
measure.

::: {#prp-cau-markov}
[Markov factorization and truncated factorization]

Under this assumption:

::: {.enumerate options="label=(\alph*)"}
1. \( (V_1,\dots,V_m) \) has density
   \[
   p(v_1,\dots,v_m)=\prod_{j=1}^mp_j\bigl(v_j\mid v_{\mathrm{pa}(j)}\bigr),
   \]{#eq-cau-markov}

   and \( p_j(\cdot\mid v_{\mathrm{pa}(j)}) \) is the conditional density of \( V_j \) given \( V_1,\dots,V_{j-1} \), which
   depends on them only through the parents;

2. under \( \operatorname{do}(V_k=x) \), the other variables have density
   \( \prod_{j\ne k}p_j\bigl(v_j\mid v_{\mathrm{pa}(j)}\bigr) \), evaluated with \( v_k=x \).
:::

:::

::: {.proof}
(a) The variables \( V_1,\dots,V_{j-1} \) are functions of \( U_1,\dots,U_{j-1} \), which are independent of \( U_j \). For
independent inputs, the conditional distribution of \( f_j(V_{\mathrm{pa}(j)},U_j) \) given \( (V_1,\dots,V_{j-1})=(v_1,\dots,v_{j-1}) \)
is the distribution of \( f_j(v_{\mathrm{pa}(j)},U_j) \) (Billingsley 1995), whose density is \( p_j(\cdot\mid v_{\mathrm{pa}(j)}) \). The
factorization follows by multiplying these conditional densities in the order \( j=1,\dots,m \).
(b) The modified model is a structural causal model for the variables other than \( V_k \), in which \( V_k \) is the constant
\( x \) and every other assignment and disturbance is unchanged. Apply (a) to it.
:::

So an intervention deletes one factor and keeps the others, which are observational conditional densities; with every
variable observed, the interventional distribution is computable from the observational one.

The linear case has an explicit form. Write a linear model as \( \mathbf{V}=\mathbf{c}+\B\mathbf{V}+\mathbf{U} \), where
\( \B=(b_{jl}) \) has \( b_{jl}\ne0 \) only if \( l\in\mathrm{pa}(j) \), so \( \B \) is strictly lower triangular, and
\( \boldsymbol{\Omega}=\Cov(\mathbf{U})=\diag(\omega_1,\dots,\omega_m) \). The coefficient \( b_{jl} \) labels the arrow
\( l\to j \).

::: {#prp-cau-path-tracing}
[Path tracing]

In a linear structural causal model:

::: {.enumerate options="label=(\alph*)"}
1. \( \I-\B \) is invertible, \( \mathbf{V}=\A(\mathbf{c}+\mathbf{U}) \) with \( \A=(\I-\B)^{-1} \), and
   \( \Cov(\mathbf{V})=\A\boldsymbol{\Omega}\A\T \);

2. for \( j\ne l \), the entry \( a_{jl} \) of \( \A \) is the sum, over all directed paths from \( l \) to \( j \), of the products of
   the coefficients on the arrows of the path (zero if there is none), and \( a_{jj}=1 \);

3. \( \E\bigl(V_j\mid\operatorname{do}(V_l=x)\bigr)=\E\bigl(V_j\mid\operatorname{do}(V_l=0)\bigr)+a_{jl}\,x \).
:::

The number \( a_{jl} \) is the **total effect** of \( V_l \) on \( V_j \): the change in the interventional mean of \( V_j \) per
unit change in the value imposed on \( V_l \).
:::

::: {.proof}
(a) \( \I-\B \) is unit lower triangular, hence invertible, and \( (\I-\B)\mathbf{V}=\mathbf{c}+\mathbf{U} \). The covariance
follows from @thm-rv-linear.
(b) \( \B \) is strictly lower triangular, so \( \B^m=\bzero \) and \( (\I-\B)(\I+\B+\dots+\B^{m-1})=\I-\B^m=\I \). The
\( (j,l) \) entry of \( \B^r \) is \( \sum b_{jk_1}b_{k_1k_2}\cdots b_{k_{r-1}l} \) over all \( k_1,\dots,k_{r-1} \); a term is
nonzero only if \( l\to k_{r-1}\to\dots\to k_1\to j \) are arrows, which is a directed path of length \( r \) (its nodes are
distinct because the graph is acyclic). Summing over \( r \) gives (b).
(c) The intervened model is linear with \( \B \) replaced by \( \tilde{\B} \), in which row \( l \) is zero, and with \( c_l+U_l \)
replaced by \( x \). A directed path out of \( l \) never uses an arrow into \( l \), so the directed paths from \( l \) in the
graph of \( \tilde{\B} \) are those of \( \B \), and by (b) the \( (j,l) \) entry of \( (\I-\tilde{\B})^{-1} \) is \( a_{jl} \). By (a)
applied to the intervened model, \( V_j \) is \( a_{jl}x \) plus terms not involving \( x \).
:::

This is Wright's method of path coefficients (Wright 1921, 1934); the covariance formula in (a) gives his rule for
correlations (@exr-cau-trek).

## Chains, forks and colliders

Three configurations of three variables are the building blocks of everything the graph says about conditional
independence
([Figure 25.3.1](#fig-cau-structures)).

\begin{center}
\begin{tikzpicture}[>=Stealth, line cap=round]
\definecolor{ink}{HTML}{1B1F27}\definecolor{accent}{HTML}{2F5F96}\definecolor{draft}{HTML}{A45C25}
\tikzset{var/.style={draw=ink, circle, minimum size=7mm, inner sep=0pt, font=\small}}
\node[var] (X1) at (0,0) {$X$}; \node[var] (M1) at (1.3,0) {$M$}; \node[var] (Y1) at (2.6,0) {$Y$};
\draw[->, thick, accent] (X1) -- node[above, font=\small] {$a$} (M1);
\draw[->, thick, accent] (M1) -- node[above, font=\small] {$b$} (Y1);
\node[font=\small] at (1.3,-0.9) {(a) chain};
\node[var] (X2) at (3.8,0) {$X$}; \node[var] (Z2) at (5.1,0.9) {$Z$}; \node[var] (Y2) at (6.4,0) {$Y$};
\draw[->, thick, draft] (Z2) -- node[above left, font=\small] {$a$} (X2);
\draw[->, thick, draft] (Z2) -- node[above right, font=\small] {$b$} (Y2);
\node[font=\small] at (5.1,-0.9) {(b) fork};
\node[var] (X3) at (7.6,0.9) {$X$}; \node[var] (C3) at (8.9,0) {$C$}; \node[var] (Y3) at (10.2,0.9) {$Y$};
\draw[->, thick, ink] (X3) -- node[below left, font=\small] {$a$} (C3);
\draw[->, thick, ink] (Y3) -- node[below right, font=\small] {$b$} (C3);
\node[font=\small] at (8.9,-0.9) {(c) collider};
\end{tikzpicture}
\end{center}

[**Figure 25.3.1.** The three structures of @prp-cau-three-structures. In (a) and (b) the middle node transmits
association between \( X \) and \( Y \), and conditioning on it removes the association. In (c) the middle node transmits
nothing until it is conditioned on.]{#fig-cau-structures}

::: {#prp-cau-three-structures}
[Chains, forks and colliders]

Consider linear models with independent disturbances, \( \Var U_X=\sigma_X^2>0 \) and so on, and write
\( \sigma_{XY\cdot W} \) for the partial covariance of \( X \) and \( Y \) given \( W \) (@def-mvn-partial-correlation).

::: {.enumerate options="label=(\alph*)"}
1. **Chain** \( X=U_X \), \( M=aX+U_M \), \( Y=bM+U_Y \): \( \Cov(X,Y)=ab\sigma_X^2 \) and \( \sigma_{XY\cdot M}=0 \).

2. **Fork** \( Z=U_Z \), \( X=aZ+U_X \), \( Y=bZ+U_Y \): \( \Cov(X,Y)=ab\sigma_Z^2 \) and \( \sigma_{XY\cdot Z}=0 \).

3. **Collider** \( X=U_X \), \( Y=U_Y \), \( C=aX+bY+U_C \): \( \Cov(X,Y)=0 \) and
   \[
   \sigma_{XY\cdot C}=-\frac{ab\,\sigma_X^2\sigma_Y^2}{a^2\sigma_X^2+b^2\sigma_Y^2+\sigma_C^2}.
   \]{#eq-cau-collider-partial}

:::

With normal disturbances, a zero partial covariance means conditional independence (@prp-mvn-partial-meaning), so in
(a) and (b) \( X \) and \( Y \) are dependent but conditionally independent given the middle node, and in (c) they are
independent but, if \( ab\ne0 \), conditionally dependent given \( C \).
:::

::: {.proof}
The covariances follow from @prp-cau-path-tracing(a). (a) \( Y-bM=U_Y \) is uncorrelated with \( X \) and \( M \), so
\( L(Y\mid M,X)=bM \) has no \( X \) term, and by @prp-cor-partial-coefficient \( \sigma_{XY\cdot M}=0 \). (b) The same argument
with \( Y-bZ=U_Y \). (c) By definition,
\[
\sigma_{XY\cdot C}=\Cov(X,Y)-\frac{\Cov(X,C)\Cov(Y,C)}{\Var(C)}=0-\frac{(a\sigma_X^2)(b\sigma_Y^2)}{\Var(C)}.
\]
:::

In the collider, independent causes become dependent once their common effect is held fixed. In the listing, with \( a=b=1 \) and unit variances, the correlations of \( X \) and \( Y \) in the chain and the fork
are \( 0.577 \) and \( 0.500 \) (simulated: \( 0.582 \) and
\( 0.505 \)), and both partial correlations given the middle node are zero. In the collider the
correlation is zero and the partial correlation given \( C \) is \( -0.500 \) (simulated
\( -0.501 \)).

```{.python .run #cell-structures-structures}
import numpy as np

rng = np.random.default_rng(2503)
n = 100_000
a, b = 1.0, 1.0                                   # the two edge coefficients

def corr_and_partial(x, y, w):
    """Correlation of x and y, and their partial correlation given w (residuals on 1, w)."""
    W = np.column_stack([np.ones(len(w)), w])
    rx = x - W @ np.linalg.lstsq(W, x, rcond=None)[0]
    ry = y - W @ np.linalg.lstsq(W, y, rcond=None)[0]
    return np.corrcoef(x, y)[0, 1], np.corrcoef(rx, ry)[0, 1]

U = rng.normal(size=(3, n))
x = U[0]; m = a * x + U[1]; y = b * m + U[2]      # chain   X -> M -> Y
chain = corr_and_partial(x, y, m)
z = U[0]; x = a * z + U[1]; y = b * z + U[2]      # fork    X <- Z -> Y
fork = corr_and_partial(x, y, z)
x = U[0]; y = U[1]; c = a * x + b * y + U[2]      # collider X -> C <- Y
collider = corr_and_partial(x, y, c)
for name, (r, rp) in [("chain", chain), ("fork", fork), ("collider", collider)]:
    print(f"{name:9s} corr {r:6.3f}   partial corr given middle {rp:6.3f}")
```

## d-separation

In a larger graph association can travel along many paths at once; d-separation combines the three structures.

::: {#def-cau-d-separation}
[d-separation]

Let \( \mathcal Z \) be a set of nodes. A path is **blocked** by \( \mathcal Z \) if

::: {.enumerate options="label=(\roman*)"}
1. some non-collider on the path is in \( \mathcal Z \), or

2. some collider on the path is not in \( \mathcal Z \) and has no descendant in \( \mathcal Z \).
:::

Otherwise the path is **open** given \( \mathcal Z \). Disjoint sets of nodes \( \mathcal A \) and \( \mathcal B \) are
**d-separated** by \( \mathcal Z \) (disjoint from both) if every path between a node of \( \mathcal A \) and a node of
\( \mathcal B \) is blocked by \( \mathcal Z \).
:::

A path is open when each non-collider is left free and each collider is conditioned on, directly or through a
descendant, a noisy measurement of it (@exr-cau-collider-descendant).

::: {#thm-cau-d-separation}
[d-separation and conditional independence]

Let \( \mathbf{V} \) follow a structural causal model with graph \( G \), and let \( \mathcal A,\mathcal B,\mathcal Z \) be disjoint
sets of nodes.

::: {.enumerate options="label=(\alph*)"}
1. If \( \mathcal Z \) d-separates \( \mathcal A \) and \( \mathcal B \) in \( G \), then \( \mathbf{V}_{\mathcal A} \) and
   \( \mathbf{V}_{\mathcal B} \) are conditionally independent given \( \mathbf{V}_{\mathcal Z} \).

2. If \( \mathcal Z \) does not d-separate them, there is a linear structural causal model on \( G \) with normal disturbances
   in which they are conditionally dependent; in fact this holds for all edge coefficients and disturbance variances
   outside a set of Lebesgue measure zero.
:::

:::

::: {.proof}
*Sketch; full proofs are in the references.* For (a) with densities (Lauritzen, Dawid, Larsen and Leimer 1990): the
nodes of \( \mathcal A\cup\mathcal B\cup\mathcal Z \) and their ancestors form a set \( \mathcal N \) whose density, by
@prp-cau-markov, is a product of functions of the sets \( \{j\}\cup\mathrm{pa}(j) \), \( j\in\mathcal N \). These sets are
complete in the *moral graph* on \( \mathcal N \), which joins each node to its parents and any two parents of a common child.
If \( \mathcal Z \) separates \( \mathcal A \) from \( \mathcal B \) there, the factors split into two groups, giving the conditional
independence; a combinatorial lemma shows that d-separation in \( G \) is equivalent to this separation. Geiger, Verma and Pearl (1990) give
another proof. For (b), partial correlations of a linear normal model are rational functions of the parameters, not
identically zero along an open path, and so vanish only on a set of measure zero (Spirtes, Glymour and Scheines 2000).
:::

Part (b) says the converse of (a) fails only by accident, as in @exr-cau-zero-slope, where the effect of visits and the
association through illness cancel exactly. Distributions with no
such cancellation are called **faithful** to the graph.

## Confounding as an open path

::: {#def-cau-confounding}
[Back-door paths and confounding]

A **back-door path** from \( X \) to \( Y \) is a path between them whose arrow at \( X \) points into \( X \),
\( X\leftarrow\cdots Y \). The effect of \( X \) on \( Y \) is **confounded** if some back-door path from \( X \) to \( Y \) is open
given the empty set, that is, contains no collider.
:::

Directed paths from \( X \) to \( Y \) carry the effect of \( X \); back-door paths carry association that reflects whatever
caused \( X \). In a linear model the split is exact.

::: {#prp-cau-confounding-bias}
[Confounding bias in a linear model]

In a linear structural causal model with disturbances of positive variance, let \( X=V_k \), \( Y=V_j \), and let
\( \tau=a_{jk} \) be the total effect of \( X \) on \( Y \). For each \( l \), let \( d_l \) be the sum over the directed paths from
\( V_l \) to \( Y \) that do not pass through \( X \) of the products of coefficients (so \( d_k=0 \), and \( d_j=1 \) from the trivial path of
\( Y \) to itself).

::: {.enumerate options="label=(\alph*)"}
1. \( Y=\tau X+R \), where \( R=\sum_ld_l(c_l+U_l) \) is the value of \( Y \) under \( \operatorname{do}(X=0) \). Under
   \( \operatorname{do}(X=x) \), \( Y=\tau x+R \).

2. The slope of the best linear predictor \( L(Y\mid X) \) is
   \[
   \tau+\frac{\Cov(R,X)}{\Var(X)},\qquad \Cov(R,X)=\sum_l\omega_l\,d_l\,a_{kl}.
   \]{#eq-cau-backdoor-bias}

3. If the effect of \( X \) on \( Y \) is not confounded, the slope of \( L(Y\mid X) \) is \( \tau \).
:::

:::

::: {.proof}
(a) A directed path from \( V_l \) to \( Y \) either avoids \( X \) or passes through it exactly once, and a path of the second
kind splits at \( X \) into a directed path from \( V_l \) to \( X \) and one from \( X \) to \( Y \). Summing products over paths
with @prp-cau-path-tracing(b), \( a_{jl}=d_l+a_{kl}\tau \) for every \( l \) (for \( l=k \), \( a_{kk}=1 \) and \( d_k=0 \); for \( l=j \), \( a_{jj}=1=d_j \) and \( a_{kj}=0 \)). Hence
\( Y=\sum_la_{jl}(c_l+U_l)=\tau\sum_la_{kl}(c_l+U_l)+R=\tau X+R \). Under \( \operatorname{do}(X=x) \), as in the proof of
@prp-cau-path-tracing(c), the arrows into \( X \) are removed, so the directed paths from \( V_l \), \( l\ne k \), to \( Y \) in
the intervened graph are exactly those that avoid \( X \). Thus \( Y=\tau x+\sum_{l\ne k}d_l(c_l+U_l)=\tau x+R \).

(b) The slope is \( \Cov(X,Y)/\Var(X)=\tau+\Cov(R,X)/\Var(X) \) by @thm-proj-blp, and
\( \Cov(R,X)=\sum_{l,l'}d_la_{kl'}\Cov(U_l,U_{l'})=\sum_l\omega_ld_la_{kl} \), since the disturbances are uncorrelated.

(c) Suppose \( \omega_ld_la_{kl}\ne0 \) for some \( l \). Then \( l\ne k \), and there are a directed path \( \pi_1 \) from \( V_l \) to
\( X \) and a directed path \( \pi_2 \) from \( V_l \) to \( Y \) that avoids \( X \). Let \( w \) be the node of \( \pi_1 \) closest to \( X \)
that also lies on \( \pi_2 \); it exists because \( V_l \) lies on both, and \( w\ne X \) because \( \pi_2 \) avoids \( X \). Follow
\( \pi_1 \) backwards from \( X \) to \( w \), then \( \pi_2 \) forwards from \( w \) to \( Y \). By the choice of \( w \) the nodes are
distinct, so this is a path. Its arrow at \( X \) is the last arrow of \( \pi_1 \), which points into \( X \), and every arrow
on it points away from \( w \), so it has no collider. The effect is therefore confounded. Contrapositively, if it is not
confounded, every term of \( \Cov(R,X) \) vanishes.
:::

Part (a) is the potential outcome \( Y(x)=\tau x+Y(0) \): linear structural models have constant effects. Part (b) writes
the regression slope as the effect plus one contribution for each common ancestor of \( X \) and \( Y \). In the model of
[Figure 25.1.1](01-prediction-intervention.html#fig-cau-confounder-dag), \( R=c_Y+\gamma Z+U_Y \) and the only term comes from
\( Z \): \( \omega d a=\sigma_Z^2\cdot\gamma\cdot\alpha \), which is @eq-cau-confounded-slope.

This is the omitted-variable bias of least squares: by @prp-lm-omitted the short regression estimates \( \theta+\pi\gamma \),
with \( \pi \) the sample slope of \( \mathbf{z} \) on \( \x \), whose limit is \( \alpha\sigma_Z^2/\Var(X) \)
(see also @thm-dep-omitted in [Chapter 19](../ch19-theory-of-departures/index.html)). The graph adds that an omitted variable
causes bias only if it lies on an open back-door path; an omitted mediator should stay omitted
([Section 25.5](05-bad-controls.html)).

## Exercises

### A. Check your understanding

::: {#exr-cau-paths}
[A1]

A DAG has arrows \( A\to X \), \( A\to B \), \( B\to Y \), \( X\to Y \), \( X\to C \), \( Y\to C \) and \( D\to C \). List every path between
\( X \) and \( Y \), mark the back-door paths, and say which paths are open given \( \emptyset \), given \( \{B\} \) and given
\( \{C\} \).
:::

::: {.solution}
There are three paths: \( X\to Y \); \( X\leftarrow A\to B\to Y \), a back-door path; and \( X\to C\leftarrow Y \), with collider
\( C \). (The arrow \( D\to C \) gives no path to \( Y \) except through \( C \), where \( D\to C\leftarrow Y \) would require revisiting
\( C \).) Given \( \emptyset \): the first two are open and the third is blocked by the collider. Given \( \{B\} \): only
\( X\to Y \) is open. Given \( \{C\} \): all three are open, since conditioning on the collider opens the third.
:::

### B. Practice

::: {#exr-cau-collider-descendant}
[B1]

Extend the collider by \( D=C+U_D \), with \( \Var U_D=\delta \) and all other variances one. Compute \( \sigma_{XY\cdot D} \) and
show that it lies between \( 0 \) and \( \sigma_{XY\cdot C} \), approaching them as \( \delta\to\infty \) and \( \delta\to0 \).
:::

::: {.solution}
\( \Cov(X,D)=a \), \( \Cov(Y,D)=b \) and \( \Var D=a^2+b^2+1+\delta \), so
\( \sigma_{XY\cdot D}=-ab/(a^2+b^2+1+\delta) \). At \( \delta=0 \) this is @eq-cau-collider-partial, and it decreases in absolute
value to zero as \( \delta \) grows. A proxy for a collider opens the path partially, the more so the better the proxy.
:::

### C. Going deeper

::: {#exr-cau-trek}
[C1]

Show from @prp-cau-path-tracing(a) that in a linear model \( \Cov(V_j,V_k)=\sum_l\omega_la_{jl}a_{kl} \). Interpret each term
as a sum over pairs of directed paths from \( V_l \) to \( V_j \) and to \( V_k \), and use this to compute \( \Cov(X,Y) \) in the chain
of @prp-cau-three-structures.
:::
