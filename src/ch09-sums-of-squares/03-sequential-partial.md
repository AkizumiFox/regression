# Sequential and partial sums of squares

How should the variation explained by a model be shared among its terms? A term is a block of
columns representing one idea: a regressor, the coding columns of a factor, or an interaction.
One answer enters the terms one at a time and credits each with what it adds to those already
present. Another credits each term with what it adds when entered last. This section defines
both, explains the three “types” of sums of squares that software prints, and identifies the
hypothesis that each of them addresses.

## Sequential sums of squares

Write the model matrix as a sequence of terms,
\[
\X=[\X_0,\X_1,\dots,\X_k],
\]
where \( \X_0 \) is usually the intercept column \( \bone \). Let \( \M_j \) be the projection onto
\( \C(\X_0,\X_1,\dots,\X_j) \), so that \( \M_k=\M \).

::: {#def-ss-sequential}
[Sequential sums of squares]

The **sequential** (or **Type I**) sum of squares of the term \( \X_j \), for the order
\( \X_0,\X_1,\dots,\X_k \), is
\[
\text{SS}(\X_j\mid\X_0,\dots,\X_{j-1})=\y\T(\M_j-\M_{j-1})\y ,\qquad j=1,\dots,k,
\]
with \( d_j=\rank(\M_j)-\rank(\M_{j-1}) \) degrees of freedom.
:::

The sequential sums of squares are the analysis of variance decomposition of
@thm-ss-decomposition(e) for the chain
\( \C(\X_0)\subseteq\C(\X_0,\X_1)\subseteq\dots\subseteq\C(\X) \), after the first piece
\( \M_0 \) has been set aside. They add up to the regression sum of squares
\( \text{SS}(\X_1,\dots,\X_k\mid\X_0) \), their degrees of freedom add up to \( r-\rank(\X_0) \), and
under normality they are independent of each other and of the residual sum of squares
(@thm-qf-orthogonal-projections). A term linearly dependent on earlier ones receives zero degrees
of freedom. The price is that each piece depends on the order: \( \X_j \) is credited only with the
part of its span orthogonal to *the terms that happen to precede it*. For the murder rates of
[Section 6.5](../ch06-projections/05-nested.html), poverty received \( 102.49 \) first and
\( 16.37 \) last (@exm-proj-sequential-order).

The order dependence can be more surprising than “entering first gives more credit”. The
listing computes \( \text{SS}(\x\mid\bone,S) \) for each of the three regressors and every set
\( S \) of the other two.

| Regressor | given \( \bone \) only | given \( \bone \) and one other (named) | given \( \bone \) and the remaining other (named) | given \( \bone \) and both others |
|---|---|---|---|---|
| poverty | 102.49 | 19.69 (single) | 125.45 (urban) | 16.37 |
| single | 162.84 | 80.04 (poverty) | 160.13 (urban) | 51.06 |
| urban | 6.30 | 29.26 (poverty) | 3.59 (single) | 0.27 |

Adjusting for urbanization *increases* the sum of squares of poverty from \( 102.49 \) to
\( 125.45 \), and adjusting for poverty more than quadruples that of urbanization. Urban states have somewhat less poverty but more murder, so holding one fixed sharpens the
association of the other; such a regressor is called a **suppressor** (@exr-ss-two-regressor-formula).
A sequential sum of squares is not a property of a regressor, and not even a monotone function of
the set adjusted for.

```{.python .run #cell-extra-ss-subsets}
from itertools import combinations
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
y = data["murder"].to_numpy()
n = len(y)
one = np.ones(n)
cols = {"poverty": data["poverty"].to_numpy(), "single": data["single"].to_numpy(),
        "urban": data["urban"].to_numpy()}

def sse(Z):
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    e = y - Z @ b
    return e @ e

table = {}
for name, x in cols.items():
    others = [k for k in cols if k != name]
    for size in range(3):
        for S in combinations(others, size):
            base = np.column_stack([one] + [cols[k] for k in S])
            table[name, S] = sse(base) - sse(np.column_stack([base, x]))
            print(f"SS({name} | 1{''.join(', ' + k for k in S)}) = {table[name, S]:7.2f}")
```

## Partial sums of squares: Types II and III

The second answer credits a term with what it adds last. For models whose terms are factors
and their interactions, “last” needs care, because a main effect is part of every interaction
that contains it. Say that a term \( \X_i \) **contains** a term \( \X_j \) if the factors of
\( \X_j \) are a proper subset of the factors of \( \X_i \). For example the interaction \( AB \)
contains the main effects \( A \) and \( B \).

::: {#def-ss-partial}
[Partial sums of squares; Types II and III]

For a model with terms \( \X_0,\X_1,\dots,\X_k \), write \( \X_{(-j)} \) for \( \X \) with the columns of
\( \X_j \) removed.

::: {.enumerate options="label=(\alph*)"}
1. The **partial** (or **last-in**) sum of squares of \( \X_j \) is
   \( \text{SS}(\X_j\mid\X_{(-j)})=\y\T(\M-\M_{(-j)})\y \), where \( \M_{(-j)} \) projects onto
   \( \C(\X_{(-j)}) \).

2. The **Type II** sum of squares of \( \X_j \) is its extra sum of squares given all terms that do
   not contain it: \( \text{SS}(\X_j\mid\X_i:\ \X_i\text{ does not contain }\X_j,\ i\neq j) \).

3. The **Type III** sum of squares of \( \X_j \) is its partial sum of squares (a) computed when
   every factor is coded by sum-to-zero contrast columns and every interaction by the products of
   the columns of its factors.
:::

:::

When no term contains another, as in a regression on numerical regressors without
interactions, the three coincide, and they equal the sequential sum of squares of \( \X_j \) when
it is entered last. They also coincide for a term that is not contained in any other, such as
the highest-order interaction. They differ for main effects in models with interactions.

Partial sums of squares are not a decomposition. The subspaces \( \C(\X)\ominus\C(\X_{(-j)}) \) are
generally not orthogonal, so the partial sums of squares can add up to more or less than the
regression sum of squares (@exr-ss-partial-sum-sign), and they need not be independent of each
other, though each is independent of the residual sum of squares under normality.

Why does the definition of Type III name a coding? Because in models with interactions the
partial sum of squares of a main effect depends on the coding. Deleting the columns of \( A \) while
keeping those of \( AB \) leaves the reduced space \( \C(\bone,\X_B,\X_{AB}) \), which depends on
which columns represent \( AB \), although the full space does not. Different codings give
different reduced models and so test different hypotheses.

## What each sum of squares tests

A sum of squares \( \y\T\bP\y \) with \( \bP \) a projection of rank \( d \) has expectation
\[
\E(\Y\T\bP\Y)=\sigma^2d+\norm{\bP\bmu}^2,\qquad \bmu=\E(\Y),
\]
by @thm-rv-quadform-mean with \( \Cov(\Y)=\sigma^2\I \). So the hypothesis a sum of squares addresses
is \( \bP\bmu=\bzero \): *the mean has no component in the subspace onto which it projects*. Under
normality this is exactly the condition for its \( F \) ratio to be central
([Section 9.6](06-expected-mean-squares.html)). We translate it into statements about parameters.

::: {#prp-ss-tested-regression}
[Hypotheses of sequential and partial sums of squares]

Let \( \X=[\X_0,\dots,\X_k] \) have full column rank and \( \bmu=\sum_l\X_l\bbeta_l \).

::: {.enumerate options="label=(\alph*)"}
1. The partial sum of squares of \( \X_j \) addresses \( \bbeta_j=\bzero \):
   \( (\M-\M_{(-j)})\bmu=\bzero \) iff \( \bbeta_j=\bzero \).

2. The sequential sum of squares of \( \X_j \) addresses
   \[
\bbeta_j+(\tilde{\X}_j\T\tilde{\X}_j)^{-1}\tilde{\X}_j\T\sum_{l>j}\X_l\bbeta_l=\bzero,
   \qquad \tilde{\X}_j=(\I-\M_{j-1})\X_j .
\]
   This is \( \bbeta_j=\bzero \) whenever the later terms have zero coefficients, or whenever
   \( \tilde{\X}_j\T\X_l=\bzero \) for all \( l>j \).
:::

:::

::: {.proof}
(a) For \( l\neq j \), \( \X_l\bbeta_l \) lies in \( \C(\X_{(-j)})\subseteq\C(\X) \), so
\( (\M-\M_{(-j)})\X_l\bbeta_l=\bzero \). Let \( \check{\X}_j=(\I-\M_{(-j)})\X_j \). By
@lem-ss-residualized, \( \M-\M_{(-j)} \) projects onto \( \C(\check{\X}_j) \), so
\( (\M-\M_{(-j)})\X_j\bbeta_j \) is the projection of \( \X_j\bbeta_j \) onto that space. Since
\( \X_j\bbeta_j-\check{\X}_j\bbeta_j=\M_{(-j)}\X_j\bbeta_j \) is orthogonal to \( \C(\check{\X}_j) \), the
projection is \( \check{\X}_j\bbeta_j \). The columns of \( \check{\X}_j \) are independent
(@lem-proj-fwl-split(b)), so it vanishes iff \( \bbeta_j=\bzero \).

(b) Terms with \( l<j \) are annihilated by \( \M_j-\M_{j-1} \), which projects onto
\( \C(\tilde{\X}_j) \) by @lem-ss-residualized. As in (a), \( \X_j\bbeta_j \) projects to
\( \tilde{\X}_j\bbeta_j \). A later term projects to
\( \tilde{\X}_j(\tilde{\X}_j\T\tilde{\X}_j)^{-1}\tilde{\X}_j\T\X_l\bbeta_l \) by @thm-proj-M-formula.
Adding, \( (\M_j-\M_{j-1})\bmu=\tilde{\X}_j\bm v \), where \( \bm v \) is the vector in the display, and
this is zero iff \( \bm v=\bzero \) because \( \tilde{\X}_j \) has independent columns.
:::

So a sequential sum of squares answers the question “does \( \X_j \) matter?” only if the terms
fitted after it do not matter, or are orthogonal to it after adjusting for the terms before it.
Otherwise the effects of later terms leak into it, through the coefficients
\( (\tilde{\X}_j\T\tilde{\X}_j)^{-1}\tilde{\X}_j\T\X_l \) of the regressions of later terms on
\( \tilde{\X}_j \).

For factors, the natural parameters are the cell means. Consider a two-way layout with factor
\( A \) at \( a \) levels and \( B \) at \( b \) levels, \( n_{ij}\ge1 \) observations in cell \( (i,j) \), row
totals \( n_{i\cdot} \), and cell means \( \mu_{ij} \). The full model with interaction has model space
\( \C(\X)= \) all vectors that are constant within cells, of dimension \( ab \). Define the
**weighted** and **unweighted** row means
\[
\bar{\mu}^{w}_{i}=\frac{1}{n_{i\cdot}}\sum_jn_{ij}\mu_{ij},
\qquad
\bar{\mu}_{i\cdot}=\frac1b\sum_j\mu_{ij}.
\]

::: {#thm-ss-two-way-hypotheses}
[What the sums of squares for a main effect test]

In the two-way layout with all \( n_{ij}\ge1 \):

::: {.enumerate options="label=(\alph*)"}
1. the sequential sum of squares \( \text{SS}(A\mid\mu) \) of \( A \) entered first addresses
   \( \bar{\mu}^w_1=\dots=\bar{\mu}^w_a \);

2. the Type III sum of squares of \( A \) addresses
   \( \bar{\mu}_{1\cdot}=\dots=\bar{\mu}_{a\cdot} \);

3. the partial sum of squares of the columns of \( A \) under reference-level coding (indicators of
   levels \( 2,\dots,a \) of \( A \) and \( 2,\dots,b \) of \( B \), and their products for \( AB \))
   addresses \( \mu_{11}=\mu_{21}=\dots=\mu_{a1} \), the equality of the \( A \) means at the first
   level of \( B \);

4. the Type II sum of squares \( \text{SS}(A\mid\mu,B) \) addresses a hypothesis that, when the
   cell means are additive (\( \mu_{ij}=\alpha_i+\beta_j \)), is equivalent to
   \( \alpha_1=\dots=\alpha_a \).
:::

Each sum of squares equals the sum of squares @eq-ss-hypothesis of the hypothesis it addresses,
in the cell-means model.
:::

::: {.proof}
Let \( \bP_0 \), \( \bP_A \), \( \bP_B \) and \( \bP_{A+B} \) be the projections onto the constant
vectors, the vectors constant within rows, those constant within columns, and their sum. Each
hypothesis below defines a subspace \( \mathcal S_0 \) of \( \C(\X) \), and each sum of squares is
\( \norm{\bP\y}^2 \) with \( \bP \) the projection onto \( \C(\X)\ominus\mathcal S_0 \). Then
\( \bP\bmu=\bzero \) iff \( \bmu\in\mathcal S_0 \), and the sum of squares equals \( \text{SSE}_0-\text{SSE} \)
for \( \mathcal S_0 \), which is the sum of squares of the hypothesis by @thm-ss-hypothesis.

(a) \( \bP_A \) replaces each entry of a vector by the average of its row, so \( \bP_A\bmu \) has entry
\( \bar{\mu}^w_i \) in row \( i \), and \( \bP_0\bmu \) has every entry equal to the overall weighted mean
\( \sum_in_{i\cdot}\bar{\mu}^w_i/n \). So \( (\bP_A-\bP_0)\bmu=\bzero \) iff all \( \bar{\mu}^w_i \) equal
their weighted average, that is, iff they are all equal. The relevant \( \mathcal S_0 \) is
\( \{\bmu\in\C(\X):(\bP_A-\bP_0)\bmu=\bzero\} \), whose orthogonal complement in \( \C(\X) \) is
\( \C(\bP_A-\bP_0) \), because \( \C(\bP_A-\bP_0)\subseteq\C(\X) \).

(b) Any \( a\times b \) array of cell means decomposes uniquely as
\[
\mu_{ij}=\bar{\mu}_{\cdot\cdot}+(\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot\cdot})+(\bar{\mu}_{\cdot j}-\bar{\mu}_{\cdot\cdot})
+(\mu_{ij}-\bar{\mu}_{i\cdot}-\bar{\mu}_{\cdot j}+\bar{\mu}_{\cdot\cdot}),
\]
where bars denote unweighted averages. The four arrays lie in the mutually orthogonal
subspaces (for the unweighted inner product on arrays) of constants, row effects summing to
zero, column effects summing to zero, and arrays with all row and column sums zero, of
dimensions \( 1,a-1,b-1,(a-1)(b-1) \). Sum-to-zero coding spans exactly these: the columns
\( \bm u_i-\bm u_1 \) (\( i=2,\dots,a \)), where \( \bm u_i \) indicates row \( i \), give all row effects
\( \alpha_i \) with \( \sum_i\alpha_i=0 \); similarly for \( B \); and the \( (a-1)(b-1) \) products are
linearly independent arrays with zero row and column sums, so they span that subspace. Deleting
the columns of \( A \) removes exactly the row-effect component, so the reduced model is
\( \mathcal S_0=\{\bmu:\bar{\mu}_{i\cdot}=\bar{\mu}_{\cdot\cdot}\text{ for all }i\} \).

(c) With reference coding, every cell-mean array is uniquely
\( \mu_{ij}=m+\alpha_i+\beta_j+\gamma_{ij} \) with \( \alpha_1=\beta_1=0 \) and
\( \gamma_{1j}=\gamma_{i1}=0 \); explicitly \( m=\mu_{11} \), \( \alpha_i=\mu_{i1}-\mu_{11} \),
\( \beta_j=\mu_{1j}-\mu_{11} \) and \( \gamma_{ij}=\mu_{ij}-\mu_{i1}-\mu_{1j}+\mu_{11} \). Deleting the
columns of \( A \) forces \( \alpha_i=0 \), that is, \( \mu_{i1}=\mu_{11} \) for all \( i \).

(d) The Type II sum of squares is \( \y\T(\bP_{A+B}-\bP_B)\y \). If \( \bmu \) is additive, then
\( \bmu\in\C(\bP_{A+B}) \) and \( (\bP_{A+B}-\bP_B)\bmu=\bmu-\bP_B\bmu \), which is zero iff \( \bmu \) is
constant within columns, that is, iff \( \alpha_i \) does not depend on \( i \).
:::

These are four different questions. Part (a) averages over \( B \) with the design's weights, so
if \( B \) matters and the rows have different mixes of \( B \) levels, it can detect “an effect of
\( A \)” that is entirely an effect of \( B \) (@exr-ss-type1-leak). Part (b) averages with equal
weights. Part (c) compares the rows at one level of \( B \): a
simple effect, not a main effect, and one that changes with the reference level. Part (d) is the
natural test of \( A \) without interaction; with interaction its hypothesis has weights that
depend on all the \( n_{ij} \).

::: {.warning}
[“Type III” with the default coding]

Many regression programs code factors by reference levels unless told otherwise. Asking them for
last-in sums of squares of main effects in a model with interactions then produces the
hypotheses of @thm-ss-two-way-hypotheses(c), which concern a single level of the other factor.
They are not Type III sums of squares, although they are sometimes reported under that name.
Use sum-to-zero coding, or better, write down the hypothesis and use @eq-ss-hypothesis.
:::

## An unbalanced example

::: {#exm-ss-anes-types}
[Party, education and political self-placement]

In the 1996 American National Election Study, \( 944 \) respondents placed themselves on a
seven-point scale from extremely liberal (1) to extremely conservative (7). Classify them by party
identification (Democrat, Independent, Republican, counting those who only lean to a party as Independents) and by
education (high school or less, college, graduate study). The cell counts and means are:

| | high school | college | graduate |
|---|---|---|---|
| Democrat | 146, 3.97 | 108, 3.43 | 126, 2.92 |
| Independent | 80, 4.39 | 64, 4.28 | 95, 4.09 |
| Republican | 87, 5.40 | 105, 5.36 | 133, 5.40 |

The design is unbalanced, and education is associated with party: \( 38\% \) of the Democrats but
only \( 27\% \) of the Republicans have at most a high school education. The pattern of means shows an
interaction. Among Democrats self-placement moves left with education by more than a scale
point, while among Republicans it does not move at all. The sums of squares for the terms are:

| Term | df | Type I, party first | Type I, education first | Type II | Type III | reference coding |
|---|---|---|---|---|---|---|
| party | 2 | 647.46 | 674.42 | 674.42 | 650.41 | 112.27 |
| education | 2 | 43.26 | 16.31 | 43.26 | 32.03 | 75.12 |
| interaction | 4 | 35.81 | 35.81 | 35.81 | 35.81 | 35.81 |
| residual | 935 | 1224.62 | 1224.62 | 1224.62 | 1224.62 | 1224.62 |

The two Type I columns each add up to the model sum of squares \( 726.54 \), up to rounding. The
other columns do not: the Type II entries add up to more. Education receives \( 16.31 \) when
entered first and \( 43.26 \) after party. Since the least educated group contains the most
Democrats, who are the most liberal, the marginal comparison mixes the leftward drift of educated
Democrats with the partisan composition of each education group, and the two partly cancel. The Type III sum of
squares \( 32.03 \) tests whether the unweighted averages over the three parties of the education
means are equal. The last column, computed with reference coding, answers different questions:
\( 75.12 \) for education measures the education gradient among *Democrats only* (the first party
level), which is the steepest, and \( 112.27 \) for party compares the parties among respondents
with a high school education only. The interaction row is the same in every column because the
interaction is the highest-order term.
:::

```{.python .run #cell-types-ss-data}
import numpy as np
import pandas as pd
import statsmodels.api as sm
d = sm.datasets.anes96.load_pandas().data
party = pd.cut(d["PID"], [-1, 1.5, 4.5, 6], labels=["Dem", "Ind", "Rep"]).rename("party")
educ = pd.cut(d["educ"], [0, 3.5, 5.5, 7], labels=["HS", "College", "Graduate"]).rename("educ")
y = d["selfLR"].to_numpy()
n = len(y)
print(pd.crosstab(party, educ))
print(d.groupby([party, educ], observed=True)["selfLR"].mean().unstack().round(3))
```

```{.python .run #cell-types-ss-projections}
def proj(*blocks):
    """Orthogonal projection onto the span of the given column blocks."""
    Z = np.column_stack(blocks)
    U, s, _ = np.linalg.svd(Z, full_matrices=False)
    U = U[:, s > 1e-10 * s[0]]                     # orthonormal basis of C(Z), any rank
    return U @ U.T

one = np.ones((n, 1))
A = pd.get_dummies(party).to_numpy(float)           # indicators of party
B = pd.get_dummies(educ).to_numpy(float)            # indicators of education
AB = np.column_stack([A[:, [i]] * B for i in range(3)])   # the nine cells
P0, PA, PB = proj(one), proj(one, A), proj(one, B)
PAB, M = proj(one, A, B), proj(AB)                  # additive model, cell-means model
I = np.eye(n)

def ss(P):
    return y @ P @ y

type1_AB = {"party": ss(PA - P0), "educ": ss(PAB - PA), "inter": ss(M - PAB)}
type1_BA = {"educ": ss(PB - P0), "party": ss(PAB - PB), "inter": ss(M - PAB)}
type2 = {"party": ss(PAB - PB), "educ": ss(PAB - PA), "inter": ss(M - PAB)}
sse = ss(I - M)
```

```{.python .run #cell-types-ss-type3}
def coded(labels, kind):
    """Columns for a factor: sum-to-zero ("sum") or reference-level ("treatment") coding."""
    D = pd.get_dummies(labels).to_numpy(float)
    return D[:, 1:] - D[:, [0]] if kind == "sum" else D[:, 1:]

type3 = {}
for kind in ["sum", "treatment"]:
    Ca, Cb = coded(party, kind), coded(educ, kind)
    Cab = np.column_stack([Ca[:, [i]] * Cb for i in range(2)])
    full = proj(one, Ca, Cb, Cab)                    # always the cell-means space
    type3[kind] = {"party": ss(full - proj(one, Cb, Cab)),
                   "educ": ss(full - proj(one, Ca, Cab))}
    assert np.allclose(full, M)
for name, table in [("Type I, party first", type1_AB), ("Type I, education first", type1_BA),
                    ("Type II", type2), ("Type III, sum coding", type3["sum"]),
                    ("'Type III', treatment coding", type3["treatment"])]:
    print(f"{name:30s}", {k: round(float(v), 2) for k, v in table.items()})
```

The script behind these cells also checks every entry against statsmodels, and it confirms
@thm-ss-two-way-hypotheses numerically: the party entries of the first, fourth and fifth
columns equal @eq-ss-hypothesis for the hypotheses of equal weighted row means, equal unweighted
row means, and equal means in the high-school column.

## Which one to use

There is no single right sum of squares for a main effect in an unbalanced design with
interaction, because there is no single main-effect hypothesis. Decide the hypothesis first, from
the four above, and compute its sum of squares with @eq-ss-hypothesis. If \( A \) and \( B \)
interact, a single “main effect of \( A \)” summarizes a heterogeneous comparison with arbitrary
weights; Nelder (1977) argued that such tests are rarely meaningful, and Type II sums of squares
respect this marginality by adjusting \( A \) for \( B \) but not for \( AB \). Without
interaction the Type II and Type III hypotheses coincide, and Type II, which does not adjust for
the interaction columns, generally gives the more powerful test. In regression without
interactions all partial types agree; the distinction matters only for terms contained in other
terms.

Chapter 17 treats empty cells, where even Type III hypotheses may be nonestimable
(@exr-ss-empty-cell).

## Exercises

### A. Check your understanding

::: {#exr-ss-list-types}
[A1]

For a two-way model with terms \( A \), \( B \) and \( AB \), write each of the Type I sums of squares
(both orders), the Type II sums of squares and the partial sums of squares in the notation
\( \text{SS}(\cdot\mid\cdot) \). Which of them are equal for every data set, whatever the design?
:::

::: {.solution}
Type I, \( A \) first: \( \text{SS}(A\mid\mu) \), \( \text{SS}(B\mid\mu,A) \), \( \text{SS}(AB\mid\mu,A,B) \).
Type I, \( B \) first: \( \text{SS}(B\mid\mu) \), \( \text{SS}(A\mid\mu,B) \), \( \text{SS}(AB\mid\mu,A,B) \).
Type II: \( \text{SS}(A\mid\mu,B) \), \( \text{SS}(B\mid\mu,A) \), \( \text{SS}(AB\mid\mu,A,B) \). Partial:
\( \text{SS}(A\mid\mu,B,AB) \), \( \text{SS}(B\mid\mu,A,AB) \), \( \text{SS}(AB\mid\mu,A,B) \), where the
first two depend on the coding. Always equal: the interaction entries of all four lists; the
Type II entry for \( A \) and the Type I entry for \( A \) with \( B \) first; the Type II entry for
\( B \) and the Type I entry for \( B \) with \( A \) first.
:::

### B. Practice

::: {#exr-ss-two-regressor-formula}
[B1]

Let \( \X=[\bone,\x_1,\x_2] \), and let \( r_{y1} \), \( r_{y2} \) and \( r_{12} \) be the sample correlations
among \( \y \), \( \x_1 \) and \( \x_2 \), with \( |r_{12}|<1 \). Show that
\[
\frac{\text{SS}(\x_2\mid\bone)}{\text{SST}}=r_{y2}^2,
\qquad
\frac{\text{SS}(\x_2\mid\bone,\x_1)}{\text{SST}}=\frac{(r_{y2}-r_{y1}r_{12})^2}{1-r_{12}^2}.
\]
Show that adjusting for \( \x_1 \) increases the sum of squares of \( \x_2 \) whenever \( r_{y2}=0 \) and
\( r_{y1}r_{12}\neq0 \), and find the condition for it in general.
:::

::: {.solution}
Work in \( \bone\perpc \) with the centred vectors \( \bu=\y-\bar{y}\bone \), \( \bv_1 \), \( \bv_2 \), scaled to
unit length, so that inner products are correlations and \( \text{SST}=\norm{\y-\bar{y}\bone}^2 \). By
@lem-ss-residualized, \( \text{SS}(\x_2\mid\bone)=\text{SST}\,(\bu\T\bv_2)^2=\text{SST}\,r_{y2}^2 \).
Adjusting for \( \x_1 \) replaces \( \bv_2 \) by \( \bw=\bv_2-r_{12}\bv_1 \), with \( \norm{\bw}^2=1-r_{12}^2 \)
and \( \bu\T\bw=r_{y2}-r_{y1}r_{12} \), so \( \text{SS}(\x_2\mid\bone,\x_1)=\text{SST}(\bu\T\bw)^2/\norm{\bw}^2 \).
If \( r_{y2}=0 \) the first ratio is zero and the second is \( r_{y1}^2r_{12}^2/(1-r_{12}^2)>0 \). In
general the sum of squares increases iff \( (r_{y2}-r_{y1}r_{12})^2>r_{y2}^2(1-r_{12}^2) \).
:::

::: {#exr-ss-weights-coincide}
[B2]

Show that the hypotheses of @thm-ss-two-way-hypotheses(a) and (b) are the same subspace of cell
means iff \( n_{ij}=n_{i\cdot}/b \) for all \( i,j \), that is, iff within every row all cells have
the same count. Show that under proportional frequencies, \( n_{ij}=n_{i\cdot}n_{\cdot j}/n \) (@thm-ss-proportional), the
hypothesis (a) compares row means weighted by the column proportions \( n_{\cdot j}/n \), the same
weights in every row.
:::

::: {.solution}
Write \( w_{ij}=n_{ij}/n_{i\cdot} \), so \( \sum_jw_{ij}=1 \), and let \( f_i(\bmu)=\sum_jw_{ij}\mu_{ij} \) and
\( g_i(\bmu)=\sum_j\mu_{ij}/b \). Hypothesis (a) is the null space of the linear functionals
\( \sum_ic_if_i \) with \( \sum_ic_i=0 \), and (b) that of the functionals \( \sum_id_ig_i \) with
\( \sum_id_i=0 \). Two subspaces defined as null spaces coincide iff the spaces of functionals
coincide. The functional \( f_i-f_1 \) has coefficient \( w_{ij} \) on cell \( (i,j) \) and \( -w_{1j} \) on
cell \( (1,j) \). If it equals some \( \sum_id_ig_i \), whose coefficients are constant across each row,
then \( w_{ij} \) does not depend on \( j \), so \( w_{ij}=1/b \); likewise \( w_{1j}=1/b \). Hence
coincidence forces \( n_{ij}=n_{i\cdot}/b \). Conversely, if all \( w_{ij}=1/b \), then \( f_i=g_i \). Under
proportional frequencies \( w_{ij}=n_{\cdot j}/n \), which does not depend on \( i \), so (a) compares
\( \sum_j(n_{\cdot j}/n)\mu_{ij} \) across rows.
:::

::: {#exr-ss-type1-leak}
[B3]

Take \( a=b=2 \), counts \( n_{11}=n_{22}=9 \) and \( n_{12}=n_{21}=1 \), and additive cell means
\( \mu_{ij}=\beta_j \) with \( \beta_1=0 \) and \( \beta_2=1 \), so that \( A \) has no effect. Compute the
weighted row means \( \bar{\mu}^w_i \) and show that the hypothesis addressed by the Type I sum of
squares of \( A \) entered first is false. What do the Type II and Type III hypotheses say here?
:::

::: {#exr-ss-empty-cell}
[B4]

Suppose that in a \( 2\times2 \) layout the cell \( (2,2) \) is empty. Show that
\( \bar{\mu}_{1\cdot}-\bar{\mu}_{2\cdot} \) is not estimable. Compute what the last-in sum of squares
of \( A \) under sum-to-zero coding tests in this case. *Hint:* with a cell missing, the interaction
columns restricted to the observed cells are no longer independent of the main-effect columns.
:::

### C. Going deeper

::: {#exr-ss-partial-sum-sign}
[C1]

With the notation of @exr-ss-two-regressor-formula, show that
\[
\frac{\text{SS}(\x_1\mid\bone,\x_2)+\text{SS}(\x_2\mid\bone,\x_1)-\text{SSR}}{\text{SST}}
=\frac{r_{12}\bigl[r_{12}(r_{y1}^2+r_{y2}^2)-2r_{y1}r_{y2}\bigr]}{1-r_{12}^2}.
\]
Deduce that the partial sums of squares add up to the regression sum of squares when
\( r_{12}=0 \), that they add up to less when \( r_{y1}=r_{y2}\neq0 \) and \( 0<r_{12}<1 \), and give
correlations for which they add up to more.
:::

