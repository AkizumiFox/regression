# Several factors, interactions and covariates

With one factor, rank deficiency is a nuisance, and a coding removes it once and for all. With two
or more factors, or a factor together with a numeric regressor, it becomes informative. The rank of the
model matrix now depends on which combinations of levels were observed. Whether a comparison is
estimable depends on how the design links the levels involved. This section works out the three cases
that occur most often: two factors without interaction, two factors with interaction, and a factor
with a numeric covariate.

## Two factors without interaction

Let factor \( A \) have levels \( i=1,\dots,a \) and factor \( B \) levels \( j=1,\dots,b \). Call a pair
\( (i,j) \) a **cell**. A cell is *occupied* if at least one observation has that combination, and we assume
every level of each factor occurs in some occupied cell. The *additive* model is
\[
\E(y)=\mu+\alpha_i+\beta_j\qquad\text{for an observation in cell }(i,j),
\]{#eq-est-additive}

with model matrix \( \X=[\bone,\Z_A,\Z_B] \) built from the two indicator matrices. The number of parameters is
\( p=1+a+b \). Since \( \Z_A\bone_a=\Z_B\bone_b=\bone \), the vectors \( (1,-\bone_a\T,\bzero\T)\T \) and
\( (1,\bzero\T,-\bone_b\T)\T \) always lie in \( \Null(\X) \), so \( \rank(\X)\le a+b-1 \). Whether equality holds, and which
differences \( \alpha_i-\alpha_k \) are estimable, depends on the pattern of occupied cells.

Picture the design as a graph. The vertices are the \( a \) levels of \( A \) and the \( b \) levels of \( B \), and each
occupied cell \( (i,j) \) is an edge joining level \( i \) of \( A \) to level \( j \) of \( B \). The design is
**connected** if this graph is connected, that is, if one can walk from any level to any other through
occupied cells.

::: {#thm-est-connected}
[Connectedness and estimability in the additive model]

In the additive model @eq-est-additive, let the design graph have \( c \) connected components. Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \rank(\X)=a+b-c \);

2. \( \alpha_i-\alpha_k \) is estimable iff levels \( i \) and \( k \) of \( A \) lie in the same component, and likewise
           for \( \beta_j-\beta_l \);

3. if the design is connected, the estimable functions are exactly the linear functions of the cell means
           \( \mu+\alpha_i+\beta_j \) over *all* cells, occupied or not. In particular every difference between levels of
           the same factor is estimable.
:::

:::

::: {.proof}
A vector \( (m,\bu\T,\bv\T)\T \) lies in \( \Null(\X) \) iff \( m+u_i+v_j=0 \) for every occupied cell \( (i,j) \). Fix
such a vector. If \( (i,j) \) and \( (k,j) \) are both occupied, then \( u_i=u_k \). If \( (i,j) \) and \( (i,l) \) are both
occupied, then \( v_j=v_l \). Following edges, \( u \) is constant on the \( A \)-levels of each component, say
\( u_i=t_C \) for all \( A \)-levels \( i \) in component \( C \), and then \( v_j=-m-t_C \) for the \( B \)-levels of \( C \).
Conversely, any choice of \( m \) and of numbers \( t_1,\dots,t_c \), one for each component, defines in this way a vector satisfying
all the equations. So \( \Null(\X) \) is parameterized linearly and injectively by \( (m,t_1,\dots,t_c) \), and it has dimension
\( c+1 \). Rank–nullity gives \( \rank(\X)=1+a+b-(c+1) \), which is (a).

For (b), by @thm-est-characterization(c), \( \alpha_i-\alpha_k \) is estimable iff \( u_i-u_k=0 \) for every null vector,
that is, iff \( t_{C(i)}=t_{C(k)} \) for all choices of the \( t \)'s. This holds iff \( C(i)=C(k) \). The argument for \( B \) is
the same.

For (c), if \( c=1 \), then (b) makes every \( \alpha_i-\alpha_k \) and \( \beta_j-\beta_l \) estimable. Pick an occupied cell
\( (i_0,j_0) \). Then \( \mu+\alpha_i+\beta_j=(\mu+\alpha_{i_0}+\beta_{j_0})+(\alpha_i-\alpha_{i_0})+(\beta_j-\beta_{j_0}) \) is estimable
for every \( (i,j) \). Conversely, every estimable function is a combination of the means of occupied cells
(@thm-est-estimable-identifiable), and these are among the cell means.
:::

Part (c) contains a point worth stating separately. In a connected additive design, the mean of an
*empty* cell is estimable. The data never observed that combination, but additivity lets it be pieced
together from cells that were observed: \( \mu+\alpha_i+\beta_j \) is the mean of an occupied cell \( (i,j') \), plus the
difference \( \beta_j-\beta_{j'} \) learned elsewhere. The estimate is only as good as the additivity assumption, which
cannot be checked at the empty cell. This is extrapolation of a subtle kind.

::: {#exm-est-connected}
[Two islands]

Take \( a=4 \), \( b=5 \) and seven occupied cells:
\( (1,1),(1,2),(2,2),(2,3) \) and \( (3,4),(4,4),(4,5) \) (levels numbered from 1). The design graph has two components:
rows 1–2 with columns 1–3, and rows 3–4 with columns 4–5. By @thm-est-connected, \( \rank(\X)=4+5-2=7 \), one
less than \( a+b-1=8 \). The difference \( \alpha_1-\alpha_2 \) is estimable, since rows 1 and 2 are linked through
column 2, but \( \alpha_1-\alpha_3 \) is not. The two islands could sit at any vertical offset from each other without
changing any fitted value. Occupying a single further cell, \( (3,3) \), joins the islands. The rank rises to \( 8 \),
and \( \alpha_1-\alpha_3 \) becomes estimable. The listing checks each statement with the rank test
@thm-est-characterization(h). The script also checks the rank formula on \( 200 \) random designs.
:::

```{.python .run #cell-two-factors-additive}
import numpy as np

def additive_matrix(cells, a, b):
    """Model matrix [1, row indicators, column indicators] for occupied (i, j) cells."""
    rows = np.array([i for i, j in cells])
    cols = np.array([j for i, j in cells])
    return np.column_stack([np.ones(len(cells)), np.eye(a)[rows], np.eye(b)[cols]])

def estimable(lam, X):
    """lambda is estimable iff appending it to the rows of X does not raise the rank."""
    return np.linalg.matrix_rank(np.vstack([X, lam])) == np.linalg.matrix_rank(X)

a, b = 4, 5
cells = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 3), (3, 3), (3, 4)]   # two separate blocks
X = additive_matrix(cells, a, b)
row_diff = lambda i, k: np.r_[0, np.eye(a)[i] - np.eye(a)[k], np.zeros(b)]
print("rank:", np.linalg.matrix_rank(X), " (a + b - 1 =", a + b - 1, ")")
print("alpha_1 - alpha_2 estimable:", estimable(row_diff(0, 1), X))
print("alpha_1 - alpha_3 estimable:", estimable(row_diff(0, 2), X))
X_linked = additive_matrix(cells + [(2, 2)], a, b)                # one more occupied cell
print("after adding cell (3, 3): rank", np.linalg.matrix_rank(X_linked),
      " alpha_1 - alpha_3 estimable:", estimable(row_diff(0, 2), X_linked))
```

A complete design, with every cell occupied, is connected. So is any design in which one level of \( B \)
appears with every level of \( A \) and one level of \( A \) with every level of \( B \). Disconnected designs arise in
practice more often than one might expect: in observational data with many sparse levels (workers and firms,
students and schools, players and teams), and in experiments where blocks were filled carelessly. The
estimable functions of a disconnected design are characterized in @exr-est-disconnected.

## Two factors with interaction

The additive model assumes that the difference between two levels of \( A \) is the same at every level of \( B \).
The *interaction model* drops that assumption:
\[
\E(y)=\mu+\alpha_i+\beta_j+\gamma_{ij}\qquad\text{for an observation in cell }(i,j),
\]{#eq-est-interaction}

with an interaction parameter \( \gamma_{ij} \) for each of the \( ab \) cells. Write \( \bW \) for the \( n\times ab \) indicator
matrix of the cells. The column of an empty cell is zero. The model matrix is \( \X=[\bone,\Z_A,\Z_B,\bW] \), with
\( p=1+a+b+ab \) columns, and \( m \) denotes the number of occupied cells.

::: {#prp-est-interaction}
[Estimability in the interaction model]

In the interaction model @eq-est-interaction with \( m \) occupied cells, write \( \mu_{ij}=\mu+\alpha_i+\beta_j+\gamma_{ij} \) for the
mean of an occupied cell.

::: {.enumerate options="label=(\alph*)"}
1. \( \C(\X)=\C(\bW) \) and \( \rank(\X)=m \): the model is the cell-means model in disguise.

2. The estimable functions are exactly the linear combinations \( \sum c_{ij}\mu_{ij} \) over occupied cells.

3. An interaction contrast \( \sum_{i,j}d_{ij}\gamma_{ij} \) whose coefficient table \( (d_{ij}) \) has zero row and column
           sums is estimable iff \( d_{ij}=0 \) for every empty cell. In particular,
           \( \gamma_{ij}-\gamma_{il}-\gamma_{kj}+\gamma_{kl} \) is estimable iff all four cells are occupied.
:::

:::

::: {.proof}
(a) Every column of \( \bone \), \( \Z_A \) and \( \Z_B \) is a sum of columns of \( \bW \). For example, the indicator of level
\( i \) of \( A \) is the sum of the indicators of the cells in row \( i \). So \( \C(\X)=\C(\bW) \). The nonzero columns of \( \bW \) are
the \( m \) indicators of occupied cells, which have disjoint supports.
(b) By @thm-est-estimable-identifiable, the estimable functions are the combinations of the expected responses,
and these are the \( \mu_{ij} \).
(c) If \( d_{ij}=0 \) off the occupied cells, then
\( \sum_{\text{occ}}d_{ij}\mu_{ij}=\sum d_{ij}\mu+\sum_i\alpha_i\sum_jd_{ij}+\sum_j\beta_j\sum_id_{ij}+\sum d_{ij}\gamma_{ij}=\sum d_{ij}\gamma_{ij} \),
because the row and column sums vanish. So the contrast is estimable. Conversely, if \( d_{i_0j_0}\ne0 \) for an empty cell
\( (i_0,j_0) \), the column of \( \gamma_{i_0j_0} \) in \( \X \) is zero, so the coordinate vector \( \bv \) of that parameter lies in
\( \Null(\X) \). The contrast's coefficient vector has inner product \( d_{i_0j_0}\ne0 \) with \( \bv \), so it is not estimable
(@thm-est-characterization(c)).
:::

In the interaction model the fitted value of each occupied cell is its own mean, and an empty cell has no
estimable mean at all. The two models treat empty cells in opposite ways. The additive model fills them by
assumption, while the interaction model declines to say anything about them. Which is appropriate is a
substantive question, and [Chapter 17](../ch17-unbalanced-data/index.html) returns to it (@thm-ub-empty-cells).

::: {#exm-est-empty-cells}
[Empty cells in a survey]

Cross-classify the \( 944 \) ANES respondents of @exm-est-party-coding by party identification (seven levels) and
education (seven levels, labelled in the statsmodels documentation from “grades 1–8” to “PhD”). Two of the \( 49 \) cells are empty: no pure independent and no
weak Republican in the sample has only an elementary-school education. Other cells are thin: only
\( 3 \) pure independents report some high school. The overparameterized interaction model has
\( p=64 \) columns and rank \( 47 \), the number of occupied cells, as
@prp-est-interaction(a) says. The interaction contrast comparing strong Democrats with pure independents across
“some high school” and “high-school graduate” involves four occupied cells and is estimable. It rests on only
\( 3 \) respondents in one cell, so it is imprecise, but estimable. The same comparison across “grades
1–8” and “high-school graduate” involves an empty cell and is not estimable. No amount of computation will produce it, and
software that prints a number for it is printing an artefact of its side condition. The additive model, with rank
\( 7+7-1=13 \), treats both comparisons as zero by assumption.
:::

```{.python .run #cell-two-factors-emptycells}
import statsmodels.api as sm

anes = sm.datasets.anes96.load_pandas().data
pid = anes["PID"].to_numpy().astype(int)              # 7 levels
educ = anes["educ"].to_numpy().astype(int) - 1        # 7 levels, now 0..6
counts = np.zeros((7, 7), dtype=int)
np.add.at(counts, (pid, educ), 1)
print("empty (PID, educ) cells:",
      [(int(i), int(j) + 1) for i, j in zip(*np.nonzero(counts == 0))])

n = len(pid)
cell = pid * 7 + educ
X_int = np.column_stack([np.ones(n), np.eye(7)[pid], np.eye(7)[educ], np.eye(49)[cell]])
print("interaction model: p =", X_int.shape[1], " rank =", np.linalg.matrix_rank(X_int),
      " occupied cells =", int((counts > 0).sum()))

def interaction_contrast(i, k, j, l):
    """gamma_ij - gamma_il - gamma_kj + gamma_kl, as a coefficient vector on X_int."""
    lam = np.zeros(X_int.shape[1])
    for (r, c), s in [((i, j), 1), ((i, l), -1), ((k, j), -1), ((k, l), 1)]:
        lam[15 + r * 7 + c] += s
    return lam

print("PID 0 vs 3, educ 2 vs 3:", estimable(interaction_contrast(0, 3, 1, 2), X_int))
print("PID 0 vs 3, educ 1 vs 3:", estimable(interaction_contrast(0, 3, 0, 2), X_int))
```

## A factor and a numeric covariate

Now let a factor with \( g \) levels, with indicator matrix \( \Z \), share the model with a numeric covariate \( \x \).
The two basic models are the **parallel lines** model
\[
\E(y)=\mu_k+\theta x\qquad\text{for an observation at level }k,
\]{#eq-est-parallel}

with model matrix \( [\Z,\x] \), and the **separate lines** model
\[
\E(y)=\mu_k+\theta_kx ,
\]{#eq-est-separate}

with model matrix \( [\Z,\Z\circ\x] \), where \( \Z\circ\x \) multiplies each indicator column entrywise by \( \x \), so its
\( k \)th column is \( \x \) on level \( k \) and zero elsewhere. The same models can be written with an intercept and any coding
of the factor. For example, the separate lines model in reference coding has columns
\( \bone,\ \Z\mathbf{C},\ \x,\ (\Z\mathbf{C})\circ\x \). The coefficients of the last block are then differences of slopes from the reference
level, and the column space is unchanged.

By @exr-est-ancova-identifiable, the common slope \( \theta \) is identifiable iff \( \x\notin\C(\Z) \), that is, iff \( x \)
varies within at least one level. In the separate lines model, the slope \( \theta_k \) is identifiable iff \( x \) varies within
level \( k \), since the \( k \)th column of \( \Z\circ\x \) is a multiple of \( \bz_k \) exactly when \( x \) is constant on level \( k \).
When the slopes are identifiable, the separate lines model is simply \( g \) simple regressions fitted at once. Its
coefficients are those of the separate fits, because the columns for different levels have disjoint supports.

The parallel lines model is where adjustment happens.

::: {#prp-est-adjusted}
[Adjusted differences]

In the parallel lines model @eq-est-parallel with \( \x\notin\C(\Z) \), let \( \bar{y}_k \) and \( \bar{x}_k \) be the level means. Then
\[
\hat{\theta}=\frac{\sum_k\sum_{i\in k}(x_i-\bar{x}_k)(y_i-\bar{y}_k)}{\sum_k\sum_{i\in k}(x_i-\bar{x}_k)^2},
\qquad
\hat{\mu}_k-\hat{\mu}_l=(\bar{y}_k-\bar{y}_l)-\hat{\theta}\,(\bar{x}_k-\bar{x}_l).
\]
:::

::: {.proof}
The model matrix \( [\Z,\x] \) has full column rank because \( \x\notin\C(\Z) \). Apply @thm-proj-fwl with \( \X_1=\x \) and
\( \X_2=\Z \). Projecting onto \( \C(\Z) \) replaces each entry by its level mean (@exm-proj-oneway-M), so \( (\I-\M_Z)\x \) and
\( (\I-\M_Z)\y \) are the within-level deviations, and \( \hat{\theta} \) is the no-intercept slope of one on the other.
This is the first formula. By @thm-proj-fwl(c), \( \hat{\bmu}=(\Z\T\Z)^{-1}\Z\T(\y-\hat{\theta}\x) \), whose \( k \)th entry is
\( \bar{y}_k-\hat{\theta}\bar{x}_k \). Subtracting gives the second formula.
:::

The difference \( \bar{y}_k-\bar{y}_l \) compares the levels as they happen to be. The adjusted difference compares them at
equal values of the covariate, correcting for the fact that the groups differ in \( x \). The correction is the slope times
the difference in covariate means. It is large when the groups are unbalanced on a covariate that matters. The
estimated slope uses only variation *within* levels, so that differences between levels in both \( x \) and \( y \) cannot be
mistaken for an effect of \( x \).

::: {#exm-est-adjusted}
[Education and television news, adjusted for age]

In the ANES data, respondents report how many times per week they watch television news. Coding education with
seven levels and “high-school graduate” as the reference, respondents with some high school but no diploma watch
\( 1.005 \) times a week more than high-school graduates. They are also older, with mean ages
\( 59.8 \) and \( 48.2 \), and news watching rises with age. Adding age to the model, the
within-level slope is \( 0.0674 \) times per week per year of age, and the adjusted difference shrinks to
\( 0.224 \). More than three quarters of the raw difference is an age difference, and
@prp-est-adjusted accounts for it exactly. Whether the adjusted figure is the “effect of education” is a
different question. Here age plausibly influences both schooling (older cohorts had less of it) and viewing habits, which
argues for adjusting, but Chapter 25 shows that adjustment is not always appropriate.

Allowing separate age slopes for three education groups (codes 1–3, high school or less, \( 313 \) respondents;
codes 4–5, \( 277 \); and the two highest levels, codes 6–7, \( 354 \)) gives slopes
\( 0.0686 \), \( 0.0729 \) and \( 0.0597 \), close to the common slope
\( 0.0672 \) of the parallel lines model with the same three groups. The residual sum of squares falls only
from \( 5623.6 \) to \( 5616.7 \) for the two extra parameters.
[Figure 8.6.1](#fig-est-slopes) shows the three fitted lines with binned means. The lines are nearly parallel, and the
highest education group starts higher and rises a little more slowly.
:::

::: {when-format="html"}
![**Figure 8.6.1.** Times per week of television news against age for three education groups in the 1996 ANES:
separate least squares lines (@eq-est-separate) and the group means within ten-year age bins (points).](separate_slopes.svg){#fig-est-slopes width=66%}
:::

::: {when-format="pdf"}
![Times per week of television news against age for three education groups in the 1996 ANES:
separate least squares lines (@eq-est-separate) and the group means within ten-year age bins (points).](separate_slopes.pdf){width=66%}
:::

```{.python .run #cell-two-factors-covariate}
tv = anes["TVnews"].to_numpy()                       # times per week watching TV news
age = anes["age"].to_numpy()
E = np.eye(7)[educ][:, [0, 1, 3, 4, 5, 6]]          # reference: high-school graduate
unadjusted = np.linalg.lstsq(np.column_stack([np.ones(n), E]), tv, rcond=None)[0]
adjusted = np.linalg.lstsq(np.column_stack([np.ones(n), E, age]), tv, rcond=None)[0]
print("some high school minus high-school graduate:",
      f"unadjusted {unadjusted[2]:.3f}, adjusted for age {adjusted[2]:.3f}")
print("mean age: some high school", age[educ == 1].mean().round(1),
      " high-school graduate", age[educ == 2].mean().round(1),
      " common age slope", adjusted[-1].round(4))
```

```{.python .run #cell-two-factors-slopes}
grp = np.select([educ <= 2, educ <= 4], [0, 1], 2)  # education codes 1-3 / 4-5 / 6-7
D = np.eye(3)[grp]
X_par = np.column_stack([D, age])                   # parallel lines: 3 intercepts, one slope
X_sep = np.column_stack([D, D * age[:, None]])      # separate lines: 3 intercepts, 3 slopes
b_par = np.linalg.lstsq(X_par, tv, rcond=None)[0]
b_sep = np.linalg.lstsq(X_sep, tv, rcond=None)[0]
sse_par = np.sum((tv - X_par @ b_par) ** 2)
sse_sep = np.sum((tv - X_sep @ b_sep) ** 2)
print("separate slopes:", b_sep[3:].round(4),
      f" SSE parallel {sse_par:.1f}, separate {sse_sep:.1f}")
```

::: {.warning}
In a separate lines model the difference between two levels is \( (\mu_k-\mu_l)+(\theta_k-\theta_l)x \), which depends on \( x \).
The coefficient that software labels as the “main effect” of level \( k \) is this difference at \( x=0 \). For a covariate like
age, that is a comparison of newborns, far outside the data. Centre the covariate at a meaningful value before reading such
coefficients, or report the difference at several values of \( x \).
:::

## Broken lines

Indicators can also switch a slope on. Let \( c \) be a known value of \( x \) and put \( (x-c)_+=\max(x-c,0) \), which
is \( x-c \) times the indicator of \( x>c \). The model
\[
\E(y)=\beta_0+\beta_1x+\beta_2(x-c)_+
\]
is a continuous line with a bend at \( c \): slope \( \beta_1 \) to the left of \( c \) and \( \beta_1+\beta_2 \) to the right. Adding the
indicator \( \mathbf{1}\{x>c\} \) itself allows a jump at \( c \) as well. With the bend point known, these are ordinary linear models, and
all of this chapter applies. For example, \( \beta_2 \) is identifiable iff the data contain values of \( x \) on both sides of
\( c \) (strictly below and strictly above) and at least three distinct values of \( x \) in all (@exr-est-broken). When \( c \) is unknown the model is nonlinear in
\( c \). Such piecewise linear terms are the simplest splines, the subject of Part IX.

## Exercises

### A. Check your understanding

::: {#exr-est-two-by-two}
[A1]

In the additive model with \( a=b=2 \), suppose cells \( (1,1) \), \( (2,1) \) and \( (2,2) \) are occupied and \( (1,2) \) is empty.
Find the rank of \( \X \) and a basis of \( \Null(\X) \). Show that \( \mu+\alpha_1+\beta_2 \), the mean of the empty cell, is estimable by
writing it as a combination of the means of occupied cells.
:::

::: {#exr-est-separate-rank}
[A2]

A factor has three levels, and the covariate \( x \) is constant on level 2 but varies on levels 1 and 3. Find the rank of the
separate lines model matrix \( [\Z,\Z\circ\x] \), and say which of \( \mu_1,\mu_2,\mu_3,\theta_1,\theta_2,\theta_3 \) are identifiable.
:::

### B. Practice

::: {#exr-est-mean-at-x0}
[B1]

In the separate lines model, show that the mean at level \( k \) and covariate value \( x_0 \), \( \mu_k+\theta_kx_0 \), is estimable iff
either \( x \) varies within level \( k \), or \( x_0 \) equals the common value of \( x \) on level \( k \). Interpret this as a statement about
extrapolation.
:::

::: {.solution}
The mean is \( \blambda\T\bbeta \) with \( \blambda \) having \( 1 \) in the position of \( \mu_k \) and \( x_0 \) in the position of \( \theta_k \). The
columns for \( \mu_k \) and \( \theta_k \) are \( \bz_k \) and \( \bz_k\circ\x \), and they involve only level \( k \). If \( x \) varies on level
\( k \), both coefficients are identifiable (@prp-est-coordinate), so \( \blambda\T\bbeta \) is estimable. If \( x\equiv x^* \) on level
\( k \), then \( \bz_k\circ\x=x^*\bz_k \), and the null space contains the vector with \( -x^* \) at \( \mu_k \) and \( 1 \) at \( \theta_k \).
Then \( \blambda\T\bv=x_0-x^* \), which vanishes iff \( x_0=x^* \). A single covariate value within a level supports no
statement about other values, however many observations share it.
:::

::: {#exr-est-broken}
[B2]

Write the broken-line model \( \E(y)=\beta_0+\beta_1x+\beta_2(x-c)_+ \) in matrix form and prove the identifiability condition
stated in the text. Show that the continuous broken line is the same as the model with two separate lines, one on each side of
\( c \), constrained to meet at \( c \), and count parameters to explain the difference of one.
:::

::: {#exr-est-general-interaction}
[B3]

In the interaction model, suppose all cells are occupied except \( (1,1) \). Show that the \( (a-1)(b-1) \) contrasts
\( \gamma_{ij}-\gamma_{i1}-\gamma_{1j}+\gamma_{11} \) are no longer estimable, but that
\( \gamma_{ij}-\gamma_{i2}-\gamma_{1j}+\gamma_{12} \) is estimable for \( i\ge2 \), \( j\ge3 \). How many linearly independent
estimable interaction contrasts remain?
:::

### C. Going deeper

::: {#exr-est-disconnected}
[C1]

In the additive model, let the design graph have components \( C_1,\dots,C_c \). Show that
\( \lambda_0\mu+\sum_i\lambda_i\alpha_i+\sum_j\kappa_j\beta_j \) is estimable iff \( \lambda_0=\sum_j\kappa_j \) and, for every
component \( C \), \( \sum_{i\in C}\lambda_i=\sum_{j\in C}\kappa_j \). Deduce @thm-est-connected(b).
:::

::: {.solution}
From the proof of @thm-est-connected, the null vectors are \( (m,\bu,\bv) \) with
\( u_i=t_{C(i)} \) and \( v_j=-m-t_{C(j)} \), for arbitrary \( m,t_1,\dots,t_c \). The inner product with the coefficient vector is
\[
\lambda_0m+\sum_i\lambda_it_{C(i)}-\sum_j\kappa_j(m+t_{C(j)})
=m\Bigl(\lambda_0-\sum_j\kappa_j\Bigr)+\sum_Ct_C\Bigl(\sum_{i\in C}\lambda_i-\sum_{j\in C}\kappa_j\Bigr).
\]
This vanishes for all \( m,t \) iff every bracket is zero, which is the stated condition
(@thm-est-characterization(c)). For \( \alpha_i-\alpha_k \), \( \lambda_0=0 \), all \( \kappa_j=0 \), and the component sums of
\( \lambda \) are zero iff \( i \) and \( k \) lie in the same component.
:::

