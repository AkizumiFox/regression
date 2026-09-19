# All subsets and stepwise search

With \( m \) candidate regressors and the intercept always in, there are \( 2^m \) subsets. This section describes exhaustive
search, branch and bound and the stepwise shortcuts, and shows that the subset they find depends on small features of
the data, far more than the fitted values do.

## Searching all subsets

The search splits into two stages, because every criterion based on the residual sum of squares agrees within a
size.

::: {#prp-sel-same-size}
[One size, one ranking]

Among candidates with the same number of columns \( p \), the orderings by \( C_p \), AIC, AICc, BIC, \( R^2 \),
\( \bar{R}^2 \) and \( \text{SSE} \) coincide.
:::

::: {.proof}
For fixed \( p \), \( n \) and \( s^2 \), each of \( C_p \), \( n\log(\text{SSE}/n) \) plus a penalty depending on \( p \),
\( 1-\text{SSE}/\text{SST} \) and \( 1-(n-1)\text{SSE}/\{(n-p)\text{SST}\} \) is a monotone function of \( \text{SSE} \), increasing for the
criteria to be minimized and decreasing for the two \( R^2 \).
:::

So it suffices to find, for each size, the subset with the smallest residual sum of squares, and then compare the
\( m+1 \) winners with the criterion. Cross-validation does not reduce in this way, since its value is not a function
of \( \text{SSE} \) and \( p \); in practice it is applied to the size winners, as in @exm-sel-state-cv.

Visiting all \( 2^m \) subsets is cheap per subset. Sweeping on one pivot (@def-cmp-sweep) adds or removes a column in
\( O(m^2) \) operations, and a Gray-code order, in which consecutive subsets differ in one column, gives every residual sum
of squares with one sweep (Schatzoff, Tsao and Fienberg 1968). That makes \( m \) up to about \( 20 \) routine. Beyond
that the number of subsets is the obstacle: finding the best subset of a given size is NP-hard in general (Natarajan 1995).

## Branch and bound

The key fact is that removing columns never reduces the residual sum of squares: if
\( \C(\X_A)\subseteq\C(\X_B) \), then \( \text{SSE}_A-\text{SSE}_B=\norm{(\M_B-\M_A)\y}^2\ge0 \) by @thm-proj-nested. A subset whose
residual sum of squares is already too large therefore cannot have a good subset below it.

Organize the subsets as a tree. The root is the full set \( \{1,\dots,m\} \). A node is a subset \( S \) together with the
index \( \ell \) of the last column deleted to reach it (\( \ell=0 \) at the root), and its children delete one more column
\( j\in S \) with \( j>\ell \). Deleting in increasing order of index reaches every subset along exactly one path. The
descendants of \( (S,\ell) \) are the subsets \( T\subseteq S \) that contain \( F=\{j\in S:j\le\ell\} \), the columns that can no
longer be deleted.

::: {#prp-sel-branch-bound}
[Branch and bound for \( C_p \)]

Let \( c^* \) be the smallest \( C_p \) found so far. If a node \( (S,\ell) \) satisfies
\[
\frac{\text{SSE}_S}{s^2}-n+2(\lvert F\rvert+1)\ \ge\ c^* ,
\]
then no descendant of the node has \( C_p<c^* \), and the node need not be expanded. A depth-first search that
evaluates every expanded node and skips such nodes' descendants returns a subset of minimal \( C_p \).
:::

::: {.proof}
A descendant \( T \) contains \( F \) and is contained in \( S \), so \( \text{SSE}_T\ge\text{SSE}_S \) and \( \lvert T\rvert\ge\lvert F\rvert \),
with the intercept counted separately. Hence its \( C_p=\text{SSE}_T/s^2-n+2(\lvert T\rvert+1) \) is at least the left side of the
display. The search only skips subsets that cannot beat a value already attained, so the minimum it reports is the
minimum over all subsets.
:::

The same bound works for AIC, BIC or the best subset of each size. Furnival and Wilson (1974) combined bounds of this
kind with an efficient order of sweeps in the "leaps and bounds" algorithm, which handles \( 30 \) to \( 40 \) regressors.
Bertsimas, King and Mazumder (2016) solved much larger problems by mixed-integer optimization, with a certificate of
optimality.

::: {#exm-sel-branch-bound}
[How much branch and bound saves]

Take \( n=100 \) cases and \( m=15 \) independent normal regressors, five with nonzero slopes
(\( 1,-0.8,0.6,0.4,0.3 \)) and ten with zero slopes. Exhaustive search evaluates \( 32768 \) subsets. The branch and bound
of @prp-sel-branch-bound finds the same \( C_p \)-best subset, with \( 6 \) regressors (the four largest slopes and two
irrelevant columns), after \( 1110 \) evaluations, a fraction \( 0.034 \) of the work.
:::

```{.python .run #cell-subset-search-bb}
import itertools
import numpy as np
rng = np.random.default_rng(2908)
n, m = 100, 15
X = rng.normal(size=(n, m))
beta = np.r_[1.0, -0.8, 0.6, 0.4, 0.3, np.zeros(m - 5)]
y = 2 + X @ beta + rng.normal(size=n)
Xc, yc = X - X.mean(axis=0), y - y.mean()                  # centring accounts for the intercept
G, g = Xc.T @ Xc, Xc.T @ yc


def sse(S):
    S = list(S)
    if not S:
        return yc @ yc
    return yc @ yc - g[S] @ np.linalg.solve(G[np.ix_(S, S)], g[S])


s2 = sse(range(m)) / (n - m - 1)


def cp(S):
    return sse(S) / s2 - n + 2 * (len(S) + 1)


def branch_and_bound():
    """Children delete a column with a larger index than the last one deleted, so every
    subset is visited at most once; a node is not expanded when no descendant can win."""
    best = [np.inf, None]
    count = [0]

    def visit(S, last):
        e = sse(S)
        count[0] += 1
        value = e / s2 - n + 2 * (len(S) + 1)
        if value < best[0]:
            best[:] = [value, S]
        forced = [j for j in S if j <= last]              # kept by every descendant
        if e / s2 - n + 2 * (len(forced) + 1) >= best[0]:
            return                                          # SSE only grows as columns go
        for j in S:
            if j > last:
                visit(tuple(k for k in S if k != j), j)

    visit(tuple(range(m)), -1)
    return best[1], best[0], count[0]


S_bb, cp_bb, fits = branch_and_bound()
print("branch and bound:", S_bb, round(cp_bb, 3), "after", fits, "fits of", 2**m)
```

## Stepwise search

Stepwise methods give up the guarantee of finding the best subset in exchange for examining only \( O(m^2) \) of them.

::: {#def-sel-stepwise}
[Stepwise selection]

Let the current subset be \( S \), with \( p_S \) columns including the intercept.

::: {.enumerate options="label=(\alph*)"}
1. **Forward selection** starts from the intercept alone. At each step it adds the column \( j\notin S \) that most reduces
   \( \text{SSE} \), provided its **\( F \)-to-enter** \( F_j=(\text{SSE}_S-\text{SSE}_{S+j})/\{\text{SSE}_{S+j}/(n-p_S-1)\} \) exceeds a
   threshold \( F_{\text{in}} \), and stops otherwise.

2. **Backward elimination** starts from the full model. At each step it removes the column \( j\in S \) with the smallest
   **\( F \)-to-remove** \( F_j=(\text{SSE}_{S-j}-\text{SSE}_S)/\{\text{SSE}_S/(n-p_S)\} \), provided this is below a threshold
   \( F_{\text{out}} \), and stops otherwise.

3. **Stepwise regression** (Efroymson 1960) alternates: after each forward step it removes any column whose
   \( F \)-to-remove has fallen below \( F_{\text{out}} \), with \( F_{\text{out}}\le F_{\text{in}} \).
:::

The thresholds may be replaced by a criterion: add or remove the column that lowers AIC (or \( C_p \), or BIC) most, and
stop when no move lowers it.
:::

The single step of forward selection can be described in several equivalent ways, all of them versions of the
Frisch–Waugh–Lovell theorem.

::: {#prp-sel-forward-step}
[What a stepwise step maximizes]

Let \( \tilde{\y}=(\I-\M_S)\y \) and \( \tilde{\x}_j=(\I-\M_S)\x_j\ne\bzero \), and let
\( r_j=\tilde{\y}\T\tilde{\x}_j/(\norm{\tilde{\y}}\norm{\tilde{\x}_j}) \) be the partial correlation of \( \y \) and \( \x_j \) given the
columns of \( S \) (@def-cor-sample-partial).

::: {.enumerate options="label=(\alph*)"}
1. \( \text{SSE}_S-\text{SSE}_{S+j}=r_j^2\,\text{SSE}_S \) and \( F_j=(n-p_S-1)\,r_j^2/(1-r_j^2) \). So the column that most reduces
   \( \text{SSE} \) is the one with the largest \( r_j^2 \), the largest \( F \)-to-enter, and the smallest \( C_p \), AIC or BIC
   of the enlarged model.

2. For \( j\in S \), \( \text{SSE}_{S-j}-\text{SSE}_S=t_j^2s_S^2 \), where \( t_j \) is the \( t \) statistic of \( \x_j \) in the model \( S \)
   and \( s_S^2=\text{SSE}_S/(n-p_S) \). So backward elimination removes the column with the smallest \( \lvert t_j\rvert \), and its
   \( F \)-to-remove is \( t_j^2 \).
:::

:::

::: {.proof}
(a) By @thm-proj-fwl(b), the residual vector of the model \( S+j \) is the residual of the regression of \( \tilde{\y} \) on
\( \tilde{\x}_j \), so \( \text{SSE}_{S+j}=\norm{\tilde{\y}}^2-(\tilde{\x}_j\T\tilde{\y})^2/\norm{\tilde{\x}_j}^2=(1-r_j^2)\,\text{SSE}_S \). Then
\( F_j=(n-p_S-1)r_j^2/(1-r_j^2) \) is increasing in \( r_j^2 \). At a fixed size the criteria rank by \( \text{SSE} \) (@prp-sel-same-size).
(b) Apply (a) to the model \( S-j \) and the column \( \x_j \): the increase is \( (\tilde{\x}_j\T\tilde{\y})^2/\norm{\tilde{\x}_j}^2 \) with
residuals now taken on \( S-j \). By @thm-proj-fwl(a) this is \( \hat{\beta}_j^2\norm{\tilde{\x}_j}^2 \), and by
@eq-proj-vif-preview \( \Var(\hat{\beta}_j)=\sigma^2/\norm{\tilde{\x}_j}^2 \), so the increase is
\( \hat{\beta}_j^2s_S^2/\widehat{\Var}(\hat{\beta}_j)=t_j^2s_S^2 \).
:::

Criterion-based stopping is therefore threshold-based stopping with a particular threshold: forward selection by AIC
adds a column iff its \( F \)-to-enter exceeds about \( 2 \), and by BIC iff it exceeds about \( \log n \) (the table in
[Section 29.2](02-cp-aic-bic.html)). The traditional thresholds \( F_{\text{in}}=F_{\text{out}}=4 \) correspond to BIC at
\( n\approx55 \).

For the state data, forward selection by AIC adds single parenthood and then urbanization and stops. Backward elimination
by AIC removes white, high school, poverty (in that order) and stops at the same two regressors, which are also the
all-subsets AIC choice. Forward and backward computed the criterion for \( 13 \) and \( 15 \) subsets, against \( 32 \) for
exhaustive search. Agreement is common but not guaranteed.

::: {#exm-sel-forward-miss}
[A pair that forward selection cannot find]

Take \( n=60 \) cases of three regressors, \( \x_1 \) and \( \x_2 \) with sample correlation \( 0.96 \) and \( \x_3 \)
independent of both, and let \( \y=3(\x_1-\x_2)+0.4\,\x_3+\be \). The signal lies mostly in the difference
\( \x_1-\x_2 \), a direction in which the regressors vary little. The sample correlations of \( \y \) with \( \x_1 \), \( \x_2 \)
and \( \x_3 \) are \( -0.02 \), \( -0.19 \) and \( 0.41 \). Forward selection enters \( \x_3 \) first and then \( \x_2 \), reaching
\( R^2=0.19 \). The best pair is \( \{\x_1,\x_2\} \), with \( R^2=0.39 \). Neither of its members looks useful alone, so a
procedure that adds one column at a time never sees it. Backward elimination starts from a model that contains the pair and can
keep it, but it has blind spots of its own (@exr-sel-not-nested) and cannot start when \( m\ge n \).
:::

```{.python .run #cell-subset-search-forward-miss}
import itertools
import numpy as np
rng = np.random.default_rng(2909)
n2 = 60
z = rng.normal(size=(n2, 3))
x1 = z[:, 0]
x2 = 0.95 * z[:, 0] + np.sqrt(1 - 0.95**2) * z[:, 1]      # corr(x1, x2) about 0.95
x3 = z[:, 2]
Xf = np.column_stack([x1, x2, x3])
yf = 3 * (x1 - x2) + 0.4 * x3 + rng.normal(size=n2)       # the signal is in x1 - x2


def rss(cols):
    A = np.column_stack([np.ones(n2)] + [Xf[:, j] for j in cols])
    b, *_ = np.linalg.lstsq(A, yf, rcond=None)
    return np.sum((yf - A @ b) ** 2)


first = min(range(3), key=lambda j: rss([j]))
second = min((j for j in range(3) if j != first), key=lambda j: rss([first, j]))
best_pair = min(itertools.combinations(range(3), 2), key=rss)
print("forward chooses", sorted([first, second]), "RSS", round(rss([first, second]), 1))
print("best pair      ", list(best_pair), "RSS", round(rss(best_pair), 1))
```

Collinear regressors are exactly where this happens: a contrast between them can be highly predictive while each
member is not ([Chapter 26](../ch26-collinearity/index.html) describes such directions through the eigenvectors of
\( \X\T\X \)).

## Comparing the procedures

::: {#exm-sel-compare}
[Selection procedures in a simulation]

Take \( n=60 \) cases of eight normal regressors with all pairwise correlations \( 0.5 \), \( \sigma=1 \), an intercept, and
three nonzero slopes, either weak (\( 1,0.5,0.25 \)) or strong (\( 1,0.8,0.6 \)). Each procedure is scored by its loss at
the design points, \( \norm{\X(\hat{\bbeta}-\bbeta)}^2/\sigma^2 \), with the intercept's expected contribution \( 1 \) included. Over
\( 2000 \) data sets:

*Weak slopes:*

| procedure | loss | Pr(true) | size |
|---|---|---|---|
| least squares, all eight | 8.98 | 0.00 | 8.00 |
| true subset (oracle) | 3.97 | 1.00 | 3.00 |
| all subsets, AIC | 7.81 | 0.18 | 3.42 |
| all subsets, BIC | 7.33 | 0.18 | 2.52 |
| all subsets, leave-one-out | 7.76 | 0.19 | 3.34 |
| forward, AIC | 7.70 | 0.19 | 3.34 |
| backward, AIC | 7.89 | 0.17 | 3.50 |

*Strong slopes:*

| procedure | loss | Pr(true) | size |
|---|---|---|---|
| least squares, all eight | 9.02 | 0.00 | 8.00 |
| true subset (oracle) | 4.06 | 1.00 | 3.00 |
| all subsets, AIC | 7.27 | 0.37 | 3.93 |
| all subsets, BIC | 6.17 | 0.70 | 3.20 |
| all subsets, leave-one-out | 7.17 | 0.39 | 3.85 |
| forward, AIC | 7.16 | 0.39 | 3.86 |
| backward, AIC | 7.35 | 0.35 | 3.98 |

Every selection rule beats least squares on all eight regressors, whose expected loss is \( 9 \) (@thm-sel-optimism(c)), but none comes close to the oracle's \( 4 \). The search method hardly matters: forward selection
by AIC returned the all-subsets AIC choice in a fraction \( 0.95 \) of the weak data sets and \( 0.96 \) of the strong ones.
The criterion matters more: with strong effects BIC finds the true subset in a fraction \( 0.70 \) of data sets and has
the smallest loss, while with a weak third slope every rule usually misses it, BIC keeping \( 2.52 \) regressors on average.
:::

The simulation illustrates @prp-sel-aic-bic in a finite sample: where a small true model exists and its effects are
large, BIC wins; with many small effects the comparison tilts towards AIC and cross-validation.

## Instability

A procedure is **unstable** if small changes in the data change its output a great deal. Subset selection is a
discontinuous function of the data: a column is in or out. Breiman (1996b) showed that this instability inflates the
variance of the predictions that follow, and that averaging over perturbed data sets can reduce it.

::: {#exm-sel-instability}
[Bootstrapping the state regressions]

Resample the \( 50 \) states with replacement (the case bootstrap of
[Section 23.1](../ch23-resampling-inference/01-cases-and-residuals.html)) \( 1000 \) times and repeat the all-subsets AIC
and BIC choices on each resample. AIC chooses \( 27 \) different subsets of the \( 32 \) and BIC \( 25 \). The AIC choice
for the original data, single parenthood and urbanization, is chosen in a fraction \( 0.15 \) of resamples, and the
most frequent AIC choice, poverty with single parenthood and urbanization, in \( 0.20 \). The BIC choice, single
parenthood alone, is its most frequent choice, at \( 0.31 \). [Figure 29.4.1](#fig-sel-instability) shows how often each
regressor is included. Even single parenthood, whose \( t \) statistic in the chosen model is about \( 5.5 \), is left out
by AIC in a fifth of the resamples. It is correlated \( -0.70 \) with the percentage of white residents, and the two can
stand in for each other.
:::

::: {when-format="html"}
![**Figure 29.4.1.** Fraction of \( 1000 \) bootstrap resamples of the \( 50 \) states in which each regressor is in the
subset chosen by AIC or BIC (log violent crime, all \( 32 \) subsets).](instability.svg){#fig-sel-instability width=72%}
:::

::: {when-format="pdf"}
![Fraction of \( 1000 \) bootstrap resamples of the \( 50 \) states in which each regressor is in the
subset chosen by AIC or BIC (log violent crime, all \( 32 \) subsets).](instability.pdf){width=72%}
:::

```{.python .run #cell-instability-bootstrap}
import itertools
from collections import Counter
import numpy as np
import statsmodels.api as sm
data = sm.datasets.statecrime.load_pandas().data.drop(index="District of Columbia")
names = ["hs_grad", "poverty", "single", "white", "urban"]
y_all = np.log(data["violent"].to_numpy())
Z_all = data[names].to_numpy()
n = len(y_all)
subsets = [S for k in range(6) for S in itertools.combinations(range(5), k)]


def select(Z, y, penalty):
    def crit(S):
        X = np.column_stack([np.ones(len(y)), Z[:, list(S)]])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        return n * np.log(np.sum((y - X @ b) ** 2) / n) + penalty * (len(S) + 2)
    return min(subsets, key=crit)


rng = np.random.default_rng(2911)
B = 200                                                       # the book's figure uses 1000
chosen = {"AIC": Counter(), "BIC": Counter()}
for _ in range(B):
    idx = rng.integers(0, n, size=n)                          # resample whole states
    chosen["AIC"][select(Z_all[idx], y_all[idx], 2.0)] += 1
    chosen["BIC"][select(Z_all[idx], y_all[idx], np.log(n))] += 1
for crit, counts in chosen.items():
    print(crit, "distinct models:", len(counts))
    for S, c in counts.most_common(3):
        print("   ", [names[j] for j in S], c / B)
```

Instability is a property of the question as much as of the method. The three models favoured by the criteria and the
bootstrap (single parenthood alone, with urbanization, and with poverty and urbanization) predict almost equally
well (their leave-one-out errors in @exm-sel-state-cv lie between \( 0.1129 \) and \( 0.1184 \)), so any rule that must
name one will name different ones for similar data. The remedies are to report the uncertainty, as in the figure, or
to stop requiring a single model: average over bootstrap resamples (bagging, Breiman 1996a) or over models weighted by
posterior probability (Hoeting, Madigan, Raftery and Volinsky 1999), or shrink instead of selecting
([Chapter 27](../ch27-shrinkage/index.html)).

## Exercises

### A. Check your understanding

::: {#exr-sel-count}
[A1]

How many subsets must exhaustive search examine for \( m=20 \) and \( m=40 \) candidate regressors? At a million subsets per
second, how long does each take?
:::

### B. Practice

::: {#exr-sel-tree}
[B1]

Show that in the tree of @prp-sel-branch-bound every subset of \( \{1,\dots,m\} \) is reached along exactly one path from the
root, and that the descendants of \( (S,\ell) \) are exactly the subsets \( T \) with \( F\subseteq T\subseteq S \).
:::

::: {.solution}
A subset \( T \) with complement \( \{j_1<\dots<j_r\} \) in \( \{1,\dots,m\} \) is reached by deleting \( j_1,\dots,j_r \) in increasing order,
and every path deletes in increasing order, so this is the only path. A descendant of \( (S,\ell) \) deletes further columns
of \( S \) with indices above \( \ell \), so it contains \( F=\{j\in S:j\le\ell\} \); conversely any \( T \) with \( F\subseteq T\subseteq S \) is
reached from \( (S,\ell) \) by deleting the columns of \( S\setminus T \), all of which exceed \( \ell \), in increasing order.
:::

::: {#exr-sel-no-cycle}
[B2]

In stepwise regression, show that the column just entered has \( F \)-to-remove equal to its \( F \)-to-enter. Deduce that
it cannot be removed at the next step if \( F_{\text{out}}\le F_{\text{in}} \), and explain what can happen if
\( F_{\text{out}}>F_{\text{in}} \).
:::

::: {.solution}
After \( j \) enters, the model is \( S+j \) with \( p_S+1 \) columns, and its \( F \)-to-remove is
\( (\text{SSE}_S-\text{SSE}_{S+j})/\{\text{SSE}_{S+j}/(n-p_S-1)\} \), which is the \( F \)-to-enter of \( j \) at the previous step. It exceeded
\( F_{\text{in}}\ge F_{\text{out}} \), so \( j \) stays. If \( F_{\text{out}}>F_{\text{in}} \), a column with \( F_{\text{in}}<F_j<F_{\text{out}} \)
is added and immediately removed, returning to \( S \), where it is added again: the procedure cycles.
:::

::: {#exr-sel-backward-t}
[B3]

Show that a column's \( F \)-to-remove in backward elimination is the square of its \( t \) statistic, and that backward elimination
with \( F_{\text{out}}=t^2_{\alpha/2}(\nu) \) is the rule "drop the least significant coefficient while any is not significant
at level \( \alpha \)". Why do the \( p \)-values of the surviving coefficients not have their nominal meaning?
:::

### C. Going deeper

::: {#exr-sel-not-nested}
[C1]

In @exm-sel-forward-miss, show that the best subset of size one is \( \{\x_3\} \) and the best of size two is \( \{\x_1,\x_2\} \), so
the best subsets of consecutive sizes need not be nested. Why does this not contradict @thm-proj-nested?
:::

::: {.solution}
Among single columns the residual sum of squares is smallest for the largest squared correlation with \( \y \), and
\( 0.41^2 \) is far above \( 0.19^2 \) and \( 0.02^2 \). The best pair was computed in the example. @thm-proj-nested compares a
subset with its own supersets, so \( \text{SSE}_{\{3\}}\ge\text{SSE}_{\{1,3\}} \) and \( \text{SSE}_{\{3\}}\ge\text{SSE}_{\{2,3\}} \); it says nothing
about \( \{1,2\} \), which does not contain \( \{3\} \). Because best subsets need not be nested, a search that only adds columns
(forward) or only removes them (backward) can miss the best subset of some size.
:::

