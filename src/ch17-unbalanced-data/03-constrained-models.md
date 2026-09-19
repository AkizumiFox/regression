# Constrained models

The cell-means model assumes nothing about the table of means. Smaller models, such as the additive model or a
three-way model without the three-factor interaction, are the cell-means model with linear constraints
\( \R\bmu=\bzero \). This section fits and tests them with the restricted least squares of
[Section 11.3](../ch11-general-linear-hypothesis/03-testable-hypotheses.html), and shows that main-effect hypotheses
can merge once the constraints hold.

## Estimation and tests under constraints

::: {#thm-ub-constrained}
[Constrained cell-means models]

In the cell-means model of @def-ub-cell-means, let \( \R \) be a \( g\times p \) matrix of rank \( g \), and let the
columns of the \( p\times(p-g) \) matrix \( \mathbf{K} \) form a basis of \( \Null(\R) \). Assume that \( \W\mathbf{K} \) has full column
rank \( p-g \). The **constrained model** is \( \E(\Y)=\W\bmu \) with \( \R\bmu=\bzero \).

::: {.enumerate options="label=(\alph*)"}
1. The constrained model is the full-rank linear model \( \E(\Y)=\W\mathbf{K}\boldsymbol{\uptheta} \), with \( \bmu=\mathbf{K}\boldsymbol{\uptheta} \).
   Every cell mean, including the means of empty cells, is estimable in it. The least squares estimate is
   \[
\hat{\bmu}_{\R}=\V_{\R}\bD\bar{\y},\qquad \V_{\R}=\mathbf{K}(\mathbf{K}\T\bD\mathbf{K})^{-1}\mathbf{K}\T,
   \qquad \Cov(\hat{\bmu}_{\R})=\sigma^2\V_{\R},
\]{#eq-ub-constrained-estimate}

   and its residual sum of squares \( \text{SSE}_{\R} \) has \( n-p+g \) degrees of freedom.

2. If every cell is occupied, then
   \( \V_{\R}=\bD^{-1}-\bD^{-1}\R\T(\R\bD^{-1}\R\T)^{-1}\R\bD^{-1} \), so
   \( \hat{\bmu}_{\R}=\bar{\y}-\bD^{-1}\R\T(\R\bD^{-1}\R\T)^{-1}\R\bar{\y} \), and
   \( \text{SSE}_{\R}-\text{SSE}=(\R\bar{\y})\T(\R\bD^{-1}\R\T)^{-1}\R\bar{\y} \) is the sum of squares, on \( g \)
   degrees of freedom, of the hypothesis \( \R\bmu=\bzero \) in the unconstrained model.

3. Under normality, for an \( \ell\times p \) matrix \( \bL \) with \( q=\rank(\bL\mathbf{K})\ge1 \), the \( F \) statistic for
   \( H:\bL\bmu=\bzero \) within the constrained model is
   \[
F=\frac{(\bL\hat{\bmu}_{\R})\T(\bL\V_{\R}\bL\T)\ginv\bL\hat{\bmu}_{\R}}{q\,s_{\R}^2},
   \qquad s_{\R}^2=\frac{\text{SSE}_{\R}}{n-p+g},
\]{#eq-ub-constrained-f}

   for any generalized inverse. It has the \( F(q,n-p+g,\gamma) \) distribution with
   \( \gamma=(\bL\bmu)\T(\bL\V_{\R}\bL\T)\ginv\bL\bmu/\sigma^2 \).

4. Within the constrained model, \( \bL_1\bmu=\bzero \) and \( \bL_2\bmu=\bzero \) are the same hypothesis iff
   \( \bL_1\mathbf{K} \) and \( \bL_2\mathbf{K} \) have the same row space, and then their \( F \) statistics coincide.
:::

:::

::: {.proof}
(a) A table satisfies \( \R\bmu=\bzero \) iff \( \bmu=\mathbf{K}\boldsymbol{\uptheta} \) for a unique \( \boldsymbol{\uptheta} \), so the constrained
model is the linear model with model matrix \( \W\mathbf{K} \), of full column rank by assumption. Every coordinate of
\( \boldsymbol{\uptheta} \) is then estimable, and so is every linear function of \( \bmu=\mathbf{K}\boldsymbol{\uptheta} \). Since
\( (\W\mathbf{K})\T\W\mathbf{K}=\mathbf{K}\T\bD\mathbf{K} \) and \( \W\T\y=\bD\bar{\y} \), the least squares estimate is
\( \hat{\boldsymbol{\uptheta}}=(\mathbf{K}\T\bD\mathbf{K})^{-1}\mathbf{K}\T\bD\bar{\y} \), with covariance \( \sigma^2(\mathbf{K}\T\bD\mathbf{K})^{-1} \) (@thm-lm-moments). Multiplying by \( \mathbf{K} \) gives @eq-ub-constrained-estimate. The model matrix has rank
\( p-g \), which gives the degrees of freedom.

(b) With every cell occupied, \( \W \) has full column rank, \( \bD^{-1} \) is the inverse of \( \W\T\W \), and
\( \bar{\y} \) is the unrestricted estimate. The constrained model is the reduced model of the testable hypothesis
\( \R\bmu=\bzero \) in the unconstrained one. So @prp-glh-restricted-ls(a) and (b), with \( \bLambda\T=\R \) and
\( \G=\bD^{-1} \), give the displayed formula for \( \hat{\bmu}_{\R} \) and
\( \text{SSE}_{\R}=\text{SSE}+(\R\bar{\y})\T(\R\bD^{-1}\R\T)^{-1}\R\bar{\y} \), the sum of squares @eq-ub-cell-f of the
hypothesis. For the covariance, let \( \V' \) be the right side of the claimed identity. On vectors \( \bD\mathbf{K}\mathbf{a} \) both
matrices give \( \mathbf{K}\mathbf{a} \), because \( \R\mathbf{K}=\bzero \). On vectors \( \R\T\mathbf{b} \) both give \( \bzero \), because
\( \mathbf{K}\T\R\T=\bzero \) and \( \V'\R\T\mathbf{b}=\bD^{-1}\R\T\mathbf{b}-\bD^{-1}\R\T\mathbf{b} \). These vectors span \( \Real^p \): the
spaces \( \C(\bD\mathbf{K}) \) and \( \C(\R\T) \) have dimensions \( p-g \) and \( g \), and if \( \bD\mathbf{K}\mathbf{a}=\R\T\mathbf{b} \), then
\( \mathbf{a}\T\mathbf{K}\T\bD\mathbf{K}\mathbf{a}=(\R\mathbf{K}\mathbf{a})\T\mathbf{b}=0 \), so \( \mathbf{K}\mathbf{a}=\bzero \) because \( \bD \) is positive definite. Two
matrices that agree on a spanning set are equal.

(c) In the full-rank model \( \W\mathbf{K}\boldsymbol{\uptheta} \) the hypothesis is \( \bL\mathbf{K}\boldsymbol{\uptheta}=\bzero \), which is testable.
@thm-glh-general-f with \( \bLambda\T=\bL\mathbf{K} \) and \( \G=(\mathbf{K}\T\bD\mathbf{K})^{-1} \) gives the statistic, since
\( \bL\mathbf{K}\hat{\boldsymbol{\uptheta}}=\bL\hat{\bmu}_{\R} \) and \( \bL\mathbf{K}(\mathbf{K}\T\bD\mathbf{K})^{-1}\mathbf{K}\T\bL\T=\bL\V_{\R}\bL\T \).

(d) Within the model, \( \bL\bmu=\bzero \) iff \( \boldsymbol{\uptheta}\in\Null(\bL\mathbf{K}) \), and two matrices have the same null
space iff they have the same row space. The hypotheses then have the same reduced model, and by
@thm-glh-general-f(b) the same \( F \) statistic.
:::

Part (d) is the new phenomenon: hypotheses that differ in the cell-means model can coincide once the constraints
hold. Part (b) pulls the cell averages onto the constraint set along the directions \( \bD^{-1}\R\T \). The rank assumption on \( \W\mathbf{K} \) holds whenever every cell is occupied; with
empty cells it is a real condition (@exr-ub-constrained-identified).

## The additive model

Take \( \R=\mathbf{C}_a\T\otimes\mathbf{C}_b\T \), the interaction hypothesis matrix of @eq-ub-kronecker-hypotheses, with
\( g=(a-1)(b-1) \). Its null space is the space of additive tables \( \mu_{ij}=\alpha_i+\beta_j \)
([Section 17.1](01-cell-means.html)), so the constrained model is the additive model
@eq-est-additive written in cell means.

::: {#cor-ub-additive}
[Main effects in the additive model]

In the two-way layout with all cells occupied, let \( \R=\mathbf{C}_a\T\otimes\mathbf{C}_b\T \). Then:

::: {.enumerate options="label=(\alph*)"}
1. \( \hat{\bmu}_{\R} \) is the table of fitted cell means of the additive model, and
   \( \text{SSE}_{\R}-\text{SSE}=\text{SS}(AB\mid\mu,A,B) \);

2. for any weights \( \bw\in\Real^b \) with \( \sum_jw_j\neq0 \), the hypothesis that \( \sum_jw_j\mu_{ij} \) is the same for all
   \( i \) is, within the additive model, the hypothesis \( \alpha_1=\dots=\alpha_a \). Its sum of squares is
   \( \text{SS}(A\mid\mu,B) \), whatever \( \bw \);

3. weights that differ from row to row do not have this property: the hypothesis that
   \( \sum_jw_{ij}\mu_{ij} \) is the same for all \( i \) depends on the column effects unless
   \( \sum_jw_{ij}\beta_j \) does not depend on \( i \).
:::

:::

::: {.proof}
(a) The constrained mean vectors \( \W\bmu \) with \( \bmu \) additive form \( \C([\bone,\Z_A,\Z_B]) \), so the
constrained fit is the additive fit. The difference of residual sums of squares is then \( \text{SS}(AB\mid\mu,A,B) \)
by definition. (b) For an additive table, \( \sum_jw_j\mu_{ij}=(\sum_jw_j)\alpha_i+\sum_jw_j\beta_j \), which is the same
for all \( i \) iff the \( \alpha_i \) are equal. So all these hypotheses define the same subspace of the constrained
model, and by @thm-ub-constrained(d) the same test. Taking \( \bw=\mathbf{e}_1 \) shows that the subspace is
\( \{\mu_{ij}=\beta_j\} \), the model with \( B \) only, and the drop in residual sum of squares from it to the additive
model is \( \text{SS}(A\mid\mu,B) \). (c) Now \( \sum_jw_{ij}\mu_{ij}=(\sum_jw_{ij})\alpha_i+\sum_jw_{ij}\beta_j \), and the second
term varies with \( i \) in general.
:::

So once additivity is accepted, the argument about types disappears for the main effects: any common weighting,
and even a single column, gives the comparison \( \alpha_i-\alpha_k \), tested by the Type II sum of squares (and by
Type III, which coincides with it there). Only the count-weighted means of Type I stay
contaminated, by part (c). The usual summaries after an additive fit are the **least squares means** (estimated
marginal means), the unweighted row means of \( \hat{\bmu}_{\R} \), with standard errors from \( \V_{\R} \).

::: {#exm-ub-additive-trial}
[The variety trial under additivity]

In @exm-ub-types-trial the interaction sum of squares \( 2.146 \) gave \( F=1.71 \) on \( 4 \) and \( 28 \) degrees of
freedom, \( p=0.176 \). Suppose we accept additivity. The constrained fit @eq-ub-constrained-estimate gives the
table of fitted means

| | S1 | S2 | S3 |
|---|---|---|---|
| V1 | 5.56 | 7.17 | 7.95 |
| V2 | 5.11 | 6.73 | 7.51 |
| V3 | 4.48 | 6.09 | 6.87 |

with \( \text{SSE}_{\R}=10.927 \) on \( 37-9+4=32 \) degrees of freedom and \( s_{\R}^2=0.3415 \). The fitted mean of
V1 at S3 is \( 7.95 \), although its two plots averaged \( 7.45 \): additivity carries V1's advantage at S1 and S2
over to S3. The test of equal unweighted variety means, @eq-ub-constrained-f, gives \( 6.184 \), the Type II sum of
squares, with \( F=9.06 \) on \( 2 \) and \( 32 \) degrees of freedom, \( p=0.00076 \). The script checks that site weights
\( 0.5,0.3,0.2 \), and the column S1 alone, give the same \( 6.184 \) (@cor-ub-additive(b)), while each variety's
own count weights give \( 0.571 \), the Type I value. The least squares means are \( 6.894 \), \( 6.449 \) and
\( 5.811 \) (standard errors \( 0.169 \), \( 0.170 \), \( 0.178 \)), and V1 exceeds V2 by \( 0.445 \) (standard error
\( 0.243 \)).

The generating means had a small interaction: V1's advantage over V2 is \( 0.8 \) at S1,
\( 0.7 \) at S2 and \( 0.1 \) at S3, and the test on \( 4 \) and \( 28 \) degrees of freedom has little power to detect it. A non-significant
interaction test is not evidence of additivity.
:::

```{.python .run #cell-constrained-constrained}
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import null_space
counts = np.array([[7, 4, 2],        # rows: varieties V1-V3; columns: sites S1-S3
                   [3, 5, 4],
                   [2, 3, 7]])
true_means = np.array([[5.8, 7.2, 7.6],
                       [5.0, 6.5, 7.5],
                       [4.4, 5.9, 7.2]])
rng = np.random.default_rng(20172)
rows = [(f"V{i + 1}", f"S{j + 1}", round(true_means[i, j] + rng.normal(0, 0.5), 1))
        for i in range(3) for j in range(3) for _ in range(counts[i, j])]
trial = pd.DataFrame(rows, columns=["variety", "site", "y"])

cells = trial.groupby(["variety", "site"])["y"]
ybar = cells.mean().to_numpy()                        # 9 cell means, row-major
n_c = cells.size().to_numpy()
Dinv = np.diag(1 / n_c)
C3 = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]])
G = np.kron(C3.T, C3.T)                               # 4 interaction contrasts: G mu = 0
K = null_space(G)                                     # 9 x 5 basis of additive tables

V_G = K @ np.linalg.inv(K.T @ np.diag(n_c) @ K) @ K.T           # Cov(mu_hat_G) / sigma^2
mu_G = V_G @ np.diag(n_c) @ ybar                                # constrained estimate
mu_G_direct = ybar - Dinv @ G.T @ np.linalg.solve(G @ Dinv @ G.T, G @ ybar)
print("additive fitted cell means:\n", mu_G.reshape(3, 3).round(3))

n, m, g = n_c.sum(), 9, G.shape[0]
sse = ((trial["y"] - cells.transform("mean")) ** 2).sum()
Gy = G @ ybar
ss_constraint = Gy @ np.linalg.solve(G @ Dinv @ G.T, Gy)       # test of additivity
sse_G = sse + ss_constraint
s2_G = sse_G / (n - m + g)                            # residual mean square of the additive model

def constrained_test(L):
    """F test of L mu = 0 within the constrained model G mu = 0."""
    u = L @ mu_G
    V = L @ V_G @ L.T
    ss = u @ np.linalg.pinv(V) @ u
    q = np.linalg.matrix_rank(L @ K)
    F = ss / q / s2_G
    return ss, q, F, stats.f.sf(F, q, n - m + g)

avg3 = np.ones((1, 3)) / 3
print("variety, unweighted:", np.round(constrained_test(np.kron(C3.T, avg3)), 4))
print("site,    unweighted:", np.round(constrained_test(np.kron(avg3, C3.T)), 4))
```

```{.python .run #cell-constrained-lsmeans}
lsmeans = (mu_G.reshape(3, 3)).mean(axis=1)           # estimated unweighted variety means
L_rows = np.kron(np.eye(3), avg3)
se_ls = np.sqrt(np.diag(L_rows @ V_G @ L_rows.T) * s2_G)
d = np.array([1.0, -1.0, 0.0]) @ L_rows
diff12, se_diff12 = d @ mu_G, np.sqrt(d @ V_G @ d * s2_G)
print("least squares means:", lsmeans.round(3), " se", se_ls.round(3))
print(f"V1 - V2: {diff12:.3f} (se {se_diff12:.3f})")
```

In the listing the constraint matrix is called `G`. The script checks the fit, the identity of
@thm-ub-constrained(b) and the tests against statsmodels.

## Several factors: choosing among constrained models

With three or more factors the useful constrained models are the hierarchical ones, named by their highest
interactions. In an \( a\times b\times d \) layout, where \( C \) has \( d \) levels, \( [AB][C] \) contains the \( AB \)
interaction, the main effect of \( C \), and everything below them. The model \( [AB][AC][BC] \) sets the
\( (a-1)(b-1)(d-1) \) contrasts \( \mathbf{C}_a\T\otimes\mathbf{C}_b\T\otimes\mathbf{C}_d\T \) to zero. Unbalanced data make these models
non-orthogonal, so we remove one term at a time, judging each step by its drop in residual sum of squares over
the error mean square of the full model.

::: {#exm-ub-threeway}
[Simplifying a three-factor curing study]

A (synthetic) study of a resin coating measured hardness after curing at three temperatures (\( T \)), with two
hardeners (\( H \)) and resin from three suppliers (\( S \)). The \( 18 \) cells received one to three panels each,
\( 37 \) in all, so the full model \( [THS] \) leaves \( 19 \) error degrees of freedom, with mean square \( 0.408 \). The
backward path:

| Step | From | To | Term removed | df | SS | \( F \) | \( p \) |
|----|--------------|------------|----------|---|-------|-------|--------|
| 1 | \( [THS] \) | \( [TH][TS][HS] \) | \( THS \) | 4 | 1.137 | 0.70 | 0.6036 |
| 2 | \( [TH][TS][HS] \) | \( [TH][HS] \) | \( TS \) | 4 | 1.845 | 1.13 | 0.3719 |
| 3 | \( [TH][HS] \) | \( [TH][S] \) | \( HS \) | 2 | 4.132 | 5.06 | 0.0173 |
| 4 | \( [TH][HS] \) | \( [T][HS] \) | \( TH \) | 2 | 13.716 | 16.81 | 0.0001 |
| 5 | \( [TH][HS] \) | \( [HS] \) | \( T \) and \( TH \) | 4 | 19.100 | 11.70 | 0.0001 |

Steps 1 and 2 cost little. Step 1 is the test of the constraints \( \mathbf{G}_3\bmu=\bzero \) with
\( \mathbf{G}_3=\mathbf{C}_3\T\otimes\mathbf{C}_2\T\otimes\mathbf{C}_3\T \), and the script recovers its sum of squares from the cell means
through @thm-ub-constrained(b). Steps 3 and 4 show that neither interaction with the hardener can go, so the path
stops at \( [TH][HS] \). In that model the effect of temperature under a given hardener is the same for all three
suppliers, so whether temperature matters at all is a single question inside the model, answered by step 5.
:::

```{.python .run #cell-constrained-threeway}
import statsmodels.formula.api as smf
rng3 = np.random.default_rng(31)
levels = [(i, j, k) for i in range(3) for j in range(2) for k in range(3)]
n3 = rng3.integers(1, 4, size=len(levels))                     # 1 to 3 runs per cell

def true_hardness(i, j, k):
    """Temperature i, hardener j, supplier k: hardener interacts with both, nothing else does."""
    return (40 + [0, 1.0, 2.0][i] + 1.5 * j + [0, 0.5, -0.5][k]
            + j * [0, -1.0, -2.0][i] + j * [0, 1.2, -0.8][k])

cure = pd.DataFrame([(f"t{i + 1}", f"h{j + 1}", f"s{k + 1}", true_hardness(i, j, k) + rng3.normal(0, 0.7))
                     for (i, j, k), r in zip(levels, n3) for _ in range(r)],
                    columns=["temp", "hardener", "supplier", "y"])

def fit(formula):
    f = smf.ols(formula, cure).fit()
    return f.ssr, int(f.df_resid)

sse_full, df_full = fit("y ~ temp * hardener * supplier")
mse_full = sse_full / df_full
path = [("[THS]", "[TH][TS][HS]", "y ~ temp * hardener * supplier", "y ~ (temp + hardener + supplier) ** 2"),
        ("[TH][TS][HS]", "[TH][HS]", "y ~ (temp + hardener + supplier) ** 2", "y ~ temp * hardener + hardener * supplier"),
        ("[TH][HS]", "[TH][S]", "y ~ temp * hardener + hardener * supplier", "y ~ temp * hardener + supplier"),
        ("[TH][HS]", "[T][HS]", "y ~ temp * hardener + hardener * supplier", "y ~ temp + hardener * supplier"),
        ("[TH][HS]", "[HS]", "y ~ temp * hardener + hardener * supplier", "y ~ hardener * supplier")]
steps = []
for big, small, f_big, f_small in path:
    (s1, d1), (s0, d0) = fit(f_big), fit(f_small)
    F = (s0 - s1) / (d0 - d1) / mse_full
    steps.append((big, small, d0 - d1, s0 - s1, F, stats.f.sf(F, d0 - d1, df_full)))
    print(f"{big:13s} -> {small:13s} df = {d0 - d1}  SS = {s0 - s1:6.3f}  F = {F:5.2f}  p = {steps[-1][5]:.4f}")
```

Such comparisons can look contradictory: a factor may matter alone but not after another. That is the unbalanced
analogue of collinearity, and it should be reported as such, not resolved by choosing a type of sum of squares.

## Exercises

### A. Check your understanding

::: {#exr-ub-constraint-count}
[A1]

In an \( a\times b\times d \) layout with every cell occupied and \( n \) observations, how many constraints define the model
\( [AB][AC][BC] \), and how many error degrees of freedom does it have? Check your answer against
step 1 of @exm-ub-threeway.
:::

### B. Practice

::: {#exr-ub-additive-2x2}
[B1]

For a \( 2\times2 \) layout with all cells occupied, take \( \R=(1,-1,-1,1) \). Use @thm-ub-constrained(b) to show that the
additive fit is \( \hat{\mu}_{ij}=\bar{y}_{ij}-(-1)^{i+j}\,\hat{\delta}/(n_{ij}H) \), where
\( \hat{\delta}=\bar{y}_{11}-\bar{y}_{12}-\bar{y}_{21}+\bar{y}_{22} \) and \( H=\sum_{ij}1/n_{ij} \). Which cells move most?
:::

::: {.solution}
Here \( \R\bar{\y}=\hat{\delta} \) and \( \R\bD^{-1}\R\T=H \), so
\( \hat{\bmu}_{\R}=\bar{\y}-\bD^{-1}\R\T\hat{\delta}/H \). The \( (i,j) \) entry of \( \bD^{-1}\R\T \) is \( (-1)^{i+j}/n_{ij} \). Each
cell average moves by \( \hat{\delta}/(n_{ij}H) \) in the direction that removes the interaction, so the least replicated
cells move most. The four moves remove \( \hat{\delta} \) exactly: \( \sum_{ij}(-1)^{i+j}\cdot(-1)^{i+j}\hat{\delta}/(n_{ij}H)=\hat{\delta} \).
:::

::: {#exr-ub-partial-constraint}
[B2]

In a \( 2\times3 \) layout with every cell occupied, impose only \( \mu_{11}-\mu_{12}-\mu_{21}+\mu_{22}=0 \). Describe the
constrained model in words. How many degrees of freedom does its residual have, and what does the \( F \) test of the
remaining interaction contrast within it test?
:::

::: {#exr-ub-lsmeans}
[B3]

In the additive model with all cells occupied, show that the least squares means satisfy
\( \hat{\bar{\mu}}_{i\cdot}-\hat{\bar{\mu}}_{k\cdot}=\hat{\alpha}_i-\hat{\alpha}_k \) for any solution of the normal equations, and that
the same difference is obtained if the fitted table is averaged over \( j \) with any fixed weights summing to one. Why does
the choice of weights still matter for the least squares means themselves?
:::

::: {.solution}
The fitted table is \( \hat{\mu}_{ij}=\hat{m}+\hat{\alpha}_i+\hat{\beta}_j \), unique as a table (@thm-proj-normal-equations). Averaging over \( j \) with weights \( w_j \) summing to one gives
\( \hat{m}+\hat{\alpha}_i+\sum_jw_j\hat{\beta}_j \), and the last term does not depend on \( i \). So the differences are
\( \hat{\alpha}_i-\hat{\alpha}_k \), an estimable function whose value does not depend on the solution chosen. The
means themselves contain \( \sum_jw_j\hat{\beta}_j \), which does depend on the weights. They are levels, not
comparisons, and a level is a mean over some population of columns.
:::

### C. Going deeper

::: {#exr-ub-unweighted-means}
[C1]

The *method of unweighted means* estimates \( \alpha_1-\alpha_2 \) in a \( 2\times b \) additive layout by
\( b^{-1}\sum_j(\bar{y}_{1j}-\bar{y}_{2j}) \), as if every cell had the same count. Show that the least squares estimate is
\( \sum_jh_j(\bar{y}_{1j}-\bar{y}_{2j})/\sum_jh_j \) with \( h_j=n_{1j}n_{2j}/n_{\cdot j} \). Deduce that both are unbiased under
additivity, and that the least squares estimate has the smaller variance unless all \( h_j \) are equal.
:::

::: {.solution}
Under additivity each \( \bar{y}_{1j}-\bar{y}_{2j} \) is unbiased for \( \delta=\alpha_1-\alpha_2 \), with variance \( \sigma^2/h_j \) (@exr-ub-type2-harmonic), and the \( b \) differences are independent. Any weighted average with weights summing to one is
unbiased. Its variance \( \sigma^2\sum_jv_j^2/h_j \) is minimized, by the Cauchy–Schwarz argument of
@prp-ub-marginal-means, at \( v_j\propto h_j \). The minimum variance is \( \sigma^2/\sum_jh_j \). The least squares estimate is the best
linear unbiased estimator (@thm-opt-gauss-markov), and it is a linear function of the cell averages, which are
sufficient here. So it is the \( h \)-weighted average. The unweighted average has \( v_j=1/b \), and it is optimal only if all
\( h_j \) are equal.
:::

::: {#exr-ub-constrained-identified}
[C2]

Allow empty cells. Show that \( \W\mathbf{K} \) in @thm-ub-constrained has full column rank iff no nonzero table satisfying
\( \R\bmu=\bzero \) vanishes on every occupied cell. Deduce that for the additive constraints this holds iff the design graph of
@thm-est-connected is connected.
:::

::: {.solution}
\( \W\mathbf{K}\boldsymbol{\uptheta}=\bzero \) iff the table \( \mathbf{K}\boldsymbol{\uptheta} \) vanishes on every occupied cell, since the rows of \( \W \) pick out the
occupied cells. The columns of \( \mathbf{K} \) are independent, so \( \W\mathbf{K} \) has a nonzero null vector iff some nonzero
table satisfying the constraints vanishes on \( O \). For additive tables \( \alpha_i+\beta_j \), vanishing on \( O \) means
\( \alpha_i=-\beta_j \) for every occupied \( (i,j) \). As in the proof of @thm-est-connected, this forces \( \alpha_i=t_C \) and
\( \beta_j=-t_C \) on each component \( C \) of the design graph. The table is nonzero iff the \( t_C \) are not all equal, which is
possible iff there are at least two components.
:::
