# Missing observations in designed experiments

Designed experiments lose observations: a plot floods, a sample is dropped, a strip of timber splits. The balanced formulas of
the last two sections then no longer apply. The correct response is simple: the model still holds for the remaining observations,
so fit it to them and test by comparing models ([Chapter 11](../ch11-general-linear-hypothesis/index.html)). Before computers,
experimenters instead *filled in the missing value so that the balanced formulas give the exact least squares answer.* The trick
shows exactly what the balanced formulas get right and wrong with a filled-in value, and it is the simplest instance of the EM
algorithm.

## Filling in, exactly

Let \( \X \) be the \( n\times p \) model matrix of the complete design, of rank \( r \), and suppose the observations in a set of \( m \) rows are
missing. Reorder so that \( \X=\begin{pmatrix}\X_o\\\X_m\end{pmatrix} \) and \( \y=\begin{pmatrix}\y_o\\\cdot\end{pmatrix} \), with \( \y_o \) the \( n-m \)
observed values. Let \( \M \) be the projection onto \( \C(\X) \) for the *complete* design, partitioned conformably into blocks
\( \M_{oo},\M_{om},\M_{mo},\M_{mm} \), and let \( \mathbf{D}=\begin{pmatrix}\bzero\\\I_m\end{pmatrix} \) be the \( n\times m \) matrix of indicators of the
missing rows. The key assumption is that the loss does not destroy estimability:
\[
\rank(\X_o)=\rank(\X)=r .
\]{#eq-dsn-missing-rank}

Since \( \Null(\X)\subseteq\Null(\X_o) \) always, @eq-dsn-missing-rank says that \( \Null(\X_o)=\Null(\X) \): every function estimable in the complete
design remains estimable. A single missing observation in a randomized block design with \( t,b\ge2 \), or in a Latin square with \( a\ge3 \), never violates it.
Losing every observation on some treatment does, and so does losing the only observation of a cell in an unreplicated factorial fitted
with all its interactions.

::: {#thm-dsn-missing}
[Missing observations]

Assume @eq-dsn-missing-rank, and let \( \hbeta_o \) be any least squares estimate from the observed data alone, that is, for the model
\( \E(\Y_o)=\X_o\bbeta \). Put \( \bz^*=\X_m\hbeta_o \), the fitted values at the missing positions, and for \( \bz\in\Real^m \) let
\( \y(\bz)=\begin{pmatrix}\y_o\\\bz\end{pmatrix} \) be the data completed with \( \bz \).

::: {.enumerate options="label=(\alph*)"}
1. *(Filling in.)* \( \hbeta_o \) is a least squares estimate for the complete data \( \y(\bz^*) \), and the residual sum of squares of the
   complete-data analysis of \( \y(\bz^*) \) equals \( \text{SSE}_o \), the residual sum of squares of the observed data.

2. *(Characterization.)* The eigenvalues of \( \M_{mm} \) lie in \( [0,1) \), and \( \bz^* \) is the unique solution of
   \[
   \bz=\M_{mo}\y_o+\M_{mm}\bz ,
   \]{#eq-dsn-fill-equation}

   that is, the unique \( \bz \) that equals its own complete-data fitted value.

3. *(Iteration.)* For any starting value \( \bz_0 \), the sequence \( \bz_{s+1}=\M_{mo}\y_o+\M_{mm}\bz_s \), which fills in the current values,
   fits the complete-data model and replaces the missing values by their fitted values, converges to \( \bz^* \). The error
   \( \bz_s-\bz^* \) equals \( \M_{mm}^s(\bz_0-\bz^*) \).

4. *(Covariate method.)* For any \( \bw\in\Real^m \), the model \( \E\,\y(\bw)=\X\bbeta+\mathbf{D}\bgamma \) has rank \( r+m \). Its least squares fit
   has \( \X_o\hbeta=\X_o\hbeta_o \), \( \hat{\bgamma}=\bw-\bz^* \), and residual sum of squares \( \text{SSE}_o \), on \( n-m-r \) degrees of
   freedom.

5. *(Tests.)* Let \( \C(\X_0)\subseteq\C(\X) \) with \( \rank(\X_{0,o})=\rank(\X_0) \). The correct sum of squares for the reduced model is
   \( \text{SSE}_{0,o}-\text{SSE}_o \), the difference of the residual sums of squares of the two models fitted to the observed data, on
   \( r-r_0 \) degrees of freedom, with error degrees of freedom \( n-m-r \). The sum of squares computed from the complete-data formulas
   applied to \( \y(\bz^*) \) is at least as large.
:::

:::

::: {.proof}
(a) For any \( \bb \), the complete-data criterion splits as
\[
\begin{aligned}
\norm{\y(\bz^*)-\X\bb}^2&=\norm{\y_o-\X_o\bb}^2+\norm{\bz^*-\X_m\bb}^2\\
&\ge\norm{\y_o-\X_o\bb}^2\ \ge\ \text{SSE}_o .
\end{aligned}
\]
At \( \bb=\hbeta_o \) both inequalities are equalities, because \( \bz^*=\X_m\hbeta_o \) and \( \hbeta_o \) minimizes the observed-data criterion. So
\( \hbeta_o \) minimizes the complete-data criterion, and the minimum is \( \text{SSE}_o \).

(b) By (a), the complete-data fitted vector of \( \y(\bz^*) \) is \( \X\hbeta_o \), whose missing part is \( \bz^* \). The fitted vector is
\( \M\y(\bz^*) \), whose missing part is \( \M_{mo}\y_o+\M_{mm}\bz^* \). So \( \bz^* \) solves @eq-dsn-fill-equation. For uniqueness, note that
\( \M_{mm}=\mathbf{D}\T\M\mathbf{D} \) is symmetric and nonnegative definite, and for \( \bu\in\Real^m \),
\[
\bu\T\M_{mm}\bu=\norm{\M\mathbf{D}\bu}^2\le\norm{\mathbf{D}\bu}^2=\norm{\bu}^2 ,
\]
with equality iff \( \mathbf{D}\bu\in\C(\X) \), that is, iff \( \mathbf{D}\bu=\X\bb \) for some \( \bb \). Such a \( \bb \) has \( \X_o\bb=\bzero \) and \( \X_m\bb=\bu \).
Under @eq-dsn-missing-rank, \( \X_o\bb=\bzero \) implies \( \X\bb=\bzero \), so \( \bu=\bzero \). Hence \( \bu\T\M_{mm}\bu<\norm{\bu}^2 \) for \( \bu\ne\bzero \),
every eigenvalue of \( \M_{mm} \) is less than \( 1 \), and \( \I-\M_{mm} \) is nonsingular. Equation @eq-dsn-fill-equation is
\( (\I-\M_{mm})\bz=\M_{mo}\y_o \), which has exactly one solution.

(c) Subtracting @eq-dsn-fill-equation for \( \bz^* \) from the recursion gives \( \bz_{s+1}-\bz^*=\M_{mm}(\bz_s-\bz^*) \). The spectral radius of the
symmetric matrix \( \M_{mm} \) is its largest eigenvalue, which is less than \( 1 \) by (b), so \( \M_{mm}^s\to\bzero \).

(d) The argument of (b) shows that \( \C(\X)\cap\C(\mathbf{D})=\{\bzero\} \), so \( \rank[\X,\mathbf{D}]=r+m \). For any \( \bb \) and \( \mathbf{g} \),
\[
\norm{\y(\bw)-\X\bb-\mathbf{D}\mathbf{g}}^2=\norm{\y_o-\X_o\bb}^2+\norm{\bw-\X_m\bb-\mathbf{g}}^2 .
\]
For each \( \bb \) the second term is minimized, at zero, by \( \mathbf{g}=\bw-\X_m\bb \). What remains is the observed-data criterion, so the minimizing
\( \bb \) are the observed-data least squares estimates, \( \X_o\hbeta=\X_o\hbeta_o \), and the minimum is \( \text{SSE}_o \). Since \( \X_m\bb \) is the same for every
least squares estimate (by @eq-dsn-missing-rank, every row of \( \X \) is an estimable function in the observed-data model), \( \hat{\bgamma}=\bw-\X_m\hbeta_o=\bw-\bz^* \).

(e) The first statement is @thm-glh-f-test for the observed-data model, of rank \( r \) with \( n-m \) observations. For the second, the complete-data
reduced-model fit of \( \y(\bz^*) \) has residual sum of squares
\( \min_{\bb}\bigl(\norm{\y_o-\X_{0,o}\bb}^2+\norm{\bz^*-\X_{0,m}\bb}^2\bigr)\ge\text{SSE}_{0,o} \), while by (a) the complete-data full-model fit has
residual sum of squares exactly \( \text{SSE}_o \). Their difference is therefore at least \( \text{SSE}_{0,o}-\text{SSE}_o \).
:::

Part (a) is the classical claim. Part (e) gives the two classical corrections. The residual sum of squares is right, but it has
\( n-m-r \) degrees of freedom, not the \( n-r \) that the complete-data formulas assign, since the filled-in values fit perfectly. And the
treatment sum of squares is biased upward, because the filled-in value suits the full model and the reduced model cannot choose its
own. Part (d) is @thm-dsn-ancova with indicator covariates: each gap gets its own parameter, the model fits it exactly and learns
nothing from it, and the indicator's coefficient measures how far the number in the gap was from the least squares value.

::: {.idea}
A missing observation can be filled in with the value that the model, fitted to the remaining data, predicts for it. Then the balanced
analysis of the completed data gives the exact least squares estimates and residual sum of squares. Subtract one error degree of freedom for
each filled-in value, and compute tests by fitting the reduced model to the observed data.
:::

## The randomized block design with one missing value

For the randomized complete block design, @eq-dsn-fill-equation can be solved by hand. The projection of the complete design is
\( \M=\bP_T+\bP_B-\bP_0 \) (@exm-proj-two-factor), so the entry of \( \M \) for two units is
\[
m_{uv}=\frac{[\text{same treatment}]}{b}+\frac{[\text{same block}]}{t}-\frac1{bt},
\]
where \( [\cdot] \) is \( 1 \) if the condition holds and \( 0 \) otherwise.

::: {#cor-dsn-yates}
[Yates's missing value formula]

In a randomized complete block design with \( t \) treatments and \( b \) blocks, suppose the observation on treatment \( i \) in block \( j \) is missing.
Let \( T' \) and \( B' \) be the totals of the observed values of treatment \( i \) and of block \( j \), and \( G' \) the total of all observed values.
The least squares estimate of the missing value is
\[
z^*=\frac{tT'+bB'-G'}{(t-1)(b-1)} .
\]{#eq-dsn-yates}

The iteration of @thm-dsn-missing(c) reduces the error by the factor \( 1/b+1/t-1/(bt) \) at every step.
:::

::: {.proof}
Here \( m=1 \), and \( \M_{mm}=1/b+1/t-1/(bt) \), so \( 1-\M_{mm}=(b-1)(t-1)/(bt) \). The row of \( \M \) belonging to the missing unit, restricted
to the observed units, sums the observed values of treatment \( i \) with weight \( 1/b-1/(bt) \), those of block \( j \) with weight \( 1/t-1/(bt) \),
and all others with weight \( -1/(bt) \). No observed unit shares both the treatment and the block. So
\( \M_{mo}\y_o=T'/b+B'/t-G'/(bt) \). Solving @eq-dsn-fill-equation,
\[
z^*=\frac{bt}{(b-1)(t-1)}\Bigl(\frac{T'}{b}+\frac{B'}{t}-\frac{G'}{bt}\Bigr)=\frac{tT'+bB'-G'}{(t-1)(b-1)} .
\]
The rate is \( \M_{mm} \) by @thm-dsn-missing(c).
:::

::: {#exm-dsn-adhesives-missing}
[A strip that split]

Suppose that in @exm-dsn-adhesives the strip from board 4 for adhesive 2 had split while being cut, so its strength (in fact
\( 8.62 \)) was never measured. The observed totals are \( T'=45.33 \) for adhesive 2,
\( B'=36.25 \) for board 4, and \( G'=263.95 \) overall. Yates's formula gives
\[
z^*=\frac{5\times45.33+6\times36.25-263.95}{4\times5}=9.010 .
\]
The iteration started from the mean of the observed values, \( 9.102 \), gives \( 9.041 \), \( 9.020 \),
\( 9.013 \), \( \dots \). The distance to \( 9.010 \) shrinks by the factor \( 1/6+1/5-1/30=1/3 \) at each step. Filling
the gap with \( 0 \) and adding an indicator column gives \( \hat{\gamma}=-9.010 \), so the covariate method recovers the same value.

The completed table gives residual sum of squares \( 1.2788 \), exactly that of the \( 29 \) observed values, but on
\( 19 \) degrees of freedom: the error mean square is \( 1.2788/19=0.0673 \), not \( 0.0639 \). Its adhesive sum of squares,
\( 4.0534 \), exceeds the correct \( 4.0512 \) by \( 0.0022 \). The correct \( F \) is \( 15.05 \) on 4 and 19 degrees
of freedom; the naive analysis would report \( 15.85 \), mostly because of the degrees of freedom.
:::

```{.python .run #cell-missing-exact}
import numpy as np
from scipy import stats

rng = np.random.default_rng(606)
t, b = 5, 6                                           # adhesives, boards
plan = np.array([rng.permutation(t) for _ in range(b)])   # which strip of each board gets which adhesive
tau = np.array([0.0, 0.40, 0.65, 0.15, 0.80])         # adhesive effects (unknown in practice)
board = rng.normal(0, 0.75, b)                        # board effects
y = np.round(9 + tau[:, None] + board[None, :] + rng.normal(0, 0.3, (t, b)), 2)

n = t * b
trt = np.repeat(np.arange(t), b)                      # y.ravel() lists adhesive 1 first
blk = np.tile(np.arange(b), t)
X = np.column_stack([np.ones(n), np.eye(t)[trt], np.eye(b)[blk]])   # complete-data design
miss = np.array([1 * b + 3])                          # adhesive 2, board 4
obs = np.setdiff1d(np.arange(n), miss)
yv = y.ravel().copy()
yv[miss] = np.nan                                     # the value is lost

def fit(Xm, yy):
    coef = np.linalg.lstsq(Xm, yy, rcond=None)[0]
    return coef, np.sum((yy - Xm @ coef) ** 2)

coef_obs, sse_obs = fit(X[obs], yv[obs])              # exact least squares, 29 observations
df_obs = len(obs) - np.linalg.matrix_rank(X[obs])     # 29 - 10 = 19
X_red = np.column_stack([np.ones(n), np.eye(b)[blk]])  # blocks only
_, sse_red = fit(X_red[obs], yv[obs])
F_exact = ((sse_red - sse_obs) / (t - 1)) / (sse_obs / df_obs)
print(f"exact: SSE {sse_obs:.3f} on {df_obs} df, adhesive SS {sse_red - sse_obs:.3f},"
      f" F = {F_exact:.2f};  fitted value for the lost cell {X[miss] @ coef_obs}")
```

```{.python .run #cell-missing-yates}
Y_obs = np.where(np.isnan(yv), 0.0, yv).reshape(t, b)
T_i = Y_obs[1].sum()                                  # observed total of adhesive 2
B_j = Y_obs[:, 3].sum()                               # observed total of board 4
G = Y_obs.sum()
z_yates = (t * T_i + b * B_j - G) / ((t - 1) * (b - 1))
print(f"Yates's value {z_yates:.3f}")

d = np.zeros(n)
d[miss] = 1.0                                         # one indicator column per missing cell
y_zero = np.where(np.isnan(yv), 0.0, yv)              # any number will do in the gap
coef_cov, sse_cov = fit(np.column_stack([X, d]), y_zero)
print(f"covariate method: gamma-hat = {coef_cov[-1]:.3f}, so the filled value is {-coef_cov[-1]:.3f};"
      f" SSE {sse_cov:.3f}")

z = np.array([np.nanmean(yv)])                        # start from the grand mean
history = [z[0]]
for step in range(30):
    y_fill = yv.copy()
    y_fill[miss] = z
    coef_fill, _ = fit(X, y_fill)                     # complete-data analysis
    z = X[miss] @ coef_fill                           # replace the gap by its fitted value
    history.append(z[0])
print("iterates", np.round(history[:6], 3))
```

With several missing values, @eq-dsn-fill-equation is an \( m\times m \) system, and the iteration converges at the rate of the largest
eigenvalue of \( \M_{mm} \). In the adhesive experiment with a second value also lost (adhesive 5 on board 1), that rate is
\( 0.367 \), a little slower than \( 1/3 \), because the two gaps are coupled through the grand mean (@exr-dsn-two-missing).

## Why the method works, and when it does not

The iteration of @thm-dsn-missing(c) is the EM algorithm for the normal linear model with missing responses: the E-step replaces each
missing value by its current fitted value, and the M-step refits. Part (b) proves convergence directly.

Everything assumes that the loss is unrelated to what the value would have been. A strip that split before gluing qualifies; a joint
that failed in handling because its glue was weak does not, and then both filling in and dropping make the adhesive look stronger than it
is. The observed data cannot reveal this. [Chapter 41](../ch41-missing-data/index.html) treats missing data in general, including the kinds of
missingness (@def-mis-mechanisms) and multiple imputation (@def-mis-multiple-imputation).

::: {.warning}
Never analyse filled-in values as if they were data without the corrections of @thm-dsn-missing(e). The residual degrees of freedom must be
reduced by one for each filled-in value, and tests must be computed by comparing models fitted to the observed data. Standard errors computed
from the completed data with the balanced formulas are also too small, because they treat \( \bz^* \) as observed.
:::

## Exercises

### A. Check your understanding

::: {#exr-dsn-yates-hand}
[A1]

Verify the arithmetic of Yates's formula in @exm-dsn-adhesives-missing, and check that the filled-in value \( 9.010 \) makes the residual of that
cell in the completed table exactly zero.
:::

::: {#exr-dsn-missing-rank-fail}
[A2]

In a randomized block design with \( t=3 \) and \( b=2 \), suppose both observations on treatment 1 are lost. Show that @eq-dsn-missing-rank fails, and
that \( \M_{mm} \) has an eigenvalue equal to \( 1 \). What happens to the iteration?
:::

::: {#exr-dsn-latin-missing}
[A3]

Derive the analogue of Yates's formula for one missing value in an \( a\times a \) Latin square:
\( z^*=\{a(R'+C'+T')-2G'\}/\{(a-1)(a-2)\} \), where \( R' \), \( C' \) and \( T' \) are the observed totals of the row, column and treatment of the
missing cell.
:::

### B. Practice

::: {#exr-dsn-yates-bias}
[B1]

In the randomized block design with one missing value, show that the treatment sum of squares computed from the completed table exceeds the
correct one by
\[
\frac{\{B'-(t-1)z^*\}^2}{t(t-1)} .
\]
Check the formula against the numbers in @exm-dsn-adhesives-missing.
:::

::: {.solution}
The correct treatment sum of squares is \( \text{SSE}_{0,o}-\text{SSE}_o \), where \( \text{SSE}_{0,o} \) is the residual sum of squares of the blocks-only
model on the observed data. The complete-data sum of squares is \( \text{SSE}_0(\y(z^*))-\text{SSE}_o \) by @thm-dsn-missing(a). So the bias is
\( \text{SSE}_0(\y(z^*))-\text{SSE}_{0,o} \). In the blocks-only model, the missing cell's fitted value from the observed data is the mean of the observed
values of its block, \( B'/(t-1) \), and by @thm-dsn-missing(a) applied to that model, filling the gap with \( B'/(t-1) \) would give
\( \text{SSE}_{0,o} \) exactly. Only the within-block sum of squares of block \( j \) depends on the filled value \( z \). It equals a constant plus
\( z^2-(B'+z)^2/t \), a quadratic in \( z \) with leading coefficient \( 1-1/t \) and minimum at \( z=B'/(t-1) \). Hence
\( \text{SSE}_0(\y(z))-\text{SSE}_{0,o}=(1-1/t)\{z-B'/(t-1)\}^2=\{(t-1)z-B'\}^2/(t(t-1)) \). With \( B'=36.25 \), \( t=5 \) and \( z^*=9.010 \): \( (36.25-36.04)^2/20=0.0022 \), as in the example.
:::

::: {#exr-dsn-two-missing}
[B2]

Suppose two values are missing in a randomized block design, in different treatments and different blocks. Write \( \M_{mm} \) explicitly and find its
eigenvalues. Check the rate \( 0.367 \) quoted for the adhesive experiment.
:::

::: {.solution}
With the two missing cells in different treatments and different blocks, they share neither, so \( \M_{mm} \) has diagonal entries
\( d=1/b+1/t-1/(bt) \) and off-diagonal entries \( -1/(bt) \). Its eigenvalues are \( d\pm1/(bt) \), and the larger is
\( 1/b+1/t=1/6+1/5=0.367 \) for \( t=5 \), \( b=6 \).
:::

### C. Going deeper

::: {#exr-dsn-missing-covariate-se}
[C1]

Show that the covariance method of @thm-dsn-missing(d) gives correct standard errors for estimable functions of \( \bbeta \), without any correction,
when the software's estimated covariance matrix of \( (\hbeta,\hat{\bgamma}) \) is used. Why does this not contradict the warning at the end of the section?
:::
