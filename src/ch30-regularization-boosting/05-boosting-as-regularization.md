# Boosting as regularization

The boosting path of [Figure 30.4.1](04-componentwise-boosting.html#fig-reg-boosting)(a) looks like a lasso path, and not by
coincidence. A close relative of boosting with tiny steps traces the lasso path exactly under a condition on the design
(Efron, Hastie, Johnstone and Tibshirani, 2004, under their positive cone condition; Hastie, Taylor, Tibshirani and Walther, 2007,
in terms of monotone paths). This section proves the simplest case, sketches the general
one, compares the procedures, and ends with boosting beyond linear least squares. Throughout, \( \y \) and the columns of \( \X \) are
centred and the columns have unit length.

## Forward stagewise regression

::: {#def-reg-stagewise}
[Incremental forward stagewise regression]

Fix \( \varepsilon>0 \). Start from \( \bb^{(0)}=\bzero \). At step \( k \), with residual \( \br^{(k)}=\y-\X\bb^{(k)} \), choose
\( j=\operatorname{arg\,max}_j\lvert\x_j\T\br^{(k)}\rvert \) and set \( b_j^{(k+1)}=b_j^{(k)}+\varepsilon\operatorname{sign}(\x_j\T\br^{(k)}) \), leaving the other
coordinates unchanged. Write FS\(_\varepsilon\) for this procedure, and FS\(_0\) for the limit of its paths as \( \varepsilon\to0 \), when the
limit exists.
:::

Stagewise and componentwise boosting choose the same coordinate; boosting moves it by \( \nu \) times its correlation, stagewise
by a fixed \( \varepsilon \), so the \( \ell_1 \) length of its path after \( k \) steps is exactly \( k\varepsilon \).

## The orthonormal case

Under an orthonormal design the lasso estimate at \( \lambda \) is \( \operatorname{sign}(z_j)(\lvert z_j\rvert-\lambda)_+ \) with
\( \mathbf{z}=\X\T\y \) (@thm-shr-lasso-orthonormal), and stagewise follows it to within one step.

::: {#prp-reg-stagewise-orthonormal}
[Stagewise under an orthonormal design]

Let \( \X\T\X=\I \), \( \mathbf{z}=\X\T\y \), and let \( \tau_k=\max_j\lvert z_j-b_j^{(k)}\rvert \) be the largest absolute residual correlation after
\( k \) steps of FS\(_\varepsilon\). Let \( K \) be the last step before \( \tau \) first falls below \( \varepsilon \), so that
\( \tau_0,\dots,\tau_K\ge\varepsilon \). For every \( k\le K \):

::: {.enumerate options="label=(\alph*)"}
1. \( \tau_{k+1}\le\tau_k \) (for \( k<K \));

2. every coordinate satisfies \( \lvert b_j^{(k)}-\operatorname{sign}(z_j)(\lvert z_j\rvert-\tau_k)_+\rvert\le\varepsilon \), that is, FS\(_\varepsilon\) is within
   \( \varepsilon \) of the lasso estimate at \( \lambda=\tau_k \).
:::

Consequently, the part of the stagewise path with \( \tau_k\ge\varepsilon \) and the lasso path are within \( 3\varepsilon \) of each other,
in the sense that every point of either is within \( 3\varepsilon \) of some point of the other in every coordinate. As \( \varepsilon\to0 \) the
stagewise path converges to the lasso path.
:::

::: {.proof}
With \( \X\T\X=\I \), the residual correlations are \( \X\T(\y-\X\bb)=\mathbf{z}-\bb \). The coordinate chosen at step \( k \) has
\( \lvert z_j-b_j^{(k)}\rvert=\tau_k\ge\varepsilon \). Moving \( b_j \) by \( \varepsilon \) towards \( z_j \) reduces this distance to \( \tau_k-\varepsilon \) without
overshooting, and the other coordinates are unchanged. This proves (a). It also shows that every move of \( b_j \) made while
\( \tau\ge\varepsilon \) goes towards \( z_j \) without passing it. Starting from \( 0 \), the coordinate therefore stays on the segment
from \( 0 \) to \( z_j \), and \( b_j^{(k)}=\operatorname{sign}(z_j)(\lvert z_j\rvert-d_j) \) with \( d_j=\lvert z_j-b_j^{(k)}\rvert\in[0,\lvert z_j\rvert] \).

Fix \( j \) and \( k \). If \( j \) has not been chosen before step \( k \), then \( b_j^{(k)}=0 \) and \( \lvert z_j\rvert=d_j\le\tau_k \), so the lasso
coordinate at \( \tau_k \) is also \( 0 \). Otherwise, let \( l<k \) be the last step at which \( j \) was chosen. Then
\( d_j=\tau_l-\varepsilon\ge\tau_k-\varepsilon \) by (a), and \( d_j\le\tau_k \) because \( \tau_k \) is a maximum. If \( \lvert z_j\rvert\ge\tau_k \), the difference from
the lasso coordinate is \( (\lvert z_j\rvert-d_j)-(\lvert z_j\rvert-\tau_k)=\tau_k-d_j\in[0,\varepsilon] \). If \( \lvert z_j\rvert<\tau_k \), the lasso coordinate is \( 0 \)
and \( \lvert b_j^{(k)}\rvert=\lvert z_j\rvert-d_j\le\tau_k-d_j\le\varepsilon \). This proves (b).

For the last statement, (b) puts every stagewise point with \( \tau_k\ge\varepsilon \) within \( \varepsilon \) of the lasso path. Conversely,
with \( K \) as in the statement, each step lowers \( \tau \) by at most \( \varepsilon \), because
\( \tau_{k+1}\ge d_j=\tau_k-\varepsilon \) for the chosen \( j \). Hence \( \tau_K<2\varepsilon \), and every \( \lambda\in[\tau_K,\tau_0] \) is within \( \varepsilon \) of
some \( \tau_k \) with \( k\le K \). Soft thresholding is \( 1 \)-Lipschitz in \( \lambda \), so for such \( \lambda \) the lasso estimate is within
\( 2\varepsilon \) of \( \bb^{(k)} \). For \( \lambda\ge\tau_0=\max_j\lvert z_j\rvert \) it is \( \bzero=\bb^{(0)} \). For \( \lambda<\tau_K \) it is within
\( \tau_K-\lambda<2\varepsilon \) of the lasso estimate at \( \tau_K \), hence within \( 3\varepsilon \) of \( \bb^{(K)} \).
:::

For componentwise boosting the same argument gives tolerance \( \nu\tau_k \) (@exr-reg-boost-soft).

## The general relation

In general the lasso path is piecewise linear with explicit directions. Part (a) below restates the segment formula of
@prp-shr-lasso-path in the notation of this chapter; what is new is the reading in terms of tied correlations.

::: {#prp-reg-stagewise}
[Stagewise and the lasso path]

::: {.enumerate options="label=(\alph*)"}
1. Suppose that for \( \lambda \) in an interval the lasso estimate is unique and has a fixed active set \( A \) and fixed signs \( \mathbf{s}_A \), and \( \X_A \) has
   full column rank. Then on that interval
   \[
   \hbeta_A(\lambda)=(\X_A\T\X_A)^{-1}(\X_A\T\y-\lambda\mathbf{s}_A),
   \]{#eq-reg-lasso-segment}

   so as \( \lambda \) decreases the active coefficients move along the direction \( \mathbf{d}_A=(\X_A\T\X_A)^{-1}\mathbf{s}_A \), and all active
   variables keep the common absolute correlation \( \lambda \) with the residual.

2. *(Sketch; Efron, Hastie, Johnstone and Tibshirani, 2004, who prove it under their positive cone condition on the design; Hastie et al., 2007, for monotone paths.)* If along the whole lasso path every direction \( \mathbf{d}_A \) has the signs
   \( \mathbf{s}_A \) coordinatewise, that is, if every lasso coefficient is monotone in \( \lambda \), then FS\(_0\) exists and coincides with the
   lasso path.

3. *(Sketch; Hastie, Taylor, Tibshirani and Walther, 2007.)* In general FS\(_0\) is the path of the *monotone lasso*. Write
   \( \bb=\bb^+-\bb^- \) with \( \bb^\pm\ge\bzero \). Along the path, every coordinate of \( \bb^+ \) and of \( \bb^- \) is nondecreasing, and the path
   lowers the residual sum of squares as fast as possible per unit of \( \ell_1 \) arc length, subject to that constraint.
:::

:::

::: {.proof}
(a) The formula @eq-reg-lasso-segment is @prp-shr-lasso-path. By @lem-reg-optimality and @lem-reg-norm-subgradient, the active
coordinates satisfy \( \X_A\T(\y-\X_A\hbeta_A)=\lambda\mathbf{s}_A \), which says that the active correlations are \( \lambda\mathbf{s}_A \).

*Sketch of (b) and (c).* Complete proofs need a careful limiting argument; see the papers cited. Choosing a variable of
maximal absolute correlation lowers its correlation relative to the others, so in the limit the active correlations are
tied, as for the lasso. The net movement is a direction \( \mathbf{d} \) with
\( d_js_j\ge0 \) on \( A \), and the only direction keeping the correlations tied is \( \mathbf{d}_A \). If \( \mathbf{d}_A \) has the signs
\( \mathbf{s}_A \), stagewise follows it; that is (b). Otherwise the lasso shrinks a coefficient against the sign of its correlation, which
stagewise cannot do. Stagewise instead holds it and moves along a nonnegative least squares fit of the residual on the
sign-adjusted active columns. Efron et al. (2004) build this into a modified least angle regression and prove it gives FS\(_0\);
Hastie et al. (2007) identify the path with the monotone lasso.
:::

The lasso at \( \ell_1 \) norm \( t \) minimizes the residual sum of squares subject to \( \norm{\bb}_1\le t \) (the constrained form of @def-shr-lasso);
stagewise budgets the \( \ell_1 \) arc length and never undoes a move. For monotone paths the two budgets coincide (@exr-reg-arc-length).

::: {#exm-reg-stagewise}
[Stagewise and lasso paths]

*Orthonormal design.* With \( 5 \) orthonormal columns and \( \varepsilon=0.01 \), the largest distance from the lasso estimate at
\( \lambda=\tau_k \) over all steps is \( 0.0092 \).

*A non-monotone lasso path.* For \( n=40 \) observations of \( p=6 \) regressors sharing a common factor (correlations near 0.64),
the lasso coefficient of regressor \( 6 \) grows to \( 0.308 \) at \( \norm{\bb}_1=4.12 \) and then shrinks back. FS\(_\varepsilon\) with
\( \varepsilon=0.002 \) agrees with the lasso to within \( 0.039 \) until then, and afterwards holds the coefficient: halfway to the
least squares norm \( 7.53 \), the two differ by \( 0.165 \) there ([Figure 30.5.1](#fig-reg-stagewise)). At every \( \ell_1 \) norm
checked, the lasso has the smaller residual sum of squares, as it must; the excess of stagewise is at most \( 0.011 \).

*Prediction.* Over \( 50 \) data sets from the design of @exm-reg-boosting, the average smallest risk along the path is
\( 2.857 \) for the lasso, \( 2.857 \) for stagewise and \( 2.856 \) for boosting with \( \nu=0.05 \), against \( 4.607 \) for least
squares. Data set by data set the points lie close to the diagonal (panel (c)): the best risks of stagewise and boosting never
differ from the lasso's by more than \( 6.6 \) per cent.
:::

::: {when-format="html"}
![**Figure 30.5.1.** (a) The lasso path and (b) the forward stagewise path (\( \varepsilon=0.002 \)) for the correlated design of
@exm-reg-stagewise, against \( \norm{\bb}_1 \). The thick purple curve is the coefficient that the lasso shrinks back after the dashed
line; stagewise holds it. (c) Smallest risk along the path for 50 simulated data sets: forward stagewise and boosting
(\( \nu=0.05 \)) against the lasso.](stagewise_lasso.svg){#fig-reg-stagewise width=100%}
:::

::: {when-format="pdf"}
![(a) The lasso path and (b) the forward stagewise path (\( \varepsilon=0.002 \)) for the correlated design of
@exm-reg-stagewise, against \( \norm{\bb}_1 \). The thick purple curve is the coefficient that the lasso shrinks back after the dashed
line; stagewise holds it. (c) Smallest risk along the path for 50 simulated data sets: forward stagewise and boosting
(\( \nu=0.05 \)) against the lasso.](stagewise_lasso.pdf){width=100%}
:::

```{.python .run #cell-stagewise-compare}
import numpy as np

def soft(z, lam):
    return np.sign(z) * np.maximum(np.abs(z) - lam, 0.0)

def lasso_cd(X, y, lam, b=None, tol=1e-12):
    G, c = X.T @ X, X.T @ y
    b = np.zeros(X.shape[1]) if b is None else b.copy()
    while True:
        largest = 0.0
        for j in range(len(b)):
            old = b[j]
            b[j] = soft(c[j] - G[j] @ b + G[j, j] * old, lam) / G[j, j]
            largest = max(largest, abs(b[j] - old))
        if largest < tol:
            return b

def forward_stagewise(X, y, eps, steps):
    """Incremental forward stagewise: move the most correlated coefficient by eps."""
    b = np.zeros(X.shape[1])
    r = y.astype(float).copy()
    path = [b.copy()]
    for _ in range(steps):
        c = X.T @ r
        j = np.argmax(np.abs(c))
        delta = eps * np.sign(c[j])
        b[j] += delta
        r -= delta * X[:, j]
        path.append(b.copy())
    return np.array(path)

rng = np.random.default_rng(4)
n, p = 40, 6
z = rng.normal(size=(n, 1))
X_raw = 0.8 * z + 0.6 * rng.normal(size=(n, p))        # a shared factor: correlations near 0.64
center = X_raw.mean(axis=0)
scale = np.linalg.norm(X_raw - center, axis=0)
X = (X_raw - center) / scale                           # centred, unit-length columns
beta = np.array([3.0, -2.0, 2.0, 0.0, 0.0, -1.0])
y = X @ beta + 0.5 * rng.normal(size=n)
y_bar = y.mean()
y = y - y_bar

lam_max = np.max(np.abs(X.T @ y))
grid = lam_max * np.logspace(0, -3, 400)
lasso_path = []
b = np.zeros(p)
for lam in grid:
    b = lasso_cd(X, y, lam, b)
    lasso_path.append(b.copy())
lasso_path = np.array(lasso_path)

eps = 0.002
fs_path = forward_stagewise(X, y, eps, 6000)
l1_fs = np.abs(fs_path).sum(axis=1)
arc_fs = eps * np.arange(len(fs_path))                  # total variation of the stagewise path
print("lasso at the end:     ", np.round(lasso_path[-1], 3))
print("stagewise at the end: ", np.round(fs_path[-1], 3))
```

Near least squares, stagewise dithers, stepping back and forth by \( \varepsilon \): its arc length grows while its \( \ell_1 \) norm does not.

## Beyond linear least squares

For a differentiable loss \( \sum_iL(y_i,f_i) \), boosting fits the negative gradient \( -\partial L(y_i,f_i)/\partial f_i \) in
place of the residual (Friedman 2001): absolute error gives residual signs, the Huber loss clipped residuals, and the
exponential loss \( e^{-yf} \), \( y=\pm1 \), recovers AdaBoost (Friedman, Hastie and Tibshirani, 2000). Small regression trees as
base learners give gradient tree boosting; penalized splines in one regressor at a time give componentwise boosting of an
additive model (Bühlmann and Yu, 2003; Bühlmann and Hothorn, 2007), a topic for Chapters [43](../ch43-smoothing/index.html) and [44](../ch44-additive-models/index.html) (@def-smo-pspline, @def-add-model).

Zhang and Yu (2005) prove consistency of early-stopped boosting for general convex losses over rich base classes; Bühlmann
(2006) proves it for componentwise \( L_2 \) boosting with \( p \) growing almost exponentially in \( n \) under sparsity. We state
this without proof; it says that some sequence of stopping iterations works, not how to find it.

::: {.warning}
Early stopping is an implicit regularizer: the estimates minimize no stated criterion. For prediction, tune \( m \) like any
smoothing parameter. For interpreting coefficients, boosting inherits every caution of
[Chapter 29](../ch29-model-selection/index.html) about selection, and inference after boosting meets the obstacles of inference
after the lasso (@thm-hd-debiased).
:::

## Exercises

### A. Check your understanding

::: {#exr-reg-fs-hand}
[A1]

Let \( \X\T\X=\I \), \( \mathbf{z}=(3,-1,0.5)\T \) and \( \varepsilon=0.5 \). Run FS\(_\varepsilon\) by hand until the largest absolute residual correlation
first falls below \( 1 \), and check @prp-reg-stagewise-orthonormal(b) at each step.
:::

::: {#exr-reg-arc-length}
[A2]

Show that the \( \ell_1 \) arc length of the FS\(_\varepsilon\) path after \( k \) steps is \( k\varepsilon \), that it is at least \( \norm{\bb^{(k)}}_1 \), and
that equality holds iff no coordinate has ever moved in both directions. Which of the two paths in [Figure 30.5.1](#fig-reg-stagewise)
has equality up to the least squares region?
:::

### B. Practice

::: {#exr-reg-boost-soft}
[B1]

Adapt the proof of @prp-reg-stagewise-orthonormal to componentwise boosting with step \( \nu\in(0,1] \) under an orthonormal design. Show that
\( \tau_k \) is nonincreasing and that \( \lvert b_j^{(k)}-\operatorname{sign}(z_j)(\lvert z_j\rvert-\tau_k)_+\rvert\le\nu\tau_k \) for all \( j \) and \( k \).
:::

::: {.solution}
The chosen coordinate's remaining correlation changes from \( \tau_k \) to \( (1-\nu)\tau_k \), with no overshoot since \( \nu\le1 \). So \( \tau \) is
nonincreasing, and each coordinate stays between \( 0 \) and \( z_j \). If \( j \) was last chosen at step \( l<k \), then
\( d_j=(1-\nu)\tau_l\ge(1-\nu)\tau_k \) and \( d_j\le\tau_k \). As in the proof, the distance from the lasso coordinate at \( \tau_k \) is at most
\( \tau_k-d_j\le\nu\tau_k \). A coordinate that has never been chosen is \( 0 \), and so is the lasso coordinate.
:::

::: {#exr-reg-lars-breakpoints}
[B2]

On a segment of the lasso path with active set \( A \) and signs \( \mathbf{s}_A \), show that the correlation of an inactive regressor is
\( c_j(\lambda)=\x_j\T(\I-\M_A)\y+\lambda\,\x_j\T\X_A(\X_A\T\X_A)^{-1}\mathbf{s}_A \), where \( \M_A \) projects onto \( \C(\X_A) \). Deduce the value of \( \lambda \) at
which \( j \) joins the active set, and the value at which an active coefficient reaches zero. This is the lasso version of least angle
regression.
:::

::: {.solution}
By @eq-reg-lasso-segment, \( \y-\X_A\hbeta_A(\lambda)=(\I-\M_A)\y+\lambda\X_A(\X_A\T\X_A)^{-1}\mathbf{s}_A \). Take the inner product with \( \x_j \). Regressor \( j \)
joins when \( \lvert c_j(\lambda)\rvert \) reaches \( \lambda \). Since \( c_j \) is affine in \( \lambda \), this happens at the largest \( \lambda \) below the current one
solving \( c_j(\lambda)=\pm\lambda \). An active coefficient \( \hat\beta_k(\lambda) \) is affine in \( \lambda \) by @eq-reg-lasso-segment, and it leaves when it
reaches \( 0 \). The next breakpoint is the largest of these candidate values. At a breakpoint, \( A \) and \( \mathbf{s}_A \) are updated, and the
path continues linearly.
:::

::: {#exr-reg-fs-two}
[B3]

Show that for \( p=2 \) (unit-length columns, correlation \( \rho\in(-1,1) \)) the direction \( \mathbf{d}_A \) always has the signs \( \mathbf{s}_A \). Conclude,
using @prp-reg-stagewise(b), that FS\(_0\) and the lasso path coincide for every data set with two regressors.
:::

::: {.solution}
For \( A=\{j\} \), \( \mathbf{d}_A=s_j \). For \( A=\{1,2\} \), \( (\X_A\T\X_A)^{-1}=(1-\rho^2)^{-1}\begin{pmatrix}1&-\rho\\-\rho&1\end{pmatrix} \). So
\( \mathbf{d}_A=(1-\rho^2)^{-1}(s_1-\rho s_2,\,s_2-\rho s_1)\T \). If \( s_1=s_2=s \), this is \( s(1-\rho)(1-\rho^2)^{-1}(1,1)\T \), with the signs of
\( s \). If \( s_2=-s_1 \), it is \( (1+\rho)(1-\rho^2)^{-1}(s_1,s_2)\T \). In both cases the signs are right, because \( \lvert\rho\rvert<1 \).
:::
