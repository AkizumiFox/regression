# Robust fitting as a diagnostic

Least squares hides bad cases by accommodating them. A fit that refuses to accommodate them leaves them with large
residuals, visible without any deletion. That is the diagnostic use of robust regression, and this section develops only
as much theory as that use needs: does a fit the bad cases cannot capture tell a different story, and which cases
account for the difference?

## M-estimators

Least squares minimizes \( \sum_i\rho(y_i-\x_{(i)}\T\bb) \) with \( \rho(u)=u^2/2 \), which lets one large residual outweigh many
small ones. An M-estimator (Huber 1964) uses a \( \rho \) that grows more slowly.

::: {#def-res-m-estimator}
[M-estimator of regression]

Let \( \rho:\Real\to[0,\infty) \) be even, nondecreasing on \( [0,\infty) \), differentiable, with \( \rho(0)=0 \), and put \( \psi=\rho' \) and
\( w(u)=\psi(u)/u \) (with \( w(0)=\lim_{u\to0}\psi(u)/u \)). Given a scale \( \hat{\sigma}>0 \), an **M-estimate** \( \hbeta_M \) minimizes
\[
Q(\bb)=\sum_{i=1}^n\rho\Bigl(\frac{y_i-\x_{(i)}\T\bb}{\hat{\sigma}}\Bigr).
\]
Two choices of \( \rho \) are standard.

::: {.enumerate options="label=(\alph*)"}
1. **Huber:** \( \rho(u)=u^2/2 \) for \( \lvert u\rvert\le k \) and \( k\lvert u\rvert-k^2/2 \) beyond, so \( \psi(u)=\max(-k,\min(k,u)) \) and
   \( w(u)=\min(1,k/\lvert u\rvert) \).

2. **Bisquare** (Tukey's biweight): \( \rho(u)=\tfrac{c^2}6\bigl[1-(1-u^2/c^2)^3\bigr] \) for \( \lvert u\rvert\le c \) and \( c^2/6 \) beyond, so
   \( \psi(u)=u(1-u^2/c^2)^2 \) and \( w(u)=(1-u^2/c^2)^2 \) for \( \lvert u\rvert<c \), and \( \psi(u)=w(u)=0 \) beyond.
:::

:::

Setting the gradient of \( Q \) to zero gives the estimating equations
\[
\sum_{i=1}^n\psi\Bigl(\frac{y_i-\x_{(i)}\T\bb}{\hat{\sigma}}\Bigr)\x_{(i)}=\bzero
\quad\Longleftrightarrow\quad
\X\T\W(\bb)(\y-\X\bb)=\bzero ,
\]{#eq-res-m-equations}

where \( \W(\bb) \) is diagonal with entries \( w\bigl((y_i-\x_{(i)}\T\bb)/\hat{\sigma}\bigr) \): the normal equations of weighted least
squares ([Section 6.9](../ch06-projections/09-inner-products.html)) with weights that depend on the fit. A residual beyond
\( k\hat{\sigma} \) has its Huber influence capped; a residual beyond \( c\hat{\sigma} \) gets bisquare weight zero. The Huber \( \rho \) is
convex; the bisquare \( \rho \) is not, and is called **redescending** because \( \psi \) falls back to zero.

A standard scale is \( \hat{\sigma}=\operatorname{median}_i\lvert y_i-\x_{(i)}\T\bb\rvert/0.6745 \), where \( 0.6745=\Phi^{-1}(0.75) \) makes it
consistent for \( \sigma \) under centred normal errors; software re-estimates it at each step. The constants \( k=1.345 \) and
\( c=4.685 \) make the location efficiency \( (\E\psi'(Z))^2/\E\psi(Z)^2 \), \( Z\sim\Normal(0,1) \), equal to \( 0.950 \) for both (computed
by numerical integration in the script). For convex \( \rho \) such as Huber's, with \( p \) fixed and under smoothness conditions on \( \psi \), that this ratio governs
the asymptotic variance of regression M-estimators when \( \max_ih_{ii}\to0 \) is a theorem of Huber (1973) that this book does not
prove; redescending estimators such as the bisquare also need a consistent starting value, as in the MM-estimator of Yohai (1987).

## Iteratively reweighted least squares

@eq-res-m-equations suggests **iteratively reweighted least squares** (IRLS), studied for robust regression by Holland and
Welsch (1977): from \( \bb^{(t)} \), compute \( w_i=w\bigl((y_i-\x_{(i)}\T\bb^{(t)})/\hat{\sigma}\bigr) \) and let \( \bb^{(t+1)} \) be the weighted least
squares estimate; repeat. With the scale fixed, it never goes uphill.

::: {#prp-res-irls}
[IRLS descends]

Suppose \( \rho(u)=g(u^2) \), where \( g:[0,\infty)\to[0,\infty) \) is concave, nondecreasing, and differentiable on \( (0,\infty) \) with a
finite limit \( g'(0^+) \). Then \( w(u)=2g'(u^2)\ge0 \). Hold \( \hat{\sigma} \) fixed, and suppose \( \X\T\W\X \) is nonsingular at every step.
Then every IRLS step satisfies \( Q(\bb^{(t+1)})\le Q(\bb^{(t)}) \). The Huber and bisquare \( \rho \) satisfy the hypotheses.
:::

::: {.proof}
Since \( \rho(u)=g(u^2) \), \( \psi(u)=2ug'(u^2) \), so \( w(u)=2g'(u^2) \), which is nonnegative because \( g \) is nondecreasing. A concave
function lies below each of its tangent lines: \( g(v)\le g(v_0)+g'(v_0)(v-v_0) \) for \( v\ge0 \), \( v_0>0 \) (and for \( v_0=0 \) with
\( g'(0^+) \), by continuity). With \( v=u^2 \) and \( v_0=u_0^2 \),
\[
\rho(u)\le\rho(u_0)+\tfrac12w(u_0)\,(u^2-u_0^2).
\]
Write \( u_i(\bb)=(y_i-\x_{(i)}\T\bb)/\hat{\sigma} \) and \( w_i=w\bigl(u_i(\bb^{(t)})\bigr) \). Summing over \( i \),
\[
Q(\bb)\le G(\bb)=Q(\bb^{(t)})+\tfrac12\sum_iw_i\bigl(u_i(\bb)^2-u_i(\bb^{(t)})^2\bigr)\quad\text{for every }\bb,
\]
with equality at \( \bb=\bb^{(t)} \). Up to terms that do not involve \( \bb \), \( G(\bb) \) is \( \frac1{2\hat{\sigma}^2}\sum_iw_i(y_i-\x_{(i)}\T\bb)^2 \), which the
weighted least squares step minimizes. Hence \( Q(\bb^{(t+1)})\le G(\bb^{(t+1)})\le G(\bb^{(t)})=Q(\bb^{(t)}) \).

For Huber, \( g(v)=v/2 \) for \( v\le k^2 \) and \( k\sqrt v-k^2/2 \) beyond, so \( g'(v)=\tfrac12 \) and then \( k/(2\sqrt v) \): continuous
at \( v=k^2 \), nonincreasing, hence \( g \) is concave. For the bisquare, \( g(v)=\tfrac{c^2}6[1-(1-v/c^2)^3] \) for \( v\le c^2 \) and constant
beyond, so \( g'(v)=\tfrac12(1-v/c^2)^2 \), falling to \( 0 \) at \( v=c^2 \), and \( 0 \) beyond: again nonincreasing and continuous.
:::

This is majorization: \( Q \) is replaced by a quadratic that lies above it and touches it at the current point. For Huber's
convex \( \rho \), a fixed point of IRLS solves @eq-res-m-equations and is a global minimizer. For the bisquare, \( Q \) can have
several local minima, and IRLS finds one near its start.

```{.python .run #cell-robust-fits-irls}
import numpy as np

def huber_weight(u, k=1.345):
    return np.minimum(1.0, k / np.maximum(np.abs(u), 1e-12))

def bisquare_weight(u, c=4.685):
    return np.where(np.abs(u) < c, (1 - (u / c) ** 2) ** 2, 0.0)

def mad_scale(r):
    return np.median(np.abs(r)) / 0.6745           # consistent for sigma at the normal

def irls(X, y, weight, beta=None, scale=None, tol=1e-10, max_iter=500):
    """M-estimate by iteratively reweighted least squares.

    Starts from least squares unless `beta` is given. If `scale` is None it is re-estimated
    from the current residuals at every step; otherwise it is held fixed."""
    if beta is None:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for _ in range(max_iter):
        r = y - X @ beta
        s = mad_scale(r) if scale is None else scale
        w = weight(r / s)
        sw = np.sqrt(w)
        new, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        if np.max(np.abs(new - beta)) < tol * (1 + np.max(np.abs(beta))):
            beta = new
            break
        beta = new
    r = y - X @ beta
    s = mad_scale(r) if scale is None else scale
    return beta, weight(r / s), s
```

## Breakdown

How many bad cases can an estimator tolerate? The finite-sample breakdown point of @def-dep-breakdown counts the cases
that must be replaced to make the estimate arbitrary, and least squares breaks down at a single case: changing \( y_i \) by
\( \Delta \) moves \( \hbeta \) by \( (\X\T\X)^{-1}\x_{(i)}\Delta \) (@prp-dep-breakdown).

Provided no single row dominates the design, a Huber estimate cannot be broken by changing \( y_i \) alone, since \( \psi \) is bounded; how much a \( y \)-outlier can do depends on the leverages. But one case bad in both \( \x \) and \( y \), a
**bad leverage point**, breaks it: the fit turns towards the case, whose residual becomes moderate and whose large row then
dominates the estimating equations (@exr-res-huber-breakdown). This is masking again: M-estimates downweight by residual,
and a high-leverage case controls its own residual. For regression M-estimates with monotone \( \psi \) the breakdown point is
\( 1/n \) (Maronna, Martin and Yohai 2006). Added to the thirty clean cases of @exm-res-masking, one case at \( x=50 \), \( y=0 \)
turns the least squares slope from \( 0.514 \) to \( -0.025 \), and the Huber slope to \( -0.026 \).

## High-breakdown fits

Rousseeuw (1984) proposed estimators that fit a majority of the data and ignore the rest.

::: {#def-res-lts}
[Least trimmed squares]

Let \( q \) be an integer with \( n/2\le q\le n \). For \( \bb\in\Real^p \), let \( r_{[1]}^2(\bb)\le\dots\le r_{[n]}^2(\bb) \) be the ordered
squared residuals \( (y_i-\x_{(i)}\T\bb)^2 \). The **least trimmed squares** (LTS) estimate minimizes
\[
T(\bb)=\sum_{k=1}^qr_{[k]}^2(\bb),
\]
the sum of the \( q \) smallest squared residuals. The **least median of squares** (LMS) estimate minimizes \( r_{[\lceil n/2\rceil]}^2(\bb) \)
instead.
:::

With \( q=\lfloor(n+p+1)/2\rfloor \) the breakdown point of LTS is close to \( 1/2 \) (see Rousseeuw and Leroy 1987): the fit follows
the best-fitting half of the data. Computing it is combinatorial, since \( \min T \) is the minimum over \( q \)-subsets of their least
squares criteria (@exr-res-lts-subsets). The practical algorithm rests on the **concentration step** of Rousseeuw and Van
Driessen (2006).

::: {#prp-res-cstep}
[Concentration step]

Let \( \bb\in\Real^p \), let \( H \) be a set of \( q \) cases with the smallest squared residuals at \( \bb \), and suppose the rows of \( \X \) in \( H \)
have full column rank. Let \( \bb' \) be the least squares estimate computed from the cases in \( H \). Then \( T(\bb')\le T(\bb) \).
:::

::: {.proof}
By the choice of \( H \), \( T(\bb)=\sum_{i\in H}(y_i-\x_{(i)}\T\bb)^2 \), and since \( \bb' \) minimizes the sum of squares over \( H \),
this is at least \( \sum_{i\in H}(y_i-\x_{(i)}\T\bb')^2 \). Any \( q \) squared residuals at \( \bb' \) sum to at least the \( q \) smallest,
\( T(\bb') \).
:::

Repeated concentration steps give a nonincreasing sequence of values of \( T \), constant after finitely many steps since there
are finitely many \( q \)-subsets. The algorithm starts from many random **elemental** fits through \( p \) cases and keeps the best
result. It needs only one start whose \( p \) cases are clean: with a fraction \( \epsilon \) of bad cases, \( N \) starts all fail with
probability about \( \bigl(1-(1-\epsilon)^p\bigr)^N \), and with \( \epsilon=0.2 \), \( p=4 \), \( N=9 \) starts push this below \( 1\% \) (@exr-res-lts-starts). LTS with \( q\approx n/2 \) is inefficient; the MM-estimator of Yohai (1987) uses a high-breakdown fit only to
start a redescending M-estimate with a fixed robust scale, keeping the high breakdown point and recovering efficiency. A
simplified version, bisquare IRLS from LTS with the scale fixed at the median absolute LTS residual over \( 0.6745 \), is used
below.

```{.python .run #cell-robust-fits-lts}
def lts(X, y, q, n_starts=500, seed=0):
    """Least trimmed squares: minimize the sum of the q smallest squared residuals.

    Random elemental starts (p cases), each followed by concentration steps: fit least
    squares to the current q cases, then keep the q cases with the smallest residuals."""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    best, best_obj = None, np.inf
    for _ in range(n_starts):
        idx = rng.choice(n, size=p, replace=False)
        if np.linalg.matrix_rank(X[idx]) < p:
            continue
        beta = np.linalg.solve(X[idx], y[idx])
        for _ in range(100):                       # concentration steps
            keep = np.argsort((y - X @ beta) ** 2)[:q]
            new, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            if np.allclose(new, beta):
                break
            beta = new
        obj = np.sort((y - X @ beta) ** 2)[:q].sum()
        if obj < best_obj:
            best, best_obj = beta, obj
    return best, best_obj
```

## Using a robust fit diagnostically

The procedure: (1) fit least squares and a high-breakdown or MM-type fit; (2) compare the coefficients, since a difference of
more than a standard error means a minority of cases drives least squares; (3) standardize the robust residuals by the robust
scale and list cases beyond about \( 2.5 \) or with small weights, plotting them against a robust leverage measure to separate
outlying responses from bad leverage points (Rousseeuw and van Zomeren 1990); (4) examine the flagged set substantively, refit
without it, and report both fits. @thm-res-group-deletion(d) measures its incompatibility with the rest, with an optimistic
\( p \)-value because the set was chosen from the data.

::: {#exm-res-robust-masking}
[Robust fits on the masking data]

For the data of @exm-res-masking ([Figure 20.6.1](#fig-res-robust)(a)), least squares has slope \( 0.007 \), and the Huber fit
\( 0.005 \), with Huber weight \( 1.000 \) for each bad record: they have captured the fit so completely that their residuals are
not even downweighted. The bisquare fit started from least squares is captured too (slope \( 0.005 \)). LTS with \( q=18 \) of the
\( 34 \) cases has slope \( 0.546 \) and robust scale \( 0.587 \); in that scale each bad record has standardized residual \( -18.6 \),
and no clean case exceeds \( 2.81 \). The bisquare fit started from LTS has slope \( 0.509 \) and intercept \( 1.930 \), close to the
clean fit (\( 0.514 \) and \( 1.862 \)), and gives the four records weight zero.
:::

```{.python .run #cell-robust-fits-masking}
rng = np.random.default_rng(20)
n0, m = 30, 4
x_clean = np.sort(rng.uniform(0, 10, n0))
y_clean = 2 + 0.5 * x_clean + rng.normal(0, 0.5, n0)
x = np.r_[x_clean, np.full(m, 20.0)]
y = np.r_[y_clean, np.full(m, 2.0)]
X = np.column_stack([np.ones(len(x)), x])
n, p = X.shape

b_ls, *_ = np.linalg.lstsq(X, y, rcond=None)
b_hub, w_hub, _ = irls(X, y, huber_weight)
b_lts, _ = lts(X, y, q=(n + p + 1) // 2)
s_lts = mad_scale(y - X @ b_lts)
b_bis, w_bis, _ = irls(X, y, bisquare_weight, beta=b_lts, scale=s_lts)
for name, b in [("least squares", b_ls), ("Huber", b_hub), ("LTS", b_lts), ("bisquare from LTS", b_bis)]:
    print(f"{name:18s} intercept {b[0]:6.3f}  slope {b[1]:6.3f}")
print("scaled LTS residuals of the cluster:", ((y - X @ b_lts) / s_lts)[-m:].round(1))
```

::: {#exm-res-engel-robust}
[Engel's budgets, robustly]

For food expenditure on income in Engel's data (@exm-res-engel-leverage), least squares gives slope \( 0.485 \), Huber
\( 0.537 \) and bisquare \( 0.553 \). The bisquare fit gives weight zero to \( 3 \) households, rows \( 59, 105, 138 \)
([Figure 20.6.1](#fig-res-robust)(b)). Row \( 138 \), the case of highest leverage (income \( 4958 \), food \( 1827 \), far below the
line through the others, Huber weight \( 0.117 \)), pulls the least squares slope down; in the least squares fit its
externally studentized residual is \( -8.39 \) and its Cook's distance \( 9.28 \).

But the spread of food expenditure grows with income. On the logarithmic scale the three fits nearly agree (slopes \( 0.856 \),
\( 0.868 \), \( 0.871 \)) and the smallest bisquare weight is \( 0.16 \): nobody is rejected. The robust fit correctly reported that a
few cases disagree with a constant-variance line; whether the data or the scale is wrong no fit can decide
([Chapter 21](../ch21-nonnormality-heteroscedasticity-serial/index.html), [Chapter 22](../ch22-transformations/index.html)).
:::

```{.python .run #cell-robust-fits-engel}
import statsmodels.api as sm
engel = sm.datasets.engel.load_pandas().data
inc = engel["income"].to_numpy()
food = engel["foodexp"].to_numpy()
XE = np.column_stack([np.ones(len(inc)), inc])
fits = {
    "least squares": sm.OLS(food, XE).fit(),
    "Huber": sm.RLM(food, XE, M=sm.robust.norms.HuberT()).fit(),
    "bisquare": sm.RLM(food, XE, M=sm.robust.norms.TukeyBiweight()).fit(),
}
for name, f in fits.items():
    print(f"{name:14s} slope {f.params[1]:.4f}")
w = fits["bisquare"].weights
print("households with bisquare weight 0 (row numbers from 1):", np.where(w == 0)[0] + 1)
```

::: {when-format="html"}
![**Figure 20.6.1.** Robust fits. (a) The masking data: least squares and Huber are captured, LTS and bisquare from LTS
are not. (b) Engel's budgets; circled points have bisquare weight zero.](robust_fits.svg){#fig-res-robust width=100%}
:::

::: {when-format="pdf"}
![Robust fits. (a) The masking data: least squares and Huber are captured, LTS and bisquare from LTS
are not. (b) Engel's budgets; circled points have bisquare weight zero.](robust_fits.pdf){width=100%}
:::

::: {.warning}
Deleting flagged cases and then reporting ordinary standard errors treats a data-driven choice as fixed in advance, and the
intervals are too short. Report both analyses; [Chapter 23](../ch23-resampling-inference/index.html) discusses resampling that
can include the selection step.
:::

## Exercises

### A. Check your understanding

::: {#exr-res-huber-weights}
[A1]

With \( k=1.345 \), compute the Huber weights of standardized residuals \( 0.5 \), \( 2 \) and \( 10 \). With \( c=4.685 \), compute the bisquare
weights of the same residuals.
:::


::: {#exr-res-ls-one-step}
[A2]

Show that least squares is the M-estimator with \( \rho(u)=u^2/2 \), that its weights are all one, and that IRLS started anywhere
reaches it in one step.
:::

### B. Practice

::: {#exr-res-lts-subsets}
[B1]

Show that \( \min_{\bb}T(\bb)=\min_H\text{SSE}_H \), where \( H \) ranges over the \( q \)-subsets of cases for which the rows of \( \X \) in \( H \)
have full column rank and \( \text{SSE}_H \) is the residual sum of squares of least squares on \( H \) (assume at least one such subset
exists). Deduce that an LTS estimate is the least squares estimate on some \( q \)-subset.
:::

::: {.solution}
For any \( \bb \) and any \( q \)-subset \( H \), \( T(\bb)\le\sum_{i\in H}(y_i-\x_{(i)}\T\bb)^2 \), with equality when \( H \) consists of cases with
the \( q \) smallest squared residuals. So \( T(\bb)=\min_H\sum_{i\in H}(y_i-\x_{(i)}\T\bb)^2 \), and
\( \min_{\bb}T(\bb)=\min_H\min_{\bb}\sum_{i\in H}(y_i-\x_{(i)}\T\bb)^2 \), the order of minimization being interchangeable. For a subset of full
rank the inner minimum is \( \text{SSE}_H \), attained at the least squares estimate on \( H \). (Subsets without full rank can be handled by a
limiting argument or excluded, as the statement assumes.) An LTS estimate is therefore the least squares estimate on a minimizing \( H \).
:::

::: {#exr-res-lts-starts}
[B2]

Find the smallest number \( N \) of elemental starts for which \( \bigl(1-(1-\epsilon)^p\bigr)^N\le0.01 \) when \( \epsilon=0.2 \) and \( p=4 \), and
when \( \epsilon=0.4 \) and \( p=10 \). What does the second answer say about high-breakdown fitting with many regressors?
:::


::: {#exr-res-robust-weights-plot}
[B3]

For the state murder-rate regression with all \( 51 \) jurisdictions, fit Huber and bisquare M-estimates with statsmodels and list the
jurisdictions with the smallest weights. Compare the poverty coefficient with the least squares values with and without the District of
Columbia (@exm-res-state-influence).
:::

### C. Going deeper

::: {#exr-res-huber-breakdown}
[C1]

Consider regression through the origin, \( \E(Y)=\beta x \), and an M-estimate with fixed scale \( \hat{\sigma}=1 \) and a continuous, nondecreasing,
bounded, odd \( \psi \) with \( \psi(u)>0 \) for \( u>0 \). Add to fixed data \( (x_i,y_i) \), \( i<n \), one case \( (x_n,y_n)=(L,aL) \). Show that as \( L\to\infty \)
every solution \( \hat{\beta}_L \) of \( \sum_i\psi(y_i-\hat{\beta}x_i)x_i=0 \) tends to \( a \). Conclude that the breakdown point is \( 1/n \).
:::

::: {.solution}
Let \( B=\sup\lvert\psi\rvert \) and \( S=\sum_{i<n}\lvert x_i\rvert \), so the first \( n-1 \) terms contribute at most \( BS \) in absolute value, whatever
\( \hat{\beta} \). Fix \( \eta>0 \). If \( \hat{\beta}\le a-\eta \), the last term is \( \psi\bigl(L(a-\hat{\beta})\bigr)L\ge\psi(L\eta)L \), which exceeds \( BS \) once
\( L \) is large, because \( \psi(L\eta)\to\sup\psi>0 \) (as \( \psi \) is nondecreasing and positive on \( (0,\infty) \)). Then the sum is positive and
\( \hat{\beta} \) is not a solution. Symmetrically for \( \hat{\beta}\ge a+\eta \). So for large \( L \) every solution lies in \( (a-\eta,a+\eta) \). Since
\( a \) is arbitrary, one case can place the estimate anywhere, and the breakdown point is \( 1/n \).
:::
